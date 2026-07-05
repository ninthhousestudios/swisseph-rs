pub mod config;
pub mod date;
pub mod eclipse;
pub mod error;
pub mod heliacal;
pub mod houses;
pub mod pheno;
pub mod util;

use std::ffi::c_char;

use swisseph::Ephemeris;
use swisseph::config::{EphemerisConfig, TopoPosition};
use swisseph::flags::CalcFlags;
use swisseph::types::{Body, SiderealMode};

use crate::config::SweConfig;
use crate::error::{SweErrorCode, error_code, ffi_guard, write_err};

/// Opaque handle wrapping an `Ephemeris`. Never `#[repr(C)]` — Dart/C sees only `*mut SweEphemeris`.
pub struct SweEphemeris(Ephemeris);

/// Per-call sidereal mode override. Nullable in all FFI signatures — NULL means
/// "use the handle's configured sidereal mode".
#[repr(C)]
pub struct SweSidMode {
    /// Raw `swe_set_sid_mode` value (bits 0-7 = mode index, upper bits = SiderealBits).
    pub sid_mode: i32,
    /// Reference epoch for user-defined sidereal (mode index 255).
    pub t0: f64,
    /// Initial ayanamsa at `t0`.
    pub ayan_t0: f64,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a config with per-call overrides applied. Returns `None` when no overrides
/// are present (callers should use the plain Ephemeris method in that case).
pub(crate) unsafe fn build_config(
    eph: &Ephemeris,
    geopos: *const f64,
    sid_mode: *const SweSidMode,
) -> Option<EphemerisConfig> {
    if geopos.is_null() && sid_mode.is_null() {
        return None;
    }
    let mut config = eph.config().clone();
    if !geopos.is_null() {
        let gp = unsafe { std::slice::from_raw_parts(geopos, 3) };
        config.topographic = Some(TopoPosition {
            longitude: gp[0],
            latitude: gp[1],
            altitude: gp[2],
        });
    }
    if !sid_mode.is_null() {
        let sm = unsafe { &*sid_mode };
        config.set_sidereal_mode(sm.sid_mode, sm.t0, sm.ayan_t0);
    }
    Some(config)
}

/// Build a config with only a sidereal override (no geopos). Returns `None` when
/// `sid_mode` is NULL.
unsafe fn build_sid_config(
    eph: &Ephemeris,
    sid_mode: *const SweSidMode,
) -> Option<EphemerisConfig> {
    unsafe { build_config(eph, std::ptr::null(), sid_mode) }
}

/// Write a CalcResult's data and flags_used to FFI out-params.
unsafe fn write_calc_result(r: &swisseph::CalcResult, xx: *mut f64, flags_used: *mut i32) {
    unsafe {
        for i in 0..6 {
            *xx.add(i) = r.data[i];
        }
        if !flags_used.is_null() {
            *flags_used = r.flags_used.bits() as i32;
        }
    }
}

/// Convert a `*const c_char` to a `&str`. Returns `Err` with a static message on
/// null or invalid UTF-8.
pub(crate) unsafe fn cstr_to_str<'a>(ptr: *const c_char) -> Result<&'a str, &'static str> {
    if ptr.is_null() {
        return Err("null star name");
    }
    let cstr = unsafe { std::ffi::CStr::from_ptr(ptr) };
    cstr.to_str().map_err(|_| "invalid UTF-8 in star name")
}

// ---------------------------------------------------------------------------
// Handle lifecycle + introspection
// ---------------------------------------------------------------------------

/// Return the library version as a static NUL-terminated UTF-8 string.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_version() -> *const c_char {
    static VERSION: &[u8] = b"0.1.0\0";
    VERSION.as_ptr() as *const c_char
}

/// Create a new ephemeris handle from a flattened config.
///
/// On success, writes the handle to `*out` and returns 0.
/// On failure, returns a negative error code and writes a message to `err_buf`.
///
/// # Safety
/// - `config` must point to a valid `SweConfig` (use `swisseph_config_default` to initialize).
/// - `out` must point to a writable `*mut SweEphemeris`.
/// - `err_buf` may be NULL; if non-NULL, `err_cap` bytes must be writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_new(
    config: *const SweConfig,
    out: *mut *mut SweEphemeris,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if config.is_null() || out.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let rust_config = match unsafe { crate::config::config_to_rust(&*config) } {
            Ok(c) => c,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        match Ephemeris::new(rust_config) {
            Ok(eph) => {
                let boxed = Box::new(SweEphemeris(eph));
                unsafe { *out = Box::into_raw(boxed) };
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

/// Free an ephemeris handle. Null-safe (no-op on NULL).
///
/// # Safety
/// `handle` must be NULL or a pointer previously returned by `swisseph_new` that has not
/// already been freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_free(handle: *mut SweEphemeris) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Return the resolved tidal acceleration (arcsec/century^2).
/// NAN when unresolved (e.g. Moshier without explicit override).
///
/// # Safety
/// `handle` must be a valid, non-NULL handle from `swisseph_new`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_tid_acc(handle: *const SweEphemeris) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    let eph = unsafe { &(*handle).0 };
    eph.config().tidal_acceleration.unwrap_or(f64::NAN)
}

/// Write the resolved astro model values into `out` (array of 8 i32s).
/// Order: prec_longterm, prec_shortterm, nutation, bias, jplhor, jplhora, sidereal_time, delta_t.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle from `swisseph_new`.
/// - `out` must point to at least 8 writable `i32` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_astro_models(
    handle: *const SweEphemeris,
    out: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || out.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }
        let eph = unsafe { &(*handle).0 };
        let m = &eph.config().astro_models;
        let vals = astro_models_to_i32s(m);
        unsafe {
            for (i, v) in vals.iter().enumerate() {
                *out.add(i) = *v;
            }
        }
        SweErrorCode::Ok as i32
    })
}

fn astro_models_to_i32s(m: &swisseph::AstroModels) -> [i32; 8] {
    [
        m.prec_longterm as i32,
        m.prec_shortterm as i32,
        m.nutation as i32,
        m.bias as i32,
        m.jplhor_mode as i32,
        m.jplhora_mode as i32,
        m.sidereal_time as i32,
        m.delta_t as i32,
    ]
}

// ---------------------------------------------------------------------------
// File data introspection (stateless swe_get_current_file_data)
// ---------------------------------------------------------------------------

/// Query which ephemeris file would serve a calculation at `jd` for the given
/// file category `ifno`.
///
/// This is the stateless equivalent of C's `swe_get_current_file_data(ifno)`.
/// Instead of reporting the file used by the last `swe_calc` call, the caller
/// provides `jd` to select the file explicitly.
///
/// `ifno` values (mirrors C):
/// - 0 = planet (`sepl*.se1` or JPL `.eph`)
/// - 1 = moon (`semo*.se1`)
/// - 2 = main asteroid (`seas*.se1`)
/// - 3 = individual asteroid (always returns "no data" — stateless, no "last used")
/// - 4 = planet moon (always returns "no data" — stateless, no "last used")
///
/// Returns 0 on success, negative error code on failure. Returns
/// `EphemerisNotAvailable` when no file covers the given `jd` (including
/// Moshier-only configs which have no files).
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle from `swisseph_new`.
/// - `path_buf`, if non-NULL, must point to at least `path_cap` writable bytes.
/// - `tfstart`, `tfend`, `denum` must each point to a writable slot (or be NULL).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_file_data(
    handle: *const SweEphemeris,
    ifno: i32,
    jd: f64,
    path_buf: *mut c_char,
    path_cap: usize,
    tfstart: *mut f64,
    tfend: *mut f64,
    denum: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let kind = match swisseph::FileDataKind::try_from(ifno) {
            Ok(k) => k,
            Err(_) => {
                let msg = format!("invalid ifno: {ifno}");
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        match eph.file_data(kind, jd) {
            Some(fd) => unsafe {
                let path_str = fd.path.to_string_lossy();
                write_err(path_buf, path_cap, &path_str);
                if !tfstart.is_null() {
                    *tfstart = fd.start_jd;
                }
                if !tfend.is_null() {
                    *tfend = fd.end_jd;
                }
                if !denum.is_null() {
                    *denum = fd.denum;
                }
                SweErrorCode::Ok as i32
            },
            None => {
                unsafe {
                    write_err(
                        err_buf,
                        err_cap,
                        "no file data available for the given ifno/jd",
                    )
                };
                SweErrorCode::EphemerisNotAvailable as i32
            }
        }
    })
}

// ---------------------------------------------------------------------------
// calc / calc_ut
// ---------------------------------------------------------------------------

/// Compute planetary position at `tjd_ut` (Julian Day, UT1).
///
/// # Parameters
/// - `handle`: ephemeris handle from `swisseph_new`
/// - `tjd_ut`: Julian Day in UT1
/// - `ipl`: body number (C `ipl` values, matching `Body` discriminant)
/// - `iflag`: calculation flags (C `SEFLG_*` bit values)
/// - `geopos`: NULL, or pointer to `[lon, lat, alt]` for a per-call topographic override
/// - `sid_mode`: NULL, or pointer to a `SweSidMode` for a per-call sidereal override
/// - `xx`: out-param, pointer to 6 `f64` slots receiving [lon, lat, dist, lon_speed, lat_speed, dist_speed]
/// - `flags_used`: out-param (may be NULL), pointer to `i32` receiving the flags actually applied
/// - `err_buf`, `err_cap`: optional error message buffer
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `xx` must point to at least 6 writable `f64` slots.
/// - `geopos`, if non-NULL, must point to 3 readable `f64` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_calc_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
    geopos: *const f64,
    sid_mode: *const SweSidMode,
    xx: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || xx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let body = match Body::try_from(ipl) {
            Ok(b) => b,
            Err(_) => {
                let msg = format!("invalid body ID: {ipl}");
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidBody as i32;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        let result = match unsafe { build_config(eph, geopos, sid_mode) } {
            Some(config) => eph.calc_ut_with_config(tjd_ut, body, calc_flags, &config),
            None => eph.calc_ut(tjd_ut, body, calc_flags),
        };

        match result {
            Ok(r) => {
                unsafe { write_calc_result(&r, xx, flags_used) };
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

/// Compute planetary position at `tjd_et` (Julian Day, TT/ET).
///
/// # Parameters
/// - `handle`: ephemeris handle from `swisseph_new`
/// - `tjd_et`: Julian Day in TT/ET
/// - `ipl`: body number (C `ipl` values, matching `Body` discriminant)
/// - `iflag`: calculation flags (C `SEFLG_*` bit values)
/// - `geopos`: NULL, or pointer to `[lon, lat, alt]` for a per-call topographic override
/// - `sid_mode`: NULL, or pointer to a `SweSidMode` for a per-call sidereal override
/// - `xx`: out-param, pointer to 6 `f64` slots receiving [lon, lat, dist, lon_speed, lat_speed, dist_speed]
/// - `flags_used`: out-param (may be NULL), pointer to `i32` receiving the flags actually applied
/// - `err_buf`, `err_cap`: optional error message buffer
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `xx` must point to at least 6 writable `f64` slots.
/// - `geopos`, if non-NULL, must point to 3 readable `f64` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_calc(
    handle: *const SweEphemeris,
    tjd_et: f64,
    ipl: i32,
    iflag: i32,
    geopos: *const f64,
    sid_mode: *const SweSidMode,
    xx: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || xx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let body = match Body::try_from(ipl) {
            Ok(b) => b,
            Err(_) => {
                let msg = format!("invalid body ID: {ipl}");
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidBody as i32;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        let result = match unsafe { build_config(eph, geopos, sid_mode) } {
            Some(config) => eph.calc_with_config(tjd_et, body, calc_flags, &config),
            None => eph.calc(tjd_et, body, calc_flags),
        };

        match result {
            Ok(r) => {
                unsafe { write_calc_result(&r, xx, flags_used) };
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// calc_pctr
// ---------------------------------------------------------------------------

/// Compute planetocentric position at `tjd_et` (TT).
/// Swiss/JPL only — Moshier returns an error.
///
/// # Parameters
/// - `ipl`: target body
/// - `iplctr`: center body (must differ from `ipl`)
/// - No `geopos` override (C's `swe_calc_pctr` has no topocentric path)
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `xx` must point to at least 6 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_calc_pctr(
    handle: *const SweEphemeris,
    tjd_et: f64,
    ipl: i32,
    iplctr: i32,
    iflag: i32,
    xx: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || xx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let body = match Body::try_from(ipl) {
            Ok(b) => b,
            Err(_) => {
                let msg = format!("invalid body ID: {ipl}");
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidBody as i32;
            }
        };
        let center = match Body::try_from(iplctr) {
            Ok(b) => b,
            Err(_) => {
                let msg = format!("invalid center body ID: {iplctr}");
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidBody as i32;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.calc_pctr(tjd_et, body, center, calc_flags) {
            Ok(r) => {
                unsafe { write_calc_result(&r, xx, flags_used) };
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// fixstar2 family
// ---------------------------------------------------------------------------

/// Compute fixed-star position at `tjd_et` (TT).
///
/// # Parameters
/// - `star`: input star name (NUL-terminated UTF-8)
/// - `star_out`: buffer receiving the resolved "name,bayer" canonical name (NUL-terminated)
/// - `star_out_cap`: capacity of `star_out` in bytes (including NUL)
/// - `geopos`: NULL, or `[lon, lat, alt]` for per-call topographic override
/// - `sid_mode`: NULL, or per-call sidereal override
/// - `xx`: out-param, 6 `f64` slots
/// - `flags_used`: out-param (may be NULL)
///
/// # Safety
/// - `star` must be a valid NUL-terminated UTF-8 string.
/// - `star_out` may be NULL; if non-NULL, `star_out_cap` bytes must be writable.
/// - `xx` must point to at least 6 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_fixstar2(
    handle: *const SweEphemeris,
    star: *const c_char,
    star_out: *mut c_char,
    star_out_cap: usize,
    tjd_et: f64,
    iflag: i32,
    geopos: *const f64,
    sid_mode: *const SweSidMode,
    xx: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || xx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let star_str = match unsafe { cstr_to_str(star) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        let result = match unsafe { build_config(eph, geopos, sid_mode) } {
            Some(config) => eph.fixstar2_with_config(star_str, tjd_et, calc_flags, &config),
            None => eph.fixstar2(star_str, tjd_et, calc_flags),
        };

        match result {
            Ok((name, r)) => {
                unsafe {
                    write_err(star_out, star_out_cap, &name);
                    write_calc_result(&r, xx, flags_used);
                }
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

/// Compute fixed-star position at `tjd_ut` (UT1).
///
/// # Parameters
/// - `star`: input star name (NUL-terminated UTF-8)
/// - `star_out`: buffer receiving the resolved "name,bayer" canonical name (NUL-terminated);
///   may be NULL if the resolved name is not needed
/// - `star_out_cap`: capacity of `star_out` in bytes (including NUL)
/// - `geopos`: NULL, or `[lon, lat, alt]` for per-call topographic override
/// - `sid_mode`: NULL, or per-call sidereal override
/// - `xx`: out-param, pointer to 6 `f64` slots
/// - `flags_used`: out-param (may be NULL)
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `star` must be a valid NUL-terminated UTF-8 string.
/// - `star_out` may be NULL; if non-NULL, `star_out_cap` bytes must be writable.
/// - `xx` must point to at least 6 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_fixstar2_ut(
    handle: *const SweEphemeris,
    star: *const c_char,
    star_out: *mut c_char,
    star_out_cap: usize,
    tjd_ut: f64,
    iflag: i32,
    geopos: *const f64,
    sid_mode: *const SweSidMode,
    xx: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || xx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let star_str = match unsafe { cstr_to_str(star) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        // fixstar2_ut computes deltaT internally, but we need the _with_config
        // path for per-call overrides. When overrides are present, manually
        // compute deltaT and route through fixstar2_with_config (TT).
        let result = match unsafe { build_config(eph, geopos, sid_mode) } {
            Some(config) => {
                let dt = swisseph::deltat::calc_deltat(tjd_ut, &config);
                eph.fixstar2_with_config(star_str, tjd_ut + dt, calc_flags, &config)
            }
            None => eph.fixstar2_ut(star_str, tjd_ut, calc_flags),
        };

        match result {
            Ok((name, r)) => {
                unsafe {
                    write_err(star_out, star_out_cap, &name);
                    write_calc_result(&r, xx, flags_used);
                }
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

/// Look up the magnitude of a fixed star by name.
///
/// # Parameters
/// - `star`: input star name (NUL-terminated UTF-8)
/// - `star_out`: buffer receiving the resolved canonical name (may be NULL)
/// - `star_out_cap`: capacity of `star_out`
/// - `mag`: out-param receiving the magnitude
///
/// # Safety
/// - `star` must be a valid NUL-terminated UTF-8 string.
/// - `mag` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_fixstar2_mag(
    handle: *const SweEphemeris,
    star: *const c_char,
    star_out: *mut c_char,
    star_out_cap: usize,
    mag: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || mag.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let star_str = match unsafe { cstr_to_str(star) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        match eph.fixstar2_mag(star_str) {
            Ok((name, m)) => {
                unsafe {
                    write_err(star_out, star_out_cap, &name);
                    *mag = m;
                }
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Ayanamsa
// ---------------------------------------------------------------------------

/// Ayanamsa at `tjd_et` (TT) with flags. Nutation added unless `NONUT` is set.
///
/// # Parameters
/// - `sid_mode`: NULL uses the handle's configured sidereal mode; non-NULL overrides per-call
/// - `daya`: out-param receiving the ayanamsa in degrees
/// - `flags_used`: out-param (may be NULL)
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `daya` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_ayanamsa_ex(
    handle: *const SweEphemeris,
    tjd_et: f64,
    iflag: i32,
    sid_mode: *const SweSidMode,
    daya: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || daya.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        let result = match unsafe { build_sid_config(eph, sid_mode) } {
            Some(config) => eph.get_ayanamsa_ex_with_config(tjd_et, calc_flags, &config),
            None => eph.get_ayanamsa_ex(tjd_et, calc_flags),
        };

        match result {
            Ok(val) => {
                unsafe {
                    *daya = val;
                    if !flags_used.is_null() {
                        *flags_used = iflag;
                    }
                }
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

/// Ayanamsa at `tjd_ut` (UT1) with flags. Nutation added unless `NONUT` is set.
///
/// # Parameters
/// - `sid_mode`: NULL uses the handle's configured sidereal mode; non-NULL overrides per-call
/// - `daya`: out-param receiving the ayanamsa in degrees
/// - `flags_used`: out-param (may be NULL)
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `daya` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_ayanamsa_ex_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    iflag: i32,
    sid_mode: *const SweSidMode,
    daya: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || daya.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        let result = match unsafe { build_sid_config(eph, sid_mode) } {
            Some(config) => {
                let dt = swisseph::deltat::calc_deltat(tjd_ut, &config);
                eph.get_ayanamsa_ex_with_config(tjd_ut + dt, calc_flags, &config)
            }
            None => eph.get_ayanamsa_ut(tjd_ut, calc_flags),
        };

        match result {
            Ok(val) => {
                unsafe {
                    *daya = val;
                    if !flags_used.is_null() {
                        *flags_used = iflag;
                    }
                }
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

/// Legacy ayanamsa at `tjd_et` (TT), no nutation, returns degrees directly.
/// Returns NAN on error (e.g. no sidereal mode configured and no per-call override).
///
/// # Safety
/// `handle` must be valid, non-NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_ayanamsa(
    handle: *const SweEphemeris,
    tjd_et: f64,
    sid_mode: *const SweSidMode,
) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    let eph = unsafe { &(*handle).0 };

    let result = match unsafe { build_sid_config(eph, sid_mode) } {
        Some(config) => swisseph::ayanamsa::get_ayanamsa_ex(
            &config,
            tjd_et,
            CalcFlags::empty(),
            &config.astro_models,
        ),
        None => eph.get_ayanamsa(tjd_et),
    };

    result.unwrap_or(f64::NAN)
}

/// Legacy ayanamsa at `tjd_ut` (UT1), no nutation, returns degrees directly.
/// Returns NAN on error.
///
/// # Safety
/// `handle` must be valid, non-NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_ayanamsa_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    sid_mode: *const SweSidMode,
) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    let eph = unsafe { &(*handle).0 };

    let config = match unsafe { build_sid_config(eph, sid_mode) } {
        Some(c) => c,
        None => eph.config().clone(),
    };
    let dt = swisseph::deltat::calc_deltat(tjd_ut, &config);
    swisseph::ayanamsa::get_ayanamsa_ex(
        &config,
        tjd_ut + dt,
        CalcFlags::empty(),
        &config.astro_models,
    )
    .unwrap_or(f64::NAN)
}

// ---------------------------------------------------------------------------
// Ayanamsa name (handle-free)
// ---------------------------------------------------------------------------

/// Write the human-readable name for a sidereal mode into `buf`.
/// User-defined mode (255) writes an empty string. Unknown mode returns an error.
///
/// Handle-free — does not require an ephemeris instance.
///
/// # Safety
/// - `buf` must point to at least `cap` writable bytes, or be NULL (returns error).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_ayanamsa_name(
    sid_mode_raw: i32,
    buf: *mut c_char,
    cap: usize,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if buf.is_null() || cap == 0 {
            unsafe { write_err(err_buf, err_cap, "null or zero-capacity buffer") };
            return SweErrorCode::InvalidArg as i32;
        }

        let mode_index = sid_mode_raw & 0xFF;
        match SiderealMode::try_from(mode_index) {
            Ok(mode) => {
                let name = mode.name().unwrap_or("");
                unsafe { write_err(buf, cap, name) };
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                SweErrorCode::InvalidSiderealMode as i32
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Planet name
// ---------------------------------------------------------------------------

/// Write the display name for a body into `buf` (e.g. "Sun", "Chiron", asteroid name).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `buf` must point to at least `cap` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_planet_name(
    handle: *const SweEphemeris,
    ipl: i32,
    buf: *mut c_char,
    cap: usize,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || buf.is_null() || cap == 0 {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let body = match Body::try_from(ipl) {
            Ok(b) => b,
            Err(_) => {
                let msg = format!("invalid body ID: {ipl}");
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidBody as i32;
            }
        };

        let name = eph.get_planet_name(body);
        unsafe { write_err(buf, cap, &name) };
        SweErrorCode::Ok as i32
    })
}
