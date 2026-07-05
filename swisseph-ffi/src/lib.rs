pub mod config;
pub mod error;

use std::ffi::c_char;

use swisseph::Ephemeris;
use swisseph::config::TopoPosition;
use swisseph::flags::CalcFlags;
use swisseph::types::Body;

use crate::config::SweConfig;
use crate::error::{SweErrorCode, error_code, ffi_guard, write_err};

/// Opaque handle wrapping an `Ephemeris`. Never `#[repr(C)]` — Dart/C sees only `*mut SweEphemeris`.
pub struct SweEphemeris(Ephemeris);

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
/// Ephemeris::new resolves the value from the file denum when unset — the caller
/// can't know it without asking.
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

/// Compute planetary position at `tjd_ut` (Julian Day, UT1).
///
/// # Parameters
/// - `handle`: ephemeris handle from `swisseph_new`
/// - `tjd_ut`: Julian Day in UT1
/// - `ipl`: body number (C `ipl` values, matching `Body` discriminant)
/// - `iflag`: calculation flags (C `SEFLG_*` bit values)
/// - `geopos`: NULL, or pointer to `[lon, lat, alt]` for a per-call topographic override
/// - `xx`: out-param, pointer to 6 `f64` slots receiving [lon, lat, dist, lon_speed, lat_speed, dist_speed]
/// - `flags_used`: out-param, pointer to `i32` receiving the flags actually applied
/// - `err_buf`, `err_cap`: optional error message buffer
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `xx` must point to at least 6 writable `f64` slots.
/// - `flags_used` may be NULL.
/// - `geopos`, if non-NULL, must point to 3 readable `f64` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_calc_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
    geopos: *const f64,
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

        let result = if !geopos.is_null() {
            let gp = unsafe { std::slice::from_raw_parts(geopos, 3) };
            let mut config = eph.config().clone();
            config.topographic = Some(TopoPosition {
                longitude: gp[0],
                latitude: gp[1],
                altitude: gp[2],
            });
            eph.calc_ut_with_config(tjd_ut, body, calc_flags, &config)
        } else {
            eph.calc_ut(tjd_ut, body, calc_flags)
        };

        match result {
            Ok(r) => {
                unsafe {
                    for i in 0..6 {
                        *xx.add(i) = r.data[i];
                    }
                    if !flags_used.is_null() {
                        *flags_used = r.flags_used.bits() as i32;
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
