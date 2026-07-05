use std::ffi::c_char;

use swisseph::flags::CalcFlags;
use swisseph::nodaps::NodApsMethod;
use swisseph::types::Body;

use crate::SweEphemeris;
use crate::error::{SweErrorCode, error_code, ffi_guard, write_err};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

unsafe fn write_phenomena(p: &swisseph::Phenomena, attr: *mut f64) {
    unsafe {
        *attr = p.phase_angle;
        *attr.add(1) = p.phase;
        *attr.add(2) = p.elongation;
        *attr.add(3) = p.apparent_diameter;
        *attr.add(4) = p.apparent_magnitude;
        *attr.add(5) = p.horizontal_parallax;
        for i in 6..20 {
            *attr.add(i) = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// swisseph_pheno — swe_pheno
// ---------------------------------------------------------------------------

/// Planetary phenomena (phase angle, elongation, magnitude, etc.) at `tjd_et` (TT).
///
/// # Parameters
/// - `ipl`: body number
/// - `iflag`: calculation flags
/// - `attr`: out-param, pointer to 20 `f64` slots. `attr[0]`=phase_angle, `[1]`=phase,
///   `[2]`=elongation, `[3]`=apparent_diameter, `[4]`=apparent_magnitude,
///   `[5]`=horizontal_parallax, `[6..19]`=0.
/// - `flags_used`: out-param (may be NULL), flags actually applied
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_pheno(
    handle: *const SweEphemeris,
    tjd_et: f64,
    ipl: i32,
    iflag: i32,
    attr: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || attr.is_null() {
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

        match eph.pheno(tjd_et, body, calc_flags) {
            Ok((pheno, used)) => {
                unsafe {
                    write_phenomena(&pheno, attr);
                    if !flags_used.is_null() {
                        *flags_used = used.bits() as i32;
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

/// Planetary phenomena at `tjd_ut` (UT1). See [`swisseph_pheno`] for details.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_pheno_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
    attr: *mut f64,
    flags_used: *mut i32,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || attr.is_null() {
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

        match eph.pheno_ut(tjd_ut, body, calc_flags) {
            Ok((pheno, used)) => {
                unsafe {
                    write_phenomena(&pheno, attr);
                    if !flags_used.is_null() {
                        *flags_used = used.bits() as i32;
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

// ---------------------------------------------------------------------------
// swisseph_nod_aps — swe_nod_aps
// ---------------------------------------------------------------------------

/// Nodes and apsides at `tjd_et` (TT).
///
/// # Parameters
/// - `ipl`: body number
/// - `iflag`: calculation flags
/// - `method`: `NodApsMethod` bits (1=MEAN, 2=OSCU, 4=OSCU_BAR, 256=FOPOINT)
/// - `xnasc`, `xndsc`, `xperi`, `xaphe`: out-params, each pointer to 6 `f64` slots
///   receiving ascending node, descending node, perihelion, aphelion respectively
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - All four output pointers must point to at least 6 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_nod_aps(
    handle: *const SweEphemeris,
    tjd_et: f64,
    ipl: i32,
    iflag: i32,
    method: i32,
    xnasc: *mut f64,
    xndsc: *mut f64,
    xperi: *mut f64,
    xaphe: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null()
            || xnasc.is_null()
            || xndsc.is_null()
            || xperi.is_null()
            || xaphe.is_null()
        {
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
        let nod_method = NodApsMethod::from_bits_retain(method as u32);

        match eph.nod_aps(tjd_et, body, calc_flags, nod_method) {
            Ok(na) => {
                unsafe {
                    for i in 0..6 {
                        *xnasc.add(i) = na.ascending[i];
                        *xndsc.add(i) = na.descending[i];
                        *xperi.add(i) = na.perihelion[i];
                        *xaphe.add(i) = na.aphelion[i];
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

/// Nodes and apsides at `tjd_ut` (UT1). See [`swisseph_nod_aps`] for details.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - All four output pointers must point to at least 6 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_nod_aps_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
    method: i32,
    xnasc: *mut f64,
    xndsc: *mut f64,
    xperi: *mut f64,
    xaphe: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null()
            || xnasc.is_null()
            || xndsc.is_null()
            || xperi.is_null()
            || xaphe.is_null()
        {
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
        let nod_method = NodApsMethod::from_bits_retain(method as u32);

        match eph.nod_aps_ut(tjd_ut, body, calc_flags, nod_method) {
            Ok(na) => {
                unsafe {
                    for i in 0..6 {
                        *xnasc.add(i) = na.ascending[i];
                        *xndsc.add(i) = na.descending[i];
                        *xperi.add(i) = na.perihelion[i];
                        *xaphe.add(i) = na.aphelion[i];
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

// ---------------------------------------------------------------------------
// swisseph_get_orbital_elements — swe_get_orbital_elements
// ---------------------------------------------------------------------------

/// Osculating orbital elements at `tjd_et` (TT).
///
/// # Parameters
/// - `ipl`: body number
/// - `iflag`: calculation flags
/// - `dret`: out-param, pointer to 50 `f64` slots. Slots `[0..16]` receive the 17 named
///   orbital element fields (see `OrbitalElements::as_array` for slot meanings):
///   `[0]`=semi_major_axis, `[1]`=eccentricity, `[2]`=inclination,
///   `[3]`=ascending_node (Ω), `[4]`=arg_perihelion (ω), `[5]`=perihelion_lon (ϖ),
///   `[6]`=mean_anomaly, `[7]`=true_anomaly, `[8]`=eccentric_anomaly,
///   `[9]`=mean_longitude, `[10]`=sidereal_period, `[11]`=mean_daily_motion,
///   `[12]`=tropical_period, `[13]`=synodic_period, `[14]`=perihelion_passage (JD TT),
///   `[15]`=perihelion_distance, `[16]`=aphelion_distance.
///   Slots `[17..49]` are zeroed.
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dret` must point to at least 50 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_get_orbital_elements(
    handle: *const SweEphemeris,
    tjd_et: f64,
    ipl: i32,
    iflag: i32,
    dret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || dret.is_null() {
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

        match eph.get_orbital_elements(tjd_et, body, calc_flags) {
            Ok(elems) => {
                let arr = elems.as_array();
                unsafe {
                    for i in 0..17 {
                        *dret.add(i) = arr[i];
                    }
                    for i in 17..50 {
                        *dret.add(i) = 0.0;
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

// ---------------------------------------------------------------------------
// swisseph_orbit_max_min_true_distance — swe_orbit_max_min_true_distance
// ---------------------------------------------------------------------------

/// Maximum, minimum, and current true distance for a body at `tjd_et` (TT).
///
/// # Parameters
/// - `dmax`, `dmin`, `dtrue`: out-params, pointers to writable `f64` values
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dmax`, `dmin`, `dtrue` must each point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_orbit_max_min_true_distance(
    handle: *const SweEphemeris,
    tjd_et: f64,
    ipl: i32,
    iflag: i32,
    dmax: *mut f64,
    dmin: *mut f64,
    dtrue: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || dmax.is_null() || dmin.is_null() || dtrue.is_null() {
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

        match eph.orbit_max_min_true_distance(tjd_et, body, calc_flags) {
            Ok((mx, mn, tr)) => {
                unsafe {
                    *dmax = mx;
                    *dmin = mn;
                    *dtrue = tr;
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
// Crossings
// ---------------------------------------------------------------------------
//
// IMPORTANT: C's swe_solcross etc. return the crossing JD as a double return
// value, with jd_cross < tjd_start meaning error. This FFI deliberately
// diverges: all crossings return i32 (0=OK, negative=error) with the crossing
// JD written to a `jx` out-param, matching our uniform error convention. The
// Dart binding must adapt accordingly.

/// Next Julian Day (TT) at which the Sun's ecliptic longitude equals `x2cross` (degrees).
///
/// **Return convention differs from C:** returns `i32` status (0=OK, negative=error),
/// crossing JD written to `*jx`. C returns the JD directly with error signaled by
/// `jd < tjd` — this FFI uses the uniform out-param convention instead.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_solcross(
    handle: *const SweEphemeris,
    x2cross: f64,
    tjd_et: f64,
    iflag: i32,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.solcross(x2cross, tjd_et, calc_flags) {
            Ok(jd) => {
                unsafe { *jx = jd };
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

/// UT-based [`swisseph_solcross`]. See that function for the return convention note.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_solcross_ut(
    handle: *const SweEphemeris,
    x2cross: f64,
    tjd_ut: f64,
    iflag: i32,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.solcross_ut(x2cross, tjd_ut, calc_flags) {
            Ok(jd) => {
                unsafe { *jx = jd };
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

/// Next Julian Day (TT) at which the Moon's ecliptic longitude equals `x2cross` (degrees).
///
/// **Return convention differs from C** — see [`swisseph_solcross`].
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_mooncross(
    handle: *const SweEphemeris,
    x2cross: f64,
    tjd_et: f64,
    iflag: i32,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.mooncross(x2cross, tjd_et, calc_flags) {
            Ok(jd) => {
                unsafe { *jx = jd };
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

/// UT-based [`swisseph_mooncross`].
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_mooncross_ut(
    handle: *const SweEphemeris,
    x2cross: f64,
    tjd_ut: f64,
    iflag: i32,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.mooncross_ut(x2cross, tjd_ut, calc_flags) {
            Ok(jd) => {
                unsafe { *jx = jd };
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

/// Next Julian Day (TT) at which the Moon crosses its node (ecliptic latitude = 0).
///
/// **Return convention differs from C** — see [`swisseph_solcross`].
///
/// # Parameters
/// - `jx`: out, crossing JD
/// - `xlon`: out, Moon's ecliptic longitude at the crossing (degrees)
/// - `xlat`: out, Moon's ecliptic latitude at the crossing (degrees, near zero)
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx`, `xlon`, `xlat` must each point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_mooncross_node(
    handle: *const SweEphemeris,
    tjd_et: f64,
    iflag: i32,
    xlon: *mut f64,
    xlat: *mut f64,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() || xlon.is_null() || xlat.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.mooncross_node(tjd_et, calc_flags) {
            Ok(mc) => {
                unsafe {
                    *jx = mc.jd;
                    *xlon = mc.longitude;
                    *xlat = mc.latitude;
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

/// UT-based [`swisseph_mooncross_node`].
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx`, `xlon`, `xlat` must each point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_mooncross_node_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    iflag: i32,
    xlon: *mut f64,
    xlat: *mut f64,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() || xlon.is_null() || xlat.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.mooncross_node_ut(tjd_ut, calc_flags) {
            Ok(mc) => {
                unsafe {
                    *jx = mc.jd;
                    *xlon = mc.longitude;
                    *xlat = mc.latitude;
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

/// Next Julian Day (TT) at which `ipl`'s heliocentric longitude equals `x2cross` (degrees).
/// `dir >= 0` searches forward, `dir < 0` searches backward.
///
/// **Return convention differs from C** — see [`swisseph_solcross`].
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx` must point to a writable `f64`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_helio_cross(
    handle: *const SweEphemeris,
    ipl: i32,
    x2cross: f64,
    tjd_et: f64,
    iflag: i32,
    dir: i32,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() {
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

        match eph.helio_cross(body, x2cross, tjd_et, calc_flags, dir) {
            Ok(jd) => {
                unsafe { *jx = jd };
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

/// UT-based [`swisseph_helio_cross`].
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `jx` must point to a writable `f64`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_helio_cross_ut(
    handle: *const SweEphemeris,
    ipl: i32,
    x2cross: f64,
    tjd_ut: f64,
    iflag: i32,
    dir: i32,
    jx: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || jx.is_null() {
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

        match eph.helio_cross_ut(body, x2cross, tjd_ut, calc_flags, dir) {
            Ok(jd) => {
                unsafe { *jx = jd };
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
