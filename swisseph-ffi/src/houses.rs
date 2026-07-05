use std::ffi::c_char;

use swisseph::azalt::{AzAltDir, HorDir, RefracDir};
use swisseph::flags::CalcFlags;
use swisseph::houses;
use swisseph::types::{Body, HouseSystem};

use crate::SweEphemeris;
use crate::cstr_to_str;
use crate::error::{SweErrorCode, error_code, ffi_guard, write_err};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn hsys_from_char(c: i32) -> Result<HouseSystem, i32> {
    HouseSystem::try_from(c as u8).map_err(|_| SweErrorCode::InvalidHouseSystem as i32)
}

unsafe fn write_ascmc(ascmc: &houses::AscMc, out: *mut f64) {
    if out.is_null() {
        return;
    }
    let arr = ascmc.as_array();
    unsafe {
        for i in 0..8 {
            *out.add(i) = arr[i];
        }
        // ascmc[8..9] zeroed per C convention
        *out.add(8) = 0.0;
        *out.add(9) = 0.0;
    }
}

unsafe fn write_cusps(result: &houses::HouseResult, hsys: HouseSystem, cusps: *mut f64) {
    if cusps.is_null() {
        return;
    }
    let n = if hsys == HouseSystem::Gauquelin {
        36
    } else {
        12
    };
    unsafe {
        for i in 0..=n {
            *cusps.add(i) = result.cusps[i];
        }
    }
}

unsafe fn write_cusp_speeds(result: &houses::HouseResult, hsys: HouseSystem, out: *mut f64) {
    if out.is_null() {
        return;
    }
    let n = if hsys == HouseSystem::Gauquelin {
        36
    } else {
        12
    };
    unsafe {
        for i in 0..=n {
            *out.add(i) = result.cusp_speeds[i];
        }
    }
}

// ---------------------------------------------------------------------------
// swisseph_houses — swe_houses
// ---------------------------------------------------------------------------

/// Compute tropical house cusps and angular points at `tjd_ut` (UT1).
///
/// # Parameters
/// - `handle`: ephemeris handle from `swisseph_new`
/// - `tjd_ut`: Julian Day in UT1
/// - `geolat`: geographic latitude (north positive), degrees
/// - `geolon`: geographic longitude (east positive), degrees
/// - `hsys`: house system as ASCII char code (e.g. `'P'` for Placidus)
/// - `cusps`: out, must point to at least 13 writable `f64` slots (37 for `'G'` Gauquelin)
/// - `ascmc`: out, must point to at least 10 writable `f64` slots
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `cusps` must point to at least 13 (or 37 for Gauquelin) writable `f64` slots.
/// - `ascmc` must point to at least 10 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_houses(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    geolat: f64,
    geolon: f64,
    hsys: i32,
    cusps: *mut f64,
    ascmc: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || cusps.is_null() || ascmc.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };

        match eph.houses(tjd_ut, geolat, geolon, hs) {
            Ok(r) => {
                unsafe {
                    write_cusps(&r, hs, cusps);
                    write_ascmc(&r.ascmc, ascmc);
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
// swisseph_houses_ex — swe_houses_ex
// ---------------------------------------------------------------------------

/// Compute house cusps and angular points at `tjd_ut` (UT1) with flags.
///
/// Sidereal mode (when `SEFLG_SIDEREAL` is set in `iflag`) uses the mode
/// configured on the handle at construction time.
///
/// # Parameters
/// - `iflag`: calculation flags (e.g. `SEFLG_SIDEREAL`)
/// - `cusps`: out, at least 13 `f64` slots (37 for Gauquelin)
/// - `ascmc`: out, at least 10 `f64` slots
///
/// ascmc layout: [0]=Asc, [1]=MC, [2]=ARMC, [3]=Vertex, [4]=equatorial Asc,
/// [5]=co-Asc (Koch), [6]=co-Asc (Munkasey), [7]=polar Asc, [8..9]=0.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `cusps` must point to at least 13 (or 37 for Gauquelin) writable `f64` slots.
/// - `ascmc` must point to at least 10 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_houses_ex(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    iflag: i32,
    geolat: f64,
    geolon: f64,
    hsys: i32,
    cusps: *mut f64,
    ascmc: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || cusps.is_null() || ascmc.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.houses_ex(tjd_ut, calc_flags, geolat, geolon, hs) {
            Ok(r) => {
                unsafe {
                    write_cusps(&r, hs, cusps);
                    write_ascmc(&r.ascmc, ascmc);
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
// swisseph_houses_ex2 — swe_houses_ex2 (with speeds)
// ---------------------------------------------------------------------------

/// Compute house cusps, angular points, and their speeds at `tjd_ut` (UT1).
///
/// Sidereal mode (when `SEFLG_SIDEREAL` is set in `iflag`) uses the mode
/// configured on the handle at construction time.
///
/// # Parameters
/// - `cusp_speed`: out, nullable — cusp speeds (degrees/day), same sizing as `cusps`
/// - `ascmc_speed`: out, nullable — angular-point speeds (10 slots, same layout as `ascmc`)
///
/// Cusps must be sized 13 for all systems except Gauquelin ('G') which needs 37.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `cusps`, `ascmc` must be valid and properly sized.
/// - `cusp_speed`, `ascmc_speed` may be NULL (speeds omitted).
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_houses_ex2(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    iflag: i32,
    geolat: f64,
    geolon: f64,
    hsys: i32,
    cusps: *mut f64,
    ascmc: *mut f64,
    cusp_speed: *mut f64,
    ascmc_speed: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || cusps.is_null() || ascmc.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.houses_ex2(tjd_ut, calc_flags, geolat, geolon, hs) {
            Ok(r) => {
                unsafe {
                    write_cusps(&r, hs, cusps);
                    write_ascmc(&r.ascmc, ascmc);
                    write_cusp_speeds(&r, hs, cusp_speed);
                    write_ascmc(&r.ascmc_speeds, ascmc_speed);
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
// swisseph_houses_armc — swe_houses_armc (free-fn, no handle)
// ---------------------------------------------------------------------------

/// Compute house cusps from ARMC, obliquity, and geographic latitude directly.
///
/// Handle-free function — does not need an Ephemeris instance.
///
/// # Safety
/// - `cusps` must point to at least 13 (or 37 for Gauquelin) writable `f64` slots.
/// - `ascmc` must point to at least 10 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_houses_armc(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: i32,
    cusps: *mut f64,
    ascmc: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if cusps.is_null() || ascmc.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };

        match houses::houses_armc(armc, geolat, eps, hs, None) {
            Ok(r) => {
                unsafe {
                    write_cusps(&r, hs, cusps);
                    write_ascmc(&r.ascmc, ascmc);
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
// swisseph_houses_armc_ex2 — swe_houses_armc_ex2 (with speeds)
// ---------------------------------------------------------------------------

/// Compute house cusps and speeds from ARMC, obliquity, and geographic latitude.
///
/// Handle-free. `sundec` is required for Sunshine house systems (`'I'`/`'i'`);
/// pass NULL for all others.
///
/// # Safety
/// - `cusps`, `ascmc` must be valid and properly sized.
/// - `cusp_speed`, `ascmc_speed` may be NULL.
/// - `sundec` may be NULL (required for Sunshine systems).
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_houses_armc_ex2(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: i32,
    sundec: *const f64,
    cusps: *mut f64,
    ascmc: *mut f64,
    cusp_speed: *mut f64,
    ascmc_speed: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if cusps.is_null() || ascmc.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };

        let sd = if sundec.is_null() {
            None
        } else {
            Some(unsafe { *sundec })
        };

        match houses::houses_armc(armc, geolat, eps, hs, sd) {
            Ok(r) => {
                unsafe {
                    write_cusps(&r, hs, cusps);
                    write_ascmc(&r.ascmc, ascmc);
                    write_cusp_speeds(&r, hs, cusp_speed);
                    write_ascmc(&r.ascmc_speeds, ascmc_speed);
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
// swisseph_house_pos — swe_house_pos
// ---------------------------------------------------------------------------

/// Compute the house position of a planet (continuous 1.0..13.0 for 12-house
/// systems, 1.0..37.0 for Gauquelin).
///
/// # Parameters
/// - `armc`: sidereal time as ARMC (degrees)
/// - `geolat`: geographic latitude (degrees)
/// - `eps`: obliquity of ecliptic (degrees)
/// - `hsys`: house system (ASCII char)
/// - `xpin`: ecliptic longitude and latitude of the planet, 2 `f64` values
/// - `sundec`: Sun declination, required for Sunshine (`'I'`/`'i'`); NULL otherwise
/// - `hpos`: out, house position (1.0..13.0 or 1.0..37.0)
///
/// # Safety
/// - `xpin` must point to 2 readable `f64` values.
/// - `hpos` must point to a writable `f64`.
/// - `sundec` may be NULL.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_house_pos(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: i32,
    xpin: *const f64,
    sundec: *const f64,
    hpos: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if xpin.is_null() || hpos.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };

        let xp = unsafe { [*xpin, *xpin.add(1)] };
        let sd = if sundec.is_null() {
            None
        } else {
            Some(unsafe { *sundec })
        };

        match houses::house_pos(armc, geolat, eps, hs, xp, sd) {
            Ok(pos) => {
                unsafe { *hpos = pos };
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
// swisseph_house_name — swe_house_name
// ---------------------------------------------------------------------------

/// Get the human-readable name of a house system.
///
/// Handle-free. Writes the name to `buf` (NUL-terminated, truncated if needed).
///
/// # Safety
/// - `buf` must point to at least `cap` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_house_name(
    hsys: i32,
    buf: *mut c_char,
    cap: usize,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if buf.is_null() || cap == 0 {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let hs = match hsys_from_char(hsys) {
            Ok(h) => h,
            Err(code) => {
                let msg = format!("invalid house system: '{}'", hsys as u8 as char);
                unsafe { write_err(err_buf, err_cap, &msg) };
                return code;
            }
        };

        let name = hs.name();
        unsafe { write_err(buf, cap, name) };
        SweErrorCode::Ok as i32
    })
}

// ---------------------------------------------------------------------------
// swisseph_gauquelin_sector — swe_gauquelin_sector
// ---------------------------------------------------------------------------

/// Compute a Gauquelin sector position at `tjd_ut` (UT1).
///
/// # Parameters
/// - `ipl`: planet body ID
/// - `starname`: fixed star name (NUL-terminated), or NULL for a planet
/// - `iflag`: calculation flags
/// - `imeth`: method (0/1 = geometric, 2–5 = rise/set-based)
/// - `geopos`: geographic position [lon, lat, height], 3 `f64` values
/// - `atpress`: atmospheric pressure (hPa), for rise/set methods
/// - `attemp`: atmospheric temperature (°C), for rise/set methods
/// - `dgsect`: out, Gauquelin sector (1.0..36.x)
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `geopos` must point to 3 readable `f64` values.
/// - `dgsect` must point to a writable `f64`.
/// - `starname` may be NULL.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_gauquelin_sector(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    starname: *const c_char,
    iflag: i32,
    imeth: i32,
    geopos: *const f64,
    atpress: f64,
    attemp: f64,
    dgsect: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || dgsect.is_null() {
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

        let star = if starname.is_null() {
            None
        } else {
            match unsafe { cstr_to_str(starname) } {
                Ok(s) if s.is_empty() => None,
                Ok(s) => Some(s),
                Err(msg) => {
                    unsafe { write_err(err_buf, err_cap, msg) };
                    return SweErrorCode::InvalidArg as i32;
                }
            }
        };

        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);

        match eph.gauquelin_sector(tjd_ut, body, star, calc_flags, imeth, gp, atpress, attemp) {
            Ok(sector) => {
                unsafe { *dgsect = sector };
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
// swisseph_azalt — swe_azalt
// ---------------------------------------------------------------------------

/// Convert ecliptic/equatorial coordinates to horizontal (azimuth + altitude).
///
/// # Parameters
/// - `calc_flag`: `SE_ECL2HOR` (0) or `SE_EQU2HOR` (1)
/// - `geopos`: [longitude, latitude, height], 3 `f64` values
/// - `atpress`: atmospheric pressure (hPa) — 0 for no refraction
/// - `attemp`: atmospheric temperature (°C)
/// - `xin`: input [longitude/RA, latitude/declination], 2 `f64` values
/// - `xaz`: out [azimuth, true altitude, apparent altitude], 3 `f64` values
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `xin` must point to 2 readable `f64` values (only [0] and [1] used).
/// - `xaz` must point to 3 writable `f64` values.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_azalt(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    calc_flag: i32,
    geopos: *const f64,
    atpress: f64,
    attemp: f64,
    xin: *const f64,
    xaz: *mut f64,
) {
    if handle.is_null() || geopos.is_null() || xin.is_null() || xaz.is_null() {
        return;
    }

    let eph = unsafe { &(*handle).0 };
    let dir = match calc_flag {
        0 => AzAltDir::EclToHor,
        _ => AzAltDir::EquToHor,
    };

    let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
    let xi = unsafe { [*xin, *xin.add(1)] };

    let result = eph.azalt(tjd_ut, dir, gp, atpress, attemp, 0.0, xi);
    unsafe {
        *xaz = result[0];
        *xaz.add(1) = result[1];
        *xaz.add(2) = result[2];
    }
}

// ---------------------------------------------------------------------------
// swisseph_azalt_rev — swe_azalt_rev
// ---------------------------------------------------------------------------

/// Convert horizontal (azimuth + altitude) back to ecliptic/equatorial.
///
/// # Parameters
/// - `calc_flag`: `SE_HOR2ECL` (0) or `SE_HOR2EQU` (1)
/// - `geopos`: [longitude, latitude, height], 3 `f64` values
/// - `xin`: [azimuth (from south, clockwise), true altitude], 2 `f64` values
/// - `xout`: out [lon/RA, lat/dec], 2 `f64` values
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `xin` must point to 2 readable `f64` values.
/// - `xout` must point to 2 writable `f64` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_azalt_rev(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    calc_flag: i32,
    geopos: *const f64,
    xin: *const f64,
    xout: *mut f64,
) {
    if handle.is_null() || geopos.is_null() || xin.is_null() || xout.is_null() {
        return;
    }

    let eph = unsafe { &(*handle).0 };
    let dir = match calc_flag {
        0 => HorDir::HorToEcl,
        _ => HorDir::HorToEqu,
    };

    let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
    let xi = unsafe { [*xin, *xin.add(1)] };

    let result = eph.azalt_rev(tjd_ut, dir, gp, xi);
    unsafe {
        *xout = result[0];
        *xout.add(1) = result[1];
    }
}

// ---------------------------------------------------------------------------
// swisseph_refrac — swe_refrac
// ---------------------------------------------------------------------------

/// Simple atmospheric refraction (sea-level, no dip).
///
/// Handle-free.
///
/// # Parameters
/// - `inalt`: input altitude (degrees)
/// - `atpress`: atmospheric pressure (hPa)
/// - `attemp`: atmospheric temperature (°C)
/// - `calc_flag`: `SE_TRUE_TO_APP` (0) or `SE_APP_TO_TRUE` (1)
///
/// Returns the refracted/de-refracted altitude (degrees).
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_refrac(inalt: f64, atpress: f64, attemp: f64, calc_flag: i32) -> f64 {
    let dir = match calc_flag {
        0 => RefracDir::TrueToApp,
        _ => RefracDir::AppToTrue,
    };
    swisseph::azalt::refrac(inalt, atpress, attemp, dir)
}

// ---------------------------------------------------------------------------
// swisseph_refrac_extended — swe_refrac_extended
// ---------------------------------------------------------------------------

/// Extended atmospheric refraction with horizon dip for an elevated observer.
///
/// Handle-free.
///
/// # Parameters
/// - `inalt`: input altitude (degrees)
/// - `geoalt`: observer height above sea level (meters)
/// - `atpress`: atmospheric pressure (hPa)
/// - `attemp`: atmospheric temperature (°C)
/// - `lapse_rate`: temperature lapse rate (K/m), typically 0.0065
/// - `calc_flag`: `SE_TRUE_TO_APP` (0) or `SE_APP_TO_TRUE` (1)
/// - `dret`: out, 4 `f64` values: [true alt, apparent alt, refraction, dip]
///
/// Returns the refracted/de-refracted altitude.
///
/// # Safety
/// - `dret` must point to at least 4 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_refrac_extended(
    inalt: f64,
    geoalt: f64,
    atpress: f64,
    attemp: f64,
    lapse_rate: f64,
    calc_flag: i32,
    dret: *mut f64,
) -> f64 {
    let dir = match calc_flag {
        0 => RefracDir::TrueToApp,
        _ => RefracDir::AppToTrue,
    };

    if dret.is_null() {
        let mut d = [0.0f64; 4];
        return swisseph::azalt::refrac_extended(
            inalt, geoalt, atpress, attemp, lapse_rate, dir, &mut d,
        );
    }

    let mut d = [0.0f64; 4];
    let result =
        swisseph::azalt::refrac_extended(inalt, geoalt, atpress, attemp, lapse_rate, dir, &mut d);
    unsafe {
        for i in 0..4 {
            *dret.add(i) = d[i];
        }
    }
    result
}
