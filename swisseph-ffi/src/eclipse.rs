use std::ffi::c_char;

use swisseph::flags::{CalcFlags, EclipseFlags, RiseSetFlags};
use swisseph::types::Body;

use crate::SweEphemeris;
use crate::cstr_to_str;
use crate::error::{SweErrorCode, error_code, ffi_guard, write_err};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

unsafe fn parse_starname_opt<'a>(starname: *const c_char) -> Result<Option<&'a str>, &'static str> {
    if starname.is_null() {
        return Ok(None);
    }
    match unsafe { cstr_to_str(starname) } {
        Ok(s) if s.is_empty() => Ok(None),
        Ok(s) => Ok(Some(s)),
        Err(msg) => Err(msg),
    }
}

unsafe fn write_eclipse_how_attr(how: &swisseph::eclipse::EclipseHow, attr: *mut f64) {
    if attr.is_null() {
        return;
    }
    unsafe {
        *attr = how.magnitude;
        *attr.add(1) = how.diameter_ratio;
        *attr.add(2) = how.obscuration;
        *attr.add(3) = how.core_diameter_km;
        *attr.add(4) = how.azimuth;
        *attr.add(5) = how.true_altitude;
        *attr.add(6) = how.apparent_altitude;
        *attr.add(7) = how.elongation;
        *attr.add(8) = how.nasa_magnitude;
        *attr.add(9) = how.saros_series;
        *attr.add(10) = how.saros_member;
        for i in 11..20 {
            *attr.add(i) = 0.0;
        }
    }
}

unsafe fn write_lun_eclipse_how_attr(how: &swisseph::eclipse::LunarEclipseHow, attr: *mut f64) {
    if attr.is_null() {
        return;
    }
    unsafe {
        *attr = how.umbral_magnitude;
        *attr.add(1) = how.penumbral_magnitude;
        *attr.add(2) = 0.0;
        *attr.add(3) = 0.0;
        *attr.add(4) = how.azimuth;
        *attr.add(5) = how.true_altitude;
        *attr.add(6) = how.apparent_altitude;
        *attr.add(7) = how.distance_from_opposition;
        *attr.add(8) = how.umbral_magnitude;
        *attr.add(9) = how.saros_series;
        *attr.add(10) = how.saros_member;
        for i in 11..20 {
            *attr.add(i) = 0.0;
        }
    }
}

unsafe fn write_eclipse_where_geopos(w: &swisseph::eclipse::EclipseWhere, geopos: *mut f64) {
    if geopos.is_null() {
        return;
    }
    unsafe {
        *geopos = w.central_longitude;
        *geopos.add(1) = w.central_latitude;
        *geopos.add(2) = w.core_diameter_km;
        *geopos.add(3) = w.penumbra_diameter_km;
        *geopos.add(4) = w.shadow_axis_distance_km;
        *geopos.add(5) = w.umbra_diameter_fundamental_km;
        *geopos.add(6) = w.penumbra_diameter_fundamental_km;
        *geopos.add(7) = w.cos_umbra_half_angle;
        *geopos.add(8) = w.cos_penumbra_half_angle;
        *geopos.add(9) = 0.0;
    }
}

unsafe fn zero_f64_array(ptr: *mut f64, n: usize) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        for i in 0..n {
            *ptr.add(i) = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// swisseph_rise_trans — swe_rise_trans
// ---------------------------------------------------------------------------

/// Rise/set/meridian-transit search at `tjd_ut` (UT1).
///
/// Returns 0 on success, **-2** for circumpolar body (no error message — this is
/// a status, not an error, matching C convention), or a negative error code with
/// a message in `err_buf`.
///
/// # Parameters
/// - `ipl`: body number
/// - `starname`: fixed star name (NUL-terminated), or NULL for a planet
/// - `epheflag`: calculation flags (SEFLG_MOSEPH etc.)
/// - `rsmi`: rise/set event selector (SE_CALC_RISE=1, SE_CALC_SET=2, etc.)
/// - `geopos`: [lon, lat, height], 3 `f64` values
/// - `atpress`: atmospheric pressure (hPa)
/// - `attemp`: atmospheric temperature (°C)
/// - `tret`: out, pointer to 1 `f64` receiving the event time (UT)
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `tret` must point to a writable `f64`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_rise_trans(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    starname: *const c_char,
    epheflag: i32,
    rsmi: i32,
    geopos: *const f64,
    atpress: f64,
    attemp: f64,
    tret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || tret.is_null() {
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

        let star = match unsafe { parse_starname_opt(starname) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(epheflag as u32);
        let rsmi_flags = RiseSetFlags::from_bits_retain(rsmi as u32);

        match eph.rise_trans(
            tjd_ut, body, star, calc_flags, rsmi_flags, gp, atpress, attemp,
        ) {
            Ok(r) => {
                unsafe { *tret = r.time };
                SweErrorCode::Ok as i32
            }
            Err(swisseph::Error::CircumpolarBody) => -2,
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// swisseph_rise_trans_true_hor — swe_rise_trans_true_hor
// ---------------------------------------------------------------------------

/// Rise/set/meridian-transit search with custom horizon height at `tjd_ut` (UT1).
///
/// Same return convention as [`swisseph_rise_trans`]: 0 = success, -2 = circumpolar,
/// negative = error.
///
/// # Parameters
/// - `horhgt`: horizon height in degrees (-100 = auto-dip from `geopos[2]`)
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `tret` must point to a writable `f64`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_rise_trans_true_hor(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    starname: *const c_char,
    epheflag: i32,
    rsmi: i32,
    geopos: *const f64,
    atpress: f64,
    attemp: f64,
    horhgt: f64,
    tret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || tret.is_null() {
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

        let star = match unsafe { parse_starname_opt(starname) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(epheflag as u32);
        let rsmi_flags = RiseSetFlags::from_bits_retain(rsmi as u32);

        match eph.rise_trans_true_hor(
            tjd_ut, body, star, calc_flags, rsmi_flags, gp, atpress, attemp, horhgt,
        ) {
            Ok(r) => {
                unsafe { *tret = r.time };
                SweErrorCode::Ok as i32
            }
            Err(swisseph::Error::CircumpolarBody) => -2,
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// swisseph_sol_eclipse_where — swe_sol_eclipse_where
// ---------------------------------------------------------------------------

/// Geographic position of greatest solar eclipse at `tjd_ut` (UT1).
///
/// On success returns **positive** EclipseFlags bits (CENTRAL/NONCENTRAL/TOTAL/ANNULAR/PARTIAL);
/// 0 (empty flags) means no eclipse at this instant. Negative = error.
///
/// `geopos[10]` out: [0]=lon, [1]=lat, [2..8]=shadow-cone geometry (core_diameter_km,
/// penumbra_diameter_km, shadow_axis_distance_km, umbra_fundamental_km, penumbra_fundamental_km,
/// cos_umbra_half, cos_penumbra_half), [9]=0.
///
/// `attr[20]` out: local circumstances at the central point (magnitude, azimuth, etc.),
/// populated via an internal `eclipse_how` call. `attr[3]` is `dcore[0]` (core shadow diameter).
/// Zeroed when no eclipse is found.
///
/// # Safety
/// - `handle`, `geopos`, `attr` must be valid, non-NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_sol_eclipse_where(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ifl: i32,
    geopos: *mut f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || attr.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        match eph.sol_eclipse_where(tjd_ut, calc_flags) {
            Ok(w) => {
                unsafe {
                    write_eclipse_where_geopos(&w, geopos);
                }
                let ifl_masked =
                    calc_flags & (CalcFlags::JPLEPH | CalcFlags::SWIEPH | CalcFlags::MOSEPH);
                if !w.flags.is_empty() {
                    if let Ok(how) = eph.eclipse_how_at(
                        tjd_ut,
                        Body::Sun,
                        None,
                        ifl_masked,
                        [w.central_longitude, w.central_latitude, 0.0],
                    ) {
                        unsafe {
                            write_eclipse_how_attr(&how, attr);
                            *attr.add(3) = w.core_diameter_km;
                        }
                    } else {
                        unsafe { zero_f64_array(attr, 20) };
                    }
                } else {
                    unsafe { zero_f64_array(attr, 20) };
                }
                w.flags.bits() as i32
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
// swisseph_sol_eclipse_how — swe_sol_eclipse_how
// ---------------------------------------------------------------------------

/// Local circumstances of a solar eclipse at observer `geopos` at `tjd_ut` (UT1).
///
/// On success returns **positive** EclipseFlags bits; 0 = no eclipse visible.
/// Negative = error.
///
/// `attr[20]` out: [0]=magnitude, [1]=diameter_ratio, [2]=obscuration,
/// [3]=core_diameter_km, [4]=azimuth, [5]=true_altitude, [6]=apparent_altitude,
/// [7]=elongation, [8]=nasa_magnitude, [9]=saros_series, [10]=saros_member, [11..19]=0.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values [lon, lat, height].
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_sol_eclipse_how(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ifl: i32,
    geopos: *const f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || attr.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        match eph.sol_eclipse_how(tjd_ut, calc_flags, gp) {
            Ok(how) => {
                unsafe { write_eclipse_how_attr(&how, attr) };
                how.flags.bits() as i32
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
// swisseph_sol_eclipse_when_glob — swe_sol_eclipse_when_glob
// ---------------------------------------------------------------------------

/// Global solar eclipse search from `tjd_start` (UT1).
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `tret[10]` out: [0]=time_maximum, [1]=time_ra_conjunction, [2]=time_begin, [3]=time_end,
/// [4]=time_totality_begin, [5]=time_totality_end, [6]=time_centerline_begin,
/// [7]=time_centerline_end, [8..9]=0.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `tret` must point to at least 10 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_sol_eclipse_when_glob(
    handle: *const SweEphemeris,
    tjd_start: f64,
    ifl: i32,
    ifltype: i32,
    backward: i32,
    tret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || tret.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);
        let ecl_type = EclipseFlags::from_bits_retain(ifltype as u32);

        match eph.sol_eclipse_when_glob(tjd_start, calc_flags, ecl_type, backward != 0) {
            Ok(g) => {
                unsafe {
                    *tret = g.time_maximum;
                    *tret.add(1) = g.time_ra_conjunction;
                    *tret.add(2) = g.time_begin;
                    *tret.add(3) = g.time_end;
                    *tret.add(4) = g.time_totality_begin;
                    *tret.add(5) = g.time_totality_end;
                    *tret.add(6) = g.time_centerline_begin;
                    *tret.add(7) = g.time_centerline_end;
                    *tret.add(8) = 0.0;
                    *tret.add(9) = 0.0;
                }
                g.flags.bits() as i32
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
// swisseph_sol_eclipse_when_loc — swe_sol_eclipse_when_loc
// ---------------------------------------------------------------------------

/// Local solar eclipse search from `tjd_start` (UT1), visible from `geopos`.
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `tret[10]` out: [0]=time_maximum, [1]=time_first_contact, [2]=time_second_contact,
/// [3]=time_third_contact, [4]=time_fourth_contact, [5]=time_sunrise, [6]=time_sunset,
/// [7..9]=0. **Different slot semantics from sol_eclipse_when_glob.**
///
/// `attr[20]` out: local circumstances at maximum (same layout as sol_eclipse_how).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `tret` must point to at least 10 writable `f64` slots.
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_sol_eclipse_when_loc(
    handle: *const SweEphemeris,
    tjd_start: f64,
    ifl: i32,
    geopos: *const f64,
    backward: i32,
    tret: *mut f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || tret.is_null() || attr.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        match eph.sol_eclipse_when_loc(tjd_start, calc_flags, gp, backward != 0) {
            Ok(loc) => {
                unsafe {
                    *tret = loc.time_maximum;
                    *tret.add(1) = loc.time_first_contact;
                    *tret.add(2) = loc.time_second_contact;
                    *tret.add(3) = loc.time_third_contact;
                    *tret.add(4) = loc.time_fourth_contact;
                    *tret.add(5) = loc.time_sunrise;
                    *tret.add(6) = loc.time_sunset;
                    for i in 7..10 {
                        *tret.add(i) = 0.0;
                    }
                    write_eclipse_how_attr(&loc.attr, attr);
                }
                loc.flags.bits() as i32
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
// swisseph_lun_eclipse_how — swe_lun_eclipse_how
// ---------------------------------------------------------------------------

/// Local circumstances of a lunar eclipse at observer `geopos` at `tjd_ut` (UT1).
///
/// On success returns **positive** EclipseFlags bits; 0 = no eclipse visible.
/// Negative = error.
///
/// `attr[20]` out: [0]=umbral_magnitude, [1]=penumbral_magnitude, [2]=0, [3]=0,
/// [4]=azimuth, [5]=true_altitude, [6]=apparent_altitude, [7]=distance_from_opposition,
/// [8]=umbral_magnitude (dup), [9]=saros_series, [10]=saros_member, [11..19]=0.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_lun_eclipse_how(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ifl: i32,
    geopos: *const f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || attr.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        match eph.lun_eclipse_how(tjd_ut, calc_flags, gp) {
            Ok(how) => {
                unsafe { write_lun_eclipse_how_attr(&how, attr) };
                how.flags.bits() as i32
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
// swisseph_lun_eclipse_when — swe_lun_eclipse_when
// ---------------------------------------------------------------------------

/// Global lunar eclipse search from `tjd_start` (UT1).
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `tret[10]` out: [0]=time_maximum, [1]=0 (unused), [2]=time_partial_begin,
/// [3]=time_partial_end, [4]=time_totality_begin, [5]=time_totality_end,
/// [6]=time_penumbral_begin, [7]=time_penumbral_end, [8..9]=0.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `tret` must point to at least 10 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_lun_eclipse_when(
    handle: *const SweEphemeris,
    tjd_start: f64,
    ifl: i32,
    ifltype: i32,
    backward: i32,
    tret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || tret.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);
        let ecl_type = EclipseFlags::from_bits_retain(ifltype as u32);

        match eph.lun_eclipse_when(tjd_start, calc_flags, ecl_type, backward != 0) {
            Ok(g) => {
                unsafe {
                    *tret = g.time_maximum;
                    *tret.add(1) = 0.0;
                    *tret.add(2) = g.time_partial_begin;
                    *tret.add(3) = g.time_partial_end;
                    *tret.add(4) = g.time_totality_begin;
                    *tret.add(5) = g.time_totality_end;
                    *tret.add(6) = g.time_penumbral_begin;
                    *tret.add(7) = g.time_penumbral_end;
                    *tret.add(8) = 0.0;
                    *tret.add(9) = 0.0;
                }
                g.flags.bits() as i32
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
// swisseph_lun_eclipse_when_loc — swe_lun_eclipse_when_loc
// ---------------------------------------------------------------------------

/// Local lunar eclipse search from `tjd_start` (UT1), visible from `geopos`.
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `tret[10]` out: [0]=time_maximum, [1]=0 (unused), [2]=time_partial_begin,
/// [3]=time_partial_end, [4]=time_totality_begin, [5]=time_totality_end,
/// [6]=time_penumbral_begin, [7]=time_penumbral_end, [8]=time_moonrise,
/// [9]=time_moonset.
///
/// `attr[20]` out: lunar eclipse circumstances (same layout as lun_eclipse_how).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `tret` must point to at least 10 writable `f64` slots.
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_lun_eclipse_when_loc(
    handle: *const SweEphemeris,
    tjd_start: f64,
    ifl: i32,
    geopos: *const f64,
    backward: i32,
    tret: *mut f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || tret.is_null() || attr.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let eph = unsafe { &(*handle).0 };
        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        match eph.lun_eclipse_when_loc(tjd_start, calc_flags, gp, backward != 0) {
            Ok(loc) => {
                unsafe {
                    *tret = loc.time_maximum;
                    *tret.add(1) = 0.0;
                    *tret.add(2) = loc.time_partial_begin;
                    *tret.add(3) = loc.time_partial_end;
                    *tret.add(4) = loc.time_totality_begin;
                    *tret.add(5) = loc.time_totality_end;
                    *tret.add(6) = loc.time_penumbral_begin;
                    *tret.add(7) = loc.time_penumbral_end;
                    *tret.add(8) = loc.time_moonrise;
                    *tret.add(9) = loc.time_moonset;
                    write_lun_eclipse_how_attr(&loc.attr, attr);
                }
                loc.flags.bits() as i32
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
// swisseph_lun_occult_where — swe_lun_occult_where
// ---------------------------------------------------------------------------

/// Geographic position of maximal lunar occultation at `tjd_ut` (UT1).
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `geopos[10]` out: same layout as [`swisseph_sol_eclipse_where`].
/// `attr[20]` out: local circumstances at the central point, same as [`swisseph_sol_eclipse_where`].
/// `attr[3]` is `dcore[0]`. Zeroed when no eclipse is found.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos`, `attr` must point to at least 10/20 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_lun_occult_where(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    ipl: i32,
    starname: *const c_char,
    ifl: i32,
    geopos: *mut f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || attr.is_null() {
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

        let star = match unsafe { parse_starname_opt(starname) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        // Normalize asteroid-134340 → Pluto, matching lun_occult_where's internal alias.
        let body_for_how = match body {
            Body::Asteroid(id) if id.mpc_number() == 134340 => Body::Pluto,
            b => b,
        };

        match eph.lun_occult_where(tjd_ut, body, star, calc_flags) {
            Ok(w) => {
                unsafe {
                    write_eclipse_where_geopos(&w, geopos);
                }
                let ifl_masked =
                    calc_flags & (CalcFlags::JPLEPH | CalcFlags::SWIEPH | CalcFlags::MOSEPH);
                if !w.flags.is_empty() {
                    if let Ok(how) = eph.eclipse_how_at(
                        tjd_ut,
                        body_for_how,
                        star,
                        ifl_masked,
                        [w.central_longitude, w.central_latitude, 0.0],
                    ) {
                        unsafe {
                            write_eclipse_how_attr(&how, attr);
                            *attr.add(3) = w.core_diameter_km;
                        }
                    } else {
                        unsafe { zero_f64_array(attr, 20) };
                    }
                } else {
                    unsafe { zero_f64_array(attr, 20) };
                }
                w.flags.bits() as i32
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
// swisseph_lun_occult_when_glob — swe_lun_occult_when_glob
// ---------------------------------------------------------------------------

/// Global occultation search from `tjd_start` (UT1).
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `tret[10]` out: same slot layout as [`swisseph_sol_eclipse_when_glob`], but `tret[1]`
/// is the occulted body's transit instant (not specifically the Sun's).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `tret` must point to at least 10 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_lun_occult_when_glob(
    handle: *const SweEphemeris,
    tjd_start: f64,
    ipl: i32,
    starname: *const c_char,
    ifl: i32,
    ifltype: i32,
    backward: i32,
    tret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || tret.is_null() {
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

        let star = match unsafe { parse_starname_opt(starname) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);
        let ecl_type = EclipseFlags::from_bits_retain(ifltype as u32);

        match eph.lun_occult_when_glob(tjd_start, body, star, calc_flags, ecl_type, backward != 0) {
            Ok(g) => {
                unsafe {
                    *tret = g.time_maximum;
                    *tret.add(1) = g.time_ra_conjunction;
                    *tret.add(2) = g.time_begin;
                    *tret.add(3) = g.time_end;
                    *tret.add(4) = g.time_totality_begin;
                    *tret.add(5) = g.time_totality_end;
                    *tret.add(6) = g.time_centerline_begin;
                    *tret.add(7) = g.time_centerline_end;
                    *tret.add(8) = 0.0;
                    *tret.add(9) = 0.0;
                }
                g.flags.bits() as i32
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
// swisseph_lun_occult_when_loc — swe_lun_occult_when_loc
// ---------------------------------------------------------------------------

/// Local occultation search from `tjd_start` (UT1), visible from `geopos`.
///
/// On success returns **positive** EclipseFlags bits. Negative = error.
///
/// `tret[10]` out: [0]=time_maximum, [1]=time_first_contact, [2]=time_second_contact,
/// [3]=time_third_contact, [4]=time_fourth_contact, [5]=time_rise, [6]=time_set,
/// [7..9]=0. Same slot semantics as [`swisseph_sol_eclipse_when_loc`]; for a fixed star,
/// contacts 1/4 alias contacts 2/3.
///
/// `attr[20]` out: local circumstances at maximum (same layout as sol_eclipse_how).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `geopos` must point to 3 readable `f64` values.
/// - `tret` must point to at least 10 writable `f64` slots.
/// - `attr` must point to at least 20 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_lun_occult_when_loc(
    handle: *const SweEphemeris,
    tjd_start: f64,
    ipl: i32,
    starname: *const c_char,
    ifl: i32,
    geopos: *const f64,
    backward: i32,
    tret: *mut f64,
    attr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || geopos.is_null() || tret.is_null() || attr.is_null() {
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

        let star = match unsafe { parse_starname_opt(starname) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let gp = unsafe { [*geopos, *geopos.add(1), *geopos.add(2)] };
        let calc_flags = CalcFlags::from_bits_retain(ifl as u32);

        match eph.lun_occult_when_loc(tjd_start, body, star, calc_flags, gp, backward != 0) {
            Ok(loc) => {
                unsafe {
                    *tret = loc.time_maximum;
                    *tret.add(1) = loc.time_first_contact;
                    *tret.add(2) = loc.time_second_contact;
                    *tret.add(3) = loc.time_third_contact;
                    *tret.add(4) = loc.time_fourth_contact;
                    *tret.add(5) = loc.time_rise;
                    *tret.add(6) = loc.time_set;
                    for i in 7..10 {
                        *tret.add(i) = 0.0;
                    }
                    write_eclipse_how_attr(&loc.attr, attr);
                }
                loc.flags.bits() as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                error_code(&e)
            }
        }
    })
}
