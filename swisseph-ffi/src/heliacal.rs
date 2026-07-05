use std::ffi::c_char;

use swisseph::flags::{CalcFlags, HeliacalFlags};
use swisseph::heliacal::HeliacalEventType;

use crate::error::{SweErrorCode, error_code, ffi_guard, write_err};
use crate::{SweEphemeris, cstr_to_str};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Split a combined C-style iflag into (CalcFlags for ephemeris source, HeliacalFlags).
/// C's heliacal API uses a single `iflag` where bits 0-2 are the ephemeris source
/// (JPLEPH=1, SWIEPH=2, MOSEPH=4) and bits 7+ are heliacal-specific flags.
fn split_helflag(iflag: i32) -> (CalcFlags, HeliacalFlags) {
    let bits = iflag as u32;
    let epheflag = CalcFlags::from_bits_retain(bits & 7);
    let helflag = HeliacalFlags::from_bits_retain(bits);
    (epheflag, helflag)
}

/// Read a `[f64; 3]` from a pointer, returning InvalidArg on null.
unsafe fn read_dgeo(dgeo: *const f64) -> Result<[f64; 3], i32> {
    if dgeo.is_null() {
        return Err(SweErrorCode::InvalidArg as i32);
    }
    Ok(unsafe { [*dgeo, *dgeo.add(1), *dgeo.add(2)] })
}

/// Read a mutable `[f64; 4]` from a pointer (copied in, modified by the lib, not written back).
unsafe fn read_datm(datm: *const f64) -> Result<[f64; 4], i32> {
    if datm.is_null() {
        return Err(SweErrorCode::InvalidArg as i32);
    }
    Ok(unsafe { [*datm, *datm.add(1), *datm.add(2), *datm.add(3)] })
}

/// Read a mutable `[f64; 6]` from a pointer.
unsafe fn read_dobs(dobs: *const f64) -> Result<[f64; 6], i32> {
    if dobs.is_null() {
        return Err(SweErrorCode::InvalidArg as i32);
    }
    Ok(unsafe {
        [
            *dobs,
            *dobs.add(1),
            *dobs.add(2),
            *dobs.add(3),
            *dobs.add(4),
            *dobs.add(5),
        ]
    })
}

// ---------------------------------------------------------------------------
// swisseph_heliacal_ut — swe_heliacal_ut
// ---------------------------------------------------------------------------

/// Find the next heliacal event for `object_name` after `tjd_start` (UT1).
///
/// # Parameters
/// - `dgeo`: `[longitude (°E+), latitude (°N+), altitude (m)]` — 3 readable `f64` values
/// - `datm`: `[pressure (hPa), temperature (°C), rel. humidity (%), extinction_coeff]` — 4 readable `f64` values
/// - `dobs`: `[age, Snellen_ratio, optic_type (0=eye/1=bino/2=tele), magnification, aperture_mm, transmission]` — 6 readable `f64` values
/// - `object_name`: NUL-terminated UTF-8 (planet name or star designation)
/// - `event_type`: `HeliacalEventType` (1=MorningFirst, 2=EveningLast, 3=EveningFirst,
///   4=MorningLast, 5=AcronymchalRising, 6=AcronymchalSetting)
/// - `helflag`: combined flag: bits 0–2 = ephemeris source (1=JPL, 2=Swiss, 4=Moshier),
///   bits 7+ = `SE_HELFLAG_*` heliacal flags
/// - `dret`: out-param, pointer to 50 `f64` slots. `dret[0]`=start_visible,
///   `dret[1]`=optimum_visibility, `dret[2]`=end_visible, `dret[3..49]`=0.
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dgeo` must point to 3 readable `f64` values.
/// - `datm` must point to 4 readable `f64` values.
/// - `dobs` must point to 6 readable `f64` values.
/// - `object_name` must be a valid NUL-terminated UTF-8 string.
/// - `dret` must point to at least 50 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_heliacal_ut(
    handle: *const SweEphemeris,
    tjd_start: f64,
    dgeo: *const f64,
    datm: *const f64,
    dobs: *const f64,
    object_name: *const c_char,
    event_type: i32,
    helflag: i32,
    dret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || dret.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let dgeo_arr = match unsafe { read_dgeo(dgeo) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dgeo pointer") };
                return code;
            }
        };
        let mut datm_arr = match unsafe { read_datm(datm) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null datm pointer") };
                return code;
            }
        };
        let mut dobs_arr = match unsafe { read_dobs(dobs) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dobs pointer") };
                return code;
            }
        };

        let name = match unsafe { cstr_to_str(object_name) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let event = match HeliacalEventType::try_from(event_type) {
            Ok(e) => e,
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let (epheflag, helflags) = split_helflag(helflag);
        let eph = unsafe { &(*handle).0 };

        match eph.heliacal_ut(
            tjd_start,
            &dgeo_arr,
            &mut datm_arr,
            &mut dobs_arr,
            name,
            event,
            epheflag,
            helflags,
        ) {
            Ok(result) => {
                unsafe {
                    *dret = result.start_visible;
                    *dret.add(1) = result.optimum_visibility;
                    *dret.add(2) = result.end_visible;
                    for i in 3..50 {
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
// swisseph_heliacal_pheno_ut — swe_heliacal_pheno_ut
// ---------------------------------------------------------------------------

/// Detailed heliacal-phenomena report at `tjd_ut` (UT1).
///
/// # Parameters
/// - `dgeo`: `[longitude (°E+), latitude (°N+), altitude (m)]` — 3 `f64`
/// - `datm`: `[pressure (hPa), temperature (°C), rel. humidity (%), extinction_coeff]` — 4 `f64`
/// - `dobs`: `[age, Snellen, optic_type, magnification, aperture_mm, transmission]` — 6 `f64`
/// - `object_name`: NUL-terminated UTF-8
/// - `event_type`: `HeliacalEventType` (1–6)
/// - `helflag`: combined ephemeris-source + heliacal flags
/// - `darr`: out-param, pointer to 50 `f64` slots. `darr[0..27]` receive the 28
///   `HeliacalPheno` fields in C's `dret[]` slot order:
///   `[0]`=topo_altitude, `[1]`=topo_apparent_altitude, `[2]`=geo_altitude,
///   `[3]`=azimuth_object, `[4]`=topo_sun_altitude, `[5]`=sun_azimuth,
///   `[6]`=TAV_actual, `[7]`=arcv_actual, `[8]`=DAZ_actual, `[9]`=arcl_actual,
///   `[10]`=extinction_coeff, `[11]`=min_TAV, `[12]`=t_first_visible,
///   `[13]`=t_best_visible, `[14]`=t_last_visible, `[15]`=t_best_yallop,
///   `[16]`=crescent_width, `[17]`=q_yallop, `[18]`=q_criterion,
///   `[19]`=parallax, `[20]`=magnitude, `[21]`=rise_object, `[22]`=rise_sun,
///   `[23]`=lag, `[24]`=visibility_duration, `[25]`=crescent_length,
///   `[26]`=elongation, `[27]`=illumination. `darr[28..49]`=0.
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dgeo` must point to 3 readable `f64` values.
/// - `datm` must point to 4 readable `f64` values.
/// - `dobs` must point to 6 readable `f64` values.
/// - `object_name` must be a valid NUL-terminated UTF-8 string.
/// - `darr` must point to at least 50 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_heliacal_pheno_ut(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    dgeo: *const f64,
    datm: *const f64,
    dobs: *const f64,
    object_name: *const c_char,
    event_type: i32,
    helflag: i32,
    darr: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || darr.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let dgeo_arr = match unsafe { read_dgeo(dgeo) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dgeo pointer") };
                return code;
            }
        };
        let mut datm_arr = match unsafe { read_datm(datm) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null datm pointer") };
                return code;
            }
        };
        let mut dobs_arr = match unsafe { read_dobs(dobs) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dobs pointer") };
                return code;
            }
        };

        let name = match unsafe { cstr_to_str(object_name) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let event = match HeliacalEventType::try_from(event_type) {
            Ok(e) => e,
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let (epheflag, helflags) = split_helflag(helflag);
        let eph = unsafe { &(*handle).0 };

        match eph.heliacal_pheno_ut(
            tjd_ut,
            &dgeo_arr,
            &mut datm_arr,
            &mut dobs_arr,
            name,
            event,
            epheflag,
            helflags,
        ) {
            Ok(pheno) => {
                let arr = pheno.as_array();
                unsafe {
                    for i in 0..28 {
                        *darr.add(i) = arr[i];
                    }
                    for i in 28..50 {
                        *darr.add(i) = 0.0;
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
// swisseph_vis_limit_mag — swe_vis_limit_mag
// ---------------------------------------------------------------------------

/// Visual limiting magnitude at `tjd_ut` (UT1).
///
/// # Parameters
/// - `dgeo`: `[longitude (°E+), latitude (°N+), altitude (m)]` — 3 `f64`
/// - `datm`: `[pressure (hPa), temperature (°C), rel. humidity (%), extinction_coeff]` — 4 `f64`
/// - `dobs`: `[age, Snellen, optic_type, magnification, aperture_mm, transmission]` — 6 `f64`
/// - `object_name`: NUL-terminated UTF-8
/// - `helflag`: combined ephemeris-source + heliacal flags
/// - `dret`: out-param, pointer to 8 `f64` slots:
///   `[0]`=limiting_magnitude, `[1]`=altitude_object, `[2]`=azimuth_object,
///   `[3]`=altitude_sun, `[4]`=azimuth_sun, `[5]`=altitude_moon,
///   `[6]`=azimuth_moon, `[7]`=magnitude_object.
///
/// Returns the vision-mode flags (non-negative: 0=photopic, 1=scotopic, 2=mixed) on success,
/// `-2` when the object is below the horizon (`dret[0]` = −100), or a negative error code
/// on failure — mirroring C's `swe_vis_limit_mag` return convention.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dgeo` must point to 3 readable `f64` values.
/// - `datm` must point to 4 readable `f64` values.
/// - `dobs` must point to 6 readable `f64` values.
/// - `object_name` must be a valid NUL-terminated UTF-8 string.
/// - `dret` must point to at least 8 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_vis_limit_mag(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    dgeo: *const f64,
    datm: *const f64,
    dobs: *const f64,
    object_name: *const c_char,
    helflag: i32,
    dret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || dret.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let dgeo_arr = match unsafe { read_dgeo(dgeo) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dgeo pointer") };
                return code;
            }
        };
        let mut datm_arr = match unsafe { read_datm(datm) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null datm pointer") };
                return code;
            }
        };
        let mut dobs_arr = match unsafe { read_dobs(dobs) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dobs pointer") };
                return code;
            }
        };

        let name = match unsafe { cstr_to_str(object_name) } {
            Ok(s) => s,
            Err(msg) => {
                unsafe { write_err(err_buf, err_cap, msg) };
                return SweErrorCode::InvalidArg as i32;
            }
        };

        let (epheflag, helflags) = split_helflag(helflag);
        let eph = unsafe { &(*handle).0 };

        match eph.vis_limit_mag(
            tjd_ut,
            &dgeo_arr,
            &mut datm_arr,
            &mut dobs_arr,
            name,
            epheflag,
            helflags,
        ) {
            Ok(result) => {
                unsafe {
                    *dret = result.limiting_magnitude;
                    *dret.add(1) = result.altitude_object;
                    *dret.add(2) = result.azimuth_object;
                    *dret.add(3) = result.altitude_sun;
                    *dret.add(4) = result.azimuth_sun;
                    *dret.add(5) = result.altitude_moon;
                    *dret.add(6) = result.azimuth_moon;
                    *dret.add(7) = result.magnitude_object;
                }
                if result.below_horizon {
                    -2
                } else {
                    result.vision.bits() as i32
                }
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
// swisseph_heliacal_angle — swe_heliacal_angle
// ---------------------------------------------------------------------------

/// Heliacal angle (optimal altitude / arcus visionis) at `tjd_ut` (UT1).
///
/// # Parameters
/// - `dgeo`: `[longitude (°E+), latitude (°N+), altitude (m)]` — 3 `f64`
/// - `datm`: `[pressure (hPa), temperature (°C), rel. humidity (%), extinction_coeff]` — 4 `f64`
/// - `dobs`: `[age, Snellen, optic_type, magnification, aperture_mm, transmission]` — 6 `f64`
/// - `helflag`: combined ephemeris-source + heliacal flags
/// - `mag`: object's visual magnitude
/// - `azi_obj`: object's azimuth (degrees)
/// - `azi_sun`: Sun's azimuth (degrees)
/// - `azi_moon`: Moon's azimuth (degrees)
/// - `alt_moon`: Moon's altitude (degrees)
/// - `dret`: out-param, pointer to 3 `f64` slots:
///   `[0]`=optimal_altitude, `[1]`=arcus_visionis, `[2]`=sun_altitude_diff.
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dgeo` must point to 3 readable `f64` values.
/// - `datm` must point to 4 readable `f64` values.
/// - `dobs` must point to 6 readable `f64` values.
/// - `dret` must point to at least 3 writable `f64` slots.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_heliacal_angle(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    dgeo: *const f64,
    datm: *const f64,
    dobs: *const f64,
    helflag: i32,
    mag: f64,
    azi_obj: f64,
    azi_sun: f64,
    azi_moon: f64,
    alt_moon: f64,
    dret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || dret.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let dgeo_arr = match unsafe { read_dgeo(dgeo) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dgeo pointer") };
                return code;
            }
        };
        let mut datm_arr = match unsafe { read_datm(datm) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null datm pointer") };
                return code;
            }
        };
        let mut dobs_arr = match unsafe { read_dobs(dobs) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dobs pointer") };
                return code;
            }
        };

        let (_, helflags) = split_helflag(helflag);
        let eph = unsafe { &(*handle).0 };

        match eph.heliacal_angle(
            tjd_ut,
            &dgeo_arr,
            &mut datm_arr,
            &mut dobs_arr,
            helflags,
            mag,
            azi_obj,
            azi_sun,
            azi_moon,
            alt_moon,
        ) {
            Ok(result) => {
                unsafe {
                    *dret = result.optimal_altitude;
                    *dret.add(1) = result.arcus_visionis;
                    *dret.add(2) = result.sun_altitude_diff;
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
// swisseph_topo_arcus_visionis — swe_topo_arcus_visionis
// ---------------------------------------------------------------------------

/// Topocentric arcus visionis at `tjd_ut` (UT1), degrees.
///
/// All geometry is caller-supplied. Angles in degrees.
///
/// # Parameters
/// - `dgeo`: `[longitude (°E+), latitude (°N+), altitude (m)]` — 3 `f64`
/// - `datm`: `[pressure (hPa), temperature (°C), rel. humidity (%), extinction_coeff]` — 4 `f64`
/// - `dobs`: `[age, Snellen, optic_type, magnification, aperture_mm, transmission]` — 6 `f64`
/// - `helflag`: combined ephemeris-source + heliacal flags
/// - `mag`: object's visual magnitude
/// - `azi_obj`: object's azimuth (degrees)
/// - `alt_obj`: object's altitude (degrees)
/// - `azi_sun`: Sun's azimuth (degrees)
/// - `azi_moon`: Moon's azimuth (degrees)
/// - `alt_moon`: Moon's altitude (degrees)
/// - `dret`: out-param, pointer to a writable `f64` receiving the arcus visionis (degrees)
///
/// Returns 0 on success, negative error code on failure.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `dgeo` must point to 3 readable `f64` values.
/// - `datm` must point to 4 readable `f64` values.
/// - `dobs` must point to 6 readable `f64` values.
/// - `dret` must point to a writable `f64`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swisseph_topo_arcus_visionis(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    dgeo: *const f64,
    datm: *const f64,
    dobs: *const f64,
    helflag: i32,
    mag: f64,
    azi_obj: f64,
    alt_obj: f64,
    azi_sun: f64,
    azi_moon: f64,
    alt_moon: f64,
    dret: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || dret.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }

        let dgeo_arr = match unsafe { read_dgeo(dgeo) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dgeo pointer") };
                return code;
            }
        };
        let mut datm_arr = match unsafe { read_datm(datm) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null datm pointer") };
                return code;
            }
        };
        let mut dobs_arr = match unsafe { read_dobs(dobs) } {
            Ok(a) => a,
            Err(code) => {
                unsafe { write_err(err_buf, err_cap, "null dobs pointer") };
                return code;
            }
        };

        let (_, helflags) = split_helflag(helflag);
        let eph = unsafe { &(*handle).0 };

        match eph.topo_arcus_visionis(
            tjd_ut,
            &dgeo_arr,
            &mut datm_arr,
            &mut dobs_arr,
            helflags,
            mag,
            azi_obj,
            alt_obj,
            azi_sun,
            azi_moon,
            alt_moon,
        ) {
            Ok(tav) => {
                unsafe { *dret = tav };
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
