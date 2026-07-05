use std::ffi::c_char;

use swisseph::date;
use swisseph::deltat::calc_deltat;
use swisseph::flags::CalcFlags;
use swisseph::format;
use swisseph::math;
use swisseph::sidereal_time;
use swisseph::types::{CalendarType, JdTt, JdUt1, UtcComponents};

use crate::SweEphemeris;
use crate::error::{SweErrorCode, ffi_guard, write_err};

// ---------------------------------------------------------------------------
// Handle-free date/time functions
// ---------------------------------------------------------------------------

/// Convert a calendar date + fractional hour (UT) to a Julian Day number.
/// `gregflag`: 0 = Julian calendar, 1 = Gregorian calendar.
///
/// # Safety
/// No pointer arguments.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_julday(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    gregflag: i32,
) -> f64 {
    let cal = if gregflag == 0 {
        CalendarType::Julian
    } else {
        CalendarType::Gregorian
    };
    date::julday(year, month, day, hour, cal)
}

/// Convert a Julian Day number to calendar date components.
/// `gregflag`: 0 = Julian calendar, 1 = Gregorian calendar.
///
/// # Safety
/// - `year`, `month`, `day` must point to writable `i32` slots.
/// - `hour` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_revjul(
    jd: f64,
    gregflag: i32,
    year: *mut i32,
    month: *mut i32,
    day: *mut i32,
    hour: *mut f64,
) {
    let cal = if gregflag == 0 {
        CalendarType::Julian
    } else {
        CalendarType::Gregorian
    };
    let (y, m, d, h) = date::revjul(jd, cal);
    unsafe {
        if !year.is_null() {
            *year = y;
        }
        if !month.is_null() {
            *month = m;
        }
        if !day.is_null() {
            *day = d;
        }
        if !hour.is_null() {
            *hour = h;
        }
    }
}

/// Validate and convert a calendar date to a Julian Day.
/// Returns 0 on success (writes JD to `*tjd`), or ERR on invalid date.
/// `cal`: 'g'/'G' = Gregorian, 'j'/'J' = Julian.
///
/// # Safety
/// - `tjd` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    cal: c_char,
    tjd: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if tjd.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }
        let cal_type = match cal as u8 {
            b'g' | b'G' => CalendarType::Gregorian,
            b'j' | b'J' => CalendarType::Julian,
            _ => {
                unsafe { write_err(err_buf, err_cap, "invalid calendar type (use 'g' or 'j')") };
                return SweErrorCode::InvalidCalendarType as i32;
            }
        };
        match date::date_conversion(year, month, day, hour, cal_type) {
            Ok(jd) => {
                unsafe { *tjd = jd };
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                crate::error::error_code(&e)
            }
        }
    })
}

/// Day of the week for a Julian Day: 0 = Monday .. 6 = Sunday.
///
/// # Safety
/// No pointer arguments.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_day_of_week(jd: f64) -> i32 {
    date::day_of_week(jd) as i32
}

/// Shift a UTC date/time by a timezone offset (hours, east positive).
/// All time components are both input and output (in-place conversion).
///
/// # Safety
/// - All pointer arguments must point to writable locations.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_utc_time_zone(
    iyear: i32,
    imonth: i32,
    iday: i32,
    ihour: i32,
    imin: i32,
    dsec: f64,
    d_timezone: f64,
    oyear: *mut i32,
    omonth: *mut i32,
    oday: *mut i32,
    ohour: *mut i32,
    omin: *mut i32,
    osec: *mut f64,
) {
    let input = UtcComponents {
        year: iyear,
        month: imonth,
        day: iday,
        hour: ihour,
        minute: imin,
        second: dsec,
    };
    let out = date::utc_time_zone(&input, d_timezone);
    unsafe {
        if !oyear.is_null() {
            *oyear = out.year;
        }
        if !omonth.is_null() {
            *omonth = out.month;
        }
        if !oday.is_null() {
            *oday = out.day;
        }
        if !ohour.is_null() {
            *ohour = out.hour;
        }
        if !omin.is_null() {
            *omin = out.minute;
        }
        if !osec.is_null() {
            *osec = out.second;
        }
    }
}

// ---------------------------------------------------------------------------
// Handle-based UTC / delta-T / sidereal time
// ---------------------------------------------------------------------------

/// Convert UTC to Julian Day (TT and UT1). Writes `dret[0]` = JD(TT), `dret[1]` = JD(UT1).
/// `gregflag`: 0 = Julian, 1 = Gregorian.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - `dret` must point to at least 2 writable `f64` slots.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_utc_to_jd(
    handle: *const SweEphemeris,
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    min: i32,
    sec: f64,
    gregflag: i32,
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
        let cal = if gregflag == 0 {
            CalendarType::Julian
        } else {
            CalendarType::Gregorian
        };
        let utc = UtcComponents {
            year,
            month,
            day,
            hour,
            minute: min,
            second: sec,
        };
        match date::utc_to_jd(&utc, cal, eph.leap_seconds(), eph) {
            Ok(result) => {
                unsafe {
                    *dret = result.tt.0;
                    *dret.add(1) = result.ut1.0;
                }
                SweErrorCode::Ok as i32
            }
            Err(e) => {
                let msg = e.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                crate::error::error_code(&e)
            }
        }
    })
}

/// Convert Julian Day (TT) to UTC calendar components.
/// `gregflag`: 0 = Julian, 1 = Gregorian.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - All out-param pointers must be writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_jdet_to_utc(
    handle: *const SweEphemeris,
    tjd_et: f64,
    gregflag: i32,
    year: *mut i32,
    month: *mut i32,
    day: *mut i32,
    hour: *mut i32,
    min: *mut i32,
    sec: *mut f64,
) {
    if handle.is_null() {
        return;
    }
    let eph = unsafe { &(*handle).0 };
    let cal = if gregflag == 0 {
        CalendarType::Julian
    } else {
        CalendarType::Gregorian
    };
    let utc = date::jdet_to_utc(JdTt(tjd_et), cal, eph.leap_seconds(), eph);
    unsafe {
        if !year.is_null() {
            *year = utc.year;
        }
        if !month.is_null() {
            *month = utc.month;
        }
        if !day.is_null() {
            *day = utc.day;
        }
        if !hour.is_null() {
            *hour = utc.hour;
        }
        if !min.is_null() {
            *min = utc.minute;
        }
        if !sec.is_null() {
            *sec = utc.second;
        }
    }
}

/// Convert Julian Day (UT1) to UTC calendar components.
/// `gregflag`: 0 = Julian, 1 = Gregorian.
///
/// # Safety
/// - `handle` must be a valid, non-NULL handle.
/// - All out-param pointers must be writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_jdut1_to_utc(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    gregflag: i32,
    year: *mut i32,
    month: *mut i32,
    day: *mut i32,
    hour: *mut i32,
    min: *mut i32,
    sec: *mut f64,
) {
    if handle.is_null() {
        return;
    }
    let eph = unsafe { &(*handle).0 };
    let cal = if gregflag == 0 {
        CalendarType::Julian
    } else {
        CalendarType::Gregorian
    };
    let utc = date::jdut1_to_utc(JdUt1(tjd_ut), cal, eph.leap_seconds(), eph);
    unsafe {
        if !year.is_null() {
            *year = utc.year;
        }
        if !month.is_null() {
            *month = utc.month;
        }
        if !day.is_null() {
            *day = utc.day;
        }
        if !hour.is_null() {
            *hour = utc.hour;
        }
        if !min.is_null() {
            *min = utc.minute;
        }
        if !sec.is_null() {
            *sec = utc.second;
        }
    }
}

/// Compute Delta T (TT - UT1) in days for `tjd_ut`.
/// Uses the handle's configured ephemeris model.
///
/// # Safety
/// `handle` must be a valid, non-NULL handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_deltat(handle: *const SweEphemeris, tjd_ut: f64) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    let eph = unsafe { &(*handle).0 };
    calc_deltat(tjd_ut, eph.config())
}

/// Compute Delta T (TT - UT1) in days for `tjd_ut`, with iflag-based ephemeris resolution.
/// Writes the result to `*deltat`. Returns 0 on success.
///
/// The `iflag` EPHMASK bits select which ephemeris backend's tidal acceleration to use.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `deltat` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_deltat_ex(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    iflag: i32,
    deltat: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || deltat.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }
        let eph = unsafe { &(*handle).0 };
        let calc_flags = CalcFlags::from_bits_retain(iflag as u32);
        let config = eph.effective_config(calc_flags, eph.config());
        let dt = calc_deltat(tjd_ut, &config);
        unsafe { *deltat = dt };
        SweErrorCode::Ok as i32
    })
}

/// Greenwich Apparent Sidereal Time (hours, 0..24) from `tjd_ut`.
/// Uses the handle's configured astro models. Forces TIDAL_DEFAULT internally
/// (matching C's `swe_sidtime` behavior).
///
/// # Safety
/// `handle` must be valid, non-NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_sidtime(handle: *const SweEphemeris, tjd_ut: f64) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    let eph = unsafe { &(*handle).0 };
    let mut deltat_config = eph.config().clone();
    deltat_config.tidal_acceleration = Some(swisseph::constants::TIDAL_DEFAULT);
    sidereal_time::sidereal_time(tjd_ut, &deltat_config)
}

/// Greenwich Apparent Sidereal Time (hours, 0..24) from pre-computed
/// true obliquity `eps` (degrees) and nutation in longitude `nut` (degrees).
/// Forces TIDAL_DEFAULT internally (matching C's `swe_sidtime0` behavior).
///
/// # Safety
/// `handle` must be valid, non-NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_sidtime0(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    eps: f64,
    nut: f64,
) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    let eph = unsafe { &(*handle).0 };
    let mut deltat_config = eph.config().clone();
    deltat_config.tidal_acceleration = Some(swisseph::constants::TIDAL_DEFAULT);
    sidereal_time::sidereal_time0(tjd_ut, eps, nut, &deltat_config)
}

/// Equation of time at `tjd_ut` (UT1). Returns `E = LAT − LMT` in **days**.
/// Positive means the Sun is ahead of the mean Sun.
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `e` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_time_equ(
    handle: *const SweEphemeris,
    tjd_ut: f64,
    e: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || e.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }
        let eph = unsafe { &(*handle).0 };
        match eph.time_equ(tjd_ut) {
            Ok(val) => {
                unsafe { *e = val };
                SweErrorCode::Ok as i32
            }
            Err(err) => {
                let msg = err.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                crate::error::error_code(&err)
            }
        }
    })
}

/// Convert Local Mean Time to Local Apparent Time.
/// `geolon` in degrees (east-positive). `tjd_lmt` and output are Julian Day (UT-scale).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `tjd_lat` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_lmt_to_lat(
    handle: *const SweEphemeris,
    tjd_lmt: f64,
    geolon: f64,
    tjd_lat: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || tjd_lat.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }
        let eph = unsafe { &(*handle).0 };
        match eph.lmt_to_lat(tjd_lmt, geolon) {
            Ok(val) => {
                unsafe { *tjd_lat = val };
                SweErrorCode::Ok as i32
            }
            Err(err) => {
                let msg = err.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                crate::error::error_code(&err)
            }
        }
    })
}

/// Convert Local Apparent Time to Local Mean Time.
/// `geolon` in degrees (east-positive). `tjd_lat` and output are Julian Day (UT-scale).
///
/// # Safety
/// - `handle` must be valid, non-NULL.
/// - `tjd_lmt` must point to a writable `f64`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_lat_to_lmt(
    handle: *const SweEphemeris,
    tjd_lat: f64,
    geolon: f64,
    tjd_lmt: *mut f64,
    err_buf: *mut c_char,
    err_cap: usize,
) -> i32 {
    ffi_guard!(err_buf, err_cap, {
        if handle.is_null() || tjd_lmt.is_null() {
            unsafe { write_err(err_buf, err_cap, "null pointer argument") };
            return SweErrorCode::InvalidArg as i32;
        }
        let eph = unsafe { &(*handle).0 };
        match eph.lat_to_lmt(tjd_lat, geolon) {
            Ok(val) => {
                unsafe { *tjd_lmt = val };
                SweErrorCode::Ok as i32
            }
            Err(err) => {
                let msg = err.to_string();
                unsafe { write_err(err_buf, err_cap, &msg) };
                crate::error::error_code(&err)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Formatting (centisecond-based, handle-free)
// ---------------------------------------------------------------------------

/// Round centiseconds to the nearest arcsecond with a 30°-boundary guard.
///
/// # Safety
/// No pointer arguments.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_csroundsec(x: i32) -> i32 {
    format::csroundsec(x)
}

/// Format time-of-day centiseconds as `"HH:MM:SS"` into `buf`.
/// `sep` is the separator character (e.g. ':').
/// `suppress_zero`: if true, omits seconds when they are zero.
///
/// # Safety
/// - `buf` must point to at least `cap` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_cs2timestr(
    t: i32,
    sep: c_char,
    suppress_zero: bool,
    buf: *mut c_char,
    cap: usize,
) {
    let sep_char = sep as u8 as char;
    let s = format::cs2timestr(t, sep_char, suppress_zero);
    unsafe { write_err(buf, cap, &s) };
}

/// Format arc centiseconds as a longitude/latitude string with direction letters.
/// `pchar` for positive, `mchar` for negative.
///
/// # Safety
/// - `buf` must point to at least `cap` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_cs2lonlatstr(
    t: i32,
    pchar: c_char,
    mchar: c_char,
    buf: *mut c_char,
    cap: usize,
) {
    let pc = pchar as u8 as char;
    let mc = mchar as u8 as char;
    let s = format::cs2lonlatstr(t, pc, mc);
    unsafe { write_err(buf, cap, &s) };
}

/// Format arc centiseconds as degrees-within-sign (`" D°MM'SS"`).
/// Truncates (no rounding) and wraps into [0, 30°).
///
/// # Safety
/// - `buf` must point to at least `cap` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_cs2degstr(t: i32, buf: *mut c_char, cap: usize) {
    let s = format::cs2degstr(t);
    unsafe { write_err(buf, cap, &s) };
}

/// Split a decimal degree value into components. Mirrors C's `swe_split_deg`.
///
/// # Parameters
/// - `ddeg`: decimal degrees to split
/// - `roundflag`: combination of SE_SPLIT_DEG_* flags
/// - `deg`, `min`, `sec`: out-params for integer components
/// - `secfr`: out-param for fractional seconds
/// - `sign`: out-param for sign (zodiacal sign index, or ±1)
///
/// # Safety
/// All out-param pointers must be writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_split_deg(
    ddeg: f64,
    roundflag: i32,
    deg: *mut i32,
    min: *mut i32,
    sec: *mut i32,
    secfr: *mut f64,
    sign: *mut i32,
) {
    use swisseph::flags::SplitDegFlags;
    let flags = SplitDegFlags::from_bits_retain(roundflag as u32);
    let parts = math::split_degrees(ddeg, flags);
    unsafe {
        if !deg.is_null() {
            *deg = parts.degrees;
        }
        if !min.is_null() {
            *min = parts.minutes;
        }
        if !sec.is_null() {
            *sec = parts.seconds;
        }
        if !secfr.is_null() {
            *secfr = parts.second_fraction;
        }
        if !sign.is_null() {
            *sign = parts.sign;
        }
    }
}
