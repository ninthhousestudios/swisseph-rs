//! Julian Day ↔ calendar conversion, delta-T dispatch, and UTC ↔ JD with
//! leap-second handling.

use crate::constants::{J1972, NLEAP_INIT};
use crate::error::Error;
use crate::types::{CalendarType, DeltaT, JdTt, JdUt1, UtcComponents, UtcToJd};

/// Built-in IERS leap-second table, as packed dates `YYYYMMDD` (the day the leap second is
/// inserted at the end of).
pub const LEAP_SECONDS: [i32; 27] = [
    19720630, 19721231, 19731231, 19741231, 19751231, 19761231, 19771231, 19781231, 19791231,
    19810630, 19820630, 19830630, 19850630, 19871231, 19891231, 19901231, 19920630, 19930630,
    19940630, 19951231, 19970630, 19981231, 20051231, 20081231, 20120630, 20150630, 20161231,
];

fn count_leap_seconds(ndat: i32, table: &[i32]) -> i32 {
    let mut count = NLEAP_INIT;
    for &entry in table {
        if ndat <= entry {
            break;
        }
        count += 1;
    }
    count
}

fn pack_date(year: i32, month: i32, day: i32) -> i32 {
    year * 10000 + month * 100 + day
}

fn split_jd_to_utc(jd: f64, cal: CalendarType) -> UtcComponents {
    let (year, month, day, d) = revjul(jd, cal);
    let hour = d as i32;
    let d = (d - hour as f64) * 60.0;
    let minute = d as i32;
    let second = (d - minute as f64) * 60.0;
    UtcComponents {
        year,
        month,
        day,
        hour,
        minute,
        second,
    }
}

// ---------------------------------------------------------------------------
// julday / revjul
// ---------------------------------------------------------------------------

/// Convert a calendar date + hour (fractional, UT) to a Julian Day. Port of `swe_julday`.
pub fn julday(year: i32, month: i32, day: i32, hour: f64, cal: CalendarType) -> f64 {
    let mut u = year as f64;
    if month < 3 {
        u -= 1.0;
    }
    let u0 = u + 4712.0;
    let mut u1 = (month + 1) as f64;
    if u1 < 4.0 {
        u1 += 12.0;
    }
    let mut jd =
        (u0 * 365.25).floor() + (30.6 * u1 + 0.000001).floor() + day as f64 + hour / 24.0 - 63.5;
    if cal == CalendarType::Gregorian {
        let mut u2 = (u.abs() / 100.0).floor() - (u.abs() / 400.0).floor();
        if u < 0.0 {
            u2 = -u2;
        }
        jd = jd - u2 + 2.0;
        if (u < 0.0) && (u / 100.0 == (u / 100.0).floor()) && (u / 400.0 != (u / 400.0).floor()) {
            jd -= 1.0;
        }
    }
    jd
}

/// Convert a Julian Day to `(year, month, day, hour)`, `hour` fractional UT. Port of `swe_revjul`.
pub fn revjul(jd: f64, cal: CalendarType) -> (i32, i32, i32, f64) {
    let mut u0 = jd + 32082.5;
    if cal == CalendarType::Gregorian {
        let mut u1 = u0 + (u0 / 36525.0).floor() - (u0 / 146100.0).floor() - 38.0;
        if jd >= 1830691.5 {
            u1 += 1.0;
        }
        u0 = u0 + (u1 / 36525.0).floor() - (u1 / 146100.0).floor() - 38.0;
    }
    let u2 = (u0 + 123.0).floor();
    let u3 = ((u2 - 122.2) / 365.25).floor();
    let u4 = ((u2 - (365.25 * u3).floor()) / 30.6001).floor();
    let mut month = (u4 - 1.0) as i32;
    if month > 12 {
        month -= 12;
    }
    let day = (u2 - (365.25 * u3).floor() - (30.6001 * u4).floor()) as i32;
    let year = (u3 + ((u4 - 2.0) / 12.0).floor() - 4800.0) as i32;
    let hour = (jd - (jd + 0.5).floor() + 0.5) * 24.0;
    (year, month, day, hour)
}

// ---------------------------------------------------------------------------
// date_conversion / day_of_week
// ---------------------------------------------------------------------------

/// Validate and convert a calendar date + hour to a Julian Day, rejecting dates that don't
/// round-trip through [`julday`]/[`revjul`] (e.g. February 30). Port of `swe_date_conversion`.
pub fn date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    cal: CalendarType,
) -> crate::Result<f64> {
    let jd = julday(year, month, day, hour, cal);
    let (ry, rm, rd, _) = revjul(jd, cal);
    if ry == year && rm == month && rd == day {
        Ok(jd)
    } else {
        Err(Error::InvalidDate {
            year,
            month,
            day: day as f64,
        })
    }
}

/// Day of the week for `jd`: `0` = Monday .. `6` = Sunday. Port of `swe_day_of_week`.
pub fn day_of_week(jd: f64) -> u8 {
    (((jd - 2433282.0 - 1.5).floor() as i64 % 7 + 7) % 7) as u8
}

// ---------------------------------------------------------------------------
// utc_time_zone
// ---------------------------------------------------------------------------

/// Shift a UTC calendar date/time by a time-zone offset (hours, east positive), handling
/// day rollover and preserving a leap-second flag. Port of `swe_utc_time_zone`.
pub fn utc_time_zone(input: &UtcComponents, tz_offset: f64) -> UtcComponents {
    let have_leapsec = input.second >= 60.0;
    let sec = if have_leapsec {
        input.second - 1.0
    } else {
        input.second
    };
    let mut dhour = input.hour as f64 + input.minute as f64 / 60.0 + sec / 3600.0;
    let mut tjd = julday(
        input.year,
        input.month,
        input.day,
        0.0,
        CalendarType::Gregorian,
    );
    dhour -= tz_offset;
    if dhour < 0.0 {
        tjd -= 1.0;
        dhour += 24.0;
    }
    if dhour >= 24.0 {
        tjd += 1.0;
        dhour -= 24.0;
    }
    let (year, month, day, _) = revjul(tjd + 0.001, CalendarType::Gregorian);
    let hour = dhour as i32;
    let d = (dhour - hour as f64) * 60.0;
    let minute = d as i32;
    let mut second = (d - minute as f64) * 60.0;
    if have_leapsec {
        second += 1.0;
    }
    UtcComponents {
        year,
        month,
        day,
        hour,
        minute,
        second,
    }
}

// ---------------------------------------------------------------------------
// utc_to_jd
// ---------------------------------------------------------------------------

/// Convert UTC calendar date/time to Julian Day in both TT and UT1, accounting for leap seconds
/// and (pre-1972) treating the input as UT1 directly. Port of `swe_utc_to_jd`.
pub fn utc_to_jd(
    utc: &UtcComponents,
    cal: CalendarType,
    leap_secs: &[i32],
    dt: &impl DeltaT,
) -> crate::Result<UtcToJd> {
    // Phase 1: Validate date via round-trip
    let tjd_ut1 = julday(utc.year, utc.month, utc.day, 0.0, cal);
    let (ry, rm, rd, _) = revjul(tjd_ut1, cal);
    if ry != utc.year || rm != utc.month || rd != utc.day {
        return Err(Error::InvalidDate {
            year: utc.year,
            month: utc.month,
            day: utc.day as f64,
        });
    }
    // Phase 2: Validate time
    if utc.hour < 0
        || utc.hour > 23
        || utc.minute < 0
        || utc.minute > 59
        || utc.second < 0.0
        || utc.second >= 61.0
        || (utc.second >= 60.0 && (utc.minute < 59 || utc.hour < 23 || tjd_ut1 < J1972))
    {
        return Err(Error::InvalidTime {
            hour: utc.hour,
            minute: utc.minute,
            second: utc.second,
        });
    }
    let dhour = utc.hour as f64 + utc.minute as f64 / 60.0 + utc.second / 3600.0;
    // Phase 3: Pre-1972 fast path — treat as UT1
    if tjd_ut1 < J1972 {
        let ut1 = JdUt1(julday(utc.year, utc.month, utc.day, dhour, cal));
        let tt = JdTt(ut1.0 + dt.delta_t(ut1));
        return Ok(UtcToJd { tt, ut1 });
    }
    // Phase 4: Convert to Gregorian for leap-sec counting
    let (gy, gm, gd) = if cal == CalendarType::Julian {
        let (y, m, d, _) = revjul(tjd_ut1, CalendarType::Gregorian);
        (y, m, d)
    } else {
        (utc.year, utc.month, utc.day)
    };
    // Phase 5: Count accumulated leap seconds
    let ndat = pack_date(gy, gm, gd);
    let nleap = count_leap_seconds(ndat, leap_secs);
    // Phase 6: Validate leap second 60 against table
    if utc.second >= 60.0 && !leap_secs.contains(&ndat) {
        return Err(Error::InvalidLeapSecond {
            year: gy,
            month: gm,
            day: gd,
        });
    }
    // Phase 7: Stale-table fallback
    let d = dt.delta_t(JdUt1(tjd_ut1)) * 86400.0;
    if d - nleap as f64 - 32.184 >= 1.0 {
        let ut1 = JdUt1(tjd_ut1 + dhour / 24.0);
        let tt = JdTt(ut1.0 + dt.delta_t(ut1));
        return Ok(UtcToJd { tt, ut1 });
    }
    // Phase 8: Compute TT and UT1
    let d = tjd_ut1 - J1972
        + utc.hour as f64 / 24.0
        + utc.minute as f64 / 1440.0
        + utc.second / 86400.0;
    let tjd_et_1972 = J1972 + (32.184 + NLEAP_INIT as f64) / 86400.0;
    let tjd_et = tjd_et_1972 + d + (nleap - NLEAP_INIT) as f64 / 86400.0;
    // Two-pass delta-T inversion: TT -> UT1
    let d1 = dt.delta_t(JdUt1(tjd_et));
    let d2 = dt.delta_t(JdUt1(tjd_et - d1));
    let tjd_ut1_final = tjd_et - dt.delta_t(JdUt1(tjd_et - d2));
    Ok(UtcToJd {
        tt: JdTt(tjd_et),
        ut1: JdUt1(tjd_ut1_final),
    })
}

// ---------------------------------------------------------------------------
// jdet_to_utc / jdut1_to_utc
// ---------------------------------------------------------------------------

/// Convert a Julian Day (TT) to UTC calendar date/time, counting accumulated leap seconds.
/// Port of `swe_jdet_to_utc`.
pub fn jdet_to_utc(
    jd_tt: JdTt,
    cal: CalendarType,
    leap_secs: &[i32],
    dt: &impl DeltaT,
) -> UtcComponents {
    let tjd_et_1972 = J1972 + (32.184 + NLEAP_INIT as f64) / 86400.0;
    // Two-pass TT -> UT1
    let d = dt.delta_t(JdUt1(jd_tt.0));
    let tjd_ut = jd_tt.0 - dt.delta_t(JdUt1(jd_tt.0 - d));
    let tjd_ut = jd_tt.0 - dt.delta_t(JdUt1(tjd_ut));
    // Pre-1972: output UT1 directly
    if jd_tt.0 < tjd_et_1972 {
        return split_jd_to_utc(tjd_ut, cal);
    }
    // Count leap seconds (conservative: use tjd_ut - 1 day)
    let (y2, m2, d2, _) = revjul(tjd_ut - 1.0, CalendarType::Gregorian);
    let ndat = pack_date(y2, m2, d2);
    let mut nleap = 0i32;
    for (i, &entry) in leap_secs.iter().enumerate() {
        if ndat <= entry {
            break;
        }
        nleap = (i + 1) as i32;
    }
    // Probe for the potentially-missing leap second
    let mut second_60 = 0;
    if (nleap as usize) < leap_secs.len() {
        let next_entry = leap_secs[nleap as usize];
        let ny = next_entry / 10000;
        let nm = (next_entry % 10000) / 100;
        let nd = next_entry % 100;
        let tjd_next = julday(ny, nm, nd, 0.0, CalendarType::Gregorian);
        let (y3, m3, d3, _) = revjul(tjd_next + 1.0, CalendarType::Gregorian);
        let boundary_utc = UtcComponents {
            year: y3,
            month: m3,
            day: d3,
            hour: 0,
            minute: 0,
            second: 0.0,
        };
        if let Ok(boundary) = utc_to_jd(&boundary_utc, CalendarType::Gregorian, leap_secs, dt) {
            let diff = jd_tt.0 - boundary.tt.0;
            if diff >= 0.0 {
                nleap += 1;
            } else if diff > -1.0 / 86400.0 {
                second_60 = 1;
            }
        }
    }
    // Convert TT -> UTC calendar
    let tjd = J1972 + (jd_tt.0 - tjd_et_1972) - (nleap + second_60) as f64 / 86400.0;
    let (mut year, mut month, mut day, d) = revjul(tjd, CalendarType::Gregorian);
    let hour = d as i32;
    let d = (d - hour as f64) * 60.0;
    let minute = d as i32;
    let second = (d - minute as f64) * 60.0 + second_60 as f64;
    // Stale-table fallback
    let d = dt.delta_t(JdUt1(jd_tt.0 - dt.delta_t(JdUt1(jd_tt.0))));
    if d * 86400.0 - (nleap + NLEAP_INIT) as f64 - 32.184 >= 1.0 {
        return split_jd_to_utc(jd_tt.0 - d, cal);
    }
    // Julian calendar output conversion
    if cal == CalendarType::Julian {
        let tjd = julday(year, month, day, 0.0, CalendarType::Gregorian);
        let (jy, jm, jd, _) = revjul(tjd, CalendarType::Julian);
        year = jy;
        month = jm;
        day = jd;
    }
    UtcComponents {
        year,
        month,
        day,
        hour,
        minute,
        second,
    }
}

/// Convert a Julian Day (UT1) to UTC calendar date/time. Port of `swe_jdut1_to_utc`.
pub fn jdut1_to_utc(
    jd_ut: JdUt1,
    cal: CalendarType,
    leap_secs: &[i32],
    dt: &impl DeltaT,
) -> UtcComponents {
    let jd_tt = JdTt(jd_ut.0 + dt.delta_t(jd_ut));
    jdet_to_utc(jd_tt, cal, leap_secs, dt)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    struct ConstDeltaT(f64);
    impl DeltaT for ConstDeltaT {
        fn delta_t(&self, _: JdUt1) -> f64 {
            self.0
        }
    }

    // -- julday --

    #[test]
    fn julday_j2000() {
        let jd = julday(2000, 1, 1, 12.0, CalendarType::Gregorian);
        assert_eq!(jd, 2451545.0);
    }

    #[test]
    fn julday_unix_epoch() {
        let jd = julday(1970, 1, 1, 0.0, CalendarType::Gregorian);
        assert_eq!(jd, 2440587.5);
    }

    #[test]
    fn julday_jd_epoch_julian() {
        let jd = julday(-4712, 1, 1, 12.0, CalendarType::Julian);
        assert!((jd - 0.0).abs() < 1e-9);
    }

    #[test]
    fn julday_gregorian_adoption() {
        let jd_before = julday(1582, 10, 4, 12.0, CalendarType::Julian);
        let jd_after = julday(1582, 10, 15, 12.0, CalendarType::Gregorian);
        assert!((jd_after - jd_before - 1.0).abs() < 1e-9);
    }

    // -- revjul round-trip --

    #[test]
    fn revjul_roundtrip_gregorian() {
        let cases = [
            (2000, 1, 1, 12.0),
            (1970, 1, 1, 0.0),
            (1582, 10, 15, 12.0),
            (2024, 6, 15, 6.5),
            (1, 1, 1, 0.0),
        ];
        for (y, m, d, h) in cases {
            let jd = julday(y, m, d, h, CalendarType::Gregorian);
            let (ry, rm, rd, rh) = revjul(jd, CalendarType::Gregorian);
            assert_eq!((ry, rm, rd), (y, m, d), "round-trip failed for {y}-{m}-{d}");
            assert!((rh - h).abs() < 1e-8, "hour mismatch for {y}-{m}-{d}");
        }
    }

    #[test]
    fn revjul_roundtrip_julian() {
        let cases = [(-4712, 1, 1, 12.0), (1582, 10, 4, 12.0), (100, 3, 1, 0.0)];
        for (y, m, d, h) in cases {
            let jd = julday(y, m, d, h, CalendarType::Julian);
            let (ry, rm, rd, rh) = revjul(jd, CalendarType::Julian);
            assert_eq!((ry, rm, rd), (y, m, d), "round-trip failed for {y}-{m}-{d}");
            assert!((rh - h).abs() < 1e-8, "hour mismatch for {y}-{m}-{d}");
        }
    }

    #[test]
    fn revjul_roundtrip_negative_years() {
        for y in [-1000, -500, -100, -1, 0] {
            let jd = julday(y, 6, 15, 12.0, CalendarType::Gregorian);
            let (ry, rm, rd, _) = revjul(jd, CalendarType::Gregorian);
            assert_eq!((ry, rm, rd), (y, 6, 15), "round-trip failed for year {y}");
        }
    }

    // -- date_conversion --

    #[test]
    fn date_conversion_valid() {
        assert!(date_conversion(2000, 1, 1, 0.0, CalendarType::Gregorian).is_ok());
    }

    #[test]
    fn date_conversion_feb30() {
        assert!(date_conversion(2000, 2, 30, 0.0, CalendarType::Gregorian).is_err());
    }

    #[test]
    fn date_conversion_month13() {
        assert!(date_conversion(2000, 13, 1, 0.0, CalendarType::Gregorian).is_err());
    }

    #[test]
    fn date_conversion_proleptic_gregorian() {
        // Proleptic Gregorian: Oct 10 1582 is mathematically valid (no historical gap check)
        assert!(date_conversion(1582, 10, 10, 0.0, CalendarType::Gregorian).is_ok());
        assert!(date_conversion(1582, 10, 4, 0.0, CalendarType::Gregorian).is_ok());
        assert!(date_conversion(1582, 10, 15, 0.0, CalendarType::Gregorian).is_ok());
    }

    // -- day_of_week --

    #[test]
    fn day_of_week_known_dates() {
        // 2000-01-01 was a Saturday (5)
        let jd = julday(2000, 1, 1, 12.0, CalendarType::Gregorian);
        assert_eq!(day_of_week(jd), 5);
        // 2024-01-01 was a Monday (0)
        let jd = julday(2024, 1, 1, 12.0, CalendarType::Gregorian);
        assert_eq!(day_of_week(jd), 0);
        // 1969-07-20 (Moon landing) was a Sunday (6)
        let jd = julday(1969, 7, 20, 12.0, CalendarType::Gregorian);
        assert_eq!(day_of_week(jd), 6);
    }

    // -- utc_time_zone --

    #[test]
    fn utc_time_zone_eastward_rollover() {
        let input = UtcComponents {
            year: 2024,
            month: 1,
            day: 1,
            hour: 22,
            minute: 0,
            second: 0.0,
        };
        let out = utc_time_zone(&input, -5.5);
        assert_eq!(out.day, 2);
        assert_eq!(out.hour, 3);
        assert_eq!(out.minute, 30);
    }

    #[test]
    fn utc_time_zone_westward_rollover() {
        let input = UtcComponents {
            year: 2024,
            month: 1,
            day: 1,
            hour: 2,
            minute: 0,
            second: 0.0,
        };
        let out = utc_time_zone(&input, 8.0);
        assert_eq!(out.year, 2023);
        assert_eq!(out.month, 12);
        assert_eq!(out.day, 31);
        assert_eq!(out.hour, 18);
    }

    #[test]
    fn utc_time_zone_preserves_leap_second() {
        let input = UtcComponents {
            year: 2024,
            month: 1,
            day: 1,
            hour: 5,
            minute: 59,
            second: 60.5,
        };
        let out = utc_time_zone(&input, 5.0);
        assert_eq!(out.hour, 0);
        assert_eq!(out.minute, 59);
        assert!(out.second >= 60.0);
    }

    // -- utc_to_jd --

    #[test]
    fn utc_to_jd_pre1972() {
        let dt = ConstDeltaT(42.184 / 86400.0);
        let utc = UtcComponents {
            year: 1960,
            month: 6,
            day: 15,
            hour: 12,
            minute: 0,
            second: 0.0,
        };
        let result = utc_to_jd(&utc, CalendarType::Gregorian, &LEAP_SECONDS, &dt).unwrap();
        let expected_ut1 = julday(1960, 6, 15, 12.0, CalendarType::Gregorian);
        assert!((result.ut1.0 - expected_ut1).abs() < 1e-12);
        assert!((result.tt.0 - result.ut1.0 - 42.184 / 86400.0).abs() < 1e-9);
    }

    #[test]
    fn utc_to_jd_invalid_date() {
        let dt = ConstDeltaT(0.0);
        let utc = UtcComponents {
            year: 2000,
            month: 2,
            day: 30,
            hour: 0,
            minute: 0,
            second: 0.0,
        };
        assert!(utc_to_jd(&utc, CalendarType::Gregorian, &LEAP_SECONDS, &dt).is_err());
    }

    #[test]
    fn utc_to_jd_invalid_leap_second() {
        let dt = ConstDeltaT(69.184 / 86400.0);
        let utc = UtcComponents {
            year: 2020,
            month: 6,
            day: 30,
            hour: 23,
            minute: 59,
            second: 60.0,
        };
        assert!(matches!(
            utc_to_jd(&utc, CalendarType::Gregorian, &LEAP_SECONDS, &dt),
            Err(Error::InvalidLeapSecond { .. })
        ));
    }

    #[test]
    fn utc_to_jd_invalid_leap_second_high_delta_t() {
        // Regression: a large delta-T must not bypass leap-second validation
        let dt = ConstDeltaT(75.184 / 86400.0);
        let utc = UtcComponents {
            year: 2020,
            month: 6,
            day: 30,
            hour: 23,
            minute: 59,
            second: 60.0,
        };
        assert!(matches!(
            utc_to_jd(&utc, CalendarType::Gregorian, &LEAP_SECONDS, &dt),
            Err(Error::InvalidLeapSecond { .. })
        ));
    }

    #[test]
    fn utc_to_jd_valid_leap_second() {
        let dt = ConstDeltaT(69.184 / 86400.0);
        let utc = UtcComponents {
            year: 2016,
            month: 12,
            day: 31,
            hour: 23,
            minute: 59,
            second: 60.0,
        };
        assert!(utc_to_jd(&utc, CalendarType::Gregorian, &LEAP_SECONDS, &dt).is_ok());
    }

    #[test]
    fn utc_to_jd_post1972_tt_ut1_relationship() {
        let delta_t_days = 69.184 / 86400.0;
        let dt = ConstDeltaT(delta_t_days);
        let utc = UtcComponents {
            year: 2020,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.0,
        };
        let result = utc_to_jd(&utc, CalendarType::Gregorian, &LEAP_SECONDS, &dt).unwrap();
        assert!(result.tt.0 > result.ut1.0);
    }

    // -- jdet_to_utc round-trip --

    #[test]
    fn jdet_to_utc_roundtrip_pre1972() {
        let dt = ConstDeltaT(42.184 / 86400.0);
        let utc_in = UtcComponents {
            year: 1960,
            month: 6,
            day: 15,
            hour: 12,
            minute: 30,
            second: 45.0,
        };
        let jd = utc_to_jd(&utc_in, CalendarType::Gregorian, &LEAP_SECONDS, &dt).unwrap();
        let utc_out = jdet_to_utc(jd.tt, CalendarType::Gregorian, &LEAP_SECONDS, &dt);
        assert_eq!(utc_out.year, utc_in.year);
        assert_eq!(utc_out.month, utc_in.month);
        assert_eq!(utc_out.day, utc_in.day);
        assert_eq!(utc_out.hour, utc_in.hour);
        assert_eq!(utc_out.minute, utc_in.minute);
        assert!((utc_out.second - utc_in.second).abs() < 0.001);
    }

    #[test]
    fn jdut1_to_utc_roundtrip_pre1972() {
        let dt = ConstDeltaT(42.184 / 86400.0);
        let utc_in = UtcComponents {
            year: 1960,
            month: 6,
            day: 15,
            hour: 12,
            minute: 30,
            second: 45.0,
        };
        let jd = utc_to_jd(&utc_in, CalendarType::Gregorian, &LEAP_SECONDS, &dt).unwrap();
        let utc_out = jdut1_to_utc(jd.ut1, CalendarType::Gregorian, &LEAP_SECONDS, &dt);
        assert_eq!(utc_out.year, utc_in.year);
        assert_eq!(utc_out.month, utc_in.month);
        assert_eq!(utc_out.day, utc_in.day);
        assert_eq!(utc_out.hour, utc_in.hour);
        assert_eq!(utc_out.minute, utc_in.minute);
        assert!((utc_out.second - utc_in.second).abs() < 0.001);
    }
}
