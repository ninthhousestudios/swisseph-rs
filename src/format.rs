//! Degree/time string formatting (centisecond precision).
//!
//! Low-level internals; exposed for golden tests and advanced use.

const DEG30: i64 = 30 * 360000;

/// Round centiseconds to the nearest arcsecond (multiple of 100), with a
/// zodiac-sign-boundary guard: if rounding up would land exactly on a 30-degree
/// boundary, round down instead. Port of `swe_csroundsec` (swephlib.c:3836–3843).
pub fn csroundsec(x: i32) -> i32 {
    let x = x as i64;
    let t = (x + 50) / 100 * 100;
    let result = if t > x && t % DEG30 == 0 {
        x / 100 * 100
    } else {
        t
    };
    result as i32
}

/// Format time-of-day centiseconds as `"HH:MM:SS"` (or `"HH:MM"` when
/// `suppress_zero` is true and seconds are zero). `sep` is the separator
/// character placed between HH/MM and MM/SS.
/// Port of `swe_cs2timestr` (swephlib.c:3864–3886).
pub fn cs2timestr(t: i32, sep: char, suppress_zero: bool) -> String {
    let t = ((t as i64 + 50) / 100) % (24 * 3600);
    let s = t % 60;
    let m = (t / 60) % 60;
    let h = t / 3600 % 100;

    if s == 0 && suppress_zero {
        format!("{:02}{}{:02}", h, sep, m)
    } else {
        format!("{:02}{}{:02}{}{:02}", h, sep, m, sep, s)
    }
}

/// Format arc centiseconds as a longitude/latitude string with a direction
/// letter. Negative values use `mchar`; positive use `pchar`. Seconds are
/// suppressed when zero. Leading degree zeros are suppressed.
/// Port of `swe_cs2lonlatstr` (swephlib.c:3888–3916).
pub fn cs2lonlatstr(t: i32, pchar: char, mchar: char) -> String {
    let dir = if t < 0 { mchar } else { pchar };
    let t = ((t as i64).unsigned_abs() as i64 + 50) / 100;
    let s = t % 60;
    let m = (t / 60) % 60;
    let h = t / 3600 % 1000;

    if s == 0 {
        format!("{}{}{:02}", h, dir, m)
    } else {
        format!("{}{}{:02}'{:02}", h, dir, m, s)
    }
}

/// Format arc centiseconds as degrees-within-sign (`" D°MM'SS"`).
/// **Truncates** (no rounding) and wraps into `[0, 30°)`.
/// Port of `swe_cs2degstr` (swephlib.c:3918–3929).
pub fn cs2degstr(t: i32) -> String {
    let t = (t as i64) / 100 % (30 * 3600);
    let s = t % 60;
    let m = (t / 60) % 60;
    let h = t / 3600 % 100;
    format!("{:2}\u{b0}{:02}'{:02}", h, m, s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csroundsec_basic() {
        assert_eq!(csroundsec(100), 100);
        assert_eq!(csroundsec(149), 100);
        assert_eq!(csroundsec(150), 200);
        assert_eq!(csroundsec(0), 0);
    }

    #[test]
    fn csroundsec_sign_boundary_guard() {
        // 29°59'59.50" = 29*360000 + 59*6000 + 59*100 + 50
        let near_30 = 29 * 360000 + 59 * 6000 + 59 * 100 + 50;
        let rounded = csroundsec(near_30);
        assert!((rounded as i64) < DEG30);
        assert_eq!(rounded, near_30 / 100 * 100);
    }

    #[test]
    fn cs2timestr_basic() {
        // 12:30:00 = 12*3600*100 + 30*60*100 = 4320000 + 180000 = 4500000
        let cs = 12 * 3600 * 100 + 30 * 60 * 100;
        assert_eq!(cs2timestr(cs, ':', false), "12:30:00");
        assert_eq!(cs2timestr(cs, ':', true), "12:30");
    }

    #[test]
    fn cs2timestr_with_seconds() {
        let cs = 12 * 3600 * 100 + 30 * 60 * 100 + 45 * 100;
        assert_eq!(cs2timestr(cs, ':', false), "12:30:45");
        assert_eq!(cs2timestr(cs, ':', true), "12:30:45");
    }

    #[test]
    fn cs2lonlatstr_positive() {
        let cs = 122 * 360000 + 30 * 6000 + 45 * 100;
        assert_eq!(cs2lonlatstr(cs, 'E', 'W'), "122E30'45");
    }

    #[test]
    fn cs2lonlatstr_negative() {
        let cs = -(5 * 360000 + 30 * 6000);
        assert_eq!(cs2lonlatstr(cs, 'N', 'S'), "5S30");
    }

    #[test]
    fn cs2degstr_basic() {
        let cs = 15 * 360000 + 30 * 6000 + 45 * 100;
        assert_eq!(cs2degstr(cs), "15\u{b0}30'45");
    }

    #[test]
    fn cs2degstr_truncates() {
        let cs = 15 * 360000 + 30 * 6000 + 45 * 100 + 99;
        assert_eq!(cs2degstr(cs), "15\u{b0}30'45");
    }

    #[test]
    fn cs2degstr_wraps_at_30() {
        let cs = 31 * 360000;
        assert_eq!(cs2degstr(cs), " 1\u{b0}00'00");
    }
}
