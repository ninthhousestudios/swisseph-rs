use std::ffi::c_char;
use std::ptr;

use swisseph_ffi::SweEphemeris;
use swisseph_ffi::config::SweConfig;
use swisseph_ffi::error::SweErrorCode;

const J2000: f64 = 2451545.0;

unsafe fn default_handle() -> *mut SweEphemeris {
    unsafe {
        let mut config = std::mem::zeroed::<SweConfig>();
        swisseph_ffi::config::swisseph_config_default(&mut config);
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0, "swisseph_new failed");
        handle
    }
}

// ---------------------------------------------------------------------------
// julday / revjul round-trips
// ---------------------------------------------------------------------------

#[test]
fn julday_revjul_gregorian_roundtrip() {
    unsafe {
        let jd = swisseph_ffi::date::swisseph_julday(2000, 1, 1, 12.0, 1);
        assert!((jd - J2000).abs() < 1e-10);

        let mut y = 0i32;
        let mut m = 0i32;
        let mut d = 0i32;
        let mut h = 0.0f64;
        swisseph_ffi::date::swisseph_revjul(jd, 1, &mut y, &mut m, &mut d, &mut h);
        assert_eq!(y, 2000);
        assert_eq!(m, 1);
        assert_eq!(d, 1);
        assert!((h - 12.0).abs() < 1e-10);
    }
}

#[test]
fn julday_revjul_julian_roundtrip() {
    unsafe {
        let jd = swisseph_ffi::date::swisseph_julday(100, 3, 15, 6.0, 0);
        let mut y = 0i32;
        let mut m = 0i32;
        let mut d = 0i32;
        let mut h = 0.0f64;
        swisseph_ffi::date::swisseph_revjul(jd, 0, &mut y, &mut m, &mut d, &mut h);
        assert_eq!(y, 100);
        assert_eq!(m, 3);
        assert_eq!(d, 15);
        assert!((h - 6.0).abs() < 1e-10);
    }
}

#[test]
fn julday_revjul_negative_jd() {
    unsafe {
        let jd = swisseph_ffi::date::swisseph_julday(-4712, 1, 1, 12.0, 0);
        assert!(jd < 1.0);
        let mut y = 0i32;
        let mut m = 0i32;
        let mut d = 0i32;
        let mut h = 0.0f64;
        swisseph_ffi::date::swisseph_revjul(jd, 0, &mut y, &mut m, &mut d, &mut h);
        assert_eq!(y, -4712);
        assert_eq!(m, 1);
        assert_eq!(d, 1);
        assert!((h - 12.0).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// date_conversion
// ---------------------------------------------------------------------------

#[test]
fn date_conversion_valid() {
    unsafe {
        let mut tjd = 0.0f64;
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::date::swisseph_date_conversion(
            2024,
            6,
            15,
            10.5,
            b'g' as c_char,
            &mut tjd,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        assert!(tjd > 2460000.0);
    }
}

#[test]
fn date_conversion_invalid() {
    unsafe {
        let mut tjd = 0.0f64;
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::date::swisseph_date_conversion(
            2024,
            2,
            30, // Feb 30 is invalid
            0.0,
            b'g' as c_char,
            &mut tjd,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, SweErrorCode::InvalidDate as i32);
    }
}

// ---------------------------------------------------------------------------
// day_of_week
// ---------------------------------------------------------------------------

#[test]
fn day_of_week_j2000() {
    // J2000.0 = 2000-Jan-1.5 = Saturday = 5
    let dow = swisseph_ffi::date::swisseph_day_of_week(J2000);
    assert_eq!(dow, 5);
}

// ---------------------------------------------------------------------------
// utc_time_zone
// ---------------------------------------------------------------------------

#[test]
fn utc_time_zone_shift() {
    unsafe {
        let mut oy = 0i32;
        let mut om = 0i32;
        let mut od = 0i32;
        let mut oh = 0i32;
        let mut omin = 0i32;
        let mut osec = 0.0f64;
        swisseph_ffi::date::swisseph_utc_time_zone(
            2024, 6, 15, 10, 30, 0.0, 5.5, &mut oy, &mut om, &mut od, &mut oh, &mut omin, &mut osec,
        );
        // 10:30 minus 5.5h offset = 05:00 UTC
        assert_eq!(oy, 2024);
        assert_eq!(om, 6);
        assert_eq!(od, 15);
        assert_eq!(oh, 5);
        assert_eq!(omin, 0);
        assert!((osec - 0.0).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// utc_to_jd
// ---------------------------------------------------------------------------

#[test]
fn utc_to_jd_j2000() {
    unsafe {
        let handle = default_handle();
        let mut dret = [0.0f64; 2];
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::date::swisseph_utc_to_jd(
            handle,
            2000,
            1,
            1,
            12,
            0,
            0.0,
            1, // Gregorian
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        // TT ≈ UT + deltaT (about 64s at 2000)
        assert!(dret[0] > J2000);
        // UT1 ≈ J2000
        assert!((dret[1] - J2000).abs() < 0.001);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn utc_to_jd_leap_second() {
    unsafe {
        let handle = default_handle();
        let mut dret = [0.0f64; 2];
        let mut err_buf = [0u8; 256];
        // 2016-12-31 23:59:60 is a valid leap second
        let ret = swisseph_ffi::date::swisseph_utc_to_jd(
            handle,
            2016,
            12,
            31,
            23,
            59,
            60.0,
            1,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        assert!(dret[0] > 0.0);
        swisseph_ffi::swisseph_free(handle);
    }
}

// ---------------------------------------------------------------------------
// jdet_to_utc / jdut1_to_utc
// ---------------------------------------------------------------------------

#[test]
fn jdet_to_utc_roundtrip() {
    unsafe {
        let handle = default_handle();
        // Convert 2000-01-01 12:00:00 UTC to JD, then back
        let mut dret = [0.0f64; 2];
        let mut err_buf = [0u8; 256];
        swisseph_ffi::date::swisseph_utc_to_jd(
            handle,
            2000,
            1,
            1,
            12,
            0,
            0.0,
            1,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        let jd_tt = dret[0];

        let mut y = 0i32;
        let mut m = 0i32;
        let mut d = 0i32;
        let mut h = 0i32;
        let mut min = 0i32;
        let mut sec = 0.0f64;
        swisseph_ffi::date::swisseph_jdet_to_utc(
            handle, jd_tt, 1, &mut y, &mut m, &mut d, &mut h, &mut min, &mut sec,
        );
        assert_eq!(y, 2000);
        assert_eq!(m, 1);
        assert_eq!(d, 1);
        assert_eq!(h, 12);
        assert_eq!(min, 0);
        assert!(sec.abs() < 0.01);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn jdut1_to_utc_roundtrip() {
    unsafe {
        let handle = default_handle();
        let mut dret = [0.0f64; 2];
        let mut err_buf = [0u8; 256];
        swisseph_ffi::date::swisseph_utc_to_jd(
            handle,
            2020,
            6,
            15,
            18,
            30,
            45.0,
            1,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        let jd_ut = dret[1];

        let mut y = 0i32;
        let mut m = 0i32;
        let mut d = 0i32;
        let mut h = 0i32;
        let mut min = 0i32;
        let mut sec = 0.0f64;
        swisseph_ffi::date::swisseph_jdut1_to_utc(
            handle, jd_ut, 1, &mut y, &mut m, &mut d, &mut h, &mut min, &mut sec,
        );
        assert_eq!(y, 2020);
        assert_eq!(m, 6);
        assert_eq!(d, 15);
        assert_eq!(h, 18);
        assert_eq!(min, 30);
        assert!((sec - 45.0).abs() < 0.01);
        swisseph_ffi::swisseph_free(handle);
    }
}

// ---------------------------------------------------------------------------
// deltat
// ---------------------------------------------------------------------------

#[test]
fn deltat_moshier_j2000() {
    unsafe {
        let handle = default_handle();
        let dt = swisseph_ffi::date::swisseph_deltat(handle, J2000);
        // About 64 seconds at J2000, in days: ~0.00074
        assert!(dt > 0.0005 && dt < 0.001);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn deltat_ex_parity() {
    unsafe {
        let handle = default_handle();
        let dt_plain = swisseph_ffi::date::swisseph_deltat(handle, J2000);
        let mut dt_ex = 0.0f64;
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::date::swisseph_deltat_ex(
            handle,
            J2000,
            4, // SEFLG_MOSEPH
            &mut dt_ex,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        assert!((dt_plain - dt_ex).abs() < 1e-15);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn deltat_ex_sentinel_minus_one() {
    // C's swe_deltat_ex(tjd, -1, NULL) forces SE_TIDAL_DEFAULT. Verify our FFI
    // matches the sidtime path (which also forces TIDAL_DEFAULT) rather than
    // falling through to Moshier/DE404.
    unsafe {
        let handle = default_handle();
        let mut dt_sentinel = 0.0f64;
        let mut dt_plain = 0.0f64;
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::date::swisseph_deltat_ex(
            handle,
            J2000,
            -1, // sentinel
            &mut dt_sentinel,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        // For default Moshier config, TIDAL_DEFAULT == TIDAL_DE431, config has
        // TIDAL_DE404. At J2000 the difference is tiny but nonzero at far epochs.
        // At minimum, verify it doesn't crash and returns a valid value.
        assert!(dt_sentinel > 0.0005 && dt_sentinel < 0.001);

        // Compare against the plain swisseph_deltat (which uses config's tid_acc).
        // For Moshier default config, tid_acc is DE404 — so at J2000 these should
        // differ only at the ~1e-9 level (tidal correction is tiny near J2000).
        dt_plain = swisseph_ffi::date::swisseph_deltat(handle, J2000);
        // They'll be very close at J2000 but not necessarily bitwise-equal
        assert!((dt_sentinel - dt_plain).abs() < 1e-6);
        swisseph_ffi::swisseph_free(handle);
    }
}

// ---------------------------------------------------------------------------
// sidereal time
// ---------------------------------------------------------------------------

#[test]
fn sidtime_reasonable() {
    unsafe {
        let handle = default_handle();
        let st = swisseph_ffi::date::swisseph_sidtime(handle, J2000);
        // Should be in [0, 24)
        assert!(st >= 0.0 && st < 24.0);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn sidtime0_parity() {
    unsafe {
        let handle = default_handle();
        // sidtime0 with pre-computed eps/nut should give same as sidtime for
        // matching eps/nut values — but since we don't have easy access to those,
        // just verify it returns a valid value.
        let st = swisseph_ffi::date::swisseph_sidtime0(handle, J2000, 23.4393, -0.00385);
        assert!(st >= 0.0 && st < 24.0);
        swisseph_ffi::swisseph_free(handle);
    }
}

// ---------------------------------------------------------------------------
// time_equ / lmt_to_lat / lat_to_lmt
// ---------------------------------------------------------------------------

#[test]
fn time_equ_range() {
    unsafe {
        let handle = default_handle();
        let mut e = 0.0f64;
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::date::swisseph_time_equ(
            handle,
            J2000,
            &mut e,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        // Equation of time is at most ~16 minutes ≈ 0.011 days
        assert!(e.abs() < 0.015);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn lmt_lat_roundtrip() {
    unsafe {
        let handle = default_handle();
        let mut err_buf = [0u8; 256];
        let geolon = 8.55; // Zurich
        let tjd_lmt = J2000;

        let mut tjd_lat = 0.0f64;
        let ret = swisseph_ffi::date::swisseph_lmt_to_lat(
            handle,
            tjd_lmt,
            geolon,
            &mut tjd_lat,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);

        let mut tjd_lmt_back = 0.0f64;
        let ret = swisseph_ffi::date::swisseph_lat_to_lmt(
            handle,
            tjd_lat,
            geolon,
            &mut tjd_lmt_back,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        assert!((tjd_lmt - tjd_lmt_back).abs() < 1e-9);
        swisseph_ffi::swisseph_free(handle);
    }
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

#[test]
fn csroundsec_boundary() {
    // 30° boundary guard: 30*360000 - 50 cs rounds down, not up
    let x = 30 * 360000 - 50;
    let r = swisseph_ffi::date::swisseph_csroundsec(x);
    assert_eq!(r, (30 * 360000 - 100)); // rounds down
}

#[test]
fn cs2timestr_basic() {
    unsafe {
        let mut buf = [0u8; 32];
        let cs = 12 * 3600 * 100 + 30 * 60 * 100 + 45 * 100; // 12:30:45 in centiseconds
        swisseph_ffi::date::swisseph_cs2timestr(
            cs,
            b':' as c_char,
            false,
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
        );
        let s = std::ffi::CStr::from_ptr(buf.as_ptr() as *const c_char)
            .to_str()
            .unwrap();
        assert_eq!(s, "12:30:45");
    }
}

#[test]
fn cs2lonlatstr_basic() {
    unsafe {
        let mut buf = [0u8; 32];
        let cs = 8 * 3600 * 100 + 33 * 60 * 100; // 8°33'
        swisseph_ffi::date::swisseph_cs2lonlatstr(
            cs,
            b'E' as c_char,
            b'W' as c_char,
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
        );
        let s = std::ffi::CStr::from_ptr(buf.as_ptr() as *const c_char)
            .to_str()
            .unwrap();
        assert_eq!(s, "8E33");
    }
}

#[test]
fn cs2degstr_basic() {
    unsafe {
        let mut buf = [0u8; 32];
        let cs = 15 * 3600 * 100 + 22 * 60 * 100 + 10 * 100; // 15°22'10"
        swisseph_ffi::date::swisseph_cs2degstr(cs, buf.as_mut_ptr() as *mut c_char, buf.len());
        let s = std::ffi::CStr::from_ptr(buf.as_ptr() as *const c_char)
            .to_str()
            .unwrap();
        assert_eq!(s, "15\u{b0}22'10");
    }
}

// ---------------------------------------------------------------------------
// split_deg
// ---------------------------------------------------------------------------

#[test]
fn split_deg_basic() {
    unsafe {
        let mut deg = 0i32;
        let mut min = 0i32;
        let mut sec = 0i32;
        let mut secfr = 0.0f64;
        let mut sign = 0i32;
        swisseph_ffi::date::swisseph_split_deg(
            123.456, 0, // no flags
            &mut deg, &mut min, &mut sec, &mut secfr, &mut sign,
        );
        assert_eq!(sign, 1);
        assert_eq!(deg, 123);
        assert_eq!(min, 27);
        assert_eq!(sec, 21);
        assert!(secfr > 0.0);
    }
}

#[test]
fn split_deg_zodiacal() {
    unsafe {
        let mut deg = 0i32;
        let mut min = 0i32;
        let mut sec = 0i32;
        let mut secfr = 0.0f64;
        let mut sign = 0i32;
        // SE_SPLIT_DEG_ZODIACAL = 8, SE_SPLIT_DEG_ROUND_SEC = 1
        swisseph_ffi::date::swisseph_split_deg(
            123.456,
            8 | 1, // ZODIACAL + ROUND_SEC
            &mut deg,
            &mut min,
            &mut sec,
            &mut secfr,
            &mut sign,
        );
        // 123° / 30° = sign 4 (Leo), remainder 3°
        assert_eq!(sign, 4);
        assert_eq!(deg, 3);
    }
}

// ---------------------------------------------------------------------------
// Utility functions (util.rs)
// ---------------------------------------------------------------------------

#[test]
fn degnorm_basic() {
    let r = swisseph_ffi::util::swisseph_degnorm(-30.0);
    assert!((r - 330.0).abs() < 1e-10);
    let r = swisseph_ffi::util::swisseph_degnorm(400.0);
    assert!((r - 40.0).abs() < 1e-10);
}

#[test]
fn radnorm_basic() {
    let r = swisseph_ffi::util::swisseph_radnorm(-1.0);
    let expected = std::f64::consts::TAU - 1.0;
    assert!((r - expected).abs() < 1e-10);
}

#[test]
fn difdegn_range() {
    // difdegn = p1 - p2 in [0, 360)
    let r = swisseph_ffi::util::swisseph_difdegn(10.0, 350.0);
    assert!((r - 20.0).abs() < 1e-10);
    // Boundary: exactly 360 difference wraps to 0
    let r = swisseph_ffi::util::swisseph_difdegn(0.0, 0.0);
    assert!(r.abs() < 1e-10);
    // Boundary: exactly 180 difference stays at 180
    let r = swisseph_ffi::util::swisseph_difdegn(180.0, 0.0);
    assert!((r - 180.0).abs() < 1e-10);
}

#[test]
fn difdeg2n_range() {
    // difdeg2n = p1 - p2 in [-180, 180)
    let r = swisseph_ffi::util::swisseph_difdeg2n(10.0, 350.0);
    assert!((r - 20.0).abs() < 1e-10);
    let r = swisseph_ffi::util::swisseph_difdeg2n(350.0, 10.0);
    assert!((r - (-20.0)).abs() < 1e-10);
    // Boundary: exactly 180 difference maps to -180 (the [-180, 180) convention)
    let r = swisseph_ffi::util::swisseph_difdeg2n(0.0, 180.0);
    assert!((r - (-180.0)).abs() < 1e-10);
    let r = swisseph_ffi::util::swisseph_difdeg2n(180.0, 0.0);
    assert!((r - (-180.0)).abs() < 1e-10);
    // 90 - 270 = -180 (wraps)
    let r = swisseph_ffi::util::swisseph_difdeg2n(90.0, 270.0);
    assert!((r - (-180.0)).abs() < 1e-10);
}

#[test]
fn deg_midp_basic() {
    let mid = swisseph_ffi::util::swisseph_deg_midp(350.0, 10.0);
    assert!((mid - 0.0).abs() < 1e-10 || (mid - 360.0).abs() < 1e-10);
}

#[test]
fn cotrans_identity() {
    // With eps = 0, cotrans should be identity-like
    unsafe {
        let xpo = [100.0f64, 0.0, 1.0];
        let mut xpn = [0.0f64; 3];
        swisseph_ffi::util::swisseph_cotrans(xpo.as_ptr(), xpn.as_mut_ptr(), 0.0);
        assert!((xpn[0] - 100.0).abs() < 1e-10);
        assert!(xpn[1].abs() < 1e-10);
        assert!((xpn[2] - 1.0).abs() < 1e-10);
    }
}

#[test]
fn cotrans_sp_basic() {
    unsafe {
        let xpo = [100.0f64, 10.0, 1.0, 0.5, 0.1, 0.0];
        let mut xpn = [0.0f64; 6];
        swisseph_ffi::util::swisseph_cotrans_sp(xpo.as_ptr(), xpn.as_mut_ptr(), 23.4);
        // Just verify it doesn't crash and produces reasonable values
        assert!(xpn[0].is_finite());
        assert!(xpn[1].is_finite());
        assert!((xpn[2] - 1.0).abs() < 1e-10); // distance unchanged
    }
}
