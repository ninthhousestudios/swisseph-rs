use std::ffi::{CStr, CString, c_char};
use std::ptr;

use swisseph_ffi::SweEphemeris;
use swisseph_ffi::SweSidMode;
use swisseph_ffi::config::SweConfig;
use swisseph_ffi::error::SweErrorCode;

const J2000: f64 = 2451545.0;
const SEFLG_SPEED: i32 = 256;
const SEFLG_EQUATORIAL: i32 = 2048;
const SEFLG_XYZ: i32 = 4096;
const SEFLG_SIDEREAL: i32 = 65536;

const EPHE_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../ephe\0");

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

unsafe fn handle_with_stars() -> *mut SweEphemeris {
    unsafe {
        let mut config = std::mem::zeroed::<SweConfig>();
        swisseph_ffi::config::swisseph_config_default(&mut config);
        config.ephe_path = EPHE_PATH.as_ptr() as *const c_char;
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0, "swisseph_new with ephe_path failed");
        handle
    }
}

unsafe fn sidereal_handle() -> *mut SweEphemeris {
    unsafe {
        let mut config = std::mem::zeroed::<SweConfig>();
        swisseph_ffi::config::swisseph_config_default(&mut config);
        config.has_sidereal = true;
        config.sid_mode = 1; // Lahiri
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0, "swisseph_new with Lahiri failed");
        handle
    }
}

fn assert_bitwise(ffi: &[f64; 6], rust: &[f64; 6], label: &str) {
    for i in 0..6 {
        assert!(
            (ffi[i] - rust[i]).abs() < 1e-15,
            "{label} xx[{i}] mismatch: ffi={} rust={}",
            ffi[i],
            rust[i]
        );
    }
}

// ---------------------------------------------------------------------------
// swisseph_calc (TT)
// ---------------------------------------------------------------------------

#[test]
fn calc_tt_planets() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let epochs = [J2000, 2460000.5, 2415020.5];
    let bodies: [(i32, swisseph::Body); 10] = [
        (0, swisseph::Body::Sun),
        (1, swisseph::Body::Moon),
        (2, swisseph::Body::Mercury),
        (3, swisseph::Body::Venus),
        (4, swisseph::Body::Mars),
        (5, swisseph::Body::Jupiter),
        (6, swisseph::Body::Saturn),
        (7, swisseph::Body::Uranus),
        (8, swisseph::Body::Neptune),
        (9, swisseph::Body::Pluto),
    ];
    let flag_combos: [(i32, swisseph::CalcFlags); 3] = [
        (SEFLG_SPEED, swisseph::CalcFlags::SPEED),
        (
            SEFLG_SPEED | SEFLG_EQUATORIAL,
            swisseph::CalcFlags::SPEED | swisseph::CalcFlags::EQUATORIAL,
        ),
        (
            SEFLG_SPEED | SEFLG_XYZ,
            swisseph::CalcFlags::SPEED | swisseph::CalcFlags::XYZ,
        ),
    ];

    for &tjd in &epochs {
        for &(ipl, body) in &bodies {
            for &(iflag, flags) in &flag_combos {
                let mut xx = [0.0f64; 6];
                let mut flags_used: i32 = 0;
                let mut err_buf = [0u8; 256];
                let ret = unsafe {
                    swisseph_ffi::swisseph_calc(
                        handle,
                        tjd,
                        ipl,
                        iflag,
                        ptr::null(),
                        ptr::null(),
                        xx.as_mut_ptr(),
                        &mut flags_used,
                        err_buf.as_mut_ptr() as *mut c_char,
                        err_buf.len(),
                    )
                };
                assert_eq!(ret, 0, "calc failed for ipl={ipl} tjd={tjd} iflag={iflag}");

                let result = eph.calc(tjd, body, flags).unwrap();
                let label = format!("calc ipl={ipl} tjd={tjd} iflag={iflag}");
                assert_bitwise(&xx, &result.data, &label);
            }
        }
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// swisseph_calc — TrueNode, MeanNode
// ---------------------------------------------------------------------------

#[test]
fn calc_tt_nodes() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    for &(ipl, body) in &[
        (10i32, swisseph::Body::MeanNode),
        (11, swisseph::Body::TrueNode),
    ] {
        let mut xx = [0.0f64; 6];
        let mut err_buf = [0u8; 256];
        let ret = unsafe {
            swisseph_ffi::swisseph_calc(
                handle,
                J2000,
                ipl,
                SEFLG_SPEED,
                ptr::null(),
                ptr::null(),
                xx.as_mut_ptr(),
                ptr::null_mut(),
                err_buf.as_mut_ptr() as *mut c_char,
                err_buf.len(),
            )
        };
        assert_eq!(ret, 0);
        let result = eph.calc(J2000, body, swisseph::CalcFlags::SPEED).unwrap();
        assert_bitwise(&xx, &result.data, &format!("node ipl={ipl}"));
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// swisseph_calc — fictitious body
// ---------------------------------------------------------------------------

#[test]
fn calc_tt_fictitious() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let ipl = 40; // Cupido
    let body = swisseph::Body::try_from(ipl).unwrap();
    let mut xx = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_calc(
            handle,
            J2000,
            ipl,
            SEFLG_SPEED,
            ptr::null(),
            ptr::null(),
            xx.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);
    let result = eph.calc(J2000, body, swisseph::CalcFlags::SPEED).unwrap();
    assert_bitwise(&xx, &result.data, "fictitious Cupido");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// fixstar2 — Aldebaran
// ---------------------------------------------------------------------------

#[test]
fn fixstar2_aldebaran() {
    let ephe_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../ephe");
    let mut config = swisseph::EphemerisConfig::default();
    config.ephe_path = Some(ephe_path);
    let eph = swisseph::Ephemeris::new(config).unwrap();
    let handle = unsafe { handle_with_stars() };

    let star = CString::new("Aldebaran").unwrap();
    let mut star_out = [0u8; 128];
    let mut xx = [0.0f64; 6];
    let mut flags_used: i32 = 0;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_fixstar2(
            handle,
            star.as_ptr(),
            star_out.as_mut_ptr() as *mut c_char,
            star_out.len(),
            J2000,
            SEFLG_SPEED,
            ptr::null(),
            ptr::null(),
            xx.as_mut_ptr(),
            &mut flags_used,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "fixstar2 Aldebaran failed");

    let (name, result) = eph
        .fixstar2("Aldebaran", J2000, swisseph::CalcFlags::SPEED)
        .unwrap();

    // Check resolved name was written back
    let star_out_str = unsafe { CStr::from_ptr(star_out.as_ptr() as *const c_char) }
        .to_str()
        .unwrap();
    assert_eq!(star_out_str, name);
    assert!(
        star_out_str.contains(","),
        "resolved name should contain comma: {star_out_str}"
    );

    assert_bitwise(&xx, &result.data, "fixstar2 Aldebaran");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// fixstar2_mag
// ---------------------------------------------------------------------------

#[test]
fn fixstar2_mag_aldebaran() {
    let ephe_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../ephe");
    let mut config = swisseph::EphemerisConfig::default();
    config.ephe_path = Some(ephe_path);
    let eph = swisseph::Ephemeris::new(config).unwrap();
    let handle = unsafe { handle_with_stars() };

    let star = CString::new("Aldebaran").unwrap();
    let mut star_out = [0u8; 128];
    let mut mag: f64 = 0.0;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_fixstar2_mag(
            handle,
            star.as_ptr(),
            star_out.as_mut_ptr() as *mut c_char,
            star_out.len(),
            &mut mag,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "fixstar2_mag Aldebaran failed");

    let (name, expected_mag) = eph.fixstar2_mag("Aldebaran").unwrap();
    assert!((mag - expected_mag).abs() < 1e-15, "mag mismatch");

    let star_out_str = unsafe { CStr::from_ptr(star_out.as_ptr() as *const c_char) }
        .to_str()
        .unwrap();
    assert_eq!(star_out_str, name);

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// fixstar2 — unknown star
// ---------------------------------------------------------------------------

#[test]
fn fixstar2_unknown_star() {
    let handle = unsafe { handle_with_stars() };

    let star = CString::new("NonexistentStar12345").unwrap();
    let mut xx = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_fixstar2(
            handle,
            star.as_ptr(),
            ptr::null_mut(),
            0,
            J2000,
            SEFLG_SPEED,
            ptr::null(),
            ptr::null(),
            xx.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert!(ret < 0, "expected error for unknown star, got {ret}");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// ayanamsa — all four entry points with per-call Lahiri override
// ---------------------------------------------------------------------------

#[test]
fn ayanamsa_ex_per_call_lahiri() {
    let handle = unsafe { default_handle() };

    let sid = SweSidMode {
        sid_mode: 1, // Lahiri
        t0: 0.0,
        ayan_t0: 0.0,
    };
    let mut daya: f64 = 0.0;
    let mut flags_used: i32 = 0;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_ex(
            handle,
            J2000,
            0, // no flags
            &sid,
            &mut daya,
            &mut flags_used,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "get_ayanamsa_ex failed");

    // Verify against Rust API
    let mut config = swisseph::EphemerisConfig::default();
    config.set_sidereal_mode(1, 0.0, 0.0);
    let eph = swisseph::Ephemeris::new(config).unwrap();
    let expected = eph
        .get_ayanamsa_ex(J2000, swisseph::CalcFlags::empty())
        .unwrap();
    assert!(
        (daya - expected).abs() < 1e-15,
        "ayanamsa_ex mismatch: ffi={daya} rust={expected}"
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn ayanamsa_ex_ut_per_call_lahiri() {
    let handle = unsafe { default_handle() };

    let sid = SweSidMode {
        sid_mode: 1,
        t0: 0.0,
        ayan_t0: 0.0,
    };
    let mut daya: f64 = 0.0;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_ex_ut(
            handle,
            J2000,
            0,
            &sid,
            &mut daya,
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "get_ayanamsa_ex_ut failed");

    let mut config = swisseph::EphemerisConfig::default();
    config.set_sidereal_mode(1, 0.0, 0.0);
    let eph = swisseph::Ephemeris::new(config).unwrap();
    let expected = eph
        .get_ayanamsa_ut(J2000, swisseph::CalcFlags::empty())
        .unwrap();
    assert!(
        (daya - expected).abs() < 1e-15,
        "ayanamsa_ex_ut mismatch: ffi={daya} rust={expected}"
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn ayanamsa_plain_per_call_lahiri() {
    let handle = unsafe { default_handle() };

    let sid = SweSidMode {
        sid_mode: 1,
        t0: 0.0,
        ayan_t0: 0.0,
    };
    let daya = unsafe { swisseph_ffi::swisseph_get_ayanamsa(handle, J2000, &sid) };
    assert!(!daya.is_nan(), "plain ayanamsa returned NAN");

    let mut config = swisseph::EphemerisConfig::default();
    config.set_sidereal_mode(1, 0.0, 0.0);
    let eph = swisseph::Ephemeris::new(config).unwrap();
    let expected = eph.get_ayanamsa(J2000).unwrap();
    assert!(
        (daya - expected).abs() < 1e-15,
        "plain ayanamsa mismatch: ffi={daya} rust={expected}"
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn ayanamsa_ut_per_call_lahiri() {
    let handle = unsafe { default_handle() };

    let sid = SweSidMode {
        sid_mode: 1,
        t0: 0.0,
        ayan_t0: 0.0,
    };
    let daya = unsafe { swisseph_ffi::swisseph_get_ayanamsa_ut(handle, J2000, &sid) };
    assert!(!daya.is_nan(), "ayanamsa_ut returned NAN");
    assert!(
        daya > 20.0 && daya < 30.0,
        "ayanamsa_ut out of range: {daya}"
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// ayanamsa — no sidereal mode configured, no per-call override → error
// ---------------------------------------------------------------------------

#[test]
fn ayanamsa_no_configured_mode_succeeds() {
    let handle = unsafe { default_handle() };

    let mut daya: f64 = 0.0;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_ex(
            handle,
            J2000,
            0,
            ptr::null(), // no override, no configured mode — uses mode 0 fallback
            &mut daya,
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "ayanamsa should not error without configured mode");
    assert!(!daya.is_nan(), "ayanamsa should not be NaN");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// ayanamsa — handle-configured sidereal mode (no per-call override)
// ---------------------------------------------------------------------------

#[test]
fn ayanamsa_ex_handle_configured() {
    let handle = unsafe { sidereal_handle() };

    let mut daya: f64 = 0.0;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_ex(
            handle,
            J2000,
            0,
            ptr::null(), // use handle's Lahiri
            &mut daya,
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "handle-configured ayanamsa_ex failed");
    assert!(daya > 20.0 && daya < 30.0, "ayanamsa out of range: {daya}");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// get_ayanamsa_name
// ---------------------------------------------------------------------------

#[test]
fn ayanamsa_name_lahiri() {
    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_name(
            1, // Lahiri
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);
    let name = unsafe { CStr::from_ptr(buf.as_ptr() as *const c_char) }
        .to_str()
        .unwrap();
    assert_eq!(name, "Lahiri");
}

#[test]
fn ayanamsa_name_user_defined() {
    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_name(
            255, // User
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);
    let name = unsafe { CStr::from_ptr(buf.as_ptr() as *const c_char) }
        .to_str()
        .unwrap();
    assert_eq!(name, "");
}

#[test]
fn ayanamsa_name_invalid() {
    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_ayanamsa_name(
            999,
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, SweErrorCode::InvalidSiderealMode as i32);
}

// ---------------------------------------------------------------------------
// get_planet_name
// ---------------------------------------------------------------------------

#[test]
fn planet_name_sun() {
    let handle = unsafe { default_handle() };

    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_planet_name(
            handle,
            0, // Sun
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);
    let name = unsafe { CStr::from_ptr(buf.as_ptr() as *const c_char) }
        .to_str()
        .unwrap();
    assert_eq!(name, "Sun");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn planet_name_fictitious() {
    let handle = unsafe { default_handle() };

    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_planet_name(
            handle,
            40, // Cupido
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);
    let name = unsafe { CStr::from_ptr(buf.as_ptr() as *const c_char) }
        .to_str()
        .unwrap();
    assert!(!name.is_empty(), "fictitious body name should not be empty");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn planet_name_invalid_body() {
    let handle = unsafe { default_handle() };

    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_get_planet_name(
            handle,
            -999,
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, SweErrorCode::InvalidBody as i32);

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// Per-call sidereal override on calc
// ---------------------------------------------------------------------------

#[test]
fn calc_ut_sidereal_override() {
    let handle = unsafe { default_handle() };

    let sid = SweSidMode {
        sid_mode: 1, // Lahiri
        t0: 0.0,
        ayan_t0: 0.0,
    };
    let mut xx = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_calc_ut(
            handle,
            J2000,
            0, // Sun
            SEFLG_SPEED | SEFLG_SIDEREAL,
            ptr::null(),
            &sid,
            xx.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "calc_ut with sidereal override failed");

    // Verify against Rust API with configured sidereal mode
    let mut config = swisseph::EphemerisConfig::default();
    config.set_sidereal_mode(1, 0.0, 0.0);
    let eph = swisseph::Ephemeris::new(config).unwrap();
    let result = eph
        .calc_ut(
            J2000,
            swisseph::Body::Sun,
            swisseph::CalcFlags::SPEED | swisseph::CalcFlags::SIDEREAL,
        )
        .unwrap();
    assert_bitwise(&xx, &result.data, "calc_ut sidereal Lahiri Sun");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}
