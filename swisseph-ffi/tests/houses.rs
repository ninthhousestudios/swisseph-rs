use std::ffi::c_char;
use std::ptr;

use swisseph_ffi::SweEphemeris;
use swisseph_ffi::config::SweConfig;
use swisseph_ffi::error::SweErrorCode;

const J2000: f64 = 2451545.0;
const SEFLG_SPEED: i32 = 256;

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

fn assert_close(a: f64, b: f64, eps: f64, label: &str) {
    assert!(
        (a - b).abs() < eps,
        "{label}: {a} vs {b}, diff={}",
        (a - b).abs()
    );
}

// ---------------------------------------------------------------------------
// houses — basic Placidus
// ---------------------------------------------------------------------------

#[test]
fn houses_placidus() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geolat = 47.37;
    let geolon = 8.55;
    let hsys = b'P' as i32;

    let mut cusps = [0.0f64; 13];
    let mut ascmc = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses(
            handle,
            J2000,
            geolat,
            geolon,
            hsys,
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "swisseph_houses failed");

    let result = eph
        .houses(J2000, geolat, geolon, swisseph::HouseSystem::Placidus)
        .unwrap();

    for i in 1..=12 {
        assert_close(cusps[i], result.cusps[i], 1e-15, &format!("cusp[{i}]"));
    }

    let ascmc_rust = result.ascmc.as_array();
    for i in 0..8 {
        assert_close(ascmc[i], ascmc_rust[i], 1e-15, &format!("ascmc[{i}]"));
    }
    assert_eq!(ascmc[8], 0.0);
    assert_eq!(ascmc[9], 0.0);

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// houses_ex2 — Placidus with speeds
// ---------------------------------------------------------------------------

#[test]
fn houses_ex2_placidus_with_speeds() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geolat = 47.37;
    let geolon = 8.55;
    let hsys = b'P' as i32;

    let mut cusps = [0.0f64; 13];
    let mut ascmc = [0.0f64; 10];
    let mut cusp_speed = [0.0f64; 13];
    let mut ascmc_speed = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses_ex2(
            handle,
            J2000,
            SEFLG_SPEED,
            geolat,
            geolon,
            hsys,
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            cusp_speed.as_mut_ptr(),
            ascmc_speed.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "swisseph_houses_ex2 failed");

    let result = eph
        .houses_ex2(
            J2000,
            swisseph::CalcFlags::SPEED,
            geolat,
            geolon,
            swisseph::HouseSystem::Placidus,
        )
        .unwrap();

    for i in 1..=12 {
        assert_close(cusps[i], result.cusps[i], 1e-15, &format!("cusp[{i}]"));
        assert_close(
            cusp_speed[i],
            result.cusp_speeds[i],
            1e-15,
            &format!("cusp_speed[{i}]"),
        );
    }

    let ascmc_rust = result.ascmc.as_array();
    let ascmc_speed_rust = result.ascmc_speeds.as_array();
    for i in 0..8 {
        assert_close(ascmc[i], ascmc_rust[i], 1e-15, &format!("ascmc[{i}]"));
        assert_close(
            ascmc_speed[i],
            ascmc_speed_rust[i],
            1e-15,
            &format!("ascmc_speed[{i}]"),
        );
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// houses_ex2 — Gauquelin with 37 slots
// ---------------------------------------------------------------------------

#[test]
fn houses_ex2_gauquelin() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geolat = 48.85;
    let geolon = 2.35;
    let hsys = b'G' as i32;

    let mut cusps = [0.0f64; 37];
    let mut ascmc = [0.0f64; 10];
    let mut cusp_speed = [0.0f64; 37];
    let mut ascmc_speed = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses_ex2(
            handle,
            J2000,
            0,
            geolat,
            geolon,
            hsys,
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            cusp_speed.as_mut_ptr(),
            ascmc_speed.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "swisseph_houses_ex2 Gauquelin failed");

    let result = eph
        .houses_ex2(
            J2000,
            swisseph::CalcFlags::empty(),
            geolat,
            geolon,
            swisseph::HouseSystem::Gauquelin,
        )
        .unwrap();

    for i in 1..=36 {
        assert_close(cusps[i], result.cusps[i], 1e-15, &format!("G cusp[{i}]"));
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// houses_ex2 — Koch, Whole Sign, Sunshine
// ---------------------------------------------------------------------------

#[test]
fn houses_ex2_multiple_systems() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geolat = 47.37;
    let geolon = 8.55;

    for (hsys_char, hsys_enum) in [
        (b'K', swisseph::HouseSystem::Koch),
        (b'W', swisseph::HouseSystem::WholeSign),
        (b'I', swisseph::HouseSystem::Sunshine),
    ] {
        let mut cusps = [0.0f64; 13];
        let mut ascmc = [0.0f64; 10];
        let mut err_buf = [0u8; 256];

        let ret = unsafe {
            swisseph_ffi::houses::swisseph_houses_ex2(
                handle,
                J2000,
                0,
                geolat,
                geolon,
                hsys_char as i32,
                cusps.as_mut_ptr(),
                ascmc.as_mut_ptr(),
                ptr::null_mut(),
                ptr::null_mut(),
                err_buf.as_mut_ptr() as *mut c_char,
                err_buf.len(),
            )
        };
        assert_eq!(ret, 0, "houses_ex2 failed for '{}'", hsys_char as char);

        let result = eph
            .houses_ex2(
                J2000,
                swisseph::CalcFlags::empty(),
                geolat,
                geolon,
                hsys_enum,
            )
            .unwrap();

        for i in 1..=12 {
            assert_close(
                cusps[i],
                result.cusps[i],
                1e-15,
                &format!("'{}' cusp[{i}]", hsys_char as char),
            );
        }
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// houses_ex2 — polar latitude (Placidus fallback)
// ---------------------------------------------------------------------------

#[test]
fn houses_ex2_polar_latitude() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geolat = 70.0;
    let geolon = 25.0;

    let mut cusps = [0.0f64; 13];
    let mut ascmc = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses_ex2(
            handle,
            J2000,
            0,
            geolat,
            geolon,
            b'P' as i32,
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);

    let result = eph
        .houses_ex2(
            J2000,
            swisseph::CalcFlags::empty(),
            geolat,
            geolon,
            swisseph::HouseSystem::Placidus,
        )
        .unwrap();

    for i in 1..=12 {
        assert_close(
            cusps[i],
            result.cusps[i],
            1e-15,
            &format!("polar cusp[{i}]"),
        );
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// houses_armc — handle-free
// ---------------------------------------------------------------------------

#[test]
fn houses_armc_basic() {
    let armc = 150.0;
    let geolat = 47.37;
    let eps = 23.44;

    let mut cusps = [0.0f64; 13];
    let mut ascmc = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses_armc(
            armc,
            geolat,
            eps,
            b'P' as i32,
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);

    let result =
        swisseph::houses::houses_armc(armc, geolat, eps, swisseph::HouseSystem::Placidus, None)
            .unwrap();

    for i in 1..=12 {
        assert_close(cusps[i], result.cusps[i], 1e-15, &format!("armc cusp[{i}]"));
    }
}

// ---------------------------------------------------------------------------
// houses_armc_ex2 — with speeds and sundec
// ---------------------------------------------------------------------------

#[test]
fn houses_armc_ex2_with_speeds() {
    let armc = 150.0;
    let geolat = 47.37;
    let eps = 23.44;

    let mut cusps = [0.0f64; 13];
    let mut ascmc = [0.0f64; 10];
    let mut cusp_speed = [0.0f64; 13];
    let mut ascmc_speed = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses_armc_ex2(
            armc,
            geolat,
            eps,
            b'P' as i32,
            ptr::null(),
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            cusp_speed.as_mut_ptr(),
            ascmc_speed.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);

    let result =
        swisseph::houses::houses_armc(armc, geolat, eps, swisseph::HouseSystem::Placidus, None)
            .unwrap();

    for i in 1..=12 {
        assert_close(
            cusps[i],
            result.cusps[i],
            1e-15,
            &format!("armc_ex2 cusp[{i}]"),
        );
        assert_close(
            cusp_speed[i],
            result.cusp_speeds[i],
            1e-15,
            &format!("armc_ex2 speed[{i}]"),
        );
    }
}

// ---------------------------------------------------------------------------
// house_pos — Placidus + Gauquelin
// ---------------------------------------------------------------------------

#[test]
fn house_pos_placidus() {
    let armc = 150.0;
    let geolat = 47.37;
    let eps = 23.44;
    let xpin = [120.0f64, 5.0];

    let mut hpos: f64 = 0.0;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_house_pos(
            armc,
            geolat,
            eps,
            b'P' as i32,
            xpin.as_ptr(),
            ptr::null(),
            &mut hpos,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);

    let rust_hpos = swisseph::houses::house_pos(
        armc,
        geolat,
        eps,
        swisseph::HouseSystem::Placidus,
        xpin,
        None,
    )
    .unwrap();

    assert_close(hpos, rust_hpos, 1e-15, "house_pos Placidus");
}

#[test]
fn house_pos_gauquelin() {
    let armc = 150.0;
    let geolat = 47.37;
    let eps = 23.44;
    let xpin = [120.0f64, 5.0];

    let mut hpos: f64 = 0.0;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_house_pos(
            armc,
            geolat,
            eps,
            b'G' as i32,
            xpin.as_ptr(),
            ptr::null(),
            &mut hpos,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);

    let rust_hpos = swisseph::houses::house_pos(
        armc,
        geolat,
        eps,
        swisseph::HouseSystem::Gauquelin,
        xpin,
        None,
    )
    .unwrap();

    assert_close(hpos, rust_hpos, 1e-15, "house_pos Gauquelin");
}

#[test]
fn house_pos_sunshine_with_sundec() {
    let armc = 150.0;
    let geolat = 47.37;
    let eps = 23.44;
    let xpin = [120.0f64, 5.0];
    let sundec: f64 = 12.5;

    let mut hpos: f64 = 0.0;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_house_pos(
            armc,
            geolat,
            eps,
            b'I' as i32,
            xpin.as_ptr(),
            &sundec,
            &mut hpos,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0);

    let rust_hpos = swisseph::houses::house_pos(
        armc,
        geolat,
        eps,
        swisseph::HouseSystem::Sunshine,
        xpin,
        Some(sundec),
    )
    .unwrap();

    assert_close(hpos, rust_hpos, 1e-15, "house_pos Sunshine");
}

// ---------------------------------------------------------------------------
// house_name
// ---------------------------------------------------------------------------

#[test]
fn house_name_all_systems() {
    let systems: [(i32, &str); 5] = [
        (b'P' as i32, "Placidus"),
        (b'K' as i32, "Koch"),
        (b'W' as i32, "equal/ whole sign"),
        (b'G' as i32, "Gauquelin sectors"),
        (b'I' as i32, "Sunshine"),
    ];

    for (hsys, expected) in &systems {
        let mut buf = [0u8; 64];
        let mut err_buf = [0u8; 256];

        let ret = unsafe {
            swisseph_ffi::houses::swisseph_house_name(
                *hsys,
                buf.as_mut_ptr() as *mut c_char,
                buf.len(),
                err_buf.as_mut_ptr() as *mut c_char,
                err_buf.len(),
            )
        };
        assert_eq!(ret, 0, "house_name failed for '{}'", *hsys as u8 as char);

        let name = unsafe { std::ffi::CStr::from_ptr(buf.as_ptr() as *const c_char) }
            .to_str()
            .unwrap();
        assert_eq!(
            name, *expected,
            "house_name mismatch for '{}'",
            *hsys as u8 as char
        );
    }
}

#[test]
fn house_name_invalid() {
    let mut buf = [0u8; 64];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_house_name(
            b'Z' as i32,
            buf.as_mut_ptr() as *mut c_char,
            buf.len(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, SweErrorCode::InvalidHouseSystem as i32);
}

// ---------------------------------------------------------------------------
// azalt / azalt_rev round-trip
// ---------------------------------------------------------------------------

#[test]
fn azalt_round_trip() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geopos = [8.55f64, 47.37, 500.0];
    let xin_ecl = [120.0f64, 5.0];

    // Forward: ecliptic -> horizontal
    let mut xaz = [0.0f64; 3];
    unsafe {
        swisseph_ffi::houses::swisseph_azalt(
            handle,
            J2000,
            0, // SE_ECL2HOR
            geopos.as_ptr(),
            1013.25,
            15.0,
            xin_ecl.as_ptr(),
            xaz.as_mut_ptr(),
        );
    }

    let rust_xaz = eph.azalt(
        J2000,
        swisseph::azalt::AzAltDir::EclToHor,
        geopos,
        1013.25,
        15.0,
        0.0,
        xin_ecl,
    );

    for i in 0..3 {
        assert_close(xaz[i], rust_xaz[i], 1e-15, &format!("azalt[{i}]"));
    }

    // Reverse: horizontal -> ecliptic
    let xin_hor = [xaz[0], xaz[1]]; // azimuth, true altitude
    let mut xout = [0.0f64; 2];
    unsafe {
        swisseph_ffi::houses::swisseph_azalt_rev(
            handle,
            J2000,
            0, // SE_HOR2ECL
            geopos.as_ptr(),
            xin_hor.as_ptr(),
            xout.as_mut_ptr(),
        );
    }

    let rust_xout = eph.azalt_rev(J2000, swisseph::azalt::HorDir::HorToEcl, geopos, xin_hor);

    for i in 0..2 {
        assert_close(xout[i], rust_xout[i], 1e-15, &format!("azalt_rev[{i}]"));
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn azalt_equ2hor() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geopos = [8.55f64, 47.37, 500.0];
    let xin_equ = [45.0f64, 20.0];

    let mut xaz = [0.0f64; 3];
    unsafe {
        swisseph_ffi::houses::swisseph_azalt(
            handle,
            J2000,
            1, // SE_EQU2HOR
            geopos.as_ptr(),
            0.0,
            0.0,
            xin_equ.as_ptr(),
            xaz.as_mut_ptr(),
        );
    }

    let rust_xaz = eph.azalt(
        J2000,
        swisseph::azalt::AzAltDir::EquToHor,
        geopos,
        0.0,
        0.0,
        0.0,
        xin_equ,
    );

    for i in 0..3 {
        assert_close(xaz[i], rust_xaz[i], 1e-15, &format!("equ2hor[{i}]"));
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// refrac — both directions
// ---------------------------------------------------------------------------

#[test]
fn refrac_true_to_app() {
    let result = swisseph_ffi::houses::swisseph_refrac(30.0, 1013.25, 15.0, 0);
    let rust_result =
        swisseph::azalt::refrac(30.0, 1013.25, 15.0, swisseph::azalt::RefracDir::TrueToApp);
    assert_close(result, rust_result, 1e-15, "refrac true->app");
}

#[test]
fn refrac_app_to_true() {
    let result = swisseph_ffi::houses::swisseph_refrac(30.0, 1013.25, 15.0, 1);
    let rust_result =
        swisseph::azalt::refrac(30.0, 1013.25, 15.0, swisseph::azalt::RefracDir::AppToTrue);
    assert_close(result, rust_result, 1e-15, "refrac app->true");
}

// ---------------------------------------------------------------------------
// refrac_extended — both directions
// ---------------------------------------------------------------------------

#[test]
fn refrac_extended_true_to_app() {
    let mut dret = [0.0f64; 4];

    let result = unsafe {
        swisseph_ffi::houses::swisseph_refrac_extended(
            30.0,
            500.0,
            1013.25,
            15.0,
            0.0065,
            0,
            dret.as_mut_ptr(),
        )
    };

    let mut rust_dret = [0.0f64; 4];
    let rust_result = swisseph::azalt::refrac_extended(
        30.0,
        500.0,
        1013.25,
        15.0,
        0.0065,
        swisseph::azalt::RefracDir::TrueToApp,
        &mut rust_dret,
    );

    assert_close(result, rust_result, 1e-15, "refrac_ext return");
    for i in 0..4 {
        assert_close(
            dret[i],
            rust_dret[i],
            1e-15,
            &format!("refrac_ext dret[{i}]"),
        );
    }
}

#[test]
fn refrac_extended_app_to_true() {
    let mut dret = [0.0f64; 4];

    let result = unsafe {
        swisseph_ffi::houses::swisseph_refrac_extended(
            30.0,
            500.0,
            1013.25,
            15.0,
            0.0065,
            1,
            dret.as_mut_ptr(),
        )
    };

    let mut rust_dret = [0.0f64; 4];
    let rust_result = swisseph::azalt::refrac_extended(
        30.0,
        500.0,
        1013.25,
        15.0,
        0.0065,
        swisseph::azalt::RefracDir::AppToTrue,
        &mut rust_dret,
    );

    assert_close(result, rust_result, 1e-15, "refrac_ext_rev return");
    for i in 0..4 {
        assert_close(
            dret[i],
            rust_dret[i],
            1e-15,
            &format!("refrac_ext_rev dret[{i}]"),
        );
    }
}

// ---------------------------------------------------------------------------
// gauquelin_sector — imeth 0 and 2
// ---------------------------------------------------------------------------

#[test]
fn gauquelin_sector_geometric() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geopos = [8.55f64, 47.37, 500.0];
    let mut dgsect: f64 = 0.0;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_gauquelin_sector(
            handle,
            J2000,
            0, // Sun
            ptr::null(),
            SEFLG_SPEED,
            0, // imeth=0 geometric
            geopos.as_ptr(),
            0.0,
            0.0,
            &mut dgsect,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "gauquelin_sector imeth=0 failed");

    let rust_sector = eph
        .gauquelin_sector(
            J2000,
            swisseph::Body::Sun,
            None,
            swisseph::CalcFlags::SPEED,
            0,
            geopos,
            0.0,
            0.0,
        )
        .unwrap();

    assert_close(dgsect, rust_sector, 1e-15, "gauquelin imeth=0");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn gauquelin_sector_risetrans() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let geopos = [8.55f64, 47.37, 500.0];
    let mut dgsect: f64 = 0.0;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_gauquelin_sector(
            handle,
            J2000,
            0, // Sun
            ptr::null(),
            SEFLG_SPEED,
            2, // imeth=2 rise/set
            geopos.as_ptr(),
            1013.25,
            15.0,
            &mut dgsect,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, 0, "gauquelin_sector imeth=2 failed");

    let rust_sector = eph
        .gauquelin_sector(
            J2000,
            swisseph::Body::Sun,
            None,
            swisseph::CalcFlags::SPEED,
            2,
            geopos,
            1013.25,
            15.0,
        )
        .unwrap();

    assert_close(dgsect, rust_sector, 1e-15, "gauquelin imeth=2");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// Error paths
// ---------------------------------------------------------------------------

#[test]
fn houses_invalid_hsys() {
    let handle = unsafe { default_handle() };

    let mut cusps = [0.0f64; 13];
    let mut ascmc = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_houses(
            handle,
            J2000,
            47.37,
            8.55,
            b'Z' as i32,
            cusps.as_mut_ptr(),
            ascmc.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(ret, SweErrorCode::InvalidHouseSystem as i32);

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn gauquelin_circumpolar() {
    let handle = unsafe { default_handle() };

    let geopos = [25.0f64, 70.0, 0.0];
    let mut dgsect: f64 = 0.0;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::houses::swisseph_gauquelin_sector(
            handle,
            J2000,
            1, // Moon
            ptr::null(),
            0,
            2, // imeth=2 rise/set
            geopos.as_ptr(),
            1013.25,
            15.0,
            &mut dgsect,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    // At polar latitudes, the Moon may be circumpolar
    // This should return an error, not crash
    if ret != 0 {
        assert_eq!(ret, SweErrorCode::CircumpolarBody as i32);
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// Two epochs, two latitudes — houses_ex2 matrix
// ---------------------------------------------------------------------------

#[test]
fn houses_ex2_matrix() {
    let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
    let handle = unsafe { default_handle() };

    let epochs = [J2000, 2460000.5];
    let locations = [(47.37, 8.55), (70.0, 25.0)];
    let systems: [(i32, swisseph::HouseSystem); 5] = [
        (b'P' as i32, swisseph::HouseSystem::Placidus),
        (b'K' as i32, swisseph::HouseSystem::Koch),
        (b'W' as i32, swisseph::HouseSystem::WholeSign),
        (b'G' as i32, swisseph::HouseSystem::Gauquelin),
        (b'I' as i32, swisseph::HouseSystem::Sunshine),
    ];

    for &tjd in &epochs {
        for &(geolat, geolon) in &locations {
            for &(hsys_char, hsys_enum) in &systems {
                let n = if hsys_enum == swisseph::HouseSystem::Gauquelin {
                    37
                } else {
                    13
                };
                let mut cusps = vec![0.0f64; n];
                let mut ascmc = [0.0f64; 10];
                let mut cusp_speed = vec![0.0f64; n];
                let mut ascmc_speed = [0.0f64; 10];
                let mut err_buf = [0u8; 256];

                let ret = unsafe {
                    swisseph_ffi::houses::swisseph_houses_ex2(
                        handle,
                        tjd,
                        0,
                        geolat,
                        geolon,
                        hsys_char,
                        cusps.as_mut_ptr(),
                        ascmc.as_mut_ptr(),
                        cusp_speed.as_mut_ptr(),
                        ascmc_speed.as_mut_ptr(),
                        err_buf.as_mut_ptr() as *mut c_char,
                        err_buf.len(),
                    )
                };
                let label = format!("hsys='{}' tjd={tjd} lat={geolat}", hsys_char as u8 as char);
                assert_eq!(ret, 0, "houses_ex2 failed: {label}");

                let result = eph
                    .houses_ex2(tjd, swisseph::CalcFlags::empty(), geolat, geolon, hsys_enum)
                    .unwrap();

                let cusp_count = if hsys_enum == swisseph::HouseSystem::Gauquelin {
                    36
                } else {
                    12
                };
                for i in 1..=cusp_count {
                    assert_close(
                        cusps[i],
                        result.cusps[i],
                        1e-15,
                        &format!("{label} cusp[{i}]"),
                    );
                    assert_close(
                        cusp_speed[i],
                        result.cusp_speeds[i],
                        1e-15,
                        &format!("{label} cusp_speed[{i}]"),
                    );
                }
            }
        }
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}
