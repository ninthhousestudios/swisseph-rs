use std::ffi::c_char;
use std::ptr;

use swisseph::Ephemeris;
use swisseph::config::EphemerisConfig;
use swisseph::flags::CalcFlags;
use swisseph::nodaps::NodApsMethod;
use swisseph::types::Body;

use swisseph_ffi::SweEphemeris;
use swisseph_ffi::config::SweConfig;

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

fn make_eph() -> Ephemeris {
    Ephemeris::new(EphemerisConfig::default()).unwrap()
}

fn assert_eps(a: f64, b: f64, eps: f64, label: &str) {
    let diff = (a - b).abs();
    assert!(diff < eps, "{label}: ffi={a} lib={b} diff={diff} eps={eps}");
}

// ---------------------------------------------------------------------------
// pheno / pheno_ut
// ---------------------------------------------------------------------------

#[test]
fn pheno_venus() {
    let eph = make_eph();
    let tjd_et = 2451545.0;
    let body = Body::Venus;

    let (lib_p, lib_flags) = eph.pheno(tjd_et, body, CalcFlags::SPEED).unwrap();

    let handle = unsafe { default_handle() };
    let mut attr = [0.0f64; 20];
    let mut flags_used = 0i32;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_pheno(
            handle,
            tjd_et,
            body.to_raw_id(),
            CalcFlags::SPEED.bits() as i32,
            ptr::null(),
            ptr::null(),
            attr.as_mut_ptr(),
            &mut flags_used,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(attr[0], lib_p.phase_angle, 1e-15, "phase_angle");
    assert_eps(attr[1], lib_p.phase, 1e-15, "phase");
    assert_eps(attr[2], lib_p.elongation, 1e-15, "elongation");
    assert_eps(attr[3], lib_p.apparent_diameter, 1e-15, "apparent_diameter");
    assert_eps(
        attr[4],
        lib_p.apparent_magnitude,
        1e-15,
        "apparent_magnitude",
    );
    assert_eps(
        attr[5],
        lib_p.horizontal_parallax,
        1e-15,
        "horizontal_parallax",
    );
    assert_eq!(flags_used, lib_flags.bits() as i32);
    for i in 6..20 {
        assert_eq!(attr[i], 0.0, "attr[{i}] should be zero");
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn pheno_moon() {
    let eph = make_eph();
    let tjd_et = 2460600.5;
    let body = Body::Moon;

    let (lib_p, _) = eph.pheno(tjd_et, body, CalcFlags::SPEED).unwrap();

    let handle = unsafe { default_handle() };
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_pheno(
            handle,
            tjd_et,
            body.to_raw_id(),
            CalcFlags::SPEED.bits() as i32,
            ptr::null(),
            ptr::null(),
            attr.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(attr[0], lib_p.phase_angle, 1e-15, "phase_angle");
    assert!(
        attr[5] > 0.0,
        "Moon should have nonzero horizontal_parallax"
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn pheno_ut_venus() {
    let eph = make_eph();
    let tjd_ut = 2451545.0;
    let body = Body::Venus;

    let (lib_p, _) = eph.pheno_ut(tjd_ut, body, CalcFlags::SPEED).unwrap();

    let handle = unsafe { default_handle() };
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_pheno_ut(
            handle,
            tjd_ut,
            body.to_raw_id(),
            CalcFlags::SPEED.bits() as i32,
            ptr::null(),
            ptr::null(),
            attr.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(attr[0], lib_p.phase_angle, 1e-15, "phase_angle_ut");
    assert_eps(attr[2], lib_p.elongation, 1e-15, "elongation_ut");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn pheno_moon_topocentric_geopos_override() {
    let geopos = [8.55, 47.37, 500.0];
    let tjd_et = 2451545.0;
    let body = Body::Moon;
    let flags = CalcFlags::SPEED | CalcFlags::TOPOCTR;

    let mut config = EphemerisConfig::default();
    config.topographic = Some(swisseph::config::TopoPosition {
        longitude: geopos[0],
        latitude: geopos[1],
        altitude: geopos[2],
    });
    let eph = Ephemeris::new(config).unwrap();
    let (lib_p, _) = eph.pheno(tjd_et, body, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_pheno(
            handle,
            tjd_et,
            body.to_raw_id(),
            flags.bits() as i32,
            geopos.as_ptr(),
            ptr::null(),
            attr.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(attr[5], lib_p.horizontal_parallax, 1e-15, "topo_hpar");
    assert!(
        attr[5] > 0.0,
        "topocentric Moon should have nonzero horizontal_parallax"
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// nod_aps / nod_aps_ut
// ---------------------------------------------------------------------------

#[test]
fn nod_aps_mean_moon() {
    let eph = make_eph();
    let tjd_et = 2451545.0;
    let body = Body::Moon;
    let flags = CalcFlags::SPEED;

    let lib_na = eph
        .nod_aps(tjd_et, body, flags, NodApsMethod::MEAN)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut xnasc = [0.0f64; 6];
    let mut xndsc = [0.0f64; 6];
    let mut xperi = [0.0f64; 6];
    let mut xaphe = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_nod_aps(
            handle,
            tjd_et,
            body.to_raw_id(),
            flags.bits() as i32,
            NodApsMethod::MEAN.bits() as i32,
            xnasc.as_mut_ptr(),
            xndsc.as_mut_ptr(),
            xperi.as_mut_ptr(),
            xaphe.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    for i in 0..6 {
        assert_eps(xnasc[i], lib_na.ascending[i], 1e-15, &format!("asc[{i}]"));
        assert_eps(xndsc[i], lib_na.descending[i], 1e-15, &format!("desc[{i}]"));
        assert_eps(xperi[i], lib_na.perihelion[i], 1e-15, &format!("peri[{i}]"));
        assert_eps(xaphe[i], lib_na.aphelion[i], 1e-15, &format!("aphe[{i}]"));
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn nod_aps_oscu_mars() {
    let eph = make_eph();
    let tjd_et = 2451545.0;
    let body = Body::Mars;
    let flags = CalcFlags::SPEED;

    let lib_na = eph
        .nod_aps(tjd_et, body, flags, NodApsMethod::OSCU)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut xnasc = [0.0f64; 6];
    let mut xndsc = [0.0f64; 6];
    let mut xperi = [0.0f64; 6];
    let mut xaphe = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_nod_aps(
            handle,
            tjd_et,
            body.to_raw_id(),
            flags.bits() as i32,
            NodApsMethod::OSCU.bits() as i32,
            xnasc.as_mut_ptr(),
            xndsc.as_mut_ptr(),
            xperi.as_mut_ptr(),
            xaphe.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    for i in 0..6 {
        assert_eps(
            xnasc[i],
            lib_na.ascending[i],
            1e-15,
            &format!("osc_asc[{i}]"),
        );
        assert_eps(
            xperi[i],
            lib_na.perihelion[i],
            1e-15,
            &format!("osc_peri[{i}]"),
        );
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn nod_aps_ut_mean_moon() {
    let eph = make_eph();
    let tjd_ut = 2451545.0;
    let body = Body::Moon;
    let flags = CalcFlags::SPEED;

    let lib_na = eph
        .nod_aps_ut(tjd_ut, body, flags, NodApsMethod::MEAN)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut xnasc = [0.0f64; 6];
    let mut xndsc = [0.0f64; 6];
    let mut xperi = [0.0f64; 6];
    let mut xaphe = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_nod_aps_ut(
            handle,
            tjd_ut,
            body.to_raw_id(),
            flags.bits() as i32,
            NodApsMethod::MEAN.bits() as i32,
            xnasc.as_mut_ptr(),
            xndsc.as_mut_ptr(),
            xperi.as_mut_ptr(),
            xaphe.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    for i in 0..6 {
        assert_eps(
            xnasc[i],
            lib_na.ascending[i],
            1e-15,
            &format!("ut_asc[{i}]"),
        );
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// orbital elements
// ---------------------------------------------------------------------------

#[test]
fn orbital_elements_jupiter() {
    let eph = make_eph();
    let tjd_et = 2451545.0;
    let body = Body::Jupiter;
    let flags = CalcFlags::SPEED;

    let lib_e = eph.get_orbital_elements(tjd_et, body, flags).unwrap();
    let lib_arr = lib_e.as_array();

    let handle = unsafe { default_handle() };
    let mut dret = [0.0f64; 50];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_get_orbital_elements(
            handle,
            tjd_et,
            body.to_raw_id(),
            flags.bits() as i32,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    for i in 0..17 {
        assert_eps(dret[i], lib_arr[i], 1e-15, &format!("dret[{i}]"));
    }
    for i in 17..50 {
        assert_eq!(dret[i], 0.0, "dret[{i}] should be zero");
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// orbit_max_min_true_distance
// ---------------------------------------------------------------------------

#[test]
fn orbit_max_min_true_distance_mars() {
    let eph = make_eph();
    let tjd_et = 2451545.0;
    let body = Body::Mars;
    let flags = CalcFlags::SPEED;

    let (lib_max, lib_min, lib_true) = eph
        .orbit_max_min_true_distance(tjd_et, body, flags)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut dmax = 0.0f64;
    let mut dmin = 0.0f64;
    let mut dtrue = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_orbit_max_min_true_distance(
            handle,
            tjd_et,
            body.to_raw_id(),
            flags.bits() as i32,
            &mut dmax,
            &mut dmin,
            &mut dtrue,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(dmax, lib_max, 1e-15, "dmax");
    assert_eps(dmin, lib_min, 1e-15, "dmin");
    assert_eps(dtrue, lib_true, 1e-15, "dtrue");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// crossings — solcross
// ---------------------------------------------------------------------------

#[test]
fn solcross_forward() {
    let eph = make_eph();
    let x2cross = 0.0;
    let tjd_et = 2451500.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.solcross(x2cross, tjd_et, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_solcross(
            handle,
            x2cross,
            tjd_et,
            flags.bits() as i32,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "solcross");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn solcross_ut_forward() {
    let eph = make_eph();
    let x2cross = 180.0;
    let tjd_ut = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.solcross_ut(x2cross, tjd_ut, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_solcross_ut(
            handle,
            x2cross,
            tjd_ut,
            flags.bits() as i32,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "solcross_ut");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// crossings — mooncross
// ---------------------------------------------------------------------------

#[test]
fn mooncross_forward() {
    let eph = make_eph();
    let x2cross = 90.0;
    let tjd_et = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.mooncross(x2cross, tjd_et, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_mooncross(
            handle,
            x2cross,
            tjd_et,
            flags.bits() as i32,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "mooncross");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn mooncross_ut_forward() {
    let eph = make_eph();
    let x2cross = 270.0;
    let tjd_ut = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.mooncross_ut(x2cross, tjd_ut, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_mooncross_ut(
            handle,
            x2cross,
            tjd_ut,
            flags.bits() as i32,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "mooncross_ut");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// crossings — mooncross_node
// ---------------------------------------------------------------------------

#[test]
fn mooncross_node_forward() {
    let eph = make_eph();
    let tjd_et = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_mc = eph.mooncross_node(tjd_et, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut xlon = 0.0f64;
    let mut xlat = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_mooncross_node(
            handle,
            tjd_et,
            flags.bits() as i32,
            &mut xlon,
            &mut xlat,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_mc.jd, 1e-15, "mooncross_node jd");
    assert_eps(xlon, lib_mc.longitude, 1e-15, "mooncross_node lon");
    assert_eps(xlat, lib_mc.latitude, 1e-15, "mooncross_node lat");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn mooncross_node_ut_forward() {
    let eph = make_eph();
    let tjd_ut = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_mc = eph.mooncross_node_ut(tjd_ut, flags).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut xlon = 0.0f64;
    let mut xlat = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_mooncross_node_ut(
            handle,
            tjd_ut,
            flags.bits() as i32,
            &mut xlon,
            &mut xlat,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_mc.jd, 1e-15, "mooncross_node_ut jd");
    assert_eps(xlon, lib_mc.longitude, 1e-15, "mooncross_node_ut lon");
    assert_eps(xlat, lib_mc.latitude, 1e-15, "mooncross_node_ut lat");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// crossings — helio_cross
// ---------------------------------------------------------------------------

#[test]
fn helio_cross_forward() {
    let eph = make_eph();
    let body = Body::Mars;
    let x2cross = 120.5;
    let tjd_et = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.helio_cross(body, x2cross, tjd_et, flags, 1).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_helio_cross(
            handle,
            body.to_raw_id(),
            x2cross,
            tjd_et,
            flags.bits() as i32,
            1,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "helio_cross forward");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn helio_cross_backward() {
    let eph = make_eph();
    let body = Body::Jupiter;
    let x2cross = 0.0;
    let tjd_et = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.helio_cross(body, x2cross, tjd_et, flags, -1).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_helio_cross(
            handle,
            body.to_raw_id(),
            x2cross,
            tjd_et,
            flags.bits() as i32,
            -1,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "helio_cross backward");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn helio_cross_ut_forward() {
    let eph = make_eph();
    let body = Body::Mercury;
    let x2cross = 0.0;
    let tjd_ut = 2451545.0;
    let flags = CalcFlags::SPEED;

    let lib_jd = eph.helio_cross_ut(body, x2cross, tjd_ut, flags, 1).unwrap();

    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_helio_cross_ut(
            handle,
            body.to_raw_id(),
            x2cross,
            tjd_ut,
            flags.bits() as i32,
            1,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(jx, lib_jd, 1e-15, "helio_cross_ut");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// Error paths
// ---------------------------------------------------------------------------

#[test]
fn pheno_null_handle() {
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_pheno(
            ptr::null(),
            2451545.0,
            Body::Venus.to_raw_id(),
            0,
            ptr::null(),
            ptr::null(),
            attr.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert!(ret < 0);
}

#[test]
fn nod_aps_null_xnasc() {
    let handle = unsafe { default_handle() };
    let mut xndsc = [0.0f64; 6];
    let mut xperi = [0.0f64; 6];
    let mut xaphe = [0.0f64; 6];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_nod_aps(
            handle,
            2451545.0,
            Body::Moon.to_raw_id(),
            CalcFlags::SPEED.bits() as i32,
            NodApsMethod::MEAN.bits() as i32,
            ptr::null_mut(),
            xndsc.as_mut_ptr(),
            xperi.as_mut_ptr(),
            xaphe.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert!(ret < 0);
    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn helio_cross_sun_rejected() {
    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_helio_cross(
            handle,
            Body::Sun.to_raw_id(),
            0.0,
            2451545.0,
            CalcFlags::SPEED.bits() as i32,
            1,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(
        ret, -2,
        "helio_cross(Sun) should return UnsupportedFlags (-2)"
    );
    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn helio_cross_moon_rejected() {
    let handle = unsafe { default_handle() };
    let mut jx = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_helio_cross(
            handle,
            Body::Moon.to_raw_id(),
            0.0,
            2451545.0,
            CalcFlags::SPEED.bits() as i32,
            1,
            &mut jx,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert_eq!(
        ret, -2,
        "helio_cross(Moon) should return UnsupportedFlags (-2)"
    );
    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn solcross_null_jx() {
    let handle = unsafe { default_handle() };
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::pheno::swisseph_solcross(
            handle,
            0.0,
            2451545.0,
            0,
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };
    assert!(ret < 0);
    unsafe { swisseph_ffi::swisseph_free(handle) };
}
