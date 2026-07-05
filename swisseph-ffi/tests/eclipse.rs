use std::ffi::c_char;
use std::ptr;

use swisseph::Ephemeris;
use swisseph::config::EphemerisConfig;
use swisseph::flags::{CalcFlags, EclipseFlags, RiseSetFlags};
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
// rise_trans
// ---------------------------------------------------------------------------

#[test]
fn rise_trans_sun_rise() {
    let eph = make_eph();
    let tjd_ut = 2451545.0;
    let geopos = [8.55, 47.37, 500.0];
    let body = Body::Sun;
    let rsmi = RiseSetFlags::RISE;

    let lib_result = eph
        .rise_trans(
            tjd_ut,
            body,
            None,
            CalcFlags::empty(),
            rsmi,
            geopos,
            0.0,
            0.0,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_rise_trans(
            handle,
            tjd_ut,
            body.to_raw_id(),
            ptr::null(),
            0,
            rsmi.bits() as i32,
            geopos.as_ptr(),
            0.0,
            0.0,
            &mut tret,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(tret, lib_result.time, 1e-15, "rise_trans Sun rise");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn rise_trans_moon_set() {
    let eph = make_eph();
    let tjd_ut = 2451545.0;
    let geopos = [8.55, 47.37, 500.0];
    let body = Body::Moon;
    let rsmi = RiseSetFlags::SET;

    let lib_result = eph
        .rise_trans(
            tjd_ut,
            body,
            None,
            CalcFlags::empty(),
            rsmi,
            geopos,
            0.0,
            0.0,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_rise_trans(
            handle,
            tjd_ut,
            body.to_raw_id(),
            ptr::null(),
            0,
            rsmi.bits() as i32,
            geopos.as_ptr(),
            0.0,
            0.0,
            &mut tret,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(tret, lib_result.time, 1e-15, "rise_trans Moon set");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn rise_trans_circumpolar() {
    let handle = unsafe { default_handle() };
    let tjd_ut = 2451545.0;
    let geopos = [18.0, 69.65, 0.0]; // Tromso, polar
    let mut tret = 0.0f64;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_rise_trans(
            handle,
            tjd_ut,
            Body::Sun.to_raw_id(),
            ptr::null(),
            0,
            (RiseSetFlags::RISE | RiseSetFlags::FORCE_SLOW).bits() as i32,
            geopos.as_ptr(),
            0.0,
            0.0,
            &mut tret,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, -2, "circumpolar should return -2");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// sol_eclipse_when_glob
// ---------------------------------------------------------------------------

#[test]
fn sol_eclipse_when_glob_forward() {
    let eph = make_eph();
    let tjd_start = 2457000.0;

    let lib_result = eph
        .sol_eclipse_when_glob(tjd_start, CalcFlags::empty(), EclipseFlags::empty(), false)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_sol_eclipse_when_glob(
            handle,
            tjd_start,
            0,
            0,
            0,
            tret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0, "should return positive eclipse flags, got {ret}");
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(tret[0], lib_result.time_maximum, 1e-15, "time_maximum");
    assert_eps(tret[2], lib_result.time_begin, 1e-15, "time_begin");
    assert_eps(tret[3], lib_result.time_end, 1e-15, "time_end");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn sol_eclipse_when_glob_backward() {
    let eph = make_eph();
    let tjd_start = 2460000.0;

    let lib_result = eph
        .sol_eclipse_when_glob(tjd_start, CalcFlags::empty(), EclipseFlags::empty(), true)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_sol_eclipse_when_glob(
            handle,
            tjd_start,
            0,
            0,
            1,
            tret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0);
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(
        tret[0],
        lib_result.time_maximum,
        1e-15,
        "backward time_maximum",
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// sol_eclipse_how
// ---------------------------------------------------------------------------

#[test]
fn sol_eclipse_how_known_eclipse() {
    let eph = make_eph();
    let tjd_ut = 2458353.03681; // 2018-08-11 near a solar eclipse
    let geopos = [8.55, 47.37, 500.0];

    let lib_result = eph
        .sol_eclipse_how(tjd_ut, CalcFlags::empty(), geopos)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_sol_eclipse_how(
            handle,
            tjd_ut,
            0,
            geopos.as_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(attr[0], lib_result.magnitude, 1e-15, "magnitude");
    assert_eps(attr[1], lib_result.diameter_ratio, 1e-15, "diameter_ratio");
    assert_eps(attr[4], lib_result.azimuth, 1e-15, "azimuth");
    assert_eps(attr[9], lib_result.saros_series, 1e-15, "saros_series");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// sol_eclipse_where
// ---------------------------------------------------------------------------

#[test]
fn sol_eclipse_where_known() {
    let eph = make_eph();
    // 1999-Aug-11 total solar eclipse at maximum
    let tjd_ut = 2451401.9604166667;

    let lib_result = eph.sol_eclipse_where(tjd_ut, CalcFlags::empty()).unwrap();
    assert!(
        !lib_result.flags.is_empty(),
        "expected a central eclipse at this epoch"
    );

    let handle = unsafe { default_handle() };
    let mut geopos = [0.0f64; 10];
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_sol_eclipse_where(
            handle,
            tjd_ut,
            0,
            geopos.as_mut_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(geopos[0], lib_result.central_longitude, 1e-15, "lon");
    assert_eps(geopos[1], lib_result.central_latitude, 1e-15, "lat");
    assert_eps(geopos[2], lib_result.core_diameter_km, 1e-15, "core_diam");

    // attr[20] populated via eclipse_how at the central point
    let how = eph
        .eclipse_how_at(
            tjd_ut,
            Body::Sun,
            None,
            CalcFlags::empty(),
            [
                lib_result.central_longitude,
                lib_result.central_latitude,
                0.0,
            ],
        )
        .unwrap();
    assert_eps(attr[0], how.magnitude, 1e-15, "attr[0] magnitude");
    assert_eps(attr[1], how.diameter_ratio, 1e-15, "attr[1] diameter_ratio");
    assert_eps(attr[2], how.obscuration, 1e-15, "attr[2] obscuration");
    assert_eps(
        attr[3],
        lib_result.core_diameter_km,
        1e-15,
        "attr[3] core_diameter_km (dcore[0])",
    );
    assert_eps(attr[4], how.azimuth, 1e-15, "attr[4] azimuth");
    assert_eps(attr[5], how.true_altitude, 1e-15, "attr[5] true_altitude");
    assert_eps(
        attr[6],
        how.apparent_altitude,
        1e-15,
        "attr[6] apparent_altitude",
    );
    assert_eps(attr[7], how.elongation, 1e-15, "attr[7] elongation");
    assert_eps(attr[8], how.nasa_magnitude, 1e-15, "attr[8] nasa_magnitude");
    assert_eps(attr[9], how.saros_series, 1e-15, "attr[9] saros_series");
    assert_eps(attr[10], how.saros_member, 1e-15, "attr[10] saros_member");
    for i in 11..20 {
        assert_eq!(attr[i], 0.0, "attr[{i}] should be 0");
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// sol_eclipse_when_loc
// ---------------------------------------------------------------------------

#[test]
fn sol_eclipse_when_loc_forward() {
    let eph = make_eph();
    let tjd_start = 2457000.0;
    let geopos = [8.55, 47.37, 500.0];

    let lib_result = eph
        .sol_eclipse_when_loc(tjd_start, CalcFlags::empty(), geopos, false)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_sol_eclipse_when_loc(
            handle,
            tjd_start,
            0,
            geopos.as_ptr(),
            0,
            tret.as_mut_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0);
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(tret[0], lib_result.time_maximum, 1e-15, "time_maximum");
    assert_eps(tret[1], lib_result.time_first_contact, 1e-15, "1st_contact");
    assert_eps(attr[0], lib_result.attr.magnitude, 1e-15, "attr magnitude");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// lun_eclipse_how
// ---------------------------------------------------------------------------

#[test]
fn lun_eclipse_how_known() {
    let eph = make_eph();
    let tjd_ut = 2452279.9283; // 2001-01-09 total lunar eclipse
    let geopos = [8.55, 47.37, 500.0];

    let lib_result = eph
        .lun_eclipse_how(tjd_ut, CalcFlags::empty(), geopos)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_lun_eclipse_how(
            handle,
            tjd_ut,
            0,
            geopos.as_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(attr[0], lib_result.umbral_magnitude, 1e-15, "umbral_mag");
    assert_eps(attr[1], lib_result.penumbral_magnitude, 1e-15, "pen_mag");
    assert_eps(attr[4], lib_result.azimuth, 1e-15, "azimuth");
    assert_eps(attr[8], lib_result.umbral_magnitude, 1e-15, "attr[8] dup");
    assert_eps(attr[9], lib_result.saros_series, 1e-15, "saros_series");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// lun_eclipse_when
// ---------------------------------------------------------------------------

#[test]
fn lun_eclipse_when_forward() {
    let eph = make_eph();
    let tjd_start = 2457000.0;

    let lib_result = eph
        .lun_eclipse_when(tjd_start, CalcFlags::empty(), EclipseFlags::empty(), false)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_lun_eclipse_when(
            handle,
            tjd_start,
            0,
            0,
            0,
            tret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0);
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(tret[0], lib_result.time_maximum, 1e-15, "time_maximum");
    assert_eq!(tret[1], 0.0, "tret[1] should be 0 for lunar");
    assert_eps(tret[6], lib_result.time_penumbral_begin, 1e-15, "pen_begin");
    assert_eps(tret[7], lib_result.time_penumbral_end, 1e-15, "pen_end");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// lun_eclipse_when_loc
// ---------------------------------------------------------------------------

#[test]
fn lun_eclipse_when_loc_forward() {
    let eph = make_eph();
    let tjd_start = 2457000.0;
    let geopos = [8.55, 47.37, 500.0];

    let lib_result = eph
        .lun_eclipse_when_loc(tjd_start, CalcFlags::empty(), geopos, false)
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_lun_eclipse_when_loc(
            handle,
            tjd_start,
            0,
            geopos.as_ptr(),
            0,
            tret.as_mut_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0);
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(tret[0], lib_result.time_maximum, 1e-15, "time_maximum");
    assert_eps(tret[8], lib_result.time_moonrise, 1e-15, "moonrise");
    assert_eps(tret[9], lib_result.time_moonset, 1e-15, "moonset");
    assert_eps(
        attr[0],
        lib_result.attr.umbral_magnitude,
        1e-15,
        "umbral_mag",
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// lun_occult_when_glob (requires ephe for star data)
// ---------------------------------------------------------------------------

#[test]
fn lun_occult_when_glob_venus() {
    let eph = make_eph();
    let tjd_start = 2451545.0;

    let lib_result = eph
        .lun_occult_when_glob(
            tjd_start,
            Body::Venus,
            None,
            CalcFlags::empty(),
            EclipseFlags::empty(),
            false,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_lun_occult_when_glob(
            handle,
            tjd_start,
            Body::Venus.to_raw_id(),
            ptr::null(),
            0,
            0,
            0,
            tret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0);
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(tret[0], lib_result.time_maximum, 1e-15, "time_maximum");
    assert_eps(tret[2], lib_result.time_begin, 1e-15, "time_begin");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// lun_occult_where
// ---------------------------------------------------------------------------

#[test]
fn lun_occult_where_venus() {
    let eph = make_eph();
    // Venus occultation maximum from occ_when_glob (CENTRAL|TOTAL, retval=5)
    let tjd_ut = 2451607.5448415945;

    let lib_result = eph
        .lun_occult_where(tjd_ut, Body::Venus, None, CalcFlags::empty())
        .unwrap();
    assert!(
        !lib_result.flags.is_empty(),
        "expected an occultation at this epoch"
    );

    let handle = unsafe { default_handle() };
    let mut geopos = [0.0f64; 10];
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_lun_occult_where(
            handle,
            tjd_ut,
            Body::Venus.to_raw_id(),
            ptr::null(),
            0,
            geopos.as_mut_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(geopos[0], lib_result.central_longitude, 1e-15, "lon");
    assert_eps(geopos[1], lib_result.central_latitude, 1e-15, "lat");

    // attr[20] populated via eclipse_how at the central point
    let how = eph
        .eclipse_how_at(
            tjd_ut,
            Body::Venus,
            None,
            CalcFlags::empty(),
            [
                lib_result.central_longitude,
                lib_result.central_latitude,
                0.0,
            ],
        )
        .unwrap();
    assert_eps(attr[0], how.magnitude, 1e-15, "attr[0] magnitude");
    assert_eps(attr[1], how.diameter_ratio, 1e-15, "attr[1] diameter_ratio");
    assert_eps(attr[2], how.obscuration, 1e-15, "attr[2] obscuration");
    assert_eps(
        attr[3],
        lib_result.core_diameter_km,
        1e-15,
        "attr[3] core_diameter_km (dcore[0])",
    );
    assert_eps(attr[4], how.azimuth, 1e-15, "attr[4] azimuth");
    assert_eps(attr[5], how.true_altitude, 1e-15, "attr[5] true_altitude");
    assert_eps(attr[7], how.elongation, 1e-15, "attr[7] elongation");
    for i in 11..20 {
        assert_eq!(attr[i], 0.0, "attr[{i}] should be 0");
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// lun_occult_when_loc
// ---------------------------------------------------------------------------

#[test]
fn lun_occult_when_loc_venus() {
    let eph = make_eph();
    let tjd_start = 2451545.0;
    let geopos = [8.55, 47.37, 500.0];

    let lib_result = eph
        .lun_occult_when_loc(
            tjd_start,
            Body::Venus,
            None,
            CalcFlags::empty(),
            geopos,
            false,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = [0.0f64; 10];
    let mut attr = [0.0f64; 20];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_lun_occult_when_loc(
            handle,
            tjd_start,
            Body::Venus.to_raw_id(),
            ptr::null(),
            0,
            geopos.as_ptr(),
            0,
            tret.as_mut_ptr(),
            attr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret > 0);
    assert_eq!(ret as u32, lib_result.flags.bits());
    assert_eps(tret[0], lib_result.time_maximum, 1e-15, "time_maximum");
    assert_eps(attr[0], lib_result.attr.magnitude, 1e-15, "magnitude");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// rise_trans_true_hor
// ---------------------------------------------------------------------------

#[test]
fn rise_trans_true_hor_with_height() {
    let eph = make_eph();
    let tjd_ut = 2451545.0;
    let geopos = [8.55, 47.37, 500.0];
    let body = Body::Sun;
    let rsmi = RiseSetFlags::RISE | RiseSetFlags::FORCE_SLOW;

    let lib_result = eph
        .rise_trans_true_hor(
            tjd_ut,
            body,
            None,
            CalcFlags::empty(),
            rsmi,
            geopos,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut tret = 0.0f64;
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::eclipse::swisseph_rise_trans_true_hor(
            handle,
            tjd_ut,
            body.to_raw_id(),
            ptr::null(),
            0,
            rsmi.bits() as i32,
            geopos.as_ptr(),
            0.0,
            0.0,
            0.0,
            &mut tret,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0);
    assert_eps(tret, lib_result.time, 1e-15, "rise_trans_true_hor Sun rise");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}
