use std::ffi::{CString, c_char};
use std::ptr;

use swisseph::Ephemeris;
use swisseph::config::EphemerisConfig;
use swisseph::flags::{CalcFlags, HeliacalFlags};
use swisseph::heliacal::HeliacalEventType;

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
// vis_limit_mag
// ---------------------------------------------------------------------------

#[test]
fn vis_limit_mag_venus() {
    let eph = make_eph();
    let tjd_ut = 2452275.5;
    let dgeo = [31.25, 30.1, 30.0];
    let mut datm = [1013.25, 15.0, 40.0, 40.0];
    let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let epheflag = CalcFlags::empty();
    let helflag = HeliacalFlags::empty();

    let lib_result = eph
        .vis_limit_mag(
            tjd_ut, &dgeo, &mut datm, &mut dobs, "venus", epheflag, helflag,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let name = CString::new("venus").unwrap();
    let mut dret = [0.0f64; 8];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_vis_limit_mag(
            handle,
            tjd_ut,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            name.as_ptr(),
            (epheflag.bits() | helflag.bits()) as i32,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret >= 0, "vis_limit_mag returned error: {ret}");
    assert_eq!(ret, lib_result.vision.bits() as i32, "vision flags");

    assert_eps(
        dret[0],
        lib_result.limiting_magnitude,
        1e-15,
        "limiting_magnitude",
    );
    assert_eps(
        dret[1],
        lib_result.altitude_object,
        1e-15,
        "altitude_object",
    );
    assert_eps(dret[2], lib_result.azimuth_object, 1e-15, "azimuth_object");
    assert_eps(dret[3], lib_result.altitude_sun, 1e-15, "altitude_sun");
    assert_eps(dret[4], lib_result.azimuth_sun, 1e-15, "azimuth_sun");
    assert_eps(dret[5], lib_result.altitude_moon, 1e-15, "altitude_moon");
    assert_eps(dret[6], lib_result.azimuth_moon, 1e-15, "azimuth_moon");
    assert_eps(
        dret[7],
        lib_result.magnitude_object,
        1e-15,
        "magnitude_object",
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// heliacal_angle
// ---------------------------------------------------------------------------

#[test]
fn heliacal_angle_basic() {
    let eph = make_eph();
    let tjd_ut = 2452275.5;
    let dgeo = [31.25, 30.1, 30.0];
    let mut datm = [1013.25, 15.0, 40.0, 40.0];
    let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let helflag = HeliacalFlags::empty();
    let mag = -3.9;
    let azi_obj = 100.0;
    let azi_sun = 90.0;
    let azi_moon = 200.0;
    let alt_moon = -10.0;

    let lib_result = eph
        .heliacal_angle(
            tjd_ut, &dgeo, &mut datm, &mut dobs, helflag, mag, azi_obj, azi_sun, azi_moon, alt_moon,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut dret = [0.0f64; 3];
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_heliacal_angle(
            handle,
            tjd_ut,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            helflag.bits() as i32,
            mag,
            azi_obj,
            azi_sun,
            azi_moon,
            alt_moon,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0, "heliacal_angle returned error: {ret}");
    assert_eps(
        dret[0],
        lib_result.optimal_altitude,
        1e-15,
        "optimal_altitude",
    );
    assert_eps(dret[1], lib_result.arcus_visionis, 1e-15, "arcus_visionis");
    assert_eps(
        dret[2],
        lib_result.sun_altitude_diff,
        1e-15,
        "sun_altitude_diff",
    );

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// topo_arcus_visionis
// ---------------------------------------------------------------------------

#[test]
fn topo_arcus_visionis_basic() {
    let eph = make_eph();
    let tjd_ut = 2452275.5;
    let dgeo = [31.25, 30.1, 30.0];
    let mut datm = [1013.25, 15.0, 40.0, 40.0];
    let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let helflag = HeliacalFlags::empty();
    let mag = -3.9;
    let azi_obj = 100.0;
    let alt_obj = 5.0;
    let azi_sun = 90.0;
    let azi_moon = 200.0;
    let alt_moon = -10.0;

    let lib_result = eph
        .topo_arcus_visionis(
            tjd_ut, &dgeo, &mut datm, &mut dobs, helflag, mag, azi_obj, alt_obj, azi_sun, azi_moon,
            alt_moon,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let mut dret_val = 0.0f64;
    let mut err_buf = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_topo_arcus_visionis(
            handle,
            tjd_ut,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            helflag.bits() as i32,
            mag,
            azi_obj,
            alt_obj,
            azi_sun,
            azi_moon,
            alt_moon,
            &mut dret_val,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0, "topo_arcus_visionis returned error: {ret}");
    assert_eps(dret_val, lib_result, 1e-15, "topo_arcus_visionis");

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// heliacal_pheno_ut
// ---------------------------------------------------------------------------

#[test]
fn heliacal_pheno_ut_venus() {
    let eph = make_eph();
    let tjd_ut = 2452275.5;
    let dgeo = [31.25, 30.1, 30.0];
    let mut datm = [1013.25, 15.0, 40.0, 40.0];
    let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let epheflag = CalcFlags::empty();
    let helflag = HeliacalFlags::empty();
    let event = HeliacalEventType::MorningFirst;

    let lib_result = eph
        .heliacal_pheno_ut(
            tjd_ut, &dgeo, &mut datm, &mut dobs, "venus", event, epheflag, helflag,
        )
        .unwrap();

    let lib_arr = lib_result.as_array();

    let handle = unsafe { default_handle() };
    let name = CString::new("venus").unwrap();
    let mut darr = [0.0f64; 50];
    let mut err_buf = [0u8; 256];
    let combined = (epheflag.bits() | helflag.bits()) as i32;
    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_heliacal_pheno_ut(
            handle,
            tjd_ut,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            name.as_ptr(),
            event as i32,
            combined,
            darr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0, "heliacal_pheno_ut returned error: {ret}");

    for i in 0..28 {
        assert_eps(darr[i], lib_arr[i], 1e-15, &format!("darr[{i}]"));
    }
    for i in 28..50 {
        assert_eq!(darr[i], 0.0, "darr[{i}] should be zeroed");
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// heliacal_ut
// ---------------------------------------------------------------------------

#[test]
fn heliacal_ut_venus_morning_first() {
    let eph = make_eph();
    let tjd_start = 2452275.5;
    let dgeo = [31.25, 30.1, 30.0];
    let mut datm = [1013.25, 15.0, 40.0, 40.0];
    let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let epheflag = CalcFlags::empty();
    let helflag = HeliacalFlags::empty();
    let event = HeliacalEventType::MorningFirst;

    let lib_result = eph
        .heliacal_ut(
            tjd_start, &dgeo, &mut datm, &mut dobs, "venus", event, epheflag, helflag,
        )
        .unwrap();

    let handle = unsafe { default_handle() };
    let name = CString::new("venus").unwrap();
    let mut dret = [0.0f64; 50];
    let mut err_buf = [0u8; 256];
    let combined = (epheflag.bits() | helflag.bits()) as i32;
    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_heliacal_ut(
            handle,
            tjd_start,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            name.as_ptr(),
            event as i32,
            combined,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert_eq!(ret, 0, "heliacal_ut returned error: {ret}");
    assert_eps(dret[0], lib_result.start_visible, 1e-15, "start_visible");
    assert_eps(
        dret[1],
        lib_result.optimum_visibility,
        1e-15,
        "optimum_visibility",
    );
    assert_eps(dret[2], lib_result.end_visible, 1e-15, "end_visible");
    for i in 3..50 {
        assert_eq!(dret[i], 0.0, "dret[{i}] should be zeroed");
    }

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

// ---------------------------------------------------------------------------
// Error paths
// ---------------------------------------------------------------------------

#[test]
fn heliacal_ut_bad_event_type() {
    let handle = unsafe { default_handle() };
    let name = CString::new("venus").unwrap();
    let dgeo = [31.25, 30.1, 30.0];
    let datm = [1013.25, 15.0, 40.0, 40.0];
    let dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let mut dret = [0.0f64; 50];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_heliacal_ut(
            handle,
            2452275.5,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            name.as_ptr(),
            99, // invalid event type
            0,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret < 0, "expected error for invalid event type, got {ret}");
    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn heliacal_ut_null_object_name() {
    let handle = unsafe { default_handle() };
    let dgeo = [31.25, 30.1, 30.0];
    let datm = [1013.25, 15.0, 40.0, 40.0];
    let dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let mut dret = [0.0f64; 50];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_heliacal_ut(
            handle,
            2452275.5,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            ptr::null(), // null object name
            1,
            0,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret < 0, "expected error for null object name, got {ret}");
    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn vis_limit_mag_null_dgeo() {
    let handle = unsafe { default_handle() };
    let name = CString::new("venus").unwrap();
    let datm = [1013.25, 15.0, 40.0, 40.0];
    let dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let mut dret = [0.0f64; 8];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_vis_limit_mag(
            handle,
            2452275.5,
            ptr::null(), // null dgeo
            datm.as_ptr(),
            dobs.as_ptr(),
            name.as_ptr(),
            0,
            dret.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret < 0, "expected error for null dgeo, got {ret}");
    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn heliacal_pheno_ut_bad_event_type() {
    let handle = unsafe { default_handle() };
    let name = CString::new("venus").unwrap();
    let dgeo = [31.25, 30.1, 30.0];
    let datm = [1013.25, 15.0, 40.0, 40.0];
    let dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let mut darr = [0.0f64; 50];
    let mut err_buf = [0u8; 256];

    let ret = unsafe {
        swisseph_ffi::heliacal::swisseph_heliacal_pheno_ut(
            handle,
            2452275.5,
            dgeo.as_ptr(),
            datm.as_ptr(),
            dobs.as_ptr(),
            name.as_ptr(),
            0, // invalid event type (must be 1-6)
            0,
            darr.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        )
    };

    assert!(ret < 0, "expected error for invalid event type, got {ret}");
    unsafe { swisseph_ffi::swisseph_free(handle) };
}
