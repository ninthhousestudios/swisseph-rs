use std::ffi::c_char;
use std::ptr;

use swisseph_ffi::SweEphemeris;
use swisseph_ffi::config::SweConfig;
use swisseph_ffi::error::SweErrorCode;

unsafe fn default_config() -> SweConfig {
    unsafe {
        let mut config = std::mem::zeroed::<SweConfig>();
        swisseph_ffi::config::swisseph_config_default(&mut config);
        config
    }
}

#[test]
fn new_free_roundtrip() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, SweErrorCode::Ok as i32);
        assert!(!handle.is_null());
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn free_null_is_noop() {
    unsafe {
        swisseph_ffi::swisseph_free(ptr::null_mut());
    }
}

#[test]
fn null_out_returns_invalid_arg() {
    unsafe {
        let config = default_config();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, SweErrorCode::InvalidArg as i32);
    }
}

#[test]
fn tiny_err_buf_truncates_with_nul() {
    unsafe {
        let config = default_config();
        let mut err_buf = [0xFFu8; 4];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, SweErrorCode::InvalidArg as i32);
        // Must be NUL-terminated within the 4-byte buffer
        let nul_pos = err_buf.iter().position(|&b| b == 0).unwrap();
        assert!(nul_pos < 4);
    }
}

#[test]
fn calc_ut_moshier_sun() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);

        // Compute Sun at J2000.0 with SPEED flag (256)
        let mut xx = [0.0f64; 6];
        let mut flags_used: i32 = 0;
        let ret = swisseph_ffi::swisseph_calc_ut(
            handle,
            2451545.0,   // J2000.0
            0,           // Sun
            256,         // SEFLG_SPEED
            ptr::null(), // no geopos override
            ptr::null(), // no sidereal override
            xx.as_mut_ptr(),
            &mut flags_used,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);

        // Verify against Ephemeris::calc_ut directly
        let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
        let result = eph
            .calc_ut(2451545.0, swisseph::Body::Sun, swisseph::CalcFlags::SPEED)
            .unwrap();
        for i in 0..6 {
            assert!(
                (xx[i] - result.data[i]).abs() < 1e-15,
                "xx[{i}] mismatch: ffi={} rust={}",
                xx[i],
                result.data[i]
            );
        }
        assert_eq!(flags_used, result.flags_used.bits() as i32);

        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn calc_ut_moshier_moon() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);

        let mut xx = [0.0f64; 6];
        let mut flags_used: i32 = 0;
        let ret = swisseph_ffi::swisseph_calc_ut(
            handle,
            2451545.0,
            1, // Moon
            256,
            ptr::null(),
            ptr::null(),
            xx.as_mut_ptr(),
            &mut flags_used,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);

        let eph = swisseph::Ephemeris::new(swisseph::EphemerisConfig::default()).unwrap();
        let result = eph
            .calc_ut(2451545.0, swisseph::Body::Moon, swisseph::CalcFlags::SPEED)
            .unwrap();
        for i in 0..6 {
            assert!(
                (xx[i] - result.data[i]).abs() < 1e-15,
                "xx[{i}] mismatch: ffi={} rust={}",
                xx[i],
                result.data[i]
            );
        }

        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn calc_ut_invalid_body() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );

        let mut xx = [0.0f64; 6];
        let ret = swisseph_ffi::swisseph_calc_ut(
            handle,
            2451545.0,
            -999, // truly invalid body (no Body variant for negative values below valid range)
            256,
            ptr::null(),
            ptr::null(),
            xx.as_mut_ptr(),
            ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert!(ret < 0, "expected negative error code, got {ret}");

        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn calc_ut_with_geopos_override() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );

        // With topographic override — TOPOCTR flag is 32768
        let geopos = [8.55f64, 47.37, 500.0];
        let mut xx = [0.0f64; 6];
        let mut flags_used: i32 = 0;
        let ret = swisseph_ffi::swisseph_calc_ut(
            handle,
            2451545.0,
            1, // Moon
            256 | 32768,
            geopos.as_ptr(),
            ptr::null(),
            xx.as_mut_ptr(),
            &mut flags_used,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);

        // Verify against Rust API with config override
        let mut rust_config = swisseph::EphemerisConfig::default();
        rust_config.topographic = Some(swisseph::TopoPosition {
            longitude: 8.55,
            latitude: 47.37,
            altitude: 500.0,
        });
        let eph = swisseph::Ephemeris::new(rust_config).unwrap();
        let result = eph
            .calc_ut(
                2451545.0,
                swisseph::Body::Moon,
                swisseph::CalcFlags::SPEED | swisseph::CalcFlags::TOPOCTR,
            )
            .unwrap();
        for i in 0..6 {
            assert!(
                (xx[i] - result.data[i]).abs() < 1e-15,
                "xx[{i}] mismatch: ffi={} rust={}",
                xx[i],
                result.data[i]
            );
        }

        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn version_not_null() {
    let ptr = swisseph_ffi::swisseph_version();
    assert!(!ptr.is_null());
    let cstr = unsafe { std::ffi::CStr::from_ptr(ptr) };
    assert_eq!(cstr.to_str().unwrap(), "0.1.0");
}

#[test]
fn get_tid_acc_moshier_unresolved() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        let tid_acc = swisseph_ffi::swisseph_get_tid_acc(handle);
        // Moshier doesn't open files, so tid_acc stays unresolved (NaN from the FFI)
        assert!(tid_acc.is_nan());
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn get_tid_acc_explicit() {
    unsafe {
        let mut config = default_config();
        config.tidal_acceleration = -25.8; // explicit override, not NaN
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        let tid_acc = swisseph_ffi::swisseph_get_tid_acc(handle);
        assert!((tid_acc - (-25.8)).abs() < 1e-10);
        swisseph_ffi::swisseph_free(handle);
    }
}

#[test]
fn get_astro_models() {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err_buf = [0u8; 256];
        swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        let mut models = [0i32; 8];
        let ret = swisseph_ffi::swisseph_get_astro_models(
            handle,
            models.as_mut_ptr(),
            err_buf.as_mut_ptr() as *mut c_char,
            err_buf.len(),
        );
        assert_eq!(ret, 0);
        // All should be non-zero (resolved defaults)
        for (i, &v) in models.iter().enumerate() {
            assert!(v > 0, "models[{i}] should be > 0, got {v}");
        }
        swisseph_ffi::swisseph_free(handle);
    }
}
