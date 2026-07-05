use std::ffi::c_char;
use std::ptr;

use swisseph_ffi::SweEphemeris;
use swisseph_ffi::config::SweConfig;
use swisseph_ffi::error::SweErrorCode;

const PLANETS: [i32; 10] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
const SEFLG_SPEED: i32 = 256;

unsafe fn default_config() -> SweConfig {
    unsafe {
        let mut config = std::mem::zeroed::<SweConfig>();
        swisseph_ffi::config::swisseph_config_default(&mut config);
        config
    }
}

unsafe fn make_handle() -> *mut SweEphemeris {
    unsafe {
        let config = default_config();
        let mut handle: *mut SweEphemeris = ptr::null_mut();
        let mut err = [0u8; 256];
        let ret = swisseph_ffi::swisseph_new(
            &config,
            &mut handle,
            err.as_mut_ptr() as *mut c_char,
            err.len(),
        );
        assert_eq!(ret, SweErrorCode::Ok as i32);
        assert!(!handle.is_null());
        handle
    }
}

fn epochs() -> impl Iterator<Item = f64> {
    (0..100).map(|i| 2415020.5 + i as f64 * 500.0)
}

fn to_bits(xx: &[f64; 6], flags_used: i32) -> ([u64; 6], i32) {
    (xx.map(f64::to_bits), flags_used)
}

unsafe fn compute_all_ffi(handle: *const SweEphemeris) -> Vec<([u64; 6], i32)> {
    let mut out = Vec::with_capacity(100 * PLANETS.len());
    for jd in epochs() {
        for &ipl in &PLANETS {
            let mut xx = [0.0f64; 6];
            let mut flags_used: i32 = 0;
            let mut err = [0u8; 256];
            let ret = unsafe {
                swisseph_ffi::swisseph_calc_ut(
                    handle,
                    jd,
                    ipl,
                    SEFLG_SPEED,
                    ptr::null(),
                    ptr::null(),
                    xx.as_mut_ptr(),
                    &mut flags_used,
                    err.as_mut_ptr() as *mut c_char,
                    err.len(),
                )
            };
            assert_eq!(
                ret,
                SweErrorCode::Ok as i32,
                "calc_ut(jd={jd}, ipl={ipl}) failed: {}",
                String::from_utf8_lossy(&err)
            );
            out.push(to_bits(&xx, flags_used));
        }
    }
    out
}

#[test]
fn concurrent_ffi_bitwise_deterministic() {
    let handle = unsafe { make_handle() };
    let handle_int = handle as usize;
    let baseline = unsafe { compute_all_ffi(handle) };

    std::thread::scope(|s| {
        let handles: Vec<_> = (0..8)
            .map(|_| {
                s.spawn(move || {
                    let ptr = handle_int as *const SweEphemeris;
                    unsafe { compute_all_ffi(ptr) }
                })
            })
            .collect();
        for (i, h) in handles.into_iter().enumerate() {
            let results = h.join().expect("thread panicked");
            assert_eq!(
                results, baseline,
                "thread {i} diverged bitwise from the serial baseline"
            );
        }
    });

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn concurrent_ffi_via_usize_handle() {
    let handle = unsafe { make_handle() };

    let handle_int = handle as usize;
    let baseline = unsafe { compute_all_ffi(handle) };

    std::thread::scope(|s| {
        let handles: Vec<_> = (0..8)
            .map(|_| {
                s.spawn(move || {
                    let ptr = handle_int as *const SweEphemeris;
                    unsafe { compute_all_ffi(ptr) }
                })
            })
            .collect();
        for (i, h) in handles.into_iter().enumerate() {
            let results = h.join().expect("thread panicked");
            assert_eq!(
                results, baseline,
                "thread {i} (via usize round-trip) diverged from baseline"
            );
        }
    });

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn err_buf_isolation_across_threads() {
    let handle = unsafe { make_handle() };
    let handle_int = handle as usize;

    std::thread::scope(|s| {
        let handles: Vec<_> = (0..8)
            .map(|tid| {
                s.spawn(move || {
                    let ptr = handle_int as *const SweEphemeris;
                    let mut err = [0u8; 256];
                    let mut xx = [0.0f64; 6];

                    let ret = unsafe {
                        swisseph_ffi::swisseph_calc_ut(
                            ptr,
                            2451545.0,
                            -999,
                            SEFLG_SPEED,
                            ptr::null(),
                            ptr::null(),
                            xx.as_mut_ptr(),
                            ptr::null_mut(),
                            err.as_mut_ptr() as *mut c_char,
                            err.len(),
                        )
                    };
                    assert!(ret < 0, "thread {tid}: expected error for invalid body");

                    let msg = unsafe { std::ffi::CStr::from_ptr(err.as_ptr() as *const c_char) };
                    let msg_str = msg.to_str().unwrap();
                    assert!(
                        msg_str.contains("invalid body"),
                        "thread {tid}: err_buf = '{msg_str}', expected 'invalid body'"
                    );
                })
            })
            .collect();

        for h in handles {
            h.join().expect("thread panicked");
        }
    });

    unsafe { swisseph_ffi::swisseph_free(handle) };
}

#[test]
fn handle_as_int_roundtrip() {
    let handle = unsafe { make_handle() };

    let int_val = handle as usize;
    let recovered = int_val as *mut SweEphemeris;
    assert_eq!(handle, recovered);

    let mut xx = [0.0f64; 6];
    let mut err = [0u8; 256];
    let ret = unsafe {
        swisseph_ffi::swisseph_calc_ut(
            recovered,
            2451545.0,
            0,
            SEFLG_SPEED,
            ptr::null(),
            ptr::null(),
            xx.as_mut_ptr(),
            ptr::null_mut(),
            err.as_mut_ptr() as *mut c_char,
            err.len(),
        )
    };
    assert_eq!(ret, SweErrorCode::Ok as i32);
    assert!(xx[0] > 270.0 && xx[0] < 290.0);

    unsafe { swisseph_ffi::swisseph_free(recovered) };
}
