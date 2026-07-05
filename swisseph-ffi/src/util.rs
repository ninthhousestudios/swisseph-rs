use swisseph::math;

// ---------------------------------------------------------------------------
// Handle-free math/coordinate utilities
// ---------------------------------------------------------------------------

/// Normalize degrees to [0, 360). Port of `swe_degnorm`.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_degnorm(x: f64) -> f64 {
    math::normalize_degrees(x)
}

/// Normalize radians to [0, 2π). Port of `swe_radnorm`.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_radnorm(x: f64) -> f64 {
    math::normalize_radians(x)
}

/// Difference `p1 - p2` normalized to [0, 360). Port of `swe_difdegn`.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_difdegn(p1: f64, p2: f64) -> f64 {
    math::diff_degrees_norm(p1, p2)
}

/// Difference `p1 - p2` normalized to [-180, 180). Port of `swe_difdeg2n`.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_difdeg2n(p1: f64, p2: f64) -> f64 {
    math::diff_degrees(p1, p2)
}

/// Midpoint of two degree values on the 360° circle. Port of `swe_deg_midp`.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_deg_midp(x1: f64, x0: f64) -> f64 {
    math::midpoint_degrees(x1, x0)
}

/// Midpoint of two radian values on the 2π circle. Port of `swe_rad_midp`.
#[unsafe(no_mangle)]
pub extern "C" fn swisseph_rad_midp(x1: f64, x0: f64) -> f64 {
    math::midpoint_radians(x1, x0)
}

/// Coordinate transformation (ecliptic ↔ equatorial).
/// `xpo[3]` in, `xpn[3]` out. `eps` in degrees.
///
/// # Safety
/// - `xpo` must point to 3 readable `f64` values.
/// - `xpn` must point to 3 writable `f64` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_cotrans(xpo: *const f64, xpn: *mut f64, eps: f64) {
    if xpo.is_null() || xpn.is_null() {
        return;
    }
    let input = unsafe { [*xpo, *xpo.add(1), *xpo.add(2)] };
    let result = math::cotrans(input, eps);
    unsafe {
        *xpn = result[0];
        *xpn.add(1) = result[1];
        *xpn.add(2) = result[2];
    }
}

/// Coordinate transformation with speed (ecliptic ↔ equatorial).
/// `xpo[6]` in (pos + speed), `xpn[6]` out. `eps` in degrees.
///
/// # Safety
/// - `xpo` must point to 6 readable `f64` values.
/// - `xpn` must point to 6 writable `f64` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_cotrans_sp(xpo: *const f64, xpn: *mut f64, eps: f64) {
    if xpo.is_null() || xpn.is_null() {
        return;
    }
    let input = unsafe {
        [
            *xpo,
            *xpo.add(1),
            *xpo.add(2),
            *xpo.add(3),
            *xpo.add(4),
            *xpo.add(5),
        ]
    };
    let result = math::cotrans_with_speed(input, eps);
    unsafe {
        for i in 0..6 {
            *xpn.add(i) = result[i];
        }
    }
}
