use std::f64::consts::PI;

use crate::constants::{DEG180, DEG360, DEGTORAD, RADTODEG, TWOPI};
use crate::flags::SplitDegFlags;
use crate::types::DegreeParts;

// ---------------------------------------------------------------------------
// Angle normalization
// ---------------------------------------------------------------------------

pub fn normalize_degrees(x: f64) -> f64 {
    let mut y = x % 360.0;
    if y.abs() < 1e-13 {
        y = 0.0;
    }
    if y < 0.0 {
        y += 360.0;
    }
    y
}

pub fn normalize_radians(x: f64) -> f64 {
    let mut y = x % TWOPI;
    if y.abs() < 1e-13 {
        y = 0.0;
    }
    if y < 0.0 {
        y += TWOPI;
    }
    y
}

pub fn mod_2pi(x: f64) -> f64 {
    let mut y = x % TWOPI;
    if y < 0.0 {
        y += TWOPI;
    }
    y
}

pub fn mods3600(x: f64) -> f64 {
    x - 1_296_000.0 * (x / 1_296_000.0).floor()
}

// ---------------------------------------------------------------------------
// Angle differences
// ---------------------------------------------------------------------------

pub fn diff_degrees_norm(p1: f64, p2: f64) -> f64 {
    normalize_degrees(p1 - p2)
}

pub fn diff_degrees(p1: f64, p2: f64) -> f64 {
    let dif = normalize_degrees(p1 - p2);
    if dif >= 180.0 { dif - 360.0 } else { dif }
}

pub fn diff_radians(p1: f64, p2: f64) -> f64 {
    let dif = normalize_radians(p1 - p2);
    if dif >= PI { dif - TWOPI } else { dif }
}

// ---------------------------------------------------------------------------
// Midpoints
// ---------------------------------------------------------------------------

pub fn midpoint_degrees(x1: f64, x0: f64) -> f64 {
    let d = diff_degrees(x1, x0);
    normalize_degrees(x0 + d / 2.0)
}

pub fn midpoint_radians(x1: f64, x0: f64) -> f64 {
    DEGTORAD * midpoint_degrees(x1 * RADTODEG, x0 * RADTODEG)
}

// ---------------------------------------------------------------------------
// Centisecond angle functions
// ---------------------------------------------------------------------------

pub fn csnorm(mut p: i32) -> i32 {
    if p < 0 {
        loop {
            p += DEG360;
            if p >= 0 {
                break;
            }
        }
    } else if p >= DEG360 {
        loop {
            p -= DEG360;
            if p < DEG360 {
                break;
            }
        }
    }
    p
}

pub fn difcsn(p1: i32, p2: i32) -> i32 {
    csnorm(p1 - p2)
}

pub fn difcs2n(p1: i32, p2: i32) -> i32 {
    let dif = csnorm(p1 - p2);
    if dif >= DEG180 { dif - DEG360 } else { dif }
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

pub fn d2l(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5) as i32
    } else {
        -((0.5 - x) as i32)
    }
}

// ---------------------------------------------------------------------------
// Chebyshev evaluation (Broucke/Clenshaw, ACM algorithm 446)
// ---------------------------------------------------------------------------

pub fn chebyshev_eval(x: f64, coeffs: &[f64]) -> f64 {
    let x2 = x * 2.0;
    let mut br = 0.0;
    let mut brp2 = 0.0;
    let mut brpp = 0.0;
    for &c in coeffs.iter().rev() {
        brp2 = brpp;
        brpp = br;
        br = x2 * brpp - brp2 + c;
    }
    (br - brp2) * 0.5
}

pub fn chebyshev_deriv(x: f64, coeffs: &[f64]) -> f64 {
    let ncf = coeffs.len();
    if ncf <= 1 {
        return 0.0;
    }
    let x2 = x * 2.0;
    let mut bf = 0.0;
    let mut bj = 0.0;
    let mut xjp2 = 0.0;
    let mut xjpl = 0.0;
    let mut bjp2 = 0.0;
    let mut bjpl = 0.0;
    for j in (1..ncf).rev() {
        let dj = (j + j) as f64;
        let xj = coeffs[j] * dj + xjp2;
        bj = x2 * bjpl - bjp2 + xj;
        bf = bjp2;
        bjp2 = bjpl;
        bjpl = bj;
        xjp2 = xjpl;
        xjpl = xj;
    }
    (bj - bf) * 0.5
}

// ---------------------------------------------------------------------------
// Coordinate transforms — basic
// ---------------------------------------------------------------------------

pub fn rotate_x(pos: [f64; 3], eps: f64) -> [f64; 3] {
    let sineps = eps.sin();
    let coseps = eps.cos();
    rotate_x_sincos(pos, sineps, coseps)
}

pub fn rotate_x_sincos(pos: [f64; 3], sineps: f64, coseps: f64) -> [f64; 3] {
    [
        pos[0],
        pos[1] * coseps + pos[2] * sineps,
        -pos[1] * sineps + pos[2] * coseps,
    ]
}

pub fn cartesian_to_polar(x: [f64; 3]) -> [f64; 3] {
    if x[0] == 0.0 && x[1] == 0.0 && x[2] == 0.0 {
        return [0.0; 3];
    }
    let rxy_sq = x[0] * x[0] + x[1] * x[1];
    let dist = (rxy_sq + x[2] * x[2]).sqrt();
    let rxy = rxy_sq.sqrt();
    let mut lon = x[1].atan2(x[0]);
    if lon < 0.0 {
        lon += TWOPI;
    }
    let lat = if rxy == 0.0 {
        if x[2] >= 0.0 { PI / 2.0 } else { -(PI / 2.0) }
    } else {
        (x[2] / rxy).atan()
    };
    [lon, lat, dist]
}

pub fn polar_to_cartesian(l: [f64; 3]) -> [f64; 3] {
    let cosl1 = l[1].cos();
    [
        l[2] * cosl1 * l[0].cos(),
        l[2] * cosl1 * l[0].sin(),
        l[2] * l[1].sin(),
    ]
}

// ---------------------------------------------------------------------------
// Coordinate transforms — with speed (Jacobian velocity transform)
// ---------------------------------------------------------------------------

pub fn cartesian_to_polar_with_speed(x: [f64; 6]) -> [f64; 6] {
    if x[0] == 0.0 && x[1] == 0.0 && x[2] == 0.0 {
        let speed = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5]).sqrt();
        let vel_dir = cartesian_to_polar([x[3], x[4], x[5]]);
        return [vel_dir[0], vel_dir[1], 0.0, 0.0, 0.0, speed];
    }
    if x[3] == 0.0 && x[4] == 0.0 && x[5] == 0.0 {
        let pos = cartesian_to_polar([x[0], x[1], x[2]]);
        return [pos[0], pos[1], pos[2], 0.0, 0.0, 0.0];
    }
    let rxy_sq = x[0] * x[0] + x[1] * x[1];
    let rxyz = (rxy_sq + x[2] * x[2]).sqrt();
    let rxy = rxy_sq.sqrt();
    let mut lon = x[1].atan2(x[0]);
    if lon < 0.0 {
        lon += TWOPI;
    }
    let lat = (x[2] / rxy).atan();
    let coslon = x[0] / rxy;
    let sinlon = x[1] / rxy;
    let coslat = rxy / rxyz;
    let sinlat = x[2] / rxyz;
    let xx3 = x[3] * coslon + x[4] * sinlon;
    let xx4 = -x[3] * sinlon + x[4] * coslon;
    let speed_lon = xx4 / rxy;
    let xx4b = -sinlat * xx3 + coslat * x[5];
    let xx5 = coslat * xx3 + sinlat * x[5];
    let speed_lat = xx4b / rxyz;
    [lon, lat, rxyz, speed_lon, speed_lat, xx5]
}

pub fn polar_to_cartesian_with_speed(l: [f64; 6]) -> [f64; 6] {
    if l[3] == 0.0 && l[4] == 0.0 && l[5] == 0.0 {
        let pos = polar_to_cartesian([l[0], l[1], l[2]]);
        return [pos[0], pos[1], pos[2], 0.0, 0.0, 0.0];
    }
    let coslon = l[0].cos();
    let sinlon = l[0].sin();
    let coslat = l[1].cos();
    let sinlat = l[1].sin();
    let xx0 = l[2] * coslat * coslon;
    let xx1 = l[2] * coslat * sinlon;
    let xx2 = l[2] * sinlat;
    let rxyz = l[2];
    let rxy = (xx0 * xx0 + xx1 * xx1).sqrt();
    let xx5 = l[5];
    let xx4 = l[4] * rxyz;
    let x5 = sinlat * xx5 + coslat * xx4;
    let xx3 = coslat * xx5 - sinlat * xx4;
    let xx4b = l[3] * rxy;
    let x3 = coslon * xx3 - sinlon * xx4b;
    let x4 = sinlon * xx3 + coslon * xx4b;
    [xx0, xx1, xx2, x3, x4, x5]
}

// ---------------------------------------------------------------------------
// Coordinate system conversion (ecliptic ↔ equatorial, degree-level wrappers)
// ---------------------------------------------------------------------------

pub fn cotrans(xpo: [f64; 3], eps: f64) -> [f64; 3] {
    let e = eps * DEGTORAD;
    let x = [xpo[0] * DEGTORAD, xpo[1] * DEGTORAD, 1.0];
    let cart = polar_to_cartesian(x);
    let rotated = rotate_x(cart, e);
    let result = cartesian_to_polar(rotated);
    [result[0] * RADTODEG, result[1] * RADTODEG, xpo[2]]
}

pub fn cotrans_with_speed(xpo: [f64; 6], eps: f64) -> [f64; 6] {
    let e = eps * DEGTORAD;
    let x = [
        xpo[0] * DEGTORAD,
        xpo[1] * DEGTORAD,
        1.0,
        xpo[3] * DEGTORAD,
        xpo[4] * DEGTORAD,
        xpo[5],
    ];
    let cart = polar_to_cartesian_with_speed(x);
    let pos = rotate_x([cart[0], cart[1], cart[2]], e);
    let vel = rotate_x([cart[3], cart[4], cart[5]], e);
    let result = cartesian_to_polar_with_speed([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]);
    [
        result[0] * RADTODEG,
        result[1] * RADTODEG,
        xpo[2],
        result[3] * RADTODEG,
        result[4] * RADTODEG,
        xpo[5],
    ]
}

// ---------------------------------------------------------------------------
// Degree splitting
// ---------------------------------------------------------------------------

pub fn split_degrees(mut ddeg: f64, flags: SplitDegFlags) -> DegreeParts {
    let mut sign = 1i32;
    if ddeg < 0.0 {
        sign = -1;
        ddeg = -ddeg;
    }
    // NAKSHATRA: silently ignored (deferred to swisseph-rs/9)
    let mut dadd = if flags.contains(SplitDegFlags::ROUND_DEG) {
        0.5
    } else if flags.contains(SplitDegFlags::ROUND_MIN) {
        0.5 / 60.0
    } else if flags.contains(SplitDegFlags::ROUND_SEC) {
        0.5 / 3600.0
    } else {
        0.0
    };
    if flags.contains(SplitDegFlags::KEEP_DEG) {
        if (ddeg + dadd) as i32 - ddeg as i32 > 0 {
            dadd = 0.0;
        }
    } else if flags.contains(SplitDegFlags::KEEP_SIGN) {
        if ddeg % 30.0 + dadd >= 30.0 {
            dadd = 0.0;
        }
    }
    ddeg += dadd;
    if flags.contains(SplitDegFlags::ZODIACAL) {
        sign = (ddeg / 30.0) as i32;
        if sign == 12 {
            sign = 0;
        }
        ddeg %= 30.0;
    }
    let degrees = ddeg as i32;
    ddeg -= degrees as f64;
    let minutes = (ddeg * 60.0) as i32;
    ddeg -= minutes as f64 / 60.0;
    let seconds = (ddeg * 3600.0) as i32;
    let second_fraction = if !flags
        .intersects(SplitDegFlags::ROUND_DEG | SplitDegFlags::ROUND_MIN | SplitDegFlags::ROUND_SEC)
    {
        ddeg * 3600.0 - seconds as f64
    } else {
        0.0
    };
    DegreeParts {
        degrees,
        minutes: if flags.contains(SplitDegFlags::ROUND_DEG) {
            0
        } else {
            minutes
        },
        seconds: if flags.intersects(SplitDegFlags::ROUND_DEG | SplitDegFlags::ROUND_MIN) {
            0
        } else {
            seconds
        },
        second_fraction,
        sign,
    }
}

pub fn poly_eval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

// ---------------------------------------------------------------------------
// Owen 1990 shared utilities
// ---------------------------------------------------------------------------

pub const OWEN_T0S: [f64; 5] = [-3392455.5, -470455.5, 2451544.5, 5373544.5, 8295544.5];

pub fn owen_t0_icof(jd: f64) -> (f64, usize) {
    let mut t0 = OWEN_T0S[0];
    let mut icof = 0;
    for i in 1..5 {
        if jd >= (OWEN_T0S[i - 1] + OWEN_T0S[i]) / 2.0 {
            t0 = OWEN_T0S[i];
            icof = i;
        }
    }
    (t0, icof)
}

pub fn owen_chebyshev_basis(jd: f64) -> (usize, [f64; 10]) {
    let (t0, icof) = owen_t0_icof(jd);
    let x = (jd - t0) / 36525.0 / 40.0;
    let mut tau = [0.0; 10];
    tau[1] = x;
    for i in 2..=9 {
        tau[i] = x * tau[i - 1];
    }
    let k = [
        1.0,
        tau[1],
        2.0 * tau[2] - 1.0,
        4.0 * tau[3] - 3.0 * tau[1],
        8.0 * tau[4] - 8.0 * tau[2] + 1.0,
        16.0 * tau[5] - 20.0 * tau[3] + 5.0 * tau[1],
        32.0 * tau[6] - 48.0 * tau[4] + 18.0 * tau[2] - 1.0,
        64.0 * tau[7] - 112.0 * tau[5] + 56.0 * tau[3] - 7.0 * tau[1],
        128.0 * tau[8] - 256.0 * tau[6] + 160.0 * tau[4] - 32.0 * tau[2] + 1.0,
        256.0 * tau[9] - 576.0 * tau[7] + 432.0 * tau[5] - 120.0 * tau[3] + 9.0 * tau[1],
    ];
    (icof, k)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_degrees_basics() {
        assert_eq!(normalize_degrees(0.0), 0.0);
        assert_eq!(normalize_degrees(360.0), 0.0);
        assert_eq!(normalize_degrees(-90.0), 270.0);
        assert_eq!(normalize_degrees(720.0), 0.0);
        assert_eq!(normalize_degrees(450.0), 90.0);
    }

    #[test]
    fn near_zero_guard() {
        assert_eq!(normalize_degrees(1e-14), 0.0);
        assert_eq!(normalize_degrees(-1e-14), 0.0);
        assert_eq!(normalize_radians(1e-14), 0.0);
        assert_eq!(normalize_radians(-1e-14), 0.0);
    }

    #[test]
    fn mod_2pi_no_guard() {
        let tiny = 1e-14;
        assert_eq!(mod_2pi(tiny), tiny);
    }

    #[test]
    fn diff_degrees_signed() {
        assert_eq!(diff_degrees(10.0, 350.0), 20.0);
        assert_eq!(diff_degrees(350.0, 10.0), -20.0);
        assert_eq!(diff_degrees(0.0, 0.0), 0.0);
    }

    #[test]
    fn midpoint_degrees_wrap() {
        let mid = midpoint_degrees(350.0, 10.0);
        assert_eq!(mid, 0.0);
    }

    #[test]
    fn csnorm_basics() {
        assert_eq!(csnorm(0), 0);
        assert_eq!(csnorm(DEG360), 0);
        assert_eq!(csnorm(-1), DEG360 - 1);
        assert_eq!(csnorm(DEG360 + 100), 100);
    }

    #[test]
    fn d2l_rounding() {
        assert_eq!(d2l(0.5), 1);
        assert_eq!(d2l(-0.5), -1);
        assert_eq!(d2l(1.5), 2);
        assert_eq!(d2l(-1.5), -2);
        assert_eq!(d2l(0.0), 0);
    }

    #[test]
    fn chebyshev_t0() {
        // T₀(t) = 1: single coeff [2.0] → value should be 1.0
        // Broucke: with coef=[c], br = x2*0 - 0 + c = c, brp2 = 0, result = (c - 0) * 0.5 = c/2
        assert_eq!(chebyshev_eval(0.5, &[2.0]), 1.0);
        assert_eq!(chebyshev_eval(-0.3, &[2.0]), 1.0);
    }

    #[test]
    fn chebyshev_t1() {
        // T₁(t) = t: coeffs [0.0, 1.0]
        // Broucke with these coeffs should give t
        let t = 0.7;
        let result = chebyshev_eval(t, &[0.0, 1.0]);
        assert!(
            (result - t).abs() < 1e-15,
            "T₁({t}) = {result}, expected {t}"
        );
    }

    #[test]
    fn chebyshev_deriv_t2() {
        // T₂(t) = 2t²-1, T₂'(t) = 4t
        // Coeffs for T₂: [-1.0, 0.0, 1.0] (in Chebyshev basis: c₀=-1, c₁=0, c₂=1 → T₂ itself)
        // Wait — actually Broucke uses the convention that sum = Σ cⱼTⱼ with T₀ halved
        // For a pure T₂ polynomial: c₀=0, c₁=0, c₂=1 → function value = T₂(t)
        // Derivative = 4t
        let t = 0.3;
        let result = chebyshev_deriv(t, &[0.0, 0.0, 1.0]);
        let expected = 4.0 * t;
        assert!(
            (result - expected).abs() < 1e-14,
            "T₂'({t}) = {result}, expected {expected}"
        );
    }

    #[test]
    fn cartpol_roundtrip() {
        let original = [1.0, 2.0, 3.0];
        let polar = cartesian_to_polar(original);
        let back = polar_to_cartesian(polar);
        for i in 0..3 {
            assert!(
                (original[i] - back[i]).abs() < 1e-15,
                "component {i}: {} vs {}",
                original[i],
                back[i]
            );
        }
    }

    #[test]
    fn cartpol_zero() {
        assert_eq!(cartesian_to_polar([0.0, 0.0, 0.0]), [0.0; 3]);
    }

    #[test]
    fn cartpol_pole() {
        let result = cartesian_to_polar([0.0, 0.0, 1.0]);
        assert!((result[1] - PI / 2.0).abs() < 1e-15);
        let result_neg = cartesian_to_polar([0.0, 0.0, -1.0]);
        assert!((result_neg[1] + PI / 2.0).abs() < 1e-15);
    }

    #[test]
    fn cartpol_sp_zero_position() {
        let input = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let result = cartesian_to_polar_with_speed(input);
        assert_eq!(result[2], 0.0);
        assert_eq!(result[3], 0.0);
        assert_eq!(result[4], 0.0);
        assert!((result[5] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn cartpol_sp_zero_speed() {
        let input = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        let result = cartesian_to_polar_with_speed(input);
        assert_eq!(result[3], 0.0);
        assert_eq!(result[4], 0.0);
        assert_eq!(result[5], 0.0);
        assert!(result[2] > 0.0);
    }

    #[test]
    fn polcart_sp_roundtrip() {
        let polar = [1.0, 0.5, 2.0, 0.1, 0.05, 0.3];
        let cart = polar_to_cartesian_with_speed(polar);
        let back = cartesian_to_polar_with_speed(cart);
        for i in 0..6 {
            assert!(
                (polar[i] - back[i]).abs() < 1e-14,
                "component {i}: {} vs {}",
                polar[i],
                back[i]
            );
        }
    }

    #[test]
    fn rotate_x_identity() {
        let pos = [1.0, 2.0, 3.0];
        let rotated = rotate_x(pos, 0.0);
        for i in 0..3 {
            assert!(
                (pos[i] - rotated[i]).abs() < 1e-15,
                "component {i}: {} vs {}",
                pos[i],
                rotated[i]
            );
        }
    }

    #[test]
    fn split_degrees_basic() {
        let result = split_degrees(123.456789, SplitDegFlags::empty());
        assert_eq!(result.degrees, 123);
        assert_eq!(result.minutes, 27);
        assert_eq!(result.sign, 1);
    }

    #[test]
    fn split_degrees_negative() {
        let result = split_degrees(-45.5, SplitDegFlags::empty());
        assert_eq!(result.sign, -1);
        assert_eq!(result.degrees, 45);
        assert_eq!(result.minutes, 30);
    }

    #[test]
    fn split_degrees_zodiacal() {
        let result = split_degrees(95.0, SplitDegFlags::ZODIACAL);
        assert_eq!(result.sign, 3);
        assert_eq!(result.degrees, 5);
    }

    #[test]
    fn split_degrees_round_sec() {
        let result = split_degrees(10.999999, SplitDegFlags::ROUND_SEC);
        assert_eq!(result.second_fraction, 0.0);
    }
}
