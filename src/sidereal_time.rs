use crate::config::EphemerisConfig;
use crate::constants::*;
use crate::deltat::calc_deltat;
use crate::flags::CalcFlags;
use crate::math::{cartesian_to_polar, normalize_degrees, polar_to_cartesian, rotate_x};
use crate::nutation;
use crate::obliquity;
use crate::precession;
use crate::types::*;

// 33 pairs (sin, cos) in microarcseconds — swephlib.c:3341–3375
static STCF: [[f64; 2]; 33] = [
    [2640.96, -0.39],
    [63.52, -0.02],
    [11.75, 0.01],
    [11.21, 0.01],
    [-4.55, 0.00],
    [2.02, 0.00],
    [1.98, 0.00],
    [-1.72, 0.00],
    [-1.41, -0.01],
    [-1.26, -0.01],
    [-0.63, 0.00],
    [-0.63, 0.00],
    [0.46, 0.00],
    [0.45, 0.00],
    [0.36, 0.00],
    [-0.24, -0.12],
    [0.32, 0.00],
    [0.28, 0.00],
    [0.27, 0.00],
    [0.26, 0.00],
    [-0.21, 0.00],
    [0.19, 0.00],
    [0.18, 0.00],
    [-0.10, 0.05],
    [0.15, 0.00],
    [-0.14, 0.00],
    [0.14, 0.00],
    [-0.14, 0.00],
    [0.14, 0.00],
    [0.13, 0.00],
    [-0.11, 0.00],
    [0.11, 0.00],
    [0.11, 0.00],
];

// 33×14 integer multipliers for [l, l', F, D, Om, L_Me..L_Ne, p_A]
// swephlib.c:3378–3412
static STFARG: [[i8; 14]; 33] = [
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 4, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, -2, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, -2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, -1],
    [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [2, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 4, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, -2, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
];

// ---------------------------------------------------------------------------
// Complementary terms — 33-term Fourier series (swephlib.c:3413–3450)
// ---------------------------------------------------------------------------

fn complementary_terms(tt: f64) -> f64 {
    let delm = nutation::planetary_args_iau2000(tt);
    let mut dadd = -0.87 * delm[4].sin() * tt;

    for i in 0..33 {
        let mut darg = 0.0_f64;
        for j in 0..14 {
            darg += f64::from(STFARG[i][j]) * delm[j];
        }
        dadd += STCF[i][0] * darg.sin() + STCF[i][1] * darg.cos();
    }

    dadd / (3600.0 * 1_000_000.0)
}

// ---------------------------------------------------------------------------
// GMST models
// ---------------------------------------------------------------------------

fn gmst_iau1976(jd0: f64, secs: f64) -> f64 {
    let tu = (jd0 - J2000) / 36525.0;
    let gmst_0h = ((-6.2e-6 * tu + 9.3104e-2) * tu + 8640184.812866) * tu + 24110.54841;
    let msday = 1.0 + ((-1.86e-5 * tu + 0.186208) * tu + 8640184.812866) / (86400.0 * 36525.0);
    gmst_0h + msday * secs
}

fn gmst_iau2006(jd0: f64, secs: f64, tu: f64, config: &EphemerisConfig) -> f64 {
    let tt = (jd0 + calc_deltat(jd0, config) - J2000) / 36525.0;

    let gmst_0h =
        (((-0.000000002454 * tt - 0.00000199708) * tt - 0.0000002926) * tt + 0.092772110) * tt * tt
            + 307.4771013 * (tt - tu)
            + 8640184.79447825 * tu
            + 24110.5493771;

    let msday = 1.0
        + ((((-0.000000012270 * tt - 0.00000798832) * tt - 0.0000008778) * tt + 0.185544220) * tt
            + 8640184.79447825)
            / (86400.0 * 36525.0);

    gmst_0h + msday * secs
}

// Verbatim ERA / GMST polynomial coefficients from the C source.
#[allow(clippy::excessive_precision)]
fn gmst_era(tjd_ut: f64, config: &EphemerisConfig) -> f64 {
    let jdrel = tjd_ut - J2000;
    let tt = (tjd_ut + calc_deltat(tjd_ut, config) - J2000) / 36525.0;

    let mut gmst = normalize_degrees((0.7790572732640 + 1.00273781191135448 * jdrel) * 360.0);

    gmst += (0.014506
        + tt * (4612.156534
            + tt * (1.3915817 + tt * (-0.00000044 + tt * (-0.000029956 + tt * -0.0000000368)))))
        / 3600.0;

    gmst = normalize_degrees(gmst + complementary_terms(tt));

    gmst / 15.0 * 3600.0
}

// ---------------------------------------------------------------------------
// Long-term model (swephlib.c:3285–3324)
// ---------------------------------------------------------------------------

// eps and nut in DEGREES (matching C convention)
fn sidtime_long_term(tjd_ut: f64, eps: f64, nut: f64, config: &EphemerisConfig) -> f64 {
    let flags = CalcFlags::empty();
    let models = &config.astro_models;
    // C's sidtime_long_term (swephlib.c:3291, 3301) resolves both of its deltaT calls via
    // swe_deltat_ex(tjd, -1, NULL) -- the `-1` sentinel forces SE_TIDAL_DEFAULT, independent of
    // whatever ephemeris source is actually configured (swi_get_tid_acc's iflag=0/denum=9999
    // path falls straight to the `default:` case). Force the same override here — same
    // deltaT-tid_acc-inconsistency pattern already documented for swe_houses_ex2 (see
    // Ephemeris::houses_ex2's `deltat_config`).
    let deltat_config = {
        let mut c = config.clone();
        c.tidal_acceleration = Some(TIDAL_DEFAULT);
        c
    };
    let tjd_et = tjd_ut + calc_deltat(tjd_ut, &deltat_config);
    let t = (tjd_et - J2000) / 365250.0;
    let t2 = t * t;
    let t3 = t * t2;
    let dlt = AUNIT / CLIGHT / 86400.0;

    let dlon = 100.46645683 + (1295977422.83429 * t - 2.04411 * t2 - 0.00523 * t3) / 3600.0;
    let dlon = normalize_degrees(dlon - dlt * 360.0 / 365.2425);

    let mut xs = polar_to_cartesian([dlon * DEGTORAD, 0.0, 1.0]);

    let eps_j2000 = obliquity::obliquity(J2000 + calc_deltat(J2000, &deltat_config), flags, models);
    xs = rotate_x(xs, -eps_j2000.eps);

    precession::precess(
        &mut xs,
        tjd_et,
        flags,
        models,
        PrecessionDirection::J2000ToDate,
    );

    let eps_date = obliquity::obliquity(tjd_et, flags, models);
    let nut_date = nutation::nutation(tjd_et, flags, models);

    xs = rotate_x(xs, eps_date.eps);
    let mut pol = cartesian_to_polar(xs);
    pol[0] *= RADTODEG;

    let dhour = ((tjd_ut - 0.5) % 1.0) * 360.0;

    if eps == 0.0 {
        let eps_true_deg = eps_date.eps * RADTODEG + nut_date.deps * RADTODEG;
        let dpsi_deg = nut_date.dpsi * RADTODEG;
        pol[0] += dpsi_deg * (eps_true_deg * DEGTORAD).cos();
    } else {
        pol[0] += nut * (eps * DEGTORAD).cos();
    }

    normalize_degrees(pol[0] + dhour) / 15.0
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Greenwich Apparent Sidereal Time from pre-computed obliquity and nutation.
///
/// `eps` = true obliquity in **degrees**, `nut` = dpsi in **degrees**.
/// Returns hours in [0, 24).
pub fn sidereal_time0(tjd_ut: f64, eps: f64, nut: f64, config: &EphemerisConfig) -> f64 {
    let sidt_model = config.astro_models.sidereal_time;

    if sidt_model == SiderealTimeModel::Longterm
        && (tjd_ut <= SIDT_LTERM_T0 || tjd_ut >= SIDT_LTERM_T1)
    {
        let mut gmst = sidtime_long_term(tjd_ut, eps, nut, config);
        if tjd_ut <= SIDT_LTERM_T0 {
            gmst -= SIDT_LTERM_OFS0;
        } else {
            gmst -= SIDT_LTERM_OFS1;
        }
        if gmst >= 24.0 {
            gmst -= 24.0;
        }
        if gmst < 0.0 {
            gmst += 24.0;
        }
        return gmst;
    }

    let mut jd0 = tjd_ut.floor();
    let mut secs = tjd_ut - jd0;
    if secs < 0.5 {
        jd0 -= 0.5;
        secs += 0.5;
    } else {
        jd0 += 0.5;
        secs -= 0.5;
    }
    secs *= 86400.0;
    let tu = (jd0 - J2000) / 36525.0;

    let mut gmst = match sidt_model {
        SiderealTimeModel::IersConv2010 | SiderealTimeModel::Longterm => gmst_era(tjd_ut, config),
        SiderealTimeModel::IAU2006 => gmst_iau2006(jd0, secs, tu, config),
        SiderealTimeModel::IAU1976 => gmst_iau1976(jd0, secs),
    };

    let eqeq = 240.0 * nut * (eps * DEGTORAD).cos();
    gmst += eqeq;

    gmst -= 86400.0 * (gmst / 86400.0).floor();
    gmst /= 3600.0;

    if gmst >= 24.0 {
        gmst = 0.0;
    }

    gmst
}

/// Greenwich Apparent Sidereal Time with automatic obliquity/nutation.
///
/// Returns hours in [0, 24).
pub fn sidereal_time(tjd_ut: f64, config: &EphemerisConfig) -> f64 {
    let flags = CalcFlags::empty();
    let models = &config.astro_models;

    let tjd_et = tjd_ut + calc_deltat(tjd_ut, config);

    let eps = obliquity::obliquity(tjd_et, flags, models);
    let nut = nutation::nutation(tjd_et, flags, models);

    let eps_deg = eps.eps * RADTODEG;
    let dpsi_deg = nut.dpsi * RADTODEG;
    let deps_deg = nut.deps * RADTODEG;

    sidereal_time0(tjd_ut, eps_deg + deps_deg, dpsi_deg, config)
}
