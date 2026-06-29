use crate::constants::{
    B1950, DEGTORAD, J1900, J2000, RADTODEG, SSY_PLANE_INCL, SSY_PLANE_NODE_E2000,
};
use crate::context::EphemerisConfig;
use crate::error::Error;
use crate::flags::{CalcFlags, SiderealBits};
use crate::math::{
    cartesian_to_polar, cartesian_to_polar_with_speed, normalize_degrees, normalize_radians,
    polar_to_cartesian, polar_to_cartesian_with_speed, rotate_x, rotate_x_sincos,
};
use crate::types::{AstroModels, PrecessionDirection, SiderealMode};

#[derive(Debug, Clone, Copy)]
pub(crate) struct AyaInit {
    pub t0: f64,
    pub ayan_t0: f64,
    pub t0_is_ut: bool,
    pub prec_offset: i32,
}

pub(crate) const AYANAMSA: [AyaInit; 47] = [
    AyaInit {
        t0: 2433282.42346,
        ayan_t0: 24.042044444,
        t0_is_ut: false,
        prec_offset: 11,
    }, //  0 FAGAN_BRADLEY
    AyaInit {
        t0: 2435553.5,
        ayan_t0: 23.250182778 - 0.004658035,
        t0_is_ut: false,
        prec_offset: 1,
    }, //  1 LAHIRI
    AyaInit {
        t0: 1721057.5,
        ayan_t0: 0.0,
        t0_is_ut: true,
        prec_offset: 0,
    }, //  2 DELUCE
    AyaInit {
        t0: J1900,
        ayan_t0: 360.0 - 338.98556,
        t0_is_ut: false,
        prec_offset: 11,
    }, //  3 RAMAN
    AyaInit {
        t0: J1900,
        ayan_t0: 360.0 - 341.33904,
        t0_is_ut: false,
        prec_offset: -1,
    }, //  4 USHASHASHI
    AyaInit {
        t0: J1900,
        ayan_t0: 360.0 - 337.636111,
        t0_is_ut: false,
        prec_offset: 11,
    }, //  5 KRISHNAMURTI
    AyaInit {
        t0: J1900,
        ayan_t0: 360.0 - 333.0369024,
        t0_is_ut: false,
        prec_offset: 0,
    }, //  6 DJWHAL_KHUL
    AyaInit {
        t0: J1900,
        ayan_t0: 360.0 - 338.917778,
        t0_is_ut: false,
        prec_offset: -1,
    }, //  7 YUKTESHWAR
    AyaInit {
        t0: J1900,
        ayan_t0: 360.0 - 338.634444,
        t0_is_ut: false,
        prec_offset: -1,
    }, //  8 JN_BHASIN
    AyaInit {
        t0: 1684532.5,
        ayan_t0: -5.66667,
        t0_is_ut: true,
        prec_offset: -1,
    }, //  9 BABYL_KUGLER1
    AyaInit {
        t0: 1684532.5,
        ayan_t0: -4.26667,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 10 BABYL_KUGLER2
    AyaInit {
        t0: 1684532.5,
        ayan_t0: -3.41667,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 11 BABYL_KUGLER3
    AyaInit {
        t0: 1684532.5,
        ayan_t0: -4.46667,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 12 BABYL_HUBER
    AyaInit {
        t0: 1673941.0,
        ayan_t0: -5.079167,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 13 BABYL_ETPSC
    AyaInit {
        t0: 1684532.5,
        ayan_t0: -4.44138598,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 14 ALDEBARAN_15TAU
    AyaInit {
        t0: 1674484.0,
        ayan_t0: -9.33333,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 15 HIPPARCHOS
    AyaInit {
        t0: 1927135.8747793,
        ayan_t0: 0.0,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 16 SASSANIAN
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 17 GALCENT_0SAG (fixed-star)
    AyaInit {
        t0: J2000,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 18 J2000
    AyaInit {
        t0: J1900,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 19 J1900
    AyaInit {
        t0: B1950,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 20 B1950
    AyaInit {
        t0: 1903396.8128654,
        ayan_t0: 0.0,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 21 SURYASIDDHANTA
    AyaInit {
        t0: 1903396.8128654,
        ayan_t0: -0.21463395,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 22 SURYASIDDHANTA_MSUN
    AyaInit {
        t0: 1903396.7895321,
        ayan_t0: 0.0,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 23 ARYABHATA
    AyaInit {
        t0: 1903396.7895321,
        ayan_t0: -0.23763238,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 24 ARYABHATA_MSUN
    AyaInit {
        t0: 1903396.8128654,
        ayan_t0: -0.79167046,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 25 SS_REVATI
    AyaInit {
        t0: 1903396.8128654,
        ayan_t0: 2.11070444,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 26 SS_CITRA
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 27 TRUE_CITRA (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 28 TRUE_REVATI (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 29 TRUE_PUSHYA (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 30 GALCENT_RGILBRAND (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 31 GALEQU_IAU1958 (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 32 GALEQU_TRUE (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 33 GALEQU_MULA (fixed-star)
    AyaInit {
        t0: 2451079.734892,
        ayan_t0: 30.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 34 GALALIGN_MARDYKS
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 35 TRUE_MULA (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 36 GALCENT_MULA_WILHELM (fixed-star)
    AyaInit {
        t0: 1911797.740782065,
        ayan_t0: 0.0,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 37 ARYABHATA_522
    AyaInit {
        t0: 1721057.5,
        ayan_t0: -3.2,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 38 BABYL_BRITTON
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 39 TRUE_SHEORAN (fixed-star)
    AyaInit {
        t0: 0.0,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 40 GALCENT_COCHRANE (fixed-star)
    AyaInit {
        t0: 2451544.5,
        ayan_t0: 25.0,
        t0_is_ut: true,
        prec_offset: 0,
    }, // 41 GALEQU_FIORENZA
    AyaInit {
        t0: 1775845.5,
        ayan_t0: -2.9422,
        t0_is_ut: true,
        prec_offset: -1,
    }, // 42 VALENS_MOON
    AyaInit {
        t0: J1900,
        ayan_t0: 22.44597222,
        t0_is_ut: false,
        prec_offset: 11,
    }, // 43 LAHIRI_1940
    AyaInit {
        t0: 1825235.2458513028,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 44 LAHIRI_VP285
    AyaInit {
        t0: 1827424.752255678,
        ayan_t0: 0.0,
        t0_is_ut: false,
        prec_offset: 0,
    }, // 45 KRISHNAMURTI_VP291
    AyaInit {
        t0: 2435553.5,
        ayan_t0: 23.25 - 0.00464207,
        t0_is_ut: false,
        prec_offset: 11,
    }, // 46 LAHIRI_ICRC
];

#[allow(dead_code)]
pub(crate) fn aya_init(index: usize) -> AyaInit {
    AYANAMSA[index]
}

// ---------------------------------------------------------------------------
// Fixed-star ayanamsa indices: require star catalog, deferred.
// ---------------------------------------------------------------------------
const FIXED_STAR_INDICES: [usize; 12] = [17, 27, 28, 29, 30, 31, 32, 33, 35, 36, 39, 40];

fn sidereal_index(config: &EphemerisConfig) -> usize {
    match config.sidereal_mode {
        Some(mode) => mode as i32 as usize,
        None => 0, // FaganBradley
    }
}

pub(crate) fn resolve_t0(config: &EphemerisConfig, flags: CalcFlags) -> f64 {
    let mut t0 = config.sidereal_t0;
    if config.sidereal_t0_is_ut {
        t0 += crate::deltat::calc_deltat(t0, config);
    }
    let _ = flags; // passed for API parity with C; doesn't affect precession output
    t0
}

/// Core ayanamsa computation — no nutation. Matches `swi_get_ayanamsa_ex`.
pub fn get_ayanamsa_ex(
    config: &EphemerisConfig,
    jd_tt: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<f64, Error> {
    let idx = sidereal_index(config);

    if FIXED_STAR_INDICES.contains(&idx) {
        let mode = config.sidereal_mode.unwrap_or(SiderealMode::FaganBradley);
        return Err(Error::SiderealModeRequiresFixedStars(mode));
    }

    let t0 = resolve_t0(config, flags);
    let ayan_t0 = config.sidereal_ayan_t0;

    let lon = if config.sidereal_bits.contains(SiderealBits::ECL_DATE) {
        // Method 2: propagate the zero-point through ecliptic-of-date
        let mut x = polar_to_cartesian([normalize_degrees(ayan_t0) * DEGTORAD, 0.0, 1.0]);
        let eps_t0 = crate::obliquity::obliquity(t0, CalcFlags::empty(), models).eps;
        x = rotate_x(x, -eps_t0); // ecliptic of t0 → equatorial of t0
        if t0 != J2000 {
            crate::precession::precess(&mut x, t0, flags, models, PrecessionDirection::DateToJ2000);
        }
        crate::precession::precess(
            &mut x,
            jd_tt,
            flags,
            models,
            PrecessionDirection::J2000ToDate,
        );
        let eps_d = crate::obliquity::obliquity(jd_tt, CalcFlags::empty(), models).eps;
        x = rotate_x(x, eps_d); // equatorial of date → ecliptic of date
        let polar = cartesian_to_polar(x);
        normalize_degrees(polar[0] * RADTODEG)
    } else {
        // Method 1: precess vernal point at date back to t0
        let mut x = [1.0_f64, 0.0, 0.0];
        if jd_tt != J2000 {
            crate::precession::precess(
                &mut x,
                jd_tt,
                flags,
                models,
                PrecessionDirection::DateToJ2000,
            );
        }
        crate::precession::precess(&mut x, t0, flags, models, PrecessionDirection::J2000ToDate);
        let eps_t0 = crate::obliquity::obliquity(t0, CalcFlags::empty(), models).eps;
        x = rotate_x(x, eps_t0); // equatorial of t0 → ecliptic of t0
        let polar = cartesian_to_polar(x);
        // FP note 1: write as negation-then-add, not subtract.
        -polar[0] * RADTODEG + ayan_t0
    };

    let corr = get_aya_correction(config, flags, models);
    Ok(normalize_degrees(lon - corr))
}

/// Precession-model correction. Returns 0 when not applicable.
pub fn get_aya_correction(config: &EphemerisConfig, flags: CalcFlags, models: &AstroModels) -> f64 {
    let idx = sidereal_index(config);
    let prec_offset = if idx == 255 {
        0
    } else {
        AYANAMSA[idx].prec_offset
    };
    let prec_model = models.prec_longterm as i32;

    if config.sidereal_t0 == J2000
        || config.sidereal_bits.contains(SiderealBits::NO_PREC_OFFSET)
        || prec_offset == 0
        || prec_offset < 0
        || prec_model == prec_offset
    {
        return 0.0;
    }

    let t0 = resolve_t0(config, flags);
    let mut x = [1.0_f64, 0.0, 0.0];

    // Precess t0→J2000 with current model
    crate::precession::precess(&mut x, t0, flags, models, PrecessionDirection::DateToJ2000);

    // Precess J2000→t0 with the ayanamsa's original model
    let orig_prec = match prec_offset {
        1 => crate::types::PrecessionModel::IAU1976,
        11 => crate::types::PrecessionModel::Newcomb,
        _ => unreachable!("prec_offset {prec_offset} not reachable after guards"),
    };
    let models_orig = AstroModels {
        prec_longterm: orig_prec,
        prec_shortterm: orig_prec,
        ..*models
    };
    crate::precession::precess(
        &mut x,
        t0,
        flags,
        &models_orig,
        PrecessionDirection::J2000ToDate,
    );

    let eps = crate::obliquity::obliquity(t0, CalcFlags::empty(), models).eps;
    x = rotate_x(x, eps); // equatorial → ecliptic
    let polar = cartesian_to_polar(x);
    let mut corr = polar[0] * RADTODEG;
    if corr > 350.0 {
        corr -= 360.0;
    }
    corr
}

/// Public ayanamsa with nutation added (unless NONUT). Matches `swe_get_ayanamsa_ex`.
pub fn get_ayanamsa_ex_nut(
    config: &EphemerisConfig,
    jd_tt: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<f64, Error> {
    let mut daya = get_ayanamsa_ex(config, jd_tt, flags, models)?;
    if !flags.contains(CalcFlags::NONUT) {
        daya += crate::nutation::nutation(jd_tt, flags, models).dpsi * RADTODEG;
    }
    Ok(daya)
}

/// Two-point numerical ayanamsa speed derivative. Matches `swi_get_ayanamsa_with_speed`.
/// Returns `[ayanamsa_deg, speed_deg_per_day]`.
pub fn get_ayanamsa_with_speed(
    config: &EphemerisConfig,
    jd_tt: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<[f64; 2], Error> {
    const TINTV: f64 = 0.001;
    let d0 = get_ayanamsa_ex(config, jd_tt, flags, models)?;
    let d2 = get_ayanamsa_ex(config, jd_tt - TINTV, flags, models)?;
    Ok([d0, (d0 - d2) / TINTV])
}

/// ECL_T0 projection: project body from J2000 equatorial onto ecliptic of t0.
///
/// Returns `(xecl, xequ)` — ecliptic-sidereal Cartesian and equatorial-sidereal Cartesian,
/// both in the frame of epoch t0. Matches C `swi_trop_ra2sid_lon`.
pub(crate) fn trop_ra2sid_lon(
    x2000: &[f64; 6],
    config: &EphemerisConfig,
    models: &AstroModels,
    flags: CalcFlags,
) -> ([f64; 6], [f64; 6]) {
    let mut x = *x2000;
    let t0 = resolve_t0(config, flags);

    // Step 1: precess J2000 → t0 (position and speed as SEPARATE calls, matching C)
    if config.sidereal_t0 != J2000 {
        let mut pos3 = [x[0], x[1], x[2]];
        crate::precession::precess(
            &mut pos3,
            t0,
            flags,
            models,
            PrecessionDirection::J2000ToDate,
        );
        x[0] = pos3[0];
        x[1] = pos3[1];
        x[2] = pos3[2];
        if flags.contains(CalcFlags::SPEED) {
            let mut vel3 = [x[3], x[4], x[5]];
            crate::precession::precess(
                &mut vel3,
                t0,
                flags,
                models,
                PrecessionDirection::J2000ToDate,
            );
            x[3] = vel3[0];
            x[4] = vel3[1];
            x[5] = vel3[2];
        }
    }

    let xequ = x; // equatorial sidereal (frame of t0)

    // Step 2: equatorial t0 → ecliptic t0 (C passes iflag here, not 0)
    let oe = crate::obliquity::obliquity(t0, flags, models);
    let pos3 = rotate_x_sincos([x[0], x[1], x[2]], oe.sin_eps, oe.cos_eps);
    x[0] = pos3[0];
    x[1] = pos3[1];
    x[2] = pos3[2];
    if flags.contains(CalcFlags::SPEED) {
        let vel3 = rotate_x_sincos([x[3], x[4], x[5]], oe.sin_eps, oe.cos_eps);
        x[3] = vel3[0];
        x[4] = vel3[1];
        x[5] = vel3[2];
    }

    // Step 3: Cartesian ecliptic → polar
    let mut pol = cartesian_to_polar_with_speed(x);

    // Step 4: subtract ayanamsa (in RADIANS) and apply correction
    let corr = get_aya_correction(config, flags, models);
    pol[0] -= config.sidereal_ayan_t0 * DEGTORAD;
    pol[0] = normalize_radians(pol[0] + corr * DEGTORAD);

    // Step 5: back to Cartesian
    let xecl = polar_to_cartesian_with_speed(pol);

    (xecl, xequ)
}

/// SSY_PLANE projection: project body onto solar-system invariable plane.
///
/// Returns ecliptic-sidereal Cartesian in the SSY plane. Matches C `swi_trop_ra2sid_lon_sosy`.
pub(crate) fn trop_ra2sid_lon_sosy(
    x2000: &[f64; 6],
    config: &EphemerisConfig,
    models: &AstroModels,
    flags: CalcFlags,
) -> [f64; 6] {
    let oe = crate::obliquity::obliquity(J2000, flags, models);

    // === Planet path ===
    let mut x = *x2000;

    // (a) equatorial J2000 → ecliptic J2000
    let pos3 = rotate_x_sincos([x[0], x[1], x[2]], oe.sin_eps, oe.cos_eps);
    x[0] = pos3[0];
    x[1] = pos3[1];
    x[2] = pos3[2];
    if flags.contains(CalcFlags::SPEED) {
        let vel3 = rotate_x_sincos([x[3], x[4], x[5]], oe.sin_eps, oe.cos_eps);
        x[3] = vel3[0];
        x[4] = vel3[1];
        x[5] = vel3[2];
    }

    // (b) ecliptic Cartesian → polar
    let mut xpol = cartesian_to_polar_with_speed(x);

    // (c) longitude shift by -plane_node
    xpol[0] -= SSY_PLANE_NODE_E2000;

    // (d) polar → Cartesian, then tilt to SSY plane
    let mut xcart = polar_to_cartesian_with_speed(xpol);
    let pos3 = rotate_x([xcart[0], xcart[1], xcart[2]], SSY_PLANE_INCL);
    xcart[0] = pos3[0];
    xcart[1] = pos3[1];
    xcart[2] = pos3[2];
    if flags.contains(CalcFlags::SPEED) {
        let vel3 = rotate_x([xcart[3], xcart[4], xcart[5]], SSY_PLANE_INCL);
        xcart[3] = vel3[0];
        xcart[4] = vel3[1];
        xcart[5] = vel3[2];
    }

    // (e) Cartesian → polar in SSY plane
    let mut x = cartesian_to_polar_with_speed(xcart);

    // === Zero-point path (vernal point of t0 in SSY plane) ===
    let mut x0 = [1.0_f64, 0.0, 0.0];
    let t0 = resolve_t0(config, flags);

    if config.sidereal_t0 != J2000 {
        crate::precession::precess(&mut x0, t0, flags, models, PrecessionDirection::DateToJ2000);
    }

    // equatorial J2000 → ecliptic J2000
    x0 = rotate_x_sincos(x0, oe.sin_eps, oe.cos_eps);
    let mut x0pol = cartesian_to_polar(x0);
    x0pol[0] -= SSY_PLANE_NODE_E2000;
    let mut x0cart = polar_to_cartesian(x0pol);
    x0cart = rotate_x(x0cart, SSY_PLANE_INCL);
    let x0pol2 = cartesian_to_polar(x0cart);

    // === Measure planet relative to zero point (work in DEGREES) ===
    x[0] -= x0pol2[0]; // angle difference in radians
    x[0] *= RADTODEG; // now in degrees

    let corr = get_aya_correction(config, flags, models);
    x[0] -= config.sidereal_ayan_t0;
    x[0] = normalize_degrees(x[0] + corr) * DEGTORAD; // normalize in degrees then back to radians

    polar_to_cartesian_with_speed(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::EphemerisConfig;
    use crate::flags::SiderealBits;

    #[test]
    fn table_length() {
        assert_eq!(AYANAMSA.len(), 47);
    }

    #[test]
    fn spot_check_values() {
        assert!((AYANAMSA[0].ayan_t0 - 24.042044444).abs() < 1e-9);
        assert_eq!(AYANAMSA[0].prec_offset, 11);
        assert!(AYANAMSA[9].t0_is_ut);
        assert!((AYANAMSA[1].ayan_t0 - 23.245524743).abs() < 1e-9);
    }

    #[test]
    fn set_sidereal_mode_j2000() {
        let mut cfg = EphemerisConfig::default();
        cfg.set_sidereal_mode(18, 0.0, 0.0);
        assert!(cfg.sidereal_bits.contains(SiderealBits::ECL_T0));
        assert_eq!(cfg.sidereal_t0, crate::constants::J2000);
    }

    #[test]
    fn set_sidereal_mode_true_citra_strips_bits() {
        let mut cfg = EphemerisConfig::default();
        cfg.set_sidereal_mode(27 | 256, 0.0, 0.0);
        assert!(cfg.sidereal_bits.is_empty());
    }

    #[test]
    fn set_sidereal_mode_user_ut() {
        let mut cfg = EphemerisConfig::default();
        cfg.set_sidereal_mode(255 | 1024, 1900000.0, 5.0);
        assert_eq!(cfg.sidereal_t0, 1900000.0);
        assert_eq!(cfg.sidereal_ayan_t0, 5.0);
        assert!(cfg.sidereal_t0_is_ut);
    }
}
