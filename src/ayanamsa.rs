use crate::constants::{B1950, J1900, J2000};

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
