use serde::Deserialize;
use swisseph::flags::{CalcFlags, HeliacalFlags, VisLimFlags};
use swisseph::heliacal::HeliacalEventType;
use swisseph::{Ephemeris, EphemerisConfig};

#[derive(Deserialize)]
struct VisLimitCase {
    tjd_ut: f64,
    object: String,
    helflag: u32,
    flag_desc: String,
    epheflag: u32,
    retval: i32,
    dret: [f64; 8],
}

#[derive(Deserialize)]
struct ArcVisCase {
    tjd_ut: f64,
    mag: f64,
    azi_obj: f64,
    alt_obj: f64,
    azi_sun: f64,
    azi_moon: f64,
    alt_moon: f64,
    helflag: u32,
    #[allow(dead_code)]
    retval: i32,
    dret: f64,
}

#[derive(Deserialize)]
struct HelAngleCase {
    tjd_ut: f64,
    mag: f64,
    azi_obj: f64,
    azi_sun: f64,
    azi_moon: f64,
    alt_moon: f64,
    helflag: u32,
    #[allow(dead_code)]
    retval: i32,
    dret: [f64; 3],
}

#[derive(Deserialize)]
struct PhenoCase {
    tjd_ut: f64,
    object: String,
    type_event: i32,
    geo: [f64; 3],
    helflag: u32,
    desc: String,
    #[allow(dead_code)]
    retval: i32,
    darr: [f64; 28],
}

#[derive(Deserialize)]
struct GoldenData {
    vis_limit: Vec<VisLimitCase>,
    arcvis: Vec<ArcVisCase>,
    helangle: Vec<HelAngleCase>,
    pheno: Vec<PhenoCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("heliacal.json");
    let data = std::fs::read_to_string(path).expect("read heliacal.json");
    serde_json::from_str(&data).expect("parse heliacal.json")
}

fn make_eph() -> Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: swisseph::EphemerisSource::Swiss,
        ephe_path: Some("../swisseph/ephe".into()),
        ..Default::default()
    };
    Ephemeris::new(config).unwrap()
}

fn make_eph_moseph() -> Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: swisseph::EphemerisSource::Moshier,
        ..Default::default()
    };
    Ephemeris::new(config).unwrap()
}

fn retval_to_vision_flags(retval: i32) -> (VisLimFlags, bool) {
    if retval == -2 {
        return (VisLimFlags::empty(), true);
    }
    let mut flags = VisLimFlags::empty();
    if retval & 1 != 0 {
        flags |= VisLimFlags::SCOTOPIC;
    }
    if retval & 2 != 0 {
        flags |= VisLimFlags::MIXED;
    }
    (flags, false)
}

#[test]
fn golden_vis_limit_mag() {
    let data = load();
    let eph_swi = make_eph();
    let eph_mos = make_eph_moseph();

    let eps = 1e-7;

    for (i, case) in data.vis_limit.iter().enumerate() {
        let epheflag = CalcFlags::from_bits_truncate(case.epheflag);
        let helflag = HeliacalFlags::from_bits_truncate(case.helflag);
        let eph = if epheflag.contains(CalcFlags::SWIEPH) {
            &eph_swi
        } else {
            &eph_mos
        };

        let dgeo = [31.25, 30.1, 30.0];
        let mut datm = [1013.25, 15.0, 40.0, 40.0];
        let mut dobs = if helflag.contains(HeliacalFlags::OPTICAL_PARAMS) {
            [36.0, 1.0, 1.0, 10.0, 50.0, 0.8]
        } else {
            [36.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        };

        let result = eph
            .vis_limit_mag(
                case.tjd_ut,
                &dgeo,
                &mut datm,
                &mut dobs,
                &case.object,
                epheflag,
                helflag,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "case {i}: {} obj={} jd={:.4} flags={}: {e}",
                    case.flag_desc, case.object, case.tjd_ut, case.helflag
                )
            });

        let (expected_vision, expected_below) = retval_to_vision_flags(case.retval);

        assert_eq!(
            result.below_horizon,
            expected_below,
            "case {i} ({} obj={} jd={:.4}): below_horizon mismatch: got {}, expected {} (retval={})",
            case.flag_desc,
            case.object,
            case.tjd_ut,
            result.below_horizon,
            expected_below,
            case.retval
        );

        if result.below_horizon {
            super::assert_f64_eps(
                &format!("case {i} below_horizon: limiting_magnitude"),
                case.dret[0],
                result.limiting_magnitude,
                1e-10,
            );
            continue;
        }

        assert_eq!(
            result.vision, expected_vision,
            "case {i} ({} obj={} jd={:.4}): vision flags mismatch: got {:?}, expected {:?} (retval={})",
            case.flag_desc, case.object, case.tjd_ut, result.vision, expected_vision, case.retval
        );

        // Limiting magnitude (dret[0]) sits at the end of the chain
        // azalt→refraction→brightness→extinction→optics→vis_lim_magn,
        // compounding FP drift up to ~3e-7. Positions are tighter (direct from object_loc).
        let eps_lim_mag = 5e-7;
        let fields: &[(&str, f64, f64, f64)] = &[
            (
                "dret[0] limiting_magnitude",
                case.dret[0],
                result.limiting_magnitude,
                eps_lim_mag,
            ),
            (
                "dret[1] altitude_object",
                case.dret[1],
                result.altitude_object,
                eps,
            ),
            (
                "dret[2] azimuth_object",
                case.dret[2],
                result.azimuth_object,
                eps,
            ),
            (
                "dret[3] altitude_sun",
                case.dret[3],
                result.altitude_sun,
                eps,
            ),
            ("dret[4] azimuth_sun", case.dret[4], result.azimuth_sun, eps),
            (
                "dret[5] altitude_moon",
                case.dret[5],
                result.altitude_moon,
                eps,
            ),
            (
                "dret[6] azimuth_moon",
                case.dret[6],
                result.azimuth_moon,
                eps,
            ),
            (
                "dret[7] magnitude_object",
                case.dret[7],
                result.magnitude_object,
                eps,
            ),
        ];
        for (name, expected, actual, field_eps) in fields {
            super::assert_f64_eps(
                &format!(
                    "case {i} ({} obj={} jd={:.4}): {name}",
                    case.flag_desc, case.object, case.tjd_ut
                ),
                *expected,
                *actual,
                *field_eps,
            );
        }
    }
}

#[test]
fn golden_topo_arcus_visionis() {
    let data = load();
    let eph = make_eph();
    let eps = 1e-6;

    for (i, case) in data.arcvis.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(case.helflag);
        let dgeo = [31.25, 30.1, 30.0];
        let mut datm = [1013.25, 15.0, 40.0, 40.0];
        let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];

        let result = eph
            .topo_arcus_visionis(
                case.tjd_ut,
                &dgeo,
                &mut datm,
                &mut dobs,
                helflag,
                case.mag,
                case.azi_obj,
                case.alt_obj,
                case.azi_sun,
                case.azi_moon,
                case.alt_moon,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "arcvis case {i}: mag={:.2} alt_obj={:.1} azi_obj={:.1} jd={:.4}: {e}",
                    case.mag, case.alt_obj, case.azi_obj, case.tjd_ut
                )
            });

        super::assert_f64_eps(
            &format!(
                "arcvis case {i}: mag={:.2} alt_obj={:.1} azi_obj={:.1} jd={:.4}",
                case.mag, case.alt_obj, case.azi_obj, case.tjd_ut
            ),
            case.dret,
            result,
            eps,
        );
    }
}

#[test]
fn golden_heliacal_angle() {
    let data = load();
    let eph = make_eph();
    let eps = 1e-6;

    for (i, case) in data.helangle.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(case.helflag);
        let dgeo = [31.25, 30.1, 30.0];
        let mut datm = [1013.25, 15.0, 40.0, 40.0];
        let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];

        let result = eph
            .heliacal_angle(
                case.tjd_ut,
                &dgeo,
                &mut datm,
                &mut dobs,
                helflag,
                case.mag,
                case.azi_obj,
                case.azi_sun,
                case.azi_moon,
                case.alt_moon,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "helangle case {i}: mag={:.2} azi_obj={:.1} jd={:.4}: {e}",
                    case.mag, case.azi_obj, case.tjd_ut
                )
            });

        let fields: &[(&str, f64, f64)] = &[
            ("optimal_altitude", case.dret[0], result.optimal_altitude),
            ("arcus_visionis", case.dret[1], result.arcus_visionis),
            ("sun_altitude_diff", case.dret[2], result.sun_altitude_diff),
        ];
        for (name, expected, actual) in fields {
            super::assert_f64_eps(
                &format!(
                    "helangle case {i}: mag={:.2} azi_obj={:.1} jd={:.4}: {name}",
                    case.mag, case.azi_obj, case.tjd_ut
                ),
                *expected,
                *actual,
                eps,
            );
        }
    }
}

#[test]
fn golden_heliacal_pheno_ut() {
    let data = load();
    let eph_swi = make_eph();
    let eph_mos = make_eph_moseph();

    let slot_names = [
        "AltO", "AppAltO", "GeoAltO", "AziO", "AltS", "AziS", "TAVact", "ARCVact", "DAZact",
        "ARCLact", "kact", "MinTAV", "TfirstVR", "TbVR", "TlastVR", "TbYallop", "WMoon", "qYal",
        "qCrit", "ParO", "MagnO", "RiseO", "RiseS", "Lag", "TvisVR", "LMoon", "elong", "illum",
    ];

    for (i, case) in data.pheno.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(case.helflag);
        let epheflag = CalcFlags::from_bits_truncate(case.helflag & 0x7);
        let eph = if epheflag.contains(CalcFlags::SWIEPH) {
            &eph_swi
        } else {
            &eph_mos
        };

        let event = HeliacalEventType::try_from(case.type_event).unwrap_or_else(|e| {
            panic!(
                "pheno case {i} ({}): bad type_event {}: {e}",
                case.desc, case.type_event
            )
        });

        let mut datm = [1013.25, 15.0, 40.0, 40.0];
        let mut dobs = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];

        let result = eph
            .heliacal_pheno_ut(
                case.tjd_ut,
                &case.geo,
                &mut datm,
                &mut dobs,
                &case.object,
                event,
                epheflag,
                helflag,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "pheno case {i} ({}): obj={} jd={:.4}: {e}",
                    case.desc, case.object, case.tjd_ut
                )
            });

        let actual = result.as_array();
        for (slot, (expected, got)) in case.darr.iter().zip(actual.iter()).enumerate() {
            let eps =
                if slot >= 12 && slot <= 15 || slot == 21 || slot == 22 || slot == 23 || slot == 24
                {
                    // Time/duration slots (from rise/set searches or crossing/parabola fits)
                    1e-5
                } else {
                    // Instantaneous geometry slots
                    1e-7
                };
            super::assert_f64_eps(
                &format!(
                    "pheno case {i} ({}) slot {slot} ({}): obj={} jd={:.4}",
                    case.desc, slot_names[slot], case.object, case.tjd_ut
                ),
                *expected,
                *got,
                eps,
            );
        }
    }
}
