use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource, TopoPosition};

#[derive(Deserialize)]
struct CalcCase {
    body: i32,
    #[serde(default)]
    #[allow(dead_code)]
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    retflag: u32,
    output: [f64; 6],
}

fn load(name: &str) -> Vec<CalcCase> {
    let path = super::golden_data_path(name);
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
}

fn body_from_c_id(id: i32) -> Body {
    match id {
        15 => Body::Chiron,
        16 => Body::Pholus,
        17 => Body::Ceres,
        18 => Body::Pallas,
        19 => Body::Juno,
        20 => Body::Vesta,
        _ => panic!("unexpected body id {id}"),
    }
}

fn check_cases<'a>(
    cases: &[CalcCase],
    eph_fn: impl Fn(&CalcCase) -> &'a Ephemeris,
    body_fn: impl Fn(i32) -> Body,
    check_retflag: bool,
) {
    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_fn(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let is_jpl = flags.contains(CalcFlags::JPLEPH);
        let eph = eph_fn(c);

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "case {i} body={} jd={:.1} {}: error: {e}",
                    c.body, c.jd, c.flag_name
                ));
                continue;
            }
        };

        let label = format!("case {i} body={} jd={:.1} {}", c.body, c.jd, c.flag_name);

        if check_retflag && !is_jpl {
            let retflag_expected = CalcFlags::from_bits_truncate(c.retflag);
            let retflag_mask =
                CalcFlags::SWIEPH | CalcFlags::MOSEPH | CalcFlags::SPEED | CalcFlags::HELCTR;
            if result.flags_used & retflag_mask != retflag_expected & retflag_mask {
                failures.push(format!(
                    "{label}: retflag mismatch: expected {:?}, got {:?}",
                    retflag_expected & retflag_mask,
                    result.flags_used & retflag_mask,
                ));
            }
        }

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            let eps = if is_jpl {
                if k >= 3 { 1e-5 } else { 2e-6 }
            } else if k >= 3 {
                1e-7
            } else {
                1e-9
            };
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                    c.output[k], result.data[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(200) {
            eprintln!("{f}");
        }
        panic!("{n} element failures (showing first 200)");
    }
}

#[test]
fn golden_asteroid() {
    let topo = TopoPosition {
        longitude: 8.55,
        latitude: 47.37,
        altitude: 500.0,
    };

    let eph_sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        topographic: Some(topo),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let eph_jpl = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Jpl,
        ephe_path: Some(ephe_path()),
        topographic: Some(topo),
        ..EphemerisConfig::default()
    })
    .expect("JPL ephemeris required for asteroid golden tests (de441.eph in ephe/)");

    let cases = load("asteroid.json");
    assert!(
        cases.len() >= 300,
        "expected 300+ cases, got {}",
        cases.len()
    );

    check_cases(
        &cases,
        |c| {
            let flags = CalcFlags::from_bits_truncate(c.flags);
            if flags.contains(CalcFlags::JPLEPH) {
                &eph_jpl
            } else {
                &eph_sweph
            }
        },
        body_from_c_id,
        true,
    );
}

#[test]
fn golden_asteroid_numbered() {
    let ast_numbers = vec![433, 7066, 136199, 2060];

    let eph_sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        asteroid_numbers: ast_numbers.clone(),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let eph_jpl = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Jpl,
        ephe_path: Some(ephe_path()),
        asteroid_numbers: ast_numbers,
        ..EphemerisConfig::default()
    })
    .expect("JPL ephemeris required for numbered asteroid golden tests (de441.eph in ephe/)");

    let cases = load("asteroid_numbered.json");
    assert_eq!(cases.len(), 72, "expected 72 numbered asteroid cases");

    check_cases(
        &cases,
        |c| {
            let flags = CalcFlags::from_bits_truncate(c.flags);
            if flags.contains(CalcFlags::JPLEPH) {
                &eph_jpl
            } else {
                &eph_sweph
            }
        },
        |id| {
            assert!(id >= 10000, "expected SE_AST_OFFSET body, got {id}");
            Body::asteroid(id - 10000).unwrap()
        },
        true,
    );
}

#[test]
fn golden_asteroid_moseph() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        asteroid_numbers: vec![433],
        ..EphemerisConfig::default()
    })
    .unwrap();

    let cases = load("asteroid_moseph.json");
    assert_eq!(cases.len(), 27, "expected 27 MOSEPH asteroid cases");

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = if c.body >= 10000 {
            Body::asteroid(c.body - 10000).unwrap()
        } else {
            body_from_c_id(c.body)
        };
        let flags = CalcFlags::from_bits_truncate(c.flags);

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "case {i} body={} jd={:.1} {}: error: {e}",
                    c.body, c.jd, c.flag_name
                ));
                continue;
            }
        };

        let label = format!("case {i} body={} jd={:.1} {}", c.body, c.jd, c.flag_name);

        if !result.flags_used.contains(CalcFlags::MOSEPH) {
            failures.push(format!(
                "{label}: flags_used missing MOSEPH: {:?}",
                result.flags_used
            ));
        }
        let retflag_expected = CalcFlags::from_bits_truncate(c.retflag);
        let retflag_mask = CalcFlags::MOSEPH | CalcFlags::SPEED;
        if result.flags_used & retflag_mask != retflag_expected & retflag_mask {
            failures.push(format!(
                "{label}: retflag mismatch: expected {:?}, got {:?}",
                retflag_expected & retflag_mask,
                result.flags_used & retflag_mask,
            ));
        }

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            // MOSEPH positions carry the stateless-vs-stateful architecture
            // tolerance (~0.4 mas, same class as CLAUDE.md <stateless_tolerance>):
            // C's pipeline reads obliquity/nutation/deflection from global caches
            // populated within the same swe_calc call; our stateless port
            // recomputes them independently. Worst observed: 3.8e-7 (Eros lon).
            let eps = if k >= 3 { 1e-7 } else { 5e-7 };
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                    c.output[k], result.data[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(200) {
            eprintln!("{f}");
        }
        panic!("{n} element failures (showing first 200)");
    }
}

#[test]
fn asteroid_error_chiron_beyond_limits() {
    let jd_before_chiron = 1967601.0;

    let eph_swiss = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let eph_moshier = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    for eph in [&eph_swiss, &eph_moshier] {
        match eph.calc(jd_before_chiron, Body::Chiron, CalcFlags::SPEED) {
            Err(swisseph::Error::BeyondEphemerisLimits { .. }) => {}
            Err(e) => panic!("expected BeyondEphemerisLimits, got: {e}"),
            Ok(_) => panic!("expected error for Chiron before CHIRON_START"),
        }
    }
}

#[test]
fn asteroid_error_pholus_beyond_limits() {
    let jd_before_pholus = 640648.0;

    let eph_swiss = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let eph_moshier = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    for eph in [&eph_swiss, &eph_moshier] {
        match eph.calc(jd_before_pholus, Body::Pholus, CalcFlags::SPEED) {
            Err(swisseph::Error::BeyondEphemerisLimits { .. }) => {}
            Err(e) => panic!("expected BeyondEphemerisLimits, got: {e}"),
            Ok(_) => panic!("expected error for Pholus before PHOLUS_START"),
        }
    }
}

#[test]
fn asteroid_error_not_in_numbers() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        asteroid_numbers: vec![433],
        ..EphemerisConfig::default()
    })
    .unwrap();

    let body = Body::asteroid(99999).unwrap();
    match eph.calc(2451545.0, body, CalcFlags::SPEED) {
        Err(swisseph::Error::EphemerisNotAvailable { .. }) => {}
        Err(e) => panic!("expected EphemerisNotAvailable, got: {e}"),
        Ok(_) => panic!("expected error for asteroid(99999) not in asteroid_numbers"),
    }
}

#[test]
fn asteroid_error_outside_file_range() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        asteroid_numbers: vec![433],
        ..EphemerisConfig::default()
    })
    .unwrap();

    let body = Body::asteroid(433).unwrap();
    let jd_before_file = 2200000.5;
    match eph.calc(jd_before_file, body, CalcFlags::SPEED) {
        Err(swisseph::Error::EphemerisNotAvailable { .. }) => {}
        Err(e) => panic!("expected EphemerisNotAvailable, got: {e}"),
        Ok(_) => panic!("expected error for asteroid 433 at jd before file range"),
    }
}

#[test]
fn asteroid_alias_identities() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let jd = 2451545.0;
    let flags = CalcFlags::SPEED;

    let pluto_named = eph.calc(jd, Body::Pluto, flags).unwrap();
    let pluto_numbered = eph
        .calc(jd, Body::asteroid(134340).unwrap(), flags)
        .unwrap();
    for k in 0..6 {
        assert_eq!(
            pluto_named.data[k].to_bits(),
            pluto_numbered.data[k].to_bits(),
            "Pluto alias [{k}]: named={:.15e}, asteroid(134340)={:.15e}",
            pluto_named.data[k],
            pluto_numbered.data[k]
        );
    }

    let ceres_named = eph.calc(jd, Body::Ceres, flags).unwrap();
    let ceres_numbered = eph.calc(jd, Body::asteroid(1).unwrap(), flags).unwrap();
    for k in 0..6 {
        assert_eq!(
            ceres_named.data[k].to_bits(),
            ceres_numbered.data[k].to_bits(),
            "Ceres alias [{k}]: named={:.15e}, asteroid(1)={:.15e}",
            ceres_named.data[k],
            ceres_numbered.data[k]
        );
    }
}
