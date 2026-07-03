use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource, TopoPosition};

#[derive(Deserialize)]
struct CalcCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    retflag: u32,
    output: [f64; 6],
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("asteroid.json");
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

    let cases = load();
    assert!(
        cases.len() >= 300,
        "expected 300+ cases, got {}",
        cases.len()
    );

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let is_jpl = flags.contains(CalcFlags::JPLEPH);

        let eph: &Ephemeris = if is_jpl { &eph_jpl } else { &eph_sweph };

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "case {i} {} jd={:.1} {}: error: {e}",
                    c.body_name, c.jd, c.flag_name
                ));
                continue;
            }
        };

        let label = format!("case {i} {} jd={:.1} {}", c.body_name, c.jd, c.flag_name);

        // Retflag: C returns SWIEPH for JPLEPH+asteroid cases (the asteroid itself
        // is always from a .se1 file). Only check non-ephemeris-source bits.
        if !is_jpl {
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
            // JPLEPH: asteroid from .se1, Earth/Sun from JPL DE — different Earth/Sun
            // source produces ~1e-8 (near J2000) to ~5e-7 (far epochs) position diffs.
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
