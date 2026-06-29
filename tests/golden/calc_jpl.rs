use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct CalcCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    output: [f64; 6],
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("calc_jpl.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
}

fn jpl_file_exists() -> bool {
    ephe_path().join("de441.eph").exists()
}

fn body_from_c_id(id: i32) -> swisseph::Body {
    use swisseph::Body;
    match id {
        0 => Body::Sun,
        1 => Body::Moon,
        2 => Body::Mercury,
        3 => Body::Venus,
        4 => Body::Mars,
        5 => Body::Jupiter,
        6 => Body::Saturn,
        7 => Body::Uranus,
        8 => Body::Neptune,
        9 => Body::Pluto,
        _ => panic!("unexpected body id {id}"),
    }
}

#[test]
fn golden_calc_jpl() {
    if !jpl_file_exists() {
        eprintln!("SKIP: ephe/de441.eph not found");
        return;
    }

    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Jpl,
        ephe_path: Some(ephe_path()),
        jpl_filename: Some("de441.eph".to_string()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let cases = load();
    assert!(
        cases.len() >= 700,
        "expected 700+ cases, got {}",
        cases.len()
    );

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
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
        let is_speed3 = flags.contains(CalcFlags::SPEED3);

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            // Speed components: stateless deflection geometry differs from C's
            // cached sun position by <1e-7 deg/day. SPEED3 has no file
            // boundaries in JPL (single continuous file), so 1e-7 is tight.
            // Branches kept distinct to document each tolerance's rationale even
            // where the values currently coincide.
            #[allow(clippy::if_same_then_else)]
            let eps = if k >= 3 {
                1e-7
            } else if is_speed3 {
                1e-8
            } else {
                1e-8
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
        for f in failures.iter().take(3000) {
            eprintln!("{f}");
        }
        panic!("{n} element failures");
    }
}
