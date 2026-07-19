use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct CalcCase {
    body: i32,
    #[allow(dead_code)]
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    retflag: u32,
    output: [f64; 6],
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("fictitious.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
}

#[test]
fn golden_fictitious() {
    let eph_moshier = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let eph_sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let cases = load();
    assert!(
        cases.len() >= 240,
        "expected 240+ cases, got {}",
        cases.len()
    );

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = Body::fictitious(c.body).unwrap_or_else(|e| {
            panic!("case {i} body={}: {e}", c.body);
        });
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = if flags.contains(CalcFlags::SWIEPH) {
            &eph_sweph
        } else {
            &eph_moshier
        };

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "case {i} body={} ({}) jd={:.1} {}: error: {e}",
                    c.body, c.body_name, c.jd, c.flag_name
                ));
                continue;
            }
        };

        let label = format!(
            "case {i} body={} ({}) jd={:.1} {}",
            c.body, c.body_name, c.jd, c.flag_name
        );

        let retflag_expected = CalcFlags::from_bits_truncate(c.retflag);
        let retflag_mask = CalcFlags::SWIEPH | CalcFlags::MOSEPH | CalcFlags::SPEED;
        if result.flags_used & retflag_mask != retflag_expected & retflag_mask {
            failures.push(format!(
                "{label}: retflag mismatch: expected {:?}, got {:?}",
                retflag_expected & retflag_mask,
                result.flags_used & retflag_mask,
            ));
        }

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            let is_swieph = flags.contains(CalcFlags::SWIEPH);
            let eps = if k >= 3 && is_swieph {
                5e-7
            } else if k >= 3 {
                1e-7
            } else if is_swieph {
                5e-9
            } else {
                1e-9
            };
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e} (eps {eps:.0e})",
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
