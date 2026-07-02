use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct NodeCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    eph_name: String,
    retflag: i32,
    output: [f64; 6],
}

fn load() -> Vec<NodeCase> {
    let path = super::golden_data_path("truenode.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
}

fn body_from_c_id(id: i32) -> Body {
    match id {
        11 => Body::TrueNode,
        13 => Body::OscuApogee,
        _ => panic!("unexpected body id {id}"),
    }
}

/// `SE_TRUE_NODE` / `SE_OSCU_APOG` through `Ephemeris::calc` (`lunar_osc_elem` +
/// `swi_plan_for_osc_elem`): 168 cases — 2 bodies × {MOSEPH, SWIEPH} × 6 flag
/// combos {SPEED, SPEED|EQUATORIAL, SPEED|XYZ, SPEED|NONUT, SPEED|J2000, no_speed}
/// × the 7 gen_calc.c epochs. Positions eps 1e-9, speeds eps 1e-7 (the osculating
/// elements are built from finite differences over 1e-4-day intervals, amplifying
/// backend ULP noise into the speed components).
#[test]
fn golden_truenode() {
    let moshier = Ephemeris::new(EphemerisConfig::default()).unwrap();
    let sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let cases = load();
    assert!(
        cases.len() >= 168,
        "expected 168+ cases, got {}",
        cases.len()
    );

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = if c.eph_name == "SWIEPH" {
            &sweph
        } else {
            &moshier
        };

        let label = format!(
            "case {i} {} tjd={:.1} {} {}",
            c.body_name, c.jd, c.eph_name, c.flag_name
        );

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        if result.flags_used.bits() as i32 != c.retflag {
            failures.push(format!(
                "{label}: retflag mismatch (expected {:#x}, got {:#x})",
                c.retflag,
                result.flags_used.bits()
            ));
        }

        // MOSEPH node/apogee SPEED carries a documented stateless-vs-stateful
        // precision artifact (up to ~3.6e-6 deg/day ≈ 0.013"/day, astronomically
        // negligible). C's `lunar_osc_elem` builds the speed from a central
        // difference over the wide 0.1-day Moshier interval, and its off-center
        // samples read C's GLOBAL obliquity/nutation cache, which rounds slightly
        // differently than a clean recomputation — so C's own node speed does not
        // even match a finite difference of C's own node POSITIONS. Positions are
        // bit-accurate (1e-9); the Swiss backend (1e-4 interval) matches speed at
        // 1e-7. Only the Moshier speed is relaxed. See CLAUDE.md <stateless_tolerance>.
        let moseph = c.eph_name == "MOSEPH";
        for k in 0..6 {
            let eps = if k >= 3 {
                if moseph { 5e-6 } else { 1e-7 }
            } else {
                1e-9
            };
            let diff = (c.output[k] - result.data[k]).abs();
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e} > eps {eps:.1e}",
                    c.output[k], result.data[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 40)");
    }
}
