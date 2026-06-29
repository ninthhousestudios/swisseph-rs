use serde::Deserialize;
use std::path::PathBuf;
use swisseph::jpl::{JplFile, jpl_pleph};

#[derive(Deserialize)]
struct PlephCase {
    ntarg: i32,
    ncent: i32,
    jd: f64,
    rrd: [f64; 6],
}

#[derive(Deserialize)]
struct GoldenData {
    cases: Vec<PlephCase>,
}

fn jpl_file() -> Option<JplFile> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("ephe")
        .join("de441.eph");
    JplFile::open(&path).ok()
}

#[test]
fn jpl_pleph_matches_c() {
    let file = match jpl_file() {
        Some(f) => f,
        None => {
            eprintln!("SKIP: ephe/de441.eph not found");
            return;
        }
    };

    let data: GoldenData = serde_json::from_str(
        &std::fs::read_to_string(super::golden_data_path("jpl_pleph.json")).unwrap(),
    )
    .unwrap();

    let mut tested = 0;
    for case in &data.cases {
        let result = jpl_pleph(&file, case.jd, case.ntarg, case.ncent, true).unwrap_or_else(|e| {
            panic!(
                "jpl_pleph failed ntarg={} ncent={} jd={}: {e}",
                case.ntarg, case.ncent, case.jd
            )
        });

        let label = format!(
            "ntarg={} ncent={} jd={:.1}",
            case.ntarg, case.ncent, case.jd
        );
        let eps = 1e-9;
        for (k, &r) in result.iter().enumerate() {
            super::assert_f64_eps(&format!("{label}:rrd[{k}]"), case.rrd[k], r, eps);
        }
        tested += 1;
    }

    assert!(
        tested >= 80,
        "expected at least 80 cases, only tested {tested}"
    );
}
