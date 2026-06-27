use serde::Deserialize;
use std::path::PathBuf;
use swisseph::sweph_file::{SwissEphFile, evaluate_body};

#[derive(Deserialize)]
struct EvalCase {
    body_id: i32,
    jd: f64,
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    #[allow(dead_code)]
    neval: usize,
}

#[derive(Deserialize)]
struct GoldenData {
    cases: Vec<EvalCase>,
}

fn ephe_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
        .join(name)
}

#[test]
fn evaluate_body_matches_c() {
    let data: GoldenData = serde_json::from_str(
        &std::fs::read_to_string(super::golden_data_path("sweph_eval.json")).unwrap(),
    )
    .unwrap();

    let planet_file = SwissEphFile::open(&ephe_path("sepl_18.se1")).unwrap();
    let moon_file = SwissEphFile::open(&ephe_path("semo_18.se1")).unwrap();

    let mut tested = 0;
    for case in &data.cases {
        let file = if case.body_id == 1 {
            &moon_file
        } else {
            &planet_file
        };
        let result = match evaluate_body(file, case.body_id, case.jd, true) {
            Ok(r) => r,
            Err(swisseph::Error::BeyondEphemerisLimits { .. }) => continue,
            Err(e) => panic!("unexpected error for body{}@{}: {e}", case.body_id, case.jd),
        };
        tested += 1;
        let label = format!("body{}@{:.1}", case.body_id, case.jd);

        super::assert_f64_exact(&format!("{label}:x"), case.x, result[0]);
        super::assert_f64_exact(&format!("{label}:y"), case.y, result[1]);
        super::assert_f64_exact(&format!("{label}:z"), case.z, result[2]);
        super::assert_f64_exact(&format!("{label}:vx"), case.vx, result[3]);
        super::assert_f64_exact(&format!("{label}:vy"), case.vy, result[4]);
        super::assert_f64_exact(&format!("{label}:vz"), case.vz, result[5]);
    }
    assert!(
        tested >= 70,
        "expected at least 70 cases, only tested {tested}"
    );
}
