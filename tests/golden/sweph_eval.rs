use serde::Deserialize;
use std::path::PathBuf;
use swisseph::sweph_file::{SwissEphFile, evaluate_body, find_file_for_jd, open_ephemeris_files};

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
    neval: usize,
}

#[derive(Deserialize)]
struct GoldenData {
    cases: Vec<EvalCase>,
}

fn ephe_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
}

#[test]
fn evaluate_body_matches_c() {
    let data: GoldenData = serde_json::from_str(
        &std::fs::read_to_string(super::golden_data_path("sweph_eval.json")).unwrap(),
    )
    .unwrap();

    let dir = ephe_dir();
    let planet_files = open_ephemeris_files(&dir, "sepl").unwrap();
    let moon_files = open_ephemeris_files(&dir, "semo").unwrap();

    let mut tested = 0;
    for case in &data.cases {
        let files: &[SwissEphFile] = if case.body_id == 1 {
            &moon_files
        } else {
            &planet_files
        };
        let file = match find_file_for_jd(files, case.body_id, case.jd) {
            Some(f) => f,
            None => continue,
        };
        let (result, neval) = match evaluate_body(file, case.body_id, case.jd, true) {
            Ok(r) => r,
            Err(swisseph::Error::BeyondEphemerisLimits { .. }) => continue,
            Err(e) => panic!("unexpected error for body{}@{}: {e}", case.body_id, case.jd),
        };
        tested += 1;
        let label = format!("body{}@{:.1}", case.body_id, case.jd);

        assert_eq!(neval, case.neval, "{label}:neval");
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
