use serde::Deserialize;
use std::path::PathBuf;
use swisseph::jpl::JplFile;

#[derive(Deserialize)]
struct GoldenData {
    ss: [f64; 3],
    au: f64,
    emrat: f64,
    denum: i32,
    ncon: i32,
    ipt: Vec<i32>,
    ksize: usize,
    ncoeffs: usize,
}

#[test]
fn jpl_header_matches_c() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("ephe")
        .join("de441.eph");
    if !path.exists() {
        return;
    }

    let data: GoldenData = serde_json::from_str(
        &std::fs::read_to_string(super::golden_data_path("jpl_header.json")).unwrap(),
    )
    .unwrap();

    let jpl = JplFile::open(&path).unwrap_or_else(|e| panic!("failed to open de441.eph: {e}"));
    let h = jpl.header();

    for (i, (expected, actual)) in data.ss.iter().zip(h.ss.iter()).enumerate() {
        super::assert_f64_exact(&format!("ss[{i}]"), *expected, *actual);
    }
    super::assert_f64_exact("au", data.au, h.au);
    super::assert_f64_exact("emrat", data.emrat, h.emrat);
    assert_eq!(h.denum, data.denum, "denum");
    assert_eq!(h.ncon, data.ncon, "ncon");
    for (i, (expected, actual)) in data.ipt.iter().zip(h.ipt.iter()).enumerate() {
        assert_eq!(*actual, *expected, "ipt[{i}]");
    }
    assert_eq!(h.ksize, data.ksize, "ksize");
    assert_eq!(h.ncoeffs, data.ncoeffs, "ncoeffs");
}
