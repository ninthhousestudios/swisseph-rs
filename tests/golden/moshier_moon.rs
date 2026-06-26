use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize)]
struct MoonCase {
    jd: f64,
    lon: f64,
    lat: f64,
    dist: f64,
}

type GoldenData = HashMap<String, Vec<MoonCase>>;

fn load() -> GoldenData {
    let path = super::golden_data_path("moshier_moon.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

#[test]
fn golden_moshmoon2() {
    use swisseph::moshier::moon::moshmoon2;

    let data = load();
    let cases = data
        .get("moon")
        .unwrap_or_else(|| panic!("missing 'moon' key in golden data"));
    for (i, c) in cases.iter().enumerate() {
        let result = moshmoon2(c.jd);
        let label = format!("moon case {i} jd={:.1}", c.jd);
        super::assert_f64_exact(&format!("{label} lon"), c.lon, result[0]);
        super::assert_f64_exact(&format!("{label} lat"), c.lat, result[1]);
        super::assert_f64_exact(&format!("{label} dist"), c.dist, result[2]);
    }
    assert!(cases.len() >= 10, "expected >= 10 test epochs");
}
