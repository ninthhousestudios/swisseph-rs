use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize)]
struct BackendCase {
    jd: f64,
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}

type GoldenData = HashMap<String, Vec<BackendCase>>;

fn load() -> GoldenData {
    let path = super::golden_data_path("moshier_backend.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

#[test]
fn golden_moshier_backend() {
    use swisseph::constants::J2000;
    use swisseph::flags::CalcFlags;
    use swisseph::moshier::backend::compute;
    use swisseph::obliquity::obliquity;
    use swisseph::types::{AstroModels, Body};

    let eps_j2000 = obliquity(J2000, CalcFlags::empty(), &AstroModels::default());

    let body_map: &[(&str, Body)] = &[
        ("sun", Body::Sun),
        ("moon", Body::Moon),
        ("mercury", Body::Mercury),
        ("venus", Body::Venus),
        ("mars", Body::Mars),
        ("jupiter", Body::Jupiter),
        ("saturn", Body::Saturn),
        ("uranus", Body::Uranus),
        ("neptune", Body::Neptune),
        ("pluto", Body::Pluto),
    ];

    let data = load();
    let mut total = 0;
    for (name, body) in body_map {
        let cases = data
            .get(*name)
            .unwrap_or_else(|| panic!("missing body {name}"));
        for (i, c) in cases.iter().enumerate() {
            let result =
                compute(c.jd, *body, &eps_j2000).unwrap_or_else(|e| panic!("{name} case {i}: {e}"));
            let label = format!("{name} case {i} jd={:.1}", c.jd);
            super::assert_f64_eps(&format!("{label} x"), c.x, result[0], 1e-15);
            super::assert_f64_eps(&format!("{label} y"), c.y, result[1], 1e-15);
            super::assert_f64_eps(&format!("{label} z"), c.z, result[2], 1e-15);
            super::assert_f64_eps(&format!("{label} vx"), c.vx, result[3], 1e-10);
            super::assert_f64_eps(&format!("{label} vy"), c.vy, result[4], 1e-10);
            super::assert_f64_eps(&format!("{label} vz"), c.vz, result[5], 1e-10);
            total += 1;
        }
    }
    assert_eq!(total, 110, "expected 10 bodies × 11 epochs = 110 cases");
}

#[test]
fn earth_geocentric_is_zero() {
    use swisseph::constants::J2000;
    use swisseph::flags::CalcFlags;
    use swisseph::moshier::backend::compute;
    use swisseph::obliquity::obliquity;
    use swisseph::types::{AstroModels, Body};

    let eps = obliquity(J2000, CalcFlags::empty(), &AstroModels::default());
    let result = compute(J2000, Body::Earth, &eps).unwrap();
    assert_eq!(result, [0.0; 6]);
}
