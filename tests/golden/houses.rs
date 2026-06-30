use serde::Deserialize;
use swisseph::houses::houses_armc;
use swisseph::types::HouseSystem;

#[derive(Deserialize)]
struct AnglesSpecialCase {
    armc: f64,
    geolat: f64,
    eps: f64,
    ascmc: [f64; 8],
    ascmc_speed: [f64; 8],
}

#[derive(Deserialize)]
struct EqualFamilyCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct GoldenData {
    angles_special: Vec<AnglesSpecialCase>,
    equal_family: Vec<EqualFamilyCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("houses.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn parse_hsys(s: &str) -> HouseSystem {
    HouseSystem::try_from(s.as_bytes()[0])
        .unwrap_or_else(|e| panic!("Unknown house system {s}: {e}"))
}

#[test]
fn angles_special() {
    let data = load();
    assert_eq!(
        data.angles_special.len(),
        30,
        "expected 30 golden cases (6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.angles_special.iter().enumerate() {
        // The special points are system-independent; 'P' was used by the C generator,
        // but only Equal-family systems are ported so far, so use Equal here.
        let result = houses_armc(c.armc, c.geolat, c.eps, HouseSystem::Equal, None)
            .unwrap_or_else(|e| panic!("case {i}: houses_armc failed: {e}"));

        let actual = result.ascmc.as_array();
        let actual_speed = result.ascmc_speeds.as_array();
        let label_base = format!(
            "case {i} (armc={:.6} geolat={:.6} eps={:.6})",
            c.armc, c.geolat, c.eps
        );
        for j in 0..8 {
            super::assert_f64_exact(&format!("{label_base} ascmc[{j}]"), c.ascmc[j], actual[j]);
            super::assert_f64_exact(
                &format!("{label_base} ascmc_speed[{j}]"),
                c.ascmc_speed[j],
                actual_speed[j],
            );
        }
    }
}

#[test]
fn equal_family() {
    let data = load();
    assert_eq!(
        data.equal_family.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.equal_family.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_exact(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
            );
            super::assert_f64_exact(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
            );
        }
    }
}
