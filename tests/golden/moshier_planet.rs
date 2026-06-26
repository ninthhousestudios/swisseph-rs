use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize)]
struct MoshierCase {
    jd: f64,
    lon: f64,
    lat: f64,
    dist: f64,
}

type GoldenData = HashMap<String, Vec<MoshierCase>>;

fn load() -> GoldenData {
    let path = super::golden_data_path("moshier_planet.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

#[test]
fn golden_moshplan2() {
    use swisseph::moshier::{PlantTbl, planets::moshplan2, tables};

    let planet_tables: &[(&str, &PlantTbl)] = &[
        ("mercury", &tables::MER404),
        ("venus", &tables::VEN404),
        ("earth", &tables::EAR404),
        ("mars", &tables::MAR404),
        ("jupiter", &tables::JUP404),
        ("saturn", &tables::SAT404),
        ("uranus", &tables::URA404),
        ("neptune", &tables::NEP404),
        ("pluto", &tables::PLU404),
    ];

    let data = load();
    let mut total = 0;
    for (name, tbl) in planet_tables {
        let cases = data
            .get(*name)
            .unwrap_or_else(|| panic!("missing planet {name}"));
        for (i, c) in cases.iter().enumerate() {
            let result = moshplan2(c.jd, tbl);
            let label = format!("{name} case {i} jd={:.1}", c.jd);
            super::assert_f64_exact(&format!("{label} lon"), c.lon, result[0]);
            super::assert_f64_exact(&format!("{label} lat"), c.lat, result[1]);
            super::assert_f64_exact(&format!("{label} dist"), c.dist, result[2]);
            total += 1;
        }
    }
    assert_eq!(total, 81, "expected 9 planets × 9 epochs = 81 cases");
}
