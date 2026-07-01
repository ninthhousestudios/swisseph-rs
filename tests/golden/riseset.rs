use serde::Deserialize;
use swisseph::{Body, CalcFlags, Ephemeris, Error, RiseSetFlags};

#[derive(Deserialize)]
struct FullCase {
    geopos: [f64; 3],
    #[allow(dead_code)]
    geopos_name: String,
    body: String,
    tjd_ut: f64,
    rsmi: String,
    retval: i32,
    tret0: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    full: Vec<FullCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("riseset.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn body_of(s: &str) -> Body {
    match s {
        "Sun" => Body::Sun,
        "Moon" => Body::Moon,
        other => panic!("Unknown body: {other}"),
    }
}

fn rsmi_of(s: &str) -> RiseSetFlags {
    match s {
        "RISE" => RiseSetFlags::RISE,
        "SET" => RiseSetFlags::SET,
        "MTRANSIT" => RiseSetFlags::MTRANSIT,
        other => panic!("Unknown rsmi: {other}"),
    }
}

#[test]
fn full() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.full.iter().enumerate() {
        let body = body_of(&c.body);
        let rsmi = rsmi_of(&c.rsmi) | RiseSetFlags::FORCE_SLOW;
        let label = format!(
            "full[{i}][geopos={:?},body={},tjd_ut={},rsmi={}]",
            c.geopos, c.body, c.tjd_ut, c.rsmi
        );
        let actual = ephe.rise_trans_true_hor(
            c.tjd_ut,
            body,
            None,
            CalcFlags::MOSEPH,
            rsmi,
            c.geopos,
            1013.25,
            15.0,
            0.0,
        );
        if c.retval == -2 {
            match actual {
                Err(Error::CircumpolarBody) => {}
                other => panic!("{label}: expected CircumpolarBody, got {other:?}"),
            }
        } else {
            let result = actual.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
            super::assert_f64_eps(&format!("{label}.time"), c.tret0, result.time, 1e-6);
        }
    }
}
