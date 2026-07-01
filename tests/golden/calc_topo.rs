use serde::Deserialize;
use swisseph::{CalcFlags, Ephemeris, EphemerisConfig, TopoPosition};

#[derive(Deserialize)]
struct CalcTopoCase {
    lon: f64,
    lat: f64,
    alt: f64,
    body: i32,
    body_name: String,
    jd: f64,
    output: [f64; 6],
}

fn load() -> Vec<CalcTopoCase> {
    let path = super::golden_data_path("calc_topo.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn body_from_c_id(id: i32) -> swisseph::Body {
    use swisseph::Body;
    match id {
        0 => Body::Sun,
        1 => Body::Moon,
        2 => Body::Mercury,
        3 => Body::Venus,
        4 => Body::Mars,
        _ => panic!("unexpected body id {id}"),
    }
}

#[test]
fn golden_calc_topo() {
    let cases = load();
    assert!(!cases.is_empty(), "expected golden cases, got none");

    let mut failures = Vec::new();
    let mut current_topo = None;
    let mut eph: Option<Ephemeris> = None;

    for (i, c) in cases.iter().enumerate() {
        let topo = (c.lon, c.lat, c.alt);
        if current_topo != Some(topo) {
            current_topo = Some(topo);
            let config = EphemerisConfig {
                topographic: Some(TopoPosition {
                    longitude: c.lon,
                    latitude: c.lat,
                    altitude: c.alt,
                }),
                ..Default::default()
            };
            eph = Some(Ephemeris::new(config).unwrap());
        }

        let body = body_from_c_id(c.body);
        let flags =
            CalcFlags::MOSEPH | CalcFlags::TOPOCTR | CalcFlags::EQUATORIAL | CalcFlags::SPEED;
        let result = eph.as_ref().unwrap().calc(c.jd, body, flags).unwrap();

        let label = format!(
            "case {i} {} jd={:.1} lon={} lat={} alt={}",
            c.body_name, c.jd, c.lon, c.lat, c.alt
        );

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            let eps = if k < 3 { 1e-9 } else { 1e-7 };
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                    c.output[k], result.data[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(30) {
            eprintln!("{f}");
        }
        panic!("{n} element failures (showing first 30)");
    }
}
