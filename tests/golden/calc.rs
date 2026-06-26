use serde::Deserialize;
use swisseph::{CalcFlags, Ephemeris, EphemerisConfig};

#[derive(Deserialize)]
struct CalcCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    output: [f64; 6],
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("calc.json");
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
        5 => Body::Jupiter,
        6 => Body::Saturn,
        _ => panic!("unexpected body id {id}"),
    }
}

#[test]
fn golden_calc() {
    let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();
    let cases = load();
    assert_eq!(cases.len(), 350);

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let result = eph.calc(c.jd, body, flags).unwrap();

        let label = format!("case {i} {} jd={:.1} {}", c.body_name, c.jd, c.flag_name);

        for k in 0..6 {
            if k >= 3 && !flags.contains(CalcFlags::SPEED) {
                continue;
            }
            let diff = (c.output[k] - result.data[k]).abs();
            let eps = 1e-10;
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
