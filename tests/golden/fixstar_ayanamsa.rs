use serde::Deserialize;
use swisseph::types::Body;
use swisseph::{CalcFlags, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct AyaCase {
    idx: i32,
    tjd: f64,
    daya: f64,
    retflag: u32,
}

#[derive(Deserialize)]
struct CalcCase {
    idx: i32,
    tjd: f64,
    lon: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    ayanamsa: Vec<AyaCase>,
    calc: Vec<CalcCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("fixstar_ayanamsa.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn make_eph(idx: i32) -> swisseph::Ephemeris {
    let mut cfg = EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some("../swisseph/ephe".into()),
        ..Default::default()
    };
    cfg.set_sidereal_mode(idx, 0.0, 0.0);
    swisseph::Ephemeris::new(cfg).unwrap()
}

#[test]
fn golden_fixstar_ayanamsa() {
    let data = load();

    assert_eq!(data.ayanamsa.len(), 12 * 4, "unexpected case count");

    let mut failures: Vec<String> = Vec::new();

    for (i, c) in data.ayanamsa.iter().enumerate() {
        let eph = make_eph(c.idx);
        let label = format!("case {i} idx={} tjd={:.1}", c.idx, c.tjd);

        match eph.get_ayanamsa_ex(c.tjd, CalcFlags::MOSEPH) {
            Ok(daya) => {
                let diff = (daya - c.daya).abs();
                if diff > 1e-8 {
                    failures.push(format!(
                        "{label}: expected {:.17e}, got {:.17e}, diff {diff:.3e}",
                        c.daya, daya
                    ));
                }
                let _ = c.retflag; // retflag for ayanamsa is epheflag; not exposed by our API
            }
            Err(e) => {
                failures.push(format!("{label}: error {e}"));
            }
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("{f}");
        }
        panic!("{} fixstar_ayanamsa failures", failures.len());
    }
}

#[test]
fn golden_fixstar_ayanamsa_nonut() {
    // Spot-check with NONUT — ensures the nutation branch is also tested.
    let data = load();
    let mut failures: Vec<String> = Vec::new();

    for c in data.ayanamsa.iter().filter(|c| c.tjd == 2451545.0) {
        let eph = make_eph(c.idx);
        // We don't have a separate no-nut golden value; just confirm it doesn't error.
        let label = format!("NONUT idx={}", c.idx);
        if let Err(e) = eph.get_ayanamsa_ex(c.tjd, CalcFlags::MOSEPH | CalcFlags::NONUT) {
            failures.push(format!("{label}: {e}"));
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("{f}");
        }
        panic!("{} fixstar_ayanamsa NONUT failures", failures.len());
    }
}

#[test]
fn golden_fixstar_ayanamsa_calc_sidereal() {
    let data = load();
    let mut failures: Vec<String> = Vec::new();

    for (i, c) in data.calc.iter().enumerate() {
        let eph = make_eph(c.idx);
        let label = format!("calc case {i} idx={} tjd={:.1}", c.idx, c.tjd);

        match eph.calc(c.tjd, Body::Sun, CalcFlags::MOSEPH | CalcFlags::SIDEREAL) {
            Ok(result) => {
                let lon = result.data[0];
                let diff = (lon - c.lon).abs();
                if diff > 1e-8 {
                    failures.push(format!(
                        "{label}: expected {:.17e}, got {:.17e}, diff {diff:.3e}",
                        c.lon, lon
                    ));
                }
            }
            Err(e) => {
                failures.push(format!("{label}: error {e}"));
            }
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("{f}");
        }
        panic!("{} fixstar_ayanamsa calc failures", failures.len());
    }
}
