use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct CalcCase {
    backend: String,
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    output: [f64; 6],
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("calc_helctr.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Swiss (.se1) files live in the sibling C repo; the JPL DE file lives in this crate's ephe/.
fn sweph_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
}

fn jpl_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
}

fn body_from_c_id(id: i32) -> Body {
    match id {
        0 => Body::Sun,
        1 => Body::Moon,
        2 => Body::Mercury,
        3 => Body::Venus,
        4 => Body::Mars,
        5 => Body::Jupiter,
        6 => Body::Saturn,
        7 => Body::Uranus,
        8 => Body::Neptune,
        9 => Body::Pluto,
        14 => Body::Earth,
        _ => panic!("unexpected body id {id}"),
    }
}

/// Dedicated `SEFLG_HELCTR`/`SEFLG_BARYCTR` golden coverage for the calc pipeline
/// (swisseph-rs/94, Earth added swisseph-rs/96).
///
/// Asserts `swe_calc(.., SEFLG_HELCTR, ..)` against C directly for Sun..Pluto + Moon + Earth
/// across the polar and XYZ frames, with/without J2000/EQUATORIAL, with/without SPEED, over the
/// Moshier / Swiss / JPL backends. BARYCTR Earth cases for Swiss/JPL (Moshier rejects BARYCTR).
/// 1760 cases; JPL rows skipped if de441.eph is absent.
///
/// Heliocentric Sun is the origin (Sun relative to itself) → all-zero output. Positions eps 1e-9;
/// speed eps 1e-7.
#[test]
fn golden_calc_helctr() {
    let moshier = Ephemeris::new(EphemerisConfig::default()).unwrap();
    let sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(sweph_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let jpl_file = jpl_path().join("de441.eph");
    let jpl = if jpl_file.exists() {
        Some(
            Ephemeris::new(EphemerisConfig {
                ephemeris_source: EphemerisSource::Jpl,
                ephe_path: Some(jpl_path()),
                jpl_filename: Some("de441.eph".to_string()),
                ..EphemerisConfig::default()
            })
            .unwrap(),
        )
    } else {
        eprintln!("SKIP jpl rows: ephe/de441.eph not found");
        None
    };

    let cases = load();
    assert!(
        cases.len() >= 1350,
        "expected 1350+ cases, got {}",
        cases.len()
    );

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let eph = match c.backend.as_str() {
            "moshier" => &moshier,
            "sweph" => &sweph,
            "jpl" => match &jpl {
                Some(e) => e,
                None => continue,
            },
            other => panic!("unexpected backend {other}"),
        };

        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let label = format!(
            "case {i} {} {} jd={:.1} {}",
            c.backend, c.body_name, c.jd, c.flag_name
        );

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            let eps = if k >= 3 { 1e-7 } else { 1e-9 };
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e} > eps {eps:.1e}",
                    c.output[k], result.data[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} element failures (showing first 40)");
    }
}
