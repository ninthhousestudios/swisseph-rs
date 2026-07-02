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
        _ => panic!("unexpected body id {id}"),
    }
}

/// Dedicated `SEFLG_HELCTR` golden coverage for the calc pipeline (swisseph-rs/94).
///
/// HELCTR was ported for phenomena (swisseph-rs/83) but until now only verified transitively
/// through the pheno battery, which never sets `SEFLG_SPEED` on its heliocentric calls and never
/// exercises the JPL Moon. This asserts `swe_calc(.., SEFLG_HELCTR, ..)` against C directly for
/// Sun..Pluto + Moon across the polar and XYZ frames, with/without J2000/EQUATORIAL, with/without
/// SPEED, over the Moshier / Swiss / JPL backends (1200 cases; JPL rows skipped if de441.eph is
/// absent).
///
/// Heliocentric Sun is the origin (Sun relative to itself) → all-zero output; C returns zeros and
/// so does the Rust dispatch (context.rs short-circuit). Positions eps 1e-9; speed eps 1e-7 (the
/// stateless deflection geometry is off, but HELCTR forces NOABERR|NOGDEFL so the geometric speed
/// is tight — 1e-7 documents the shared calc tolerance).
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
        cases.len() >= 1200,
        "expected 1200+ cases, got {}",
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
