use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct PhenoCase {
    tjd_et: f64,
    ipl: i32,
    body_name: String,
    iflag: u32,
    flag_name: String,
    retflag: i32,
    attr: [f64; 6],
}

fn load() -> Vec<PhenoCase> {
    let path = super::golden_data_path("pheno.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
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

/// `swe_pheno` (`Ephemeris::pheno`): phase angle, illuminated fraction, elongation, apparent
/// diameter, apparent magnitude, and Moon horizontal parallax over Sun..Pluto × 4 epochs × the
/// {MOSEPH, MOSEPH|TRUEPOS, SWIEPH} flag battery (120 cases). Exercises magnitude branches 5a-5j;
/// the Bowell §5k asteroid branch is unreachable via golden data (no backend computes asteroid
/// positions yet -- see gen_pheno.c).
#[test]
fn golden_pheno() {
    let moshier = Ephemeris::new(EphemerisConfig::default()).unwrap();
    let sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let cases = load();
    assert!(
        cases.len() >= 120,
        "expected 120+ cases, got {}",
        cases.len()
    );

    let attr_labels = [
        "phase_angle",
        "phase",
        "elongation",
        "apparent_diameter",
        "apparent_magnitude",
        "horizontal_parallax",
    ];

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.ipl);
        let flags = CalcFlags::from_bits_truncate(c.iflag);
        let eph = if flags.contains(CalcFlags::SWIEPH) {
            &sweph
        } else {
            &moshier
        };

        let label = format!(
            "case {i} {} tjd={:.1} {}",
            c.body_name, c.tjd_et, c.flag_name
        );
        let (result, retflag) = match eph.pheno(c.tjd_et, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        if retflag.bits() as i32 != c.retflag {
            failures.push(format!(
                "{label}: retflag mismatch (expected {:#x}, got {:#x})",
                c.retflag,
                retflag.bits()
            ));
        }

        let got = [
            result.phase_angle,
            result.phase,
            result.elongation,
            result.apparent_diameter,
            result.apparent_magnitude,
            result.horizontal_parallax,
        ];
        for k in 0..6 {
            // Moon apparent magnitude (attr[4]) compounds the geocentric AND heliocentric lunar
            // distances (lbr[2]*lbr2[2]); the ~1e-9 residual on each pushes it just past 1e-9, so
            // relax that one field to 1e-8 (per the task's escape note). Everything else is 1e-9.
            let eps = if k == 4 && body == Body::Moon {
                1e-8
            } else {
                1e-9
            };
            let diff = (c.attr[k] - got[k]).abs();
            if diff > eps {
                failures.push(format!(
                    "{label} [{}]: expected {:.15e}, got {:.15e}, diff {diff:.3e} > eps {eps:.1e}",
                    attr_labels[k], c.attr[k], got[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(30) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 30)");
    }
}
