use serde::Deserialize;
use std::collections::HashMap;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig};

#[derive(Deserialize)]
struct MeanElementCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    /// C `sid_mode` passed to `swe_set_sid_mode` (0 = tropical, no SEFLG_SIDEREAL).
    #[serde(default)]
    sid_mode: i32,
    output: [f64; 6],
}

fn load() -> Vec<MeanElementCase> {
    let path = super::golden_data_path("mean_elements.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn body_from_c_id(id: i32) -> Body {
    match id {
        10 => Body::MeanNode,
        12 => Body::MeanApogee,
        -1 => Body::EclipticNutation,
        _ => panic!("unexpected body id {id}"),
    }
}

/// Build a Moshier `Ephemeris` for a given `sid_mode` (0 = tropical; non-zero is
/// the raw C `swe_set_sid_mode` argument, ayanamsa index OR'd with SE_SIDBIT_*).
fn make_eph(sid_mode: i32) -> Ephemeris {
    let mut cfg = EphemerisConfig::default();
    if sid_mode != 0 {
        cfg.set_sidereal_mode(sid_mode, 0.0, 0.0);
    }
    Ephemeris::new(cfg).expect("Ephemeris::new")
}

/// Mean node / mean apogee / ecliptic-nutation via `Ephemeris::calc`: 231 cases.
///
/// - 165 tropical: 3 bodies × 11 epochs × 5 flag combos.
/// - 66 sidereal (SEFLG_SIDEREAL|SEFLG_SPEED): 2 lunar bodies × 11 epochs × 3
///   sid_modes {Lahiri traditional, Lahiri|ECL_T0, Lahiri|SSY_PLANE}. Regression
///   guard for the swisseph-rs/84 review follow-up: `mean_element_pipeline` now
///   threads the J2000 equatorial `x2000` so the ECL_T0/SSY_PLANE rigorous
///   projections no longer silently degrade to traditional ayanamsa subtraction.
///   The x2000 precession is guarded `jd != J2000` (C's app_pos_etc_mean skips it
///   at exactly J2000; otherwise `precess_speed` adds a spurious ~3.8e-5 deg/day
///   rate term to the sidereal speed).
#[test]
fn golden_mean_elements() {
    let cases = load();
    assert_eq!(cases.len(), 231);

    let mut ephemerides: HashMap<i32, Ephemeris> = HashMap::new();

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemerides
            .entry(c.sid_mode)
            .or_insert_with(|| make_eph(c.sid_mode));
        let result = eph.calc(c.jd, body, flags).unwrap();

        let sid_label = if c.sid_mode == 0 {
            String::new()
        } else {
            format!(" sid={}", c.sid_mode)
        };
        let label = format!(
            "case {i} {} jd={:.1} {}{sid_label}",
            c.body_name, c.jd, c.flag_name
        );

        // Tropical cases are bitwise-tight (1e-10). Sidereal cases route through the
        // extra x2000 precession + projection, so allow 1e-9 positions / 1e-7 speed.
        let sidereal = c.sid_mode != 0;
        for k in 0..6 {
            if k >= 3 && !flags.contains(CalcFlags::SPEED) {
                continue;
            }
            let eps = if !sidereal {
                1e-10
            } else if k >= 3 {
                1e-7
            } else {
                1e-9
            };
            let diff = (c.output[k] - result.data[k]).abs();
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
        for f in failures.iter().take(30) {
            eprintln!("{f}");
        }
        panic!("{n} element failures (showing first 30)");
    }
}
