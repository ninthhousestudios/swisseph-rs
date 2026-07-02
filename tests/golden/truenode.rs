use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct NodeCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    eph_name: String,
    /// C `sid_mode` passed to `swe_set_sid_mode` (0 = tropical, no SEFLG_SIDEREAL).
    #[serde(default)]
    sid_mode: i32,
    retflag: i32,
    output: [f64; 6],
}

fn load() -> Vec<NodeCase> {
    let path = super::golden_data_path("truenode.json");
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
        11 => Body::TrueNode,
        13 => Body::OscuApogee,
        _ => panic!("unexpected body id {id}"),
    }
}

/// Build an `Ephemeris` for a `(backend, sid_mode)` pair. `sid_mode == 0` is
/// tropical; a non-zero value is the raw C `swe_set_sid_mode` argument (ayanamsa
/// index OR'd with the SE_SIDBIT_* projection bits).
fn make_eph(eph_name: &str, sid_mode: i32) -> Ephemeris {
    let mut cfg = if eph_name == "SWIEPH" {
        EphemerisConfig {
            ephemeris_source: EphemerisSource::Swiss,
            ephe_path: Some(ephe_path()),
            ..EphemerisConfig::default()
        }
    } else {
        EphemerisConfig::default()
    };
    if sid_mode != 0 {
        cfg.set_sidereal_mode(sid_mode, 0.0, 0.0);
    }
    Ephemeris::new(cfg).expect("Ephemeris::new")
}

/// `SE_TRUE_NODE` / `SE_OSCU_APOG` through `Ephemeris::calc` (`lunar_osc_elem` +
/// `swi_plan_for_osc_elem`): 252 cases.
///
/// - 168 tropical: 2 bodies × {MOSEPH, SWIEPH} × 6 flag combos {SPEED,
///   SPEED|EQUATORIAL, SPEED|XYZ, SPEED|NONUT, SPEED|J2000, no_speed} × 7 epochs.
/// - 84 sidereal (SEFLG_SIDEREAL|SEFLG_SPEED): 2 bodies × {MOSEPH, SWIEPH} × 3
///   sid_modes {Lahiri traditional, Lahiri|ECL_T0, Lahiri|SSY_PLANE} × 7 epochs.
///   The ECL_T0 / SSY_PLANE "rigorous" projections read the J2000 equatorial
///   vector `lunar_osc_elem` now threads through as `x2000` (swisseph-rs/84 review
///   follow-up); before the fix they silently fell back to traditional ayanamsa
///   subtraction because `calc_inner` returned an all-zero `x2000`.
///
/// Positions eps 1e-9, speeds eps 1e-7 — except MOSEPH speeds, relaxed to 5e-6:
/// the osculating elements are built from a central difference over the wide
/// 0.1-day Moshier interval whose off-center samples read C's global
/// obliquity/nutation cache, a stateless-vs-stateful artifact (~0.013"/day) that
/// the sidereal projection carries through from the tropical speed. See
/// CLAUDE.md <stateless_tolerance> §3.
#[test]
fn golden_truenode() {
    let cases = load();
    assert!(
        cases.len() >= 252,
        "expected 252+ cases, got {}",
        cases.len()
    );

    let mut ephemerides: HashMap<(String, i32), Ephemeris> = HashMap::new();

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = body_from_c_id(c.body);
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemerides
            .entry((c.eph_name.clone(), c.sid_mode))
            .or_insert_with(|| make_eph(&c.eph_name, c.sid_mode));

        let sid_label = if c.sid_mode == 0 {
            String::new()
        } else {
            format!(" sid={}", c.sid_mode)
        };
        let label = format!(
            "case {i} {} tjd={:.1} {} {}{sid_label}",
            c.body_name, c.jd, c.eph_name, c.flag_name
        );

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        if result.flags_used.bits() as i32 != c.retflag {
            failures.push(format!(
                "{label}: retflag mismatch (expected {:#x}, got {:#x})",
                c.retflag,
                result.flags_used.bits()
            ));
        }

        let moseph = c.eph_name == "MOSEPH";
        for k in 0..6 {
            let eps = if k >= 3 {
                if moseph { 5e-6 } else { 1e-7 }
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
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 40)");
    }
}
