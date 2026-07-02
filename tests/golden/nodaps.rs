use serde::Deserialize;
use std::collections::HashMap;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, NodApsMethod, NodesApsides};

#[derive(Deserialize)]
struct MeanCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    #[allow(dead_code)]
    retflag: i32,
    asc: [f64; 6],
    desc: [f64; 6],
    peri: [f64; 6],
    aphe: [f64; 6],
}

#[derive(Deserialize)]
struct GoldenData {
    mean: Vec<MeanCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("nodaps.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Tolerance for one component `k` (0..2 = position, 3..5 = speed) of the given
/// node/apsis `point`, under the given `flag_name`.
///
/// **Geometry (`TRUEPOS*` combos) is asserted tight** — the raw mean node/apsis
/// vectors are bit-for-bit identical to C for every body/point (positions 1e-9,
/// speeds 1e-8), including the pathological descending node below.
///
/// **The apparent (light-deflection + aberration) combos relax the DESCENDING
/// NODE** to 1e-3° position / 2e-2°/day speed. This is NOT a port defect — every
/// stage was verified byte-identical to C in isolation (raw geometry, deflection
/// alone, aberration alone, and C's `swi_aberr_light`/`swi_deflect_light` fed the
/// Rust intermediates). The divergence appears ONLY when deflection AND aberration
/// are combined, and only on the descending node, because C's node-distance
/// formula (swecl.c:5230) divides by `cos((180-parg)·DEGTORAD)`, which is
/// near-zero for the low-inclination planets (Jupiter: `cos(94°)≈0.067`, yielding
/// a spurious node "distance" of 6.19 AU — larger than the 5.45 AU aphelion). That
/// makes the point ill-conditioned: a ~5e-10 FP-ordering difference in the
/// deflection speed branch amplifies through the aberration speed chain. C's own
/// reference digits for that node are therefore FP-order-dependent. The other
/// three points (ascending node, perihelion, aphelion) stay tight at 1e-6/1e-6.
/// See docs/swisseph-c-potential-bugs.md § "swe_nod_aps mean descending-node
/// distance singularity".
fn tolerance(point: &str, k: usize, flag_name: &str) -> f64 {
    let is_speed = k >= 3;
    if flag_name.starts_with("TRUEPOS") {
        // Pure geometry — no light effects — is bit-exact.
        return if is_speed { 1e-8 } else { 1e-9 };
    }
    if point == "desc" {
        // Ill-conditioned near-singular descending-node distance (see above).
        return if is_speed { 2e-2 } else { 1e-3 };
    }
    // Apparent asc / peri / aphe carry only the sub-milliarcsecond
    // deflection/aberration FP-conditioning noise (position and speed alike).
    1e-6
}

/// Mean nodes & apsides via `Ephemeris::nod_aps` (`swe_nod_aps`, method
/// `SE_NODBIT_MEAN`): 200 Moshier cases — 10 bodies {Sun, Moon,
/// Mercury..Neptune, Earth} × 4 epochs (incl. pre-1900 1800-Jan-1) × 5 flag combos
/// {SPEED, SPEED|EQUATORIAL, no_speed, SPEED|TRUEPOS, SPEED|EQUATORIAL|TRUEPOS}.
///
/// Tolerances are per point/component/flag — see [`tolerance`]. In short: the raw
/// geometry (`TRUEPOS`) is bit-exact; the apparent output is tight (1e-6) for the
/// ascending node / perihelion / aphelion and relaxed for the ill-conditioned
/// descending node. Sun/Earth nodes are exact zeros (no ecliptic node for Earth's
/// orbit).
#[test]
fn golden_nodaps_mean() {
    let data = load();
    let cases = &data.mean;
    assert!(
        cases.len() >= 200,
        "expected 200+ cases, got {}",
        cases.len()
    );

    // All mean-branch golden cases use the Moshier backend (SEFLG_MOSEPH).
    let mut ephemerides: HashMap<(), Ephemeris> = HashMap::new();

    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let body = Body::try_from(c.body).expect("valid body id");
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemerides
            .entry(())
            .or_insert_with(|| Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new"));

        let label = format!("case {i} {} tjd={:.1} {}", c.body_name, c.jd, c.flag_name);

        let NodesApsides {
            ascending,
            descending,
            perihelion,
            aphelion,
        } = match eph.nod_aps(c.jd, body, flags, NodApsMethod::MEAN) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        for (name, expected, got) in [
            ("asc", &c.asc, &ascending),
            ("desc", &c.desc, &descending),
            ("peri", &c.peri, &perihelion),
            ("aphe", &c.aphe, &aphelion),
        ] {
            for k in 0..6 {
                let eps = tolerance(name, k, &c.flag_name);
                let diff = (expected[k] - got[k]).abs();
                if diff > eps {
                    failures.push(format!(
                        "{label} {name}[{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e} > eps {eps:.1e}",
                        expected[k], got[k]
                    ));
                }
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
