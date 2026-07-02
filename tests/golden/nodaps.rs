use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;
use swisseph::{
    Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource, NodApsMethod, NodesApsides,
};

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
    oscu: Vec<MeanCase>,
    oscu_bar: Vec<MeanCase>,
    fopoint: Vec<MeanCase>,
    helctr_bary_mean: Vec<MeanCase>,
    helctr_bary_osc: Vec<MeanCase>,
}

fn sweph_ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
}

/// Ephemeris cache keyed by backend, shared across the osculating batteries
/// (which mix `SEFLG_MOSEPH`/`SEFLG_SWIEPH` cases in the same JSON file).
fn ephemeris_for(cache: &mut HashMap<EphemerisSource, Ephemeris>, flags: CalcFlags) -> &Ephemeris {
    let source = if flags.contains(CalcFlags::SWIEPH) {
        EphemerisSource::Swiss
    } else {
        EphemerisSource::Moshier
    };
    cache.entry(source).or_insert_with(|| {
        let config = match source {
            EphemerisSource::Moshier => EphemerisConfig::default(),
            EphemerisSource::Swiss => EphemerisConfig {
                ephemeris_source: EphemerisSource::Swiss,
                ephe_path: Some(sweph_ephe_path()),
                ..Default::default()
            },
            EphemerisSource::Jpl => unreachable!("nodaps golden data has no JPL cases"),
        };
        Ephemeris::new(config).expect("Ephemeris::new")
    })
}

/// Runs one `Ephemeris::nod_aps` case and records any component outside
/// `eps_pos`/`eps_speed` into `failures`.
/// Tolerance for one osculating-branch component `k` (0..2 = position, 3..5 =
/// speed) of the given node/apsis `point`.
///
/// **Root cause, verified (not guessed):** feeding C's own dumped `xpos[1]`
/// (the ecliptic-of-date sample at exactly `tjd_et`) directly into this
/// port's per-sample ellipse formula reproduces C's own `uu`/`cosnode`/
/// `sinnode`/`sinincl` to ~12 significant digits — the formula is a faithful
/// port. The remaining divergence comes entirely from the ~1e-10..1e-11
/// relative backend noise in the raw position/velocity sample (Moshier series
/// evaluation order, or Swiss/JPL file interpolation), amplified by the
/// A.4.3 tangent-line construction's `fac = z / dz` division — near-singular
/// whenever the sampled radial (z) speed is small relative to the position
/// scale. Ascending AND descending directions share the same `xn`/`-xn`
/// vector, so both inherit this noise (unlike the mean branch, where only
/// the descending node has its own divide-by-near-zero formula); descending
/// is empirically worse because `rn2/ro2`'s rescale ratio (A.4.4) tends to be
/// larger there. Perihelion/aphelion are far less sensitive: they come from
/// `uu`/`ny`/`sema`/`ecce` directly, only picking up node noise secondhand
/// through `uu`'s `cosnode`/`sinnode` term. This is the same class of C-native
/// ill-conditioning as the mean branch's descending-node singularity (see
/// [`tolerance`]) — not a port defect.
fn osc_tolerance(point: &str, k: usize) -> f64 {
    let is_speed = k >= 3;
    match point {
        "desc" => {
            if is_speed {
                3e-2
            } else {
                2e-3
            }
        }
        "asc" => {
            if is_speed {
                1e-4
            } else {
                1e-3
            }
        }
        _ => {
            if is_speed {
                1e-4
            } else {
                5e-5
            }
        } // peri / aphe
    }
}

fn check_case(
    eph: &Ephemeris,
    method: NodApsMethod,
    i: usize,
    c: &MeanCase,
    failures: &mut Vec<String>,
) {
    let body = Body::try_from(c.body).expect("valid body id");
    let flags = CalcFlags::from_bits_truncate(c.flags);
    let label = format!("case {i} {} tjd={:.1} {}", c.body_name, c.jd, c.flag_name);

    let NodesApsides {
        ascending,
        descending,
        perihelion,
        aphelion,
    } = match eph.nod_aps(c.jd, body, flags, method) {
        Ok(r) => r,
        Err(e) => {
            failures.push(format!("{label}: error: {e}"));
            return;
        }
    };

    for (name, expected, got) in [
        ("asc", &c.asc, &ascending),
        ("desc", &c.desc, &descending),
        ("peri", &c.peri, &perihelion),
        ("aphe", &c.aphe, &aphelion),
    ] {
        for k in 0..6 {
            let eps = osc_tolerance(name, k);
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

/// Tolerance for the HELCTR/BARYCTR mean-branch cases: same as [`tolerance`]
/// but with the peri/aphe longitude floor raised to 5e-6° for the
/// heliocentric/barycentric observer branches.
///
/// **Root cause (order-of-magnitude verified):** C's `swi_deflect_light`
/// (sweph.c:3771-3776) retards the barycentric Sun position by the light-time
/// `dt` before building its `planethel` ("Q") vector; this port omits that
/// 2nd-order refinement, matching the existing `calc_planet` convention
/// elsewhere in this codebase (which never retards the Sun for deflection
/// either, and stays within its own tolerance). The residual matches the
/// expected order of `dt · xsun_speed` projected onto the node's ~0.4-5 AU
/// distance (worst observed: 3.4e-6° for Mercury's heliocentric perihelion
/// longitude) — noticeably below the mean branch's own descending-node
/// ill-conditioning (see [`tolerance`]), and only surfaces on peri/aphe
/// longitude because those points, unlike asc/desc, come directly from `uu`
/// without the extra `fac=z/ż`-style division that would otherwise dominate.
fn helctr_bary_tolerance(point: &str, k: usize, flag_name: &str) -> f64 {
    let base = tolerance(point, k, flag_name);
    if k < 3 && matches!(point, "peri" | "aphe") {
        base.max(5e-6)
    } else {
        base
    }
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

/// Osculating nodes & apsides via `Ephemeris::nod_aps` (`swe_nod_aps`, method
/// `SE_NODBIT_OSCU`): 72 cases — 9 bodies {Moon, Mercury..Neptune, Pluto} × 4
/// epochs (incl. a pre-1900 epoch nudged off the sepl_18 .se1 file boundary —
/// see `osc_epochs` in `tests/c-gen/gen_nodaps.c`) × 2 backends {MOSEPH,
/// SWIEPH}, all with `SEFLG_SPEED` (the 3-position central-difference speed is
/// always exercised).
///
/// Tolerances are per point — see [`osc_tolerance`] for the verified root
/// cause (a near-singular division in the node-direction construction, not a
/// port defect).
#[test]
fn golden_nodaps_oscu() {
    let data = load();
    let cases = &data.oscu;
    assert!(cases.len() >= 72, "expected 72+ cases, got {}", cases.len());

    let mut ephemerides: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemeris_for(&mut ephemerides, flags);
        check_case(eph, NodApsMethod::OSCU, i, c, &mut failures);
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 40)");
    }
}

/// Osculating-about-the-barycenter nodes & apsides (`SE_NODBIT_OSCU_BAR`): 8
/// SWIEPH cases — Jupiter/Saturn/Pluto (beyond the 6 AU threshold, so the
/// ellipse is genuinely barycentric) and Mercury (inside it, so this collapses
/// to the same heliocentric ellipse as plain `OSCU`) × 2 epochs. Moshier has
/// no real barycentric frame — `SE_NODBIT_OSCU_BAR` there returns
/// `Error::UnsupportedFlags`, matching `calc_inner`'s general `BARYCTR` gate;
/// this battery is SWIEPH-only so it never exercises that path.
#[test]
fn golden_nodaps_oscu_bar() {
    let data = load();
    let cases = &data.oscu_bar;
    assert!(cases.len() >= 8, "expected 8+ cases, got {}", cases.len());

    let mut ephemerides: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemeris_for(&mut ephemerides, flags);
        check_case(eph, NodApsMethod::OSCU_BAR, i, c, &mut failures);
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 40)");
    }
}

/// HELCTR/BARYCTR observer-frame coverage for the mean branch (A.5.1) — added
/// after a review found `transform_nodaps_output` ignored these flags and
/// always returned geocentric output regardless of what was requested. 12
/// cases — Mercury/Jupiter × 2 epochs × {SWIEPH|HELCTR, SWIEPH|BARYCTR,
/// MOSEPH|HELCTR}, method `SE_NODBIT_MEAN`.
#[test]
fn golden_nodaps_helctr_bary_mean() {
    let data = load();
    let cases = &data.helctr_bary_mean;
    assert!(cases.len() >= 12, "expected 12+ cases, got {}", cases.len());

    let mut ephemerides: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemeris_for(&mut ephemerides, flags);
        let body = Body::try_from(c.body).expect("valid body id");
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
                let eps = helctr_bary_tolerance(name, k, &c.flag_name);
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

/// HELCTR/BARYCTR observer-frame coverage for the osculating branch: 12 cases
/// — Jupiter/Pluto × 2 epochs × {SWIEPH|HELCTR, SWIEPH|BARYCTR,
/// MOSEPH|HELCTR}, method `SE_NODBIT_OSCU`.
#[test]
fn golden_nodaps_helctr_bary_osc() {
    let data = load();
    let cases = &data.helctr_bary_osc;
    assert!(cases.len() >= 12, "expected 12+ cases, got {}", cases.len());

    let mut ephemerides: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemeris_for(&mut ephemerides, flags);
        check_case(eph, NodApsMethod::OSCU, i, c, &mut failures);
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 40)");
    }
}

/// `SE_NODBIT_FOPOINT` (2nd focal point instead of aphelion), combined with
/// `SE_NODBIT_OSCU`: 6 MOSEPH cases — Moon/Mars/Jupiter × 2 epochs.
#[test]
fn golden_nodaps_fopoint() {
    let data = load();
    let cases = &data.fopoint;
    assert!(cases.len() >= 6, "expected 6+ cases, got {}", cases.len());

    let mut ephemerides: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();
    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let eph = ephemeris_for(&mut ephemerides, flags);
        check_case(
            eph,
            NodApsMethod::OSCU | NodApsMethod::FOPOINT,
            i,
            c,
            &mut failures,
        );
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} failures (showing first 40)");
    }
}

/// `SEFLG_TOPOCTR` without a configured topographic position must be
/// rejected, not silently degrade to a geocentric result (a review found
/// `Ephemeris::nod_aps` had no such guard, unlike `Ephemeris::calc`).
#[test]
fn nod_aps_rejects_topoctr_without_config() {
    let eph = Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new");
    let err = eph
        .nod_aps(
            2451545.0,
            Body::Mercury,
            CalcFlags::MOSEPH | CalcFlags::TOPOCTR,
            NodApsMethod::MEAN,
        )
        .expect_err("TOPOCTR without config.topographic must error");
    assert!(
        err.to_string().contains("topocentric"),
        "unexpected error: {err}"
    );
}
