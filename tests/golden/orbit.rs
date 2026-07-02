use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct ElementsCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    #[allow(dead_code)]
    retflag: i32,
    dret: [f64; 17],
}

#[derive(Deserialize)]
struct MaxMinCase {
    body: i32,
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    #[allow(dead_code)]
    retflag: i32,
    dmax: f64,
    dmin: f64,
    dtrue: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    elements: Vec<ElementsCase>,
    maxmin: Vec<MaxMinCase>,
}

fn sweph_ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
}

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
            EphemerisSource::Jpl => unreachable!("orbit golden data has no JPL cases"),
        };
        Ephemeris::new(config).expect("Ephemeris::new")
    })
}

fn load() -> GoldenData {
    let path = super::golden_data_path("orbit.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Slots 2..=9 (inclination, node, arg-peri, peri-lon, mean/true/eccentric
/// anomaly, mean longitude) are degree-valued and may straddle the 0/360 wrap.
fn is_angle_slot(k: usize) -> bool {
    (2..=9).contains(&k)
}

/// Signed shortest angular difference (deg), handling the 0/360 wrap.
fn angle_diff(a: f64, b: f64) -> f64 {
    let mut d = (a - b).abs() % 360.0;
    if d > 180.0 {
        d = 360.0 - d;
    }
    d
}

/// Per-slot tolerance for the orbital elements. The internal derivation always
/// runs on TRUEPOS/J2000/NONUT geometry, so the only divergence source is the
/// backend state-vector FP noise (calc golden tests: 1e-9 position / 1e-7
/// speed). Velocity-derived quantities (a, e, and the anomalies computed from
/// `dot_prod(xpos, v)`) therefore carry ~1e-7-relative noise; period/distance
/// slots inherit it too.
///
/// `incl` (the case's inclination, deg) gates a relaxation of the ascending
/// node (slot 3) and argument of perihelion (slot 4): for a near-planar orbit
/// both are ill-conditioned — the node line where the orbit meets the reference
/// ecliptic is barely defined, so a ~1e-9 state-vector difference swings each by
/// up to ~4e-4°. Their SUM, the longitude of perihelion (slot 5), stays
/// well-conditioned and is asserted tight. Earth (inclination ~1e-4° to the
/// J2000 ecliptic) is the only body that hits this; the inclined planets
/// (1°..17°) keep the 1e-6 bound. Same class of C-native ill-conditioning as
/// `nodaps`'s descending-node singularity.
fn elements_tolerance(k: usize, incl: f64) -> f64 {
    // Empirically, every planet (Mercury..Pluto) matches C bit-for-bit on all
    // 17 slots; the only non-zero residual is Earth's, whose heliocentric
    // position this port derives as -(geocentric Sun) rather than through a
    // native SE_EARTH path (max ~4e-8° on the anomalies, ~1e-11 on the metric
    // slots) — so a uniform 1e-6 both accommodates Earth and stays a real
    // regression guard for the planets.
    match k {
        // eccentricity (dimensionless) — bit-exact for planets, ~6e-12 Earth.
        1 => 1e-8,
        // mean daily motion (deg/day) — bit-exact for planets, ~9e-12 Earth.
        11 => 1e-9,
        // ascending node / argument of perihelion — ill-conditioned for a
        // near-planar orbit (see doc comment). Earth's inclination to the J2000
        // ecliptic stays under 0.03° across all epochs; the next-lowest body
        // (Uranus) is 0.77°, so this threshold isolates the near-planar case.
        3 | 4 if incl < 0.1 => 1e-3,
        _ => 1e-6,
    }
}

/// Osculating orbital elements via `Ephemeris::get_orbital_elements`
/// (`swe_get_orbital_elements`): 130 cases — Mercury..Pluto + Earth across
/// Moshier {default, HELCTR, ORBEL_AA} + Swiss {default, HELCTR}, plus the
/// inside-6-AU BARYCTR-gate fallback. Asserts all 17 `dret` slots (angle slots
/// wrap-aware). See [`elements_tolerance`].
#[test]
fn golden_orbit_elements() {
    let data = load();
    assert!(
        data.elements.len() >= 130,
        "expected 130+ cases, got {}",
        data.elements.len()
    );

    let mut cache: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();

    for (i, c) in data.elements.iter().enumerate() {
        let body = Body::try_from(c.body).expect("valid body id");
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let label = format!("case {i} {} tjd={:.1} {}", c.body_name, c.jd, c.flag_name);

        let eph = ephemeris_for(&mut cache, flags);
        let got = match eph.get_orbital_elements(c.jd, body, flags) {
            Ok(e) => e.as_array(),
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        for (k, &g) in got.iter().enumerate() {
            let expected = c.dret[k];
            let eps = elements_tolerance(k, c.dret[2]);
            let diff = if is_angle_slot(k) {
                angle_diff(g, expected)
            } else {
                (g - expected).abs()
            };
            if diff > eps {
                failures.push(format!(
                    "{label} dret[{k}]: expected {expected:.15e}, got {g:.15e}, diff {diff:.3e} > eps {eps:.1e}"
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} orbit-elements failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Max/min/true distance via `Ephemeris::orbit_max_min_true_distance`
/// (`swe_orbit_max_min_true_distance`): 30 Moshier cases — Mercury/Venus/Mars/
/// Jupiter/Pluto × 3 epochs × {geocentric two-ellipse, heliocentric}. `dtrue`
/// is bit-tight; `dmax`/`dmin` come from the 300-iteration refinement.
#[test]
fn golden_orbit_maxmin() {
    let data = load();
    assert!(
        data.maxmin.len() >= 30,
        "expected 30+ cases, got {}",
        data.maxmin.len()
    );

    let mut cache: HashMap<EphemerisSource, Ephemeris> = HashMap::new();
    let mut failures = Vec::new();

    for (i, c) in data.maxmin.iter().enumerate() {
        let body = Body::try_from(c.body).expect("valid body id");
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let label = format!("case {i} {} tjd={:.1} {}", c.body_name, c.jd, c.flag_name);

        let eph = ephemeris_for(&mut cache, flags);
        let (dmax, dmin, dtrue) = match eph.orbit_max_min_true_distance(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error: {e}"));
                continue;
            }
        };

        for (name, expected, got, eps) in [
            // dmax/dmin come from the 300-iteration refinement (actual residual
            // ~1e-10); dtrue is a direct ellipse evaluation (~1e-12).
            ("dmax", c.dmax, dmax, 1e-8),
            ("dmin", c.dmin, dmin, 1e-8),
            ("dtrue", c.dtrue, dtrue, 1e-9),
        ] {
            let diff = (expected - got).abs();
            if diff > eps {
                failures.push(format!(
                    "{label} {name}: expected {expected:.15e}, got {got:.15e}, diff {diff:.3e} > eps {eps:.1e}"
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} orbit-maxmin failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
