use serde::Deserialize;
use swisseph::{CalcFlags, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct FixstarCase {
    star: String,
    tjd: f64,
    iflag: u32,
    flag_name: String,
    xx: [f64; 6],
    retflag: u32,
    star_out: String,
}

#[derive(Deserialize)]
struct MagCase {
    star: String,
    mag: f64,
    star_out: String,
}

#[derive(Deserialize)]
struct GoldenData {
    fixstar: Vec<FixstarCase>,
    mag: Vec<MagCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("fixstar.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn make_eph() -> swisseph::Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some("../swisseph/ephe".into()),
        ..Default::default()
    };
    swisseph::Ephemeris::new(config).unwrap()
}

fn make_eph_sweph() -> swisseph::Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some("../swisseph/ephe".into()),
        ..Default::default()
    };
    swisseph::Ephemeris::new(config).unwrap()
}

fn make_eph_jpl() -> swisseph::Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Jpl,
        ephe_path: Some("../swisseph/ephe".into()),
        jpl_filename: Some("de441.eph".into()),
        ..Default::default()
    };
    swisseph::Ephemeris::new(config).unwrap()
}

/// Shared tolerance logic and assertion runner for fixstar position tests.
/// `cases` is a pre-filtered slice; `label_prefix` is used for failure messages.
fn run_fixstar_cases(eph: &swisseph::Ephemeris, cases: &[FixstarCase]) {
    // position tolerance: 1e-8 degrees (or 1e-10 radians via RADIANS flag)
    // speed tolerance: 1e-6 degrees/day (aberration-speed approximation differs vs C)
    // distance tolerance: 1e-6 relative (very large distances for parallax-less stars)
    // XYZ: 1e-3 AU absolute (relative ~1e-9 for near stars, ~1e-12 for distant)

    let mut failures: Vec<String> = Vec::new();

    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.iflag);
        let has_speed = flags.contains(CalcFlags::SPEED);
        let is_xyz = flags.contains(CalcFlags::XYZ);
        let is_radians = flags.contains(CalcFlags::RADIANS);

        let label = format!("case {i} star={} tjd={:.1} {}", c.star, c.tjd, c.flag_name);

        let result = match eph.fixstar2(&c.star, c.tjd, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{label}: error {e}"));
                continue;
            }
        };

        let (name, calc) = result;

        if name != c.star_out {
            failures.push(format!(
                "{label}: star_out mismatch: got '{name}', expected '{}'",
                c.star_out
            ));
        }

        if calc.flags_used.bits() != c.retflag {
            failures.push(format!(
                "{label}: retflag mismatch: got {}, expected {}",
                calc.flags_used.bits(),
                c.retflag
            ));
        }

        for k in 0..6 {
            let expected = c.xx[k];
            let actual = calc.data[k];

            if k >= 3 && !has_speed {
                continue;
            }

            let (eps, kind) = if is_xyz {
                (1e-3, "AU")
            } else if k == 2 {
                let rel = expected.abs() * 1e-6;
                (rel.max(1e-6), "AU_rel")
            } else if k == 5 {
                // Distance speed: C says "speed is incorrect !!!" and stateless
                // deflection adds additional error. Not meaningful to < ~0.1.
                (0.1, "AU/day_dist")
            } else if k >= 3 {
                if is_radians {
                    (1e-8, "rad/day")
                } else {
                    (1e-6, "deg/day")
                }
            } else if is_radians {
                (1e-10, "rad")
            } else {
                (1e-8, "deg")
            };

            let diff = (expected - actual).abs();
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}] {kind}: expected {expected:.17e}, got {actual:.17e}, diff {diff:.3e} > eps {eps:.3e}"
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(40) {
            eprintln!("{f}");
        }
        panic!("{n} fixstar position failures (first 40 shown)");
    }
}

#[test]
fn golden_fixstar_positions() {
    let eph = make_eph();
    let data = load();

    let cases: Vec<FixstarCase> = data
        .fixstar
        .into_iter()
        .filter(|c| !c.flag_name.starts_with("swieph") && !c.flag_name.starts_with("jpl"))
        .collect();
    let expected_count = 7 * 4 * 7; // NSTARS × NEPOCHS × NFLAGS_MOSHIER
    assert_eq!(
        cases.len(),
        expected_count,
        "unexpected number of Moshier fixstar cases"
    );
    run_fixstar_cases(&eph, &cases);
}

#[test]
fn golden_fixstar_positions_sweph() {
    let eph = make_eph_sweph();
    let data = load();

    let cases: Vec<FixstarCase> = data
        .fixstar
        .into_iter()
        .filter(|c| c.flag_name.starts_with("swieph"))
        .collect();
    let expected_count = 7 * 4 * 3; // NSTARS × NEPOCHS × NFLAGS_SWIEPH
    assert_eq!(
        cases.len(),
        expected_count,
        "unexpected number of SWIEPH fixstar cases"
    );
    run_fixstar_cases(&eph, &cases);
}

#[test]
fn golden_fixstar_positions_jpl() {
    let eph = make_eph_jpl();
    let data = load();

    let cases: Vec<FixstarCase> = data
        .fixstar
        .into_iter()
        .filter(|c| c.flag_name.starts_with("jpl"))
        .collect();
    let expected_count = 7 * 4 * 3; // NSTARS × NEPOCHS × NFLAGS_JPL
    assert_eq!(
        cases.len(),
        expected_count,
        "unexpected number of JPL fixstar cases"
    );
    run_fixstar_cases(&eph, &cases);
}

#[test]
fn golden_fixstar_mag() {
    let eph = make_eph();
    let data = load();

    assert_eq!(data.mag.len(), 4);

    let mut failures: Vec<String> = Vec::new();

    for c in &data.mag {
        match eph.fixstar2_mag(&c.star) {
            Ok((name, mag)) => {
                if name != c.star_out {
                    failures.push(format!(
                        "mag {} star_out: got '{name}', expected '{}'",
                        c.star, c.star_out
                    ));
                }
                let diff = (mag - c.mag).abs();
                if diff > 1e-14 {
                    failures.push(format!(
                        "mag {}: got {mag:.17e}, expected {:.17e}, diff {diff:.3e}",
                        c.star, c.mag
                    ));
                }
            }
            Err(e) => {
                failures.push(format!("mag {}: error {e}", c.star));
            }
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("{f}");
        }
        panic!("{} mag failures", failures.len());
    }
}
