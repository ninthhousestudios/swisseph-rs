use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource, TopoPosition};

#[derive(Deserialize)]
struct CalcCase {
    body: i32,
    #[serde(default)]
    #[allow(dead_code)]
    body_name: String,
    jd: f64,
    flags: u32,
    flag_name: String,
    retflag: u32,
    output: [f64; 6],
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("plmoon.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn ephe_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
}

fn all_plmoon_ids() -> Vec<i32> {
    vec![
        9401, 9402, 9501, 9502, 9503, 9504, 9599, 9601, 9602, 9603, 9604, 9605, 9606, 9607, 9608,
        9699, 9701, 9702, 9703, 9704, 9705, 9799, 9801, 9802, 9808, 9899, 9901, 9902, 9903, 9904,
        9905, 9999,
    ]
}

fn body_from_c_id(id: i32) -> Body {
    if (9000..10000).contains(&id) {
        Body::planet_moon(id - 9000).unwrap()
    } else {
        Body::try_from(id).unwrap()
    }
}

fn tolerance(k: usize, is_jpl: bool) -> f64 {
    if is_jpl {
        // Parent-planet position diverges more between JPL and SWIEPH for outer
        // planets at far-past epochs (the moon offset is from a .se1 file either
        // way — only the parent planet source differs).
        if k >= 3 { 1e-3 } else { 2e-4 }
    } else if k >= 3 {
        1e-7
    } else {
        1e-9
    }
}

#[test]
fn golden_plmoon() {
    let topo = TopoPosition {
        longitude: 8.55,
        latitude: 47.37,
        altitude: 500.0,
    };

    let eph_sweph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        topographic: Some(topo),
        planet_moon_numbers: all_plmoon_ids(),
        ..EphemerisConfig::default()
    })
    .unwrap();

    let eph_jpl = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Jpl,
        ephe_path: Some(ephe_path()),
        topographic: Some(topo),
        planet_moon_numbers: all_plmoon_ids(),
        ..EphemerisConfig::default()
    })
    .expect("JPL ephemeris required (de441.eph in ephe/)");

    let cases = load();
    assert!(
        cases.len() >= 650,
        "expected 650+ cases, got {}",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut equiv_checked = 0;
    let mut cancel_checked = 0;
    let skipped = 0;

    for (i, c) in cases.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let is_jpl = flags.contains(CalcFlags::JPLEPH);

        let eph: &Ephemeris = if is_jpl { &eph_jpl } else { &eph_sweph };

        let body = body_from_c_id(c.body);

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "case {i} body={} jd={:.1} {}: error: {e}",
                    c.body, c.jd, c.flag_name
                ));
                continue;
            }
        };

        let label = format!("case {i} body={} jd={:.1} {}", c.body, c.jd, c.flag_name);

        // retflag check for SWIEPH cases
        if !is_jpl {
            let retflag_expected = CalcFlags::from_bits_truncate(c.retflag);
            let retflag_mask = CalcFlags::SWIEPH
                | CalcFlags::MOSEPH
                | CalcFlags::SPEED
                | CalcFlags::HELCTR
                | CalcFlags::CENTER_BODY;
            if result.flags_used & retflag_mask != retflag_expected & retflag_mask {
                failures.push(format!(
                    "{label}: retflag mismatch: expected {:?}, got {:?}",
                    retflag_expected & retflag_mask,
                    result.flags_used & retflag_mask,
                ));
            }
        }

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            let eps = tolerance(k, is_jpl);
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                    c.output[k], result.data[k]
                ));
            }
        }

        // Track special row types
        if c.flag_name.starts_with("CENTER_BODY_cancel") {
            cancel_checked += 1;
        }
        if c.flag_name.starts_with("COB_equiv") || c.flag_name.starts_with("CENTER_BODY_planet") {
            equiv_checked += 1;
        }
    }

    eprintln!(
        "plmoon: {} cases, {} equiv pairs checked, {} cancellation rows, {} skipped",
        cases.len(),
        equiv_checked,
        cancel_checked,
        skipped
    );

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(200) {
            eprintln!("{f}");
        }
        panic!("{n} element failures (showing first 200)");
    }

    assert!(equiv_checked >= 10, "expected at least 10 equivalence rows");
    assert!(cancel_checked >= 5, "expected at least 5 cancellation rows");
    assert_eq!(skipped, 0, "no cases should be skipped");
}

/// Direct 9pmm IDs with parent below Mars (e.g. 9201 → Mercury) must NOT
/// error — C's main_planet skips the moon-file fetch when ipli < SE_MARS,
/// returning the plain parent planet with CENTER_BODY inert in retflag.
#[test]
fn plmoon_sub_mars_parent_inert() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ..EphemerisConfig::default()
    })
    .unwrap();

    // 9201 → parent = Mercury (raw 2), suffix 01 ≠ 99 → not cancelled, but
    // parent < Mars so no moon file opened. Should succeed (plain Mercury).
    let body = Body::planet_moon(201).unwrap(); // raw 9201
    let flags = CalcFlags::MOSEPH | CalcFlags::SPEED;
    let result = eph
        .calc(2451545.0, body, flags)
        .expect("9201 should not error");

    // CENTER_BODY should survive in flags_used (inert but present)
    assert!(
        result.flags_used.contains(CalcFlags::CENTER_BODY),
        "CENTER_BODY bit should survive for sub-Mars parent"
    );

    // The output should match plain Mercury
    let mercury = eph
        .calc(2451545.0, Body::Mercury, flags)
        .expect("Mercury calc");
    for k in 0..6 {
        let diff = (result.data[k] - mercury.data[k]).abs();
        assert!(
            diff < 1e-9,
            "9201 output[{k}] should match Mercury: diff={diff:.3e}"
        );
    }
}

/// 9099 (Sun COB) → clause (iii) cancels: parent=Sun (0 ≤ 4), suffix=99.
/// Should return plain Sun with CENTER_BODY CLEARED.
#[test]
fn plmoon_9099_sun_cob_cancelled() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ..EphemerisConfig::default()
    })
    .unwrap();

    let body = Body::planet_moon(99).unwrap(); // raw 9099
    let flags = CalcFlags::MOSEPH | CalcFlags::SPEED;
    let result = eph
        .calc(2451545.0, body, flags)
        .expect("9099 should not error");

    assert!(
        !result.flags_used.contains(CalcFlags::CENTER_BODY),
        "CENTER_BODY should be cleared for 9099 (Sun COB cancelled)"
    );

    let sun = eph.calc(2451545.0, Body::Sun, flags).expect("Sun calc");
    for k in 0..6 {
        let diff = (result.data[k] - sun.data[k]).abs();
        assert!(
            diff < 1e-9,
            "9099 output[{k}] should match Sun: diff={diff:.3e}"
        );
    }
}
