use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{
    Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource, NodApsMethod, TopoPosition,
};

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

#[derive(Deserialize)]
struct PhenoCase {
    ipl: i32,
    jd: f64,
    flags: u32,
    retflag: u32,
    attr: [f64; 6],
}

#[derive(Deserialize)]
struct OrbitCase {
    ipl: i32,
    jd: f64,
    flags: u32,
    #[allow(dead_code)]
    retflag: i32,
    dret: [Option<f64>; 17],
}

#[derive(Deserialize)]
struct MosephGoldenData {
    moseph: Vec<CalcCase>,
    quirks: Vec<CalcCase>,
    pheno: Vec<PhenoCase>,
    orbit: Vec<OrbitCase>,
}

fn load() -> Vec<CalcCase> {
    let path = super::golden_data_path("plmoon.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn load_moseph() -> MosephGoldenData {
    let path = super::golden_data_path("plmoon_moseph.json");
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

// ---------------------------------------------------------------------------
// MOSEPH golden tests (swisseph-rs/127 step 1)
// ---------------------------------------------------------------------------

// MOSEPH plmoon routes through apparent_planet() (the generic Swiss/JPL pipeline),
// while C's MOSEPH plmoon goes through app_pos_etc_plan() which in Rust is the
// structurally different calc_planet + calc_inner path. The two pipelines handle
// light-time, aberration, deflection, and frame conversion differently when fed
// Moshier positions (equatorial-of-date vs J2000, niter=0 vs niter=1 light-time,
// analytic vs re-evaluated speed). Worst case is Io (large planetocentric offset):
// ~2.5e-2° position, ~7e-5 deg/day speed. All SWIEPH/JPLEPH paths match C at
// 1e-9 / 1e-7 — this tolerance only applies to MOSEPH.
fn moseph_tolerance(k: usize) -> f64 {
    if k >= 3 { 1e-4 } else { 5e-2 }
}

#[test]
fn golden_plmoon_moseph() {
    let data = load_moseph();
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9401, 9501, 9599, 9699, 9901],
        ..EphemerisConfig::default()
    })
    .unwrap();

    assert_eq!(data.moseph.len(), 45, "expected 45 MOSEPH cases");
    let mut failures = Vec::new();

    for (i, c) in data.moseph.iter().enumerate() {
        let flags = CalcFlags::from_bits_truncate(c.flags);
        let body = body_from_c_id(c.body);

        let result = match eph.calc(c.jd, body, flags) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "moseph case {i} body={} jd={:.1} {}: error: {e}",
                    c.body, c.jd, c.flag_name
                ));
                continue;
            }
        };

        let label = format!(
            "moseph case {i} body={} jd={:.1} {}",
            c.body, c.jd, c.flag_name
        );

        // retflag: MOSEPH + CENTER_BODY should be set
        let retflag_expected = CalcFlags::from_bits_truncate(c.retflag);
        let retflag_mask =
            CalcFlags::SWIEPH | CalcFlags::MOSEPH | CalcFlags::SPEED | CalcFlags::CENTER_BODY;
        if result.flags_used & retflag_mask != retflag_expected & retflag_mask {
            failures.push(format!(
                "{label}: retflag mismatch: expected {:?}, got {:?}",
                retflag_expected & retflag_mask,
                result.flags_used & retflag_mask,
            ));
        }

        for k in 0..6 {
            let diff = (c.output[k] - result.data[k]).abs();
            let eps = moseph_tolerance(k);
            if diff > eps {
                failures.push(format!(
                    "{label} [{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                    c.output[k], result.data[k]
                ));
            }
        }
    }

    if !failures.is_empty() {
        let n = failures.len();
        for f in failures.iter().take(100) {
            eprintln!("{f}");
        }
        panic!("{n} MOSEPH plmoon failures");
    }
}

// ---------------------------------------------------------------------------
// Quirk golden tests (swisseph-rs/127 step 2): §2 normalization edge cases
// verified against C golden data (bitwise comparison with the plain planet).
// ---------------------------------------------------------------------------

#[test]
fn golden_plmoon_quirks() {
    let data = load_moseph();
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ..EphemerisConfig::default()
    })
    .unwrap();

    assert_eq!(data.quirks.len(), 6, "expected 6 quirk rows (3 pairs)");

    // Process in pairs: [quirk_ipl, plain_planet]
    for pair in data.quirks.chunks(2) {
        let quirk = &pair[0];
        let plain = &pair[1];

        let flags = CalcFlags::from_bits_truncate(quirk.flags);
        let body = body_from_c_id(quirk.body);

        let result = eph
            .calc(2451545.0, body, flags)
            .unwrap_or_else(|e| panic!("quirk {} should not error: {e}", quirk.flag_name));

        // C golden data: quirk output should bitwise-match the plain planet
        for k in 0..6 {
            let diff = (quirk.output[k] - plain.output[k]).abs();
            assert!(
                diff == 0.0,
                "C golden: {} output[{k}] should match plain {}: C diff={diff:.3e}",
                quirk.flag_name,
                plain.flag_name
            );
        }

        // Rust: our output should match the C golden plain planet
        for k in 0..6 {
            let diff = (result.data[k] - plain.output[k]).abs();
            let eps = if k >= 3 { 1e-7 } else { 1e-9 };
            assert!(
                diff <= eps,
                "Rust vs C: {} output[{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                quirk.flag_name,
                plain.output[k],
                result.data[k]
            );
        }

        // retflag CENTER_BODY bit check
        let c_retflag = CalcFlags::from_bits_truncate(quirk.retflag);
        let c_has_cob = c_retflag.contains(CalcFlags::CENTER_BODY);
        let rust_has_cob = result.flags_used.contains(CalcFlags::CENTER_BODY);
        assert_eq!(
            rust_has_cob, c_has_cob,
            "{}: CENTER_BODY bit mismatch: C={c_has_cob}, Rust={rust_has_cob}",
            quirk.flag_name
        );
    }
}

// ---------------------------------------------------------------------------
// Error-path tests (swisseph-rs/127 step 3): Rust-only, no C golden data.
// ---------------------------------------------------------------------------

#[test]
fn plmoon_outside_time_range_errors() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9599],
        ..EphemerisConfig::default()
    })
    .unwrap();

    // Jupiter COB (9599) at jd 2200000.5 (~1050 AD) — outside sepm file range
    // (2378491.5–2524599.5). Should error.
    let body = Body::planet_moon(599).unwrap();
    let flags = CalcFlags::SWIEPH | CalcFlags::SPEED;
    let err = eph
        .calc(2200000.5, body, flags)
        .err()
        .expect("9599 at 1050 AD should error");

    // Plain Jupiter at the same epoch should succeed (planet files span millennia).
    eph.calc(2200000.5, Body::Jupiter, flags)
        .expect("plain Jupiter at 1050 AD should succeed");

    assert!(
        format!("{err:?}").contains("NotAvailable")
            || format!("{err:?}").contains("EphemerisNotAvailable")
            || format!("{err:?}").contains("BeyondEphemerisLimits"),
        "expected an ephemeris-not-available error, got: {err:?}"
    );
}

#[test]
fn plmoon_phobos_outside_mars_range_errors() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9401],
        ..EphemerisConfig::default()
    })
    .unwrap();

    // Phobos (9401) at jd 2380000.5 (~1804 AD) — inside the Jupiter-moon range
    // (2378491.5–2524599.5) but OUTSIDE the Mars-moon range (2415015.5–2469082.5).
    let body = Body::planet_moon(401).unwrap();
    let flags = CalcFlags::SWIEPH | CalcFlags::SPEED;
    eph.calc(2380000.5, body, flags)
        .err()
        .expect("Phobos at 1804 AD should error (outside Mars-moon file range)");
}

#[test]
fn plmoon_unlisted_moon_id_errors() {
    // Moon id not in planet_moon_numbers → should error.
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9501], // only Io
        ..EphemerisConfig::default()
    })
    .unwrap();

    // 9502 (Europa) not listed
    let body = Body::planet_moon(502).unwrap();
    let flags = CalcFlags::SWIEPH | CalcFlags::SPEED;
    eph.calc(2451545.0, body, flags)
        .err()
        .expect("9502 not in planet_moon_numbers should error");
}

#[test]
fn plmoon_center_body_on_unlisted_planet_errors() {
    // CENTER_BODY on Jupiter without 9599 configured
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9501], // Io but NOT 9599
        ..EphemerisConfig::default()
    })
    .unwrap();

    let flags = CalcFlags::SWIEPH | CalcFlags::SPEED | CalcFlags::CENTER_BODY;
    eph.calc(2451545.0, Body::Jupiter, flags)
        .err()
        .expect("Jupiter+CENTER_BODY without 9599 in config should error");
}

// ---------------------------------------------------------------------------
// nod_aps rejection (swisseph-rs/127 step 4)
// ---------------------------------------------------------------------------

#[test]
fn plmoon_nod_aps_rejected() {
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ..EphemerisConfig::default()
    })
    .unwrap();

    let body = Body::planet_moon(501).unwrap(); // 9501 = Io
    let flags = CalcFlags::MOSEPH | CalcFlags::SPEED;
    let err = eph
        .nod_aps(2451545.0, body, flags, NodApsMethod::OSCU)
        .expect_err("nod_aps should reject PlanetMoon");
    assert!(
        format!("{err:?}").contains("InvalidBody"),
        "expected InvalidBody error, got: {err:?}"
    );
}

// ---------------------------------------------------------------------------
// pheno golden tests (swisseph-rs/127 step 5)
// ---------------------------------------------------------------------------

#[test]
fn golden_plmoon_pheno() {
    let data = load_moseph();
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9501, 9599],
        ..EphemerisConfig::default()
    })
    .unwrap();

    assert_eq!(data.pheno.len(), 2, "expected 2 pheno cases");

    for c in &data.pheno {
        let body = body_from_c_id(c.ipl);
        let flags = CalcFlags::from_bits_truncate(c.flags);

        let (result, _retflag) = eph
            .pheno(c.jd, body, flags)
            .unwrap_or_else(|e| panic!("pheno ipl={} should not error: {e}", c.ipl));

        let label = format!("pheno ipl={}", c.ipl);
        let attrs = [
            result.phase_angle,
            result.phase,
            result.elongation,
            result.apparent_diameter,
            result.apparent_magnitude,
            result.horizontal_parallax,
        ];

        for k in 0..6 {
            let diff = (c.attr[k] - attrs[k]).abs();
            // Position-derived attrs (phase_angle, phase, elongation) carry the same
            // MOSEPH pipeline divergence as the calc path (see moseph_tolerance comment).
            // Magnitude/diameter/parallax are zero for plmoon (no table entry) — exact.
            let eps = if c.attr[k] == 0.0 { 0.0 } else { 5e-2 };
            assert!(
                diff <= eps,
                "{label} attr[{k}]: expected {:.15e}, got {:.15e}, diff {diff:.3e}",
                c.attr[k],
                attrs[k]
            );
        }
    }
}

// ---------------------------------------------------------------------------
// orbital elements golden test (swisseph-rs/127 step 6)
// Preserved C quirk: plmoon produces parent-planet-equivalent heliocentric
// elements (zero-mass two-body around the Sun). Not an endorsed feature.
// ---------------------------------------------------------------------------

#[test]
fn golden_plmoon_orbital_elements() {
    let data = load_moseph();
    let eph = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Moshier,
        ephe_path: Some(ephe_path()),
        planet_moon_numbers: vec![9501],
        ..EphemerisConfig::default()
    })
    .unwrap();

    assert_eq!(data.orbit.len(), 1, "expected 1 orbit case");
    let c = &data.orbit[0];

    let body = body_from_c_id(c.ipl);
    let flags = CalcFlags::from_bits_truncate(c.flags);

    let result = eph
        .get_orbital_elements(c.jd, body, flags)
        .unwrap_or_else(|e| panic!("orbital elements ipl={} should not error: {e}", c.ipl));

    let dret = result.as_array();
    for k in 0..17 {
        match c.dret[k] {
            Some(expected) => {
                let diff = (expected - dret[k]).abs();
                assert!(
                    diff <= 1e-6,
                    "orbit dret[{k}]: expected {expected:.15e}, got {:.15e}, diff {diff:.3e}",
                    dret[k]
                );
            }
            None => {
                // C produced NaN — verify Rust also produces NaN.
                assert!(
                    dret[k].is_nan(),
                    "orbit dret[{k}]: expected NaN, got {:.15e}",
                    dret[k]
                );
            }
        }
    }
}
