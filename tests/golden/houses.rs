use serde::Deserialize;
use swisseph::houses::houses_armc;
use swisseph::types::HouseSystem;

#[derive(Deserialize)]
struct AnglesSpecialCase {
    armc: f64,
    geolat: f64,
    eps: f64,
    ascmc: [f64; 8],
    ascmc_speed: [f64; 8],
}

#[derive(Deserialize)]
struct EqualFamilyCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct QuadArithCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct GreatCircleCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct IterativeCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct Gauquelin36Case {
    armc: f64,
    geolat: f64,
    eps: f64,
    // serde's array impl tops out at 32 elements; 36 cusps need a Vec.
    cusps: Vec<f64>,
    cusp_speed: Vec<f64>,
}

#[derive(Deserialize)]
struct ClosedFormMiscCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct SunshineCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    sundec: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct GoldenData {
    angles_special: Vec<AnglesSpecialCase>,
    equal_family: Vec<EqualFamilyCase>,
    quad_arith: Vec<QuadArithCase>,
    great_circle: Vec<GreatCircleCase>,
    iterative: Vec<IterativeCase>,
    gauquelin36: Vec<Gauquelin36Case>,
    closed_form_misc: Vec<ClosedFormMiscCase>,
    sunshine: Vec<SunshineCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("houses.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn parse_hsys(s: &str) -> HouseSystem {
    HouseSystem::try_from(s.as_bytes()[0])
        .unwrap_or_else(|e| panic!("Unknown house system {s}: {e}"))
}

#[test]
fn angles_special() {
    let data = load();
    assert_eq!(
        data.angles_special.len(),
        30,
        "expected 30 golden cases (6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.angles_special.iter().enumerate() {
        // The special points are system-independent; 'P' was used by the C generator,
        // but only Equal-family systems are ported so far, so use Equal here.
        let result = houses_armc(c.armc, c.geolat, c.eps, HouseSystem::Equal, None)
            .unwrap_or_else(|e| panic!("case {i}: houses_armc failed: {e}"));

        let actual = result.ascmc.as_array();
        let actual_speed = result.ascmc_speeds.as_array();
        let label_base = format!(
            "case {i} (armc={:.6} geolat={:.6} eps={:.6})",
            c.armc, c.geolat, c.eps
        );
        for j in 0..8 {
            super::assert_f64_exact(&format!("{label_base} ascmc[{j}]"), c.ascmc[j], actual[j]);
            super::assert_f64_exact(
                &format!("{label_base} ascmc_speed[{j}]"),
                c.ascmc_speed[j],
                actual_speed[j],
            );
        }
    }
}

#[test]
fn quad_arith() {
    let data = load();
    assert_eq!(
        data.quad_arith.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.quad_arith.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        // Porphyry cusp speeds are analytical (linear quadrant-rate interpolation);
        // S/X/M/F use the driver-level finite-difference path, which is not
        // bitwise-exact against C's central difference.
        let speed_eps = if c.hsys == "O" { 1e-9 } else { 1e-7 };

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                speed_eps,
            );
        }
    }
}

#[test]
fn great_circle() {
    let data = load();
    assert_eq!(
        data.great_circle.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.great_circle.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-9,
            );
        }
    }
}

#[test]
fn iterative() {
    let data = load();
    assert_eq!(
        data.iterative.len(),
        84,
        "expected 84 golden cases (2 systems x 6 armc x 7 geolat incl. polar x 1 eps)"
    );
    for (i, c) in data.iterative.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-7,
            );
        }
    }
}

#[test]
fn gauquelin36() {
    let data = load();
    assert_eq!(
        data.gauquelin36.len(),
        42,
        "expected 42 golden cases (6 armc x 7 geolat incl. polar x 1 eps)"
    );
    for (i, c) in data.gauquelin36.iter().enumerate() {
        let result = houses_armc(c.armc, c.geolat, c.eps, HouseSystem::Gauquelin, None)
            .unwrap_or_else(|e| panic!("case {i}: houses_armc failed: {e}"));

        let label_base = format!(
            "case {i} (G armc={:.6} geolat={:.6} eps={:.6})",
            c.armc, c.geolat, c.eps
        );
        for h in 1..=36usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-7,
            );
        }
    }
}

#[test]
fn equal_family() {
    let data = load();
    assert_eq!(
        data.equal_family.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.equal_family.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_exact(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
            );
            super::assert_f64_exact(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
            );
        }
    }
}

#[test]
fn closed_form_misc() {
    let data = load();
    assert_eq!(
        data.closed_form_misc.len(),
        120,
        "expected 120 golden cases (4 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.closed_form_misc.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            // U's cusp speeds are stale pre-switch values (not analytical or finite-diff,
            // see c-ref-houses.md §4.2e) — assert them exactly as C produces, including zeros.
            let speed_eps = if c.hsys == "U" { 0.0 } else { 1e-7 };
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                speed_eps,
            );
        }
    }
}

#[test]
fn sunshine() {
    let data = load();
    assert_eq!(
        data.sunshine.len(),
        60,
        "expected 60 golden cases (2 systems x 6 armc x 5 geolat, 1 sundec per case)"
    );
    for (i, c) in data.sunshine.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, Some(c.sundec))
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6} sundec={:.6})",
            c.hsys, c.armc, c.geolat, c.eps, c.sundec
        );
        // Sunshine is closed-form per house (Treindl directly, Makransky via a quadrant case
        // split); Makransky's case split may need the looser tolerance.
        let cusp_eps = if c.hsys == "i" { 1e-8 } else { 1e-9 };
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                cusp_eps,
            );
            // I/i use the driver-level finite-difference cusp speed path (do_interpol).
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-7,
            );
        }
    }
}

#[test]
fn sunshine_requires_sundec() {
    let err = houses_armc(0.0, 51.5, 23.4392911, HouseSystem::Sunshine, None)
        .expect_err("Sunshine without sundec must error");
    assert!(matches!(err, swisseph::error::Error::CError(_)));
}
