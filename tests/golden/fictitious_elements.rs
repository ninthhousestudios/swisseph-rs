use serde::Deserialize;
use std::path::PathBuf;
use swisseph::fictitious::{
    FictitiousCatalog, kepler, load_fictitious_catalog, osc_el_plan, resolve_elements,
};
use swisseph::types::AstroModels;

#[derive(Deserialize)]
struct Elements {
    tjd0: f64,
    tequ: f64,
    mano: f64,
    sema: f64,
    ecce: f64,
    parg: f64,
    node: f64,
    incl: f64,
}

#[derive(Deserialize)]
struct Case {
    row: usize,
    ipl: i32,
    tjd: f64,
    name: String,
    is_geo: i32,
    elem: Elements,
    xearth: [f64; 6],
    xsun: [f64; 6],
    xp: [f64; 6],
}

fn load() -> Vec<Case> {
    let path = super::golden_data_path("fictitious_elements.json");
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

/// Golden test for resolved elements (bitwise-exact for pure arithmetic)
/// and osc_el_plan output (1e-9 position, 1e-7 velocity).
#[test]
fn golden_fictitious_elements() {
    let cases = load();
    assert!(
        cases.len() >= 114,
        "expected 114+ cases, got {}",
        cases.len()
    );

    let path = ephe_path();
    if !path.join("seorbel.txt").exists() {
        eprintln!("Skipping: seorbel.txt not found at {}", path.display());
        return;
    }
    let catalog = load_fictitious_catalog(Some(&path));
    let models = AstroModels::default();

    let mut elem_pass = 0;
    let mut osc_pass = 0;

    for case in &cases {
        let label = format!("ipl={} ({}) tjd={:.1}", case.ipl, case.name, case.tjd);

        // ── Element resolution ──────────────────────────────────────
        let elem = resolve_elements(&catalog, case.row, case.tjd)
            .unwrap_or_else(|e| panic!("{label}: resolve_elements failed: {e}"));

        // Name
        assert_eq!(
            elem.name, case.name,
            "{label}: name mismatch: got '{}', expected '{}'",
            elem.name, case.name
        );

        // is_geo
        assert_eq!(elem.is_geo, case.is_geo != 0, "{label}: is_geo mismatch");

        // Element values — bitwise exact for pure polynomial arithmetic
        super::assert_f64_exact(&format!("{label} tjd0"), case.elem.tjd0, elem.tjd0);
        super::assert_f64_exact(&format!("{label} tequ"), case.elem.tequ, elem.tequ);
        super::assert_f64_exact(&format!("{label} sema"), case.elem.sema, elem.sema);
        super::assert_f64_exact(&format!("{label} ecce"), case.elem.ecce, elem.ecce);
        super::assert_f64_exact(&format!("{label} mano"), case.elem.mano, elem.mano);
        super::assert_f64_exact(&format!("{label} parg"), case.elem.parg, elem.parg);
        super::assert_f64_exact(&format!("{label} node"), case.elem.node, elem.node);
        super::assert_f64_exact(&format!("{label} incl"), case.elem.incl, elem.incl);

        elem_pass += 1;

        // ── osc_el_plan ─────────────────────────────────────────────
        let xp = osc_el_plan(
            case.tjd,
            &catalog,
            case.row,
            &case.xearth,
            &case.xsun,
            &models,
        )
        .unwrap_or_else(|e| panic!("{label}: osc_el_plan failed: {e}"));

        let pos_eps = 1e-9;
        let vel_eps = 1e-7;

        for i in 0..3 {
            let comp = ["x", "y", "z"][i];
            super::assert_f64_eps(&format!("{label} xp[{comp}]"), case.xp[i], xp[i], pos_eps);
        }
        for i in 3..6 {
            let comp = ["vx", "vy", "vz"][i - 3];
            super::assert_f64_eps(&format!("{label} xp[{comp}]"), case.xp[i], xp[i], vel_eps);
        }

        osc_pass += 1;
    }

    eprintln!(
        "fictitious_elements: {elem_pass} element-resolution + {osc_pass} osc_el_plan cases passed"
    );
}

/// Bodies 55-58 (Vulcan, White Moon, Proserpina, Waldemath) must fail without seorbel.txt.
#[test]
fn fictitious_bodies_55_58_need_file() {
    let catalog = FictitiousCatalog::builtin();
    for row in 15..=18 {
        let result = resolve_elements(&catalog, row, 2451545.0);
        assert!(
            result.is_err(),
            "Expected error for row {row} (ipl {}) with builtin catalog",
            row + 40
        );
    }
}
