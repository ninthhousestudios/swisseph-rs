use serde::Deserialize;
use swisseph::ayanamsa;
use swisseph::config::EphemerisConfig;
use swisseph::flags::CalcFlags;
use swisseph::types::AstroModels;

#[derive(Deserialize)]
struct AyaCase {
    index: i32,
    tjd: f64,
    with_nut: f64,
    no_nut: f64,
    speed: f64,
}

#[derive(Deserialize)]
struct UserCase {
    t0: f64,
    ayan_t0: f64,
    tjd: f64,
    no_nut: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    cases: Vec<AyaCase>,
    user: Vec<UserCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("ayanamsa.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn config_for_index(index: i32) -> EphemerisConfig {
    let mut cfg = EphemerisConfig::default();
    cfg.set_sidereal_mode(index, 0.0, 0.0);
    cfg
}

#[test]
fn golden_ayanamsa_no_nut() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::NONUT;
    for (i, c) in data.cases.iter().enumerate() {
        let cfg = config_for_index(c.index);
        let models = AstroModels::default();
        let result = ayanamsa::get_ayanamsa_ex(&cfg, c.tjd, flags, &models)
            .unwrap_or_else(|e| panic!("case {i} (idx={} tjd={}): {e}", c.index, c.tjd));
        let label = format!("case {i} (idx={} tjd={:.6})", c.index, c.tjd);
        super::assert_f64_eps(&label, c.no_nut, result, 1e-8);
    }
}

#[test]
fn golden_ayanamsa_with_nut() {
    let data = load();
    let flags = CalcFlags::MOSEPH;
    for (i, c) in data.cases.iter().enumerate() {
        let cfg = config_for_index(c.index);
        let models = AstroModels::default();
        let result = ayanamsa::get_ayanamsa_ex_nut(&cfg, c.tjd, flags, &models)
            .unwrap_or_else(|e| panic!("case {i} (idx={} tjd={}): {e}", c.index, c.tjd));
        let label = format!("case {i} (idx={} tjd={:.6})", c.index, c.tjd);
        super::assert_f64_eps(&label, c.with_nut, result, 1e-8);
    }
}

#[test]
fn golden_ayanamsa_speed() {
    let data = load();
    let flags = CalcFlags::MOSEPH;
    for (i, c) in data.cases.iter().enumerate() {
        let cfg = config_for_index(c.index);
        let models = AstroModels::default();
        let result = ayanamsa::get_ayanamsa_with_speed(&cfg, c.tjd, flags, &models)
            .unwrap_or_else(|e| panic!("case {i} (idx={} tjd={}): {e}", c.index, c.tjd));
        let label = format!("case {i} (idx={} tjd={:.6}) speed", c.index, c.tjd);
        super::assert_f64_eps(&label, c.speed, result[1], 1e-7);
    }
}

#[test]
fn golden_ayanamsa_user_mode() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::NONUT;
    for (i, c) in data.user.iter().enumerate() {
        let mut cfg = EphemerisConfig::default();
        cfg.set_sidereal_mode(255, c.t0, c.ayan_t0);
        let models = AstroModels::default();
        let result = ayanamsa::get_ayanamsa_ex(&cfg, c.tjd, flags, &models)
            .unwrap_or_else(|e| panic!("user case {i}: {e}"));
        let label = format!(
            "user case {i} (t0={} ayan_t0={} tjd={})",
            c.t0, c.ayan_t0, c.tjd
        );
        super::assert_f64_eps(&label, c.no_nut, result, 1e-8);
    }
}
