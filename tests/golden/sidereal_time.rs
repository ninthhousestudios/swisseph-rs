use serde::Deserialize;
use swisseph::config::EphemerisConfig;
use swisseph::constants::TIDAL_DEFAULT;
use swisseph::sidereal_time;
use swisseph::types::*;

#[derive(Deserialize)]
struct SidtimeCase {
    model: String,
    tjd_ut: f64,
    expected: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    sidtime: Vec<SidtimeCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("sidereal_time.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn parse_model(name: &str) -> SiderealTimeModel {
    match name {
        "IAU1976" => SiderealTimeModel::IAU1976,
        "IAU2006" => SiderealTimeModel::IAU2006,
        "IersConv2010" => SiderealTimeModel::IersConv2010,
        "Longterm" => SiderealTimeModel::Longterm,
        _ => panic!("Unknown sidereal time model: {name}"),
    }
}

// C's swe_sidtime uses swe_deltat_ex(tjd, -1, NULL) which defaults to
// SEFLG_SWIEPH tidal acceleration (TIDAL_DE431). Our EphemerisConfig defaults
// to Moshier (TIDAL_DE404). Pin tidal_acceleration to match C's golden data.
fn config_for_model(model: SiderealTimeModel) -> EphemerisConfig {
    EphemerisConfig {
        tidal_acceleration: Some(TIDAL_DEFAULT),
        astro_models: AstroModels {
            sidereal_time: model,
            ..AstroModels::default()
        },
        ..EphemerisConfig::default()
    }
}

#[test]
fn golden_sidtime() {
    let data = load();
    assert_eq!(
        data.sidtime.len(),
        128,
        "expected 128 golden cases (4 models x 32 epochs)"
    );
    for (i, c) in data.sidtime.iter().enumerate() {
        let model = parse_model(&c.model);
        let config = config_for_model(model);
        let result = sidereal_time::sidereal_time(c.tjd_ut, &config);
        let label = format!("case {i} ({} tjd={:.6})", c.model, c.tjd_ut);
        super::assert_f64_exact(&label, c.expected, result);
    }
}
