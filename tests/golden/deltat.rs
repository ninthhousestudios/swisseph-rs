use serde::Deserialize;
use swisseph::context::EphemerisConfig;
use swisseph::deltat;
use swisseph::types::*;

#[derive(Deserialize)]
struct DeltaTCase {
    model: String,
    tjd: f64,
    expected: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    deltat: Vec<DeltaTCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("deltat.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn parse_model(name: &str) -> DeltaTModel {
    match name {
        "StephensonMorrison1984" => DeltaTModel::StephensonMorrison1984,
        "Stephenson1997" => DeltaTModel::Stephenson1997,
        "StephensonMorrison2004" => DeltaTModel::StephensonMorrison2004,
        "EspenakMeeus2006" => DeltaTModel::EspenakMeeus2006,
        "StephensonEtc2016" => DeltaTModel::StephensonEtc2016,
        _ => panic!("Unknown delta-T model: {name}"),
    }
}

fn config_for_model(model: DeltaTModel) -> EphemerisConfig {
    EphemerisConfig {
        astro_models: AstroModels {
            delta_t: model,
            ..AstroModels::default()
        },
        ..EphemerisConfig::default()
    }
}

#[test]
fn golden_deltat() {
    let data = load();
    for (i, c) in data.deltat.iter().enumerate() {
        let model = parse_model(&c.model);
        let config = config_for_model(model);
        let result = deltat::calc_deltat(c.tjd, &config);
        let label = format!("case {i} ({} tjd={:.6})", c.model, c.tjd);
        super::assert_f64_exact(&label, c.expected, result);
    }
}
