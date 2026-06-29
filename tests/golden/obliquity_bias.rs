use serde::Deserialize;
use swisseph::flags::CalcFlags;
use swisseph::types::*;
use swisseph::{bias, obliquity};

// ---------------------------------------------------------------------------
// Deserialization types
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct ObliquityCase {
    model: String,
    jd: f64,
    eps: f64,
}

#[derive(Deserialize)]
struct BiasCase {
    bias_model: String,
    direction: String,
    jd: f64,
    flags: u32,
    input: [f64; 6],
    output: [f64; 6],
}

#[derive(Deserialize)]
struct GoldenData {
    obliquity: Vec<ObliquityCase>,
    bias: Vec<BiasCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("obliquity_bias.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn models_for_prec(name: &str) -> AstroModels {
    let prec = match name {
        "IAU1976" => PrecessionModel::IAU1976,
        "Laskar1986" => PrecessionModel::Laskar1986,
        "WillEpsLask" => PrecessionModel::WillEpsLask,
        "Williams1994" => PrecessionModel::Williams1994,
        "Simon1994" => PrecessionModel::Simon1994,
        "IAU2000" => PrecessionModel::IAU2000,
        "Bretagnon2003" => PrecessionModel::Bretagnon2003,
        "IAU2006" => PrecessionModel::IAU2006,
        "Vondrak2011" => PrecessionModel::Vondrak2011,
        "Owen1990" => PrecessionModel::Owen1990,
        "Newcomb" => PrecessionModel::Newcomb,
        _ => return AstroModels::default(),
    };
    AstroModels {
        prec_longterm: prec,
        prec_shortterm: prec,
        ..AstroModels::default()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn golden_obliquity() {
    let data = load();
    for (i, c) in data.obliquity.iter().enumerate() {
        let (flags, models) = match c.model.as_str() {
            "shortterm_IAU2006" => (
                CalcFlags::empty(),
                AstroModels {
                    prec_longterm: PrecessionModel::Vondrak2011,
                    prec_shortterm: PrecessionModel::IAU2006,
                    ..AstroModels::default()
                },
            ),
            "fallback_Vondrak" => (
                CalcFlags::empty(),
                AstroModels {
                    prec_longterm: PrecessionModel::Vondrak2011,
                    prec_shortterm: PrecessionModel::IAU2006,
                    ..AstroModels::default()
                },
            ),
            "Vondrak_JPLHOR_APPROX" => (CalcFlags::JPLHOR_APPROX, AstroModels::default()),
            "JPLHOR_IAU1976" => (CalcFlags::DPSIDEPS_1980, AstroModels::default()),
            "JPLHOR_Owen" => (CalcFlags::DPSIDEPS_1980, AstroModels::default()),
            name => (CalcFlags::empty(), models_for_prec(name)),
        };

        let actual = obliquity::obliquity(c.jd, flags, &models);
        let label = format!("obliquity[{}][{}]({})", i, c.model, c.jd);

        let uses_trig = matches!(
            c.model.as_str(),
            "Vondrak2011" | "Vondrak_JPLHOR_APPROX" | "shortterm_IAU2006" | "fallback_Vondrak"
        );

        if uses_trig {
            super::assert_f64_eps(&label, c.eps, actual.eps, 1e-15);
        } else {
            super::assert_f64_exact(&label, c.eps, actual.eps);
        }
    }
}

#[test]
fn golden_bias() {
    let data = load();
    for (i, c) in data.bias.iter().enumerate() {
        let bias_model = match c.bias_model.as_str() {
            "None" => BiasModel::None,
            "IAU2000" => BiasModel::IAU2000,
            "IAU2006" => BiasModel::IAU2006,
            other => panic!("Unknown bias model: {other}"),
        };
        let direction = match c.direction.as_str() {
            "GcrsToJ2000" => FrameTransform::GcrsToJ2000,
            "J2000ToGcrs" => FrameTransform::J2000ToGcrs,
            other => panic!("Unknown direction: {other}"),
        };

        let flags = CalcFlags::from_bits_truncate(c.flags);
        let models = AstroModels {
            bias: bias_model,
            ..AstroModels::default()
        };

        let mut pos = c.input;
        bias::frame_bias(&mut pos, c.jd, flags, &models, direction);

        let label = format!(
            "bias[{}][{},{},jd={},flags={}]",
            i, c.bias_model, c.direction, c.jd, c.flags
        );
        for (j, &p) in pos.iter().enumerate() {
            super::assert_f64_exact(&format!("{label}[{j}]"), c.output[j], p);
        }
    }
}
