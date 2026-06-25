use serde::Deserialize;
use swisseph::flags::CalcFlags;
use swisseph::precession;
use swisseph::types::*;

#[derive(Deserialize)]
struct PrecessionCase {
    model: String,
    direction: String,
    jd: f64,
    flags: u32,
    input: [f64; 3],
    output: [f64; 3],
}

#[derive(Deserialize)]
struct GoldenData {
    precession: Vec<PrecessionCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("precession.json");
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

#[test]
fn golden_precession() {
    let data = load();
    for (i, c) in data.precession.iter().enumerate() {
        if c.direction == "roundtrip" {
            continue;
        }

        let base_model = c
            .model
            .strip_prefix("roundtrip_")
            .unwrap_or(c.model.as_str());

        let (flags, models) = match c.model.as_str() {
            "JPLHOR" => (CalcFlags::DPSIDEPS_1980, AstroModels::default()),
            "JPLHOR_APPROX" => (CalcFlags::JPLHOR_APPROX, AstroModels::default()),
            "mixed_IAU2006short_Vondrak" => (
                CalcFlags::empty(),
                AstroModels {
                    prec_longterm: PrecessionModel::Vondrak2011,
                    prec_shortterm: PrecessionModel::IAU2006,
                    ..AstroModels::default()
                },
            ),
            "mixed_IAU1976short_Laskar" => (
                CalcFlags::empty(),
                AstroModels {
                    prec_longterm: PrecessionModel::Laskar1986,
                    prec_shortterm: PrecessionModel::IAU1976,
                    ..AstroModels::default()
                },
            ),
            "mixed_IAU2000short_Owen" => (
                CalcFlags::empty(),
                AstroModels {
                    prec_longterm: PrecessionModel::Owen1990,
                    prec_shortterm: PrecessionModel::IAU2000,
                    ..AstroModels::default()
                },
            ),
            _ => (CalcFlags::empty(), models_for_prec(base_model)),
        };

        let direction = match c.direction.as_str() {
            "J2000ToDate" => PrecessionDirection::J2000ToDate,
            "DateToJ2000" => PrecessionDirection::DateToJ2000,
            other => panic!("Unknown direction: {other}"),
        };

        let flags = flags | CalcFlags::from_bits_truncate(c.flags);
        let mut pos = c.input;
        precession::precess(&mut pos, c.jd, flags, &models, direction);

        let label = format!(
            "precession[{}][{},{},jd={},flags={}]",
            i, c.model, c.direction, c.jd, c.flags
        );

        let uses_trig = matches!(
            base_model,
            "Vondrak2011"
                | "Owen1990"
                | "Laskar1986"
                | "Simon1994"
                | "Williams1994"
                | "WillEpsLask"
        ) || matches!(
            c.model.as_str(),
            "JPLHOR"
                | "JPLHOR_APPROX"
                | "mixed_IAU2006short_Vondrak"
                | "mixed_IAU1976short_Laskar"
                | "mixed_IAU2000short_Owen"
        );

        for j in 0..3 {
            if uses_trig {
                super::assert_f64_eps(&format!("{label}[{j}]"), c.output[j], pos[j], 1e-15);
            } else {
                super::assert_f64_exact(&format!("{label}[{j}]"), c.output[j], pos[j]);
            }
        }
    }
}

#[test]
fn golden_precession_roundtrip() {
    let data = load();
    for (i, c) in data.precession.iter().enumerate() {
        if c.direction != "roundtrip" {
            continue;
        }

        let base_model = c.model.strip_prefix("roundtrip_").unwrap();
        let models = models_for_prec(base_model);

        let mut pos = c.input;
        precession::precess(
            &mut pos,
            c.jd,
            CalcFlags::empty(),
            &models,
            PrecessionDirection::J2000ToDate,
        );
        precession::precess(
            &mut pos,
            c.jd,
            CalcFlags::empty(),
            &models,
            PrecessionDirection::DateToJ2000,
        );

        let label = format!("roundtrip[{}][{},jd={}]", i, base_model, c.jd);
        for j in 0..3 {
            super::assert_f64_eps(&format!("{label}[{j}]"), c.output[j], pos[j], 1e-14);
        }
    }
}
