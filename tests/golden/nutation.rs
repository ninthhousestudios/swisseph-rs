use serde::Deserialize;
use swisseph::flags::CalcFlags;
use swisseph::nutation;
use swisseph::types::*;

#[derive(Deserialize)]
struct NutationCase {
    model: String,
    jd: f64,
    #[allow(dead_code)]
    flags: u32,
    dpsi: f64,
    deps: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    nutation: Vec<NutationCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("nutation.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn models_for_nut(name: &str) -> AstroModels {
    let nut = match name {
        "IAU1980" => NutationModel::IAU1980,
        "IAUCorr1987" => NutationModel::IAUCorr1987,
        "IAU2000A" => NutationModel::IAU2000A,
        "IAU2000B" => NutationModel::IAU2000B,
        "Woolard" => NutationModel::Woolard,
        _ => return AstroModels::default(),
    };
    AstroModels {
        nutation: nut,
        ..AstroModels::default()
    }
}

#[test]
fn golden_nutation() {
    let data = load();
    for (i, c) in data.nutation.iter().enumerate() {
        let models = models_for_nut(&c.model);
        let nut = nutation::nutation(c.jd, CalcFlags::empty(), &models);
        let label_dpsi = format!("case {i} ({} jd={}) dpsi", c.model, c.jd);
        let label_deps = format!("case {i} ({} jd={}) deps", c.model, c.jd);
        super::assert_f64_exact(&label_dpsi, c.dpsi, nut.dpsi);
        super::assert_f64_exact(&label_deps, c.deps, nut.deps);
    }
}

#[test]
fn jplhor_approx_v3_offset() {
    use swisseph::constants::*;
    let jd = 2415020.0; // below HORIZONS TJD0
    let flags = CalcFlags::JPLHOR_APPROX;
    let models = AstroModels {
        jplhora_mode: JplHoraMode::V3,
        ..AstroModels::default()
    };
    let nut_approx = nutation::nutation(jd, flags, &models);

    let nut_base = nutation::nutation(
        jd,
        CalcFlags::empty(),
        &AstroModels {
            nutation: NutationModel::IAU1980,
            ..AstroModels::default()
        },
    );
    let expected_dpsi = nut_base.dpsi + DPSI_IAU1980_TJD0 / 3600.0 * DEGTORAD;
    let expected_deps = nut_base.deps + DEPS_IAU1980_TJD0 / 3600.0 * DEGTORAD;
    super::assert_f64_exact("JPLHOR_APPROX_V3 dpsi", expected_dpsi, nut_approx.dpsi);
    super::assert_f64_exact("JPLHOR_APPROX_V3 deps", expected_deps, nut_approx.deps);
}

#[test]
fn jplhor_approx_v2_offset() {
    use swisseph::constants::DEGTORAD;
    let jd = 2451545.0;
    let flags = CalcFlags::JPLHOR_APPROX;
    let models_v2 = AstroModels {
        nutation: NutationModel::IAU2000A,
        jplhora_mode: JplHoraMode::V2,
        ..AstroModels::default()
    };
    let nut_v2 = nutation::nutation(jd, flags, &models_v2);

    let nut_base = nutation::nutation(
        jd,
        CalcFlags::empty(),
        &AstroModels {
            nutation: NutationModel::IAU2000A,
            ..AstroModels::default()
        },
    );
    let expected_dpsi = nut_base.dpsi + (-41.7750 / 3600.0 / 1000.0 * DEGTORAD);
    let expected_deps = nut_base.deps + (-6.8192 / 3600.0 / 1000.0 * DEGTORAD);
    super::assert_f64_exact("JPLHOR_APPROX_V2 dpsi", expected_dpsi, nut_v2.dpsi);
    super::assert_f64_exact("JPLHOR_APPROX_V2 deps", expected_deps, nut_v2.deps);
}

#[test]
fn jplhor_uses_iau1980() {
    let jd = 2451545.0;
    let nut_jplhor = nutation::nutation(jd, CalcFlags::DPSIDEPS_1980, &AstroModels::default());
    let nut_iau1980 = nutation::nutation(
        jd,
        CalcFlags::empty(),
        &AstroModels {
            nutation: NutationModel::IAU1980,
            ..AstroModels::default()
        },
    );
    super::assert_f64_exact("JPLHOR == IAU1980 dpsi", nut_iau1980.dpsi, nut_jplhor.dpsi);
    super::assert_f64_exact("JPLHOR == IAU1980 deps", nut_iau1980.deps, nut_jplhor.deps);
}
