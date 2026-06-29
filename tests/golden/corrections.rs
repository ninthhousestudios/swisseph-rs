use serde::Deserialize;

#[derive(Deserialize)]
struct MeffCase {
    r: f64,
    result: f64,
}

#[derive(Deserialize)]
struct AberrCase {
    input: [f64; 6],
    earth: [f64; 6],
    output: [f64; 6],
}

#[derive(Deserialize)]
struct DeflCase {
    label: String,
    input: [f64; 6],
    earth_helio: [f64; 3],
    planet_helio: [f64; 3],
    output: [f64; 3],
}

#[derive(Deserialize)]
struct PipelineCase {
    tjd: f64,
    body: i32,
    body_name: String,
    true_pos: [f64; 6],
    aberr_pos: [f64; 6],
    defl_pos: [f64; 6],
    both_pos: [f64; 6],
}

#[derive(Deserialize)]
struct GoldenData {
    meff: Vec<MeffCase>,
    aberr_light: Vec<AberrCase>,
    deflect_light: Vec<DeflCase>,
    pipeline: Vec<PipelineCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("corrections.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

#[test]
fn golden_meff() {
    use swisseph::corrections::meff;

    let data = load();
    for (i, c) in data.meff.iter().enumerate() {
        let actual = meff(c.r);
        let label = format!("meff case {i} r={:.3}", c.r);
        super::assert_f64_exact(&label, c.result, actual);
    }
    assert_eq!(data.meff.len(), 30);
}

#[test]
fn golden_aberr_light() {
    use swisseph::corrections::aberr_light;

    let data = load();
    for (i, c) in data.aberr_light.iter().enumerate() {
        let mut xx = c.input;
        let earth_vel = [c.earth[3], c.earth[4], c.earth[5]];
        aberr_light(&mut xx, &earth_vel, false);
        let label = format!("aberr case {i}");
        for (k, &x) in xx.iter().enumerate().take(3) {
            super::assert_f64_exact(&format!("{label} [{k}]"), c.output[k], x);
        }
    }
    assert_eq!(data.aberr_light.len(), 40);
}

#[test]
fn golden_deflect_light() {
    use swisseph::corrections::deflect_light;

    let data = load();
    for (i, c) in data.deflect_light.iter().enumerate() {
        let mut xx = c.input;
        let earth6 = [
            c.earth_helio[0],
            c.earth_helio[1],
            c.earth_helio[2],
            0.0,
            0.0,
            0.0,
        ];
        let planet6 = [
            c.planet_helio[0],
            c.planet_helio[1],
            c.planet_helio[2],
            0.0,
            0.0,
            0.0,
        ];
        deflect_light(&mut xx, &earth6, &planet6, false);
        let label = format!("defl case {i} {}", c.label);
        for (k, &x) in xx.iter().enumerate().take(3) {
            super::assert_f64_exact(&format!("{label} [{k}]"), c.output[k], x);
        }
    }
    assert_eq!(data.deflect_light.len(), 12);
}

#[test]
fn golden_pipeline() {
    let data = load();
    assert_eq!(data.pipeline.len(), 15);
    // Pipeline tests will be validated when the full calc pipeline is
    // implemented (swisseph-rs/30). For now, verify the data loads and
    // has the expected structure.
    for (i, c) in data.pipeline.iter().enumerate() {
        assert!(
            c.body >= 2 && c.body <= 6,
            "pipeline case {i}: unexpected body {}",
            c.body
        );
        assert!(c.tjd > 2400000.0, "pipeline case {i}: bad tjd {}", c.tjd);
        let _ = &c.body_name;
        let _ = c.true_pos;
        let _ = c.aberr_pos;
        let _ = c.defl_pos;
        let _ = c.both_pos;
    }
}
