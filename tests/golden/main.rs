mod calc;
mod calc_sweph;
mod corrections;
mod date;
mod deltat;
mod jpl_header;
mod jpl_pleph;
mod math;
mod mean_elements;
mod moshier_backend;
mod moshier_moon;
mod moshier_planet;
mod nutation;
mod obliquity_bias;
mod precession;
mod se1_header;
mod sidereal_time;
mod sweph_eval;

use std::path::PathBuf;

fn golden_data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("golden-data")
        .join(name)
}

fn assert_f64_exact(label: &str, expected: f64, actual: f64) {
    assert!(
        expected.to_bits() == actual.to_bits(),
        "{label}: expected {expected:.20} (bits {:016x}), got {actual:.20} (bits {:016x})",
        expected.to_bits(),
        actual.to_bits(),
    );
}

#[allow(dead_code)]
fn assert_f64_eps(label: &str, expected: f64, actual: f64, eps: f64) {
    let diff = (expected - actual).abs();
    assert!(
        diff <= eps,
        "{label}: expected {expected:.20}, got {actual:.20}, diff {diff:.20e} > eps {eps:.20e}",
    );
}
mod debug_sun;
