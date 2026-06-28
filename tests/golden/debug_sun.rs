// Quick debug: run from project root with cargo test
#[test]
fn debug_sun_positions() {
    use std::path::PathBuf;
    use swisseph::constants::EARTH_MOON_MRAT;
    use swisseph::sweph_file::types::SEI_MOON;
    use swisseph::sweph_file::{SEI_SUNBARY, SwissEphFile, evaluate_body};

    let ephe = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe");
    let pf = SwissEphFile::open(&ephe.join("sepl_18.se1")).unwrap();
    let mf = SwissEphFile::open(&ephe.join("semo_18.se1")).unwrap();

    let jd = 2451545.0; // J2000

    let (emb, _) = evaluate_body(&pf, 0, jd, true).unwrap();
    let (moon, _) = evaluate_body(&mf, SEI_MOON, jd, true).unwrap();
    let (helio_earth_raw, _) = evaluate_body(&pf, SEI_SUNBARY, jd, true).unwrap();

    eprintln!(
        "EMB(body0):     {:+.15e} {:+.15e} {:+.15e}",
        emb[0], emb[1], emb[2]
    );
    eprintln!(
        "Moon(body1):    {:+.15e} {:+.15e} {:+.15e}",
        moon[0], moon[1], moon[2]
    );
    eprintln!(
        "HelioE(body10): {:+.15e} {:+.15e} {:+.15e}",
        helio_earth_raw[0], helio_earth_raw[1], helio_earth_raw[2]
    );

    let mut earth_bary = [0.0f64; 6];
    for i in 0..6 {
        earth_bary[i] = emb[i] - moon[i] / (EARTH_MOON_MRAT + 1.0);
    }

    let mut sun_bary = [0.0f64; 6];
    for i in 0..6 {
        sun_bary[i] = emb[i] - helio_earth_raw[i];
    }

    let mut earth_helio = [0.0f64; 6];
    for i in 0..6 {
        earth_helio[i] = earth_bary[i] - sun_bary[i];
    }

    eprintln!(
        "earth_bary:  {:+.15e} {:+.15e} {:+.15e}",
        earth_bary[0], earth_bary[1], earth_bary[2]
    );
    eprintln!(
        "sun_bary:    {:+.15e} {:+.15e} {:+.15e}",
        sun_bary[0], sun_bary[1], sun_bary[2]
    );
    eprintln!(
        "earth_helio: {:+.15e} {:+.15e} {:+.15e}",
        earth_helio[0], earth_helio[1], earth_helio[2]
    );
    eprintln!(
        "geo_sun(-eh): {:+.15e} {:+.15e} {:+.15e}",
        -earth_helio[0], -earth_helio[1], -earth_helio[2]
    );

    // Compare with Moshier
    use swisseph::flags::CalcFlags;
    use swisseph::moshier::backend::compute_pipeline;
    use swisseph::obliquity::obliquity;
    use swisseph::types::Body;
    use swisseph::types::{AstroModels, Epsilon};
    let eps_j2000 = obliquity(2451545.0, CalcFlags::empty(), &AstroModels::default());
    let pp = compute_pipeline(jd, Body::Sun, &eps_j2000).unwrap();
    eprintln!(
        "Moshier earth_helio: {:+.15e} {:+.15e} {:+.15e}",
        pp.earth_helio[0], pp.earth_helio[1], pp.earth_helio[2]
    );
    eprintln!(
        "Moshier geo_sun:     {:+.15e} {:+.15e} {:+.15e}",
        -pp.earth_helio[0], -pp.earth_helio[1], -pp.earth_helio[2]
    );
}
