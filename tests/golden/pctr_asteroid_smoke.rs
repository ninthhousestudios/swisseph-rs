use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

fn ephe_path() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("swisseph")
        .join("ephe")
}

fn eph() -> Ephemeris {
    Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    })
    .expect("Ephemeris::new")
}

#[test]
fn pctr_main_asteroids_from_sun() {
    let eph = eph();
    let jd = 2460000.5;
    let flags = CalcFlags::SPEED | CalcFlags::SWIEPH;

    for body in [
        Body::Chiron,
        Body::Pholus,
        Body::Ceres,
        Body::Pallas,
        Body::Juno,
        Body::Vesta,
    ] {
        let r = eph
            .calc_pctr(jd, body, Body::Sun, flags)
            .unwrap_or_else(|e| panic!("calc_pctr({body:?}, Sun) failed: {e}"));
        assert!(
            (0.0..360.0).contains(&r.data[0]),
            "{body:?}: longitude {:.4} out of range",
            r.data[0]
        );
        assert!(r.data[2] > 0.0, "{body:?}: distance should be positive");
    }
}

#[test]
fn pctr_main_asteroids_from_jupiter() {
    let eph = eph();
    let jd = 2460000.5;
    let flags = CalcFlags::SPEED | CalcFlags::SWIEPH;

    for body in [Body::Chiron, Body::Ceres, Body::Vesta] {
        eph.calc_pctr(jd, body, Body::Jupiter, flags)
            .unwrap_or_else(|e| panic!("calc_pctr({body:?}, Jupiter) failed: {e}"));
    }
}

#[test]
fn pctr_asteroid_as_center() {
    let eph = eph();
    let jd = 2460000.5;
    let flags = CalcFlags::SPEED | CalcFlags::SWIEPH;

    eph.calc_pctr(jd, Body::Mars, Body::Ceres, flags)
        .unwrap_or_else(|e| panic!("calc_pctr(Mars, Ceres) failed: {e}"));
}
