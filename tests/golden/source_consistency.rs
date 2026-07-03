use swisseph::Ephemeris;
use swisseph::config::EphemerisConfig;
use swisseph::flags::CalcFlags;
use swisseph::types::{Body, EphemerisSource};

fn make_moshier() -> Ephemeris {
    Ephemeris::new(EphemerisConfig::default()).unwrap()
}

fn make_swiss() -> Ephemeris {
    let ephe = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("swisseph")
        .join("ephe");
    Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe),
        ..Default::default()
    })
    .unwrap()
}

#[test]
fn swieph_flags_on_moshier_config_clamps_to_moshier() {
    let eph = make_moshier();
    let flags = CalcFlags::SWIEPH | CalcFlags::SPEED;
    let result = eph.calc(2451545.0, Body::Moon, flags).unwrap();
    // No Swiss files loaded → clamps to Moshier, signaled via flags_used
    assert!(
        result.flags_used.contains(CalcFlags::MOSEPH),
        "expected MOSEPH in flags_used, got {:?}",
        result.flags_used
    );
    assert!(
        !result.flags_used.contains(CalcFlags::SWIEPH),
        "SWIEPH should not be in flags_used when clamped to Moshier"
    );
}

#[test]
fn moseph_flags_on_swiss_config_uses_moshier() {
    let eph = make_swiss();
    let flags = CalcFlags::MOSEPH | CalcFlags::SPEED;
    let result = eph.calc(2451545.0, Body::Moon, flags).unwrap();
    assert!(
        result.flags_used.contains(CalcFlags::MOSEPH),
        "expected MOSEPH in flags_used, got {:?}",
        result.flags_used
    );

    // Verify position matches a pure Moshier Ephemeris
    let mosh = make_moshier();
    let mosh_result = mosh
        .calc(2451545.0, Body::Moon, CalcFlags::MOSEPH | CalcFlags::SPEED)
        .unwrap();
    for i in 0..6 {
        super::assert_f64_exact(
            &format!("Moon[{i}] MOSEPH-on-Swiss vs pure Moshier"),
            mosh_result.data[i],
            result.data[i],
        );
    }
}

#[test]
fn moseph_flags_on_swiss_config_deltat_uses_de404() {
    let eph_swiss = make_swiss();
    let eph_mosh = make_moshier();

    // At 1800 AD the DE441 vs DE404 tid_acc gap produces a measurable deltaT
    // difference (~0.3s). calc_ut with MOSEPH flags on a Swiss config should use
    // the Moshier tid_acc (DE404), producing the same result as a pure Moshier
    // Ephemeris.
    let jd_ut_1800 = 2378496.5;
    let flags = CalcFlags::MOSEPH | CalcFlags::SPEED;

    let swiss_result = eph_swiss.calc_ut(jd_ut_1800, Body::Moon, flags).unwrap();
    let mosh_result = eph_mosh.calc_ut(jd_ut_1800, Body::Moon, flags).unwrap();

    for i in 0..6 {
        super::assert_f64_exact(
            &format!("Moon[{i}] calc_ut MOSEPH-on-Swiss vs pure Moshier at 1800"),
            mosh_result.data[i],
            swiss_result.data[i],
        );
    }
}

#[test]
fn matched_flags_identity() {
    // When flags match config, behavior is unchanged — Swiss flags on Swiss config
    let eph = make_swiss();
    let flags = CalcFlags::SWIEPH | CalcFlags::SPEED;
    let result = eph.calc(2451545.0, Body::Moon, flags).unwrap();
    assert!(
        result.flags_used.contains(CalcFlags::SWIEPH),
        "expected SWIEPH in flags_used, got {:?}",
        result.flags_used
    );
}

#[test]
fn no_ephmask_uses_config_default() {
    // No EPHMASK bits → uses config.ephemeris_source
    let eph = make_moshier();
    let flags = CalcFlags::SPEED;
    let result = eph.calc(2451545.0, Body::Moon, flags).unwrap();
    assert!(
        result.flags_used.contains(CalcFlags::MOSEPH),
        "expected MOSEPH (config default) in flags_used, got {:?}",
        result.flags_used
    );
}
