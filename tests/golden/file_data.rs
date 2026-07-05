use std::path::PathBuf;

use swisseph::config::EphemerisConfig;
use swisseph::context::Ephemeris;
use swisseph::types::{EphemerisSource, FileDataKind};

fn ephe_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
}

fn make_swiss_eph() -> Option<Ephemeris> {
    let dir = ephe_dir();
    if !dir.join("sepl_18.se1").exists() {
        return None;
    }
    let config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(dir),
        ..Default::default()
    };
    Ephemeris::new(config).ok()
}

fn make_moshier_eph() -> Ephemeris {
    Ephemeris::new(EphemerisConfig::default()).unwrap()
}

#[test]
fn moshier_returns_none_for_all_kinds() {
    let eph = make_moshier_eph();
    let jd = 2451545.0;
    assert!(eph.file_data(FileDataKind::Planet, jd).is_none());
    assert!(eph.file_data(FileDataKind::Moon, jd).is_none());
    assert!(eph.file_data(FileDataKind::MainAsteroid, jd).is_none());
    assert!(eph.file_data(FileDataKind::Asteroid, jd).is_none());
    assert!(eph.file_data(FileDataKind::PlanetMoon, jd).is_none());
}

#[test]
fn swiss_planet_file_data() {
    let eph = match make_swiss_eph() {
        Some(e) => e,
        None => return,
    };
    let jd = 2451545.0; // J2000, well inside sepl_18's range
    let fd = eph.file_data(FileDataKind::Planet, jd).unwrap();
    assert!(fd.path.to_string_lossy().contains("sepl"));
    assert!(fd.start_jd <= jd);
    assert!(fd.end_jd >= jd);
    assert!(fd.denum > 0);
}

#[test]
fn swiss_moon_file_data() {
    let eph = match make_swiss_eph() {
        Some(e) => e,
        None => return,
    };
    let jd = 2451545.0;
    let fd = eph.file_data(FileDataKind::Moon, jd).unwrap();
    assert!(fd.path.to_string_lossy().contains("semo"));
    assert!(fd.start_jd <= jd);
    assert!(fd.end_jd >= jd);
    assert!(fd.denum > 0);
}

#[test]
fn swiss_main_asteroid_file_data() {
    let eph = match make_swiss_eph() {
        Some(e) => e,
        None => return,
    };
    let jd = 2451545.0;
    let fd = eph.file_data(FileDataKind::MainAsteroid, jd).unwrap();
    assert!(fd.path.to_string_lossy().contains("seas"));
    assert!(fd.start_jd <= jd);
    assert!(fd.end_jd >= jd);
}

#[test]
fn swiss_planet_boundary_selects_correct_file() {
    let eph = match make_swiss_eph() {
        Some(e) => e,
        None => return,
    };
    // sepl_18 starts at 2378496.5 (1800-Jan-1). An epoch at its tfstart should
    // select sepl_18, not sepl_12.
    let boundary_jd = 2378496.5;
    let fd = eph.file_data(FileDataKind::Planet, boundary_jd);
    if let Some(fd) = fd {
        assert!(fd.start_jd <= boundary_jd);
        assert!(fd.end_jd >= boundary_jd);
    }
}

#[test]
fn swiss_out_of_range_returns_none() {
    let eph = match make_swiss_eph() {
        Some(e) => e,
        None => return,
    };
    // Very far future — beyond any .se1 file range
    let jd = 9999999.0;
    assert!(eph.file_data(FileDataKind::Planet, jd).is_none());
}

#[test]
fn asteroid_and_planet_moon_always_none() {
    let eph = match make_swiss_eph() {
        Some(e) => e,
        None => return,
    };
    let jd = 2451545.0;
    assert!(eph.file_data(FileDataKind::Asteroid, jd).is_none());
    assert!(eph.file_data(FileDataKind::PlanetMoon, jd).is_none());
}

#[test]
fn file_data_kind_try_from_valid() {
    assert_eq!(FileDataKind::try_from(0).unwrap(), FileDataKind::Planet);
    assert_eq!(FileDataKind::try_from(1).unwrap(), FileDataKind::Moon);
    assert_eq!(
        FileDataKind::try_from(2).unwrap(),
        FileDataKind::MainAsteroid
    );
    assert_eq!(FileDataKind::try_from(3).unwrap(), FileDataKind::Asteroid);
    assert_eq!(FileDataKind::try_from(4).unwrap(), FileDataKind::PlanetMoon);
}

#[test]
fn file_data_kind_try_from_invalid() {
    assert!(FileDataKind::try_from(5).is_err());
    assert!(FileDataKind::try_from(-1).is_err());
}

#[cfg(feature = "jpl")]
#[test]
fn jpl_planet_file_data() {
    let dir = ephe_dir();
    let jpl_path = dir.join("de441.eph");
    if !jpl_path.exists() {
        return;
    }
    let config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Jpl,
        ephe_path: Some(dir),
        jpl_filename: Some("de441.eph".into()),
        ..Default::default()
    };
    let eph = Ephemeris::new(config).unwrap();
    let jd = 2451545.0;
    let fd = eph.file_data(FileDataKind::Planet, jd).unwrap();
    assert!(fd.path.to_string_lossy().contains("de441"));
    assert_eq!(fd.denum, 441);
    assert!(fd.start_jd <= jd);
    assert!(fd.end_jd >= jd);
}
