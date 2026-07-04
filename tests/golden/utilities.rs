use serde::Deserialize;
use swisseph::Ephemeris;
use swisseph::config::EphemerisConfig;
use swisseph::types::Body;

#[derive(Deserialize)]
struct TimeEquCase {
    jd_ut: f64,
    #[serde(rename = "E")]
    e: f64,
}

#[derive(Deserialize)]
struct LmtToLatCase {
    jd_lmt: f64,
    geolon: f64,
    tjd_lat: f64,
}

#[derive(Deserialize)]
struct LatToLmtCase {
    jd_lat: f64,
    geolon: f64,
    tjd_lmt: f64,
}

#[derive(Deserialize)]
struct PlanetNameCase {
    ipl: i32,
    name: String,
}

#[derive(Deserialize)]
struct HouseNameCase {
    hsys: String,
    name: String,
}

#[derive(Deserialize)]
struct AyanamsaNameCase {
    sidm: i32,
    name: Option<String>,
}

#[derive(Deserialize)]
struct CsRoundCase {
    input: i32,
    output: i32,
}

#[derive(Deserialize)]
struct CsTimeCase {
    input: i32,
    suppress_zero: bool,
    output: String,
}

#[derive(Deserialize)]
struct CsLonlatCase {
    input: i32,
    output: String,
}

#[derive(Deserialize)]
struct CsDegCase {
    input: i32,
    output: String,
}

#[derive(Deserialize)]
struct GoldenData {
    time_equ: Vec<TimeEquCase>,
    lmt_to_lat: Vec<LmtToLatCase>,
    lat_to_lmt: Vec<LatToLmtCase>,
    planet_name: Vec<PlanetNameCase>,
    house_name: Vec<HouseNameCase>,
    ayanamsa_name: Vec<AyanamsaNameCase>,
    csroundsec: Vec<CsRoundCase>,
    cs2timestr: Vec<CsTimeCase>,
    cs2lonlatstr: Vec<CsLonlatCase>,
    cs2degstr: Vec<CsDegCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("utilities.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn make_eph() -> Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: swisseph::types::EphemerisSource::Swiss,
        ephe_path: Some(std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")),
        ..EphemerisConfig::default()
    };
    Ephemeris::new(config).expect("Ephemeris::new")
}

#[test]
fn golden_time_equ() {
    let data = load();
    let eph = make_eph();
    assert!(!data.time_equ.is_empty());
    for (i, c) in data.time_equ.iter().enumerate() {
        let result = eph.time_equ(c.jd_ut).expect("time_equ");
        let label = format!("time_equ case {i} jd={:.6}", c.jd_ut);
        super::assert_f64_eps(&label, c.e, result, 1e-9);
    }
}

#[test]
fn golden_lmt_to_lat() {
    let data = load();
    let eph = make_eph();
    assert!(!data.lmt_to_lat.is_empty());
    for (i, c) in data.lmt_to_lat.iter().enumerate() {
        let result = eph.lmt_to_lat(c.jd_lmt, c.geolon).expect("lmt_to_lat");
        let label = format!("lmt_to_lat case {i} jd={:.6} geolon={}", c.jd_lmt, c.geolon);
        super::assert_f64_eps(&label, c.tjd_lat, result, 1e-9);
    }
}

#[test]
fn golden_lat_to_lmt() {
    let data = load();
    let eph = make_eph();
    assert!(!data.lat_to_lmt.is_empty());
    for (i, c) in data.lat_to_lmt.iter().enumerate() {
        let result = eph.lat_to_lmt(c.jd_lat, c.geolon).expect("lat_to_lmt");
        let label = format!("lat_to_lmt case {i} jd={:.6} geolon={}", c.jd_lat, c.geolon);
        super::assert_f64_eps(&label, c.tjd_lmt, result, 1e-9);
    }
}

#[test]
fn golden_planet_name() {
    let data = load();
    let eph = make_eph();
    assert!(!data.planet_name.is_empty());
    for (i, c) in data.planet_name.iter().enumerate() {
        let body = Body::try_from(c.ipl).unwrap_or_else(|_| panic!("invalid ipl {}", c.ipl));
        let result = eph.get_planet_name(body);
        assert_eq!(
            result, c.name,
            "planet_name case {i} ipl={}: expected {:?}, got {:?}",
            c.ipl, c.name, result
        );
    }
}

#[test]
fn golden_house_name() {
    let data = load();
    assert!(!data.house_name.is_empty());
    for (i, c) in data.house_name.iter().enumerate() {
        let ch = c.hsys.as_bytes()[0];
        let hsys = swisseph::types::HouseSystem::try_from(ch)
            .unwrap_or_else(|_| panic!("invalid hsys {:?}", c.hsys));
        let result = hsys.name();
        assert_eq!(
            result, c.name,
            "house_name case {i} hsys={}: expected {:?}, got {:?}",
            c.hsys, c.name, result
        );
    }
}

#[test]
fn golden_ayanamsa_name() {
    let data = load();
    assert!(!data.ayanamsa_name.is_empty());
    for (i, c) in data.ayanamsa_name.iter().enumerate() {
        let mode = swisseph::types::SiderealMode::try_from(c.sidm)
            .unwrap_or_else(|_| panic!("invalid sidm {}", c.sidm));
        let result = mode.name();
        let expected = c.name.as_deref();
        assert_eq!(
            result, expected,
            "ayanamsa_name case {i} sidm={}: expected {:?}, got {:?}",
            c.sidm, expected, result
        );
    }
}

#[test]
fn golden_csroundsec() {
    let data = load();
    assert!(!data.csroundsec.is_empty());
    for (i, c) in data.csroundsec.iter().enumerate() {
        let result = swisseph::format::csroundsec(c.input);
        assert_eq!(
            result, c.output,
            "csroundsec case {i} input={}: expected {}, got {}",
            c.input, c.output, result
        );
    }
}

#[test]
fn golden_cs2timestr() {
    let data = load();
    assert!(!data.cs2timestr.is_empty());
    for (i, c) in data.cs2timestr.iter().enumerate() {
        let result = swisseph::format::cs2timestr(c.input, ':', c.suppress_zero);
        assert_eq!(
            result, c.output,
            "cs2timestr case {i} input={} suppress={}: expected {:?}, got {:?}",
            c.input, c.suppress_zero, c.output, result
        );
    }
}

#[test]
fn golden_cs2lonlatstr() {
    let data = load();
    assert!(!data.cs2lonlatstr.is_empty());
    for (i, c) in data.cs2lonlatstr.iter().enumerate() {
        let result = swisseph::format::cs2lonlatstr(c.input, 'E', 'W');
        assert_eq!(
            result, c.output,
            "cs2lonlatstr case {i} input={}: expected {:?}, got {:?}",
            c.input, c.output, result
        );
    }
}

#[test]
fn golden_cs2degstr() {
    let data = load();
    assert!(!data.cs2degstr.is_empty());
    for (i, c) in data.cs2degstr.iter().enumerate() {
        let result = swisseph::format::cs2degstr(c.input);
        assert_eq!(
            result, c.output,
            "cs2degstr case {i} input={}: expected {:?}, got {:?}",
            c.input, c.output, result
        );
    }
}

#[test]
fn delta_t_userdef() {
    let config = EphemerisConfig {
        delta_t_userdef: Some(0.123),
        ..EphemerisConfig::default()
    };
    let result = swisseph::deltat::calc_deltat(2451545.0, &config);
    assert_eq!(result, 0.123, "delta_t_userdef should short-circuit");

    let config_auto = EphemerisConfig::default();
    let result_auto = swisseph::deltat::calc_deltat(2451545.0, &config_auto);
    assert_ne!(
        result_auto, 0.123,
        "default config should compute deltaT normally"
    );
}
