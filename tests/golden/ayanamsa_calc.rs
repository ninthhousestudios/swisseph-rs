use serde::Deserialize;
use swisseph::config::EphemerisConfig;
use swisseph::context::Ephemeris;
use swisseph::flags::CalcFlags;
use swisseph::types::Body;

#[derive(Deserialize)]
struct Case {
    index: i32,
    body: String,
    tjd: f64,
    lon: f64,
    lat: f64,
    dist: f64,
    lon_speed: f64,
}

#[derive(Deserialize)]
struct EquCase {
    tjd: f64,
    sid_ra: f64,
    trop_ra: f64,
    sid_dec: f64,
    trop_dec: f64,
    dist: f64,
}

#[derive(Deserialize)]
struct NoIndexCase {
    body: String,
    tjd: f64,
    lon: f64,
    lat: f64,
    dist: f64,
    lon_speed: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    cases: Vec<Case>,
    equ: Vec<EquCase>,
    ecl_t0: Vec<Case>,
    user_ecl_t0: Vec<NoIndexCase>,
    ssy: Vec<NoIndexCase>,
    speed3: Vec<Case>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("ayanamsa_calc.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn body_from_str(s: &str) -> Body {
    match s {
        "Sun" => Body::Sun,
        "Moon" => Body::Moon,
        "Mars" => Body::Mars,
        "Jupiter" => Body::Jupiter,
        "MeanNode" => Body::MeanNode,
        other => panic!("unknown body: {other}"),
    }
}

fn eph_for_index(index: i32) -> Ephemeris {
    let mut cfg = EphemerisConfig::default();
    cfg.set_sidereal_mode(index, 0.0, 0.0);
    Ephemeris::new(cfg).expect("Ephemeris::new")
}

#[test]
fn golden_ayanamsa_calc_sidereal() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::SIDEREAL | CalcFlags::SPEED;

    for (i, c) in data.cases.iter().enumerate() {
        let eph = eph_for_index(c.index);
        let body = body_from_str(&c.body);
        let result = eph.calc(c.tjd, body, flags).unwrap_or_else(|e| {
            panic!(
                "case {i} (idx={} body={} tjd={}): {e}",
                c.index, c.body, c.tjd
            )
        });

        let label = |field: &str| {
            format!(
                "case {i} idx={} body={} tjd={:.1} {field}",
                c.index, c.body, c.tjd
            )
        };
        super::assert_f64_eps(&label("lon"), c.lon, result.data[0], 1e-9);
        super::assert_f64_eps(&label("lat"), c.lat, result.data[1], 1e-9);
        super::assert_f64_eps(&label("dist"), c.dist, result.data[2], 1e-9);
        super::assert_f64_eps(&label("lon_speed"), c.lon_speed, result.data[3], 1e-7);
    }
}

/// Equatorial output with SIDEREAL set matches C swe_calc golden data.
/// SIDEREAL forces NONUT (via plaus_iflag), so sidereal equatorial differs from
/// fully-nutated tropical equatorial — but equals tropical-with-NONUT equatorial,
/// confirming the projection leaves xreturn[12..24] untouched.
#[test]
fn golden_ayanamsa_calc_equatorial_tropical() {
    let data = load();
    let sid_flags =
        CalcFlags::MOSEPH | CalcFlags::SIDEREAL | CalcFlags::SPEED | CalcFlags::EQUATORIAL;
    let trop_flags = CalcFlags::MOSEPH | CalcFlags::SPEED | CalcFlags::EQUATORIAL;

    for (i, c) in data.equ.iter().enumerate() {
        let eph_sid = eph_for_index(1); // Lahiri
        let eph_trop = Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new");

        let sid = eph_sid
            .calc(c.tjd, Body::Sun, sid_flags)
            .unwrap_or_else(|e| panic!("equ case {i} sidereal: {e}"));
        let trop = eph_trop
            .calc(c.tjd, Body::Sun, trop_flags)
            .unwrap_or_else(|e| panic!("equ case {i} tropical: {e}"));

        let label = |s: &str| format!("equ case {i} tjd={:.1} {s}", c.tjd);

        // Verify against golden C values
        super::assert_f64_eps(&label("sid_ra"), c.sid_ra, sid.data[0], 1e-9);
        super::assert_f64_eps(&label("trop_ra"), c.trop_ra, trop.data[0], 1e-9);
        super::assert_f64_eps(&label("sid_dec"), c.sid_dec, sid.data[1], 1e-9);
        super::assert_f64_eps(&label("trop_dec"), c.trop_dec, trop.data[1], 1e-9);

        // Sidereal equatorial must equal tropical-with-NONUT equatorial.
        // The projection is ecliptic-only; xreturn[12..24] is untouched.
        // (Difference from trop is due to NONUT, not the projection.)
        let nonut_flags =
            CalcFlags::MOSEPH | CalcFlags::SPEED | CalcFlags::EQUATORIAL | CalcFlags::NONUT;
        let nonut = eph_trop
            .calc(c.tjd, Body::Sun, nonut_flags)
            .unwrap_or_else(|e| panic!("equ case {i} nonut: {e}"));
        super::assert_f64_eps(
            &label("equ_ra_proj_unchanged"),
            nonut.data[0],
            sid.data[0],
            1e-12,
        );
        super::assert_f64_eps(
            &label("equ_dec_proj_unchanged"),
            nonut.data[1],
            sid.data[1],
            1e-12,
        );
        super::assert_f64_eps(&label("dist"), c.dist, sid.data[2], 1e-9);
    }
}

/// ECL_T0 modes (18/19/20/34) project onto the ecliptic of epoch t0.
/// swe_set_sid_mode auto-sets SE_SIDBIT_ECL_T0 for these indices.
#[test]
fn golden_ayanamsa_calc_ecl_t0() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::SIDEREAL | CalcFlags::SPEED;

    for (i, c) in data.ecl_t0.iter().enumerate() {
        let mut cfg = EphemerisConfig::default();
        cfg.set_sidereal_mode(c.index, 0.0, 0.0);
        let eph = Ephemeris::new(cfg).expect("Ephemeris::new");
        let body = body_from_str(&c.body);
        let result = eph.calc(c.tjd, body, flags).unwrap_or_else(|e| {
            panic!(
                "ecl_t0 case {i} (idx={} body={} tjd={}): {e}",
                c.index, c.body, c.tjd
            )
        });

        let label = |field: &str| {
            format!(
                "ecl_t0 case {i} idx={} body={} tjd={:.1} {field}",
                c.index, c.body, c.tjd
            )
        };
        super::assert_f64_eps(&label("lon"), c.lon, result.data[0], 1e-9);
        super::assert_f64_eps(&label("lat"), c.lat, result.data[1], 1e-9);
        super::assert_f64_eps(&label("dist"), c.dist, result.data[2], 1e-9);
        super::assert_f64_eps(&label("lon_speed"), c.lon_speed, result.data[3], 1e-7);
    }
}

/// USER mode with ECL_T0 bit: SE_SIDM_USER|SE_SIDBIT_ECL_T0 = 511, t0=J2000, ayan=25°.
#[test]
fn golden_ayanamsa_calc_user_ecl_t0() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::SIDEREAL | CalcFlags::SPEED;

    let mut cfg = EphemerisConfig::default();
    cfg.set_sidereal_mode(255 | 256, 2451545.0, 25.0);
    let eph = Ephemeris::new(cfg).expect("Ephemeris::new");

    for (i, c) in data.user_ecl_t0.iter().enumerate() {
        let body = body_from_str(&c.body);
        let result = eph.calc(c.tjd, body, flags).unwrap_or_else(|e| {
            panic!("user_ecl_t0 case {i} (body={} tjd={}): {e}", c.body, c.tjd)
        });

        let label = |field: &str| {
            format!(
                "user_ecl_t0 case {i} body={} tjd={:.1} {field}",
                c.body, c.tjd
            )
        };
        super::assert_f64_eps(&label("lon"), c.lon, result.data[0], 1e-9);
        super::assert_f64_eps(&label("lat"), c.lat, result.data[1], 1e-9);
        super::assert_f64_eps(&label("dist"), c.dist, result.data[2], 1e-9);
        super::assert_f64_eps(&label("lon_speed"), c.lon_speed, result.data[3], 1e-7);
    }
}

/// SIDEREAL with SEFLG_SPEED3 (no SEFLG_SPEED) must match C's `use_speed3`, which
/// projects each of the three points and differences the projected positions.
/// Covers both a default-branch mode (Lahiri) and an ECL_T0 mode (J2000) — the
/// latter is a regression guard for swisseph-rs/53, where the projection branch
/// previously discarded the 3-point speed and returned lon_speed≈0.
#[test]
fn golden_ayanamsa_calc_speed3() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::SIDEREAL | CalcFlags::SPEED3;

    for (i, c) in data.speed3.iter().enumerate() {
        let eph = eph_for_index(c.index);
        let body = body_from_str(&c.body);
        let result = eph.calc(c.tjd, body, flags).unwrap_or_else(|e| {
            panic!(
                "speed3 case {i} (idx={} body={} tjd={}): {e}",
                c.index, c.body, c.tjd
            )
        });

        let label = |field: &str| {
            format!(
                "speed3 case {i} idx={} body={} tjd={:.1} {field}",
                c.index, c.body, c.tjd
            )
        };
        super::assert_f64_eps(&label("lon"), c.lon, result.data[0], 1e-9);
        super::assert_f64_eps(&label("lat"), c.lat, result.data[1], 1e-9);
        super::assert_f64_eps(&label("dist"), c.dist, result.data[2], 1e-9);
        super::assert_f64_eps(&label("lon_speed"), c.lon_speed, result.data[3], 1e-7);
    }
}

/// SSY_PLANE: SE_SIDM_LAHIRI|SE_SIDBIT_SSY_PLANE = 513, Sun at J2000.
#[test]
fn golden_ayanamsa_calc_ssy_plane() {
    let data = load();
    let flags = CalcFlags::MOSEPH | CalcFlags::SIDEREAL | CalcFlags::SPEED;

    let mut cfg = EphemerisConfig::default();
    cfg.set_sidereal_mode(1 | 512, 0.0, 0.0);
    let eph = Ephemeris::new(cfg).expect("Ephemeris::new");

    for (i, c) in data.ssy.iter().enumerate() {
        let body = body_from_str(&c.body);
        let result = eph
            .calc(c.tjd, body, flags)
            .unwrap_or_else(|e| panic!("ssy case {i} (body={} tjd={}): {e}", c.body, c.tjd));

        let label = |field: &str| format!("ssy case {i} body={} tjd={:.1} {field}", c.body, c.tjd);
        super::assert_f64_eps(&label("lon"), c.lon, result.data[0], 1e-9);
        super::assert_f64_eps(&label("lat"), c.lat, result.data[1], 1e-9);
        super::assert_f64_eps(&label("dist"), c.dist, result.data[2], 1e-9);
        super::assert_f64_eps(&label("lon_speed"), c.lon_speed, result.data[3], 1e-7);
    }
}
