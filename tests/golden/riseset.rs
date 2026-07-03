use serde::Deserialize;
use std::path::PathBuf;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource, Error, RiseSetFlags};

#[derive(Deserialize)]
struct FullCase {
    geopos: [f64; 3],
    #[allow(dead_code)]
    geopos_name: String,
    body: String,
    tjd_ut: f64,
    rsmi: String,
    retval: i32,
    tret0: f64,
}

#[derive(Deserialize)]
struct DipCase {
    geopos: [f64; 3],
    #[allow(dead_code)]
    geopos_name: String,
    tjd_ut: f64,
    atpress: f64,
    retval: i32,
    tret0: f64,
}

#[derive(Deserialize)]
struct MtransFlagsCase {
    geopos: [f64; 3],
    #[allow(dead_code)]
    geopos_name: String,
    body: String,
    tjd_ut: f64,
    rsmi: String,
    retval: i32,
    tret0: f64,
}

#[derive(Deserialize)]
struct FastCase {
    geopos: [f64; 3],
    #[allow(dead_code)]
    geopos_name: String,
    body: String,
    tjd_ut: f64,
    rsmi: String,
    retval: i32,
    tret0: f64,
}

#[derive(Deserialize)]
struct AsteroidCase {
    geopos: [f64; 3],
    #[allow(dead_code)]
    geopos_name: String,
    body: String,
    tjd_ut: f64,
    rsmi: String,
    iflag: u32,
    retval: i32,
    tret0: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    full: Vec<FullCase>,
    dip: Vec<DipCase>,
    mtrans_flags: Vec<MtransFlagsCase>,
    fast: Vec<FastCase>,
    asteroid: Vec<AsteroidCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("riseset.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn body_of(s: &str) -> Body {
    match s {
        "Sun" => Body::Sun,
        "Moon" => Body::Moon,
        other => panic!("Unknown body: {other}"),
    }
}

fn rsmi_of(s: &str) -> RiseSetFlags {
    match s {
        "RISE" => RiseSetFlags::RISE,
        "SET" => RiseSetFlags::SET,
        "MTRANSIT" => RiseSetFlags::MTRANSIT,
        "ITRANSIT" => RiseSetFlags::ITRANSIT,
        other => panic!("Unknown rsmi: {other}"),
    }
}

#[test]
fn full() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.full.iter().enumerate() {
        let body = body_of(&c.body);
        let rsmi = rsmi_of(&c.rsmi) | RiseSetFlags::FORCE_SLOW;
        let label = format!(
            "full[{i}][geopos={:?},body={},tjd_ut={},rsmi={}]",
            c.geopos, c.body, c.tjd_ut, c.rsmi
        );
        let actual = ephe.rise_trans_true_hor(
            c.tjd_ut,
            body,
            None,
            CalcFlags::MOSEPH,
            rsmi,
            c.geopos,
            1013.25,
            15.0,
            0.0,
        );
        if c.retval == -2 {
            match actual {
                Err(Error::CircumpolarBody) => {}
                other => panic!("{label}: expected CircumpolarBody, got {other:?}"),
            }
        } else {
            let result = actual.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
            super::assert_f64_eps(&format!("{label}.time"), c.tret0, result.time, 1e-6);
        }
    }
}

/// Covers the `horhgt == -100` auto-dip sentinel combined with `atpress == 0`: calc_dip must
/// receive `atpress` unmodified (not routed through the atpress-auto-estimate used elsewhere),
/// per swecl.c:4415-4416.
#[test]
fn dip() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.dip.iter().enumerate() {
        let label = format!(
            "dip[{i}][geopos={:?},tjd_ut={},atpress={}]",
            c.geopos, c.tjd_ut, c.atpress
        );
        let actual = ephe.rise_trans_true_hor(
            c.tjd_ut,
            Body::Sun,
            None,
            CalcFlags::MOSEPH,
            RiseSetFlags::RISE
                | RiseSetFlags::NO_REFRACTION
                | RiseSetFlags::DISC_CENTER
                | RiseSetFlags::FORCE_SLOW,
            c.geopos,
            c.atpress,
            15.0,
            -100.0,
        );
        if c.retval == -2 {
            match actual {
                Err(Error::CircumpolarBody) => {}
                other => panic!("{label}: expected CircumpolarBody, got {other:?}"),
            }
        } else {
            let result = actual.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
            super::assert_f64_eps(&format!("{label}.time"), c.tret0, result.time, 1e-6);
        }
    }
}

/// `swe_rise_trans`'s fast path (`rise_set_fast`), dispatched via `Ephemeris::rise_trans`.
/// Also sanity-checks that the fast and full algorithms agree to ~1e-5 day for the same inputs.
#[test]
fn fast() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.fast.iter().enumerate() {
        let body = body_of(&c.body);
        let rsmi = rsmi_of(&c.rsmi);
        let label = format!(
            "fast[{i}][geopos={:?},body={},tjd_ut={},rsmi={}]",
            c.geopos, c.body, c.tjd_ut, c.rsmi
        );
        let actual = ephe.rise_trans(
            c.tjd_ut,
            body,
            None,
            CalcFlags::MOSEPH,
            rsmi,
            c.geopos,
            1013.25,
            15.0,
        );
        assert_eq!(c.retval, 0, "{label}: unexpected non-OK C retval");
        let result = actual.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
        super::assert_f64_eps(&format!("{label}.time"), c.tret0, result.time, 1e-6);

        let full = ephe
            .rise_trans_true_hor(
                c.tjd_ut,
                body,
                None,
                CalcFlags::MOSEPH,
                rsmi,
                c.geopos,
                1013.25,
                15.0,
                0.0,
            )
            .unwrap_or_else(|e| panic!("{label}: full algorithm unexpected error {e}"));
        super::assert_f64_eps(
            &format!("{label}.fast_vs_full"),
            full.time,
            result.time,
            1e-5,
        );
    }
}

/// Covers `calc_mer_trans` with `SEFLG_NONUT | SEFLG_TRUEPOS` set on `epheflag`: C masks
/// `epheflag` down to `SEFLG_EPHMASK` only for meridian transits (swecl.c:4701), dropping
/// NONUT/TRUEPOS -- unlike the rise/set branch, which keeps them (swecl.c:4425).
#[test]
fn mtrans_flags() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.mtrans_flags.iter().enumerate() {
        let body = body_of(&c.body);
        let rsmi = rsmi_of(&c.rsmi);
        let label = format!(
            "mtrans_flags[{i}][geopos={:?},body={},tjd_ut={},rsmi={}]",
            c.geopos, c.body, c.tjd_ut, c.rsmi
        );
        let actual = ephe.rise_trans_true_hor(
            c.tjd_ut,
            body,
            None,
            CalcFlags::MOSEPH | CalcFlags::NONUT | CalcFlags::TRUEPOS,
            rsmi,
            c.geopos,
            1013.25,
            15.0,
            0.0,
        );
        if c.retval == -2 {
            match actual {
                Err(Error::CircumpolarBody) => {}
                other => panic!("{label}: expected CircumpolarBody, got {other:?}"),
            }
        } else {
            let result = actual.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
            super::assert_f64_eps(&format!("{label}.time"), c.tret0, result.time, 1e-6);
        }
    }
}

/// Rise/set/transit for numbered asteroid Eros (433) via SWIEPH, exercising `disc_diameter_m`'s
/// asteroid-metadata branch.
#[test]
fn asteroid() {
    let data = load();
    let ephe = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")),
        asteroid_numbers: vec![433],
        ..EphemerisConfig::default()
    })
    .unwrap();
    for (i, c) in data.asteroid.iter().enumerate() {
        let body = Body::asteroid(433).unwrap();
        let rsmi = rsmi_of(&c.rsmi) | RiseSetFlags::FORCE_SLOW;
        let flags = CalcFlags::from_bits_truncate(c.iflag);
        let label = format!(
            "asteroid[{i}][body={},tjd_ut={},rsmi={}]",
            c.body, c.tjd_ut, c.rsmi
        );
        let actual = ephe.rise_trans_true_hor(
            c.tjd_ut, body, None, flags, rsmi, c.geopos, 1013.25, 15.0, 0.0,
        );
        if c.retval == -2 {
            match actual {
                Err(Error::CircumpolarBody) => {}
                other => panic!("{label}: expected CircumpolarBody, got {other:?}"),
            }
        } else {
            let result = actual.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
            super::assert_f64_eps(&format!("{label}.time"), c.tret0, result.time, 1e-6);
        }
    }
}
