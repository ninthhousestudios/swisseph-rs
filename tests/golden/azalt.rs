use serde::Deserialize;
use swisseph::Ephemeris;
use swisseph::azalt::{self, AzAltDir, HorDir, RefracDir};

#[derive(Deserialize)]
struct RefracCase {
    inalt: f64,
    atpress: f64,
    attemp: f64,
    dir: String,
    out: f64,
}

#[derive(Deserialize)]
struct RefracExtCase {
    inalt: f64,
    geoalt: f64,
    atpress: f64,
    attemp: f64,
    lapse_rate: f64,
    dir: String,
    out: f64,
    dret: [f64; 4],
}

#[derive(Deserialize)]
struct AzaltCase {
    tjd_ut: f64,
    geopos: [f64; 3],
    dir: String,
    xin: [f64; 2],
    xaz: [f64; 3],
}

#[derive(Deserialize)]
struct AzaltRevCase {
    tjd_ut: f64,
    geopos: [f64; 3],
    dir: String,
    xin: [f64; 2],
    xout: [f64; 2],
}

#[derive(Deserialize)]
struct GoldenData {
    refrac: Vec<RefracCase>,
    refrac_ext: Vec<RefracExtCase>,
    azalt: Vec<AzaltCase>,
    azalt_rev: Vec<AzaltRevCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("azalt.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn refrac_dir(s: &str) -> RefracDir {
    match s {
        "TrueToApp" => RefracDir::TrueToApp,
        "AppToTrue" => RefracDir::AppToTrue,
        other => panic!("Unknown refrac dir: {other}"),
    }
}

#[test]
fn golden_refrac() {
    let data = load();
    for (i, c) in data.refrac.iter().enumerate() {
        let actual = azalt::refrac(c.inalt, c.atpress, c.attemp, refrac_dir(&c.dir));
        let label = format!(
            "refrac[{i}][inalt={},atpress={},dir={}]",
            c.inalt, c.atpress, c.dir
        );
        if actual.to_bits() != c.out.to_bits() {
            super::assert_f64_eps(&label, c.out, actual, 1e-9);
        }
    }
}

#[test]
fn golden_refrac_extended() {
    let data = load();
    for (i, c) in data.refrac_ext.iter().enumerate() {
        let mut dret = [0.0; 4];
        let actual = azalt::refrac_extended(
            c.inalt,
            c.geoalt,
            c.atpress,
            c.attemp,
            c.lapse_rate,
            refrac_dir(&c.dir),
            &mut dret,
        );
        let label = format!(
            "refrac_ext[{i}][inalt={},geoalt={},atpress={},dir={}]",
            c.inalt, c.geoalt, c.atpress, c.dir
        );
        if actual.to_bits() != c.out.to_bits() {
            super::assert_f64_eps(&label, c.out, actual, 1e-9);
        }
        for (j, &expected) in c.dret.iter().enumerate() {
            let got = dret[j];
            if got.to_bits() != expected.to_bits() {
                super::assert_f64_eps(&format!("{label}.dret[{j}]"), expected, got, 1e-9);
            }
        }
    }
}

#[test]
fn golden_azalt() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.azalt.iter().enumerate() {
        let dir = match c.dir.as_str() {
            "EclToHor" => AzAltDir::EclToHor,
            "EquToHor" => AzAltDir::EquToHor,
            other => panic!("Unknown azalt dir: {other}"),
        };
        let actual = ephe.azalt(c.tjd_ut, dir, c.geopos, 0.0, 15.0, 0.0065, c.xin);
        let label = format!(
            "azalt[{i}][tjd_ut={},dir={},geopos={:?}]",
            c.tjd_ut, c.dir, c.geopos
        );
        for (j, &expected) in c.xaz.iter().enumerate() {
            super::assert_f64_eps(&format!("{label}.xaz[{j}]"), expected, actual[j], 1e-7);
        }
    }
}

#[test]
fn golden_azalt_rev() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.azalt_rev.iter().enumerate() {
        let dir = match c.dir.as_str() {
            "HorToEcl" => HorDir::HorToEcl,
            "HorToEqu" => HorDir::HorToEqu,
            other => panic!("Unknown azalt_rev dir: {other}"),
        };
        let actual = ephe.azalt_rev(c.tjd_ut, dir, c.geopos, c.xin);
        let label = format!(
            "azalt_rev[{i}][tjd_ut={},dir={},geopos={:?}]",
            c.tjd_ut, c.dir, c.geopos
        );
        for (j, &expected) in c.xout.iter().enumerate() {
            super::assert_f64_eps(&format!("{label}.xout[{j}]"), expected, actual[j], 1e-7);
        }
    }
}
