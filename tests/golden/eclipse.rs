use std::path::PathBuf;

use serde::Deserialize;
use swisseph::{Body, CalcFlags, EclipseFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct SolWhereCase {
    tjd_ut: f64,
    nonut: bool,
    retval: i32,
    geopos: [f64; 10],
    attr: [f64; 11],
    dcore: [f64; 7],
}

#[derive(Deserialize)]
struct SolHowCase {
    tjd_ut: f64,
    geopos: [f64; 3],
    retval: i32,
    attr: [f64; 11],
}

#[derive(Deserialize)]
struct SolWhenGlobCase {
    tjd_start: f64,
    backward: bool,
    retval: i32,
    tret: [f64; 10],
}

#[derive(Deserialize)]
struct SolWhenLocCase {
    geopos: [f64; 3],
    tjd_start: f64,
    backward: bool,
    retval: i32,
    tret: [f64; 10],
    attr: [f64; 11],
}

#[derive(Deserialize)]
struct LunHowCase {
    tjd_ut: f64,
    geopos: [f64; 3],
    retval: i32,
    attr: [f64; 11],
}

#[derive(Deserialize)]
struct LunWhenCase {
    tjd_start: f64,
    backward: bool,
    retval: i32,
    tret: [f64; 8],
}

#[derive(Deserialize)]
struct LunWhenLocCase {
    geopos: [f64; 3],
    tjd_start: f64,
    backward: bool,
    retval: i32,
    tret: [f64; 10],
    attr: [f64; 11],
}

#[derive(Deserialize)]
struct OccWhereCase {
    tjd_ut: f64,
    ipl: i32,
    starname: Option<String>,
    retval: i32,
    geopos: [f64; 10],
    dcore: [f64; 7],
}

#[derive(Deserialize)]
struct OccWhenGlobCase {
    tjd_start: f64,
    ipl: i32,
    starname: Option<String>,
    backward: bool,
    retval: i32,
    tret: [f64; 10],
}

#[derive(Deserialize)]
struct OccWhenLocCase {
    geopos: [f64; 3],
    tjd_start: f64,
    ipl: i32,
    starname: Option<String>,
    backward: bool,
    retval: i32,
    tret: [f64; 10],
    attr: [f64; 11],
}

#[derive(Deserialize)]
struct OccWhenGlobIfltypeCase {
    tjd_start: f64,
    ipl: i32,
    starname: Option<String>,
    ifltype: i32,
    retval: i32,
    tret: [f64; 10],
}

#[derive(Deserialize)]
struct GoldenData {
    sol_where: Vec<SolWhereCase>,
    sol_how: Vec<SolHowCase>,
    sol_when_glob: Vec<SolWhenGlobCase>,
    sol_when_loc: Vec<SolWhenLocCase>,
    lun_how: Vec<LunHowCase>,
    lun_when: Vec<LunWhenCase>,
    lun_when_loc: Vec<LunWhenLocCase>,
    occ_where: Vec<OccWhereCase>,
    occ_when_glob: Vec<OccWhenGlobCase>,
    occ_when_loc: Vec<OccWhenLocCase>,
    occ_when_glob_ifltype: Vec<OccWhenGlobIfltypeCase>,
    occ_where_asteroid: Vec<OccWhereCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("eclipse.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Ephemeris with `ephe_path` set (unlike this file's other tests' `Default::default()`) so the
/// occ_where/occ_when_glob Aldebaran cases can load the fixed-star catalog (same setup as
/// `tests/golden/fixstar.rs`).
fn make_eph() -> Ephemeris {
    let config = EphemerisConfig {
        ephe_path: Some("../swisseph/ephe".into()),
        ..Default::default()
    };
    Ephemeris::new(config).unwrap()
}

/// `swe_sol_eclipse_where` (`Ephemeris::sol_eclipse_where`): shadow-cone geometry only. C's
/// `attr[]` beyond index 3 (`dcore[0]`, core shadow diameter) comes from `eclipse_how`, which is
/// a later task (RSE 6, swisseph-rs/73) -- not asserted here.
///
/// `dcore[1..6]` (penumbra diameter, shadow-axis distance, fundamental-plane diameters, cone
/// half-angle cosines) come from `swi_test_eclipse_where_dcore`, a non-static test-only hook
/// added to `../swisseph/swecl.c` (see `tests/c-gen/gen_eclipse.c`) -- C's public
/// `swe_sol_eclipse_where` never exposes these beyond `dcore[0]`, so they have no other oracle
/// (Codex review, swisseph-rs/72).
#[test]
fn sol_where() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.sol_where.iter().enumerate() {
        let label = format!("sol_where[{i}][tjd_ut={},nonut={}]", c.tjd_ut, c.nonut);
        let ifl = CalcFlags::MOSEPH
            | if c.nonut {
                CalcFlags::NONUT
            } else {
                CalcFlags::empty()
            };
        let result = ephe
            .sol_eclipse_where(c.tjd_ut, ifl)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.central_longitude"),
            c.geopos[0],
            result.central_longitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.central_latitude"),
            c.geopos[1],
            result.central_latitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.core_diameter_km"),
            c.attr[3],
            result.core_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.core_diameter_km (dcore[0])"),
            c.dcore[0],
            result.core_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_km"),
            c.dcore[1],
            result.penumbra_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.shadow_axis_distance_km"),
            c.dcore[2],
            result.shadow_axis_distance_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.umbra_diameter_fundamental_km"),
            c.dcore[3],
            result.umbra_diameter_fundamental_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_fundamental_km"),
            c.dcore[4],
            result.penumbra_diameter_fundamental_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.cos_umbra_half_angle"),
            c.dcore[5],
            result.cos_umbra_half_angle,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.cos_penumbra_half_angle"),
            c.dcore[6],
            result.cos_penumbra_half_angle,
            1e-7,
        );
    }
}

/// `swe_sol_eclipse_how` (`Ephemeris::sol_eclipse_how`): local circumstances (magnitude,
/// obscuration, contact geometry, az/alt, NASA magnitude, Saros series/member) at an observer.
/// Same epochs as `sol_where` (incl. the no-eclipse epoch, exercising the horizon-visibility /
/// "no eclipse here" clearing path) crossed with a near-central and an off-track observer.
#[test]
fn sol_how() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.sol_how.iter().enumerate() {
        let label = format!("sol_how[{i}][tjd_ut={},geopos={:?}]", c.tjd_ut, c.geopos);
        let result = ephe
            .sol_eclipse_how(c.tjd_ut, CalcFlags::MOSEPH, c.geopos)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.magnitude"),
            c.attr[0],
            result.magnitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.diameter_ratio"),
            c.attr[1],
            result.diameter_ratio,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.obscuration"),
            c.attr[2],
            result.obscuration,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.core_diameter_km"),
            c.attr[3],
            result.core_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(&format!("{label}.azimuth"), c.attr[4], result.azimuth, 1e-7);
        super::assert_f64_eps(
            &format!("{label}.true_altitude"),
            c.attr[5],
            result.true_altitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.apparent_altitude"),
            c.attr[6],
            result.apparent_altitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.elongation"),
            c.attr[7],
            result.elongation,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.nasa_magnitude"),
            c.attr[8],
            result.nasa_magnitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.saros_series"),
            c.attr[9],
            result.saros_series,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.saros_member"),
            c.attr[10],
            result.saros_member,
            1e-7,
        );
    }
}

/// `swe_sol_eclipse_when_glob` (`Ephemeris::sol_eclipse_when_glob`): global next/previous solar
/// eclipse search. `ifltype` fixed at "all types" (0) per the golden-data generator; `tjd_start`
/// x `backward` battery exercises both search directions from two epochs, landing on a mix of
/// total/annular/partial-noncentral eclipses.
#[test]
fn sol_when_glob() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.sol_when_glob.iter().enumerate() {
        let label = format!(
            "sol_when_glob[{i}][tjd_start={},backward={}]",
            c.tjd_start, c.backward
        );
        let result = ephe
            .sol_eclipse_when_glob(
                c.tjd_start,
                CalcFlags::MOSEPH,
                EclipseFlags::empty(),
                c.backward,
            )
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.tret[0] (time_maximum)"),
            c.tret[0],
            result.time_maximum,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[1] (time_ra_conjunction)"),
            c.tret[1],
            result.time_ra_conjunction,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[2] (time_begin)"),
            c.tret[2],
            result.time_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[3] (time_end)"),
            c.tret[3],
            result.time_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[4] (time_totality_begin)"),
            c.tret[4],
            result.time_totality_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[5] (time_totality_end)"),
            c.tret[5],
            result.time_totality_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[6] (time_centerline_begin)"),
            c.tret[6],
            result.time_centerline_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[7] (time_centerline_end)"),
            c.tret[7],
            result.time_centerline_end,
            1e-5,
        );
    }
}

/// `swe_sol_eclipse_when_loc` (`Ephemeris::sol_eclipse_when_loc`): local next/previous solar
/// eclipse search, visible-from-`geopos` (topocentric, unlike `sol_when_glob`'s geocentric
/// search). `tret[]` index semantics differ from `sol_when_glob`'s: `tret[1]/[4]` = 1st/4th
/// (penumbra) contact, `tret[2]/[3]` = 2nd/3rd (umbra) contact -- see `SolarEclipseLocal`'s doc
/// comments. Two observers (near-central for the sol_where set; Chile, near the 2019/2020
/// tracks) x 2 start epochs x 2 search directions.
#[test]
fn sol_when_loc() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.sol_when_loc.iter().enumerate() {
        let label = format!(
            "sol_when_loc[{i}][geopos={:?},tjd_start={},backward={}]",
            c.geopos, c.tjd_start, c.backward
        );
        let result = ephe
            .sol_eclipse_when_loc(c.tjd_start, CalcFlags::MOSEPH, c.geopos, c.backward)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.tret[0] (time_maximum)"),
            c.tret[0],
            result.time_maximum,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[1] (time_first_contact)"),
            c.tret[1],
            result.time_first_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[2] (time_second_contact)"),
            c.tret[2],
            result.time_second_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[3] (time_third_contact)"),
            c.tret[3],
            result.time_third_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[4] (time_fourth_contact)"),
            c.tret[4],
            result.time_fourth_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[5] (time_sunrise)"),
            c.tret[5],
            result.time_sunrise,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[6] (time_sunset)"),
            c.tret[6],
            result.time_sunset,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.magnitude"),
            c.attr[0],
            result.attr.magnitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.diameter_ratio"),
            c.attr[1],
            result.attr.diameter_ratio,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.obscuration"),
            c.attr[2],
            result.attr.obscuration,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.core_diameter_km"),
            c.attr[3],
            result.attr.core_diameter_km,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.azimuth"),
            c.attr[4],
            result.attr.azimuth,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.true_altitude"),
            c.attr[5],
            result.attr.true_altitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.apparent_altitude"),
            c.attr[6],
            result.attr.apparent_altitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.elongation"),
            c.attr[7],
            result.attr.elongation,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.nasa_magnitude"),
            c.attr[8],
            result.attr.nasa_magnitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.saros_series"),
            c.attr[9],
            result.attr.saros_series,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.saros_member"),
            c.attr[10],
            result.attr.saros_member,
            1e-5,
        );
    }
}

/// `swe_lun_eclipse_how` (`Ephemeris::lun_eclipse_how`): geocentric shadow-cone geometry plus the
/// Moon's azimuth/altitude at a single observer (Zurich). Three epochs: a total eclipse visible
/// from Zurich, a total eclipse geocentrically but with the Moon below Zurich's horizon
/// (exercises the horizon-visibility "no eclipse here" clearing path -- `retval == 0` while
/// `umbral_magnitude`/`saros_*` stay populated from the geocentric geometry), and a small partial
/// eclipse.
#[test]
fn lun_how() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.lun_how.iter().enumerate() {
        let label = format!("lun_how[{i}][tjd_ut={},geopos={:?}]", c.tjd_ut, c.geopos);
        let result = ephe
            .lun_eclipse_how(c.tjd_ut, CalcFlags::MOSEPH, c.geopos)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.attr[0] (umbral_magnitude)"),
            c.attr[0],
            result.umbral_magnitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[1] (penumbral_magnitude)"),
            c.attr[1],
            result.penumbral_magnitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[4] (azimuth)"),
            c.attr[4],
            result.azimuth,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[5] (true_altitude)"),
            c.attr[5],
            result.true_altitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[6] (apparent_altitude)"),
            c.attr[6],
            result.apparent_altitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[7] (distance_from_opposition)"),
            c.attr[7],
            result.distance_from_opposition,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[8] (umbral_magnitude duplicate)"),
            c.attr[8],
            result.umbral_magnitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[9] (saros_series)"),
            c.attr[9],
            result.saros_series,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[10] (saros_member)"),
            c.attr[10],
            result.saros_member,
            1e-7,
        );
    }
}

/// `swe_lun_eclipse_when` (`Ephemeris::lun_eclipse_when`): global next/previous lunar eclipse
/// search, purely geocentric. `ifltype` fixed at "all types" (0) per the golden-data generator;
/// `tjd_start` x `backward` battery lands on a mix of total/partial/penumbral eclipses (incl. a
/// partial-only case where `tret[4]`/`tret[5]` (totality) stay `0.0`, and a penumbral-only case
/// where `tret[2..=5]` all stay `0.0`).
#[test]
fn lun_when() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.lun_when.iter().enumerate() {
        let label = format!(
            "lun_when[{i}][tjd_start={},backward={}]",
            c.tjd_start, c.backward
        );
        let result = ephe
            .lun_eclipse_when(
                c.tjd_start,
                CalcFlags::MOSEPH,
                EclipseFlags::empty(),
                c.backward,
            )
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.tret[0] (time_maximum)"),
            c.tret[0],
            result.time_maximum,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[2] (time_partial_begin)"),
            c.tret[2],
            result.time_partial_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[3] (time_partial_end)"),
            c.tret[3],
            result.time_partial_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[4] (time_totality_begin)"),
            c.tret[4],
            result.time_totality_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[5] (time_totality_end)"),
            c.tret[5],
            result.time_totality_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[6] (time_penumbral_begin)"),
            c.tret[6],
            result.time_penumbral_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[7] (time_penumbral_end)"),
            c.tret[7],
            result.time_penumbral_end,
            1e-5,
        );
    }
}

/// `swe_lun_eclipse_when_loc` (`Ephemeris::lun_eclipse_when_loc`): local next/previous lunar
/// eclipse search, visible-from-`geopos`, with moonrise/moonset clipping. Same `tret[]` index
/// semantics as `lun_when` (unlike solar's differing global/local layouts), plus `tret[8]`/
/// `tret[9]` for moonrise/moonset.
#[test]
fn lun_when_loc() {
    let data = load();
    let ephe = Ephemeris::new(Default::default()).unwrap();
    for (i, c) in data.lun_when_loc.iter().enumerate() {
        let label = format!(
            "lun_when_loc[{i}][geopos={:?},tjd_start={},backward={}]",
            c.geopos, c.tjd_start, c.backward
        );
        let result = ephe
            .lun_eclipse_when_loc(c.tjd_start, CalcFlags::MOSEPH, c.geopos, c.backward)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.tret[0] (time_maximum)"),
            c.tret[0],
            result.time_maximum,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[2] (time_partial_begin)"),
            c.tret[2],
            result.time_partial_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[3] (time_partial_end)"),
            c.tret[3],
            result.time_partial_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[4] (time_totality_begin)"),
            c.tret[4],
            result.time_totality_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[5] (time_totality_end)"),
            c.tret[5],
            result.time_totality_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[6] (time_penumbral_begin)"),
            c.tret[6],
            result.time_penumbral_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[7] (time_penumbral_end)"),
            c.tret[7],
            result.time_penumbral_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[8] (time_moonrise)"),
            c.tret[8],
            result.time_moonrise,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[9] (time_moonset)"),
            c.tret[9],
            result.time_moonset,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[0] (umbral_magnitude)"),
            c.attr[0],
            result.attr.umbral_magnitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[1] (penumbral_magnitude)"),
            c.attr[1],
            result.attr.penumbral_magnitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[4] (azimuth)"),
            c.attr[4],
            result.attr.azimuth,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[5] (true_altitude)"),
            c.attr[5],
            result.attr.true_altitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[6] (apparent_altitude)"),
            c.attr[6],
            result.attr.apparent_altitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[7] (distance_from_opposition)"),
            c.attr[7],
            result.attr.distance_from_opposition,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[9] (saros_series)"),
            c.attr[9],
            result.attr.saros_series,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr[10] (saros_member)"),
            c.attr[10],
            result.attr.saros_member,
            1e-5,
        );
    }
}

/// `swe_lun_occult_where` (`Ephemeris::lun_occult_where`): shadow-cone geometry for two planets
/// (Venus, Mars), one fixed star (Aldebaran, via `starname`; `ipl` is then just C's `-1` sentinel,
/// translated verbatim through `Body::try_from` -- unused once `starname` is set), and
/// numbered-asteroid Pluto (`SE_AST_OFFSET + 134340` -> `Body::Asteroid`), which
/// `lun_occult_where` aliases to `Body::Pluto`; the golden value is C's own aliased-to-Pluto
/// result, so a broken alias fails here (loudly, since MOSEPH has no bare-asteroid ephemeris).
/// `dcore[1..6]` come from the same `swi_test_eclipse_where_dcore` test-only hook `sol_where`
/// uses (generic over `ipl`/`starname` already).
#[test]
fn occ_where() {
    let data = load();
    let ephe = make_eph();
    for (i, c) in data.occ_where.iter().enumerate() {
        let label = format!(
            "occ_where[{i}][tjd_ut={},ipl={},starname={:?}]",
            c.tjd_ut, c.ipl, c.starname
        );
        let body = Body::try_from(c.ipl).unwrap();
        let result = ephe
            .lun_occult_where(c.tjd_ut, body, c.starname.as_deref(), CalcFlags::MOSEPH)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.central_longitude"),
            c.geopos[0],
            result.central_longitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.central_latitude"),
            c.geopos[1],
            result.central_latitude,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.core_diameter_km (dcore[0])"),
            c.dcore[0],
            result.core_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_km"),
            c.dcore[1],
            result.penumbra_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.shadow_axis_distance_km"),
            c.dcore[2],
            result.shadow_axis_distance_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.umbra_diameter_fundamental_km"),
            c.dcore[3],
            result.umbra_diameter_fundamental_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_fundamental_km"),
            c.dcore[4],
            result.penumbra_diameter_fundamental_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.cos_umbra_half_angle"),
            c.dcore[5],
            result.cos_umbra_half_angle,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.cos_penumbra_half_angle"),
            c.dcore[6],
            result.cos_penumbra_half_angle,
            1e-7,
        );
    }
}

/// `swe_lun_occult_when_glob` (`Ephemeris::lun_occult_when_glob`): global occultation search for
/// the same four occulted bodies as `occ_where`, `ifltype = 0` (all types valid for the occulted
/// body), both search directions.
#[test]
fn occ_when_glob() {
    let data = load();
    let ephe = make_eph();
    for (i, c) in data.occ_when_glob.iter().enumerate() {
        let label = format!(
            "occ_when_glob[{i}][tjd_start={},ipl={},starname={:?},backward={}]",
            c.tjd_start, c.ipl, c.starname, c.backward
        );
        let body = Body::try_from(c.ipl).unwrap();
        let result = ephe
            .lun_occult_when_glob(
                c.tjd_start,
                body,
                c.starname.as_deref(),
                CalcFlags::MOSEPH,
                EclipseFlags::empty(),
                c.backward,
            )
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.tret[0] (time_maximum)"),
            c.tret[0],
            result.time_maximum,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[1] (time_ra_conjunction)"),
            c.tret[1],
            result.time_ra_conjunction,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[2] (time_begin)"),
            c.tret[2],
            result.time_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[3] (time_end)"),
            c.tret[3],
            result.time_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[4] (time_totality_begin)"),
            c.tret[4],
            result.time_totality_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[5] (time_totality_end)"),
            c.tret[5],
            result.time_totality_end,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[6] (time_centerline_begin)"),
            c.tret[6],
            result.time_centerline_begin,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[7] (time_centerline_end)"),
            c.tret[7],
            result.time_centerline_end,
            1e-5,
        );
    }
}

/// `swe_lun_occult_when_glob` with a non-empty `ifltype` (swisseph-rs/92 follow-up): exercises the
/// type-filter logic the `ifltype = 0` `occ_when_glob` battery never reaches. Covers the retry that
/// skips non-matching types (`PARTIAL` searching past the nearer total occultations of Venus), the
/// pass-through when the first occultation already matches (`TOTAL`), and the two hard-error returns
/// (`ANNULAR` for a planet, `PARTIAL | CENTRAL`). C signals an error with `retval < 0`; the Rust
/// port returns `Err` for those, so those cases assert an error instead of comparing `tret`.
#[test]
fn occ_when_glob_ifltype() {
    let data = load();
    let ephe = make_eph();
    for (i, c) in data.occ_when_glob_ifltype.iter().enumerate() {
        let label = format!(
            "occ_when_glob_ifltype[{i}][ipl={},starname={:?},ifltype={:#x}]",
            c.ipl, c.starname, c.ifltype
        );
        let body = Body::try_from(c.ipl).unwrap();
        let ifltype = EclipseFlags::from_bits(c.ifltype as u32)
            .unwrap_or_else(|| panic!("{label}: bad ifltype bits {:#x}", c.ifltype));
        let result = ephe.lun_occult_when_glob(
            c.tjd_start,
            body,
            c.starname.as_deref(),
            CalcFlags::MOSEPH,
            ifltype,
            false,
        );

        if c.retval < 0 {
            assert!(
                result.is_err(),
                "{label}: expected error (C retval {}), got {result:?}",
                c.retval
            );
            continue;
        }

        let result = result.unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));
        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        let got = [
            result.time_maximum,
            result.time_ra_conjunction,
            result.time_begin,
            result.time_end,
            result.time_totality_begin,
            result.time_totality_end,
            result.time_centerline_begin,
            result.time_centerline_end,
        ];
        for (slot, &g) in got.iter().enumerate() {
            super::assert_f64_eps(&format!("{label}.tret[{slot}]"), c.tret[slot], g, 1e-5);
        }
    }
}

/// `swe_lun_occult_when_loc` (`Ephemeris::lun_occult_when_loc`): local occultation search visible
/// from Zurich, for the same Venus (finite disc) and Aldebaran (point-source star) bodies as
/// `occ_where`/`occ_when_glob`, both search directions. The star cases exercise the point-source
/// contact-1/4-aliased-from-2/3 branch (swecl.c:2696-2699, c-ref-occultation.md §3 step 9); the
/// Venus/forward case exercises `OCC_BEG_DAYLIGHT`/`OCC_END_DAYLIGHT` (retval bits 8192/16384).
#[test]
fn occ_when_loc() {
    let data = load();
    let ephe = make_eph();
    for (i, c) in data.occ_when_loc.iter().enumerate() {
        let label = format!(
            "occ_when_loc[{i}][geopos={:?},tjd_start={},ipl={},starname={:?},backward={}]",
            c.geopos, c.tjd_start, c.ipl, c.starname, c.backward
        );
        let body = Body::try_from(c.ipl).unwrap();
        let result = ephe
            .lun_occult_when_loc(
                c.tjd_start,
                body,
                c.starname.as_deref(),
                CalcFlags::MOSEPH,
                c.geopos,
                c.backward,
            )
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.tret[0] (time_maximum)"),
            c.tret[0],
            result.time_maximum,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[1] (time_first_contact)"),
            c.tret[1],
            result.time_first_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[2] (time_second_contact)"),
            c.tret[2],
            result.time_second_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[3] (time_third_contact)"),
            c.tret[3],
            result.time_third_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[4] (time_fourth_contact)"),
            c.tret[4],
            result.time_fourth_contact,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[5] (time_rise)"),
            c.tret[5],
            result.time_rise,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.tret[6] (time_set)"),
            c.tret[6],
            result.time_set,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.magnitude"),
            c.attr[0],
            result.attr.magnitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.diameter_ratio"),
            c.attr[1],
            result.attr.diameter_ratio,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.obscuration"),
            c.attr[2],
            result.attr.obscuration,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.core_diameter_km"),
            c.attr[3],
            result.attr.core_diameter_km,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.azimuth"),
            c.attr[4],
            result.attr.azimuth,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.true_altitude"),
            c.attr[5],
            result.attr.true_altitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.apparent_altitude"),
            c.attr[6],
            result.attr.apparent_altitude,
            1e-5,
        );
        super::assert_f64_eps(
            &format!("{label}.attr.elongation"),
            c.attr[7],
            result.attr.elongation,
            1e-5,
        );
    }
}

/// `swe_lun_occult_where` with numbered asteroid Eros (433) via SWIEPH, exercising
/// `body_radius_au`'s asteroid-metadata branch.
#[test]
fn occ_where_asteroid() {
    let data = load();
    let ephe = Ephemeris::new(EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")),
        asteroid_numbers: vec![433],
        ..Default::default()
    })
    .unwrap();
    for (i, c) in data.occ_where_asteroid.iter().enumerate() {
        let label = format!("occ_where_asteroid[{i}][tjd_ut={},ipl={}]", c.tjd_ut, c.ipl);
        let body = Body::try_from(c.ipl).unwrap();
        let result = ephe
            .lun_occult_where(c.tjd_ut, body, c.starname.as_deref(), CalcFlags::SWIEPH)
            .unwrap_or_else(|e| panic!("{label}: unexpected error {e}"));

        assert_eq!(
            c.retval,
            result.flags.bits() as i32,
            "{label}: retval mismatch (expected {:#x}, got {:#x})",
            c.retval,
            result.flags.bits()
        );
        super::assert_f64_eps(
            &format!("{label}.central_longitude"),
            c.geopos[0],
            result.central_longitude,
            2e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.central_latitude"),
            c.geopos[1],
            result.central_latitude,
            2e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.core_diameter_km (dcore[0])"),
            c.dcore[0],
            result.core_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_km"),
            c.dcore[1],
            result.penumbra_diameter_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.shadow_axis_distance_km"),
            c.dcore[2],
            result.shadow_axis_distance_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.umbra_diameter_fundamental_km"),
            c.dcore[3],
            result.umbra_diameter_fundamental_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_fundamental_km"),
            c.dcore[4],
            result.penumbra_diameter_fundamental_km,
            1e-2,
        );
        super::assert_f64_eps(
            &format!("{label}.cos_umbra_half_angle"),
            c.dcore[5],
            result.cos_umbra_half_angle,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.cos_penumbra_half_angle"),
            c.dcore[6],
            result.cos_penumbra_half_angle,
            1e-7,
        );
    }
}
