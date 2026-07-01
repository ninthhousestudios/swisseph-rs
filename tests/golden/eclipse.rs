use serde::Deserialize;
use swisseph::{CalcFlags, EclipseFlags, Ephemeris};

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
struct GoldenData {
    sol_where: Vec<SolWhereCase>,
    sol_how: Vec<SolHowCase>,
    sol_when_glob: Vec<SolWhenGlobCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("eclipse.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
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
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.umbra_diameter_fundamental_km"),
            c.dcore[3],
            result.umbra_diameter_fundamental_km,
            1e-7,
        );
        super::assert_f64_eps(
            &format!("{label}.penumbra_diameter_fundamental_km"),
            c.dcore[4],
            result.penumbra_diameter_fundamental_km,
            1e-7,
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
