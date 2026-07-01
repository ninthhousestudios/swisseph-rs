use serde::Deserialize;
use swisseph::{CalcFlags, Ephemeris};

#[derive(Deserialize)]
struct SolWhereCase {
    tjd_ut: f64,
    nonut: bool,
    retval: i32,
    geopos: [f64; 10],
    attr: [f64; 11],
}

#[derive(Deserialize)]
struct GoldenData {
    sol_where: Vec<SolWhereCase>,
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
    }
}
