use serde::Deserialize;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct PctrCase {
    ipl: i32,
    iplctr: i32,
    tjd: f64,
    iflag: i32,
    retflag: i32,
    xx: [f64; 6],
    ok: i32,
}

#[derive(Deserialize)]
struct GoldenData {
    pctr: Vec<PctrCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("pctr.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn make_eph(source: EphemerisSource) -> Option<Ephemeris> {
    match source {
        EphemerisSource::Moshier => {
            Some(Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new moshier"))
        }
        EphemerisSource::Swiss => {
            let ephe_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .join("swisseph")
                .join("ephe");
            if !ephe_path.exists() {
                return None;
            }
            let mut config = EphemerisConfig::default();
            config.ephemeris_source = EphemerisSource::Swiss;
            config.ephe_path = Some(ephe_path);
            Some(Ephemeris::new(config).expect("Ephemeris::new sweph"))
        }
        _ => None,
    }
}

// Positions 5e-8, speeds 1e-7. The position tolerance is wider than the
// standard pipeline's 1e-9 because pctr's deflection (§5) uses earth_helio
// while C uses earth_bary (pedp->x). The ~0.005 AU difference propagates
// into a deflection correction error up to ~1.4e-8 for distant body pairs
// (Saturn-Jupiter worst case). Same architectural cause as the documented
// 1e-7 speed tolerance (CLAUDE.md <stateless_tolerance> §1).
const POS_EPS: f64 = 5e-8;
const SPEED_EPS: f64 = 1e-7;

#[test]
fn pctr() {
    let data = load();
    let eph_moshier = make_eph(EphemerisSource::Moshier).unwrap();
    let eph_sweph = make_eph(EphemerisSource::Swiss);

    let mut ok_count = 0;
    let mut err_count = 0;

    for (i, c) in data.pctr.iter().enumerate() {
        let body =
            Body::try_from(c.ipl).unwrap_or_else(|e| panic!("case {i}: bad ipl {}: {e}", c.ipl));
        let center = Body::try_from(c.iplctr)
            .unwrap_or_else(|e| panic!("case {i}: bad iplctr {}: {e}", c.iplctr));

        let raw_flags = CalcFlags::from_bits_truncate(c.iflag as u32);
        let is_sweph = raw_flags.contains(CalcFlags::SWIEPH);

        let eph = if is_sweph {
            match &eph_sweph {
                Some(e) => e,
                None => continue,
            }
        } else {
            &eph_moshier
        };

        let result = eph.calc_pctr(c.tjd, body, center, raw_flags);

        if c.ok == 0 {
            assert!(
                result.is_err(),
                "case {i}: ipl={}, iplctr={}, iflag={:#x}: expected Err but got Ok",
                c.ipl,
                c.iplctr,
                c.iflag
            );
            err_count += 1;
        } else {
            let r = result.unwrap_or_else(|e| {
                panic!(
                    "case {i}: ipl={}, iplctr={}, tjd={}, iflag={:#x}: error: {e}",
                    c.ipl, c.iplctr, c.tjd, c.iflag
                )
            });

            for j in 0..3 {
                let diff = (r.data[j] - c.xx[j]).abs();
                assert!(
                    diff <= POS_EPS,
                    "case {i}: ipl={}, iplctr={}, tjd={}, iflag={:#x}: xx[{j}] diff {diff:.2e} > {POS_EPS:.0e}",
                    c.ipl,
                    c.iplctr,
                    c.tjd,
                    c.iflag
                );
            }
            for j in 3..6 {
                let diff = (r.data[j] - c.xx[j]).abs();
                assert!(
                    diff <= SPEED_EPS,
                    "case {i}: ipl={}, iplctr={}, tjd={}, iflag={:#x}: xx[{j}] diff {diff:.2e} > {SPEED_EPS:.0e}",
                    c.ipl,
                    c.iplctr,
                    c.tjd,
                    c.iflag
                );
            }

            let expected_flags = CalcFlags::from_bits_truncate(c.retflag as u32);
            assert_eq!(
                r.flags_used, expected_flags,
                "case {i}: ipl={}, iplctr={}, tjd={}: flags_used mismatch",
                c.ipl, c.iplctr, c.tjd
            );

            ok_count += 1;
        }
    }

    println!(
        "All {} pctr cases passed ({ok_count} ok, {err_count} err).",
        data.pctr.len()
    );
}
