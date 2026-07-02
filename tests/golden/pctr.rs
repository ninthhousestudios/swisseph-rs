use serde::Deserialize;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, EphemerisSource};

#[derive(Deserialize)]
struct PctrCase {
    ipl: i32,
    iplctr: i32,
    tjd: f64,
    iflag: i32,
    retflag: i32,
    /// Sidereal mode index (with SE_SIDBIT bits) to feed set_sidereal_mode, or
    /// -1 for tropical. Sidereal cases are always SWIEPH.
    sid_mode: i32,
    sid_t0: f64,
    sid_ayan: f64,
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

fn ephe_path() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("swisseph")
        .join("ephe")
}

fn eph_moshier() -> Ephemeris {
    Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new moshier")
}

/// Swiss ephemeris, optionally with a sidereal mode set. Built unconditionally
/// (like tests/golden/calc_sweph.rs) — a missing ../swisseph/ephe fixture is a
/// hard failure, not a silent skip, so the numeric cases can never pass vacuously.
fn eph_sweph(sidereal: Option<(i32, f64, f64)>) -> Ephemeris {
    let mut config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path()),
        ..EphemerisConfig::default()
    };
    if let Some((sid_mode, t0, ayan)) = sidereal {
        config.set_sidereal_mode(sid_mode, t0, ayan);
    }
    Ephemeris::new(config).expect("Ephemeris::new sweph")
}

// Positions and speeds both to 3e-8 (worst observed: pos 1.44e-8, speed
// 1.21e-8). See swisseph-rs/99. The residual is NOT deflection: it is present
// under NOGDEFL and NOABERR alike and vanishes entirely under TRUEPOS, so it
// originates in the §3d retarded-time (t = tjd − dt) re-evaluation of the ipl
// body — a stateless-vs-stateful precision difference between C's evaluation at
// the retarded epoch and a fresh recompute. It shows up as a latitude offset
// (longitude stays ~4e-10), largest for the most widely-separated pair
// (Saturn↔Jupiter, ~10 AU): ≈0.05 mas — astronomically negligible.
// (The §5 deflection geometry itself matches C's swi_deflect_light exactly:
//  e = earth−sun, q = xx + earth−sun — verified during the 90 review.)
const POS_EPS: f64 = 3e-8;
const SPEED_EPS: f64 = 3e-8;

#[test]
fn pctr() {
    let data = load();
    let eph_m = eph_moshier();
    let eph_s = eph_sweph(None);

    let expected_ok = data.pctr.iter().filter(|c| c.ok == 1).count();
    assert!(expected_ok > 0, "golden data has no success cases");

    let mut ok_count = 0;
    let mut err_count = 0;
    let mut max_pos = 0.0_f64;
    let mut max_speed = 0.0_f64;

    for (i, c) in data.pctr.iter().enumerate() {
        let body =
            Body::try_from(c.ipl).unwrap_or_else(|e| panic!("case {i}: bad ipl {}: {e}", c.ipl));
        let center = Body::try_from(c.iplctr)
            .unwrap_or_else(|e| panic!("case {i}: bad iplctr {}: {e}", c.iplctr));
        let raw_flags = CalcFlags::from_bits_truncate(c.iflag as u32);

        // Sidereal cases need a per-case sidereal-mode config; tropical SWIEPH
        // and Moshier reuse the shared instances.
        let sid_eph;
        let eph: &Ephemeris = if c.sid_mode >= 0 {
            sid_eph = eph_sweph(Some((c.sid_mode, c.sid_t0, c.sid_ayan)));
            &sid_eph
        } else if raw_flags.contains(CalcFlags::SWIEPH) {
            &eph_s
        } else {
            &eph_m
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
            continue;
        }

        let r = result.unwrap_or_else(|e| {
            panic!(
                "case {i}: ipl={}, iplctr={}, tjd={}, iflag={:#x} sid={}: error: {e}",
                c.ipl, c.iplctr, c.tjd, c.iflag, c.sid_mode
            )
        });

        for j in 0..3 {
            let diff = (r.data[j] - c.xx[j]).abs();
            max_pos = max_pos.max(diff);
            assert!(
                diff <= POS_EPS,
                "case {i}: ipl={}, iplctr={}, tjd={}, iflag={:#x} sid={}: xx[{j}] diff {diff:.2e} > {POS_EPS:.0e}",
                c.ipl,
                c.iplctr,
                c.tjd,
                c.iflag,
                c.sid_mode
            );
        }
        for j in 3..6 {
            let diff = (r.data[j] - c.xx[j]).abs();
            max_speed = max_speed.max(diff);
            assert!(
                diff <= SPEED_EPS,
                "case {i}: ipl={}, iplctr={}, tjd={}, iflag={:#x} sid={}: xx[{j}] diff {diff:.2e} > {SPEED_EPS:.0e}",
                c.ipl,
                c.iplctr,
                c.tjd,
                c.iflag,
                c.sid_mode
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

    assert_eq!(
        ok_count, expected_ok,
        "expected {expected_ok} numeric successes, got {ok_count}"
    );

    println!(
        "All {} pctr cases passed ({ok_count} ok, {err_count} err). max pos diff {max_pos:.2e}, max speed diff {max_speed:.2e}.",
        data.pctr.len()
    );
}
