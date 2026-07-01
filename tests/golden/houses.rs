use serde::Deserialize;
use swisseph::houses::houses_armc;
use swisseph::types::HouseSystem;
use swisseph::{CalcFlags, Ephemeris, EphemerisConfig};

#[derive(Deserialize)]
struct AnglesSpecialCase {
    armc: f64,
    geolat: f64,
    eps: f64,
    ascmc: [f64; 8],
    ascmc_speed: [f64; 8],
}

#[derive(Deserialize)]
struct EqualFamilyCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct QuadArithCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct GreatCircleCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct IterativeCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct Gauquelin36Case {
    armc: f64,
    geolat: f64,
    eps: f64,
    // serde's array impl tops out at 32 elements; 36 cusps need a Vec.
    cusps: Vec<f64>,
    cusp_speed: Vec<f64>,
}

#[derive(Deserialize)]
struct ClosedFormMiscCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct SunshineCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    sundec: f64,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
}

#[derive(Deserialize)]
struct UtWrapperCase {
    tjd_ut: f64,
    geolat: f64,
    geolon: f64,
    hsys: String,
    nonut: bool,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
    ascmc: [f64; 8],
    ascmc_speed: [f64; 8],
}

#[derive(Deserialize)]
struct SiderealTradCase {
    tjd_ut: f64,
    geolat: f64,
    geolon: f64,
    hsys: String,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
    ascmc: [f64; 8],
    ascmc_speed: [f64; 8],
}

#[derive(Deserialize)]
struct SiderealGeomCase {
    tjd_ut: f64,
    geolat: f64,
    geolon: f64,
    hsys: String,
    sid_mode: i32,
    cusps: [f64; 12],
    cusp_speed: [f64; 12],
    ascmc: [f64; 8],
    ascmc_speed: [f64; 8],
}

#[derive(Deserialize)]
struct HousePosCase {
    hsys: String,
    armc: f64,
    geolat: f64,
    eps: f64,
    xpin: [f64; 2],
    sundec: f64,
    hpos: f64,
    err: bool,
}

#[derive(Deserialize)]
struct GauquelinSectorCase {
    tjd_ut: f64,
    ipl: i32,
    imeth: i32,
    geolon: f64,
    geolat: f64,
    dgsect: f64,
}

#[derive(Deserialize)]
struct GoldenData {
    angles_special: Vec<AnglesSpecialCase>,
    equal_family: Vec<EqualFamilyCase>,
    quad_arith: Vec<QuadArithCase>,
    great_circle: Vec<GreatCircleCase>,
    iterative: Vec<IterativeCase>,
    gauquelin36: Vec<Gauquelin36Case>,
    closed_form_misc: Vec<ClosedFormMiscCase>,
    sunshine: Vec<SunshineCase>,
    ut_wrapper: Vec<UtWrapperCase>,
    sidereal_trad: Vec<SiderealTradCase>,
    sidereal_geom: Vec<SiderealGeomCase>,
    house_pos: Vec<HousePosCase>,
    gauquelin_sector: Vec<GauquelinSectorCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("houses.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn parse_hsys(s: &str) -> HouseSystem {
    HouseSystem::try_from(s.as_bytes()[0])
        .unwrap_or_else(|e| panic!("Unknown house system {s}: {e}"))
}

#[test]
fn angles_special() {
    let data = load();
    assert_eq!(
        data.angles_special.len(),
        30,
        "expected 30 golden cases (6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.angles_special.iter().enumerate() {
        // The special points are system-independent; 'P' was used by the C generator,
        // but only Equal-family systems are ported so far, so use Equal here.
        let result = houses_armc(c.armc, c.geolat, c.eps, HouseSystem::Equal, None)
            .unwrap_or_else(|e| panic!("case {i}: houses_armc failed: {e}"));

        let actual = result.ascmc.as_array();
        let actual_speed = result.ascmc_speeds.as_array();
        let label_base = format!(
            "case {i} (armc={:.6} geolat={:.6} eps={:.6})",
            c.armc, c.geolat, c.eps
        );
        for j in 0..8 {
            super::assert_f64_exact(&format!("{label_base} ascmc[{j}]"), c.ascmc[j], actual[j]);
            super::assert_f64_exact(
                &format!("{label_base} ascmc_speed[{j}]"),
                c.ascmc_speed[j],
                actual_speed[j],
            );
        }
    }
}

#[test]
fn quad_arith() {
    let data = load();
    assert_eq!(
        data.quad_arith.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.quad_arith.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        // Porphyry cusp speeds are analytical (linear quadrant-rate interpolation);
        // S/X/M/F use the driver-level finite-difference path, which is not
        // bitwise-exact against C's central difference.
        let speed_eps = if c.hsys == "O" { 1e-9 } else { 1e-7 };

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                speed_eps,
            );
        }
    }
}

#[test]
fn great_circle() {
    let data = load();
    assert_eq!(
        data.great_circle.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.great_circle.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-9,
            );
        }
    }
}

#[test]
fn iterative() {
    let data = load();
    assert_eq!(
        data.iterative.len(),
        84,
        "expected 84 golden cases (2 systems x 6 armc x 7 geolat incl. polar x 1 eps)"
    );
    for (i, c) in data.iterative.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-7,
            );
        }
    }
}

#[test]
fn gauquelin36() {
    let data = load();
    assert_eq!(
        data.gauquelin36.len(),
        42,
        "expected 42 golden cases (6 armc x 7 geolat incl. polar x 1 eps)"
    );
    for (i, c) in data.gauquelin36.iter().enumerate() {
        let result = houses_armc(c.armc, c.geolat, c.eps, HouseSystem::Gauquelin, None)
            .unwrap_or_else(|e| panic!("case {i}: houses_armc failed: {e}"));

        let label_base = format!(
            "case {i} (G armc={:.6} geolat={:.6} eps={:.6})",
            c.armc, c.geolat, c.eps
        );
        for h in 1..=36usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-7,
            );
        }
    }
}

#[test]
fn equal_family() {
    let data = load();
    assert_eq!(
        data.equal_family.len(),
        150,
        "expected 150 golden cases (5 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.equal_family.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_exact(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
            );
            super::assert_f64_exact(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
            );
        }
    }
}

#[test]
fn closed_form_misc() {
    let data = load();
    assert_eq!(
        data.closed_form_misc.len(),
        120,
        "expected 120 golden cases (4 systems x 6 armc x 5 geolat x 1 eps)"
    );
    for (i, c) in data.closed_form_misc.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, None)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6})",
            c.hsys, c.armc, c.geolat, c.eps
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-9,
            );
            // U's cusp speeds are stale pre-switch values (not analytical or finite-diff,
            // see c-ref-houses.md §4.2e) — assert them exactly as C produces, including zeros.
            let speed_eps = if c.hsys == "U" { 0.0 } else { 1e-7 };
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                speed_eps,
            );
        }
    }
}

#[test]
fn sunshine() {
    let data = load();
    assert_eq!(
        data.sunshine.len(),
        76,
        "expected 76 golden cases (60: 2 systems x 6 armc x 5 geolat, 1 sundec per case; \
         + 16: 2 systems x 2 armc x geolat in {{70,-70}} x sundec in {{23,-23}}, all four of \
         which satisfy |tand(geolat)*tand(sundec)|>=1 -- circumpolar-Sun combinations that \
         exercise Makransky's sunshine_init ERR -> Porphyry fallback (do_interpol stays false \
         on that path; Treindl never short-circuits on it and is included at the same \
         combinations for contrast)"
    );
    for (i, c) in data.sunshine.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = houses_armc(c.armc, c.geolat, c.eps, hsys, Some(c.sundec))
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_armc failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6} sundec={:.6})",
            c.hsys, c.armc, c.geolat, c.eps, c.sundec
        );
        // The polar-battery cases (|geolat|=70, outside the standard battery's max of 64) are
        // the only ones that can trigger Makransky's circumpolar ERR; at those combinations 'i'
        // falls back to fill_porphyry (closed-form, already bitwise-exact elsewhere — see
        // quad_arith's 'O' cases) rather than the quadrant case-split, so it gets the tighter
        // tolerance too.
        let fallback = hsys == HouseSystem::SunshineAlt && c.geolat.abs() > 64.0;
        // Sunshine is closed-form per house (Treindl directly, Makransky via a quadrant case
        // split); Makransky's case split may need the looser tolerance.
        let cusp_eps = if hsys == HouseSystem::SunshineAlt && !fallback {
            1e-8
        } else {
            1e-9
        };
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                cusp_eps,
            );
            // I/i use the driver-level finite-difference cusp speed path (do_interpol) --
            // except the Porphyry-fallback cases, which use fill_porphyry's analytical speeds
            // (do_interpol is never set on that path) and so get the tighter tolerance.
            let speed_eps = if fallback { 1e-9 } else { 1e-7 };
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                speed_eps,
            );
        }
    }
}

#[test]
fn sunshine_requires_sundec() {
    let err = houses_armc(0.0, 51.5, 23.4392911, HouseSystem::Sunshine, None)
        .expect_err("Sunshine without sundec must error");
    assert!(matches!(err, swisseph::error::Error::CError(_)));
}

#[test]
fn ut_wrapper() {
    let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();
    let data = load();
    assert_eq!(
        data.ut_wrapper.len(),
        42,
        "expected 42 golden cases (36: 6 triples x 6 systems + 6: 1 triple x 6 systems NONUT)"
    );
    for (i, c) in data.ut_wrapper.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let flags = if c.nonut {
            CalcFlags::NONUT
        } else {
            CalcFlags::empty()
        };
        let result = eph
            .houses_ex2(c.tjd_ut, flags, c.geolat, c.geolon, hsys)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_ex2 failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} tjd_ut={:.6} geolat={:.6} geolon={:.6} nonut={})",
            c.hsys, c.tjd_ut, c.geolat, c.geolon, c.nonut
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-7,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-6,
            );
        }
        let actual_ascmc = result.ascmc.as_array();
        let actual_ascmc_speed = result.ascmc_speeds.as_array();
        for j in 0..8 {
            super::assert_f64_eps(
                &format!("{label_base} ascmc[{j}]"),
                c.ascmc[j],
                actual_ascmc[j],
                1e-7,
            );
            super::assert_f64_eps(
                &format!("{label_base} ascmc_speed[{j}]"),
                c.ascmc_speed[j],
                actual_ascmc_speed[j],
                1e-6,
            );
        }
    }
}

#[test]
fn sidereal_trad() {
    let mut config = EphemerisConfig::default();
    config.set_sidereal_mode(1 /* Lahiri */, 0.0, 0.0);
    let eph = Ephemeris::new(config).unwrap();
    let data = load();
    assert_eq!(
        data.sidereal_trad.len(),
        9,
        "expected 9 golden cases (3 systems x 3 triples)"
    );
    for (i, c) in data.sidereal_trad.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let result = eph
            .houses_ex2(c.tjd_ut, CalcFlags::SIDEREAL, c.geolat, c.geolon, hsys)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_ex2 failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} tjd_ut={:.6} geolat={:.6} geolon={:.6})",
            c.hsys, c.tjd_ut, c.geolat, c.geolon
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-7,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-6,
            );
        }
        let actual_ascmc = result.ascmc.as_array();
        let actual_ascmc_speed = result.ascmc_speeds.as_array();
        for j in 0..8 {
            super::assert_f64_eps(
                &format!("{label_base} ascmc[{j}]"),
                c.ascmc[j],
                actual_ascmc[j],
                1e-7,
            );
            super::assert_f64_eps(
                &format!("{label_base} ascmc_speed[{j}]"),
                c.ascmc_speed[j],
                actual_ascmc_speed[j],
                1e-6,
            );
        }
    }
}

#[test]
fn sidereal_geom() {
    let data = load();
    assert_eq!(
        data.sidereal_geom.len(),
        18,
        "expected 18 golden cases (2 sid_modes x 3 systems x 3 triples)"
    );
    for (i, c) in data.sidereal_geom.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let mut config = EphemerisConfig::default();
        config.set_sidereal_mode(c.sid_mode, 0.0, 0.0);
        let eph = Ephemeris::new(config).unwrap();
        let result = eph
            .houses_ex2(c.tjd_ut, CalcFlags::SIDEREAL, c.geolat, c.geolon, hsys)
            .unwrap_or_else(|e| panic!("case {i} ({}): houses_ex2 failed: {e}", c.hsys));

        let label_base = format!(
            "case {i} ({} tjd_ut={:.6} geolat={:.6} geolon={:.6} sid_mode={})",
            c.hsys, c.tjd_ut, c.geolat, c.geolon, c.sid_mode
        );
        for h in 1..=12usize {
            super::assert_f64_eps(
                &format!("{label_base} cusp[{h}]"),
                c.cusps[h - 1],
                result.cusps[h],
                1e-7,
            );
            super::assert_f64_eps(
                &format!("{label_base} cusp_speed[{h}]"),
                c.cusp_speed[h - 1],
                result.cusp_speeds[h],
                1e-6,
            );
        }
        let actual_ascmc = result.ascmc.as_array();
        let actual_ascmc_speed = result.ascmc_speeds.as_array();
        for j in 0..8 {
            super::assert_f64_eps(
                &format!("{label_base} ascmc[{j}]"),
                c.ascmc[j],
                actual_ascmc[j],
                1e-7,
            );
            super::assert_f64_eps(
                &format!("{label_base} ascmc_speed[{j}]"),
                c.ascmc_speed[j],
                actual_ascmc_speed[j],
                1e-6,
            );
        }
    }
}

#[test]
fn house_pos() {
    let data = load();
    assert_eq!(
        data.house_pos.len(),
        150,
        "expected 150 golden cases (25 systems x 2 armc/geolat/eps triples x 3 xpin)"
    );
    for (i, c) in data.house_pos.iter().enumerate() {
        let hsys = parse_hsys(&c.hsys);
        let label = format!(
            "case {i} ({} armc={:.6} geolat={:.6} eps={:.6} xpin=[{:.6},{:.6}])",
            c.hsys, c.armc, c.geolat, c.eps, c.xpin[0], c.xpin[1]
        );
        let result =
            swisseph::houses::house_pos(c.armc, c.geolat, c.eps, hsys, c.xpin, Some(c.sundec));
        if c.err {
            result
                .err()
                .unwrap_or_else(|| panic!("{label}: expected Err, got Ok"));
        } else {
            let hpos = result.unwrap_or_else(|e| panic!("{label}: house_pos failed: {e}"));
            super::assert_f64_eps(&format!("{label} hpos"), c.hpos, hpos, 1e-7);
        }
    }
}

#[test]
fn gauquelin_sector() {
    let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();
    let data = load();
    assert_eq!(
        data.gauquelin_sector.len(),
        36,
        "expected 36 golden cases (6 tjd_ut/geolat/geolon triples x 3 bodies x 2 imeth)"
    );
    for (i, c) in data.gauquelin_sector.iter().enumerate() {
        let body = swisseph::types::Body::try_from(c.ipl)
            .unwrap_or_else(|e| panic!("case {i}: unknown body {}: {e}", c.ipl));
        let label = format!(
            "case {i} (ipl={} imeth={} tjd_ut={:.6} geolon={:.6} geolat={:.6})",
            c.ipl, c.imeth, c.tjd_ut, c.geolon, c.geolat
        );
        let dgsect = eph
            .gauquelin_sector_geometric(
                c.tjd_ut,
                body,
                c.imeth,
                CalcFlags::empty(),
                c.geolon,
                c.geolat,
            )
            .unwrap_or_else(|e| panic!("{label}: gauquelin_sector_geometric failed: {e}"));
        super::assert_f64_eps(&format!("{label} dgsect"), c.dgsect, dgsect, 1e-6);
    }
}
