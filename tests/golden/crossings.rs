use serde::Deserialize;
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig};

#[derive(Deserialize)]
struct LonCrossCase {
    x2cross: f64,
    jd_start: f64,
    variant: String,
    jd_result: f64,
    ok: i32,
}

#[derive(Deserialize)]
struct NodeCrossCase {
    jd_start: f64,
    variant: String,
    jd_result: f64,
    xlon: f64,
    xlat: f64,
    ok: i32,
}

#[derive(Deserialize)]
struct HelioCrossCase {
    ipl: i32,
    x2cross: f64,
    jd_start: f64,
    dir: i32,
    variant: String,
    jd_result: f64,
    rc: i32,
}

#[derive(Deserialize)]
struct GoldenData {
    solcross: Vec<LonCrossCase>,
    mooncross: Vec<LonCrossCase>,
    mooncross_node: Vec<NodeCrossCase>,
    helio_cross: Vec<HelioCrossCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("crossings.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn make_eph() -> Ephemeris {
    Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new")
}

const TIME_EPS: f64 = 1e-6;
const LON_EPS: f64 = 1e-7;
const LAT_EPS: f64 = 5e-9;

#[test]
fn crossings() {
    let data = load();
    let eph = make_eph();
    let flags = CalcFlags::MOSEPH;

    println!("--- solcross ({} cases) ---", data.solcross.len());
    for (i, c) in data.solcross.iter().enumerate() {
        assert_eq!(c.ok, 1, "solcross case {i}: C reported failure");
        let result = match c.variant.as_str() {
            "et" => eph.solcross(c.x2cross, c.jd_start, flags),
            "ut" => eph.solcross_ut(c.x2cross, c.jd_start, flags),
            v => panic!("unknown variant: {v}"),
        };
        let jd = result.unwrap_or_else(|e| panic!("solcross case {i} ({}): error: {e}", c.variant));
        let diff = (jd - c.jd_result).abs();
        assert!(
            diff < TIME_EPS,
            "solcross case {i} ({}): x2cross={}, jd_start={}: jd diff {diff:.2e} (got {jd}, expected {})",
            c.variant,
            c.x2cross,
            c.jd_start,
            c.jd_result
        );
    }

    println!("--- mooncross ({} cases) ---", data.mooncross.len());
    for (i, c) in data.mooncross.iter().enumerate() {
        assert_eq!(c.ok, 1, "mooncross case {i}: C reported failure");
        let result = match c.variant.as_str() {
            "et" => eph.mooncross(c.x2cross, c.jd_start, flags),
            "ut" => eph.mooncross_ut(c.x2cross, c.jd_start, flags),
            v => panic!("unknown variant: {v}"),
        };
        let jd =
            result.unwrap_or_else(|e| panic!("mooncross case {i} ({}): error: {e}", c.variant));
        let diff = (jd - c.jd_result).abs();
        assert!(
            diff < TIME_EPS,
            "mooncross case {i} ({}): x2cross={}, jd_start={}: jd diff {diff:.2e} (got {jd}, expected {})",
            c.variant,
            c.x2cross,
            c.jd_start,
            c.jd_result
        );
    }

    println!(
        "--- mooncross_node ({} cases) ---",
        data.mooncross_node.len()
    );
    for (i, c) in data.mooncross_node.iter().enumerate() {
        assert_eq!(c.ok, 1, "mooncross_node case {i}: C reported failure");
        let result = match c.variant.as_str() {
            "et" => eph.mooncross_node(c.jd_start, flags),
            "ut" => eph.mooncross_node_ut(c.jd_start, flags),
            v => panic!("unknown variant: {v}"),
        };
        let mc = result
            .unwrap_or_else(|e| panic!("mooncross_node case {i} ({}): error: {e}", c.variant));
        let jd_diff = (mc.jd - c.jd_result).abs();
        let lon_diff = (mc.longitude - c.xlon).abs();
        let lat_diff = (mc.latitude - c.xlat).abs();
        assert!(
            jd_diff < TIME_EPS,
            "mooncross_node case {i} ({}): jd diff {jd_diff:.2e} (got {}, expected {})",
            c.variant,
            mc.jd,
            c.jd_result
        );
        assert!(
            lon_diff < LON_EPS,
            "mooncross_node case {i} ({}): lon diff {lon_diff:.2e}",
            c.variant
        );
        assert!(
            mc.latitude.abs() < LAT_EPS,
            "mooncross_node case {i} ({}): |lat| = {:.2e} (should be ~0)",
            c.variant,
            mc.latitude.abs()
        );
        assert!(
            lat_diff < LAT_EPS,
            "mooncross_node case {i} ({}): lat diff {lat_diff:.2e}",
            c.variant
        );
    }

    println!("--- helio_cross ({} cases) ---", data.helio_cross.len());
    for (i, c) in data.helio_cross.iter().enumerate() {
        assert_eq!(c.rc, 0, "helio_cross case {i}: C returned error");
        let body = Body::try_from(c.ipl)
            .unwrap_or_else(|e| panic!("helio_cross case {i}: bad ipl {}: {e}", c.ipl));
        let result = match c.variant.as_str() {
            "et" => eph.helio_cross(body, c.x2cross, c.jd_start, flags, c.dir),
            "ut" => eph.helio_cross_ut(body, c.x2cross, c.jd_start, flags, c.dir),
            v => panic!("unknown variant: {v}"),
        };
        let jd = result.unwrap_or_else(|e| {
            panic!(
                "helio_cross case {i} ({}, ipl={}, dir={}): error: {e}",
                c.variant, c.ipl, c.dir
            )
        });
        let diff = (jd - c.jd_result).abs();
        assert!(
            diff < TIME_EPS,
            "helio_cross case {i} ({}, ipl={}, dir={}, x2cross={}): jd diff {diff:.2e} (got {jd}, expected {})",
            c.variant,
            c.ipl,
            c.dir,
            c.x2cross,
            c.jd_result
        );
    }

    println!(
        "All {} crossings cases passed.",
        data.solcross.len()
            + data.mooncross.len()
            + data.mooncross_node.len()
            + data.helio_cross.len()
    );
}
