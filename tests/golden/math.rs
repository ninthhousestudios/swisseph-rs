use serde::Deserialize;
use swisseph::flags::SplitDegFlags;
use swisseph::math;

// ---------------------------------------------------------------------------
// Deserialization types
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct NormCase {
    input: f64,
    output: f64,
}

#[derive(Deserialize)]
struct DiffCase {
    p1: f64,
    p2: f64,
    output: f64,
}

#[derive(Deserialize)]
struct MidpCase {
    x1: f64,
    x0: f64,
    output: f64,
}

#[derive(Deserialize)]
struct CsnormCase {
    input: i32,
    output: i32,
}

#[derive(Deserialize)]
struct DifcsCase {
    p1: i32,
    p2: i32,
    output: i32,
}

#[derive(Deserialize)]
struct D2lCase {
    input: f64,
    output: i32,
}

#[derive(Deserialize)]
struct ChebyshevCase {
    t: f64,
    label: String,
    coef: Vec<f64>,
    value: f64,
    deriv: f64,
}

#[derive(Deserialize)]
struct CartpolCase {
    x: f64,
    y: f64,
    z: f64,
    lon: f64,
    lat: f64,
    dist: f64,
}

#[derive(Deserialize)]
struct PolcartCase {
    lon: f64,
    lat: f64,
    dist: f64,
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]
struct CoordTransCase {
    xi: Vec<f64>,
    eps: f64,
    xo: Vec<f64>,
}

#[derive(Deserialize)]
struct SpCase {
    xi: Vec<f64>,
    xo: Vec<f64>,
}

#[derive(Deserialize)]
struct SpCaseL {
    li: Vec<f64>,
    xo: Vec<f64>,
}

#[derive(Deserialize)]
struct SplitDegCase {
    ddeg: f64,
    flags: u32,
    deg: i32,
    min: i32,
    sec: i32,
    secfr: f64,
    sign: i32,
}

#[derive(Deserialize)]
struct MathGolden {
    degnorm: Vec<NormCase>,
    radnorm: Vec<NormCase>,
    difdeg2n: Vec<DiffCase>,
    difdegn: Vec<DiffCase>,
    difrad2n: Vec<DiffCase>,
    midp_deg: Vec<MidpCase>,
    midp_rad: Vec<MidpCase>,
    csnorm: Vec<CsnormCase>,
    difcsn: Vec<DifcsCase>,
    difcs2n: Vec<DifcsCase>,
    d2l: Vec<D2lCase>,
    chebyshev: Vec<ChebyshevCase>,
    cartpol: Vec<CartpolCase>,
    polcart: Vec<PolcartCase>,
    coortrf: Vec<CoordTransCase>,
    cartpol_sp: Vec<SpCase>,
    polcart_sp: Vec<SpCaseL>,
    cotrans: Vec<CoordTransCase>,
    cotrans_sp: Vec<CoordTransCase>,
    split_deg: Vec<SplitDegCase>,
}

fn load() -> MathGolden {
    let path = super::golden_data_path("math.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn golden_degnorm() {
    let data = load();
    for (i, c) in data.degnorm.iter().enumerate() {
        let actual = math::normalize_degrees(c.input);
        super::assert_f64_exact(&format!("degnorm[{i}]({:.20})", c.input), c.output, actual);
    }
}

#[test]
fn golden_radnorm() {
    let data = load();
    for (i, c) in data.radnorm.iter().enumerate() {
        let actual = math::normalize_radians(c.input);
        super::assert_f64_exact(&format!("radnorm[{i}]({:.20})", c.input), c.output, actual);
    }
}

#[test]
fn golden_difdeg2n() {
    let data = load();
    for (i, c) in data.difdeg2n.iter().enumerate() {
        let actual = math::diff_degrees(c.p1, c.p2);
        super::assert_f64_exact(
            &format!("difdeg2n[{i}]({:.20}, {:.20})", c.p1, c.p2),
            c.output,
            actual,
        );
    }
}

#[test]
fn golden_difdegn() {
    let data = load();
    for (i, c) in data.difdegn.iter().enumerate() {
        let actual = math::diff_degrees_norm(c.p1, c.p2);
        super::assert_f64_exact(
            &format!("difdegn[{i}]({:.20}, {:.20})", c.p1, c.p2),
            c.output,
            actual,
        );
    }
}

#[test]
fn golden_difrad2n() {
    let data = load();
    for (i, c) in data.difrad2n.iter().enumerate() {
        let actual = math::diff_radians(c.p1, c.p2);
        super::assert_f64_exact(
            &format!("difrad2n[{i}]({:.20}, {:.20})", c.p1, c.p2),
            c.output,
            actual,
        );
    }
}

#[test]
fn golden_midp_deg() {
    let data = load();
    for (i, c) in data.midp_deg.iter().enumerate() {
        let actual = math::midpoint_degrees(c.x1, c.x0);
        super::assert_f64_exact(
            &format!("midp_deg[{i}]({:.20}, {:.20})", c.x1, c.x0),
            c.output,
            actual,
        );
    }
}

#[test]
fn golden_midp_rad() {
    let data = load();
    for (i, c) in data.midp_rad.iter().enumerate() {
        let actual = math::midpoint_radians(c.x1, c.x0);
        super::assert_f64_exact(
            &format!("midp_rad[{i}]({:.20}, {:.20})", c.x1, c.x0),
            c.output,
            actual,
        );
    }
}

#[test]
fn golden_csnorm() {
    let data = load();
    for (i, c) in data.csnorm.iter().enumerate() {
        let actual = math::csnorm(c.input);
        assert_eq!(
            actual, c.output,
            "csnorm[{i}]({}) expected {}, got {actual}",
            c.input, c.output
        );
    }
}

#[test]
fn golden_difcsn() {
    let data = load();
    for (i, c) in data.difcsn.iter().enumerate() {
        let actual = math::difcsn(c.p1, c.p2);
        assert_eq!(
            actual, c.output,
            "difcsn[{i}]({}, {}) expected {}, got {actual}",
            c.p1, c.p2, c.output
        );
    }
}

#[test]
fn golden_difcs2n() {
    let data = load();
    for (i, c) in data.difcs2n.iter().enumerate() {
        let actual = math::difcs2n(c.p1, c.p2);
        assert_eq!(
            actual, c.output,
            "difcs2n[{i}]({}, {}) expected {}, got {actual}",
            c.p1, c.p2, c.output
        );
    }
}

#[test]
fn golden_d2l() {
    let data = load();
    for (i, c) in data.d2l.iter().enumerate() {
        let actual = math::d2l(c.input);
        assert_eq!(
            actual, c.output,
            "d2l[{i}]({}) expected {}, got {actual}",
            c.input, c.output
        );
    }
}

#[test]
fn golden_chebyshev_eval() {
    let data = load();
    for (i, c) in data.chebyshev.iter().enumerate() {
        let actual = math::chebyshev_eval(c.t, &c.coef);
        super::assert_f64_exact(
            &format!("echeb[{i}] {}(t={:.20})", c.label, c.t),
            c.value,
            actual,
        );
    }
}

#[test]
fn golden_chebyshev_deriv() {
    let data = load();
    for (i, c) in data.chebyshev.iter().enumerate() {
        let actual = math::chebyshev_deriv(c.t, &c.coef);
        super::assert_f64_exact(
            &format!("edcheb[{i}] {}(t={:.20})", c.label, c.t),
            c.deriv,
            actual,
        );
    }
}

#[test]
fn golden_cartpol() {
    let data = load();
    for (i, c) in data.cartpol.iter().enumerate() {
        let actual = math::cartesian_to_polar([c.x, c.y, c.z]);
        let label = format!("cartpol[{i}]({}, {}, {})", c.x, c.y, c.z);
        super::assert_f64_exact(&format!("{label} lon"), c.lon, actual[0]);
        super::assert_f64_exact(&format!("{label} lat"), c.lat, actual[1]);
        super::assert_f64_exact(&format!("{label} dist"), c.dist, actual[2]);
    }
}

#[test]
fn golden_polcart() {
    let data = load();
    for (i, c) in data.polcart.iter().enumerate() {
        let actual = math::polar_to_cartesian([c.lon, c.lat, c.dist]);
        let label = format!("polcart[{i}]({}, {}, {})", c.lon, c.lat, c.dist);
        super::assert_f64_exact(&format!("{label} x"), c.x, actual[0]);
        super::assert_f64_exact(&format!("{label} y"), c.y, actual[1]);
        super::assert_f64_exact(&format!("{label} z"), c.z, actual[2]);
    }
}

#[test]
fn golden_coortrf() {
    let data = load();
    for (i, c) in data.coortrf.iter().enumerate() {
        let actual = math::rotate_x([c.xi[0], c.xi[1], c.xi[2]], c.eps);
        let label = format!("coortrf[{i}]");
        super::assert_f64_exact(&format!("{label} [0]"), c.xo[0], actual[0]);
        super::assert_f64_exact(&format!("{label} [1]"), c.xo[1], actual[1]);
        super::assert_f64_exact(&format!("{label} [2]"), c.xo[2], actual[2]);
    }
}

#[test]
fn golden_cartpol_sp() {
    let data = load();
    for (i, c) in data.cartpol_sp.iter().enumerate() {
        let actual = math::cartesian_to_polar_with_speed([
            c.xi[0], c.xi[1], c.xi[2], c.xi[3], c.xi[4], c.xi[5],
        ]);
        let label = format!("cartpol_sp[{i}]");
        for (j, &a) in actual.iter().enumerate() {
            super::assert_f64_exact(&format!("{label} [{j}]"), c.xo[j], a);
        }
    }
}

#[test]
fn golden_polcart_sp() {
    let data = load();
    for (i, c) in data.polcart_sp.iter().enumerate() {
        let actual = math::polar_to_cartesian_with_speed([
            c.li[0], c.li[1], c.li[2], c.li[3], c.li[4], c.li[5],
        ]);
        let label = format!("polcart_sp[{i}]");
        for (j, &a) in actual.iter().enumerate() {
            super::assert_f64_exact(&format!("{label} [{j}]"), c.xo[j], a);
        }
    }
}

#[test]
fn golden_cotrans() {
    let data = load();
    for (i, c) in data.cotrans.iter().enumerate() {
        let actual = math::cotrans([c.xi[0], c.xi[1], c.xi[2]], c.eps);
        let label = format!("cotrans[{i}]");
        super::assert_f64_exact(&format!("{label} [0]"), c.xo[0], actual[0]);
        super::assert_f64_exact(&format!("{label} [1]"), c.xo[1], actual[1]);
        super::assert_f64_exact(&format!("{label} [2]"), c.xo[2], actual[2]);
    }
}

#[test]
fn golden_cotrans_sp() {
    let data = load();
    for (i, c) in data.cotrans_sp.iter().enumerate() {
        let actual = math::cotrans_with_speed(
            [c.xi[0], c.xi[1], c.xi[2], c.xi[3], c.xi[4], c.xi[5]],
            c.eps,
        );
        let label = format!("cotrans_sp[{i}]");
        for (j, &a) in actual.iter().enumerate() {
            super::assert_f64_exact(&format!("{label} [{j}]"), c.xo[j], a);
        }
    }
}

#[test]
fn golden_split_deg() {
    let data = load();
    for (i, c) in data.split_deg.iter().enumerate() {
        let flags = SplitDegFlags::from_bits_truncate(c.flags);
        let actual = math::split_degrees(c.ddeg, flags);
        let label = format!("split_deg[{i}](ddeg={}, flags={:#x})", c.ddeg, c.flags);
        assert_eq!(actual.degrees, c.deg, "{label}: degrees");
        assert_eq!(actual.minutes, c.min, "{label}: minutes");
        assert_eq!(actual.seconds, c.sec, "{label}: seconds");
        super::assert_f64_exact(&format!("{label}: secfr"), c.secfr, actual.second_fraction);
        assert_eq!(actual.sign, c.sign, "{label}: sign");
    }
}
