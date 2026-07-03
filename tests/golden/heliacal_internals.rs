use serde::Deserialize;
use swisseph::flags::{CalcFlags, HeliacalFlags};
use swisseph::heliacal;
use swisseph::{Ephemeris, EphemerisConfig};

#[derive(Deserialize)]
struct ExtinctionCase {
    #[serde(rename = "AltO")]
    alt_o: f64,
    #[serde(rename = "AltS")]
    alt_s: f64,
    sunra: f64,
    #[serde(rename = "Lat")]
    lat: f64,
    #[serde(rename = "HeightEye")]
    height_eye: f64,
    datm: [f64; 4],
    #[serde(rename = "Deltam")]
    deltam: f64,
    kt: f64,
    #[serde(rename = "kR")]
    k_r: f64,
    #[serde(rename = "kOZ")]
    k_oz: f64,
    #[serde(rename = "kW")]
    k_w: f64,
    ka: f64,
}

#[derive(Deserialize)]
struct AirmassCase {
    #[serde(rename = "AppAltO")]
    app_alt_o: f64,
    #[serde(rename = "Press")]
    press: f64,
    #[serde(rename = "Airmass")]
    airmass: f64,
    #[serde(rename = "Xext_rayleigh")]
    xext_rayleigh: f64,
    #[serde(rename = "Xext_water")]
    xext_water: f64,
    #[serde(rename = "Xext_aerosol")]
    xext_aerosol: f64,
    #[serde(rename = "Xlay_ozone")]
    xlay_ozone: f64,
}

#[derive(Deserialize)]
struct AppAltCase {
    alt: f64,
    #[serde(rename = "TempE")]
    temp_e: f64,
    #[serde(rename = "PresE")]
    pres_e: f64,
    #[serde(rename = "AppAltfromTopoAlt")]
    app_alt_from_topo_alt: f64,
    #[serde(rename = "TopoAltfromAppAlt")]
    topo_alt_from_app_alt: f64,
}

#[derive(Deserialize)]
struct OpticCase {
    #[serde(rename = "B")]
    b: f64,
    config: String,
    dobs: [f64; 6],
    helflag: u32,
    #[serde(rename = "kX")]
    k_x: f64,
    #[serde(rename = "CVA")]
    cva: f64,
    #[serde(rename = "PupilDia")]
    pupil_dia: f64,
    #[serde(rename = "OpticFactor_intensity")]
    optic_factor_intensity: f64,
    #[serde(rename = "OpticFactor_background")]
    optic_factor_background: f64,
}

#[derive(Deserialize)]
struct BrightnessCase {
    #[serde(rename = "AltO")]
    alt_o: f64,
    #[serde(rename = "AziO")]
    azi_o: f64,
    #[serde(rename = "AltM")]
    alt_m: f64,
    #[serde(rename = "AziM")]
    azi_m: f64,
    #[serde(rename = "AltS")]
    alt_s: f64,
    #[serde(rename = "AziS")]
    azi_s: f64,
    sunra: f64,
    #[serde(rename = "Lat")]
    lat: f64,
    #[serde(rename = "HeightEye")]
    height_eye: f64,
    datm: [f64; 4],
    #[serde(rename = "JDNDaysUT")]
    jdn_days_ut: f64,
    #[serde(rename = "Bn")]
    bn: f64,
    #[serde(rename = "Bm")]
    bm: f64,
    #[serde(rename = "Btwi")]
    btwi: f64,
    #[serde(rename = "Bday")]
    bday: f64,
    #[serde(rename = "Bcity")]
    bcity: f64,
    #[serde(rename = "Bsky")]
    bsky: f64,
}

#[derive(Deserialize)]
struct ObjectLocCase {
    object: String,
    #[serde(rename = "Angle")]
    angle: i32,
    jd_ut: f64,
    dgeo: [f64; 3],
    datm: [f64; 4],
    helflag: u32,
    dret: f64,
}

#[derive(Deserialize)]
struct MagnitudeCase {
    object: String,
    jd_ut: f64,
    dgeo: [f64; 3],
    helflag: u32,
    dmag: f64,
}

#[derive(Deserialize)]
struct AzaltCartCase {
    object: String,
    jd_ut: f64,
    dgeo: [f64; 3],
    datm: [f64; 4],
    helflag: u32,
    dret: [f64; 6],
}

#[derive(Deserialize)]
struct GoldenData {
    extinction: Vec<ExtinctionCase>,
    airmass: Vec<AirmassCase>,
    app_alt: Vec<AppAltCase>,
    optic: Vec<OpticCase>,
    brightness: Vec<BrightnessCase>,
    objectloc: Vec<ObjectLocCase>,
    magnitude: Vec<MagnitudeCase>,
    azaltcart: Vec<AzaltCartCase>,
}

fn load() -> GoldenData {
    let path = super::golden_data_path("heliacal_internals.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn assert_exact_or_eps(label: &str, expected: f64, actual: f64, eps: f64) {
    if actual.to_bits() == expected.to_bits() {
        return;
    }
    super::assert_f64_eps(label, expected, actual, eps);
}

#[test]
fn golden_extinction() {
    let data = load();
    for (i, c) in data.extinction.iter().enumerate() {
        let helflag = HeliacalFlags::empty();
        let label_base = format!(
            "extinction[{i}][AltO={},AltS={},sunra={},Lat={},H={}]",
            c.alt_o, c.alt_s, c.sunra, c.lat, c.height_eye
        );

        let actual_kr = heliacal::kr(c.alt_s, c.height_eye);
        assert_exact_or_eps(&format!("{label_base}/kR"), c.k_r, actual_kr, 1e-12);

        let actual_koz = heliacal::koz(c.alt_s, c.sunra, c.lat);
        assert_exact_or_eps(&format!("{label_base}/kOZ"), c.k_oz, actual_koz, 1e-12);

        let actual_kw = heliacal::kw(c.height_eye, c.datm[1], c.datm[2]);
        assert_exact_or_eps(&format!("{label_base}/kW"), c.k_w, actual_kw, 1e-12);

        let actual_ka = heliacal::ka(
            c.alt_s,
            c.sunra,
            c.lat,
            c.height_eye,
            c.datm[1],
            c.datm[2],
            c.datm[3],
        );
        assert_exact_or_eps(&format!("{label_base}/ka"), c.ka, actual_ka, 1e-12);

        let actual_kt = heliacal::kt(
            c.alt_s,
            c.sunra,
            c.lat,
            c.height_eye,
            c.datm[1],
            c.datm[2],
            c.datm[3],
            4,
        );
        assert_exact_or_eps(&format!("{label_base}/kt"), c.kt, actual_kt, 1e-12);

        let actual_dm = heliacal::deltam(
            c.alt_o,
            c.alt_s,
            c.sunra,
            c.lat,
            c.height_eye,
            &c.datm,
            helflag,
        );
        assert_exact_or_eps(&format!("{label_base}/Deltam"), c.deltam, actual_dm, 1e-12);
    }
}

#[test]
fn golden_airmass() {
    let data = load();
    for (i, c) in data.airmass.iter().enumerate() {
        let label_base = format!("airmass[{i}][AppAltO={},Press={}]", c.app_alt_o, c.press);

        let actual_airm = heliacal::airmass(c.app_alt_o, c.press);
        assert_exact_or_eps(
            &format!("{label_base}/Airmass"),
            c.airmass,
            actual_airm,
            1e-12,
        );

        let zend = {
            let mut z = (90.0 - c.app_alt_o) * std::f64::consts::PI / 180.0;
            if z > std::f64::consts::FRAC_PI_2 {
                z = std::f64::consts::FRAC_PI_2;
            }
            z
        };

        let actual_xr = heliacal::xext(8515.0, zend, c.press);
        assert_exact_or_eps(
            &format!("{label_base}/Xext_rayleigh"),
            c.xext_rayleigh,
            actual_xr,
            1e-12,
        );

        let actual_xw = heliacal::xext(3000.0, zend, c.press);
        assert_exact_or_eps(
            &format!("{label_base}/Xext_water"),
            c.xext_water,
            actual_xw,
            1e-12,
        );

        let actual_xa = heliacal::xext(3745.0, zend, c.press);
        assert_exact_or_eps(
            &format!("{label_base}/Xext_aerosol"),
            c.xext_aerosol,
            actual_xa,
            1e-12,
        );

        let actual_xoz = heliacal::xlay(20000.0, zend, c.press);
        assert_exact_or_eps(
            &format!("{label_base}/Xlay_ozone"),
            c.xlay_ozone,
            actual_xoz,
            1e-12,
        );
    }
}

#[test]
fn golden_app_alt() {
    let data = load();
    for (i, c) in data.app_alt.iter().enumerate() {
        let label_base = format!(
            "app_alt[{i}][alt={},TempE={},PresE={}]",
            c.alt, c.temp_e, c.pres_e
        );

        let actual_app =
            heliacal::app_alt_from_topo_alt(c.alt, c.temp_e, c.pres_e, HeliacalFlags::empty());
        assert_exact_or_eps(
            &format!("{label_base}/AppAltfromTopoAlt"),
            c.app_alt_from_topo_alt,
            actual_app,
            1e-12,
        );

        let actual_topo = heliacal::topo_alt_from_app_alt(c.alt, c.temp_e, c.pres_e);
        assert_exact_or_eps(
            &format!("{label_base}/TopoAltfromAppAlt"),
            c.topo_alt_from_app_alt,
            actual_topo,
            1e-12,
        );
    }
}

#[test]
fn golden_optic() {
    let data = load();
    for (i, c) in data.optic.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(c.helflag);
        let label_base = format!("optic[{i}][B={},config={}]", c.b, c.config);

        let actual_cva = heliacal::cva(c.b, c.dobs[1], helflag);
        assert_exact_or_eps(&format!("{label_base}/CVA"), c.cva, actual_cva, 1e-12);

        let actual_pd = heliacal::pupil_dia(c.dobs[0], c.b);
        assert_exact_or_eps(
            &format!("{label_base}/PupilDia"),
            c.pupil_dia,
            actual_pd,
            1e-12,
        );

        let actual_ofi = heliacal::optic_factor(c.b, c.k_x, &c.dobs, false, 0, helflag);
        assert_exact_or_eps(
            &format!("{label_base}/OpticFactor_intensity"),
            c.optic_factor_intensity,
            actual_ofi,
            1e-12,
        );

        let actual_ofb = heliacal::optic_factor(c.b, c.k_x, &c.dobs, false, 1, helflag);
        assert_exact_or_eps(
            &format!("{label_base}/OpticFactor_background"),
            c.optic_factor_background,
            actual_ofb,
            1e-12,
        );
    }
}

fn assert_rel_or_abs(label: &str, expected: f64, actual: f64, eps: f64) {
    if actual.to_bits() == expected.to_bits() {
        return;
    }
    let diff = (actual - expected).abs();
    if expected.abs() > 1.0 {
        let rel = diff / expected.abs();
        assert!(
            rel < eps,
            "{label}: expected {expected:.17e}, got {actual:.17e} (rel err {rel:.3e} > {eps:.1e})"
        );
    } else {
        assert!(
            diff < eps,
            "{label}: expected {expected:.17e}, got {actual:.17e} (abs err {diff:.3e} > {eps:.1e})"
        );
    }
}

#[test]
fn golden_brightness() {
    let data = load();
    let helflag = HeliacalFlags::empty();
    for (i, c) in data.brightness.iter().enumerate() {
        let label_base = format!(
            "brightness[{i}][AltO={},AltS={},AltM={},sunra={}]",
            c.alt_o, c.alt_s, c.alt_m, c.sunra
        );

        let actual_bn = heliacal::bn(
            c.alt_o,
            c.jdn_days_ut,
            c.alt_s,
            c.sunra,
            c.lat,
            c.height_eye,
            &c.datm,
            helflag,
        );
        assert_rel_or_abs(&format!("{label_base}/Bn"), c.bn, actual_bn, 1e-12);

        let actual_bm = heliacal::bm(
            c.alt_o,
            c.azi_o,
            c.alt_m,
            c.azi_m,
            c.alt_s,
            c.azi_s,
            c.sunra,
            c.lat,
            c.height_eye,
            &c.datm,
            helflag,
        );
        assert_rel_or_abs(&format!("{label_base}/Bm"), c.bm, actual_bm, 1e-12);

        let actual_btwi = heliacal::btwi(
            c.alt_o,
            c.azi_o,
            c.alt_s,
            c.azi_s,
            c.sunra,
            c.lat,
            c.height_eye,
            &c.datm,
            helflag,
        );
        assert_rel_or_abs(&format!("{label_base}/Btwi"), c.btwi, actual_btwi, 1e-12);

        let actual_bday = heliacal::bday(
            c.alt_o,
            c.azi_o,
            c.alt_s,
            c.azi_s,
            c.sunra,
            c.lat,
            c.height_eye,
            &c.datm,
            helflag,
        );
        assert_rel_or_abs(&format!("{label_base}/Bday"), c.bday, actual_bday, 1e-12);

        let actual_bcity = heliacal::bcity(0.0);
        assert_rel_or_abs(&format!("{label_base}/Bcity"), c.bcity, actual_bcity, 1e-12);

        let actual_bsky = heliacal::bsky(
            c.alt_o,
            c.azi_o,
            c.alt_m,
            c.azi_m,
            c.jdn_days_ut,
            c.alt_s,
            c.azi_s,
            c.sunra,
            c.lat,
            c.height_eye,
            &c.datm,
            helflag,
        );
        assert_rel_or_abs(&format!("{label_base}/Bsky"), c.bsky, actual_bsky, 1e-12);
    }
}

fn make_eph() -> Ephemeris {
    let config = EphemerisConfig {
        ephemeris_source: swisseph::EphemerisSource::Swiss,
        ephe_path: Some("../swisseph/ephe".into()),
        topographic: Some(swisseph::config::TopoPosition {
            longitude: 31.25,
            latitude: 30.1,
            altitude: 30.0,
        }),
        ..Default::default()
    };
    Ephemeris::new(config).unwrap()
}

#[test]
fn golden_objectloc() {
    let data = load();
    let eph = make_eph();
    for (i, c) in data.objectloc.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(c.helflag);
        let epheflag = CalcFlags::from_bits_truncate(c.helflag);
        let label_base = format!(
            "objectloc[{i}][obj={},Angle={},jd={}]",
            c.object, c.angle, c.jd_ut
        );

        let actual = heliacal::object_loc(
            &eph, c.jd_ut, &c.dgeo, &c.datm, &c.object, c.angle, epheflag, helflag,
        );
        match actual {
            Ok(val) => {
                super::assert_f64_eps(&label_base, c.dret, val, 1e-7);
            }
            Err(e) => panic!("{label_base}: unexpected error: {e}"),
        }
    }
}

#[test]
fn golden_magnitude() {
    let data = load();
    let eph = make_eph();
    for (i, c) in data.magnitude.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(c.helflag);
        let epheflag = CalcFlags::from_bits_truncate(c.helflag);
        let label_base = format!("magnitude[{i}][obj={},jd={}]", c.object, c.jd_ut);

        let actual = heliacal::magnitude(&eph, c.jd_ut, &c.dgeo, &c.object, epheflag, helflag);
        match actual {
            Ok(val) => {
                super::assert_f64_eps(&label_base, c.dmag, val, 1e-8);
            }
            Err(e) => panic!("{label_base}: unexpected error: {e}"),
        }
    }
}

#[test]
fn golden_azaltcart() {
    let data = load();
    let eph = make_eph();
    for (i, c) in data.azaltcart.iter().enumerate() {
        let helflag = HeliacalFlags::from_bits_truncate(c.helflag);
        let epheflag = CalcFlags::from_bits_truncate(c.helflag);
        let label_base = format!("azaltcart[{i}][obj={},jd={}]", c.object, c.jd_ut);

        let actual = heliacal::azalt_cart(
            &eph, c.jd_ut, &c.dgeo, &c.datm, &c.object, epheflag, helflag,
        );
        match actual {
            Ok(val) => {
                let labels = ["az", "topo_alt", "app_alt", "cart_x", "cart_y", "cart_z"];
                for (j, &lbl) in labels.iter().enumerate() {
                    super::assert_f64_eps(&format!("{label_base}/{lbl}"), c.dret[j], val[j], 1e-7);
                }
            }
            Err(e) => panic!("{label_base}: unexpected error: {e}"),
        }
    }
}
