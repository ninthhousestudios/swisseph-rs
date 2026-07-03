#![allow(clippy::too_many_arguments)]

use crate::constants::{AST_OFFSET, DEGTORAD};
use crate::date::revjul;
use crate::flags::HeliacalFlags;
use crate::types::{Body, CalendarType};

// ── Heliacal event types ───────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum HeliacalEventType {
    MorningFirst = 1,
    EveningLast = 2,
    EveningFirst = 3,
    MorningLast = 4,
    AcronymchalRising = 5,
    AcronymchalSetting = 6,
}

impl TryFrom<i32> for HeliacalEventType {
    type Error = crate::Error;
    fn try_from(v: i32) -> crate::Result<Self> {
        match v {
            1 => Ok(Self::MorningFirst),
            2 => Ok(Self::EveningLast),
            3 => Ok(Self::EveningFirst),
            4 => Ok(Self::MorningLast),
            5 => Ok(Self::AcronymchalRising),
            6 => Ok(Self::AcronymchalSetting),
            _ => Err(crate::Error::InvalidBody(v)),
        }
    }
}

// ── Constants (live only — swehel.c:76–200) ────────────────────────
// Staged for sub-task 3/8 (VisLimMagn)
#[allow(dead_code)]
const BNIGHT: f64 = 1479.0; // [nL]
#[allow(dead_code)]
const BNIGHT_FACTOR: f64 = 1.0;
const NL2ERG: f64 = 1.02e-15;
const ERG2NL: f64 = 1.0 / NL2ERG;
const MOON_DISTANCE: f64 = 384410.4978; // [km]

const SCALE_H_WATER: f64 = 3000.0; // [m]
const SCALE_H_RAYLEIGH: f64 = 8515.0; // [m]
const SCALE_H_AEROSOL: f64 = 3745.0; // [m]
const SCALE_H_OZONE: f64 = 20000.0; // [m]

const ASTR2TAU: f64 = 0.921034037197618; // ln(10^0.4)
const TAU2ASTR: f64 = 1.0 / ASTR2TAU;

const C2K: f64 = 273.15; // [K]
const LAPSE_SA: f64 = 0.0065; // [K/m] standard atmosphere

const LOWEST_APP_ALT: f64 = -3.5; // [deg]
const MIN2DEG: f64 = 1.0 / 60.0;

const RA: f64 = 6378136.6; // [m] WGS84 equatorial radius

// ── Object resolution (swehel.c:305–336) ───────────────────────────

pub fn object_to_body(name: &str) -> Option<Body> {
    let lower = name.to_ascii_lowercase();
    if lower.starts_with("sun") {
        return Some(Body::Sun);
    }
    if lower.starts_with("venus") {
        return Some(Body::Venus);
    }
    if lower.starts_with("mars") {
        return Some(Body::Mars);
    }
    if lower.starts_with("mercur") {
        return Some(Body::Mercury);
    }
    if lower.starts_with("jupiter") {
        return Some(Body::Jupiter);
    }
    if lower.starts_with("saturn") {
        return Some(Body::Saturn);
    }
    if lower.starts_with("uranus") {
        return Some(Body::Uranus);
    }
    if lower.starts_with("neptun") {
        return Some(Body::Neptune);
    }
    if lower.starts_with("moon") {
        return Some(Body::Moon);
    }
    // C uses atoi(s) which parses leading digits, ignoring trailing text
    let leading: String = lower.chars().take_while(|c| c.is_ascii_digit()).collect();
    if !leading.is_empty() {
        if let Ok(n) = leading.parse::<i32>() {
            if n > 0 {
                return Body::try_from(n + AST_OFFSET).ok();
            }
        }
    }
    None
}

pub fn tolower_string_star(name: &str) -> String {
    if let Some(comma_pos) = name.find(',') {
        let mut result = name[..comma_pos].to_ascii_lowercase();
        result.push_str(&name[comma_pos..]);
        result
    } else {
        name.to_ascii_lowercase()
    }
}

// ── Default heliacal parameters (swehel.c:1324–1361) ───────────────

pub fn default_heliacal_parameters(
    datm: &mut [f64; 4],
    dgeo: &[f64; 3],
    dobs: &mut [f64; 6],
    helflag: HeliacalFlags,
) {
    if datm[0] <= 0.0 {
        // ISA pressure estimate
        datm[0] = 1013.25 * (1.0 - 0.0065 * dgeo[2] / 288.0).powf(5.255);
        if datm[1] == 0.0 {
            datm[1] = 15.0 - 0.0065 * dgeo[2];
        }
        if datm[2] == 0.0 {
            datm[2] = 40.0;
        }
    }
    // SIMULATE_VICTORVB always defined → the #ifndef block (RH clamp in else branch) is DEAD

    if dobs[0] == 0.0 {
        dobs[0] = 36.0;
    }
    if dobs[1] == 0.0 {
        dobs[1] = 1.0;
    }
    if !helflag.contains(HeliacalFlags::OPTICAL_PARAMS) {
        for i in 2..=5 {
            dobs[i] = 0.0;
        }
    }
    if dobs[3] == 0.0 {
        dobs[2] = 1.0; // Binocular = 1
        dobs[3] = 1.0; // OpticMagn = 1: use eye
    }
}

// ── Meteorological / coordinate helpers (swehel.c §3) ──────────────

fn mymin(a: f64, b: f64) -> f64 {
    if a <= b { a } else { b }
}

fn mymax(a: f64, b: f64) -> f64 {
    if a >= b { a } else { b }
}

pub fn tanh_manual(x: f64) -> f64 {
    (x.exp() - (-x).exp()) / (x.exp() + (-x).exp())
}

pub fn kelvin(temp: f64) -> f64 {
    temp + C2K
}

pub fn topo_alt_from_app_alt(app_alt: f64, temp_e: f64, pres_e: f64) -> f64 {
    if app_alt >= LOWEST_APP_ALT {
        let r = if app_alt > 17.904104638432 {
            0.97 / (app_alt * DEGTORAD).tan()
        } else {
            (34.46 + 4.23 * app_alt + 0.004 * app_alt * app_alt)
                / (1.0 + 0.505 * app_alt + 0.0845 * app_alt * app_alt)
        };
        let r = (pres_e - 80.0) / 930.0 / (1.0 + 0.00008 * (r + 39.0) * (temp_e - 10.0)) * r;
        app_alt - r * MIN2DEG
    } else {
        app_alt
    }
}

pub fn app_alt_from_topo_alt(
    topo_alt: f64,
    temp_e: f64,
    pres_e: f64,
    helflag: HeliacalFlags,
) -> f64 {
    let nloop: i32 = if helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        5
    } else {
        2
    };
    let mut new_app_alt = topo_alt;
    let mut new_topo_alt = 0.0;
    let mut oud_app_alt = new_app_alt;
    let mut oud_topo_alt = new_topo_alt;

    for _i in 0..=nloop {
        new_topo_alt = new_app_alt - topo_alt_from_app_alt(new_app_alt, temp_e, pres_e);
        let mut verschil = new_app_alt - oud_app_alt;
        oud_app_alt = new_topo_alt - oud_topo_alt - verschil;
        if verschil != 0.0 && oud_app_alt != 0.0 {
            verschil =
                new_app_alt - verschil * (topo_alt + new_topo_alt - new_app_alt) / oud_app_alt;
        } else {
            verschil = topo_alt + new_topo_alt;
        }
        oud_app_alt = new_app_alt;
        oud_topo_alt = new_topo_alt;
        new_app_alt = verschil;
    }

    let retalt = topo_alt + new_topo_alt;
    if retalt < LOWEST_APP_ALT {
        topo_alt
    } else {
        retalt
    }
}

pub fn hour_angle(topo_alt: f64, topo_decl: f64, lat: f64) -> f64 {
    let alti = topo_alt * DEGTORAD;
    let decli = topo_decl * DEGTORAD;
    let lati = lat * DEGTORAD;
    let mut ha = (alti.sin() - lati.sin() * decli.sin()) / lati.cos() / decli.cos();
    if ha < -1.0 {
        ha = -1.0;
    }
    if ha > 1.0 {
        ha = 1.0;
    }
    ha.acos() / DEGTORAD / 15.0
}

pub fn distance_angle(lat_a: f64, long_a: f64, lat_b: f64, long_b: f64) -> f64 {
    let dlon = long_b - long_a;
    let dlat = lat_b - lat_a;
    let sindlat2 = (dlat / 2.0).sin();
    let sindlon2 = (dlon / 2.0).sin();
    let mut corde = sindlat2 * sindlat2 + lat_a.cos() * lat_b.cos() * sindlon2 * sindlon2;
    if corde > 1.0 {
        corde = 1.0;
    }
    2.0 * corde.sqrt().asin()
}

pub fn temp_e_from_temp_s(temp_s: f64, height_eye: f64, lapse: f64) -> f64 {
    temp_s - lapse * height_eye
}

pub fn pres_e_from_pres_s(temp_s: f64, press: f64, height_eye: f64) -> f64 {
    press
        * (-9.80665 * 0.0289644 / (kelvin(temp_s) + 3.25 * height_eye / 1000.0) / 8.31441
            * height_eye)
            .exp()
}

fn sgn(x: f64) -> f64 {
    if x < 0.0 { -1.0 } else { 1.0 }
}

// ── Extinction layer (swehel.c §4) ─────────────────────────────────

pub fn kw(height_eye: f64, temp_s: f64, rh: f64) -> f64 {
    let mut wt = 0.031;
    wt *= 0.94 * (rh / 100.0) * (temp_s / 15.0).exp() * (-1.0 * height_eye / SCALE_H_WATER).exp();
    wt
}

pub fn koz(alt_s: f64, sunra: f64, lat: f64) -> f64 {
    let oz = 0.031;
    let lt = lat * DEGTORAD;
    let koz_ret = oz * (3.0 + 0.4 * (lt * (sunra * DEGTORAD).cos() - (3.0 * lt).cos())) / 3.0;
    let mut altslim = -alt_s - 12.0;
    if altslim < 0.0 {
        altslim = 0.0;
    }
    let changeko = (100.0 - 11.6 * mymin(6.0, altslim)) / 100.0;
    koz_ret * changeko
}

pub fn kr(alt_s: f64, height_eye: f64) -> f64 {
    let mut val = -alt_s - 12.0;
    if val < 0.0 {
        val = 0.0;
    }
    if val > 6.0 {
        val = 6.0;
    }
    let changek = 1.0 - 0.166667 * val;
    let lambda = 0.55 + (changek - 1.0) * 0.04;
    0.1066 * (-1.0 * height_eye / SCALE_H_RAYLEIGH).exp() * (lambda / 0.55_f64).powf(-4.0)
}

pub fn ka(alt_s: f64, sunra: f64, lat: f64, height_eye: f64, temp_s: f64, rh: f64, vr: f64) -> f64 {
    let sl = sgn(lat);
    let changeka = 1.0 - 0.166667 * mymin(6.0, mymax(-alt_s - 12.0, 0.0));
    let lambda = 0.55 + (changeka - 1.0) * 0.04;

    let kaact;
    if vr != 0.0 {
        if vr >= 1.0 {
            let beta_vr = 3.912 / vr;
            let betaa = beta_vr
                - (kw(height_eye, temp_s, rh) / SCALE_H_WATER
                    + kr(alt_s, height_eye) / SCALE_H_RAYLEIGH)
                    * 1000.0
                    * ASTR2TAU;
            kaact = betaa * SCALE_H_AEROSOL / 1000.0 * TAU2ASTR;
        } else {
            kaact =
                vr - kw(height_eye, temp_s, rh) - kr(alt_s, height_eye) - koz(alt_s, sunra, lat);
        }
    } else {
        // SIMULATE_VICTORVB is always active — clamp RH
        let mut rh_clamped = rh;
        if rh_clamped <= 0.00000001 {
            rh_clamped = 0.00000001;
        }
        if rh_clamped >= 99.99999999 {
            rh_clamped = 99.99999999;
        }
        let base = 0.1
            * (-1.0 * height_eye / SCALE_H_AEROSOL).exp()
            * (1.0 - 0.32 / (rh_clamped / 100.0).ln()).powf(1.33)
            * (1.0 + 0.33 * sl * (sunra * DEGTORAD).sin());
        kaact = base * (lambda / 0.55_f64).powf(-1.3);
    }
    kaact
}

pub fn kt(
    alt_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    temp_s: f64,
    rh: f64,
    vr: f64,
    ext_type: i32,
) -> f64 {
    let mut kw_act = 0.0;
    let mut kr_act = 0.0;
    let mut koz_act = 0.0;
    let mut ka_act = 0.0;
    match ext_type {
        0 => ka_act = ka(alt_s, sunra, lat, height_eye, temp_s, rh, vr),
        1 => kw_act = kw(height_eye, temp_s, rh),
        2 => kr_act = kr(alt_s, height_eye),
        3 => koz_act = koz(alt_s, sunra, lat),
        4 => {
            ka_act = ka(alt_s, sunra, lat, height_eye, temp_s, rh, vr);
            kw_act = kw(height_eye, temp_s, rh);
            kr_act = kr(alt_s, height_eye);
            koz_act = koz(alt_s, sunra, lat);
        }
        _ => {}
    }
    if ka_act < 0.0 {
        ka_act = 0.0;
    }
    kw_act + kr_act + koz_act + ka_act
}

pub fn airmass(app_alt_o: f64, press: f64) -> f64 {
    let mut zend = (90.0 - app_alt_o) * DEGTORAD;
    if zend > std::f64::consts::FRAC_PI_2 {
        zend = std::f64::consts::FRAC_PI_2;
    }
    let airm = 1.0 / (zend.cos() + 0.025 * (-11.0 * zend.cos()).exp());
    press / 1013.0 * airm
}

pub fn xext(scale_h: f64, zend: f64, press: f64) -> f64 {
    press
        / 1013.0
        / (zend.cos()
            + 0.01
                * (scale_h / 1000.0).sqrt()
                * (-30.0 / (scale_h / 1000.0).sqrt() * zend.cos()).exp())
}

pub fn xlay(scale_h: f64, zend: f64, press: f64) -> f64 {
    let a = zend.sin() / (1.0 + scale_h / RA);
    press / 1013.0 / (1.0 - a * a).sqrt()
}

pub fn deltam(
    alt_o: f64,
    alt_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let pres_e = pres_e_from_pres_s(datm[1], datm[0], height_eye);
    let temp_e = temp_e_from_temp_s(datm[1], height_eye, LAPSE_SA);
    let app_alt_o = app_alt_from_topo_alt(alt_o, temp_e, pres_e, helflag);

    // staticAirmass == 0 → always take this branch
    let mut zend = (90.0 - app_alt_o) * DEGTORAD;
    if zend > std::f64::consts::FRAC_PI_2 {
        zend = std::f64::consts::FRAC_PI_2;
    }
    let x_r = xext(SCALE_H_RAYLEIGH, zend, datm[0]);
    let x_w = xext(SCALE_H_WATER, zend, datm[0]);
    let x_a = xext(SCALE_H_AEROSOL, zend, datm[0]);
    let x_oz = xlay(SCALE_H_OZONE, zend, datm[0]);

    kr(alt_s, height_eye) * x_r
        + kt(alt_s, sunra, lat, height_eye, datm[1], datm[2], datm[3], 0) * x_a
        + koz(alt_s, sunra, lat) * x_oz
        + kw(height_eye, datm[1], datm[2]) * x_w
}

// ── Optics & vision helpers (swehel.c §5) ──────────────────────────

pub fn cva(b: f64, sn: f64, helflag: HeliacalFlags) -> f64 {
    let mut is_scotopic = b < 1394.0;
    if helflag.contains(HeliacalFlags::VISLIM_PHOTOPIC) {
        is_scotopic = false;
    }
    if helflag.contains(HeliacalFlags::VISLIM_SCOTOPIC) {
        is_scotopic = true;
    }
    if is_scotopic {
        mymin(900.0, 380.0 / sn * 10.0_f64.powf(0.3 * b.powf(-0.29))) / 60.0 / 60.0
    } else {
        (40.0 / sn) * 10.0_f64.powf(8.28 * b.powf(-0.29)) / 60.0 / 60.0
    }
}

pub fn pupil_dia(age: f64, b: f64) -> f64 {
    (0.534
        - 0.00211 * age
        - (0.236 - 0.00127 * age) * tanh_manual(0.4 * b.ln() / 10.0_f64.ln() - 2.2))
        * 10.0
}

pub fn optic_factor(
    bback: f64,
    k_x: f64,
    dobs: &[f64; 6],
    is_moon: bool,
    type_factor: i32,
    helflag: HeliacalFlags,
) -> f64 {
    let age = dobs[0];
    let sn = dobs[1];
    let sni = if sn < 1e-8 { 1e-8 } else { sn };
    let binocular = dobs[2];
    let optic_mag = dobs[3];
    let mut optic_dia = dobs[4];
    let mut optic_trans = dobs[5];
    let _ = is_moon; // ObjectName "moon" check is a no-op in C

    let pst = pupil_dia(23.0, bback);

    if optic_mag == 1.0 {
        optic_trans = 1.0;
        optic_dia = pst;
    }

    let cib = 0.7;
    let _cii = 0.5;
    let object_size = 0.0;

    let fb = if binocular == 0.0 { 1.41 } else { 1.0 };

    let mut is_scotopic = bback < 1645.0;
    if helflag.contains(HeliacalFlags::VISLIM_PHOTOPIC) {
        is_scotopic = false;
    }
    if helflag.contains(HeliacalFlags::VISLIM_SCOTOPIC) {
        is_scotopic = true;
    }

    let (fe, fsc, fci, fcb);
    if is_scotopic {
        fe = 10.0_f64.powf(0.48 * k_x);
        fsc = mymin(
            1.0,
            (1.0 - (pst / 124.4_f64).powf(4.0)) / (1.0 - (optic_dia / optic_mag / 124.4).powf(4.0)),
        );
        fci = 10.0_f64.powf(-0.4 * (1.0 - _cii / 2.0));
        fcb = 10.0_f64.powf(-0.4 * (1.0 - cib / 2.0));
    } else {
        fe = 10.0_f64.powf(0.4 * k_x);
        fsc = mymin(
            1.0,
            (optic_dia / optic_mag / pst).powf(2.0) * (1.0 - (-(pst / 6.2_f64).powf(2.0)).exp())
                / (1.0 - (-(optic_dia / optic_mag / 6.2).powf(2.0)).exp()),
        );
        fci = 1.0;
        fcb = 1.0;
    }

    let ft = 1.0 / optic_trans;
    let fp = mymax(1.0, (pst / (optic_mag * pupil_dia(age, bback))).powf(2.0));
    let fa = (pst / optic_dia).powf(2.0);
    let fr = (1.0 + 0.03 * (optic_mag * object_size / cva(bback, sni, helflag)).powf(2.0))
        / sni.powf(2.0);
    let fm = optic_mag.powf(2.0);

    if type_factor == 0 {
        fb * fe * ft * fp * fa * fr * fsc * fci
    } else {
        fb * ft * fp * fa * fm * fsc * fcb
    }
}

// ── Sky brightness model (swehel.c §6) ────────────────────────────

#[allow(clippy::approx_constant)]
const LN10: f64 = 2.302585092994;

pub fn moons_brightness(dist: f64, phasemoon: f64) -> f64 {
    -21.62
        + 5.0 * (dist / (RA / 1000.0)).ln() / LN10
        + 0.026 * phasemoon.abs()
        + 0.000000004 * phasemoon.powi(4)
}

pub fn moon_phase(alt_m: f64, azi_m: f64, alt_s: f64, azi_s: f64) -> f64 {
    let moon_avg_par = 0.95;
    let alt_mi = (alt_m + moon_avg_par) * DEGTORAD;
    let alt_si = alt_s * DEGTORAD;
    let azi_mi = azi_m * DEGTORAD;
    let azi_si = azi_s * DEGTORAD;
    180.0
        - ((azi_si - azi_mi - moon_avg_par * DEGTORAD).cos() * alt_mi.cos() * alt_si.cos()
            + alt_si.sin() * alt_mi.sin())
        .acos()
            / DEGTORAD
}

#[allow(clippy::approx_constant)]
pub fn bn(
    alt_o: f64,
    jdn_days_ut: f64,
    alt_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let pres_e = pres_e_from_pres_s(datm[1], datm[0], height_eye);
    let temp_e = temp_e_from_temp_s(datm[1], height_eye, LAPSE_SA);
    let mut app_alt_o = app_alt_from_topo_alt(alt_o, temp_e, pres_e, helflag);
    if app_alt_o < 10.0 {
        app_alt_o = 10.0;
    }
    let zend = (90.0 - app_alt_o) * DEGTORAD;

    let (iyar, imon, iday, _dut) = revjul(jdn_days_ut, CalendarType::Gregorian);
    let year_b = iyar as f64;
    let month_b = imon as f64;
    let day_b = iday as f64;

    let b0: f64 = 0.0000000000001;
    let bna = b0
        * (1.0
            + 0.3
                * (6.283 * (year_b + ((day_b - 1.0) / 30.4 + month_b - 1.0) / 12.0 - 1990.33)
                    / 11.1)
                    .cos());

    let k_x = deltam(alt_o, alt_s, sunra, lat, height_eye, datm, helflag);
    let bnb =
        bna * (0.4 + 0.6 / (1.0 - 0.96 * zend.sin().powi(2)).sqrt()) * 10.0_f64.powf(-0.4 * k_x);

    mymax(bnb, 0.0) * ERG2NL
}

pub fn bm(
    alt_o: f64,
    azi_o: f64,
    alt_m: f64,
    azi_m: f64,
    alt_s: f64,
    azi_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let m0 = -11.05;
    let lunar_radius = 0.25 * DEGTORAD;
    let object_is_moon = alt_o == alt_m && azi_o == azi_m;
    let mut bm_val = 0.0;

    if alt_m > -0.26 && !object_is_moon {
        let mut rm = distance_angle(
            alt_o * DEGTORAD,
            azi_o * DEGTORAD,
            alt_m * DEGTORAD,
            azi_m * DEGTORAD,
        ) / DEGTORAD;
        if rm <= lunar_radius {
            rm = lunar_radius;
        }

        let k_xm = deltam(alt_m, alt_s, sunra, lat, height_eye, datm, helflag);
        let k_x = deltam(alt_o, alt_s, sunra, lat, height_eye, datm, helflag);
        let c3 = 10.0_f64.powf(-0.4 * k_xm);
        let fm = 62000000.0 / rm / rm
            + 10.0_f64.powf(6.15 - rm / 40.0)
            + 10.0_f64.powf(5.36) * (1.06 + (rm * DEGTORAD).cos().powi(2));
        bm_val = fm * c3 + 440000.0 * (1.0 - c3);

        let phasemoon = moon_phase(alt_m, azi_m, alt_s, azi_s);
        let mm = moons_brightness(MOON_DISTANCE, phasemoon);
        bm_val *= 10.0_f64.powf(-0.4 * (mm - m0 + 43.27));
        bm_val *= 1.0 - 10.0_f64.powf(-0.4 * k_x);
    }

    mymax(bm_val, 0.0) * ERG2NL
}

pub fn btwi(
    alt_o: f64,
    azi_o: f64,
    alt_s: f64,
    azi_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let m0 = -11.05;
    let ms = -26.74;

    let pres_e = pres_e_from_pres_s(datm[1], datm[0], height_eye);
    let temp_e = temp_e_from_temp_s(datm[1], height_eye, LAPSE_SA);
    let app_alt_o = app_alt_from_topo_alt(alt_o, temp_e, pres_e, helflag);
    let zend_o = 90.0 - app_alt_o;

    let rs = distance_angle(
        alt_o * DEGTORAD,
        azi_o * DEGTORAD,
        alt_s * DEGTORAD,
        azi_s * DEGTORAD,
    ) / DEGTORAD;

    let k_x = deltam(alt_o, alt_s, sunra, lat, height_eye, datm, helflag);
    let k = kt(alt_s, sunra, lat, height_eye, datm[1], datm[2], datm[3], 4);

    let mut btwi_val = 10.0_f64.powf(-0.4 * (ms - m0 + 32.5 - alt_s - (zend_o / (360.0 * k))));
    btwi_val = btwi_val * (100.0 / rs) * (1.0 - 10.0_f64.powf(-0.4 * k_x));

    mymax(btwi_val, 0.0) * ERG2NL
}

pub fn bday(
    alt_o: f64,
    azi_o: f64,
    alt_s: f64,
    azi_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let m0 = -11.05;
    let ms = -26.74;

    let rs = distance_angle(
        alt_o * DEGTORAD,
        azi_o * DEGTORAD,
        alt_s * DEGTORAD,
        azi_s * DEGTORAD,
    ) / DEGTORAD;

    let k_xs = deltam(alt_s, alt_s, sunra, lat, height_eye, datm, helflag);
    let k_x = deltam(alt_o, alt_s, sunra, lat, height_eye, datm, helflag);
    let c4 = 10.0_f64.powf(-0.4 * k_xs);
    let fs = 62000000.0 / rs / rs
        + 10.0_f64.powf(6.15 - rs / 40.0)
        + 10.0_f64.powf(5.36) * (1.06 + (rs * DEGTORAD).cos().powi(2));
    let mut bday_val = fs * c4 + 440000.0 * (1.0 - c4);
    bday_val *= 10.0_f64.powf(-0.4 * (ms - m0 + 43.27));
    bday_val *= 1.0 - 10.0_f64.powf(-0.4 * k_x);

    mymax(bday_val, 0.0) * ERG2NL
}

pub fn bcity(value: f64) -> f64 {
    mymax(value, 0.0)
}

pub fn bsky(
    alt_o: f64,
    azi_o: f64,
    alt_m: f64,
    azi_m: f64,
    jdn_days_ut: f64,
    alt_s: f64,
    azi_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let mut bsky_val = 0.0;

    if alt_s < -3.0 {
        bsky_val += btwi(
            alt_o, azi_o, alt_s, azi_s, sunra, lat, height_eye, datm, helflag,
        );
    } else if alt_s > 4.0 {
        bsky_val += bday(
            alt_o, azi_o, alt_s, azi_s, sunra, lat, height_eye, datm, helflag,
        );
    } else {
        bsky_val += mymin(
            bday(
                alt_o, azi_o, alt_s, azi_s, sunra, lat, height_eye, datm, helflag,
            ),
            btwi(
                alt_o, azi_o, alt_s, azi_s, sunra, lat, height_eye, datm, helflag,
            ),
        );
    }

    if bsky_val < 200000000.0 {
        bsky_val += bm(
            alt_o, azi_o, alt_m, azi_m, alt_s, azi_s, sunra, lat, height_eye, datm, helflag,
        );
    }

    if alt_s <= 0.0 {
        bsky_val += bcity(0.0);
    }

    if bsky_val < 5000.0 {
        bsky_val += bn(
            alt_o,
            jdn_days_ut,
            alt_s,
            sunra,
            lat,
            height_eye,
            datm,
            helflag,
        );
    }

    bsky_val
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn object_to_body_planets() {
        assert_eq!(object_to_body("Sun"), Some(Body::Sun));
        assert_eq!(object_to_body("VENUS"), Some(Body::Venus));
        assert_eq!(object_to_body("mercury"), Some(Body::Mercury));
        assert_eq!(object_to_body("Mercurius"), Some(Body::Mercury));
        assert_eq!(object_to_body("neptune"), Some(Body::Neptune));
        assert_eq!(object_to_body("neptunus"), Some(Body::Neptune));
        assert_eq!(object_to_body("Moon"), Some(Body::Moon));
        assert_eq!(object_to_body("Jupiter"), Some(Body::Jupiter));
        assert_eq!(object_to_body("Saturn"), Some(Body::Saturn));
        assert_eq!(object_to_body("Uranus"), Some(Body::Uranus));
        assert_eq!(object_to_body("Mars"), Some(Body::Mars));
    }

    #[test]
    fn object_to_body_asteroid() {
        let b = object_to_body("433").unwrap();
        assert_eq!(b.to_raw_id(), 433 + AST_OFFSET);
        // C atoi semantics: leading digits with trailing text
        let b2 = object_to_body("433 Eros").unwrap();
        assert_eq!(b2.to_raw_id(), 433 + AST_OFFSET);
        let b3 = object_to_body("433, Eros").unwrap();
        assert_eq!(b3.to_raw_id(), 433 + AST_OFFSET);
    }

    #[test]
    fn object_to_body_star() {
        assert_eq!(object_to_body("Aldebaran"), None);
        assert_eq!(object_to_body("0"), None);
        assert_eq!(object_to_body(""), None);
        // Non-numeric leading text → star
        assert_eq!(object_to_body("alpha Tau"), None);
    }

    #[test]
    fn tolower_preserves_bayer() {
        assert_eq!(tolower_string_star("SIRIUS,alCMa"), "sirius,alCMa");
        assert_eq!(tolower_string_star("Aldebaran"), "aldebaran");
    }

    #[test]
    fn sgn_zero_is_positive() {
        assert_eq!(sgn(0.0), 1.0);
        assert_eq!(sgn(-1.0), -1.0);
        assert_eq!(sgn(5.0), 1.0);
    }

    #[test]
    fn kelvin_basic() {
        assert_eq!(kelvin(0.0), 273.15);
        assert_eq!(kelvin(15.0), 288.15);
    }
}
