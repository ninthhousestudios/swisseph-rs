//! Heliacal visibility: first/last sightings of a body near the Sun.
//!
//! Pure-Rust port of the Swiss Ephemeris heliacal module (`swehel.c`): visual
//! limiting-magnitude modeling (atmospheric extinction, sky background
//! brightness, and visual physiology) plus the search machinery for heliacal
//! risings/settings and acronychal (cosmical) risings/settings of planets,
//! the Moon, and fixed stars.

#![allow(clippy::too_many_arguments)]

use crate::azalt::AzAltDir;
use crate::calc;
use crate::config::TopoPosition;
use crate::constants::{AST_OFFSET, DEGTORAD};
use crate::context::Ephemeris;
use crate::date::revjul;
use crate::error::Error;
use crate::flags::{CalcFlags, HeliacalFlags, RiseSetFlags, VisLimFlags};
use crate::math::{normalize_degrees, polar_to_cartesian};
use crate::types::{Body, CalendarType};

// ── Heliacal event types ───────────────────────────────────────────

/// Classifies which heliacal or acronychal event a search targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum HeliacalEventType {
    /// First visibility in the morning sky before sunrise.
    MorningFirst = 1,
    /// Last visibility in the evening sky after sunset.
    EveningLast = 2,
    /// First visibility in the evening sky after sunset.
    EveningFirst = 3,
    /// Last visibility in the morning sky before sunrise.
    MorningLast = 4,
    /// Acronychal (cosmical) rising: the body rises as the Sun sets.
    AcronymchalRising = 5,
    /// Acronychal (cosmical) setting: the body sets as the Sun rises.
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

/// Resolves an object name (planet name, "moon", or a leading asteroid number)
/// to a `Body`; returns `None` for fixed stars.
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
    if !leading.is_empty()
        && let Ok(n) = leading.parse::<i32>()
        && n > 0
    {
        return Body::try_from(n + AST_OFFSET).ok();
    }
    None
}

/// Lowercases an object name, preserving the case of any comma-separated
/// star-designation suffix (e.g. component letter).
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

/// Fills in unset atmospheric (`datm`) and observer (`dobs`) parameters with
/// their standard-atmosphere / naked-eye defaults.
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

/// Internal: hyperbolic tangent computed from `exp`, matching the C
/// implementation's expression order for bit-for-bit fidelity.
pub fn tanh_manual(x: f64) -> f64 {
    (x.exp() - (-x).exp()) / (x.exp() + (-x).exp())
}

/// Converts a Celsius temperature to Kelvin.
pub fn kelvin(temp: f64) -> f64 {
    temp + C2K
}

/// Converts an apparent (refracted) altitude to a topocentric (unrefracted)
/// altitude, given the effective temperature and pressure at the eye.
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

/// Converts a topocentric (unrefracted) altitude to an apparent (refracted)
/// altitude via iterative refinement (2 iterations, or 5 under
/// `HIGH_PRECISION`).
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

/// Computes the hour angle (in hours) at which an object of given
/// declination reaches a target topocentric altitude, for an observer at
/// the given latitude.
pub fn hour_angle(topo_alt: f64, topo_decl: f64, lat: f64) -> f64 {
    let alti = topo_alt * DEGTORAD;
    let decli = topo_decl * DEGTORAD;
    let lati = lat * DEGTORAD;
    let mut ha = (alti.sin() - lati.sin() * decli.sin()) / lati.cos() / decli.cos();
    ha = ha.clamp(-1.0, 1.0);
    ha.acos() / DEGTORAD / 15.0
}

/// Computes the great-circle angular distance (radians) between two points
/// given as (latitude, longitude) pairs in radians.
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

/// Internal: computes the effective temperature at the observer's eye from
/// surface temperature, eye height, and the atmospheric lapse rate.
pub fn temp_e_from_temp_s(temp_s: f64, height_eye: f64, lapse: f64) -> f64 {
    temp_s - lapse * height_eye
}

/// Internal: computes the effective pressure at the observer's eye from
/// surface temperature, surface pressure, and eye height.
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

/// Internal: computes the water-vapor extinction coefficient.
pub fn kw(height_eye: f64, temp_s: f64, rh: f64) -> f64 {
    let mut wt = 0.031;
    wt *= 0.94 * (rh / 100.0) * (temp_s / 15.0).exp() * (-height_eye / SCALE_H_WATER).exp();
    wt
}

/// Internal: computes the ozone extinction coefficient.
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

/// Internal: computes the Rayleigh (molecular) extinction coefficient.
pub fn kr(alt_s: f64, height_eye: f64) -> f64 {
    let mut val = -alt_s - 12.0;
    val = val.clamp(0.0, 6.0);
    let changek = 1.0 - 0.166667 * val;
    let lambda = 0.55 + (changek - 1.0) * 0.04;
    0.1066 * (-height_eye / SCALE_H_RAYLEIGH).exp() * (lambda / 0.55_f64).powf(-4.0)
}

/// Internal: computes the aerosol extinction coefficient, optionally derived
/// from a supplied visual range `vr` (in km) instead of the empirical model.
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
        rh_clamped = rh_clamped.clamp(0.00000001, 99.99999999);
        let base = 0.1
            * (-height_eye / SCALE_H_AEROSOL).exp()
            * (1.0 - 0.32 / (rh_clamped / 100.0).ln()).powf(1.33)
            * (1.0 + 0.33 * sl * (sunra * DEGTORAD).sin());
        kaact = base * (lambda / 0.55_f64).powf(-1.3);
    }
    kaact
}

/// Internal: computes the total atmospheric extinction coefficient by
/// selecting/summing the aerosol (`ka`), water-vapor (`kw`), Rayleigh (`kr`),
/// and ozone (`koz`) components per `ext_type` (0=aerosol, 1=water, 2=Rayleigh,
/// 3=ozone, 4=all combined).
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

/// Internal: computes the relative optical air mass for a given apparent
/// altitude and surface pressure.
pub fn airmass(app_alt_o: f64, press: f64) -> f64 {
    let mut zend = (90.0 - app_alt_o) * DEGTORAD;
    if zend > std::f64::consts::FRAC_PI_2 {
        zend = std::f64::consts::FRAC_PI_2;
    }
    let airm = 1.0 / (zend.cos() + 0.025 * (-11.0 * zend.cos()).exp());
    press / 1013.0 * airm
}

/// Internal: computes the along-path column depth of an atmospheric layer of
/// given scale height at a given zenith distance (radians), scaled by pressure.
pub fn xext(scale_h: f64, zend: f64, press: f64) -> f64 {
    press
        / 1013.0
        / (zend.cos()
            + 0.01
                * (scale_h / 1000.0).sqrt()
                * (-30.0 / (scale_h / 1000.0).sqrt() * zend.cos()).exp())
}

/// Internal: computes the airmass path length through a spherical-shell
/// atmospheric layer of given scale height, using the Earth-radius geometry
/// term instead of a flat-atmosphere approximation.
pub fn xlay(scale_h: f64, zend: f64, press: f64) -> f64 {
    let a = zend.sin() / (1.0 + scale_h / RA);
    press / 1013.0 / (1.0 - a * a).sqrt()
}

/// Internal: computes the total magnitude correction due to atmospheric
/// extinction along the line of sight to the object, combining Rayleigh,
/// aerosol, ozone, and water-vapor path depths (`xext`/`xlay`) with their
/// respective extinction coefficients.
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

/// Internal: computes the critical visual acuity (minimum resolvable
/// separation, degrees) at background brightness `b` for a viewer with
/// Snellen fraction `sn`, choosing the scotopic or photopic model.
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

/// Internal: computes the pupil diameter (mm) as a function of observer age
/// and background brightness `b`.
pub fn pupil_dia(age: f64, b: f64) -> f64 {
    (0.534
        - 0.00211 * age
        - (0.236 - 0.00127 * age) * tanh_manual(0.4 * b.ln() / 10.0_f64.ln() - 2.2))
        * 10.0
}

/// Internal: computes the optical-aid brightness/contrast correction factor
/// (naked eye, binoculars, or telescope per `dobs`) used to scale either the
/// limiting-magnitude term (`type_factor == 0`) or the background-brightness
/// term (`type_factor != 0`).
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

/// Internal: computes the Moon's apparent visual magnitude at the given
/// Earth-Moon distance and phase angle (degrees).
pub fn moons_brightness(dist: f64, phasemoon: f64) -> f64 {
    -21.62
        + 5.0 * (dist / (RA / 1000.0)).ln() / LN10
        + 0.026 * phasemoon.abs()
        + 0.000000004 * phasemoon.powi(4)
}

/// Internal: computes the Moon's phase angle (degrees) from the Moon's and
/// Sun's topocentric altitude/azimuth.
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
/// Internal: computes the natural night-sky background brightness (nL),
/// including its slow 11-year (solar-cycle) modulation.
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

/// Internal: computes the sky background brightness (nL) contributed by
/// moonlight scattered near the line of sight to the object.
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

/// Internal: computes the twilight sky background brightness (nL)
/// contributed by sunlight scattered near the line of sight, for Sun
/// altitudes between civil/nautical twilight and day.
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

/// Internal: computes the daytime sky background brightness (nL)
/// contributed by sunlight scattered near the line of sight to the object.
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

/// Internal: clamps an artificial (city-light) sky-brightness contribution
/// to be non-negative.
pub fn bcity(value: f64) -> f64 {
    mymax(value, 0.0)
}

/// Internal: computes the total sky background brightness (nL) by combining
/// the twilight/daytime, moonlight, city-light, and night-sky components
/// (`btwi`/`bday`, `bm`, `bcity`, `bn`) applicable at the given Sun altitude.
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

// ── Object location & magnitude (c-ref-heliacal-vision.md §7) ─────

const LAPSE_RATE_DEFAULT: f64 = 0.0065;

fn topo_config(
    eph: &Ephemeris,
    dgeo: &[f64; 3],
    epheflag: CalcFlags,
) -> crate::config::EphemerisConfig {
    let mut config = eph.effective_config(epheflag, eph.config()).into_owned();
    config.topographic = Some(TopoPosition {
        longitude: dgeo[0],
        latitude: dgeo[1],
        altitude: dgeo[2],
    });
    config
}

/// Internal: crude calendar approximation of the Sun's right ascension
/// (degrees) for a given Julian day (UT), used only for extinction-model
/// seasonal terms.
pub fn sun_ra(jdn_days_ut: f64) -> f64 {
    let (_, imon, iday, _) = revjul(jdn_days_ut, CalendarType::Gregorian);
    normalize_degrees((imon as f64 + (iday as f64 - 1.0) / 30.4 - 3.69) * 30.0)
}

/// Internal: computes one of several angle values for `object_name` at
/// `jd_ut` — selected via `angle` (0=apparent altitude, 1=azimuth, 2=right
/// ascension, 3=ecliptic longitude, 4=apparent-refracted altitude,
/// 5=declination, 6=ecliptic latitude, 7=geocentric altitude).
#[allow(clippy::collapsible_else_if)]
pub fn object_loc(
    eph: &Ephemeris,
    jd_ut: f64,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
    object_name: &str,
    angle: i32,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<f64, Error> {
    // C checks Angle<5 for TOPOCTR BEFORE rewriting Angle==7→0,
    // so angle=7 (geocentric altitude) does NOT get TOPOCTR.
    let mut iflag = CalcFlags::EQUATORIAL | (epheflag & calc::EPHMASK);
    if !helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        iflag |= CalcFlags::NONUT | CalcFlags::TRUEPOS;
    }
    if angle < 5 {
        iflag |= CalcFlags::TOPOCTR;
    }
    let angle = if angle == 7 { 0 } else { angle };

    let tjd_tt =
        jd_ut + crate::deltat::calc_deltat(jd_ut, &eph.effective_config(epheflag, eph.config()));
    let planet = object_to_body(object_name);

    let x = if let Some(body) = planet {
        if iflag.contains(CalcFlags::TOPOCTR) {
            let config = topo_config(eph, dgeo, epheflag);
            eph.calc_with_config(tjd_tt, body, iflag, &config)?
        } else {
            eph.calc(tjd_tt, body, iflag)?
        }
    } else {
        if iflag.contains(CalcFlags::TOPOCTR) {
            let config = topo_config(eph, dgeo, epheflag);
            eph.fixstar2_with_config(object_name, tjd_tt, iflag, &config)?
                .1
        } else {
            eph.fixstar2(object_name, tjd_tt, iflag)?.1
        }
    };

    if angle == 2 || angle == 5 {
        Ok(x.data[1])
    } else if angle == 3 || angle == 6 {
        Ok(x.data[0])
    } else {
        let xin = [x.data[0], x.data[1]];
        let xaz = eph.azalt(
            jd_ut,
            AzAltDir::EquToHor,
            [dgeo[0], dgeo[1], dgeo[2]],
            datm[0],
            datm[1],
            LAPSE_RATE_DEFAULT,
            xin,
        );
        if angle == 0 {
            Ok(xaz[1])
        } else if angle == 4 {
            // C's argument-order quirk: datm[0] (pressure) in TempE slot, datm[1] (temp) in PresE slot
            Ok(app_alt_from_topo_alt(xaz[1], datm[0], datm[1], helflag))
        } else {
            // angle == 1: azimuth, flipped 180° from swe_azalt's south-origin convention
            let mut azi = xaz[0] + 180.0;
            if azi >= 360.0 {
                azi -= 360.0;
            }
            Ok(azi)
        }
    }
}

/// Internal: computes the object's raw azimuth, true altitude, apparent
/// altitude, and the Cartesian unit vector of the apparent-altitude
/// direction (matching C's `azalt_cart`, `dret[0..5]`).
pub fn azalt_cart(
    eph: &Ephemeris,
    jd_ut: f64,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<[f64; 6], Error> {
    let mut iflag = CalcFlags::EQUATORIAL | CalcFlags::TOPOCTR | (epheflag & calc::EPHMASK);
    if !helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        iflag |= CalcFlags::NONUT | CalcFlags::TRUEPOS;
    }

    let tjd_tt =
        jd_ut + crate::deltat::calc_deltat(jd_ut, &eph.effective_config(epheflag, eph.config()));
    let planet = object_to_body(object_name);
    let config = topo_config(eph, dgeo, epheflag);

    let x = if let Some(body) = planet {
        eph.calc_with_config(tjd_tt, body, iflag, &config)?
    } else {
        eph.fixstar2_with_config(object_name, tjd_tt, iflag, &config)?
            .1
    };

    let xin = [x.data[0], x.data[1]];
    let xaz = eph.azalt(
        jd_ut,
        AzAltDir::EquToHor,
        [dgeo[0], dgeo[1], dgeo[2]],
        datm[0],
        datm[1],
        LAPSE_RATE_DEFAULT,
        xin,
    );

    // C feeds degree-valued az/alt directly to swi_polcart (which calls cos/sin
    // expecting radians) — intentional reproduction for golden parity.
    let cart = polar_to_cartesian([xaz[0], xaz[2], 1.0]);

    Ok([xaz[0], xaz[1], xaz[2], cart[0], cart[1], cart[2]])
}

/// Internal: computes the apparent visual magnitude of a planet, the Moon,
/// or a fixed star at the given Julian day (UT).
pub fn magnitude(
    eph: &Ephemeris,
    jd_ut: f64,
    dgeo: &[f64; 3],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<f64, Error> {
    let planet = object_to_body(object_name);

    if let Some(body) = planet {
        let mut iflag = CalcFlags::TOPOCTR | CalcFlags::EQUATORIAL | (epheflag & calc::EPHMASK);
        if !helflag.contains(HeliacalFlags::HIGH_PRECISION) {
            iflag |= CalcFlags::NONUT | CalcFlags::TRUEPOS;
        }
        let config = topo_config(eph, dgeo, epheflag);
        let (pheno, _) = crate::phenomena::pheno_ut_with_config(eph, jd_ut, body, iflag, &config)?;
        Ok(pheno.apparent_magnitude)
    } else {
        let (_, mag) = eph.fixstar2_mag(object_name)?;
        Ok(mag)
    }
}

// ── Rise/set wrappers (c-ref-heliacal-vision.md §2/§7) ───────────

const HELIACAL_AU: f64 = 1.49597870691e+11;
const SUN_RADIUS_M: f64 = 696000000.0;
const MOON_RADIUS_M: f64 = 1737000.0;
const LAT_THRESHOLD_FAST: f64 = 63.0;

/// Internal: computes the rise or set Julian day (UT) of `ipl` using the
/// fast direct-search algorithm (valid for `|latitude| < 63°`); see
/// `my_rise_trans` for the dispatcher that falls back to the general
/// rise/transit search at high latitudes.
pub fn calc_rise_and_set(
    eph: &Ephemeris,
    tjd_start: f64,
    ipl: Body,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
    eventflag: RiseSetFlags,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<f64, Error> {
    let mut iflag = epheflag & calc::EPHMASK;
    if !helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        iflag |= CalcFlags::NONUT | CalcFlags::TRUEPOS;
    }

    let tjd0 = tjd_start;
    let geopos = [dgeo[0], dgeo[1], dgeo[2]];

    // Step 2: local-noon estimate
    let mut tjdnoon = (tjd0 as i64) as f64 - dgeo[0] / 15.0 / 24.0;

    // Step 3: compute Sun and object RA at tjd0
    let sun_iflag = iflag | CalcFlags::EQUATORIAL;
    let xs = eph.calc_ut(tjd0, Body::Sun, sun_iflag)?.data;
    let xx_init = eph.calc_ut(tjd0, ipl, sun_iflag)?.data;
    tjdnoon -= normalize_degrees(xs[0] - xx_init[0]) / 360.0;

    // Step 4: is the object currently above/below horizon?
    let xin = [xx_init[0], xx_init[1]];
    let xaz = eph.azalt(
        tjd0,
        AzAltDir::EquToHor,
        geopos,
        datm[0],
        datm[1],
        LAPSE_RATE_DEFAULT,
        xin,
    );
    let above = xaz[2] > 0.0;

    // Step 5: day-anchoring
    let is_rise = eventflag.contains(RiseSetFlags::RISE);
    if is_rise {
        if above {
            while tjdnoon <= tjd0 + 0.5 {
                tjdnoon += 1.0;
            }
            while tjdnoon > tjd0 + 1.5 {
                tjdnoon -= 1.0;
            }
        } else {
            while tjdnoon < tjd0 {
                tjdnoon += 1.0;
            }
            while tjdnoon > tjd0 + 1.0 {
                tjdnoon -= 1.0;
            }
        }
    } else {
        if above {
            while tjdnoon < tjd0 - 0.5 {
                tjdnoon += 1.0;
            }
            while tjdnoon > tjd0 + 0.5 {
                tjdnoon -= 1.0;
            }
        } else {
            while tjdnoon < tjd0 {
                tjdnoon += 1.0;
            }
            while tjdnoon > tjd0 + 1.0 {
                tjdnoon -= 1.0;
            }
        }
    }

    // Step 6: recompute position at tjdnoon for declination
    let xx_noon = eph.calc_ut(tjdnoon, ipl, sun_iflag)?.data;

    // Step 7: disc radius
    let rdi = if eventflag.contains(RiseSetFlags::DISC_CENTER) {
        0.0
    } else if ipl == Body::Sun {
        (SUN_RADIUS_M / HELIACAL_AU / xx_noon[2]).asin() / DEGTORAD
    } else if ipl == Body::Moon {
        (MOON_RADIUS_M / HELIACAL_AU / xx_noon[2]).asin() / DEGTORAD
    } else {
        0.0
    };

    // Step 8: target altitude
    let rh = -(34.5 / 60.0 + rdi);

    // Step 9: semi-diurnal arc
    let sda = (-dgeo[1].to_radians().tan() * xx_noon[1].to_radians().tan())
        .acos()
        .to_degrees();

    // Step 10: initial estimate
    let mut tjdrise = if is_rise {
        tjdnoon - sda / 360.0
    } else {
        tjdnoon + sda / 360.0
    };

    // Step 11: refinement loop (2 iterations)
    let config = topo_config(eph, dgeo, epheflag);
    let mut refine_iflag = (epheflag & calc::EPHMASK) | CalcFlags::SPEED | CalcFlags::EQUATORIAL;
    if ipl == Body::Moon {
        refine_iflag |= CalcFlags::TOPOCTR;
    }
    if !helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        refine_iflag |= CalcFlags::NONUT | CalcFlags::TRUEPOS;
    }

    let dfac: f64 = 1.0 / 365.25;

    for _ in 0..2 {
        let xx = if refine_iflag.contains(CalcFlags::TOPOCTR) {
            eph.calc_ut_with_config(tjdrise, ipl, refine_iflag, &config)?
                .data
        } else {
            eph.calc_ut(tjdrise, ipl, refine_iflag)?.data
        };

        let xin1 = [xx[0], xx[1]];
        let xaz1 = eph.azalt(
            tjdrise,
            AzAltDir::EquToHor,
            geopos,
            datm[0],
            datm[1],
            LAPSE_RATE_DEFAULT,
            xin1,
        );

        // Back-propagate RA/decl by dfac using speed
        let xin2 = [xx[0] - xx[3] * dfac, xx[1] - xx[4] * dfac];
        let xaz2 = eph.azalt(
            tjdrise - dfac,
            AzAltDir::EquToHor,
            geopos,
            datm[0],
            datm[1],
            LAPSE_RATE_DEFAULT,
            xin2,
        );

        // Secant-style update
        let dalt = xaz1[1] - xaz2[1];
        if dalt != 0.0 {
            tjdrise -= (xaz1[1] - rh) / dalt * dfac;
        }
    }

    Ok(tjdrise)
}

/// Internal: rise/set dispatcher for a planet or (via `starname`) a fixed
/// star; uses the fast `calc_rise_and_set` path for recognized planets at
/// `|latitude| < 63°`, otherwise falls back to the general rise/transit search.
pub fn my_rise_trans(
    eph: &Ephemeris,
    tjd: f64,
    ipl: Body,
    starname: Option<&str>,
    eventtype: RiseSetFlags,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
) -> Result<f64, Error> {
    // If starname is provided, resolve to body first
    let resolved_ipl = if let Some(name) = starname {
        if !name.is_empty() {
            if let Some(body) = object_to_body(name) {
                body
            } else {
                // Fixed star: always use full rise_trans
                let rsmi = eventtype;
                let atpress = datm[0];
                let attemp = datm[1];
                let result = eph.rise_trans(
                    tjd,
                    Body::Sun, // unused for stars
                    Some(name),
                    epheflag & calc::EPHMASK,
                    rsmi,
                    [dgeo[0], dgeo[1], dgeo[2]],
                    atpress,
                    attemp,
                )?;
                return Ok(result.time);
            }
        } else {
            ipl
        }
    } else {
        ipl
    };

    // Fast path: recognized planet AND |lat| < 63°
    if dgeo[1].abs() < LAT_THRESHOLD_FAST {
        calc_rise_and_set(
            eph,
            tjd,
            resolved_ipl,
            dgeo,
            datm,
            eventtype,
            epheflag,
            helflag,
        )
    } else {
        let rsmi = eventtype;
        let atpress = datm[0];
        let attemp = datm[1];
        let result = eph.rise_trans(
            tjd,
            resolved_ipl,
            None,
            epheflag & calc::EPHMASK,
            rsmi,
            [dgeo[0], dgeo[1], dgeo[2]],
            atpress,
            attemp,
        )?;
        Ok(result.time)
    }
}

/// Internal: resolves `object_name` to a body or star and computes its rise
/// or set time (`rim == 0` forces disc-center timing).
pub fn rise_set(
    eph: &Ephemeris,
    jdn_days_ut: f64,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
    object_name: &str,
    rs_event: RiseSetFlags,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    rim: i32,
) -> Result<f64, Error> {
    let mut eventflags = rs_event;
    if rim == 0 {
        eventflags |= RiseSetFlags::DISC_CENTER;
    }

    let planet = object_to_body(object_name);
    if let Some(body) = planet {
        my_rise_trans(
            eph,
            jdn_days_ut,
            body,
            Some(""),
            eventflags,
            epheflag,
            helflag,
            dgeo,
            datm,
        )
    } else {
        my_rise_trans(
            eph,
            jdn_days_ut,
            Body::Sun, // placeholder, unused for stars
            Some(object_name),
            eventflags,
            epheflag,
            helflag,
            dgeo,
            datm,
        )
    }
}

// ── VisLimMagn & swe_vis_limit_mag (c-ref-heliacal-vision.md §8) ──

/// Output of `vis_limit_mag`: the visual limiting magnitude at the observer's
/// location and time, alongside the geometry used to compute it. Port of
/// `swe_vis_limit_mag`'s `dret[0..7]` output array (swehel.c).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct VisLimitResult {
    /// Limiting visual magnitude the object would need to reach to be just
    /// visible (from `vis_lim_magn`). C `dret[0]`.
    pub limiting_magnitude: f64,
    /// Object's true (topocentric) altitude, degrees. C `dret[1]`.
    pub altitude_object: f64,
    /// Object's azimuth, degrees. C `dret[2]`.
    pub azimuth_object: f64,
    /// Sun's altitude, degrees (`-90` under `VISLIM_DARK`). C `dret[3]`.
    pub altitude_sun: f64,
    /// Sun's azimuth, degrees (`0` under `VISLIM_DARK`). C `dret[4]`.
    pub azimuth_sun: f64,
    /// Moon's altitude, degrees (`-90` if the object is the Moon, or under
    /// `VISLIM_DARK`/`VISLIM_NOMOON`). C `dret[5]`.
    pub altitude_moon: f64,
    /// Moon's azimuth, degrees (`0` under the same conditions as
    /// `altitude_moon`). C `dret[6]`.
    pub azimuth_moon: f64,
    /// Object's actual apparent visual magnitude (from `magnitude`). C `dret[7]`.
    pub magnitude_object: f64,
    /// Scotopic/mixed vision-mode flags reported alongside the limiting
    /// magnitude.
    pub vision: VisLimFlags,
    /// `true` if the object is below the local horizon, in which case
    /// `limiting_magnitude` is the `-100.0` sentinel and the other fields are
    /// left at their default (zero) values.
    pub below_horizon: bool,
}

/// Computes the visual limiting magnitude and scotopic/mixed vision flags
/// for an object at the given topocentric geometry, from the combined sky
/// background brightness (`bsky`) and extinction (`deltam`). Port of
/// `VisLimMagn` (swehel.c).
#[allow(clippy::approx_constant, clippy::impossible_comparisons)]
pub fn vis_lim_magn(
    dobs: &[f64; 6],
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
) -> (f64, VisLimFlags) {
    let bsk = bsky(
        alt_o,
        azi_o,
        alt_m,
        azi_m,
        jdn_days_ut,
        alt_s,
        azi_s,
        sunra,
        lat,
        height_eye,
        datm,
        helflag,
    );
    let k_x = deltam(alt_o, alt_s, sunra, lat, height_eye, datm, helflag);
    let corr_factor1 = optic_factor(bsk, k_x, dobs, false, 1, helflag);
    let corr_factor2 = optic_factor(bsk, k_x, dobs, false, 0, helflag);

    let mut is_scotopic = bsk < 1645.0;
    if helflag.contains(HeliacalFlags::VISLIM_PHOTOPIC) {
        is_scotopic = false;
    }
    if helflag.contains(HeliacalFlags::VISLIM_SCOTOPIC) {
        is_scotopic = true;
    }

    let (c1, c2) = if is_scotopic {
        (1.5848931924611e-10, 0.012589254117942)
    } else {
        (4.4668359215096e-9, 1.2589254117942e-6)
    };

    let mut scotopic_flag = if is_scotopic {
        VisLimFlags::SCOTOPIC
    } else {
        VisLimFlags::empty()
    };

    if BNIGHT * BNIGHT_FACTOR > bsk && BNIGHT / BNIGHT_FACTOR < bsk {
        scotopic_flag |= VisLimFlags::MIXED;
    }

    let bsk_corr = bsk * corr_factor1;
    let th = c1 * (1.0 + (c2 * bsk_corr).sqrt()).powi(2) * corr_factor2;
    let log10_val: f64 = 2.302585092994;
    let mag = -16.57 - 2.5 * (th.ln() / log10_val);

    (mag, scotopic_flag)
}

/// Computes the visual limiting magnitude for `object_name` at `tjd_ut` and
/// the given geographic/atmospheric/observer parameters, including the
/// object's, Sun's, and Moon's positions used in the computation. Port of
/// `swe_vis_limit_mag` (swehel.c).
pub fn vis_limit_mag(
    eph: &Ephemeris,
    tjd_ut: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<VisLimitResult, Error> {
    let name = tolower_string_star(object_name);

    if object_to_body(&name) == Some(Body::Sun) {
        return Err(Error::CError("object name is Sun for vis_limit_mag".into()));
    }

    let sunra = sun_ra(tjd_ut);
    default_heliacal_parameters(datm, dgeo, dobs, helflag);

    let alt_o = object_loc(eph, tjd_ut, dgeo, datm, &name, 0, epheflag, helflag)?;

    if alt_o < 0.0 {
        return Ok(VisLimitResult {
            limiting_magnitude: -100.0,
            altitude_object: 0.0,
            azimuth_object: 0.0,
            altitude_sun: 0.0,
            azimuth_sun: 0.0,
            altitude_moon: 0.0,
            azimuth_moon: 0.0,
            magnitude_object: 0.0,
            vision: VisLimFlags::empty(),
            below_horizon: true,
        });
    }

    let azi_o = object_loc(eph, tjd_ut, dgeo, datm, &name, 1, epheflag, helflag)?;

    let (alt_s, azi_s) = if helflag.contains(HeliacalFlags::VISLIM_DARK) {
        (-90.0, 0.0)
    } else {
        let a = object_loc(eph, tjd_ut, dgeo, datm, "sun", 0, epheflag, helflag)?;
        let z = object_loc(eph, tjd_ut, dgeo, datm, "sun", 1, epheflag, helflag)?;
        (a, z)
    };

    let is_moon_object = name.starts_with("moon");
    let (alt_m, azi_m) = if is_moon_object
        || helflag.contains(HeliacalFlags::VISLIM_DARK)
        || helflag.contains(HeliacalFlags::VISLIM_NOMOON)
    {
        (-90.0, 0.0)
    } else {
        let a = object_loc(eph, tjd_ut, dgeo, datm, "moon", 0, epheflag, helflag)?;
        let z = object_loc(eph, tjd_ut, dgeo, datm, "moon", 1, epheflag, helflag)?;
        (a, z)
    };

    let (lim_mag, scotopic_flag) = vis_lim_magn(
        dobs, alt_o, azi_o, alt_m, azi_m, tjd_ut, alt_s, azi_s, sunra, dgeo[1], dgeo[2], datm,
        helflag,
    );

    let mag_obj = magnitude(eph, tjd_ut, dgeo, &name, epheflag, helflag)?;

    Ok(VisLimitResult {
        limiting_magnitude: lim_mag,
        altitude_object: alt_o,
        azimuth_object: azi_o,
        altitude_sun: alt_s,
        azimuth_sun: azi_s,
        altitude_moon: alt_m,
        azimuth_moon: azi_m,
        magnitude_object: mag_obj,
        vision: scotopic_flag,
        below_horizon: false,
    })
}

// ── TopoArcVisionis: bisection for arcus visionis ──────────────────

/// Bisection search for the Sun-depression angle at which an object of known
/// magnitude becomes exactly visible. Port of `TopoArcVisionis` (swehel.c:1562-1599).
pub fn topo_arc_visionis(
    magn: f64,
    dobs: &[f64; 6],
    alt_o: f64,
    azi_o: f64,
    alt_m: f64,
    azi_m: f64,
    jdn_days_ut: f64,
    azi_s: f64,
    sunra: f64,
    lat: f64,
    height_eye: f64,
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> f64 {
    let epsilon = 0.001;
    let mut xl = 45.0_f64;
    let mut xr = 0.0_f64;

    let yl_vlm = vis_lim_magn(
        dobs,
        alt_o,
        azi_o,
        alt_m,
        azi_m,
        jdn_days_ut,
        alt_o - xl,
        azi_s,
        sunra,
        lat,
        height_eye,
        datm,
        helflag,
    );
    let mut yl = magn - yl_vlm.0;

    let yr_vlm = vis_lim_magn(
        dobs,
        alt_o,
        azi_o,
        alt_m,
        azi_m,
        jdn_days_ut,
        alt_o - xr,
        azi_s,
        sunra,
        lat,
        height_eye,
        datm,
        helflag,
    );
    let yr = magn - yr_vlm.0;

    let mut xm;
    if yl * yr <= 0.0 {
        while (xr - xl).abs() > epsilon {
            xm = (xr + xl) / 2.0;
            let alt_si = alt_o - xm;
            let ym_vlm = vis_lim_magn(
                dobs,
                alt_o,
                azi_o,
                alt_m,
                azi_m,
                jdn_days_ut,
                alt_si,
                azi_s,
                sunra,
                lat,
                height_eye,
                datm,
                helflag,
            );
            let ym = magn - ym_vlm.0;
            if yl * ym > 0.0 {
                xl = xm;
                yl = ym;
            } else {
                xr = xm;
            }
        }
        xm = (xr + xl) / 2.0;
    } else {
        xm = 99.0;
    }
    if xm < alt_o {
        xm = alt_o;
    }
    xm
}

/// Public wrapper for `topo_arc_visionis`. Port of `swe_topo_arcus_visionis`
/// (swehel.c:1601-1610). All geometry is caller-supplied; no ephemeris lookups
/// beyond the crude `sun_ra` calendar approximation.
#[allow(clippy::too_many_arguments)]
pub fn topo_arcus_visionis(
    tjd_ut: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    helflag: HeliacalFlags,
    mag: f64,
    azi_obj: f64,
    alt_obj: f64,
    azi_sun: f64,
    azi_moon: f64,
    alt_moon: f64,
) -> Result<f64, Error> {
    let sunra = sun_ra(tjd_ut);
    default_heliacal_parameters(datm, dgeo, dobs, helflag);
    Ok(topo_arc_visionis(
        mag, dobs, alt_obj, azi_obj, alt_moon, azi_moon, tjd_ut, azi_sun, sunra, dgeo[1], dgeo[2],
        datm, helflag,
    ))
}

// ── HeliacalAngle: optimum-altitude / arcus-visionis search ────────

/// Output of `heliacal_angle`: the optimal object altitude, arcus visionis,
/// and implied Sun altitude for first/last visibility.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HeliacalAngleResult {
    /// Object's altitude at the optimum (degrees).
    pub optimal_altitude: f64,
    /// Arcus visionis at the optimum (degrees).
    pub arcus_visionis: f64,
    /// Implied Sun altitude at the optimum (degrees) = optimal_altitude - arcus_visionis.
    pub sun_altitude_diff: f64,
}

/// 2-D search for the optimal object-altitude / arcus-visionis pair. Port of
/// `HeliacalAngle` (swehel.c:1636-1693). For each candidate altitude 2..20°,
/// calls `topo_arc_visionis`; then refines the minimum via bisection.
pub fn heliacal_angle_core(
    magn: f64,
    dobs: &[f64; 6],
    azi_o: f64,
    alt_m: f64,
    azi_m: f64,
    jdn_days_ut: f64,
    azi_s: f64,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
    helflag: HeliacalFlags,
) -> HeliacalAngleResult {
    let sunra = sun_ra(jdn_days_ut);
    let lat = dgeo[1];
    let height_eye = dgeo[2];

    // Coarse scan: integer altitudes 2..=20
    let mut xmin = 0.0_f64;
    let mut ymin = 10000.0_f64;
    for ix in 2..=20 {
        let x = ix as f64;
        let arc = topo_arc_visionis(
            magn,
            dobs,
            x,
            azi_o,
            alt_m,
            azi_m,
            jdn_days_ut,
            azi_s,
            sunra,
            lat,
            height_eye,
            datm,
            helflag,
        );
        if arc < ymin {
            ymin = arc;
            xmin = x;
        }
    }

    // Bracket the coarse minimum by ±1°
    let mut xl = xmin - 1.0;
    let mut xr = xmin + 1.0;
    let mut _yl = topo_arc_visionis(
        magn,
        dobs,
        xl,
        azi_o,
        alt_m,
        azi_m,
        jdn_days_ut,
        azi_s,
        sunra,
        lat,
        height_eye,
        datm,
        helflag,
    );
    let mut _yr = topo_arc_visionis(
        magn,
        dobs,
        xr,
        azi_o,
        alt_m,
        azi_m,
        jdn_days_ut,
        azi_s,
        sunra,
        lat,
        height_eye,
        datm,
        helflag,
    );

    // Minimum-finding bisection (one-sided finite-difference slope check)
    while (xr - xl).abs() > 0.1 {
        let xm = (xr + xl) / 2.0;
        let delta_x = 0.025;
        let xmd = xm + delta_x;
        let ym = topo_arc_visionis(
            magn,
            dobs,
            xm,
            azi_o,
            alt_m,
            azi_m,
            jdn_days_ut,
            azi_s,
            sunra,
            lat,
            height_eye,
            datm,
            helflag,
        );
        let ymd = topo_arc_visionis(
            magn,
            dobs,
            xmd,
            azi_o,
            alt_m,
            azi_m,
            jdn_days_ut,
            azi_s,
            sunra,
            lat,
            height_eye,
            datm,
            helflag,
        );
        if ym >= ymd {
            xl = xm;
            _yl = ym;
        } else {
            xr = xm;
            _yr = ym;
        }
    }

    let xm = (xr + xl) / 2.0;
    // C averages the last iteration's Yl/Yr, not a fresh evaluation at Xm
    let ym = (_yr + _yl) / 2.0;

    HeliacalAngleResult {
        optimal_altitude: xm,
        arcus_visionis: ym,
        sun_altitude_diff: xm - ym,
    }
}

/// Public wrapper for `heliacal_angle_core`. Port of `swe_heliacal_angle`
/// (swehel.c:1695-1705). Validates observer altitude, applies defaults, delegates.
#[allow(clippy::too_many_arguments)]
pub fn heliacal_angle(
    tjd_ut: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    helflag: HeliacalFlags,
    mag: f64,
    azi_obj: f64,
    azi_sun: f64,
    azi_moon: f64,
    alt_moon: f64,
) -> Result<HeliacalAngleResult, Error> {
    if !(crate::constants::RISE_SET_GEOALT_MIN..=crate::constants::RISE_SET_GEOALT_MAX)
        .contains(&dgeo[2])
    {
        return Err(Error::CError(format!(
            "location for heliacal events must be between {} and {} m above sea",
            crate::constants::RISE_SET_GEOALT_MIN,
            crate::constants::RISE_SET_GEOALT_MAX,
        )));
    }
    default_heliacal_parameters(datm, dgeo, dobs, helflag);
    Ok(heliacal_angle_core(
        mag, dobs, azi_obj, alt_moon, azi_moon, tjd_ut, azi_sun, dgeo, datm, helflag,
    ))
}

// ── Moon crescent geometry (swehel.c:1715-1745) ───────────────────

const AVG_RADIUS_MOON: f64 = 15.541 / 60.0;

fn width_moon(alt_o: f64, azi_o: f64, alt_s: f64, azi_s: f64, parallax: f64) -> f64 {
    let geo_alt_o = alt_o + parallax;
    0.27245
        * parallax
        * (1.0 + (geo_alt_o * DEGTORAD).sin() * (parallax * DEGTORAD).sin())
        * (1.0 - ((alt_s - geo_alt_o) * DEGTORAD).cos() * ((azi_s - azi_o) * DEGTORAD).cos())
}

fn length_moon(w: f64, mut diamoon: f64) -> f64 {
    if diamoon == 0.0 {
        diamoon = AVG_RADIUS_MOON * 2.0;
    }
    let wi = w * 60.0;
    let d = diamoon * 60.0;
    (d - 0.3 * (d + wi) / 2.0 / wi) / 60.0
}

fn q_yallop(w: f64, geo_arcv_act: f64) -> f64 {
    let wi = w * 60.0;
    (geo_arcv_act - (11.8371 - 6.3226 * wi + 0.7319 * wi * wi - 0.1018 * wi * wi * wi)) / 10.0
}

fn yallop_grade(q: f64) -> i32 {
    if q > 0.216 {
        return 1;
    }
    if q < 0.216 && q > -0.014 {
        return 2;
    }
    if q < -0.014 && q > -0.16 {
        return 3;
    }
    if q < -0.16 && q > -0.232 {
        return 4;
    }
    if q < -0.232 && q > -0.293 {
        return 5;
    }
    if q < -0.293 {
        return 6;
    }
    0
}

// ── Interpolation helpers (swehel.c:1753-1810) ────────────────────

fn crossing(a: f64, b: f64, c: f64, d: f64) -> f64 {
    (c - a) / ((b - a) - (d - c))
}

fn x2min(a: f64, b: f64, c: f64) -> f64 {
    let term = a + c - 2.0 * b;
    if term == 0.0 {
        return 0.0;
    }
    -(a - c) / 2.0 / term
}

fn funct2(a: f64, b: f64, c: f64, x: f64) -> f64 {
    (a + c - 2.0 * b) / 2.0 * x * x + (a - c) / 2.0 * x + b
}

// ── DeterTAV (swehel.c:1759-1783) ─────────────────────────────────

fn deter_tav(
    eph: &Ephemeris,
    dobs: &[f64; 6],
    jdn_days_ut: f64,
    dgeo: &[f64; 3],
    datm: &[f64; 4],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<f64, Error> {
    let sunra = sun_ra(jdn_days_ut);
    let magn = magnitude(eph, jdn_days_ut, dgeo, object_name, epheflag, helflag)?;
    let alt_o = object_loc(
        eph,
        jdn_days_ut,
        dgeo,
        datm,
        object_name,
        0,
        epheflag,
        helflag,
    )?;
    let azi_o = object_loc(
        eph,
        jdn_days_ut,
        dgeo,
        datm,
        object_name,
        1,
        epheflag,
        helflag,
    )?;

    let (alt_m, azi_m) = if object_name.starts_with("moon") {
        (-90.0, 0.0)
    } else {
        let a = object_loc(eph, jdn_days_ut, dgeo, datm, "moon", 0, epheflag, helflag)?;
        let z = object_loc(eph, jdn_days_ut, dgeo, datm, "moon", 1, epheflag, helflag)?;
        (a, z)
    };

    let azi_s = object_loc(eph, jdn_days_ut, dgeo, datm, "sun", 1, epheflag, helflag)?;

    let result = topo_arc_visionis(
        magn,
        dobs,
        alt_o,
        azi_o,
        alt_m,
        azi_m,
        jdn_days_ut,
        azi_s,
        sunra,
        dgeo[1],
        dgeo[2],
        datm,
        helflag,
    );
    Ok(result)
}

// ── HeliacalPheno output struct ───────────────────────────────────

const TJD_INVALID: f64 = 99999999.0;
const MAX_TRY_HOURS: f64 = 4.0;
const TIME_STEP_DEFAULT: f64 = 1.0;
const LOCAL_MIN_STEP: f64 = 8.0;

/// Detailed heliacal-phenomena report for an object at a specific instant.
/// Port of `swe_heliacal_pheno_ut`'s 28-element `dret[]` output array
/// (swehel.c:1862-2074).
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HeliacalPheno {
    /// Object's topocentric altitude, degrees. C `dret[0]`.
    pub tc_altitude: f64,
    /// Object's topocentric apparent (refracted) altitude, degrees. C `dret[1]`.
    pub tc_apparent_altitude: f64,
    /// Object's geocentric altitude, degrees. C `dret[2]`.
    pub gc_altitude: f64,
    /// Object's azimuth, degrees. C `dret[3]`.
    pub azimuth_object: f64,
    /// Sun's topocentric altitude, degrees. C `dret[4]`.
    pub tc_sun_altitude: f64,
    /// Sun's azimuth, degrees. C `dret[5]`.
    pub sun_azimuth: f64,
    /// Actual topocentric arcus visionis (object altitude minus Sun
    /// altitude), degrees. C `dret[6]`.
    pub tav_act: f64,
    /// Actual (parallax-corrected) arcus visionis, degrees. C `dret[7]`.
    pub arcv_act: f64,
    /// Actual azimuth difference between Sun and object, degrees. C `dret[8]`.
    pub daz_act: f64,
    /// Actual longitude difference between Sun and object along the great
    /// circle, degrees. C `dret[9]`.
    pub arcl_act: f64,
    /// Extinction coefficient at the Sun's altitude. C `dret[10]`.
    pub kact: f64,
    /// Smallest topocentric arcus visionis found in the visibility-window
    /// search. C `dret[11]`.
    pub min_tav: f64,
    /// First time (Julian day, UT) the object is possibly visible. C `dret[12]`.
    pub t_first_vr: f64,
    /// Best (optimum) time (Julian day, UT) the object is visible. C `dret[13]`.
    pub t_best_vr: f64,
    /// Last time (Julian day, UT) the object is possibly visible. C `dret[14]`.
    pub t_last_vr: f64,
    /// Optimum time of visibility (Julian day, UT) per Yallop's criterion,
    /// Moon only. C `dret[15]`.
    pub t_best_yallop: f64,
    /// Crescent width of the Moon, Moon only. C `dret[16]`.
    pub w_moon: f64,
    /// Yallop's q-test value, Moon only. C `dret[17]`.
    pub q_yallop: f64,
    /// Yallop visibility-grade classification derived from `q_yallop`, Moon
    /// only. C `dret[18]`.
    pub q_crit: f64,
    /// Object's parallax, degrees. C `dret[19]`.
    pub par_o: f64,
    /// Object's actual apparent visual magnitude. C `dret[20]`.
    pub magn_o: f64,
    /// Object's rise/set time (Julian day, UT). C `dret[21]`.
    pub rise_o: f64,
    /// Sun's rise/set time (Julian day, UT). C `dret[22]`.
    pub rise_s: f64,
    /// Time lag between the object's and the Sun's rise/set, days. C `dret[23]`.
    pub lag: f64,
    /// Duration of the visibility window, days. C `dret[24]`.
    pub t_vis_vr: f64,
    /// Crescent length of the Moon, Moon only. C `dret[25]`.
    pub l_moon: f64,
    /// Elongation of the object from the Sun, degrees. C `dret[26]`.
    pub elongation: f64,
    /// Illuminated fraction of the object's disc, percent. C `dret[27]`.
    pub illumination: f64,
}

impl HeliacalPheno {
    /// Flattens the fields into the 28-element array matching C's `dret[]`
    /// layout, in field-declaration order.
    pub fn as_array(&self) -> [f64; 28] {
        [
            self.tc_altitude,
            self.tc_apparent_altitude,
            self.gc_altitude,
            self.azimuth_object,
            self.tc_sun_altitude,
            self.sun_azimuth,
            self.tav_act,
            self.arcv_act,
            self.daz_act,
            self.arcl_act,
            self.kact,
            self.min_tav,
            self.t_first_vr,
            self.t_best_vr,
            self.t_last_vr,
            self.t_best_yallop,
            self.w_moon,
            self.q_yallop,
            self.q_crit,
            self.par_o,
            self.magn_o,
            self.rise_o,
            self.rise_s,
            self.lag,
            self.t_vis_vr,
            self.l_moon,
            self.elongation,
            self.illumination,
        ]
    }
}

// ── swe_heliacal_pheno_ut (swehel.c:1862-2074) ───────────────────

/// Computes the full `HeliacalPheno` report for `object_name` at `tjd_ut`
/// and the given event type (morning/evening first/last). Port of
/// `swe_heliacal_pheno_ut` (swehel.c:1862-2074).
pub fn heliacal_pheno_ut(
    eph: &Ephemeris,
    tjd_ut: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    event: HeliacalEventType,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<HeliacalPheno, Error> {
    // Step 1: altitude bound check
    if !(crate::constants::RISE_SET_GEOALT_MIN..=crate::constants::RISE_SET_GEOALT_MAX)
        .contains(&dgeo[2])
    {
        return Err(Error::CError(format!(
            "location for heliacal events must be between {} and {} m above sea",
            crate::constants::RISE_SET_GEOALT_MIN,
            crate::constants::RISE_SET_GEOALT_MAX,
        )));
    }

    // Step 3: sun RA
    let sunra = sun_ra(tjd_ut);

    let type_event = event as i32;
    if !(1..=4).contains(&type_event) {
        return Err(Error::CError(format!(
            "heliacal_pheno_ut: TypeEvent {} not supported (must be 1-4)",
            type_event,
        )));
    }

    // Step 4: name normalization
    let name = tolower_string_star(object_name);

    // Step 5: default parameters
    default_heliacal_parameters(datm, dgeo, dobs, helflag);

    let iflag = epheflag & calc::EPHMASK;

    // Step 7: Sun and object geometry
    let azi_s = object_loc(eph, tjd_ut, dgeo, datm, "sun", 1, iflag, helflag)?;
    let alt_s = object_loc(eph, tjd_ut, dgeo, datm, "sun", 0, iflag, helflag)?;
    let azi_o = object_loc(eph, tjd_ut, dgeo, datm, &name, 1, iflag, helflag)?;
    let alt_o = object_loc(eph, tjd_ut, dgeo, datm, &name, 0, iflag, helflag)?;
    let geo_alt_o = object_loc(eph, tjd_ut, dgeo, datm, &name, 7, iflag, helflag)?;

    // Step 8: derived quantities
    let app_alt_o = app_alt_from_topo_alt(alt_o, datm[1], datm[0], helflag);
    let daz_act = azi_s - azi_o;
    let tav_act = alt_o - alt_s;
    let par_o = geo_alt_o - alt_o;
    let magn_o = magnitude(eph, tjd_ut, dgeo, &name, iflag, helflag)?;
    let arcv_act = tav_act + par_o;
    let arcl_act = (arcv_act * DEGTORAD).cos() * (daz_act * DEGTORAD).cos();
    let arcl_act = arcl_act.acos() / DEGTORAD;

    // Step 9: elongation & illumination
    let planet = object_to_body(&name);
    let (elong, illum) = if let Some(body) = planet {
        let pheno_iflag = iflag | CalcFlags::TOPOCTR | CalcFlags::EQUATORIAL;
        let config = topo_config(eph, dgeo, epheflag);
        let (pheno, _) =
            crate::phenomena::pheno_ut_with_config(eph, tjd_ut, body, pheno_iflag, &config)?;
        (pheno.elongation, pheno.phase * 100.0)
    } else {
        (arcl_act, 100.0)
    };

    // Step 10: extinction coefficient
    let kact = kt(alt_s, sunra, dgeo[1], dgeo[2], datm[1], datm[2], datm[3], 4);

    // Step 12: Moon-only Yallop crescent block
    let mut w_moon_val = 0.0_f64;
    let mut q_yal = 0.0_f64;
    let mut q_crit_val = 0.0_f64;
    let mut l_moon_val = 0.0_f64;

    let is_moon = planet == Some(Body::Moon);

    if is_moon {
        w_moon_val = width_moon(alt_o, azi_o, alt_s, azi_s, par_o);
        l_moon_val = length_moon(w_moon_val, 0.0);
        q_yal = q_yallop(w_moon_val, arcv_act);
        q_crit_val = yallop_grade(q_yal) as f64;
    }

    // Step 13: rise/set of Sun and object
    let rs: RiseSetFlags = if type_event == 1 || type_event == 4 {
        RiseSetFlags::RISE | RiseSetFlags::DISC_CENTER
    } else {
        RiseSetFlags::SET | RiseSetFlags::DISC_CENTER
    };

    let rise_set_s = rise_set(
        eph,
        tjd_ut - 4.0 / 24.0,
        dgeo,
        datm,
        "sun",
        rs,
        iflag,
        helflag,
        0,
    )?;

    let mut tb_yallop = TJD_INVALID;
    let lag;
    let mut norise_o = false;
    let mut rise_set_o = 0.0;

    match rise_set(
        eph,
        tjd_ut - 4.0 / 24.0,
        dgeo,
        datm,
        &name,
        rs,
        iflag,
        helflag,
        0,
    ) {
        Ok(t) => {
            rise_set_o = t;
            lag = rise_set_o - rise_set_s;
            if is_moon {
                tb_yallop = (rise_set_o * 4.0 + rise_set_s * 5.0) / 9.0;
            }
        }
        Err(Error::CircumpolarBody) => {
            lag = 0.0;
            norise_o = true;
        }
        Err(e) => return Err(e),
    }

    // Step 14: early-exit guard for TypeEvent 3/4 and non-inner planets
    let is_star = planet.is_none();
    let is_mars_or_beyond = match planet {
        Some(b) => b.to_raw_id() >= Body::Mars.to_raw_id(),
        None => false,
    };

    let mut t_first_vr;
    let mut t_best_vr;
    let mut t_last_vr;
    let t_vis_vr;
    let mut min_tav;

    if (type_event == 3 || type_event == 4) && (is_star || is_mars_or_beyond) {
        t_first_vr = TJD_INVALID;
        t_best_vr = TJD_INVALID;
        t_last_vr = TJD_INVALID;
        t_vis_vr = 0.0;
        min_tav = 0.0;
    } else {
        // Step 15: visibility-window search loop
        let mut min_tav_act = 199.0_f64;
        let mut delta_alt = 0.0_f64;
        #[allow(unused_assignments)]
        let mut oldest_min_tav = 0.0_f64;
        let mut ta = 0.0_f64;
        let mut tc = 0.0_f64;
        t_best_vr = 0.0;
        min_tav = 0.0;

        let mut time_step = -TIME_STEP_DEFAULT / 24.0 / 60.0;
        let is_setting = type_event == 2 || type_event == 3;
        if is_setting {
            time_step = -time_step;
        }
        let mut time_pointer = rise_set_s - time_step;

        #[allow(unused_assignments)]
        let mut min_tav_oud = 0.0_f64;
        #[allow(unused_assignments)]
        let mut delta_alt_oud = 0.0_f64;

        loop {
            time_pointer += time_step;
            oldest_min_tav = min_tav_oud;
            min_tav_oud = min_tav_act;
            delta_alt_oud = delta_alt;

            let alt_s2 = object_loc(eph, time_pointer, dgeo, datm, "sun", 0, iflag, helflag)?;
            let alt_o2 = object_loc(eph, time_pointer, dgeo, datm, &name, 0, iflag, helflag)?;
            delta_alt = alt_o2 - alt_s2;
            min_tav_act = deter_tav(eph, dobs, time_pointer, dgeo, datm, &name, iflag, helflag)?;

            // Local-minimum detection
            if min_tav_oud < min_tav_act && t_best_vr == 0.0 {
                let mut time_check = time_pointer + sgn(time_step) * LOCAL_MIN_STEP / 24.0 / 60.0;
                if !norise_o {
                    if time_step > 0.0 {
                        time_check = mymin(time_check, rise_set_o);
                    } else {
                        time_check = mymax(time_check, rise_set_o);
                    }
                }
                let local_min_check =
                    deter_tav(eph, dobs, time_check, dgeo, datm, &name, iflag, helflag)?;
                if local_min_check > min_tav_act {
                    let extrax = x2min(min_tav_act, min_tav_oud, oldest_min_tav);
                    t_best_vr = time_pointer - (1.0 - extrax) * time_step;
                    min_tav = funct2(min_tav_act, min_tav_oud, oldest_min_tav, extrax);
                }
            }

            // Visibility-start crossing
            if delta_alt > min_tav_act && tc == 0.0 && t_best_vr == 0.0 {
                let crosspoint = crossing(delta_alt_oud, delta_alt, min_tav_oud, min_tav_act);
                tc = time_pointer - time_step * (1.0 - crosspoint);
            }

            // Visibility-end crossing
            if delta_alt < min_tav_act && ta == 0.0 && tc != 0.0 {
                let crosspoint = crossing(delta_alt_oud, delta_alt, min_tav_oud, min_tav_act);
                ta = time_pointer - time_step * (1.0 - crosspoint);
            }

            // Loop termination
            let within_range = (time_pointer - rise_set_s).abs() <= MAX_TRY_HOURS / 24.0;
            let found_end = ta != 0.0;
            let early_exit_34 = t_best_vr != 0.0
                && (type_event == 3 || type_event == 4)
                && !name.starts_with("moon")
                && !name.starts_with("venus")
                && !name.starts_with("mercury");

            if !within_range || found_end || early_exit_34 {
                break;
            }
        }

        // Step 16: post-loop assembly
        if is_setting {
            t_first_vr = tc;
            t_last_vr = ta;
        } else {
            t_first_vr = ta;
            t_last_vr = tc;
        }

        if t_first_vr == 0.0 && t_last_vr == 0.0 {
            if !is_setting {
                t_first_vr = t_best_vr - 0.000001;
            } else {
                t_last_vr = t_best_vr + 0.000001;
            }
        }

        if !norise_o {
            if !is_setting {
                t_first_vr = mymax(t_first_vr, rise_set_o);
            } else {
                t_last_vr = mymin(t_last_vr, rise_set_o);
            }
        }

        t_vis_vr = if t_last_vr != 0.0 && t_first_vr != 0.0 {
            t_last_vr - t_first_vr
        } else {
            TJD_INVALID
        };

        if t_last_vr == 0.0 {
            t_last_vr = TJD_INVALID;
        }
        if t_best_vr == 0.0 {
            t_best_vr = TJD_INVALID;
        }
        if t_first_vr == 0.0 {
            t_first_vr = TJD_INVALID;
        }
    }

    // Step 18: output
    Ok(HeliacalPheno {
        tc_altitude: alt_o,
        tc_apparent_altitude: app_alt_o,
        gc_altitude: geo_alt_o,
        azimuth_object: azi_o,
        tc_sun_altitude: alt_s,
        sun_azimuth: azi_s,
        tav_act,
        arcv_act,
        daz_act,
        arcl_act,
        kact,
        min_tav,
        t_first_vr,
        t_best_vr,
        t_last_vr,
        t_best_yallop: tb_yallop,
        w_moon: w_moon_val,
        q_yallop: q_yal,
        q_crit: q_crit_val,
        par_o,
        magn_o,
        rise_o: rise_set_o,
        rise_s: rise_set_s,
        lag,
        t_vis_vr,
        l_moon: l_moon_val,
        elongation: elong,
        illumination: illum,
    })
}

// --- event search infrastructure (c-ref-heliacal-search.md §1–§4) ---

// §1: Synodic periods & conjunction table

fn get_synodic_period(body: Body) -> f64 {
    match body {
        Body::Moon => 29.530588853,
        Body::Mercury => 115.8775,
        Body::Venus => 583.9214,
        Body::Mars => 779.9361,
        Body::Jupiter => 398.8840,
        Body::Saturn => 378.0919,
        Body::Uranus => 369.6560,
        Body::Neptune => 367.4867,
        Body::Pluto => 366.7207,
        _ => 366.0,
    }
}

const TCON: [f64; 18] = [
    0.0, 0.0, // Sun (placeholder)
    2451550.0, 2451550.0, // Moon
    2451604.0, 2451670.0, // Mercury
    2451980.0, 2452280.0, // Venus
    2451727.0, 2452074.0, // Mars
    2451673.0, 2451877.0, // Jupiter
    2451675.0, 2451868.0, // Saturn
    2451581.0, 2451768.0, // Uranus
    2451568.0, 2451753.0, // Neptune
];

/// Internal: finds the Julian day of the conjunction (or opposition, for
/// acronychal events) between `ipl` and the Sun nearest after `tjd_start`,
/// seeded from a lookup table (`TCON`) of known conjunction epochs.
pub fn find_conjunct_sun(
    eph: &Ephemeris,
    tjd_start: f64,
    ipl: Body,
    epheflag: CalcFlags,
    type_event: i32,
) -> Result<f64, Error> {
    let ipl_idx = ipl.to_raw_id() as usize;

    // Pluto (index 9) has no TCON row — latent C out-of-bounds bug.
    if ipl_idx > 8 {
        return Err(Error::CError(format!(
            "find_conjunct_sun: body {} not supported (no TCON entry beyond Neptune)",
            ipl_idx,
        )));
    }

    let daspect = if ipl_idx >= 4 && type_event >= 3 {
        180.0
    } else {
        0.0
    };

    let i = (type_event - 1) / 2 + (ipl_idx as i32) * 2;
    let tjd0 = TCON[i as usize];
    let dsynperiod = get_synodic_period(ipl);

    let mut tjdcon = tjd0 + (((tjd_start - tjd0) / dsynperiod).floor() + 1.0) * dsynperiod;

    let efl = epheflag & calc::EPHMASK | CalcFlags::SPEED;
    let mut ds = 100.0_f64;
    let mut niter = 0;
    while ds > 0.5 {
        let x = eph.calc(tjdcon, ipl, efl)?;
        let xs = eph.calc(tjdcon, Body::Sun, efl)?;
        ds = normalize_degrees(x.data[0] - xs.data[0] - daspect);
        if ds > 180.0 {
            ds -= 360.0;
        }
        tjdcon -= ds / (x.data[3] - xs.data[3]);
        niter += 1;
        if niter > 10_000 {
            return Err(Error::CError(
                "find_conjunct_sun: convergence failure".into(),
            ));
        }
    }
    Ok(tjdcon)
}

// §2: Oblique-ascension machinery

fn get_asc_obl(
    eph: &Ephemeris,
    tjd: f64,
    ipl: Body,
    starname: Option<&str>,
    epheflag: CalcFlags,
    dgeo: &[f64; 3],
    desc_obl: bool,
) -> Result<f64, Error> {
    let efl = (epheflag & calc::EPHMASK) | CalcFlags::EQUATORIAL;

    let (ra, decl) = match starname.filter(|s| !s.is_empty()) {
        Some(star) => {
            let (_, r) = eph.fixstar2(star, tjd, efl)?;
            (r.data[0], r.data[1])
        }
        None => {
            let r = eph.calc(tjd, ipl, efl)?;
            (r.data[0], r.data[1])
        }
    };

    let adp = (dgeo[1] * DEGTORAD).tan() * (decl * DEGTORAD).tan();
    if adp.abs() > 1.0 {
        return Err(Error::CircumpolarBody);
    }
    let adp_deg = adp.asin() / DEGTORAD;

    let daop = if desc_obl { ra + adp_deg } else { ra - adp_deg };
    Ok(normalize_degrees(daop))
}

fn get_asc_obl_diff(
    eph: &Ephemeris,
    tjd: f64,
    ipl: Body,
    starname: Option<&str>,
    epheflag: CalcFlags,
    dgeo: &[f64; 3],
    desc_obl: bool,
    is_acronychal: bool,
) -> Result<f64, Error> {
    let aosun = get_asc_obl(eph, tjd, Body::Sun, None, epheflag, dgeo, desc_obl)?;

    let body_desc = if is_acronychal { !desc_obl } else { desc_obl };
    let aopl = get_asc_obl(eph, tjd, ipl, starname, epheflag, dgeo, body_desc)?;

    let mut dsunpl = normalize_degrees(aosun - aopl);
    if is_acronychal {
        dsunpl = normalize_degrees(dsunpl - 180.0);
    }
    if dsunpl > 180.0 {
        dsunpl -= 360.0;
    }
    Ok(dsunpl)
}

/// Internal: finds the Julian day at which the oblique-ascension difference
/// between `ipl` (or `starname`) and the Sun crosses zero for the given
/// event type, via coarse forward search then bisection.
pub fn get_asc_obl_with_sun(
    eph: &Ephemeris,
    tjd_start: f64,
    ipl: Body,
    starname: Option<&str>,
    epheflag: CalcFlags,
    evtyp: i32,
    dperiod: f64,
    dgeo: &[f64; 3],
) -> Result<f64, Error> {
    let desc_obl = evtyp == 2 || evtyp == 3 || evtyp == 5;
    let is_acronychal = evtyp == 5 || evtyp == 6;

    let mut retro = evtyp == 1 || evtyp == 2;
    if is_acronychal && ipl != Body::Moon {
        retro = true;
    }

    let efl = epheflag & calc::EPHMASK;
    let mut tjd = tjd_start;
    let mut dsunpl = get_asc_obl_diff(eph, tjd, ipl, starname, efl, dgeo, desc_obl, is_acronychal)?;
    let mut dsunpl_save = -999999999.0_f64;

    // Coarse forward search with 10-day steps
    let mut i = 0;
    while dsunpl_save == -999999999.0
        || (dsunpl.abs() + dsunpl_save.abs() > 180.0)
        || (retro && !(dsunpl_save < 0.0 && dsunpl >= 0.0))
        || !(retro || dsunpl_save >= 0.0 && dsunpl < 0.0)
    {
        dsunpl_save = dsunpl;
        tjd += 10.0;
        if dperiod > 0.0 && tjd - tjd_start > dperiod {
            return Err(Error::CircumpolarBody);
        }
        dsunpl = get_asc_obl_diff(eph, tjd, ipl, starname, efl, dgeo, desc_obl, is_acronychal)?;
        i += 1;
        if i > 5000 {
            return Err(Error::CError("loop in get_asc_obl_with_sun() (1)".into()));
        }
    }

    // Bisection with 20-day initial bracket
    let mut daystep = 20.0_f64;
    let mut tjd_lo = tjd - daystep;
    daystep /= 2.0;
    tjd = tjd_lo + daystep;
    let mut dsunpl_test =
        get_asc_obl_diff(eph, tjd, ipl, starname, efl, dgeo, desc_obl, is_acronychal)?;

    i = 0;
    while dsunpl.abs() > 0.00001 {
        if dsunpl_save * dsunpl_test >= 0.0 {
            dsunpl_save = dsunpl_test;
            tjd_lo = tjd;
        } else {
            dsunpl = dsunpl_test;
        }
        daystep /= 2.0;
        tjd = tjd_lo + daystep;
        dsunpl_test =
            get_asc_obl_diff(eph, tjd, ipl, starname, efl, dgeo, desc_obl, is_acronychal)?;
        i += 1;
        if i > 5000 {
            return Err(Error::CError("loop in get_asc_obl_with_sun() (2)".into()));
        }
    }

    Ok(tjd)
}

// §3: Day-level search

/// Internal: searches day-by-day (then minute-by-minute) from `tjd` for the
/// Julian day (UT) on which `object_name`'s heliacal event of `type_event`
/// first occurs.
pub fn get_heliacal_day(
    eph: &Ephemeris,
    tjd: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<f64, Error> {
    let (is_rise_or_set, direct_day, direct_time) = match type_event {
        1 => (RiseSetFlags::RISE, 1.0_f64, -1.0_f64),
        2 => (RiseSetFlags::SET, -1.0, 1.0),
        3 => (RiseSetFlags::SET, 1.0, 1.0),
        4 => (RiseSetFlags::RISE, -1.0, -1.0),
        _ => {
            return Err(Error::CError(format!(
                "get_heliacal_day: invalid TypeEvent {}",
                type_event,
            )));
        }
    };

    let ipl = object_to_body(object_name);
    let ipl_id = ipl.map(|b| b.to_raw_id()).unwrap_or(-1);

    let (ndays, tjd_adj, mut daystep, tfac): (i32, f64, f64, f64) = match ipl {
        Some(Body::Moon) => (16, tjd, 1.0, 1.0),
        Some(Body::Mercury) => (60, tjd, 5.0, 5.0),
        Some(Body::Venus) => {
            let ds = if type_event >= 3 { 15.0 } else { 5.0 };
            let tf = if type_event >= 3 { 3.0 } else { 1.0 };
            (300, tjd - 30.0 * direct_day, ds, tf)
        }
        Some(Body::Mars) => (400, tjd, 15.0, 5.0),
        Some(Body::Saturn) => (300, tjd, 20.0, 5.0),
        None => {
            // Fixed star: get magnitude to adjust tfac
            let mag = magnitude(eph, tjd, dgeo, object_name, epheflag, helflag)?;
            let tf = if mag < 0.0 { 3.0 } else { 10.0 };
            (300, tjd, 15.0, tf)
        }
        _ => (300, tjd, 15.0, 3.0), // Jupiter, Uranus, Neptune, Pluto, etc.
    };

    let tend = tjd_adj + (ndays as f64) * direct_day;
    let mut retval_old: i32 = -2; // sentinel
    let div = 1440.0_f64; // minutes per day

    let mut tday = tjd_adj;
    let mut iter_count = 0;
    loop {
        if (direct_day > 0.0 && tday >= tend) || (direct_day < 0.0 && tday <= tend) {
            break;
        }

        if iter_count > 0 {
            tday -= 0.3 * direct_day;
        }

        // Sun rise/set for this day
        let sun_result = my_rise_trans(
            eph,
            tday,
            Body::Sun,
            Some(""),
            is_rise_or_set,
            epheflag,
            helflag,
            dgeo,
            datm,
        );
        match sun_result {
            Err(Error::CircumpolarBody) => {
                retval_old = -2;
                tday += daystep * direct_day;
                iter_count += 1;
                continue;
            }
            Err(e) => return Err(e),
            Ok(_) => {}
        }
        let mut tret = sun_result.unwrap();

        let vlm = vis_limit_mag(eph, tret, dgeo, datm, dobs, object_name, epheflag, helflag)?;

        let retval: i32 = if vlm.below_horizon { -2 } else { 0 };

        // Daystep-shrink-on-first-appearance
        if retval_old == -2 && retval >= 0 && daystep > 1.0 {
            retval_old = retval;
            tday -= daystep * direct_day;
            daystep = if ipl_id >= 4 || ipl_id == -1 {
                5.0
            } else {
                1.0
            };
            tday += daystep * direct_day;
            iter_count += 1;
            continue;
        }
        retval_old = retval;

        if retval == -2 {
            tday += daystep * direct_day;
            iter_count += 1;
            continue;
        }

        // Minute-level refinement within the day
        let mut visible_at_sunsetrise = true;
        let mut vd = vlm.limiting_magnitude - vlm.magnitude_object;
        let mut minute_iters = 0;
        while !vlm.below_horizon && vd < 0.0 {
            visible_at_sunsetrise = false;
            let step = if vd < -1.0 {
                5.0 / div * direct_time * tfac
            } else if vd < -0.5 {
                2.0 / div * direct_time * tfac
            } else if vd < -0.1 {
                1.0 / div * direct_time * tfac
            } else {
                1.0 / div * direct_time // no tfac for finest bracket
            };
            tret += step;
            let vlm2 = vis_limit_mag(eph, tret, dgeo, datm, dobs, object_name, epheflag, helflag)?;
            if vlm2.below_horizon {
                break;
            }
            vd = vlm2.limiting_magnitude - vlm2.magnitude_object;
            minute_iters += 1;
            if minute_iters > 10_000 {
                break;
            }
        }

        // Sunset/sunrise-instant edge nudge
        if visible_at_sunsetrise {
            for _ in 0..10 {
                let vlm2 = vis_limit_mag(
                    eph,
                    tret + 1.0 / div * direct_time,
                    dgeo,
                    datm,
                    dobs,
                    object_name,
                    epheflag,
                    helflag,
                )?;
                if !vlm2.below_horizon {
                    let vd_new = vlm2.limiting_magnitude - vlm2.magnitude_object;
                    if vd_new > vd {
                        vd = vd_new;
                        tret += 1.0 / div * direct_time;
                    }
                }
            }
        }

        // Acceptance
        if vd > 0.0 {
            if (ipl_id >= 4 || ipl_id == -1) && daystep > 1.0 {
                tday -= daystep * direct_day;
                daystep = 1.0;
                tday += daystep * direct_day;
                iter_count += 1;
                continue;
            }
            return Ok(tret);
        }

        tday += daystep * direct_day;
        iter_count += 1;
    }

    Err(Error::CError("heliacal event does not happen".into()))
}

/// Internal: searches from `tjd` for the Julian day (UT) of `object_name`'s
/// acronychal (cosmical) rising or setting of `type_event`, by iterating
/// rise/set plus the dark- and no-moon visibility-limit boundaries until
/// they converge.
pub fn get_acronychal_day(
    eph: &Ephemeris,
    tjd: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<f64, Error> {
    let helflag2 = helflag | HeliacalFlags::VISLIM_PHOTOPIC;

    let (is_rise_or_set, direct) = if type_event == 3 || type_event == 5 {
        (RiseSetFlags::RISE, -1.0_f64)
    } else {
        (RiseSetFlags::SET, 1.0_f64)
    };

    let ipl = object_to_body(object_name);
    let body = ipl.unwrap_or(Body::Sun);
    let star = if ipl.is_none() {
        Some(object_name)
    } else {
        Some("")
    };

    let mut tjd_work = tjd;
    let mut dtret = 999.0_f64;
    let mut tret;
    let mut niter = 0;
    while dtret.abs() > 0.5 / 1440.0 {
        tjd_work += 0.7 * direct;
        if direct < 0.0 {
            tjd_work -= 1.0;
        }

        tjd_work = my_rise_trans(
            eph,
            tjd_work,
            body,
            star,
            is_rise_or_set,
            epheflag,
            helflag2,
            dgeo,
            datm,
        )?;

        // Walk until visible
        let mut vlm = vis_limit_mag(
            eph,
            tjd_work,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag2,
        )?;
        let mut walk_iters = 0;
        while !vlm.below_horizon && vlm.limiting_magnitude < vlm.magnitude_object {
            tjd_work += 10.0 / 1440.0 * (-direct);
            vlm = vis_limit_mag(
                eph,
                tjd_work,
                dgeo,
                datm,
                dobs,
                object_name,
                epheflag,
                helflag2,
            )?;
            walk_iters += 1;
            if walk_iters > 10_000 {
                break;
            }
        }

        let mut tret_dark = tjd_work;
        time_limit_invisible(
            eph,
            tjd_work,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag2 | HeliacalFlags::VISLIM_DARK,
            direct,
            &mut tret_dark,
        )?;

        tret = tjd_work;
        time_limit_invisible(
            eph,
            tjd_work,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag2 | HeliacalFlags::VISLIM_NOMOON,
            direct,
            &mut tret,
        )?;

        dtret = (tret - tret_dark).abs();
        tjd_work = tret;
        niter += 1;
        if niter > 10_000 {
            return Err(Error::CError(
                "get_acronychal_day: convergence failure".into(),
            ));
        }
    }

    Ok(tjd_work)
}

// §4: Visibility timing

/// Internal: hill-climbs forward and backward from `tjd` to find the time of
/// maximum visibility margin (limiting magnitude minus actual magnitude) for
/// `object_name`; returns the optimum time and whether it is uncertain due
/// to a scotopic/photopic transition.
pub fn time_optimum_visibility(
    eph: &Ephemeris,
    tjd: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<(f64, bool), Error> {
    let vlm = vis_limit_mag(eph, tjd, dgeo, datm, dobs, object_name, epheflag, helflag)?;
    let mut retval_sv = vlm.vision;
    let mut phot_scot_sv = vlm.vision.contains(VisLimFlags::SCOTOPIC);

    let mut t1 = tjd;
    let mut t2 = tjd;
    let mut vl1 = -1.0_f64;
    let mut vl2 = -1.0_f64;

    // Forward hill-climb
    let mut d = 100.0 / 86400.0;
    for _ in 0..3 {
        t1 += d;
        let mut t_has_changed = false;
        loop {
            let vlm2 = vis_limit_mag(
                eph,
                t1 - d,
                dgeo,
                datm,
                dobs,
                object_name,
                epheflag,
                helflag,
            )?;
            if vlm2.below_horizon {
                break;
            }
            let margin = vlm2.limiting_magnitude - vlm2.magnitude_object;
            if margin <= 0.0 || margin <= vl1 {
                break;
            }
            t1 -= d;
            vl1 = margin;
            retval_sv = vlm2.vision;
            phot_scot_sv = vlm2.vision.contains(VisLimFlags::SCOTOPIC);
            t_has_changed = true;
        }
        if !t_has_changed {
            t1 -= d;
        }
        d /= 10.0;
    }

    // Backward hill-climb
    d = 100.0 / 86400.0;
    for _ in 0..3 {
        t2 -= d;
        let mut t_has_changed = false;
        loop {
            let vlm2 = vis_limit_mag(
                eph,
                t2 + d,
                dgeo,
                datm,
                dobs,
                object_name,
                epheflag,
                helflag,
            )?;
            if vlm2.below_horizon {
                break;
            }
            let margin = vlm2.limiting_magnitude - vlm2.magnitude_object;
            if margin <= 0.0 || margin <= vl2 {
                break;
            }
            t2 += d;
            vl2 = margin;
            retval_sv = vlm2.vision;
            phot_scot_sv = vlm2.vision.contains(VisLimFlags::SCOTOPIC);
            t_has_changed = true;
        }
        if !t_has_changed {
            t2 += d;
        }
        d /= 10.0;
    }

    let result_tjd = if vl2 > vl1 { t2 } else { t1 };

    // Scotopic/photopic transition check
    let vlm_final = vis_limit_mag(
        eph,
        result_tjd,
        dgeo,
        datm,
        dobs,
        object_name,
        epheflag,
        helflag,
    )?;
    if !vlm_final.below_horizon {
        let phot_scot_final = vlm_final.vision.contains(VisLimFlags::SCOTOPIC);
        if phot_scot_final != phot_scot_sv {
            return Ok((result_tjd, true)); // uncertain due to vision mode change
        }
        if retval_sv.contains(VisLimFlags::MIXED) {
            return Ok((result_tjd, true));
        }
    }

    Ok((result_tjd, false))
}

/// Internal: walks from `tjd` in direction `direct` to find the boundary
/// time (written to `*tret`) beyond which `object_name` becomes invisible;
/// returns whether the result is uncertain due to a scotopic/photopic
/// transition.
pub fn time_limit_invisible(
    eph: &Ephemeris,
    tjd: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    direct: f64,
    tret: &mut f64,
) -> Result<bool, Error> {
    let is_moon = object_name == "moon";
    let d0 = if is_moon {
        100.0 / 86400.0 * 10.0
    } else {
        100.0 / 86400.0
    };
    let ncnt = if is_moon { 4 } else { 3 };

    let vlm = vis_limit_mag(eph, tjd, dgeo, datm, dobs, object_name, epheflag, helflag)?;
    let mut retval_sv = vlm.vision;
    let phot_scot_sv = vlm.vision.contains(VisLimFlags::SCOTOPIC);

    let mut tjd_work = tjd;
    let mut d = d0;
    for _ in 0..ncnt {
        loop {
            let vlm2 = vis_limit_mag(
                eph,
                tjd_work + d * direct,
                dgeo,
                datm,
                dobs,
                object_name,
                epheflag,
                helflag,
            )?;
            if vlm2.below_horizon {
                break;
            }
            if vlm2.limiting_magnitude <= vlm2.magnitude_object {
                break;
            }
            tjd_work += d * direct;
            retval_sv = vlm2.vision;
        }
        d /= 10.0;
    }

    *tret = tjd_work;

    // Scotopic/photopic transition check
    let phot_scot_final = retval_sv.contains(VisLimFlags::SCOTOPIC);
    if phot_scot_final != phot_scot_sv {
        return Ok(true); // uncertain
    }
    if retval_sv.contains(VisLimFlags::MIXED) {
        return Ok(true);
    }

    Ok(false)
}

/// Internal: given a day already known to contain the heliacal event
/// (`tday`), computes the `[start_visible, optimum, end_visible]` triple via
/// `time_optimum_visibility` and `time_limit_invisible`, reordering the
/// boundaries chronologically for evening events. Port of
/// `get_heliacal_details` (swehel.c:3107-3161).
pub fn get_heliacal_details(
    eph: &Ephemeris,
    tday: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<[f64; 3], Error> {
    let mut dret = [0.0_f64; 3];

    // dret[1] = optimum
    let (optimum, _optimum_uncertain) =
        time_optimum_visibility(eph, tday, dgeo, datm, dobs, object_name, epheflag, helflag)?;
    dret[1] = optimum;

    // dret[0] = first boundary
    let mut direct: f64 = if type_event == 1 || type_event == 4 {
        -1.0
    } else {
        1.0
    };
    let mut _limit_1_uncertain = false;
    let mut t0 = tday;
    _limit_1_uncertain = time_limit_invisible(
        eph,
        tday,
        dgeo,
        datm,
        dobs,
        object_name,
        epheflag,
        helflag,
        direct,
        &mut t0,
    )?;
    dret[0] = t0;

    // dret[2] = second boundary (seeded from optimum)
    direct *= -1.0;
    let mut t2 = dret[1];
    let mut _limit_2_uncertain = false;
    _limit_2_uncertain = time_limit_invisible(
        eph,
        dret[1],
        dgeo,
        datm,
        dobs,
        object_name,
        epheflag,
        helflag,
        direct,
        &mut t2,
    )?;
    dret[2] = t2;

    // Reorder for evening events
    if type_event == 2 || type_event == 3 {
        dret.swap(0, 2);
    }

    Ok(dret)
}

// ── §5 Event drivers — vis_lim path ───────────────────────────────

/// Output of `heliacal_ut`: the three Julian-day (UT) instants bracketing
/// a heliacal event.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HeliacalEvent {
    /// Beginning of visibility (or, for arc_vis path, the single event instant).
    pub start_visible: f64,
    /// Optimum visibility (0.0 if arc_vis path or NO_DETAILS).
    pub optimum_visibility: f64,
    /// End of visibility (0.0 if arc_vis path or NO_DETAILS).
    pub end_visible: f64,
}

const MAX_COUNT_SYNPER: i32 = 5;
const MAX_COUNT_SYNPER_MAX: i32 = 1_000_000;

fn heliacal_ut_vis_lim(
    eph: &Ephemeris,
    tjd_start: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<[f64; 3], Error> {
    let mut dret = [0.0_f64; 3];
    let ipl = object_to_body(object_name);

    let tjd = if matches!(ipl, Some(Body::Mercury)) {
        tjd_start - 30.0
    } else {
        tjd_start - 50.0
    };

    let is_heliacal = matches!(ipl, Some(Body::Mercury) | Some(Body::Venus)) || type_event <= 2;

    let star_arg: Option<&str> = if ipl.is_none() {
        Some(object_name)
    } else {
        None
    };
    let body_arg = ipl.unwrap_or(Body::Sun);

    let tjd_seed;
    if is_heliacal {
        if let Some(body) = ipl {
            // Planet: conjunction search
            tjd_seed = find_conjunct_sun(eph, tjd, body, epheflag, type_event)?;
        } else {
            // Fixed star: oblique ascension search
            tjd_seed = get_asc_obl_with_sun(
                eph, tjd, body_arg, star_arg, epheflag, type_event, 0.0, dgeo,
            )?;
        }

        let tday = get_heliacal_day(
            eph,
            tjd_seed,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag,
            type_event,
        )?;
        dret[0] = tday;
    } else {
        // Acronychal branch: outer planets/stars with TypeEvent 3 or 4
        tjd_seed = get_asc_obl_with_sun(
            eph, tjd, body_arg, star_arg, epheflag, type_event, 0.0, dgeo,
        )?;

        let tday = get_acronychal_day(
            eph,
            tjd_seed,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag,
            type_event,
        )?;
        dret[0] = tday;
    }

    // Details refinement (unless NO_DETAILS)
    if !helflag.contains(HeliacalFlags::NO_DETAILS) && is_heliacal {
        dret = get_heliacal_details(
            eph,
            dret[0],
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag,
            type_event,
        )?;
    }
    // Acronychal branch: dret[1]/dret[2] stay 0.0 (details refinement is dead code in C)

    Ok(dret)
}

fn moon_event_vis_lim(
    eph: &Ephemeris,
    tjd_start: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<[f64; 3], Error> {
    if type_event == 1 || type_event == 2 {
        return Err(Error::CError(
            "the moon has no morning first or evening last".into(),
        ));
    }

    // HIGH_PRECISION stripped for day search (unlike planet path)
    let helflag2 = helflag & !HeliacalFlags::HIGH_PRECISION;

    let tjd = tjd_start - 30.0;
    let tjd_conj = find_conjunct_sun(eph, tjd, Body::Moon, epheflag, type_event)?;

    let mut tjd_work = get_heliacal_day(
        eph, tjd_conj, dgeo, datm, dobs, "moon", epheflag, helflag2, type_event,
    )?;

    let mut dret = [0.0_f64; 3];
    dret[0] = tjd_work;

    // Optimum (full helflag)
    let (optimum, _) =
        time_optimum_visibility(eph, tjd_work, dgeo, datm, dobs, "moon", epheflag, helflag)?;
    dret[1] = optimum;
    tjd_work = optimum;

    // End boundary
    let direct: f64 = if type_event == 4 { -1.0 } else { 1.0 };
    let mut t_end = tjd_work;
    time_limit_invisible(
        eph, tjd_work, dgeo, datm, dobs, "moon", epheflag, helflag, direct, &mut t_end,
    )?;
    dret[2] = t_end;

    // Start boundary (seeded from optimum, opposite direction)
    let mut t_start = dret[1];
    time_limit_invisible(
        eph,
        dret[1],
        dgeo,
        datm,
        dobs,
        "moon",
        epheflag,
        helflag,
        -direct,
        &mut t_start,
    )?;
    dret[0] = t_start;

    // Sunset/sunrise clamp
    if type_event == 3 {
        // Evening first: clamp start to sunset
        let trise = my_rise_trans(
            eph,
            t_start,
            Body::Sun,
            Some(""),
            RiseSetFlags::SET,
            epheflag,
            helflag,
            dgeo,
            datm,
        )?;
        if trise < dret[1] {
            dret[0] = trise;
        }
    } else {
        // Morning last (TypeEvent==4): clamp end to sunrise
        let trise = my_rise_trans(
            eph,
            dret[1],
            Body::Sun,
            Some(""),
            RiseSetFlags::RISE,
            epheflag,
            helflag,
            dgeo,
            datm,
        )?;
        if dret[0] > trise {
            dret[0] = trise;
        }
    }

    // Reorder for TypeEvent==4
    if type_event == 4 {
        dret.swap(0, 2);
    }

    Ok(dret)
}

fn moon_event_jd_ut(
    eph: &Ephemeris,
    tjd_start: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<[f64; 3], Error> {
    if helflag.intersects(HeliacalFlags::AVKIND) {
        moon_event_arc_vis(
            eph, tjd_start, dgeo, datm, dobs, epheflag, helflag, type_event,
        )
    } else {
        moon_event_vis_lim(
            eph, tjd_start, dgeo, datm, dobs, epheflag, helflag, type_event,
        )
    }
}

// ── §6 Event drivers — arc_vis path ──────────────────────────────

fn moon_event_arc_vis(
    eph: &Ephemeris,
    tjd_start: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<[f64; 3], Error> {
    let avkind = helflag & HeliacalFlags::AVKIND;
    let avkind = if avkind.is_empty() {
        HeliacalFlags::AVKIND_VR
    } else {
        avkind
    };
    if avkind != HeliacalFlags::AVKIND_VR {
        return Err(Error::CError("error: invalid AV kind for the moon".into()));
    }
    if type_event == 1 || type_event == 2 {
        return Err(Error::CError(
            "error: the moon has no morning first or evening last".into(),
        ));
    }

    // C sets TOPOCTR|EQUATORIAL but pheno_ut only needs the phase angle — strip
    // TOPOCTR to avoid requiring a topographic config (C relies on global swe_set_topo).
    let efl = (epheflag & calc::EPHMASK) | CalcFlags::EQUATORIAL;
    let efl = if helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        efl
    } else {
        efl | CalcFlags::NONUT | CalcFlags::TRUEPOS
    };

    let (remapped_event, daystep): (RiseSetFlags, f64) = if type_event == 3 {
        (RiseSetFlags::SET, 1.0)
    } else {
        (RiseSetFlags::RISE, -1.0)
    };

    // New-moon-date determination via pheno_ut phase angle
    let mut jdn_days_ut = tjd_start;
    if type_event == 3 {
        jdn_days_ut += 30.0;
    }
    let ph = eph.pheno_ut(jdn_days_ut, Body::Moon, efl)?;
    let mut phase2 = ph.0.phase_angle;
    let mut goingup = false;
    let mut nm_iters = 0;
    loop {
        let phase1 = phase2;
        jdn_days_ut += daystep;
        let ph = eph.pheno_ut(jdn_days_ut, Body::Moon, efl)?;
        phase2 = ph.0.phase_angle;
        if phase2 > phase1 {
            goingup = true;
        }
        if goingup && phase2 <= phase1 {
            break;
        }
        nm_iters += 1;
        if nm_iters > 10_000 {
            return Err(Error::CError(
                "moon_event_arc_vis: new moon search failure".into(),
            ));
        }
    }
    jdn_days_ut -= daystep; // back to day with smallest phase

    let jdn_days_ut_i = jdn_days_ut;
    jdn_days_ut -= daystep;

    let mut min_tav_oud = 199.0_f64;
    let mut min_tav;
    let mut oldest_min_tav;
    let mut delta_alt = 0.0_f64;
    let mut delta_alt_oud;
    let mut tjd_moonevent;

    loop {
        jdn_days_ut += daystep;

        tjd_moonevent = rise_set(
            eph,
            jdn_days_ut,
            dgeo,
            datm,
            "moon",
            remapped_event,
            epheflag,
            helflag,
            0, // Rim=0 (disc-center)
        )?;

        let tjd_moonevent_start = tjd_moonevent;

        // Inner per-minute loop
        min_tav = 199.0;
        #[allow(unused_assignments)]
        {
            oldest_min_tav = min_tav_oud;
            delta_alt_oud = delta_alt;
        }

        let mut inner_iters = 0;
        loop {
            oldest_min_tav = min_tav_oud;
            min_tav_oud = min_tav;
            delta_alt_oud = delta_alt;

            tjd_moonevent -= (1.0 / 60.0 / 24.0) * sgn(daystep);

            let alt_s = object_loc(eph, tjd_moonevent, dgeo, datm, "sun", 0, epheflag, helflag)?;
            let alt_o = object_loc(eph, tjd_moonevent, dgeo, datm, "moon", 0, epheflag, helflag)?;
            delta_alt = alt_o - alt_s;

            min_tav = deter_tav(
                eph,
                dobs,
                tjd_moonevent,
                dgeo,
                datm,
                "moon",
                epheflag,
                helflag,
            )?;

            let time_check = tjd_moonevent - (8.0 / 60.0 / 24.0) * sgn(daystep);
            let localmin_check =
                deter_tav(eph, dobs, time_check, dgeo, datm, "moon", epheflag, helflag)?;

            inner_iters += 1;
            if !((min_tav <= min_tav_oud || localmin_check < min_tav)
                && (tjd_moonevent - tjd_moonevent_start).abs() < 120.0 / 60.0 / 24.0)
            {
                break;
            }
            if inner_iters > 10_000 {
                break;
            }
        }

        if !(delta_alt_oud < min_tav_oud && (jdn_days_ut - jdn_days_ut_i).abs() < 15.0) {
            break;
        }
    }

    if (jdn_days_ut - jdn_days_ut_i).abs() < 15.0 {
        let extrax = x2min(min_tav, min_tav_oud, oldest_min_tav);
        tjd_moonevent += (1.0 - extrax) * sgn(daystep) / 60.0 / 24.0;
    } else {
        return Err(Error::CError("no date found for lunar event".into()));
    }

    Ok([tjd_moonevent, 0.0, 0.0])
}

fn heliacal_ut_arc_vis(
    eph: &Ephemeris,
    tjd_start: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event_in: i32,
) -> Result<[f64; 3], Error> {
    let planet = object_to_body(object_name);

    let mut objectmagn = magnitude(eph, tjd_start, dgeo, object_name, epheflag, helflag)?;

    let efl = (epheflag & calc::EPHMASK) | CalcFlags::TOPOCTR | CalcFlags::EQUATORIAL;
    let efl = if helflag.contains(HeliacalFlags::HIGH_PRECISION) {
        efl
    } else {
        efl | CalcFlags::NONUT | CalcFlags::TRUEPOS
    };
    let tc = topo_config(eph, dgeo, epheflag);

    let (mut day_step, maxlength): (f64, f64) = match planet {
        Some(Body::Mercury) => (1.0, 100.0),
        Some(Body::Venus) => (64.0, 384.0),
        Some(Body::Mars) => (128.0, 640.0),
        Some(Body::Jupiter) => (64.0, 384.0),
        Some(Body::Saturn) => (64.0, 256.0),
        _ => (64.0, 256.0),
    };

    let mut eventtype = type_event_in;
    if eventtype == 2 {
        day_step = -day_step;
    }
    if eventtype == 4 {
        eventtype = 1;
        day_step = -day_step;
    }
    if eventtype == 3 {
        eventtype = 2;
    }
    let eventtype_with_disc =
        RiseSetFlags::from_bits_truncate(eventtype as u32) | RiseSetFlags::DISC_CENTER;

    // Outer adaptive day-stepping search
    let mut jdn_days_ut_final = tjd_start + maxlength;
    let mut jdn_days_ut_step = tjd_start - 1.0;
    if day_step < 0.0 {
        std::mem::swap(&mut jdn_days_ut_step, &mut jdn_days_ut_final);
        jdn_days_ut_step = tjd_start + 1.0;
        jdn_days_ut_final = tjd_start - maxlength;
    }
    jdn_days_ut_step -= day_step;

    let mut arcus_vis_delta = 199.0_f64;
    let mut arcus_vis_pto = -5.55_f64;
    #[allow(unused_assignments)]
    let mut jdn_arc_vis_ut = tjd_start;
    let mut doneoneday = false;

    loop {
        if day_step.abs() == 1.0 {
            doneoneday = true;
        }

        // Inner loop: step through days
        let mut jdn_days_ut_step_oud;
        let mut arcus_vis_delta_oud;
        loop {
            jdn_days_ut_step_oud = jdn_days_ut_step;
            arcus_vis_delta_oud = arcus_vis_delta;
            jdn_days_ut_step += day_step;

            // Sun rise/set at this candidate day
            let tret = my_rise_trans(
                eph,
                jdn_days_ut_step,
                Body::Sun,
                Some(""),
                eventtype_with_disc,
                epheflag,
                helflag,
                dgeo,
                datm,
            )?;

            // Sun equatorial position at its rise/set
            let dt =
                crate::deltat::calc_deltat(tret, &eph.effective_config(epheflag, eph.config()));
            let tjd_tt = tret + dt;
            let xs = eph.calc_with_config(tjd_tt, Body::Sun, efl, &tc)?;
            let xaz_sun = eph.azalt(
                tret,
                AzAltDir::EquToHor,
                [dgeo[0], dgeo[1], dgeo[2]],
                datm[0],
                datm[1],
                LAPSE_RATE_DEFAULT,
                [xs.data[0], xs.data[1]],
            );

            let t_rise = hour_angle(xaz_sun[1], xs.data[1], dgeo[1]);

            let mut sunsangle = arcus_vis_pto;
            if helflag.contains(HeliacalFlags::AVKIND_MIN7) {
                sunsangle = -7.0;
            }
            if helflag.contains(HeliacalFlags::AVKIND_MIN9) {
                sunsangle = -9.0;
            }

            let t_heliacal = hour_angle(sunsangle, xs.data[1], dgeo[1]);
            let mut t_delta = t_heliacal - t_rise;
            if type_event_in == 2 || type_event_in == 3 {
                t_delta = -t_delta;
            }

            jdn_arc_vis_ut = tret - t_delta / 24.0;

            // Recompute Sun position at candidate instant
            let dt2 = crate::deltat::calc_deltat(
                jdn_arc_vis_ut,
                &eph.effective_config(epheflag, eph.config()),
            );
            let xs2 = eph.calc_with_config(jdn_arc_vis_ut + dt2, Body::Sun, efl, &tc)?;
            let xaz_sun2 = eph.azalt(
                jdn_arc_vis_ut,
                AzAltDir::EquToHor,
                [dgeo[0], dgeo[1], dgeo[2]],
                datm[0],
                datm[1],
                LAPSE_RATE_DEFAULT,
                [xs2.data[0], xs2.data[1]],
            );
            let azi_s = normalize_degrees(xaz_sun2[0] + 180.0);
            let alt_s = xaz_sun2[1];

            // Object/star position at candidate instant
            let eff = eph.effective_config(epheflag, eph.config());
            let (azi_o, alt_o) = if let Some(body) = planet {
                let dt3 = crate::deltat::calc_deltat(jdn_arc_vis_ut, &eff);
                let xp = eph.calc_with_config(jdn_arc_vis_ut + dt3, body, efl, &tc)?;
                objectmagn = magnitude(eph, jdn_arc_vis_ut, dgeo, object_name, epheflag, helflag)?;
                let xaz_p = eph.azalt(
                    jdn_arc_vis_ut,
                    AzAltDir::EquToHor,
                    [dgeo[0], dgeo[1], dgeo[2]],
                    datm[0],
                    datm[1],
                    LAPSE_RATE_DEFAULT,
                    [xp.data[0], xp.data[1]],
                );
                (normalize_degrees(xaz_p[0] + 180.0), xaz_p[1])
            } else {
                // Fixed star — magnitude NOT refreshed per §6 step 5
                let dt3 = crate::deltat::calc_deltat(jdn_arc_vis_ut, &eff);
                let star_result =
                    eph.fixstar2_with_config(object_name, jdn_arc_vis_ut + dt3, efl, &tc)?;
                let xaz_s = eph.azalt(
                    jdn_arc_vis_ut,
                    AzAltDir::EquToHor,
                    [dgeo[0], dgeo[1], dgeo[2]],
                    datm[0],
                    datm[1],
                    LAPSE_RATE_DEFAULT,
                    [star_result.1.data[0], star_result.1.data[1]],
                );
                (normalize_degrees(xaz_s[0] + 180.0), xaz_s[1])
            };

            let delta_alt = alt_o - alt_s;

            // HeliacalAngle with AltM=-1, AziM=0 (Moon interference never factored)
            let dang = heliacal_angle_core(
                objectmagn,
                dobs,
                azi_o,
                -1.0,
                0.0,
                jdn_arc_vis_ut,
                azi_s,
                dgeo,
                datm,
                helflag,
            );
            let arcus_vis = dang.arcus_visionis;
            arcus_vis_pto = dang.sun_altitude_diff;
            arcus_vis_delta = delta_alt - arcus_vis;

            if !((arcus_vis_delta_oud > 0.0 || arcus_vis_delta < 0.0)
                && (jdn_days_ut_final - jdn_days_ut_step) * sgn(day_step) > 0.0)
            {
                break;
            }
        }

        // Backoff-on-first-bracket
        if !doneoneday && (jdn_days_ut_final - jdn_days_ut_step) * sgn(day_step) > 0.0 {
            arcus_vis_delta = arcus_vis_delta_oud;
            jdn_days_ut_step = jdn_days_ut_step_oud;
            day_step = ((day_step.abs() / 2.0) as i32) as f64 * sgn(day_step);
        } else {
            break;
        }
    }

    // Window-exhaustion check
    let d = (jdn_days_ut_final - jdn_days_ut_step) * sgn(day_step);
    if d <= 0.0 || d >= maxlength {
        return Err(Error::CError(format!(
            "heliacal event not found within maxlength {}",
            maxlength
        )));
    }

    let time_step_default = 1.0 / 24.0 / 60.0; // 1 minute in days
    let mut direct = time_step_default;
    if day_step < 0.0 {
        direct = -direct;
    }

    // AVKIND_VR per-minute walkthrough
    if helflag.contains(HeliacalFlags::AVKIND_VR) {
        let time_step = direct;
        let mut tb_vr = 0.0_f64;
        let mut time_pointer = jdn_arc_vis_ut;

        let mut oldest_min_tav = deter_tav(
            eph,
            dobs,
            time_pointer,
            dgeo,
            datm,
            object_name,
            epheflag,
            helflag,
        )?;
        time_pointer += time_step;
        let mut min_tav_oud = deter_tav(
            eph,
            dobs,
            time_pointer,
            dgeo,
            datm,
            object_name,
            epheflag,
            helflag,
        )?;

        let mut min_tav_act;
        if min_tav_oud > oldest_min_tav {
            time_pointer = jdn_arc_vis_ut;
            direct = -direct;
            min_tav_act = oldest_min_tav;
        } else {
            min_tav_act = min_tav_oud;
            min_tav_oud = oldest_min_tav;
        }

        let mut vr_iters = 0;
        loop {
            time_pointer += direct;
            oldest_min_tav = min_tav_oud;
            min_tav_oud = min_tav_act;
            min_tav_act = deter_tav(
                eph,
                dobs,
                time_pointer,
                dgeo,
                datm,
                object_name,
                epheflag,
                helflag,
            )?;

            if min_tav_oud < min_tav_act {
                let extrax = x2min(min_tav_act, min_tav_oud, oldest_min_tav);
                tb_vr = time_pointer - (1.0 - extrax) * direct;
            }

            vr_iters += 1;
            if tb_vr != 0.0 || vr_iters > 10_000 {
                break;
            }
        }

        if tb_vr == 0.0 {
            return Err(Error::CError(
                "heliacal_ut_arc_vis: AVKIND_VR refinement did not converge".into(),
            ));
        }
        jdn_arc_vis_ut = tb_vr;
    }

    // AVKIND_PTO symmetric-crossing averaging
    if helflag.contains(HeliacalFlags::AVKIND_PTO) {
        let mut oude_datum;
        let mut angle;
        let mut pto_iters = 0;
        loop {
            oude_datum = jdn_arc_vis_ut;
            jdn_arc_vis_ut -= direct;

            // Object altitude at this instant
            let dt = crate::deltat::calc_deltat(
                jdn_arc_vis_ut,
                &eph.effective_config(epheflag, eph.config()),
            );
            let tjd_tt = jdn_arc_vis_ut + dt;
            if let Some(body) = planet {
                let xp = eph.calc_with_config(tjd_tt, body, efl, &tc)?;
                let xaz_p = eph.azalt(
                    jdn_arc_vis_ut,
                    AzAltDir::EquToHor,
                    [dgeo[0], dgeo[1], dgeo[2]],
                    datm[0],
                    datm[1],
                    LAPSE_RATE_DEFAULT,
                    [xp.data[0], xp.data[1]],
                );
                angle = xaz_p[1];
            } else {
                let star_result = eph.fixstar2_with_config(object_name, tjd_tt, efl, &tc)?;
                let xaz_s = eph.azalt(
                    jdn_arc_vis_ut,
                    AzAltDir::EquToHor,
                    [dgeo[0], dgeo[1], dgeo[2]],
                    datm[0],
                    datm[1],
                    LAPSE_RATE_DEFAULT,
                    [star_result.1.data[0], star_result.1.data[1]],
                );
                angle = xaz_s[1];
            }

            pto_iters += 1;
            if angle <= 0.0 || pto_iters > 10_000 {
                break;
            }
        }
        jdn_arc_vis_ut = (jdn_arc_vis_ut + oude_datum) / 2.0;
    }

    // Sanity bound — C uses stale JDNDaysUT variable; we use tjd_start (intentional fix)
    if !(-9999999.0..=9999999.0).contains(&jdn_arc_vis_ut) {
        return Err(Error::CError("no heliacal date found".into()));
    }

    Ok([jdn_arc_vis_ut, 0.0, 0.0])
}

// ── §7 Top-level dispatch ────────────────────────────────────────

fn heliacal_ut_dispatch(
    eph: &Ephemeris,
    tjd_start: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
    type_event: i32,
) -> Result<[f64; 3], Error> {
    if helflag.intersects(HeliacalFlags::AVKIND) {
        heliacal_ut_arc_vis(
            eph,
            tjd_start,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag,
            type_event,
        )
    } else {
        heliacal_ut_vis_lim(
            eph,
            tjd_start,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag,
            type_event,
        )
    }
}

/// Port of `swe_heliacal_ut` (swehel.c:3385-3511). Finds the next heliacal
/// event (rising, setting, evening first, morning last, acronychal) for the
/// named celestial object after `tjd_start_ut`.
pub fn heliacal_ut(
    eph: &Ephemeris,
    tjd_start_ut: f64,
    dgeo: &[f64; 3],
    datm: &mut [f64; 4],
    dobs: &mut [f64; 6],
    object_name: &str,
    event: HeliacalEventType,
    epheflag: CalcFlags,
    helflag: HeliacalFlags,
) -> Result<HeliacalEvent, Error> {
    // Altitude validation
    if !(crate::constants::RISE_SET_GEOALT_MIN..=crate::constants::RISE_SET_GEOALT_MAX)
        .contains(&dgeo[2])
    {
        return Err(Error::CError(format!(
            "location for heliacal events must be between {} and {} m above sea",
            crate::constants::RISE_SET_GEOALT_MIN,
            crate::constants::RISE_SET_GEOALT_MAX,
        )));
    }

    let max_count = if helflag.contains(HeliacalFlags::LONG_SEARCH) {
        MAX_COUNT_SYNPER_MAX
    } else {
        MAX_COUNT_SYNPER
    };

    let object_lower = tolower_string_star(object_name);
    default_heliacal_parameters(datm, dgeo, dobs, helflag);

    // C retries on retval==-2 ("not found this period") but propagates ERR (-1)
    // immediately. Our CError conflates both. Classify: "not found" / "does not
    // happen" / "no date found" messages are retryable; everything else
    // (convergence failure, unsupported body, etc.) is fatal.
    let is_retryable = |e: &Error| -> bool {
        match e {
            Error::CircumpolarBody => true,
            Error::CError(msg) => {
                msg.contains("not found")
                    || msg.contains("does not happen")
                    || msg.contains("no date found")
            }
            _ => false,
        }
    };

    let planet = object_to_body(&object_lower);
    let mut type_event = event as i32;

    // Sun rejection
    if matches!(planet, Some(Body::Sun)) {
        return Err(Error::CError(
            "the sun has no heliacal rising or setting".into(),
        ));
    }

    let tjd0 = tjd_start_ut;

    // Moon branch
    if matches!(planet, Some(Body::Moon)) {
        if type_event == 1 || type_event == 2 {
            return Err(Error::CError(
                "the moon has no morning first or evening last".into(),
            ));
        }

        let mut tjd = tjd0;
        let mut dret = moon_event_jd_ut(eph, tjd, dgeo, datm, dobs, epheflag, helflag, type_event)?;
        let mut retval_ok = true;

        // Retry-forward loop until result >= tjd0
        while retval_ok && dret[0] < tjd0 {
            tjd += 15.0;
            match moon_event_jd_ut(eph, tjd, dgeo, datm, dobs, epheflag, helflag, type_event) {
                Ok(d) => dret = d,
                Err(ref e) if is_retryable(e) => {
                    retval_ok = false;
                }
                Err(e) => return Err(e),
            }
        }

        if !retval_ok {
            return Err(Error::CError("no heliacal date found for the moon".into()));
        }

        return Ok(HeliacalEvent {
            start_visible: dret[0],
            optimum_visibility: dret[1],
            end_visible: dret[2],
        });
    }

    // Planets/stars branch
    let planet_id = planet.map(|b| b.to_raw_id()).unwrap_or(-1);

    // Event-type applicability gate
    if !helflag.intersects(HeliacalFlags::AVKIND) {
        // vis_lim path
        if (planet.is_none() || planet_id >= 4) && (type_event == 3 || type_event == 4) {
            return Err(Error::CError(
                "evening first/morning last not provided for outer planets/stars via vis_lim path"
                    .into(),
            ));
        }
        if type_event == 5 || type_event == 6 {
            return Err(Error::CError(
                "acronychal rising/setting is not provided for the vis_lim path".into(),
            ));
        }
    } else {
        // arc_vis path: remap acronychal TypeEvent 5/6 to 3/4
        if planet.is_none() || planet_id >= 4 {
            if type_event == 5 {
                type_event = 3;
            } else if type_event == 6 {
                type_event = 4;
            }
        }
    }

    let body_for_period = planet.unwrap_or(Body::Sun);
    let dsynperiod = get_synodic_period(body_for_period);
    let tjdmax = tjd0 + dsynperiod * max_count as f64;
    let tadd = if matches!(planet, Some(Body::Mercury)) {
        30.0
    } else {
        dsynperiod * 0.6
    };

    // Outer synodic-period loop
    let mut tjd = tjd0;
    let mut last_dret: Option<[f64; 3]> = None;
    let mut found = false;

    while tjd < tjdmax && !found {
        match heliacal_ut_dispatch(
            eph,
            tjd,
            dgeo,
            datm,
            dobs,
            &object_lower,
            epheflag,
            helflag,
            type_event,
        ) {
            Ok(mut dret) => {
                // Inner retry-forward loop until result >= tjd0
                while dret[0] < tjd0 {
                    tjd += tadd;
                    if tjd >= tjdmax {
                        break;
                    }
                    match heliacal_ut_dispatch(
                        eph,
                        tjd,
                        dgeo,
                        datm,
                        dobs,
                        &object_lower,
                        epheflag,
                        helflag,
                        type_event,
                    ) {
                        Ok(d) => dret = d,
                        Err(ref e) if is_retryable(e) => break,
                        Err(e) => return Err(e),
                    }
                }
                if dret[0] >= tjd0 {
                    last_dret = Some(dret);
                    found = true;
                }
            }
            Err(ref e) if is_retryable(e) => {}
            Err(e) => return Err(e),
        }
        tjd += tadd;
    }

    // Final result classification
    if let Some(dret) = last_dret {
        if helflag.contains(HeliacalFlags::SEARCH_1_PERIOD) && dret[0] > tjd0 + dsynperiod * 1.5 {
            return Err(Error::CError(
                "no heliacal date found within this synodic period".into(),
            ));
        }
        Ok(HeliacalEvent {
            start_visible: dret[0],
            optimum_visibility: dret[1],
            end_visible: dret[2],
        })
    } else {
        Err(Error::CError(format!(
            "no heliacal date found within {} synodic periods",
            max_count
        )))
    }
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
