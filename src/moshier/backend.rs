use crate::constants::{
    DEGTORAD, EARTH_MOON_MRAT, J1900, MOON_SPEED_INTV, MOSHLUEPH_END, MOSHLUEPH_START,
    MOSHPLEPH_END, MOSHPLEPH_START, PLAN_SPEED_INTV,
};

pub struct PipelinePositions {
    pub planet_helio: [f64; 6],
    pub earth_helio: [f64; 6],
}
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{normalize_degrees, polar_to_cartesian, rotate_x_sincos};
use crate::obliquity::obliquity;
use crate::precession::precess;
use crate::types::{AstroModels, Body, EphemerisSource, Epsilon, PrecessionDirection};

use super::moon::moshmoon2;
use super::planets::moshplan2;
use super::tables;

fn embofs_mosh(jd: f64, earth: &mut [f64; 3], eps_date: &Epsilon) {
    let t = (jd - J1900) / 36525.0;

    let mp = normalize_degrees(((1.44e-5 * t + 0.009192) * t + 477198.8491) * t + 296.104608);
    let d_half = normalize_degrees(((1.9e-6 * t - 0.001436) * t + 445267.1142) * t + 350.737486);
    let f = normalize_degrees(((-3e-7 * t - 0.003211) * t + 483202.0251) * t + 11.250889);

    let mp_rad = mp * DEGTORAD;
    let (smp, cmp) = mp_rad.sin_cos();
    let d2_rad = 2.0 * DEGTORAD * d_half;
    let (s2d, c2d) = d2_rad.sin_cos();
    let f_rad = f * DEGTORAD;
    let (sf, cf) = f_rad.sin_cos();

    let sx = s2d * cmp - c2d * smp; // sin(2D - MP)
    let cx = c2d * cmp + s2d * smp; // cos(2D - MP)
    let s2mp = 2.0 * smp * cmp;
    let c2mp = cmp * cmp - smp * smp;
    let s2f = 2.0 * sf * cf;

    let mut l = ((1.9e-6 * t - 0.001133) * t + 481267.8831) * t + 270.434164;
    let m = normalize_degrees(((-3.3e-6 * t - 1.50e-4) * t + 35999.0498) * t + 358.475833);
    let sm = (m * DEGTORAD).sin();

    l = l + 6.288750 * smp + 1.274018 * sx + 0.658309 * s2d + 0.213616 * s2mp
        - 0.185596 * sm
        - 0.114336 * s2f;

    let b = 5.128189 * sf
        + 0.280606 * (smp * cf + cmp * sf)  // sin(MP + F)
        + 0.277693 * (smp * cf - cmp * sf)  // sin(MP - F)
        + 0.173238 * (s2d * cf - c2d * sf); // sin(2D - F)

    let p = 0.950724 + 0.051818 * cmp + 0.009531 * cx + 0.007843 * c2d + 0.002824 * c2mp;
    let dist = 4.263523e-5 / (p * DEGTORAD).sin();

    let l = normalize_degrees(l);
    let l_rad = l * DEGTORAD;
    let b_rad = b * DEGTORAD;
    let mut moon_xyz = polar_to_cartesian([l_rad, b_rad, dist]);
    moon_xyz = rotate_x_sincos(moon_xyz, -eps_date.sin_eps, eps_date.cos_eps);
    precess(
        &mut moon_xyz,
        jd,
        CalcFlags::empty(),
        &AstroModels::default(),
        PrecessionDirection::DateToJ2000,
    );
    let factor = 1.0 / (EARTH_MOON_MRAT + 1.0);
    earth[0] -= moon_xyz[0] * factor;
    earth[1] -= moon_xyz[1] * factor;
    earth[2] -= moon_xyz[2] * factor;
}

fn helio_to_equatorial_j2000(jd: f64, table: &super::PlantTbl, eps_j2000: &Epsilon) -> [f64; 3] {
    let pol = moshplan2(jd, table);
    let cart = polar_to_cartesian(pol);
    rotate_x_sincos(cart, -eps_j2000.sin_eps, eps_j2000.cos_eps)
}

fn earth_position(jd: f64, eps_j2000: &Epsilon, eps_date: &Epsilon) -> [f64; 3] {
    let mut pos = helio_to_equatorial_j2000(jd, &tables::EAR404, eps_j2000);
    embofs_mosh(jd, &mut pos, eps_date);
    pos
}

fn compute_planet(
    jd: f64,
    table: Option<&super::PlantTbl>,
    eps_j2000: &Epsilon,
    eps_date: &Epsilon,
) -> [f64; 6] {
    let earth = earth_position(jd, eps_j2000, eps_date);
    let earth_prev = earth_position(jd - PLAN_SPEED_INTV, eps_j2000, eps_date);

    let (geo, geo_prev) = match table {
        Some(tbl) => {
            let planet = helio_to_equatorial_j2000(jd, tbl, eps_j2000);
            let planet_prev = helio_to_equatorial_j2000(jd - PLAN_SPEED_INTV, tbl, eps_j2000);
            (
                [
                    planet[0] - earth[0],
                    planet[1] - earth[1],
                    planet[2] - earth[2],
                ],
                [
                    planet_prev[0] - earth_prev[0],
                    planet_prev[1] - earth_prev[1],
                    planet_prev[2] - earth_prev[2],
                ],
            )
        }
        None => (
            [-earth[0], -earth[1], -earth[2]],
            [-earth_prev[0], -earth_prev[1], -earth_prev[2]],
        ),
    };

    [
        geo[0],
        geo[1],
        geo[2],
        (geo[0] - geo_prev[0]) / PLAN_SPEED_INTV,
        (geo[1] - geo_prev[1]) / PLAN_SPEED_INTV,
        (geo[2] - geo_prev[2]) / PLAN_SPEED_INTV,
    ]
}

fn moon_equatorial_j2000(eval_jd: f64, eps_date: &Epsilon) -> [f64; 3] {
    let pol = moshmoon2(eval_jd);
    let cart = polar_to_cartesian(pol);
    let mut pos = rotate_x_sincos(cart, -eps_date.sin_eps, eps_date.cos_eps);
    precess(
        &mut pos,
        eval_jd,
        CalcFlags::empty(),
        &AstroModels::default(),
        PrecessionDirection::DateToJ2000,
    );
    pos
}

fn compute_moon(jd: f64, eps_date: &Epsilon) -> [f64; 6] {
    let pos = moon_equatorial_j2000(jd, eps_date);
    let pos_plus = moon_equatorial_j2000(jd + MOON_SPEED_INTV, eps_date);
    let pos_minus = moon_equatorial_j2000(jd - MOON_SPEED_INTV, eps_date);

    let mut result = [0.0; 6];
    for i in 0..3 {
        result[i] = pos[i];
        let b = (pos_plus[i] - pos_minus[i]) / 2.0;
        let a = (pos_plus[i] + pos_minus[i]) / 2.0 - pos[i];
        result[i + 3] = (2.0 * a + b) / MOON_SPEED_INTV;
    }
    result
}

fn planet_table(body: Body) -> Result<&'static super::PlantTbl, Error> {
    match body {
        Body::Mercury => Ok(&tables::MER404),
        Body::Venus => Ok(&tables::VEN404),
        Body::Mars => Ok(&tables::MAR404),
        Body::Jupiter => Ok(&tables::JUP404),
        Body::Saturn => Ok(&tables::SAT404),
        Body::Uranus => Ok(&tables::URA404),
        Body::Neptune => Ok(&tables::NEP404),
        Body::Pluto => Ok(&tables::PLU404),
        _ => Err(Error::EphemerisNotAvailable {
            body,
            source: EphemerisSource::Moshier,
        }),
    }
}

pub fn planet_helio_velocity_at(
    jd: f64,
    body: Body,
    eps_j2000: &Epsilon,
) -> Result<[f64; 3], Error> {
    let table = planet_table(body)?;
    let pos = helio_to_equatorial_j2000(jd, table, eps_j2000);
    let pos_prev = helio_to_equatorial_j2000(jd - PLAN_SPEED_INTV, table, eps_j2000);
    Ok([
        (pos[0] - pos_prev[0]) / PLAN_SPEED_INTV,
        (pos[1] - pos_prev[1]) / PLAN_SPEED_INTV,
        (pos[2] - pos_prev[2]) / PLAN_SPEED_INTV,
    ])
}

pub fn earth_helio_velocity_at(jd: f64, eps_j2000: &Epsilon) -> [f64; 3] {
    let eps_date = obliquity(jd, CalcFlags::empty(), &AstroModels::default());
    let eps_date_prev = obliquity(
        jd - PLAN_SPEED_INTV,
        CalcFlags::empty(),
        &AstroModels::default(),
    );
    let earth = earth_position(jd, eps_j2000, &eps_date);
    let earth_prev = earth_position(jd - PLAN_SPEED_INTV, eps_j2000, &eps_date_prev);
    [
        (earth[0] - earth_prev[0]) / PLAN_SPEED_INTV,
        (earth[1] - earth_prev[1]) / PLAN_SPEED_INTV,
        (earth[2] - earth_prev[2]) / PLAN_SPEED_INTV,
    ]
}

pub fn compute_pipeline(
    jd: f64,
    body: Body,
    eps_j2000: &Epsilon,
) -> Result<PipelinePositions, Error> {
    let eps_date = obliquity(jd, CalcFlags::empty(), &AstroModels::default());
    let eps_date_prev = obliquity(
        jd - PLAN_SPEED_INTV,
        CalcFlags::empty(),
        &AstroModels::default(),
    );

    let earth = earth_position(jd, eps_j2000, &eps_date);
    let earth_prev = earth_position(jd - PLAN_SPEED_INTV, eps_j2000, &eps_date_prev);
    let earth_helio = [
        earth[0],
        earth[1],
        earth[2],
        (earth[0] - earth_prev[0]) / PLAN_SPEED_INTV,
        (earth[1] - earth_prev[1]) / PLAN_SPEED_INTV,
        (earth[2] - earth_prev[2]) / PLAN_SPEED_INTV,
    ];

    let planet_helio = match body {
        Body::Sun => [0.0; 6],
        Body::Moon => {
            let moon = moon_equatorial_j2000(jd, &eps_date);
            let moon_plus = moon_equatorial_j2000(jd + MOON_SPEED_INTV, &eps_date);
            let moon_minus = moon_equatorial_j2000(jd - MOON_SPEED_INTV, &eps_date);
            let mut result = [0.0; 6];
            for i in 0..3 {
                result[i] = moon[i];
                let b = (moon_plus[i] - moon_minus[i]) / 2.0;
                let a = (moon_plus[i] + moon_minus[i]) / 2.0 - moon[i];
                result[i + 3] = (2.0 * a + b) / MOON_SPEED_INTV;
            }
            result
        }
        _ => {
            let table = planet_table(body)?;
            let pos = helio_to_equatorial_j2000(jd, table, eps_j2000);
            let pos_prev = helio_to_equatorial_j2000(jd - PLAN_SPEED_INTV, table, eps_j2000);
            [
                pos[0],
                pos[1],
                pos[2],
                (pos[0] - pos_prev[0]) / PLAN_SPEED_INTV,
                (pos[1] - pos_prev[1]) / PLAN_SPEED_INTV,
                (pos[2] - pos_prev[2]) / PLAN_SPEED_INTV,
            ]
        }
    };

    Ok(PipelinePositions {
        planet_helio,
        earth_helio,
    })
}

pub fn compute(jd: f64, body: Body, eps_j2000: &Epsilon) -> Result<[f64; 6], Error> {
    match body {
        Body::Moon => {
            if jd < MOSHLUEPH_START - 0.2 || jd > MOSHLUEPH_END + 0.2 {
                return Err(Error::BeyondEphemerisLimits {
                    jd_tt: jd,
                    start: MOSHLUEPH_START,
                    end: MOSHLUEPH_END,
                });
            }
            let eps_date = obliquity(jd, CalcFlags::empty(), &AstroModels::default());
            Ok(compute_moon(jd, &eps_date))
        }
        Body::Earth => Ok([0.0; 6]),
        body => {
            let table = match body {
                Body::Sun => None,
                Body::Mercury => Some(&tables::MER404),
                Body::Venus => Some(&tables::VEN404),
                Body::Mars => Some(&tables::MAR404),
                Body::Jupiter => Some(&tables::JUP404),
                Body::Saturn => Some(&tables::SAT404),
                Body::Uranus => Some(&tables::URA404),
                Body::Neptune => Some(&tables::NEP404),
                Body::Pluto => Some(&tables::PLU404),
                _ => {
                    return Err(Error::EphemerisNotAvailable {
                        body,
                        source: EphemerisSource::Moshier,
                    });
                }
            };
            if jd < MOSHPLEPH_START - 0.3 || jd > MOSHPLEPH_END + 0.3 {
                return Err(Error::BeyondEphemerisLimits {
                    jd_tt: jd,
                    start: MOSHPLEPH_START,
                    end: MOSHPLEPH_END,
                });
            }
            let eps_date = obliquity(jd, CalcFlags::empty(), &AstroModels::default());
            Ok(compute_planet(jd, table, eps_j2000, &eps_date))
        }
    }
}
