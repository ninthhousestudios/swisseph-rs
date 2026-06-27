use crate::bias::frame_bias;
use crate::constants::*;
use crate::corrections::{aberr_light, deflect_light};
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{
    cartesian_to_polar_with_speed, diff_radians, polar_to_cartesian_with_speed, rotate_x_sincos,
};
use crate::moshier::backend::{
    PipelinePositions, compute_pipeline, earth_helio_velocity_at, planet_helio_velocity_at,
};
use crate::nutation::nutation;
use crate::obliquity::obliquity;
use crate::precession::{ldp_peps, precess};
use crate::types::{
    AstroModels, Body, Epsilon, FrameTransform, Nutation as NutationType, PrecessionDirection,
    PrecessionModel,
};

const EPHMASK: CalcFlags = CalcFlags::MOSEPH
    .union(CalcFlags::SWIEPH)
    .union(CalcFlags::JPLEPH);

pub fn plaus_iflag(mut flags: CalcFlags) -> CalcFlags {
    if flags.contains(CalcFlags::DPSIDEPS_1980) && flags.contains(CalcFlags::JPLHOR_APPROX) {
        flags.remove(CalcFlags::JPLHOR_APPROX);
    }
    if flags.contains(CalcFlags::TOPOCTR) {
        flags.remove(CalcFlags::HELCTR | CalcFlags::BARYCTR);
    }
    if flags.contains(CalcFlags::BARYCTR) {
        flags.remove(CalcFlags::HELCTR);
    }
    if flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
        flags.insert(CalcFlags::NOABERR | CalcFlags::NOGDEFL);
    }
    if flags.contains(CalcFlags::J2000) {
        flags.insert(CalcFlags::NONUT);
    }
    if flags.contains(CalcFlags::SIDEREAL) {
        flags.insert(CalcFlags::NONUT);
        flags.remove(CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX);
    }
    if flags.contains(CalcFlags::TRUEPOS) {
        flags.insert(CalcFlags::NOGDEFL | CalcFlags::NOABERR);
    }
    if flags.contains(CalcFlags::XYZ) {
        flags.remove(CalcFlags::RADIANS);
    }
    if flags.contains(CalcFlags::SPEED) && flags.contains(CalcFlags::SPEED3) {
        flags.remove(CalcFlags::SPEED3);
    }

    // Ephemeris selection: force MOSEPH for now (only backend available).
    // Clear Horizons flags — they only apply to JPL ephemeris.
    flags.remove(EPHMASK | CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX);
    flags.insert(CalcFlags::MOSEPH);

    flags
}

fn precess_and_ephem(
    xx: &mut [f64; 6],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> (Epsilon, NutationType, Option<NutationType>) {
    if !flags.contains(CalcFlags::J2000) {
        let mut pos3 = [xx[0], xx[1], xx[2]];
        precess(
            &mut pos3,
            jd,
            flags,
            models,
            PrecessionDirection::J2000ToDate,
        );
        xx[0] = pos3[0];
        xx[1] = pos3[1];
        xx[2] = pos3[2];
        if flags.contains(CalcFlags::SPEED) {
            precess_speed(xx, jd, flags, models, PrecessionDirection::J2000ToDate);
        }
        let eps = obliquity(jd, flags, models);
        let nut_val = nutation(jd, flags, models);
        let nutv = if flags.contains(CalcFlags::SPEED) {
            Some(nutation(jd - NUT_SPEED_INTV, flags, models))
        } else {
            None
        };
        (eps, nut_val, nutv)
    } else {
        (
            obliquity(J2000, flags, models),
            nutation(jd, flags, models),
            None,
        )
    }
}

fn precess_speed(
    xx: &mut [f64; 6],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    direction: PrecessionDirection,
) {
    let fac = match direction {
        PrecessionDirection::J2000ToDate => 1.0,
        PrecessionDirection::DateToJ2000 => -1.0,
    };

    // Rotate velocity through the precession matrix
    let mut vel3 = [xx[3], xx[4], xx[5]];
    precess(&mut vel3, jd, flags, models, direction);
    xx[3] = vel3[0];
    xx[4] = vel3[1];
    xx[5] = vel3[2];

    // C uses oec (of date) for J2000→Date, oec2000 for Date→J2000
    let oe = match direction {
        PrecessionDirection::J2000ToDate => obliquity(jd, flags, models),
        PrecessionDirection::DateToJ2000 => obliquity(J2000, flags, models),
    };

    // Equatorial → ecliptic (position and velocity)
    let pos_ecl = rotate_x_sincos([xx[0], xx[1], xx[2]], oe.sin_eps, oe.cos_eps);
    let vel_ecl = rotate_x_sincos([xx[3], xx[4], xx[5]], oe.sin_eps, oe.cos_eps);
    let mut ecl = [
        pos_ecl[0], pos_ecl[1], pos_ecl[2], vel_ecl[0], vel_ecl[1], vel_ecl[2],
    ];

    // Cartesian → polar
    ecl = cartesian_to_polar_with_speed(ecl);

    // Add precession rate to ecliptic longitude speed
    if models.prec_longterm == PrecessionModel::Vondrak2011 {
        let dpre = ldp_peps(jd);
        let dpre2 = ldp_peps(jd + 1.0);
        ecl[3] += (dpre2 - dpre) * fac;
    } else {
        let tprec = (jd - J2000) / 36525.0;
        ecl[3] += (50.290966 + 0.0222226 * tprec) / 3600.0 * DEGTORAD / 365.25 * fac;
    }

    // Polar → cartesian
    ecl = polar_to_cartesian_with_speed(ecl);

    // Ecliptic → equatorial (negate sin for inverse rotation)
    let pos_eq = rotate_x_sincos([ecl[0], ecl[1], ecl[2]], -oe.sin_eps, oe.cos_eps);
    let vel_eq = rotate_x_sincos([ecl[3], ecl[4], ecl[5]], -oe.sin_eps, oe.cos_eps);
    xx[0] = pos_eq[0];
    xx[1] = pos_eq[1];
    xx[2] = pos_eq[2];
    xx[3] = vel_eq[0];
    xx[4] = vel_eq[1];
    xx[5] = vel_eq[2];
}

fn nut_matrix(eps: &Epsilon, nut: &NutationType) -> [[f64; 3]; 3] {
    let (sin_eps, cos_eps) = (eps.eps.sin(), eps.eps.cos());
    let sin_nut_eps = (eps.eps + nut.deps).sin();
    let cos_nut_eps = (eps.eps + nut.deps).cos();
    let (sin_dpsi, cos_dpsi) = (nut.dpsi.sin(), nut.dpsi.cos());
    [
        [cos_dpsi, -sin_dpsi * cos_eps, -sin_dpsi * sin_eps],
        [
            sin_dpsi * cos_nut_eps,
            cos_dpsi * cos_nut_eps * cos_eps + sin_nut_eps * sin_eps,
            cos_dpsi * cos_nut_eps * sin_eps - sin_nut_eps * cos_eps,
        ],
        [
            sin_dpsi * sin_nut_eps,
            cos_dpsi * sin_nut_eps * cos_eps - cos_nut_eps * sin_eps,
            cos_dpsi * sin_nut_eps * sin_eps + cos_nut_eps * cos_eps,
        ],
    ]
}

fn nutate(
    pos: &mut [f64; 6],
    eps: &Epsilon,
    nut: &NutationType,
    nutv: Option<&NutationType>,
    has_speed: bool,
) {
    let matrix = nut_matrix(eps, nut);

    let x = pos[0];
    let y = pos[1];
    let z = pos[2];
    pos[0] = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
    pos[1] = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
    pos[2] = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;

    if has_speed {
        let vx = pos[3];
        let vy = pos[4];
        let vz = pos[5];
        pos[3] = matrix[0][0] * vx + matrix[0][1] * vy + matrix[0][2] * vz;
        pos[4] = matrix[1][0] * vx + matrix[1][1] * vy + matrix[1][2] * vz;
        pos[5] = matrix[2][0] * vx + matrix[2][1] * vy + matrix[2][2] * vz;

        // Apparent motion from nutation rate change (same obliquity, earlier nutation)
        if let Some(nv) = nutv {
            let matv = nut_matrix(eps, nv);
            let xv0 = matv[0][0] * x + matv[0][1] * y + matv[0][2] * z;
            let xv1 = matv[1][0] * x + matv[1][1] * y + matv[1][2] * z;
            let xv2 = matv[2][0] * x + matv[2][1] * y + matv[2][2] * z;
            pos[3] += (pos[0] - xv0) / NUT_SPEED_INTV;
            pos[4] += (pos[1] - xv1) / NUT_SPEED_INTV;
            pos[5] += (pos[2] - xv2) / NUT_SPEED_INTV;
        }
    }
}

fn app_pos_rest(
    xx: &mut [f64; 6],
    flags: CalcFlags,
    eps: &Epsilon,
    nut: &NutationType,
    nutv: Option<&NutationType>,
) -> [f64; 24] {
    let mut xreturn = [0.0; 24];
    let has_speed = flags.contains(CalcFlags::SPEED);

    // Step 1: Nutation (equatorial cartesian)
    if !flags.contains(CalcFlags::NONUT) {
        nutate(xx, eps, nut, nutv, has_speed);
    }
    xreturn[18..24].copy_from_slice(xx);

    // Step 2: Equatorial → ecliptic via obliquity rotation
    let (sin_eps, cos_eps) = (eps.sin_eps, eps.cos_eps);
    let pos3 = rotate_x_sincos([xx[0], xx[1], xx[2]], sin_eps, cos_eps);
    xx[0] = pos3[0];
    xx[1] = pos3[1];
    xx[2] = pos3[2];
    if has_speed {
        let vel3 = rotate_x_sincos([xx[3], xx[4], xx[5]], sin_eps, cos_eps);
        xx[3] = vel3[0];
        xx[4] = vel3[1];
        xx[5] = vel3[2];
    }

    // Step 3: Nutation obliquity rotation (ecliptic nutation)
    if !flags.contains(CalcFlags::NONUT) {
        let sin_nut = nut.deps.sin();
        let cos_nut = nut.deps.cos();
        let pos3 = rotate_x_sincos([xx[0], xx[1], xx[2]], sin_nut, cos_nut);
        xx[0] = pos3[0];
        xx[1] = pos3[1];
        xx[2] = pos3[2];
        if has_speed {
            let vel3 = rotate_x_sincos([xx[3], xx[4], xx[5]], sin_nut, cos_nut);
            xx[3] = vel3[0];
            xx[4] = vel3[1];
            xx[5] = vel3[2];
        }
    }
    xreturn[6..12].copy_from_slice(xx);

    // Step 4: Polar coordinates
    let eq_polar = cartesian_to_polar_with_speed([
        xreturn[18],
        xreturn[19],
        xreturn[20],
        xreturn[21],
        xreturn[22],
        xreturn[23],
    ]);
    xreturn[12..18].copy_from_slice(&eq_polar);

    let ecl_polar = cartesian_to_polar_with_speed([
        xreturn[6],
        xreturn[7],
        xreturn[8],
        xreturn[9],
        xreturn[10],
        xreturn[11],
    ]);
    xreturn[0..6].copy_from_slice(&ecl_polar);

    // Step 5: Radians → degrees (angles only, not distances)
    for i in 0..2 {
        xreturn[i] *= RADTODEG;
        xreturn[i + 3] *= RADTODEG;
        xreturn[i + 12] *= RADTODEG;
        xreturn[i + 15] *= RADTODEG;
    }

    xreturn
}

pub fn calc_planet(
    jd: f64,
    body: Body,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<[f64; 24], Error> {
    let pp = compute_pipeline(jd, body, eps_j2000)?;
    let PipelinePositions {
        planet_helio,
        earth_helio,
    } = pp;

    // Geocentric position
    let mut xx = [0.0; 6];
    for i in 0..6 {
        xx[i] = planet_helio[i] - earth_helio[i];
    }

    // Light-time (C gates entire block on !TRUEPOS)
    let mut dt = 0.0;
    if !flags.contains(CalcFlags::TRUEPOS) {
        let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
        dt = dist * AUNIT / CLIGHT / 86400.0;
        // Moshier niter=0: linear approximation only
        for i in 0..3 {
            xx[i] = planet_helio[i] - dt * planet_helio[i + 3] - earth_helio[i];
        }
        if flags.contains(CalcFlags::SPEED) {
            // Velocity at apparent time: C calls swi_moshplan at retarded time t
            // and takes only the velocity, subtracts Earth velocity at teval
            let vel_at_t = planet_helio_velocity_at(jd - dt, body, eps_j2000)?;
            for i in 0..3 {
                xx[i + 3] = vel_at_t[i] - earth_helio[i + 3];
            }
            // xxsp: change-of-dt speed correction. Light-time changes as the
            // planet moves, affecting apparent speed. Correction =
            // (dt - dt_prev) * planet_helio_vel, where dt_prev is light-time
            // at t-1 day.
            let geo_prev = [
                planet_helio[0] - earth_helio[0] - (planet_helio[3] - earth_helio[3]),
                planet_helio[1] - earth_helio[1] - (planet_helio[4] - earth_helio[4]),
                planet_helio[2] - earth_helio[2] - (planet_helio[5] - earth_helio[5]),
            ];
            let dist_prev =
                (geo_prev[0].powi(2) + geo_prev[1].powi(2) + geo_prev[2].powi(2)).sqrt();
            let dt_sp = dist_prev * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xx[i + 3] -= (dt - dt_sp) * planet_helio[i + 3];
            }
        }
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // Planet heliocentric at retarded time (for deflection geometry)
    let mut planet_helio_retarded = [0.0; 6];
    for i in 0..3 {
        planet_helio_retarded[i] = xx[i] + earth_helio[i];
        planet_helio_retarded[i + 3] = planet_helio[i + 3];
    }

    // Gravitational deflection
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOGDEFL) {
        deflect_light(
            &mut xx,
            &earth_helio,
            &planet_helio_retarded,
            flags.contains(CalcFlags::SPEED),
        );
    }

    // Aberration
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(
            &mut xx,
            &[earth_helio[3], earth_helio[4], earth_helio[5]],
            flags.contains(CalcFlags::SPEED),
        );
        // Earth velocity correction: observer velocity changed between
        // emission (retarded time t) and reception (teval)
        if flags.contains(CalcFlags::SPEED) {
            let earth_vel_t = earth_helio_velocity_at(jd - dt, eps_j2000);
            for i in 0..3 {
                xx[i + 3] += earth_helio[i + 3] - earth_vel_t[i];
            }
        }
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // Frame bias (ICRS → J2000)
    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }

    // Precession + ephemeris data
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);

    Ok(app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()))
}

pub fn calc_sun(
    jd: f64,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<[f64; 24], Error> {
    let pp = compute_pipeline(jd, Body::Sun, eps_j2000)?;
    let earth_helio = pp.earth_helio;

    // Geocentric Sun = -Earth heliocentric
    // For Moshier, Sun is at heliocentric origin — light-time retardation of
    // the Sun gives zero (it doesn't move). Earth position stays at time t.
    let mut xx = [0.0; 6];
    for i in 0..3 {
        xx[i] = -earth_helio[i];
        xx[i + 3] = -earth_helio[i + 3];
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // No deflection for Sun

    // Aberration
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(
            &mut xx,
            &[earth_helio[3], earth_helio[4], earth_helio[5]],
            flags.contains(CalcFlags::SPEED),
        );
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // Frame bias
    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }

    // Precession + ephemeris data
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);

    Ok(app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()))
}

pub fn calc_moon(
    jd: f64,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<[f64; 24], Error> {
    let pp = compute_pipeline(jd, Body::Moon, eps_j2000)?;
    let earth_helio = pp.earth_helio;

    // Moon is already geocentric from backend (planet_helio is geocentric for Moon)
    let mut xx = pp.planet_helio;

    // Light-time (C gates entire light-time on !TRUEPOS for Moon)
    if !flags.contains(CalcFlags::TRUEPOS) {
        let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
        let dt = dist * AUNIT / CLIGHT / 86400.0;
        // C does a barycentric detour: converts geocentric→barycentric, retards
        // with barycentric velocity, then subtracts unretarded Earth. Net effect
        // on geocentric position: subtract dt * (geo_vel + earth_vel).
        for i in 0..3 {
            xx[i] -= dt * (xx[i + 3] + earth_helio[i + 3]);
        }
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // No deflection for Moon (too close for GR bending)

    // Aberration
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(
            &mut xx,
            &[earth_helio[3], earth_helio[4], earth_helio[5]],
            flags.contains(CalcFlags::SPEED),
        );
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // Frame bias
    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }

    // Precession + ephemeris data
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);

    Ok(app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()))
}

pub fn extract_output(xreturn: &[f64; 24], flags: CalcFlags) -> [f64; 6] {
    let base = if flags.contains(CalcFlags::EQUATORIAL) {
        12
    } else {
        0
    };
    let cart = if flags.contains(CalcFlags::XYZ) { 6 } else { 0 };
    let src = base + cart;

    let mut data = [
        xreturn[src],
        xreturn[src + 1],
        xreturn[src + 2],
        0.0,
        0.0,
        0.0,
    ];
    if flags.contains(CalcFlags::SPEED) {
        data[3] = xreturn[src + 3];
        data[4] = xreturn[src + 4];
        data[5] = xreturn[src + 5];
    }
    if flags.contains(CalcFlags::RADIANS) && !flags.contains(CalcFlags::XYZ) {
        data[0] *= DEGTORAD;
        data[1] *= DEGTORAD;
        data[3] *= DEGTORAD;
        data[4] *= DEGTORAD;
    }
    data
}

pub fn extract_ecl_nut(ecl_nut: &[f64; 6], flags: CalcFlags) -> [f64; 6] {
    if flags.intersects(CalcFlags::EQUATORIAL | CalcFlags::XYZ) {
        return [0.0; 6];
    }
    let mut data = *ecl_nut;
    if flags.contains(CalcFlags::RADIANS) {
        for v in &mut data[..4] {
            *v *= DEGTORAD;
        }
    }
    data
}

pub fn speed3_interval(body: Body) -> f64 {
    match body {
        Body::Moon => MOON_SPEED_INTV,
        Body::OscuApogee | Body::TrueNode => 0.1,
        _ => PLAN_SPEED_INTV,
    }
}

pub fn denormalize_positions(x0: &mut [f64; 24], x1: &[f64; 24], x2: &mut [f64; 24]) {
    // Only ecliptic longitude [0] and right ascension [12] can wrap ±360°.
    for i in [0, 12] {
        if x1[i] - x0[i] < -180.0 {
            x0[i] -= 360.0;
        }
        if x1[i] - x0[i] > 180.0 {
            x0[i] += 360.0;
        }
        if x1[i] - x2[i] < -180.0 {
            x2[i] -= 360.0;
        }
        if x1[i] - x2[i] > 180.0 {
            x2[i] += 360.0;
        }
    }
}

pub fn calc_speed_3point(x1: &mut [f64; 24], x0: &[f64; 24], x2: &[f64; 24], dt: f64) {
    // Quadratic interpolation derivative at t+dt (matches C's calc_speed).
    for base in [0, 6, 12, 18] {
        for j in 0..3 {
            let k = base + j;
            let b = (x2[k] - x0[k]) / 2.0;
            let a = (x2[k] + x0[k]) / 2.0 - x1[k];
            x1[k + 3] = (2.0 * a + b) / dt;
        }
    }
}

fn mean_element_pipeline(
    xx: &mut [f64; 6],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> [f64; 24] {
    // Ecliptic polar → ecliptic cartesian (with speed)
    let cart = polar_to_cartesian_with_speed(*xx);
    *xx = cart;

    // Ecliptic → equatorial: rotate by -obliquity of date
    let eps_date = obliquity(jd, flags, models);
    let pos_eq = rotate_x_sincos([xx[0], xx[1], xx[2]], -eps_date.sin_eps, eps_date.cos_eps);
    let vel_eq = rotate_x_sincos([xx[3], xx[4], xx[5]], -eps_date.sin_eps, eps_date.cos_eps);
    *xx = [
        pos_eq[0], pos_eq[1], pos_eq[2], vel_eq[0], vel_eq[1], vel_eq[2],
    ];

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // J2000: precess equatorial back to J2000
    let eps = if flags.contains(CalcFlags::J2000) {
        let mut pos3 = [xx[0], xx[1], xx[2]];
        precess(
            &mut pos3,
            jd,
            flags,
            models,
            PrecessionDirection::DateToJ2000,
        );
        xx[0] = pos3[0];
        xx[1] = pos3[1];
        xx[2] = pos3[2];
        if flags.contains(CalcFlags::SPEED) {
            precess_speed(xx, jd, flags, models, PrecessionDirection::DateToJ2000);
        }
        obliquity(J2000, flags, models)
    } else {
        eps_date
    };

    let nut_val = nutation(jd, flags, models);
    let nutv = if flags.contains(CalcFlags::SPEED) {
        Some(nutation(jd - NUT_SPEED_INTV, flags, models))
    } else {
        None
    };

    app_pos_rest(xx, flags, &eps, &nut_val, nutv.as_ref())
}

pub fn calc_mean_node(jd: f64, flags: CalcFlags, models: &AstroModels) -> Result<[f64; 24], Error> {
    let pos = crate::moshier::moon::mean_node(jd)?;
    let pos_prev = crate::moshier::moon::mean_node(jd - MEAN_NODE_SPEED_INTV)?;

    let mut xx = [0.0; 6];
    xx[0] = pos[0];
    xx[1] = pos[1];
    xx[2] = pos[2];
    xx[3] = diff_radians(pos[0], pos_prev[0]) / MEAN_NODE_SPEED_INTV;
    xx[4] = 0.0;
    xx[5] = 0.0;

    let mut xreturn = mean_element_pipeline(&mut xx, jd, flags, models);

    if !flags.contains(CalcFlags::SIDEREAL) && !flags.contains(CalcFlags::J2000) {
        xreturn[1] = 0.0;
        xreturn[4] = 0.0;
        xreturn[5] = 0.0;
        xreturn[8] = 0.0;
        xreturn[11] = 0.0;
    }

    Ok(xreturn)
}

pub fn calc_mean_apogee(
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<[f64; 24], Error> {
    let pos = crate::moshier::moon::mean_apogee(jd)?;
    let pos_prev = crate::moshier::moon::mean_apogee(jd - MEAN_NODE_SPEED_INTV)?;

    let mut xx = [0.0; 6];
    xx[0] = pos[0];
    xx[1] = pos[1];
    xx[2] = pos[2];
    xx[3] = diff_radians(pos[0], pos_prev[0]) / MEAN_NODE_SPEED_INTV;
    xx[4] = diff_radians(pos[1], pos_prev[1]) / MEAN_NODE_SPEED_INTV;
    xx[5] = 0.0;

    let mut xreturn = mean_element_pipeline(&mut xx, jd, flags, models);

    xreturn[5] = 0.0;

    Ok(xreturn)
}

pub fn calc_ecl_nut(jd: f64, flags: CalcFlags, models: &AstroModels) -> [f64; 6] {
    let eps = obliquity(jd, flags, models);
    let nut_val = nutation(jd, flags, models);
    [
        (eps.eps + nut_val.deps) * RADTODEG,
        eps.eps * RADTODEG,
        nut_val.dpsi * RADTODEG,
        nut_val.deps * RADTODEG,
        0.0,
        0.0,
    ]
}
