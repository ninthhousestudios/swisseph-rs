//! Calculation pipeline internals: light-time iteration, aberration, deflection,
//! frame transforms, and SPEED3 numerical differentiation.
//!
//! Low-level internals; exposed for golden tests and advanced use. The primary
//! entry points are [`Ephemeris::calc`](crate::Ephemeris::calc) and
//! [`Ephemeris::calc_ut`](crate::Ephemeris::calc_ut).

use crate::bias::frame_bias;
use crate::config::EphemerisConfig;
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
#[cfg(feature = "swisseph-files")]
use crate::sweph_file::types::SEI_MOON;
#[cfg(feature = "swisseph-files")]
use crate::sweph_file::{
    SEI_FLG_HELIO, SEI_SUNBARY, SwissEphFile, body_file_id, evaluate_body, find_file_for_jd,
};
use crate::types::{
    AstroModels, Body, EphemerisSource, Epsilon, FrameTransform, Nutation as NutationType,
    PrecessionDirection, PrecessionModel,
};

pub(crate) const EPHMASK: CalcFlags = CalcFlags::MOSEPH
    .union(CalcFlags::SWIEPH)
    .union(CalcFlags::JPLEPH);

/// Extract the caller's requested ephemeris source from the EPHMASK bits in
/// `flags`, using C's precedence (sweph.c:375-381): MOSEPH > JPLEPH > SWIEPH.
/// Returns `None` when no EPHMASK bit is set (caller accepts the config default).
pub(crate) fn requested_source(flags: CalcFlags) -> Option<EphemerisSource> {
    if flags.contains(CalcFlags::MOSEPH) {
        Some(EphemerisSource::Moshier)
    } else if flags.contains(CalcFlags::JPLEPH) {
        Some(EphemerisSource::Jpl)
    } else if flags.contains(CalcFlags::SWIEPH) {
        Some(EphemerisSource::Swiss)
    } else {
        None
    }
}

/// Internal: sanitizes and normalizes the caller-supplied `CalcFlags` (plausibilization),
/// resolving mutually exclusive bits and stamping in the resolved ephemeris `source`.
pub fn plaus_iflag(mut flags: CalcFlags, source: EphemerisSource) -> CalcFlags {
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
    // Topocentric + aberration doesn't trust analytic SPEED; C silently forces
    // the 3-point numerical differentiation instead (sweph.c:402-410). Runs
    // after the rule above, so SPEED3 wins here even though SPEED is kept.
    if flags.contains(CalcFlags::SPEED)
        && flags.contains(CalcFlags::TOPOCTR)
        && !flags.contains(CalcFlags::NOABERR)
    {
        flags.insert(CalcFlags::SPEED3);
    }

    flags.remove(EPHMASK | CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX);
    match source {
        EphemerisSource::Swiss => flags.insert(CalcFlags::SWIEPH),
        EphemerisSource::Jpl => flags.insert(CalcFlags::JPLEPH),
        EphemerisSource::Moshier => flags.insert(CalcFlags::MOSEPH),
    }

    flags
}

pub(crate) fn normalize_asteroid_aliases(body: Body) -> Body {
    match body {
        Body::Asteroid(id) if id.mpc_number() == 134340 => Body::Pluto,
        Body::Asteroid(id) if (1..=4).contains(&id.mpc_number()) => {
            Body::try_from(17 + id.mpc_number() - 1).expect("Ceres..Vesta are valid body ids")
        }
        other => other,
    }
}

/// Ports sweph.c:416-437's three-clause CENTER_BODY / PlanetMoon normalization.
/// Returns `(body, moon_raw, flags)` where `moon_raw > 0` means a plmoon
/// offset must be fetched and added, and `flags` may have CENTER_BODY set or
/// cleared relative to the input.
pub(crate) fn normalize_center_body(
    body: Body,
    flags: CalcFlags,
) -> (Body, Option<i32>, CalcFlags) {
    let mut flags = flags;
    let mut moon_raw: Option<i32> = None;

    // Clause (i): planet ipl <= SE_PLUTO + CENTER_BODY → synthesize COB number
    if flags.contains(CalcFlags::CENTER_BODY) {
        let raw = body.to_raw_id();
        if (0..=9).contains(&raw) {
            moon_raw = Some(raw * 100 + 9099);
        }
    }

    // Clause (ii): direct 9pmm ipl → extract parent, force CENTER_BODY
    if let Body::PlanetMoon(id) = body {
        let encoded = id.encoded();
        let raw = PLMOON_OFFSET + encoded;
        moon_raw = Some(raw);
        let parent_raw = (raw - 9000) / 100;
        // parent_raw is 0..=9 by construction (encoded 0..=999)
        let parent = Body::try_from(parent_raw).expect("parent planet is valid");
        flags |= CalcFlags::CENTER_BODY;
        return normalize_center_body_cancel(parent, moon_raw, flags);
    }

    normalize_center_body_cancel(body, moon_raw, flags)
}

fn normalize_center_body_cancel(
    body: Body,
    moon_raw: Option<i32>,
    mut flags: CalcFlags,
) -> (Body, Option<i32>, CalcFlags) {
    // Clause (iii): parent <= SE_MARS && suffix == 99 → cancel
    if flags.contains(CalcFlags::CENTER_BODY) {
        let raw = body.to_raw_id();
        if (0..=4).contains(&raw)
            && let Some(mr) = moon_raw
            && mr % 100 == 99
        {
            flags -= CalcFlags::CENTER_BODY;
            return (body, None, flags);
        }
    }
    (body, moon_raw, flags)
}

/// Observer offset (position + velocity, AU / AU-day, J2000 mean equatorial)
/// at `jd`, or the zero vector when TOPOCTR isn't requested. Stateless
/// equivalent of C's `swed.topd.xobs` cache: recomputed fresh every call
/// (docs/c-ref-topocentric.md §2). `pub(crate)`: also reused by
/// `context::Ephemeris::calc_fixstar_*` (docs/c-ref-fixstar.md step 6).
pub(crate) fn topo_offset(
    jd: f64,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> [f64; 6] {
    if !flags.contains(CalcFlags::TOPOCTR) {
        return [0.0; 6];
    }
    match &config.topographic {
        Some(topo) => crate::topocentric::get_observer(jd, topo, flags, config, models),
        None => [0.0; 6],
    }
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

pub(crate) fn precess_speed(
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

pub(crate) fn nutate(
    pos: &mut [f64; 6],
    eps: &Epsilon,
    nut: &NutationType,
    nutv: Option<&NutationType>,
    has_speed: bool,
    backward: bool,
) {
    // C `swi_nutate`: `backward=false` applies the nutation matrix as-is (mean ->
    // true), `backward=true` applies its transpose (true -> mean, i.e. removes
    // nutation). Our `nut_matrix` is the transpose of C's `swed.nut.matrix`, so the
    // index orders below are the mirror of the C code's `backward` flag.
    let matrix = nut_matrix(eps, nut);

    let x = pos[0];
    let y = pos[1];
    let z = pos[2];
    if backward {
        pos[0] = matrix[0][0] * x + matrix[1][0] * y + matrix[2][0] * z;
        pos[1] = matrix[0][1] * x + matrix[1][1] * y + matrix[2][1] * z;
        pos[2] = matrix[0][2] * x + matrix[1][2] * y + matrix[2][2] * z;
    } else {
        pos[0] = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
        pos[1] = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
        pos[2] = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;
    }

    if has_speed {
        let vx = pos[3];
        let vy = pos[4];
        let vz = pos[5];
        if backward {
            pos[3] = matrix[0][0] * vx + matrix[1][0] * vy + matrix[2][0] * vz;
            pos[4] = matrix[0][1] * vx + matrix[1][1] * vy + matrix[2][1] * vz;
            pos[5] = matrix[0][2] * vx + matrix[1][2] * vy + matrix[2][2] * vz;
        } else {
            pos[3] = matrix[0][0] * vx + matrix[0][1] * vy + matrix[0][2] * vz;
            pos[4] = matrix[1][0] * vx + matrix[1][1] * vy + matrix[1][2] * vz;
            pos[5] = matrix[2][0] * vx + matrix[2][1] * vy + matrix[2][2] * vz;
        }

        // Apparent motion from nutation rate change (same obliquity, earlier nutation)
        if let Some(nv) = nutv {
            let matv = nut_matrix(eps, nv);
            let (xv0, xv1, xv2) = if backward {
                (
                    matv[0][0] * x + matv[1][0] * y + matv[2][0] * z,
                    matv[0][1] * x + matv[1][1] * y + matv[2][1] * z,
                    matv[0][2] * x + matv[1][2] * y + matv[2][2] * z,
                )
            } else {
                (
                    matv[0][0] * x + matv[0][1] * y + matv[0][2] * z,
                    matv[1][0] * x + matv[1][1] * y + matv[1][2] * z,
                    matv[2][0] * x + matv[2][1] * y + matv[2][2] * z,
                )
            };
            pos[3] += (pos[0] - xv0) / NUT_SPEED_INTV;
            pos[4] += (pos[1] - xv1) / NUT_SPEED_INTV;
            pos[5] += (pos[2] - xv2) / NUT_SPEED_INTV;
        }
    }
}

pub(crate) fn app_pos_rest(
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
        nutate(xx, eps, nut, nutv, has_speed, false);
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

/// Common tail for every `SEFLG_HELCTR` path: `xx` is the heliocentric position (light-time
/// corrected, in the ICRS/J2000-equatorial frame — no aberration/deflection, which plaus_iflag
/// forced off, and no geocenter conversion). Applies the shared bias → precession/nutation →
/// `app_pos_rest` finish, exactly like the geocentric paths' tail. Zeroes velocity when SPEED is
/// off. Returns `(xreturn[24], x2000[6])`.
fn finish_helctr(
    mut xx: [f64; 6],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> ([f64; 24], [f64; 6]) {
    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }
    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }
    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);
    (
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    )
}

/// Internal: computes `body`'s apparent position (24-slot `xreturn` plus J2000
/// equatorial `x2000`) using the Moshier backend, running light-time, deflection,
/// aberration, frame bias, and precession/nutation per `flags`.
pub fn calc_planet(
    jd: f64,
    body: Body,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let pp = compute_pipeline(jd, body, eps_j2000)?;
    let PipelinePositions {
        planet_helio,
        earth_helio,
    } = pp;

    // Heliocentric (Moshier): observer is the Sun (xobs = 0 in the Moshier heliocentric frame),
    // so the position is just `planet_helio`, light-time retarded with niter=0 (single analytic
    // pass, no re-evaluation — sweph.c:2659-2664). No aberration/deflection/geocenter conversion.
    if flags.contains(CalcFlags::HELCTR) {
        let mut xx = planet_helio;
        if !flags.contains(CalcFlags::TRUEPOS) {
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xx[i] = planet_helio[i] - dt * planet_helio[i + 3];
            }
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    // Observer: geocenter, or geocenter + topocentric offset (§1 "xobs replaces
    // the geocenter"). Moshier has no separate bary/helio frame, so the offset
    // is added directly to earth_helio.
    let offset = topo_offset(jd, flags, config, models);
    let mut xobs = earth_helio;
    for i in 0..6 {
        xobs[i] += offset[i];
    }

    // Geocentric (or topocentric) position
    let mut xx = [0.0; 6];
    for i in 0..6 {
        xx[i] = planet_helio[i] - xobs[i];
    }

    // Light-time (C gates entire block on !TRUEPOS)
    let mut dt = 0.0;
    if !flags.contains(CalcFlags::TRUEPOS) {
        let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
        dt = dist * AUNIT / CLIGHT / 86400.0;
        // Moshier niter=0: linear approximation only
        for i in 0..3 {
            xx[i] = planet_helio[i] - dt * planet_helio[i + 3] - xobs[i];
        }
        if flags.contains(CalcFlags::SPEED) {
            // Velocity at apparent time: C calls swi_moshplan at retarded time t
            // and takes only the velocity, subtracts observer velocity at teval
            let vel_at_t = planet_helio_velocity_at(jd - dt, body, eps_j2000)?;
            for i in 0..3 {
                xx[i + 3] = vel_at_t[i] - xobs[i + 3];
            }
            // xxsp: change-of-dt speed correction. Light-time changes as the
            // planet moves, affecting apparent speed. Correction =
            // (dt - dt_prev) * planet_helio_vel, where dt_prev is light-time
            // at t-1 day.
            let geo_prev = [
                planet_helio[0] - xobs[0] - (planet_helio[3] - xobs[3]),
                planet_helio[1] - xobs[1] - (planet_helio[4] - xobs[4]),
                planet_helio[2] - xobs[2] - (planet_helio[5] - xobs[5]),
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
        planet_helio_retarded[i] = xx[i] + xobs[i];
        planet_helio_retarded[i + 3] = planet_helio[i + 3];
    }

    // Gravitational deflection
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOGDEFL) {
        deflect_light(
            &mut xx,
            &xobs,
            &planet_helio_retarded,
            flags.contains(CalcFlags::SPEED),
        );
    }

    // Aberration
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(
            &mut xx,
            &[xobs[3], xobs[4], xobs[5]],
            flags.contains(CalcFlags::SPEED),
        );
        // Observer velocity correction: observer velocity changed between
        // emission (retarded time t) and reception (teval). The retarded-epoch
        // observer offset is independently recomputed (§4) — not derived from
        // `offset` at the current epoch.
        if flags.contains(CalcFlags::SPEED) {
            let earth_vel_t = earth_helio_velocity_at(jd - dt, eps_j2000);
            let offset_ret = topo_offset(jd - dt, flags, config, models);
            for i in 0..3 {
                xx[i + 3] += xobs[i + 3] - (earth_vel_t[i] + offset_ret[i + 3]);
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
    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);

    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

/// Internal: computes the Sun's (or, with `is_earth`, Earth's) apparent position
/// using the Moshier backend, following the same light-time/aberration/precession
/// pipeline as [`calc_planet`] with the Sun-specific heliocentric shortcuts.
pub fn calc_sun(
    jd: f64,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
    is_earth: bool,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let pp = compute_pipeline(jd, Body::Sun, eps_j2000)?;
    let earth_helio = pp.earth_helio;

    let offset = topo_offset(jd, flags, config, models);
    let mut xobs = earth_helio;
    for i in 0..6 {
        xobs[i] += offset[i];
    }

    let is_hb = flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR);

    let mut xx = if is_earth && is_hb {
        // HELCTR/BARYCTR Earth: xx = xobs (Moshier's earth_helio IS the
        // heliocentric Earth; Sun ≈ barycenter in Moshier, so helio ≈ bary).
        xobs
    } else {
        // Geocentric Sun = -observer heliocentric.
        let mut v = [0.0; 6];
        for i in 0..3 {
            v[i] = -xobs[i];
            v[i + 3] = -xobs[i + 3];
        }
        v
    };

    // Light-time for HELCTR/BARYCTR Earth (Moshier): C enters the retardation
    // block when HELCTR|BARYCTR even for MOSEPH (sweph.c:3969 gate), re-evaluating
    // Earth at the retarded time via swi_moshplan (niter=1, loop 0..=1).
    // Geocentric Moshier Sun has no light-time (Sun at origin, block skipped).
    if is_earth && is_hb && !flags.contains(CalcFlags::TRUEPOS) {
        for _ in 0..=1 {
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            let pp_ret = compute_pipeline(jd - dt, Body::Sun, eps_j2000)?;
            for i in 0..6 {
                xx[i] = pp_ret.earth_helio[i] + offset[i];
            }
        }
    }

    if !flags.contains(CalcFlags::SPEED) {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // No deflection for Sun (or Earth through this path)

    // Aberration (skipped for HELCTR/BARYCTR via plaus_iflag's NOABERR)
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(
            &mut xx,
            &[xobs[3], xobs[4], xobs[5]],
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
    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);

    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

/// Internal: computes the Moon's apparent position using the Moshier backend,
/// following the same light-time/aberration/precession pipeline as
/// [`calc_planet`] with the Moon-specific geocentric shortcuts.
pub fn calc_moon(
    jd: f64,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let pp = compute_pipeline(jd, Body::Moon, eps_j2000)?;
    let earth_helio = pp.earth_helio;

    // Heliocentric Moon (Moshier): heliocentric Moon = moon_geo (pp.planet_helio) + earth_helio.
    // The Moon path computes light-time dt ONCE from the heliocentric distance (no iteration loop
    // — sweph.c:4147-4152, `xxm`), then retards analytically with the heliocentric velocity
    // (Moshier's barycentric frame == heliocentric; observer Sun == 0, so no final subtraction).
    if flags.contains(CalcFlags::HELCTR) {
        let mut moon_helio = [0.0; 6];
        for i in 0..6 {
            moon_helio[i] = pp.planet_helio[i] + earth_helio[i];
        }
        let mut xx = moon_helio;
        if !flags.contains(CalcFlags::TRUEPOS) {
            let dist = (moon_helio[0] * moon_helio[0]
                + moon_helio[1] * moon_helio[1]
                + moon_helio[2] * moon_helio[2])
                .sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xx[i] = moon_helio[i] - dt * moon_helio[i + 3];
            }
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    let offset = topo_offset(jd, flags, config, models);
    let mut xobs = earth_helio;
    for i in 0..6 {
        xobs[i] += offset[i];
    }

    // Moon is already geocentric from backend (planet_helio is geocentric for
    // Moon); shift to topocentric by subtracting the observer offset directly.
    let mut xx = pp.planet_helio;
    for i in 0..6 {
        xx[i] -= offset[i];
    }

    // Light-time (C gates entire light-time on !TRUEPOS for Moon)
    if !flags.contains(CalcFlags::TRUEPOS) {
        let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
        let dt = dist * AUNIT / CLIGHT / 86400.0;
        // C does a barycentric detour: converts geocentric→barycentric, retards
        // with barycentric velocity, then subtracts unretarded observer. Net
        // effect on geocentric position: subtract dt * (geo_vel + observer_vel).
        for i in 0..3 {
            xx[i] -= dt * (xx[i + 3] + xobs[i + 3]);
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
            &[xobs[3], xobs[4], xobs[5]],
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
    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);

    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

/// Internal: extracts the caller-facing 6-element position/speed vector from
/// the internal 24-slot `xreturn` buffer, selecting ecliptic/equatorial and
/// polar/cartesian variants and applying radians conversion per `flags`.
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

/// Internal: extracts the ecliptic-and-nutation auxiliary output (true/mean
/// obliquity, nutation in longitude/obliquity) for `swe_calc`'s side-channel
/// output, zeroed when `flags` requests equatorial or cartesian output.
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

/// Internal: returns the finite-difference step (days) used for `SPEED3`
/// numerical differentiation of `body`'s position.
pub fn speed3_interval(body: Body) -> f64 {
    match body {
        Body::Moon => MOON_SPEED_INTV,
        Body::OscuApogee | Body::TrueNode => 0.1,
        _ => PLAN_SPEED_INTV,
    }
}

/// Internal: un-wraps the ±360° longitude/right-ascension discontinuity between
/// the three `SPEED3` sample epochs so their finite difference isn't corrupted
/// by a wraparound.
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

/// Internal: fills in `x1`'s speed components (indices 3..6, 9..12, 15..18,
/// 21..24) via central-difference quadratic interpolation from the bracketing
/// samples `x0`/`x2` taken `dt` days apart (C's `calc_speed`).
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
) -> OscOutput {
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

    // J2000 equatorial vector for the SEFLG_SIDEREAL ECL_T0 / SSY_PLANE rigorous
    // branches (applied by the caller's `apply_sidereal`). Mirrors C's
    // app_pos_etc_mean sweph.c:4335-4347: precess the equatorial-of-date vector to
    // J2000 — no nutation removal here, because the mean element carries no nutation
    // yet (app_pos_rest adds it below). Captured before the J2000 re-projection.
    //
    // CRITICAL: the precession is guarded `teval != J2000` (unlike lunar_osc_elem,
    // which precesses unconditionally). At exactly J2000 the position rotation is a
    // no-op but `precess_speed` would still add its precession-rate term (~the
    // ayanamsa rate, 3.8e-5 deg/day) to the velocity — C skips the whole block, so
    // we must too, or the sidereal SPEED comes out short by that rate.
    let x2000 = if flags.contains(CalcFlags::SIDEREAL) {
        let mut x = *xx;
        if jd != J2000 {
            let mut pos3 = [x[0], x[1], x[2]];
            precess(
                &mut pos3,
                jd,
                flags,
                models,
                PrecessionDirection::DateToJ2000,
            );
            x[0..3].copy_from_slice(&pos3);
            if flags.contains(CalcFlags::SPEED) {
                precess_speed(&mut x, jd, flags, models, PrecessionDirection::DateToJ2000);
            }
        }
        x
    } else {
        [0.0; 6]
    };

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

    (
        app_pos_rest(xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    )
}

/// Mean lunar node & perigee longitudes + speeds (degrees, degrees/day) for
/// `swe_nod_aps`'s Moon/MEAN branch. Thin re-export of
/// [`crate::moshier::moon::mean_lunar_elements`] so the `nodaps` app module
/// reaches it through `calc` rather than depending on a backend directly
/// (`app-uses-calc-not-backends`).
pub fn mean_lunar_elements(tjd: f64) -> (f64, f64, f64, f64) {
    crate::moshier::moon::mean_lunar_elements(tjd)
}

/// Internal: computes the mean lunar node's apparent position for `SE_MEAN_NODE`,
/// running the shared mean-element pipeline (ecliptic-to-equatorial rotation,
/// optional sidereal/J2000 precession, nutation) on top of the raw mean-node ephemeris.
pub fn calc_mean_node(jd: f64, flags: CalcFlags, models: &AstroModels) -> Result<OscOutput, Error> {
    let pos = crate::moshier::moon::mean_node(jd)?;
    let pos_prev = crate::moshier::moon::mean_node(jd - MEAN_NODE_SPEED_INTV)?;

    let mut xx = [0.0; 6];
    xx[0] = pos[0];
    xx[1] = pos[1];
    xx[2] = pos[2];
    xx[3] = diff_radians(pos[0], pos_prev[0]) / MEAN_NODE_SPEED_INTV;
    xx[4] = 0.0;
    xx[5] = 0.0;

    let (mut xreturn, x2000) = mean_element_pipeline(&mut xx, jd, flags, models);

    if !flags.contains(CalcFlags::SIDEREAL) && !flags.contains(CalcFlags::J2000) {
        xreturn[1] = 0.0;
        xreturn[4] = 0.0;
        xreturn[5] = 0.0;
        xreturn[8] = 0.0;
        xreturn[11] = 0.0;
    }

    Ok((xreturn, x2000))
}

/// Internal: computes the mean lunar apogee's apparent position for
/// `SE_MEAN_APOG`, running the shared mean-element pipeline on top of the raw
/// mean-apogee ephemeris.
pub fn calc_mean_apogee(
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
) -> Result<OscOutput, Error> {
    let pos = crate::moshier::moon::mean_apogee(jd)?;
    let pos_prev = crate::moshier::moon::mean_apogee(jd - MEAN_NODE_SPEED_INTV)?;

    let mut xx = [0.0; 6];
    xx[0] = pos[0];
    xx[1] = pos[1];
    xx[2] = pos[2];
    xx[3] = diff_radians(pos[0], pos_prev[0]) / MEAN_NODE_SPEED_INTV;
    xx[4] = diff_radians(pos[1], pos_prev[1]) / MEAN_NODE_SPEED_INTV;
    xx[5] = 0.0;

    let (mut xreturn, x2000) = mean_element_pipeline(&mut xx, jd, flags, models);

    xreturn[5] = 0.0;

    Ok((xreturn, x2000))
}

/// Internal: computes the ecliptic obliquity and nutation side-channel output
/// (true obliquity, mean obliquity, nutation in longitude, nutation in
/// obliquity; degrees) returned alongside `swe_calc`'s main position.
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

// ---------------------------------------------------------------------------
// Osculating lunar node / apogee (SE_TRUE_NODE / SE_OSCU_APOG)
// Port of `lunar_osc_elem` + `swi_plan_for_osc_elem` (sweph.c:5168, 5758).
// See docs/c-ref-nodaps.md Parts C, D.
// ---------------------------------------------------------------------------

/// Port of C `swi_plan_for_osc_elem` (sweph.c:5758), as it actually compiles in
/// the default build. Rotates a raw (pre-bias, GCRS-equatorial) geocentric moon
/// position+speed 6-vector into the ecliptic-of-date frame that the osculating
/// ellipse computation needs.
///
/// CRITICAL: the `SEFLG_SIDEREAL`/`SEFLG_J2000` short-circuits that skip
/// precession live inside `#ifdef SID_TNODE_FROM_ECL_T0`, which is NOT defined
/// anywhere in the C tree. So this ALWAYS precesses J2000->date and ALWAYS uses
/// obliquity-of-date, regardless of the `J2000` flag (the ref doc's Part C
/// pseudocode is wrong on this point — verified against the C source). The
/// speed vector is precessed and nutated as a PURE rotation — no precession-rate
/// or nutation-rate term (`swi_precess`, not `swi_precess_speed`; matrix rotation
/// with `nutv=None`) — matching the C comments "daily precession 0.137\" may not
/// be added" and "again: speed vector must be rotated, but not added 'speed' of
/// nutation".
pub(crate) fn plan_for_osc_elem(
    flags: CalcFlags,
    tjd: f64,
    xx: &mut [f64; 6],
    models: &AstroModels,
) {
    // ICRS -> J2000 bias (unconditional on !ICRS, matching the Moon pipeline;
    // C gates on swi_get_denum>=403, which holds for every backend here).
    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(xx, tjd, flags, models, FrameTransform::GcrsToJ2000);
    }

    // Precession J2000 -> equator of date. Position AND speed via the SAME
    // rotation (`precess`, NOT `precess_speed`).
    let mut pos3 = [xx[0], xx[1], xx[2]];
    let mut vel3 = [xx[3], xx[4], xx[5]];
    precess(
        &mut pos3,
        tjd,
        flags,
        models,
        PrecessionDirection::J2000ToDate,
    );
    precess(
        &mut vel3,
        tjd,
        flags,
        models,
        PrecessionDirection::J2000ToDate,
    );
    xx[0..3].copy_from_slice(&pos3);
    xx[3..6].copy_from_slice(&vel3);

    let oe = obliquity(tjd, flags, models);
    let nut_opt = if flags.contains(CalcFlags::NONUT) {
        None
    } else {
        Some(nutation(tjd, flags, models))
    };

    // Equatorial nutation matrix (pure rotation of pos + speed, no rate term).
    if let Some(nut_val) = nut_opt.as_ref() {
        nutate(xx, &oe, nut_val, None, true, false);
    }

    // Equatorial -> ecliptic (rotate by +obliquity), pos + speed.
    let p = rotate_x_sincos([xx[0], xx[1], xx[2]], oe.sin_eps, oe.cos_eps);
    let v = rotate_x_sincos([xx[3], xx[4], xx[5]], oe.sin_eps, oe.cos_eps);
    xx[0..3].copy_from_slice(&p);
    xx[3..6].copy_from_slice(&v);

    // Ecliptic nutation (rotate by nutation-in-obliquity), pos + speed.
    if let Some(nut_val) = nut_opt.as_ref() {
        let (sn, cn) = (nut_val.deps.sin(), nut_val.deps.cos());
        let p = rotate_x_sincos([xx[0], xx[1], xx[2]], sn, cn);
        let v = rotate_x_sincos([xx[3], xx[4], xx[5]], sn, cn);
        xx[0..3].copy_from_slice(&p);
        xx[3..6].copy_from_slice(&v);
    }
}

/// A `lunar_osc_elem` output half: the `xreturn[24]` output layout plus the
/// `x2000` J2000 equatorial vector the SIDEREAL rigorous branches consume.
type OscOutput = ([f64; 24], [f64; 6]);

/// Assemble the `xreturn[24]` output layout from an ecliptic-of-date cartesian
/// 6-vector `x` (node or apogee). Port of `lunar_osc_elem`'s D.4 output-frame
/// stage (sweph.c:5487-5591). `oe` is obliquity-of-date. No physical effects
/// here — light-time/precession/nutation already happened per-sample in D.1.
///
/// Returns `(xreturn, x2000)`. `x2000` is the J2000 equatorial cartesian (pos +
/// speed) that the `SEFLG_SIDEREAL` ECL_T0 / SSY_PLANE rigorous branches consume
/// via the caller's `apply_sidereal` — C's `lunar_osc_elem` builds it by removing
/// the full nutation matrix from the equatorial-of-date vector and precessing to
/// J2000 (sweph.c:5527-5540). It is only populated when `SEFLG_SIDEREAL` is set
/// (else `[0.0; 6]`; the traditional-ayanamsa and non-sidereal paths never read
/// it, and `apply_sidereal` treats an all-zero `x2000` as "not available").
fn osc_output_frame(
    x: &[f64; 6],
    flags: CalcFlags,
    oe: &Epsilon,
    tjd: f64,
    models: &AstroModels,
) -> OscOutput {
    let has_speed = flags.contains(CalcFlags::SPEED);
    let mut xr = [0.0; 24];

    // Ecliptic cartesian + polar.
    xr[6..12].copy_from_slice(x);
    let ecl_pol = cartesian_to_polar_with_speed([xr[6], xr[7], xr[8], xr[9], xr[10], xr[11]]);
    xr[0..6].copy_from_slice(&ecl_pol);

    // Ecliptic -> equatorial cartesian (rotate by -obliquity).
    let p = rotate_x_sincos([xr[6], xr[7], xr[8]], -oe.sin_eps, oe.cos_eps);
    xr[18..21].copy_from_slice(&p);
    if has_speed {
        let v = rotate_x_sincos([xr[9], xr[10], xr[11]], -oe.sin_eps, oe.cos_eps);
        xr[21..24].copy_from_slice(&v);
    }

    // Remove ecliptic-nutation rotation from the equatorial vector (unless NONUT).
    if !flags.contains(CalcFlags::NONUT) {
        let nut_val = nutation(tjd, flags, models);
        let (sn, cn) = (nut_val.deps.sin(), nut_val.deps.cos());
        let p = rotate_x_sincos([xr[18], xr[19], xr[20]], -sn, cn);
        xr[18..21].copy_from_slice(&p);
        if has_speed {
            let v = rotate_x_sincos([xr[21], xr[22], xr[23]], -sn, cn);
            xr[21..24].copy_from_slice(&v);
        }
    }

    // Equatorial polar.
    let eq_pol = cartesian_to_polar_with_speed([xr[18], xr[19], xr[20], xr[21], xr[22], xr[23]]);
    xr[12..18].copy_from_slice(&eq_pol);

    // J2000 equatorial vector for the SEFLG_SIDEREAL ECL_T0 / SSY_PLANE rigorous
    // branches (applied by the caller's `apply_sidereal`). Mirrors C's
    // lunar_osc_elem sweph.c:5527-5540: take the equatorial-of-date vector (its
    // nutation-in-obliquity already removed above), remove the full nutation
    // matrix, then precess to J2000. Captured here — before the J2000 re-projection
    // below overwrites xr[18..24] — and only when SIDEREAL is requested.
    let x2000 = if flags.contains(CalcFlags::SIDEREAL) {
        let mut x = [xr[18], xr[19], xr[20], xr[21], xr[22], xr[23]];
        if !flags.contains(CalcFlags::NONUT) {
            let nut_val = nutation(tjd, flags, models);
            let nutv = if has_speed {
                Some(nutation(tjd - NUT_SPEED_INTV, flags, models))
            } else {
                None
            };
            nutate(&mut x, oe, &nut_val, nutv.as_ref(), has_speed, true);
        }
        let mut pos3 = [x[0], x[1], x[2]];
        precess(
            &mut pos3,
            tjd,
            flags,
            models,
            PrecessionDirection::DateToJ2000,
        );
        x[0..3].copy_from_slice(&pos3);
        if has_speed {
            precess_speed(&mut x, tjd, flags, models, PrecessionDirection::DateToJ2000);
        }
        x
    } else {
        [0.0; 6]
    };

    // SEFLG_J2000 re-projection (sweph.c:5561-5577): the node/apogee are referred
    // to the equator/ecliptic of date; transform the equatorial vector to J2000.
    // Position via `precess`, speed via `precess_speed` (WITH the rate term —
    // unlike plan_for_osc_elem's pure rotation; this matches C's swi_precess +
    // swi_precess_speed here). Skipped under SIDEREAL — C's D.4 makes the sidereal
    // and plain-J2000 handling mutually exclusive (`if SIDEREAL ... else if J2000`).
    if flags.contains(CalcFlags::J2000) && !flags.contains(CalcFlags::SIDEREAL) {
        let mut x6 = [xr[18], xr[19], xr[20], xr[21], xr[22], xr[23]];
        let mut pos3 = [x6[0], x6[1], x6[2]];
        precess(
            &mut pos3,
            tjd,
            flags,
            models,
            PrecessionDirection::DateToJ2000,
        );
        x6[0..3].copy_from_slice(&pos3);
        if has_speed {
            precess_speed(
                &mut x6,
                tjd,
                flags,
                models,
                PrecessionDirection::DateToJ2000,
            );
        }
        xr[18..24].copy_from_slice(&x6);
        let eq_pol = cartesian_to_polar_with_speed(x6);
        xr[12..18].copy_from_slice(&eq_pol);

        // Equatorial-J2000 -> ecliptic-J2000 (rotate by +obliquity of J2000).
        let oe2000 = obliquity(J2000, flags, models);
        let p = rotate_x_sincos([x6[0], x6[1], x6[2]], oe2000.sin_eps, oe2000.cos_eps);
        xr[6..9].copy_from_slice(&p);
        if has_speed {
            let v = rotate_x_sincos([x6[3], x6[4], x6[5]], oe2000.sin_eps, oe2000.cos_eps);
            xr[9..12].copy_from_slice(&v);
        }
        let ecl_pol = cartesian_to_polar_with_speed([xr[6], xr[7], xr[8], xr[9], xr[10], xr[11]]);
        xr[0..6].copy_from_slice(&ecl_pol);
    }

    // Radians -> degrees (angles only) + degnorm on the two longitudes.
    for i in 0..2 {
        xr[i] *= RADTODEG;
        xr[i + 3] *= RADTODEG;
        xr[i + 12] *= RADTODEG;
        xr[i + 15] *= RADTODEG;
    }
    xr[0] = crate::math::normalize_degrees(xr[0]);
    xr[12] = crate::math::normalize_degrees(xr[12]);
    (xr, x2000)
}

/// Core of `lunar_osc_elem` (sweph.c:5360-5592, D.2-D.4). Given the three raw
/// geocentric equatorial-J2000 moon samples (pos+vel, light-time corrected;
/// only `[istart..=2]` valid), computes BOTH the osculating node and the
/// osculating apogee together (the node is always needed to derive the apogee
/// and vice versa), and returns `((node_xreturn, node_x2000), (apogee_xreturn,
/// apogee_x2000))`. The `x2000` companions carry the J2000 equatorial vector the
/// `SEFLG_SIDEREAL` ECL_T0 / SSY_PLANE branches need (see `osc_output_frame`).
///
/// Sample epochs: `raw_samples[0]` at `tjd - speed_intv`, `[1]` at
/// `tjd + speed_intv`, `[2]` at `tjd`. `speed_intv` is the backend-specific
/// central-difference interval (`NODE_CALC_INTV` / `NODE_CALC_INTV_MOSH`).
pub(crate) fn lunar_osc_elem(
    tjd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    raw_samples: &[[f64; 6]; 3],
    istart: usize,
    speed_intv: f64,
) -> (OscOutput, OscOutput) {
    let has_speed = flags.contains(CalcFlags::SPEED);
    let plan_flags = flags | CalcFlags::SPEED; // C always passes iflag|SEFLG_SPEED

    // D.1 tail: rotate each raw sample into ecliptic-of-date via plan_for_osc_elem.
    let mut xpos = [[0.0f64; 6]; 3];
    for i in istart..=2 {
        let t = match i {
            0 => tjd - speed_intv,
            1 => tjd + speed_intv,
            _ => tjd,
        };
        let mut s = raw_samples[i];
        plan_for_osc_elem(plan_flags, t, &mut s, models);
        xpos[i] = s;
    }

    // D.2: node direction (tangent-line intersection with the ecliptic).
    // xx[i][0..3] is the node position 3-vector for sample i.
    let mut xx = [[0.0f64; 6]; 3];
    for i in istart..=2 {
        if xpos[i][5].abs() < 1e-15 {
            xpos[i][5] = 1e-15; // clamp persists into D.3 (cross_prod reads xpos[i][5])
        }
        let fac = xpos[i][2] / xpos[i][5];
        let sgn = xpos[i][5] / xpos[i][5].abs();
        for j in 0..3 {
            xx[i][j] = (xpos[i][j] - fac * xpos[i][j + 3]) * sgn;
        }
    }

    // D.3: apogee (osculating ellipse) + ellipse-corrected node distance.
    let gmsm =
        GEOGCONST * (1.0 + 1.0 / EARTH_MOON_MRAT) / AUNIT / AUNIT / AUNIT * 86400.0 * 86400.0;
    let mut xxa = [[0.0f64; 6]; 3];
    for i in istart..=2 {
        // node direction, unit angle
        let mut rxy = (xx[i][0] * xx[i][0] + xx[i][1] * xx[i][1]).sqrt();
        let cosnode = xx[i][0] / rxy;
        let sinnode = xx[i][1] / rxy;
        // inclination from the orbital angular-momentum vector
        let xnorm = crate::math::cross_prod(
            [xpos[i][0], xpos[i][1], xpos[i][2]],
            [xpos[i][3], xpos[i][4], xpos[i][5]],
        );
        rxy = xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
        let c2 = rxy + xnorm[2] * xnorm[2];
        let mut rxyz = c2.sqrt();
        rxy = rxy.sqrt();
        let sinincl = rxy / rxyz;
        let cosincl = (1.0 - sinincl * sinincl).sqrt();
        // argument of latitude
        let cosu = xpos[i][0] * cosnode + xpos[i][1] * sinnode;
        let sinu = xpos[i][2] / sinincl;
        let uu = sinu.atan2(cosu);
        // semi-major axis
        rxyz = (xpos[i][0] * xpos[i][0] + xpos[i][1] * xpos[i][1] + xpos[i][2] * xpos[i][2]).sqrt();
        let v2 = xpos[i][3] * xpos[i][3] + xpos[i][4] * xpos[i][4] + xpos[i][5] * xpos[i][5];
        let sema = 1.0 / (2.0 / rxyz - v2 / gmsm);
        // eccentricity
        let pp = c2 / gmsm;
        let ecce = (1.0 - pp / sema).sqrt();
        // eccentric anomaly
        let mut cos_e = 1.0 / ecce * (1.0 - rxyz / sema);
        let dot = xpos[i][0] * xpos[i][3] + xpos[i][1] * xpos[i][4] + xpos[i][2] * xpos[i][5];
        let sin_e = 1.0 / ecce / (sema * gmsm).sqrt() * dot;
        // true anomaly
        let mut ny = 2.0 * (((1.0 + ecce) / (1.0 - ecce)).sqrt() * sin_e / (1.0 + cos_e)).atan();
        // apogee = perihelion + PI, distance a(1+e) unconditionally
        let mut apg = [
            crate::math::normalize_radians(uu - ny + std::f64::consts::PI),
            0.0,
            sema * (1.0 + ecce),
        ];
        apg = crate::math::polar_to_cartesian(apg);
        apg = rotate_x_sincos(apg, -sinincl, cosincl);
        apg = crate::math::cartesian_to_polar(apg);
        apg[0] += sinnode.atan2(cosnode);
        let apg_cart = crate::math::polar_to_cartesian(apg);
        xxa[i][0..3].copy_from_slice(&apg_cart);

        // ellipse-corrected node distance (reusing this sample's ecce/sema/uu)
        ny = crate::math::normalize_radians(ny - uu);
        cos_e = (2.0 * ((ny / 2.0).tan() / ((1.0 + ecce) / (1.0 - ecce)).sqrt()).atan()).cos();
        let r0 = sema * (1.0 - ecce * cos_e);
        let r1 = (xx[i][0] * xx[i][0] + xx[i][1] * xx[i][1] + xx[i][2] * xx[i][2]).sqrt();
        for v in &mut xx[i][..3] {
            *v *= r0 / r1;
        }
    }

    // Save node/apogee position + central-difference speed (both plain central
    // differences — the D.2 quadratic node speed is dead code, overwritten here).
    let mut node = [0.0f64; 6];
    let mut apog = [0.0f64; 6];
    for i in 0..3 {
        apog[i] = xxa[2][i];
        node[i] = xx[2][i];
        if has_speed {
            apog[i + 3] = (xxa[1][i] - xxa[0][i]) / speed_intv / 2.0;
            node[i + 3] = (xx[1][i] - xx[0][i]) / speed_intv / 2.0;
        }
    }

    // D.4: output frame for both bodies. oe = obliquity of date.
    let oe = obliquity(tjd, flags, models);
    (
        osc_output_frame(&node, flags, &oe, tjd, models),
        osc_output_frame(&apog, flags, &oe, tjd, models),
    )
}

/// Raw geocentric equatorial-J2000 (pre-bias) Moshier moon (pos+vel) at `t`,
/// as C's `swi_moshmoon` returns it before `swi_plan_for_osc_elem`.
pub(crate) fn raw_osc_moon_moshier(t: f64, eps_j2000: &Epsilon) -> Result<[f64; 6], Error> {
    crate::moshier::backend::compute(t, Body::Moon, eps_j2000)
}

/// Raw geocentric equatorial-J2000 (pre-bias) Swiss-ephemeris moon (pos+vel) at `t`.
#[cfg(feature = "swisseph-files")]
pub(crate) fn raw_osc_moon_sweph(moon_files: &[SwissEphFile], t: f64) -> Result<[f64; 6], Error> {
    SwephProvider {
        planet_files: &[],
        moon_files,
    }
    .moon_geo(t, true)
}

/// Raw geocentric equatorial-J2000 (pre-bias) JPL moon (pos+vel) at `t`.
#[cfg(feature = "jpl")]
pub(crate) fn raw_osc_moon_jpl(file: &crate::jpl::JplFile, t: f64) -> Result<[f64; 6], Error> {
    JplProvider { file }.moon_geo(t, true)
}

/// Apply default-branch sidereal projection (Branch 3) to `xreturn`.
///
/// `daya[0]` is the ayanamsa in degrees; `daya[1]` is the ayanamsa speed in
/// degrees/day. Adjusts ecliptic polar `[0..6]` and recomputes ecliptic
/// Cartesian `[6..12]`. Leaves equatorial `[12..24]` untouched (matches C).
pub fn apply_sidereal_default(xreturn: &mut [f64; 24], daya: [f64; 2], has_speed: bool) {
    xreturn[0] = crate::math::normalize_degrees(xreturn[0] - daya[0]);
    if has_speed {
        xreturn[3] -= daya[1];
    }
    let polar = [
        xreturn[0] * DEGTORAD,
        xreturn[1] * DEGTORAD,
        xreturn[2],
        xreturn[3] * DEGTORAD,
        xreturn[4] * DEGTORAD,
        xreturn[5],
    ];
    let cart = crate::math::polar_to_cartesian_with_speed(polar);
    xreturn[6..12].copy_from_slice(&cart);
}

// ---------------------------------------------------------------------------
// SwissEph (.se1) backend
// ---------------------------------------------------------------------------

#[cfg_attr(
    not(any(feature = "swisseph-files", feature = "jpl")),
    allow(dead_code)
)]
pub(crate) struct SwephPositions {
    pub(crate) planet_bary: [f64; 6],
    pub(crate) earth_bary: [f64; 6],
    pub(crate) earth_helio: [f64; 6],
    pub(crate) sun_bary: [f64; 6],
}

#[cfg(feature = "swisseph-files")]
pub(crate) fn sweph_positions(
    planet_file: &SwissEphFile,
    moon_file: &SwissEphFile,
    body_id: i32,
    jd: f64,
    need_speed: bool,
) -> Result<SwephPositions, Error> {
    let n = if need_speed { 6 } else { 3 };

    let (emb, _) = evaluate_body(planet_file, 0, jd, need_speed)?;
    let (moon_geo, _) = evaluate_body(moon_file, SEI_MOON, jd, need_speed)?;

    let mut earth_bary = [0.0; 6];
    for i in 0..n {
        earth_bary[i] = emb[i] - moon_geo[i] / (EARTH_MOON_MRAT + 1.0);
    }

    let (helio_earth, _) = evaluate_body(planet_file, SEI_SUNBARY, jd, need_speed)?;
    let mut sun_bary = [0.0; 6];
    for i in 0..n {
        sun_bary[i] = emb[i] - helio_earth[i];
    }

    let mut earth_helio = [0.0; 6];
    for i in 0..n {
        earth_helio[i] = earth_bary[i] - sun_bary[i];
    }

    let (mut planet, _) = evaluate_body(planet_file, body_id, jd, need_speed)?;
    if let Some(pd) = planet_file.planet_data(body_id)
        && pd.iflg & SEI_FLG_HELIO != 0
    {
        for i in 0..n {
            planet[i] += sun_bary[i];
        }
    }

    Ok(SwephPositions {
        planet_bary: planet,
        earth_bary,
        earth_helio,
        sun_bary,
    })
}

#[cfg_attr(
    not(any(feature = "swisseph-files", feature = "jpl")),
    allow(dead_code)
)]
pub(crate) trait PositionProvider {
    /// Barycentric equatorial-J2000 positions of `body`, Earth, and Sun at `jd`.
    fn positions(&self, body: Body, jd: f64, need_speed: bool) -> Result<SwephPositions, Error>;
    /// Geocentric equatorial-J2000 Moon at `jd`.
    fn moon_geo(&self, jd: f64, need_speed: bool) -> Result<[f64; 6], Error>;
    /// Whether `positions()` returns sun_bary freshly evaluated at the requested
    /// epoch. Swiss (sweplan) batch-fetches sun_bary alongside earth; JPL
    /// (swi_pleph) only fetches the target body. Affects Earth HELCTR light-time.
    fn updates_sun_in_light_time(&self) -> bool {
        true
    }
}

#[cfg(feature = "swisseph-files")]
pub(crate) struct SwephProvider<'a> {
    pub(crate) planet_files: &'a [SwissEphFile],
    pub(crate) moon_files: &'a [SwissEphFile],
}

/// Days of slop allowed past a file's per-body coverage before declaring a jd
/// out of range. Larger than the maximum one-way light-time in the solar system
/// (Pluto ~0.23 d) but far smaller than a file's century-scale span.
#[cfg(feature = "swisseph-files")]
const BOUNDARY_SLOP: f64 = 0.3;

/// Select the `.se1` file for `(body_id, jd)`. Falls back to the nearest covering
/// file when no file's `file_start <= jd` window matches but `jd` lies within
/// `BOUNDARY_SLOP` of a file's per-body range.
///
/// This restores the pre-refactor `find_file_for_jd(..).unwrap_or(primary_file)`
/// behavior used for retarded light-time epochs: when `jd - dt` falls just past
/// the absolute earliest ephemeris boundary, C extrapolates the cached segment
/// rather than erroring. Genuine out-of-range epochs (more than a light-time off)
/// still return `None`, preserving the `BeyondEphemerisLimits` contract.
#[cfg(feature = "swisseph-files")]
pub(crate) fn find_file_or_nearest(
    files: &[SwissEphFile],
    body_id: i32,
    jd: f64,
) -> Option<&SwissEphFile> {
    if let Some(f) = find_file_for_jd(files, body_id, jd) {
        return Some(f);
    }
    let mut best: Option<(&SwissEphFile, f64)> = None;
    for f in files {
        let Some(pd) = f.planet_data(body_id) else {
            continue;
        };
        let dist = if jd < pd.tfstart {
            pd.tfstart - jd
        } else if jd > pd.tfend {
            jd - pd.tfend
        } else {
            0.0
        };
        if dist <= BOUNDARY_SLOP && best.is_none_or(|(_, d)| dist < d) {
            best = Some((f, dist));
        }
    }
    best.map(|(f, _)| f)
}

#[cfg(feature = "swisseph-files")]
impl<'a> PositionProvider for SwephProvider<'a> {
    fn positions(&self, body: Body, jd: f64, need_speed: bool) -> Result<SwephPositions, Error> {
        let body_id = body_file_id(body).ok_or(Error::EphemerisNotAvailable {
            body,
            source: EphemerisSource::Swiss,
        })?;
        let planet_file = find_file_or_nearest(self.planet_files, body_id, jd).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd,
                start: 0.0,
                end: 0.0,
            },
        )?;
        let moon_file = find_file_or_nearest(self.moon_files, SEI_MOON, jd).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd,
                start: 0.0,
                end: 0.0,
            },
        )?;
        sweph_positions(planet_file, moon_file, body_id, jd, need_speed)
    }

    fn moon_geo(&self, jd: f64, need_speed: bool) -> Result<[f64; 6], Error> {
        let moon_file = find_file_or_nearest(self.moon_files, SEI_MOON, jd).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd,
                start: 0.0,
                end: 0.0,
            },
        )?;
        let (pos, _) = evaluate_body(moon_file, SEI_MOON, jd, need_speed)?;
        Ok(pos)
    }
}

#[cfg_attr(
    not(any(feature = "swisseph-files", feature = "jpl")),
    allow(dead_code)
)]
fn apparent_planet<P: PositionProvider>(
    p: &P,
    jd: f64,
    body: Body,
    _eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let need_speed = flags.contains(CalcFlags::SPEED);

    let pos = p.positions(body, jd, true)?;

    // Heliocentric (Swiss/JPL): heliocentric = planet_bary - sun_bary, light-time retarded with
    // niter=1 (two analytic passes to converge dt, then re-evaluate the ephemeris at t-dt —
    // sweph.c:2648-2696). No aberration/deflection/geocenter conversion.
    if flags.contains(CalcFlags::HELCTR) {
        let mut planet_helio = [0.0; 6];
        for (i, ph) in planet_helio.iter_mut().enumerate() {
            *ph = pos.planet_bary[i] - pos.sun_bary[i];
        }
        let mut xx = planet_helio;
        if !flags.contains(CalcFlags::TRUEPOS) {
            // C quirk (sweph.c:2513-2594): the light-time loop's initial dx is the *heliocentric*
            // position (so the first dt uses the Sun-planet distance), but the analytic
            // extrapolation base `xx0` is the *barycentric* position/velocity (xx0 is saved
            // before the SUNBARY subtraction). So the converged dt comes from the barycentric
            // distance, not the heliocentric one -- a ~1% difference that matters at 1e-9.
            let xx0 = pos.planet_bary;
            let mut dt = 0.0;
            for _ in 0..=1 {
                let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
                dt = dist * AUNIT / CLIGHT / 86400.0;
                for i in 0..3 {
                    xx[i] = xx0[i] - dt * xx0[i + 3];
                }
            }
            // C re-evaluates the planet at t-dt with NO_SAVE, so the cached SEI_SUNBARY it then
            // subtracts is still the ORIGINAL-epoch Sun (sweph.c:2692-2695) -- not the retarded
            // one. Mixed-epoch heliocentric: planet_bary(t-dt) - sun_bary(t).
            let pos_ret = p.positions(body, jd - dt, true)?;
            for (i, x) in xx.iter_mut().enumerate() {
                *x = pos_ret.planet_bary[i] - pos.sun_bary[i];
            }
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    // Barycentric (Swiss/JPL): planet_bary directly, light-time retarded with niter=1.
    // Structurally identical to HELCTR but without the sun_bary subtraction — C's
    // app_pos_etc_plan gates the subtraction on HELCTR only (sweph.c:2516-2520, 2692-2696).
    if flags.contains(CalcFlags::BARYCTR) {
        let mut xx = pos.planet_bary;
        if !flags.contains(CalcFlags::TRUEPOS) {
            let xx0 = pos.planet_bary;
            let mut dt = 0.0;
            for _ in 0..=1 {
                let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
                dt = dist * AUNIT / CLIGHT / 86400.0;
                for i in 0..3 {
                    xx[i] = xx0[i] - dt * xx0[i + 3];
                }
            }
            let pos_ret = p.positions(body, jd - dt, true)?;
            xx = pos_ret.planet_bary;
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    // Observer: barycentric Earth, or Earth + topocentric offset (§1 "xobs
    // replaces the geocenter"), evaluated once at the current (un-retarded)
    // epoch. `xobs_helio` is the same offset applied to the heliocentric frame
    // for the deflection geometry (§7).
    let offset = topo_offset(jd, flags, config, models);
    let mut xobs = pos.earth_bary;
    for i in 0..6 {
        xobs[i] += offset[i];
    }
    let mut xobs_helio = pos.earth_helio;
    for i in 0..6 {
        xobs_helio[i] += offset[i];
    }

    // Geocentric (or topocentric)
    let mut xx = [0.0; 6];
    for (i, x) in xx.iter_mut().enumerate() {
        *x = pos.planet_bary[i] - xobs[i];
    }

    // Light-time with niter=1
    let mut dt = 0.0;
    if !flags.contains(CalcFlags::TRUEPOS) {
        let xxsv = pos.planet_bary;

        // Speed correction: compute light-time at t-1 day
        let mut xxsp = [xxsv[0] - xxsv[3], xxsv[1] - xxsv[4], xxsv[2] - xxsv[5]];
        let xxsv_sp = xxsp;

        for _ in 0..=1 {
            let dx = [
                xxsp[0] - (xobs[0] - xobs[3]),
                xxsp[1] - (xobs[1] - xobs[4]),
                xxsp[2] - (xobs[2] - xobs[5]),
            ];
            let dist_sp = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
            let dt_sp = dist_sp * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xxsp[i] = xxsv_sp[i] - dt_sp * xxsv[i + 3];
            }
        }
        // Δ@(t-1) = true@(t-1) - apparent@(t-1)
        for i in 0..3 {
            xxsp[i] = xxsv_sp[i] - xxsp[i];
        }

        // Main light-time loop
        for _ in 0..=1 {
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            dt = dist * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xx[i] = xxsv[i] - dt * xxsv[i + 3] - xobs[i];
            }
        }
        // Change-of-dt correction = Δ@t - Δ@(t-1)
        for i in 0..3 {
            xxsp[i] = dt * xxsv[i + 3] - xxsp[i];
        }

        // Re-evaluate planet at retarded time t' = t - dt; file selection is
        // inside the provider (handles file-boundary retarded-epoch lookup).
        let pos_ret = p.positions(body, jd - dt, true)?;

        // Geocentric from re-evaluated position
        for (i, x) in xx.iter_mut().enumerate() {
            *x = pos_ret.planet_bary[i] - xobs[i];
        }

        // Apply change-of-dt speed correction
        for i in 0..3 {
            xx[i + 3] -= xxsp[i];
        }
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // Deflection
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOGDEFL) {
        let mut planet_helio_retarded = [0.0; 6];
        for i in 0..3 {
            planet_helio_retarded[i] = xx[i] + xobs_helio[i];
            planet_helio_retarded[i + 3] = pos.planet_bary[i + 3];
        }
        deflect_light(&mut xx, &xobs_helio, &planet_helio_retarded, need_speed);
    }

    // Aberration
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(&mut xx, &[xobs[3], xobs[4], xobs[5]], need_speed);
        if need_speed {
            // Observer at the retarded epoch — an independent offset (§4), not
            // derived from `offset` at the current epoch.
            let pos_ret = p.positions(body, jd - dt, true)?;
            let offset_ret = topo_offset(jd - dt, flags, config, models);
            let mut xobs2 = pos_ret.earth_bary;
            for i in 0..6 {
                xobs2[i] += offset_ret[i];
            }
            for i in 0..3 {
                xx[i + 3] += xobs[i + 3] - xobs2[i + 3];
            }
        }
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }

    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);
    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

#[cfg_attr(
    not(any(feature = "swisseph-files", feature = "jpl")),
    allow(dead_code)
)]
fn apparent_sun<P: PositionProvider>(
    p: &P,
    jd: f64,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
    is_earth: bool,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let need_speed = flags.contains(CalcFlags::SPEED);

    // Always pass need_speed=true internally — velocity needed for aberration
    // even when the caller doesn't request speed output.
    let pos = p.positions(Body::Sun, jd, true)?;

    // Observer, evaluated once at the current epoch (§6 — no retarded-time
    // xobs2 term for the Sun, unlike the planet/Moon paths).
    let offset = topo_offset(jd, flags, config, models);
    let mut xobs = pos.earth_bary;
    for i in 0..6 {
        xobs[i] += offset[i];
    }

    let is_hb = flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR);
    let is_bary = flags.contains(CalcFlags::BARYCTR);

    // Barycentric Sun (C's app_pos_etc_sbar, sweph.c:4254-4293): returns sun_bary directly.
    // C handles this in a completely separate swecalc branch (sweph.c:733-815) because
    // SEI_SUN == SEI_EARTH in C's internal indices — Rust doesn't have that aliasing.
    // Single analytic light-time retardation (no loop, no re-eval), then finish_helctr.
    if !is_earth && is_bary {
        let mut xx = pos.sun_bary;
        if !flags.contains(CalcFlags::TRUEPOS) {
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xx[i] -= dt * xx[i + 3];
            }
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    // Frame construction (sweph.c:3944-3950):
    //   BARYCTR      → xx = xobs (earth_bary, the barycentric Earth directly)
    //   HELCTR       → xx = xobs - sun_bary (earth_bary - sun_bary = helio_earth)
    //   geocentric   → xx = -(xobs - sun_bary) (geo_sun = -helio_earth)
    let mut xx = [0.0; 6];
    if is_earth && is_bary {
        for (i, x) in xx.iter_mut().enumerate() {
            *x = xobs[i];
        }
    } else if is_earth && is_hb {
        for (i, x) in xx.iter_mut().enumerate() {
            *x = xobs[i] - pos.sun_bary[i];
        }
    } else {
        for (i, x) in xx.iter_mut().enumerate() {
            *x = -(xobs[i] - pos.sun_bary[i]);
        }
    }

    // Light-time (sweph.c:3968-4030). The loop always uses niter=1 (0..=1 = 2
    // iterations). For HELCTR/BARYCTR Earth, re-evaluate Earth at retarded time;
    // for geocentric Sun, re-evaluate Sun at retarded time.
    // For HELCTR Earth: C's Swiss path (sweplan) batch-fetches both xearth and
    // xsun at retarded epoch; C's JPL path (swi_pleph) only re-fetches xearth,
    // leaving xsun at the original epoch. We match per-backend via
    // updates_sun_in_light_time.
    if !flags.contains(CalcFlags::TRUEPOS) {
        let orig_sun_bary = pos.sun_bary;
        for _ in 0..=1 {
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            let pos_ret = p.positions(Body::Sun, jd - dt, true)?;
            if is_earth && is_hb {
                if is_bary {
                    for (i, x) in xx.iter_mut().enumerate() {
                        *x = pos_ret.earth_bary[i] + offset[i];
                    }
                } else {
                    let sun_bary = if p.updates_sun_in_light_time() {
                        &pos_ret.sun_bary
                    } else {
                        &orig_sun_bary
                    };
                    for (i, x) in xx.iter_mut().enumerate() {
                        *x = pos_ret.earth_bary[i] + offset[i] - sun_bary[i];
                    }
                }
            } else {
                for (i, x) in xx.iter_mut().enumerate() {
                    *x = -(xobs[i] - pos_ret.sun_bary[i]);
                }
            }
        }
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // No deflection for Sun (or Earth through this path)

    // Aberration (skipped for HELCTR/BARYCTR via plaus_iflag's NOABERR)
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(&mut xx, &[xobs[3], xobs[4], xobs[5]], need_speed);
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }

    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);
    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

#[cfg_attr(
    not(any(feature = "swisseph-files", feature = "jpl")),
    allow(dead_code)
)]
fn apparent_moon<P: PositionProvider>(
    p: &P,
    jd: f64,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let need_speed = flags.contains(CalcFlags::SPEED);

    // Earth context for aberration (body=Earth maps to EMB file entry, same as
    // the original sweph_positions call with body_id=0)
    let pos = p.positions(Body::Earth, jd, true)?;

    // Heliocentric Moon (Swiss/JPL): heliocentric Moon at jd = moon_geo + earth_bary - sun_bary.
    // Light-time dt is computed ONCE from this heliocentric distance (sweph.c:4138-4152, `xxm`),
    // then the ephemeris is re-evaluated at t-dt to barycentric moon(t-dt) = moon_geo(t-dt) +
    // earth_bary(t-dt), minus the ORIGINAL-epoch Sun (xobs = SUNBARY(teval), cached across the
    // NO_SAVE re-eval — sweph.c:4168-4206). No deflection/aberration.
    if flags.contains(CalcFlags::HELCTR) {
        let moon_geo = p.moon_geo(jd, true)?;
        let mut moon_helio = [0.0; 6];
        for i in 0..6 {
            moon_helio[i] = moon_geo[i] + pos.earth_bary[i] - pos.sun_bary[i];
        }
        let mut xx = moon_helio;
        if !flags.contains(CalcFlags::TRUEPOS) {
            let dist = (moon_helio[0] * moon_helio[0]
                + moon_helio[1] * moon_helio[1]
                + moon_helio[2] * moon_helio[2])
                .sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            let pos_ret = p.positions(Body::Earth, jd - dt, true)?;
            let moon_geo_ret = p.moon_geo(jd - dt, true)?;
            for i in 0..6 {
                xx[i] = moon_geo_ret[i] + pos_ret.earth_bary[i] - pos.sun_bary[i];
            }
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    // Barycentric Moon (Swiss/JPL): moon_geo + earth_bary. No sun_bary subtraction — C's
    // app_pos_etc_moon sets xobs=0 for BARYCTR (sweph.c:4133-4137), so the "to correct
    // center" step at sweph.c:4205 is a no-op. Single dt from |moon_bary|, same as HELCTR.
    if flags.contains(CalcFlags::BARYCTR) {
        let moon_geo = p.moon_geo(jd, true)?;
        let mut moon_bary = [0.0; 6];
        for i in 0..6 {
            moon_bary[i] = moon_geo[i] + pos.earth_bary[i];
        }
        let mut xx = moon_bary;
        if !flags.contains(CalcFlags::TRUEPOS) {
            let dist = (moon_bary[0] * moon_bary[0]
                + moon_bary[1] * moon_bary[1]
                + moon_bary[2] * moon_bary[2])
                .sqrt();
            let dt = dist * AUNIT / CLIGHT / 86400.0;
            let pos_ret = p.positions(Body::Earth, jd - dt, true)?;
            let moon_geo_ret = p.moon_geo(jd - dt, true)?;
            for i in 0..6 {
                xx[i] = moon_geo_ret[i] + pos_ret.earth_bary[i];
            }
        }
        return Ok(finish_helctr(xx, jd, flags, models));
    }

    // Observer, evaluated once at the current (un-retarded) epoch.
    let offset = topo_offset(jd, flags, config, models);
    let mut xobs = pos.earth_bary;
    for i in 0..6 {
        xobs[i] += offset[i];
    }

    // Moon geocentric from provider, shifted to topocentric (§5 — the extra
    // `xxm[i] -= xobs[i]` line that doesn't appear in the planet/Sun cases).
    let mut xx = p.moon_geo(jd, true)?;
    for i in 0..6 {
        xx[i] -= offset[i];
    }

    // Light-time
    let mut xobs2 = xobs;
    if !flags.contains(CalcFlags::TRUEPOS) {
        let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
        let dt = dist * AUNIT / CLIGHT / 86400.0;
        let pos_ret = p.positions(Body::Earth, jd - dt, true)?;
        let offset_ret = topo_offset(jd - dt, flags, config, models);
        xobs2 = pos_ret.earth_bary;
        for i in 0..6 {
            xobs2[i] += offset_ret[i];
        }
        let moon_geo_ret = p.moon_geo(jd - dt, true)?;

        // moon_bary(t') = moon_geo(t') + earth_bary(t')
        // topocentric at retarded time: moon_bary(t') - xobs(t)
        for i in 0..6 {
            xx[i] = moon_geo_ret[i] + pos_ret.earth_bary[i] - xobs[i];
        }
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // No deflection for Moon

    // Aberration
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(&mut xx, &[xobs[3], xobs[4], xobs[5]], need_speed);
        if need_speed {
            for i in 0..3 {
                xx[i + 3] += xobs[i + 3] - xobs2[i + 3];
            }
        }
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
    }

    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);
    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

// ---------------------------------------------------------------------------
// JPL DE backend
// ---------------------------------------------------------------------------

#[cfg(feature = "jpl")]
fn body_to_jpl_index(body: Body) -> Option<i32> {
    use crate::jpl::{
        J_EARTH, J_JUPITER, J_MARS, J_MERCURY, J_MOON, J_NEPTUNE, J_PLUTO, J_SATURN, J_SUN,
        J_URANUS, J_VENUS,
    };
    match body {
        Body::Sun => Some(J_SUN),
        Body::Moon => Some(J_MOON),
        Body::Earth => Some(J_EARTH),
        Body::Mercury => Some(J_MERCURY),
        Body::Venus => Some(J_VENUS),
        Body::Mars => Some(J_MARS),
        Body::Jupiter => Some(J_JUPITER),
        Body::Saturn => Some(J_SATURN),
        Body::Uranus => Some(J_URANUS),
        Body::Neptune => Some(J_NEPTUNE),
        Body::Pluto => Some(J_PLUTO),
        _ => None,
    }
}

#[cfg(feature = "jpl")]
pub(crate) struct JplProvider<'a> {
    pub(crate) file: &'a crate::jpl::JplFile,
}

#[cfg(feature = "jpl")]
impl<'a> PositionProvider for JplProvider<'a> {
    fn positions(&self, body: Body, jd: f64, need_speed: bool) -> Result<SwephPositions, Error> {
        use crate::jpl::{J_EARTH, J_SBARY, J_SUN, jpl_pleph};
        let j_target = body_to_jpl_index(body).ok_or(Error::EphemerisNotAvailable {
            body,
            source: EphemerisSource::Jpl,
        })?;
        let planet_bary = jpl_pleph(self.file, jd, j_target, J_SBARY, need_speed)?;
        let earth_bary = jpl_pleph(self.file, jd, J_EARTH, J_SBARY, need_speed)?;
        let sun_bary = jpl_pleph(self.file, jd, J_SUN, J_SBARY, need_speed)?;
        let mut earth_helio = [0.0f64; 6];
        for i in 0..6 {
            earth_helio[i] = earth_bary[i] - sun_bary[i];
        }
        Ok(SwephPositions {
            planet_bary,
            earth_bary,
            earth_helio,
            sun_bary,
        })
    }

    fn moon_geo(&self, jd: f64, need_speed: bool) -> Result<[f64; 6], Error> {
        use crate::jpl::{J_EARTH, J_MOON, jpl_pleph};
        jpl_pleph(self.file, jd, J_MOON, J_EARTH, need_speed)
    }

    fn updates_sun_in_light_time(&self) -> bool {
        false
    }
}

/// Computes `body`'s apparent position using the JPL DE ephemeris `file`,
/// running the shared light-time/aberration/precession pipeline.
#[cfg(feature = "jpl")]
pub fn calc_planet_jpl(
    jd: f64,
    body: Body,
    file: &crate::jpl::JplFile,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = JplProvider { file };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

/// Computes the Sun's (or, with `is_earth`, Earth's) apparent position using
/// the JPL DE ephemeris `file`.
#[cfg(feature = "jpl")]
pub fn calc_sun_jpl(
    jd: f64,
    file: &crate::jpl::JplFile,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
    is_earth: bool,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = JplProvider { file };
    apparent_sun(&p, jd, flags, config, models, is_earth)
}

/// Computes the Moon's apparent position using the JPL DE ephemeris `file`.
#[cfg(feature = "jpl")]
pub fn calc_moon_jpl(
    jd: f64,
    file: &crate::jpl::JplFile,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = JplProvider { file };
    apparent_moon(&p, jd, flags, config, models)
}

/// Computes `body`'s apparent position using the Swiss Ephemeris (.se1) file
/// backend, selecting the appropriate `planet_files`/`moon_files` entries.
#[cfg(feature = "swisseph-files")]
#[allow(clippy::too_many_arguments)]
pub fn calc_planet_sweph(
    jd: f64,
    body: Body,
    planet_files: &[SwissEphFile],
    moon_files: &[SwissEphFile],
    _eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = SwephProvider {
        planet_files,
        moon_files,
    };
    apparent_planet(&p, jd, body, _eps_j2000, flags, config, models)
}

/// Computes the Sun's (or, with `is_earth`, Earth's) apparent position using
/// the Swiss Ephemeris (.se1) file backend.
#[cfg(feature = "swisseph-files")]
pub fn calc_sun_sweph(
    jd: f64,
    planet_files: &[SwissEphFile],
    moon_files: &[SwissEphFile],
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
    is_earth: bool,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = SwephProvider {
        planet_files,
        moon_files,
    };
    apparent_sun(&p, jd, flags, config, models, is_earth)
}

/// Computes the Moon's apparent position using the Swiss Ephemeris (.se1) file
/// backend.
#[cfg(feature = "swisseph-files")]
pub fn calc_moon_sweph(
    jd: f64,
    planet_files: &[SwissEphFile],
    moon_files: &[SwissEphFile],
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = SwephProvider {
        planet_files,
        moon_files,
    };
    apparent_moon(&p, jd, flags, config, models)
}

/// Light-time iteration for swe_calc_pctr (c-ref-pctr §3a–§3c).
/// Returns `(retarded_time, dt, xxsp)`.
pub(crate) fn pctr_light_time(
    tjd: f64,
    xx0: &[f64; 6],
    xxctr: &[f64; 6],
    has_speed: bool,
) -> (f64, f64, [f64; 3]) {
    let niter = 1;
    let mut xxsp = [0.0; 3];

    // §3a: SPEED pre-pass — "change of dt" correction seed
    if has_speed {
        let mut xxsv = [0.0; 3];
        for i in 0..3 {
            xxsv[i] = xx0[i] - xx0[i + 3];
            xxsp[i] = xxsv[i];
        }
        for _ in 0..=niter {
            let mut dx = [0.0; 3];
            for i in 0..3 {
                dx[i] = xxsp[i] - (xxctr[i] - xxctr[i + 3]);
            }
            let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
            let dt_sp = dist * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xxsp[i] = xxsv[i] - dt_sp * xx0[i + 3];
            }
        }
        for i in 0..3 {
            xxsp[i] = xxsv[i] - xxsp[i];
        }
    }

    // §3b: Main light-time loop
    let mut xx_pos = [xx0[0], xx0[1], xx0[2]];
    let mut dt = 0.0;
    for _ in 0..=niter {
        let mut dx = [0.0; 3];
        for i in 0..3 {
            dx[i] = xx_pos[i] - xxctr[i];
        }
        let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        dt = dist * AUNIT / CLIGHT / 86400.0;
        for i in 0..3 {
            xx_pos[i] = xx0[i] - dt * xx0[i + 3];
        }
    }
    let t = tjd - dt;

    // §3c: Finalize speed correction
    if has_speed {
        for i in 0..3 {
            xxsp[i] = xx0[i] - xx_pos[i] - xxsp[i];
        }
    }

    (t, dt, xxsp)
}

/// Pipeline for swe_calc_pctr §4–§9: planetocentric subtraction, deflection,
/// aberration, frame bias, precession, app_pos_rest.
///
/// `nut_epoch` is the §1 priming epoch (tjd + Δt(tjd)) at which `eps`/`nut`
/// were computed; needed for the nutation speed derivative (nutv).
#[allow(clippy::too_many_arguments)]
pub(crate) fn pctr_pipeline(
    xx_ipl: &[f64; 6],
    xxctr: &[f64; 6],
    xxctr2: &[f64; 6],
    xxsp: &[f64; 3],
    t: f64,
    tjd: f64,
    nut_epoch: f64,
    earth_bary: &[f64; 6],
    sun_bary: &[f64; 6],
    flags: CalcFlags,
    eps: &Epsilon,
    nut: &NutationType,
    models: &AstroModels,
) -> ([f64; 24], [f64; 6]) {
    let has_speed = flags.contains(CalcFlags::SPEED);

    // §4: Planetocentric subtraction (unconditional — HELCTR/BARYCTR stripped)
    let mut xx = [0.0; 6];
    for i in 0..6 {
        xx[i] = xx_ipl[i] - xxctr[i];
    }
    if !flags.contains(CalcFlags::TRUEPOS) && has_speed {
        for i in 0..3 {
            xx[i + 3] -= xxsp[i];
        }
    }
    if !has_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // §5: Gravitational deflection (Earth-observer geometry regardless of center)
    // C reads pedp->x (earth_bary) and computes q = xx + earth_bary. Our
    // deflect_light is designed around earth_helio (standard pipeline convention),
    // so we compute earth_helio = earth_bary - sun_bary, and planet_for_defl =
    // xx + earth_helio (self-consistent: in the standard pipeline planet_helio =
    // geocentric + earth_helio). This avoids NaN when ipl=Sun (planet_helio would
    // be [0;6]) while staying consistent with deflect_light's speed perturbation.
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOGDEFL) {
        let mut earth_helio = [0.0; 6];
        let mut planet_for_defl = [0.0; 6];
        for i in 0..6 {
            earth_helio[i] = earth_bary[i] - sun_bary[i];
            planet_for_defl[i] = xx[i] + earth_helio[i];
        }
        deflect_light(&mut xx, &earth_helio, &planet_for_defl, has_speed);
    }

    // §6: Annual aberration (center body velocity as observer)
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        let ctr_vel = [xxctr[3], xxctr[4], xxctr[5]];
        aberr_light(&mut xx, &ctr_vel, has_speed);
        if has_speed {
            for i in 3..6 {
                xx[i] += xxctr[i] - xxctr2[i];
            }
        }
    }
    if !has_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // §7: ICRS → J2000 frame bias (at retarded time t)
    if !flags.contains(CalcFlags::ICRS) {
        frame_bias(&mut xx, t, flags, models, FrameTransform::GcrsToJ2000);
    }

    // Save J2000 coordinates for sidereal projection
    let x2000 = xx;

    // §8–§9: Precession (at tjd) + nutation/ecliptic/polar/degrees
    if !flags.contains(CalcFlags::J2000) {
        let mut pos3 = [xx[0], xx[1], xx[2]];
        precess(
            &mut pos3,
            tjd,
            flags,
            models,
            PrecessionDirection::J2000ToDate,
        );
        xx[0] = pos3[0];
        xx[1] = pos3[1];
        xx[2] = pos3[2];
        if has_speed {
            precess_speed(
                &mut xx,
                tjd,
                flags,
                models,
                PrecessionDirection::J2000ToDate,
            );
        }
    }

    let nutv = if has_speed && !flags.contains(CalcFlags::J2000) {
        Some(nutation(nut_epoch - NUT_SPEED_INTV, flags, models))
    } else {
        None
    };

    (app_pos_rest(&mut xx, flags, eps, nut, nutv.as_ref()), x2000)
}

// ---------------------------------------------------------------------------
// Asteroid calc pipeline
// ---------------------------------------------------------------------------

#[cfg(feature = "swisseph-files")]
pub(crate) struct AsteroidProvider<'a, P: PositionProvider> {
    inner: &'a P,
    ast_file: &'a SwissEphFile,
    ast_id: i32,
    body: Body,
    source: EphemerisSource,
}

#[cfg(feature = "swisseph-files")]
impl<'a, P: PositionProvider> PositionProvider for AsteroidProvider<'a, P> {
    fn positions(&self, _body: Body, jd: f64, need_speed: bool) -> Result<SwephPositions, Error> {
        let pos = self.inner.positions(Body::Sun, jd, need_speed)?;
        let n = if need_speed { 6 } else { 3 };
        let (mut ast, _) =
            evaluate_body(self.ast_file, self.ast_id, jd, need_speed).map_err(|e| match e {
                Error::BeyondEphemerisLimits { .. } | Error::InvalidBody(_) => {
                    Error::EphemerisNotAvailable {
                        body: self.body,
                        source: self.source,
                    }
                }
                other => other,
            })?;
        // C's sweph() adds sun_bary unconditionally for ipl >= SEI_ANYBODY (slot-index
        // check, sweph.c:2332-2343) — the file's SEI_FLG_HELIO flag is NOT checked for
        // asteroids (seas files don't set it, even though their data is heliocentric).
        // MoshierEarthProvider returns sun_bary=[0;6], making this a no-op under MOSEPH
        // — matching C's flag guard that skips the add when !(SWIEPH|JPLEPH).
        for (a, &s) in ast.iter_mut().zip(pos.sun_bary.iter()).take(n) {
            *a += s;
        }
        Ok(SwephPositions {
            planet_bary: ast,
            earth_bary: pos.earth_bary,
            earth_helio: pos.earth_helio,
            sun_bary: pos.sun_bary,
        })
    }

    fn moon_geo(&self, jd: f64, need_speed: bool) -> Result<[f64; 6], Error> {
        self.inner.moon_geo(jd, need_speed)
    }

    fn updates_sun_in_light_time(&self) -> bool {
        self.inner.updates_sun_in_light_time()
    }
}

pub(crate) struct MoshierEarthProvider<'a> {
    pub(crate) eps_j2000: &'a Epsilon,
}

impl<'a> PositionProvider for MoshierEarthProvider<'a> {
    fn positions(&self, _body: Body, jd: f64, _need_speed: bool) -> Result<SwephPositions, Error> {
        let pp = compute_pipeline(jd, Body::Sun, self.eps_j2000)?;
        Ok(SwephPositions {
            planet_bary: [0.0; 6],
            earth_bary: pp.earth_helio,
            earth_helio: pp.earth_helio,
            sun_bary: [0.0; 6],
        })
    }

    fn moon_geo(&self, jd: f64, _need_speed: bool) -> Result<[f64; 6], Error> {
        raw_osc_moon_moshier(jd, self.eps_j2000)
    }
}

#[cfg(feature = "swisseph-files")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_asteroid_sweph(
    jd: f64,
    body: Body,
    ast_file: &SwissEphFile,
    ast_id: i32,
    planet_files: &[SwissEphFile],
    moon_files: &[SwissEphFile],
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let inner = SwephProvider {
        planet_files,
        moon_files,
    };
    let p = AsteroidProvider {
        inner: &inner,
        ast_file,
        ast_id,
        body,
        source: EphemerisSource::Swiss,
    };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

#[cfg(all(feature = "swisseph-files", feature = "jpl"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_asteroid_jpl(
    jd: f64,
    body: Body,
    ast_file: &SwissEphFile,
    ast_id: i32,
    jpl_file: &crate::jpl::JplFile,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let inner = JplProvider { file: jpl_file };
    let p = AsteroidProvider {
        inner: &inner,
        ast_file,
        ast_id,
        body,
        source: EphemerisSource::Jpl,
    };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

#[cfg(feature = "swisseph-files")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_asteroid_moshier(
    jd: f64,
    body: Body,
    ast_file: &SwissEphFile,
    ast_id: i32,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let inner = MoshierEarthProvider { eps_j2000 };
    let p = AsteroidProvider {
        inner: &inner,
        ast_file,
        ast_id,
        body,
        source: EphemerisSource::Moshier,
    };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

// ---------------------------------------------------------------------------
// Planet-moon calc pipeline — ports calc_center_body (sweph.c:2445)
// ---------------------------------------------------------------------------

#[cfg(feature = "swisseph-files")]
pub(crate) struct PlanetMoonProvider<'a, P: PositionProvider> {
    inner: &'a P,
    moon_file: &'a SwissEphFile,
    moon_id: i32,
    parent: Body,
}

#[cfg(feature = "swisseph-files")]
impl<P: PositionProvider> PositionProvider for PlanetMoonProvider<'_, P> {
    fn positions(&self, _body: Body, jd: f64, need_speed: bool) -> Result<SwephPositions, Error> {
        let mut pos = self.inner.positions(self.parent, jd, need_speed)?;
        let n = if need_speed { 6 } else { 3 };
        let (offset, _) = evaluate_body(self.moon_file, self.moon_id, jd, need_speed)?;
        for i in 0..n {
            pos.planet_bary[i] += offset[i];
        }
        Ok(pos)
    }

    fn moon_geo(&self, jd: f64, need_speed: bool) -> Result<[f64; 6], Error> {
        self.inner.moon_geo(jd, need_speed)
    }

    fn updates_sun_in_light_time(&self) -> bool {
        self.inner.updates_sun_in_light_time()
    }
}

#[cfg(feature = "swisseph-files")]
pub(crate) struct MoshierPlanetProvider<'a> {
    pub(crate) eps_j2000: &'a Epsilon,
}

#[cfg(feature = "swisseph-files")]
impl PositionProvider for MoshierPlanetProvider<'_> {
    fn positions(&self, body: Body, jd: f64, _need_speed: bool) -> Result<SwephPositions, Error> {
        let pp = compute_pipeline(jd, body, self.eps_j2000)?;
        let earth_pp = compute_pipeline(jd, Body::Sun, self.eps_j2000)?;
        Ok(SwephPositions {
            planet_bary: pp.planet_helio,
            earth_bary: earth_pp.earth_helio,
            earth_helio: earth_pp.earth_helio,
            sun_bary: [0.0; 6],
        })
    }

    fn moon_geo(&self, jd: f64, _need_speed: bool) -> Result<[f64; 6], Error> {
        raw_osc_moon_moshier(jd, self.eps_j2000)
    }
}

#[cfg(feature = "swisseph-files")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_plmoon_sweph(
    jd: f64,
    body: Body,
    moon_file: &SwissEphFile,
    moon_id: i32,
    parent: Body,
    planet_files: &[SwissEphFile],
    moon_files: &[SwissEphFile],
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let inner = SwephProvider {
        planet_files,
        moon_files,
    };
    let p = PlanetMoonProvider {
        inner: &inner,
        moon_file,
        moon_id,
        parent,
    };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

#[cfg(all(feature = "swisseph-files", feature = "jpl"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_plmoon_jpl(
    jd: f64,
    body: Body,
    moon_file: &SwissEphFile,
    moon_id: i32,
    parent: Body,
    jpl_file: &crate::jpl::JplFile,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let inner = JplProvider { file: jpl_file };
    let p = PlanetMoonProvider {
        inner: &inner,
        moon_file,
        moon_id,
        parent,
    };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

#[cfg(feature = "swisseph-files")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_plmoon_moshier(
    jd: f64,
    body: Body,
    moon_file: &SwissEphFile,
    moon_id: i32,
    parent: Body,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let inner = MoshierPlanetProvider { eps_j2000 };
    let p = PlanetMoonProvider {
        inner: &inner,
        moon_file,
        moon_id,
        parent,
    };
    apparent_planet(&p, jd, body, eps_j2000, flags, config, models)
}

// ---------------------------------------------------------------------------
// Fictitious planet calc pipeline — ports app_pos_etc_plan_osc (sweph.c:3365)
// ---------------------------------------------------------------------------
//
// Structurally distinct from apparent_planet (which ports app_pos_etc_plan):
// - Light-time loop works on barycentric/heliocentric positions; observer
//   subtraction happens AFTER the loop (sweph.c:3497).
// - HELCTR/BARYCTR handled via observer-zeroing, not early returns.
// - niter=1 always (2 passes), regardless of backend.
// - Speed refinement re-evaluates osc_el_plan at t-dt.

#[allow(clippy::too_many_arguments)]
fn apparent_fictitious<P: PositionProvider>(
    p: &P,
    jd: f64,
    catalog: &crate::fictitious::FictitiousCatalog,
    ipl: usize,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let need_speed = flags.contains(CalcFlags::SPEED);

    let pos = p.positions(Body::Sun, jd, true)?;
    let pdp_x =
        crate::fictitious::osc_el_plan(jd, catalog, ipl, &pos.earth_bary, &pos.sun_bary, models)?;
    let mut xx = pdp_x;

    // Observer: geocenter, topocenter, heliocenter, or barycenter (sweph.c:3396-3422)
    let offset = topo_offset(jd, flags, config, models);
    let xobs = if flags.contains(CalcFlags::BARYCTR) {
        [0.0; 6]
    } else if flags.contains(CalcFlags::HELCTR) {
        pos.sun_bary
    } else {
        let mut o = pos.earth_bary;
        for i in 0..6 {
            o[i] += offset[i];
        }
        o
    };

    let is_geo = !flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR);

    // Light-time (sweph.c:3426-3493)
    let mut dt = 0.0;
    let mut xxsp = [0.0; 3];
    let mut xobs2 = [0.0; 6];
    if !flags.contains(CalcFlags::TRUEPOS) {
        // Speed pre-pass: estimate dt at t-1 (sweph.c:3428-3453)
        if need_speed {
            let xxsv_sp = [
                pdp_x[0] - pdp_x[3],
                pdp_x[1] - pdp_x[4],
                pdp_x[2] - pdp_x[5],
            ];
            let mut xxsp_tmp = xxsv_sp;
            for _ in 0..=1 {
                let mut dx = xxsp_tmp;
                if is_geo {
                    for i in 0..3 {
                        dx[i] -= xobs[i] - xobs[i + 3];
                    }
                }
                let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                let dt_sp = dist * AUNIT / CLIGHT / 86400.0;
                for i in 0..3 {
                    xxsp_tmp[i] = xxsv_sp[i] - dt_sp * pdp_x[i + 3];
                }
            }
            for i in 0..3 {
                xxsp[i] = xxsv_sp[i] - xxsp_tmp[i];
            }
        }

        // Main light-time loop (sweph.c:3456-3471)
        for _ in 0..=1 {
            let mut dx = [xx[0], xx[1], xx[2]];
            if is_geo {
                for i in 0..3 {
                    dx[i] -= xobs[i];
                }
            }
            dt = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt() * AUNIT / CLIGHT / 86400.0;
            for i in 0..3 {
                xx[i] = pdp_x[i] - dt * pdp_x[i + 3];
                xx[i + 3] = pdp_x[i + 3];
            }
        }

        // Speed refinement: re-evaluate at t-dt (sweph.c:3472-3492)
        if need_speed {
            for i in 0..3 {
                xxsp[i] = pdp_x[i] - xx[i] - xxsp[i];
            }
            let t = jd - dt;
            let pos_ret = p.positions(Body::Sun, t, true)?;
            xx = crate::fictitious::osc_el_plan(
                t,
                catalog,
                ipl,
                &pos_ret.earth_bary,
                &pos_ret.sun_bary,
                models,
            )?;
            if flags.contains(CalcFlags::TOPOCTR) {
                let offset_ret = topo_offset(t, flags, config, models);
                xobs2 = pos_ret.earth_bary;
                for i in 0..6 {
                    xobs2[i] += offset_ret[i];
                }
            } else {
                xobs2 = pos_ret.earth_bary;
            }
        }
    }

    // Geocentric conversion — uses original-epoch xobs (sweph.c:3497-3498)
    for i in 0..6 {
        xx[i] -= xobs[i];
    }
    // Speed dt-change correction (sweph.c:3505-3507)
    if !flags.contains(CalcFlags::TRUEPOS) && need_speed {
        for i in 0..3 {
            xx[i + 3] -= xxsp[i];
        }
    }
    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // Deflection (sweph.c:3515-3517) — C's swi_deflect_light adds swed.topd.xobs
    // to the Earth position when TOPOCTR (sweph.c:3758-3760).
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOGDEFL) {
        let mut xobs_helio = [0.0; 6];
        for i in 0..6 {
            xobs_helio[i] = pos.earth_bary[i] - pos.sun_bary[i] + offset[i];
        }
        let mut planet_helio = [0.0; 6];
        for i in 0..3 {
            planet_helio[i] = xx[i] + xobs_helio[i];
            planet_helio[i + 3] = pdp_x[i + 3];
        }
        deflect_light(&mut xx, &xobs_helio, &planet_helio, need_speed);
    }

    // Aberration (sweph.c:3521-3531) — xobs is original epoch, xobs2 is t-dt
    if !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOABERR) {
        aberr_light(&mut xx, &[xobs[3], xobs[4], xobs[5]], need_speed);
        if need_speed {
            for i in 0..3 {
                xx[i + 3] += xobs[i + 3] - xobs2[i + 3];
            }
        }
    }

    if !need_speed {
        xx[3] = 0.0;
        xx[4] = 0.0;
        xx[5] = 0.0;
    }

    // C's app_pos_etc_plan_osc does NOT call swi_bias — frame bias is omitted
    // for fictitious bodies (unlike the general app_pos_etc_plan at sweph.c:2758).

    let x2000 = xx;
    let (eps, nut_val, nutv) = precess_and_ephem(&mut xx, jd, flags, models);
    Ok((
        app_pos_rest(&mut xx, flags, &eps, &nut_val, nutv.as_ref()),
        x2000,
    ))
}

#[cfg(feature = "swisseph-files")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_fictitious_sweph(
    jd: f64,
    _body: Body,
    catalog: &crate::fictitious::FictitiousCatalog,
    ipl: usize,
    planet_files: &[SwissEphFile],
    moon_files: &[SwissEphFile],
    _eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = SwephProvider {
        planet_files,
        moon_files,
    };
    apparent_fictitious(&p, jd, catalog, ipl, flags, config, models)
}

#[cfg(feature = "jpl")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_fictitious_jpl(
    jd: f64,
    _body: Body,
    catalog: &crate::fictitious::FictitiousCatalog,
    ipl: usize,
    jpl_file: &crate::jpl::JplFile,
    _eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = JplProvider { file: jpl_file };
    apparent_fictitious(&p, jd, catalog, ipl, flags, config, models)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_fictitious_moshier(
    jd: f64,
    _body: Body,
    catalog: &crate::fictitious::FictitiousCatalog,
    ipl: usize,
    eps_j2000: &Epsilon,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> Result<([f64; 24], [f64; 6]), Error> {
    let p = MoshierEarthProvider { eps_j2000 };
    apparent_fictitious(&p, jd, catalog, ipl, flags, config, models)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::CalcFlags;
    use crate::types::EphemerisSource;

    #[test]
    fn requested_source_precedence() {
        assert_eq!(
            requested_source(CalcFlags::MOSEPH),
            Some(EphemerisSource::Moshier)
        );
        assert_eq!(
            requested_source(CalcFlags::JPLEPH),
            Some(EphemerisSource::Jpl)
        );
        assert_eq!(
            requested_source(CalcFlags::SWIEPH),
            Some(EphemerisSource::Swiss)
        );
        assert_eq!(requested_source(CalcFlags::empty()), None);
        // MOSEPH wins over JPLEPH (C precedence sweph.c:375-381)
        assert_eq!(
            requested_source(CalcFlags::MOSEPH | CalcFlags::JPLEPH),
            Some(EphemerisSource::Moshier)
        );
        // MOSEPH wins over SWIEPH
        assert_eq!(
            requested_source(CalcFlags::MOSEPH | CalcFlags::SWIEPH),
            Some(EphemerisSource::Moshier)
        );
        // JPLEPH wins over SWIEPH
        assert_eq!(
            requested_source(CalcFlags::JPLEPH | CalcFlags::SWIEPH),
            Some(EphemerisSource::Jpl)
        );
        // All three set → MOSEPH wins
        assert_eq!(
            requested_source(CalcFlags::MOSEPH | CalcFlags::JPLEPH | CalcFlags::SWIEPH),
            Some(EphemerisSource::Moshier)
        );
    }
}
