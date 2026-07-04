//! Osculating (Keplerian) orbital elements and geocentric distance extrema —
//! `swe_get_orbital_elements` / `swe_orbit_max_min_true_distance`.
//!
//! * [`Ephemeris::get_orbital_elements`](crate::Ephemeris::get_orbital_elements) derives the momentary two-body Kepler elements of a
//!   planet, asteroid, or the Earth-Moon barycentre from its J2000 state vector
//!   (a full port of swecl.c:5783-5971, all 17 `dret` slots).
//! * [`Ephemeris::orbit_max_min_true_distance`](crate::Ephemeris::orbit_max_min_true_distance) returns the maximum, minimum, and current
//!   true distance of a body: heliocentric bodies read the extrema straight off
//!   their own Kepler ellipse, geocentric ones run the two-ellipse
//!   coordinate-descent search of swecl.c:6170-6287.
//!
//! Reference: `docs/c-ref-orbital-elements.md`; C `swecl.c:5687-6287`.
//!
//! ## `SEFLG_ORBEL_AA` ≡ `SEFLG_TOPOCTR`
//! Orbital elements have no topocentric variant, so C repurposes the `TOPOCTR`
//! bit as `SEFLG_ORBEL_AA` (swephexp.h:207). We replicate the bit-aliasing:
//! passing [`CalcFlags::TOPOCTR`] to this API means "sum masses inside the
//! orbit" (Astronomical Almanac method), never a topocentric request — the bit
//! never reaches `eph.calc` (it is not in `EPHMASK`), it only gates
//! `get_gmsm`'s mass summation.
//!
//! ## `SEFLG_BARYCTR` limitation
//! The barycentric branch (`r > 6` AU) issues an `eph.calc(.. | BARYCTR)` which
//! neither this codebase's `calc` pipeline nor C's Moshier backend supports —
//! both return an error there (verified against C: `swe_get_orbital_elements`
//! with `SEFLG_MOSEPH | SEFLG_BARYCTR` for Pluto returns "barycentric Moshier
//! positions are not supported"). BARYCTR requests for bodies inside 6 AU fall
//! back to HELCTR and work normally; the golden battery covers only that
//! reachable path (`tests/golden/orbit.rs`).

use std::f64::consts::PI;

use crate::calc::EPHMASK;
use crate::constants::{
    AUNIT, DEGTORAD, EARTH_MOON_MRAT, GEOGCONST, HELGRAVCONST, IPL_TO_ELEM, J2000, PLMASS, RADTODEG,
};
use crate::context::Ephemeris;
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{
    cartesian_to_polar, cross_prod, mod_2pi, normalize_degrees, polar_to_cartesian, rotate_x_sincos,
};
use crate::types::Body;

/// `SEFLG_ORBEL_AA` — bit-aliased onto [`CalcFlags::TOPOCTR`] (see module docs).
const ORBEL_AA: CalcFlags = CalcFlags::TOPOCTR;

/// Reference distance (metres) the Moon's semimajor axis is rescaled against to
/// obtain its sidereal period in months (C literal, swecl.c:5928).
const MOON_SEMA_REF_M: f64 = 383397772.5;

/// Osculating (Keplerian) orbital elements — the 17 `dret` slots of
/// `swe_get_orbital_elements` (swecl.c:5949-5965), one named field each.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OrbitalElements {
    /// `dret[0]` — semimajor axis (AU).
    pub semi_major_axis: f64,
    /// `dret[1]` — eccentricity.
    pub eccentricity: f64,
    /// `dret[2]` — inclination (deg).
    pub inclination: f64,
    /// `dret[3]` — longitude of ascending node Ω (deg).
    pub ascending_node: f64,
    /// `dret[4]` — argument of perihelion ω (deg).
    pub arg_perihelion: f64,
    /// `dret[5]` — longitude of perihelion ϖ = Ω + ω (deg).
    pub perihelion_lon: f64,
    /// `dret[6]` — mean anomaly at epoch M₀ (deg).
    pub mean_anomaly: f64,
    /// `dret[7]` — true anomaly at epoch (deg).
    pub true_anomaly: f64,
    /// `dret[8]` — eccentric anomaly at epoch (deg).
    pub eccentric_anomaly: f64,
    /// `dret[9]` — mean longitude at epoch (deg).
    pub mean_longitude: f64,
    /// `dret[10]` — sidereal orbital period (tropical years, J2000).
    pub sidereal_period: f64,
    /// `dret[11]` — mean daily motion (deg/day).
    pub mean_daily_motion: f64,
    /// `dret[12]` — tropical period (years).
    pub tropical_period: f64,
    /// `dret[13]` — synodic period (days); negative for inner planets / Moon.
    pub synodic_period: f64,
    /// `dret[14]` — JD (TT) of perihelion passage.
    pub perihelion_passage: f64,
    /// `dret[15]` — perihelion distance (AU).
    pub perihelion_distance: f64,
    /// `dret[16]` — aphelion distance (AU).
    pub aphelion_distance: f64,
}

impl OrbitalElements {
    fn from_dret(d: &[f64; 17]) -> Self {
        Self {
            semi_major_axis: d[0],
            eccentricity: d[1],
            inclination: d[2],
            ascending_node: d[3],
            arg_perihelion: d[4],
            perihelion_lon: d[5],
            mean_anomaly: d[6],
            true_anomaly: d[7],
            eccentric_anomaly: d[8],
            mean_longitude: d[9],
            sidereal_period: d[10],
            mean_daily_motion: d[11],
            tropical_period: d[12],
            synodic_period: d[13],
            perihelion_passage: d[14],
            perihelion_distance: d[15],
            aphelion_distance: d[16],
        }
    }

    /// The 17 fields as a flat `dret`-ordered array (slots 0..16), for
    /// differential comparison against C's `dret[]`.
    pub fn as_array(&self) -> [f64; 17] {
        [
            self.semi_major_axis,
            self.eccentricity,
            self.inclination,
            self.ascending_node,
            self.arg_perihelion,
            self.perihelion_lon,
            self.mean_anomaly,
            self.true_anomaly,
            self.eccentric_anomaly,
            self.mean_longitude,
            self.sidereal_period,
            self.mean_daily_motion,
            self.tropical_period,
            self.synodic_period,
            self.perihelion_passage,
            self.perihelion_distance,
            self.aphelion_distance,
        ]
    }
}

// ---------------------------------------------------------------------------
// Small vector helpers (C `square_sum`/`dot_prod` macros, sweph.h:308-309)
// ---------------------------------------------------------------------------

#[inline]
fn square_sum(x: &[f64]) -> f64 {
    x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
}

#[inline]
fn dot_prod(x: &[f64], y: &[f64]) -> f64 {
    x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
}

// ---------------------------------------------------------------------------
// get_gmsm — GM of the orbit's central body (swecl.c:5687-5742)
// ---------------------------------------------------------------------------

/// Central-body GM (AU³/day²) for `ipl`'s orbit. `r` is the body's already-known
/// heliocentric distance (AU), used only by the asteroid/`ORBEL_AA` branch.
///
/// PORTS TWO C QUIRKS LITERALLY (see `docs/c-ref-orbital-elements.md`):
/// `IPL_TO_ELEM[SE_PLUTO] == 0` (Pluto reuses Mercury's mass row), which in the
/// `ORBEL_AA` summation loop for Pluto also double-counts Mercury. ~1.6e-7
/// relative error — replicated for golden fidelity, not fixed.
fn get_gmsm(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    iflag: CalcFlags,
    r: f64,
) -> Result<f64, Error> {
    // Flags for the asteroid/AA sub-branch's per-planet re-queries.
    let mut iflj2000p = (iflag & (EPHMASK | CalcFlags::HELCTR | CalcFlags::BARYCTR))
        | CalcFlags::J2000
        | CalcFlags::TRUEPOS
        | CalcFlags::NONUT;
    if !iflj2000p.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
        iflj2000p |= CalcFlags::HELCTR;
    }
    let aa = iflag.contains(ORBEL_AA);

    if ipl == Body::Moon {
        return Ok(
            GEOGCONST * (1.0 + 1.0 / EARTH_MOON_MRAT) / AUNIT / AUNIT / AUNIT * 86400.0 * 86400.0,
        );
    }

    let raw = ipl.to_raw_id();
    let mut plm = 0.0;
    if (2..=9).contains(&raw) || raw == 14 {
        // Mercury..Pluto or Earth.
        if aa {
            if raw == 14 {
                // Earth: explicit Earth + Venus + Mercury (swecl.c:5703-5705).
                plm = 1.0 / PLMASS[IPL_TO_ELEM[14]];
                plm += 1.0 / PLMASS[IPL_TO_ELEM[3]]; // Venus
                plm += 1.0 / PLMASS[IPL_TO_ELEM[2]]; // Mercury
            } else {
                // Sum masses inside the orbit (swecl.c:5707-5711).
                let mut j = raw;
                while j >= 2 {
                    plm += 1.0 / PLMASS[IPL_TO_ELEM[j as usize]];
                    j -= 1;
                }
                if raw >= 4 {
                    // ipl >= SE_MARS: fold in Earth (the descending loop stops
                    // before Earth's id 14).
                    plm += 1.0 / PLMASS[IPL_TO_ELEM[14]];
                }
            }
        } else {
            // Two-body: single term (subject to the Pluto quirk above).
            plm = 1.0 / PLMASS[IPL_TO_ELEM[raw as usize]];
        }
        Ok(HELGRAVCONST * (1.0 + plm) / AUNIT / AUNIT / AUNIT * 86400.0 * 86400.0)
    } else {
        // Asteroid / fictitious body.
        if aa {
            for j in 2..=9 {
                let body = Body::try_from(j)?;
                let x = eph.calc(tjd_et, body, iflj2000p)?;
                if r > x.data[2] {
                    plm += 1.0 / PLMASS[IPL_TO_ELEM[j as usize]];
                }
            }
            // calc(Body::Earth) returns [0;6] in this stateless port (Earth
            // is the observer origin).  Geocentric Sun distance == Earth's
            // heliocentric distance, so query Sun without HELCTR/BARYCTR.
            let earth_r = eph
                .calc(
                    tjd_et,
                    Body::Sun,
                    iflj2000p & !(CalcFlags::HELCTR | CalcFlags::BARYCTR),
                )?
                .data[2];
            if r > earth_r {
                plm += 1.0 / PLMASS[IPL_TO_ELEM[14]];
            }
        }
        Ok(HELGRAVCONST * (1.0 + plm) / AUNIT / AUNIT / AUNIT * 86400.0 * 86400.0)
    }
}

// ---------------------------------------------------------------------------
// swe_get_orbital_elements (swecl.c:5783-5971)
// ---------------------------------------------------------------------------

/// Osculating Kepler elements of `ipl` at `tjd_et` (TT). Port of
/// `swe_get_orbital_elements` (swecl.c:5783-5971). Rejects the Sun, the lunar
/// nodes, and the apsides ([`Error::InvalidBody`]).
pub(crate) fn get_orbital_elements(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    iflag: CalcFlags,
) -> Result<OrbitalElements, Error> {
    let dret = get_orbital_elements_dret(eph, tjd_et, ipl, iflag)?;
    Ok(OrbitalElements::from_dret(&dret))
}

/// Core derivation, returning the raw 17-slot `dret` array (reused internally by
/// the distance-search routines).
fn get_orbital_elements_dret(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    iflag: CalcFlags,
) -> Result<[f64; 17], Error> {
    // 2.1 — reject Sun (id 0 / <=0), nodes, apsides.
    let raw = ipl.to_raw_id();
    if raw <= 0
        || matches!(
            ipl,
            Body::MeanNode
                | Body::TrueNode
                | Body::MeanApogee
                | Body::OscuApogee
                | Body::IntpApogee
                | Body::IntpPerigee
        )
    {
        return Err(Error::InvalidBody(raw));
    }

    // Cartesian J2000, geometric, no nutation, with speed; center chosen below.
    let mut iflj2000 = (iflag & EPHMASK)
        | CalcFlags::J2000
        | CalcFlags::XYZ
        | CalcFlags::TRUEPOS
        | CalcFlags::NONUT
        | CalcFlags::SPEED;
    // Distance probe (polar; no center flag => default center, no XYZ).
    let iflj2000p = (iflag & EPHMASK)
        | CalcFlags::J2000
        | CalcFlags::TRUEPOS
        | CalcFlags::NONUT
        | CalcFlags::SPEED;

    // 2.2 — heliocentric distance probe + center-flag decision.
    let r = eph.calc(tjd_et, ipl, iflj2000p)?.data[2];
    if ipl != Body::Moon {
        if iflag.contains(CalcFlags::BARYCTR) && r > 6.0 {
            iflj2000 |= CalcFlags::BARYCTR; // only planets beyond ~Jupiter
        } else {
            iflj2000 |= CalcFlags::HELCTR;
        }
    }

    // 2.3 — GM and final position query.
    let gmsm = get_gmsm(eph, tjd_et, ipl, iflag, r)?;
    let mut xpos = eph.calc(tjd_et, ipl, iflj2000)?.data;
    if ipl == Body::Earth {
        // "Earth" elements are actually the Earth-Moon barycentre.
        let xposm = eph
            .calc(
                tjd_et,
                Body::Moon,
                iflj2000 & !(CalcFlags::BARYCTR | CalcFlags::HELCTR),
            )?
            .data;
        for j in 0..6 {
            xpos[j] += xposm[j] / (EARTH_MOON_MRAT + 1.0);
        }
    }

    // 2.4 — first-pass node vector via r_z / v_z projection.
    let fac = xpos[2] / xpos[5];
    let sgn = xpos[5] / xpos[5].abs();
    let mut xn = [0.0f64; 3];
    let mut xs = [0.0f64; 3];
    for j in 0..3 {
        xn[j] = (xpos[j] - fac * xpos[j + 3]) * sgn;
        xs[j] = -xn[j];
    }
    let mut rxy = (xn[0] * xn[0] + xn[1] * xn[1]).sqrt();
    let cosnode = xn[0] / rxy;
    let sinnode = xn[1] / rxy;

    // 2.5 — inclination via r × v.
    let xnorm = cross_prod([xpos[0], xpos[1], xpos[2]], [xpos[3], xpos[4], xpos[5]]);
    rxy = xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
    let c2 = rxy + xnorm[2] * xnorm[2];
    let mut rxyz = c2.sqrt();
    rxy = rxy.sqrt();
    let sinincl = rxy / rxyz;
    let mut cosincl = (1.0 - sinincl * sinincl).sqrt();
    if xnorm[2] < 0.0 {
        cosincl = -cosincl; // retrograde orbit (e.g. 20461 Dioretsa)
    }
    let incl = cosincl.acos() * RADTODEG;

    // 2.6 — argument of latitude, semimajor axis, eccentricity.
    let cosu = xpos[0] * cosnode + xpos[1] * sinnode;
    let sinu = xpos[2] / sinincl;
    let uu = sinu.atan2(cosu);
    rxyz = square_sum(&xpos[0..3]).sqrt();
    let v2 = square_sum(&xpos[3..6]);
    let sema = 1.0 / (2.0 / rxyz - v2 / gmsm);
    let pp = c2 / gmsm;
    let mut ecce = pp / sema;
    if ecce > 1.0 {
        ecce = 1.0;
    }
    ecce = (1.0 - ecce).sqrt();

    // 2.7 — eccentric and true anomaly.
    let mut ecce2 = ecce;
    if ecce2 == 0.0 {
        ecce2 = 0.0000000001;
    }
    let cos_ea = 1.0 / ecce2 * (1.0 - rxyz / sema);
    let sin_ea = 1.0 / ecce2 / (sema * gmsm).sqrt() * dot_prod(&xpos[0..3], &xpos[3..6]);
    let eanom = normalize_degrees(sin_ea.atan2(cos_ea) * RADTODEG);
    let mut ny = 2.0 * (((1.0 + ecce) / (1.0 - ecce)).sqrt() * sin_ea / (1.0 + cos_ea)).atan();
    let mut tanom = normalize_degrees(ny * RADTODEG);
    if eanom > 180.0 && tanom < 180.0 {
        tanom += 180.0;
    }
    if eanom < 180.0 && tanom > 180.0 {
        tanom -= 180.0;
    }

    // 2.8 — mean anomaly.
    let manom = normalize_degrees(eanom - ecce * RADTODEG * (eanom * DEGTORAD).sin());

    // 2.9 — perihelion / aphelion direction; node/apsis refinement.
    let mut xq = [0.0f64; 3];
    xq[0] = mod_2pi(uu - ny);
    let parg = xq[0] * RADTODEG;
    xq[1] = 0.0;
    xq[2] = sema * (1.0 - ecce);
    xq = polar_to_cartesian(xq);
    xq = rotate_x_sincos(xq, -sinincl, cosincl);
    xq = cartesian_to_polar(xq);
    xq[0] += sinnode.atan2(cosnode);
    let mut xa = [0.0f64; 3];
    xa[0] = mod_2pi(xq[0] + PI);
    xa[1] = -xq[1];
    xa[2] = sema * (1.0 + ecce);
    xq = polar_to_cartesian(xq);
    xa = polar_to_cartesian(xa);
    let _ = xa; // aphelion vector computed for parity; not exposed via dret.

    ny = mod_2pi(ny - uu);
    let ny2 = mod_2pi(ny + PI);
    let cos_en = (2.0 * ((ny / 2.0).tan() / ((1.0 + ecce) / (1.0 - ecce)).sqrt()).atan()).cos();
    let cos_en2 = (2.0 * ((ny2 / 2.0).tan() / ((1.0 + ecce) / (1.0 - ecce)).sqrt()).atan()).cos();
    let rn = sema * (1.0 - ecce * cos_en);
    let rn2 = sema * (1.0 - ecce * cos_en2);
    let ro = square_sum(&xn).sqrt();
    let ro2 = square_sum(&xs).sqrt();
    for j in 0..3 {
        xn[j] *= rn / ro;
        xs[j] *= rn2 / ro2;
    }
    xn = cartesian_to_polar(xn);
    xq = cartesian_to_polar(xq);
    let _ = xq;

    // 2.10 — final angle assembly.
    let node = xn[0] * RADTODEG;
    let peri = normalize_degrees(node + parg);
    let mlon = normalize_degrees(manom + peri);

    // 2.11 — period and daily motion.
    let mut csid = sema * sema.sqrt();
    if ipl == Body::Moon {
        let semam = sema * AUNIT / MOON_SEMA_REF_M;
        csid = semam * semam.sqrt();
        csid *= 27.32166 / 365.25636300;
    }
    let dmot = 0.9856076686 / csid;
    csid *= 365.25636 / 365.242189;
    let t = (tjd_et - J2000) / 365250.0;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    let pa = (50288.200 + 222.4045 * t + 0.2095 * t2 - 0.9408 * t3 - 0.0090 * t4 + 0.0010 * t5)
        / 3600.0
        / 365250.0;
    let mut ysid =
        (1295977422.83429 - 2.0 * 2.0441 * t - 3.0 * 0.00523 * t * t) / 3600.0 / 365250.0;
    ysid = 360.0 / ysid;
    let mut ytrop = (1296027711.03429 + 2.0 * 109.15809 * t + 3.0 * 0.07207 * t2
        - 4.0 * 0.23530 * t3
        - 5.0 * 0.00180 * t4
        + 6.0 * 0.00020 * t5)
        / 3600.0
        / 365250.0;
    ytrop = 360.0 / ytrop;
    let mut ctro = 360.0 / (dmot + pa) / 365.242189;
    ctro *= ysid / ytrop;
    let csyn = if ipl == Body::Earth {
        0.0
    } else {
        360.0 / (0.9856076686 - dmot)
    };

    // 2.12 — dret assembly.
    let mut dret = [0.0f64; 17];
    dret[0] = sema;
    dret[1] = ecce;
    dret[2] = incl;
    dret[3] = node;
    dret[4] = parg;
    dret[5] = peri;
    dret[6] = manom;
    dret[7] = tanom;
    dret[8] = eanom;
    dret[9] = mlon;
    dret[10] = csid;
    dret[11] = dmot;
    dret[12] = ctro;
    dret[13] = csyn;
    dret[14] = tjd_et - dret[6] / dmot;
    dret[15] = sema * (1.0 - ecce);
    dret[16] = sema * (1.0 + ecce);
    Ok(dret)
}

// ---------------------------------------------------------------------------
// Ellipse sampling helpers (swecl.c:5973-6096)
// ---------------------------------------------------------------------------

/// Precompute the 12-element Gauss P/Q rotation + shape block from Kepler
/// elements `dp[0..5]` = a, e, i, Ω, ω. Port of `osc_get_orbit_constants`.
fn osc_get_orbit_constants(dp: &[f64]) -> [f64; 12] {
    let sema = dp[0];
    let ecce = dp[1];
    let incl = dp[2];
    let node = dp[3];
    let parg = dp[4];
    let cosnode = (node * DEGTORAD).cos();
    let sinnode = (node * DEGTORAD).sin();
    let cosincl = (incl * DEGTORAD).cos();
    let sinincl = (incl * DEGTORAD).sin();
    let cosparg = (parg * DEGTORAD).cos();
    let sinparg = (parg * DEGTORAD).sin();
    let fac = ((1.0 - ecce) * (1.0 + ecce)).sqrt();
    [
        cosparg * cosnode - sinparg * cosincl * sinnode,
        -sinparg * cosnode - cosparg * cosincl * sinnode,
        sinincl * sinnode,
        cosparg * sinnode + sinparg * cosincl * cosnode,
        -sinparg * sinnode + cosparg * cosincl * cosnode,
        -sinincl * cosnode,
        sinparg * sinincl,
        cosparg * sinincl,
        cosincl,
        sema,
        ecce,
        fac,
    ]
}

/// Ecliptic-cartesian position on the ellipse at eccentric anomaly `ean` (deg).
/// Port of `osc_get_ecl_pos`.
fn osc_get_ecl_pos(ean: f64, pqr: &[f64; 12]) -> [f64; 3] {
    let cose = (ean * DEGTORAD).cos();
    let sine = (ean * DEGTORAD).sin();
    let sema = pqr[9];
    let ecce = pqr[10];
    let fac = pqr[11];
    let x0 = sema * (cose - ecce);
    let x1 = sema * fac * sine;
    [
        pqr[0] * x0 + pqr[1] * x1,
        pqr[3] * x0 + pqr[4] * x1,
        pqr[6] * x0 + pqr[7] * x1,
    ]
}

/// Euclidean distance between two 3-vectors. Port of `get_dist_from_2_vectors`.
fn get_dist_from_2_vectors(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
    let r0 = x1[0] - x2[0];
    let r1 = x1[1] - x2[1];
    let r2 = x1[2] - x2[2];
    (r0 * r0 + r1 * r1 + r2 * r2).sqrt()
}

/// Coordinate-descent hill-climb for the local MAX distance from the fixed body
/// `xb`, over one ellipse's eccentric anomaly. Port of `osc_iterate_max_dist`
/// (swecl.c:6026-6060).
///
/// NOTE (verified against C): the search always **restarts from `ean = 0`** —
/// the caller's rough-scan anomaly is discarded (swecl.c:6032). On return `xa`
/// holds the last-evaluated (overshoot) position, NOT the position at the
/// returned optimum; the caller feeds that overshoot buffer as the fixed `xb`
/// of the alternating call, so this side effect must be preserved exactly.
fn osc_iterate_max_dist(pqr: &[f64; 12], xa: &mut [f64; 3], xb: &[f64; 3], high_prec: bool) -> f64 {
    let dstep_min = if high_prec { 0.000001 } else { 1.0 };
    let mut ean = 0.0;
    let mut eansv = 0.0;
    *xa = osc_get_ecl_pos(ean, pqr);
    let mut r = get_dist_from_2_vectors(xb, xa);
    let mut rmax = r;
    let mut dstep = 1.0;
    while dstep >= dstep_min {
        for i in 0..2 {
            while r >= rmax {
                eansv = ean;
                if i == 0 {
                    ean += dstep;
                } else {
                    ean -= dstep;
                }
                *xa = osc_get_ecl_pos(ean, pqr);
                r = get_dist_from_2_vectors(xb, xa);
                if r > rmax {
                    rmax = r;
                }
            }
            ean = eansv;
            r = rmax;
        }
        ean = eansv;
        r = rmax;
        dstep /= 10.0;
    }
    rmax
}

/// Coordinate-descent hill-climb for the local MIN distance. Port of
/// `osc_iterate_min_dist` (swecl.c:6062-6096); see [`osc_iterate_max_dist`] for
/// the shared `ean = 0` restart / `xa` overshoot semantics.
fn osc_iterate_min_dist(pqr: &[f64; 12], xa: &mut [f64; 3], xb: &[f64; 3], high_prec: bool) -> f64 {
    let dstep_min = if high_prec { 0.000001 } else { 1.0 };
    let mut ean = 0.0;
    let mut eansv = 0.0;
    *xa = osc_get_ecl_pos(ean, pqr);
    let mut r = get_dist_from_2_vectors(xb, xa);
    let mut rmin = r;
    let mut dstep = 1.0;
    while dstep >= dstep_min {
        for i in 0..2 {
            while r <= rmin {
                eansv = ean;
                if i == 0 {
                    ean += dstep;
                } else {
                    ean -= dstep;
                }
                *xa = osc_get_ecl_pos(ean, pqr);
                r = get_dist_from_2_vectors(xb, xa);
                if r < rmin {
                    rmin = r;
                }
            }
            ean = eansv;
            r = rmin;
        }
        ean = eansv;
        r = rmin;
        dstep /= 10.0;
    }
    rmin
}

// ---------------------------------------------------------------------------
// swe_orbit_max_min_true_distance (swecl.c:6101-6287)
// ---------------------------------------------------------------------------

/// Heliocentric-only branch (Sun, Moon, or an explicit HELCTR/BARYCTR request):
/// max/min come straight off the Kepler ellipse, true distance from the ellipse
/// evaluated at the body's own eccentric anomaly. Port of
/// `orbit_max_min_true_distance_helio` (swecl.c:6101-6128).
fn orbit_max_min_true_distance_helio(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    iflag: CalcFlags,
) -> Result<(f64, f64, f64), Error> {
    let iflagi = iflag & (EPHMASK | CalcFlags::HELCTR | CalcFlags::BARYCTR);
    let ipli = if ipl == Body::Sun { Body::Earth } else { ipl };
    let de = get_orbital_elements_dret(eph, tjd_et, ipli, iflagi)?;
    let dmax = de[16];
    let dmin = de[15];
    let pqri = osc_get_orbit_constants(&de);
    let xinner = osc_get_ecl_pos(de[8], &pqri);
    let dtrue = (xinner[0] * xinner[0] + xinner[1] * xinner[1] + xinner[2] * xinner[2]).sqrt();
    Ok((dmax, dmin, dtrue))
}

/// Maximum, minimum, and current true distance of `ipl` (AU). Port of
/// `swe_orbit_max_min_true_distance` (swecl.c:6170-6287). Returns
/// `(dmax, dmin, dtrue)`.
pub(crate) fn orbit_max_min_true_distance(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    iflag: CalcFlags,
) -> Result<(f64, f64, f64), Error> {
    let iflagi = iflag & (EPHMASK | CalcFlags::HELCTR | CalcFlags::BARYCTR);
    if ipl == Body::Sun
        || ipl == Body::Moon
        || iflagi.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR)
    {
        return orbit_max_min_true_distance_helio(eph, tjd_et, ipl, iflagi);
    }

    // Geocentric two-ellipse search: the target's ellipse vs the EMB's.
    let dp = get_orbital_elements_dret(eph, tjd_et, ipl, iflagi)?;
    let de = get_orbital_elements_dret(eph, tjd_et, Body::Earth, iflagi)?;
    let (douter, dinner) = if de[0] > dp[0] { (de, dp) } else { (dp, de) };
    let pqro = osc_get_orbit_constants(&douter);
    let pqri = osc_get_orbit_constants(&dinner);
    let mut xouter = osc_get_ecl_pos(douter[8], &pqro);
    let mut xinner = osc_get_ecl_pos(dinner[8], &pqri);
    let rtrue = get_dist_from_2_vectors(&xouter, &xinner);

    // Rough grid scan. LOOP-BOUND QUIRK ported literally (swecl.c:6226-6231):
    // outer scans j*2° over 0..362°, inner scans i*1° over only 0..181° — half
    // the inner ellipse is never sampled. This is almost certainly a C bug, but
    // it selects which local extremum the refinement converges to, so it must be
    // reproduced (docs/c-ref-orbital-elements.md §8.2).
    let ncnt = 182;
    let dstep = 2.0;
    let mut rmax = 0.0;
    let mut rmin = 100000000.0;
    let mut max_xouter = [0.0f64; 3];
    let mut max_xinner = [0.0f64; 3];
    let mut min_xouter = [0.0f64; 3];
    let mut min_xinner = [0.0f64; 3];
    for j in 0..ncnt {
        let eano = j as f64 * dstep;
        xouter = osc_get_ecl_pos(eano, &pqro);
        for i in 0..ncnt {
            let eani = i as f64;
            xinner = osc_get_ecl_pos(eani, &pqri);
            let r = get_dist_from_2_vectors(&xouter, &xinner);
            if r > rmax {
                rmax = r;
                max_xouter = xouter;
                max_xinner = xinner;
            }
            if r < rmin {
                rmin = r;
                min_xouter = xouter;
                min_xinner = xinner;
            }
        }
    }

    // Refine maximum: block coordinate ascent, alternating ellipses. `rmax`
    // after each pass is the outer-ellipse refine's result (C's shared `rmax`,
    // reset inside each call, ends holding the second call's value).
    xouter = max_xouter;
    xinner = max_xinner;
    let mut rmaxsv = 0.0;
    for k in 0..=300 {
        osc_iterate_max_dist(&pqri, &mut xinner, &xouter, true);
        rmax = osc_iterate_max_dist(&pqro, &mut xouter, &xinner, true);
        if k > 0 && (rmax - rmaxsv).abs() < 0.00000001 {
            break;
        }
        rmaxsv = rmax;
    }

    // Refine minimum.
    xouter = min_xouter;
    xinner = min_xinner;
    let mut rminsv = 0.0;
    for k in 0..=300 {
        osc_iterate_min_dist(&pqri, &mut xinner, &xouter, true);
        rmin = osc_iterate_min_dist(&pqro, &mut xouter, &xinner, true);
        if k > 0 && (rmin - rminsv).abs() < 0.00000001 {
            break;
        }
        rminsv = rmin;
    }

    Ok((rmax, rmin, rtrue))
}
