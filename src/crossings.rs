// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

//! Longitude and latitude crossing search — find the instant a body reaches a
//! given ecliptic longitude or zero latitude (node crossing).

use crate::context::Ephemeris;
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{diff_degrees, normalize_degrees};
use crate::types::Body;

const CROSS_PRECISION: f64 = 1.0 / 3_600_000.0;

/// Newton normally converges in 3-5 iterations; hitting this cap means the iteration is
/// stuck in a cycle (see swisseph-rs/156) and we fall back to bisection.
const MAX_NEWTON_ITERATIONS: u32 = 50;
const MAX_BISECT_ITERATIONS: u32 = 100;

/// Newton refinement of a longitude crossing, shared by the sol/moon/helio crossers.
///
/// `eval(jd)` returns the body's (ecliptic longitude, longitude speed) at `jd`. On
/// convergence, returns the Julian Day after the final Newton step (matching C, which
/// applies the step before testing the tolerance).
///
/// With DE441 Moon files, Newton can enter a stable 2-cycle straddling the crossing and
/// never converge — the C library hangs on such inputs (swisseph-rs/156). After
/// [`MAX_NEWTON_ITERATIONS`] we bisect the sign-changing bracket accumulated during the
/// Newton phase. The cycle is caused by a longitude jump discontinuity at a Chebyshev
/// segment boundary in the compressed .se1 data, so bisection may terminate by bracket
/// collapse rather than by tolerance — see [`bisect_cross`].
fn refine_lon_cross<F>(x2cross: f64, mut jd: f64, mut eval: F) -> Result<f64, Error>
where
    F: FnMut(f64) -> Result<(f64, f64), Error>,
{
    let mut pos: Option<(f64, f64)> = None;
    let mut neg: Option<(f64, f64)> = None;
    for _ in 0..MAX_NEWTON_ITERATIONS {
        let (lon, speed) = eval(jd)?;
        let dist = diff_degrees(x2cross, lon);
        let jd_next = jd + dist / speed;
        if dist.abs() < CROSS_PRECISION {
            return Ok(jd_next);
        }
        if dist > 0.0 {
            pos = Some((jd, dist));
        } else {
            neg = Some((jd, dist));
        }
        jd = jd_next;
    }
    let (Some((jd_pos, dist_pos)), Some((jd_neg, _))) = (pos, neg) else {
        return Err(Error::NoConvergence);
    };
    bisect_cross(jd_pos, dist_pos, jd_neg, |jd| {
        Ok(diff_degrees(x2cross, eval(jd)?.0))
    })
}

/// Bisect between `jd_a` and `jd_b`, whose signed distances to the crossing have opposite
/// signs, until the distance falls below [`CROSS_PRECISION`].
fn bisect_cross<F>(
    mut jd_a: f64,
    mut dist_a: f64,
    mut jd_b: f64,
    mut eval_dist: F,
) -> Result<f64, Error>
where
    F: FnMut(f64) -> Result<f64, Error>,
{
    for _ in 0..MAX_BISECT_ITERATIONS {
        let jd_mid = 0.5 * (jd_a + jd_b);
        if jd_mid == jd_a || jd_mid == jd_b {
            // The bracket has collapsed to adjacent f64s without meeting the tolerance:
            // the compressed ephemeris longitude has a jump discontinuity at the crossing
            // (DE441 Chebyshev segment boundary) larger than CROSS_PRECISION. The collapsed
            // point locates the crossing time itself to full f64 precision.
            return Ok(jd_mid);
        }
        let dist_mid = eval_dist(jd_mid)?;
        if dist_mid.abs() < CROSS_PRECISION {
            return Ok(jd_mid);
        }
        if (dist_mid > 0.0) == (dist_a > 0.0) {
            jd_a = jd_mid;
            dist_a = dist_mid;
        } else {
            jd_b = jd_mid;
        }
    }
    Err(Error::NoConvergence)
}

/// Newton refinement of a node (zero-latitude) crossing, shared by the TT/UT variants.
///
/// `jd`/`x` are the day and position that bracketed the latitude sign change in the daily
/// scan; `eval(jd)` returns the Moon's position array at `jd`. Same cycle guard and
/// bisection fallback as [`refine_lon_cross`], keyed on latitude instead of longitude.
fn refine_node_cross<F>(mut jd: f64, mut x: [f64; 6], mut eval: F) -> Result<MoonCrossing, Error>
where
    F: FnMut(f64) -> Result<[f64; 6], Error>,
{
    let mut pos: Option<(f64, f64)> = None;
    let mut neg: Option<(f64, f64)> = None;
    for _ in 0..MAX_NEWTON_ITERATIONS {
        let dist = x[1];
        if dist > 0.0 {
            pos = Some((jd, dist));
        } else {
            neg = Some((jd, dist));
        }
        jd -= dist / x[4];
        x = eval(jd)?;
        if x[1].abs() < CROSS_PRECISION {
            return Ok(MoonCrossing {
                jd,
                longitude: x[0],
                latitude: x[1],
            });
        }
    }
    let (Some((jd_pos, dist_pos)), Some((jd_neg, _))) = (pos, neg) else {
        return Err(Error::NoConvergence);
    };
    let jd = bisect_cross(jd_pos, dist_pos, jd_neg, |jd| Ok(eval(jd)?[1]))?;
    let x = eval(jd)?;
    Ok(MoonCrossing {
        jd,
        longitude: x[0],
        latitude: x[1],
    })
}

/// Result of a Moon node-crossing search: the Julian Day and the Moon's ecliptic position at
/// that instant.
#[derive(Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MoonCrossing {
    /// Julian Day of the crossing (same time scale as the search, TT or UT).
    pub jd: f64,
    /// Moon's ecliptic longitude at the crossing, degrees.
    pub longitude: f64,
    /// Moon's ecliptic latitude at the crossing, degrees (near zero by construction).
    pub latitude: f64,
}

// ---------------------------------------------------------------------------
// solcross / solcross_ut
// ---------------------------------------------------------------------------

/// Next Julian Day (TT) at or after `jd_et` at which the Sun's ecliptic longitude equals
/// `x2cross` (degrees).
pub fn solcross(eph: &Ephemeris, x2cross: f64, jd_et: f64, flags: CalcFlags) -> Result<f64, Error> {
    let flags = flags | CalcFlags::SPEED;
    let body = Body::Sun;
    let r = eph.calc(jd_et, body, flags)?;
    let dist = normalize_degrees(x2cross - r.data[0]);
    let jd = jd_et + dist / (360.0 / 365.24);
    refine_lon_cross(x2cross, jd, |jd| {
        let r = eph.calc(jd, body, flags)?;
        Ok((r.data[0], r.data[3]))
    })
}

/// UT-based [`solcross`]: next Julian Day (UT) at which the Sun's ecliptic longitude equals
/// `x2cross` (degrees).
pub fn solcross_ut(
    eph: &Ephemeris,
    x2cross: f64,
    jd_ut: f64,
    flags: CalcFlags,
) -> Result<f64, Error> {
    let flags = flags | CalcFlags::SPEED;
    let body = Body::Sun;
    let r = eph.calc_ut(jd_ut, body, flags)?;
    let dist = normalize_degrees(x2cross - r.data[0]);
    let jd = jd_ut + dist / (360.0 / 365.24);
    refine_lon_cross(x2cross, jd, |jd| {
        let r = eph.calc_ut(jd, body, flags)?;
        Ok((r.data[0], r.data[3]))
    })
}

// ---------------------------------------------------------------------------
// mooncross / mooncross_ut
// ---------------------------------------------------------------------------

/// Next Julian Day (TT) at or after `jd_et` at which the Moon's ecliptic longitude equals
/// `x2cross` (degrees).
pub fn mooncross(
    eph: &Ephemeris,
    x2cross: f64,
    jd_et: f64,
    flags: CalcFlags,
) -> Result<f64, Error> {
    let flags = flags | CalcFlags::SPEED;
    let body = Body::Moon;
    let r = eph.calc(jd_et, body, flags)?;
    let dist = normalize_degrees(x2cross - r.data[0]);
    let jd = jd_et + dist / (360.0 / 27.32);
    refine_lon_cross(x2cross, jd, |jd| {
        let r = eph.calc(jd, body, flags)?;
        Ok((r.data[0], r.data[3]))
    })
}

/// UT-based [`mooncross`]: next Julian Day (UT) at which the Moon's ecliptic longitude equals
/// `x2cross` (degrees).
pub fn mooncross_ut(
    eph: &Ephemeris,
    x2cross: f64,
    jd_ut: f64,
    flags: CalcFlags,
) -> Result<f64, Error> {
    let flags = flags | CalcFlags::SPEED;
    let body = Body::Moon;
    let r = eph.calc_ut(jd_ut, body, flags)?;
    let dist = normalize_degrees(x2cross - r.data[0]);
    let jd = jd_ut + dist / (360.0 / 27.32);
    refine_lon_cross(x2cross, jd, |jd| {
        let r = eph.calc_ut(jd, body, flags)?;
        Ok((r.data[0], r.data[3]))
    })
}

// ---------------------------------------------------------------------------
// mooncross_node / mooncross_node_ut
// ---------------------------------------------------------------------------

/// Next Julian Day (TT) at or after `jd_et` at which the Moon crosses its (mean) orbital node
/// (ecliptic latitude passes through zero).
pub fn mooncross_node(
    eph: &Ephemeris,
    jd_et: f64,
    flags: CalcFlags,
) -> Result<MoonCrossing, Error> {
    let flags = flags | CalcFlags::SPEED;
    let body = Body::Moon;
    let r = eph.calc(jd_et, body, flags)?;
    let xlat = r.data[1];
    let mut jd = jd_et + 1.0;
    loop {
        let x = eph.calc(jd, body, flags)?.data;
        if (x[1] >= 0.0 && xlat < 0.0) || (x[1] < 0.0 && xlat > 0.0) {
            return refine_node_cross(jd, x, |jd| Ok(eph.calc(jd, body, flags)?.data));
        }
        jd += 1.0;
    }
}

/// UT-based [`mooncross_node`]: next Julian Day (UT) at which the Moon crosses its orbital node.
pub fn mooncross_node_ut(
    eph: &Ephemeris,
    jd_ut: f64,
    flags: CalcFlags,
) -> Result<MoonCrossing, Error> {
    let flags = flags | CalcFlags::SPEED;
    let body = Body::Moon;
    let r = eph.calc_ut(jd_ut, body, flags)?;
    let xlat = r.data[1];
    let mut jd = jd_ut + 1.0;
    loop {
        let x = eph.calc_ut(jd, body, flags)?.data;
        if (x[1] >= 0.0 && xlat < 0.0) || (x[1] < 0.0 && xlat > 0.0) {
            return refine_node_cross(jd, x, |jd| Ok(eph.calc_ut(jd, body, flags)?.data));
        }
        jd += 1.0;
    }
}

// ---------------------------------------------------------------------------
// helio_cross / helio_cross_ut
// ---------------------------------------------------------------------------

fn reject_helio_body(body: Body) -> bool {
    matches!(
        body,
        Body::Sun
            | Body::Moon
            | Body::MeanNode
            | Body::TrueNode
            | Body::MeanApogee
            | Body::OscuApogee
            | Body::IntpApogee
            | Body::IntpPerigee
    )
}

/// Next Julian Day (TT) at which `body`'s heliocentric ecliptic longitude equals `x2cross`
/// (degrees), starting from `jd_et`. `dir >= 0` searches forward in time, `dir < 0` backward.
pub fn helio_cross(
    eph: &Ephemeris,
    body: Body,
    x2cross: f64,
    jd_et: f64,
    flags: CalcFlags,
    dir: i32,
) -> Result<f64, Error> {
    if reject_helio_body(body) {
        return Err(Error::UnsupportedFlags(CalcFlags::HELCTR));
    }
    let flags = flags | CalcFlags::SPEED | CalcFlags::HELCTR;
    let r = eph.calc(jd_et, body, flags)?;
    let xlp = if body == Body::Chiron {
        0.01971
    } else {
        r.data[3]
    };
    let jd = if dir >= 0 {
        let dist = normalize_degrees(x2cross - r.data[0]);
        jd_et + dist / xlp
    } else {
        let dist = 360.0 - normalize_degrees(x2cross - r.data[0]);
        jd_et - dist / xlp
    };
    refine_lon_cross(x2cross, jd, |jd| {
        let r = eph.calc(jd, body, flags)?;
        Ok((r.data[0], r.data[3]))
    })
}

/// UT-based [`helio_cross`]: next Julian Day (UT) at which `body`'s heliocentric ecliptic
/// longitude equals `x2cross` (degrees).
pub fn helio_cross_ut(
    eph: &Ephemeris,
    body: Body,
    x2cross: f64,
    jd_ut: f64,
    flags: CalcFlags,
    dir: i32,
) -> Result<f64, Error> {
    if reject_helio_body(body) {
        return Err(Error::UnsupportedFlags(CalcFlags::HELCTR));
    }
    let flags = flags | CalcFlags::SPEED | CalcFlags::HELCTR;
    let r = eph.calc_ut(jd_ut, body, flags)?;
    let xlp = if body == Body::Chiron {
        0.01971
    } else {
        r.data[3]
    };
    let jd = if dir >= 0 {
        let dist = normalize_degrees(x2cross - r.data[0]);
        jd_ut + dist / xlp
    } else {
        let dist = 360.0 - normalize_degrees(x2cross - r.data[0]);
        jd_ut - dist / xlp
    };
    refine_lon_cross(x2cross, jd, |jd| {
        let r = eph.calc_ut(jd, body, flags)?;
        Ok((r.data[0], r.data[3]))
    })
}
