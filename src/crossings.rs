use crate::context::Ephemeris;
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{diff_degrees, normalize_degrees};
use crate::types::Body;

const CROSS_PRECISION: f64 = 1.0 / 3_600_000.0;

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
    let mut jd = jd_et + dist / (360.0 / 365.24);
    loop {
        let r = eph.calc(jd, body, flags)?;
        let dist = diff_degrees(x2cross, r.data[0]);
        jd += dist / r.data[3];
        if dist.abs() < CROSS_PRECISION {
            break;
        }
    }
    Ok(jd)
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
    let mut jd = jd_ut + dist / (360.0 / 365.24);
    loop {
        let r = eph.calc_ut(jd, body, flags)?;
        let dist = diff_degrees(x2cross, r.data[0]);
        jd += dist / r.data[3];
        if dist.abs() < CROSS_PRECISION {
            break;
        }
    }
    Ok(jd)
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
    let mut jd = jd_et + dist / (360.0 / 27.32);
    loop {
        let r = eph.calc(jd, body, flags)?;
        let dist = diff_degrees(x2cross, r.data[0]);
        jd += dist / r.data[3];
        if dist.abs() < CROSS_PRECISION {
            break;
        }
    }
    Ok(jd)
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
    let mut jd = jd_ut + dist / (360.0 / 27.32);
    loop {
        let r = eph.calc_ut(jd, body, flags)?;
        let dist = diff_degrees(x2cross, r.data[0]);
        jd += dist / r.data[3];
        if dist.abs() < CROSS_PRECISION {
            break;
        }
    }
    Ok(jd)
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
        let mut x = eph.calc(jd, body, flags)?.data;
        if (x[1] >= 0.0 && xlat < 0.0) || (x[1] < 0.0 && xlat > 0.0) {
            let mut dist = x[1];
            loop {
                jd -= dist / x[4];
                x = eph.calc(jd, body, flags)?.data;
                dist = x[1];
                if dist.abs() < CROSS_PRECISION {
                    return Ok(MoonCrossing {
                        jd,
                        longitude: x[0],
                        latitude: x[1],
                    });
                }
            }
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
        let mut x = eph.calc_ut(jd, body, flags)?.data;
        if (x[1] >= 0.0 && xlat < 0.0) || (x[1] < 0.0 && xlat > 0.0) {
            let mut dist = x[1];
            loop {
                jd -= dist / x[4];
                x = eph.calc_ut(jd, body, flags)?.data;
                dist = x[1];
                if dist.abs() < CROSS_PRECISION {
                    return Ok(MoonCrossing {
                        jd,
                        longitude: x[0],
                        latitude: x[1],
                    });
                }
            }
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
    let mut jd = if dir >= 0 {
        let dist = normalize_degrees(x2cross - r.data[0]);
        jd_et + dist / xlp
    } else {
        let dist = 360.0 - normalize_degrees(x2cross - r.data[0]);
        jd_et - dist / xlp
    };
    loop {
        let r = eph.calc(jd, body, flags)?;
        let dist = diff_degrees(x2cross, r.data[0]);
        jd += dist / r.data[3];
        if dist.abs() < CROSS_PRECISION {
            break;
        }
    }
    Ok(jd)
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
    let mut jd = if dir >= 0 {
        let dist = normalize_degrees(x2cross - r.data[0]);
        jd_ut + dist / xlp
    } else {
        let dist = 360.0 - normalize_degrees(x2cross - r.data[0]);
        jd_ut - dist / xlp
    };
    loop {
        let r = eph.calc_ut(jd, body, flags)?;
        let dist = diff_degrees(x2cross, r.data[0]);
        jd += dist / r.data[3];
        if dist.abs() < CROSS_PRECISION {
            break;
        }
    }
    Ok(jd)
}
