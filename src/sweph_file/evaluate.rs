use std::f64::consts::TAU;

use crate::error::Error;
use crate::math::{chebyshev_deriv, chebyshev_eval};

use super::SwissEphFile;
use super::segment::unpack_segment;
use super::types::{SEI_FLG_ELLIPSE, SEI_FLG_ROTATE, SEI_MOON};

// Verbatim sin/cos of the J2000 obliquity from the C source; full digits preserved.
#[allow(clippy::excessive_precision)]
const SEPS2000: f64 = 0.39777715572793088;
#[allow(clippy::excessive_precision)]
const CEPS2000: f64 = 0.91748206215761929;

pub fn evaluate_body(
    file: &SwissEphFile,
    body_id: i32,
    jd: f64,
    need_speed: bool,
) -> Result<([f64; 6], usize), Error> {
    let planet = file
        .planet_data(body_id)
        .ok_or(Error::InvalidBody(body_id))?;
    if jd < planet.tfstart || jd > planet.tfend {
        return Err(Error::BeyondEphemerisLimits {
            jd_tt: jd,
            start: planet.tfstart,
            end: planet.tfend,
        });
    }

    let mut seg = unpack_segment(file.bytes(), planet, jd, file.header().byte_order)?;

    let is_moon = body_id == SEI_MOON;
    let ncoe = planet.ncoe;

    let neval = if planet.iflg & SEI_FLG_ROTATE != 0 {
        rot_back(&mut seg.coeffs, planet, ncoe, seg.tseg0, is_moon)
    } else {
        ncoe
    };

    let t = ((jd - seg.tseg0) / planet.dseg * 2.0 - 1.0).clamp(-1.0, 1.0);

    let mut result = [0.0; 6];
    for i in 0..3 {
        let offset = i * ncoe;
        let slice = &seg.coeffs[offset..offset + neval];
        result[i] = chebyshev_eval(t, slice);
        if need_speed {
            result[i + 3] = chebyshev_deriv(t, slice) / planet.dseg * 2.0;
        }
    }

    Ok((result, neval))
}

fn rot_back(
    coeffs: &mut [f64],
    planet: &super::types::PlanetFileData,
    ncoe: usize,
    tseg0: f64,
    is_moon: bool,
) -> usize {
    let t = tseg0 + planet.dseg / 2.0;
    let tdiff = (t - planet.telem) / 365250.0;

    let (qav, pav) = if is_moon {
        let mut dn = planet.prot + tdiff * planet.dprot;
        // C uses (int) cast (truncation toward zero), not floor
        dn -= (dn / TAU) as i32 as f64 * TAU;
        let qr = planet.qrot + tdiff * planet.dqrot;
        (qr * dn.cos(), qr * dn.sin())
    } else {
        (
            planet.qrot + tdiff * planet.dqrot,
            planet.prot + planet.dprot * tdiff,
        )
    };

    if planet.iflg & SEI_FLG_ELLIPSE != 0
        && let Some(ref refep) = planet.refep
    {
        let mut omtild = planet.peri + tdiff * planet.dperi;
        omtild -= (omtild / TAU) as i32 as f64 * TAU;
        let (som, com) = omtild.sin_cos();
        for i in 0..ncoe {
            let rx = refep[i];
            let ry = refep[ncoe + i];
            coeffs[i] = coeffs[i] + com * rx - som * ry;
            coeffs[ncoe + i] = coeffs[ncoe + i] + com * ry + som * rx;
        }
    }

    let cosih2 = 1.0 / (1.0 + qav * qav + pav * pav);
    let uix = [
        (1.0 + qav * qav - pav * pav) * cosih2,
        2.0 * qav * pav * cosih2,
        -2.0 * pav * cosih2,
    ];
    let uiy = [
        2.0 * qav * pav * cosih2,
        (1.0 - qav * qav + pav * pav) * cosih2,
        2.0 * qav * cosih2,
    ];
    let uiz = [
        2.0 * pav * cosih2,
        -2.0 * qav * cosih2,
        (1.0 - qav * qav - pav * pav) * cosih2,
    ];

    let mut neval = 0usize;
    for i in 0..ncoe {
        let cx = coeffs[i];
        let cy = coeffs[ncoe + i];
        let cz = coeffs[2 * ncoe + i];
        let xrot = cx * uix[0] + cy * uiy[0] + cz * uiz[0];
        let yrot = cx * uix[1] + cy * uiy[1] + cz * uiz[1];
        let zrot = cx * uix[2] + cy * uiy[2] + cz * uiz[2];
        if xrot.abs() + yrot.abs() + zrot.abs() >= 1e-14 {
            neval = i;
        }
        if is_moon {
            coeffs[i] = xrot;
            coeffs[ncoe + i] = CEPS2000 * yrot - SEPS2000 * zrot;
            coeffs[2 * ncoe + i] = SEPS2000 * yrot + CEPS2000 * zrot;
        } else {
            coeffs[i] = xrot;
            coeffs[ncoe + i] = yrot;
            coeffs[2 * ncoe + i] = zrot;
        }
    }

    neval
}
