//! Observer geocentric offset for topocentric calculations.
//!
//! Low-level internals; exposed for golden tests and advanced use.

use crate::config::{EphemerisConfig, TopoPosition};
use crate::constants::*;
use crate::flags::CalcFlags;
use crate::math::{cartesian_to_polar, polar_to_cartesian_with_speed};
use crate::precession::precess;
use crate::types::{AstroModels, PrecessionDirection};

/// Observer's geocentric offset (position + velocity, AU / AU-day) in the
/// J2000 mean-equatorial frame. Port of `swi_get_observer` restricted to the
/// NONUT-forced mean-frame path exercised by `swe_calc`
/// (docs/c-ref-topocentric.md §3.1–§3.2) — nutation is applied once, later,
/// together with the celestial body, so it must not be applied here too.
pub(crate) fn get_observer(
    jd_tt: f64,
    topo: &TopoPosition,
    flags: CalcFlags,
    config: &EphemerisConfig,
    models: &AstroModels,
) -> [f64; 6] {
    let tjd_ut = jd_tt - crate::deltat::calc_deltat(jd_tt, config);

    let eps = crate::obliquity::obliquity(jd_tt, flags, models);
    let sidt = crate::sidereal_time::sidereal_time0(tjd_ut, eps.eps * RADTODEG, 0.0, config) * 15.0;

    let cosfi = (topo.latitude * DEGTORAD).cos();
    let sinfi = (topo.latitude * DEGTORAD).sin();
    let f = EARTH_OBLATENESS;
    let cc = 1.0 / (cosfi * cosfi + (1.0 - f) * (1.0 - f) * sinfi * sinfi).sqrt();
    let ss = (1.0 - f) * (1.0 - f) * cc;
    let cosl = ((topo.longitude + sidt) * DEGTORAD).cos();
    let sinl = ((topo.longitude + sidt) * DEGTORAD).sin();
    let h = topo.altitude;

    let mut xobs = [
        (EARTH_RADIUS * cc + h) * cosfi * cosl,
        (EARTH_RADIUS * cc + h) * cosfi * sinl,
        (EARTH_RADIUS * ss + h) * sinfi,
        0.0,
        0.0,
        0.0,
    ];

    // swi_cartpol (position only) then swi_polcart_sp (position+speed) — not a
    // single _sp round trip, matching C's asymmetric conversion exactly.
    let polar = cartesian_to_polar([xobs[0], xobs[1], xobs[2]]);
    xobs = polar_to_cartesian_with_speed([polar[0], polar[1], polar[2], EARTH_ROT_SPEED, 0.0, 0.0]);

    for v in xobs.iter_mut() {
        *v /= AUNIT;
    }

    // Frame bias is deliberately neglected for the observer offset (see §3.2 step 9).
    let mut pos3 = [xobs[0], xobs[1], xobs[2]];
    precess(
        &mut pos3,
        jd_tt,
        flags,
        models,
        PrecessionDirection::DateToJ2000,
    );
    xobs[0] = pos3[0];
    xobs[1] = pos3[1];
    xobs[2] = pos3[2];
    crate::calc::precess_speed(
        &mut xobs,
        jd_tt,
        flags,
        models,
        PrecessionDirection::DateToJ2000,
    );

    xobs
}
