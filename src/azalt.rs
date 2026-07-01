//! Atmospheric refraction and horizontal-coordinate transforms.
//!
//! Port of `swe_refrac`, `swe_refrac_extended`, `swe_azalt`, `swe_azalt_rev`, and their shared
//! static helpers `calc_astronomical_refr` / `calc_dip` (all `swecl.c`). See
//! `docs/c-ref-refraction-azalt.md`.
//!
//! `azalt`/`azalt_rev` here are the pure geometry cores: they take a precomputed ARMC and true
//! obliquity rather than a UT Julian day. [`crate::context::Ephemeris::azalt`] /
//! [`crate::context::Ephemeris::azalt_rev`] resolve ARMC/eps/deltaT and delegate here, matching
//! how `Ephemeris::houses_ex2` delegates to `houses::houses_armc`.

use crate::constants::{DEGTORAD, EARTH_RADIUS};
use std::f64::consts::PI;

/// Direction for [`refrac`] / [`refrac_extended`]. `SE_TRUE_TO_APP` / `SE_APP_TO_TRUE`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefracDir {
    TrueToApp,
    AppToTrue,
}

/// Input-coordinate direction for [`azalt`]. `SE_ECL2HOR` / `SE_EQU2HOR`. Kept distinct from
/// [`HorDir`] since C reuses the same two integers (0/1) with different meaning per function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AzAltDir {
    EclToHor,
    EquToHor,
}

/// Output-coordinate direction for [`azalt_rev`]. `SE_HOR2ECL` / `SE_HOR2EQU`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HorDir {
    HorToEcl,
    HorToEqu,
}

/// Port of `swe_refrac` (swecl.c:2887-2984). Simple true<->apparent altitude refraction
/// (Meeus formula), sea-level observer with an ideal horizon — no dip, no elevated-observer
/// geometry (see [`refrac_extended`] for that). `atpress` in hPa, `attemp` in deg C.
pub fn refrac(inalt: f64, atpress: f64, attemp: f64, dir: RefracDir) -> f64 {
    let pt_factor = atpress / 1010.0 * 283.0 / (273.0 + attemp);

    match dir {
        RefracDir::TrueToApp => {
            let trualt = inalt;
            let refr = if trualt > 15.0 {
                let a = ((90.0 - trualt) * DEGTORAD).tan();
                (58.276 * a - 0.0824 * a * a * a) * pt_factor / 3600.0
            } else if trualt > -5.0 {
                let a = trualt + 10.3 / (trualt + 5.11);
                let r = if a + 1e-10 >= 90.0 {
                    0.0
                } else {
                    1.02 / (a * DEGTORAD).tan()
                };
                r * pt_factor / 60.0
            } else {
                0.0
            };
            let mut appalt = trualt;
            if appalt + refr > 0.0 {
                appalt += refr;
            }
            appalt
        }
        RefracDir::AppToTrue => {
            let appalt = inalt;
            let a = appalt + 7.31 / (appalt + 4.4);
            let mut refr = if a + 1e-10 >= 90.0 {
                0.0
            } else {
                let r = 1.00 / (a * DEGTORAD).tan();
                r - 0.06 * (14.7 * r + 13.0).sin()
            };
            refr *= pt_factor / 60.0;
            let mut trualt = appalt;
            if appalt - refr > 0.0 {
                trualt = appalt - refr;
            }
            trualt
        }
    }
}

/// Port of `swe_refrac_extended` (swecl.c:3035-3115). Elevated-observer true<->apparent
/// altitude, with horizon dip. `geoalt` = observer height above sea level, meters. `dret[0]` =
/// true altitude, `[1]` = apparent altitude, `[2]` = refraction, `[3]` = dip of horizon, all
/// degrees. The body is above the horizon iff `dret[0] != dret[1]`.
pub fn refrac_extended(
    inalt: f64,
    geoalt: f64,
    atpress: f64,
    attemp: f64,
    lapse_rate: f64,
    dir: RefracDir,
    dret: &mut [f64; 4],
) -> f64 {
    let dip = calc_dip(geoalt, atpress, attemp, lapse_rate);
    let mut inalt = inalt;
    if inalt > 90.0 {
        inalt = 180.0 - inalt;
    }

    match dir {
        RefracDir::TrueToApp => {
            if inalt < -10.0 {
                *dret = [inalt, inalt, 0.0, dip];
                return inalt;
            }

            // 5 fixed Newton iterations inverting calc_astronomical_refr, reusing consecutive
            // evaluations as a secant-like derivative -- "sic !!! code by Moshier" (swecl.c:3064).
            // No convergence check; replicate the loop and variable reuse exactly.
            let mut y = inalt;
            let mut yy0 = 0.0;
            let mut d0 = 0.0;
            let mut refr = 0.0;
            for _ in 0..5 {
                let d = calc_astronomical_refr(y, atpress, attemp);
                refr = d;
                let n = y - yy0;
                let denom = d - d0 - n;
                let n = if n != 0.0 && denom != 0.0 {
                    y - n * (inalt + d - y) / denom
                } else {
                    inalt + d
                };
                yy0 = y;
                d0 = d;
                y = n;
            }

            if inalt + refr < dip {
                *dret = [inalt, inalt, 0.0, dip];
                return inalt;
            }

            *dret = [inalt, inalt + refr, refr, dip];
            inalt + refr
        }
        RefracDir::AppToTrue => {
            let refr = calc_astronomical_refr(inalt, atpress, attemp);
            let trualt = inalt - refr;
            if inalt > dip {
                *dret = [trualt, inalt, refr, dip];
            } else {
                *dret = [inalt, inalt, 0.0, dip];
            }
            // dret-fill uses `inalt > dip` (strict), return uses `inalt >= dip` (inclusive) --
            // intentional asymmetry, "bug fix dieter, 4 feb 20" (swecl.c:3111).
            if inalt >= dip { trualt } else { inalt }
        }
    }
}

/// Sinclair's astronomical refraction formula (swecl.c:3124-3148). `inalt` is an *apparent*
/// altitude.
fn calc_astronomical_refr(inalt: f64, atpress: f64, attemp: f64) -> f64 {
    // 17.904104638432: chosen so the two branches are C0-continuous, not a rounded "15".
    let r = if inalt > 17.904104638432 {
        0.97 / (inalt * DEGTORAD).tan()
    } else {
        (34.46 + 4.23 * inalt + 0.004 * inalt * inalt)
            / (1.0 + 0.505 * inalt + 0.0845 * inalt * inalt)
    };
    ((atpress - 80.0) / 930.0 / (1.0 + 0.00008 * (r + 39.0) * (attemp - 10.0)) * r) / 60.0
}

/// Geometric + refractive dip of the horizon (swecl.c:3158-3169), degrees, negative for
/// `geoalt > 0`. Does NOT auto-estimate `atpress` when `0` -- that estimate is the caller's
/// (`azalt`'s) responsibility, not this helper's (see module docs / c-ref §7).
fn calc_dip(geoalt: f64, atpress: f64, attemp: f64, lapse_rate: f64) -> f64 {
    let krefr = (0.0342 + lapse_rate) / (0.154 * 0.0238);
    let d = 1.0 - 1.8480 * krefr * atpress / (273.15 + attemp) / (273.15 + attemp);
    -180.0 / PI * (1.0 / (1.0 + geoalt / EARTH_RADIUS)).acos() * d.sqrt()
}

/// Ecliptic/equatorial -> azimuth + true/apparent altitude. Pure geometry core of `swe_azalt`
/// (swecl.c:2788-2825 steps 2-8), given a precomputed ARMC and true obliquity. `geopos` =
/// [longitude (unused here), latitude, height above sea (m)]; `xin` = [lon/RA, lat/dec], degrees.
/// Returns `[azimuth (from south, positive clockwise via west), true altitude, apparent
/// altitude]`, degrees.
#[allow(clippy::too_many_arguments)]
pub fn azalt(
    dir: AzAltDir,
    armc: f64,
    eps_true: f64,
    geopos: [f64; 3],
    atpress: f64,
    attemp: f64,
    lapse_rate: f64,
    xin: [f64; 2],
) -> [f64; 3] {
    let mut xra = [xin[0], xin[1], 1.0];
    if dir == AzAltDir::EclToHor {
        xra = crate::math::cotrans(xra, -eps_true);
    }

    let mdd = crate::math::normalize_degrees(xra[0] - armc);
    let mut x = crate::math::cotrans(
        [crate::math::normalize_degrees(mdd - 90.0), xra[1], 1.0],
        90.0 - geopos[1],
    );

    x[0] = crate::math::normalize_degrees(x[0] + 90.0);
    let azimuth = 360.0 - x[0];
    let true_alt = x[1];

    let atpress = if atpress == 0.0 {
        1013.25 * (1.0 - 0.0065 * geopos[2] / 288.0).powf(5.255)
    } else {
        atpress
    };

    let mut dret = [0.0; 4];
    let app_alt = refrac_extended(
        true_alt,
        geopos[2],
        atpress,
        attemp,
        lapse_rate,
        RefracDir::TrueToApp,
        &mut dret,
    );

    [azimuth, true_alt, app_alt]
}

/// Azimuth + true altitude -> equatorial (and optionally ecliptic) coordinates. Pure geometry
/// core of `swe_azalt_rev` (swecl.c:2839-2873), given a precomputed ARMC and true obliquity.
/// Inverse of [`azalt`]'s geometric transform only -- does NOT de-refract; `xin[1]` must already
/// be a true altitude. `xin` = [azimuth (from south, clockwise), true altitude], degrees.
/// Returns [lon/RA, lat/dec], degrees.
pub fn azalt_rev(dir: HorDir, armc: f64, eps_true: f64, geolat: f64, xin: [f64; 2]) -> [f64; 2] {
    let mut xaz = [xin[0], xin[1], 1.0];
    xaz[0] = 360.0 - xaz[0];
    xaz[0] = crate::math::normalize_degrees(xaz[0] - 90.0);

    xaz = crate::math::cotrans(xaz, geolat - 90.0);
    xaz[0] = crate::math::normalize_degrees(xaz[0] + armc + 90.0);

    if dir == HorDir::HorToEqu {
        return [xaz[0], xaz[1]];
    }

    let x = crate::math::cotrans(xaz, eps_true);
    [x[0], x[1]]
}
