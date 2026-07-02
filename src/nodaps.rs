//! Nodes & apsides — `swe_nod_aps` / `swe_nod_aps_ut`.
//!
//! Standalone public API for the ascending/descending nodes and the
//! perihelion/aphelion (apogee) of any body. Two families:
//!
//! * **Mean** elements (this module, PNOC 4): VSOP-style mean-equinox-of-date
//!   polynomials for Sun..Neptune / Earth, and `swi_mean_lunar_elements` for the
//!   Moon. Implemented here.
//! * **Osculating** elements (`SE_NODBIT_OSCU` / `SE_NODBIT_OSCU_BAR`): the true
//!   instantaneous two-body ellipse — not yet implemented (PNOC 5, swisseph-rs/86).
//!
//! Both families share the [`transform_nodaps_output`] pipeline (C `swe_nod_aps`
//! A.5), which takes the four raw node/apsis vectors in heliocentric
//! ecliptic-of-date cartesian and produces the observer-relative apparent output
//! (light deflection, aberration, precession, nutation, sidereal, frame/units).
//!
//! Reference: `docs/c-ref-nodaps.md` Parts A, B; C `swecl.c:5075-5665`.

use bitflags::bitflags;

use crate::calc::{app_pos_rest, extract_output, precess_speed};
use crate::constants::{
    AUNIT, CLIGHT, DEGTORAD, IPL_TO_ELEM, J2000, MOON_MEAN_DIST, MOON_MEAN_ECC, MOON_MEAN_INCL,
    NUT_SPEED_INTV,
};
use crate::context::Ephemeris;
use crate::corrections::{aberr_light, deflect_light};
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{cotrans, polar_to_cartesian_with_speed, rotate_x_sincos};
use crate::nutation::nutation;
use crate::obliquity::obliquity;
use crate::precession::precess;
use crate::types::{AstroModels, Body, PrecessionDirection};

bitflags! {
    /// Method selector for [`Ephemeris::nod_aps`](crate::Ephemeris::nod_aps),
    /// mirroring C's `SE_NODBIT_*` (swephexp.h:291-294). Combine bits with `|`.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct NodApsMethod: u32 {
        /// Mean nodes/apsides (VSOP mean elements). Also the behaviour when the
        /// method is empty.
        const MEAN     = 1;
        /// Osculating nodes/apsides about the Sun (or geocentric for the Moon).
        const OSCU     = 2;
        /// Osculating about the barycenter (bodies beyond ~6 AU only; heliocentric
        /// otherwise).
        const OSCU_BAR = 4;
        /// Return the ellipse's second focal point instead of the aphelion.
        const FOPOINT  = 256;
    }
}

/// Nodes & apsides output — four apparent state vectors, each
/// `[lon, lat, dist, dlon, dlat, ddist]` (or equatorial / cartesian per the
/// request flags), matching C's `xnasc`/`xndsc`/`xperi`/`xaphe`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodesApsides {
    /// Ascending node.
    pub ascending: [f64; 6],
    /// Descending node.
    pub descending: [f64; 6],
    /// Perihelion (or perigee).
    pub perihelion: [f64; 6],
    /// Aphelion (or apogee), or the ellipse's 2nd focal point if
    /// [`NodApsMethod::FOPOINT`] was requested.
    pub aphelion: [f64; 6],
}

/// Observer / origin geometry at one epoch, in equatorial-J2000 cartesian
/// (pos + speed, AU / AU-day). Built per backend by
/// [`Ephemeris::nodaps_observer`](crate::Ephemeris) and consumed by
/// [`transform_nodaps_output`]. Mirrors C's `xsun`/`xear`/`xobs` globals
/// (swecl.c A.5.1).
pub(crate) struct ObsFrame {
    /// Barycentric Sun (`swed.pldat[SEI_SUNBARY].x`). All-zero for Moshier
    /// (no barycenter).
    pub sun_bary: [f64; 6],
    /// Earth in the node's heliocentric/barycentric frame, WITHOUT the
    /// topocentric offset (`swed.pldat[SEI_EARTH].x`).
    pub xear: [f64; 6],
    /// Observer position = `xear` + topocentric offset. Also serves as the
    /// heliocentric-observer vector for the deflection geometry (matching
    /// `calc_planet`, which passes its `xobs` to `deflect_light`).
    pub xobs: [f64; 6],
}

// ---------------------------------------------------------------------------
// A.0 — VSOP mean-equinox-of-date element tables (swecl.c:5012-5040)
//
// Each row is `[c0, c1, c2, c3]`, evaluated as `c0 + c1·t + c2·t² + c3·t³`
// with `t = (tjd_et − J2000) / 36525`. Rows 0..7 = Mercury, Venus, Earth,
// Mars, Jupiter, Saturn, Uranus, Neptune. Earth's node/incl rows are all-zero
// (no sensible ecliptic node for Earth's own orbit).
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const EL_NODE: [[f64; 4]; 8] = [
    [ 48.330893,  1.1861890,  0.00017587,  0.000000211], // Mercury
    [ 76.679920,  0.9011190,  0.00040665, -0.000000080], // Venus
    [  0.0,       0.0,        0.0,         0.0],          // Earth
    [ 49.558093,  0.7720923,  0.00001605,  0.000002325], // Mars
    [100.464441,  1.0209550,  0.00040117,  0.000000569], // Jupiter
    [113.665524,  0.8770970, -0.00012067, -0.000002380], // Saturn
    [ 74.005947,  0.5211258,  0.00133982,  0.000018516], // Uranus
    [131.784057,  1.1022057,  0.00026006, -0.000000636], // Neptune
];

#[rustfmt::skip]
const EL_PERI: [[f64; 4]; 8] = [
    [ 77.456119,  1.5564775,  0.00029589,  0.000000056], // Mercury
    [131.563707,  1.4022188, -0.00107337, -0.000005315], // Venus
    [102.937348,  1.7195269,  0.00045962,  0.000000499], // Earth
    [336.060234,  1.8410331,  0.00013515,  0.000000318], // Mars
    [ 14.331309,  1.6126668,  0.00103127, -0.000004569], // Jupiter
    [ 93.056787,  1.9637694,  0.00083757,  0.000004899], // Saturn
    [173.005159,  1.4863784,  0.00021450,  0.000000433], // Uranus
    [ 48.123691,  1.4262677,  0.00037918, -0.000000003], // Neptune
];

#[rustfmt::skip]
const EL_INCL: [[f64; 4]; 8] = [
    [  7.004986,  0.0018215, -0.00001809,  0.000000053], // Mercury
    [  3.394662,  0.0010037, -0.00000088, -0.000000007], // Venus
    [  0.0,       0.0,        0.0,         0.0],          // Earth
    [  1.849726, -0.0006010,  0.00001276, -0.000000006], // Mars
    [  1.303270, -0.0054966,  0.00000465, -0.000000004], // Jupiter
    [  2.488878, -0.0037363, -0.00001516,  0.000000089], // Saturn
    [  0.773196,  0.0007744,  0.00003749, -0.000000092], // Uranus
    [  1.769952, -0.0093082, -0.00000708,  0.000000028], // Neptune
];

#[rustfmt::skip]
const EL_ECCE: [[f64; 4]; 8] = [
    [ 0.20563175,  0.000020406, -0.0000000284, -0.00000000017], // Mercury
    [ 0.00677188, -0.000047766,  0.0000000975,  0.00000000044], // Venus
    [ 0.01670862, -0.000042037, -0.0000001236,  0.00000000004], // Earth
    [ 0.09340062,  0.000090483, -0.0000000806, -0.00000000035], // Mars
    [ 0.04849485,  0.000163244, -0.0000004719, -0.00000000197], // Jupiter
    [ 0.05550862, -0.000346818, -0.0000006456,  0.00000000338], // Saturn
    [ 0.04629590, -0.000027337,  0.0000000790,  0.00000000025], // Uranus
    [ 0.00898809,  0.000006408, -0.0000000008, -0.00000000005], // Neptune
];

#[rustfmt::skip]
const EL_SEMA: [[f64; 4]; 8] = [
    [  0.387098310, 0.0,          0.0,           0.0], // Mercury
    [  0.723329820, 0.0,          0.0,           0.0], // Venus
    [  1.000001018, 0.0,          0.0,           0.0], // Earth
    [  1.523679342, 0.0,          0.0,           0.0], // Mars
    [  5.202603191, 0.0000001913, 0.0,           0.0], // Jupiter
    [  9.554909596, 0.0000021389, 0.0,           0.0], // Saturn
    [ 19.218446062, -0.0000000372, 0.00000000098, 0.0], // Uranus
    [ 30.110386869, -0.0000001663, 0.00000000069, 0.0], // Neptune
];

/// Evaluate a `[c0, c1, c2, c3]` element row as `c0 + c1·t + c2·t² + c3·t³`,
/// matching C's explicit `ep[0] + ep[1]*t + ep[2]*t*t + ep[3]*t*t*t` order.
#[inline]
fn el_poly(ep: &[f64; 4], t: f64) -> f64 {
    ep[0] + ep[1] * t + ep[2] * t * t + ep[3] * t * t * t
}

/// `swe_nod_aps` (swecl.c:5075-5654) — mean branch only (PNOC 4). Osculating
/// (`OSCU`/`OSCU_BAR`) returns a not-yet-implemented error until PNOC 5.
pub(crate) fn nod_aps(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    flags: CalcFlags,
    method: NodApsMethod,
) -> Result<NodesApsides, Error> {
    // A.1 — special-case remap: asteroid-number Pluto -> SE_PLUTO (swecl.c:5116).
    let ipl = match ipl {
        Body::Asteroid(id) if id.mpc_number() == 134340 => Body::Pluto,
        other => other,
    };

    // A.1 — reject the node/apsis point bodies themselves and reserved ids.
    let raw = ipl.to_raw_id();
    if matches!(
        ipl,
        Body::MeanNode | Body::TrueNode | Body::MeanApogee | Body::OscuApogee
    ) || raw < 0
    {
        return Err(Error::InvalidBody(raw));
    }

    // A.2 — setup. Strip the JPL-Horizons approximation bits (swecl.c:5121).
    let flags = flags & !(CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX);
    let t = (tjd_et - J2000) / 36525.0;
    let do_focal_point = method.contains(NodApsMethod::FOPOINT);
    let method = method - NodApsMethod::FOPOINT; // strip FOPOINT for dispatch
    let ipli = if ipl == Body::Sun { Body::Earth } else { ipl };

    let mut do_aberr = !flags.intersects(CalcFlags::TRUEPOS | CalcFlags::NOABERR);
    let mut do_defl = !flags.contains(CalcFlags::TRUEPOS) && !flags.contains(CalcFlags::NOGDEFL);
    if ipl == Body::Moon {
        do_defl = false;
        if !flags.contains(CalcFlags::HELCTR) {
            do_aberr = false;
        }
    }

    // A.3 — mean-branch eligibility: Sun..Neptune (raw 0..8) or Earth, and the
    // method is MEAN or unspecified. Pluto (raw 9) and everything else fall to
    // the osculating branch regardless of method.
    let mean_eligible = (0..=8).contains(&raw) || raw == 14;
    let use_mean = (method.is_empty() || method.contains(NodApsMethod::MEAN)) && mean_eligible;

    let models = &eph.config().astro_models;

    if !use_mean {
        // OSCU / OSCU_BAR (and Pluto/asteroids/fictitious) — PNOC 5.
        return Err(Error::CError(
            "swe_nod_aps: osculating nodes/apsides (SE_NODBIT_OSCU / SE_NODBIT_OSCU_BAR) not yet \
             implemented — PNOC 5 (swisseph-rs/86)"
                .to_string(),
        ));
    }

    // A.3 — build the four raw node/apsis vectors (heliocentric ecliptic-of-date
    // cartesian, pos + speed).
    let mut points = mean_branch(ipl, t, do_focal_point);

    // A.5 — shared observer/apparent-position + output pipeline.
    let outputs = transform_nodaps_output(
        eph,
        &mut points,
        /* is_true_nodaps = */ false,
        ipl,
        ipli,
        flags,
        do_defl,
        do_aberr,
        /* ellipse_is_bary = */ false,
        models,
        tjd_et,
    )?;

    Ok(NodesApsides {
        ascending: outputs[0],
        descending: outputs[1],
        perihelion: outputs[2],
        aphelion: outputs[3],
    })
}

/// A.3 — mean node/apsis vectors for a mean-eligible body. Returns the four
/// six-vectors `[ascending, descending, perihelion, aphelion]` in heliocentric
/// ecliptic-of-date cartesian (pos + speed). Faithful port of swecl.c:5161-5245.
fn mean_branch(ipl: Body, t: f64, do_focal_point: bool) -> [[f64; 6]; 4] {
    let mut xna = [0.0f64; 6];
    let mut xnd = [0.0f64; 6];
    let mut xpe = [0.0f64; 6];
    let mut xap = [0.0f64; 6];

    let (incl, vincl, ecce, vecce, sema, vsema);

    if ipl == Body::Moon {
        // A.3.1 — Moon: numerical mean-element longitudes + speeds (degrees).
        let (node, dnode, peri, dperi) = crate::calc::mean_lunar_elements(
            // tjd_et back out of t: t = (tjd - J2000)/36525.
            t * 36525.0 + J2000,
        );
        xna[0] = node;
        xna[3] = dnode;
        xpe[0] = peri;
        xpe[3] = dperi;
        incl = MOON_MEAN_INCL;
        vincl = 0.0;
        ecce = MOON_MEAN_ECC;
        vecce = 0.0;
        sema = MOON_MEAN_DIST / AUNIT;
        vsema = 0.0;
    } else {
        // A.3.2 — planets (and the Sun via Earth's row): 4-term polynomials.
        let iplx = IPL_TO_ELEM[ipl.to_raw_id() as usize];
        incl = el_poly(&EL_INCL[iplx], t);
        vincl = EL_INCL[iplx][1] / 36525.0;
        sema = el_poly(&EL_SEMA[iplx], t);
        vsema = EL_SEMA[iplx][1] / 36525.0;
        ecce = el_poly(&EL_ECCE[iplx], t);
        vecce = EL_ECCE[iplx][1] / 36525.0;
        xna[0] = el_poly(&EL_NODE[iplx], t);
        xna[3] = EL_NODE[iplx][1] / 36525.0;
        xpe[0] = el_poly(&EL_PERI[iplx], t);
        xpe[3] = EL_PERI[iplx][1] / 36525.0;
    }

    // A.3.3 — shared post-processing (degrees throughout until the final
    // polar->cartesian). `parg`/`pargx` are the arg-of-perihelion (and its
    // speed) FROM the node; kept as locals because xpe[0]/xpe[3] are overwritten.
    xnd[0] = normalize_deg(xna[0] + 180.0);
    xnd[3] = xna[3];
    let parg = normalize_deg(xpe[0] - xna[0]);
    xpe[0] = parg;
    let pargx = normalize_deg(xpe[0] + xpe[3] - xna[3]);
    xpe[3] = pargx;

    // Rotate the perihelion direction from the orbital plane to the mean ecliptic
    // of date (xpe[0..3] as position, xpe[3..6] as an AUXILIARY position — "not a
    // speed" per the C comment). cotrans preserves the distance component.
    let r0 = cotrans([xpe[0], xpe[1], xpe[2]], -incl);
    xpe[0] = r0[0];
    xpe[1] = r0[1];
    let r1 = cotrans([xpe[3], xpe[4], xpe[5]], -incl - vincl);
    xpe[3] = r1[0];
    xpe[4] = r1[1];

    xpe[0] = normalize_deg(xpe[0] + xna[0]);
    xpe[3] = normalize_deg(xpe[3] + xna[0] + xna[3]);
    xpe[3] = normalize_deg(xpe[3] - xpe[0]);

    xpe[2] = sema * (1.0 - ecce);
    xpe[5] = (sema + vsema) * (1.0 - ecce - vecce) - xpe[2];

    xap[0] = normalize_deg(xpe[0] + 180.0);
    xap[1] = -xpe[1];
    xap[3] = xpe[3];
    xap[4] = -xpe[4];
    if do_focal_point {
        xap[2] = sema * ecce * 2.0;
        xap[5] = (sema + vsema) * (ecce + vecce) * 2.0 - xap[2];
    } else {
        xap[2] = sema * (1.0 + ecce);
        xap[5] = (sema + vsema) * (1.0 + ecce + vecce) - xap[2];
    }

    // Node / descending-node distances from the osculating ellipse (the ellipse's
    // radius in the node direction, not just `sema`). swecl.c:5223-5240.
    let ea = (((-parg) * DEGTORAD / 2.0).tan() * ((1.0 - ecce) / (1.0 + ecce)).sqrt()).atan() * 2.0;
    let eax = (((-pargx) * DEGTORAD / 2.0).tan()
        * ((1.0 - ecce - vecce) / (1.0 + ecce + vecce)).sqrt())
    .atan()
        * 2.0;
    xna[2] = sema * (ea.cos() - ecce) / (parg * DEGTORAD).cos();
    xna[5] = (sema + vsema) * (eax.cos() - ecce - vecce) / (pargx * DEGTORAD).cos();
    xna[5] -= xna[2];

    let ea = (((180.0 - parg) * DEGTORAD / 2.0).tan() * ((1.0 - ecce) / (1.0 + ecce)).sqrt())
        .atan()
        * 2.0;
    let eax = (((180.0 - pargx) * DEGTORAD / 2.0).tan()
        * ((1.0 - ecce - vecce) / (1.0 + ecce + vecce)).sqrt())
    .atan()
        * 2.0;
    xnd[2] = sema * (ea.cos() - ecce) / ((180.0 - parg) * DEGTORAD).cos();
    xnd[5] = (sema + vsema) * (eax.cos() - ecce - vecce) / ((180.0 - pargx) * DEGTORAD).cos();
    xnd[5] -= xnd[2];

    // Degrees -> radians -> cartesian (with speed), for all four points.
    let mut points = [xna, xnd, xpe, xap];
    for xp in &mut points {
        xp[0] *= DEGTORAD;
        xp[1] *= DEGTORAD;
        xp[3] *= DEGTORAD;
        xp[4] *= DEGTORAD;
        *xp = polar_to_cartesian_with_speed(*xp);
    }
    points
}

/// A.5 — shared output-transform pipeline (swecl.c:5401-5652). Takes the four raw
/// node/apsis vectors in heliocentric ecliptic-of-date cartesian and produces the
/// four apparent output state vectors (`[lon, lat, dist, dlon, dlat, ddist]` or
/// equatorial / cartesian per `flags`). `is_true_nodaps` gates the osculating-only
/// nutation-to-equator steps; the mean branch passes `false`.
///
/// Reused unchanged by PNOC 5's osculating branch (which produces its four raw
/// vectors differently and passes `is_true_nodaps = true`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn transform_nodaps_output(
    eph: &Ephemeris,
    points: &mut [[f64; 6]; 4],
    is_true_nodaps: bool,
    ipl: Body,
    ipli: Body,
    flags: CalcFlags,
    do_defl: bool,
    do_aberr: bool,
    ellipse_is_bary: bool,
    models: &AstroModels,
    tjd_et: f64,
) -> Result<[[f64; 6]; 4], Error> {
    let has_speed = flags.contains(CalcFlags::SPEED);
    let is_moseph = eph.config().ephemeris_source == crate::types::EphemerisSource::Moshier;

    // A.5.1 — obliquity frame: J2000 if requested, else of-date.
    let oe = if flags.contains(CalcFlags::J2000) {
        obliquity(J2000, flags, models)
    } else {
        obliquity(tjd_et, flags, models)
    };

    // Observer frame at tjd_et (xsun / xear / xobs / earth_helio).
    let frame0 = eph.nodaps_observer(tjd_et, flags)?;

    // Nutation of date, needed by is_true_nodaps steps and by the app_pos_rest tail.
    let nut_val = nutation(tjd_et, flags, models);

    let mut out = [[0.0f64; 6]; 4];

    for (ij, xp) in points.iter_mut().enumerate() {
        // Earth itself has no ascending/descending node.
        if ipli == Body::Earth && ij <= 1 {
            *xp = [0.0; 6];
            out[ij] = extract_output(&[0.0; 24], flags);
            continue;
        }

        // --- to equator ---
        if is_true_nodaps && !flags.contains(CalcFlags::NONUT) {
            // Remove the ecliptic-nutation rotation (about x by -Δε).
            let (sn, cn) = (nut_val.deps.sin(), nut_val.deps.cos());
            let p = rotate_x_sincos([xp[0], xp[1], xp[2]], -sn, cn);
            let v = rotate_x_sincos([xp[3], xp[4], xp[5]], -sn, cn);
            xp[0..3].copy_from_slice(&p);
            xp[3..6].copy_from_slice(&v);
        }
        // Ecliptic -> equatorial (rotate by -obliquity), pos + speed.
        let p = rotate_x_sincos([xp[0], xp[1], xp[2]], -oe.sin_eps, oe.cos_eps);
        let v = rotate_x_sincos([xp[3], xp[4], xp[5]], -oe.sin_eps, oe.cos_eps);
        xp[0..3].copy_from_slice(&p);
        xp[3..6].copy_from_slice(&v);
        if is_true_nodaps && !flags.contains(CalcFlags::NONUT) {
            // Remove the full nutation matrix (swi_nutate backward=true).
            crate::calc::nutate(xp, &oe, &nut_val, None, has_speed, true);
        }

        // --- to J2000 (always) ---
        precess_vec(
            xp,
            tjd_et,
            flags,
            models,
            PrecessionDirection::DateToJ2000,
            has_speed,
        );

        // --- to barycenter ---
        if ipli == Body::Moon {
            for (v, a) in xp.iter_mut().zip(frame0.xear.iter()) {
                *v += a;
            }
        } else if !is_moseph && !ellipse_is_bary {
            for (v, a) in xp.iter_mut().zip(frame0.sun_bary.iter()) {
                *v += a;
            }
        }

        // --- to observer (geocenter / topocenter) ---
        for (v, o) in xp.iter_mut().zip(frame0.xobs.iter()) {
            *v -= o;
        }
        if ipl == Body::Sun && !flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
            for v in xp.iter_mut() {
                *v = -*v;
            }
        }

        // --- light deflection ---
        let r = (xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]).sqrt();
        let dt = r * AUNIT / CLIGHT / 86400.0;
        if do_defl {
            // The node's heliocentric position = geocentric node + heliocentric
            // earth (mirrors calc_planet's planet_helio_retarded).
            let mut planet_helio = [0.0f64; 6];
            for i in 0..6 {
                planet_helio[i] = xp[i] + frame0.xobs[i];
            }
            deflect_light(xp, &frame0.xobs, &planet_helio, has_speed);
        }

        // --- aberration ---
        if do_aberr {
            aberr_light(
                xp,
                &[frame0.xobs[3], frame0.xobs[4], frame0.xobs[5]],
                has_speed,
            );
            if has_speed {
                // Observer-velocity change between emission (t-dt) and reception.
                let frame_ret = eph.nodaps_observer(tjd_et - dt, flags)?;
                for i in 0..3 {
                    xp[i + 3] += frame0.xobs[i + 3] - frame_ret.xobs[i + 3];
                }
            }
        }

        if !has_speed {
            xp[3] = 0.0;
            xp[4] = 0.0;
            xp[5] = 0.0;
        }

        // --- save the J2000-frame copy for the sidereal rigorous branches ---
        let x2000 = *xp;

        // --- precession back to date (unless J2000 requested) ---
        let eps_tail = if !flags.contains(CalcFlags::J2000) {
            precess_vec(
                xp,
                tjd_et,
                flags,
                models,
                PrecessionDirection::J2000ToDate,
                has_speed,
            );
            oe // obliquity of date
        } else {
            oe // obliquity of J2000
        };

        // --- app_pos_rest tail: nutation, equ save, ecliptic, polar, degrees ---
        let nutv = if has_speed {
            Some(nutation(tjd_et - NUT_SPEED_INTV, flags, models))
        } else {
            None
        };
        let mut xreturn = app_pos_rest(xp, flags, &eps_tail, &nut_val, nutv.as_ref());

        // --- sidereal projection (ECL_T0 / SSY_PLANE rigorous, or ayanamsa) ---
        if flags.contains(CalcFlags::SIDEREAL) {
            eph.apply_sidereal(&mut xreturn, &x2000, tjd_et, flags)?;
        }

        out[ij] = extract_output(&xreturn, flags);
    }

    Ok(out)
}

/// Precess a 6-vector (pos always, speed only when `has_speed`) in the given
/// direction, matching the calc pipeline's `precess` + `precess_speed` pairing.
#[inline]
fn precess_vec(
    xp: &mut [f64; 6],
    tjd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    direction: PrecessionDirection,
    has_speed: bool,
) {
    let mut pos3 = [xp[0], xp[1], xp[2]];
    precess(&mut pos3, tjd, flags, models, direction);
    xp[0..3].copy_from_slice(&pos3);
    if has_speed {
        precess_speed(xp, tjd, flags, models, direction);
    }
}

/// `swe_degnorm`-equivalent — normalize degrees to `[0, 360)`.
#[inline]
fn normalize_deg(x: f64) -> f64 {
    crate::math::normalize_degrees(x)
}
