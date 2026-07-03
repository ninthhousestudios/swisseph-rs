//! Nodes & apsides — `swe_nod_aps` / `swe_nod_aps_ut`.
//!
//! Standalone public API for the ascending/descending nodes and the
//! perihelion/aphelion (apogee) of any body. Two families:
//!
//! * **Mean** elements (PNOC 4): VSOP-style mean-equinox-of-date polynomials
//!   for Sun..Neptune / Earth, and `swi_mean_lunar_elements` for the Moon.
//! * **Osculating** elements (PNOC 5, `SE_NODBIT_OSCU` / `SE_NODBIT_OSCU_BAR` /
//!   `SE_NODBIT_FOPOINT`): the true instantaneous two-body (angular-momentum)
//!   ellipse, sampled at up to 3 epochs for a central-difference speed.
//!
//! Both families share the [`transform_nodaps_output`] pipeline (C `swe_nod_aps`
//! A.5), which takes the four raw node/apsis vectors in heliocentric
//! ecliptic-of-date cartesian and produces the observer-relative apparent output
//! (light deflection, aberration, precession, nutation, sidereal, frame/units).
//!
//! Reference: `docs/c-ref-nodaps.md` Parts A, B; C `swecl.c:5075-5665`.

use bitflags::bitflags;

use crate::calc::{app_pos_rest, extract_output, plan_for_osc_elem, precess_speed};
use crate::constants::{
    AUNIT, CLIGHT, DEGTORAD, EARTH_MOON_MRAT, GEOGCONST, HELGRAVCONST, IPL_TO_ELEM, J2000,
    MOON_MEAN_DIST, MOON_MEAN_ECC, MOON_MEAN_INCL, NODE_CALC_INTV, NUT_SPEED_INTV,
    OSCU_BAR_DISTANCE_THRESHOLD_AU, PLMASS,
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
    /// Topocentric offset alone (`swi_get_observer`'s output), zero unless
    /// `SEFLG_TOPOCTR` is set. Combined with `sun_bary`/`xear` by
    /// [`select_xobs`] per A.5.1 — NOT pre-added to `xear` here, since the
    /// HELCTR/BARYCTR branches of A.5.1 discard the Earth term entirely.
    pub topo: [f64; 6],
}

/// A.5.1 observer-frame selection (swecl.c:5401-5436). HELCTR (real
/// ephemerides only) and BARYCTR requests observe from the barycentric Sun
/// (HELCTR) or the barycenter/origin (BARYCTR, or HELCTR on Moshier which has
/// no true barycenter) — bypassing the Earth/topocentric offset entirely. A
/// bare `SE_SUN` request (real ephemerides only) also observes from the
/// barycentric Sun, since "the Sun's node/apsis" means Earth's orbital
/// node/apsis mirrored through heliocentric space (swecl.c:5524-5525's sign
/// flip). Every other request observes from Earth (+ topocentric offset).
fn select_xobs(frame: &ObsFrame, flags: CalcFlags, ipl: Body, is_moseph: bool) -> [f64; 6] {
    let mut xobs = if flags.contains(CalcFlags::TOPOCTR) {
        frame.topo
    } else {
        [0.0; 6]
    };
    if flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
        if flags.contains(CalcFlags::HELCTR) && !is_moseph {
            xobs = frame.sun_bary;
        }
    } else if ipl == Body::Sun && !is_moseph {
        xobs = frame.sun_bary;
    } else {
        for (v, a) in xobs.iter_mut().zip(frame.xear.iter()) {
            *v += a;
        }
    }
    xobs
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

    // A.3 (mean) — heliocentric ecliptic-of-date cartesian, pos + speed — or
    // A.4 (osculating): the instantaneous two-body ellipse. Pluto/asteroids/
    // fictitious bodies always fall here (mean_eligible is false for them).
    let (mut points, ellipse_is_bary, is_true_nodaps) = if use_mean {
        (mean_branch(ipl, t, do_focal_point), false, false)
    } else {
        let (points, ellipse_is_bary) =
            osculating_branch(eph, tjd_et, ipli, flags, method, do_focal_point)?;
        (points, ellipse_is_bary, true)
    };

    // A.5 — shared observer/apparent-position + output pipeline.
    let outputs = transform_nodaps_output(
        eph,
        &mut points,
        is_true_nodaps,
        ipl,
        ipli,
        flags,
        do_defl,
        do_aberr,
        ellipse_is_bary,
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

/// A.4 — osculating branch (swecl.c:5249-5400): the instantaneous two-body
/// (angular-momentum) ellipse, sampled at up to 3 epochs (`istart..=iend`) for
/// a central-difference speed. Returns the four raw `[ascending, descending,
/// perihelion, aphelion]` vectors (ecliptic-of-date cartesian, pos + speed)
/// plus `ellipse_is_bary` (for A.5's barycenter-add gate). Faithful port of
/// swecl.c:5249-5399; reference `docs/c-ref-nodaps.md` §A.4.
fn osculating_branch(
    eph: &Ephemeris,
    tjd_et: f64,
    ipli: Body,
    flags: CalcFlags,
    method: NodApsMethod,
    do_focal_point: bool,
) -> Result<([[f64; 6]; 4], bool), Error> {
    let has_speed = flags.contains(CalcFlags::SPEED);

    // A.4.1 — reference (heliocentric) distance, Gmsm, dt, ellipse_is_bary.
    let mut ellipse_is_bary = false;
    let (dt, dzmin, gmsm) = if ipli == Body::Moon {
        (
            NODE_CALC_INTV,
            1e-15,
            GEOGCONST * (1.0 + 1.0 / EARTH_MOON_MRAT) / AUNIT / AUNIT / AUNIT * 86400.0 * 86400.0,
        )
    } else {
        let raw0 = eph.nodaps_osc_body_j2000(tjd_et, ipli, false, flags)?;
        let dist = (raw0[0] * raw0[0] + raw0[1] * raw0[1] + raw0[2] * raw0[2]).sqrt();
        if method.contains(NodApsMethod::OSCU_BAR) && dist > OSCU_BAR_DISTANCE_THRESHOLD_AU {
            ellipse_is_bary = true;
        }
        let raw_id = ipli.to_raw_id();
        let plm = if (2..=9).contains(&raw_id) || raw_id == 14 {
            1.0 / PLMASS[IPL_TO_ELEM[raw_id as usize]]
        } else {
            0.0
        };
        let dt = NODE_CALC_INTV * 10.0 * dist;
        (
            dt,
            1e-15 * dt / NODE_CALC_INTV,
            HELGRAVCONST * (1.0 + plm) / AUNIT / AUNIT / AUNIT * 86400.0 * 86400.0,
        )
    };

    // A.4.2 — up to 3 samples (heliocentric/barycentric J2000 equatorial
    // cartesian, TRUEPOS), rotated into ecliptic-of-date via plan_for_osc_elem.
    let (istart, iend) = if has_speed {
        (0usize, 2usize)
    } else {
        (0usize, 0usize)
    };
    let mut xpos = [[0.0f64; 6]; 3];
    for (i, slot) in xpos.iter_mut().enumerate().take(iend + 1).skip(istart) {
        let t = if istart == iend {
            tjd_et
        } else {
            match i {
                0 => tjd_et - dt,
                2 => tjd_et + dt,
                _ => tjd_et,
            }
        };
        let mut raw = eph.nodaps_osc_body_j2000(t, ipli, ellipse_is_bary, flags)?;
        plan_for_osc_elem(flags, t, &mut raw, &eph.config().astro_models);
        *slot = raw;
    }

    // A.4.3-A.4.4 — per-sample ellipse elements: perihelion/aphelion(-or-2nd-
    // focal-point) direction + ellipse-corrected ascending/descending node
    // distance (replacing A.4.3's tangent-line approximation).
    let mut xq = [[0.0f64; 3]; 3]; // perihelion
    let mut xa = [[0.0f64; 3]; 3]; // aphelion / 2nd focal point
    let mut xn = [[0.0f64; 3]; 3]; // ascending node
    let mut xs = [[0.0f64; 3]; 3]; // descending node
    for i in istart..=iend {
        // A.4.3 — tangent-line node/antinode direction.
        if xpos[i][5].abs() < dzmin {
            xpos[i][5] = dzmin;
        }
        let fac = xpos[i][2] / xpos[i][5];
        let sgn = xpos[i][5] / xpos[i][5].abs();
        let mut xn_tan = [0.0f64; 3];
        for j in 0..3 {
            xn_tan[j] = (xpos[i][j] - fac * xpos[i][j + 3]) * sgn;
        }
        let xs_tan = [-xn_tan[0], -xn_tan[1], -xn_tan[2]];

        // A.4.4 — node longitude direction.
        let rxy0 = (xn_tan[0] * xn_tan[0] + xn_tan[1] * xn_tan[1]).sqrt();
        let cosnode = xn_tan[0] / rxy0;
        let sinnode = xn_tan[1] / rxy0;

        // Inclination from the orbital angular-momentum vector.
        let xnorm = crate::math::cross_prod(
            [xpos[i][0], xpos[i][1], xpos[i][2]],
            [xpos[i][3], xpos[i][4], xpos[i][5]],
        );
        let mut rxy = xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
        let c2 = rxy + xnorm[2] * xnorm[2];
        let mut rxyz = c2.sqrt();
        rxy = rxy.sqrt();
        let sinincl = rxy / rxyz;
        let mut cosincl = (1.0 - sinincl * sinincl).sqrt();
        if xnorm[2] < 0.0 {
            // Retrograde (e.g. 20461 Dioretsa) — A.4.4 only; lunar_osc_elem's
            // D.3 never flips (the Moon's inclination is never retrograde).
            cosincl = -cosincl;
        }

        // Argument of latitude.
        let cosu = xpos[i][0] * cosnode + xpos[i][1] * sinnode;
        let sinu = xpos[i][2] / sinincl;
        let uu = sinu.atan2(cosu);

        // Vis-viva semi-major axis.
        rxyz = (xpos[i][0] * xpos[i][0] + xpos[i][1] * xpos[i][1] + xpos[i][2] * xpos[i][2]).sqrt();
        let v2 = xpos[i][3] * xpos[i][3] + xpos[i][4] * xpos[i][4] + xpos[i][5] * xpos[i][5];
        let sema = 1.0 / (2.0 / rxyz - v2 / gmsm);

        // Eccentricity from specific angular momentum.
        let pp = c2 / gmsm;
        let ecce = (1.0 - pp / sema).sqrt();

        // Eccentric/true anomaly of the body.
        let cos_e = 1.0 / ecce * (1.0 - rxyz / sema);
        let dot = xpos[i][0] * xpos[i][3] + xpos[i][1] * xpos[i][4] + xpos[i][2] * xpos[i][5];
        let sin_e = 1.0 / ecce / (sema * gmsm).sqrt() * dot;
        let ny0 = 2.0 * (((1.0 + ecce) / (1.0 - ecce)).sqrt() * sin_e / (1.0 + cos_e)).atan();

        // Perihelion direction: distance of perihelion from the ascending node.
        let mut q = [
            crate::math::normalize_radians(uu - ny0),
            0.0,
            sema * (1.0 - ecce),
        ];
        q = crate::math::polar_to_cartesian(q);
        q = rotate_x_sincos(q, -sinincl, cosincl);
        q = crate::math::cartesian_to_polar(q);
        q[0] += sinnode.atan2(cosnode);

        // Aphelion, or the ellipse's 2nd focal point (SE_NODBIT_FOPOINT).
        let a = [
            crate::math::normalize_radians(q[0] + std::f64::consts::PI),
            -q[1],
            if do_focal_point {
                sema * ecce * 2.0
            } else {
                sema * (1.0 + ecce)
            },
        ];
        xq[i] = crate::math::polar_to_cartesian(q);
        xa[i] = crate::math::polar_to_cartesian(a);

        // Ellipse-corrected ascending/descending node distance (reusing this
        // sample's ecce/sema/uu), replacing A.4.3's tangent-line approximation.
        let ny_node = crate::math::normalize_radians(ny0 - uu);
        let ny_desc = crate::math::normalize_radians(ny_node + std::f64::consts::PI);
        let cos_e_node =
            (2.0 * ((ny_node / 2.0).tan() / ((1.0 + ecce) / (1.0 - ecce)).sqrt()).atan()).cos();
        let cos_e_desc =
            (2.0 * ((ny_desc / 2.0).tan() / ((1.0 + ecce) / (1.0 - ecce)).sqrt()).atan()).cos();
        let rn = sema * (1.0 - ecce * cos_e_node);
        let rn2 = sema * (1.0 - ecce * cos_e_desc);
        let ro = (xn_tan[0] * xn_tan[0] + xn_tan[1] * xn_tan[1] + xn_tan[2] * xn_tan[2]).sqrt();
        let ro2 = (xs_tan[0] * xs_tan[0] + xs_tan[1] * xs_tan[1] + xs_tan[2] * xs_tan[2]).sqrt();
        for j in 0..3 {
            xn[i][j] = xn_tan[j] * rn / ro;
            xs[i][j] = xs_tan[j] * rn2 / ro2;
        }
    }

    // A.4.5 — assemble output + (central-difference) speed.
    let mut xna = [0.0f64; 6];
    let mut xnd = [0.0f64; 6];
    let mut xpe = [0.0f64; 6];
    let mut xap = [0.0f64; 6];
    for i in 0..3 {
        if has_speed {
            xpe[i] = xq[1][i];
            xpe[i + 3] = (xq[2][i] - xq[0][i]) / dt / 2.0;
            xap[i] = xa[1][i];
            xap[i + 3] = (xa[2][i] - xa[0][i]) / dt / 2.0;
            xna[i] = xn[1][i];
            xna[i + 3] = (xn[2][i] - xn[0][i]) / dt / 2.0;
            xnd[i] = xs[1][i];
            xnd[i + 3] = (xs[2][i] - xs[0][i]) / dt / 2.0;
        } else {
            xpe[i] = xq[0][i];
            xap[i] = xa[0][i];
            xna[i] = xn[0][i];
            xnd[i] = xs[0][i];
        }
    }

    Ok(([xna, xnd, xpe, xap], ellipse_is_bary))
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
    let is_moseph = eph.effective_config(flags, eph.config()).ephemeris_source
        == crate::types::EphemerisSource::Moshier;

    // A.5.1 — obliquity frame: J2000 if requested, else of-date.
    let oe = if flags.contains(CalcFlags::J2000) {
        obliquity(J2000, flags, models)
    } else {
        obliquity(tjd_et, flags, models)
    };

    // Observer frame at tjd_et (xsun / xear / topo), then A.5.1's xobs selection.
    let frame0 = eph.nodaps_observer(tjd_et, flags)?;
    let xobs = select_xobs(&frame0, flags, ipl, is_moseph);
    // `swi_deflect_light` (sweph.c:3743) always reads the TRUE Earth/Sun
    // globals for its geometry, independent of whatever `xobs` was
    // reassigned to for HELCTR/BARYCTR output framing — only `swi_aberr_light`
    // (and the position-shift step) use the reassigned `xobs`. `xear` is
    // already heliocentric for Moshier (`sun_bary` is zero there).
    let mut earth_helio_true = frame0.xear;
    for (v, s) in earth_helio_true.iter_mut().zip(frame0.sun_bary.iter()) {
        *v -= s;
    }

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

        // --- to observer (geocenter / topocenter / heliocenter / barycenter) ---
        for (v, o) in xp.iter_mut().zip(xobs.iter()) {
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
            // earth (mirrors calc_planet's planet_helio_retarded). Uses the
            // TRUE Earth (`earth_helio_true`), not the HELCTR/BARYCTR-selected
            // `xobs` (sweph.c:3743's `swi_deflect_light` ignores the caller's
            // observer-frame choice and always reads the real Earth/Sun-bary
            // globals).
            let mut planet_helio = [0.0f64; 6];
            for i in 0..6 {
                planet_helio[i] = xp[i] + earth_helio_true[i];
            }
            deflect_light(xp, &earth_helio_true, &planet_helio, has_speed);
        }

        // --- aberration ---
        if do_aberr {
            aberr_light(xp, &[xobs[3], xobs[4], xobs[5]], has_speed);
            if has_speed {
                // Observer-velocity change between emission (t-dt) and reception.
                let frame_ret = eph.nodaps_observer(tjd_et - dt, flags)?;
                let xobs_ret = select_xobs(&frame_ret, flags, ipl, is_moseph);
                for i in 0..3 {
                    xp[i + 3] += xobs[i + 3] - xobs_ret[i + 3];
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
