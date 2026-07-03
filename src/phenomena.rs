//! Planetary phenomena: phase angle, illuminated fraction, elongation, apparent diameter,
//! apparent magnitude, and (Moon only) horizontal parallax. Port of `swe_pheno` / `swe_pheno_ut`
//! (swecl.c:3802-4142). See `docs/c-ref-phenomena.md`.
//!
//! Everything routes through [`Ephemeris::calc`](crate::context::Ephemeris::calc) — never a
//! backend directly (enforced by the sutra constraint `app-uses-calc-not-backends:phenomena→*`).

use crate::calc;
use crate::constants::{
    AST_OFFSET, AUNIT, CLIGHT, DEGTORAD, EARTH_RADIUS, J2000, PLANETARY_DIAMETERS, RADTODEG,
};
use crate::context::Ephemeris;
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::dot_prod_unit;
use crate::types::Body;

/// Truncated Euler literal used by the Bowell H-G phase functions (swecl.c:3758, `EULER`). Kept
/// distinct from [`f64::consts::E`] and from Saturn's inline `2.7182818` for bit-fidelity — see
/// the C ref doc §5g/§5k. (The deliberately-truncated value trips `clippy::approx_constant`.)
#[allow(clippy::approx_constant)]
const EULER: f64 = 2.718281828459;

/// Saturn's own shorter (7-sig-fig) inline Euler literal (swecl.c:3986), distinct from [`EULER`]
/// — must not be conflated with it or with [`f64::consts::E`] (C ref doc §5g).
#[allow(clippy::approx_constant)]
const EULER_SATURN: f64 = 2.7182818;

/// Number of built-in bodies with a `mag_elem` row (`NMAG_ELEM = SE_VESTA + 1`, swecl.c:3759).
const NMAG_ELEM: i32 = 21;
/// Boundary between the "planet" and Bowell-H-G magnitude branches (`SE_CHIRON`, swephexp.h:116).
const SE_CHIRON: i32 = 15;

/// `mag_elem[NMAG_ELEM][4]` (swecl.c:3773-3801), transcribed verbatim. Columns are
/// `[H_or_base, G_or_c1, c2, c3]`; a `99` in column 0 is the "no simple formula" sentinel.
/// Rows 2/3 (Mercury/Venus) are inert placeholders — the live code always routes those through
/// the dedicated Mallama polynomial branches — but transcribed for table parity (the C comment
/// says not to delete them).
const MAG_ELEM: [[f64; 4]; 21] = [
    [-26.86, 0.0, 0.0, 0.0],     // 0  Sun
    [-12.55, 0.0, 0.0, 0.0],     // 1  Moon
    [-0.42, 3.80, -2.73, 2.00],  // 2  Mercury (inert placeholder)
    [-4.40, 0.09, 2.39, -0.65],  // 3  Venus (inert placeholder)
    [-1.52, 1.60, 0.0, 0.0],     // 4  Mars
    [-9.40, 0.5, 0.0, 0.0],      // 5  Jupiter
    [-8.88, -2.60, 1.25, 0.044], // 6  Saturn
    [-7.19, 0.0, 0.0, 0.0],      // 7  Uranus
    [-6.87, 0.0, 0.0, 0.0],      // 8  Neptune
    [-1.00, 0.0, 0.0, 0.0],      // 9  Pluto
    [99.0, 0.0, 0.0, 0.0],       // 10 Mean Node
    [99.0, 0.0, 0.0, 0.0],       // 11 True Node
    [99.0, 0.0, 0.0, 0.0],       // 12 Mean Apogee
    [99.0, 0.0, 0.0, 0.0],       // 13 Oscu Apogee
    [99.0, 0.0, 0.0, 0.0],       // 14 Earth
    [6.5, 0.15, 0.0, 0.0],       // 15 Chiron
    [7.0, 0.15, 0.0, 0.0],       // 16 Pholus
    [3.34, 0.12, 0.0, 0.0],      // 17 Ceres
    [4.13, 0.11, 0.0, 0.0],      // 18 Pallas
    [5.33, 0.32, 0.0, 0.0],      // 19 Juno
    [3.20, 0.32, 0.0, 0.0],      // 20 Vesta
];

/// Output of [`pheno`] — `attr[0..5]` of C's `swe_pheno` (swecl.c:3744-3750).
#[derive(Debug, Clone, Copy)]
pub struct Phenomena {
    /// Phase angle (Sun-planet-Earth), degrees. `attr[0]`.
    pub phase_angle: f64,
    /// Illuminated fraction of the disc, 0..1. `attr[1]`.
    pub phase: f64,
    /// Elongation (Sun-Earth-planet angle), degrees. `attr[2]`.
    pub elongation: f64,
    /// Apparent diameter of the disc, degrees. `attr[3]`.
    pub apparent_diameter: f64,
    /// Apparent magnitude. `attr[4]`.
    pub apparent_magnitude: f64,
    /// Geocentric (or topocentric, with `TOPOCTR`) horizontal parallax — Moon only, degrees.
    /// `attr[5]`. Zero for every other body.
    pub horizontal_parallax: f64,
}

/// §1 body remapping (swecl.c:3820-3835): delegates to the shared
/// `normalize_asteroid_aliases` (identical mapping: 134340→Pluto, 1..4→Ceres..Vesta).
fn normalize_pheno_body(body: Body) -> Body {
    crate::calc::normalize_asteroid_aliases(body)
}

/// Compute planetary phenomena at `tjd_et` (Ephemeris/Dynamical Time). Port of `swe_pheno`
/// (swecl.c:3802-4123). Returns the [`Phenomena`] plus the flags actually used (C's `return
/// iflag`, the masked/patched flag copy signalling any ephemeris fallback).
///
/// `Body::Asteroid` (numbered asteroids other than 1566 Icarus) needs the SE1 orbital-element H/G
/// globals (`swed.ast_H`/`ast_G`), which aren't threaded through a stateless config yet — the same
/// gap `eclipse::body_radius_au` stubs — so it returns [`Error::EphemerisNotAvailable`]. Main
/// planets plus Ceres/Pallas/Juno/Vesta (via `MAG_ELEM`) work normally.
pub fn pheno(
    eph: &Ephemeris,
    tjd_et: f64,
    body: Body,
    flags: CalcFlags,
) -> Result<(Phenomena, CalcFlags), Error> {
    // §1 — input sanitization and the two independently-masked flag copies.
    let ipl = normalize_pheno_body(body);
    let raw = ipl.to_raw_id();

    let iflag_mask = calc::EPHMASK
        | CalcFlags::TRUEPOS
        | CalcFlags::J2000
        | CalcFlags::NONUT
        | CalcFlags::NOGDEFL
        | CalcFlags::NOABERR
        | CalcFlags::TOPOCTR;
    let mut iflag = flags & iflag_mask;
    // `iflagp` derives from the already-masked `iflag`, drops NOGDEFL/TOPOCTR, and forces HELCTR.
    let iflagp_mask = calc::EPHMASK
        | CalcFlags::TRUEPOS
        | CalcFlags::J2000
        | CalcFlags::NONUT
        | CalcFlags::NOABERR;
    let mut iflagp = (iflag & iflagp_mask) | CalcFlags::HELCTR;
    let mut epheflag = iflag & calc::EPHMASK;

    let mut attr = [0.0_f64; 6];

    // §2 — geocentric position: cartesian (for dot products) then polar (for distances/lon/lat).
    let geo_xyz = eph.calc(tjd_et, ipl, iflag | CalcFlags::XYZ)?;
    let epheflag2 = geo_xyz.flags_used & calc::EPHMASK;
    if epheflag != epheflag2 {
        // Ephemeris fallback: patch both flag copies to the ephemeris actually used.
        iflag = (iflag & !epheflag) | epheflag2;
        iflagp = (iflagp & !epheflag) | epheflag2;
        epheflag = epheflag2;
    }
    let xx = geo_xyz.data; // geocentric-apparent cartesian at tjd
    let lbr = eph.calc(tjd_et, ipl, iflag)?.data; // geocentric polar (lon, lat, dist AU) at tjd

    // §3 — light-time-corrected heliocentric position → phase angle + illuminated fraction.
    // Skipped for the Sun/Earth/nodes/apogees (attr[0], attr[1], dt stay 0).
    let mut dt = 0.0;
    let mut lbr2 = [0.0_f64; 6]; // heliocentric polar at tjd-dt ("planet-Sun distance" via lbr2[2])
    let skip_phase = matches!(
        ipl,
        Body::Sun
            | Body::Earth
            | Body::MeanNode
            | Body::TrueNode
            | Body::MeanApogee
            | Body::OscuApogee
    );
    if !skip_phase {
        dt = lbr[2] * AUNIT / CLIGHT / 86400.0;
        if iflag.contains(CalcFlags::TRUEPOS) {
            dt = 0.0;
        }
        let xx2 = eph.calc(tjd_et - dt, ipl, iflagp | CalcFlags::XYZ)?.data; // helio cartesian
        lbr2 = eph.calc(tjd_et - dt, ipl, iflagp)?.data; // helio polar at tjd-dt
        attr[0] = dot_prod_unit([xx[0], xx[1], xx[2]], [xx2[0], xx2[1], xx2[2]]).acos() * RADTODEG;
        attr[1] = (1.0 + (attr[0] * DEGTORAD).cos()) / 2.0;
    }

    // §4 — apparent diameter of the disc (uses the geocentric distance lbr[2], not light-time
    // corrected).
    let dd = if raw < NMAG_ELEM {
        PLANETARY_DIAMETERS[raw as usize]
    } else {
        // ipl > SE_AST_OFFSET would read swed.ast_diam (named-asteroid data) — not available in
        // the stateless port; same gap as eclipse::body_radius_au. See §5 for the magnitude side.
        0.0
    };
    if lbr[2] < dd / 2.0 / AUNIT {
        attr[3] = 180.0; // observer inside the body — assume on the surface
    } else {
        attr[3] = (dd / 2.0 / AUNIT / lbr[2]).asin() * 2.0 * RADTODEG;
    }

    // §5 — apparent magnitude. Guard: numbered asteroids beyond the offset, or any built-in body
    // whose mag_elem row isn't the 99 sentinel. Branches are the live MAG_MALLAMA_2018 /
    // MAG_MOON_VREIJS forms (swecl.c:3899-4068); the `#else` halves are dead code.
    let compute_mag =
        raw > AST_OFFSET || ((0..NMAG_ELEM).contains(&raw) && MAG_ELEM[raw as usize][0] < 99.0);
    if compute_mag {
        let a = attr[0];
        if ipl == Body::Sun {
            // §5a
            let mut fac =
                attr[3] / ((PLANETARY_DIAMETERS[0] / 2.0 / AUNIT).asin() * 2.0 * RADTODEG);
            fac *= fac;
            attr[4] = MAG_ELEM[0][0] - 2.5 * fac.log10();
        } else if ipl == Body::Moon {
            // §5b — Allen 1976 below the 147.1385465° stitch, Samaha cube-phase above it.
            if a <= 147.1385465 {
                attr[4] = -21.62 + 0.026 * a.abs() + 0.000000004 * a.powf(4.0);
            } else {
                attr[4] = -4.5444 - (2.5 * (180.0 - a).powf(3.0).log10());
            }
            attr[4] += 5.0 * (lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS).log10();
        } else if ipl == Body::Mercury {
            // §5c — Mallama 2018, powers by repeated multiplication (FP fidelity).
            let a2 = a * a;
            let a3 = a2 * a;
            let a4 = a3 * a;
            let a5 = a4 * a;
            let a6 = a5 * a;
            attr[4] = -0.613 + a * 6.3280E-02 - a2 * 1.6336E-03 + a3 * 3.3644E-05 - a4 * 3.4265E-07
                + a5 * 1.6893E-09
                - a6 * 3.0334E-12;
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
        } else if ipl == Body::Venus {
            // §5d — Mallama 2018, two regimes. The out-of-range (>179°) advisory warning is
            // dropped: this port carries no serr channel and the magnitude is still returned.
            let a2 = a * a;
            let a3 = a2 * a;
            let a4 = a3 * a;
            if a <= 163.7 {
                attr[4] = -4.384 - a * 1.044E-03 + a2 * 3.687E-04 - a3 * 2.814E-06 + a4 * 8.938E-09;
            } else {
                attr[4] = 236.05828 - a * 2.81914E+00 + a2 * 8.39034E-03;
            }
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
        } else if ipl == Body::Mars {
            // §5e — Mallama 2018, two regimes.
            let a2 = a * a;
            if a <= 50.0 {
                attr[4] = -1.601 + a * 0.02267 - a2 * 0.0001302;
            } else {
                attr[4] = -0.367 - a * 0.02573 + a2 * 0.0003445;
            }
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
        } else if ipl == Body::Jupiter {
            // §5f — Mallama 2018, single regime.
            let a2 = a * a;
            attr[4] = -9.395 - a * 3.7E-04 + a2 * 6.16E-04;
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
        } else if ipl == Body::Saturn {
            // §5g — Mallama 2018 + Meeus ring geometry. `T` uses tjd-dt (light-time corrected).
            let t = (tjd_et - dt - J2000) / 36525.0;
            let inc = (28.075216 - 0.012998 * t + 0.000004 * t * t) * DEGTORAD;
            let om = (169.508470 + 1.394681 * t + 0.000412 * t * t) * DEGTORAD;
            let sin_b = inc.sin() * (lbr[1] * DEGTORAD).cos() * (lbr[0] * DEGTORAD - om).sin()
                - inc.cos() * (lbr[1] * DEGTORAD).sin();
            let sin_b2 = inc.sin() * (lbr2[1] * DEGTORAD).cos() * (lbr2[0] * DEGTORAD - om).sin()
                - inc.cos() * (lbr2[1] * DEGTORAD).sin();
            // Angle-averaged (asin → mean → sin), then |·|, per Meeus.
            let sin_b = ((sin_b.asin() + sin_b2.asin()) / 2.0).sin().abs();
            // EULER_SATURN: Saturn's own shorter inline Euler literal, NOT the file-level EULER.
            attr[4] =
                -8.914 - 1.825 * sin_b + 0.026 * a - 0.378 * sin_b * EULER_SATURN.powf(-2.25 * a);
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
        } else if ipl == Body::Uranus {
            // §5h — Mallama 2018, simplified (sub-Earth latitude fi_ ignored, then the empirical
            // -0.05 compensation). Both `fi_` lines kept for literal ordering fidelity.
            let a2 = a * a;
            let fi_ = 0.0;
            attr[4] = -7.110 - 8.4E-04 * fi_ + a * 6.587E-3 + a2 * 1.045E-4;
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
            attr[4] -= 0.05;
        } else if ipl == Body::Neptune {
            // §5i — Mallama 2018, three time regimes (raw tjd for the boundaries). The 0.0055
            // slope is deliberate (keeps the piecewise function continuous), not a typo for 0.0054.
            if tjd_et < 2444239.5 {
                attr[4] = -6.89;
            } else if tjd_et <= 2451544.5 {
                attr[4] = -6.89 - 0.0055 * (tjd_et - 2444239.5) / 365.25;
            } else {
                attr[4] = -7.00;
            }
            attr[4] += 5.0 * (lbr2[2] * lbr[2]).log10();
        } else if raw < SE_CHIRON {
            // §5j — generic old-style polynomial (reached by Pluto: mag_elem[9] = [-1, 0, 0, 0],
            // so this evaluates to 5*log10(r*Δ) - 1.00). Powers by repeated multiplication.
            let me = MAG_ELEM[raw as usize];
            attr[4] = 5.0 * (lbr2[2] * lbr[2]).log10()
                + me[1] * attr[0] / 100.0
                + me[2] * attr[0] * attr[0] / 10000.0
                + me[3] * attr[0] * attr[0] * attr[0] / 1000000.0
                + me[0];
        } else if !(NMAG_ELEM..=AST_OFFSET).contains(&raw) {
            // §5k — Bowell H-G system (Chiron/Pholus/Ceres/Pallas/Juno/Vesta and numbered
            // asteroids). EULER here is the file-level 13-sig-fig literal.
            let ph1 = EULER.powf(-3.33 * (attr[0] * DEGTORAD / 2.0).tan().powf(0.63));
            let ph2 = EULER.powf(-1.87 * (attr[0] * DEGTORAD / 2.0).tan().powf(1.22));
            let me = if raw < NMAG_ELEM {
                [MAG_ELEM[raw as usize][0], MAG_ELEM[raw as usize][1]]
            } else if matches!(ipl, Body::Asteroid(id) if id.mpc_number() == 1566) {
                [16.9, 0.15] // 1566 Icarus: JPL-database H/G override
            } else {
                // Other numbered asteroids need swed.ast_H/ast_G from the SE1 orbital-element
                // file — the same stateless gap as eclipse::body_radius_au's asteroid diameter.
                // TODO(asteroid SE1 metadata): same gap as body_radius_au.
                return Err(Error::EphemerisNotAvailable {
                    body,
                    source: eph.config().ephemeris_source,
                });
            };
            attr[4] = 5.0 * (lbr2[2] * lbr[2]).log10() + me[0]
                - 2.5 * ((1.0 - me[1]) * ph1 + me[1] * ph2).log10();
        } else {
            // §5l — fictitious bodies fallback (unreachable given the guard, present for parity).
            attr[4] = 0.0;
        }
    }

    // §6 — elongation (Sun-Earth-planet angle). Skipped for the Sun/Earth. Uses a fresh Sun
    // cartesian; the C buffer reuse of xx2/lbr2 here is a hazard we avoid with a named binding.
    if ipl != Body::Sun && ipl != Body::Earth {
        let sun_xyz = eph.calc(tjd_et, Body::Sun, iflag | CalcFlags::XYZ)?.data;
        attr[2] = dot_prod_unit([xx[0], xx[1], xx[2]], [sun_xyz[0], sun_xyz[1], sun_xyz[2]]).acos()
            * RADTODEG;
    }

    // §7 — horizontal parallax, Moon only. Uses just the ephemeris-source bits (epheflag).
    if ipl == Body::Moon {
        let xm = eph
            .calc(
                tjd_et,
                Body::Moon,
                epheflag | CalcFlags::TRUEPOS | CalcFlags::EQUATORIAL | CalcFlags::RADIANS,
            )?
            .data;
        let sinhp = EARTH_RADIUS / xm[2] / AUNIT; // xm[2]: true geocentric distance, AU
        attr[5] = sinhp.asin() / DEGTORAD;
        if iflag.contains(CalcFlags::TOPOCTR) {
            // Topocentric: actual angular displacement between geocentric and topocentric apparent
            // directions (uses the Ephemeris's configured topographic position, like C's global).
            let xm_topo = eph
                .calc(
                    tjd_et,
                    Body::Moon,
                    epheflag | CalcFlags::XYZ | CalcFlags::TOPOCTR,
                )?
                .data;
            let xm_geo = eph
                .calc(tjd_et, Body::Moon, epheflag | CalcFlags::XYZ)?
                .data;
            attr[5] = dot_prod_unit(
                [xm_topo[0], xm_topo[1], xm_topo[2]],
                [xm_geo[0], xm_geo[1], xm_geo[2]],
            )
            .acos()
                / DEGTORAD;
        }
    }

    // §8 — return the masked/patched flag copy (C's "flags actually used").
    Ok((
        Phenomena {
            phase_angle: attr[0],
            phase: attr[1],
            elongation: attr[2],
            apparent_diameter: attr[3],
            apparent_magnitude: attr[4],
            horizontal_parallax: attr[5],
        },
        iflag,
    ))
}

/// UT-based wrapper around [`pheno`]. Port of `swe_pheno_ut` (swecl.c:4125-4142): defaults to
/// SWIEPH when no ephemeris is requested, converts UT→TT via deltaT, and re-calls after an
/// ephemeris fallback.
pub fn pheno_ut(
    eph: &Ephemeris,
    tjd_ut: f64,
    body: Body,
    flags: CalcFlags,
) -> Result<(Phenomena, CalcFlags), Error> {
    let mut epheflag = flags & calc::EPHMASK;
    let mut iflag = flags;
    if epheflag.is_empty() {
        epheflag = CalcFlags::SWIEPH;
        iflag |= CalcFlags::SWIEPH;
    }
    let deltat = crate::deltat::calc_deltat(tjd_ut, eph.config());
    let (attr, retflag) = pheno(eph, tjd_ut + deltat, body, iflag)?;
    if (retflag & calc::EPHMASK) != epheflag {
        // Ephemeris fallback: C recomputes deltaT with the actually-used flags. In this stateless
        // port deltaT is derived from the Ephemeris config (not the flags), so the recomputed
        // value is identical and the re-call is idempotent — kept for structural fidelity.
        let deltat = crate::deltat::calc_deltat(tjd_ut, eph.config());
        return pheno(eph, tjd_ut + deltat, body, iflag);
    }
    Ok((attr, retflag))
}
