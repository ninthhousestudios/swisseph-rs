//! Port of `swehouse.c`: house systems (`CalcH`) and the ARMC-based driver
//! (`swe_houses_armc_ex2`). See `docs/c-ref-houses.md`.

use crate::error::Error;
use crate::math::{cotrans, diff_degrees, normalize_degrees};
use crate::types::HouseSystem;

// ---------------------------------------------------------------------------
// Constants (swehouse.h:87, swehouse.c:68-70, swehouse.c:940)
// ---------------------------------------------------------------------------

const VERY_SMALL: f64 = 1e-10;
const VERY_SMALL_PLAC_ITER: f64 = 1.0 / 360000.0;
#[allow(dead_code)] // used by swe_house_pos (later sub-tasks)
const MILLIARCSEC: f64 = 1.0 / 3600000.0;
const SOLAR_YEAR: f64 = 365.242_198_93;
const ARMCS: f64 = (SOLAR_YEAR + 1.0) / SOLAR_YEAR * 360.0;
const NITER_MAX: i32 = 100;

// ---------------------------------------------------------------------------
// Degree-wrapped trig macros (swehouse.h:89-98)
// ---------------------------------------------------------------------------

fn sind(x: f64) -> f64 {
    (x * crate::constants::DEGTORAD).sin()
}

fn cosd(x: f64) -> f64 {
    (x * crate::constants::DEGTORAD).cos()
}

fn tand(x: f64) -> f64 {
    (x * crate::constants::DEGTORAD).tan()
}

fn atand(x: f64) -> f64 {
    x.atan() * crate::constants::RADTODEG
}

fn asind(x: f64) -> f64 {
    x.asin() * crate::constants::RADTODEG
}

// ---------------------------------------------------------------------------
// Public output types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AscMc {
    pub ascendant: f64,
    pub mc: f64,
    pub armc: f64,
    pub vertex: f64,
    pub equatorial_ascendant: f64,
    pub coascendant_koch: f64,
    pub coascendant_munkasey: f64,
    pub polar_ascendant: f64,
}

impl AscMc {
    pub fn as_array(&self) -> [f64; 8] {
        [
            self.ascendant,
            self.mc,
            self.armc,
            self.vertex,
            self.equatorial_ascendant,
            self.coascendant_koch,
            self.coascendant_munkasey,
            self.polar_ascendant,
        ]
    }
}

#[derive(Debug, Clone)]
pub struct HouseResult {
    pub cusps: [f64; 37],
    pub cusp_speeds: [f64; 37],
    pub ascmc: AscMc,
    pub ascmc_speeds: AscMc,
}

// ---------------------------------------------------------------------------
// Core trig (swehouse.c:2058-2177)
// ---------------------------------------------------------------------------

/// Raw oblique-ascension spherical trig, `x` ∈ [0,90]. Port of `Asc2` (swehouse.c:2100-2129).
///
/// The degenerate-branch order (`sinx==0` checked before `ass==0`) is load-bearing for
/// FP fidelity near `x=0`/`x=180` — see c-ref-houses.md §12.4.
fn asc2(x: f64, f: f64, sine: f64, cose: f64) -> f64 {
    let mut ass = -tand(f) * sine + cose * cosd(x);
    if ass.abs() < VERY_SMALL {
        ass = 0.0;
    }
    let mut sinx = sind(x);
    if sinx.abs() < VERY_SMALL {
        sinx = 0.0;
    }
    if sinx == 0.0 {
        ass = if ass < 0.0 { -VERY_SMALL } else { VERY_SMALL };
    } else if ass == 0.0 {
        ass = if sinx < 0.0 { -90.0 } else { 90.0 };
    } else {
        ass = atand(sinx / ass);
    }
    if ass < 0.0 {
        ass += 180.0;
    }
    ass
}

/// Quadrant-normalized oblique ascension. Port of `Asc1` (swehouse.c:2058-2088).
fn asc1(x1: f64, f: f64, sine: f64, cose: f64) -> f64 {
    let x1 = normalize_degrees(x1);
    let n = (x1 / 90.0 + 1.0) as i32;
    if (90.0 - f).abs() < VERY_SMALL {
        return 180.0;
    }
    if (90.0 + f).abs() < VERY_SMALL {
        return 0.0;
    }
    let mut ass = match n {
        1 => asc2(x1, f, sine, cose),
        2 => 180.0 - asc2(180.0 - x1, -f, sine, cose),
        3 => 180.0 + asc2(x1 - 180.0, -f, sine, cose),
        _ => 360.0 - asc2(360.0 - x1, f, sine, cose),
    };
    ass = normalize_degrees(ass);
    if (ass - 90.0).abs() < VERY_SMALL {
        ass = 90.0;
    }
    if (ass - 180.0).abs() < VERY_SMALL {
        ass = 180.0;
    }
    if (ass - 270.0).abs() < VERY_SMALL {
        ass = 270.0;
    }
    if (ass - 360.0).abs() < VERY_SMALL {
        ass = 0.0;
    }
    ass
}

/// Analytical derivative of `Asc1` w.r.t. armc, scaled to degrees/day. Port of `AscDash`
/// (swehouse.c:2131-2147). Must be called with the exact same `(x, f)` pair used for the
/// corresponding `asc1` position call.
fn asc_dash(x: f64, f: f64, sine: f64, cose: f64) -> f64 {
    let cosx = cosd(x);
    let sinx = sind(x);
    let sinx2 = sinx * sinx;
    let c = cose * cosx - tand(f) * sine;
    let d = sinx2 + c * c;
    let dudt = if d > VERY_SMALL {
        (cosx * c + cose * sinx2) / d
    } else {
        0.0
    };
    dudt * ARMCS
}

/// Keeps the Ascendant on the eastern hemisphere near the poles. Port of `fix_asc_polar`
/// (swehouse.c:2169-2177). Used by `swe_house_pos`, ported in a later sub-task.
#[allow(dead_code)]
fn fix_asc_polar(asc: f64, armc: f64, eps: f64, geolat: f64) -> f64 {
    let demc = atand(sind(armc) * tand(eps));
    let mut asc = asc;
    if geolat >= 0.0 && 90.0 - geolat + demc < 0.0 {
        asc = normalize_degrees(asc + 180.0);
    }
    if geolat < 0.0 && -90.0 - geolat + demc > 0.0 {
        asc = normalize_degrees(asc + 180.0);
    }
    asc
}

/// The `tand(th)/cose` ecliptic projection shared by `CalcH`'s inline MC computation and the
/// equatorial-ascendant special point, both of which unconditionally `swe_degnorm` the final
/// result (unlike the standalone `swi_armc_to_mc`/`crate::math::armc_to_mc`, which only
/// normalizes inside the `+180` branch — see c-ref-houses.md §12.1, FP-fidelity hazard #1).
fn mc_like(th: f64, cose: f64) -> f64 {
    let mc = if (th - 90.0).abs() > VERY_SMALL && (th - 270.0).abs() > VERY_SMALL {
        let tant = tand(th);
        let mut mc = atand(tant / cose);
        if th > 90.0 && th <= 270.0 {
            mc = normalize_degrees(mc + 180.0);
        }
        mc
    } else if (th - 90.0).abs() <= VERY_SMALL {
        90.0
    } else {
        270.0
    };
    normalize_degrees(mc)
}

/// Porphyry (`'O'`) cusp fill — the universal polar-circle fallback target for the iterative
/// systems (Placidus/Koch/Gauquelin, later sub-tasks) as well as a house system in its own
/// right. Port of the `porphyry:` label body (swehouse.c:1310-1335). Re-asserts `cusps[1]`/
/// `cusps[10]` unconditionally, repairing any partial writes from a failed iterative attempt.
/// Returns the (possibly polar-swapped) ascendant.
fn fill_porphyry(
    cusps: &mut [f64; 37],
    cusp_speeds: &mut [f64; 37],
    mut ac: f64,
    mc: f64,
    ac_speed: f64,
    mc_speed: f64,
    do_speed: bool,
) -> f64 {
    let mut acmc = diff_degrees(ac, mc);
    if acmc < 0.0 {
        ac = normalize_degrees(ac + 180.0);
        cusps[1] = ac;
        acmc = diff_degrees(ac, mc);
    }
    cusps[1] = ac;
    cusps[10] = mc;
    cusps[2] = normalize_degrees(ac + (180.0 - acmc) / 3.0);
    cusps[3] = normalize_degrees(ac + (180.0 - acmc) / 3.0 * 2.0);
    cusps[11] = normalize_degrees(mc + acmc / 3.0);
    cusps[12] = normalize_degrees(mc + acmc / 3.0 * 2.0);
    if do_speed {
        let q1_speed = ac_speed - mc_speed;
        cusp_speeds[1] = ac_speed;
        cusp_speeds[10] = mc_speed;
        cusp_speeds[2] = ac_speed - q1_speed / 3.0;
        cusp_speeds[3] = ac_speed - q1_speed / 3.0 * 2.0;
        cusp_speeds[11] = ac_speed + q1_speed / 3.0;
        cusp_speeds[12] = ac_speed + q1_speed / 3.0 * 2.0;
    }
    ac
}

/// Polar-circle 180° shift shared by the quadrant-trisection systems (Campanus, Horizon,
/// Savard-A, Regiomontanus): when the (co-)latitude falls within the polar circle relative to
/// the ecliptic obliquity and the ascendant has landed on the wrong side of the meridian, flip
/// `ac`, `mc`, and the four cusps adjacent to them. Cusps 4-9 are filled later by the
/// post-switch opposite-cusp mirror, so they're intentionally excluded here. Port of the shared
/// polar-handling tail in the `'C'`/`'H'`/`'J'`/`'R'` switch cases (e.g. swehouse.c:1071-1081).
fn polar_shift_subset(cusps: &mut [f64; 37], ac: &mut f64, mc: &mut f64, lat: f64, eps: f64) {
    if lat.abs() >= 90.0 - eps && diff_degrees(*ac, *mc) < 0.0 {
        *ac = normalize_degrees(*ac + 180.0);
        *mc = normalize_degrees(*mc + 180.0);
        for i in [1usize, 2, 3, 10, 11, 12] {
            cusps[i] = normalize_degrees(cusps[i] + 180.0);
        }
    }
}

/// Outcome of one Placidus/Gauquelin Newton-iteration cusp solve. Both house systems share an
/// identical per-cusp iteration skeleton (swehouse.c:1623-1730 "G", 1830-1983 default) — only
/// the pole-height seed and fractional divisor differ.
enum NewtonCusp {
    /// Converged; `f` is the pole height AT the converged cusp, required for the analytical
    /// `AscDash` speed call (speed is evaluated at the converged point, not finite-differenced).
    Converged { cusp: f64, f: f64 },
    /// `|tant| < VERY_SMALL`: the cusp coincides with the AC/DC axis. Caller uses `rectasc` as
    /// the cusp and `ARMCS` as the speed.
    DegenerateAxis,
    /// Hit `NITER_MAX` without converging, OR converged/degenerated exactly on the `NITER_MAX`th
    /// iteration. C's post-loop check is `i >= niter_max` (swehouse.c:1667 et al.), which rejects
    /// the cap iteration's result even if it satisfied the convergence/degeneracy test on that
    /// exact step. Caller falls back to Porphyry for the whole system.
    NonConverged,
}

/// Shared Newton-iteration skeleton for Placidus (§5 "Default") and Gauquelin (§5 "G") cusps.
/// `rectasc`/`fh_init`/`divisor` are precomputed per-cusp by the caller (Placidus's fixed `3`/
/// `1.5` divisors, Gauquelin's `9/ih2`); `tane`/`geolat` feed only into that precomputation, not
/// the loop itself, so they aren't parameters here.
fn placidus_newton_cusp(
    rectasc: f64,
    fh_init: f64,
    divisor: f64,
    sine: f64,
    cose: f64,
    tanfi: f64,
) -> NewtonCusp {
    let seed = asc1(rectasc, fh_init, sine, cose);
    let mut tant = tand(asind(sine * sind(seed)));
    if tant.abs() < VERY_SMALL {
        return NewtonCusp::DegenerateAxis;
    }
    let mut f = atand(sind(asind(tanfi * tant) / divisor) / tant);
    let mut cusp = asc1(rectasc, f, sine, cose);
    let mut cuspsv = 0.0;
    for i in 1..=NITER_MAX {
        tant = tand(asind(sine * sind(cusp)));
        if tant.abs() < VERY_SMALL {
            if i >= NITER_MAX {
                break;
            }
            return NewtonCusp::DegenerateAxis;
        }
        f = atand(sind(asind(tanfi * tant) / divisor) / tant);
        cusp = asc1(rectasc, f, sine, cose);
        if i > 1 && diff_degrees(cusp, cuspsv).abs() < VERY_SMALL_PLAC_ITER {
            if i >= NITER_MAX {
                break;
            }
            return NewtonCusp::Converged { cusp, f };
        }
        cuspsv = cusp;
    }
    NewtonCusp::NonConverged
}

// ---------------------------------------------------------------------------
// CalcH — THE CORE (swehouse.c:892-2050)
// ---------------------------------------------------------------------------

struct CalcH {
    cusps: [f64; 37],
    cusp_speeds: [f64; 37],
    ascmc: AscMc,
    ascmc_speeds: AscMc,
    do_interpol: bool,
}

fn calc_h(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: HouseSystem,
    sundec: Option<f64>,
    do_speed: bool,
) -> Result<CalcH, Error> {
    // Consumed by the Sunshine ('I'/'i') branch, added in a later sub-task.
    let _ = sundec;

    let th = armc;
    let cose = cosd(eps);
    let sine = sind(eps);
    let tane = tand(eps);

    let mut geolat = geolat;
    if (geolat.abs() - 90.0).abs() < VERY_SMALL {
        geolat = if geolat < 0.0 {
            -90.0 + VERY_SMALL
        } else {
            90.0 - VERY_SMALL
        };
    }
    let tanfi = tand(geolat);

    let mut mc = mc_like(th, cose);
    let mc_speed = if do_speed {
        asc_dash(th, 0.0, sine, cose)
    } else {
        0.0
    };

    // The horizon's pole height equals geographic latitude; it crosses the equator 90° east
    // of the meridian.
    let mut ac = asc1(th + 90.0, geolat, sine, cose);
    let ac_speed = if do_speed {
        asc_dash(th + 90.0, geolat, sine, cose)
    } else {
        0.0
    };

    let armc_speed = ARMCS;

    let mut cusps = [0.0; 37];
    let mut cusp_speeds = [0.0; 37];
    cusps[1] = ac;
    cusps[10] = mc;
    if do_speed {
        cusp_speeds[1] = ac_speed;
        cusp_speeds[10] = mc_speed;
    }

    let mut do_interpol = false;
    // Set only by a fully-converged Gauquelin (36 independently-filled cusps) — the one system
    // excluded, along with 'Y'/'I'/'i' (not yet ported), from the post-switch opposite-cusp
    // mirror (swehouse.c:1985-2000, c-ref-houses.md §3 step 3).
    let mut skip_mirror = false;

    match hsys {
        HouseSystem::Equal => {
            // A / E — equal houses (swehouse.c:994-1010)
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
            }
            for i in 2..=12usize {
                cusps[i] = normalize_degrees(cusps[1] + (i as f64 - 1.0) * 30.0);
            }
            if do_speed {
                for cs in cusp_speeds.iter_mut().take(13).skip(1) {
                    *cs = ac_speed;
                }
            }
        }
        HouseSystem::EqualMC => {
            // D — equal, begin at MC (swehouse.c:1011-1027)
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
            cusps[10] = mc;
            for i in 11..=12usize {
                cusps[i] = normalize_degrees(cusps[10] + (i as f64 - 10.0) * 30.0);
            }
            for i in 1..=9usize {
                cusps[i] = normalize_degrees(cusps[10] + (i as f64 + 2.0) * 30.0);
            }
            if do_speed {
                for cs in cusp_speeds.iter_mut().take(13).skip(1) {
                    *cs = mc_speed;
                }
            }
        }
        HouseSystem::EqualAries => {
            // N — equal, begin at 0° Aries (whole-sign zodiac) (swehouse.c:1301-1309)
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
            for (i, cs) in cusps.iter_mut().enumerate().take(13).skip(1) {
                *cs = (i as f64 - 1.0) * 30.0;
            }
            // No cusp_speed handling — see c-ref-houses.md §4.2(e).
        }
        HouseSystem::Vehlow => {
            // V — equal houses after Vehlow (swehouse.c:1459-1473)
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
            cusps[1] = normalize_degrees(ac - 15.0);
            for i in 2..=12usize {
                cusps[i] = normalize_degrees(cusps[1] + (i as f64 - 1.0) * 30.0);
            }
            if do_speed {
                for cs in cusp_speeds.iter_mut().take(13).skip(1) {
                    *cs = ac_speed;
                }
            }
        }
        HouseSystem::WholeSign => {
            // W — whole sign (swehouse.c:1474-1484)
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
            }
            cusps[1] = ac - (ac % 30.0);
            for i in 2..=12usize {
                cusps[i] = normalize_degrees(cusps[1] + (i as f64 - 1.0) * 30.0);
            }
            // No cusp_speed handling — see c-ref-houses.md §4.2(e).
        }
        HouseSystem::Porphyry => {
            // O — Porphyry (swehouse.c:1310-1335, label `porphyry:`)
            ac = fill_porphyry(
                &mut cusps,
                &mut cusp_speeds,
                ac,
                mc,
                ac_speed,
                mc_speed,
                do_speed,
            );
        }
        HouseSystem::Sripati => {
            // S — Sripati (swehouse.c:1410-1431): Porphyry sector midpoints.
            let mut acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                acmc = diff_degrees(ac, mc);
            }
            let q1 = 180.0 - acmc;
            let s1 = q1 / 3.0;
            let s4 = acmc / 3.0;
            cusps[1] = normalize_degrees(ac - s4 * 0.5);
            cusps[2] = normalize_degrees(ac + s1 * 0.5);
            cusps[3] = normalize_degrees(ac + s1 * 1.5);
            cusps[10] = normalize_degrees(mc - s1 * 0.5);
            cusps[11] = normalize_degrees(mc + s4 * 0.5);
            cusps[12] = normalize_degrees(mc + s4 * 1.5);
            do_interpol = do_speed;
        }
        HouseSystem::Meridian => {
            // X — Meridian / axial rotation (swehouse.c:1485-1516)
            let mut a = th;
            for i in 1..=12usize {
                let mut j = i + 10;
                if j > 12 {
                    j -= 12;
                }
                a = normalize_degrees(a + 30.0);
                cusps[j] = mc_like(a, cose);
            }
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
            do_interpol = do_speed;
        }
        HouseSystem::Morinus => {
            // M — Morinus (swehouse.c:1517-1540): same equatorial points as X, projected via
            // a full cotrans (equatorial → ecliptic, +eps) instead of the tand/cose shortcut.
            let mut a = th;
            for i in 1..=12usize {
                let mut j = i + 10;
                if j > 12 {
                    j -= 12;
                }
                a = normalize_degrees(a + 30.0);
                let x = cotrans([a, 0.0, 1.0], eps);
                cusps[j] = x[0];
            }
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
            do_interpol = do_speed;
        }
        HouseSystem::Carter => {
            // F — Carter "poli-equatorial" (swehouse.c:1541-1580)
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
            }
            let x = cotrans([ac, 0.0, 1.0], -eps);
            let a = x[0];
            for i in [2usize, 3, 10, 11, 12] {
                let ra = normalize_degrees(a + (i as f64 - 1.0) * 30.0);
                cusps[i] = mc_like(ra, cose);
            }
            do_interpol = do_speed;
        }
        HouseSystem::Regiomontanus => {
            // R — Regiomontanus (swehouse.c:1381-1409)
            let fh1 = atand(tanfi * 0.5);
            let fh2 = atand(tanfi * cosd(30.0));
            let (x11, x12, x2, x3) = (30.0 + th, 60.0 + th, 120.0 + th, 150.0 + th);
            cusps[11] = asc1(x11, fh1, sine, cose);
            cusps[12] = asc1(x12, fh2, sine, cose);
            cusps[2] = asc1(x2, fh2, sine, cose);
            cusps[3] = asc1(x3, fh1, sine, cose);
            if do_speed {
                cusp_speeds[11] = asc_dash(x11, fh1, sine, cose);
                cusp_speeds[12] = asc_dash(x12, fh2, sine, cose);
                cusp_speeds[2] = asc_dash(x2, fh2, sine, cose);
                cusp_speeds[3] = asc_dash(x3, fh1, sine, cose);
            }
            polar_shift_subset(&mut cusps, &mut ac, &mut mc, geolat, eps);
        }
        HouseSystem::PolichPage => {
            // T — Polich/Page "topocentric" (swehouse.c:1432-1458). Structurally identical to
            // Regiomontanus, but with tanfi/3 and tanfi*2/3 pole heights.
            let fh1 = atand(tanfi / 3.0);
            let fh2 = atand(tanfi * 2.0 / 3.0);
            let (x11, x12, x2, x3) = (30.0 + th, 60.0 + th, 120.0 + th, 150.0 + th);
            cusps[11] = asc1(x11, fh1, sine, cose);
            cusps[12] = asc1(x12, fh2, sine, cose);
            cusps[2] = asc1(x2, fh2, sine, cose);
            cusps[3] = asc1(x3, fh1, sine, cose);
            if do_speed {
                cusp_speeds[11] = asc_dash(x11, fh1, sine, cose);
                cusp_speeds[12] = asc_dash(x12, fh2, sine, cose);
                cusp_speeds[2] = asc_dash(x2, fh2, sine, cose);
                cusp_speeds[3] = asc_dash(x3, fh1, sine, cose);
            }
            // Polar shift on ALL 12 cusps (not the {1,2,3,10,11,12}-only subset used by
            // C/H/J/R) — cusps 4-9 aren't meaningfully set yet, so this is harmless but
            // structurally different; replicate literally (swehouse.c §"T").
            if geolat.abs() >= 90.0 - eps && diff_degrees(ac, mc) < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                mc = normalize_degrees(mc + 180.0);
                for cs in cusps.iter_mut().take(13).skip(1) {
                    *cs = normalize_degrees(*cs + 180.0);
                }
            }
        }
        HouseSystem::Campanus => {
            // C — Campanus (swehouse.c:1028-1082)
            let fh1 = asind(sind(geolat) / 2.0);
            let fh2 = asind(3.0_f64.sqrt() / 2.0 * sind(geolat));
            let cosfi = cosd(geolat);
            let (xh1, xh2) = if cosfi == 0.0 {
                let v = if geolat > 0.0 { 90.0 } else { 270.0 };
                (v, v)
            } else {
                (
                    atand(3.0_f64.sqrt() / cosfi),
                    atand((1.0 / 3.0_f64.sqrt()) / cosfi),
                )
            };
            cusps[11] = asc1(th + 90.0 - xh1, fh1, sine, cose);
            cusps[12] = asc1(th + 90.0 - xh2, fh2, sine, cose);
            cusps[2] = asc1(th + 90.0 + xh2, fh2, sine, cose);
            cusps[3] = asc1(th + 90.0 + xh1, fh1, sine, cose);
            if do_speed {
                cusp_speeds[11] = asc_dash(th + 90.0 - xh1, fh1, sine, cose);
                cusp_speeds[12] = asc_dash(th + 90.0 - xh2, fh2, sine, cose);
                cusp_speeds[2] = asc_dash(th + 90.0 + xh2, fh2, sine, cose);
                cusp_speeds[3] = asc_dash(th + 90.0 + xh1, fh1, sine, cose);
            }
            polar_shift_subset(&mut cusps, &mut ac, &mut mc, geolat, eps);
        }
        HouseSystem::Horizon => {
            // H — Horizon/Azimuth (swehouse.c:1083-1155): Campanus-style trisection of the
            // prime vertical, rotated 180° in th and with fi mapped to its co-latitude first.
            let mut fi2 = if geolat > 0.0 {
                90.0 - geolat
            } else {
                -90.0 - geolat
            };
            if (fi2.abs() - 90.0).abs() < VERY_SMALL {
                fi2 = if fi2 < 0.0 {
                    -90.0 + VERY_SMALL
                } else {
                    90.0 - VERY_SMALL
                };
            }
            let th2 = normalize_degrees(th + 180.0);
            let cosfi2 = cosd(fi2);
            let fh1 = asind(sind(fi2) / 2.0);
            let fh2 = asind(3.0_f64.sqrt() / 2.0 * sind(fi2));
            let (xh1, xh2) = if cosfi2 == 0.0 {
                let v = if fi2 > 0.0 { 90.0 } else { 270.0 };
                (v, v)
            } else {
                (
                    atand(3.0_f64.sqrt() / cosfi2),
                    atand((1.0 / 3.0_f64.sqrt()) / cosfi2),
                )
            };
            cusps[11] = asc1(th2 + 90.0 - xh1, fh1, sine, cose);
            cusps[12] = asc1(th2 + 90.0 - xh2, fh2, sine, cose);
            cusps[1] = asc1(th2 + 90.0, fi2, sine, cose);
            cusps[2] = asc1(th2 + 90.0 + xh2, fh2, sine, cose);
            cusps[3] = asc1(th2 + 90.0 + xh1, fh1, sine, cose);
            if do_speed {
                cusp_speeds[11] = asc_dash(th2 + 90.0 - xh1, fh1, sine, cose);
                cusp_speeds[12] = asc_dash(th2 + 90.0 - xh2, fh2, sine, cose);
                cusp_speeds[1] = asc_dash(th2 + 90.0, fi2, sine, cose);
                cusp_speeds[2] = asc_dash(th2 + 90.0 + xh2, fh2, sine, cose);
                cusp_speeds[3] = asc_dash(th2 + 90.0 + xh1, fh1, sine, cose);
            }
            // Polar-circle shift exactly as Campanus, evaluated against the co-latitude fi2.
            polar_shift_subset(&mut cusps, &mut ac, &mut mc, fi2, eps);
            // Unconditional re-orientation into ecliptic-house ordering (swehouse.c:1141-1144).
            cusps[1] = normalize_degrees(cusps[1] + 180.0);
            cusps[2] = normalize_degrees(cusps[2] + 180.0);
            cusps[3] = normalize_degrees(cusps[3] + 180.0);
            cusps[11] = normalize_degrees(cusps[11] + 180.0);
            cusps[12] = normalize_degrees(cusps[12] + 180.0);
            // Final AC/DC sanity check (without an MC shift this time).
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
        }
        HouseSystem::SavardA => {
            // J — Savard-A (swehouse.c:1182-1249)
            let sinfi = sind(geolat);
            let cosfi = cosd(geolat);
            let (xs1, xs2) = if geolat.abs() < VERY_SMALL {
                (asind(2.0 / 3.0), asind(1.0 / 3.0))
            } else {
                (
                    asind(sind(2.0 * geolat / 3.0) / sinfi),
                    asind(sind(geolat / 3.0) / sinfi),
                )
            };
            let (xh1, xh2) = if cosfi == 0.0 {
                let v = if geolat > 0.0 { 90.0 } else { 270.0 };
                (v, v)
            } else {
                (atand(tand(xs1) / cosfi), atand(tand(xs2) / cosfi))
            };
            let fh1 = asind(sind(geolat) * sind(90.0 - xs1));
            let fh2 = asind(sind(geolat) * sind(90.0 - xs2));
            cusps[11] = asc1(th + 90.0 - xh1, fh1, sine, cose);
            cusps[12] = asc1(th + 90.0 - xh2, fh2, sine, cose);
            cusps[2] = asc1(th + 90.0 + xh2, fh2, sine, cose);
            cusps[3] = asc1(th + 90.0 + xh1, fh1, sine, cose);
            if do_speed {
                cusp_speeds[11] = asc_dash(th + 90.0 - xh1, fh1, sine, cose);
                cusp_speeds[12] = asc_dash(th + 90.0 - xh2, fh2, sine, cose);
                cusp_speeds[2] = asc_dash(th + 90.0 + xh2, fh2, sine, cose);
                cusp_speeds[3] = asc_dash(th + 90.0 + xh1, fh1, sine, cose);
            }
            polar_shift_subset(&mut cusps, &mut ac, &mut mc, geolat, eps);
        }
        HouseSystem::Koch => {
            // K — Koch (swehouse.c:1250-1272): closed-form, no iteration. Fails outright (no
            // Newton attempt) in the polar circle, unlike the great-circle quadrant systems.
            if geolat.abs() >= 90.0 - eps {
                ac = fill_porphyry(
                    &mut cusps,
                    &mut cusp_speeds,
                    ac,
                    mc,
                    ac_speed,
                    mc_speed,
                    do_speed,
                );
            } else {
                let sina = (sind(mc) * sine / cosd(geolat)).clamp(-1.0, 1.0);
                let cosa = (1.0 - sina * sina).sqrt();
                let c = atand(tanfi / cosa);
                let ad3 = asind(sind(c) * sina) / 3.0;
                let x11 = th + 30.0 - 2.0 * ad3;
                let x12 = th + 60.0 - ad3;
                let x2 = th + 120.0 + ad3;
                let x3 = th + 150.0 + 2.0 * ad3;
                cusps[11] = asc1(x11, geolat, sine, cose);
                cusps[12] = asc1(x12, geolat, sine, cose);
                cusps[2] = asc1(x2, geolat, sine, cose);
                cusps[3] = asc1(x3, geolat, sine, cose);
                if do_speed {
                    cusp_speeds[11] = asc_dash(x11, geolat, sine, cose);
                    cusp_speeds[12] = asc_dash(x12, geolat, sine, cose);
                    cusp_speeds[2] = asc_dash(x2, geolat, sine, cose);
                    cusp_speeds[3] = asc_dash(x3, geolat, sine, cose);
                }
            }
        }
        HouseSystem::Placidus => {
            // Default — Placidus (swehouse.c:1830-1983): four independent Newton loops (cusps
            // 11, 12, 2, 3), each the same skeleton as a single Gauquelin sector but with fixed
            // fractional divisors instead of `ih2/9`.
            if geolat.abs() >= 90.0 - eps {
                ac = fill_porphyry(
                    &mut cusps,
                    &mut cusp_speeds,
                    ac,
                    mc,
                    ac_speed,
                    mc_speed,
                    do_speed,
                );
            } else {
                let a = asind(tanfi * tane);
                let fh1 = atand(sind(a / 3.0) / tane);
                let fh2 = atand(sind(a * 2.0 / 3.0) / tane);
                let specs = [
                    (11usize, fh1, normalize_degrees(30.0 + th), 3.0),
                    (12usize, fh2, normalize_degrees(60.0 + th), 1.5),
                    (2usize, fh2, normalize_degrees(120.0 + th), 1.5),
                    (3usize, fh1, normalize_degrees(150.0 + th), 3.0),
                ];
                let mut fell_back = false;
                for (idx, fh_init, rectasc, divisor) in specs {
                    match placidus_newton_cusp(rectasc, fh_init, divisor, sine, cose, tanfi) {
                        NewtonCusp::Converged { cusp, f } => {
                            cusps[idx] = cusp;
                            if do_speed {
                                cusp_speeds[idx] = asc_dash(rectasc, f, sine, cose);
                            }
                        }
                        NewtonCusp::DegenerateAxis => {
                            cusps[idx] = rectasc;
                            if do_speed {
                                cusp_speeds[idx] = ARMCS;
                            }
                        }
                        NewtonCusp::NonConverged => {
                            fell_back = true;
                            break;
                        }
                    }
                }
                if fell_back {
                    ac = fill_porphyry(
                        &mut cusps,
                        &mut cusp_speeds,
                        ac,
                        mc,
                        ac_speed,
                        mc_speed,
                        do_speed,
                    );
                }
            }
        }
        HouseSystem::Gauquelin => {
            // G — 36 Gauquelin sectors (swehouse.c:1623-1730): two mirrored Newton-iteration
            // loops (4th/2nd quarter, then 1st/3rd quarter), each filling 8 sectors plus their
            // 180°-opposite partners. Counted clockwise. Excluded from the post-switch mirror —
            // it fills all 36 cusps itself.
            if geolat.abs() >= 90.0 - eps {
                ac = fill_porphyry(
                    &mut cusps,
                    &mut cusp_speeds,
                    ac,
                    mc,
                    ac_speed,
                    mc_speed,
                    do_speed,
                );
            } else {
                let a = asind(tanfi * tane);
                let mut fell_back = false;

                // 4th/2nd quarter: ih = 2..9, ih2 = 10-ih.
                for ih in 2..=9usize {
                    let ih2 = (10 - ih) as f64;
                    let fh_init = atand(sind(a * ih2 / 9.0) / tane);
                    let rectasc = normalize_degrees(90.0 / 9.0 * ih2 + th);
                    let divisor = 9.0 / ih2;
                    match placidus_newton_cusp(rectasc, fh_init, divisor, sine, cose, tanfi) {
                        NewtonCusp::Converged { cusp, f } => {
                            cusps[ih] = cusp;
                            cusps[ih + 18] = normalize_degrees(cusp + 180.0);
                            if do_speed {
                                let sp = asc_dash(rectasc, f, sine, cose);
                                cusp_speeds[ih] = sp;
                                cusp_speeds[ih + 18] = sp;
                            }
                        }
                        NewtonCusp::DegenerateAxis => {
                            cusps[ih] = rectasc;
                            cusps[ih + 18] = normalize_degrees(rectasc + 180.0);
                            if do_speed {
                                cusp_speeds[ih] = ARMCS;
                                cusp_speeds[ih + 18] = ARMCS;
                            }
                        }
                        NewtonCusp::NonConverged => {
                            fell_back = true;
                            break;
                        }
                    }
                }

                // 1st/3rd quarter: ih = 29..36, ih2 = ih-28 — mirror-image formulas.
                if !fell_back {
                    for ih in 29..=36usize {
                        let ih2 = (ih - 28) as f64;
                        let fh_init = atand(sind(a * ih2 / 9.0) / tane);
                        let rectasc = normalize_degrees(180.0 - ih2 * 90.0 / 9.0 + th);
                        let divisor = 9.0 / ih2;
                        match placidus_newton_cusp(rectasc, fh_init, divisor, sine, cose, tanfi) {
                            NewtonCusp::Converged { cusp, f } => {
                                cusps[ih] = cusp;
                                cusps[ih - 18] = normalize_degrees(cusp + 180.0);
                                if do_speed {
                                    let sp = asc_dash(rectasc, f, sine, cose);
                                    cusp_speeds[ih] = sp;
                                    cusp_speeds[ih - 18] = sp;
                                }
                            }
                            NewtonCusp::DegenerateAxis => {
                                cusps[ih] = rectasc;
                                cusps[ih - 18] = normalize_degrees(rectasc + 180.0);
                                if do_speed {
                                    cusp_speeds[ih] = ARMCS;
                                    cusp_speeds[ih - 18] = ARMCS;
                                }
                            }
                            NewtonCusp::NonConverged => {
                                fell_back = true;
                                break;
                            }
                        }
                    }
                }

                if fell_back {
                    ac = fill_porphyry(
                        &mut cusps,
                        &mut cusp_speeds,
                        ac,
                        mc,
                        ac_speed,
                        mc_speed,
                        do_speed,
                    );
                } else {
                    cusps[1] = ac;
                    cusps[10] = mc;
                    cusps[19] = normalize_degrees(ac + 180.0);
                    cusps[28] = normalize_degrees(mc + 180.0);
                    if do_speed {
                        cusp_speeds[1] = ac_speed;
                        cusp_speeds[10] = mc_speed;
                        cusp_speeds[19] = ac_speed;
                        cusp_speeds[28] = mc_speed;
                    }
                    skip_mirror = true;
                }
            }
        }
        _ => {
            return Err(Error::CError(format!(
                "house system {hsys:?} not yet implemented"
            )));
        }
    }

    // Post-switch opposite-cusp mirror (swehouse.c:1985-2000) — skipped for G (fills all 36
    // cusps itself) on a fully-converged path, and (not yet reachable) Y/I/i.
    if !skip_mirror {
        cusps[4] = normalize_degrees(cusps[10] + 180.0);
        cusps[5] = normalize_degrees(cusps[11] + 180.0);
        cusps[6] = normalize_degrees(cusps[12] + 180.0);
        cusps[7] = normalize_degrees(cusps[1] + 180.0);
        cusps[8] = normalize_degrees(cusps[2] + 180.0);
        cusps[9] = normalize_degrees(cusps[3] + 180.0);
        if do_speed && !do_interpol {
            cusp_speeds[4] = cusp_speeds[10];
            cusp_speeds[5] = cusp_speeds[11];
            cusp_speeds[6] = cusp_speeds[12];
            cusp_speeds[7] = cusp_speeds[1];
            cusp_speeds[8] = cusp_speeds[2];
            cusp_speeds[9] = cusp_speeds[3];
        }
    }

    // Special points (swehouse.c:2001-2049), always computed.
    let f_vertex = if geolat >= 0.0 {
        90.0 - geolat
    } else {
        -90.0 - geolat
    };
    let mut vertex = asc1(th - 90.0, f_vertex, sine, cose);
    let vertex_speed = if do_speed {
        asc_dash(th - 90.0, f_vertex, sine, cose)
    } else {
        0.0
    };
    // With tropical latitudes the vertex behaves like the ascendant within the polar
    // circle; keep it on the western hemisphere.
    if geolat.abs() <= eps {
        let vemc = diff_degrees(vertex, mc);
        if vemc > 0.0 {
            vertex = normalize_degrees(vertex + 180.0);
        }
    }

    let equasc = mc_like(normalize_degrees(th + 90.0), cose);
    let equasc_speed = if do_speed {
        asc_dash(th + 90.0, 0.0, sine, cose)
    } else {
        0.0
    };

    let coasc1 = normalize_degrees(asc1(th - 90.0, geolat, sine, cose) + 180.0);
    let coasc1_speed = if do_speed {
        asc_dash(th - 90.0, geolat, sine, cose)
    } else {
        0.0
    };

    let (coasc2, coasc2_speed) = if geolat >= 0.0 {
        let f = 90.0 - geolat;
        (
            asc1(th + 90.0, f, sine, cose),
            if do_speed {
                asc_dash(th + 90.0, f, sine, cose)
            } else {
                0.0
            },
        )
    } else {
        let f = -90.0 - geolat;
        (
            asc1(th + 90.0, f, sine, cose),
            if do_speed {
                asc_dash(th + 90.0, f, sine, cose)
            } else {
                0.0
            },
        )
    };

    let polasc = asc1(th - 90.0, geolat, sine, cose);
    let polasc_speed = if do_speed {
        asc_dash(th - 90.0, geolat, sine, cose)
    } else {
        0.0
    };

    let ascmc = AscMc {
        ascendant: ac,
        mc,
        armc: th,
        vertex,
        equatorial_ascendant: equasc,
        coascendant_koch: coasc1,
        coascendant_munkasey: coasc2,
        polar_ascendant: polasc,
    };
    let ascmc_speeds = AscMc {
        ascendant: ac_speed,
        mc: mc_speed,
        armc: armc_speed,
        vertex: vertex_speed,
        equatorial_ascendant: equasc_speed,
        coascendant_koch: coasc1_speed,
        coascendant_munkasey: coasc2_speed,
        polar_ascendant: polasc_speed,
    };

    Ok(CalcH {
        cusps,
        cusp_speeds,
        ascmc,
        ascmc_speeds,
        do_interpol,
    })
}

// ---------------------------------------------------------------------------
// Driver — swe_houses_armc_ex2 (swehouse.c:622-774)
// ---------------------------------------------------------------------------

pub fn houses_armc(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: HouseSystem,
    sundec: Option<f64>,
) -> Result<HouseResult, Error> {
    let armc = normalize_degrees(armc);
    let h = calc_h(armc, geolat, eps, hsys, sundec, true)?;

    let mut cusp_speeds = h.cusp_speeds;

    if h.do_interpol {
        let dt = 1.0 / 86400.0;
        let darmc = dt * ARMCS;
        let hm = calc_h(armc - darmc, geolat, eps, hsys, sundec, false);
        let hp = calc_h(armc + darmc, geolat, eps, hsys, sundec, false);

        // Matches swe_houses_armc_ex2 (swehouse.c:704-716): if either side probe
        // fails to converge, keep the already-computed main cusp_speeds instead
        // of propagating the error.
        if let (Ok(hm), Ok(hp)) = (hm, hp) {
            let mut dt = dt;
            let mut hm_cusps = hm.cusps;
            let mut hp_cusps = hp.cusps;
            if diff_degrees(hp.ascmc.ascendant, h.ascmc.ascendant).abs() > 90.0 {
                hp_cusps = h.cusps;
                dt /= 2.0;
            } else if diff_degrees(hm.ascmc.ascendant, h.ascmc.ascendant).abs() > 90.0 {
                hm_cusps = h.cusps;
                dt /= 2.0;
            }
            for i in 1..=12usize {
                cusp_speeds[i] = diff_degrees(hp_cusps[i], hm_cusps[i]) / 2.0 / dt;
            }
        }
    }

    Ok(HouseResult {
        cusps: h.cusps,
        cusp_speeds,
        ascmc: h.ascmc,
        ascmc_speeds: h.ascmc_speeds,
    })
}
