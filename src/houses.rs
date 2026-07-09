// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

//! Port of `swehouse.c`: house systems (`CalcH`) and the ARMC-based driver
//! (`swe_houses_armc_ex2`). See `docs/c-ref-houses.md`.

use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{
    cartesian_to_polar, cartesian_to_polar_with_speed, cotrans, cross_prod, diff_degrees,
    dot_prod_unit, normalize_degrees, normalize_radians, polar_to_cartesian,
    polar_to_cartesian_with_speed, rotate_x,
};
use crate::types::{AstroModels, HouseSystem, PrecessionDirection};

// ---------------------------------------------------------------------------
// Constants (swehouse.h:87, swehouse.c:68-70, swehouse.c:940)
// ---------------------------------------------------------------------------

const VERY_SMALL: f64 = 1e-10;
const VERY_SMALL_PLAC_ITER: f64 = 1.0 / 360000.0;
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

fn acosd(x: f64) -> f64 {
    x.acos() * crate::constants::RADTODEG
}

// ---------------------------------------------------------------------------
// Public output types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Ascendant, midheaven, and related angular house-system points (degrees).
pub struct AscMc {
    /// Ascendant (rising degree of the ecliptic on the eastern horizon).
    pub ascendant: f64,
    /// Midheaven (Medium Coeli): ecliptic point on the meridian.
    pub mc: f64,
    /// Sidereal time expressed as the right ascension of the meridian (ARMC).
    pub armc: f64,
    /// Vertex: the point where the prime vertical intersects the ecliptic.
    pub vertex: f64,
    /// Equatorial ascendant (East Point).
    pub equatorial_ascendant: f64,
    /// Co-ascendant per Walter Koch's method.
    pub coascendant_koch: f64,
    /// Co-ascendant per Michael Munkasey's method.
    pub coascendant_munkasey: f64,
    /// Polar ascendant (used at extreme geographic latitudes).
    pub polar_ascendant: f64,
}

impl AscMc {
    /// Returns the eight angular points as a flat array, in the same order as
    /// C's `ascmc[]` output.
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

/// Full house-system computation result: cusps, angular points, and their speeds.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HouseResult {
    /// House cusp longitudes (degrees), indices 1..=36 populated per house
    /// system (index 0 unused, matching C's `cusp[37]`).
    #[cfg_attr(feature = "serde", serde(with = "serde_37_array"))]
    pub cusps: [f64; 37],
    /// House cusp speeds (degrees/day), same indexing as `cusps`.
    #[cfg_attr(feature = "serde", serde(with = "serde_37_array"))]
    pub cusp_speeds: [f64; 37],
    /// Ascendant, MC, and related angular points.
    pub ascmc: AscMc,
    /// Speeds (degrees/day) of the angular points in `ascmc`.
    pub ascmc_speeds: AscMc,
}

/// `serde`'s built-in array impls only cover lengths 0..=32 (`serde_core::ser::impls`), so
/// `[f64; 37]` (cusps 1..=36 plus the unused index-0 slot, matching C's `cusp[37]`) needs a
/// manual tuple-of-37 encoding here.
#[cfg(feature = "serde")]
mod serde_37_array {
    use serde::de::{Error as DeError, SeqAccess, Visitor};
    use serde::ser::SerializeTuple;
    use serde::{Deserializer, Serializer};
    use std::fmt;

    /// Serializes a `[f64; 37]` as a fixed-size tuple (serde's built-in array
    /// impls only cover lengths up to 32).
    pub fn serialize<S>(arr: &[f64; 37], serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut tup = serializer.serialize_tuple(37)?;
        for v in arr {
            tup.serialize_element(v)?;
        }
        tup.end()
    }

    struct ArrayVisitor;

    impl<'de> Visitor<'de> for ArrayVisitor {
        type Value = [f64; 37];

        fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "an array of 37 f64 values")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let mut arr = [0.0f64; 37];
            for (i, slot) in arr.iter_mut().enumerate() {
                *slot = seq
                    .next_element()?
                    .ok_or_else(|| DeError::invalid_length(i, &self))?;
            }
            Ok(arr)
        }
    }

    /// Deserializes a `[f64; 37]` from the fixed-size tuple encoding produced
    /// by [`serialize`].
    pub fn deserialize<'de, D>(deserializer: D) -> Result<[f64; 37], D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_tuple(37, ArrayVisitor)
    }
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

/// `apc_sector(n, ph, e, az)` — radians in, degrees out. Port of the static `apc_sector` helper
/// (swehouse.c:782-825) used by the `'Y'` (APC) house system. Works in radians throughout via
/// plain `tan`/`atan`/`atan2`/`sin`/`cos`, unlike the degree-macro style used everywhere else in
/// this module — see c-ref-houses.md §12.3. `ph`=geolat, `e`=obliquity, `az`=armc.
fn apc_sector(n: i32, ph: f64, e: f64, az: f64) -> f64 {
    let (kv, dasc) = if (ph * crate::constants::RADTODEG).abs() > 90.0 - VERY_SMALL {
        (0.0, 0.0)
    } else {
        let kv = (ph.tan() * e.tan() * az.cos() / (1.0 + ph.tan() * e.tan() * az.sin())).atan();
        let dasc = if (ph * crate::constants::RADTODEG).abs() < VERY_SMALL {
            let d = (90.0 - VERY_SMALL) * crate::constants::DEGTORAD;
            if ph < 0.0 { -d } else { d }
        } else {
            (kv.sin() / ph.tan()).atan()
        };
        (kv, dasc)
    };
    let is_below_hor = n < 8;
    let k = if is_below_hor {
        (n - 1) as f64
    } else {
        (n - 13) as f64
    };
    let a = if is_below_hor {
        kv + az + std::f64::consts::FRAC_PI_2 + k * (std::f64::consts::FRAC_PI_2 - kv) / 3.0
    } else {
        kv + az + std::f64::consts::FRAC_PI_2 + k * (std::f64::consts::FRAC_PI_2 + kv) / 3.0
    };
    let a = normalize_radians(a);
    let dret = (dasc.tan() * ph.tan() * az.sin() + a.sin()).atan2(
        e.cos() * (dasc.tan() * ph.tan() * az.cos() + a.cos())
            + e.sin() * ph.tan() * (az - a).sin(),
    );
    normalize_degrees(dret * crate::constants::RADTODEG)
}

/// Keeps the Ascendant on the eastern hemisphere near the poles. Port of `fix_asc_polar`
/// (swehouse.c:2169-2177).
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

/// Ascensional-difference setup shared by both Sunshine solutions. Port of `sunshine_init`
/// (swehouse.c:2878-2904). `xh[1,4,7,10]` (cardinal cusps) are left at `0.0` — those are filled
/// separately by the I/i dispatcher. Returns `(xh, ok)`; `ok=false` means the Sun is exactly
/// circumpolar at this lat/dec (`|arg|>=1`).
fn sunshine_init(lat: f64, dec: f64) -> ([f64; 13], bool) {
    let arg = tand(dec) * tand(lat);
    let ad = if arg >= 1.0 {
        90.0 - VERY_SMALL
    } else if arg <= -1.0 {
        -90.0 + VERY_SMALL
    } else {
        asind(arg)
    };
    let nsa = 90.0 - ad;
    let dsa = 90.0 + ad;
    let mut xh = [0.0; 13];
    xh[2] = -2.0 * nsa / 3.0;
    xh[3] = -nsa / 3.0;
    xh[5] = nsa / 3.0;
    xh[6] = 2.0 * nsa / 3.0;
    xh[8] = -2.0 * dsa / 3.0;
    xh[9] = -dsa / 3.0;
    xh[11] = dsa / 3.0;
    xh[12] = 2.0 * dsa / 3.0;
    (xh, arg.abs() < 1.0)
}

/// Sunshine houses, Makransky solution (`'i'`). Port of `sunshine_solution_makransky`
/// (swehouse.c:2906-3046) — read directly from C, not from the ref doc abstraction, per the
/// sub-task's escape hatch: this is the most structurally intricate function in the file. The
/// 4-to-8-way case split on the quadrant of `w` and the sign of `z-90` is ported verbatim,
/// including the C author's own uncertainty about the `z>90` remap (their comment is preserved
/// below). Returns `false` (ERR) if the Sun is exactly circumpolar at this lat/dec.
fn sunshine_solution_makransky(
    th: f64,
    fi: f64,
    ekl: f64,
    dec: f64,
    cusps: &mut [f64; 37],
) -> bool {
    let (xh, ok) = sunshine_init(fi, dec);
    if !ok {
        return false;
    }
    let sinlat = sind(fi);
    let coslat = cosd(fi);
    let tanlat = tand(fi);
    let tandec = tand(dec);
    let sinecl = sind(ekl);

    for ih in 1..=12usize {
        if (ih - 1) % 3 == 0 {
            continue; // skip 1, 4, 7, 10
        }
        let md = xh[ih].abs();
        let mut rah = if ih <= 6 {
            normalize_degrees(th + 180.0 + xh[ih])
        } else {
            normalize_degrees(th + xh[ih])
        };
        if fi < 0.0 {
            // Makransky deals with southern latitude this way.
            rah = normalize_degrees(180.0 + rah);
        }
        let zd = if md == 90.0 {
            90.0 - atand(sinlat * tandec)
        } else {
            let a = if md < 90.0 {
                atand(coslat * tand(md))
            } else {
                atand(tand(md - 90.0) / coslat)
            };
            let b = atand(tanlat * cosd(md));
            let c = if ih <= 6 { b + dec } else { b - dec };
            let f = atand(sinlat * sind(md) * tand(c));
            a + f
        };
        let pole = asind(sind(zd) * sinlat);
        let q = asind(tandec * tand(pole));
        let in_dc_quadrant = ih <= 3 || ih >= 11;
        let w = if in_dc_quadrant {
            normalize_degrees(rah - q)
        } else {
            normalize_degrees(rah + q)
        };

        let cu = if w == 90.0 {
            let r = atand(sind(ekl) * tand(pole));
            if in_dc_quadrant { 90.0 + r } else { 90.0 - r }
        } else if w == 270.0 {
            let r = atand(sinecl * tand(pole));
            if in_dc_quadrant { 270.0 - r } else { 270.0 + r }
        } else {
            let m = atand((tand(pole) / cosd(w)).abs());
            let z = if in_dc_quadrant {
                if w > 90.0 && w < 270.0 {
                    m - ekl
                } else {
                    m + ekl
                }
            } else if w > 90.0 && w < 270.0 {
                m + ekl
            } else {
                m - ekl
            };
            let mut r = 0.0;
            let mut cu = if z == 90.0 {
                if w < 180.0 { 90.0 } else { 270.0 }
            } else {
                r = atand((cosd(m) * tand(w) / cosd(z)).abs());
                if w < 90.0 {
                    r
                } else if w > 90.0 && w < 180.0 {
                    180.0 - r
                } else if w > 180.0 && w < 270.0 {
                    180.0 + r
                } else {
                    360.0 - r
                }
            };
            if z > 90.0 {
                // "I am not sure if I understood the remark 'value will fall away from
                // cancer..' on page 146 correctly." — C author's comment, swehouse.c:3037.
                // Replicated verbatim, not "fixed".
                cu = if w < 90.0 {
                    180.0 - r
                } else if w > 90.0 && w < 180.0 {
                    r
                } else if w > 180.0 && w < 270.0 {
                    360.0 - r
                } else {
                    180.0 + r
                };
            }
            if fi < 0.0 {
                // Makransky deals with southern latitude this way. Note: unlike the rah
                // adjustment above, this only applies in the general (w != 90, 270) branch.
                cu = normalize_degrees(cu + 180.0);
            }
            cu
        };
        cusps[ih] = cu;
    }
    true
}

/// Sunshine houses, Treindl solution (`'I'`). Port of `sunshine_solution_treindl`
/// (swehouse.c:3048-3143). `SUNSHINE_KEEP_MC_SOUTH` is `#define`d to `0` in C (a compile-time
/// switch always built with the `0` branch) — only that behavior (MC kept north) is ported; the
/// dead `1` branch (negating `xh[2..12]`) is omitted. `sunshine_init`'s `ok=false` return is
/// ignored here (unlike Makransky, which short-circuits on it) — Treindl proceeds with the
/// clamped `±(90-VERY_SMALL)` ascensional difference even when the Sun is exactly circumpolar.
/// Returns `false` if any house hit the `c<1e-6` degeneracy; cusps are still filled in that case
/// (the C loop does not early-return) but the caller discards them and falls back to Porphyry,
/// matching the shared `retc==ERR` check in the `'I'`/`'i'` dispatcher (swehouse.c:1176-1180).
fn sunshine_solution_treindl(th: f64, fi: f64, ekl: f64, dec: f64, cusps: &mut [f64; 37]) -> bool {
    let (xh, _) = sunshine_init(fi, dec);
    let sinlat = sind(fi);
    let coslat = cosd(fi);
    let cosdec = cosd(dec);
    let tandec = tand(dec);
    let sinecl = sind(ekl);
    let cosecl = cosd(ekl);

    let mcdec = atand(sind(th) * tand(ekl));
    let mc_under_horizon = (fi - mcdec).abs() > 90.0;

    let mut ok = true;
    for ih in 1..=12usize {
        if (ih - 1) % 3 == 0 {
            continue; // skip 1, 4, 7, 10
        }
        let xhs = 2.0 * asind(cosdec * sind(xh[ih] / 2.0));
        let cosa = tandec * tand(xhs / 2.0);
        let alph = acosd(cosa);
        let (alpha2, b) = if ih > 7 {
            (180.0 - alph, 90.0 - fi + dec) // diurnal side
        } else {
            (alph, 90.0 - fi - dec) // nocturnal side
        };
        let cosc = cosd(xhs) * cosd(b) + sind(xhs) * sind(b) * cosd(alpha2);
        let c = acosd(cosc);
        if c < 1e-6 {
            ok = false;
        }
        let sinzd = sind(xhs) * sind(alpha2) / sind(c);
        let zd = asind(sinzd);
        let rax = atand(coslat * tand(zd));
        let mut pole = asind(sinzd * sinlat);
        let a = if ih <= 6 {
            pole = -pole;
            normalize_degrees(rax + th + 180.0)
        } else {
            normalize_degrees(th + rax)
        };
        cusps[ih] = asc1(a, pole, sinecl, cosecl);
    }
    // `mc_under_horizon && !SUNSHINE_KEEP_MC_SOUTH` simplifies to `mc_under_horizon`, since the
    // compile-time constant is always 0 in C.
    if mc_under_horizon {
        for (ih, cs) in cusps.iter_mut().enumerate().take(13).skip(1) {
            if (ih - 1) % 3 == 0 {
                continue;
            }
            *cs = normalize_degrees(*cs + 180.0);
        }
    }
    ok
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
        HouseSystem::Alcabitius => {
            // B — Alcabitius semiarc (swehouse.c:1581-1622)
            if diff_degrees(ac, mc) < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
            }
            let dek = asind(sind(ac) * sine);
            let r = (-tanfi * tand(dek)).clamp(-1.0, 1.0);
            let sda = acosd(r);
            let sna = 180.0 - sda;
            let sd3 = sda / 3.0;
            let sn3 = sna / 3.0;
            cusps[11] = asc1(normalize_degrees(th + sd3), 0.0, sine, cose);
            cusps[12] = asc1(normalize_degrees(th + 2.0 * sd3), 0.0, sine, cose);
            cusps[2] = asc1(normalize_degrees(th + 180.0 - 2.0 * sn3), 0.0, sine, cose);
            cusps[3] = asc1(normalize_degrees(th + 180.0 - sn3), 0.0, sine, cose);
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
        HouseSystem::PullenSD => {
            // L — Pullen SD "sinusoidal delta", ex Neo-Porphyry (swehouse.c:1273-1300)
            let mut acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
                acmc = diff_degrees(ac, mc);
            }
            let q1 = 180.0 - acmc;
            let mut d = (acmc - 90.0) / 4.0;
            if acmc <= 30.0 {
                cusps[11] = normalize_degrees(mc + acmc / 2.0);
                cusps[12] = cusps[11];
            } else {
                cusps[11] = normalize_degrees(mc + 30.0 + d);
                cusps[12] = normalize_degrees(mc + 60.0 + 3.0 * d);
            }
            d = (q1 - 90.0) / 4.0;
            if q1 <= 30.0 {
                cusps[2] = normalize_degrees(ac + q1 / 2.0);
                cusps[3] = cusps[2];
            } else {
                cusps[2] = normalize_degrees(ac + 30.0 + d);
                cusps[3] = normalize_degrees(ac + 60.0 + 3.0 * d);
            }
            do_interpol = do_speed;
        }
        HouseSystem::PullenSR => {
            // Q — Pullen SR "sinusoidal ratio" (swehouse.c:1336-1380)
            let third = 1.0 / 3.0;
            let two23 = (2.0_f64 * 2.0).powf(third);
            let mut acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
                acmc = diff_degrees(ac, mc);
            }
            let mut q = acmc;
            if q > 90.0 {
                q = 180.0 - q;
            }
            let (x, xr, xr3, xr4) = if q < 1e-30 {
                (0.0, 0.0, 0.0, 180.0)
            } else {
                let c = (180.0 - q) / q;
                let csq = c * c;
                let ccr = (csq - c).powf(third);
                let cqx = (two23 * ccr + 1.0).sqrt();
                let r1 = 0.5 * cqx;
                let r2 = 0.5 * (-2.0 * (1.0 - 2.0 * c) / cqx - two23 * ccr + 2.0).sqrt();
                let r = r1 + r2 - 0.5;
                let x = q / (2.0 * r + 1.0);
                let xr = r * x;
                let xr3 = xr * r * r;
                let xr4 = xr3 * r;
                (x, xr, xr3, xr4)
            };
            if acmc > 90.0 {
                cusps[11] = normalize_degrees(mc + xr3);
                cusps[12] = normalize_degrees(cusps[11] + xr4);
                cusps[2] = normalize_degrees(ac + xr);
                cusps[3] = normalize_degrees(cusps[2] + x);
            } else {
                cusps[11] = normalize_degrees(mc + xr);
                cusps[12] = normalize_degrees(cusps[11] + x);
                cusps[2] = normalize_degrees(ac + xr3);
                cusps[3] = normalize_degrees(cusps[2] + xr4);
            }
            do_interpol = do_speed;
        }
        HouseSystem::KrusinskiPisaGoelzer => {
            // U — Krusinski-Pisa (swehouse.c:1731-1805): great circle through Asc and zenith,
            // divided into 12 equal 30° arcs, projected back onto the ecliptic via meridian
            // circles. A sequence of `swe_cotrans` rotations, not a closed-form trig formula.
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
            }
            // A0-A5: rotate the Asc into the asc-zenith house-circle frame.
            let mut x = cotrans([ac, 0.0, 1.0], -eps); // A1: ecliptic -> equatorial
            x[0] -= th - 90.0; // A2: rotate by RA of east point
            x = cotrans(x, -(90.0 - geolat)); // A3: equatorial -> horizontal
            let kr_horizon_lon = x[0];
            x[0] -= x[0]; // A4: rotate to 0
            x = cotrans(x, -90.0); // A5: horizontal -> asc-zenith house-circle frame
            for i in 0..6usize {
                let mut xi = [30.0 * i as f64, 0.0, x[2]];
                xi = cotrans(xi, 90.0); // B1: house-circle -> horizontal
                xi[0] += kr_horizon_lon; // B2: rotate back
                xi = cotrans(xi, 90.0 - geolat); // B3: horizontal -> equatorial
                xi[0] = normalize_degrees(xi[0] + (th - 90.0)); // B4: RA of house cusp
                let mut cusp = atand(tand(xi[0]) / cosd(eps)); // B5: equatorial -> ecliptic
                if xi[0] > 90.0 && xi[0] <= 270.0 {
                    cusp = normalize_degrees(cusp + 180.0);
                }
                cusp = normalize_degrees(cusp);
                cusps[i + 1] = cusp;
                cusps[i + 7] = normalize_degrees(cusp + 180.0);
            }
            // No cusp_speed handling: 'U' is not in the do_interpol set — cusp_speed[1,4,7,10]
            // carry the stale pre-switch ac_speed/mc_speed, the rest stay zero (§4.2e).
        }
        HouseSystem::APC => {
            // Y — APC houses (swehouse.c:1806-1829), via the radians-domain `apc_sector` helper.
            for (i, cs) in cusps.iter_mut().enumerate().take(13).skip(1) {
                *cs = apc_sector(
                    i as i32,
                    geolat * crate::constants::DEGTORAD,
                    eps * crate::constants::DEGTORAD,
                    th * crate::constants::DEGTORAD,
                );
            }
            // MC near latitude 90 from apc_sector() is not accurate; use the real MC instead.
            cusps[10] = mc;
            cusps[4] = normalize_degrees(mc + 180.0);
            if geolat.abs() >= 90.0 - eps && diff_degrees(ac, mc) < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                mc = normalize_degrees(mc + 180.0);
                for cs in cusps.iter_mut().take(13).skip(1) {
                    *cs = normalize_degrees(*cs + 180.0);
                }
            }
            do_interpol = do_speed;
            // Y fills all 12 cusps itself via independent geometry — excluded from the
            // post-switch opposite-cusp mirror (swehouse.c:1985-2000, §3 step 3).
            skip_mirror = true;
        }
        HouseSystem::Sunshine | HouseSystem::SunshineAlt => {
            // I / i — Sunshine houses, Treindl / Makransky (swehouse.c:1156-1181). Stateless:
            // sundec is a required explicit parameter, not the C `static saved_sundec` cache
            // (c-ref-houses.md §11) — no global state, no `ascmc[9]==99` sentinel.
            let dec = match sundec {
                Some(d) if (-24.0..=24.0).contains(&d) => d,
                _ => {
                    return Err(Error::CError(
                        "House system Sunshine needs valid Sun declination".into(),
                    ));
                }
            };
            let acmc = diff_degrees(ac, mc);
            if acmc < 0.0 {
                ac = normalize_degrees(ac + 180.0);
                cusps[1] = ac;
                if hsys == HouseSystem::Sunshine {
                    mc = normalize_degrees(mc + 180.0);
                    cusps[10] = mc;
                }
            }
            cusps[4] = normalize_degrees(cusps[10] + 180.0);
            cusps[7] = normalize_degrees(cusps[1] + 180.0);
            let ok = if hsys == HouseSystem::Sunshine {
                sunshine_solution_treindl(th, geolat, eps, dec, &mut cusps)
            } else {
                sunshine_solution_makransky(th, geolat, eps, dec, &mut cusps)
            };
            if ok {
                do_interpol = do_speed;
                // I/i fill all 12 cusps themselves via independent geometry — excluded from
                // the post-switch opposite-cusp mirror (swehouse.c:1985-2000, §3 step 3).
                skip_mirror = true;
            } else {
                // retc==ERR (c-1e-6 degeneracy for Treindl, circumpolar Sun for Makransky):
                // fall back to Porphyry. `hsy` becomes 'O' in C, so the post-switch mirror
                // DOES run on this path (skip_mirror stays false) — fill_porphyry only fills
                // cusps[1,2,3,10,11,12]; the mirror fills [4,5,6,7,8,9].
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

/// Computes house cusps and angular points (`swe_houses_armc_ex2`) for house
/// system `hsys` from sidereal time `armc`, geographic latitude `geolat`, and
/// obliquity `eps`, including cusp/angle speeds via ARMC finite-differencing.
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

// ---------------------------------------------------------------------------
// Traditional sidereal houses — sidereal_houses_trad (swehouse.c:535-587)
// ---------------------------------------------------------------------------

/// Traditional sidereal house computation: compute tropical houses via
/// [`houses_armc`], then shift every cusp/ascmc entry (except `armc`) by the
/// already-computed ayanamsa. Kept pure — the caller (`Ephemeris::houses_ex2`)
/// resolves `ayanamsa` via `get_ayanamsa_ex` and passes it in.
pub fn sidereal_houses_trad(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: HouseSystem,
    sundec: Option<f64>,
    ayanamsa: f64,
) -> Result<HouseResult, Error> {
    // Whole Sign is computed as Equal, then re-fixed to 30-degree boundaries below.
    let ihs2 = if hsys == HouseSystem::WholeSign {
        HouseSystem::Equal
    } else {
        hsys
    };
    let mut result = houses_armc(armc, geolat, eps, ihs2, sundec)?;

    let ito = if hsys == HouseSystem::Gauquelin {
        36
    } else {
        12
    };
    for cusp in result.cusps.iter_mut().take(ito + 1).skip(1) {
        *cusp = normalize_degrees(*cusp - ayanamsa);
        if hsys == HouseSystem::WholeSign {
            *cusp -= *cusp % 30.0;
        }
    }
    // Equal-Aries: re-fixed to exact 0,30,60... boundaries (the ayanamsa shift above
    // would otherwise leave them off-multiple).
    if hsys == HouseSystem::EqualAries {
        for (i, cusp) in result.cusps.iter_mut().enumerate().take(13).skip(1) {
            *cusp = (i as f64 - 1.0) * 30.0;
        }
    }

    // ascmc[2] is armc — skipped, per swehouse.c:892 (loop excludes it).
    result.ascmc.ascendant = normalize_degrees(result.ascmc.ascendant - ayanamsa);
    result.ascmc.mc = normalize_degrees(result.ascmc.mc - ayanamsa);
    result.ascmc.vertex = normalize_degrees(result.ascmc.vertex - ayanamsa);
    result.ascmc.equatorial_ascendant =
        normalize_degrees(result.ascmc.equatorial_ascendant - ayanamsa);
    result.ascmc.coascendant_koch = normalize_degrees(result.ascmc.coascendant_koch - ayanamsa);
    result.ascmc.coascendant_munkasey =
        normalize_degrees(result.ascmc.coascendant_munkasey - ayanamsa);
    result.ascmc.polar_ascendant = normalize_degrees(result.ascmc.polar_ascendant - ayanamsa);

    Ok(result)
}

// ---------------------------------------------------------------------------
// Sidereal houses — ecl_t0 / ssypl geometric projections
// (swehouse.c:318-403, swehouse.c:425-532)
// ---------------------------------------------------------------------------

/// Shared tail of `sidereal_houses_ecl_t0`/`sidereal_houses_ssypl`: rotate the moving point
/// `pos`/`vel` (Cartesian, mean-equator-of-`tjde` frame) onto the *true* equator of `tjde` by
/// applying the mean-obliquity rotation, the nutation-in-longitude shift, then the true-obliquity
/// rotation back. Identical in both callers (swehouse.c:349-357, swehouse.c:498-506).
fn rotate_to_true_equator(pos: [f64; 3], vel: [f64; 3], eps: f64, nutlo: [f64; 2]) -> [f64; 6] {
    let eps_mean_rad = (eps - nutlo[1]) * crate::constants::DEGTORAD;
    let pos = rotate_x(pos, eps_mean_rad);
    let vel = rotate_x(vel, eps_mean_rad);
    let mut polsp = cartesian_to_polar_with_speed([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]);
    polsp[0] += nutlo[0] * crate::constants::DEGTORAD;
    let cart = polar_to_cartesian_with_speed(polsp);
    let eps_true_rad = eps * crate::constants::DEGTORAD;
    let pos = rotate_x([cart[0], cart[1], cart[2]], -eps_true_rad);
    let vel = rotate_x([cart[3], cart[4], cart[5]], -eps_true_rad);
    [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]
}

/// Shared core of `sidereal_houses_ecl_t0`/`sidereal_houses_ssypl`: given the fully-precessed
/// moving point `x` (Cartesian position+velocity, true-equator-of-`tjde` frame), compute the
/// auxiliary obliquity `epsx` and auxiliary vernal point `xvpx` from the point's instantaneous
/// orbital-plane normal, derive `armcx`, drive `houses_armc`, and return the raw (unsigned,
/// un-offset) `dvpxe` distance for the caller to finish combining with its own ayanamsa terms.
/// Ports the shared steps 1a/1b/2/3/4/5 (swehouse.c:358-388, swehouse.c:507-517).
fn sidereal_houses_geom_core(
    x: [f64; 6],
    armc: f64,
    lat: f64,
    hsys: HouseSystem,
    sundec: Option<f64>,
) -> Result<(HouseResult, f64), Error> {
    let xnorm = cross_prod([x[0], x[1], x[2]], [x[3], x[4], x[5]]);
    let rxy_sq = xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
    let rxyz = (rxy_sq + xnorm[2] * xnorm[2]).sqrt();
    let rxy = rxy_sq.sqrt();
    let epsx = asind(rxy / rxyz);

    let mut vz = x[5];
    if vz.abs() < 1e-15 {
        vz = 1e-15;
    }
    let fac = x[2] / vz;
    let sgn = vz / vz.abs();
    let xvpx = [
        (x[0] - fac * x[3]) * sgn,
        (x[1] - fac * x[4]) * sgn,
        (x[2] - fac * vz) * sgn,
    ];

    let x2 = cartesian_to_polar(xvpx);
    let dvpx = x2[0] * crate::constants::RADTODEG;
    let armcx = normalize_degrees(armc - dvpx);

    let result = houses_armc(armcx, lat, epsx, hsys, sundec)?;
    let dvpxe = acosd(dot_prod_unit([x[0], x[1], x[2]], xvpx));

    Ok((result, dvpxe))
}

/// Subtract `shift` from every cusp and `ascmc` entry (except `armc`), matching
/// swehouse.c:391-401/519-529's tail (both callers, shared with [`sidereal_houses_trad`]'s
/// equivalent loop). Re-fixes Equal-Aries cusps to exact 30-degree multiples afterward.
fn apply_sidereal_shift(mut result: HouseResult, hsys: HouseSystem, shift: f64) -> HouseResult {
    let ito = if hsys == HouseSystem::Gauquelin {
        36
    } else {
        12
    };
    for cusp in result.cusps.iter_mut().take(ito + 1).skip(1) {
        *cusp = normalize_degrees(*cusp - shift);
    }

    result.ascmc.ascendant = normalize_degrees(result.ascmc.ascendant - shift);
    result.ascmc.mc = normalize_degrees(result.ascmc.mc - shift);
    result.ascmc.vertex = normalize_degrees(result.ascmc.vertex - shift);
    result.ascmc.equatorial_ascendant =
        normalize_degrees(result.ascmc.equatorial_ascendant - shift);
    result.ascmc.coascendant_koch = normalize_degrees(result.ascmc.coascendant_koch - shift);
    result.ascmc.coascendant_munkasey =
        normalize_degrees(result.ascmc.coascendant_munkasey - shift);
    result.ascmc.polar_ascendant = normalize_degrees(result.ascmc.polar_ascendant - shift);

    if hsys == HouseSystem::EqualAries {
        for (i, cusp) in result.cusps.iter_mut().enumerate().take(13).skip(1) {
            *cusp = (i as f64 - 1.0) * 30.0;
        }
    }

    result
}

/// Sidereal houses projected onto the ecliptic of the ayanamsa epoch `t0`. Port of
/// `sidereal_houses_ecl_t0` (swehouse.c:318-403). `t0`/`ayan_t0` are
/// `EphemerisConfig::sidereal_t0`/`sidereal_ayan_t0` passed by the caller *unresolved* — unlike
/// `swi_get_ayanamsa_ex`'s callers, C's own `sidereal_houses_ecl_t0` reads `sip->t0` raw, with no
/// `t0_is_UT` deltaT adjustment (swehouse.c:341,350,388), so the Rust port must not apply
/// `ayanamsa::resolve_t0` here either.
#[allow(clippy::too_many_arguments)]
pub fn sidereal_houses_ecl_t0(
    tjde: f64,
    armc: f64,
    eps: f64,
    nutlo: [f64; 2],
    lat: f64,
    hsys: HouseSystem,
    sundec: Option<f64>,
    t0: f64,
    ayan_t0: f64,
    models: &AstroModels,
) -> Result<HouseResult, Error> {
    let epst0 = crate::obliquity::obliquity(t0, CalcFlags::empty(), models).eps;

    // Vernal point as a unit vector on the mean ecliptic of t0 (pos=[1,0,0], unit angular
    // velocity vel=[0,1,0]), rotated to the equator of t0.
    let mut pos = rotate_x([1.0, 0.0, 0.0], -epst0);
    let mut vel = rotate_x([0.0, 1.0, 0.0], -epst0);

    // t0 -> J2000 -> tjde (position and velocity precessed as separate 3-vectors, matching C's
    // two-call-per-leg structure).
    crate::precession::precess(
        &mut pos,
        t0,
        CalcFlags::empty(),
        models,
        PrecessionDirection::DateToJ2000,
    );
    crate::precession::precess(
        &mut pos,
        tjde,
        CalcFlags::empty(),
        models,
        PrecessionDirection::J2000ToDate,
    );
    crate::precession::precess(
        &mut vel,
        t0,
        CalcFlags::empty(),
        models,
        PrecessionDirection::DateToJ2000,
    );
    crate::precession::precess(
        &mut vel,
        tjde,
        CalcFlags::empty(),
        models,
        PrecessionDirection::J2000ToDate,
    );

    let x = rotate_to_true_equator(pos, vel, eps, nutlo);
    let (result, dvpxe) = sidereal_houses_geom_core(x, armc, lat, hsys, sundec)?;
    let dvpxe = if tjde < t0 { -dvpxe } else { dvpxe };

    Ok(apply_sidereal_shift(result, hsys, dvpxe + ayan_t0))
}

/// Sidereal houses projected onto the solar-system invariable plane. Port of
/// `sidereal_houses_ssypl` (swehouse.c:425-532). Same `t0`/`ayan_t0` raw-usage caveat as
/// [`sidereal_houses_ecl_t0`] applies.
#[allow(clippy::too_many_arguments)]
pub fn sidereal_houses_ssypl(
    tjde: f64,
    armc: f64,
    eps: f64,
    nutlo: [f64; 2],
    lat: f64,
    hsys: HouseSystem,
    sundec: Option<f64>,
    t0: f64,
    ayan_t0: f64,
    models: &AstroModels,
) -> Result<HouseResult, Error> {
    use crate::constants::{J2000, RADTODEG, SSY_PLANE_INCL, SSY_PLANE_NODE, SSY_PLANE_NODE_E2000};

    let eps2000 = crate::obliquity::obliquity(J2000, CalcFlags::empty(), models).eps;

    // Zero point on the solar-system rotation plane -> ecliptic 2000.
    let pos = rotate_x([1.0, 0.0, 0.0], -SSY_PLANE_INCL);
    let vel = rotate_x([0.0, 1.0, 0.0], -SSY_PLANE_INCL);
    let mut polsp = cartesian_to_polar_with_speed([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]);
    polsp[0] += SSY_PLANE_NODE_E2000;
    let cart = polar_to_cartesian_with_speed(polsp);

    // To equator 2000.
    let mut pos = rotate_x([cart[0], cart[1], cart[2]], -eps2000);
    let mut vel = rotate_x([cart[3], cart[4], cart[5]], -eps2000);

    // To mean equator of t.
    crate::precession::precess(
        &mut pos,
        tjde,
        CalcFlags::empty(),
        models,
        PrecessionDirection::J2000ToDate,
    );
    crate::precession::precess(
        &mut vel,
        tjde,
        CalcFlags::empty(),
        models,
        PrecessionDirection::J2000ToDate,
    );

    // To true equator of t (identical tail to ecl_t0).
    let x = rotate_to_true_equator(pos, vel, eps, nutlo);
    let (result, dvpxe) = sidereal_houses_geom_core(x, armc, lat, hsys, sundec)?;
    // Always positive for dates after 5400 BC (swehouse.c:516) -- no tjde<t0 sign flip here,
    // unlike ecl_t0.
    let dvpxe = dvpxe - SSY_PLANE_NODE * RADTODEG;

    // Ayanamsa between t0 and J2000, measured on the solar-system plane (swehouse.c:518-531).
    let mut x0 = [1.0, 0.0, 0.0];
    if t0 != J2000 {
        crate::precession::precess(
            &mut x0,
            t0,
            CalcFlags::empty(),
            models,
            PrecessionDirection::DateToJ2000,
        );
    }
    let x0 = rotate_x(x0, eps2000);
    let mut x0pol = cartesian_to_polar(x0);
    x0pol[0] -= SSY_PLANE_NODE_E2000;
    let x0cart = rotate_x(polar_to_cartesian(x0pol), SSY_PLANE_INCL);
    let x00 = (cartesian_to_polar(x0cart)[0] + SSY_PLANE_NODE) * RADTODEG;

    Ok(apply_sidereal_shift(result, hsys, dvpxe + ayan_t0 + x00))
}

// ---------------------------------------------------------------------------
// swe_house_pos — planet house position (swehouse.c:2216-2876)
// ---------------------------------------------------------------------------

/// Static `armc_to_mc` used by [`house_pos`]'s A/E/D/V/W and O/B/S branches
/// (swehouse.c:2149-2166): normalizes the raw `atand` result *before* the conditional `+180`,
/// then again after. A distinct bit pattern from `crate::math::armc_to_mc` (`swi_armc_to_mc`,
/// which normalizes only once, unconditionally, at the end) and from `mc_like` (`CalcH`'s inline
/// MC, which normalizes once unconditionally without the pre-normalize step) — see
/// c-ref-houses.md §12.1.
fn armc_to_mc_house_pos(armc: f64, eps: f64) -> f64 {
    let cose = cosd(eps);
    if (armc - 90.0).abs() > VERY_SMALL && (armc - 270.0).abs() > VERY_SMALL {
        let tant = tand(armc);
        let mut mc = normalize_degrees(atand(tant / cose));
        if armc > 90.0 && armc <= 270.0 {
            mc = normalize_degrees(mc + 180.0);
        }
        mc
    } else if (armc - 90.0).abs() <= VERY_SMALL {
        90.0
    } else {
        270.0
    }
}

/// Morinus's inline `tand(a)/cose` transform on the planet's own ecliptic longitude
/// (swehouse.c:2378-2390) — normalizes only inside the `+180` branch, unlike
/// [`armc_to_mc_house_pos`]'s pre-normalize step. Same shape as `swi_armc_to_mc`, just inlined
/// under a different case label with a different input.
fn mc_transform_raw(a: f64, cose: f64) -> f64 {
    if (a - 90.0).abs() > VERY_SMALL && (a - 270.0).abs() > VERY_SMALL {
        let tant = tand(a);
        let mut hpos = atand(tant / cose);
        if a > 90.0 && a <= 270.0 {
            hpos = normalize_degrees(hpos + 180.0);
        }
        hpos
    } else if (a - 90.0).abs() <= VERY_SMALL {
        90.0
    } else {
        270.0
    }
}

/// Bracket-and-interpolate a value against a 12-entry house-cusp table, handling both prograde
/// and retrograde cusp ordering (the `swe_difdeg2n(hcusp[6], hcusp[1])` sign test). Shared by the
/// `'J'` (Savard-A) branch and the `default:` fallback (`'L'`/`'Q'`) of `swe_house_pos`
/// (swehouse.c:2503-2531, 2842-2864). Returns the final house position directly — unlike most
/// other branches, this result is NOT divided by 30 afterward; the interpolation already accounts
/// for each house's (possibly non-uniform) angular span.
fn bracket_interpolate_12(hcusp: &[f64; 37], value: f64) -> f64 {
    let retrograde = diff_degrees(hcusp[6], hcusp[1]) <= 0.0;
    let mut i = 1usize;
    let mut c2 = 360.0;
    let d = if retrograde {
        normalize_degrees(hcusp[1] - value)
    } else {
        normalize_degrees(value - hcusp[1])
    };
    for k in 1..=12usize {
        let j = k + 1;
        c2 = if j > 12 {
            360.0
        } else if retrograde {
            normalize_degrees(hcusp[1] - hcusp[j])
        } else {
            normalize_degrees(hcusp[j] - hcusp[1])
        };
        i = k;
        if d < c2 {
            break;
        }
    }
    let c1 = if retrograde {
        normalize_degrees(hcusp[1] - hcusp[i])
    } else {
        normalize_degrees(hcusp[i] - hcusp[1])
    };
    let hsize = c2 - c1;
    if hsize == 0.0 {
        i as f64
    } else {
        i as f64 + (d - c1) / hsize
    }
}

/// Koch (`'K'`) inverse (swehouse.c:2398-2460): circumpolar-aware closed form. Can fail — returns
/// `Err` (matching C's `hpos=0`/`serr="Koch house position failed in circumpolar area"` path)
/// when the diurnal-arc factor falls outside `[0,2]`. The "doubtful result in circumpolar area"
/// warning path (object or MC circumpolar but `dfac` still in range) is not an error in C either;
/// we don't surface it since `house_pos` has no separate warning channel.
fn koch_house_pos(armc: f64, geolat: f64, eps: f64, mdd: f64, de: f64) -> Result<f64, Error> {
    let adp = if 90.0 - geolat < de || -90.0 - geolat > de {
        90.0
    } else if geolat - 90.0 > de || geolat + 90.0 < de {
        -90.0
    } else {
        asind(tand(geolat) * tand(de))
    };
    let mut admc = tand(eps) * tand(geolat) * sind(armc);
    if admc.abs() > 1.0 {
        admc = if admc > 1.0 { 1.0 } else { -1.0 };
    }
    admc = asind(admc);
    let samc = 90.0 + admc;
    let mut is_invalid = samc == 0.0;
    let mut xp0 = 0.0;
    if samc.abs() > 0.0 {
        let dfac;
        if mdd >= 0.0 {
            dfac = (mdd - adp + admc) / samc;
            xp0 = normalize_degrees((dfac - 1.0) * 90.0);
        } else {
            dfac = (mdd + 180.0 + adp + admc) / samc;
            xp0 = normalize_degrees((dfac + 1.0) * 90.0);
        }
        xp0 = normalize_degrees(xp0 + MILLIARCSEC);
        if !(0.0..=2.0).contains(&dfac) {
            is_invalid = true;
        }
    }
    if is_invalid {
        return Err(Error::CError(
            "Koch house position failed in circumpolar area".into(),
        ));
    }
    Ok(xp0 / 30.0 + 1.0)
}

/// Polich/Page "topocentric" (`'T'`) inverse (swehouse.c:2745-2801): binary search on the pole
/// height, distinct from `CalcH`'s Newton loops for the same system. Mirrors below-horizon and
/// western-hemisphere points into the canonical quadrant first, then mirrors the result back.
fn topocentric_house_pos(armc: f64, geolat: f64, ra: f64, de: f64, mdd: f64) -> f64 {
    let mut fh = geolat.clamp(-89.999, 89.999);
    let mut mdd = normalize_degrees(mdd);
    let mut de = de.clamp(-90.0 + VERY_SMALL, 90.0 - VERY_SMALL);
    let sinad = (tand(de) * tand(fh)).clamp(-1.0, 1.0);
    let is_above_hor = sinad + cosd(mdd) >= 0.0;
    let mut ra = ra;
    if !is_above_hor {
        ra = normalize_degrees(ra + 180.0);
        de = -de;
        mdd = normalize_degrees(mdd + 180.0);
    }
    if mdd > 180.0 {
        ra = normalize_degrees(armc - mdd);
    }
    let tanfi = tand(fh);
    let mut ra0 = normalize_degrees(armc + 90.0);
    let mut xp1 = 1.0f64;
    let mut fac = 2.0f64;
    let mut nloop = 0;
    while xp1.abs() > 0.000_001 && nloop < 1000 {
        if xp1 > 0.0 {
            fh = atand(tand(fh) - tanfi / fac);
            ra0 -= 90.0 / fac;
        } else {
            fh = atand(tand(fh) + tanfi / fac);
            ra0 += 90.0 / fac;
        }
        let xeq0 = normalize_degrees(ra - ra0);
        let xp = cotrans([xeq0, de, 1.0], 90.0 - fh);
        xp1 = xp[1];
        fac *= 2.0;
        nloop += 1;
    }
    let mut hpos = normalize_degrees(ra0 - armc);
    if mdd > 180.0 {
        hpos = normalize_degrees(-hpos);
    }
    if !is_above_hor {
        hpos = normalize_degrees(hpos + 180.0);
    }
    normalize_degrees(hpos - 90.0) / 30.0 + 1.0
}

/// Shared Sunshine (`'I'`/`'i'`)/APC (`'Y'`) geometric-approximation inverse
/// (swehouse.c:2650-2744). `dsun` is the Sun's declination for Sunshine, or the ascendant's
/// declination for APC. Seeds a Regiomontanus-style position line, then solves for where it
/// crosses the relevant body's (Sun or ascendant) diurnal/nocturnal semiarc.
fn sunshine_apc_house_pos(mdd: f64, de: f64, geolat: f64, dsun: f64) -> f64 {
    let geolat = geolat.clamp(-90.0 + MILLIARCSEC, 90.0 - MILLIARCSEC);
    let mut de = de;
    if 90.0 - de.abs() < VERY_SMALL {
        de = if de > 0.0 {
            90.0 - VERY_SMALL
        } else {
            -90.0 + VERY_SMALL
        };
    }
    let a = tand(geolat) * tand(de) + cosd(mdd);
    let mut xp0 = normalize_degrees(atand(-a / sind(mdd)));
    if mdd < 0.0 {
        xp0 += 180.0;
    }
    xp0 = normalize_degrees(xp0);
    let sinad = tand(de) * tand(geolat);
    let is_above_hor = sinad + cosd(mdd) >= 0.0;
    let harmc = if geolat < 0.0 {
        90.0 + geolat
    } else {
        90.0 - geolat
    };
    let mut darmc = normalize_degrees(xp0 - 270.0);
    let mut is_western_half = false;
    if darmc > 180.0 {
        is_western_half = true;
        darmc = 360.0 - darmc;
    }
    let sinad2 = tand(dsun) * tand(geolat);
    let ad = if sinad2 >= 1.0 {
        90.0
    } else if sinad2 <= -1.0 {
        -90.0
    } else {
        asind(sinad2)
    };
    let sad = 90.0 + ad;
    let san = 90.0 - ad;
    if sad == 0.0 && is_above_hor {
        xp0 = 270.0;
    } else if san == 0.0 && !is_above_hor {
        xp0 = 90.0;
    } else {
        let mut dsun = dsun;
        let mut sa = sad;
        if !is_above_hor {
            dsun = -dsun;
            sa = san;
            darmc = 180.0 - darmc;
            is_western_half = !is_western_half;
        }
        let mut a = acosd(cosd(harmc) * cosd(darmc));
        if a < VERY_SMALL {
            a = VERY_SMALL;
        }
        let sinpsi = (sind(harmc) / sind(a)).clamp(-1.0, 1.0);
        let y = sind(dsun) / sinpsi;
        let y = if y > 1.0 {
            90.0 - VERY_SMALL
        } else if y < -1.0 {
            -(90.0 - VERY_SMALL)
        } else {
            asind(y)
        };
        let mut d = acosd(cosd(y) / cosd(dsun));
        if dsun < 0.0 {
            d = -d;
        }
        if geolat < 0.0 {
            d = -d;
        }
        darmc += d;
        xp0 = if is_western_half {
            270.0 - (darmc / sa) * 90.0
        } else {
            270.0 + (darmc / sa) * 90.0
        };
        if !is_above_hor {
            xp0 = normalize_degrees(xp0 + 180.0);
        }
    }
    xp0 = normalize_degrees(xp0 + MILLIARCSEC);
    xp0 / 30.0 + 1.0
}

/// Inverse problem: given `(armc, geolat, eps, hsys, xpin=[ecl.lon, ecl.lat])`, return a
/// continuous house position `1.0..13.0`. Port of `swe_house_pos` (swehouse.c:2216-2876). See
/// c-ref-houses.md §8.
///
/// Stateless departure from C: Sunshine's Sun declination and APC's ascendant declination are
/// resolved explicitly (`sundec` parameter / the pre-check's own `houses_armc` call) rather than
/// via C's `ascmc[9]==99` sentinel + `static double saved_sundec` cache (c-ref-houses.md §11).
/// Sunshine requires a valid `sundec` and returns `Err` rather than silently falling back to a
/// stale or zero value.
pub fn house_pos(
    armc: f64,
    geolat: f64,
    eps: f64,
    hsys: HouseSystem,
    xpin: [f64; 2],
    sundec: Option<f64>,
) -> Result<f64, Error> {
    let sine = sind(eps);
    let cose = cosd(eps);

    // Pre-check (swehouse.c:2231-2266): does xpin exactly match a cusp? Also the only source of
    // the ascendant needed for APC's declination-of-ascendant `dsun`. Best-effort: if `hsys` has
    // no `houses_armc` support yet (currently only Alcabitius), the shortcut/APC path is simply
    // unavailable here — matches C, which doesn't abort the whole call on this step's failure
    // either (it just skips the shortcut and proceeds into the switch with dsun=0).
    let precheck = houses_armc(armc, geolat, eps, hsys, sundec).ok();
    if let Some(ref r) = precheck {
        for i in 1..=12usize {
            if diff_degrees(xpin[0], r.cusps[i]).abs() < MILLIARCSEC && xpin[1] == 0.0 {
                return Ok(i as f64);
            }
        }
    }

    let xeq = cotrans([xpin[0], xpin[1], 1.0], -eps);
    let ra = xeq[0];
    let de = xeq[1];
    let mdd_raw = normalize_degrees(ra - armc);
    let mdn_raw = normalize_degrees(mdd_raw + 180.0);
    let mdd = if mdd_raw >= 180.0 {
        mdd_raw - 360.0
    } else {
        mdd_raw
    };
    let mdn = if mdn_raw >= 180.0 {
        mdn_raw - 360.0
    } else {
        mdn_raw
    };

    let dsun = match hsys {
        HouseSystem::Sunshine | HouseSystem::SunshineAlt => match sundec {
            Some(d) if (-24.0..=24.0).contains(&d) => d,
            _ => {
                return Err(Error::CError(
                    "House system Sunshine needs valid Sun declination".into(),
                ));
            }
        },
        HouseSystem::APC => {
            let asc = precheck
                .as_ref()
                .map(|r| r.ascmc.ascendant)
                .ok_or_else(|| {
                    Error::CError("APC house position requires a computable ascendant".into())
                })?;
            cotrans([asc, 0.0, 1.0], -eps)[1]
        }
        _ => 0.0,
    };

    let hpos = match hsys {
        HouseSystem::EqualAries => xpin[0] / 30.0 + 1.0,
        HouseSystem::Equal
        | HouseSystem::EqualMC
        | HouseSystem::Vehlow
        | HouseSystem::WholeSign => {
            let mut asc = asc1(normalize_degrees(armc + 90.0), geolat, sine, cose);
            let mc = armc_to_mc_house_pos(armc, eps);
            asc = fix_asc_polar(asc, armc, eps, geolat);
            let mut xp0 = normalize_degrees(xpin[0] - asc);
            if hsys == HouseSystem::Vehlow {
                xp0 = normalize_degrees(xp0 + 15.0);
            }
            if hsys == HouseSystem::WholeSign {
                xp0 = normalize_degrees(xp0 + asc % 30.0);
            }
            if hsys == HouseSystem::EqualMC {
                xp0 = normalize_degrees(xpin[0] - mc - 90.0);
            }
            xp0 = normalize_degrees(xp0 + MILLIARCSEC);
            xp0 / 30.0 + 1.0
        }
        HouseSystem::Porphyry | HouseSystem::Sripati | HouseSystem::Alcabitius => {
            let asc = fix_asc_polar(
                asc1(normalize_degrees(armc + 90.0), geolat, sine, cose),
                armc,
                eps,
                geolat,
            );
            if hsys == HouseSystem::Alcabitius {
                let dek = asind(sind(asc) * sine);
                let tanfi = tand(geolat);
                let r = -tanfi * tand(dek);
                let sda = acosd(r);
                let sna = 180.0 - sda;
                let mut hpos = if mdd > 0.0 {
                    if mdd < sda {
                        mdd * 90.0 / sda
                    } else {
                        90.0 + (mdd - sda) * 90.0 / sna
                    }
                } else if mdd > -sna {
                    360.0 + mdd * 90.0 / sna
                } else {
                    270.0 + (mdd + sna) * 90.0 / sda
                };
                hpos = normalize_degrees(hpos - 90.0) / 30.0 + 1.0;
                if hpos >= 13.0 {
                    hpos -= 12.0;
                }
                hpos
            } else {
                let mc = armc_to_mc_house_pos(armc, eps);
                let mut xp0 = normalize_degrees(xpin[0] - asc);
                xp0 = normalize_degrees(xp0 + MILLIARCSEC);
                let mut hpos;
                if xp0 < 180.0 {
                    hpos = 1.0;
                } else {
                    hpos = 7.0;
                    xp0 -= 180.0;
                }
                let acmc = diff_degrees(asc, mc);
                if xp0 < 180.0 - acmc {
                    hpos += xp0 * 3.0 / (180.0 - acmc);
                } else {
                    hpos += 3.0 + (xp0 - 180.0 + acmc) * 3.0 / acmc;
                }
                if hsys == HouseSystem::Sripati {
                    hpos += 0.5;
                    if hpos > 12.0 {
                        hpos = 1.0;
                    }
                }
                hpos
            }
        }
        HouseSystem::Meridian => normalize_degrees(mdd - 90.0) / 30.0 + 1.0,
        HouseSystem::Carter => {
            let mut x0 = asc1(normalize_degrees(armc + 90.0), geolat, sine, cose);
            x0 = fix_asc_polar(x0, armc, eps, geolat);
            let xeq0 = cotrans([x0, 0.0, 1.0], -eps);
            normalize_degrees(ra - xeq0[0]) / 30.0 + 1.0
        }
        HouseSystem::Morinus => {
            let hpos_mc = mc_transform_raw(xpin[0], cose);
            normalize_degrees(hpos_mc - armc - 90.0) / 30.0 + 1.0
        }
        HouseSystem::Koch => koch_house_pos(armc, geolat, eps, mdd, de)?,
        HouseSystem::Campanus => {
            let xeq0 = normalize_degrees(mdd - 90.0);
            let mut xp = cotrans([xeq0, de, 1.0], -geolat);
            xp[0] = normalize_degrees(xp[0] + MILLIARCSEC);
            xp[0] / 30.0 + 1.0
        }
        HouseSystem::Horizon => {
            let xeq0 = normalize_degrees(mdd - 90.0);
            let mut xp = cotrans([xeq0, de, 1.0], 90.0 - geolat);
            xp[0] = normalize_degrees(xp[0] + MILLIARCSEC);
            xp[0] / 30.0 + 1.0
        }
        HouseSystem::SavardA => {
            let sinfi = sind(geolat);
            let (mut xs2, mut xs1) = if geolat.abs() < VERY_SMALL {
                (1.0 / 3.0, 2.0 / 3.0)
            } else {
                (sind(geolat / 3.0) / sinfi, sind(2.0 * geolat / 3.0) / sinfi)
            };
            xs2 = asind(xs2);
            xs1 = asind(xs1);
            let mut hcusp = [0.0; 37];
            hcusp[1] = 0.0;
            hcusp[2] = xs2;
            hcusp[3] = xs1;
            hcusp[4] = 90.0;
            hcusp[5] = 180.0 - xs1;
            hcusp[6] = 180.0 - xs2;
            hcusp[7] = 180.0;
            hcusp[8] = 180.0 + xs2;
            hcusp[9] = 180.0 + xs1;
            hcusp[10] = 270.0;
            hcusp[11] = 360.0 - xs1;
            hcusp[12] = 360.0 - xs2;
            let xeq0 = normalize_degrees(mdd - 90.0);
            let xp = cotrans([xeq0, de, 1.0], -geolat);
            bracket_interpolate_12(&hcusp, xp[0])
        }
        HouseSystem::KrusinskiPisaGoelzer => {
            let mut geolat = geolat;
            if geolat.abs() < VERY_SMALL {
                geolat = if geolat >= 0.0 {
                    VERY_SMALL
                } else {
                    -VERY_SMALL
                };
            }
            let mut asc = asc1(normalize_degrees(armc + 90.0), geolat, sine, cose);
            asc = fix_asc_polar(asc, armc, eps, geolat);
            // I. plane of the 'asc-zenith' great circle relative to the equator.
            let mut x = cotrans([asc, 0.0, 1.0], -eps);
            let raep = normalize_degrees(armc + 90.0);
            x[0] = normalize_degrees(raep - x[0]);
            x = cotrans(x, -(90.0 - geolat));
            let tanx = tand(x[0]);
            let mut xtemp = if geolat == 0.0 {
                if tanx >= 0.0 { 90.0 } else { -90.0 }
            } else {
                atand(tanx / cosd(90.0 - geolat))
            };
            if x[0] > 90.0 && x[0] <= 270.0 {
                xtemp = normalize_degrees(xtemp + 180.0);
            }
            x[0] = normalize_degrees(xtemp);
            let raaz = normalize_degrees(raep - x[0]);
            // Ib. obliquity to the equator of the house plane.
            let mut xb = [raaz, 0.0, 1.0];
            xb[0] = normalize_degrees(raep - xb[0]);
            xb = cotrans(xb, -(90.0 - geolat));
            xb[1] += 90.0;
            xb = cotrans(xb, 90.0 - geolat);
            let oblaz = xb[1];
            // IIa. Asc on the house plane, relative to raaz.
            let mut xasc = cotrans([asc, 0.0, 1.0], -eps);
            xasc[0] = normalize_degrees(xasc[0] - raaz);
            let mut xtemp2 = atand(tand(xasc[0]) / cosd(oblaz));
            if xasc[0] > 90.0 && xasc[0] <= 270.0 {
                xtemp2 = normalize_degrees(xtemp2 + 180.0);
            }
            xasc[0] = normalize_degrees(xtemp2);
            // IIb. planet on the house plane, relative to raaz. (The declination-circle offset
            // xp[1] that C computes next is vestigial — never read — so it's omitted here.)
            let mut xp0 = normalize_degrees(ra - raaz);
            let mut xtemp3 = atand(tand(xp0) / cosd(oblaz));
            if xp0 > 90.0 && xp0 <= 270.0 {
                xtemp3 = normalize_degrees(xtemp3 + 180.0);
            }
            xp0 = normalize_degrees(xtemp3);
            xp0 = normalize_degrees(xp0 - xasc[0]);
            xp0 = normalize_degrees(xp0 + MILLIARCSEC);
            xp0 / 30.0 + 1.0
        }
        HouseSystem::Regiomontanus => {
            let xp0 = if mdd.abs() < VERY_SMALL {
                270.0
            } else if 180.0 - mdd.abs() < VERY_SMALL {
                90.0
            } else {
                let mut geolat = geolat;
                if 90.0 - geolat.abs() < VERY_SMALL {
                    geolat = if geolat > 0.0 {
                        90.0 - VERY_SMALL
                    } else {
                        -90.0 + VERY_SMALL
                    };
                }
                let mut de = de;
                if 90.0 - de.abs() < VERY_SMALL {
                    de = if de > 0.0 {
                        90.0 - VERY_SMALL
                    } else {
                        -90.0 + VERY_SMALL
                    };
                }
                let a = tand(geolat) * tand(de) + cosd(mdd);
                let mut xp0 = normalize_degrees(atand(-a / sind(mdd)));
                if mdd < 0.0 {
                    xp0 += 180.0;
                }
                xp0 = normalize_degrees(xp0);
                normalize_degrees(xp0 + MILLIARCSEC)
            };
            xp0 / 30.0 + 1.0
        }
        HouseSystem::Sunshine | HouseSystem::SunshineAlt | HouseSystem::APC => {
            sunshine_apc_house_pos(mdd, de, geolat, dsun)
        }
        HouseSystem::PolichPage => topocentric_house_pos(armc, geolat, ra, de, mdd),
        HouseSystem::Placidus | HouseSystem::Gauquelin => {
            let mut xp0 = if 90.0 - de.abs() <= geolat.abs() {
                if de * geolat < 0.0 {
                    normalize_degrees(90.0 + mdn / 2.0)
                } else {
                    normalize_degrees(270.0 + mdd / 2.0)
                }
            } else {
                let sinad = tand(de) * tand(geolat);
                let ad = asind(sinad);
                let is_above_hor = sinad + cosd(mdd) >= 0.0;
                let sad = 90.0 + ad;
                let san = 90.0 - ad;
                let xp0 = if is_above_hor {
                    (mdd / sad + 3.0) * 90.0
                } else {
                    (mdn / san + 1.0) * 90.0
                };
                normalize_degrees(xp0 + MILLIARCSEC)
            };
            if hsys == HouseSystem::Gauquelin {
                xp0 = 360.0 - xp0;
                xp0 / 10.0 + 1.0
            } else {
                xp0 / 30.0 + 1.0
            }
        }
        HouseSystem::PullenSD | HouseSystem::PullenSR => {
            let hcusp = match &precheck {
                Some(r) => r.cusps,
                None => houses_armc(armc, geolat, eps, hsys, sundec)?.cusps,
            };
            bracket_interpolate_12(&hcusp, xpin[0])
        }
    };

    Ok(hpos)
}
