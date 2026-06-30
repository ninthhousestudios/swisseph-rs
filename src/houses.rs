//! Port of `swehouse.c`: house systems (`CalcH`) and the ARMC-based driver
//! (`swe_houses_armc_ex2`). See `docs/c-ref-houses.md`.

use crate::error::Error;
use crate::math::{diff_degrees, normalize_degrees};
use crate::types::HouseSystem;

// ---------------------------------------------------------------------------
// Constants (swehouse.h:87, swehouse.c:68-70, swehouse.c:940)
// ---------------------------------------------------------------------------

const VERY_SMALL: f64 = 1e-10;
#[allow(dead_code)] // used by Placidus/Gauquelin Newton iteration (later sub-tasks)
const VERY_SMALL_PLAC_ITER: f64 = 1.0 / 360000.0;
#[allow(dead_code)] // used by swe_house_pos (later sub-tasks)
const MILLIARCSEC: f64 = 1.0 / 3600000.0;
const SOLAR_YEAR: f64 = 365.242_198_93;
const ARMCS: f64 = (SOLAR_YEAR + 1.0) / SOLAR_YEAR * 360.0;
#[allow(dead_code)] // used by Placidus/Gauquelin Newton iteration (later sub-tasks)
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
    // Used by the iterative house systems (Placidus, Koch, Gauquelin, ...) added later.
    let _tane = tand(eps);

    let mut geolat = geolat;
    if (geolat.abs() - 90.0).abs() < VERY_SMALL {
        geolat = if geolat < 0.0 {
            -90.0 + VERY_SMALL
        } else {
            90.0 - VERY_SMALL
        };
    }
    let _tanfi = tand(geolat);

    let mc = mc_like(th, cose);
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

    let do_interpol = false;

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
        _ => {
            return Err(Error::CError(format!(
                "house system {hsys:?} not yet implemented"
            )));
        }
    }

    // Post-switch opposite-cusp mirror (swehouse.c:1985-2000) — skipped only for G/Y/I,
    // none of which are reachable here yet.
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
