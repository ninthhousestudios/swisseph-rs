//! Rise / set / meridian transit. Port of `swe_rise_trans_true_hor` and `calc_mer_trans`
//! (`swecl.c`). See `docs/c-ref-riseset.md`. The fast-path optimization (`rise_set_fast`) and
//! the top-level `swe_rise_trans` dispatcher are a separate module (RSE 4).
//!
//! [`Ephemeris::rise_trans_true_hor`](crate::context::Ephemeris::rise_trans_true_hor) is the
//! public entry point (in `context.rs`, matching how [`crate::azalt`]'s geometry cores are
//! wrapped by `Ephemeris::azalt`/`azalt_rev`); the functions here take `&Ephemeris` explicitly
//! since the algorithm interleaves many `calc`/`azalt`/`azalt_rev` calls.

use crate::azalt::{AzAltDir, HorDir};
use crate::constants::{AUNIT, LAPSE_RATE, RADTODEG, RISE_SET_GEOALT_MAX, RISE_SET_GEOALT_MIN};
use crate::context::{Ephemeris, EphemerisConfig, TopoPosition};
use crate::error::Error;
use crate::flags::{CalcFlags, RiseSetFlags};
use crate::types::Body;

/// Result of a rise/set/transit search: a single Julian Day, UT. Port of C's `tret[0]`
/// (`swe_rise_trans_true_hor` only ever fills a single time slot for this family of events).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RiseSetResult {
    pub time: f64,
}

/// 15-point culmination pre-pass sampling window: `jmax = 14` (15 points), spaced `twohrs`
/// apart, spanning `[tjd_ut - 2h, tjd_ut + 26h]`. See docs/c-ref-riseset.md §5.2.
const JMAX: usize = 14;
const TWOHRS: f64 = 1.0 / 12.0;

/// `iflag` selector bits kept from the caller's `epheflag` in the rise/set branch
/// (docs/c-ref-riseset.md §5.1, swecl.c:4425): the ephemeris-source bits are vestigial in this
/// stateless port (the backend is `Ephemeris`'s own configured `ephemeris_source`, same as every
/// other `calc`/`calc_ut` call), but the mask is replicated for structural fidelity.
const IFLAG_KEEP_MASK: CalcFlags = CalcFlags::JPLEPH
    .union(CalcFlags::SWIEPH)
    .union(CalcFlags::MOSEPH)
    .union(CalcFlags::NONUT)
    .union(CalcFlags::TRUEPOS);

/// `iflag` selector bits kept from the caller's `epheflag` in `calc_mer_trans` (swecl.c:4701:
/// `iflag &= SEFLG_EPHMASK`) -- deliberately narrower than [`IFLAG_KEEP_MASK`]: C drops
/// `NONUT`/`TRUEPOS` here, unlike the rise/set branch.
const MER_TRANS_IFLAG_KEEP_MASK: CalcFlags = CalcFlags::JPLEPH
    .union(CalcFlags::SWIEPH)
    .union(CalcFlags::MOSEPH);

/// Full rise/set/transit algorithm. Port of `swe_rise_trans_true_hor` (swecl.c:4387-4686).
#[allow(clippy::too_many_arguments)]
pub(crate) fn rise_trans_true_hor(
    eph: &Ephemeris,
    tjd_ut: f64,
    body: Body,
    starname: Option<&str>,
    epheflag: CalcFlags,
    mut rsmi: RiseSetFlags,
    geopos: [f64; 3],
    atpress: f64,
    attemp: f64,
    mut horhgt: f64,
) -> Result<RiseSetResult, Error> {
    if !(RISE_SET_GEOALT_MIN..=RISE_SET_GEOALT_MAX).contains(&geopos[2]) {
        return Err(Error::CError(format!(
            "observer altitude {} outside valid range [{RISE_SET_GEOALT_MIN}, {RISE_SET_GEOALT_MAX}]",
            geopos[2]
        )));
    }

    // Unlike `azalt`/`refrac_extended`, `calc_dip` does not auto-estimate atpress from height
    // when atpress == 0 -- C passes the caller's raw atpress straight through (swecl.c:4415-4416,
    // calc_dip itself at swecl.c:3158-3168). Do not route this through `resolve_atpress`.
    if horhgt == -100.0 {
        horhgt = 0.0001 + crate::azalt::calc_dip(geopos[2], atpress, attemp, LAPSE_RATE);
    }

    // Pluto asteroid-number alias (swecl.c:4402-4404): ipl == SE_AST_OFFSET + 134340.
    let body = match body {
        Body::Asteroid(id) if id.mpc_number() == 134340 => Body::Pluto,
        b => b,
    };

    let mut iflag = epheflag & IFLAG_KEEP_MASK;
    let geoctr_no_ecl_lat = rsmi.contains(RiseSetFlags::GEOCTR_NO_ECL_LAT);
    let tohor_flag = if geoctr_no_ecl_lat {
        AzAltDir::EclToHor
    } else {
        iflag |= CalcFlags::EQUATORIAL | CalcFlags::TOPOCTR;
        AzAltDir::EquToHor
    };

    // TOPOCTR needs a topographic position -- thread `geopos` through a local config override
    // rather than requiring it to match `eph`'s own configured topographic position (mirrors
    // C's per-call `swe_set_topo`, but stateless; see Ephemeris::calc_with_config's doc comment).
    let topo_config = {
        let mut c = eph.config().clone();
        c.topographic = Some(TopoPosition {
            longitude: geopos[0],
            latitude: geopos[1],
            altitude: geopos[2],
        });
        c
    };

    if rsmi.intersects(RiseSetFlags::MTRANSIT | RiseSetFlags::ITRANSIT) {
        return calc_mer_trans(
            eph,
            tjd_ut,
            body,
            starname,
            epheflag,
            rsmi,
            geopos,
            &topo_config,
        );
    }

    if !rsmi.intersects(RiseSetFlags::RISE | RiseSetFlags::SET) {
        rsmi |= RiseSetFlags::RISE;
    }

    let is_fixstar = starname.is_some_and(|s| !s.is_empty());

    if !is_fixstar
        && body == Body::Sun
        && let Some(depression) = twilight_depression(rsmi)
    {
        rsmi |= RiseSetFlags::NO_REFRACTION | RiseSetFlags::DISC_CENTER;
        horhgt = -depression;
    }

    // Fixed-star position computed ONCE and reused for the whole search (swecl.c:4457-4463,
    // §5.2 step 3) -- proper motion is negligible over the ~28h sampling window. NOTE: unlike
    // C, this does not thread TOPOCTR through the fixed-star pipeline (`calc_fixstar` has no
    // topocentric parallax support yet) -- parallax at stellar distances is negligible, and
    // this path has no golden coverage (gen_riseset.c only exercises SE_SUN/SE_MOON).
    let star_pos: Option<[f64; 3]> = if is_fixstar {
        let name = starname.unwrap();
        let (_, result) = eph.fixstar2_ut(name, tjd_ut, iflag)?;
        Some([result.data[0], result.data[1], result.data[2]])
    } else {
        None
    };

    let dd_m = disc_diameter_m(body, is_fixstar, rsmi);

    // Position at `t` (UT): reused fixed-star position, or a fresh topocentric `calc_ut`.
    let resolve_xc = |t: f64| -> Result<[f64; 3], Error> {
        let mut xc = match star_pos {
            Some(p) => p,
            None => {
                let r = eph.calc_ut_with_config(t, body, iflag, &topo_config)?;
                [r.data[0], r.data[1], r.data[2]]
            }
        };
        if geoctr_no_ecl_lat {
            xc[1] = 0.0;
        }
        Ok(xc)
    };

    // Signed disc-limb offset (degrees) at a resolved position.
    let rdi_of = |xc: [f64; 3]| -> f64 {
        let radius = disc_radius_deg(rsmi, body, dd_m, xc[2]);
        if rsmi.contains(RiseSetFlags::DISC_BOTTOM) {
            -radius
        } else {
            radius
        }
    };

    // Mesh sample at `t`: (h, detect_alt) -- port of §5.2 steps 4-10 / §5.3 / §5.4's per-point
    // recompute.
    let sample = |t: f64| -> Result<(f64, f64), Error> {
        let xc = resolve_xc(t)?;
        let rdi = rdi_of(xc);
        Ok(mesh_h_and_detect(
            eph,
            t,
            [xc[0], xc[1]],
            rdi,
            rsmi,
            tohor_flag,
            geopos,
            atpress,
            attemp,
            horhgt,
        ))
    };

    // Culmination-refinement sample at `t`: always true-limb-altitude minus `horhgt`,
    // regardless of the `NO_REFRACTION` bit (§5.2 step 11, second bullet).
    let refine_sample = |t: f64| -> Result<f64, Error> {
        let xc = resolve_xc(t)?;
        let rdi = rdi_of(xc);
        let [_, true_alt] = azimuth_true_limb_alt(
            eph,
            t,
            [xc[0], xc[1]],
            rdi,
            tohor_flag,
            geopos,
            atpress,
            attemp,
        );
        Ok(true_alt - horhgt)
    };

    // --- §5.2: 15-point culmination pre-pass ------------------------------------------------
    let mut tc: Vec<f64> = Vec::with_capacity(JMAX + 5);
    let mut hv: Vec<f64> = Vec::with_capacity(JMAX + 5);
    let mut detect: Vec<f64> = Vec::with_capacity(JMAX + 1);
    let mut tculm: Vec<f64> = Vec::new();

    for ii in 0..=JMAX {
        let t = tjd_ut - TWOHRS + (ii as f64) * TWOHRS;
        let (h, detect_alt) = sample(t)?;
        tc.push(t);
        hv.push(h);
        detect.push(detect_alt);

        if ii > 1 {
            let dc = [detect[ii - 2], detect[ii - 1], detect[ii]];
            let is_max = dc[1] > dc[0] && dc[1] > dc[2];
            let is_min = dc[1] < dc[0] && dc[1] < dc[2];
            if is_max || is_min {
                let mut dt = TWOHRS;
                let mut tcu = t - dt;
                let (dtint, _) = crate::math::find_maximum(dc[0], dc[1], dc[2], dt);
                tcu += dtint + dt;
                dt /= 3.0;
                while dt > 0.0001 {
                    let a0 = refine_sample(tcu - dt)?;
                    let a1 = refine_sample(tcu)?;
                    let a2 = refine_sample(tcu + dt)?;
                    let (dtint, _) = crate::math::find_maximum(a0, a1, a2, dt);
                    tcu += dtint + dt;
                    dt /= 3.0;
                }
                tculm.push(tcu);
            }
        }
    }

    // --- §5.3: insert culminations into the mesh --------------------------------------------
    for &tcu in &tculm {
        let (h, _) = sample(tcu)?;
        let slot = tc.partition_point(|&x| x <= tcu);
        tc.insert(slot, tcu);
        hv.insert(slot, h);
    }

    // --- §5.4: zero-crossing search, sign change + 20-iteration bisection -------------------
    for ii in 1..tc.len() {
        if hv[ii - 1] * hv[ii] >= 0.0 {
            continue;
        }
        let rising = hv[ii - 1] < hv[ii];
        if rising && !rsmi.contains(RiseSetFlags::RISE) {
            continue;
        }
        if !rising && !rsmi.contains(RiseSetFlags::SET) {
            continue;
        }

        let mut t2 = [tc[ii - 1], tc[ii]];
        let mut dc = [hv[ii - 1], hv[ii]];
        let mut t = t2[0];
        for _ in 0..20 {
            t = (t2[0] + t2[1]) / 2.0;
            let (aha, _) = sample(t)?;
            if aha * dc[0] <= 0.0 {
                dc[1] = aha;
                t2[1] = t;
            } else {
                dc[0] = aha;
                t2[0] = t;
            }
        }

        if t > tjd_ut {
            return Ok(RiseSetResult { time: t });
        }
    }

    Err(Error::CircumpolarBody)
}

/// Disc angular radius (degrees) at a given distance. `dd_m` is the resolved disc diameter
/// (meters; `0` for fixed stars / `DISC_CENTER`). Port of §5.2 steps 5-6.
fn disc_radius_deg(rsmi: RiseSetFlags, body: Body, dd_m: f64, dist_au: f64) -> f64 {
    if dd_m == 0.0 {
        return 0.0;
    }
    let curdist = if rsmi.contains(RiseSetFlags::FIXED_DISC_SIZE) {
        match body {
            Body::Sun => 1.0,
            Body::Moon => 0.00257,
            _ => dist_au,
        }
    } else {
        dist_au
    };
    (dd_m / 2.0 / AUNIT / curdist).asin() * RADTODEG
}

/// Disc diameter (meters), resolved once per search. `0` for fixed stars or `DISC_CENTER`;
/// `PLANETARY_DIAMETERS[raw_id]` for `raw_id` in range; else `0` (asteroid `ast_diam` is not
/// ported -- untested path, `gen_riseset.c` only covers SE_SUN/SE_MOON). Port of §5.2 step 5.
fn disc_diameter_m(body: Body, is_fixstar: bool, rsmi: RiseSetFlags) -> f64 {
    if is_fixstar || rsmi.contains(RiseSetFlags::DISC_CENTER) {
        return 0.0;
    }
    let raw = body.to_raw_id();
    if (0..crate::constants::PLANETARY_DIAMETERS.len() as i32).contains(&raw) {
        crate::constants::PLANETARY_DIAMETERS[raw as usize]
    } else {
        0.0
    }
}

/// Civil/nautical/astronomical twilight depression (degrees). Priority: each bit-check is
/// unconditional (not `else if`), so if multiple are set simultaneously astronomical
/// (checked last) wins. Port of `rdi_twilight` (swecl.c:4164-4174).
fn twilight_depression(rsmi: RiseSetFlags) -> Option<f64> {
    let mut depression = None;
    if rsmi.contains(RiseSetFlags::CIVIL_TWILIGHT) {
        depression = Some(6.0);
    }
    if rsmi.contains(RiseSetFlags::NAUTIC_TWILIGHT) {
        depression = Some(12.0);
    }
    if rsmi.contains(RiseSetFlags::ASTRO_TWILIGHT) {
        depression = Some(18.0);
    }
    depression
}

/// `[azimuth, true-altitude-with-limb-adjustment]` at `t`, before horhgt/refraction handling.
/// Shared first step of every altitude sample (docs/c-ref-riseset.md §5.2 steps 7-8).
#[allow(clippy::too_many_arguments)]
fn azimuth_true_limb_alt(
    eph: &Ephemeris,
    t: f64,
    xc: [f64; 2],
    rdi_signed: f64,
    tohor_flag: AzAltDir,
    geopos: [f64; 3],
    atpress: f64,
    attemp: f64,
) -> [f64; 2] {
    let xh = eph.azalt(t, tohor_flag, geopos, atpress, attemp, LAPSE_RATE, xc);
    [xh[0], xh[1] + rdi_signed]
}

/// Mesh altitude sample: returns `(h, detect_alt)`. `h` is the horhgt/refraction-adjusted
/// value compared against zero for rise/set; `detect_alt` is the value used for culmination
/// detection (the un-refracted limb altitude, except when `NO_REFRACTION` is set, in which case
/// it equals `h` -- both are the same mutated `xh[ii][1]` in C). Port of §5.2 steps 7-10.
#[allow(clippy::too_many_arguments)]
fn mesh_h_and_detect(
    eph: &Ephemeris,
    t: f64,
    xc: [f64; 2],
    rdi_signed: f64,
    rsmi: RiseSetFlags,
    tohor_flag: AzAltDir,
    geopos: [f64; 3],
    atpress: f64,
    attemp: f64,
    horhgt: f64,
) -> (f64, f64) {
    let [az, true_alt] =
        azimuth_true_limb_alt(eph, t, xc, rdi_signed, tohor_flag, geopos, atpress, attemp);
    if rsmi.contains(RiseSetFlags::NO_REFRACTION) {
        let v = true_alt - horhgt;
        (v, v)
    } else {
        let xc2 = eph.azalt_rev(t, HorDir::HorToEqu, geopos, [az, true_alt]);
        let xh2 = eph.azalt(
            t,
            AzAltDir::EquToHor,
            geopos,
            atpress,
            attemp,
            LAPSE_RATE,
            xc2,
        );
        (xh2[2] - horhgt, true_alt)
    }
}

/// Meridian/anti-meridian transit. Port of `calc_mer_trans` (swecl.c:4688-4748).
#[allow(clippy::too_many_arguments)]
fn calc_mer_trans(
    eph: &Ephemeris,
    tjd_ut: f64,
    body: Body,
    starname: Option<&str>,
    epheflag: CalcFlags,
    rsmi: RiseSetFlags,
    geopos: [f64; 3],
    topo_config: &EphemerisConfig,
) -> Result<RiseSetResult, Error> {
    let iflag = (epheflag & MER_TRANS_IFLAG_KEEP_MASK) | CalcFlags::EQUATORIAL | CalcFlags::TOPOCTR;
    let is_fixstar = starname.is_some_and(|s| !s.is_empty());

    let (armc0, _) = eph.azalt_armc_eps(tjd_ut, geopos[0]);

    let x0 = if is_fixstar {
        let (_, result) = eph.fixstar2_ut(starname.unwrap(), tjd_ut, iflag)?;
        [result.data[0], result.data[1]]
    } else {
        let r = eph.calc_ut_with_config(tjd_ut, body, iflag, topo_config)?;
        [r.data[0], r.data[1]]
    };

    let mut x = x0;
    let mut t = tjd_ut;
    let mut arxc = if rsmi.contains(RiseSetFlags::ITRANSIT) {
        crate::math::normalize_degrees(armc0 + 180.0)
    } else {
        armc0
    };

    for i in 0..4 {
        let mut mdd = crate::math::normalize_degrees(x[0] - arxc);
        if i > 0 && mdd > 180.0 {
            mdd -= 360.0;
        }
        t += mdd / 361.0;

        let (armc, _) = eph.azalt_armc_eps(t, geopos[0]);
        arxc = if rsmi.contains(RiseSetFlags::ITRANSIT) {
            crate::math::normalize_degrees(armc + 180.0)
        } else {
            armc
        };

        if !is_fixstar {
            let r = eph.calc_ut_with_config(t, body, iflag, topo_config)?;
            x = [r.data[0], r.data[1]];
        }
    }

    Ok(RiseSetResult { time: t })
}
