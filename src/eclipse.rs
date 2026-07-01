//! Solar eclipse geometry and shared eclipse/occultation helpers.
//!
//! Port of `swecl.c`. `eclipse_where` and `calc_planet_star` are shared verbatim by the
//! lunar-eclipse and occultation modules (see `docs/c-ref-eclipse-solar.md` §0-3).

use crate::constants::{
    AUNIT, DEGTORAD, EARTH_OBLATENESS, PLANETARY_DIAMETERS, RADTODEG, REARTH, RMOON,
};
use crate::context::Ephemeris;
use crate::error::Error;
use crate::flags::{CalcFlags, EclipseFlags};
use crate::math::{cartesian_to_polar, normalize_degrees, polar_to_cartesian};
use crate::types::Body;

/// Geographic + shadow-cone geometry of a solar eclipse (or lunar occultation) at a given
/// geocentric instant. Mirrors C's `geopos[0..1]` + `dcore[0..6]` outputs of `eclipse_where`
/// (swecl.c:640-886, §3.6). `geopos[2..9]` are documented in C as "not implemented so far" and
/// `dcore[7..9]` are always zero there; both are omitted from this struct.
#[derive(Debug, Clone, Copy)]
pub struct EclipseWhere {
    /// Geographic longitude of the point of greatest eclipse, degrees east positive.
    pub central_longitude: f64,
    /// Geographic latitude of the point of greatest eclipse, degrees north positive.
    pub central_latitude: f64,
    /// Core (umbra) shadow diameter at the point of maximum eclipse, km. Signed: positive means
    /// annular (the antumbra, not the umbra, reaches the ground), negative/zero means total.
    pub core_diameter_km: f64,
    /// Penumbra diameter at the point of maximum eclipse, km.
    pub penumbra_diameter_km: f64,
    /// Distance of the shadow axis from the geocenter, km.
    pub shadow_axis_distance_km: f64,
    /// Umbra (core shadow) diameter on the fundamental plane, km.
    pub umbra_diameter_fundamental_km: f64,
    /// Penumbra diameter on the fundamental plane, km.
    pub penumbra_diameter_fundamental_km: f64,
    /// Cosine of the umbra cone's half-angle.
    pub cos_umbra_half_angle: f64,
    /// Cosine of the penumbra cone's half-angle.
    pub cos_penumbra_half_angle: f64,
    /// Eclipse-type classification (CENTRAL/NONCENTRAL/TOTAL/ANNULAR/PARTIAL); `0` (empty) means
    /// no eclipse anywhere on Earth at this instant.
    pub flags: EclipseFlags,
}

/// Physical radius of `ipl` (or a star, if `starname` is given) in AU. Shared `drad` resolution
/// pattern used by `eclipse_where` (swecl.c:697-704) and `eclipse_how` (swecl.c:1004-1011).
pub(crate) fn body_radius_au(ipl: Body, starname: Option<&str>) -> f64 {
    if starname.is_some_and(|s| !s.is_empty()) {
        return 0.0;
    }
    let raw = ipl.to_raw_id();
    if (0..PLANETARY_DIAMETERS.len() as i32).contains(&raw) {
        PLANETARY_DIAMETERS[raw as usize] / 2.0 / AUNIT
    } else {
        // Named-asteroid diameter (C: swed.ast_diam, populated by the SE1 orbital-element
        // loader) isn't threaded through a stateless config yet -- out of scope until asteroid
        // orbital-element support lands.
        0.0
    }
}

/// Shared body/star position dispatch (`calc_planet_star`, swecl.c:888-897). Every
/// eclipse/occultation function that needs "the position of the eclipsed/occulted body" routes
/// through this. `tjd_et` is always Ephemeris/Dynamical Time.
pub(crate) fn calc_planet_star(
    eph: &Ephemeris,
    tjd_et: f64,
    ipl: Body,
    starname: Option<&str>,
    flags: CalcFlags,
) -> Result<[f64; 6], Error> {
    match starname {
        Some(name) if !name.is_empty() => {
            let (_, result) = eph.fixstar2(name, tjd_et, flags)?;
            Ok(result.data)
        }
        _ => Ok(eph.calc(tjd_et, ipl, flags)?.data),
    }
}

/// Shadow-cone geometry core (`eclipse_where`, swecl.c:640-886). Computes, for the geocentric
/// instant `tjd_ut` (UT), the shadow-cone geometry of the Moon with respect to the
/// eclipsed/occulted body `ipl`/`starname`, and — if the shadow axis touches the Earth — the
/// geographic position of the point of greatest eclipse.
///
/// Shared between solar eclipses (`ipl=Body::Sun`, `starname=None`) and lunar occultations of
/// any planet/asteroid/star ("Sun" below always means "the eclipsed/occulted body").
pub(crate) fn eclipse_where(
    eph: &Ephemeris,
    tjd_ut: f64,
    ipl: Body,
    starname: Option<&str>,
    ifl: CalcFlags,
) -> Result<EclipseWhere, Error> {
    let config = eph.config();
    let tjd = tjd_ut + crate::deltat::calc_deltat(tjd_ut, config);

    // `iflag2` (polar, radians) is derived from the pre-XYZ flavor of `iflag`; `iflag` is then
    // reassigned to the cartesian flavor (swecl.c:667-669) -- named distinctly here to avoid
    // replicating that in-place reassignment hazard.
    let iflag_polar_deg = CalcFlags::SPEED | CalcFlags::EQUATORIAL | ifl;
    let iflag_polar_rad = iflag_polar_deg | CalcFlags::RADIANS;
    let iflag_cart = iflag_polar_deg | CalcFlags::XYZ;

    let lm = eph.calc(tjd, Body::Moon, iflag_polar_rad)?.data;
    let ls = calc_planet_star(eph, tjd, ipl, starname, iflag_polar_rad)?;
    let rmt = {
        let d = eph.calc(tjd, Body::Moon, iflag_cart)?.data;
        [d[0], d[1], d[2]]
    };
    let rst = {
        let d = calc_planet_star(eph, tjd, ipl, starname, iflag_cart)?;
        [d[0], d[1], d[2]]
    };

    let sidt = if ifl.contains(CalcFlags::NONUT) {
        let eps_mean = crate::obliquity::obliquity(tjd, ifl, &config.astro_models).eps * RADTODEG;
        crate::sidereal_time::sidereal_time0(tjd_ut, eps_mean, 0.0, config) * 15.0 * DEGTORAD
    } else {
        crate::sidereal_time::sidereal_time(tjd_ut, config) * 15.0 * DEGTORAD
    };

    let drad = body_radius_au(ipl, starname);
    let rmoon = RMOON;
    let dmoon = 2.0 * rmoon;
    let de = REARTH;

    let mut earthobl = 1.0 - EARTH_OBLATENESS;
    let mut dcore = [0.0f64; 10];
    let mut retc = EclipseFlags::empty();
    let mut no_eclipse = false;
    let mut geopos = [0.0f64; 2];
    let mut dsmt = 0.0f64;
    let mut xst_cart = [0.0f64; 3];

    // Earth-oblateness substitution (swecl.c:705-787, label `iter_where`): executes this whole
    // block exactly twice, refining `earthobl` from an ellipsoid-normal factor at the latitude
    // found on pass 1. Only pass 2's `geopos`/`xst_cart` are used; C reaches them by falling
    // through the `goto` only when `niter > 0`.
    for niter in 0..2 {
        let mut rm = polar_to_cartesian([lm[0], lm[1], lm[2]]);
        rm[2] /= earthobl;
        let dm = (rm[0] * rm[0] + rm[1] * rm[1] + rm[2] * rm[2]).sqrt();

        let mut rs = polar_to_cartesian([ls[0], ls[1], ls[2]]);
        rs[2] /= earthobl;

        let mut e = [rm[0] - rs[0], rm[1] - rs[1], rm[2] - rs[2]];
        let et = [rmt[0] - rst[0], rmt[1] - rst[1], rmt[2] - rst[2]];
        let dsm = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
        dsmt = (et[0] * et[0] + et[1] * et[1] + et[2] * et[2]).sqrt();
        for v in &mut e {
            *v /= dsm;
        }

        let sinf1 = (drad - rmoon) / dsm;
        let cosf1 = (1.0 - sinf1 * sinf1).sqrt();
        let sinf2 = (drad + rmoon) / dsm;
        let cosf2 = (1.0 - sinf2 * sinf2).sqrt();

        // Distance of moon from fundamental plane / shadow axis from geocenter / shadow
        // diameters on the fundamental plane.
        let s0 = -(rm[0] * e[0] + rm[1] * e[1] + rm[2] * e[2]);
        let r0 = (dm * dm - s0 * s0).sqrt();
        let d0 = (s0 / dsm * (drad * 2.0 - dmoon) - dmoon) / cosf1;
        let cap_d0 = (s0 / dsm * (drad * 2.0 + dmoon) + dmoon) / cosf2;

        dcore[2] = r0 * (AUNIT / 1000.0);
        dcore[3] = d0 * (AUNIT / 1000.0);
        dcore[4] = cap_d0 * (AUNIT / 1000.0);
        dcore[5] = cosf1;
        dcore[6] = cosf2;

        retc = EclipseFlags::empty();
        no_eclipse = false;
        if de * cosf1 >= r0 {
            retc |= EclipseFlags::CENTRAL;
        } else if r0 <= de * cosf1 + d0.abs() / 2.0 {
            retc |= EclipseFlags::NONCENTRAL;
        } else if r0 <= de * cosf2 + cap_d0 / 2.0 {
            retc |= EclipseFlags::PARTIAL | EclipseFlags::NONCENTRAL;
        } else {
            no_eclipse = true;
        }

        // Distance of shadow point from fundamental plane / moon -> shadow point on earth.
        // Computed regardless of `no_eclipse` -- C does not early-return here (swecl.c:780).
        let d_fp = s0 * s0 + de * de - dm * dm;
        let d_fp = if d_fp > 0.0 { d_fp.sqrt() } else { 0.0 };
        let s = s0 - d_fp;

        let xs = [rm[0] + s * e[0], rm[1] + s * e[1], rm[2] + s * e[2]];
        let mut xst = xs;
        xst[2] *= earthobl;
        let xst_polar = cartesian_to_polar(xst);

        if niter == 0 {
            let cosfi = xst_polar[1].cos();
            let sinfi = xst_polar[1].sin();
            let eobl = EARTH_OBLATENESS;
            let cc = 1.0 / (cosfi * cosfi + (1.0 - eobl) * (1.0 - eobl) * sinfi * sinfi).sqrt();
            earthobl = (1.0 - eobl) * (1.0 - eobl) * cc;
            continue;
        }

        xst_cart = polar_to_cartesian(xst_polar);
        let mut xs_polar = cartesian_to_polar(xs);
        xs_polar[0] -= sidt;
        let mut lon_deg = normalize_degrees(xs_polar[0] * RADTODEG);
        if lon_deg > 180.0 {
            lon_deg -= 360.0;
        }
        geopos = [lon_deg, xs_polar[1] * RADTODEG];
    }

    // Core-shadow diameter at the point of maximum eclipse (swecl.c:865-875), using the raw
    // (non-oblateness-adjusted) saved Moon vector and pass-2's true-z shadow point.
    let x = [
        rmt[0] - xst_cart[0],
        rmt[1] - xst_cart[1],
        rmt[2] - xst_cart[2],
    ];
    let s = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let cosf1 = dcore[5];
    let cosf2 = dcore[6];
    dcore[0] = (s / dsmt * (drad * 2.0 - dmoon) - dmoon) * cosf1 * (AUNIT / 1000.0);
    dcore[1] = (s / dsmt * (drad * 2.0 + dmoon) + dmoon) * cosf2 * (AUNIT / 1000.0);

    if !retc.contains(EclipseFlags::PARTIAL) && !no_eclipse {
        if dcore[0] > 0.0 {
            retc |= EclipseFlags::ANNULAR;
        } else {
            retc |= EclipseFlags::TOTAL;
        }
    }

    Ok(EclipseWhere {
        central_longitude: geopos[0],
        central_latitude: geopos[1],
        core_diameter_km: dcore[0],
        penumbra_diameter_km: dcore[1],
        shadow_axis_distance_km: dcore[2],
        umbra_diameter_fundamental_km: dcore[3],
        penumbra_diameter_fundamental_km: dcore[4],
        cos_umbra_half_angle: dcore[5],
        cos_penumbra_half_angle: dcore[6],
        flags: retc,
    })
}

/// Public wrapper pinning `ipl=Sun`, `starname=None` (`swe_sol_eclipse_where`, swecl.c:565-582).
/// Local-circumstance attributes (`attr[]`, via `eclipse_how`) are added by a later task
/// (RSE 6, swisseph-rs/73); this returns the shadow-geometry (`geopos`/`dcore`) portion only,
/// which is exactly what C's `retflag` (this function's return value) reflects too -- C's
/// `eclipse_how` call only fills `attr[]`, it never changes the returned classification bitmask.
///
/// C masks `ifl &= SEFLG_EPHMASK` before calling `eclipse_where` (swecl.c:568) -- this strips
/// `NONUT`/`TOPOCTR`/etc. before they ever reach the geometry core, so e.g. `eclipse_where`'s own
/// `NONUT` branch is unreachable through this entry point. Replicated here rather than passing
/// `ifl` straight through.
pub(crate) fn sol_eclipse_where(
    eph: &Ephemeris,
    tjd_ut: f64,
    ifl: CalcFlags,
) -> Result<EclipseWhere, Error> {
    eclipse_where(eph, tjd_ut, Body::Sun, None, ifl & crate::calc::EPHMASK)
}
