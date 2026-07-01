//! Solar eclipse geometry and shared eclipse/occultation helpers.
//!
//! Port of `swecl.c`. `eclipse_where`, `eclipse_how`, and `calc_planet_star` are shared verbatim
//! by the lunar-eclipse and occultation modules (see `docs/c-ref-eclipse-solar.md` §0-4).

use std::f64::consts::PI;

use crate::constants::{
    AUNIT, DEGTORAD, EARTH_OBLATENESS, LAPSE_RATE, PLANETARY_DIAMETERS, RADTODEG, REARTH, RMOON,
};
use crate::context::{Ephemeris, EphemerisConfig, TopoPosition};
use crate::error::Error;
use crate::flags::{CalcFlags, EclipseFlags};
use crate::math::{cartesian_to_polar, dot_prod_unit, normalize_degrees, polar_to_cartesian};
use crate::types::Body;

/// Saros cycle length, days (swecl.c:114).
const SAROS_CYCLE: f64 = 6585.3213;

/// `(series_no, tstart)` pairs for the 181 solar-eclipse Saros series (swecl.c:107-298), derived
/// from NASA's eclipse Saros catalogue (<https://eclipse.gsfc.nasa.gov/SEsaros/SEsaros0-180.html>).
/// `tstart` is the JD (UT) of each series' initial eclipse.
#[rustfmt::skip]
const SAROS_DATA_SOLAR: [(i32, f64); 181] = [
    (0, 641886.5), (1, 672214.5), (2, 676200.5), (3, 693357.5), (4, 723685.5),
    (5, 727671.5), (6, 744829.5), (7, 775157.5), (8, 779143.5), (9, 783131.5),
    (10, 820044.5), (11, 810859.5), (12, 748993.5), (13, 792492.5), (14, 789892.5),
    (15, 787294.5), (16, 824207.5), (17, 834779.5), (18, 838766.5), (19, 869094.5),
    (20, 886251.5), (21, 890238.5), (22, 927151.5), (23, 937722.5), (24, 941709.5),
    (25, 978623.5), (26, 989194.5), (27, 993181.5), (28, 1023510.5), (29, 1034081.5),
    (30, 972214.5), (31, 1061811.5), (32, 1006529.5), (33, 997345.5), (34, 1021088.5),
    (35, 1038245.5), (36, 1042231.5), (37, 1065974.5), (38, 1089716.5), (39, 1093703.5),
    (40, 1117446.5), (41, 1141188.5), (42, 1145175.5), (43, 1168918.5), (44, 1192660.5),
    (45, 1196647.5), (46, 1220390.5), (47, 1244132.5), (48, 1234948.5), (49, 1265277.5),
    (50, 1282433.5), (51, 1207395.5), (52, 1217968.5), (53, 1254881.5), (54, 1252282.5),
    (55, 1262855.5), (56, 1293182.5), (57, 1297169.5), (58, 1314326.5), (59, 1344654.5),
    (60, 1348640.5), (61, 1365798.5), (62, 1396126.5), (63, 1400112.5), (64, 1417270.5),
    (65, 1447598.5), (66, 1444999.5), (67, 1462157.5), (68, 1492485.5), (69, 1456959.5),
    (70, 1421434.5), (71, 1471518.5), (72, 1455748.5), (73, 1466320.5), (74, 1496648.5),
    (75, 1500634.5), (76, 1511207.5), (77, 1548120.5), (78, 1552106.5), (79, 1562679.5),
    (80, 1599592.5), (81, 1603578.5), (82, 1614150.5), (83, 1644479.5), (84, 1655050.5),
    (85, 1659037.5), (86, 1695950.5), (87, 1693351.5), (88, 1631484.5), (89, 1727666.5),
    (90, 1672384.5), (91, 1663200.5), (92, 1693529.5), (93, 1710685.5), (94, 1714672.5),
    (95, 1738415.5), (96, 1755572.5), (97, 1766144.5), (98, 1789887.5), (99, 1807044.5),
    (100, 1817616.5), (101, 1841359.5), (102, 1858516.5), (103, 1862502.5), (104, 1892831.5),
    (105, 1903402.5), (106, 1887633.5), (107, 1924547.5), (108, 1921948.5), (109, 1873251.5),
    (110, 1890409.5), (111, 1914151.5), (112, 1918138.5), (113, 1935296.5), (114, 1959038.5),
    (115, 1963024.5), (116, 1986767.5), (117, 2010510.5), (118, 2014496.5), (119, 2031654.5),
    (120, 2061982.5), (121, 2065968.5), (122, 2083126.5), (123, 2113454.5), (124, 2104269.5),
    (125, 2108256.5), (126, 2151755.5), (127, 2083302.5), (128, 2080704.5), (129, 2124203.5),
    (130, 2121603.5), (131, 2132176.5), (132, 2162504.5), (133, 2166490.5), (134, 2177062.5),
    (135, 2207390.5), (136, 2217962.5), (137, 2228534.5), (138, 2258862.5), (139, 2269434.5),
    (140, 2273421.5), (141, 2310334.5), (142, 2314320.5), (143, 2311722.5), (144, 2355221.5),
    (145, 2319695.5), (146, 2284169.5), (147, 2314498.5), (148, 2325069.5), (149, 2329056.5),
    (150, 2352799.5), (151, 2369956.5), (152, 2380528.5), (153, 2404271.5), (154, 2421428.5),
    (155, 2425414.5), (156, 2455743.5), (157, 2472900.5), (158, 2476886.5), (159, 2500629.5),
    (160, 2517786.5), (161, 2515187.5), (162, 2545516.5), (163, 2556087.5), (164, 2487635.5),
    (165, 2504793.5), (166, 2535121.5), (167, 2525936.5), (168, 2543094.5), (169, 2573422.5),
    (170, 2577408.5), (171, 2594566.5), (172, 2624894.5), (173, 2628880.5), (174, 2646038.5),
    (175, 2669780.5), (176, 2673766.5), (177, 2690924.5), (178, 2721252.5), (179, 2718653.5),
    (180, 2729226.5),
];

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
/// Local-circumstance attributes (`attr[]`) live in [`EclipseHow`]/`sol_eclipse_how` instead;
/// this returns the shadow-geometry (`geopos`/`dcore`) portion only, which is exactly what C's
/// `retflag` (this function's return value) reflects too -- C's `eclipse_how` call only fills
/// `attr[]`, it never changes the returned classification bitmask.
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

/// Local circumstances of a solar eclipse (or lunar occultation) as seen from a specific
/// observer. Mirrors C's `attr[0..10]` output of `eclipse_how` (swecl.c:967-1152, §4.9).
/// `core_diameter_km` (`attr[3]`) is left at `0.0` by `eclipse_how` itself -- callers combining
/// this with `eclipse_where`'s geocentric shadow geometry fill it from `EclipseWhere::
/// core_diameter_km` (`dcore[0]`), same as C's `swe_sol_eclipse_how` (§4.11 step 6).
#[derive(Debug, Clone, Copy)]
pub struct EclipseHow {
    /// Magnitude: fraction of the eclipsed body's diameter covered by the Moon (IMCCE convention).
    pub magnitude: f64,
    /// Ratio of the Moon's angular diameter to the eclipsed body's angular diameter.
    pub diameter_ratio: f64,
    /// Obscuration: fraction of the eclipsed body's disc area covered by the Moon.
    pub obscuration: f64,
    /// Core (umbra) shadow diameter, km -- `0.0` unless filled by the caller from
    /// [`EclipseWhere::core_diameter_km`].
    pub core_diameter_km: f64,
    /// Azimuth of the eclipsed body, degrees, measured from south, clockwise via west.
    pub azimuth: f64,
    /// True (geometric) altitude of the eclipsed body above the horizon, degrees.
    pub true_altitude: f64,
    /// Apparent (refraction-corrected) altitude of the eclipsed body above the horizon, degrees.
    pub apparent_altitude: f64,
    /// Angular separation ("elongation") of the Moon from the eclipsed body's center, degrees.
    pub elongation: f64,
    /// Magnitude per the NASA convention (`= magnitude` for partial; `= diameter_ratio` for
    /// total/annular).
    pub nasa_magnitude: f64,
    /// Saros series number (solar eclipses of the Sun only; `-99999999.0` if none found).
    pub saros_series: f64,
    /// Saros series member number, 1-based (solar eclipses of the Sun only).
    pub saros_member: f64,
    /// Eclipse-type classification (TOTAL/ANNULAR/PARTIAL, possibly OR'd with VISIBLE); empty
    /// means no eclipse visible from this location at this instant.
    pub flags: EclipseFlags,
}

/// [`calc_planet_star`] with a topocentric config override threaded through the planet branch
/// (shared by [`eclipse_how`]). Stars don't yet have a per-call topographic override (`fixstar2`
/// has no TOPOCTR path, matching riseset.rs's fixed-star note) -- only the planet branch uses
/// `config`.
fn calc_planet_star_topo(
    eph: &Ephemeris,
    config: &EphemerisConfig,
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
        _ => Ok(eph.calc_with_config(tjd_et, ipl, flags, config)?.data),
    }
}

/// Local circumstances at an observer (`eclipse_how`, swecl.c:967-1152). `geolon`/`geolat`
/// degrees, `geohgt` meters above sea. C relies on `swe_set_topo` populating global `swed.topd`
/// before its `SEFLG_TOPOCTR` calc calls; the stateless port threads a per-call config override
/// instead (mirrors `Ephemeris::calc_with_config`'s doc comment / riseset.rs's identical
/// pattern), never mutating `eph`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn eclipse_how(
    eph: &Ephemeris,
    tjd_ut: f64,
    ipl: Body,
    starname: Option<&str>,
    ifl: CalcFlags,
    geolon: f64,
    geolat: f64,
    geohgt: f64,
) -> Result<EclipseHow, Error> {
    let config = eph.config();
    let te = tjd_ut + crate::deltat::calc_deltat(tjd_ut, config);

    let topo_config = {
        let mut c = config.clone();
        c.topographic = Some(TopoPosition {
            longitude: geolon,
            latitude: geolat,
            altitude: geohgt,
        });
        c
    };

    let iflag = CalcFlags::EQUATORIAL | CalcFlags::TOPOCTR | ifl;
    let iflag_cart = iflag | CalcFlags::XYZ;

    let ls = calc_planet_star_topo(eph, &topo_config, te, ipl, starname, iflag)?;
    let lm = eph
        .calc_with_config(te, Body::Moon, iflag, &topo_config)?
        .data;
    let xs = calc_planet_star_topo(eph, &topo_config, te, ipl, starname, iflag_cart)?;
    let xm = eph
        .calc_with_config(te, Body::Moon, iflag_cart, &topo_config)?
        .data;

    let drad = body_radius_au(ipl, starname);

    let geopos = [geolon, geolat, geohgt];
    let xh = eph.azalt(
        tjd_ut,
        crate::azalt::AzAltDir::EquToHor,
        geopos,
        0.0,
        10.0,
        LAPSE_RATE,
        [ls[0], ls[1]],
    );

    let rmoon = (RMOON / lm[2]).asin() * RADTODEG;
    let rsun = (drad / ls[2]).asin() * RADTODEG;
    let rsplusrm = rsun + rmoon;
    let rsminusrm = rsun - rmoon;
    let x1 = [xs[0] / ls[2], xs[1] / ls[2], xs[2] / ls[2]];
    let x2 = [xm[0] / lm[2], xm[1] / lm[2], xm[2] / lm[2]];
    let dctr = dot_prod_unit(x1, x2).acos() * RADTODEG;

    let mut retc = EclipseFlags::empty();
    if dctr < rsminusrm {
        retc = EclipseFlags::ANNULAR;
    } else if dctr < rsminusrm.abs() {
        retc = EclipseFlags::TOTAL;
    } else if dctr < rsplusrm {
        retc = EclipseFlags::PARTIAL;
    }

    let diameter_ratio = if rsun > 0.0 { rmoon / rsun } else { 0.0 };

    let lsunleft = -dctr + rsun + rmoon;
    let magnitude = if rsun > 0.0 {
        lsunleft / rsun / 2.0
    } else {
        1.0
    };

    // Obscuration: fraction of the eclipsed body's disc area covered by the Moon
    // (circular-segment lens-area formula, swecl.c:1075-1107, §4.6).
    let lsun = rsun;
    let lmoon = rmoon;
    let lctr = dctr;
    let obscuration = if retc.is_empty() || lsun == 0.0 {
        1.0
    } else if retc == EclipseFlags::TOTAL || retc == EclipseFlags::ANNULAR {
        lmoon * lmoon / lsun / lsun
    } else {
        let a_denom = 2.0 * lctr * lmoon;
        let b_denom = 2.0 * lctr * lsun;
        if a_denom < 1e-9 {
            lmoon * lmoon / lsun / lsun
        } else {
            let a = ((lctr * lctr + lmoon * lmoon - lsun * lsun) / a_denom).clamp(-1.0, 1.0);
            let b = ((lctr * lctr + lsun * lsun - lmoon * lmoon) / b_denom).clamp(-1.0, 1.0);
            let a = a.acos();
            let b = b.acos();
            let mut sc1 = a * lmoon * lmoon / 2.0;
            let mut sc2 = b * lsun * lsun / 2.0;
            sc1 -= (a.cos() * a.sin()) * lmoon * lmoon / 2.0;
            sc2 -= (b.cos() * b.sin()) * lsun * lsun / 2.0;
            (sc1 + sc2) * 2.0 / PI / lsun / lsun
        }
    };

    // Visibility threshold (swecl.c:1108-1123): 34.4556' horizon refraction (Bennett) +
    // 1.75'/sqrt(h) horizon dip + 0.37'/sqrt(h) observer-to-horizon refraction.
    let hmin_appr = -(34.4556 + (1.75 + 0.37) * geohgt.sqrt()) / 60.0;
    if xh[1] + rsun + hmin_appr.abs() >= 0.0 && !retc.is_empty() {
        retc |= EclipseFlags::VISIBLE;
    }

    // NASA magnitude + Saros series/member: only for a genuine solar eclipse (Sun, not a star).
    let mut nasa_magnitude = 0.0;
    let mut saros_series = 0.0;
    let mut saros_member = 0.0;
    if ipl == Body::Sun && starname.unwrap_or("").is_empty() {
        nasa_magnitude = if retc.intersects(EclipseFlags::TOTAL | EclipseFlags::ANNULAR) {
            diameter_ratio
        } else {
            magnitude
        };

        let mut found = false;
        for &(series_no, tstart) in SAROS_DATA_SOLAR.iter() {
            let mut d = (tjd_ut - tstart) / SAROS_CYCLE;
            if d < 0.0 && d * SAROS_CYCLE > -2.0 {
                d = 0.0000001;
            }
            if d < 0.0 {
                continue;
            }
            let j = d as i32;
            if (d - j as f64) * SAROS_CYCLE < 2.0 {
                saros_series = series_no as f64;
                saros_member = (j + 1) as f64;
                found = true;
                break;
            }
            let k = j + 1;
            if (k as f64 - d) * SAROS_CYCLE < 2.0 {
                saros_series = series_no as f64;
                saros_member = (k + 1) as f64;
                found = true;
                break;
            }
        }
        if !found {
            saros_series = -99999999.0;
            saros_member = -99999999.0;
        }
    }

    Ok(EclipseHow {
        magnitude,
        diameter_ratio,
        obscuration,
        core_diameter_km: 0.0,
        azimuth: xh[0],
        true_altitude: xh[1],
        apparent_altitude: xh[2],
        elongation: dctr,
        nasa_magnitude,
        saros_series,
        saros_member,
        flags: retc,
    })
}

/// Public wrapper pinning `ipl=Sun`, `starname=None` (`swe_sol_eclipse_how`, swecl.c:922-964).
/// Layers a horizon-visibility gate on top of [`eclipse_how`]'s own geometric result: a second,
/// independent topocentric az/alt pass (matching C's redundant `swe_calc_ut` + `swe_azalt`) can
/// zero out the returned classification -- and with it `magnitude`/`diameter_ratio`/
/// `obscuration`/`core_diameter_km`/`nasa_magnitude`/`saros_series`/`saros_member` -- purely
/// because the Sun's apparent altitude is `<= 0`, even when `eclipse_how` found a real eclipse in
/// progress. `azimuth`/`true_altitude`/`apparent_altitude`/`elongation` stay populated either way
/// (§4.11 step 8).
pub(crate) fn sol_eclipse_how(
    eph: &Ephemeris,
    tjd_ut: f64,
    ifl: CalcFlags,
    geopos: [f64; 3],
) -> Result<EclipseHow, Error> {
    if !(crate::constants::RISE_SET_GEOALT_MIN..=crate::constants::RISE_SET_GEOALT_MAX)
        .contains(&geopos[2])
    {
        return Err(Error::CError(format!(
            "location for eclipses must be between {:.0} and {:.0} m above sea",
            crate::constants::RISE_SET_GEOALT_MIN,
            crate::constants::RISE_SET_GEOALT_MAX
        )));
    }
    let ifl = ifl & crate::calc::EPHMASK;

    let mut how = eclipse_how(
        eph,
        tjd_ut,
        Body::Sun,
        None,
        ifl,
        geopos[0],
        geopos[1],
        geopos[2],
    )?;
    let mut retflag = how.flags;

    let where_result = eclipse_where(eph, tjd_ut, Body::Sun, None, ifl)?;
    if !retflag.is_empty() {
        retflag |= where_result.flags & (EclipseFlags::CENTRAL | EclipseFlags::NONCENTRAL);
    }
    how.core_diameter_km = where_result.core_diameter_km;

    let topo_config = {
        let mut c = eph.config().clone();
        c.topographic = Some(TopoPosition {
            longitude: geopos[0],
            latitude: geopos[1],
            altitude: geopos[2],
        });
        c
    };
    let ls = eph
        .calc_ut_with_config(
            tjd_ut,
            Body::Sun,
            ifl | CalcFlags::TOPOCTR | CalcFlags::EQUATORIAL,
            &topo_config,
        )?
        .data;
    let xaz = eph.azalt(
        tjd_ut,
        crate::azalt::AzAltDir::EquToHor,
        geopos,
        0.0,
        10.0,
        LAPSE_RATE,
        [ls[0], ls[1]],
    );
    how.azimuth = xaz[0];
    how.true_altitude = xaz[1];
    how.apparent_altitude = xaz[2];

    if xaz[2] <= 0.0 {
        retflag = EclipseFlags::empty();
    }
    if retflag.is_empty() {
        how.magnitude = 0.0;
        how.diameter_ratio = 0.0;
        how.obscuration = 0.0;
        how.core_diameter_km = 0.0;
        how.nasa_magnitude = 0.0;
        how.saros_series = 0.0;
        how.saros_member = 0.0;
    }
    how.flags = retflag;

    Ok(how)
}
