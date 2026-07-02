//! Solar eclipse geometry and shared eclipse/occultation helpers.
//!
//! Port of `swecl.c`. `eclipse_where`, `eclipse_how`, and `calc_planet_star` are shared verbatim
//! by the lunar-eclipse and occultation modules (see `docs/c-ref-eclipse-solar.md` §0-4).

use std::f64::consts::PI;

use crate::config::{EphemerisConfig, TopoPosition};
use crate::constants::{
    AUNIT, DEARTH, DEGTORAD, DSUN, EARTH_OBLATENESS, J2000, LAPSE_RATE, PLANETARY_DIAMETERS,
    RADTODEG, REARTH, RMOON, RSUN,
};
use crate::context::Ephemeris;
use crate::error::Error;
use crate::flags::{CalcFlags, EclipseFlags, RiseSetFlags};
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

/// `(series_no, tstart)` pairs for the 180 lunar-eclipse Saros series (swecl.c:306-486), parallel
/// to but distinct from [`SAROS_DATA_SOLAR`]. `tstart` is the JD (UT) of each series' initial
/// eclipse.
#[rustfmt::skip]
const SAROS_DATA_LUNAR: [(i32, f64); 180] = [
    (1, 782437.5), (2, 799593.5), (3, 783824.5), (4, 754884.5), (5, 824724.5),
    (6, 762857.5), (7, 773430.5), (8, 810343.5), (9, 807743.5), (10, 824901.5),
    (11, 855229.5), (12, 859215.5), (13, 876373.5), (14, 906701.5), (15, 910687.5),
    (16, 927845.5), (17, 958173.5), (18, 962159.5), (19, 979317.5), (20, 1009645.5),
    (21, 1007046.5), (22, 1017618.5), (23, 1054531.5), (24, 979493.5), (25, 976895.5),
    (26, 1020394.5), (27, 1017794.5), (28, 1028367.5), (29, 1058695.5), (30, 1062681.5),
    (31, 1073253.5), (32, 1110167.5), (33, 1114153.5), (34, 1131311.5), (35, 1161639.5),
    (36, 1165625.5), (37, 1176197.5), (38, 1213111.5), (39, 1217097.5), (40, 1221084.5),
    (41, 1257997.5), (42, 1255398.5), (43, 1186946.5), (44, 1283128.5), (45, 1227845.5),
    (46, 1225247.5), (47, 1255575.5), (48, 1272732.5), (49, 1276719.5), (50, 1307047.5),
    (51, 1317619.5), (52, 1328191.5), (53, 1358519.5), (54, 1375676.5), (55, 1379663.5),
    (56, 1409991.5), (57, 1420562.5), (58, 1424549.5), (59, 1461463.5), (60, 1465449.5),
    (61, 1436509.5), (62, 1493179.5), (63, 1457653.5), (64, 1435298.5), (65, 1452456.5),
    (66, 1476198.5), (67, 1480184.5), (68, 1503928.5), (69, 1527670.5), (70, 1531656.5),
    (71, 1548814.5), (72, 1579142.5), (73, 1583128.5), (74, 1600286.5), (75, 1624028.5),
    (76, 1628015.5), (77, 1651758.5), (78, 1675500.5), (79, 1672901.5), (80, 1683474.5),
    (81, 1713801.5), (82, 1645349.5), (83, 1649336.5), (84, 1686249.5), (85, 1683650.5),
    (86, 1694222.5), (87, 1731136.5), (88, 1735122.5), (89, 1745694.5), (90, 1776022.5),
    (91, 1786594.5), (92, 1797166.5), (93, 1827494.5), (94, 1838066.5), (95, 1848638.5),
    (96, 1878966.5), (97, 1882952.5), (98, 1880354.5), (99, 1923853.5), (100, 1881741.5),
    (101, 1852801.5), (102, 1889715.5), (103, 1893701.5), (104, 1897688.5), (105, 1928016.5),
    (106, 1938588.5), (107, 1942575.5), (108, 1972903.5), (109, 1990059.5), (110, 1994046.5),
    (111, 2024375.5), (112, 2034946.5), (113, 2045518.5), (114, 2075847.5), (115, 2086418.5),
    (116, 2083820.5), (117, 2120733.5), (118, 2124719.5), (119, 2062852.5), (120, 2086596.5),
    (121, 2103752.5), (122, 2094568.5), (123, 2118311.5), (124, 2142054.5), (125, 2146040.5),
    (126, 2169783.5), (127, 2186940.5), (128, 2197512.5), (129, 2214670.5), (130, 2238412.5),
    (131, 2242398.5), (132, 2266142.5), (133, 2289884.5), (134, 2287285.5), (135, 2311028.5),
    (136, 2334770.5), (137, 2292659.5), (138, 2276890.5), (139, 2326974.5), (140, 2304619.5),
    (141, 2308606.5), (142, 2345520.5), (143, 2349506.5), (144, 2360078.5), (145, 2390406.5),
    (146, 2394392.5), (147, 2411550.5), (148, 2441878.5), (149, 2445864.5), (150, 2456437.5),
    (151, 2486765.5), (152, 2490751.5), (153, 2501323.5), (154, 2538236.5), (155, 2529052.5),
    (156, 2473771.5), (157, 2563367.5), (158, 2508085.5), (159, 2505486.5), (160, 2542400.5),
    (161, 2546386.5), (162, 2556958.5), (163, 2587287.5), (164, 2597858.5), (165, 2601845.5),
    (166, 2632173.5), (167, 2649330.5), (168, 2653317.5), (169, 2683645.5), (170, 2694217.5),
    (171, 2698203.5), (172, 2728532.5), (173, 2739103.5), (174, 2683822.5), (175, 2740492.5),
    (176, 2724722.5), (177, 2708952.5), (178, 2732695.5), (179, 2749852.5), (180, 2753839.5),
];

/// Saros series/member lookup shared by the solar (`eclipse_how`) and lunar (`lun_eclipse_how`)
/// eclipse cores -- identical scan algorithm over two distinct Saros tables
/// ([`SAROS_DATA_SOLAR`]/[`SAROS_DATA_LUNAR`]), swecl.c:1112-1144 / 3352-3372. Returns
/// `(series_no, member_no)`, or `(-99999999.0, -99999999.0)` if no series matches.
fn saros_lookup(tjd_ut: f64, table: &[(i32, f64)]) -> (f64, f64) {
    for &(series_no, tstart) in table {
        let mut d = (tjd_ut - tstart) / SAROS_CYCLE;
        if d < 0.0 && d * SAROS_CYCLE > -2.0 {
            d = 0.0000001;
        }
        if d < 0.0 {
            continue;
        }
        let j = d as i32;
        if (d - j as f64) * SAROS_CYCLE < 2.0 {
            return (series_no as f64, (j + 1) as f64);
        }
        let k = j + 1;
        if (k as f64 - d) * SAROS_CYCLE < 2.0 {
            return (series_no as f64, (k + 1) as f64);
        }
    }
    (-99999999.0, -99999999.0)
}

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
    /// means no eclipse visible from this location at this instant. [`sol_eclipse_how`]
    /// additionally merges CENTRAL/NONCENTRAL in from the geocentric shadow geometry
    /// (`eclipse_where`) -- `eclipse_how` itself never sets those bits.
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

        let (s, m) = saros_lookup(tjd_ut, &SAROS_DATA_SOLAR);
        saros_series = s;
        saros_member = m;
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

/// Global solar-eclipse search result: `tret[0..10]` per `swe_sol_eclipse_when_glob`
/// (swecl.c:1185-1515, §5.8). `tret[8]`/`tret[9]` (annular-total transition times) are not
/// implemented upstream and always `0.0` in C -- omitted here.
#[derive(Debug, Clone, Copy)]
pub struct SolarEclipseGlobal {
    /// Time (UT) of maximum eclipse: geocentric minimum Sun-Moon angular separation. `tret[0]`.
    pub time_maximum: f64,
    /// Time (UT) when the eclipse's RA-conjunction instant occurs (geocentric RA(Sun) ==
    /// RA(Moon)), or `0.0` if no such instant falls within the eclipse window. `tret[1]`.
    pub time_ra_conjunction: f64,
    /// Time (UT) of eclipse begin, first contact anywhere on Earth. `tret[2]`.
    pub time_begin: f64,
    /// Time (UT) of eclipse end, last contact anywhere on Earth. `tret[3]`.
    pub time_end: f64,
    /// Time (UT) of totality/annularity begin, `0.0` if partial. `tret[4]`.
    pub time_totality_begin: f64,
    /// Time (UT) of totality/annularity end, `0.0` if partial. `tret[5]`.
    pub time_totality_end: f64,
    /// Time (UT) of center-line begin, `0.0` if noncentral. `tret[6]`.
    pub time_centerline_begin: f64,
    /// Time (UT) of center-line end, `0.0` if noncentral. `tret[7]`.
    pub time_centerline_end: f64,
    /// Eclipse-type classification (CENTRAL/NONCENTRAL combined with
    /// TOTAL/ANNULAR/HYBRID/PARTIAL). Never empty -- the search retries indefinitely (bounded
    /// only by ephemeris range) until a matching eclipse is found.
    pub flags: EclipseFlags,
}

/// Meeus (German ed., p.381) synodic-month lunation estimate for signed month index `k`: the
/// approximate mean-conjunction instant (ET/TT, not yet delta-T-corrected to UT) plus the
/// F-argument node-proximity pre-filter. Textually identical in `swe_sol_eclipse_when_glob`
/// (swecl.c:1227-1265, §5.2) and `eclipse_when_loc` (swecl.c:2129-2172, §6.2) -- each C function
/// re-derives it independently; factored into one function here per the C ref doc's porting note.
/// Returns `None` if `k`'s F-argument (Moon's distance from its node at the mean conjunction)
/// falls in the `(21,159)` degree band, meaning no solar eclipse is geometrically possible for
/// this lunation -- the caller should advance `k` and retry.
fn meeus_new_moon_estimate(k: f64) -> Option<f64> {
    let tt_ = k / 1236.85;
    let t2 = tt_ * tt_;
    let t3 = t2 * tt_;
    let t4 = t3 * tt_;
    let mut ff = normalize_degrees(
        160.7108 + 390.67050274 * k - 0.0016341 * t2 - 0.00000227 * t3 + 0.000000011 * t4,
    );
    if ff > 180.0 {
        ff -= 180.0;
    }
    if ff > 21.0 && ff < 159.0 {
        return None;
    }

    // Approximate time of geocentric maximum eclipse (Meeus, German ed., p. 381).
    let mut tjd =
        2451550.09765 + 29.530588853 * k + 0.0001337 * t2 - 0.000000150 * t3 + 0.00000000073 * t4;
    let m = normalize_degrees(2.5534 + 29.10535669 * k - 0.0000218 * t2 - 0.00000011 * t3);
    let mm = normalize_degrees(
        201.5643 + 385.81693528 * k + 0.1017438 * t2 + 0.00001239 * t3 + 0.000000058 * t4,
    );
    let e = 1.0 - 0.002516 * tt_ - 0.0000074 * t2;
    let m_rad = m * DEGTORAD;
    let mm_rad = mm * DEGTORAD;
    tjd = tjd - 0.4075 * mm_rad.sin() + 0.1721 * e * m_rad.sin();
    Some(tjd)
}

/// Contact-time sample formula (swecl.c:1394-1424, §5.5), shared by the coarse `find_zero` pass
/// and the 3-pass Newton refinement. `n`: 0 = eclipse begin/end (penumbra boundary, but divides
/// by `cosf1`/umbra half-angle -- literal C quirk, not a typo, negligible impact since both
/// half-angles are well under 1°); 1 = totality/annularity begin/end (umbra boundary); 2 =
/// center-line begin/end.
fn contact_dc(n: u32, w: &EclipseWhere, de_km: f64) -> f64 {
    match n {
        0 => {
            w.penumbra_diameter_fundamental_km / 2.0 + de_km / w.cos_umbra_half_angle
                - w.shadow_axis_distance_km
        }
        1 => {
            w.umbra_diameter_fundamental_km.abs() / 2.0 + de_km / w.cos_penumbra_half_angle
                - w.shadow_axis_distance_km
        }
        _ => de_km / w.cos_penumbra_half_angle - w.shadow_axis_distance_km,
    }
}

/// Global eclipse search: find the next (or, if `backward`, previous) solar eclipse anywhere on
/// Earth after/before `tjd_start` (UT), restricted to eclipse types in `ifltype`. Port of
/// `swe_sol_eclipse_when_glob` (swecl.c:1185-1515, §5). `ifltype = EclipseFlags::empty()` means
/// all types.
pub(crate) fn sol_eclipse_when_glob(
    eph: &Ephemeris,
    tjd_start: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    backward: bool,
) -> Result<SolarEclipseGlobal, Error> {
    let ifl = ifl & crate::calc::EPHMASK;
    let config = eph.config();

    if ifltype == (EclipseFlags::PARTIAL | EclipseFlags::CENTRAL) {
        return Err(Error::CError(
            "central partial eclipses do not exist".to_string(),
        ));
    }
    if ifltype == (EclipseFlags::HYBRID | EclipseFlags::NONCENTRAL) {
        return Err(Error::CError(
            "non-central hybrid (annular-total) eclipses do not exist".to_string(),
        ));
    }

    let mut ifltype = ifltype;
    if ifltype.is_empty() {
        ifltype = EclipseFlags::ALLTYPES_SOLAR;
    }
    if ifltype == EclipseFlags::TOTAL
        || ifltype == EclipseFlags::ANNULAR
        || ifltype == EclipseFlags::HYBRID
    {
        ifltype |= EclipseFlags::NONCENTRAL | EclipseFlags::CENTRAL;
    }
    if ifltype == EclipseFlags::PARTIAL {
        ifltype |= EclipseFlags::NONCENTRAL;
    }

    let direction = if backward { -1.0 } else { 1.0 };
    let iflag = CalcFlags::EQUATORIAL | ifl;
    let iflag_cart = iflag | CalcFlags::XYZ;
    let de_km = 6378.140;

    let mut k = ((tjd_start - J2000) / 365.2425 * 12.3685).trunc();
    k -= direction;

    'next_try: loop {
        let mut tret = [0.0f64; 8];

        let mut tjd = match meeus_new_moon_estimate(k) {
            Some(tjd) => tjd,
            None => {
                k += direction;
                continue 'next_try;
            }
        };

        // Iterative refinement to the instant of minimum geocentric Sun-Moon angular separation
        // (§5.3). `tjd` is treated as ET/TT throughout this refinement (Meeus's formula is
        // dynamical time); UT conversion happens once, after convergence.
        let dtstart = if !(2_000_000.0..=2_500_000.0).contains(&tjd) {
            5.0
        } else {
            1.0
        };
        let mut dt = dtstart;
        while dt > 0.0001 {
            let mut dc = [0.0f64; 3];
            let mut t = tjd - dt;
            for dc_i in dc.iter_mut() {
                let ls = eph.calc(t, Body::Sun, iflag)?.data;
                let lm = eph.calc(t, Body::Moon, iflag)?.data;
                let xs = eph.calc(t, Body::Sun, iflag_cart)?.data;
                let xm = eph.calc(t, Body::Moon, iflag_cart)?.data;
                let xa = [xs[0] / ls[2], xs[1] / ls[2], xs[2] / ls[2]];
                let xb = [xm[0] / lm[2], xm[1] / lm[2], xm[2] / lm[2]];
                let rmoon = (RMOON / lm[2]).asin() * RADTODEG;
                let rsun = (RSUN / ls[2]).asin() * RADTODEG;
                *dc_i = dot_prod_unit(xa, xb).acos() * RADTODEG - (rmoon + rsun);
                t += dt;
            }
            let (dtint, _) = crate::math::find_maximum(dc[0], dc[1], dc[2], dt);
            tjd += dtint + dt;
            dt /= 4.0;
        }

        // 3-pass fixed-point ET->UT conversion (swecl.c:1310-1312).
        let tjds1 = tjd - crate::deltat::calc_deltat(tjd, config);
        let tjds2 = tjd - crate::deltat::calc_deltat(tjds1, config);
        let tjd = tjd - crate::deltat::calc_deltat(tjds2, config);

        let where_result = eclipse_where(eph, tjd, Body::Sun, None, ifl)?;
        // In extreme cases `eclipse_where` under-detects a tiny eclipse -- confirm via
        // `eclipse_how` with the coordinates `eclipse_where` returned.
        let how_result = eclipse_how(
            eph,
            tjd,
            Body::Sun,
            None,
            ifl,
            where_result.central_longitude,
            where_result.central_latitude,
            0.0,
        )?;
        if how_result.flags.is_empty() {
            k += direction;
            continue 'next_try;
        }
        tret[0] = tjd;
        if (backward && tret[0] >= tjd_start - 0.0001)
            || (!backward && tret[0] <= tjd_start + 0.0001)
        {
            k += direction;
            continue 'next_try;
        }

        // Eclipse type (TOTAL/ANNULAR/etc.); ANNULAR_TOTAL (hybrid) is discovered later (§5.6).
        let mut retflag = where_result.flags;
        let mut dont_times = false;
        if retflag.is_empty() {
            // Can happen with an extremely small percentage (C: "fix this????").
            retflag = EclipseFlags::PARTIAL | EclipseFlags::NONCENTRAL;
            tret[4] = tjd;
            tret[5] = tjd;
            dont_times = true;
        }

        if !ifltype.contains(EclipseFlags::NONCENTRAL) && retflag.contains(EclipseFlags::NONCENTRAL)
        {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::CENTRAL) && retflag.contains(EclipseFlags::CENTRAL) {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::ANNULAR) && retflag.contains(EclipseFlags::ANNULAR) {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::PARTIAL) && retflag.contains(EclipseFlags::PARTIAL) {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.intersects(EclipseFlags::TOTAL | EclipseFlags::HYBRID)
            && retflag.contains(EclipseFlags::TOTAL)
        {
            k += direction;
            continue 'next_try;
        }

        if dont_times {
            return Ok(SolarEclipseGlobal {
                time_maximum: tret[0],
                time_ra_conjunction: 0.0,
                time_begin: 0.0,
                time_end: 0.0,
                time_totality_begin: tret[4],
                time_totality_end: tret[5],
                time_centerline_begin: 0.0,
                time_centerline_end: 0.0,
                flags: retflag,
            });
        }

        // Contact-time refinement (§5.5): n=0 eclipse begin/end (always), n=1 totality/
        // annularity begin/end (skip if PARTIAL), n=2 center-line begin/end (skip if
        // NONCENTRAL).
        let o = if retflag.contains(EclipseFlags::PARTIAL) {
            0
        } else if retflag.contains(EclipseFlags::NONCENTRAL) {
            1
        } else {
            2
        };
        let dta = 2.0 / 24.0;
        let dtb = 10.0 / 24.0 / 60.0 / 3.0;

        for n in 0..=o {
            let (i1, i2) = match n {
                0 => (2usize, 3usize),
                1 => (4usize, 5usize),
                _ => (6usize, 7usize),
            };

            let mut dc = [0.0f64; 3];
            let mut t = tjd - dta;
            for dc_i in dc.iter_mut() {
                let w = eclipse_where(eph, t, Body::Sun, None, ifl)?;
                *dc_i = contact_dc(n, &w, de_km);
                t += dta;
            }
            // Divergence from C on `find_zero` failure (no real parabola roots): C ignores the
            // failure return and proceeds with zero-initialized `dt1`/`dt2`, i.e. assigns
            // `tret[i1] = tret[i2] = tjd + dta` and Newton-refines from there. We leave the
            // slots at 0.0 and skip refinement instead -- refining around 0.0 would evaluate
            // `eclipse_where` at JD ~0 (4713 BC). Unreachable for well-conditioned real
            // eclipses: the confirmed eclipse guarantees a sign change in `dc`.
            if let Some((dt1, dt2)) = crate::math::find_zero(dc[0], dc[1], dc[2], dta) {
                tret[i1] = tjd + dt1 + dta;
                tret[i2] = tjd + dt2 + dta;

                let mut dt = dtb;
                for _ in 0..3 {
                    for &j in &[i1, i2] {
                        let mut dc2 = [0.0f64; 2];
                        let mut t = tret[j] - dt;
                        for dc_i in dc2.iter_mut() {
                            let w = eclipse_where(eph, t, Body::Sun, None, ifl)?;
                            *dc_i = contact_dc(n, &w, de_km);
                            t += dt;
                        }
                        let dt1 = dc2[1] / ((dc2[1] - dc2[0]) / dt);
                        tret[j] -= dt1;
                    }
                    dt /= 3.0;
                }
            }
        }

        // Annular-total (hybrid) detection (§5.6).
        if retflag.contains(EclipseFlags::TOTAL) {
            let dc0 = eclipse_where(eph, tret[0], Body::Sun, None, ifl)?.core_diameter_km;
            let dc1 = eclipse_where(eph, tret[4], Body::Sun, None, ifl)?.core_diameter_km;
            let dc2 = eclipse_where(eph, tret[5], Body::Sun, None, ifl)?.core_diameter_km;
            if dc0 * dc1 < 0.0 || dc0 * dc2 < 0.0 {
                retflag |= EclipseFlags::HYBRID;
                retflag.remove(EclipseFlags::TOTAL);
            }
        }
        if !ifltype.contains(EclipseFlags::TOTAL) && retflag.contains(EclipseFlags::TOTAL) {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::HYBRID) && retflag.contains(EclipseFlags::HYBRID) {
            k += direction;
            continue 'next_try;
        }

        // Time of maximum eclipse at local apparent noon (§5.7): first check for a solar
        // RA-transit sign change between eclipse begin/end; if found, secant-iterate to the
        // instant of exact geocentric RA(Sun) == RA(Moon).
        let mut dc_transit = [0.0f64; 2];
        for (i, dc_i) in dc_transit.iter_mut().enumerate() {
            let tt = tret[2 + i] + crate::deltat::calc_deltat(tret[2 + i], config);
            let ls = eph.calc(tt, Body::Sun, iflag)?.data;
            let lm = eph.calc(tt, Body::Moon, iflag)?.data;
            let mut d = normalize_degrees(ls[0] - lm[0]);
            if d > 180.0 {
                d -= 360.0;
            }
            *dc_i = d;
        }
        if dc_transit[0] * dc_transit[1] >= 0.0 {
            tret[1] = 0.0;
        } else {
            let mut tjd_ra = tjd;
            let mut dt = 0.1;
            let dt1_init = (tret[3] - tret[2]) / 2.0;
            if dt1_init < dt {
                dt = dt1_init / 2.0;
            }
            while dt > 0.01 {
                let mut dc2 = [0.0f64; 2];
                let mut t = tjd_ra;
                for dc_i in dc2.iter_mut() {
                    let tt = t + crate::deltat::calc_deltat(t, config);
                    let ls = eph.calc(tt, Body::Sun, iflag)?.data;
                    let lm = eph.calc(tt, Body::Moon, iflag)?.data;
                    let mut d = normalize_degrees(ls[0] - lm[0]);
                    if d > 180.0 {
                        d -= 360.0;
                    }
                    if d > 180.0 {
                        d -= 360.0;
                    }
                    *dc_i = d;
                    t -= dt;
                }
                let a = (dc2[1] - dc2[0]) / dt;
                if a < 1e-10 {
                    break;
                }
                let dt1 = dc2[0] / a;
                tjd_ra += dt1;
                dt /= 3.0;
            }
            tret[1] = tjd_ra;
        }

        return Ok(SolarEclipseGlobal {
            time_maximum: tret[0],
            time_ra_conjunction: tret[1],
            time_begin: tret[2],
            time_end: tret[3],
            time_totality_begin: tret[4],
            time_totality_end: tret[5],
            time_centerline_begin: tret[6],
            time_centerline_end: tret[7],
            flags: retflag,
        });
    }
}

/// Local solar-eclipse search result: `tret[0..7]` per `swe_sol_eclipse_when_loc` +
/// `eclipse_when_loc` (swecl.c:2019-2410, §6). **Index semantics differ from
/// [`SolarEclipseGlobal`]'s `tret[]`** (§6.3) -- do not conflate the two.
#[derive(Debug, Clone, Copy)]
pub struct SolarEclipseLocal {
    /// Time (UT) of maximum eclipse as seen from this location -- re-anchored to sunrise/sunset
    /// if the true geocentric maximum wasn't visible here. `tret[0]`.
    pub time_maximum: f64,
    /// Time (UT) of first contact (penumbra ingress, i.e. visible eclipse begin here). `tret[1]`.
    pub time_first_contact: f64,
    /// Time (UT) of second contact (umbra/antumbra ingress -- totality/annularity begin here),
    /// `0.0` if only a partial eclipse is visible from this location. `tret[2]`.
    pub time_second_contact: f64,
    /// Time (UT) of third contact (umbra/antumbra egress -- totality/annularity end here), `0.0`
    /// if only a partial eclipse is visible from this location. `tret[3]`.
    pub time_third_contact: f64,
    /// Time (UT) of fourth contact (penumbra egress, i.e. visible eclipse end here). `tret[4]`.
    pub time_fourth_contact: f64,
    /// Time (UT) of sunrise between first and fourth contact, `0.0` if none (or circumpolar).
    /// `tret[5]`.
    pub time_sunrise: f64,
    /// Time (UT) of sunset between first and fourth contact, `0.0` if none (or circumpolar).
    /// `tret[6]`.
    pub time_sunset: f64,
    /// Local circumstances (`attr[]`) at whichever instant was written last: the moment of
    /// maximum eclipse, unless a sunrise/sunset re-anchor (below) overwrote it. `core_diameter_km`
    /// is filled in by [`sol_eclipse_when_loc`] from a geocentric `eclipse_where` call at
    /// `time_maximum` (§6.1 step 4), same "geocentric, not observer-specific" caveat as
    /// `sol_eclipse_how`.
    pub attr: EclipseHow,
    /// Eclipse-type classification (TOTAL/ANNULAR/PARTIAL) OR'd with VISIBLE and whichever of
    /// MAX/1ST/2ND/3RD/4TH_VISIBLE applied at some contact; `NONCENTRAL` merged in by
    /// [`sol_eclipse_when_loc`]. The search loop retries internally until a visible eclipse is
    /// found -- never empty.
    pub flags: EclipseFlags,
}

/// Overlap-gap sample used throughout [`eclipse_when_loc`]'s convergence/contact refinement:
/// `acos(dot_prod_unit(xs/|xs|, xm/|xm|)) * RADTODEG`, i.e. the angular separation between the
/// (already-topocentric) cartesian Sun/Moon position vectors `xs`/`xm`. Distances are computed
/// manually via `sqrt(sum of squares)` rather than reusing a polar-array distance component,
/// matching the C source's own `/*ls[2]*/`/`/*lm[2]*/` comments flagging this as deliberate
/// (§6.2.1) -- port literally, do not "simplify" to a polar re-fetch.
fn topo_angular_separation(xs: [f64; 3], xm: [f64; 3]) -> f64 {
    let ds = (xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2]).sqrt();
    let dm = (xm[0] * xm[0] + xm[1] * xm[1] + xm[2] * xm[2]).sqrt();
    let x1 = [xs[0] / ds, xs[1] / ds, xs[2] / ds];
    let x2 = [xm[0] / dm, xm[1] / dm, xm[2] / dm];
    dot_prod_unit(x1, x2).acos() * RADTODEG
}

/// Local eclipse search: find the next (or, if `backward`, previous) solar eclipse **visible
/// from** `geopos` (topocentric, unlike [`sol_eclipse_when_glob`]'s purely geocentric search).
/// Port of `eclipse_when_loc` (swecl.c:2100-2410, §6.2). `geopos` = [longitude, latitude, height
/// above sea (m)], degrees/degrees/meters.
#[allow(clippy::too_many_arguments)]
pub(crate) fn eclipse_when_loc(
    eph: &Ephemeris,
    tjd_start: f64,
    ifl: CalcFlags,
    geopos: [f64; 3],
    backward: bool,
) -> Result<SolarEclipseLocal, Error> {
    let config = eph.config();
    let direction = if backward { -1.0 } else { 1.0 };

    let topo_config = {
        let mut c = config.clone();
        c.topographic = Some(TopoPosition {
            longitude: geopos[0],
            latitude: geopos[1],
            altitude: geopos[2],
        });
        c
    };

    // Computed once, reused unchanged for the whole function (main loop, contacts 2/3, contacts
    // 1/4) except the visibility-scan/sunrise-sunset `eclipse_how` calls, which build their own
    // topocentric flags internally from raw `ifl` (§6.2, swecl.c:2126-2127).
    let iflag = CalcFlags::EQUATORIAL | CalcFlags::TOPOCTR | ifl;
    let iflag_cart = iflag | CalcFlags::XYZ;
    let iflag_cart_speed = iflag_cart | CalcFlags::SPEED;

    let mut k = ((tjd_start - J2000) / 365.2425 * 12.3685).trunc();
    k -= direction;

    'next_try: loop {
        let mut tjd = match meeus_new_moon_estimate(k) {
            Some(tjd) => tjd,
            None => {
                k += direction;
                continue 'next_try;
            }
        };

        // §6.2.1 main convergence loop: refine `tjd` (still ET/TT) to the instant of minimum
        // topocentric Sun-Moon angular separation. `dtstart`/`dt<=1e-5` boundary and the
        // 2-then-3 `dtdiv` schedule are intentionally different constants than §5.2's global
        // search -- do not unify.
        let mut dtdiv = 2.0;
        let dtstart = if (1_900_000.0..=2_500_000.0).contains(&tjd) {
            0.5
        } else {
            2.0
        };
        let mut dt = dtstart;
        while dt > 0.00001 {
            if dt < 0.1 {
                dtdiv = 3.0;
            }
            let mut dc = [0.0f64; 3];
            let mut t = tjd - dt;
            for dc_i in dc.iter_mut() {
                let xs = eph
                    .calc_with_config(t, Body::Sun, iflag_cart, &topo_config)?
                    .data;
                let xm = eph
                    .calc_with_config(t, Body::Moon, iflag_cart, &topo_config)?
                    .data;
                *dc_i = topo_angular_separation([xs[0], xs[1], xs[2]], [xm[0], xm[1], xm[2]]);
                t += dt;
            }
            let (dtint, _) = crate::math::find_maximum(dc[0], dc[1], dc[2], dt);
            tjd += dtint + dt;
            dt /= dtdiv;
        }

        // §6.2.2 post-convergence confirmation: a fresh set of calc calls at the converged
        // `tjd`, not a reuse of the loop's last sample.
        let xs = eph
            .calc_with_config(tjd, Body::Sun, iflag_cart, &topo_config)?
            .data;
        let ls = eph
            .calc_with_config(tjd, Body::Sun, iflag, &topo_config)?
            .data;
        let xm = eph
            .calc_with_config(tjd, Body::Moon, iflag_cart, &topo_config)?
            .data;
        let lm = eph
            .calc_with_config(tjd, Body::Moon, iflag, &topo_config)?
            .data;
        let dctr = dot_prod_unit([xs[0], xs[1], xs[2]], [xm[0], xm[1], xm[2]]).acos() * RADTODEG;
        let rmoon = (RMOON / lm[2]).asin() * RADTODEG;
        let rsun = (RSUN / ls[2]).asin() * RADTODEG;
        let rsplusrm = rsun + rmoon;
        let rsminusrm = rsun - rmoon;
        if dctr > rsplusrm {
            k += direction;
            continue 'next_try;
        }

        let mut tret = [0.0f64; 7];
        // 2-pass fixed-point ET->UT conversion (fewer passes than §5.3's 3-pass global-search
        // version).
        let t0_pass1 = tjd - crate::deltat::calc_deltat(tjd, config);
        tret[0] = tjd - crate::deltat::calc_deltat(t0_pass1, config);
        if (backward && tret[0] >= tjd_start - 0.0001)
            || (!backward && tret[0] <= tjd_start + 0.0001)
        {
            k += direction;
            continue 'next_try;
        }

        // Phase classification (swecl.c:2235-2240): mirrors §4.4's thresholds, but the partial
        // branch tests `dctr <= rsplusrm` (not strict `<` like `eclipse_how`) -- preserve the
        // exact operator (measure-zero difference given the rejection guard above).
        let mut retflag = if dctr < rsminusrm {
            EclipseFlags::ANNULAR
        } else if dctr < rsminusrm.abs() {
            EclipseFlags::TOTAL
        } else if dctr <= rsplusrm {
            EclipseFlags::PARTIAL
        } else {
            EclipseFlags::empty()
        };
        let dctrmin = dctr;

        // Contacts 2/3 (swecl.c:2242-2300): skipped (tret[2]=tret[3]=0) if only a partial
        // eclipse is visible from this location (umbra never reaches here).
        if dctrmin <= rsminusrm.abs() {
            let twomin = 2.0 / 24.0 / 60.0;
            let mut dc = [0.0f64; 3];
            dc[1] = rsminusrm.abs() - dctrmin;
            for &(i, t) in &[(0usize, tjd - twomin), (2usize, tjd + twomin)] {
                let xs = eph
                    .calc_with_config(t, Body::Sun, iflag_cart, &topo_config)?
                    .data;
                let xm = eph
                    .calc_with_config(t, Body::Moon, iflag_cart, &topo_config)?
                    .data;
                let ds = (xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2]).sqrt();
                let dm = (xm[0] * xm[0] + xm[1] * xm[1] + xm[2] * xm[2]).sqrt();
                let mut rmoon = (RMOON / dm).asin() * RADTODEG;
                rmoon *= 0.99916; // "gives better accuracy for 2nd/3rd contacts" (swecl.c)
                let rsun = (RSUN / ds).asin() * RADTODEG;
                let dctr_i = topo_angular_separation([xs[0], xs[1], xs[2]], [xm[0], xm[1], xm[2]]);
                dc[i] = (rsun - rmoon).abs() - dctr_i;
            }
            if let Some((dt1, dt2)) = crate::math::find_zero(dc[0], dc[1], dc[2], twomin) {
                tret[2] = tjd + dt1 + twomin;
                tret[3] = tjd + dt2 + twomin;
                let tensec = 10.0 / 24.0 / 60.0 / 60.0;
                let mut dt = tensec;
                for _ in 0..2 {
                    for &j in &[2usize, 3usize] {
                        let mut xs = eph
                            .calc_with_config(tret[j], Body::Sun, iflag_cart_speed, &topo_config)?
                            .data;
                        let mut xm = eph
                            .calc_with_config(tret[j], Body::Moon, iflag_cart_speed, &topo_config)?
                            .data;
                        let mut dc2 = [0.0f64; 2];
                        for (i2, dc2_i) in dc2.iter_mut().enumerate() {
                            if i2 == 1 {
                                for k in 0..3 {
                                    xs[k] -= xs[k + 3] * dt;
                                    xm[k] -= xm[k + 3] * dt;
                                }
                            }
                            let ds = (xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2]).sqrt();
                            let dm = (xm[0] * xm[0] + xm[1] * xm[1] + xm[2] * xm[2]).sqrt();
                            let mut rmoon = (RMOON / dm).asin() * RADTODEG;
                            rmoon *= 0.99916;
                            let rsun = (RSUN / ds).asin() * RADTODEG;
                            let dctr_i = topo_angular_separation(
                                [xs[0], xs[1], xs[2]],
                                [xm[0], xm[1], xm[2]],
                            );
                            *dc2_i = (rsun - rmoon).abs() - dctr_i;
                        }
                        let dt1 = -dc2[0] / ((dc2[0] - dc2[1]) / dt);
                        tret[j] += dt1;
                    }
                    dt /= 10.0;
                }
                tret[2] -= crate::deltat::calc_deltat(tret[2], config);
                tret[3] -= crate::deltat::calc_deltat(tret[3], config);
            }
        }

        // Contacts 1/4 (swecl.c:2301-2353): structurally identical to contacts 2/3 but always
        // computed (every eclipse, even partial, has a 1st/4th contact), no `0.99916`
        // correction, uses `rsplusrm` (outer penumbra boundary) instead of `fabs(rsminusrm)`,
        // and the secant-refinement `dc` formula takes `fabs(rsplusrm)` even though `rsplusrm`
        // is never negative -- an asymmetry vs. the initial-sample formula (no `fabs` there);
        // preserve both exactly.
        {
            let twohr = 2.0 / 24.0;
            let mut dc = [0.0f64; 3];
            dc[1] = rsplusrm - dctrmin;
            for &(i, t) in &[(0usize, tjd - twohr), (2usize, tjd + twohr)] {
                let xs = eph
                    .calc_with_config(t, Body::Sun, iflag_cart, &topo_config)?
                    .data;
                let xm = eph
                    .calc_with_config(t, Body::Moon, iflag_cart, &topo_config)?
                    .data;
                let ds = (xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2]).sqrt();
                let dm = (xm[0] * xm[0] + xm[1] * xm[1] + xm[2] * xm[2]).sqrt();
                let rmoon = (RMOON / dm).asin() * RADTODEG;
                let rsun = (RSUN / ds).asin() * RADTODEG;
                let dctr_i = topo_angular_separation([xs[0], xs[1], xs[2]], [xm[0], xm[1], xm[2]]);
                dc[i] = (rsun + rmoon) - dctr_i;
            }
            if let Some((dt1, dt2)) = crate::math::find_zero(dc[0], dc[1], dc[2], twohr) {
                tret[1] = tjd + dt1 + twohr;
                tret[4] = tjd + dt2 + twohr;
                let tenmin = 10.0 / 24.0 / 60.0;
                let mut dt = tenmin;
                for _ in 0..3 {
                    for &j in &[1usize, 4usize] {
                        let mut xs = eph
                            .calc_with_config(tret[j], Body::Sun, iflag_cart_speed, &topo_config)?
                            .data;
                        let mut xm = eph
                            .calc_with_config(tret[j], Body::Moon, iflag_cart_speed, &topo_config)?
                            .data;
                        let mut dc2 = [0.0f64; 2];
                        for (i2, dc2_i) in dc2.iter_mut().enumerate() {
                            if i2 == 1 {
                                for k in 0..3 {
                                    xs[k] -= xs[k + 3] * dt;
                                    xm[k] -= xm[k + 3] * dt;
                                }
                            }
                            let ds = (xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2]).sqrt();
                            let dm = (xm[0] * xm[0] + xm[1] * xm[1] + xm[2] * xm[2]).sqrt();
                            let rmoon = (RMOON / dm).asin() * RADTODEG;
                            let rsun = (RSUN / ds).asin() * RADTODEG;
                            let dctr_i = topo_angular_separation(
                                [xs[0], xs[1], xs[2]],
                                [xm[0], xm[1], xm[2]],
                            );
                            *dc2_i = (rsun + rmoon).abs() - dctr_i;
                        }
                        let dt1 = -dc2[0] / ((dc2[0] - dc2[1]) / dt);
                        tret[j] += dt1;
                    }
                    dt /= 10.0;
                }
                tret[1] -= crate::deltat::calc_deltat(tret[1], config);
                tret[4] -= crate::deltat::calc_deltat(tret[4], config);
            } else {
                // Divergence from C on `find_zero` failure: C proceeds with zero-initialized
                // offsets (`tret[1] = tret[4] = tjd + twohr`) and refines from there. Leaving
                // the slots at 0.0 is not an option here -- `tret[1]` feeds the
                // `rise_trans(tret[1] - 0.001, ...)` calls below, which would then evaluate at
                // JD ~0 (4713 BC). Unreachable for real eclipses (`dc[1] = rsplusrm - dctrmin
                // >= 0` by the rejection guard above, so the parabola has real roots); if it
                // ever fires, treat the lunation as "no eclipse found" and advance the search.
                k += direction;
                continue 'next_try;
            }
        }

        // Visibility scan (swecl.c:2354-2384): DESCENDING order so the i=0 (max) write survives
        // last in `how`, matching C's shared-`attr[]`-clobbering behavior. Every non-skipped
        // call in this scan uses the raw `ifl` (not `iflag`/`iflag_cart`), since `eclipse_how`
        // adds its own EQUATORIAL|TOPOCTR internally.
        let mut how: Option<EclipseHow> = None;
        for i in (0..=4).rev() {
            if tret[i] == 0.0 {
                continue;
            }
            let h = eclipse_how(
                eph,
                tret[i],
                Body::Sun,
                None,
                ifl,
                geopos[0],
                geopos[1],
                geopos[2],
            )?;
            if h.apparent_altitude > 0.0 {
                retflag |= EclipseFlags::VISIBLE;
                retflag |= match i {
                    0 => EclipseFlags::MAX_VISIBLE,
                    1 => EclipseFlags::PARTBEG_VISIBLE,
                    2 => EclipseFlags::TOTBEG_VISIBLE,
                    3 => EclipseFlags::TOTEND_VISIBLE,
                    4 => EclipseFlags::PARTEND_VISIBLE,
                    _ => unreachable!(),
                };
            }
            how = Some(h);
        }
        if !retflag.contains(EclipseFlags::VISIBLE) {
            k += direction;
            continue 'next_try;
        }
        let mut how = how.expect("tret[0] is always set, so the i=0 scan iteration always runs");

        // Sunrise/sunset interaction (swecl.c:2385-2420). Literal quirk: both `swe_rise_trans`
        // calls pass `iflag` (the TOPOCTR-augmented, function-top value), unlike every
        // `eclipse_how` call in this function which uses raw `ifl` -- an inconsistency in the
        // original C, preserved here rather than "normalized" to `ifl`.
        let rise = eph.rise_trans(
            tret[1] - 0.001,
            Body::Sun,
            None,
            iflag,
            RiseSetFlags::RISE | RiseSetFlags::DISC_BOTTOM,
            geopos,
            0.0,
            0.0,
        );
        let tjdr = match rise {
            Ok(r) => r.time,
            Err(Error::CircumpolarBody) => {
                return Ok(SolarEclipseLocal {
                    time_maximum: tret[0],
                    time_first_contact: tret[1],
                    time_second_contact: tret[2],
                    time_third_contact: tret[3],
                    time_fourth_contact: tret[4],
                    time_sunrise: tret[5],
                    time_sunset: tret[6],
                    attr: how,
                    flags: retflag,
                });
            }
            Err(e) => return Err(e),
        };
        let set = eph.rise_trans(
            tret[1] - 0.001,
            Body::Sun,
            None,
            iflag,
            RiseSetFlags::SET | RiseSetFlags::DISC_BOTTOM,
            geopos,
            0.0,
            0.0,
        );
        let tjds = match set {
            Ok(r) => r.time,
            Err(Error::CircumpolarBody) => {
                return Ok(SolarEclipseLocal {
                    time_maximum: tret[0],
                    time_first_contact: tret[1],
                    time_second_contact: tret[2],
                    time_third_contact: tret[3],
                    time_fourth_contact: tret[4],
                    time_sunrise: tret[5],
                    time_sunset: tret[6],
                    attr: how,
                    flags: retflag,
                });
            }
            Err(e) => return Err(e),
        };

        if tjds < tret[1] || (tjds > tjdr && tjdr > tret[4]) {
            k += direction;
            continue 'next_try;
        }
        if tjdr > tret[1] && tjdr < tret[4] {
            tret[5] = tjdr;
            if !retflag.contains(EclipseFlags::MAX_VISIBLE) {
                tret[0] = tjdr;
                let h = eclipse_how(
                    eph,
                    tret[5],
                    Body::Sun,
                    None,
                    ifl,
                    geopos[0],
                    geopos[1],
                    geopos[2],
                )?;
                retflag.remove(EclipseFlags::TOTAL | EclipseFlags::ANNULAR | EclipseFlags::PARTIAL);
                retflag |=
                    h.flags & (EclipseFlags::TOTAL | EclipseFlags::ANNULAR | EclipseFlags::PARTIAL);
                how = h;
            }
        }
        if tjds > tret[1] && tjds < tret[4] {
            tret[6] = tjds;
            if !retflag.contains(EclipseFlags::MAX_VISIBLE) {
                tret[0] = tjds;
                let h = eclipse_how(
                    eph,
                    tret[6],
                    Body::Sun,
                    None,
                    ifl,
                    geopos[0],
                    geopos[1],
                    geopos[2],
                )?;
                retflag.remove(EclipseFlags::TOTAL | EclipseFlags::ANNULAR | EclipseFlags::PARTIAL);
                retflag |=
                    h.flags & (EclipseFlags::TOTAL | EclipseFlags::ANNULAR | EclipseFlags::PARTIAL);
                how = h;
            }
        }

        return Ok(SolarEclipseLocal {
            time_maximum: tret[0],
            time_first_contact: tret[1],
            time_second_contact: tret[2],
            time_third_contact: tret[3],
            time_fourth_contact: tret[4],
            time_sunrise: tret[5],
            time_sunset: tret[6],
            attr: how,
            flags: retflag,
        });
    }
}

/// Public wrapper (`swe_sol_eclipse_when_loc`, swecl.c:2019-2041, §6.1). `geopos` = [longitude,
/// latitude, height above sea (m)], degrees/degrees/meters.
///
/// Unlike [`sol_eclipse_how`], only the `NONCENTRAL` bit (not `CENTRAL`) is merged in from the
/// geocentric `eclipse_where` call at `time_maximum` -- a genuine literal difference between the
/// two wrappers, preserved exactly (§6.1 step 4).
pub(crate) fn sol_eclipse_when_loc(
    eph: &Ephemeris,
    tjd_start: f64,
    ifl: CalcFlags,
    geopos: [f64; 3],
    backward: bool,
) -> Result<SolarEclipseLocal, Error> {
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

    let mut result = eclipse_when_loc(eph, tjd_start, ifl, geopos, backward)?;
    let where_result = eclipse_where(eph, result.time_maximum, Body::Sun, None, ifl)?;
    result.flags |= where_result.flags & EclipseFlags::NONCENTRAL;
    result.attr.core_diameter_km = where_result.core_diameter_km;

    Ok(result)
}

/// Result of `lun_eclipse_how`'s static core geometry pass (swecl.c:3248-3372,
/// docs/c-ref-eclipse-lunar.md §2): the `attr[]` subset it computes directly -- everything except
/// azimuth/altitude, which only the public wrapper adds. Also carries the `dcore[0..4]`
/// fundamental-plane shadow-cone geometry (`r0`/`d0`/`D0`/`cosf1`/`cosf2`, ref doc §2.5 table)
/// that `swe_lun_eclipse_when`'s contact-time refinement (ref doc §4.7) needs -- `swe_lun_eclipse_
/// how`'s public wrapper output doesn't carry these, only this internal core does.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LunarEclipseCore {
    /// `attr[0]`/`attr[8]` -- umbral magnitude: fraction of the Moon's diameter covered by the
    /// umbra. `0.0` for a penumbral-only eclipse or no eclipse.
    pub umbral_magnitude: f64,
    /// `attr[1]` -- penumbral magnitude, always computed regardless of eclipse type.
    pub penumbral_magnitude: f64,
    /// `attr[7]` -- Moon's distance from exact opposition, degrees; `0.0` unless an eclipse is
    /// occurring (`flags` non-empty).
    pub distance_from_opposition: f64,
    /// `attr[9]` -- Saros series number, or `-99999999.0` if none found.
    pub saros_series: f64,
    /// `attr[10]` -- Saros series member number (1-based), or `-99999999.0` if none found.
    pub saros_member: f64,
    /// `dcore[0]` -- distance of the shadow axis from the selenocenter, AU.
    pub r0: f64,
    /// `dcore[1]` -- diameter of the umbra (core shadow) on the fundamental plane, AU.
    pub d0: f64,
    /// `dcore[2]` -- diameter of the penumbra (half-shadow) on the fundamental plane, AU.
    pub cap_d0: f64,
    /// `dcore[3]` -- cosine of the umbra cone's half-angle.
    pub cosf1: f64,
    /// `dcore[4]` -- cosine of the penumbra cone's half-angle.
    pub cosf2: f64,
    /// Eclipse-type classification: empty (no eclipse), or exactly one of TOTAL/PARTIAL/
    /// PENUMBRAL (`retc`, §2.10).
    pub flags: EclipseFlags,
}

/// Selenocentric shadow-cone geometry core (`lun_eclipse_how`, swecl.c:3248-3372, ref doc §2).
/// Computes the Earth-shadow-cone geometry and eclipse magnitudes for one geocentric instant
/// `tjd_ut` (UT). This is the mirror image of [`eclipse_where`]'s geometry: here the Moon is the
/// "screen" and the Earth casts the shadow, whereas `eclipse_where` casts the Moon's shadow onto
/// the Earth.
pub(crate) fn lun_eclipse_how(
    eph: &Ephemeris,
    tjd_ut: f64,
    ifl: CalcFlags,
) -> Result<LunarEclipseCore, Error> {
    let config = eph.config();
    let tjd = tjd_ut + crate::deltat::calc_deltat(tjd_ut, config);

    let iflag = CalcFlags::SPEED | CalcFlags::EQUATORIAL | ifl | CalcFlags::XYZ;
    let rm0 = {
        let d = eph.calc(tjd, Body::Moon, iflag)?.data;
        [d[0], d[1], d[2]]
    };
    let dm = (rm0[0] * rm0[0] + rm0[1] * rm0[1] + rm0[2] * rm0[2]).sqrt();
    let rs0 = {
        let d = eph.calc(tjd, Body::Sun, iflag)?.data;
        [d[0], d[1], d[2]]
    };
    let ds = (rs0[0] * rs0[0] + rs0[1] * rs0[1] + rs0[2] * rs0[2]).sqrt();

    let x1 = [rs0[0] / ds, rs0[1] / ds, rs0[2] / ds];
    let x2 = [rm0[0] / dm, rm0[1] / dm, rm0[2] / dm];
    let dctr = dot_prod_unit(x1, x2).acos() * RADTODEG;

    // Change of origin: selenocentric frame (§2.2).
    let rs = [rs0[0] - rm0[0], rs0[1] - rm0[1], rs0[2] - rm0[2]];
    let rm = [-rm0[0], -rm0[1], -rm0[2]];

    let e_raw = [rm[0] - rs[0], rm[1] - rs[1], rm[2] - rs[2]];
    let dsm = (e_raw[0] * e_raw[0] + e_raw[1] * e_raw[1] + e_raw[2] * e_raw[2]).sqrt();
    let e = [e_raw[0] / dsm, e_raw[1] / dsm, e_raw[2] / dsm];

    let f1 = (RSUN - REARTH) / dsm;
    let cosf1 = (1.0 - f1 * f1).sqrt();
    let f2 = (RSUN + REARTH) / dsm;
    let cosf2 = (1.0 - f2 * f2).sqrt();

    // Position of the Moon relative to the shadow axis (§2.4). `dm` here is the Moon's
    // geocentric distance from before the frame change (magnitude unchanged by negation) --
    // reused directly rather than recomputed, matching the C source's own reuse.
    let s0 = -(rm[0] * e[0] + rm[1] * e[1] + rm[2] * e[2]);
    let r0 = (dm * dm - s0 * s0).sqrt();

    // Shadow diameters on the fundamental plane, with atmospheric enlargement (§2.5). Ordering
    // and the doubled cosf1/cosf2 division are literal C, not simplified -- see ref doc §2.5.
    let mut d0 = (s0 / dsm * (DSUN - DEARTH) - DEARTH).abs() * (1.0 + 1.0 / 50.0) / cosf1;
    let mut cap_d0 = (s0 / dsm * (DSUN + DEARTH) + DEARTH) * (1.0 + 1.0 / 50.0) / cosf2;
    d0 /= cosf1;
    cap_d0 /= cosf2;
    d0 *= 0.99405;
    cap_d0 *= 0.98813;

    let rmoon = RMOON;
    let dmoon = 2.0 * rmoon;

    // Phase / umbral magnitude (§2.6).
    let mut flags = EclipseFlags::empty();
    let mut umbral_magnitude = 0.0;
    if d0 / 2.0 >= r0 + rmoon / cosf1 {
        flags = EclipseFlags::TOTAL;
        umbral_magnitude = (d0 / 2.0 - r0 + rmoon) / dmoon;
    } else if d0 / 2.0 >= r0 - rmoon / cosf1 {
        flags = EclipseFlags::PARTIAL;
        umbral_magnitude = (d0 / 2.0 - r0 + rmoon) / dmoon;
    } else if cap_d0 / 2.0 >= r0 - rmoon / cosf2 {
        flags = EclipseFlags::PENUMBRAL;
    }

    // Penumbral magnitude, always computed regardless of `flags` (§2.7).
    let penumbral_magnitude = (cap_d0 / 2.0 - r0 + rmoon) / dmoon;

    // Distance from opposition (§2.8).
    let distance_from_opposition = if !flags.is_empty() {
        180.0 - dctr.abs()
    } else {
        0.0
    };

    // Saros series lookup (§2.9).
    let (saros_series, saros_member) = saros_lookup(tjd_ut, &SAROS_DATA_LUNAR);

    Ok(LunarEclipseCore {
        umbral_magnitude,
        penumbral_magnitude,
        distance_from_opposition,
        saros_series,
        saros_member,
        r0,
        d0,
        cap_d0,
        cosf1,
        cosf2,
        flags,
    })
}

/// Contact-time sample formula shared by `swe_lun_eclipse_when`'s coarse `find_zero` bracket and
/// its 3-round Newton refinement (ref doc §4.7, swecl.c:3583-3608). `n`: 0 = penumbral begin/end
/// (`dcore[2]`/`D0` boundary), 1 = partial (umbra) begin/end (`dcore[1]`/`d0` boundary, `+
/// RMOON/cosf1`), 2 = totality begin/end (same `d0` boundary, `- RMOON/cosf1` -- far limb exits
/// the umbra entirely).
fn lun_contact_dc(n: u32, core: &LunarEclipseCore) -> f64 {
    match n {
        0 => core.cap_d0 / 2.0 + RMOON / core.cosf2 - core.r0,
        1 => core.d0 / 2.0 + RMOON / core.cosf1 - core.r0,
        _ => core.d0 / 2.0 - RMOON / core.cosf1 - core.r0,
    }
}

/// Local circumstances of a lunar eclipse at a given instant, plus the Moon's azimuth/altitude
/// at an observer. Mirrors C's `attr[0..10]` output of `swe_lun_eclipse_how`/`lun_eclipse_how`
/// (swecl.c:3172-3239, ref doc §1/§3). `attr[2]`/`attr[3]` are solar-only slots (vestigial
/// index-parity padding in C), omitted here; `attr[8]` (documented as a duplicate of `attr[0]`)
/// is likewise omitted -- callers needing exact `attr[]` index parity should reuse
/// `umbral_magnitude` for both `attr[0]` and `attr[8]`.
#[derive(Debug, Clone, Copy)]
pub struct LunarEclipseHow {
    /// `attr[0]`/`attr[8]` -- umbral magnitude.
    pub umbral_magnitude: f64,
    /// `attr[1]` -- penumbral magnitude.
    pub penumbral_magnitude: f64,
    /// `attr[4]` -- azimuth of the Moon, degrees, measured from south, clockwise via west.
    pub azimuth: f64,
    /// `attr[5]` -- true (geometric) altitude of the Moon above the horizon, degrees.
    pub true_altitude: f64,
    /// `attr[6]` -- apparent (refraction-corrected) altitude of the Moon above the horizon,
    /// degrees.
    pub apparent_altitude: f64,
    /// `attr[7]` -- Moon's distance from exact opposition, degrees.
    pub distance_from_opposition: f64,
    /// `attr[9]` -- Saros series number, or `-99999999.0` if none found.
    pub saros_series: f64,
    /// `attr[10]` -- Saros series member number, 1-based, or `-99999999.0` if none found.
    pub saros_member: f64,
    /// Eclipse-type classification: empty (no eclipse, or Moon below the horizon at `geopos`),
    /// or exactly one of TOTAL/PARTIAL/PENUMBRAL.
    pub flags: EclipseFlags,
}

/// Public wrapper (`swe_lun_eclipse_how`, swecl.c:3190-3239, ref doc §3): geocentric eclipse
/// attributes at `tjd_ut` (UT) plus the Moon's topocentric azimuth/altitude at `geopos`. Unlike C
/// (where `geopos` may be `NULL` for a geocentric-only query), this port always requires
/// `geopos` -- callers needing the geocentric-only core directly should call [`lun_eclipse_how`]
/// (used internally by the eclipse-search module, a later task).
///
/// If the Moon's apparent altitude at `geopos` is `<= 0` (below the horizon), `flags` is forced
/// empty even if the geocentric geometry found a real eclipse in progress --
/// `umbral_magnitude`/`penumbral_magnitude`/`saros_series`/`saros_member` stay populated from the
/// geocentric calculation regardless (§3 step 6, matching C exactly).
pub(crate) fn swe_lun_eclipse_how(
    eph: &Ephemeris,
    tjd_ut: f64,
    ifl: CalcFlags,
    geopos: [f64; 3],
) -> Result<LunarEclipseHow, Error> {
    if !(crate::constants::RISE_SET_GEOALT_MIN..=crate::constants::RISE_SET_GEOALT_MAX)
        .contains(&geopos[2])
    {
        return Err(Error::CError(format!(
            "location for eclipses must be between {:.0} and {:.0} m above sea",
            crate::constants::RISE_SET_GEOALT_MIN,
            crate::constants::RISE_SET_GEOALT_MAX
        )));
    }

    // Strip TOPOCTR/JPLHOR/JPLHOR_APPROX before the geocentric core runs (swecl.c:3219-3220) --
    // those refinements apply only to the separate topocentric az/alt call below.
    let ifl = ifl & !CalcFlags::TOPOCTR & !(CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX);

    let core = lun_eclipse_how(eph, tjd_ut, ifl)?;
    let mut flags = core.flags;

    let topo_config = {
        let mut c = eph.config().clone();
        c.topographic = Some(TopoPosition {
            longitude: geopos[0],
            latitude: geopos[1],
            altitude: geopos[2],
        });
        c
    };
    let lm = eph
        .calc_ut_with_config(
            tjd_ut,
            Body::Moon,
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
        [lm[0], lm[1]],
    );

    if xaz[2] <= 0.0 {
        flags = EclipseFlags::empty();
    }

    Ok(LunarEclipseHow {
        umbral_magnitude: core.umbral_magnitude,
        penumbral_magnitude: core.penumbral_magnitude,
        azimuth: xaz[0],
        true_altitude: xaz[1],
        apparent_altitude: xaz[2],
        distance_from_opposition: core.distance_from_opposition,
        saros_series: core.saros_series,
        saros_member: core.saros_member,
        flags,
    })
}

/// Global lunar-eclipse search result: `tret[0..8]` per `swe_lun_eclipse_when` (swecl.c:3389-3616,
/// ref doc §4). `tret[1]` is unused for lunar eclipses (index-parity padding with the solar
/// `tret[]` layout) and omitted here.
#[derive(Debug, Clone, Copy)]
pub struct LunarEclipseGlobal {
    /// Time (UT) of maximum eclipse: minimum selenocentric Sun/Earth-shadow angular separation.
    /// `tret[0]`.
    pub time_maximum: f64,
    /// Time (UT) of partial (umbra) phase begin, `0.0` if only penumbral. `tret[2]`.
    pub time_partial_begin: f64,
    /// Time (UT) of partial (umbra) phase end, `0.0` if only penumbral. `tret[3]`.
    pub time_partial_end: f64,
    /// Time (UT) of totality begin, `0.0` unless the eclipse is total. `tret[4]`.
    pub time_totality_begin: f64,
    /// Time (UT) of totality end, `0.0` unless the eclipse is total. `tret[5]`.
    pub time_totality_end: f64,
    /// Time (UT) of penumbral phase begin (always set for any eclipse type). `tret[6]`.
    pub time_penumbral_begin: f64,
    /// Time (UT) of penumbral phase end (always set for any eclipse type). `tret[7]`.
    pub time_penumbral_end: f64,
    /// Eclipse-type classification: exactly one of TOTAL/PARTIAL/PENUMBRAL. Never empty -- the
    /// search retries indefinitely (bounded only by ephemeris range) until a matching eclipse of
    /// a requested type is found.
    pub flags: EclipseFlags,
}

/// Global lunar-eclipse search: find the next (or, if `backward`, previous) lunar eclipse of a
/// type in `ifltype` (empty = any of TOTAL/PARTIAL/PENUMBRAL) starting from `tjd_start` (UT). No
/// geographic position -- purely geocentric. Port of `swe_lun_eclipse_when` (swecl.c:3389-3616,
/// ref doc §4).
pub(crate) fn lun_eclipse_when(
    eph: &Ephemeris,
    tjd_start: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    backward: bool,
) -> Result<LunarEclipseGlobal, Error> {
    let ifl = ifl & crate::calc::EPHMASK;
    let config = eph.config();

    // `ifltype` normalization (§4.1): solar-only CENTRAL/NONCENTRAL bits are meaningless here and
    // stripped unconditionally; ANNULAR/HYBRID (annular-total) don't exist for lunar eclipses --
    // stripped, and if nothing else survives, that's an unsatisfiable request (would infinite-loop
    // the search below).
    let mut ifltype = ifltype & !(EclipseFlags::CENTRAL | EclipseFlags::NONCENTRAL);
    if ifltype.intersects(EclipseFlags::ANNULAR | EclipseFlags::HYBRID) {
        ifltype &= !(EclipseFlags::ANNULAR | EclipseFlags::HYBRID);
        if ifltype.is_empty() {
            return Err(Error::CError(
                "annular lunar eclipses don't exist".to_string(),
            ));
        }
    }
    if ifltype.is_empty() {
        ifltype = EclipseFlags::ALLTYPES_LUNAR;
    }

    let direction = if backward { -1.0 } else { 1.0 };
    let iflag_cart = CalcFlags::EQUATORIAL | ifl | CalcFlags::XYZ;

    let mut k = ((tjd_start - J2000) / 365.2425 * 12.3685).trunc();
    k -= direction;

    'next_try: loop {
        let mut tret = [0.0f64; 8];

        // Full-moon (synodic-month) stepping via Meeus's lunation number K, plus the F-argument
        // node-proximity pre-filter (§4.2/§4.3).
        let kk = k + 0.5;
        let tt_ = kk / 1236.85;
        let t2 = tt_ * tt_;
        let t3 = t2 * tt_;
        let t4 = t3 * tt_;
        let f_deg = normalize_degrees(
            160.7108 + 390.67050274 * kk - 0.0016341 * t2 - 0.00000227 * t3 + 0.000000011 * t4,
        );
        let mut ff = f_deg;
        if ff > 180.0 {
            ff -= 180.0;
        }
        if ff > 21.0 && ff < 159.0 {
            k += direction;
            continue 'next_try;
        }

        // Approximate time of geocentric maximum eclipse (Meeus, German ed., p.381, §4.4).
        let mut tjd = 2451550.09765 + 29.530588853 * kk + 0.0001337 * t2 - 0.000000150 * t3
            + 0.00000000073 * t4;
        let m = normalize_degrees(2.5534 + 29.10535669 * kk - 0.0000218 * t2 - 0.00000011 * t3);
        let mm = normalize_degrees(
            201.5643 + 385.81693528 * kk + 0.1017438 * t2 + 0.00001239 * t3 + 0.000000058 * t4,
        );
        let om = normalize_degrees(124.7746 - 1.56375580 * kk + 0.0020691 * t2 + 0.00000215 * t3);
        let e = 1.0 - 0.002516 * tt_ - 0.0000074 * t2;
        let a1 = normalize_degrees(299.77 + 0.107408 * kk - 0.009173 * t2);
        let m_rad = m * DEGTORAD;
        let mm_rad = mm * DEGTORAD;
        let f_rad = f_deg * DEGTORAD;
        let om_rad = om * DEGTORAD;
        // Literal C quirk (swecl.c:3469): `Om` is already in radians here, so `sin(Om)` is
        // dimensionless -- multiplying by `DEGTORAD` again is not a unit conversion, it's part of
        // the tabulated Meeus coefficient. Preserve exactly, do not "fix".
        let f1_rad = f_rad - 0.02665 * om_rad.sin() * DEGTORAD;
        let a1_rad = a1 * DEGTORAD;
        tjd =
            tjd - 0.4075 * mm_rad.sin() + 0.1721 * e * m_rad.sin() + 0.0161 * (2.0 * mm_rad).sin()
                - 0.0097 * (2.0 * f1_rad).sin()
                + 0.0073 * e * (mm_rad - m_rad).sin()
                - 0.0050 * e * (mm_rad + m_rad).sin()
                - 0.0023 * (mm_rad - 2.0 * f1_rad).sin()
                + 0.0021 * e * (2.0 * m_rad).sin()
                + 0.0012 * (mm_rad + 2.0 * f1_rad).sin()
                + 0.0006 * e * (2.0 * mm_rad + m_rad).sin()
                - 0.0004 * (3.0 * mm_rad).sin()
                - 0.0003 * e * (m_rad + 2.0 * f1_rad).sin()
                + 0.0003 * a1_rad.sin()
                - 0.0002 * e * (m_rad - 2.0 * f1_rad).sin()
                - 0.0002 * e * (2.0 * mm_rad - m_rad).sin()
                - 0.0002 * om_rad.sin();

        // Precise refinement to the instant of minimum selenocentric Sun/Earth-shadow angular
        // separation (§4.5). `tjd` is ET/TT throughout; UT conversion happens once, after
        // convergence.
        let dtstart = if !(2_100_000.0..=2_500_000.0).contains(&tjd) {
            5.0
        } else {
            0.1
        };
        let mut dt = dtstart;
        while dt > 0.001 {
            let mut dc = [0.0f64; 3];
            let mut t = tjd - dt;
            for dc_i in dc.iter_mut() {
                let xs = eph.calc(t, Body::Sun, iflag_cart)?.data;
                let xm = eph.calc(t, Body::Moon, iflag_cart)?.data;
                let xs_sel = [xs[0] - xm[0], xs[1] - xm[1], xs[2] - xm[2]];
                let xm_sel = [-xm[0], -xm[1], -xm[2]];
                let ds =
                    (xs_sel[0] * xs_sel[0] + xs_sel[1] * xs_sel[1] + xs_sel[2] * xs_sel[2]).sqrt();
                let dm =
                    (xm_sel[0] * xm_sel[0] + xm_sel[1] * xm_sel[1] + xm_sel[2] * xm_sel[2]).sqrt();
                let xa = [xs_sel[0] / ds, xs_sel[1] / ds, xs_sel[2] / ds];
                let xb = [xm_sel[0] / dm, xm_sel[1] / dm, xm_sel[2] / dm];
                let rearth = (REARTH / dm).asin() * RADTODEG;
                let rsun = (RSUN / ds).asin() * RADTODEG;
                *dc_i = dot_prod_unit(xa, xb).acos() * RADTODEG;
                *dc_i -= rearth + rsun;
                t += dt;
            }
            let (dtint, _) = crate::math::find_maximum(dc[0], dc[1], dc[2], dt);
            tjd += dtint + dt;
            dt /= 4.0;
        }

        // 3-pass fixed-point ET->UT conversion (swecl.c:3525-3527).
        let tjds1 = tjd - crate::deltat::calc_deltat(tjd, config);
        let tjds2 = tjd - crate::deltat::calc_deltat(tjds1, config);
        let tjd = tjd - crate::deltat::calc_deltat(tjds2, config);

        // Confirm eclipse and reject wrong types (§4.6). Uses the geocentric core directly
        // (equivalent to calling the public wrapper with `geopos = NULL`, which skips the
        // topocentric az/alt gate entirely).
        let core = lun_eclipse_how(eph, tjd, ifl)?;
        if core.flags.is_empty() {
            k += direction;
            continue 'next_try;
        }
        tret[0] = tjd;
        if (backward && tret[0] >= tjd_start - 0.0001)
            || (!backward && tret[0] <= tjd_start + 0.0001)
        {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::PENUMBRAL)
            && core.flags.contains(EclipseFlags::PENUMBRAL)
        {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::PARTIAL) && core.flags.contains(EclipseFlags::PARTIAL) {
            k += direction;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::TOTAL) && core.flags.contains(EclipseFlags::TOTAL) {
            k += direction;
            continue 'next_try;
        }
        let retflag = core.flags;

        // Contact-time computation via `dcore`-based zero search (§4.7). `o` controls how many
        // contact-pairs to compute based on the eclipse's actual type -- a total eclipse also
        // gets partial and penumbral contacts (physical containment).
        let o = if retflag.contains(EclipseFlags::PENUMBRAL) {
            0
        } else if retflag.contains(EclipseFlags::PARTIAL) {
            1
        } else {
            2
        };
        let dta = 2.0 / 24.0;

        for n in 0..=o {
            let (i1, i2) = match n {
                0 => (6usize, 7usize),
                1 => (2usize, 3usize),
                _ => (4usize, 5usize),
            };

            // Stage A: coarse bracket, sampling `dcore` at `tjd - dta, tjd, tjd + dta`.
            let mut dc = [0.0f64; 3];
            let mut t = tjd - dta;
            for dc_i in dc.iter_mut() {
                let c = lun_eclipse_how(eph, t, ifl)?;
                *dc_i = lun_contact_dc(n, &c);
                t += dta;
            }
            // Divergence from C on `find_zero` failure: C ignores the failure return and
            // proceeds with a stale/zero-initialized `dt1`/`dt2` (see docs/c-ref-eclipse-
            // lunar.md §4.7 and the analogous documented divergence in sol_eclipse_when_glob).
            // We leave the slots at 0.0 and skip refinement instead -- unreachable for a
            // confirmed eclipse (the type check above guarantees a sign change in `dc`).
            if let Some((dt1, dt2)) = crate::math::find_zero(dc[0], dc[1], dc[2], dta) {
                let dtb = (dt1 + dta) / 2.0;
                tret[i1] = tjd + dt1 + dta;
                tret[i2] = tjd + dt2 + dta;

                // Stage B: 3 rounds of 2-point secant/Newton refinement, halving `dt` each round.
                let mut dt = dtb / 2.0;
                for _ in 0..3 {
                    for &j in &[i1, i2] {
                        let mut dc2 = [0.0f64; 2];
                        let mut t = tret[j] - dt;
                        for dc2_i in dc2.iter_mut() {
                            let c = lun_eclipse_how(eph, t, ifl)?;
                            *dc2_i = lun_contact_dc(n, &c);
                            t += dt;
                        }
                        let dt1 = dc2[1] / ((dc2[1] - dc2[0]) / dt);
                        tret[j] -= dt1;
                    }
                    dt /= 2.0;
                }
            }
        }

        return Ok(LunarEclipseGlobal {
            time_maximum: tret[0],
            time_partial_begin: tret[2],
            time_partial_end: tret[3],
            time_totality_begin: tret[4],
            time_totality_end: tret[5],
            time_penumbral_begin: tret[6],
            time_penumbral_end: tret[7],
            flags: retflag,
        });
    }
}

/// Local lunar-eclipse search result: `tret[0..10]` + `attr[0..11]` per `swe_lun_eclipse_when_loc`
/// (swecl.c:3644-3739, ref doc §5). `tret[]` index semantics match [`LunarEclipseGlobal`]'s (not
/// [`SolarEclipseLocal`]'s different layout) plus two new slots for moonrise/moonset.
#[derive(Debug, Clone, Copy)]
pub struct LunarEclipseLocal {
    /// Time (UT) of maximum eclipse as visible from this location -- re-anchored to
    /// moonrise/moonset if the true geocentric maximum wasn't visible here. `tret[0]`.
    pub time_maximum: f64,
    /// Time (UT) of partial (umbra) phase begin, `0.0` if not applicable or clipped away by
    /// moonrise/moonset. `tret[2]`.
    pub time_partial_begin: f64,
    /// Time (UT) of partial (umbra) phase end, `0.0` if not applicable or clipped away. `tret[3]`.
    pub time_partial_end: f64,
    /// Time (UT) of totality begin, `0.0` if not applicable or clipped away. `tret[4]`.
    pub time_totality_begin: f64,
    /// Time (UT) of totality end, `0.0` if not applicable or clipped away. `tret[5]`.
    pub time_totality_end: f64,
    /// Time (UT) of penumbral phase begin, `0.0` if clipped away by moonrise. `tret[6]`.
    pub time_penumbral_begin: f64,
    /// Time (UT) of penumbral phase end, `0.0` if clipped away by moonset. `tret[7]`.
    pub time_penumbral_end: f64,
    /// Time (UT) of moonrise, if it occurs during the eclipse; `0.0` otherwise. `tret[8]`.
    pub time_moonrise: f64,
    /// Time (UT) of moonset, if it occurs during the eclipse; `0.0` otherwise. `tret[9]`.
    pub time_moonset: f64,
    /// Local circumstances (`attr[]`) at whichever instant was written last: the moment of
    /// maximum eclipse, unless a moonrise/moonset re-anchor overwrote it.
    pub attr: LunarEclipseHow,
    /// Eclipse-type classification (TOTAL/PARTIAL/PENUMBRAL) OR'd with VISIBLE and whichever of
    /// MAX/PARTBEG/PARTEND/TOTBEG/TOTEND/PENUMBBEG/PENUMBEND_VISIBLE applied at some contact. The
    /// search loop retries internally until a visible occurrence is found -- never empty.
    pub flags: EclipseFlags,
}

/// Local lunar-eclipse search: find the next (or, if `backward`, previous) lunar eclipse that is
/// at least partly visible (Moon above the horizon during some phase) from `geopos`, clipping
/// contact times to moonrise/moonset as needed. Port of `swe_lun_eclipse_when_loc`
/// (swecl.c:3644-3739, ref doc §5). `geopos` = [longitude, latitude, height above sea (m)],
/// degrees/degrees/meters.
pub(crate) fn lun_eclipse_when_loc(
    eph: &Ephemeris,
    tjd_start: f64,
    ifl: CalcFlags,
    geopos: [f64; 3],
    backward: bool,
) -> Result<LunarEclipseLocal, Error> {
    if !(crate::constants::RISE_SET_GEOALT_MIN..=crate::constants::RISE_SET_GEOALT_MAX)
        .contains(&geopos[2])
    {
        return Err(Error::CError(format!(
            "location for eclipses must be between {:.0} and {:.0} m above sea",
            crate::constants::RISE_SET_GEOALT_MIN,
            crate::constants::RISE_SET_GEOALT_MAX
        )));
    }
    let ifl = ifl & !(CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX);

    let mut tjd_start = tjd_start;

    'next_lun_ecl: loop {
        let glob = lun_eclipse_when(eph, tjd_start, ifl, EclipseFlags::empty(), backward)?;
        let mut tret = [
            glob.time_maximum,
            0.0,
            glob.time_partial_begin,
            glob.time_partial_end,
            glob.time_totality_begin,
            glob.time_totality_end,
            glob.time_penumbral_begin,
            glob.time_penumbral_end,
            0.0,
            0.0,
        ];

        // Visibility scan (§5 step 4): descending i=7..=0 (skip i==1, unused slot; skip any
        // tret[i]==0, not-applicable contact) -- order doesn't affect the OR'd result, only
        // evaluation order.
        let mut retflag = EclipseFlags::empty();
        for i in (0..=7).rev() {
            if i == 1 || tret[i] == 0.0 {
                continue;
            }
            let h = swe_lun_eclipse_how(eph, tret[i], ifl, geopos)?;
            if h.apparent_altitude > 0.0 {
                retflag |= EclipseFlags::VISIBLE;
                retflag |= match i {
                    0 => EclipseFlags::MAX_VISIBLE,
                    2 => EclipseFlags::PARTBEG_VISIBLE,
                    3 => EclipseFlags::PARTEND_VISIBLE,
                    4 => EclipseFlags::TOTBEG_VISIBLE,
                    5 => EclipseFlags::TOTEND_VISIBLE,
                    6 => EclipseFlags::PENUMBBEG_VISIBLE,
                    7 => EclipseFlags::PENUMBEND_VISIBLE,
                    _ => unreachable!(),
                };
            }
        }
        if !retflag.contains(EclipseFlags::VISIBLE) {
            tjd_start = if backward {
                tret[0] - 25.0
            } else {
                tret[0] + 25.0
            };
            continue 'next_lun_ecl;
        }

        // Moonrise/moonset clipping (§5 step 6). Both searches start just before penumbral
        // begin (`tret[6] - 0.001`), matching C exactly. `Error::CircumpolarBody` from either
        // call means "no rise/set found in the window" -- skip clipping entirely, matching C's
        // `retc >= 0` guard (a genuine `ERR` still propagates).
        let mut tjd_max = tret[0];
        let rise = eph.rise_trans(
            tret[6] - 0.001,
            Body::Moon,
            None,
            ifl,
            RiseSetFlags::RISE | RiseSetFlags::DISC_BOTTOM,
            geopos,
            0.0,
            0.0,
        );
        let clip = match rise {
            Ok(r) => {
                let tjdr = r.time;
                match eph.rise_trans(
                    tret[6] - 0.001,
                    Body::Moon,
                    None,
                    ifl,
                    RiseSetFlags::SET | RiseSetFlags::DISC_BOTTOM,
                    geopos,
                    0.0,
                    0.0,
                ) {
                    Ok(s) => Some((tjdr, s.time)),
                    Err(Error::CircumpolarBody) => None,
                    Err(e) => return Err(e),
                }
            }
            Err(Error::CircumpolarBody) => None,
            Err(e) => return Err(e),
        };

        if let Some((tjdr, tjds)) = clip {
            if tjds < tret[6] || (tjds > tjdr && tjdr > tret[7]) {
                tjd_start = if backward {
                    tret[0] - 25.0
                } else {
                    tret[0] + 25.0
                };
                continue 'next_lun_ecl;
            }
            // HAZARD (ref doc §5, FP/logic hazard note): the second block below reads
            // `tret[6]`/`tret[7]` which may have just been mutated by the first block -- this is
            // the C source's own behavior (sequential, not independent), preserved exactly.
            if tjdr > tret[6] && tjdr < tret[7] {
                tret[6] = 0.0;
                for t in &mut tret[2..=5] {
                    if tjdr > *t {
                        *t = 0.0;
                    }
                }
                tret[8] = tjdr;
                if tjdr > tret[0] {
                    tjd_max = tjdr;
                }
            }
            if tjds > tret[6] && tjds < tret[7] {
                tret[7] = 0.0;
                for t in &mut tret[2..=5] {
                    if tjds < *t {
                        *t = 0.0;
                    }
                }
                tret[9] = tjds;
                if tjds < tret[0] {
                    tjd_max = tjds;
                }
            }
        }

        tret[0] = tjd_max;
        let how = swe_lun_eclipse_how(eph, tjd_max, ifl, geopos)?;
        if how.flags.is_empty() {
            tjd_start = if backward {
                tret[0] - 25.0
            } else {
                tret[0] + 25.0
            };
            continue 'next_lun_ecl;
        }
        retflag |= how.flags & EclipseFlags::ALLTYPES_LUNAR;

        return Ok(LunarEclipseLocal {
            time_maximum: tret[0],
            time_partial_begin: tret[2],
            time_partial_end: tret[3],
            time_totality_begin: tret[4],
            time_totality_end: tret[5],
            time_penumbral_begin: tret[6],
            time_penumbral_end: tret[7],
            time_moonrise: tret[8],
            time_moonset: tret[9],
            attr: how,
            flags: retflag,
        });
    }
}

// === Occultations (Moon occults a planet/asteroid/fixed star) ===
//
// `swe_lun_occult_where` and `swe_lun_occult_when_glob` (swecl.c:606-630, 1572-1984,
// docs/c-ref-occultation.md). Both reuse `eclipse_where`/`calc_planet_star`/`body_radius_au`
// verbatim -- the Sun is just the `ipl=Body::Sun, starname=None` special case those already
// handle generically. The delta from the solar port is confined to: asteroid-134340->Pluto
// aliasing, and (for `_when_glob`) a generic Moon-body elongation bracketing search in place of
// solar's Meeus lunation-number estimate (occultation search must work for any sidereal period,
// including a fixed star's zero proper motion).

/// Asteroid-134340 (numbered-asteroid Pluto) aliasing to `Body::Pluto`, applied by all three
/// `swe_lun_occult_*` entry points (swecl.c:620-623, 1599-1600, 2084-2085).
fn normalize_occulted_body(ipl: Body) -> Body {
    match ipl {
        Body::Asteroid(id) if id.mpc_number() == 134340 => Body::Pluto,
        other => other,
    }
}

/// Geographic position of maximal occultation of `ipl`/`starname` by the Moon at `tjd_ut` (UT).
/// Port of `swe_lun_occult_where` (swecl.c:606-630, §1) -- a thin wrapper over [`eclipse_where`]
/// threading the occulted body through in place of the Sun; same shape/masking as
/// [`sol_eclipse_where`]. C additionally calls `eclipse_how` here purely to catch an error from
/// it and to fill `attr[3]`, both of which the Rust port omits for the same reason
/// `sol_eclipse_where` does: local-circumstance attributes live in a separate function
/// ([`eclipse_how`], exposed for occultations by a later task).
pub(crate) fn lun_occult_where(
    eph: &Ephemeris,
    tjd_ut: f64,
    ipl: Body,
    starname: Option<&str>,
    ifl: CalcFlags,
) -> Result<EclipseWhere, Error> {
    let ipl = normalize_occulted_body(ipl);
    eclipse_where(eph, tjd_ut, ipl, starname, ifl & crate::calc::EPHMASK)
}

/// Global occultation search result: `tret[0..10]` per `swe_lun_occult_when_glob`
/// (swecl.c:1572-1984, §2). Same slot layout as [`SolarEclipseGlobal`], but `time_ra_conjunction`
/// (`tret[1]`) is the transit instant of the *occulted body* (not necessarily the Sun), and the
/// search never produces `ANNULAR`/`HYBRID` for `ipl != Body::Sun` (rejected/stripped from
/// `ifltype` up front -- see [`lun_occult_when_glob`]).
#[derive(Debug, Clone, Copy)]
pub struct OccultGlobal {
    /// Time (UT) of maximum occultation: geocentric minimum Moon-body angular separation.
    /// `tret[0]`.
    pub time_maximum: f64,
    /// Time (UT) when the occulted body transits the meridian relative to the Moon (geocentric
    /// ecliptic-longitude conjunction), or `0.0` if no such instant falls within the occultation
    /// window. `tret[1]`.
    pub time_ra_conjunction: f64,
    /// Time (UT) of occultation begin, first contact anywhere on Earth. `tret[2]`.
    pub time_begin: f64,
    /// Time (UT) of occultation end, last contact anywhere on Earth. `tret[3]`.
    pub time_end: f64,
    /// Time (UT) of totality/annularity begin, `0.0` if partial. `tret[4]`.
    pub time_totality_begin: f64,
    /// Time (UT) of totality/annularity end, `0.0` if partial. `tret[5]`.
    pub time_totality_end: f64,
    /// Time (UT) of center-line begin, `0.0` if noncentral. `tret[6]`.
    pub time_centerline_begin: f64,
    /// Time (UT) of center-line end, `0.0` if noncentral. `tret[7]`.
    pub time_centerline_end: f64,
    /// Occultation-type classification (CENTRAL/NONCENTRAL combined with TOTAL/PARTIAL, plus
    /// ANNULAR/HYBRID when `ipl == Body::Sun`). Never empty -- the search retries indefinitely
    /// (bounded only by ephemeris range) until a matching occultation is found.
    pub flags: EclipseFlags,
}

/// Global occultation search: find the next (or, if `backward`, previous) occultation of
/// `ipl`/`starname` by the Moon anywhere on Earth after/before `tjd_start` (UT), restricted to
/// types in `ifltype`. Port of `swe_lun_occult_when_glob` (swecl.c:1572-1984, §2).
/// `ifltype = EclipseFlags::empty()` means all types.
///
/// Structurally a near-duplicate of [`sol_eclipse_when_glob`] -- same contact-time refinement
/// (`contact_dc`/`find_zero`), annular-total detection, and transit computation -- but the rough
/// initial estimate uses a generic Newton-style Moon-body elongation bracket (`dl/13` per
/// swecl.c:1640-1666) instead of solar's Meeus lunation-number formula, since occultation search
/// must work for any occulted body's sidereal period (including a fixed star's zero proper
/// motion). C's `SE_ECL_ONE_TRY` early-return optimization (a single-conjunction-check mode for
/// callers willing to resume the search themselves) has no equivalent here -- this always
/// searches until a matching occultation is found, matching the exposed `backward: bool`-only
/// signature (same choice already made for [`sol_eclipse_when_glob`], which has no one-try mode
/// at all).
pub(crate) fn lun_occult_when_glob(
    eph: &Ephemeris,
    tjd_start: f64,
    ipl: Body,
    starname: Option<&str>,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    backward: bool,
) -> Result<OccultGlobal, Error> {
    let ipl = normalize_occulted_body(ipl);
    let ifl = ifl & crate::calc::EPHMASK;
    let config = eph.config();

    if ifltype == (EclipseFlags::PARTIAL | EclipseFlags::CENTRAL) {
        return Err(Error::CError(
            "central partial eclipses do not exist".to_string(),
        ));
    }

    // `ipl == SE_SUN` (as opposed to a real occultation of a planet/asteroid/star) is the only
    // case where annular/annular-total ("hybrid") occultations are geometrically meaningful --
    // C's `ipl != SE_SUN` predicate, ported directly here rather than replicating C's separate
    // `if (ipl < 0) ipl = 0` clamp (a raw-int sentinel for "star lookup, ignore ipl" with no
    // equivalent needed once `starname: Option<&str>` already carries that distinction).
    let is_sun = ipl == Body::Sun && starname.unwrap_or("").is_empty();
    let mut ifltype = ifltype;
    if !is_sun {
        // C tests this with `&`/`==` differently for the hard-error vs silent-strip cases
        // (swecl.c:1626-1634) -- ported literally, not unified.
        let stripped = ifltype & !(EclipseFlags::NONCENTRAL | EclipseFlags::CENTRAL);
        if stripped == EclipseFlags::ANNULAR || ifltype == EclipseFlags::HYBRID {
            return Err(Error::CError(format!(
                "annular occulation do not exist for object {} {}",
                ipl.to_raw_id(),
                starname.unwrap_or("")
            )));
        }
        if ifltype.intersects(EclipseFlags::ANNULAR | EclipseFlags::HYBRID) {
            ifltype &= !(EclipseFlags::ANNULAR | EclipseFlags::HYBRID);
        }
    }
    if ifltype.is_empty() {
        ifltype = EclipseFlags::TOTAL
            | EclipseFlags::PARTIAL
            | EclipseFlags::NONCENTRAL
            | EclipseFlags::CENTRAL;
        if is_sun {
            ifltype |= EclipseFlags::ANNULAR | EclipseFlags::HYBRID;
        }
    }
    // C tests these two with a bitwise `&` (any-bit-set), unlike solar's `==` (bare-type-only)
    // expansion (swecl.c:1640-1642 vs swecl.c:1234-1236) -- literal divergence, not a typo.
    if ifltype.intersects(EclipseFlags::TOTAL | EclipseFlags::ANNULAR | EclipseFlags::HYBRID) {
        ifltype |= EclipseFlags::NONCENTRAL | EclipseFlags::CENTRAL;
    }
    if ifltype.contains(EclipseFlags::PARTIAL) {
        ifltype |= EclipseFlags::NONCENTRAL;
    }

    let direction = if backward { -1.0 } else { 1.0 };
    let iflag = CalcFlags::EQUATORIAL | ifl;
    let iflag_cart = iflag | CalcFlags::XYZ;
    let de_km = 6378.140;

    let mut t = tjd_start;

    'next_try: loop {
        // §2 step 1: rough conjunction in ecliptic longitude (swecl.c:1640-1666). Plain `ifl`
        // (no EQUATORIAL/SPEED) -- polar ecliptic lon/lat/dist in degrees.
        let ls0 = calc_planet_star(eph, t, ipl, starname, ifl)?;
        if let Some(name) = starname
            && !name.is_empty()
            && ls0[1].abs() > 7.0
        {
            return Err(Error::CError(format!(
                "occultation never occurs: star {name} has ecl. lat. {:.1}",
                ls0[1]
            )));
        }
        let mut ls = ls0;
        let mut lm = eph.calc(t, Body::Moon, ifl)?.data;
        let mut dl = normalize_degrees(ls[0] - lm[0]);
        if backward {
            dl -= 360.0;
        }
        while dl.abs() > 0.1 {
            t += dl / 13.0;
            ls = calc_planet_star(eph, t, ipl, starname, ifl)?;
            lm = eph.calc(t, Body::Moon, ifl)?.data;
            dl = normalize_degrees(ls[0] - lm[0]);
            if dl > 180.0 {
                dl -= 360.0;
            }
        }
        let mut tjd = t;

        // §2 step 2: latitude-difference gate.
        if (ls[1] - lm[1]).abs() > 2.0 {
            t += direction * 20.0;
            continue 'next_try;
        }

        // §2 step 3: occulted-body angular radius (reuses the same helper as eclipse_where).
        let body_radius = body_radius_au(ipl, starname);

        // §2 step 4: refine time of maximum occultation (parabola-vertex bracketing,
        // dtstart=1, dtdiv=3 -- fixed, unlike solar's conditional dtstart).
        let mut dt = 1.0;
        while dt > 0.0001 {
            let mut dc = [0.0f64; 3];
            let mut tt = tjd - dt;
            for dc_i in dc.iter_mut() {
                let ls2 = calc_planet_star(eph, tt, ipl, starname, iflag)?;
                let lm2 = eph.calc(tt, Body::Moon, iflag)?.data;
                let xs = calc_planet_star(eph, tt, ipl, starname, iflag_cart)?;
                let xm = eph.calc(tt, Body::Moon, iflag_cart)?.data;
                let rmoon = (RMOON / lm2[2]).asin() * RADTODEG;
                let rsun = (body_radius / ls2[2]).asin() * RADTODEG;
                *dc_i = dot_prod_unit([xs[0], xs[1], xs[2]], [xm[0], xm[1], xm[2]]).acos()
                    * RADTODEG
                    - (rmoon + rsun);
                tt += dt;
            }
            let (dtint, _) = crate::math::find_maximum(dc[0], dc[1], dc[2], dt);
            tjd += dtint + dt;
            dt /= 3.0;
        }

        // §2 step 5: single ET->UT deltaT subtraction (not solar's 3-pass fixed point).
        let tjd = tjd - crate::deltat::calc_deltat(tjd, config);

        // §2 step 6: C calls `eclipse_where` twice here with identical arguments (once to
        // confirm an occultation exists anywhere, once more to fetch the classification) --
        // since `tjd` is unchanged between the two calls and `eclipse_where` is a pure function
        // of its arguments, the second call is provably redundant, and reusing one result here
        // makes C's dead "retflag == 0, extremely small percentage" fallback (swecl.c:1762-1766)
        // genuinely unreachable, since that branch only fires when the *second* (identical) call
        // returns empty despite the first not being empty.
        let where_result = eclipse_where(eph, tjd, ipl, starname, ifl)?;
        if where_result.flags.is_empty() {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }

        let mut tret = [0.0f64; 8];
        tret[0] = tjd;
        if (backward && tret[0] >= tjd_start - 0.0001)
            || (!backward && tret[0] <= tjd_start + 0.0001)
        {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }

        let mut retflag = where_result.flags;

        if !ifltype.contains(EclipseFlags::NONCENTRAL) && retflag.contains(EclipseFlags::NONCENTRAL)
        {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::CENTRAL) && retflag.contains(EclipseFlags::CENTRAL) {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::ANNULAR) && retflag.contains(EclipseFlags::ANNULAR) {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::PARTIAL) && retflag.contains(EclipseFlags::PARTIAL) {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }
        if !ifltype.intersects(EclipseFlags::TOTAL | EclipseFlags::HYBRID)
            && retflag.contains(EclipseFlags::TOTAL)
        {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }

        // Contact-time refinement (§2 step 6 cont'd, reusing solar's shared `contact_dc`): n=0
        // occultation begin/end (always), n=1 totality/annularity begin/end (skip if PARTIAL),
        // n=2 center-line begin/end (skip if NONCENTRAL). `dtb` is NOT divided by 3 here, unlike
        // solar's `dtb` -- literal C divergence (swecl.c:1842-1843 vs swecl.c:1385-1386).
        let o = if retflag.contains(EclipseFlags::PARTIAL) {
            0
        } else if retflag.contains(EclipseFlags::NONCENTRAL) {
            1
        } else {
            2
        };
        let dta = 2.0 / 24.0;
        let dtb = 10.0 / 24.0 / 60.0;

        for n in 0..=o {
            let (i1, i2) = match n {
                0 => (2usize, 3usize),
                1 => (4usize, 5usize),
                _ => (6usize, 7usize),
            };

            let mut dc = [0.0f64; 3];
            let mut t2 = tjd - dta;
            for dc_i in dc.iter_mut() {
                let w = eclipse_where(eph, t2, ipl, starname, ifl)?;
                *dc_i = contact_dc(n, &w, de_km);
                t2 += dta;
            }
            // Same intentional divergence from C on `find_zero` failure as
            // `sol_eclipse_when_glob`: leave the slots at 0.0 and skip refinement rather than
            // refining around a stale/zero value -- unreachable for a confirmed occultation.
            if let Some((dt1, dt2)) = crate::math::find_zero(dc[0], dc[1], dc[2], dta) {
                tret[i1] = tjd + dt1 + dta;
                tret[i2] = tjd + dt2 + dta;

                let mut dt = dtb;
                for _ in 0..3 {
                    for &j in &[i1, i2] {
                        let mut dc2 = [0.0f64; 2];
                        let mut t3 = tret[j] - dt;
                        for dc_i in dc2.iter_mut() {
                            let w = eclipse_where(eph, t3, ipl, starname, ifl)?;
                            *dc_i = contact_dc(n, &w, de_km);
                            t3 += dt;
                        }
                        let dt1 = dc2[1] / ((dc2[1] - dc2[0]) / dt);
                        tret[j] -= dt1;
                    }
                    dt /= 3.0;
                }
            }
        }

        // Annular-total (hybrid) detection -- unreachable in practice for `ipl != Body::Sun`
        // since `body_radius_au` returns a small (or zero, for a star) disc, but ported
        // faithfully rather than short-circuited (§ "Radius handling", c-ref-occultation.md).
        if retflag.contains(EclipseFlags::TOTAL) {
            let dc0 = eclipse_where(eph, tret[0], ipl, starname, ifl)?.core_diameter_km;
            let dc1 = eclipse_where(eph, tret[4], ipl, starname, ifl)?.core_diameter_km;
            let dc2 = eclipse_where(eph, tret[5], ipl, starname, ifl)?.core_diameter_km;
            if dc0 * dc1 < 0.0 || dc0 * dc2 < 0.0 {
                retflag |= EclipseFlags::HYBRID;
                retflag.remove(EclipseFlags::TOTAL);
            }
        }
        if !ifltype.contains(EclipseFlags::TOTAL) && retflag.contains(EclipseFlags::TOTAL) {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }
        if !ifltype.contains(EclipseFlags::HYBRID) && retflag.contains(EclipseFlags::HYBRID) {
            t = tjd + direction * 20.0;
            continue 'next_try;
        }

        // Time of maximum occultation at local apparent noon (transit of the occulted body, not
        // necessarily the Sun): check for a sign change between occultation begin/end, then
        // secant-iterate to the exact geocentric ecliptic-longitude conjunction instant.
        let mut dc_transit = [0.0f64; 2];
        for (i, dc_i) in dc_transit.iter_mut().enumerate() {
            let tt = tret[2 + i] + crate::deltat::calc_deltat(tret[2 + i], config);
            let ls = calc_planet_star(eph, tt, ipl, starname, iflag)?;
            let lm = eph.calc(tt, Body::Moon, iflag)?.data;
            let mut d = normalize_degrees(ls[0] - lm[0]);
            if d > 180.0 {
                d -= 360.0;
            }
            *dc_i = d;
        }
        if dc_transit[0] * dc_transit[1] >= 0.0 {
            tret[1] = 0.0;
        } else {
            let mut tjd_ra = tjd;
            let mut dt = 0.1;
            let dt1_init = (tret[3] - tret[2]) / 2.0;
            if dt1_init < dt {
                dt = dt1_init / 2.0;
            }
            while dt > 0.01 {
                let mut dc2 = [0.0f64; 2];
                let mut t2 = tjd_ra;
                for dc_i in dc2.iter_mut() {
                    let tt = t2 + crate::deltat::calc_deltat(t2, config);
                    let ls = calc_planet_star(eph, tt, ipl, starname, iflag)?;
                    let lm = eph.calc(tt, Body::Moon, iflag)?.data;
                    let mut d = normalize_degrees(ls[0] - lm[0]);
                    if d > 180.0 {
                        d -= 360.0;
                    }
                    if d > 180.0 {
                        d -= 360.0;
                    }
                    *dc_i = d;
                    t2 -= dt;
                }
                let a = (dc2[1] - dc2[0]) / dt;
                if a < 1e-10 {
                    break;
                }
                let dt1 = dc2[0] / a;
                tjd_ra += dt1;
                dt /= 3.0;
            }
            tret[1] = tjd_ra;
        }

        return Ok(OccultGlobal {
            time_maximum: tret[0],
            time_ra_conjunction: tret[1],
            time_begin: tret[2],
            time_end: tret[3],
            time_totality_begin: tret[4],
            time_totality_end: tret[5],
            time_centerline_begin: tret[6],
            time_centerline_end: tret[7],
            flags: retflag,
        });
    }
}
