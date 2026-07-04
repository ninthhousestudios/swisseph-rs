//! Physical, astronomical, and conversion constants used throughout the crate.
//!
//! Mirrors the constants defined in the C Swiss Ephemeris headers (`sweph.h`,
//! `swephexp.h`, `swephlib.h`) and their associated source files, grouped by
//! purpose (physical constants, unit conversions, reference epochs, ephemeris
//! validity ranges, etc.).

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Astronomical unit (meters).
pub const AUNIT: f64 = 1.49597870700e11;
/// Speed of light (m/s).
pub const CLIGHT: f64 = 2.99792458e8;
/// Heliocentric gravitational constant (m^3/s^2).
pub const HELGRAVCONST: f64 = 1.32712440017987e20;
/// Geocentric gravitational constant (m^3/s^2).
pub const GEOGCONST: f64 = 3.98600448e14;
/// Gaussian gravitational constant (radians/day).
pub const KGAUSS: f64 = 0.01720209895;
/// Earth's equatorial radius (meters).
pub const EARTH_RADIUS: f64 = 6378136.6;
/// Earth's flattening (oblateness), dimensionless.
pub const EARTH_OBLATENESS: f64 = 1.0 / 298.25642;
/// Earth's rotation speed (radians/day).
pub const EARTH_ROT_SPEED: f64 = 7.2921151467e-5 * 86400.0;
/// Sun-to-Earth mass ratio.
pub const SUN_EARTH_MRAT: f64 = 332946.050895;
/// Earth-to-Moon mass ratio.
pub const EARTH_MOON_MRAT: f64 = 1.0 / 0.0123000383;
/// Moon's mean distance from Earth (meters).
pub const MOON_MEAN_DIST: f64 = 384400000.0;
/// Moon's mean orbital inclination (degrees).
pub const MOON_MEAN_INCL: f64 = 5.1453964;
/// Moon's mean orbital eccentricity.
pub const MOON_MEAN_ECC: f64 = 0.054900489;
/// Light travel time across one astronomical unit (days).
pub const LIGHTTIME_AUNIT: f64 = 499.0047838362 / 3600.0 / 24.0;
/// One parsec expressed in astronomical units.
pub const PARSEC_TO_AUNIT: f64 = 206264.8062471;
/// Standard atmospheric temperature lapse rate (K/meter), used for refraction.
pub const LAPSE_RATE: f64 = 0.0065;
/// Conversion factor from km/s to AU/century.
pub const KM_S_TO_AU_CTY: f64 = 21.095;

// Valid observer-altitude range for rise/set (sweph.h:198-199, SEI_ECL_GEOALT_MIN/_MAX), meters.
// Eclipse local-circumstance functions (`eclipse_how`/`eclipse_when_loc`) reuse the same range.
/// Minimum valid observer altitude for rise/set and eclipse local-circumstance calculations (meters).
pub const RISE_SET_GEOALT_MIN: f64 = -500.0;
/// Maximum valid observer altitude for rise/set and eclipse local-circumstance calculations (meters).
pub const RISE_SET_GEOALT_MAX: f64 = 25000.0;

// Eclipse/occultation shadow-geometry body diameters, AU (swecl.c:80-84). `DSUN`'s numerator
// matches `PLANETARY_DIAMETERS[0]` (Sun) exactly, so `RSUN` and the general `drad` lookup agree
// for solar eclipses; `DMOON`'s 3476300.0 is the mean lunar radius used for shadow-cone geometry,
// distinct from `PLANETARY_DIAMETERS[1]`'s 3475000.0 (a different, general-purpose figure).
/// Sun's diameter used for eclipse/occultation shadow geometry (AU).
pub const DSUN: f64 = 1392000000.0 / AUNIT;
/// Moon's mean diameter used for eclipse shadow-cone geometry (AU).
pub const DMOON: f64 = 3476300.0 / AUNIT;
/// Earth's diameter used for eclipse shadow geometry (AU).
pub const DEARTH: f64 = 6378140.0 * 2.0 / AUNIT;
/// Sun's radius used for eclipse/occultation shadow geometry (AU).
pub const RSUN: f64 = DSUN / 2.0;
/// Moon's radius used for eclipse shadow-cone geometry (AU).
pub const RMOON: f64 = DMOON / 2.0;
/// Earth's radius used for eclipse shadow geometry (AU).
pub const REARTH: f64 = DEARTH / 2.0;

// ---------------------------------------------------------------------------
// Unit conversions
// ---------------------------------------------------------------------------

/// Multiplier to convert degrees to radians.
pub const DEGTORAD: f64 = PI / 180.0;
/// Multiplier to convert radians to degrees.
pub const RADTODEG: f64 = 180.0 / PI;
/// Multiplier to convert centiseconds of arc to radians.
pub const CSTORAD: f64 = DEGTORAD / 360000.0;
/// Multiplier to convert radians to centiseconds of arc.
pub const RADTOCS: f64 = RADTODEG * 360000.0;
/// Multiplier to convert centiseconds of arc to degrees.
pub const CS2DEG: f64 = 1.0 / 360000.0;
/// Full circle in radians (2 * pi).
pub const TWOPI: f64 = 2.0 * PI;
/// One astronomical unit expressed in kilometers.
pub const AUNIT_TO_KM: f64 = 149597870.700;
/// Multiplier to convert astronomical units to light-years.
pub const AUNIT_TO_LIGHTYEAR: f64 = 1.0 / 63241.07708427;
/// Multiplier to convert astronomical units to parsecs.
pub const AUNIT_TO_PARSEC: f64 = 1.0 / 206264.8062471;
// Verbatim arcsec-to-radian constant from the C source; full digits preserved.
/// One arcsecond expressed in radians.
#[allow(clippy::excessive_precision)]
pub const STR: f64 = 4.8481368110953599359e-6;

// ---------------------------------------------------------------------------
// Centisecond angle constants
// ---------------------------------------------------------------------------

/// One degree expressed in centiseconds of arc.
pub const DEG: i32 = 360000;
/// 15 degrees expressed in centiseconds of arc.
pub const DEG15: i32 = 15 * DEG;
/// 30 degrees expressed in centiseconds of arc.
pub const DEG30: i32 = 30 * DEG;
/// 60 degrees expressed in centiseconds of arc.
pub const DEG60: i32 = 60 * DEG;
/// 90 degrees expressed in centiseconds of arc.
pub const DEG90: i32 = 90 * DEG;
/// 120 degrees expressed in centiseconds of arc.
pub const DEG120: i32 = 120 * DEG;
/// 150 degrees expressed in centiseconds of arc.
pub const DEG150: i32 = 150 * DEG;
/// 180 degrees expressed in centiseconds of arc.
pub const DEG180: i32 = 180 * DEG;
/// 270 degrees expressed in centiseconds of arc.
pub const DEG270: i32 = 270 * DEG;
/// 360 degrees expressed in centiseconds of arc.
pub const DEG360: i32 = 360 * DEG;

// ---------------------------------------------------------------------------
// Reference epochs (Julian Day numbers)
// ---------------------------------------------------------------------------

/// Julian day number for epoch J2000.0 (2000 Jan 1, 12:00 TT).
pub const J2000: f64 = 2451545.0;
/// Julian day number for epoch B1950.0.
pub const B1950: f64 = 2433282.42345905;
/// Julian day number for epoch J1900.0.
pub const J1900: f64 = 2415020.0;
/// Julian day number for epoch B1850.0.
pub const B1850: f64 = 2396758.2035810;

// ---------------------------------------------------------------------------
// Precession model validity ranges (Julian centuries from J2000)
// ---------------------------------------------------------------------------

/// Validity range of the IAU 1976 precession model (Julian centuries from J2000).
pub const PREC_IAU_1976_CTIES: f64 = 2.0;
/// Validity range of the IAU 2000 precession model (Julian centuries from J2000).
pub const PREC_IAU_2000_CTIES: f64 = 2.0;
/// Validity range of the IAU 2006 precession model (Julian centuries from J2000).
pub const PREC_IAU_2006_CTIES: f64 = 75.0;

// ---------------------------------------------------------------------------
// JPL Horizons epoch constants
// ---------------------------------------------------------------------------

/// Reference epoch for the IAU 1980 nutation offsets used by JPL Horizons (Julian day number).
pub const DPSI_DEPS_IAU1980_TJD0_HORIZONS: f64 = 2437684.5;
/// Nutation-in-longitude offset (dpsi) at the IAU 1980 reference epoch (arcsec).
pub const DPSI_IAU1980_TJD0: f64 = 0.064284;
/// Nutation-in-obliquity offset (deps) at the IAU 1980 reference epoch (arcsec).
pub const DEPS_IAU1980_TJD0: f64 = 0.006151;

// ---------------------------------------------------------------------------
// Sidereal time long-term model boundaries (swephlib.c:3460–3463)
// ---------------------------------------------------------------------------

/// Start boundary of the long-term sidereal time model (Julian day number).
pub const SIDT_LTERM_T0: f64 = 2396758.5;
/// End boundary of the long-term sidereal time model (Julian day number).
pub const SIDT_LTERM_T1: f64 = 2469807.5;
/// Sidereal time offset at [`SIDT_LTERM_T0`] (hours).
pub const SIDT_LTERM_OFS0: f64 = 0.000378172 / 15.0;
/// Sidereal time offset at [`SIDT_LTERM_T1`] (hours).
pub const SIDT_LTERM_OFS1: f64 = 0.001385646 / 15.0;

// ---------------------------------------------------------------------------
// Body ID offsets
// ---------------------------------------------------------------------------

/// Number of built-in main planet bodies.
pub const NPLANETS: i32 = 23;
/// Body-ID offset added to asteroid numbers to form an SE body ID.
pub const AST_OFFSET: i32 = 10000;
/// Body-ID offset added to planetary-moon numbers to form an SE body ID.
pub const PLMOON_OFFSET: i32 = 9000;
/// Body-ID offset for fictitious/hypothetical bodies.
pub const FICT_OFFSET: i32 = 40;
/// Alternate body-ID offset for fictitious bodies (one less than [`FICT_OFFSET`]).
pub const FICT_OFFSET_1: i32 = 39;
/// Maximum allowed fictitious body number.
pub const FICT_MAX: i32 = 999;
/// Number of orbital elements stored per fictitious body.
pub const NFICT_ELEM: i32 = 15;

// ---------------------------------------------------------------------------
// Ascendant/MC array indices (for house calculation results)
// ---------------------------------------------------------------------------

/// Array index of the Ascendant in house-calculation results.
pub const ASC: usize = 0;
/// Array index of the Midheaven (MC) in house-calculation results.
pub const MC: usize = 1;
/// Array index of the sidereal time (ARMC) in house-calculation results.
pub const ARMC: usize = 2;
/// Array index of the Vertex in house-calculation results.
pub const VERTEX: usize = 3;
/// Array index of the equatorial ascendant in house-calculation results.
pub const EQUASC: usize = 4;
/// Array index of the first co-ascendant (Walter Koch method) in house-calculation results.
pub const COASC1: usize = 5;
/// Array index of the second co-ascendant (Munkasey method) in house-calculation results.
pub const COASC2: usize = 6;
/// Array index of the polar/horizon ascendant in house-calculation results.
pub const POLASC: usize = 7;
/// Total number of entries in the Ascendant/MC results array.
pub const NASCMC: usize = 8;

// ---------------------------------------------------------------------------
// Tidal acceleration constants (arcsec/cy^2)
// ---------------------------------------------------------------------------

/// Tidal acceleration of the Moon per the JPL DE200 ephemeris (arcsec/century^2).
pub const TIDAL_DE200: f64 = -23.8946;
/// Tidal acceleration of the Moon per the JPL DE403 ephemeris (arcsec/century^2).
pub const TIDAL_DE403: f64 = -25.580;
/// Tidal acceleration of the Moon per the JPL DE404 ephemeris (arcsec/century^2).
pub const TIDAL_DE404: f64 = -25.580;
/// Tidal acceleration of the Moon per the JPL DE405 ephemeris (arcsec/century^2).
pub const TIDAL_DE405: f64 = -25.826;
/// Tidal acceleration of the Moon per the JPL DE406 ephemeris (arcsec/century^2).
pub const TIDAL_DE406: f64 = -25.826;
/// Tidal acceleration of the Moon per the JPL DE421 ephemeris (arcsec/century^2).
pub const TIDAL_DE421: f64 = -25.85;
/// Tidal acceleration of the Moon per the JPL DE422 ephemeris (arcsec/century^2).
pub const TIDAL_DE422: f64 = -25.85;
/// Tidal acceleration of the Moon per the JPL DE430 ephemeris (arcsec/century^2).
pub const TIDAL_DE430: f64 = -25.82;
/// Tidal acceleration of the Moon per the JPL DE431 ephemeris (arcsec/century^2).
pub const TIDAL_DE431: f64 = -25.80;
/// Tidal acceleration of the Moon per the JPL DE441 ephemeris (arcsec/century^2).
pub const TIDAL_DE441: f64 = -25.936;
/// Default tidal acceleration value, currently [`TIDAL_DE431`] (arcsec/century^2).
pub const TIDAL_DEFAULT: f64 = TIDAL_DE431;
/// Sentinel requesting automatic tidal-acceleration selection based on the active ephemeris.
pub const TIDAL_AUTOMATIC: f64 = 999999.0;
/// Fixed tidal acceleration value of -26.0 arcsec/century^2 (older convention).
pub const TIDAL_26: f64 = -26.0;
/// Tidal acceleration of the Moon per Stephenson et al. (2016) (arcsec/century^2).
pub const TIDAL_STEPHENSON_2016: f64 = -25.85;
/// Sentinel requesting automatic Delta T computation.
pub const DELTAT_AUTOMATIC: f64 = -1e-10;

// ---------------------------------------------------------------------------
// Sentinel values
// ---------------------------------------------------------------------------

/// Sentinel Julian day value marking an invalid/uninitialized time.
pub const TJD_INVALID: f64 = 99999999.0;
/// Maximum length of a star name buffer (bytes).
pub const MAX_STNAME: usize = 256;

// ---------------------------------------------------------------------------
// Ephemeris validity ranges (Julian Day numbers)
// ---------------------------------------------------------------------------

/// Start of the Moshier planetary ephemeris validity range (Julian day number).
pub const MOSHPLEPH_START: f64 = 625000.5;
/// End of the Moshier planetary ephemeris validity range (Julian day number).
pub const MOSHPLEPH_END: f64 = 2818000.5;
/// Start of the Moshier lunar ephemeris validity range (Julian day number).
pub const MOSHLUEPH_START: f64 = 625000.5;
/// End of the Moshier lunar ephemeris validity range (Julian day number).
pub const MOSHLUEPH_END: f64 = 2818000.5;
/// Start of the Moshier lunar node/apogee ephemeris validity range (Julian day number).
pub const MOSHNDEPH_START: f64 = -3100015.5;
/// End of the Moshier lunar node/apogee ephemeris validity range (Julian day number).
pub const MOSHNDEPH_END: f64 = 8000016.5;
/// Central-difference interval for planet speed computation (days).
pub const PLAN_SPEED_INTV: f64 = 0.0001;
/// Central-difference interval for fixed-star speed computation (days).
pub const FIXSTAR_DT: f64 = PLAN_SPEED_INTV * 0.1; // = 0.00001 days
/// Central-difference interval for Moon speed computation (days).
pub const MOON_SPEED_INTV: f64 = 0.00005;
/// Central-difference interval for light-deflection speed computation (days).
pub const DEFL_SPEED_INTV: f64 = 0.0000005;
/// Central-difference interval for nutation speed computation (days).
pub const NUT_SPEED_INTV: f64 = 0.0001;
/// Central-difference interval for mean node/apogee speed computation (days).
pub const MEAN_NODE_SPEED_INTV: f64 = 0.001;
/// Central-difference interval for the osculating node/apogee with the JPL/Swiss
/// moon (C `NODE_CALC_INTV`, sweph.c). Small because these backends' moon is smooth.
pub const NODE_CALC_INTV: f64 = 0.0001;
/// Wider interval for the Moshier moon (C `NODE_CALC_INTV_MOSH`): the Moshier
/// moon's short-period terms make the node/apogee oscillate wildly within small
/// intervals, so a coarser finite difference is used.
pub const NODE_CALC_INTV_MOSH: f64 = 0.1;
/// Distance threshold (AU) above which `SE_NODBIT_OSCU_BAR` computes the
/// osculating ellipse about the barycenter rather than the Sun (C's hardcoded
/// `x[2] > 6` test in `swe_nod_aps`, swecl.c:5256 — sits between Jupiter ~5.2 AU
/// and Saturn ~9.5 AU). Used by `orbit.rs` (PNOC 6) and `nodaps.rs`.
pub const OSCU_BAR_DISTANCE_THRESHOLD_AU: f64 = 6.0;

/// Sun-mass / planet-mass ratios, indexed 0..8 = Mercury, Venus, EMB, Mars,
/// Jupiter, Saturn, Uranus, Neptune, Pluto (C `plmass[9]`, swecl.c:5040-5050).
/// Used by `swe_nod_aps`'s osculating branch and `swe_get_orbital_elements`
/// (`orbit.rs`, PNOC 6) — lives here because both app modules need it and they
/// must not import each other.
pub const PLMASS: [f64; 9] = [
    6023600.0,   // Mercury
    408523.719,  // Venus
    328900.5,    // Earth and Moon (EMB)
    3098703.59,  // Mars
    1047.348644, // Jupiter
    3497.9018,   // Saturn
    22902.98,    // Uranus
    19412.26,    // Neptune
    136566000.0, // Pluto
];

/// Maps a body number (indexed directly by the raw `SE_*` id, NOT a contiguous
/// planet count) to a row index into the VSOP mean-element tables
/// (`EL_NODE`/`EL_PERI`/… rows 0..7 = Mercury..Neptune) and, reused, into
/// [`PLMASS`] (rows 0..8, same order + Pluto=8). C `ipl_to_elem[15]`
/// (swecl.c:5052). Transcribed verbatim including its quirks:
/// `IPL_TO_ELEM[0]`=2 (Sun→Earth's row), and `IPL_TO_ELEM[9]`=0 (Pluto→Mercury's
/// row, a stale C mapping used only for a negligible `plm` perturbation — MUST be
/// preserved bit-for-bit, see c-ref-orbital-elements.md's Pluto note).
/// `SE_EARTH`=14 sits far outside the planet range.
pub const IPL_TO_ELEM: [usize; 15] = [2, 0, 0, 1, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 2];

/// Reference epoch for the mean lunar node Gregorian-date correction (Julian day number).
pub const CORR_MNODE_JD_T0GREG: f64 = -3063616.5;
/// Start of the JPL DE431 ephemeris validity range (Julian day number).
pub const JPL_DE431_START: f64 = -3027215.5;
/// End of the JPL DE431 ephemeris validity range (Julian day number).
pub const JPL_DE431_END: f64 = 7930192.5;
/// Start of Chiron's valid ephemeris range (Julian day number).
pub const CHIRON_START: f64 = 1967601.5;
/// End of Chiron's valid ephemeris range (Julian day number).
pub const CHIRON_END: f64 = 3419437.5;
/// Start of Pholus's valid ephemeris range (Julian day number).
pub const PHOLUS_START: f64 = 640648.5;
/// End of Pholus's valid ephemeris range (Julian day number).
pub const PHOLUS_END: f64 = 4390617.5;

// ---------------------------------------------------------------------------
// Solar system plane constants
// ---------------------------------------------------------------------------

/// Ascending node of the solar system invariable plane on the ecliptic of J2000 (radians).
pub const SSY_PLANE_NODE_E2000: f64 = 107.582569 * DEGTORAD;
/// Ascending node of the solar system invariable plane on the mean ecliptic of date (radians).
pub const SSY_PLANE_NODE: f64 = 107.58883388 * DEGTORAD;
/// Inclination of the solar system invariable plane to the ecliptic (radians).
pub const SSY_PLANE_INCL: f64 = 1.578701 * DEGTORAD;

// ---------------------------------------------------------------------------
// Sun radius (angular, in radians)
// ---------------------------------------------------------------------------

/// Mean angular radius of the Sun as seen from Earth (radians).
pub const SUN_RADIUS: f64 = 959.63 / 3600.0 * DEGTORAD;

// ---------------------------------------------------------------------------
// Default ephemeris file names
// ---------------------------------------------------------------------------

/// Default JPL DE431 ephemeris file name.
pub const FNAME_DE431: &str = "de431.eph";
/// Default JPL DE406 ephemeris file name.
pub const FNAME_DE406: &str = "de406.eph";
/// Default ephemeris file name, currently [`FNAME_DE431`].
pub const FNAME_DFT: &str = FNAME_DE431;
/// Default fixed-star catalog file name.
pub const STARFILE: &str = "sefstars.txt";
/// Default asteroid name catalog file name.
pub const ASTNAMFILE: &str = "seasnam.txt";
/// Default fictitious/hypothetical body orbital-elements file name.
pub const FICTFILE: &str = "seorbel.txt";
/// Default Swiss Ephemeris data file suffix (without the leading dot).
pub const FILE_SUFFIX: &str = "se1";

// ---------------------------------------------------------------------------
// Planetary diameters (meters), indexed by raw body ID 0-20
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Date/time constants
// ---------------------------------------------------------------------------

/// Julian day number for 1972 Jan 1, 0h UTC (start of the leap-second era).
pub const J1972: f64 = 2441317.5;
/// Initial number of leap seconds at [`J1972`].
pub const NLEAP_INIT: i32 = 10;

// ---------------------------------------------------------------------------
// Planetary diameters (meters), indexed by raw body ID 0-20
// ---------------------------------------------------------------------------

/// Planetary (and Sun/Moon/asteroid) diameters (meters), indexed by raw body ID 0-20.
pub const PLANETARY_DIAMETERS: [f64; 21] = [
    1392000000.0,     // 0  Sun
    3475000.0,        // 1  Moon
    2439400.0 * 2.0,  // 2  Mercury
    6051800.0 * 2.0,  // 3  Venus
    3389500.0 * 2.0,  // 4  Mars
    69911000.0 * 2.0, // 5  Jupiter
    58232000.0 * 2.0, // 6  Saturn
    25362000.0 * 2.0, // 7  Uranus
    24622000.0 * 2.0, // 8  Neptune
    1188300.0 * 2.0,  // 9  Pluto
    0.0,              // 10 Mean Node
    0.0,              // 11 True Node
    0.0,              // 12 Mean Apogee
    0.0,              // 13 Oscu Apogee
    6371008.4 * 2.0,  // 14 Earth
    271370.0,         // 15 Chiron (irregular)
    290000.0,         // 16 Pholus
    939400.0,         // 17 Ceres
    545000.0,         // 18 Pallas
    246596.0,         // 19 Juno
    525400.0,         // 20 Vesta
];
