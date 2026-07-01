use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

pub const AUNIT: f64 = 1.49597870700e11;
pub const CLIGHT: f64 = 2.99792458e8;
pub const HELGRAVCONST: f64 = 1.32712440017987e20;
pub const GEOGCONST: f64 = 3.98600448e14;
pub const KGAUSS: f64 = 0.01720209895;
pub const EARTH_RADIUS: f64 = 6378136.6;
pub const EARTH_OBLATENESS: f64 = 1.0 / 298.25642;
pub const EARTH_ROT_SPEED: f64 = 7.2921151467e-5 * 86400.0;
pub const SUN_EARTH_MRAT: f64 = 332946.050895;
pub const EARTH_MOON_MRAT: f64 = 1.0 / 0.0123000383;
pub const MOON_MEAN_DIST: f64 = 384400000.0;
pub const MOON_MEAN_INCL: f64 = 5.1453964;
pub const MOON_MEAN_ECC: f64 = 0.054900489;
pub const LIGHTTIME_AUNIT: f64 = 499.0047838362 / 3600.0 / 24.0;
pub const PARSEC_TO_AUNIT: f64 = 206264.8062471;
pub const LAPSE_RATE: f64 = 0.0065;
pub const KM_S_TO_AU_CTY: f64 = 21.095;

// Valid observer-altitude range for rise/set (sweph.h:198-199, SEI_ECL_GEOALT_MIN/_MAX), meters.
pub const RISE_SET_GEOALT_MIN: f64 = -500.0;
pub const RISE_SET_GEOALT_MAX: f64 = 25000.0;

// ---------------------------------------------------------------------------
// Unit conversions
// ---------------------------------------------------------------------------

pub const DEGTORAD: f64 = PI / 180.0;
pub const RADTODEG: f64 = 180.0 / PI;
pub const CSTORAD: f64 = DEGTORAD / 360000.0;
pub const RADTOCS: f64 = RADTODEG * 360000.0;
pub const CS2DEG: f64 = 1.0 / 360000.0;
pub const TWOPI: f64 = 2.0 * PI;
pub const AUNIT_TO_KM: f64 = 149597870.700;
pub const AUNIT_TO_LIGHTYEAR: f64 = 1.0 / 63241.07708427;
pub const AUNIT_TO_PARSEC: f64 = 1.0 / 206264.8062471;
// Verbatim arcsec-to-radian constant from the C source; full digits preserved.
#[allow(clippy::excessive_precision)]
pub const STR: f64 = 4.8481368110953599359e-6;

// ---------------------------------------------------------------------------
// Centisecond angle constants
// ---------------------------------------------------------------------------

pub const DEG: i32 = 360000;
pub const DEG15: i32 = 15 * DEG;
pub const DEG30: i32 = 30 * DEG;
pub const DEG60: i32 = 60 * DEG;
pub const DEG90: i32 = 90 * DEG;
pub const DEG120: i32 = 120 * DEG;
pub const DEG150: i32 = 150 * DEG;
pub const DEG180: i32 = 180 * DEG;
pub const DEG270: i32 = 270 * DEG;
pub const DEG360: i32 = 360 * DEG;

// ---------------------------------------------------------------------------
// Reference epochs (Julian Day numbers)
// ---------------------------------------------------------------------------

pub const J2000: f64 = 2451545.0;
pub const B1950: f64 = 2433282.42345905;
pub const J1900: f64 = 2415020.0;
pub const B1850: f64 = 2396758.2035810;

// ---------------------------------------------------------------------------
// Precession model validity ranges (Julian centuries from J2000)
// ---------------------------------------------------------------------------

pub const PREC_IAU_1976_CTIES: f64 = 2.0;
pub const PREC_IAU_2000_CTIES: f64 = 2.0;
pub const PREC_IAU_2006_CTIES: f64 = 75.0;

// ---------------------------------------------------------------------------
// JPL Horizons epoch constants
// ---------------------------------------------------------------------------

pub const DPSI_DEPS_IAU1980_TJD0_HORIZONS: f64 = 2437684.5;
pub const DPSI_IAU1980_TJD0: f64 = 0.064284;
pub const DEPS_IAU1980_TJD0: f64 = 0.006151;

// ---------------------------------------------------------------------------
// Sidereal time long-term model boundaries (swephlib.c:3460–3463)
// ---------------------------------------------------------------------------

pub const SIDT_LTERM_T0: f64 = 2396758.5;
pub const SIDT_LTERM_T1: f64 = 2469807.5;
pub const SIDT_LTERM_OFS0: f64 = 0.000378172 / 15.0;
pub const SIDT_LTERM_OFS1: f64 = 0.001385646 / 15.0;

// ---------------------------------------------------------------------------
// Body ID offsets
// ---------------------------------------------------------------------------

pub const NPLANETS: i32 = 23;
pub const AST_OFFSET: i32 = 10000;
pub const PLMOON_OFFSET: i32 = 9000;
pub const COMET_OFFSET: i32 = 1000;
pub const FICT_OFFSET: i32 = 40;
pub const FICT_OFFSET_1: i32 = 39;
pub const FICT_MAX: i32 = 999;
pub const NFICT_ELEM: i32 = 15;

// ---------------------------------------------------------------------------
// Ascendant/MC array indices (for house calculation results)
// ---------------------------------------------------------------------------

pub const ASC: usize = 0;
pub const MC: usize = 1;
pub const ARMC: usize = 2;
pub const VERTEX: usize = 3;
pub const EQUASC: usize = 4;
pub const COASC1: usize = 5;
pub const COASC2: usize = 6;
pub const POLASC: usize = 7;
pub const NASCMC: usize = 8;

// ---------------------------------------------------------------------------
// Tidal acceleration constants (arcsec/cy^2)
// ---------------------------------------------------------------------------

pub const TIDAL_DE200: f64 = -23.8946;
pub const TIDAL_DE403: f64 = -25.580;
pub const TIDAL_DE404: f64 = -25.580;
pub const TIDAL_DE405: f64 = -25.826;
pub const TIDAL_DE406: f64 = -25.826;
pub const TIDAL_DE421: f64 = -25.85;
pub const TIDAL_DE422: f64 = -25.85;
pub const TIDAL_DE430: f64 = -25.82;
pub const TIDAL_DE431: f64 = -25.80;
pub const TIDAL_DE441: f64 = -25.936;
pub const TIDAL_DEFAULT: f64 = TIDAL_DE431;
pub const TIDAL_AUTOMATIC: f64 = 999999.0;
pub const TIDAL_26: f64 = -26.0;
pub const TIDAL_STEPHENSON_2016: f64 = -25.85;
pub const DELTAT_AUTOMATIC: f64 = -1e-10;

// ---------------------------------------------------------------------------
// Sentinel values
// ---------------------------------------------------------------------------

pub const TJD_INVALID: f64 = 99999999.0;
pub const MAX_STNAME: usize = 256;

// ---------------------------------------------------------------------------
// Ephemeris validity ranges (Julian Day numbers)
// ---------------------------------------------------------------------------

pub const MOSHPLEPH_START: f64 = 625000.5;
pub const MOSHPLEPH_END: f64 = 2818000.5;
pub const MOSHLUEPH_START: f64 = 625000.5;
pub const MOSHLUEPH_END: f64 = 2818000.5;
pub const MOSHNDEPH_START: f64 = -3100015.5;
pub const MOSHNDEPH_END: f64 = 8000016.5;
pub const PLAN_SPEED_INTV: f64 = 0.0001;
pub const FIXSTAR_DT: f64 = PLAN_SPEED_INTV * 0.1; // = 0.00001 days
pub const MOON_SPEED_INTV: f64 = 0.00005;
pub const DEFL_SPEED_INTV: f64 = 0.0000005;
pub const NUT_SPEED_INTV: f64 = 0.0001;
pub const MEAN_NODE_SPEED_INTV: f64 = 0.001;
pub const CORR_MNODE_JD_T0GREG: f64 = -3063616.5;
pub const JPL_DE431_START: f64 = -3027215.5;
pub const JPL_DE431_END: f64 = 7930192.5;
pub const CHIRON_START: f64 = 1967601.5;
pub const CHIRON_END: f64 = 3419437.5;
pub const PHOLUS_START: f64 = 640648.5;
pub const PHOLUS_END: f64 = 4390617.5;

// ---------------------------------------------------------------------------
// Solar system plane constants
// ---------------------------------------------------------------------------

pub const SSY_PLANE_NODE_E2000: f64 = 107.582569 * DEGTORAD;
pub const SSY_PLANE_NODE: f64 = 107.58883388 * DEGTORAD;
pub const SSY_PLANE_INCL: f64 = 1.578701 * DEGTORAD;

// ---------------------------------------------------------------------------
// Sun radius (angular, in radians)
// ---------------------------------------------------------------------------

pub const SUN_RADIUS: f64 = 959.63 / 3600.0 * DEGTORAD;

// ---------------------------------------------------------------------------
// Default ephemeris file names
// ---------------------------------------------------------------------------

pub const FNAME_DE431: &str = "de431.eph";
pub const FNAME_DE406: &str = "de406.eph";
pub const FNAME_DFT: &str = FNAME_DE431;
pub const STARFILE: &str = "sefstars.txt";
pub const ASTNAMFILE: &str = "seasnam.txt";
pub const FICTFILE: &str = "seorbel.txt";
pub const FILE_SUFFIX: &str = "se1";

// ---------------------------------------------------------------------------
// Planetary diameters (meters), indexed by raw body ID 0-20
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Date/time constants
// ---------------------------------------------------------------------------

pub const J1972: f64 = 2441317.5;
pub const NLEAP_INIT: i32 = 10;

// ---------------------------------------------------------------------------
// Planetary diameters (meters), indexed by raw body ID 0-20
// ---------------------------------------------------------------------------

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
