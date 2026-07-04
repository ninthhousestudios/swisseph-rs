//! Read-only ephemeris configuration types (`EphemerisConfig`, `TopoPosition`). Leaf module —
//! depends only on `flags`/`types`, so `calc`/`ayanamsa`/etc. can import config without cycling
//! back through `context`.

use std::path::PathBuf;

use crate::flags::SiderealBits;
use crate::types::{AstroModels, EphemerisSource, SiderealMode};

/// Geographic position for topocentric calculations. Replaces C's `swe_set_topo`.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TopoPosition {
    /// Geographic longitude in degrees, east-positive.
    pub longitude: f64,
    /// Geographic latitude in degrees, north-positive.
    pub latitude: f64,
    /// Altitude above sea level in meters.
    pub altitude: f64,
}

/// Read-only configuration for an [`Ephemeris`](crate::Ephemeris) instance.
///
/// Replaces the stateful `swe_set_*` family of C functions with a single immutable struct.
/// `set_sidereal_mode` is implemented in `ayanamsa.rs`, next to the `AYANAMSA` table it resolves.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EphemerisConfig {
    /// Which ephemeris backend to use. Replaces `swe_set_ephe_path` (Swiss/JPL)
    /// or the `SEFLG_MOSEPH` flag. Default: [`Moshier`](EphemerisSource::Moshier).
    pub ephemeris_source: EphemerisSource,
    /// Directory containing `.se1` ephemeris files. Required for `Swiss`/`Jpl` sources.
    /// Replaces `swe_set_ephe_path`.
    pub ephe_path: Option<PathBuf>,
    /// JPL DE filename (e.g. `"de441.eph"`). If `None`, uses the default in `ephe_path`.
    /// Replaces `swe_set_jpl_file`.
    pub jpl_filename: Option<String>,
    /// Sidereal mode for `SEFLG_SIDEREAL` calculations. Replaces `swe_set_sid_mode(mode, ...)`.
    pub sidereal_mode: Option<SiderealMode>,
    /// Reference epoch for user-defined sidereal mode (Julian Day, TT unless `sidereal_t0_is_ut`).
    pub sidereal_t0: f64,
    /// Initial ayanamsa value at `sidereal_t0`, degrees.
    pub sidereal_ayan_t0: f64,
    /// Sidereal projection modifier bits. See [`SiderealBits`].
    pub sidereal_bits: SiderealBits,
    /// If `true`, `sidereal_t0` is UT rather than TT.
    pub sidereal_t0_is_ut: bool,
    /// Observer position for `SEFLG_TOPOCTR`. Replaces `swe_set_topo`.
    pub topographic: Option<TopoPosition>,
    /// Astronomical model overrides. Replaces `swe_set_astro_models`.
    pub astro_models: AstroModels,
    /// Override tidal acceleration (arcsec/century^2). `None` = auto-derive from ephemeris.
    /// Replaces `swe_set_tid_acc`.
    pub tidal_acceleration: Option<f64>,
    /// User-defined Delta T value (days). When `Some`, bypasses all Delta T models.
    /// Replaces `swe_set_delta_t_userdef`.
    pub delta_t_userdef: Option<f64>,
    /// Extra leap-second years beyond the built-in table.
    pub extra_leap_seconds: Vec<i32>,
    /// Path to a custom leap-seconds file (overrides built-in table entirely).
    pub leap_seconds_file: Option<PathBuf>,
    /// MPC numbers of asteroids whose `.se1` files should be opened. Stateless design
    /// requires declaring which asteroids are needed upfront.
    pub asteroid_numbers: Vec<i32>,
    /// Raw planetary-moon/COB ids per `ephe/sat/plmolist.txt` (9401–9999;
    /// COB entries are 9n99). Both `Body::PlanetMoon(id)` calc AND
    /// `CalcFlags::CENTER_BODY` on Jupiter..Pluto (which resolves to the 9n99
    /// COB id) require the id listed here — stateless design, no lazy file opening.
    pub planet_moon_numbers: Vec<i32>,
}

impl Default for EphemerisConfig {
    fn default() -> Self {
        Self {
            ephemeris_source: EphemerisSource::Moshier,
            ephe_path: None,
            jpl_filename: None,
            sidereal_mode: None,
            sidereal_t0: 0.0,
            sidereal_ayan_t0: 0.0,
            sidereal_bits: SiderealBits::empty(),
            sidereal_t0_is_ut: false,
            topographic: None,
            astro_models: AstroModels::default(),
            tidal_acceleration: None,
            delta_t_userdef: None,
            extra_leap_seconds: Vec::new(),
            leap_seconds_file: None,
            asteroid_numbers: Vec::new(),
            planet_moon_numbers: Vec::new(),
        }
    }
}
