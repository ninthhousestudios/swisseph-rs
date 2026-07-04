//! Read-only ephemeris configuration types (`EphemerisConfig`, `TopoPosition`). Leaf module —
//! depends only on `flags`/`types`, so `calc`/`ayanamsa`/etc. can import config without cycling
//! back through `context`.

use std::path::PathBuf;

use crate::flags::SiderealBits;
use crate::types::{AstroModels, EphemerisSource, SiderealMode};

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TopoPosition {
    pub longitude: f64,
    pub latitude: f64,
    pub altitude: f64,
}

/// `set_sidereal_mode` is implemented in `ayanamsa.rs`, next to the `AYANAMSA` table it resolves.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EphemerisConfig {
    pub ephemeris_source: EphemerisSource,
    pub ephe_path: Option<PathBuf>,
    pub jpl_filename: Option<String>,
    pub sidereal_mode: Option<SiderealMode>,
    pub sidereal_t0: f64,
    pub sidereal_ayan_t0: f64,
    pub sidereal_bits: SiderealBits,
    pub sidereal_t0_is_ut: bool,
    pub topographic: Option<TopoPosition>,
    pub astro_models: AstroModels,
    pub tidal_acceleration: Option<f64>,
    pub delta_t_userdef: Option<f64>,
    pub extra_leap_seconds: Vec<i32>,
    pub leap_seconds_file: Option<PathBuf>,
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
