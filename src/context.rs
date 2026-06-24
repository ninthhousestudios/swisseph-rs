use std::path::PathBuf;

use crate::flags::CalcFlags;
use crate::types::{AstroModels, SiderealMode};

#[derive(Debug, Clone, Copy)]
pub struct TopoPosition {
    pub longitude: f64,
    pub latitude: f64,
    pub altitude: f64,
}

#[derive(Debug, Clone)]
pub struct EphemerisConfig {
    pub ephe_path: Option<PathBuf>,
    pub jpl_filename: Option<String>,
    pub sidereal_mode: Option<SiderealMode>,
    pub sidereal_t0: f64,
    pub sidereal_ayan_t0: f64,
    pub topographic: Option<TopoPosition>,
    pub astro_models: AstroModels,
    pub tidal_acceleration: Option<f64>,
}

impl Default for EphemerisConfig {
    fn default() -> Self {
        Self {
            ephe_path: None,
            jpl_filename: None,
            sidereal_mode: None,
            sidereal_t0: 0.0,
            sidereal_ayan_t0: 0.0,
            topographic: None,
            astro_models: AstroModels::default(),
            tidal_acceleration: None,
        }
    }
}

pub struct Ephemeris {
    config: EphemerisConfig,
}

impl Ephemeris {
    pub fn new(config: EphemerisConfig) -> crate::Result<Self> {
        Ok(Self { config })
    }

    pub fn config(&self) -> &EphemerisConfig {
        &self.config
    }
}

pub struct CalcResult {
    pub data: [f64; 6],
    pub flags_used: CalcFlags,
}
