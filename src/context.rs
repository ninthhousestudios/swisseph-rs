use std::fs;
use std::path::PathBuf;

use crate::date::LEAP_SECONDS;
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::types::{AstroModels, Body, DeltaT, EphemerisSource, JdUt1, SiderealMode};

#[derive(Debug, Clone, Copy)]
pub struct TopoPosition {
    pub longitude: f64,
    pub latitude: f64,
    pub altitude: f64,
}

#[derive(Debug, Clone)]
pub struct EphemerisConfig {
    pub ephemeris_source: EphemerisSource,
    pub ephe_path: Option<PathBuf>,
    pub jpl_filename: Option<String>,
    pub sidereal_mode: Option<SiderealMode>,
    pub sidereal_t0: f64,
    pub sidereal_ayan_t0: f64,
    pub topographic: Option<TopoPosition>,
    pub astro_models: AstroModels,
    pub tidal_acceleration: Option<f64>,
    pub extra_leap_seconds: Vec<i32>,
    pub leap_seconds_file: Option<PathBuf>,
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
            topographic: None,
            astro_models: AstroModels::default(),
            tidal_acceleration: None,
            extra_leap_seconds: Vec::new(),
            leap_seconds_file: None,
        }
    }
}

pub struct Ephemeris {
    config: EphemerisConfig,
    leap_seconds: Vec<i32>,
}

impl Ephemeris {
    pub fn new(config: EphemerisConfig) -> crate::Result<Self> {
        let leap_seconds = Self::build_leap_seconds(&config)?;
        Ok(Self {
            config,
            leap_seconds,
        })
    }

    pub fn config(&self) -> &EphemerisConfig {
        &self.config
    }

    pub fn leap_seconds(&self) -> &[i32] {
        &self.leap_seconds
    }

    pub fn calc(&self, jd_tt: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        let flags = crate::calc::plaus_iflag(flags);
        let models = &self.config.astro_models;
        let eps_j2000 =
            crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);

        if body == Body::Earth {
            return Ok(CalcResult {
                data: [0.0; 6],
                flags_used: flags,
            });
        }

        let xreturn = match body {
            Body::Sun => crate::calc::calc_sun(jd_tt, &eps_j2000, flags, models)?,
            Body::Moon => crate::calc::calc_moon(jd_tt, &eps_j2000, flags, models)?,
            _ => crate::calc::calc_planet(jd_tt, body, &eps_j2000, flags, models)?,
        };

        Ok(CalcResult {
            data: crate::calc::extract_output(&xreturn, flags),
            flags_used: flags,
        })
    }

    fn build_leap_seconds(config: &EphemerisConfig) -> crate::Result<Vec<i32>> {
        let last_hardcoded = *LEAP_SECONDS.last().unwrap();
        let mut table: Vec<i32> = LEAP_SECONDS.to_vec();
        // Merge extra entries from config
        for &entry in &config.extra_leap_seconds {
            if entry > last_hardcoded && !table.contains(&entry) {
                table.push(entry);
            }
        }
        // Parse file if provided
        if let Some(path) = &config.leap_seconds_file {
            match fs::read_to_string(path) {
                Ok(contents) => {
                    for line in contents.lines() {
                        let trimmed = line.trim();
                        if trimmed.is_empty() || trimmed.starts_with('#') {
                            continue;
                        }
                        if let Ok(ndat) = trimmed.parse::<i32>() {
                            if ndat > last_hardcoded && !table.contains(&ndat) {
                                table.push(ndat);
                            }
                        }
                    }
                }
                Err(_) if !path.exists() => {} // silently ignore missing file, matching C behavior
                Err(_) => return Err(Error::FileNotFound(path.clone())),
            }
        }
        table.sort_unstable();
        Ok(table)
    }
}

impl DeltaT for Ephemeris {
    fn delta_t(&self, jd_ut: JdUt1) -> f64 {
        crate::deltat::calc_deltat(jd_ut.0, &self.config)
    }
}

pub struct CalcResult {
    pub data: [f64; 6],
    pub flags_used: CalcFlags,
}
