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
    planet_files: Vec<crate::sweph_file::SwissEphFile>,
    moon_files: Vec<crate::sweph_file::SwissEphFile>,
    jpl_file: Option<crate::jpl::JplFile>,
}

impl Ephemeris {
    pub fn new(config: EphemerisConfig) -> crate::Result<Self> {
        let leap_seconds = Self::build_leap_seconds(&config)?;
        let mut jpl_file: Option<crate::jpl::JplFile> = None;
        let (planet_files, moon_files) = match config.ephemeris_source {
            EphemerisSource::Swiss => {
                let dir = config.ephe_path.as_ref().ok_or_else(|| {
                    Error::FileFormat("ephe_path required for Swisseph".to_string())
                })?;
                let planet = crate::sweph_file::open_ephemeris_files(dir, "sepl")?;
                let moon = crate::sweph_file::open_ephemeris_files(dir, "semo")?;
                if planet.is_empty() || moon.is_empty() {
                    return Err(Error::FileFormat(
                        "no planet or moon ephemeris files found".to_string(),
                    ));
                }
                (planet, moon)
            }
            EphemerisSource::Jpl => {
                let dir = config
                    .ephe_path
                    .as_ref()
                    .ok_or_else(|| Error::FileFormat("ephe_path required for Jpl".to_string()))?;
                let filename = config.jpl_filename.as_deref().unwrap_or("de441.eph");
                let path = dir.join(filename);
                jpl_file = Some(crate::jpl::JplFile::open(&path)?);
                (Vec::new(), Vec::new())
            }
            EphemerisSource::Moshier => (Vec::new(), Vec::new()),
        };
        Ok(Self {
            config,
            leap_seconds,
            planet_files,
            moon_files,
            jpl_file,
        })
    }

    pub fn config(&self) -> &EphemerisConfig {
        &self.config
    }

    pub fn leap_seconds(&self) -> &[i32] {
        &self.leap_seconds
    }

    /// Compute planetary position.
    ///
    /// Unlike the C library, this implementation does not cache computed
    /// positions. Moshier evaluations are sub-microsecond; callers needing
    /// deduplication for repeated same-JD queries should cache externally.
    pub fn calc(&self, jd_tt: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        let flags = crate::calc::plaus_iflag(flags, self.config.ephemeris_source);
        let unsupported = flags & (CalcFlags::TOPOCTR | CalcFlags::SIDEREAL);
        if !unsupported.is_empty() {
            return Err(Error::UnsupportedFlags(unsupported));
        }

        if body == Body::Earth {
            return Ok(CalcResult {
                data: [0.0; 6],
                flags_used: flags,
            });
        }

        if flags.contains(CalcFlags::SPEED3) {
            return self.calc_speed3(jd_tt, body, flags);
        }

        let (xreturn, flags_used) = self.calc_inner(jd_tt, body, flags)?;
        Ok(CalcResult {
            data: Self::extract_for_body(&xreturn, body, flags_used),
            flags_used,
        })
    }

    pub fn calc_ut(&self, jd_ut: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        let dt = crate::deltat::calc_deltat(jd_ut, &self.config);
        self.calc(jd_ut + dt, body, flags)
    }

    fn extract_for_body(xreturn: &[f64; 24], body: Body, flags: CalcFlags) -> [f64; 6] {
        if body == Body::EclipticNutation {
            crate::calc::extract_ecl_nut(
                &[
                    xreturn[0], xreturn[1], xreturn[2], xreturn[3], xreturn[4], xreturn[5],
                ],
                flags,
            )
        } else {
            crate::calc::extract_output(xreturn, flags)
        }
    }

    fn calc_inner(
        &self,
        jd_tt: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<([f64; 24], CalcFlags), Error> {
        let models = &self.config.astro_models;

        if body == Body::EclipticNutation {
            let ecl_nut = crate::calc::calc_ecl_nut(jd_tt, flags, models);
            let mut xreturn = [0.0; 24];
            xreturn[..6].copy_from_slice(&ecl_nut);
            return Ok((xreturn, flags));
        }

        if matches!(body, Body::MeanNode | Body::MeanApogee) {
            if flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
                return Ok(([0.0; 24], flags));
            }
            let xr = match body {
                Body::MeanNode => crate::calc::calc_mean_node(jd_tt, flags, models)?,
                Body::MeanApogee => crate::calc::calc_mean_apogee(jd_tt, flags, models)?,
                _ => unreachable!(),
            };
            return Ok((xr, flags));
        }

        if flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
            return Err(Error::UnsupportedFlags(
                flags & (CalcFlags::HELCTR | CalcFlags::BARYCTR),
            ));
        }

        let eps_j2000 =
            crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);

        match self.config.ephemeris_source {
            EphemerisSource::Swiss => {
                match self.calc_body_sweph(jd_tt, body, &eps_j2000, flags, models) {
                    Ok(xr) => Ok((xr, flags)),
                    Err(Error::BeyondEphemerisLimits { .. }) => {
                        let fallback_flags = (flags & !CalcFlags::SWIEPH) | CalcFlags::MOSEPH;
                        let xr = self.calc_body_moshier(
                            jd_tt,
                            body,
                            &eps_j2000,
                            fallback_flags,
                            models,
                        )?;
                        Ok((xr, fallback_flags))
                    }
                    Err(e) => Err(e),
                }
            }
            EphemerisSource::Jpl => {
                let xr = self.calc_body_jpl(jd_tt, body, &eps_j2000, flags, models)?;
                Ok((xr, flags))
            }
            EphemerisSource::Moshier => {
                let xr = self.calc_body_moshier(jd_tt, body, &eps_j2000, flags, models)?;
                Ok((xr, flags))
            }
        }
    }

    fn calc_body_moshier(
        &self,
        jd_tt: f64,
        body: Body,
        eps_j2000: &crate::types::Epsilon,
        flags: CalcFlags,
        models: &crate::types::AstroModels,
    ) -> Result<[f64; 24], Error> {
        match body {
            Body::Sun => crate::calc::calc_sun(jd_tt, eps_j2000, flags, models),
            Body::Moon => crate::calc::calc_moon(jd_tt, eps_j2000, flags, models),
            Body::Mercury
            | Body::Venus
            | Body::Mars
            | Body::Jupiter
            | Body::Saturn
            | Body::Uranus
            | Body::Neptune
            | Body::Pluto => crate::calc::calc_planet(jd_tt, body, eps_j2000, flags, models),
            _ => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Moshier,
            }),
        }
    }

    fn calc_body_sweph(
        &self,
        jd_tt: f64,
        body: Body,
        eps_j2000: &crate::types::Epsilon,
        flags: CalcFlags,
        models: &crate::types::AstroModels,
    ) -> Result<[f64; 24], Error> {
        match body {
            Body::Sun => crate::calc::calc_sun_sweph(
                jd_tt,
                &self.planet_files,
                &self.moon_files,
                flags,
                models,
            ),
            Body::Moon => crate::calc::calc_moon_sweph(
                jd_tt,
                &self.planet_files,
                &self.moon_files,
                flags,
                models,
            ),
            Body::Mercury
            | Body::Venus
            | Body::Mars
            | Body::Jupiter
            | Body::Saturn
            | Body::Uranus
            | Body::Neptune
            | Body::Pluto => crate::calc::calc_planet_sweph(
                jd_tt,
                body,
                &self.planet_files,
                &self.moon_files,
                eps_j2000,
                flags,
                models,
            ),
            _ => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Swiss,
            }),
        }
    }

    fn calc_body_jpl(
        &self,
        jd_tt: f64,
        body: Body,
        eps_j2000: &crate::types::Epsilon,
        flags: CalcFlags,
        models: &crate::types::AstroModels,
    ) -> Result<[f64; 24], Error> {
        let file = self.jpl_file.as_ref().unwrap();
        match body {
            Body::Sun => crate::calc::calc_sun_jpl(jd_tt, file, flags, models),
            Body::Moon => crate::calc::calc_moon_jpl(jd_tt, file, flags, models),
            Body::Mercury
            | Body::Venus
            | Body::Mars
            | Body::Jupiter
            | Body::Saturn
            | Body::Uranus
            | Body::Neptune
            | Body::Pluto => {
                crate::calc::calc_planet_jpl(jd_tt, body, file, eps_j2000, flags, models)
            }
            _ => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Jpl,
            }),
        }
    }

    fn calc_speed3(&self, jd_tt: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        let dt = crate::calc::speed3_interval(body);
        let inner_flags = flags & !CalcFlags::SPEED3;

        let (mut x0, _) = self.calc_inner(jd_tt - dt, body, inner_flags)?;
        let (mut x2, _) = self.calc_inner(jd_tt + dt, body, inner_flags)?;
        let (mut x1, flags_used) = self.calc_inner(jd_tt, body, inner_flags)?;

        crate::calc::denormalize_positions(&mut x0, &x1, &mut x2);
        crate::calc::calc_speed_3point(&mut x1, &x0, &x2, dt);

        Ok(CalcResult {
            data: Self::extract_for_body(&x1, body, flags | CalcFlags::SPEED),
            flags_used,
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
