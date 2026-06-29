use std::fs;
use std::path::PathBuf;

use crate::date::LEAP_SECONDS;
use crate::error::Error;
use crate::flags::{CalcFlags, SiderealBits};
use crate::types::{
    AstroModels, Body, DeltaT, EphemerisSource, JdUt1, NutationModel, PrecessionModel, SiderealMode,
};

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
    pub sidereal_bits: SiderealBits,
    pub sidereal_t0_is_ut: bool,
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
            sidereal_bits: SiderealBits::empty(),
            sidereal_t0_is_ut: false,
            topographic: None,
            astro_models: AstroModels::default(),
            tidal_acceleration: None,
            extra_leap_seconds: Vec::new(),
            leap_seconds_file: None,
        }
    }
}

impl EphemerisConfig {
    pub fn set_sidereal_mode(&mut self, mut sid_mode: i32, t0: f64, ayan_t0: f64) {
        if sid_mode < 0 {
            sid_mode = 0;
        }
        let mut index = (sid_mode % 256) as usize;

        let mut bits = if matches!(index, 18 | 19 | 20 | 34) {
            SiderealBits::ECL_T0
        } else if matches!(
            index,
            17 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 35 | 36 | 39 | 40
        ) {
            SiderealBits::empty()
        } else {
            SiderealBits::from_bits_truncate((sid_mode as u32) & !0xFF_u32)
        };

        if index >= 47 && index != 255 {
            index = 0;
            bits = SiderealBits::empty();
        }

        if index == 255 {
            self.sidereal_t0 = t0;
            self.sidereal_ayan_t0 = ayan_t0;
            self.sidereal_t0_is_ut = bits.contains(SiderealBits::USER_UT);
        } else {
            let a = crate::ayanamsa::AYANAMSA[index];
            self.sidereal_t0 = a.t0;
            self.sidereal_ayan_t0 = a.ayan_t0;
            self.sidereal_t0_is_ut = a.t0_is_ut;
        }

        if bits.contains(SiderealBits::PREC_ORIG) && index != 255 {
            let prec_offset = crate::ayanamsa::AYANAMSA[index].prec_offset;
            if prec_offset > 0 {
                let prec_model = match prec_offset {
                    1 => PrecessionModel::IAU1976,
                    11 => PrecessionModel::Newcomb,
                    _ => unreachable!(),
                };
                self.astro_models.prec_longterm = prec_model;
                self.astro_models.prec_shortterm = prec_model;
                self.astro_models.nutation = match prec_offset {
                    11 => NutationModel::Woolard,
                    1 => NutationModel::IAU1980,
                    _ => unreachable!(),
                };
            }
        }

        self.sidereal_mode =
            Some(SiderealMode::try_from(index as i32).unwrap_or(SiderealMode::FaganBradley));
        self.sidereal_bits = bits;
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
        let unsupported = flags & CalcFlags::TOPOCTR;
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

        let (mut xreturn, x2000, flags_used) = self.calc_inner(jd_tt, body, flags)?;
        if flags.contains(CalcFlags::SIDEREAL) && body != Body::EclipticNutation {
            self.apply_sidereal(&mut xreturn, &x2000, jd_tt, flags_used)?;
        }
        Ok(CalcResult {
            data: Self::extract_for_body(&xreturn, body, flags_used),
            flags_used,
        })
    }

    pub fn calc_ut(&self, jd_ut: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        let dt = crate::deltat::calc_deltat(jd_ut, &self.config);
        self.calc(jd_ut + dt, body, flags)
    }

    /// Ayanamsa at `jd_tt` (TT), with nutation added unless `NONUT` is set.
    pub fn get_ayanamsa_ex(&self, jd_tt: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::ayanamsa::get_ayanamsa_ex_nut(&self.config, jd_tt, flags, &self.config.astro_models)
    }

    /// Ayanamsa at `jd_ut` (UT), with nutation added unless `NONUT` is set.
    pub fn get_ayanamsa_ut(&self, jd_ut: f64, flags: CalcFlags) -> Result<f64, Error> {
        let dt = crate::deltat::calc_deltat(jd_ut, &self.config);
        self.get_ayanamsa_ex(jd_ut + dt, flags)
    }

    /// Legacy ayanamsa accessor (no nutation) at `jd_tt` (TT).
    pub fn get_ayanamsa(&self, jd_tt: f64) -> Result<f64, Error> {
        crate::ayanamsa::get_ayanamsa_ex(
            &self.config,
            jd_tt,
            CalcFlags::empty(),
            &self.config.astro_models,
        )
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
    ) -> Result<([f64; 24], [f64; 6], CalcFlags), Error> {
        let models = &self.config.astro_models;

        if body == Body::EclipticNutation {
            let ecl_nut = crate::calc::calc_ecl_nut(jd_tt, flags, models);
            let mut xreturn = [0.0; 24];
            xreturn[..6].copy_from_slice(&ecl_nut);
            return Ok((xreturn, [0.0; 6], flags));
        }

        if matches!(body, Body::MeanNode | Body::MeanApogee) {
            if flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
                return Ok(([0.0; 24], [0.0; 6], flags));
            }
            let xr = match body {
                Body::MeanNode => crate::calc::calc_mean_node(jd_tt, flags, models)?,
                Body::MeanApogee => crate::calc::calc_mean_apogee(jd_tt, flags, models)?,
                _ => unreachable!(),
            };
            return Ok((xr, [0.0; 6], flags));
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
                    Ok((xr, x2000)) => Ok((xr, x2000, flags)),
                    Err(Error::BeyondEphemerisLimits { .. }) => {
                        let fallback_flags = (flags & !CalcFlags::SWIEPH) | CalcFlags::MOSEPH;
                        let (xr, x2000) = self.calc_body_moshier(
                            jd_tt,
                            body,
                            &eps_j2000,
                            fallback_flags,
                            models,
                        )?;
                        Ok((xr, x2000, fallback_flags))
                    }
                    Err(e) => Err(e),
                }
            }
            EphemerisSource::Jpl => {
                let (xr, x2000) = self.calc_body_jpl(jd_tt, body, &eps_j2000, flags, models)?;
                Ok((xr, x2000, flags))
            }
            EphemerisSource::Moshier => {
                let (xr, x2000) = self.calc_body_moshier(jd_tt, body, &eps_j2000, flags, models)?;
                Ok((xr, x2000, flags))
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
    ) -> Result<([f64; 24], [f64; 6]), Error> {
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
    ) -> Result<([f64; 24], [f64; 6]), Error> {
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
    ) -> Result<([f64; 24], [f64; 6]), Error> {
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

        let (mut x0, _, _) = self.calc_inner(jd_tt - dt, body, inner_flags)?;
        let (mut x2, _, _) = self.calc_inner(jd_tt + dt, body, inner_flags)?;
        let (mut x1, x2000, flags_used) = self.calc_inner(jd_tt, body, inner_flags)?;

        crate::calc::denormalize_positions(&mut x0, &x1, &mut x2);
        crate::calc::calc_speed_3point(&mut x1, &x0, &x2, dt);

        if flags.contains(CalcFlags::SIDEREAL) && body != Body::EclipticNutation {
            self.apply_sidereal(&mut x1, &x2000, jd_tt, flags_used | CalcFlags::SPEED)?;
        }

        Ok(CalcResult {
            data: Self::extract_for_body(&x1, body, flags | CalcFlags::SPEED),
            flags_used,
        })
    }

    fn apply_sidereal(
        &self,
        xreturn: &mut [f64; 24],
        x2000: &[f64; 6],
        jd_tt: f64,
        flags: CalcFlags,
    ) -> Result<(), Error> {
        use crate::constants::RADTODEG;
        use crate::math::cartesian_to_polar_with_speed;

        let bits = self.config.sidereal_bits;
        let models = &self.config.astro_models;
        let has_speed = flags.contains(CalcFlags::SPEED);
        let has_meaningful_x2000 = *x2000 != [0.0f64; 6];

        if has_meaningful_x2000 && bits.contains(SiderealBits::ECL_T0) {
            let (xecl, xequ) = crate::ayanamsa::trop_ra2sid_lon(x2000, &self.config, models, flags);

            xreturn[6..12].copy_from_slice(&xecl);
            xreturn[18..24].copy_from_slice(&xequ);

            // Recompute ecliptic polar [0..6] from new Cartesian [6..12]
            let ecl_pol = cartesian_to_polar_with_speed(xecl);
            xreturn[0] = ecl_pol[0] * RADTODEG;
            xreturn[1] = ecl_pol[1] * RADTODEG;
            xreturn[2] = ecl_pol[2];
            xreturn[3] = if has_speed {
                ecl_pol[3] * RADTODEG
            } else {
                0.0
            };
            xreturn[4] = if has_speed {
                ecl_pol[4] * RADTODEG
            } else {
                0.0
            };
            xreturn[5] = if has_speed { ecl_pol[5] } else { 0.0 };

            // Recompute equatorial polar [12..18] from new Cartesian [18..24]
            let equ_pol = cartesian_to_polar_with_speed(xequ);
            xreturn[12] = equ_pol[0] * RADTODEG;
            xreturn[13] = equ_pol[1] * RADTODEG;
            xreturn[14] = equ_pol[2];
            xreturn[15] = if has_speed {
                equ_pol[3] * RADTODEG
            } else {
                0.0
            };
            xreturn[16] = if has_speed {
                equ_pol[4] * RADTODEG
            } else {
                0.0
            };
            xreturn[17] = if has_speed { equ_pol[5] } else { 0.0 };
        } else if has_meaningful_x2000 && bits.contains(SiderealBits::SSY_PLANE) {
            let xecl = crate::ayanamsa::trop_ra2sid_lon_sosy(x2000, &self.config, models, flags);

            xreturn[6..12].copy_from_slice(&xecl);

            // Recompute ecliptic polar [0..6] from new Cartesian [6..12]
            let ecl_pol = cartesian_to_polar_with_speed(xecl);
            xreturn[0] = ecl_pol[0] * RADTODEG;
            xreturn[1] = ecl_pol[1] * RADTODEG;
            xreturn[2] = ecl_pol[2];
            xreturn[3] = if has_speed {
                ecl_pol[3] * RADTODEG
            } else {
                0.0
            };
            xreturn[4] = if has_speed {
                ecl_pol[4] * RADTODEG
            } else {
                0.0
            };
            xreturn[5] = if has_speed { ecl_pol[5] } else { 0.0 };

            // Recompute ecliptic Cartesian [6..12] from polar (already done above)
            // Leave equatorial [12..24] untouched (matches C Branch 2)
        } else {
            // Default branch: ayanamsa subtraction on ecliptic polar
            let daya =
                crate::ayanamsa::get_ayanamsa_with_speed(&self.config, jd_tt, flags, models)?;
            crate::calc::apply_sidereal_default(xreturn, daya, has_speed);
        }

        Ok(())
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
