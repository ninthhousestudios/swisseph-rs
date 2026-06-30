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
    stars: crate::stars::StarCatalog,
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
        let stars = crate::stars::load_catalog(config.ephe_path.as_deref());
        Ok(Self {
            config,
            leap_seconds,
            planet_files,
            moon_files,
            jpl_file,
            stars,
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
        let idx = self
            .config
            .sidereal_mode
            .map(|m| m as i32 as usize)
            .unwrap_or(0);
        if crate::ayanamsa::FIXED_STAR_INDICES.contains(&idx) {
            let (daya, _) = self.fixstar_ayanamsa(jd_tt, flags)?;
            if !flags.contains(CalcFlags::NONUT) {
                let dpsi = crate::nutation::nutation(jd_tt, flags, &self.config.astro_models).dpsi;
                return Ok(daya + dpsi * crate::constants::RADTODEG);
            }
            return Ok(daya);
        }
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

        let (mut x0, x2000_0, _) = self.calc_inner(jd_tt - dt, body, inner_flags)?;
        let (mut x2, x2000_2, _) = self.calc_inner(jd_tt + dt, body, inner_flags)?;
        let (mut x1, x2000_1, flags_used) = self.calc_inner(jd_tt, body, inner_flags)?;

        // Sidereal projection must be applied to each of the three points
        // BEFORE the 3-point derivative, matching C's `use_speed3`, which calls
        // swecalc three times with SEFLG_SIDEREAL set and then differences the
        // already-projected positions (sweph.c:495-519). Project positions only
        // (no SPEED) — the speed is what the 3-point derivative produces. Applying
        // the projection after differencing would discard the 3-point speed and,
        // for the ECL_T0/SSY branches, read a zero velocity from x2000 (the inner
        // evals run without SPEED), collapsing the longitude speed to ~0.
        if flags.contains(CalcFlags::SIDEREAL) && body != Body::EclipticNutation {
            let pos_flags = flags_used & !CalcFlags::SPEED;
            self.apply_sidereal(&mut x0, &x2000_0, jd_tt - dt, pos_flags)?;
            self.apply_sidereal(&mut x2, &x2000_2, jd_tt + dt, pos_flags)?;
            self.apply_sidereal(&mut x1, &x2000_1, jd_tt, pos_flags)?;
        }

        crate::calc::denormalize_positions(&mut x0, &x1, &mut x2);
        crate::calc::calc_speed_3point(&mut x1, &x0, &x2, dt);

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
            let idx = self
                .config
                .sidereal_mode
                .map(|m| m as i32 as usize)
                .unwrap_or(0);
            let (daya_val, daya_sp) = if crate::ayanamsa::FIXED_STAR_INDICES.contains(&idx) {
                self.fixstar_ayanamsa(jd_tt, flags)?
            } else {
                let a =
                    crate::ayanamsa::get_ayanamsa_with_speed(&self.config, jd_tt, flags, models)?;
                (a[0], a[1])
            };
            crate::calc::apply_sidereal_default(xreturn, [daya_val, daya_sp], has_speed);
        }

        Ok(())
    }

    // -----------------------------------------------------------------------
    // Fixed-star API
    // -----------------------------------------------------------------------

    /// Compute apparent position of a fixed star at `jd_tt` (TT).
    ///
    /// Returns `(canonical_name, CalcResult)` where the name is
    /// `"traditional,bayer"` matching `swe_fixstar2` output.
    pub fn fixstar2(
        &self,
        star: &str,
        jd_tt: f64,
        flags: CalcFlags,
    ) -> Result<(String, CalcResult), Error> {
        // C's swe_fixstar2 returns the original input iflag unchanged (it passes
        // iflag by value to fixstar_calc_from_struct and ignores the return).
        let orig_flags = flags;
        let flags = crate::calc::plaus_iflag(flags, self.config.ephemeris_source);
        let resolved = if let Some(s) = crate::stars::builtin_star(star) {
            s
        } else {
            self.stars.search(star)?
        };
        let data = self.calc_fixstar(&resolved, jd_tt, flags)?;
        let name = format!("{},{}", resolved.name, resolved.bayer);
        Ok((
            name,
            CalcResult {
                data,
                flags_used: orig_flags,
            },
        ))
    }

    /// UT variant of `fixstar2`.
    pub fn fixstar2_ut(
        &self,
        star: &str,
        jd_ut: f64,
        flags: CalcFlags,
    ) -> Result<(String, CalcResult), Error> {
        let dt = crate::deltat::calc_deltat(jd_ut, &self.config);
        self.fixstar2(star, jd_ut + dt, flags)
    }

    /// Magnitude lookup for a star by name. Searches catalog only — builtin
    /// stars are not available via this function, matching C `swe_fixstar2_mag`.
    pub fn fixstar2_mag(&self, star: &str) -> Result<(String, f64), Error> {
        let resolved = self.stars.search(star)?;
        let name = format!("{},{}", resolved.name, resolved.bayer);
        Ok((name, resolved.mag))
    }

    /// Core fixed-star position pipeline (port of `fixstar_calc_from_struct`).
    /// Always computes speed internally; zeros it in the output if the caller
    /// did not request `SPEED`.
    fn calc_fixstar(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
    ) -> Result<[f64; 6], Error> {
        use crate::bias::{fk4_fk5, frame_bias, icrs2fk5};
        use crate::calc::{nutate, precess_speed};
        use crate::constants::*;
        use crate::corrections::{aberr_light, deflect_light};
        use crate::math::{
            cartesian_to_polar_with_speed, polar_to_cartesian_with_speed, rotate_x_sincos,
        };
        use crate::nutation::nutation;
        use crate::obliquity::obliquity;
        use crate::precession::precess;
        use crate::types::FrameTransform;
        use crate::types::PrecessionDirection;

        let models = &self.config.astro_models;
        // Force speed internally; honor caller's SPEED for output (step 18).
        let iflag = flags | CalcFlags::SPEED;

        // Step 1: Elapsed days since catalog epoch.
        let t = if star.epoch == 1950.0 {
            jd - B1950
        } else {
            jd - J2000
        };

        // Step 2: Initial polar+speed vector (radians / AU / day).
        let rdist = if star.parall == 0.0 {
            1e9
        } else {
            (1.0 / (star.parall * RADTODEG * 3600.0)) * PARSEC_TO_AUNIT
        };
        let mut x: [f64; 6] = [
            star.ra,
            star.de,
            rdist,
            star.ramot / 36525.0,
            star.demot / 36525.0,
            star.radvel / 36525.0,
        ];

        // Step 3: Polar → Cartesian with full space-motion.
        x = polar_to_cartesian_with_speed(x);

        // Step 4: FK4/FK5/ICRS frame corrections.
        if star.epoch == 1950.0 {
            fk4_fk5(&mut x, B1950);
            let mut pos3 = [x[0], x[1], x[2]];
            precess(
                &mut pos3,
                B1950,
                CalcFlags::empty(),
                models,
                PrecessionDirection::DateToJ2000,
            );
            x[0] = pos3[0];
            x[1] = pos3[1];
            x[2] = pos3[2];
            let mut vel3 = [x[3], x[4], x[5]];
            precess(
                &mut vel3,
                B1950,
                CalcFlags::empty(),
                models,
                PrecessionDirection::DateToJ2000,
            );
            x[3] = vel3[0];
            x[4] = vel3[1];
            x[5] = vel3[2];
        }
        if star.epoch != 0.0 {
            // FK5 → ICRS.
            icrs2fk5(&mut x, true, true);
            // ICRS → J2000 frame bias. C's swi_get_denum returns 403 for Moshier,
            // so this is applied unconditionally here (403 >= 403 is always true).
            frame_bias(
                &mut x,
                J2000,
                CalcFlags::SPEED,
                models,
                FrameTransform::GcrsToJ2000,
            );
        }

        // Steps 5–6: Earth heliocentric position for parallax / deflection /
        // aberration. Moshier returns heliocentric, matching C's xearth for MOSEPH.
        let eps_j2000 = obliquity(J2000, CalcFlags::empty(), models);
        let pp =
            crate::moshier::backend::compute_pipeline(jd, crate::types::Body::Sun, &eps_j2000)?;
        let xobs = pp.earth_helio; // heliocentric Earth = geocenter reference

        // Step 7: Proper motion + parallax (geocentric).
        for i in 0..3 {
            x[i] += t * x[i + 3]; // proper motion over elapsed days
            x[i] -= xobs[i]; // subtract observer (parallax)
            x[i + 3] -= xobs[i + 3]; // subtract observer velocity
        }

        // Step 8: Gravitational deflection (dt=0 for stars, matching C).
        if !iflag.contains(CalcFlags::TRUEPOS) && !iflag.contains(CalcFlags::NOGDEFL) {
            let mut planet_helio = [0.0f64; 6];
            for i in 0..3 {
                planet_helio[i] = x[i] + xobs[i]; // heliocentric star direction
                planet_helio[i + 3] = x[i + 3];
            }
            deflect_light(&mut x, &xobs, &planet_helio, true);
        }

        // Step 9: Annual aberration — swi_aberr_light_ex pattern.
        // C computes Earth state at both t and t-dt; speed = (pos_t - pos_t-dt) / FIXSTAR_DT.
        // This replaces (not adds to) x[3..6], matching C's swi_aberr_light_ex.
        if !iflag.contains(CalcFlags::TRUEPOS) && !iflag.contains(CalcFlags::NOABERR) {
            let orig = [x[0], x[1], x[2]];
            let orig_vel = [x[3], x[4], x[5]];
            let ev: [f64; 3] = [xobs[3], xobs[4], xobs[5]];
            aberr_light(&mut x, &ev, false);
            // Earth state at t-dt for speed via finite difference.
            let pp_dt = crate::moshier::backend::compute_pipeline(
                jd - FIXSTAR_DT,
                crate::types::Body::Sun,
                &eps_j2000,
            )?;
            let xobs_dt = pp_dt.earth_helio;
            let ev_dt: [f64; 3] = [xobs_dt[3], xobs_dt[4], xobs_dt[5]];
            let mut xx2: [f64; 6] = [
                orig[0] - FIXSTAR_DT * orig_vel[0],
                orig[1] - FIXSTAR_DT * orig_vel[1],
                orig[2] - FIXSTAR_DT * orig_vel[2],
                orig_vel[0],
                orig_vel[1],
                orig_vel[2],
            ];
            aberr_light(&mut xx2, &ev_dt, false);
            for i in 0..3 {
                x[i + 3] = (x[i] - xx2[i]) / FIXSTAR_DT;
            }
        }

        // Step 10: ICRS → J2000 frame bias.
        // C condition: !(iflag & SEFLG_ICRS) && (denum >= 403 || BARYCTR).
        // Moshier: swi_get_denum returns 403, so this is always applied.
        if !iflag.contains(CalcFlags::ICRS) {
            frame_bias(&mut x, jd, iflag, models, FrameTransform::GcrsToJ2000);
        }

        // Step 11: Save J2000 equatorial Cartesian for sidereal branch.
        let xxsv = x;

        // Step 12: Precession J2000 → equinox of date.
        let eps_date = if !iflag.contains(CalcFlags::J2000) {
            let mut pos3 = [x[0], x[1], x[2]];
            precess(
                &mut pos3,
                jd,
                iflag,
                models,
                PrecessionDirection::J2000ToDate,
            );
            x[0] = pos3[0];
            x[1] = pos3[1];
            x[2] = pos3[2];
            precess_speed(&mut x, jd, iflag, models, PrecessionDirection::J2000ToDate);
            obliquity(jd, iflag, models)
        } else {
            obliquity(J2000, iflag, models)
        };

        // Step 13: Nutation.
        let nut_val = nutation(jd, iflag, models);
        let nutv = Some(nutation(jd - NUT_SPEED_INTV, iflag, models));
        if !iflag.contains(CalcFlags::NONUT) {
            nutate(&mut x, &eps_date, &nut_val, nutv.as_ref(), true);
        }

        // Step 14: Equatorial → ecliptic (skip when EQUATORIAL requested).
        if !iflag.contains(CalcFlags::EQUATORIAL) {
            let pos3 = rotate_x_sincos([x[0], x[1], x[2]], eps_date.sin_eps, eps_date.cos_eps);
            x[0] = pos3[0];
            x[1] = pos3[1];
            x[2] = pos3[2];
            let vel3 = rotate_x_sincos([x[3], x[4], x[5]], eps_date.sin_eps, eps_date.cos_eps);
            x[3] = vel3[0];
            x[4] = vel3[1];
            x[5] = vel3[2];
            if !iflag.contains(CalcFlags::NONUT) {
                let snut = nut_val.deps.sin();
                let cnut = nut_val.deps.cos();
                let pos3 = rotate_x_sincos([x[0], x[1], x[2]], snut, cnut);
                x[0] = pos3[0];
                x[1] = pos3[1];
                x[2] = pos3[2];
                let vel3 = rotate_x_sincos([x[3], x[4], x[5]], snut, cnut);
                x[3] = vel3[0];
                x[4] = vel3[1];
                x[5] = vel3[2];
            }
        }

        // Step 15: Sidereal transform.
        if iflag.contains(CalcFlags::SIDEREAL) {
            let bits = self.config.sidereal_bits;
            if bits.contains(crate::flags::SiderealBits::ECL_T0) {
                let (xecl, xequ) =
                    crate::ayanamsa::trop_ra2sid_lon(&xxsv, &self.config, models, iflag);
                x = if iflag.contains(CalcFlags::EQUATORIAL) {
                    xequ
                } else {
                    xecl
                };
            } else if bits.contains(crate::flags::SiderealBits::SSY_PLANE) {
                let xecl =
                    crate::ayanamsa::trop_ra2sid_lon_sosy(&xxsv, &self.config, models, iflag);
                x = xecl;
            } else {
                // Default: subtract ayanamsa from ecliptic (or equatorial) longitude.
                x = cartesian_to_polar_with_speed(x);
                let idx = self
                    .config
                    .sidereal_mode
                    .map(|m| m as i32 as usize)
                    .unwrap_or(0);
                let (daya_val, daya_sp) = if crate::ayanamsa::FIXED_STAR_INDICES.contains(&idx) {
                    self.fixstar_ayanamsa(jd, iflag)?
                } else {
                    let a =
                        crate::ayanamsa::get_ayanamsa_with_speed(&self.config, jd, iflag, models)?;
                    (a[0], a[1])
                };
                x[0] -= daya_val * DEGTORAD;
                x[3] -= daya_sp * DEGTORAD;
                x = polar_to_cartesian_with_speed(x);
            }
        }

        // Step 16: Cartesian → polar.
        if !iflag.contains(CalcFlags::XYZ) {
            x = cartesian_to_polar_with_speed(x);
        }

        // Step 17: Radians → degrees (angles only, not distances).
        if !iflag.contains(CalcFlags::RADIANS) && !iflag.contains(CalcFlags::XYZ) {
            x[0] *= RADTODEG;
            x[1] *= RADTODEG;
            x[3] *= RADTODEG;
            x[4] *= RADTODEG;
        }

        // Step 18: Zero speeds if caller did not request them.
        if !flags.contains(CalcFlags::SPEED) {
            x[3] = 0.0;
            x[4] = 0.0;
            x[5] = 0.0;
        }

        Ok(x)
    }

    // -----------------------------------------------------------------------
    // Fixed-star ayanamsa (port of swi_get_ayanamsa_ex fixed-star branches)
    // -----------------------------------------------------------------------

    /// Ayanamsa value (degrees) for one of the 12 fixed-star modes at `jd_tt`.
    /// Mirrors C's early-return block in `swi_get_ayanamsa_ex` (sweph.c:3049–3142).
    /// Does NOT add nutation; caller adds it when appropriate.
    fn fixstar_ayanamsa_single(&self, jd_tt: f64, flags: CalcFlags) -> Result<f64, Error> {
        use crate::math::{armc_to_mc, normalize_degrees};

        // Flag construction mirrors C's swi_get_ayanamsa_ex entry (sweph.c:3007–3028).
        let ephmask = CalcFlags::JPLEPH | CalcFlags::SWIEPH | CalcFlags::MOSEPH;
        let epheflag = flags & ephmask;
        let iflag_base = epheflag | CalcFlags::NONUT;
        let iflag_galequ = iflag_base | CalcFlags::TRUEPOS;
        let mut iflag_true = iflag_base;
        if flags.contains(CalcFlags::TRUEPOS) {
            iflag_true |= CalcFlags::TRUEPOS;
        }
        if flags.contains(CalcFlags::NOABERR) {
            iflag_true |= CalcFlags::NOABERR;
        }
        if flags.contains(CalcFlags::NOGDEFL) {
            iflag_true |= CalcFlags::NOGDEFL;
        }

        let idx = self
            .config
            .sidereal_mode
            .map(|m| m as i32 as usize)
            .unwrap_or(0);

        let daya = match idx {
            17 => {
                let (_, r) = self.fixstar2(",SgrA*", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 240.0)
            }
            27 => {
                let (_, r) = self.fixstar2("Spica", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 180.0)
            }
            28 => {
                let (_, r) = self.fixstar2(",zePsc", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 359.833_333_333_3)
            }
            29 => {
                let (_, r) = self.fixstar2(",deCnc", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 106.0)
            }
            30 => {
                let (_, r) = self.fixstar2(",SgrA*", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 210.0 - 90.0 * 0.381_966_011_3)
            }
            31 => {
                let (_, r) = self.fixstar2(",GP1958", jd_tt, iflag_galequ)?;
                normalize_degrees(r.data[0] - 150.0)
            }
            32 => {
                let (_, r) = self.fixstar2(",GPol", jd_tt, iflag_galequ)?;
                normalize_degrees(r.data[0] - 150.0)
            }
            33 => {
                let (_, r) = self.fixstar2(",GPol", jd_tt, iflag_galequ)?;
                normalize_degrees(r.data[0] - 150.0 - 6.666_666_666_7)
            }
            35 => {
                let (_, r) = self.fixstar2(",laSco", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 240.0)
            }
            36 => {
                // GALCENT_MULA_WILHELM: SgrA* in equatorial, project RA → MC longitude.
                // obliquity uses iflag_base (= ephe|NONUT), matching C's `iflag` at that point.
                let (_, r) = self.fixstar2(",SgrA*", jd_tt, iflag_true | CalcFlags::EQUATORIAL)?;
                let ra = r.data[0];
                let eps_deg =
                    crate::obliquity::obliquity(jd_tt, iflag_base, &self.config.astro_models)
                        .eps
                        .to_degrees();
                normalize_degrees(armc_to_mc(ra, eps_deg) - 246.666_666_666_7)
            }
            39 => {
                let (_, r) = self.fixstar2(",deCnc", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 103.492_642_216_25)
            }
            40 => {
                let (_, r) = self.fixstar2(",SgrA*", jd_tt, iflag_true)?;
                normalize_degrees(r.data[0] - 270.0)
            }
            _ => unreachable!("fixstar_ayanamsa_single: non-fixed-star index {idx}"),
        };

        Ok(daya)
    }

    /// Returns `(ayanamsa_deg, speed_deg_per_day)` for fixed-star ayanamsa modes.
    /// Speed via two-point derivative (matches C's `swi_get_ayanamsa_with_speed` pattern).
    fn fixstar_ayanamsa(&self, jd_tt: f64, flags: CalcFlags) -> Result<(f64, f64), Error> {
        const TINTV: f64 = 0.001;
        let d0 = self.fixstar_ayanamsa_single(jd_tt, flags)?;
        let d2 = self.fixstar_ayanamsa_single(jd_tt - TINTV, flags)?;
        Ok((d0, (d0 - d2) / TINTV))
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
                        if let Ok(ndat) = trimmed.parse::<i32>()
                            && ndat > last_hardcoded
                            && !table.contains(&ndat)
                        {
                            table.push(ndat);
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
