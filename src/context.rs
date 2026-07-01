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
        let idx = crate::ayanamsa::sidereal_index(&self.config);
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

    /// Tropical houses at `tjd_ut` (UT), no flags. Port of `swe_houses` (swehouse.c:120-160).
    pub fn houses(
        &self,
        tjd_ut: f64,
        geolat: f64,
        geolon: f64,
        hsys: crate::types::HouseSystem,
    ) -> Result<crate::houses::HouseResult, Error> {
        self.houses_ex2(tjd_ut, CalcFlags::empty(), geolat, geolon, hsys)
    }

    /// Houses at `tjd_ut` (UT) with flags, no speeds requested by the caller (speeds are
    /// always computed — `HouseResult` carries them unconditionally). Port of `swe_houses_ex`
    /// (swehouse.c:200-236).
    pub fn houses_ex(
        &self,
        tjd_ut: f64,
        flags: CalcFlags,
        geolat: f64,
        geolon: f64,
        hsys: crate::types::HouseSystem,
    ) -> Result<crate::houses::HouseResult, Error> {
        self.houses_ex2(tjd_ut, flags, geolat, geolon, hsys)
    }

    /// Houses at `tjd_ut` (UT) with flags and speeds. Port of `swe_houses_ex2`
    /// (swehouse.c:238-270). Computes ARMC + true obliquity from `tjd_ut`, resolves the Sun's
    /// declination for Sunshine house systems, and dispatches to the traditional-sidereal or
    /// tropical driver. See docs/c-ref-houses.md §3, §6, §11.
    pub fn houses_ex2(
        &self,
        tjd_ut: f64,
        flags: CalcFlags,
        geolat: f64,
        geolon: f64,
        hsys: crate::types::HouseSystem,
    ) -> Result<crate::houses::HouseResult, Error> {
        use crate::constants::{DEGTORAD, RADTODEG};
        use crate::types::HouseSystem;

        // C's swe_houses_ex2 resolves deltaT tidal acceleration via swe_deltat_ex(tjd_ut, iflag,
        // NULL) where iflag is the caller's flags (never carrying an ephemeris-source bit for a
        // typical house call) -- this falls through to SE_TIDAL_DEFAULT regardless of the
        // actually-configured ephemeris backend (swephlib.c:2545-2568, ~2701). Our stateless
        // Ephemeris::calc_deltat normally resolves tid_acc from config.ephemeris_source, which
        // would silently pick up e.g. Moshier's DE404 here -- diverging from C by several
        // microdegrees at pre-1900 epochs. Force TIDAL_DEFAULT to match.
        let deltat_config = {
            let mut c = self.config.clone();
            c.tidal_acceleration = Some(crate::constants::TIDAL_DEFAULT);
            c
        };
        let tjde = tjd_ut + crate::deltat::calc_deltat(tjd_ut, &deltat_config);
        let models = &self.config.astro_models;

        // eps is always computed with iflag=0, regardless of the caller's flags
        // (swehouse.c:245, c-ref-houses.md §3).
        let eps_mean = crate::obliquity::obliquity(tjde, CalcFlags::empty(), models).eps * RADTODEG;
        let nut = crate::nutation::nutation(tjde, CalcFlags::empty(), models);
        let mut dpsi_deg = nut.dpsi * RADTODEG;
        let mut deps_deg = nut.deps * RADTODEG;
        if flags.contains(CalcFlags::NONUT) {
            dpsi_deg = 0.0;
            deps_deg = 0.0;
        }
        let eps_true = eps_mean + deps_deg;
        let armc = crate::math::normalize_degrees(
            crate::sidereal_time::sidereal_time0(tjd_ut, eps_true, dpsi_deg, &deltat_config) * 15.0
                + geolon,
        );

        let sundec = if matches!(hsys, HouseSystem::Sunshine | HouseSystem::SunshineAlt) {
            let xp = self.calc_ut(tjd_ut, Body::Sun, CalcFlags::SPEED | CalcFlags::EQUATORIAL)?;
            Some(xp.data[1])
        } else {
            None
        };

        let mut result = if flags.contains(CalcFlags::SIDEREAL) {
            if self
                .config
                .sidereal_bits
                .intersects(SiderealBits::ECL_T0 | SiderealBits::SSY_PLANE)
            {
                return Err(Error::CError(
                    "sidereal house mode ECL_T0/SSY_PLANE not yet implemented".into(),
                ));
            }
            let ayanamsa = self.get_ayanamsa_ex(tjde, flags)?;
            crate::houses::sidereal_houses_trad(armc, geolat, eps_true, hsys, sundec, ayanamsa)?
        } else {
            crate::houses::houses_armc(armc, geolat, eps_true, hsys, sundec)?
        };

        if flags.contains(CalcFlags::RADIANS) {
            let ito = if hsys == HouseSystem::Gauquelin {
                36
            } else {
                12
            };
            for cusp in result.cusps.iter_mut().take(ito + 1).skip(1) {
                *cusp *= DEGTORAD;
            }
            let ascmc = &mut result.ascmc;
            ascmc.ascendant *= DEGTORAD;
            ascmc.mc *= DEGTORAD;
            ascmc.armc *= DEGTORAD;
            ascmc.vertex *= DEGTORAD;
            ascmc.equatorial_ascendant *= DEGTORAD;
            ascmc.coascendant_koch *= DEGTORAD;
            ascmc.coascendant_munkasey *= DEGTORAD;
            ascmc.polar_ascendant *= DEGTORAD;
        }

        Ok(result)
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
            let idx = crate::ayanamsa::sidereal_index(&self.config);
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

    /// Dispatcher: routes fixed-star computation to the correct backend.
    fn calc_fixstar(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
    ) -> Result<[f64; 6], Error> {
        match self.config.ephemeris_source {
            crate::types::EphemerisSource::Swiss => self.calc_fixstar_sweph(star, jd, flags),
            crate::types::EphemerisSource::Jpl => self.calc_fixstar_jpl(star, jd, flags),
            crate::types::EphemerisSource::Moshier => self.calc_fixstar_moshier(star, jd, flags),
        }
    }

    /// Moshier backend: computes heliocentric Earth via Moshier pipeline.
    fn calc_fixstar_moshier(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
    ) -> Result<[f64; 6], Error> {
        use crate::constants::{FIXSTAR_DT, J2000};
        use crate::obliquity::obliquity;

        let models = &self.config.astro_models;
        // Moshier returns heliocentric Earth, matching C's xearth for MOSEPH.
        let eps_j2000 = obliquity(J2000, CalcFlags::empty(), models);
        let pp =
            crate::moshier::backend::compute_pipeline(jd, crate::types::Body::Sun, &eps_j2000)?;
        let xobs = pp.earth_helio;
        let pp_dt = crate::moshier::backend::compute_pipeline(
            jd - FIXSTAR_DT,
            crate::types::Body::Sun,
            &eps_j2000,
        )?;
        let xobs_dt = pp_dt.earth_helio;
        // Moshier is heliocentric; Sun is at the origin, so sun_bary = 0.
        let sun_bary = [0.0f64; 6];
        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary)
    }

    /// SWIEPH backend: barycentric Earth for parallax/aberration, sun_bary for deflection.
    fn calc_fixstar_sweph(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
    ) -> Result<[f64; 6], Error> {
        use crate::calc::{find_file_or_nearest, sweph_positions};
        use crate::constants::FIXSTAR_DT;
        use crate::sweph_file::types::{SEI_EMB, SEI_MOON};

        let planet_file = find_file_or_nearest(&self.planet_files, SEI_EMB, jd).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd,
                start: 0.0,
                end: 0.0,
            },
        )?;
        let moon_file = find_file_or_nearest(&self.moon_files, SEI_MOON, jd).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd,
                start: 0.0,
                end: 0.0,
            },
        )?;
        let pp = sweph_positions(planet_file, moon_file, SEI_EMB, jd, true)?;
        // C uses barycentric Earth (xearth) for parallax and aberration.
        let xobs = pp.earth_bary;
        let sun_bary = pp.sun_bary;

        let jd_dt = jd - FIXSTAR_DT;
        let planet_file_dt = find_file_or_nearest(&self.planet_files, SEI_EMB, jd_dt).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd_dt,
                start: 0.0,
                end: 0.0,
            },
        )?;
        let moon_file_dt = find_file_or_nearest(&self.moon_files, SEI_MOON, jd_dt).ok_or(
            Error::BeyondEphemerisLimits {
                jd_tt: jd_dt,
                start: 0.0,
                end: 0.0,
            },
        )?;
        let pp_dt = sweph_positions(planet_file_dt, moon_file_dt, SEI_EMB, jd_dt, true)?;
        let xobs_dt = pp_dt.earth_bary;

        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary)
    }

    /// JPL backend: barycentric Earth for parallax/aberration, sun_bary for deflection.
    fn calc_fixstar_jpl(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
    ) -> Result<[f64; 6], Error> {
        use crate::constants::FIXSTAR_DT;
        use crate::jpl::{J_EARTH, J_SBARY, J_SUN, jpl_pleph};

        let file = self.jpl_file.as_ref().ok_or(Error::EphemerisNotAvailable {
            body: crate::types::Body::Sun,
            source: crate::types::EphemerisSource::Jpl,
        })?;

        // C uses barycentric Earth for parallax/aberration; deflection uses earth_helio
        // computed inside swi_deflect_light as earth_bary - sun_bary.
        let xobs = jpl_pleph(file, jd, J_EARTH, J_SBARY, true)?;
        let sun_bary = jpl_pleph(file, jd, J_SUN, J_SBARY, true)?;
        let xobs_dt = jpl_pleph(file, jd - FIXSTAR_DT, J_EARTH, J_SBARY, true)?;

        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary)
    }

    /// Core fixed-star position pipeline (port of `fixstar_calc_from_struct`).
    /// Steps 1–4 (catalog→Cartesian) and 7–18 (corrections→output).
    ///
    /// `xobs`/`xobs_dt`: Earth position for parallax (step 7) and aberration velocity (step 9).
    ///   Moshier: heliocentric Earth. SWIEPH/JPL: barycentric Earth (matching C's xearth).
    /// `sun_bary`: Sun's barycentric position at `jd` (all 6 components including velocity).
    ///   Moshier passes zero (Sun at origin). SWIEPH/JPL pass the actual barycentric Sun.
    ///   Used to compute earth_helio = xobs - sun_bary for step 8 (deflection), replicating
    ///   C's swi_deflect_light which internally computes e = earth_bary - sun_bary.
    fn calc_fixstar_inner(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
        xobs: [f64; 6],
        xobs_dt: [f64; 6],
        sun_bary: [f64; 6],
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
            // ICRS → J2000 frame bias. Applied unconditionally (denum >= 403 for all
            // modern backends: Moshier returns 403, SWIEPH/JPL use DE430+).
            frame_bias(
                &mut x,
                J2000,
                CalcFlags::SPEED,
                models,
                FrameTransform::GcrsToJ2000,
            );
        }

        // Step 7: Proper motion + parallax.
        // xobs is the Earth's position in the backend's native frame
        // (heliocentric for Moshier, barycentric for SWIEPH/JPL).
        for i in 0..3 {
            x[i] += t * x[i + 3]; // proper motion over elapsed days
            x[i] -= xobs[i]; // subtract observer (parallax)
            x[i + 3] -= xobs[i + 3]; // subtract observer velocity
        }

        // Step 8: Gravitational deflection (dt=0 for stars, matching C).
        // C's swi_deflect_light internally computes e = earth_bary - sun_bary = earth_helio,
        // regardless of what frame xobs is in. Replicate that: use earth_helio for deflection.
        // For Moshier sun_bary=[0;6] so earth_helio = xobs (already heliocentric).
        if !iflag.contains(CalcFlags::TRUEPOS) && !iflag.contains(CalcFlags::NOGDEFL) {
            let mut earth_helio_defl = [0.0f64; 6];
            for i in 0..6 {
                earth_helio_defl[i] = xobs[i] - sun_bary[i];
            }
            let mut planet_ref = [0.0f64; 6];
            for i in 0..3 {
                // x = star - xobs (geocentric); x + earth_helio = star - sun_bary = heliocentric star
                planet_ref[i] = x[i] + earth_helio_defl[i];
                planet_ref[i + 3] = x[i + 3];
            }
            deflect_light(&mut x, &earth_helio_defl, &planet_ref, true);
        }

        // Step 9: Annual aberration — swi_aberr_light_ex pattern.
        // C computes Earth state at both t and t-dt; speed = (pos_t - pos_t-dt) / FIXSTAR_DT.
        // This replaces (not adds to) x[3..6], matching C's swi_aberr_light_ex.
        if !iflag.contains(CalcFlags::TRUEPOS) && !iflag.contains(CalcFlags::NOABERR) {
            let orig = [x[0], x[1], x[2]];
            let orig_vel = [x[3], x[4], x[5]];
            let ev: [f64; 3] = [xobs[3], xobs[4], xobs[5]];
            aberr_light(&mut x, &ev, false);
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
        // Applied unconditionally for all modern backends (denum >= 403 always true).
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

        // Step 13: Nutation. Only computed when NONUT is unset; the value is
        // reused by step 14 (nutation-in-ecliptic rotation).
        let nut_val = if !iflag.contains(CalcFlags::NONUT) {
            let nv = nutation(jd, iflag, models);
            let nutv = nutation(jd - NUT_SPEED_INTV, iflag, models);
            nutate(&mut x, &eps_date, &nv, Some(&nutv), true);
            Some(nv)
        } else {
            None
        };

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
            if let Some(ref nv) = nut_val {
                let snut = nv.deps.sin();
                let cnut = nv.deps.cos();
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
                let idx = crate::ayanamsa::sidereal_index(&self.config);
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

        let idx = crate::ayanamsa::sidereal_index(&self.config);

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
        // Both samples are independently normalized to [0,360); use the signed
        // shortest difference so a 360° wrap between samples doesn't blow up the
        // speed (~3.6e5 deg/day spike). diff_degrees returns a value in (-180,180].
        Ok((d0, crate::math::diff_degrees(d0, d2) / TINTV))
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
