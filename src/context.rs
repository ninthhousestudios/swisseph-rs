// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

//! Primary API — [`Ephemeris`] construction, configuration dispatch, and all
//! public calculation entry points ([`calc`](Ephemeris::calc),
//! [`houses`](Ephemeris::houses), eclipses, rise/set, heliacal, etc.).

use std::fs;

#[cfg(feature = "jpl")]
use crate::calc::JplProvider;
#[cfg(any(feature = "swisseph-files", feature = "jpl"))]
use crate::calc::PositionProvider;
#[cfg(any(feature = "swisseph-files", feature = "jpl"))]
use crate::calc::SwephPositions;
#[cfg(feature = "swisseph-files")]
use crate::calc::SwephProvider;
use crate::config::EphemerisConfig;
use crate::date::LEAP_SECONDS;
use crate::error::Error;
use crate::flags::{CalcFlags, EclipseFlags, SiderealBits};
use crate::types::{Body, DeltaT, EphemerisSource, JdUt1};

/// Represents barycentric state of body, sun and earth
type BarycentricState = ([f64; 6], [f64; 6], [f64; 6]);

#[cfg(not(feature = "swisseph-files"))]
#[allow(dead_code)]
pub(crate) struct AsteroidMetaStub {
    pub h: f64,
    pub g: f64,
    pub diameter_km: f64,
    pub name: String,
}

/// Three raw geocentric moon samples for the osculating node/apogee, plus the
/// `istart`, backend-specific central-difference interval, and backend used.
type OscMoonSamples = ([[f64; 6]; 3], usize, f64, EphemerisSource);

/// Selects `ipli`'s position from a `SwephPositions` bundle in the frame
/// `swe_nod_aps`'s osculating branch needs: barycentric (`SE_NODBIT_OSCU_BAR`)
/// or heliocentric, with `Body::Earth` reading the always-populated
/// `earth_bary`/`earth_helio` fields rather than `planet_bary` (which, for the
/// `query = Body::Sun` dummy call `nodaps_osc_body_j2000` makes for Earth, is
/// the Sun's own position, not Earth's — see `docs/c-ref-nodaps.md` §A.4.1).
#[cfg(any(feature = "swisseph-files", feature = "jpl"))]
fn nodaps_osc_frame(pos: &SwephPositions, ipli: Body, want_bary: bool) -> [f64; 6] {
    if ipli == Body::Earth {
        if want_bary {
            pos.earth_bary
        } else {
            pos.earth_helio
        }
    } else if want_bary {
        pos.planet_bary
    } else {
        let mut helio = [0.0; 6];
        for (h, (p, s)) in helio
            .iter_mut()
            .zip(pos.planet_bary.iter().zip(pos.sun_bary.iter()))
        {
            *h = p - s;
        }
        helio
    }
}

/// The main entry point for all Swiss Ephemeris calculations.
///
/// Holds read-only configuration and any opened ephemeris data files. All methods take
/// `&self` — there is no mutable state or internal caching. Thread-safe by construction.
pub struct Ephemeris {
    config: EphemerisConfig,
    /// The caller's original `tidal_acceleration` before `Ephemeris::new`
    /// resolved it from the open file's DE number. Needed by `effective_config`
    /// to re-derive the correct tid_acc when per-call flags select a different
    /// ephemeris source than the configured one (e.g. MOSEPH on a Swiss config).
    user_tidal_acceleration: Option<f64>,
    leap_seconds: Vec<i32>,
    #[cfg(feature = "swisseph-files")]
    planet_files: Vec<crate::sweph_file::SwissEphFile>,
    #[cfg(feature = "swisseph-files")]
    moon_files: Vec<crate::sweph_file::SwissEphFile>,
    #[cfg(feature = "swisseph-files")]
    main_asteroid_files: Vec<crate::sweph_file::SwissEphFile>,
    #[cfg(feature = "swisseph-files")]
    asteroid_files: Vec<crate::sweph_file::SwissEphFile>,
    #[cfg(feature = "swisseph-files")]
    planet_moon_files: Vec<crate::sweph_file::SwissEphFile>,
    #[cfg(feature = "jpl")]
    jpl_file: Option<crate::jpl::JplFile>,
    stars: crate::stars::StarCatalog,
    fictitious_catalog: crate::fictitious::FictitiousCatalog,
}

impl Ephemeris {
    /// Construct an `Ephemeris` from `config`, opening any configured ephemeris data files
    /// (planet/moon/asteroid `.se1` or JPL) and loading the fixed-star and fictitious-planet
    /// catalogs. Resolves `tidal_acceleration` from the opened file's DE number when not
    /// explicitly set.
    pub fn new(mut config: EphemerisConfig) -> crate::Result<Self> {
        let user_tidal_acceleration = config.tidal_acceleration;
        let leap_seconds = Self::build_leap_seconds(&config)?;

        #[cfg(feature = "jpl")]
        let mut jpl_file: Option<crate::jpl::JplFile> = None;

        #[cfg(feature = "swisseph-files")]
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
            #[cfg(feature = "jpl")]
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
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => return Err(Error::UnsupportedEphemeris(EphemerisSource::Jpl)),
            EphemerisSource::Moshier => (Vec::new(), Vec::new()),
        };

        #[cfg(all(not(feature = "swisseph-files"), feature = "jpl"))]
        match config.ephemeris_source {
            EphemerisSource::Swiss => {
                return Err(Error::UnsupportedEphemeris(EphemerisSource::Swiss));
            }
            EphemerisSource::Jpl => {
                let dir = config
                    .ephe_path
                    .as_ref()
                    .ok_or_else(|| Error::FileFormat("ephe_path required for Jpl".to_string()))?;
                let filename = config.jpl_filename.as_deref().unwrap_or("de441.eph");
                let path = dir.join(filename);
                jpl_file = Some(crate::jpl::JplFile::open(&path)?);
            }
            EphemerisSource::Moshier => {}
        }

        #[cfg(not(any(feature = "swisseph-files", feature = "jpl")))]
        match config.ephemeris_source {
            EphemerisSource::Swiss => {
                return Err(Error::UnsupportedEphemeris(EphemerisSource::Swiss));
            }
            EphemerisSource::Jpl => {
                return Err(Error::UnsupportedEphemeris(EphemerisSource::Jpl));
            }
            EphemerisSource::Moshier => {}
        }

        // Asteroid files load when ephe_path is set, regardless of ephemeris source
        // (C reads asteroids from .se1 even under MOSEPH/JPLEPH — c-ref-asteroid.md §1.4).
        #[cfg(feature = "swisseph-files")]
        let (main_asteroid_files, asteroid_files) = if let Some(dir) = config.ephe_path.as_ref() {
            let main_ast = match crate::sweph_file::open_ephemeris_files(dir, "seas") {
                Ok(files) => files,
                Err(Error::FileNotFound(_)) => Vec::new(),
                Err(e) => return Err(e),
            };

            if !config.asteroid_numbers.is_empty() {
                let mut nums = config.asteroid_numbers.clone();
                nums.sort_unstable();
                nums.dedup();
                let mut ast = Vec::with_capacity(nums.len());
                for &n in &nums {
                    ast.push(crate::sweph_file::open_asteroid_file(dir, n)?);
                }
                (main_ast, ast)
            } else {
                (main_ast, Vec::new())
            }
        } else if !config.asteroid_numbers.is_empty() {
            return Err(Error::FileFormat(
                "ephe_path required for asteroid files".to_string(),
            ));
        } else {
            (Vec::new(), Vec::new())
        };

        // Planet-moon files load when ephe_path is set, regardless of ephemeris source
        // (the moon/COB offset always comes from the sepm file regardless of epheflag).
        #[cfg(feature = "swisseph-files")]
        let planet_moon_files = if let Some(dir) = config.ephe_path.as_ref() {
            if !config.planet_moon_numbers.is_empty() {
                let mut nums = config.planet_moon_numbers.clone();
                nums.sort_unstable();
                nums.dedup();
                let mut files = Vec::with_capacity(nums.len());
                for &n in &nums {
                    if !(9000..=9999).contains(&n) {
                        return Err(Error::InvalidBody(n));
                    }
                    files.push(crate::sweph_file::open_planet_moon_file(dir, n)?);
                }
                files
            } else {
                Vec::new()
            }
        } else if !config.planet_moon_numbers.is_empty() {
            return Err(Error::FileFormat(
                "ephe_path required for planetary moon files".to_string(),
            ));
        } else {
            Vec::new()
        };

        // Resolve the ephemeris-specific tidal acceleration from the open file's
        // DE number, mirroring C's `swi_get_tid_acc` (swephlib.c:3211–3221): JPL
        // uses the JPL file's denum, SWIEPH the moon (SEI_FILE_MOON) file's. This
        // is what makes ΔT — and therefore the topocentric observer offset — match
        // C away from J2000 (DE441 tid_acc differs from the DE431 default). Only
        // fill in when the caller hasn't pinned tid_acc explicitly (C's
        // `is_tid_acc_manual` short-circuit).
        if config.tidal_acceleration.is_none() {
            let denum = match config.ephemeris_source {
                #[cfg(feature = "swisseph-files")]
                EphemerisSource::Swiss => moon_files.first().map(|f| f.header().denum),
                #[cfg(not(feature = "swisseph-files"))]
                EphemerisSource::Swiss => None,
                #[cfg(feature = "jpl")]
                EphemerisSource::Jpl => jpl_file.as_ref().map(|f| f.header().denum),
                #[cfg(not(feature = "jpl"))]
                EphemerisSource::Jpl => None,
                EphemerisSource::Moshier => None,
            };
            if let Some(denum) = denum {
                config.tidal_acceleration = Some(crate::deltat::denum_to_tid_acc(denum));
            }
        }
        let stars = crate::stars::load_catalog(config.ephe_path.as_deref());
        let fictitious_catalog =
            crate::fictitious::load_fictitious_catalog(config.ephe_path.as_deref())?;
        Ok(Self {
            config,
            user_tidal_acceleration,
            leap_seconds,
            #[cfg(feature = "swisseph-files")]
            planet_files,
            #[cfg(feature = "swisseph-files")]
            moon_files,
            #[cfg(feature = "swisseph-files")]
            main_asteroid_files,
            #[cfg(feature = "swisseph-files")]
            asteroid_files,
            #[cfg(feature = "swisseph-files")]
            planet_moon_files,
            #[cfg(feature = "jpl")]
            jpl_file,
            stars,
            fictitious_catalog,
        })
    }

    /// Returns a reference to the configuration this `Ephemeris` was created with.
    pub fn config(&self) -> &EphemerisConfig {
        &self.config
    }

    /// Resolve the per-call effective config: if `flags` requests a different
    /// ephemeris source than `config.ephemeris_source`, clamp it to what this
    /// `Ephemeris` can actually serve (loaded backends), adjust
    /// `tidal_acceleration` to match, and return the modified config. When the
    /// effective source matches the config's, returns a zero-cost borrow.
    ///
    /// Capability: Swiss requires `!planet_files.is_empty()`, Jpl requires
    /// `jpl_file.is_some()`. Unavailable → C's fallback cascade
    /// (Jpl→Swiss→Moshier), signaled via the returned source.
    pub fn effective_config<'a>(
        &self,
        flags: CalcFlags,
        config: &'a EphemerisConfig,
    ) -> std::borrow::Cow<'a, EphemerisConfig> {
        let requested = crate::calc::requested_source(flags);
        let effective = match requested {
            Some(src) => self.clamp_source(src),
            None => config.ephemeris_source,
        };
        if effective == config.ephemeris_source {
            std::borrow::Cow::Borrowed(config)
        } else {
            let mut c = config.clone();
            c.ephemeris_source = effective;
            c.tidal_acceleration = self.user_tidal_acceleration;
            std::borrow::Cow::Owned(c)
        }
    }

    /// Clamp a requested source to what this `Ephemeris` actually loaded.
    /// Cascade: Jpl→Swiss→Moshier (C never hard-errors on missing files in
    /// the calc path — it falls back and signals via flags_used).
    fn clamp_source(&self, requested: EphemerisSource) -> EphemerisSource {
        match requested {
            EphemerisSource::Jpl => {
                #[cfg(feature = "jpl")]
                if self.jpl_file.is_some() {
                    return EphemerisSource::Jpl;
                }
                #[cfg(feature = "swisseph-files")]
                if !self.planet_files.is_empty() {
                    return EphemerisSource::Swiss;
                }
                EphemerisSource::Moshier
            }
            EphemerisSource::Swiss => {
                #[cfg(feature = "swisseph-files")]
                if !self.planet_files.is_empty() {
                    return EphemerisSource::Swiss;
                }
                EphemerisSource::Moshier
            }
            EphemerisSource::Moshier => EphemerisSource::Moshier,
        }
    }

    /// Returns the effective leap-second table (built-in + any extras from config).
    pub fn leap_seconds(&self) -> &[i32] {
        &self.leap_seconds
    }

    /// Look up the SE1 orbital-element metadata (H, G, diameter) for a numbered asteroid.
    /// Returns `None` for non-`Asteroid` bodies (main asteroids Chiron..Vesta deliberately
    /// return `None` — C never populates the globals from seas files for those; their
    /// magnitude/diameter data live in `MAG_ELEM` / `PLANETARY_DIAMETERS`).
    #[cfg(feature = "swisseph-files")]
    pub(crate) fn asteroid_meta(&self, body: Body) -> Option<&crate::sweph_file::AsteroidMeta> {
        let n = match body {
            Body::Asteroid(id) => id.mpc_number(),
            _ => return None,
        };
        let target_id = crate::constants::AST_OFFSET + n;
        self.asteroid_files
            .iter()
            .find(|f| f.planets().first().is_some_and(|p| p.body_id == target_id))
            .and_then(|f| f.header().asteroid.as_ref())
    }

    #[cfg(not(feature = "swisseph-files"))]
    pub(crate) fn asteroid_meta(&self, _body: Body) -> Option<&AsteroidMetaStub> {
        None
    }

    /// Compute planetary position at `jd_tt` (Julian Day, TT time scale).
    ///
    /// Returns [`CalcResult`] with `data[0..3]` = position (lon/lat/dist or RA/dec/dist or
    /// x/y/z depending on flags), `data[3..6]` = speed (degrees/day or AU/day) when `SPEED`
    /// is set. `flags_used` indicates which flags were actually applied (may differ from
    /// requested if the backend forced a fallback).
    ///
    /// Unlike the C library, this implementation does not cache computed positions.
    #[doc(alias = "swe_calc")]
    pub fn calc(&self, jd_tt: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        self.calc_with_config(jd_tt, body, flags, &self.config)
    }

    /// Same as [`calc`](Self::calc) but with an explicit config override.
    ///
    /// Callers pass a clone of [`config()`](Self::config) with per-call fields changed
    /// (e.g. `topographic` for a different observer position). This threads the override
    /// through the entire pipeline without requiring a new `Ephemeris` instance.
    pub fn calc_with_config(
        &self,
        jd_tt: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<CalcResult, Error> {
        let config = self.effective_config(flags, config);
        let flags = crate::calc::plaus_iflag(flags, config.ephemeris_source);
        if flags.contains(CalcFlags::TOPOCTR) && config.topographic.is_none() {
            return Err(Error::CError(
                "topocentric requires topographic position".to_string(),
            ));
        }

        if body == Body::Earth && !flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
            return Ok(CalcResult {
                data: [0.0; 6],
                flags_used: flags,
            });
        }

        if flags.contains(CalcFlags::SPEED3) {
            return self.calc_speed3(jd_tt, body, flags, &config);
        }

        let (mut xreturn, x2000, flags_used) = self.calc_inner(jd_tt, body, flags, &config)?;
        if flags.contains(CalcFlags::SIDEREAL) && body != Body::EclipticNutation {
            self.apply_sidereal(&mut xreturn, &x2000, jd_tt, flags_used, &config)?;
        }
        Ok(CalcResult {
            data: Self::extract_for_body(&xreturn, body, flags_used),
            flags_used,
        })
    }

    /// Compute planetary position at `jd_ut` (Julian Day, UT1 time scale).
    ///
    /// Converts UT → TT via Delta T, then delegates to [`calc`](Self::calc).
    #[doc(alias = "swe_calc_ut")]
    pub fn calc_ut(&self, jd_ut: f64, body: Body, flags: CalcFlags) -> Result<CalcResult, Error> {
        let config = self.effective_config(flags, &self.config);
        let dt = crate::deltat::calc_deltat(jd_ut, &config);
        self.calc(jd_ut + dt, body, flags)
    }

    /// Same as [`calc_ut`](Self::calc_ut) but with an explicit config override; see
    /// [`calc_with_config`](Self::calc_with_config).
    pub fn calc_ut_with_config(
        &self,
        jd_ut: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<CalcResult, Error> {
        let config = self.effective_config(flags, config);
        let dt = crate::deltat::calc_deltat(jd_ut, &config);
        self.calc_with_config(jd_ut + dt, body, flags, &config)
    }

    /// Ayanamsa (precession-corrected sidereal offset) at `jd_tt` (TT), degrees.
    /// Nutation is added unless `NONUT` is set in `flags`.
    #[doc(alias = "swe_get_ayanamsa_ex")]
    pub fn get_ayanamsa_ex(&self, jd_tt: f64, flags: CalcFlags) -> Result<f64, Error> {
        self.get_ayanamsa_ex_with_config(jd_tt, flags, &self.config)
    }

    /// Same as [`get_ayanamsa_ex`](Self::get_ayanamsa_ex) but with an explicit config override;
    /// see [`calc_with_config`](Self::calc_with_config).
    pub fn get_ayanamsa_ex_with_config(
        &self,
        jd_tt: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<f64, Error> {
        let idx = crate::ayanamsa::sidereal_index(config);
        if crate::ayanamsa::FIXED_STAR_INDICES.contains(&idx) {
            let (daya, _) = self.fixstar_ayanamsa(jd_tt, flags, config)?;
            if !flags.contains(CalcFlags::NONUT) {
                let dpsi = crate::nutation::nutation(jd_tt, flags, &config.astro_models).dpsi;
                return Ok(daya + dpsi * crate::constants::RADTODEG);
            }
            return Ok(daya);
        }
        crate::ayanamsa::get_ayanamsa_ex_nut(config, jd_tt, flags, &config.astro_models)
    }

    /// Ayanamsa at `jd_ut` (UT1), degrees. Nutation added unless `NONUT` is set.
    #[doc(alias = "swe_get_ayanamsa_ex_ut")]
    pub fn get_ayanamsa_ut(&self, jd_ut: f64, flags: CalcFlags) -> Result<f64, Error> {
        let config = self.effective_config(flags, &self.config);
        let dt = crate::deltat::calc_deltat(jd_ut, &config);
        self.get_ayanamsa_ex(jd_ut + dt, flags)
    }

    /// Legacy ayanamsa accessor (no nutation) at `jd_tt` (TT), degrees.
    #[doc(alias = "swe_get_ayanamsa")]
    pub fn get_ayanamsa(&self, jd_tt: f64) -> Result<f64, Error> {
        crate::ayanamsa::get_ayanamsa_ex(
            &self.config,
            jd_tt,
            CalcFlags::empty(),
            &self.config.astro_models,
        )
    }

    /// Tropical houses at `tjd_ut` (UT1), no flags.
    ///
    /// `geolat`/`geolon` in degrees (north/east positive). Returns cusps (degrees) and
    /// special points (Asc, MC, ARMC, Vertex, etc.).
    #[doc(alias = "swe_houses")]
    pub fn houses(
        &self,
        tjd_ut: f64,
        geolat: f64,
        geolon: f64,
        hsys: crate::types::HouseSystem,
    ) -> Result<crate::houses::HouseResult, Error> {
        self.houses_ex2(tjd_ut, CalcFlags::empty(), geolat, geolon, hsys)
    }

    /// Houses at `tjd_ut` (UT1) with flags. Speeds are always included in the result.
    #[doc(alias = "swe_houses_ex")]
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

    /// Houses at `tjd_ut` (UT1) with flags and speeds.
    ///
    /// Computes ARMC + true obliquity, resolves the Sun's declination for Sunshine systems,
    /// and dispatches to the sidereal or tropical driver. `geolat`/`geolon` in degrees
    /// (north/east positive). Supports `SIDEREAL`, `NONUT`, `RADIANS` flags.
    #[doc(alias = "swe_houses_ex2")]
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
            let bits = self.config.sidereal_bits;
            if bits.contains(SiderealBits::ECL_T0) {
                crate::houses::sidereal_houses_ecl_t0(
                    tjde,
                    armc,
                    eps_true,
                    [dpsi_deg, deps_deg],
                    geolat,
                    hsys,
                    sundec,
                    self.config.sidereal_t0,
                    self.config.sidereal_ayan_t0,
                    models,
                )?
            } else if bits.contains(SiderealBits::SSY_PLANE) {
                crate::houses::sidereal_houses_ssypl(
                    tjde,
                    armc,
                    eps_true,
                    [dpsi_deg, deps_deg],
                    geolat,
                    hsys,
                    sundec,
                    self.config.sidereal_t0,
                    self.config.sidereal_ayan_t0,
                    models,
                )?
            } else {
                let ayanamsa = self.get_ayanamsa_ex(tjde, flags)?;
                crate::houses::sidereal_houses_trad(armc, geolat, eps_true, hsys, sundec, ayanamsa)?
            }
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

    /// ARMC + true obliquity at `tjd_ut`, shared setup for [`Ephemeris::azalt`] /
    /// [`Ephemeris::azalt_rev`]. Port of `swe_azalt`/`swe_azalt_rev`'s ARMC construction
    /// (`swe_sidtime(tjd_ut)*15 + geolon`) and their `SE_ECL_NUT` true-obliquity lookup, which
    /// both resolve deltaT via `swe_deltat_ex(tjd_ut, -1, NULL)` -- the `-1` sentinel forces
    /// `SE_TIDAL_DEFAULT` regardless of the configured ephemeris backend, same pattern as
    /// [`Ephemeris::houses_ex2`]. See docs/c-ref-refraction-azalt.md §1 step 3, §8.
    pub(crate) fn azalt_armc_eps(&self, tjd_ut: f64, geolon: f64) -> (f64, f64) {
        let deltat_config = {
            let mut c = self.config.clone();
            c.tidal_acceleration = Some(crate::constants::TIDAL_DEFAULT);
            c
        };
        let tjde = tjd_ut + crate::deltat::calc_deltat(tjd_ut, &deltat_config);
        let models = &self.config.astro_models;
        let eps_true = crate::calc::calc_ecl_nut(tjde, CalcFlags::empty(), models)[0];

        let armc = crate::math::normalize_degrees(
            crate::sidereal_time::sidereal_time(tjd_ut, &deltat_config) * 15.0 + geolon,
        );
        (armc, eps_true)
    }

    /// Ecliptic/equatorial → azimuth + true/apparent altitude at `tjd_ut` (UT1).
    ///
    /// `geopos` = \[longitude (east+), latitude (north+), height (meters)\].
    /// `atpress` in hPa (0 = auto from height), `attemp` in °C, `lapse_rate` in K/m (C default
    /// 0.0065). `xin` = \[lon/RA, lat/dec\] in degrees.
    ///
    /// Returns \[azimuth (from south, clockwise via west), true altitude, apparent altitude\],
    /// degrees.
    #[doc(alias = "swe_azalt")]
    #[allow(clippy::too_many_arguments)]
    pub fn azalt(
        &self,
        tjd_ut: f64,
        dir: crate::azalt::AzAltDir,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
        lapse_rate: f64,
        xin: [f64; 2],
    ) -> [f64; 3] {
        let (armc, eps_true) = self.azalt_armc_eps(tjd_ut, geopos[0]);
        crate::azalt::azalt(
            dir, armc, eps_true, geopos, atpress, attemp, lapse_rate, xin,
        )
    }

    /// Azimuth + true altitude → ecliptic/equatorial coordinates at `tjd_ut` (UT1).
    ///
    /// Inverse of [`azalt`](Self::azalt)'s geometric transform (does NOT de-refract).
    /// `geopos` = \[longitude, latitude, height (unused)\]. `xin` = \[azimuth (from south,
    /// clockwise), true altitude\], degrees. Returns \[lon/RA, lat/dec\], degrees.
    #[doc(alias = "swe_azalt_rev")]
    pub fn azalt_rev(
        &self,
        tjd_ut: f64,
        dir: crate::azalt::HorDir,
        geopos: [f64; 3],
        xin: [f64; 2],
    ) -> [f64; 2] {
        let (armc, eps_true) = self.azalt_armc_eps(tjd_ut, geopos[0]);
        crate::azalt::azalt_rev(dir, armc, eps_true, geopos[1], xin)
    }

    /// Rise/set/meridian-transit search (full precision algorithm).
    ///
    /// `starname` selects a fixed star (ignoring `body`); `horhgt` is the local horizon
    /// height in degrees (-100 = auto-dip from `geopos[2]`). `geopos` = \[longitude (east+),
    /// latitude (north+), height (m)\]. `atpress` in hPa, `attemp` in °C.
    /// Returns the UT of the next event, or [`Error::CircumpolarBody`].
    #[doc(alias = "swe_rise_trans_true_hor")]
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_true_hor(
        &self,
        tjd_ut: f64,
        body: Body,
        starname: Option<&str>,
        epheflag: CalcFlags,
        rsmi: crate::flags::RiseSetFlags,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
        horhgt: f64,
    ) -> Result<crate::riseset::RiseSetResult, Error> {
        crate::riseset::rise_trans_true_hor(
            self, tjd_ut, body, starname, epheflag, rsmi, geopos, atpress, attemp, horhgt,
        )
    }

    /// Rise/set/meridian-transit search, dispatching to the fast algorithm when eligible.
    ///
    /// Fast path: not a fixed star, RISE/SET only, no FORCE_SLOW/twilight, body in
    /// Sun..TrueNode, |lat| <= 60 (65 for Sun). Otherwise falls back to
    /// [`rise_trans_true_hor`](Self::rise_trans_true_hor) with `horhgt = 0.0`.
    #[doc(alias = "swe_rise_trans")]
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        tjd_ut: f64,
        body: Body,
        starname: Option<&str>,
        epheflag: CalcFlags,
        rsmi: crate::flags::RiseSetFlags,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
    ) -> Result<crate::riseset::RiseSetResult, Error> {
        use crate::flags::RiseSetFlags;

        let is_fixstar = starname.is_some_and(|s| !s.is_empty());
        let no_twilight = !rsmi.intersects(
            RiseSetFlags::CIVIL_TWILIGHT
                | RiseSetFlags::NAUTIC_TWILIGHT
                | RiseSetFlags::ASTRO_TWILIGHT,
        );
        let classic_body =
            (Body::Sun.to_raw_id()..=Body::TrueNode.to_raw_id()).contains(&body.to_raw_id());
        let lat_ok = geopos[1].abs() <= 60.0 || (body == Body::Sun && geopos[1].abs() <= 65.0);

        let fast_eligible = !is_fixstar
            && rsmi.intersects(RiseSetFlags::RISE | RiseSetFlags::SET)
            && !rsmi.contains(RiseSetFlags::FORCE_SLOW)
            && no_twilight
            && classic_body
            && lat_ok;

        if fast_eligible {
            crate::riseset::rise_set_fast(
                self, tjd_ut, body, epheflag, rsmi, geopos, atpress, attemp,
            )
        } else {
            self.rise_trans_true_hor(
                tjd_ut, body, starname, epheflag, rsmi, geopos, atpress, attemp, 0.0,
            )
        }
    }

    /// Solar eclipse shadow geometry at `tjd_ut` (UT1): geographic position of greatest eclipse
    /// + core/penumbra shadow diameters, geocentric. Local circumstances come from
    ///   [`sol_eclipse_how`](Self::sol_eclipse_how).
    #[doc(alias = "swe_sol_eclipse_where")]
    pub fn sol_eclipse_where(
        &self,
        tjd_ut: f64,
        ifl: CalcFlags,
    ) -> Result<crate::eclipse::EclipseWhere, Error> {
        crate::eclipse::sol_eclipse_where(self, tjd_ut, ifl)
    }

    /// Raw local eclipse/occultation circumstances at a point. Unlike
    /// [`sol_eclipse_how`](Self::sol_eclipse_how) this is the bare internal `eclipse_how`:
    /// no CENTRAL/NONCENTRAL merge, no redundant az/alt recompute, no horizon-visibility gate.
    pub fn eclipse_how_at(
        &self,
        tjd_ut: f64,
        ipl: Body,
        starname: Option<&str>,
        ifl: CalcFlags,
        geopos: [f64; 3],
    ) -> Result<crate::eclipse::EclipseHow, Error> {
        crate::eclipse::eclipse_how(
            self, tjd_ut, ipl, starname, ifl, geopos[0], geopos[1], geopos[2],
        )
    }

    /// Local circumstances of a solar eclipse at a specific observer at `tjd_ut` (UT1):
    /// magnitude, obscuration, contact geometry, azimuth/altitude.
    /// `geopos` = \[longitude (east+), latitude (north+), height (m)\].
    #[doc(alias = "swe_sol_eclipse_how")]
    pub fn sol_eclipse_how(
        &self,
        tjd_ut: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
    ) -> Result<crate::eclipse::EclipseHow, Error> {
        crate::eclipse::sol_eclipse_how(self, tjd_ut, ifl, geopos)
    }

    /// Planetary phenomena at `tjd_et` (TT): phase angle (deg), illuminated fraction,
    /// elongation (deg), apparent diameter (deg), apparent magnitude, horizontal parallax (deg
    /// Moon only). Returns [`Phenomena`](crate::phenomena::Phenomena) + flags actually used.
    #[doc(alias = "swe_pheno")]
    pub fn pheno(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<(crate::phenomena::Phenomena, CalcFlags), Error> {
        crate::phenomena::pheno(self, tjd_et, body, flags)
    }

    /// [`pheno`](Self::pheno) with a per-call config override (topographic position,
    /// sidereal mode).
    pub fn pheno_with_config(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<(crate::phenomena::Phenomena, CalcFlags), Error> {
        crate::phenomena::pheno_with_config(self, tjd_et, body, flags, config)
    }

    /// UT-based [`pheno`](Self::pheno).
    #[doc(alias = "swe_pheno_ut")]
    pub fn pheno_ut(
        &self,
        tjd_ut: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<(crate::phenomena::Phenomena, CalcFlags), Error> {
        crate::phenomena::pheno_ut(self, tjd_ut, body, flags)
    }

    /// UT-based [`pheno_with_config`](Self::pheno_with_config).
    pub fn pheno_ut_with_config(
        &self,
        tjd_ut: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<(crate::phenomena::Phenomena, CalcFlags), Error> {
        crate::phenomena::pheno_ut_with_config(self, tjd_ut, body, flags, config)
    }

    /// Limiting visual magnitude of an object at `tjd_ut` (UT1).
    ///
    /// `dgeo` = \[longitude (east+), latitude (north+), altitude (m)\].
    /// `datm` = \[atmospheric pressure hPa, temperature °C, humidity 0��1, extinction coeff\].
    /// `dobs` = \[age, Snellen ratio, 0=naked-eye/1=binocular/2=telescope, aperture, magnification, 0\].
    #[doc(alias = "swe_vis_limit_mag")]
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag(
        &self,
        tjd_ut: f64,
        dgeo: &[f64; 3],
        datm: &mut [f64; 4],
        dobs: &mut [f64; 6],
        object_name: &str,
        epheflag: crate::flags::CalcFlags,
        helflag: crate::flags::HeliacalFlags,
    ) -> Result<crate::heliacal::VisLimitResult, Error> {
        crate::heliacal::vis_limit_mag(
            self,
            tjd_ut,
            dgeo,
            datm,
            dobs,
            object_name,
            epheflag,
            helflag,
        )
    }

    /// Topocentric arcus visionis at `tjd_ut` (UT1), degrees.
    ///
    /// All geometry is caller-supplied. Angles in degrees.
    #[doc(alias = "swe_topo_arcus_visionis")]
    #[allow(clippy::too_many_arguments)]
    pub fn topo_arcus_visionis(
        &self,
        tjd_ut: f64,
        dgeo: &[f64; 3],
        datm: &mut [f64; 4],
        dobs: &mut [f64; 6],
        helflag: crate::flags::HeliacalFlags,
        mag: f64,
        azi_obj: f64,
        alt_obj: f64,
        azi_sun: f64,
        azi_moon: f64,
        alt_moon: f64,
    ) -> Result<f64, Error> {
        crate::heliacal::topo_arcus_visionis(
            tjd_ut, dgeo, datm, dobs, helflag, mag, azi_obj, alt_obj, azi_sun, azi_moon, alt_moon,
        )
    }

    /// Heliacal angle (optimal altitude / arcus visionis) at `tjd_ut` (UT1).
    #[doc(alias = "swe_heliacal_angle")]
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_angle(
        &self,
        tjd_ut: f64,
        dgeo: &[f64; 3],
        datm: &mut [f64; 4],
        dobs: &mut [f64; 6],
        helflag: crate::flags::HeliacalFlags,
        mag: f64,
        azi_obj: f64,
        azi_sun: f64,
        azi_moon: f64,
        alt_moon: f64,
    ) -> Result<crate::heliacal::HeliacalAngleResult, Error> {
        crate::heliacal::heliacal_angle(
            tjd_ut, dgeo, datm, dobs, helflag, mag, azi_obj, azi_sun, azi_moon, alt_moon,
        )
    }

    /// Heliacal phenomena (visibility window, geometry, Yallop criteria) at `tjd_ut` (UT1).
    #[doc(alias = "swe_heliacal_pheno_ut")]
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_pheno_ut(
        &self,
        tjd_ut: f64,
        dgeo: &[f64; 3],
        datm: &mut [f64; 4],
        dobs: &mut [f64; 6],
        object_name: &str,
        event: crate::heliacal::HeliacalEventType,
        epheflag: crate::flags::CalcFlags,
        helflag: crate::flags::HeliacalFlags,
    ) -> Result<crate::heliacal::HeliacalPheno, Error> {
        crate::heliacal::heliacal_pheno_ut(
            self,
            tjd_ut,
            dgeo,
            datm,
            dobs,
            object_name,
            event,
            epheflag,
            helflag,
        )
    }

    /// Find the next heliacal event (morning first, evening last, evening first, morning last,
    /// acronychal rising/setting) for `object_name` after `tjd_start_ut` (UT1).
    #[doc(alias = "swe_heliacal_ut")]
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_ut(
        &self,
        tjd_start_ut: f64,
        dgeo: &[f64; 3],
        datm: &mut [f64; 4],
        dobs: &mut [f64; 6],
        object_name: &str,
        event: crate::heliacal::HeliacalEventType,
        epheflag: crate::flags::CalcFlags,
        helflag: crate::flags::HeliacalFlags,
    ) -> Result<crate::heliacal::HeliacalEvent, Error> {
        crate::heliacal::heliacal_ut(
            self,
            tjd_start_ut,
            dgeo,
            datm,
            dobs,
            object_name,
            event,
            epheflag,
            helflag,
        )
    }

    /// Planetary nodes and apsides of `body` at `tjd_et` (TT).
    ///
    /// Returns ascending/descending nodes and perihelion/aphelion as ecliptic positions (degrees).
    /// `method` selects mean or osculating elements.
    #[doc(alias = "swe_nod_aps")]
    pub fn nod_aps(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
        method: crate::nodaps::NodApsMethod,
    ) -> Result<crate::nodaps::NodesApsides, Error> {
        let eff = self.effective_config(flags, &self.config);
        if flags.contains(CalcFlags::TOPOCTR) && eff.topographic.is_none() {
            return Err(Error::CError(
                "topocentric requires topographic position".to_string(),
            ));
        }
        crate::nodaps::nod_aps(self, tjd_et, body, flags, method)
    }

    /// UT-based [`nod_aps`](Self::nod_aps): converts UT→TT via Delta T, then delegates.
    #[doc(alias = "swe_nod_aps_ut")]
    pub fn nod_aps_ut(
        &self,
        tjd_ut: f64,
        body: Body,
        flags: CalcFlags,
        method: crate::nodaps::NodApsMethod,
    ) -> Result<crate::nodaps::NodesApsides, Error> {
        let eff = self.effective_config(flags, &self.config);
        let tjde = tjd_ut + crate::deltat::calc_deltat(tjd_ut, &eff);
        self.nod_aps(tjde, body, flags, method)
    }

    /// Osculating (Keplerian) orbital elements of `body` at `tjd_et` (TT).
    ///
    /// Rejects Sun, lunar nodes, and apsides. Note: `TOPOCTR` flag is bit-aliased as
    /// `ORBEL_AA` here (sum masses inside the orbit), not a topocentric request.
    #[doc(alias = "swe_get_orbital_elements")]
    pub fn get_orbital_elements(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<crate::orbit::OrbitalElements, Error> {
        crate::orbit::get_orbital_elements(self, tjd_et, body, flags)
    }

    /// Maximum, minimum, and current true distance of `body` at `tjd_et` (TT), in AU.
    /// Returns `(dmax, dmin, dtrue)`.
    #[doc(alias = "swe_orbit_max_min_true_distance")]
    pub fn orbit_max_min_true_distance(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<(f64, f64, f64), Error> {
        crate::orbit::orbit_max_min_true_distance(self, tjd_et, body, flags)
    }

    // -----------------------------------------------------------------------
    // Crossings (swe_solcross / mooncross / mooncross_node / helio_cross)
    // -----------------------------------------------------------------------

    /// Next JD (TT) at which the Sun's ecliptic longitude equals `x2cross` (degrees).
    #[doc(alias = "swe_solcross")]
    pub fn solcross(&self, x2cross: f64, jd_et: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::solcross(self, x2cross, jd_et, flags)
    }

    /// UT-based [`solcross`](Self::solcross).
    #[doc(alias = "swe_solcross_ut")]
    pub fn solcross_ut(&self, x2cross: f64, jd_ut: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::solcross_ut(self, x2cross, jd_ut, flags)
    }

    /// Next JD (TT) at which the Moon's ecliptic longitude equals `x2cross` (degrees).
    #[doc(alias = "swe_mooncross")]
    pub fn mooncross(&self, x2cross: f64, jd_et: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::mooncross(self, x2cross, jd_et, flags)
    }

    /// UT-based [`mooncross`](Self::mooncross).
    #[doc(alias = "swe_mooncross_ut")]
    pub fn mooncross_ut(&self, x2cross: f64, jd_ut: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::mooncross_ut(self, x2cross, jd_ut, flags)
    }

    /// Next JD (TT) at which the Moon crosses its node (ecliptic latitude = 0).
    #[doc(alias = "swe_mooncross_node")]
    pub fn mooncross_node(
        &self,
        jd_et: f64,
        flags: CalcFlags,
    ) -> Result<crate::crossings::MoonCrossing, Error> {
        crate::crossings::mooncross_node(self, jd_et, flags)
    }

    /// UT-based [`mooncross_node`](Self::mooncross_node).
    #[doc(alias = "swe_mooncross_node_ut")]
    pub fn mooncross_node_ut(
        &self,
        jd_ut: f64,
        flags: CalcFlags,
    ) -> Result<crate::crossings::MoonCrossing, Error> {
        crate::crossings::mooncross_node_ut(self, jd_ut, flags)
    }

    /// Next JD (TT) at which `body`'s heliocentric longitude equals `x2cross` (degrees).
    /// `dir >= 0` searches forward, `dir < 0` searches backward.
    #[doc(alias = "swe_helio_cross")]
    pub fn helio_cross(
        &self,
        body: Body,
        x2cross: f64,
        jd_et: f64,
        flags: CalcFlags,
        dir: i32,
    ) -> Result<f64, Error> {
        crate::crossings::helio_cross(self, body, x2cross, jd_et, flags, dir)
    }

    /// UT-based [`helio_cross`](Self::helio_cross).
    #[doc(alias = "swe_helio_cross_ut")]
    pub fn helio_cross_ut(
        &self,
        body: Body,
        x2cross: f64,
        jd_ut: f64,
        flags: CalcFlags,
        dir: i32,
    ) -> Result<f64, Error> {
        crate::crossings::helio_cross_ut(self, body, x2cross, jd_ut, flags, dir)
    }

    /// Position of `body` as seen from `center` (planetocentric coordinates) at `jd_tt` (TT).
    /// Swiss/JPL only (Moshier returns `Err`).
    #[doc(alias = "swe_calc_pctr")]
    pub fn calc_pctr(
        &self,
        jd_tt: f64,
        body: Body,
        center: Body,
        flags: CalcFlags,
    ) -> Result<CalcResult, Error> {
        let body = crate::calc::normalize_asteroid_aliases(body);
        let center = crate::calc::normalize_asteroid_aliases(center);
        if body == center {
            return Err(Error::CError(
                "ipl and iplctr must not be identical".to_string(),
            ));
        }
        let config = self.effective_config(flags, &self.config);
        let flags = crate::calc::plaus_iflag(flags, config.ephemeris_source);

        // C's swe_calc_pctr internally calls swe_calc with SEFLG_BARYCTR.
        // Moshier doesn't support barycentric positions → propagate the error.
        if config.ephemeris_source == EphemerisSource::Moshier {
            return Err(Error::CError(
                "barycentric Moshier positions are not supported".to_string(),
            ));
        }

        // Strip HELCTR/BARYCTR from user flags (sweph.c:8059)
        let flags = flags & !(CalcFlags::HELCTR | CalcFlags::BARYCTR);
        let models = &config.astro_models;
        let has_speed = flags.contains(CalcFlags::SPEED);

        // §1: Prime obliquity/nutation at tjd + Δt(tjd) — a third, distinct epoch.
        // For J2000 output the §9 ecliptic rotation uses the J2000 mean obliquity
        // (oec2000), not obliquity-of-date, mirroring calc::precess_and_ephem's
        // J2000 branch. (nut_val is unused when J2000 forces NONUT, but priming it
        // is harmless and keeps the non-J2000 path unchanged.)
        let dt_prime = crate::deltat::calc_deltat(jd_tt, &config);
        let eps = if flags.contains(CalcFlags::J2000) {
            crate::obliquity::obliquity(crate::constants::J2000, flags, models)
        } else {
            crate::obliquity::obliquity(jd_tt + dt_prime, flags, models)
        };
        let nut_val = crate::nutation::nutation(jd_tt + dt_prime, flags, models);

        // §2: Barycentric J2000-equatorial states of both bodies at tjd
        let (xx0, _, _) = self.pctr_bary_state(jd_tt, body)?;
        let (xxctr, _, _) = self.pctr_bary_state(jd_tt, center)?;

        // §3: Light-time iteration + re-eval at retarded time
        let (t, xxsp, xx, xxctr2, eb_defl, sb_defl) = if flags.contains(CalcFlags::TRUEPOS) {
            // No light-time; deflection/aberration also gated on !TRUEPOS, so
            // earth_bary/sun_bary values are unused — zero placeholders.
            (jd_tt, [0.0; 3], xx0, xxctr, [0.0; 6], [0.0; 6])
        } else {
            let (t, _dt, xxsp) = crate::calc::pctr_light_time(jd_tt, &xx0, &xxctr, has_speed);

            // §3d: Re-evaluate both bodies at retarded time
            let (xx_t, eb_t, sb_t) = self.pctr_bary_state(t, body)?;
            let (xxctr2_t, _, _) = self.pctr_bary_state(t, center)?;

            (t, xxsp, xx_t, xxctr2_t, eb_t, sb_t)
        };

        // §4–§9: Planetocentric subtraction, deflection, aberration, bias, precess, output
        // Note: §4 subtracts xxctr (center at tjd, NOT xxctr2 at t) — see c-ref-pctr §4.
        let nut_epoch = jd_tt + dt_prime;
        let (mut xreturn, x2000) = crate::calc::pctr_pipeline(
            &xx, &xxctr, &xxctr2, &xxsp, t, jd_tt, nut_epoch, &eb_defl, &sb_defl, flags, &eps,
            &nut_val, models,
        );

        // §9 sidereal tail
        if flags.contains(CalcFlags::SIDEREAL) {
            self.apply_sidereal(&mut xreturn, &x2000, jd_tt, flags, &config)?;
        }

        Ok(CalcResult {
            data: crate::calc::extract_output(&xreturn, flags),
            flags_used: flags,
        })
    }

    /// Barycentric J2000-equatorial state of `body` at epoch `t`, plus Earth and
    /// Sun barycentric states (always populated by the provider regardless of the
    /// queried body). Swiss/JPL only — Moshier is rejected before reaching here.
    #[cfg_attr(
        not(any(feature = "swisseph-files", feature = "jpl")),
        allow(unused_variables)
    )]
    fn pctr_bary_state(&self, t: f64, body: Body) -> Result<BarycentricState, Error> {
        let eps_j2000 = crate::obliquity::obliquity(
            crate::constants::J2000,
            CalcFlags::empty(),
            &self.config.astro_models,
        );

        match self.config.ephemeris_source {
            EphemerisSource::Moshier => Err(Error::CError(
                "barycentric Moshier positions are not supported".to_string(),
            )),
            #[cfg(feature = "swisseph-files")]
            EphemerisSource::Swiss => {
                let provider = SwephProvider {
                    planet_files: &self.planet_files,
                    moon_files: &self.moon_files,
                    main_asteroid_files: &self.main_asteroid_files,
                    asteroid_files: &self.asteroid_files,
                    planet_moon_files: &self.planet_moon_files,
                };
                self.pctr_bary_from_provider(&provider, t, body, &eps_j2000)
            }
            #[cfg(not(feature = "swisseph-files"))]
            EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            EphemerisSource::Jpl => {
                let provider = JplProvider {
                    file: self.jpl_file.as_ref().expect(
                        "JPL file must be initialized into Ephemeris upon calling Ephemeris::new()",
                    ),
                };
                self.pctr_bary_from_provider(&provider, t, body, &eps_j2000)
            }
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
        }
    }

    #[cfg(any(feature = "swisseph-files", feature = "jpl"))]
    fn pctr_bary_from_provider<P: crate::calc::PositionProvider>(
        &self,
        provider: &P,
        t: f64,
        body: Body,
        _eps_j2000: &crate::types::Epsilon,
    ) -> Result<BarycentricState, Error> {
        match body {
            Body::Moon => {
                let moon_geo = provider.moon_geo(t, true)?;
                let pos = provider.positions(Body::Sun, t, true)?;
                let mut body_bary = [0.0; 6];
                for i in 0..6 {
                    body_bary[i] = moon_geo[i] + pos.earth_bary[i];
                }
                Ok((body_bary, pos.earth_bary, pos.sun_bary))
            }
            Body::Earth => {
                let pos = provider.positions(Body::Sun, t, true)?;
                Ok((pos.earth_bary, pos.earth_bary, pos.sun_bary))
            }
            Body::Sun => {
                let pos = provider.positions(Body::Sun, t, true)?;
                Ok((pos.sun_bary, pos.earth_bary, pos.sun_bary))
            }
            _ => {
                let pos = provider.positions(body, t, true)?;
                Ok((pos.planet_bary, pos.earth_bary, pos.sun_bary))
            }
        }
    }

    /// Observer / origin geometry for the nodes-&-apsides pipeline at epoch `t`
    /// (TT), in equatorial-J2000 cartesian. Replaces C's `xsun`/`xear`/`xobs`
    /// globals (swecl.c A.5.1) with an explicit per-epoch computation. The
    /// final `xobs` (which of `sun_bary`/`xear`/`topo` to actually observe
    /// from, per HELCTR/BARYCTR/SE_SUN) is resolved by `nodaps::select_xobs`,
    /// not here.
    ///
    /// `xear` is Earth's position in the node's native frame: heliocentric
    /// (≈barycentric, Moshier has no true barycenter) for Moshier, real
    /// barycentric for Swiss/JPL (`swed.pldat[SEI_EARTH].x` is barycentric in
    /// C's own shared frame). `sun_bary` is `[0.0; 6]` for Moshier (matching
    /// `transform_nodaps_output`'s `is_moseph` gate, which never adds it).
    pub(crate) fn nodaps_observer(
        &self,
        t: f64,
        flags: CalcFlags,
    ) -> Result<crate::nodaps::ObsFrame, Error> {
        let _eff = self.effective_config(flags, &self.config);
        let config = &*_eff;
        let models = &config.astro_models;
        let topo = crate::calc::topo_offset(t, flags, config, models);
        match config.ephemeris_source {
            EphemerisSource::Moshier => {
                let eps_j2000 = crate::obliquity::obliquity(
                    crate::constants::J2000,
                    CalcFlags::empty(),
                    models,
                );
                let pp = crate::moshier::backend::compute_pipeline(t, Body::Sun, &eps_j2000)?;
                Ok(crate::nodaps::ObsFrame {
                    sun_bary: [0.0; 6],
                    xear: pp.earth_helio,
                    topo,
                })
            }
            #[cfg(feature = "swisseph-files")]
            EphemerisSource::Swiss => {
                let provider = SwephProvider {
                    planet_files: &self.planet_files,
                    moon_files: &self.moon_files,
                    main_asteroid_files: &self.main_asteroid_files,
                    asteroid_files: &self.asteroid_files,
                    planet_moon_files: &self.planet_moon_files,
                };
                let pos = provider.positions(Body::Sun, t, true)?;
                Ok(crate::nodaps::ObsFrame {
                    sun_bary: pos.sun_bary,
                    xear: pos.earth_bary,
                    topo,
                })
            }
            #[cfg(not(feature = "swisseph-files"))]
            EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            EphemerisSource::Jpl => {
                let provider = JplProvider {
                    file: self.jpl_file.as_ref().expect(
                        "JPL file must be initialized into Ephemeris upon calling Ephemeris::new()",
                    ),
                };
                let pos = provider.positions(Body::Sun, t, true)?;
                Ok(crate::nodaps::ObsFrame {
                    sun_bary: pos.sun_bary,
                    xear: pos.earth_bary,
                    topo,
                })
            }
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
        }
    }

    /// `SE_NODBIT_OSCU`/`SE_NODBIT_OSCU_BAR` (`docs/c-ref-nodaps.md` §A.4.1,
    /// §A.4.2) — heliocentric or barycentric J2000-equatorial cartesian
    /// position+speed of `ipli` at `t`, across all three backends, with no
    /// light-time/aberration/deflection (matching `swe_nod_aps`'s
    /// `SEFLG_TRUEPOS`-forced `iflJ2000`).
    ///
    /// The Moon is always geocentric — `want_bary` is ignored, matching C
    /// never setting HELCTR/BARYCTR for `ipli == SE_MOON`. `Body::Earth` gets
    /// the Earth->EMB correction added back in (swecl.c:5293-5296: `swe_calc`
    /// on `SE_EARTH` returns bare Earth with the Moon's contribution already
    /// subtracted out; `swe_nod_aps` wants the EMB `ipli` actually orbits
    /// around).
    ///
    /// `want_bary` (`SE_NODBIT_OSCU_BAR`) is only meaningful for Swiss/JPL,
    /// which carry a real solar-system barycenter; Moshier has no such frame
    /// and rejects it the same way `calc_inner` rejects a bare `SEFLG_BARYCTR`
    /// request (`Error::UnsupportedFlags`).
    pub(crate) fn nodaps_osc_body_j2000(
        &self,
        t: f64,
        ipli: Body,
        want_bary: bool,
        flags: CalcFlags,
    ) -> Result<[f64; 6], Error> {
        let _eff = self.effective_config(flags, &self.config);
        let source = _eff.ephemeris_source;
        let models = &_eff.astro_models;
        let eps_j2000 =
            crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);

        if ipli == Body::Moon {
            return self.raw_moon_at(source, t, &eps_j2000);
        }

        let mut xx = match source {
            EphemerisSource::Moshier => {
                if want_bary {
                    return Err(Error::UnsupportedFlags(CalcFlags::BARYCTR));
                }
                if ipli == Body::Earth {
                    crate::moshier::backend::compute_pipeline(t, Body::Sun, &eps_j2000)?.earth_helio
                } else {
                    crate::moshier::backend::compute_pipeline(t, ipli, &eps_j2000)?.planet_helio
                }
            }
            #[cfg(feature = "swisseph-files")]
            EphemerisSource::Swiss => {
                let provider = SwephProvider {
                    planet_files: &self.planet_files,
                    moon_files: &self.moon_files,
                    main_asteroid_files: &self.main_asteroid_files,
                    asteroid_files: &self.asteroid_files,
                    planet_moon_files: &self.planet_moon_files,
                };
                let query = if ipli == Body::Earth { Body::Sun } else { ipli };
                let pos = provider.positions(query, t, true)?;
                nodaps_osc_frame(&pos, ipli, want_bary)
            }
            #[cfg(not(feature = "swisseph-files"))]
            EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            EphemerisSource::Jpl => {
                let provider = JplProvider {
                    file: self.jpl_file.as_ref().expect(
                        "JPL file must be initialized into Ephemeris upon calling Ephemeris::new()",
                    ),
                };
                let query = if ipli == Body::Earth { Body::Sun } else { ipli };
                let pos = provider.positions(query, t, true)?;
                nodaps_osc_frame(&pos, ipli, want_bary)
            }
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
        };

        if ipli == Body::Earth {
            let moon = self.raw_moon_at(source, t, &eps_j2000)?;
            for i in 0..6 {
                xx[i] += moon[i] / (crate::constants::EARTH_MOON_MRAT + 1.0);
            }
        }

        Ok(xx)
    }

    /// Global solar eclipse search: next/previous eclipse anywhere on Earth from `tjd_start`
    /// (UT1). `ifltype` filters eclipse types (empty = all). `backward` searches past.
    #[doc(alias = "swe_sol_eclipse_when_glob")]
    pub fn sol_eclipse_when_glob(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        ifltype: EclipseFlags,
        backward: bool,
    ) -> Result<crate::eclipse::SolarEclipseGlobal, Error> {
        crate::eclipse::sol_eclipse_when_glob(self, tjd_start, ifl, ifltype, backward)
    }

    /// Local solar eclipse search: next/previous eclipse visible from `geopos`, with contact
    /// times and local circumstances. `geopos` = \[lon (east+), lat (north+), height (m)\].
    #[doc(alias = "swe_sol_eclipse_when_loc")]
    pub fn sol_eclipse_when_loc(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
        backward: bool,
    ) -> Result<crate::eclipse::SolarEclipseLocal, Error> {
        crate::eclipse::sol_eclipse_when_loc(self, tjd_start, ifl, geopos, backward)
    }

    /// Local circumstances of a lunar eclipse at `tjd_ut` (UT1): magnitude, Saros, Moon
    /// azimuth/altitude. `geopos` = \[lon (east+), lat (north+), height (m)\].
    #[doc(alias = "swe_lun_eclipse_how")]
    pub fn lun_eclipse_how(
        &self,
        tjd_ut: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
    ) -> Result<crate::eclipse::LunarEclipseHow, Error> {
        crate::eclipse::swe_lun_eclipse_how(self, tjd_ut, ifl, geopos)
    }

    /// Global lunar eclipse search from `tjd_start` (UT1). Geocentric — no observer position.
    /// `ifltype` filters types (empty = any).
    #[doc(alias = "swe_lun_eclipse_when")]
    pub fn lun_eclipse_when(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        ifltype: EclipseFlags,
        backward: bool,
    ) -> Result<crate::eclipse::LunarEclipseGlobal, Error> {
        crate::eclipse::lun_eclipse_when(self, tjd_start, ifl, ifltype, backward)
    }

    /// Local lunar eclipse search: visible from `geopos` (Moon above horizon), contact times
    /// clipped to moonrise/moonset. `geopos` = \[lon (east+), lat (north+), height (m)\].
    #[doc(alias = "swe_lun_eclipse_when_loc")]
    pub fn lun_eclipse_when_loc(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
        backward: bool,
    ) -> Result<crate::eclipse::LunarEclipseLocal, Error> {
        crate::eclipse::lun_eclipse_when_loc(self, tjd_start, ifl, geopos, backward)
    }

    /// Geographic position of maximal lunar occultation of `body`/`starname` at `tjd_ut` (UT1).
    /// `starname` (if non-empty) takes precedence over `body`.
    #[doc(alias = "swe_lun_occult_where")]
    pub fn lun_occult_where(
        &self,
        tjd_ut: f64,
        body: Body,
        starname: Option<&str>,
        ifl: CalcFlags,
    ) -> Result<crate::eclipse::EclipseWhere, Error> {
        crate::eclipse::lun_occult_where(self, tjd_ut, body, starname, ifl)
    }

    /// Global occultation search: next/previous occultation of `body`/`starname` by the Moon
    /// anywhere on Earth from `tjd_start` (UT1). `starname` (if non-empty) takes precedence.
    #[doc(alias = "swe_lun_occult_when_glob")]
    pub fn lun_occult_when_glob(
        &self,
        tjd_start: f64,
        body: Body,
        starname: Option<&str>,
        ifl: CalcFlags,
        ifltype: EclipseFlags,
        backward: bool,
    ) -> Result<crate::eclipse::OccultGlobal, Error> {
        crate::eclipse::lun_occult_when_glob(
            self, tjd_start, body, starname, ifl, ifltype, backward,
        )
    }

    /// Local occultation search: visible from `geopos`, with contact times and circumstances.
    /// `starname` (if non-empty) takes precedence over `body`.
    /// `geopos` = \[lon (east+), lat (north+), height (m)\].
    #[doc(alias = "swe_lun_occult_when_loc")]
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_loc(
        &self,
        tjd_start: f64,
        body: Body,
        starname: Option<&str>,
        ifl: CalcFlags,
        geopos: [f64; 3],
        backward: bool,
    ) -> Result<crate::eclipse::OccultLocal, Error> {
        crate::eclipse::lun_occult_when_loc(self, tjd_start, body, starname, ifl, geopos, backward)
    }

    /// Gauquelin sector position (geometric method, `imeth` 0 or 1) at `t_ut` (UT1).
    ///
    /// `imeth` 0 = with ecliptic latitude, 1 = without. Returns sector 1.0–36.0.
    /// `starname` (if non-empty) uses the fixed-star position instead of `body`.
    #[doc(alias = "swe_gauquelin_sector")]
    #[allow(clippy::too_many_arguments)]
    pub fn gauquelin_sector_geometric(
        &self,
        t_ut: f64,
        body: Body,
        starname: Option<&str>,
        imeth: i32,
        flags: CalcFlags,
        geolon: f64,
        geolat: f64,
    ) -> Result<f64, Error> {
        use crate::constants::RADTODEG;
        use crate::types::HouseSystem;

        if !(0..=1).contains(&imeth) {
            return Err(Error::CError(format!(
                "invalid imeth for geometric gauquelin: {imeth}"
            )));
        }

        let eff = self.effective_config(flags, &self.config);
        let models = &eff.astro_models;
        let t_et = t_ut + crate::deltat::calc_deltat(t_ut, &eff);
        let eps = crate::obliquity::obliquity(t_et, flags, models).eps * RADTODEG;
        let nut = crate::nutation::nutation(t_et, flags, models);
        let dpsi_deg = nut.dpsi * RADTODEG;
        let deps_deg = nut.deps * RADTODEG;
        let eps_true = eps + deps_deg;
        let armc = crate::math::normalize_degrees(
            crate::sidereal_time::sidereal_time0(t_ut, eps_true, dpsi_deg, &eff) * 15.0 + geolon,
        );

        let x0 = if starname.is_some_and(|s| !s.is_empty()) {
            self.fixstar2(starname.unwrap(), t_et, flags)?.1
        } else {
            self.calc(t_et, body, flags)?
        };
        let lat = if imeth == 1 { 0.0 } else { x0.data[1] };

        crate::houses::house_pos(
            armc,
            geolat,
            eps_true,
            HouseSystem::Gauquelin,
            [x0.data[0], lat],
            None,
        )
    }

    /// Full Gauquelin sector dispatcher at `t_ut` (UT1). Routes imeth 0/1 to geometric,
    /// imeth 2–5 to rise/set-based method. `geopos` = \[lon, lat, height\].
    /// `atpress` in hPa, `attemp` in °C.
    #[doc(alias = "swe_gauquelin_sector")]
    #[allow(clippy::too_many_arguments)]
    pub fn gauquelin_sector(
        &self,
        t_ut: f64,
        body: Body,
        starname: Option<&str>,
        flags: CalcFlags,
        imeth: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
    ) -> Result<f64, Error> {
        if !(0..=5).contains(&imeth) {
            return Err(Error::CError(format!("invalid method: {imeth}")));
        }

        // §1: asteroid-numbered Pluto → Body::Pluto (swecl.c:6344-6345)
        let body = match body {
            Body::Asteroid(id) if id.mpc_number() == 134340 => Body::Pluto,
            _ => body,
        };

        if imeth <= 1 {
            self.gauquelin_sector_geometric(
                t_ut, body, starname, imeth, flags, geopos[0], geopos[1],
            )
        } else {
            self.gauquelin_sector_risetrans(
                t_ut, body, starname, flags, imeth, geopos, atpress, attemp,
            )
        }
    }

    /// Rise/set-based Gauquelin sector (imeth 2–5). Finds the bracketing rise and set
    /// times around `t_ut`, then linearly interpolates into sectors 1–36.
    /// Port of swecl.c:6370-6438, docs/c-ref-gauquelin-riseset.md §2/§7.
    #[allow(clippy::too_many_arguments)]
    fn gauquelin_sector_risetrans(
        &self,
        t_ut: f64,
        body: Body,
        starname: Option<&str>,
        flags: CalcFlags,
        imeth: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
    ) -> Result<f64, Error> {
        use crate::flags::RiseSetFlags;

        let epheflag = flags & crate::calc::EPHMASK;

        // §2: derive rise method flags from imeth
        let mut risemeth = RiseSetFlags::empty();
        if imeth == 2 || imeth == 4 {
            risemeth |= RiseSetFlags::NO_REFRACTION;
        }
        if imeth == 2 || imeth == 3 {
            risemeth |= RiseSetFlags::DISC_CENTER;
        }

        let mut tret = [0.0f64; 2]; // [0] = rise, [1] = set
        let mut rise_found = true;
        let mut set_found = true;
        let above_horizon;

        // §7.1: find next rising
        match self.rise_trans(
            t_ut,
            body,
            starname,
            epheflag,
            RiseSetFlags::RISE | risemeth,
            geopos,
            atpress,
            attemp,
        ) {
            Ok(r) => tret[0] = r.time,
            Err(Error::CircumpolarBody) => rise_found = false,
            Err(e) => return Err(e),
        }

        // §7.2: find next setting
        match self.rise_trans(
            t_ut,
            body,
            starname,
            epheflag,
            RiseSetFlags::SET | risemeth,
            geopos,
            atpress,
            attemp,
        ) {
            Ok(r) => tret[1] = r.time,
            Err(Error::CircumpolarBody) => set_found = false,
            Err(e) => return Err(e),
        }

        // §7.3: bracket determination + one re-search
        if tret[0] < tret[1] && rise_found {
            above_horizon = false;
            let t = if set_found { tret[1] - 1.2 } else { t_ut - 1.2 };
            set_found = true;
            match self.rise_trans(
                t,
                body,
                starname,
                epheflag,
                RiseSetFlags::SET | risemeth,
                geopos,
                atpress,
                attemp,
            ) {
                Ok(r) => tret[1] = r.time,
                Err(Error::CircumpolarBody) => set_found = false,
                Err(e) => return Err(e),
            }
        } else if tret[0] >= tret[1] && set_found {
            above_horizon = true;
            let t = if rise_found {
                tret[0] - 1.2
            } else {
                t_ut - 1.2
            };
            rise_found = true;
            match self.rise_trans(
                t,
                body,
                starname,
                epheflag,
                RiseSetFlags::RISE | risemeth,
                geopos,
                atpress,
                attemp,
            ) {
                Ok(r) => tret[0] = r.time,
                Err(Error::CircumpolarBody) => rise_found = false,
                Err(e) => return Err(e),
            }
        } else {
            above_horizon = false;
        }

        // §7.4: sector interpolation or failure
        if rise_found && set_found {
            if above_horizon {
                Ok((t_ut - tret[0]) / (tret[1] - tret[0]) * 18.0 + 1.0)
            } else {
                Ok((t_ut - tret[1]) / (tret[0] - tret[1]) * 18.0 + 19.0)
            }
        } else {
            Err(Error::CError(format!(
                "rise or set not found for planet {}",
                body.to_raw_id()
            )))
        }
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
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6], CalcFlags), Error> {
        let body = crate::calc::normalize_asteroid_aliases(body);
        #[cfg_attr(not(feature = "swisseph-files"), allow(unused_variables))]
        let (body, moon_raw, flags) = crate::calc::normalize_center_body(body, flags);
        let models = &config.astro_models;

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
            // `x2000` carries the J2000 equatorial vector for the SIDEREAL ECL_T0 /
            // SSY_PLANE rigorous branches (see `mean_element_pipeline`); all-zero for
            // non-sidereal calls, which `apply_sidereal` never reads.
            let (xr, x2000) = match body {
                Body::MeanNode => crate::calc::calc_mean_node(jd_tt, flags, models)?,
                Body::MeanApogee => crate::calc::calc_mean_apogee(jd_tt, flags, models)?,
                _ => unreachable!(),
            };
            return Ok((xr, x2000, flags));
        }

        if matches!(body, Body::TrueNode | Body::OscuApogee) {
            // Heliocentric/barycentric node/apogee is meaningless — C returns a
            // zeroed output (sweph.c:931-967, the HELCTR|BARYCTR guard).
            if flags.intersects(CalcFlags::HELCTR | CalcFlags::BARYCTR) {
                return Ok(([0.0; 24], [0.0; 6], flags));
            }
            // D.1: three raw geocentric moon samples (with Swiss->Moshier fallback).
            let (samples, istart, speed_intv, source) =
                self.osc_moon_samples(jd_tt, flags, config)?;
            // D.2-D.4: node and apogee are computed together; keep the requested half.
            let (node_out, apog_out) =
                crate::calc::lunar_osc_elem(jd_tt, flags, models, &samples, istart, speed_intv);
            // `x2000` carries the J2000 equatorial vector for the SIDEREAL ECL_T0 /
            // SSY_PLANE rigorous branches; `apply_sidereal` consumes it (all-zero
            // for non-sidereal calls, which never read it).
            let (mut xr, x2000) = match body {
                Body::TrueNode => node_out,
                Body::OscuApogee => apog_out,
                _ => unreachable!(),
            };
            // D.5: the true node is on the ecliptic by definition — force exact
            // zero latitude/z when neither SIDEREAL nor J2000 (suppress FP noise).
            // No zeroing for the apogee (not constrained to the ecliptic).
            if body == Body::TrueNode
                && !flags.contains(CalcFlags::SIDEREAL)
                && !flags.contains(CalcFlags::J2000)
            {
                xr[1] = 0.0;
                xr[4] = 0.0;
                xr[8] = 0.0;
                xr[11] = 0.0;
            }
            let flags_used = if source == config.ephemeris_source {
                flags
            } else {
                (flags & !crate::calc::EPHMASK) | CalcFlags::MOSEPH
            };
            return Ok((xr, x2000, flags_used));
        }

        if let Body::Fictitious(fid) = body {
            let ipl = (fid.raw_id() - crate::constants::FICT_OFFSET) as usize;
            let eps_j2000 =
                crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);
            let catalog = &self.fictitious_catalog;
            return match config.ephemeris_source {
                #[cfg(feature = "swisseph-files")]
                EphemerisSource::Swiss => {
                    match crate::calc::calc_fictitious_sweph(
                        jd_tt,
                        body,
                        catalog,
                        ipl,
                        &self.planet_files,
                        &self.moon_files,
                        &eps_j2000,
                        flags,
                        config,
                        models,
                    ) {
                        Ok((xr, x2000)) => Ok((xr, x2000, flags)),
                        Err(Error::BeyondEphemerisLimits { .. }) => {
                            let fallback_flags = (flags & !CalcFlags::SWIEPH) | CalcFlags::MOSEPH;
                            let (xr, x2000) = crate::calc::calc_fictitious_moshier(
                                jd_tt,
                                body,
                                catalog,
                                ipl,
                                &eps_j2000,
                                fallback_flags,
                                config,
                                models,
                            )?;
                            Ok((xr, x2000, fallback_flags))
                        }
                        Err(e) => Err(e),
                    }
                }
                #[cfg(not(feature = "swisseph-files"))]
                EphemerisSource::Swiss => {
                    unreachable!(".se1 file support not compiled into swisseph-rs")
                }
                #[cfg(feature = "jpl")]
                EphemerisSource::Jpl => {
                    let jpl = self.jpl_file.as_ref().ok_or(Error::EphemerisNotAvailable {
                        body,
                        source: EphemerisSource::Jpl,
                    })?;
                    let (xr, x2000) = crate::calc::calc_fictitious_jpl(
                        jd_tt, body, catalog, ipl, jpl, &eps_j2000, flags, config, models,
                    )?;
                    Ok((xr, x2000, flags))
                }
                #[cfg(not(feature = "jpl"))]
                EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
                EphemerisSource::Moshier => {
                    if flags.contains(CalcFlags::BARYCTR) {
                        return Err(Error::UnsupportedFlags(CalcFlags::BARYCTR));
                    }
                    let (xr, x2000) = crate::calc::calc_fictitious_moshier(
                        jd_tt, body, catalog, ipl, &eps_j2000, flags, config, models,
                    )?;
                    Ok((xr, x2000, flags))
                }
            };
        }

        if body == Body::Chiron
            && !(crate::constants::CHIRON_START..=crate::constants::CHIRON_END).contains(&jd_tt)
        {
            return Err(Error::BeyondEphemerisLimits {
                jd_tt,
                start: crate::constants::CHIRON_START,
                end: crate::constants::CHIRON_END,
            });
        }
        if body == Body::Pholus
            && !(crate::constants::PHOLUS_START..=crate::constants::PHOLUS_END).contains(&jd_tt)
        {
            return Err(Error::BeyondEphemerisLimits {
                jd_tt,
                start: crate::constants::PHOLUS_START,
                end: crate::constants::PHOLUS_END,
            });
        }

        // Heliocentric Sun is the origin (the Sun relative to itself) — C's swe_calc returns an
        // all-zero xx (position and speed) across every output frame, since the zero vector is
        // invariant under bias/precession/nutation and polar-converts to zeros. calc_sun has no
        // HELCTR branch for the Sun itself, so short-circuit here, matching the MeanNode/MeanApogee
        // HELCTR-zeros handling above.
        if body == Body::Sun && flags.contains(CalcFlags::HELCTR) {
            return Ok(([0.0; 24], [0.0; 6], flags));
        }

        let eps_j2000 =
            crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);

        // C's main_planet only opens the moon file when ipli >= SE_MARS (4) &&
        // ipli <= SE_PLUTO (9). For parents Sun..Venus with CENTER_BODY set but
        // suffix != 99, the flag survives inertly and the ordinary planet path runs.
        #[cfg_attr(not(feature = "swisseph-files"), allow(unused_variables))]
        let parent_raw = body.to_raw_id();
        #[cfg(feature = "swisseph-files")]
        if let Some(moon_id) = moon_raw.filter(|_| (4..=9).contains(&parent_raw)) {
            let (moon_file, parent) = self.planet_moon_file_for(body, moon_id, jd_tt)?;
            return match config.ephemeris_source {
                EphemerisSource::Swiss => {
                    let (xr, x2000) = crate::calc::calc_plmoon_sweph(
                        jd_tt,
                        body,
                        moon_file,
                        moon_id,
                        parent,
                        &self.planet_files,
                        &self.moon_files,
                        &eps_j2000,
                        flags,
                        config,
                        models,
                    )?;
                    Ok((xr, x2000, flags))
                }
                #[cfg(feature = "jpl")]
                EphemerisSource::Jpl => {
                    let jpl = self.jpl_file.as_ref().ok_or(Error::EphemerisNotAvailable {
                        body,
                        source: EphemerisSource::Jpl,
                    })?;
                    let (xr, x2000) = crate::calc::calc_plmoon_jpl(
                        jd_tt, body, moon_file, moon_id, parent, jpl, &eps_j2000, flags, config,
                        models,
                    )?;
                    Ok((xr, x2000, flags))
                }
                #[cfg(not(feature = "jpl"))]
                EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
                EphemerisSource::Moshier => {
                    if flags.contains(CalcFlags::BARYCTR) {
                        return Err(Error::UnsupportedFlags(CalcFlags::BARYCTR));
                    }
                    let (xr, x2000) = crate::calc::calc_plmoon_moshier(
                        jd_tt, body, moon_file, moon_id, parent, &eps_j2000, flags, config, models,
                    )?;
                    Ok((xr, x2000, flags))
                }
            };
        }

        match config.ephemeris_source {
            #[cfg(feature = "swisseph-files")]
            EphemerisSource::Swiss => {
                match self.calc_body_sweph(jd_tt, body, &eps_j2000, flags, models, config) {
                    Ok((xr, x2000)) => Ok((xr, x2000, flags)),
                    Err(Error::BeyondEphemerisLimits { .. }) => {
                        let fallback_flags = (flags & !CalcFlags::SWIEPH) | CalcFlags::MOSEPH;
                        let (xr, x2000) = self.calc_body_moshier(
                            jd_tt,
                            body,
                            &eps_j2000,
                            fallback_flags,
                            models,
                            config,
                        )?;
                        Ok((xr, x2000, fallback_flags))
                    }
                    Err(e) => Err(e),
                }
            }
            #[cfg(not(feature = "swisseph-files"))]
            EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            EphemerisSource::Jpl => {
                let (xr, x2000) =
                    self.calc_body_jpl(jd_tt, body, &eps_j2000, flags, models, config)?;
                Ok((xr, x2000, flags))
            }
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
            EphemerisSource::Moshier => {
                let (xr, x2000) =
                    self.calc_body_moshier(jd_tt, body, &eps_j2000, flags, models, config)?;
                Ok((xr, x2000, flags))
            }
        }
    }

    #[cfg(feature = "swisseph-files")]
    fn planet_moon_file_for(
        &self,
        body: Body,
        moon_id: i32,
        jd: f64,
    ) -> Result<(&crate::sweph_file::SwissEphFile, Body), Error> {
        let file = crate::sweph_file::find_file_for_jd(&self.planet_moon_files, moon_id, jd)
            .ok_or(Error::EphemerisNotAvailable {
                body,
                source: self.config.ephemeris_source,
            })?;
        let parent_raw = (moon_id - 9000) / 100;
        let parent = Body::try_from(parent_raw).expect("parent planet is valid");
        Ok((file, parent))
    }

    #[cfg(feature = "swisseph-files")]
    fn asteroid_file_for(
        &self,
        body: Body,
        jd: f64,
    ) -> Result<(&crate::sweph_file::SwissEphFile, i32), Error> {
        let sei_id = crate::sweph_file::body_file_id(body).ok_or(Error::EphemerisNotAvailable {
            body,
            source: self.config.ephemeris_source,
        })?;
        let files = match body {
            Body::Chiron | Body::Pholus | Body::Ceres | Body::Pallas | Body::Juno | Body::Vesta => {
                &self.main_asteroid_files
            }
            Body::Asteroid(_) => &self.asteroid_files,
            _ => {
                return Err(Error::EphemerisNotAvailable {
                    body,
                    source: self.config.ephemeris_source,
                });
            }
        };
        let file = crate::sweph_file::find_file_for_jd(files, sei_id, jd).ok_or(
            Error::EphemerisNotAvailable {
                body,
                source: self.config.ephemeris_source,
            },
        )?;
        Ok((file, sei_id))
    }

    fn calc_body_moshier(
        &self,
        jd_tt: f64,
        body: Body,
        eps_j2000: &crate::types::Epsilon,
        flags: CalcFlags,
        models: &crate::types::AstroModels,
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6]), Error> {
        if flags.contains(CalcFlags::BARYCTR) {
            return Err(Error::UnsupportedFlags(CalcFlags::BARYCTR));
        }
        match body {
            Body::Sun | Body::Earth => {
                crate::calc::calc_sun(jd_tt, eps_j2000, flags, config, models, body == Body::Earth)
            }
            Body::Moon => crate::calc::calc_moon(jd_tt, eps_j2000, flags, config, models),
            Body::Mercury
            | Body::Venus
            | Body::Mars
            | Body::Jupiter
            | Body::Saturn
            | Body::Uranus
            | Body::Neptune
            | Body::Pluto => {
                crate::calc::calc_planet(jd_tt, body, eps_j2000, flags, config, models)
            }
            #[cfg(feature = "swisseph-files")]
            Body::Chiron
            | Body::Pholus
            | Body::Ceres
            | Body::Pallas
            | Body::Juno
            | Body::Vesta
            | Body::Asteroid(_) => {
                let (f, id) = self.asteroid_file_for(body, jd_tt)?;
                crate::calc::calc_asteroid_moshier(
                    jd_tt, body, f, id, eps_j2000, flags, config, models,
                )
            }
            #[cfg(not(feature = "swisseph-files"))]
            Body::Chiron
            | Body::Pholus
            | Body::Ceres
            | Body::Pallas
            | Body::Juno
            | Body::Vesta
            | Body::Asteroid(_) => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Moshier,
            }),
            _ => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Moshier,
            }),
        }
    }

    #[cfg(feature = "swisseph-files")]
    fn calc_body_sweph(
        &self,
        jd_tt: f64,
        body: Body,
        eps_j2000: &crate::types::Epsilon,
        flags: CalcFlags,
        models: &crate::types::AstroModels,
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6]), Error> {
        match body {
            Body::Sun | Body::Earth => crate::calc::calc_sun_sweph(
                jd_tt,
                &self.planet_files,
                &self.moon_files,
                flags,
                config,
                models,
                body == Body::Earth,
            ),
            Body::Moon => crate::calc::calc_moon_sweph(
                jd_tt,
                &self.planet_files,
                &self.moon_files,
                flags,
                config,
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
                config,
                models,
            ),
            Body::Chiron
            | Body::Pholus
            | Body::Ceres
            | Body::Pallas
            | Body::Juno
            | Body::Vesta
            | Body::Asteroid(_) => {
                let (f, id) = self.asteroid_file_for(body, jd_tt)?;
                crate::calc::calc_asteroid_sweph(
                    jd_tt,
                    body,
                    f,
                    id,
                    &self.planet_files,
                    &self.moon_files,
                    eps_j2000,
                    flags,
                    config,
                    models,
                )
            }
            _ => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Swiss,
            }),
        }
    }

    #[cfg(feature = "jpl")]
    fn calc_body_jpl(
        &self,
        jd_tt: f64,
        body: Body,
        eps_j2000: &crate::types::Epsilon,
        flags: CalcFlags,
        models: &crate::types::AstroModels,
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6]), Error> {
        let file = self
            .jpl_file
            .as_ref()
            .expect("JPL file must be initialized into Ephemeris upon calling Ephemeris::new()");
        match body {
            Body::Sun | Body::Earth => {
                crate::calc::calc_sun_jpl(jd_tt, file, flags, config, models, body == Body::Earth)
            }
            Body::Moon => crate::calc::calc_moon_jpl(jd_tt, file, flags, config, models),
            Body::Mercury
            | Body::Venus
            | Body::Mars
            | Body::Jupiter
            | Body::Saturn
            | Body::Uranus
            | Body::Neptune
            | Body::Pluto => {
                crate::calc::calc_planet_jpl(jd_tt, body, file, eps_j2000, flags, config, models)
            }
            #[cfg(feature = "swisseph-files")]
            Body::Chiron
            | Body::Pholus
            | Body::Ceres
            | Body::Pallas
            | Body::Juno
            | Body::Vesta
            | Body::Asteroid(_) => {
                let (f, id) = self.asteroid_file_for(body, jd_tt)?;
                crate::calc::calc_asteroid_jpl(
                    jd_tt, body, f, id, file, eps_j2000, flags, config, models,
                )
            }
            #[cfg(not(feature = "swisseph-files"))]
            Body::Chiron
            | Body::Pholus
            | Body::Ceres
            | Body::Pallas
            | Body::Juno
            | Body::Vesta
            | Body::Asteroid(_) => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Jpl,
            }),
            _ => Err(Error::EphemerisNotAvailable {
                body,
                source: EphemerisSource::Jpl,
            }),
        }
    }

    /// Raw geocentric equatorial-J2000 (pre-bias) moon (pos+vel) at `t` from the
    /// given backend, for the osculating-node/apogee D.1 sample loop.
    fn raw_moon_at(
        &self,
        source: EphemerisSource,
        t: f64,
        eps_j2000: &crate::types::Epsilon,
    ) -> Result<[f64; 6], Error> {
        match source {
            EphemerisSource::Moshier => crate::calc::raw_osc_moon_moshier(t, eps_j2000),
            #[cfg(feature = "swisseph-files")]
            EphemerisSource::Swiss => crate::calc::raw_osc_moon_sweph(&self.moon_files, t),
            #[cfg(not(feature = "swisseph-files"))]
            EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            EphemerisSource::Jpl => crate::calc::raw_osc_moon_jpl(
                self.jpl_file.as_ref().expect(
                    "JPL file must be initialized into Ephemeris upon calling Ephemeris::new()",
                ),
                t,
            ),
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
        }
    }

    /// Fetch the three light-time-corrected moon samples for one backend
    /// (D.1 loop body, sweph.c:5251-5357). Only `[istart..=2]` are populated.
    fn fetch_osc_samples(
        &self,
        source: EphemerisSource,
        tjd: f64,
        istart: usize,
        speed_intv: f64,
        truepos: bool,
        eps_j2000: &crate::types::Epsilon,
    ) -> Result<[[f64; 6]; 3], Error> {
        let mut samples = [[0.0f64; 6]; 3];
        for (i, slot) in samples.iter_mut().enumerate().skip(istart) {
            let t = match i {
                0 => tjd - speed_intv,
                1 => tjd + speed_intv,
                _ => tjd,
            };
            let mut m = self.raw_moon_at(source, t, eps_j2000)?;
            // Light-time-corrected moon for the apparent node (full re-evaluation
            // at t-dt, NOT the cheap `x -= dt*speed`; C insists on this). C ONLY
            // does this for the JPL and Swiss branches — the Moshier branch
            // (sweph.c:5336-5354) has NO light-time block, so the Moshier moon is
            // used geometrically. Match that exactly.
            if !truepos && source != EphemerisSource::Moshier {
                let dist = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
                let dt = dist * crate::constants::AUNIT / crate::constants::CLIGHT / 86400.0;
                m = self.raw_moon_at(source, t - dt, eps_j2000)?;
            }
            *slot = m;
        }
        Ok(samples)
    }

    /// Obtain the three backend moon samples for the osculating node/apogee,
    /// with the Swiss->Moshier fallback (mirrors `calc_inner`; JPL does not fall
    /// back, matching the rest of this port). Returns the samples, `istart`, the
    /// backend-specific `speed_intv`, and the backend actually used.
    fn osc_moon_samples(
        &self,
        tjd: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<OscMoonSamples, Error> {
        let istart = if flags.contains(CalcFlags::SPEED) {
            0
        } else {
            2
        };
        let truepos = flags.contains(CalcFlags::TRUEPOS);
        let models = &config.astro_models;
        let eps_j2000 =
            crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);

        match config.ephemeris_source {
            #[cfg(feature = "swisseph-files")]
            EphemerisSource::Swiss => {
                let si = crate::constants::NODE_CALC_INTV;
                match self.fetch_osc_samples(
                    EphemerisSource::Swiss,
                    tjd,
                    istart,
                    si,
                    truepos,
                    &eps_j2000,
                ) {
                    Ok(s) => Ok((s, istart, si, EphemerisSource::Swiss)),
                    Err(Error::BeyondEphemerisLimits { .. }) => {
                        let sim = crate::constants::NODE_CALC_INTV_MOSH;
                        let s = self.fetch_osc_samples(
                            EphemerisSource::Moshier,
                            tjd,
                            istart,
                            sim,
                            truepos,
                            &eps_j2000,
                        )?;
                        Ok((s, istart, sim, EphemerisSource::Moshier))
                    }
                    Err(e) => Err(e),
                }
            }
            #[cfg(not(feature = "swisseph-files"))]
            EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            EphemerisSource::Jpl => {
                let si = crate::constants::NODE_CALC_INTV;
                let s = self.fetch_osc_samples(
                    EphemerisSource::Jpl,
                    tjd,
                    istart,
                    si,
                    truepos,
                    &eps_j2000,
                )?;
                Ok((s, istart, si, EphemerisSource::Jpl))
            }
            #[cfg(not(feature = "jpl"))]
            EphemerisSource::Jpl => unreachable!("JPL support not compiled into swisseph-rs"),
            EphemerisSource::Moshier => {
                let si = crate::constants::NODE_CALC_INTV_MOSH;
                let s = self.fetch_osc_samples(
                    EphemerisSource::Moshier,
                    tjd,
                    istart,
                    si,
                    truepos,
                    &eps_j2000,
                )?;
                Ok((s, istart, si, EphemerisSource::Moshier))
            }
        }
    }

    fn calc_speed3(
        &self,
        jd_tt: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<CalcResult, Error> {
        let dt = crate::calc::speed3_interval(body);
        let inner_flags = flags & !CalcFlags::SPEED3;

        let (mut x0, x2000_0, _) = self.calc_inner(jd_tt - dt, body, inner_flags, config)?;
        let (mut x2, x2000_2, _) = self.calc_inner(jd_tt + dt, body, inner_flags, config)?;
        let (mut x1, x2000_1, flags_used) = self.calc_inner(jd_tt, body, inner_flags, config)?;

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
            self.apply_sidereal(&mut x0, &x2000_0, jd_tt - dt, pos_flags, config)?;
            self.apply_sidereal(&mut x2, &x2000_2, jd_tt + dt, pos_flags, config)?;
            self.apply_sidereal(&mut x1, &x2000_1, jd_tt, pos_flags, config)?;
        }

        crate::calc::denormalize_positions(&mut x0, &x1, &mut x2);
        crate::calc::calc_speed_3point(&mut x1, &x0, &x2, dt);

        Ok(CalcResult {
            data: Self::extract_for_body(&x1, body, flags | CalcFlags::SPEED),
            flags_used,
        })
    }

    pub(crate) fn apply_sidereal(
        &self,
        xreturn: &mut [f64; 24],
        x2000: &[f64; 6],
        jd_tt: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<(), Error> {
        use crate::constants::RADTODEG;
        use crate::math::cartesian_to_polar_with_speed;

        let bits = config.sidereal_bits;
        let models = &config.astro_models;
        let has_speed = flags.contains(CalcFlags::SPEED);
        let has_meaningful_x2000 = *x2000 != [0.0f64; 6];

        if has_meaningful_x2000 && bits.contains(SiderealBits::ECL_T0) {
            let (xecl, xequ) = crate::ayanamsa::trop_ra2sid_lon(x2000, config, models, flags);

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
            let xecl = crate::ayanamsa::trop_ra2sid_lon_sosy(x2000, config, models, flags);

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
            let idx = crate::ayanamsa::sidereal_index(config);
            let (daya_val, daya_sp) = if crate::ayanamsa::FIXED_STAR_INDICES.contains(&idx) {
                self.fixstar_ayanamsa(jd_tt, flags, config)?
            } else {
                let a = crate::ayanamsa::get_ayanamsa_with_speed(config, jd_tt, flags, models)?;
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
    /// `star` is searched case-insensitively in the star catalog. Returns
    /// `(canonical_name, CalcResult)` where the name is `"traditional,bayer"`.
    #[doc(alias = "swe_fixstar2")]
    pub fn fixstar2(
        &self,
        star: &str,
        jd_tt: f64,
        flags: CalcFlags,
    ) -> Result<(String, CalcResult), Error> {
        self.fixstar2_with_config(star, jd_tt, flags, &self.config)
    }

    /// Same as [`fixstar2`](Self::fixstar2) but with an explicit config override -- see
    /// [`calc_with_config`](Self::calc_with_config). Lets callers (e.g. `eclipse.rs`'s
    /// `eclipse_how`/`occult_when_loc`) get a topocentric fixed-star position at a caller-supplied
    /// `geopos` without requiring it to match the `Ephemeris`'s own configured topographic
    /// position.
    pub fn fixstar2_with_config(
        &self,
        star: &str,
        jd_tt: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<(String, CalcResult), Error> {
        let config = self.effective_config(flags, config);
        // C's swe_fixstar2 returns the original input iflag unchanged (it passes
        // iflag by value to fixstar_calc_from_struct and ignores the return).
        let orig_flags = flags;
        let flags = crate::calc::plaus_iflag(flags, config.ephemeris_source);
        let resolved = if let Some(s) = crate::stars::builtin_star(star) {
            s
        } else {
            self.stars.search(star)?
        };
        let data = self.calc_fixstar(&resolved, jd_tt, flags, &config)?;
        let name = format!("{},{}", resolved.name, resolved.bayer);
        Ok((
            name,
            CalcResult {
                data,
                flags_used: orig_flags,
            },
        ))
    }

    /// UT-based [`fixstar2`](Self::fixstar2).
    #[doc(alias = "swe_fixstar2_ut")]
    pub fn fixstar2_ut(
        &self,
        star: &str,
        jd_ut: f64,
        flags: CalcFlags,
    ) -> Result<(String, CalcResult), Error> {
        let config = self.effective_config(flags, &self.config);
        let dt = crate::deltat::calc_deltat(jd_ut, &config);
        self.fixstar2(star, jd_ut + dt, flags)
    }

    /// Magnitude lookup for a star by name. Searches the catalog file only (built-in reference
    /// stars are not available via this function). Returns `(canonical_name, magnitude)`.
    #[doc(alias = "swe_fixstar2_mag")]
    pub fn fixstar2_mag(&self, star: &str) -> Result<(String, f64), Error> {
        let resolved = self.stars.search(star)?;
        let name = format!("{},{}", resolved.name, resolved.bayer);
        Ok((name, resolved.mag))
    }

    /// Dispatcher: routes fixed-star computation to the correct backend. `config` is threaded
    /// explicitly (rather than always reading `self.config`) so callers needing a topocentric
    /// position (e.g. `eclipse.rs`'s `eclipse_how`/`occult_when_loc`, which build a per-call
    /// `topo_config` override) can get one -- mirrors `calc_with_config`'s pattern.
    fn calc_fixstar(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<[f64; 6], Error> {
        match config.ephemeris_source {
            #[cfg(feature = "swisseph-files")]
            crate::types::EphemerisSource::Swiss => {
                self.calc_fixstar_sweph(star, jd, flags, config)
            }
            #[cfg(not(feature = "swisseph-files"))]
            crate::types::EphemerisSource::Swiss => {
                unreachable!(".se1 file support not compiled into swisseph-rs")
            }
            #[cfg(feature = "jpl")]
            crate::types::EphemerisSource::Jpl => self.calc_fixstar_jpl(star, jd, flags, config),
            #[cfg(not(feature = "jpl"))]
            crate::types::EphemerisSource::Jpl => {
                unreachable!("JPL support not compiled into swisseph-rs")
            }
            crate::types::EphemerisSource::Moshier => {
                self.calc_fixstar_moshier(star, jd, flags, config)
            }
        }
    }

    /// Moshier backend: computes heliocentric Earth via Moshier pipeline. `xobs`/`xobs_dt` get a
    /// topocentric offset added when `TOPOCTR` is set (docs/c-ref-fixstar.md step 6) -- previously
    /// silently ignored here (no golden coverage exercised it), matching the same offset-addition
    /// pattern as `calc.rs`'s `apparent_planet`/`apparent_sun`/`apparent_moon`.
    fn calc_fixstar_moshier(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<[f64; 6], Error> {
        use crate::constants::{FIXSTAR_DT, J2000};
        use crate::obliquity::obliquity;

        let models = &config.astro_models;
        // Moshier returns heliocentric Earth, matching C's xearth for MOSEPH.
        let eps_j2000 = obliquity(J2000, CalcFlags::empty(), models);
        let pp =
            crate::moshier::backend::compute_pipeline(jd, crate::types::Body::Sun, &eps_j2000)?;
        let mut xobs = pp.earth_helio;
        let pp_dt = crate::moshier::backend::compute_pipeline(
            jd - FIXSTAR_DT,
            crate::types::Body::Sun,
            &eps_j2000,
        )?;
        let mut xobs_dt = pp_dt.earth_helio;
        let offset = crate::calc::topo_offset(jd, flags, config, models);
        let offset_dt = crate::calc::topo_offset(jd - FIXSTAR_DT, flags, config, models);
        for i in 0..6 {
            xobs[i] += offset[i];
            xobs_dt[i] += offset_dt[i];
        }
        // Moshier is heliocentric; Sun is at the origin, so sun_bary = 0.
        let sun_bary = [0.0f64; 6];
        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary, config)
    }

    /// SWIEPH backend: barycentric Earth for parallax/aberration, sun_bary for deflection.
    #[cfg(feature = "swisseph-files")]
    fn calc_fixstar_sweph(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
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
        let mut xobs = pp.earth_bary;
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
        let mut xobs_dt = pp_dt.earth_bary;

        let models = &config.astro_models;
        let offset = crate::calc::topo_offset(jd, flags, config, models);
        let offset_dt = crate::calc::topo_offset(jd_dt, flags, config, models);
        for i in 0..6 {
            xobs[i] += offset[i];
            xobs_dt[i] += offset_dt[i];
        }

        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary, config)
    }

    /// JPL backend: barycentric Earth for parallax/aberration, sun_bary for deflection.
    #[cfg(feature = "jpl")]
    fn calc_fixstar_jpl(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<[f64; 6], Error> {
        use crate::constants::FIXSTAR_DT;
        use crate::jpl::{J_EARTH, J_SBARY, J_SUN, jpl_pleph};

        let file = self.jpl_file.as_ref().ok_or(Error::EphemerisNotAvailable {
            body: crate::types::Body::Sun,
            source: crate::types::EphemerisSource::Jpl,
        })?;

        // C uses barycentric Earth for parallax/aberration; deflection uses earth_helio
        // computed inside swi_deflect_light as earth_bary - sun_bary.
        let mut xobs = jpl_pleph(file, jd, J_EARTH, J_SBARY, true)?;
        let sun_bary = jpl_pleph(file, jd, J_SUN, J_SBARY, true)?;
        let mut xobs_dt = jpl_pleph(file, jd - FIXSTAR_DT, J_EARTH, J_SBARY, true)?;
        let models = &config.astro_models;
        let offset = crate::calc::topo_offset(jd, flags, config, models);
        let offset_dt = crate::calc::topo_offset(jd - FIXSTAR_DT, flags, config, models);
        for i in 0..6 {
            xobs[i] += offset[i];
            xobs_dt[i] += offset_dt[i];
        }

        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary, config)
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
    #[allow(clippy::too_many_arguments)]
    fn calc_fixstar_inner(
        &self,
        star: &crate::stars::Star,
        jd: f64,
        flags: CalcFlags,
        xobs: [f64; 6],
        xobs_dt: [f64; 6],
        sun_bary: [f64; 6],
        config: &EphemerisConfig,
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

        let models = &config.astro_models;
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
            nutate(&mut x, &eps_date, &nv, Some(&nutv), true, false);
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
            let bits = config.sidereal_bits;
            if bits.contains(crate::flags::SiderealBits::ECL_T0) {
                let (xecl, xequ) = crate::ayanamsa::trop_ra2sid_lon(&xxsv, config, models, iflag);
                x = if iflag.contains(CalcFlags::EQUATORIAL) {
                    xequ
                } else {
                    xecl
                };
            } else if bits.contains(crate::flags::SiderealBits::SSY_PLANE) {
                let xecl = crate::ayanamsa::trop_ra2sid_lon_sosy(&xxsv, config, models, iflag);
                x = xecl;
            } else {
                // Default: subtract ayanamsa from ecliptic (or equatorial) longitude.
                x = cartesian_to_polar_with_speed(x);
                let idx = crate::ayanamsa::sidereal_index(config);
                let (daya_val, daya_sp) = if crate::ayanamsa::FIXED_STAR_INDICES.contains(&idx) {
                    self.fixstar_ayanamsa(jd, iflag, config)?
                } else {
                    let a = crate::ayanamsa::get_ayanamsa_with_speed(config, jd, iflag, models)?;
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
    fn fixstar_ayanamsa_single(
        &self,
        jd_tt: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<f64, Error> {
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

        let idx = crate::ayanamsa::sidereal_index(config);

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
                let eps_deg = crate::obliquity::obliquity(jd_tt, iflag_base, &config.astro_models)
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
    fn fixstar_ayanamsa(
        &self,
        jd_tt: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<(f64, f64), Error> {
        const TINTV: f64 = 0.001;
        let d0 = self.fixstar_ayanamsa_single(jd_tt, flags, config)?;
        let d2 = self.fixstar_ayanamsa_single(jd_tt - TINTV, flags, config)?;
        // Both samples are independently normalized to [0,360); use the signed
        // shortest difference so a 360° wrap between samples doesn't blow up the
        // speed (~3.6e5 deg/day spike). diff_degrees returns a value in (-180,180].
        Ok((d0, crate::math::diff_degrees(d0, d2) / TINTV))
    }

    fn build_leap_seconds(config: &EphemerisConfig) -> crate::Result<Vec<i32>> {
        let last_hardcoded = *LEAP_SECONDS.last().expect(
            "LEAP_SECONDS constant array is empty; the codebase has been critically misconfigured",
        );
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

    // -----------------------------------------------------------------------
    // Equation of time / LMT↔LAT (sweph.c:7387–7436)
    // -----------------------------------------------------------------------

    /// Equation of time: `E = LAT − LMT`, returned in **days**.
    /// Equation of time at `tjd_ut` (UT1), returned in days (positive = Sun ahead of mean).
    #[doc(alias = "swe_time_equ")]
    pub fn time_equ(&self, tjd_ut: f64) -> crate::Result<f64> {
        let deltat_config = {
            let mut c = self.config.clone();
            c.tidal_acceleration = Some(crate::constants::TIDAL_DEFAULT);
            c
        };
        let sidt = crate::sidereal_time::sidereal_time(tjd_ut, &deltat_config);
        let t = tjd_ut + 0.5;
        let dt_day = t - t.floor();
        let sidt_deg = (sidt - dt_day * 24.0) * 15.0;
        let result = self.calc_ut(tjd_ut, Body::Sun, CalcFlags::EQUATORIAL)?;
        let sun_ra = result.data[0];
        let mut dt = crate::math::normalize_degrees(sidt_deg - sun_ra - 180.0);
        if dt > 180.0 {
            dt -= 360.0;
        }
        dt *= 4.0;
        Ok(dt / 1440.0)
    }

    /// Convert Local Mean Time to Local Apparent Time.
    /// `geolon` in degrees (east-positive). Both input and output are Julian Day (UT-scale).
    #[doc(alias = "swe_lmt_to_lat")]
    pub fn lmt_to_lat(&self, tjd_lmt: f64, geolon: f64) -> crate::Result<f64> {
        let tjd_lmt0 = tjd_lmt - geolon / 360.0;
        let e = self.time_equ(tjd_lmt0)?;
        Ok(tjd_lmt + e)
    }

    /// Convert Local Apparent Time to Local Mean Time.
    /// `geolon` in degrees (east-positive). Both input and output are Julian Day (UT-scale).
    #[doc(alias = "swe_lat_to_lmt")]
    pub fn lat_to_lmt(&self, tjd_lat: f64, geolon: f64) -> crate::Result<f64> {
        let tjd_lmt0 = tjd_lat - geolon / 360.0;
        let mut e = self.time_equ(tjd_lmt0)?;
        e = self.time_equ(tjd_lmt0 - e)?;
        e = self.time_equ(tjd_lmt0 - e)?;
        Ok(tjd_lat - e)
    }

    // -----------------------------------------------------------------------
    // Ephemeris file introspection (stateless swe_get_current_file_data)
    // -----------------------------------------------------------------------

    /// Return metadata about the ephemeris file that would serve a calculation at
    /// the given Julian Day for the specified file kind.
    ///
    /// This is the stateless equivalent of C's `swe_get_current_file_data(ifno)`.
    /// C reports whichever file was used by the *last* `swe_calc` call (global
    /// state). Here, `jd` selects the file explicitly — the same selection logic
    /// that `calc` uses internally.
    ///
    /// Returns `None` when:
    /// - The ephemeris source is Moshier (no files)
    /// - No file covers the given `jd`
    /// - `kind` is `Asteroid` or `PlanetMoon` (stateless: no "last-used" concept)
    #[doc(alias = "swe_get_current_file_data")]
    pub fn file_data(
        &self,
        kind: crate::types::FileDataKind,
        jd: f64,
    ) -> Option<crate::types::FileData> {
        use crate::types::FileDataKind;

        match kind {
            FileDataKind::Planet => self.file_data_planet(jd),
            FileDataKind::Moon => self.file_data_sweph(self.sweph_moon_files(), jd),
            FileDataKind::MainAsteroid => {
                self.file_data_sweph(self.sweph_main_asteroid_files(), jd)
            }
            FileDataKind::Asteroid | FileDataKind::PlanetMoon => None,
        }
    }

    /// File data for a specific body at `jd`. Unlike [`file_data`](Self::file_data),
    /// this handles individual asteroids and planet-moons by looking up the
    /// body's own `.se1` file.
    pub fn file_data_for_body(&self, body: Body, jd: f64) -> Option<crate::types::FileData> {
        match body {
            Body::Sun
            | Body::Mercury
            | Body::Venus
            | Body::Earth
            | Body::Mars
            | Body::Jupiter
            | Body::Saturn
            | Body::Uranus
            | Body::Neptune
            | Body::Pluto
            | Body::MeanNode
            | Body::TrueNode
            | Body::MeanApogee
            | Body::OscuApogee
            | Body::EclipticNutation => self.file_data_planet(jd),
            Body::Moon => self.file_data(crate::types::FileDataKind::Moon, jd),
            Body::Chiron | Body::Pholus | Body::Ceres | Body::Pallas | Body::Juno | Body::Vesta => {
                self.file_data(crate::types::FileDataKind::MainAsteroid, jd)
            }
            Body::Asteroid(_) => self.file_data_asteroid(body, jd),
            Body::PlanetMoon(_) => self.file_data_planet_moon(body, jd),
            Body::IntpApogee | Body::IntpPerigee => self.file_data_planet(jd),
            Body::Fictitious(_) => None,
        }
    }

    #[cfg(feature = "swisseph-files")]
    fn file_data_asteroid(&self, body: Body, jd: f64) -> Option<crate::types::FileData> {
        self.asteroid_file_for(body, jd).ok().map(|(f, _)| {
            let h = f.header();
            crate::types::FileData {
                path: f.path().to_path_buf(),
                start_jd: h.time_range.0,
                end_jd: h.time_range.1,
                denum: h.denum,
            }
        })
    }

    #[cfg(not(feature = "swisseph-files"))]
    fn file_data_asteroid(&self, _body: Body, _jd: f64) -> Option<crate::types::FileData> {
        None
    }

    #[cfg(feature = "swisseph-files")]
    fn file_data_planet_moon(&self, body: Body, jd: f64) -> Option<crate::types::FileData> {
        let id = match body {
            Body::PlanetMoon(id) => id,
            _ => return None,
        };
        let moon_id = crate::constants::PLMOON_OFFSET + id.encoded();
        self.planet_moon_file_for(body, moon_id, jd)
            .ok()
            .map(|(f, _)| {
                let h = f.header();
                crate::types::FileData {
                    path: f.path().to_path_buf(),
                    start_jd: h.time_range.0,
                    end_jd: h.time_range.1,
                    denum: h.denum,
                }
            })
    }

    #[cfg(not(feature = "swisseph-files"))]
    fn file_data_planet_moon(&self, _body: Body, _jd: f64) -> Option<crate::types::FileData> {
        None
    }

    fn file_data_planet(&self, jd: f64) -> Option<crate::types::FileData> {
        #[cfg(feature = "jpl")]
        if self.config.ephemeris_source == EphemerisSource::Jpl
            && let Some(ref jf) = self.jpl_file
        {
            let h = jf.header();
            if jd >= h.ss[0] && jd <= h.ss[1] {
                return Some(crate::types::FileData {
                    path: jf.path().to_path_buf(),
                    start_jd: h.ss[0],
                    end_jd: h.ss[1],
                    denum: h.denum,
                });
            }
        }
        self.file_data_sweph(self.sweph_planet_files(), jd)
    }

    #[cfg(feature = "swisseph-files")]
    fn file_data_sweph(
        &self,
        files: &[crate::sweph_file::SwissEphFile],
        jd: f64,
    ) -> Option<crate::types::FileData> {
        files
            .iter()
            .rev()
            .find(|f| {
                let (start, end) = f.header().time_range;
                start <= jd && jd <= end
            })
            .map(|f| {
                let h = f.header();
                crate::types::FileData {
                    path: f.path().to_path_buf(),
                    start_jd: h.time_range.0,
                    end_jd: h.time_range.1,
                    denum: h.denum,
                }
            })
    }

    #[cfg(not(feature = "swisseph-files"))]
    fn file_data_sweph(&self, _files: &[()], _jd: f64) -> Option<crate::types::FileData> {
        None
    }

    #[cfg(feature = "swisseph-files")]
    fn sweph_planet_files(&self) -> &[crate::sweph_file::SwissEphFile] {
        &self.planet_files
    }

    #[cfg(not(feature = "swisseph-files"))]
    fn sweph_planet_files(&self) -> &[()] {
        &[]
    }

    #[cfg(feature = "swisseph-files")]
    fn sweph_moon_files(&self) -> &[crate::sweph_file::SwissEphFile] {
        &self.moon_files
    }

    #[cfg(not(feature = "swisseph-files"))]
    fn sweph_moon_files(&self) -> &[()] {
        &[]
    }

    #[cfg(feature = "swisseph-files")]
    fn sweph_main_asteroid_files(&self) -> &[crate::sweph_file::SwissEphFile] {
        &self.main_asteroid_files
    }

    #[cfg(not(feature = "swisseph-files"))]
    fn sweph_main_asteroid_files(&self) -> &[()] {
        &[]
    }

    // -----------------------------------------------------------------------
    // Body name lookup (sweph.c:6946–7125)
    // -----------------------------------------------------------------------

    /// Resolve a body to its display name (e.g. "Sun", "Chiron", asteroid name from file).
    #[doc(alias = "swe_get_planet_name")]
    pub fn get_planet_name(&self, body: Body) -> String {
        let body = crate::calc::normalize_asteroid_aliases(body);
        match body {
            Body::Sun => "Sun".into(),
            Body::Moon => "Moon".into(),
            Body::Mercury => "Mercury".into(),
            Body::Venus => "Venus".into(),
            Body::Mars => "Mars".into(),
            Body::Jupiter => "Jupiter".into(),
            Body::Saturn => "Saturn".into(),
            Body::Uranus => "Uranus".into(),
            Body::Neptune => "Neptune".into(),
            Body::Pluto => "Pluto".into(),
            Body::MeanNode => "mean Node".into(),
            Body::TrueNode => "true Node".into(),
            Body::MeanApogee => "mean Apogee".into(),
            Body::OscuApogee => "osc. Apogee".into(),
            Body::IntpApogee => "intp. Apogee".into(),
            Body::IntpPerigee => "intp. Perigee".into(),
            Body::Earth => "Earth".into(),
            Body::Chiron => "Chiron".into(),
            Body::Pholus => "Pholus".into(),
            Body::Ceres => "Ceres".into(),
            Body::Pallas => "Pallas".into(),
            Body::Juno => "Juno".into(),
            Body::Vesta => "Vesta".into(),
            Body::EclipticNutation => "Ecl. Nut.".into(),
            Body::Fictitious(id) => crate::fictitious::fictitious_name(
                &self.fictitious_catalog,
                (id.raw_id() - crate::constants::FICT_OFFSET) as usize,
            ),
            Body::Asteroid(id) => match id.mpc_number() {
                2060 => "Chiron".into(),
                5145 => "Pholus".into(),
                mpc => self.asteroid_name(mpc),
            },
            Body::PlanetMoon(_) => format!("{}", body.to_raw_id()),
        }
    }

    fn asteroid_name(&self, mpc: i32) -> String {
        let name_from_file = self.try_asteroid_name_from_file(mpc);
        let name = match &name_from_file {
            Some(n) if !n.is_empty() => n.as_str(),
            _ => return format!("{}: not found (asteroid)", mpc),
        };

        if (name.starts_with('?') || (name.len() > 1 && name.as_bytes()[1].is_ascii_digit()))
            && let Some(override_name) = self.seasnam_lookup(mpc)
        {
            return override_name;
        }
        name.to_string()
    }

    #[cfg(feature = "swisseph-files")]
    fn try_asteroid_name_from_file(&self, mpc: i32) -> Option<String> {
        let dir = self.config.ephe_path.as_ref()?;
        match crate::sweph_file::open_asteroid_file(dir, mpc) {
            Ok(f) => f.header().asteroid.as_ref().map(|a| a.name.clone()),
            Err(_) => None,
        }
    }

    /// Without `swisseph-files`, asteroid names (which live in .se1 file headers)
    /// are never available.
    #[cfg(not(feature = "swisseph-files"))]
    fn try_asteroid_name_from_file(&self, _mpc: i32) -> Option<String> {
        None
    }

    fn seasnam_lookup(&self, mpc: i32) -> Option<String> {
        let dir = self.config.ephe_path.as_ref()?;
        let path = dir.join("seasnam.txt");
        let contents = fs::read_to_string(&path).ok()?;
        for line in contents.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let stripped = trimmed.trim_start_matches(['(', '[', '{']);
            let num_str: String = stripped
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect();
            if num_str.is_empty() {
                continue;
            }
            let file_mpc: i32 = match num_str.parse() {
                Ok(n) => n,
                Err(_) => continue,
            };
            if file_mpc != mpc {
                continue;
            }
            let rest = &stripped[num_str.len()..];
            let name_part = rest.trim_start();
            if name_part.is_empty() {
                continue;
            }
            let name = name_part
                .split(['#', '\r', '\n'])
                .next()
                .unwrap_or("")
                .trim_end();
            if !name.is_empty() {
                return Some(name.to_string());
            }
        }
        None
    }
}

impl DeltaT for Ephemeris {
    fn delta_t(&self, jd_ut: JdUt1) -> f64 {
        crate::deltat::calc_deltat(jd_ut.0, &self.config)
    }
}

/// Result of a planetary position calculation.
///
/// `data[0..3]` = position (ecliptic longitude/latitude/distance in degrees/degrees/AU by
/// default; RA/dec/dist with `EQUATORIAL`; x/y/z AU with `XYZ`; radians with `RADIANS`).
/// `data[3..6]` = speed (degrees/day or AU/day) when `SPEED` is set in the input flags.
///
/// `flags_used` reports the flags actually applied — compare with requested flags to detect
/// fallbacks (e.g. `SWIEPH` requested but `MOSEPH` used because files were unavailable).
#[derive(Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CalcResult {
    /// Position and speed array. Layout depends on flags (see struct-level docs).
    pub data: [f64; 6],
    /// Flags that were actually applied (may differ from requested).
    pub flags_used: CalcFlags,
}
