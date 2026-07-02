use std::fs;

use crate::calc::{JplProvider, PositionProvider, SwephPositions, SwephProvider};
use crate::config::EphemerisConfig;
use crate::date::LEAP_SECONDS;
use crate::error::Error;
use crate::flags::{CalcFlags, EclipseFlags, SiderealBits};
use crate::types::{Body, DeltaT, EphemerisSource, JdUt1};

/// Three raw geocentric moon samples for the osculating node/apogee, plus the
/// `istart`, backend-specific central-difference interval, and backend used.
type OscMoonSamples = ([[f64; 6]; 3], usize, f64, EphemerisSource);

/// Selects `ipli`'s position from a `SwephPositions` bundle in the frame
/// `swe_nod_aps`'s osculating branch needs: barycentric (`SE_NODBIT_OSCU_BAR`)
/// or heliocentric, with `Body::Earth` reading the always-populated
/// `earth_bary`/`earth_helio` fields rather than `planet_bary` (which, for the
/// `query = Body::Sun` dummy call `nodaps_osc_body_j2000` makes for Earth, is
/// the Sun's own position, not Earth's — see `docs/c-ref-nodaps.md` §A.4.1).
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

pub struct Ephemeris {
    config: EphemerisConfig,
    leap_seconds: Vec<i32>,
    planet_files: Vec<crate::sweph_file::SwissEphFile>,
    moon_files: Vec<crate::sweph_file::SwissEphFile>,
    jpl_file: Option<crate::jpl::JplFile>,
    stars: crate::stars::StarCatalog,
}

impl Ephemeris {
    pub fn new(mut config: EphemerisConfig) -> crate::Result<Self> {
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
        // Resolve the ephemeris-specific tidal acceleration from the open file's
        // DE number, mirroring C's `swi_get_tid_acc` (swephlib.c:3211–3221): JPL
        // uses the JPL file's denum, SWIEPH the moon (SEI_FILE_MOON) file's. This
        // is what makes ΔT — and therefore the topocentric observer offset — match
        // C away from J2000 (DE441 tid_acc differs from the DE431 default). Only
        // fill in when the caller hasn't pinned tid_acc explicitly (C's
        // `is_tid_acc_manual` short-circuit).
        if config.tidal_acceleration.is_none() {
            let denum = match config.ephemeris_source {
                EphemerisSource::Swiss => moon_files.first().map(|f| f.header().denum),
                EphemerisSource::Jpl => jpl_file.as_ref().map(|f| f.header().denum),
                EphemerisSource::Moshier => None,
            };
            if let Some(denum) = denum {
                config.tidal_acceleration = Some(crate::deltat::denum_to_tid_acc(denum));
            }
        }
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
        self.calc_with_config(jd_tt, body, flags, &self.config)
    }

    /// Same as [`calc`](Self::calc) but with an explicit config override. Used by the rise/set
    /// module (`riseset.rs`) to thread a caller-supplied `geopos` into the TOPOCTR pipeline
    /// without requiring it to match the `Ephemeris`'s own configured topographic position
    /// (mirrors C's per-call `swe_set_topo` before `swe_calc_ut`, but stateless).
    pub(crate) fn calc_with_config(
        &self,
        jd_tt: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<CalcResult, Error> {
        let flags = crate::calc::plaus_iflag(flags, config.ephemeris_source);
        if flags.contains(CalcFlags::TOPOCTR) && config.topographic.is_none() {
            return Err(Error::CError(
                "topocentric requires topographic position".to_string(),
            ));
        }

        if body == Body::Earth {
            return Ok(CalcResult {
                data: [0.0; 6],
                flags_used: flags,
            });
        }

        if flags.contains(CalcFlags::SPEED3) {
            return self.calc_speed3(jd_tt, body, flags, config);
        }

        let (mut xreturn, x2000, flags_used) = self.calc_inner(jd_tt, body, flags, config)?;
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

    /// Same as [`calc_ut`](Self::calc_ut) but with an explicit config override; see
    /// [`calc_with_config`](Self::calc_with_config).
    pub(crate) fn calc_ut_with_config(
        &self,
        jd_ut: f64,
        body: Body,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<CalcResult, Error> {
        let dt = crate::deltat::calc_deltat(jd_ut, config);
        self.calc_with_config(jd_ut + dt, body, flags, config)
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

    /// Ecliptic/equatorial -> azimuth + true/apparent altitude at `tjd_ut` (UT). Port of
    /// `swe_azalt` (swecl.c:2788-2825). `geopos` = [longitude, latitude, height above sea (m)].
    /// `atpress` in hPa (`0` => standard-atmosphere estimate from `geopos[2]`), `attemp` in deg
    /// C. `lapse_rate` in deg K/m -- C's default is `0.0065`; the stateless port has no
    /// `swe_set_lapse_rate` equivalent, so callers pass it explicitly (see
    /// docs/c-ref-refraction-azalt.md §9). `xin` = [lon/RA, lat/dec], degrees. Returns
    /// `[azimuth (from south, positive clockwise via west), true altitude, apparent altitude]`,
    /// degrees.
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

    /// Azimuth + true altitude -> ecliptic/equatorial coordinates at `tjd_ut` (UT). Port of
    /// `swe_azalt_rev` (swecl.c:2839-2873). Inverse of [`Ephemeris::azalt`]'s geometric
    /// transform only -- does NOT de-refract; `xin[1]` must already be a true altitude. `geopos`
    /// = [longitude, latitude, height (unused)]. `xin` = [azimuth (from south, clockwise), true
    /// altitude], degrees. Returns [lon/RA, lat/dec], degrees.
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

    /// Rise/set/meridian-transit search (full algorithm). Port of `swe_rise_trans_true_hor`
    /// (swecl.c:4387-4686); dispatches to `calc_mer_trans` when `rsmi` requests
    /// `MTRANSIT`/`ITRANSIT`. `starname` selects a fixed star (ignoring `body`); `horhgt` is the
    /// local horizon height above/below the sea-level horizon, degrees (`-100` = auto dip from
    /// `geopos[2]`). See docs/c-ref-riseset.md. The fast-path optimization and the
    /// `swe_rise_trans` dispatcher are a separate module (RSE 4).
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
    /// Port of `swe_rise_trans` (swecl.c:4355-4383, docs/c-ref-riseset.md §4). Fast-path
    /// eligible iff: not a fixed star, requesting RISE/SET (not a transit), `FORCE_SLOW` not
    /// set, no twilight bit, `body` in `Sun..=TrueNode`, and `|geopos[1]| <= 60` (`<= 65` for
    /// the Sun). Otherwise delegates to [`Ephemeris::rise_trans_true_hor`] with `horhgt = 0.0`.
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

    /// Solar eclipse shadow geometry: geographic position of greatest eclipse + core/penumbra
    /// shadow diameters, geocentric. Port of `swe_sol_eclipse_where`'s shadow-geometry pass
    /// (`eclipse_where`, swecl.c:565-582, 640-886); local-circumstance attributes (`attr[]`) come
    /// from [`Ephemeris::sol_eclipse_how`] instead.
    pub fn sol_eclipse_where(
        &self,
        tjd_ut: f64,
        ifl: CalcFlags,
    ) -> Result<crate::eclipse::EclipseWhere, Error> {
        crate::eclipse::sol_eclipse_where(self, tjd_ut, ifl)
    }

    /// Local circumstances of a solar eclipse at a specific observer: magnitude, obscuration,
    /// contact geometry, azimuth/altitude. Port of `swe_sol_eclipse_how` (swecl.c:922-964).
    /// `geopos` = [longitude, latitude, height above sea (m)], degrees/degrees/meters.
    pub fn sol_eclipse_how(
        &self,
        tjd_ut: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
    ) -> Result<crate::eclipse::EclipseHow, Error> {
        crate::eclipse::sol_eclipse_how(self, tjd_ut, ifl, geopos)
    }

    /// Planetary phenomena (phase angle, illuminated fraction, elongation, apparent diameter,
    /// apparent magnitude, Moon horizontal parallax) at `tjd_et` (ET). Port of `swe_pheno`
    /// (swecl.c:3802-4123). Returns the [`Phenomena`](crate::phenomena::Phenomena) plus the flags
    /// actually used.
    pub fn pheno(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<(crate::phenomena::Phenomena, CalcFlags), Error> {
        crate::phenomena::pheno(self, tjd_et, body, flags)
    }

    /// UT-based [`pheno`](Self::pheno). Port of `swe_pheno_ut` (swecl.c:4125-4142).
    pub fn pheno_ut(
        &self,
        tjd_ut: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<(crate::phenomena::Phenomena, CalcFlags), Error> {
        crate::phenomena::pheno_ut(self, tjd_ut, body, flags)
    }

    /// Nodes & apsides of `body` at `tjd_et` (TT). Port of `swe_nod_aps`
    /// (swecl.c:5075-5654). `method` selects mean vs osculating elements
    /// ([`NodApsMethod`](crate::NodApsMethod)); the mean branch (Sun..Neptune,
    /// Earth, Moon) is implemented, osculating is not yet (PNOC 5).
    pub fn nod_aps(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
        method: crate::nodaps::NodApsMethod,
    ) -> Result<crate::nodaps::NodesApsides, Error> {
        if flags.contains(CalcFlags::TOPOCTR) && self.config.topographic.is_none() {
            return Err(Error::CError(
                "topocentric requires topographic position".to_string(),
            ));
        }
        crate::nodaps::nod_aps(self, tjd_et, body, flags, method)
    }

    /// UT-based [`nod_aps`](Self::nod_aps). Port of `swe_nod_aps_ut`
    /// (swecl.c:5656-5665): converts UT→TT via deltaT, then delegates.
    pub fn nod_aps_ut(
        &self,
        tjd_ut: f64,
        body: Body,
        flags: CalcFlags,
        method: crate::nodaps::NodApsMethod,
    ) -> Result<crate::nodaps::NodesApsides, Error> {
        let tjde = tjd_ut + crate::deltat::calc_deltat(tjd_ut, &self.config);
        self.nod_aps(tjde, body, flags, method)
    }

    /// Osculating (Keplerian) orbital elements of `body` at `tjd_et` (TT). Port
    /// of `swe_get_orbital_elements` (swecl.c:5783-5971). Rejects the Sun, the
    /// lunar nodes, and the apsides. Note: `SEFLG_TOPOCTR` is bit-aliased onto
    /// `SEFLG_ORBEL_AA` here ("sum masses inside the orbit"), NOT a topocentric
    /// request — see [`crate::orbit`].
    pub fn get_orbital_elements(
        &self,
        tjd_et: f64,
        body: Body,
        flags: CalcFlags,
    ) -> Result<crate::orbit::OrbitalElements, Error> {
        crate::orbit::get_orbital_elements(self, tjd_et, body, flags)
    }

    /// Maximum, minimum, and current true distance of `body` (AU) at `tjd_et`
    /// (TT), returned as `(dmax, dmin, dtrue)`. Port of
    /// `swe_orbit_max_min_true_distance` (swecl.c:6170-6287).
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

    /// Next JD (ET) at which the Sun's ecliptic longitude equals `x2cross`.
    pub fn solcross(&self, x2cross: f64, jd_et: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::solcross(self, x2cross, jd_et, flags)
    }

    /// UT-based [`solcross`](Self::solcross).
    pub fn solcross_ut(&self, x2cross: f64, jd_ut: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::solcross_ut(self, x2cross, jd_ut, flags)
    }

    /// Next JD (ET) at which the Moon's ecliptic longitude equals `x2cross`.
    pub fn mooncross(&self, x2cross: f64, jd_et: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::mooncross(self, x2cross, jd_et, flags)
    }

    /// UT-based [`mooncross`](Self::mooncross).
    pub fn mooncross_ut(&self, x2cross: f64, jd_ut: f64, flags: CalcFlags) -> Result<f64, Error> {
        crate::crossings::mooncross_ut(self, x2cross, jd_ut, flags)
    }

    /// Next JD (ET) at which the Moon crosses its node (latitude = 0).
    pub fn mooncross_node(
        &self,
        jd_et: f64,
        flags: CalcFlags,
    ) -> Result<crate::crossings::MoonCrossing, Error> {
        crate::crossings::mooncross_node(self, jd_et, flags)
    }

    /// UT-based [`mooncross_node`](Self::mooncross_node).
    pub fn mooncross_node_ut(
        &self,
        jd_ut: f64,
        flags: CalcFlags,
    ) -> Result<crate::crossings::MoonCrossing, Error> {
        crate::crossings::mooncross_node_ut(self, jd_ut, flags)
    }

    /// Next JD (ET) at which `body`'s heliocentric longitude equals `x2cross`.
    /// `dir >= 0` searches forward, `dir < 0` searches backward.
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
        let config = &self.config;
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
            EphemerisSource::Swiss => {
                let provider = SwephProvider {
                    planet_files: &self.planet_files,
                    moon_files: &self.moon_files,
                };
                let pos = provider.positions(Body::Sun, t, true)?;
                Ok(crate::nodaps::ObsFrame {
                    sun_bary: pos.sun_bary,
                    xear: pos.earth_bary,
                    topo,
                })
            }
            EphemerisSource::Jpl => {
                let provider = JplProvider {
                    file: self.jpl_file.as_ref().unwrap(),
                };
                let pos = provider.positions(Body::Sun, t, true)?;
                Ok(crate::nodaps::ObsFrame {
                    sun_bary: pos.sun_bary,
                    xear: pos.earth_bary,
                    topo,
                })
            }
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
    ) -> Result<[f64; 6], Error> {
        let source = self.config.ephemeris_source;
        let models = &self.config.astro_models;
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
            EphemerisSource::Swiss => {
                let provider = SwephProvider {
                    planet_files: &self.planet_files,
                    moon_files: &self.moon_files,
                };
                let query = if ipli == Body::Earth { Body::Sun } else { ipli };
                let pos = provider.positions(query, t, true)?;
                nodaps_osc_frame(&pos, ipli, want_bary)
            }
            EphemerisSource::Jpl => {
                let provider = JplProvider {
                    file: self.jpl_file.as_ref().unwrap(),
                };
                let query = if ipli == Body::Earth { Body::Sun } else { ipli };
                let pos = provider.positions(query, t, true)?;
                nodaps_osc_frame(&pos, ipli, want_bary)
            }
        };

        if ipli == Body::Earth {
            let moon = self.raw_moon_at(source, t, &eps_j2000)?;
            for i in 0..6 {
                xx[i] += moon[i] / (crate::constants::EARTH_MOON_MRAT + 1.0);
            }
        }

        Ok(xx)
    }

    /// Global eclipse search: next/previous solar eclipse anywhere on Earth from `tjd_start`
    /// (UT), restricted to eclipse types in `ifltype` (empty = all types). Port of
    /// `swe_sol_eclipse_when_glob` (swecl.c:1185-1515).
    pub fn sol_eclipse_when_glob(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        ifltype: EclipseFlags,
        backward: bool,
    ) -> Result<crate::eclipse::SolarEclipseGlobal, Error> {
        crate::eclipse::sol_eclipse_when_glob(self, tjd_start, ifl, ifltype, backward)
    }

    /// Local eclipse search: next/previous solar eclipse *visible from* `geopos` (topocentric),
    /// with local contact times C1-C4 and full local circumstances. Port of
    /// `swe_sol_eclipse_when_loc` (swecl.c:2019-2041, 2100-2410). `geopos` = [longitude,
    /// latitude, height above sea (m)], degrees/degrees/meters.
    pub fn sol_eclipse_when_loc(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
        backward: bool,
    ) -> Result<crate::eclipse::SolarEclipseLocal, Error> {
        crate::eclipse::sol_eclipse_when_loc(self, tjd_start, ifl, geopos, backward)
    }

    /// Local circumstances of a lunar eclipse at `geopos`: umbral/penumbral magnitude, Saros
    /// series/member, and the Moon's azimuth/true/apparent altitude. Port of
    /// `swe_lun_eclipse_how` (swecl.c:3190-3239). `geopos` = [longitude, latitude, height above
    /// sea (m)], degrees/degrees/meters.
    pub fn lun_eclipse_how(
        &self,
        tjd_ut: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
    ) -> Result<crate::eclipse::LunarEclipseHow, Error> {
        crate::eclipse::swe_lun_eclipse_how(self, tjd_ut, ifl, geopos)
    }

    /// Global lunar-eclipse search: next/previous lunar eclipse from `tjd_start` (UT), restricted
    /// to eclipse types in `ifltype` (empty = any of TOTAL/PARTIAL/PENUMBRAL). Purely geocentric,
    /// no geographic position. Port of `swe_lun_eclipse_when` (swecl.c:3389-3616).
    pub fn lun_eclipse_when(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        ifltype: EclipseFlags,
        backward: bool,
    ) -> Result<crate::eclipse::LunarEclipseGlobal, Error> {
        crate::eclipse::lun_eclipse_when(self, tjd_start, ifl, ifltype, backward)
    }

    /// Local lunar-eclipse search: next/previous lunar eclipse visible from `geopos` (Moon above
    /// the horizon during some phase), with contact times clipped to moonrise/moonset. Port of
    /// `swe_lun_eclipse_when_loc` (swecl.c:3644-3739). `geopos` = [longitude, latitude, height
    /// above sea (m)], degrees/degrees/meters.
    pub fn lun_eclipse_when_loc(
        &self,
        tjd_start: f64,
        ifl: CalcFlags,
        geopos: [f64; 3],
        backward: bool,
    ) -> Result<crate::eclipse::LunarEclipseLocal, Error> {
        crate::eclipse::lun_eclipse_when_loc(self, tjd_start, ifl, geopos, backward)
    }

    /// Geographic position of maximal occultation of `body`/`starname` by the Moon at `tjd_ut`
    /// (UT). `starname` (if given, non-empty) takes precedence over `body`. Port of
    /// `swe_lun_occult_where` (swecl.c:606-630).
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
    /// anywhere on Earth from `tjd_start` (UT), restricted to types in `ifltype` (empty = all
    /// types valid for the occulted body). `starname` (if given, non-empty) takes precedence
    /// over `body`. Port of `swe_lun_occult_when_glob` (swecl.c:1572-1984).
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

    /// Local occultation search: next/previous occultation of `body`/`starname` by the Moon
    /// *visible from* `geopos` (topocentric), with local contact times and circumstances.
    /// `starname` (if given, non-empty) takes precedence over `body`. Port of
    /// `swe_lun_occult_when_loc` (swecl.c:2071-2098, 2412-2764).
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

    /// Gauquelin sector position of a body, geometric method (`imeth` 0 = with ecliptic
    /// latitude, 1 = without). Port of `swe_gauquelin_sector`'s `imeth ∈ {0,1}` branch
    /// (swecl.c:6338-6356) — reuses `swe_house_pos`'s `'G'` branch directly. `imeth ∈ {2,3,4,5}`
    /// Gauquelin sector position via geometric house position (imeth 0/1).
    /// See docs/c-ref-houses.md §10, docs/c-ref-gauquelin-riseset.md §9.
    ///
    /// When `starname` is `Some(non_empty)`, uses `fixstar2` instead of `calc` to resolve
    /// the body position (swecl.c:6356-6362). Unlike [`Ephemeris::houses_ex2`], the
    /// deltaT/obliquity/nutation resolution here uses the caller's `flags` directly, not a
    /// forced `TIDAL_DEFAULT` override.
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

        let models = &self.config.astro_models;
        let t_et = t_ut + crate::deltat::calc_deltat(t_ut, &self.config);
        let eps = crate::obliquity::obliquity(t_et, flags, models).eps * RADTODEG;
        let nut = crate::nutation::nutation(t_et, flags, models);
        let dpsi_deg = nut.dpsi * RADTODEG;
        let deps_deg = nut.deps * RADTODEG;
        let eps_true = eps + deps_deg;
        let armc = crate::math::normalize_degrees(
            crate::sidereal_time::sidereal_time0(t_ut, eps_true, dpsi_deg, &self.config) * 15.0
                + geolon,
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

    /// Full Gauquelin sector dispatcher. Routes imeth 0/1 to the geometric method
    /// and imeth 2–5 to the rise/set method. Port of `swe_gauquelin_sector`
    /// (swecl.c:6309-6439, docs/c-ref-gauquelin-riseset.md).
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

        // Heliocentric (SEFLG_HELCTR) is supported below (per-backend/-body branches in calc.rs);
        // plaus_iflag has already forced NOABERR|NOGDEFL for it. Barycentric is still unported.
        if flags.contains(CalcFlags::BARYCTR) {
            return Err(Error::UnsupportedFlags(flags & CalcFlags::BARYCTR));
        }

        // Heliocentric Sun is the origin (the Sun relative to itself) — C's swe_calc returns an
        // all-zero xx (position and speed) across every output frame, since the zero vector is
        // invariant under bias/precession/nutation and polar-converts to zeros. calc_sun has no
        // HELCTR branch (it always builds the geocentric -observer vector), so short-circuit here,
        // matching the MeanNode/MeanApogee HELCTR handling above.
        if body == Body::Sun && flags.contains(CalcFlags::HELCTR) {
            return Ok(([0.0; 24], [0.0; 6], flags));
        }

        let eps_j2000 =
            crate::obliquity::obliquity(crate::constants::J2000, CalcFlags::empty(), models);

        match config.ephemeris_source {
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
            EphemerisSource::Jpl => {
                let (xr, x2000) =
                    self.calc_body_jpl(jd_tt, body, &eps_j2000, flags, models, config)?;
                Ok((xr, x2000, flags))
            }
            EphemerisSource::Moshier => {
                let (xr, x2000) =
                    self.calc_body_moshier(jd_tt, body, &eps_j2000, flags, models, config)?;
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
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6]), Error> {
        match body {
            Body::Sun => crate::calc::calc_sun(jd_tt, eps_j2000, flags, config, models),
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
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6]), Error> {
        match body {
            Body::Sun => crate::calc::calc_sun_sweph(
                jd_tt,
                &self.planet_files,
                &self.moon_files,
                flags,
                config,
                models,
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
        config: &EphemerisConfig,
    ) -> Result<([f64; 24], [f64; 6]), Error> {
        let file = self.jpl_file.as_ref().unwrap();
        match body {
            Body::Sun => crate::calc::calc_sun_jpl(jd_tt, file, flags, config, models),
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
            EphemerisSource::Swiss => crate::calc::raw_osc_moon_sweph(&self.moon_files, t),
            EphemerisSource::Jpl => {
                crate::calc::raw_osc_moon_jpl(self.jpl_file.as_ref().unwrap(), t)
            }
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

    pub(crate) fn apply_sidereal(
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
        self.fixstar2_with_config(star, jd_tt, flags, &self.config)
    }

    /// Same as [`fixstar2`](Self::fixstar2) but with an explicit config override -- see
    /// [`calc_with_config`](Self::calc_with_config). Lets callers (e.g. `eclipse.rs`'s
    /// `eclipse_how`/`occult_when_loc`) get a topocentric fixed-star position at a caller-supplied
    /// `geopos` without requiring it to match the `Ephemeris`'s own configured topographic
    /// position.
    pub(crate) fn fixstar2_with_config(
        &self,
        star: &str,
        jd_tt: f64,
        flags: CalcFlags,
        config: &EphemerisConfig,
    ) -> Result<(String, CalcResult), Error> {
        // C's swe_fixstar2 returns the original input iflag unchanged (it passes
        // iflag by value to fixstar_calc_from_struct and ignores the return).
        let orig_flags = flags;
        let flags = crate::calc::plaus_iflag(flags, config.ephemeris_source);
        let resolved = if let Some(s) = crate::stars::builtin_star(star) {
            s
        } else {
            self.stars.search(star)?
        };
        let data = self.calc_fixstar(&resolved, jd_tt, flags, config)?;
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
            crate::types::EphemerisSource::Swiss => {
                self.calc_fixstar_sweph(star, jd, flags, config)
            }
            crate::types::EphemerisSource::Jpl => self.calc_fixstar_jpl(star, jd, flags, config),
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
        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary)
    }

    /// SWIEPH backend: barycentric Earth for parallax/aberration, sun_bary for deflection.
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

        self.calc_fixstar_inner(star, jd, flags, xobs, xobs_dt, sun_bary)
    }

    /// JPL backend: barycentric Earth for parallax/aberration, sun_bary for deflection.
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
