use std::ffi::c_char;
use std::path::PathBuf;

use swisseph::config::{EphemerisConfig, TopoPosition};
use swisseph::types::EphemerisSource;

/// Flattened C-ABI configuration struct. Mirrors `EphemerisConfig` with C-compatible types.
///
/// Use `swisseph_config_default` to initialize — a zeroed struct has WRONG defaults
/// (NAN sentinels for tidal_acceleration/delta_t_userdef must be set explicitly).
#[repr(C)]
pub struct SweConfig {
    /// 0 = Moshier, 1 = Swiss, 2 = Jpl.
    pub ephemeris_source: i32,
    /// Path to .se1 files. NULL = none.
    pub ephe_path: *const c_char,
    /// JPL filename (e.g. "de441.eph"). NULL = default.
    pub jpl_filename: *const c_char,
    /// Path to a custom leap-seconds file. NULL = use built-in table.
    pub leap_seconds_file: *const c_char,

    // -- Sidereal --
    /// If true, sidereal mode is active.
    pub has_sidereal: bool,
    /// Raw swe_set_sid_mode value (bits 0-7 = mode index, upper bits = SiderealBits).
    pub sid_mode: i32,
    /// Reference epoch for user-defined sidereal (sid_mode & 0xFF == 255).
    pub sid_t0: f64,
    /// Initial ayanamsa at sid_t0.
    pub sid_ayan_t0: f64,

    // -- Topographic --
    /// If true, topographic position is set.
    pub has_topo: bool,
    /// Geographic longitude, degrees east-positive.
    pub geolon: f64,
    /// Geographic latitude, degrees north-positive.
    pub geolat: f64,
    /// Altitude above sea level, meters.
    pub altitude: f64,

    // -- Scalars with NAN = unset --
    /// Tidal acceleration override (arcsec/century^2). NAN = auto-derive from ephemeris.
    pub tidal_acceleration: f64,
    /// User-defined Delta T (days). NAN = use models.
    pub delta_t_userdef: f64,

    // -- Variable-length arrays (null + 0 = empty) --
    /// Pointer to asteroid MPC numbers. NULL + len=0 means empty.
    pub asteroid_numbers: *const i32,
    /// Length of asteroid_numbers array.
    pub asteroid_numbers_len: usize,
    /// Pointer to planet-moon ids. NULL + len=0 means empty.
    pub planet_moon_numbers: *const i32,
    /// Length of planet_moon_numbers array.
    pub planet_moon_numbers_len: usize,
    /// Pointer to extra leap-second years. NULL + len=0 means empty.
    pub extra_leap_seconds: *const i32,
    /// Length of extra_leap_seconds array.
    pub extra_leap_seconds_len: usize,

    // -- Astro models (i32 each, 0 = library default) --
    /// Precession long-term model. 0 = default.
    pub astro_model_prec_longterm: i32,
    /// Precession short-term model. 0 = default.
    pub astro_model_prec_shortterm: i32,
    /// Nutation model. 0 = default.
    pub astro_model_nutation: i32,
    /// Frame bias model. 0 = default.
    pub astro_model_bias: i32,
    /// JPL Horizons mode. 0 = default.
    pub astro_model_jplhor: i32,
    /// JPL Horizons approx mode. 0 = default.
    pub astro_model_jplhora: i32,
    /// Sidereal time model. 0 = default.
    pub astro_model_sidereal_time: i32,
    /// Delta T model. 0 = default.
    pub astro_model_delta_t: i32,
}

/// Fill a SweConfig with sane defaults (Moshier, no paths, NAN for unset scalars).
///
/// # Safety
/// `config` must point to a valid, writable `SweConfig`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swisseph_config_default(config: *mut SweConfig) {
    if config.is_null() {
        return;
    }
    unsafe {
        (*config) = SweConfig {
            ephemeris_source: 0,
            ephe_path: std::ptr::null(),
            jpl_filename: std::ptr::null(),
            leap_seconds_file: std::ptr::null(),
            has_sidereal: false,
            sid_mode: 0,
            sid_t0: 0.0,
            sid_ayan_t0: 0.0,
            has_topo: false,
            geolon: 0.0,
            geolat: 0.0,
            altitude: 0.0,
            tidal_acceleration: f64::NAN,
            delta_t_userdef: f64::NAN,
            asteroid_numbers: std::ptr::null(),
            asteroid_numbers_len: 0,
            planet_moon_numbers: std::ptr::null(),
            planet_moon_numbers_len: 0,
            extra_leap_seconds: std::ptr::null(),
            extra_leap_seconds_len: 0,
            astro_model_prec_longterm: 0,
            astro_model_prec_shortterm: 0,
            astro_model_nutation: 0,
            astro_model_bias: 0,
            astro_model_jplhor: 0,
            astro_model_jplhora: 0,
            astro_model_sidereal_time: 0,
            astro_model_delta_t: 0,
        };
    }
}

/// Convert a C SweConfig into a Rust EphemerisConfig.
/// Returns Err(&str) on invalid arguments (null pointers where non-null expected, bad UTF-8, etc.)
pub(crate) unsafe fn config_to_rust(c: &SweConfig) -> Result<EphemerisConfig, &'static str> {
    let ephemeris_source = match c.ephemeris_source {
        0 => EphemerisSource::Moshier,
        1 => EphemerisSource::Swiss,
        2 => EphemerisSource::Jpl,
        _ => return Err("invalid ephemeris_source (must be 0, 1, or 2)"),
    };

    let ephe_path = unsafe { nullable_cstr_to_pathbuf(c.ephe_path)? };
    let jpl_filename = unsafe { nullable_cstr_to_string(c.jpl_filename)? };
    let leap_seconds_file = unsafe { nullable_cstr_to_pathbuf(c.leap_seconds_file)? };

    let mut config = EphemerisConfig {
        ephemeris_source,
        ephe_path,
        jpl_filename,
        leap_seconds_file,
        tidal_acceleration: if c.tidal_acceleration.is_nan() {
            None
        } else {
            Some(c.tidal_acceleration)
        },
        delta_t_userdef: if c.delta_t_userdef.is_nan() {
            None
        } else {
            Some(c.delta_t_userdef)
        },
        ..EphemerisConfig::default()
    };

    if c.has_topo {
        config.topographic = Some(TopoPosition {
            longitude: c.geolon,
            latitude: c.geolat,
            altitude: c.altitude,
        });
    }

    if c.has_sidereal {
        config.set_sidereal_mode(c.sid_mode, c.sid_t0, c.sid_ayan_t0);
    }

    // Variable-length arrays
    config.asteroid_numbers =
        unsafe { slice_from_ptr(c.asteroid_numbers, c.asteroid_numbers_len)? };
    config.planet_moon_numbers =
        unsafe { slice_from_ptr(c.planet_moon_numbers, c.planet_moon_numbers_len)? };
    config.extra_leap_seconds =
        unsafe { slice_from_ptr(c.extra_leap_seconds, c.extra_leap_seconds_len)? };

    // Astro models — 0 means "default", non-zero is the raw C enum value
    apply_astro_models(&mut config, c)?;

    Ok(config)
}

fn apply_astro_models(config: &mut EphemerisConfig, c: &SweConfig) -> Result<(), &'static str> {
    if c.astro_model_prec_longterm != 0 {
        config.astro_models.prec_longterm = prec_model_from_i32(c.astro_model_prec_longterm)
            .ok_or("invalid astro_model_prec_longterm")?;
    }
    if c.astro_model_prec_shortterm != 0 {
        config.astro_models.prec_shortterm = prec_model_from_i32(c.astro_model_prec_shortterm)
            .ok_or("invalid astro_model_prec_shortterm")?;
    }
    if c.astro_model_nutation != 0 {
        config.astro_models.nutation = nutation_model_from_i32(c.astro_model_nutation)
            .ok_or("invalid astro_model_nutation")?;
    }
    if c.astro_model_bias != 0 {
        config.astro_models.bias =
            bias_model_from_i32(c.astro_model_bias).ok_or("invalid astro_model_bias")?;
    }
    if c.astro_model_jplhor != 0 {
        config.astro_models.jplhor_mode =
            jplhor_mode_from_i32(c.astro_model_jplhor).ok_or("invalid astro_model_jplhor")?;
    }
    if c.astro_model_jplhora != 0 {
        config.astro_models.jplhora_mode =
            jplhora_mode_from_i32(c.astro_model_jplhora).ok_or("invalid astro_model_jplhora")?;
    }
    if c.astro_model_sidereal_time != 0 {
        config.astro_models.sidereal_time =
            sidereal_time_model_from_i32(c.astro_model_sidereal_time)
                .ok_or("invalid astro_model_sidereal_time")?;
    }
    if c.astro_model_delta_t != 0 {
        config.astro_models.delta_t =
            delta_t_model_from_i32(c.astro_model_delta_t).ok_or("invalid astro_model_delta_t")?;
    }
    Ok(())
}

fn prec_model_from_i32(v: i32) -> Option<swisseph::PrecessionModel> {
    use swisseph::PrecessionModel;
    match v {
        1 => Some(PrecessionModel::IAU1976),
        2 => Some(PrecessionModel::Laskar1986),
        3 => Some(PrecessionModel::WillEpsLask),
        4 => Some(PrecessionModel::Williams1994),
        5 => Some(PrecessionModel::Simon1994),
        6 => Some(PrecessionModel::IAU2000),
        7 => Some(PrecessionModel::Bretagnon2003),
        8 => Some(PrecessionModel::IAU2006),
        9 => Some(PrecessionModel::Vondrak2011),
        10 => Some(PrecessionModel::Owen1990),
        11 => Some(PrecessionModel::Newcomb),
        _ => None,
    }
}

fn nutation_model_from_i32(v: i32) -> Option<swisseph::NutationModel> {
    use swisseph::NutationModel;
    match v {
        1 => Some(NutationModel::IAU1980),
        2 => Some(NutationModel::IAUCorr1987),
        3 => Some(NutationModel::IAU2000A),
        4 => Some(NutationModel::IAU2000B),
        5 => Some(NutationModel::Woolard),
        _ => None,
    }
}

fn bias_model_from_i32(v: i32) -> Option<swisseph::BiasModel> {
    use swisseph::BiasModel;
    match v {
        1 => Some(BiasModel::None),
        2 => Some(BiasModel::IAU2000),
        3 => Some(BiasModel::IAU2006),
        _ => std::option::Option::None,
    }
}

fn jplhor_mode_from_i32(v: i32) -> Option<swisseph::JplHorMode> {
    use swisseph::JplHorMode;
    match v {
        1 => Some(JplHorMode::LongAgreement),
        _ => None,
    }
}

fn jplhora_mode_from_i32(v: i32) -> Option<swisseph::JplHoraMode> {
    use swisseph::JplHoraMode;
    match v {
        1 => Some(JplHoraMode::V1),
        2 => Some(JplHoraMode::V2),
        3 => Some(JplHoraMode::V3),
        _ => None,
    }
}

fn sidereal_time_model_from_i32(v: i32) -> Option<swisseph::SiderealTimeModel> {
    use swisseph::SiderealTimeModel;
    match v {
        1 => Some(SiderealTimeModel::IAU1976),
        2 => Some(SiderealTimeModel::IAU2006),
        3 => Some(SiderealTimeModel::IersConv2010),
        4 => Some(SiderealTimeModel::Longterm),
        _ => None,
    }
}

fn delta_t_model_from_i32(v: i32) -> Option<swisseph::DeltaTModel> {
    use swisseph::DeltaTModel;
    match v {
        1 => Some(DeltaTModel::StephensonMorrison1984),
        2 => Some(DeltaTModel::Stephenson1997),
        3 => Some(DeltaTModel::StephensonMorrison2004),
        4 => Some(DeltaTModel::EspenakMeeus2006),
        5 => Some(DeltaTModel::StephensonEtc2016),
        _ => None,
    }
}

unsafe fn nullable_cstr_to_pathbuf(ptr: *const c_char) -> Result<Option<PathBuf>, &'static str> {
    if ptr.is_null() {
        return Ok(None);
    }
    let cstr = unsafe { std::ffi::CStr::from_ptr(ptr) };
    let s = cstr.to_str().map_err(|_| "invalid UTF-8 in path")?;
    Ok(Some(PathBuf::from(s)))
}

unsafe fn nullable_cstr_to_string(ptr: *const c_char) -> Result<Option<String>, &'static str> {
    if ptr.is_null() {
        return Ok(None);
    }
    let cstr = unsafe { std::ffi::CStr::from_ptr(ptr) };
    let s = cstr.to_str().map_err(|_| "invalid UTF-8 in string")?;
    Ok(Some(s.to_owned()))
}

unsafe fn slice_from_ptr(ptr: *const i32, len: usize) -> Result<Vec<i32>, &'static str> {
    if len == 0 {
        return Ok(Vec::new());
    }
    if ptr.is_null() {
        return Err("null pointer with nonzero length");
    }
    Ok(unsafe { std::slice::from_raw_parts(ptr, len).to_vec() })
}
