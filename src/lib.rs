#![warn(missing_docs)]
//! Pure-Rust, stateless port of the [Swiss Ephemeris](https://www.astro.com/swisseph/)
//! astronomical calculation library (C version 2.10.03).
//!
//! # Quick start
//!
//! The Moshier backend is self-contained — no data files needed:
//!
//! ```
//! use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig};
//!
//! let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();
//! let jd_ut = 2451545.0; // J2000.0
//! let result = eph.calc_ut(jd_ut, Body::Sun, CalcFlags::SPEED).unwrap();
//! let longitude = result.data[0];
//! let latitude  = result.data[1];
//! let distance  = result.data[2];
//! assert!((280.0..=281.0).contains(&longitude));
//! ```
//!
//! # Stateless design
//!
//! [`Ephemeris`] holds only read-only configuration — no mutable cache, no
//! internal state. All methods take `&self`, never `&mut self`. The calculation
//! pipeline is pure: inputs → math → output, no side effects. `Ephemeris` is
//! `Send + Sync`.
//!
//! Configuration is set once via [`EphemerisConfig`] (which implements
//! `Default`). The ephemeris source for a single call can be overridden via
//! [`CalcFlags`] flag bits (`MOSEPH`, `SWIEPH`, `JPLEPH`); topographic
//! position, sidereal mode, and other settings are fixed at construction.
//!
//! # Modules
//!
//! | Module | Role |
//! |--------|------|
//! | [`context`] | **Primary API** — [`Ephemeris`] and [`CalcResult`] |
//! | [`config`] | [`EphemerisConfig`] and [`TopoPosition`] |
//! | [`types`] | Domain types: [`Body`], [`HouseSystem`], [`SiderealMode`], etc. |
//! | [`flags`] | Bitflag structs: [`CalcFlags`], [`EclipseFlags`], [`RiseSetFlags`], etc. |
//! | [`error`] | [`Error`] enum |
//! | [`houses`] | House systems and cusps |
//! | [`eclipse`] | Solar/lunar eclipses and occultations |
//! | [`heliacal`] | Heliacal visibility (rising/setting phenomena) |
//! | [`riseset`] | Rise, set, and meridian transit |
//! | [`azalt`] | Horizontal coordinates and atmospheric refraction |
//! | [`phenomena`] | Phase angle, elongation, magnitude |
//! | [`nodaps`] | Planetary nodes and apsides |
//! | [`orbit`] | Keplerian orbital elements |
//! | [`crossings`] | Longitude/latitude crossing search |
//! | [`stars`] | Fixed-star catalog |
//! | [`date`] | Julian Day ↔ calendar, UTC conversion |
//! | [`mod@format`] | Degree/time string formatting |
//! | [`math`] | Coordinate transforms, Chebyshev evaluation, degree splitting |
//! | [`constants`] | Physical constants, epochs, unit conversions |
//! | [`calc`] | Calculation pipeline internals (light-time, aberration, frame transforms) |
//! | [`moshier`] | Moshier analytical backend (always available) |
//! | [`sweph_file`] | Swiss Ephemeris `.se1` file reader (requires feature `swisseph-files`) |
//! | [`jpl`] | JPL DE ephemeris reader (requires feature `jpl`) |
//! | [`precession`], [`nutation`], [`obliquity`], [`bias`], [`sidereal_time`], [`deltat`] | Low-level positional astronomy — prefer `Ephemeris` methods |
//! | [`ayanamsa`] | Sidereal ayanamsa computation |
//! | [`corrections`] | Aberration, light deflection, relativistic mass effect |
//! | [`fictitious`] | Fictitious/hypothetical planet elements |
//! | [`topocentric`] | Observer geocentric offset |

/// Ephemeris construction and configuration.
pub mod config;
/// Physical constants, epochs, and unit conversions.
pub mod constants;
/// Primary API — [`Ephemeris`] and [`CalcResult`].
pub mod context;
/// Error types.
pub mod error;
/// Bitflag structs for calculation, eclipse, rise/set, and formatting options.
pub mod flags;
/// Domain types: [`Body`], [`HouseSystem`], [`SiderealMode`], newtypes, enums.
pub mod types;

/// Sidereal ayanamsa computation.
pub mod ayanamsa;
/// Horizontal coordinates and atmospheric refraction.
pub mod azalt;
/// GCRS ↔ J2000 frame rotation (IERS 2006 bias matrix).
pub mod bias;
/// Calculation pipeline internals: light-time, aberration, frame transforms.
pub mod calc;
/// Relativistic corrections: aberration, light deflection, mass effect.
pub mod corrections;
/// Longitude / latitude crossing search.
pub mod crossings;
/// Julian Day ↔ calendar conversion, delta-T, UTC.
pub mod date;
/// Delta-T models and historical tables.
pub mod deltat;
/// Solar/lunar eclipses and stellar/planetary occultations.
pub mod eclipse;
/// Fictitious / hypothetical planet elements and orbital mechanics.
pub mod fictitious;
/// Degree/time string formatting (centisecond precision).
pub mod format;
/// Heliacal visibility: atmospheric extinction, optics, event search.
pub mod heliacal;
/// House systems and cusps.
pub mod houses;
/// JPL Development Ephemeris (DE) file reader.
#[cfg(feature = "jpl")]
pub mod jpl;
/// Coordinate transforms, Chebyshev evaluation, degree arithmetic.
pub mod math;
/// Moshier analytical ephemeris backend (always available, no data files).
pub mod moshier;
/// Planetary nodes and apsides (mean and osculating).
pub mod nodaps;
/// Nutation models and term tables.
pub mod nutation;
/// Obliquity of the ecliptic (11 models).
pub mod obliquity;
/// Osculating (Keplerian) orbital elements and distance extrema.
pub mod orbit;
/// Phase angle, elongation, apparent magnitude.
pub mod phenomena;
/// Precession models (3 algorithm families, 11 models).
pub mod precession;
/// Rise, set, and meridian transit.
pub mod riseset;
/// Greenwich Mean Sidereal Time and Equation of the Equinoxes.
pub mod sidereal_time;
/// Fixed-star catalog loading and search.
pub mod stars;
/// Swiss Ephemeris `.se1` file reader.
#[cfg(feature = "swisseph-files")]
pub mod sweph_file;
/// Observer geocentric offset for topocentric calculations.
pub mod topocentric;

pub use azalt::{AzAltDir, HorDir, RefracDir};
pub use config::{EphemerisConfig, TopoPosition};
pub use context::{CalcResult, Ephemeris};
pub use crossings::MoonCrossing;
pub use eclipse::{
    EclipseHow, EclipseWhere, LunarEclipseGlobal, LunarEclipseHow, LunarEclipseLocal, OccultGlobal,
    OccultLocal, SolarEclipseGlobal, SolarEclipseLocal,
};
pub use error::Error;
pub use flags::{
    CalcFlags, EclipseFlags, HeliacalFlags, RiseSetFlags, SiderealBits, SplitDegFlags, VisLimFlags,
};
pub use heliacal::{
    HeliacalAngleResult, HeliacalEvent, HeliacalEventType, HeliacalPheno, VisLimitResult,
};
pub use houses::{AscMc, HouseResult};
pub use nodaps::{NodApsMethod, NodesApsides};
pub use orbit::OrbitalElements;
pub use phenomena::Phenomena;
pub use riseset::RiseSetResult;
pub use stars::{Star, StarCatalog};
pub use types::{
    AsteroidId, AstroModels, BiasModel, Body, CalendarType, DegreeParts, DeltaT, DeltaTModel,
    EphemerisSource, Epsilon, FictitiousBody, FictitiousId, FileData, FileDataKind, FrameTransform,
    HouseSystem, JdTt, JdUt1, JplHorMode, JplHoraMode, Nutation, NutationModel, PlanetMoonId,
    PrecessionDirection, PrecessionModel, SiderealMode, SiderealTimeModel, UtcComponents, UtcToJd,
};

/// Convenience alias for `Result<T, swisseph::Error>`.
pub type Result<T> = std::result::Result<T, Error>;
