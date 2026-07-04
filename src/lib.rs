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
//! `Default`). Per-call overrides for ephemeris source, topographic position,
//! and sidereal mode are available through flag bits and the `*_with_config`
//! internal variants.
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
//! | [`moshier`] | Moshier analytical backend (always available) |
//! | [`precession`], [`nutation`], [`obliquity`], [`bias`], [`sidereal_time`], [`deltat`] | Low-level positional astronomy — prefer `Ephemeris` methods |
//! | [`ayanamsa`] | Sidereal ayanamsa computation |
//! | [`corrections`] | Aberration, light deflection, relativistic mass effect |
//! | [`fictitious`] | Fictitious/hypothetical planet elements |
//! | [`topocentric`] | Observer geocentric offset |

pub mod config;
pub mod constants;
pub mod context;
pub mod error;
pub mod flags;
pub mod types;

pub mod ayanamsa;
pub mod azalt;
pub mod bias;
pub mod calc;
pub mod corrections;
pub mod crossings;
pub mod date;
pub mod deltat;
pub mod eclipse;
pub mod fictitious;
pub mod format;
pub mod heliacal;
pub mod houses;
#[cfg(feature = "jpl")]
pub mod jpl;
pub mod math;
pub mod moshier;
pub mod nodaps;
pub mod nutation;
pub mod obliquity;
pub mod orbit;
pub mod phenomena;
pub mod precession;
pub mod riseset;
pub mod sidereal_time;
pub mod stars;
#[cfg(feature = "swisseph-files")]
pub mod sweph_file;
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
    EphemerisSource, Epsilon, FictitiousBody, FictitiousId, FrameTransform, HouseSystem, JdTt,
    JdUt1, JplHorMode, JplHoraMode, Nutation, NutationModel, PlanetMoonId, PrecessionDirection,
    PrecessionModel, SiderealMode, SiderealTimeModel, UtcComponents, UtcToJd,
};

pub type Result<T> = std::result::Result<T, Error>;
