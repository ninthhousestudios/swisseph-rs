pub mod constants;
pub mod context;
pub mod error;
pub mod flags;
pub mod types;

pub mod ayanamsa;
pub mod calc;
pub mod date;
pub mod eclipse;
pub mod heliacal;
pub mod houses;
pub mod jpl;
pub mod math;
pub mod moshier;
pub mod phenomena;
pub mod precession;
pub mod stars;
pub mod sweph_file;

pub use context::{CalcResult, Ephemeris, EphemerisConfig};
pub use error::Error;
pub use flags::CalcFlags;
pub use types::{
    Body, CalendarType, DegreeParts, DeltaT, EphemerisSource, HouseSystem, JdTt, JdUt1,
    SiderealMode, UtcComponents, UtcToJd,
};

pub type Result<T> = std::result::Result<T, Error>;
