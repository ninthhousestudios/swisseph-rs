use std::fmt;
use std::path::PathBuf;

use crate::flags::CalcFlags;
use crate::types::{Body, EphemerisSource, SiderealMode};

/// Errors returned by ephemeris calculations and configuration validation.
#[derive(Debug)]
pub enum Error {
    /// A raw integer body ID has no corresponding [`Body`] variant.
    InvalidBody(i32),
    /// The requested [`CalcFlags`] combination is not supported for the given operation.
    UnsupportedFlags(CalcFlags),
    /// The byte value does not map to a known [`HouseSystem`](crate::types::HouseSystem).
    InvalidHouseSystem(u8),
    /// The integer does not map to a known [`SiderealMode`].
    InvalidSiderealMode(i32),
    /// The integer does not map to a known [`CalendarType`](crate::types::CalendarType).
    InvalidCalendarType(i32),
    /// The year/month/day combination is not a valid calendar date.
    InvalidDate {
        /// Calendar year.
        year: i32,
        /// Calendar month (1-12).
        month: i32,
        /// Day of month (fractional, may include a time-of-day component).
        day: f64,
    },
    /// The requested body is not available from the configured ephemeris source.
    EphemerisNotAvailable {
        /// The body that was requested.
        body: Body,
        /// The ephemeris source it was requested from.
        source: EphemerisSource,
    },
    /// The Julian Day (TT) falls outside the loaded ephemeris file's time range.
    BeyondEphemerisLimits {
        /// The requested Julian Day, TT.
        jd_tt: f64,
        /// Start of the available range (Julian Day, TT).
        start: f64,
        /// End of the available range (Julian Day, TT).
        end: f64,
    },
    /// An ephemeris data file could not be found at the expected path.
    FileNotFound(PathBuf),
    /// An ephemeris data file is malformed or has an unexpected binary layout.
    FileFormat(String),
    /// The body is circumpolar at the given geographic latitude (never rises or never sets).
    CircumpolarBody,
    /// The hour/minute/second combination is not a valid time of day.
    InvalidTime {
        /// Hour of day (0-23).
        hour: i32,
        /// Minute of hour (0-59).
        minute: i32,
        /// Second of minute (fractional; may be up to 61 to allow for leap seconds).
        second: f64,
    },
    /// A leap second was requested on a date that has none in the IERS table.
    InvalidLeapSecond {
        /// Calendar year.
        year: i32,
        /// Calendar month (1-12).
        month: i32,
        /// Day of month.
        day: i32,
    },
    /// The ephemeris backend is compiled out or not yet implemented.
    UnsupportedEphemeris(EphemerisSource),
    /// The sidereal mode requires fixed-star data that is not yet available.
    SiderealModeRequiresFixedStars(SiderealMode),
    /// A catch-all for error messages ported from C's string-buffer error reporting.
    CError(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBody(id) => write!(f, "invalid body ID: {id}"),
            Self::UnsupportedFlags(flags) => write!(f, "unsupported flags: {flags:?}"),
            Self::InvalidHouseSystem(c) => write!(f, "invalid house system: '{}'", *c as char),
            Self::InvalidSiderealMode(id) => write!(f, "invalid sidereal mode: {id}"),
            Self::InvalidCalendarType(id) => write!(f, "invalid calendar type: {id}"),
            Self::InvalidDate { year, month, day } => {
                write!(f, "invalid date: {year}-{month}-{day}")
            }
            Self::EphemerisNotAvailable { body, source } => {
                write!(f, "ephemeris not available for {body:?} from {source:?}")
            }
            Self::BeyondEphemerisLimits { jd_tt, start, end } => {
                write!(f, "JD {jd_tt} outside ephemeris range [{start}, {end}]")
            }
            Self::FileNotFound(path) => write!(f, "file not found: {}", path.display()),
            Self::FileFormat(msg) => write!(f, "file format error: {msg}"),
            Self::CircumpolarBody => write!(f, "body is circumpolar (no rise/set)"),
            Self::InvalidTime {
                hour,
                minute,
                second,
            } => {
                write!(f, "invalid time: {hour}:{minute}:{second}")
            }
            Self::InvalidLeapSecond { year, month, day } => {
                write!(f, "no leap second on {year}-{month:02}-{day:02}")
            }
            Self::UnsupportedEphemeris(source) => {
                write!(f, "ephemeris source {source:?} is not yet supported")
            }
            Self::SiderealModeRequiresFixedStars(mode) => {
                write!(
                    f,
                    "sidereal mode {mode:?} requires the fixed-star subsystem (not yet implemented)"
                )
            }
            Self::CError(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for Error {}
