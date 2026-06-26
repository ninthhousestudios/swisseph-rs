use std::fmt;
use std::path::PathBuf;

use crate::flags::CalcFlags;
use crate::types::{Body, EphemerisSource};

#[derive(Debug)]
pub enum Error {
    InvalidBody(i32),
    UnsupportedFlags(CalcFlags),
    InvalidHouseSystem(u8),
    InvalidSiderealMode(i32),
    InvalidCalendarType(i32),
    InvalidDate { year: i32, month: i32, day: f64 },
    EphemerisNotAvailable { body: Body, source: EphemerisSource },
    BeyondEphemerisLimits { jd_tt: f64, start: f64, end: f64 },
    FileNotFound(PathBuf),
    FileFormat(String),
    CircumpolarBody,
    InvalidTime { hour: i32, minute: i32, second: f64 },
    InvalidLeapSecond { year: i32, month: i32, day: i32 },
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
            Self::CError(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for Error {}
