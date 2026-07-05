use std::ffi::c_char;

use swisseph::Error;

/// Error codes returned by all FFI functions. 0 = success, negative = error.
/// Codes are append-only — never reorder or reassign existing values.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SweErrorCode {
    /// Success.
    Ok = 0,
    /// Unknown or out-of-range body number (`ipl`).
    InvalidBody = -1,
    /// Requested flag combination not supported (e.g. BARYCTR on Moshier).
    UnsupportedFlags = -2,
    /// Unknown house system character.
    InvalidHouseSystem = -3,
    /// Unknown sidereal mode index.
    InvalidSiderealMode = -4,
    /// Calendar type not 'g'/'j' (or gregflag not 0/1).
    InvalidCalendarType = -5,
    /// Calendar date does not exist.
    InvalidDate = -6,
    /// No ephemeris file covers the requested epoch or body.
    EphemerisNotAvailable = -7,
    /// Julian Day is outside the valid range for the loaded ephemeris.
    BeyondEphemerisLimits = -8,
    /// Ephemeris file not found on disk.
    FileNotFound = -9,
    /// Ephemeris file is corrupt or has an unsupported format.
    FileFormat = -10,
    /// Body never rises or sets at this location (circumpolar or below horizon).
    CircumpolarBody = -11,
    /// UTC time is invalid (e.g. 25:00 or 61 seconds).
    InvalidTime = -12,
    /// Leap-second specification is invalid.
    InvalidLeapSecond = -13,
    /// Requested ephemeris backend is not compiled in or not configured.
    UnsupportedEphemeris = -14,
    /// Sidereal mode requires fixed-star data that is not loaded.
    SiderealModeRequiresFixedStars = -15,
    /// Error originating from C interop.
    CError = -16,
    /// A Rust panic was caught at the FFI boundary.
    Panic = -90,
    /// A required argument was NULL or otherwise invalid.
    InvalidArg = -91,
    /// Unexpected internal error.
    Internal = -99,
}

pub fn error_code(err: &Error) -> i32 {
    let code = match err {
        Error::InvalidBody(_) => SweErrorCode::InvalidBody,
        Error::UnsupportedFlags(_) => SweErrorCode::UnsupportedFlags,
        Error::InvalidHouseSystem(_) => SweErrorCode::InvalidHouseSystem,
        Error::InvalidSiderealMode(_) => SweErrorCode::InvalidSiderealMode,
        Error::InvalidCalendarType(_) => SweErrorCode::InvalidCalendarType,
        Error::InvalidDate { .. } => SweErrorCode::InvalidDate,
        Error::EphemerisNotAvailable { .. } => SweErrorCode::EphemerisNotAvailable,
        Error::BeyondEphemerisLimits { .. } => SweErrorCode::BeyondEphemerisLimits,
        Error::FileNotFound(_) => SweErrorCode::FileNotFound,
        Error::FileFormat(_) => SweErrorCode::FileFormat,
        Error::CircumpolarBody => SweErrorCode::CircumpolarBody,
        Error::InvalidTime { .. } => SweErrorCode::InvalidTime,
        Error::InvalidLeapSecond { .. } => SweErrorCode::InvalidLeapSecond,
        Error::UnsupportedEphemeris(_) => SweErrorCode::UnsupportedEphemeris,
        Error::SiderealModeRequiresFixedStars(_) => SweErrorCode::SiderealModeRequiresFixedStars,
        Error::CError(_) => SweErrorCode::CError,
    };
    code as i32
}

/// Write a UTF-8 error message into a caller-provided buffer.
/// Always NUL-terminates. Truncates at a char boundary when the message exceeds capacity.
/// No-op when buf is null or cap == 0.
pub unsafe fn write_err(buf: *mut c_char, cap: usize, msg: &str) {
    if buf.is_null() || cap == 0 {
        return;
    }
    let max_bytes = cap - 1; // reserve space for NUL
    let truncated = if msg.len() <= max_bytes {
        msg
    } else {
        // Find the last valid UTF-8 char boundary at or before max_bytes
        let mut end = max_bytes;
        while end > 0 && !msg.is_char_boundary(end) {
            end -= 1;
        }
        &msg[..end]
    };
    unsafe {
        std::ptr::copy_nonoverlapping(truncated.as_ptr(), buf as *mut u8, truncated.len());
        *buf.add(truncated.len()) = 0;
    }
}

/// Wrap an extern "C" fn body in catch_unwind. On panic, writes "panic" to err_buf
/// and returns SweErrorCode::Panic.
macro_rules! ffi_guard {
    ($err_buf:expr, $err_cap:expr, $body:expr) => {{
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| $body));
        match result {
            Ok(val) => val,
            Err(_) => {
                unsafe {
                    $crate::error::write_err($err_buf, $err_cap, "internal panic");
                }
                $crate::error::SweErrorCode::Panic as i32
            }
        }
    }};
}

pub(crate) use ffi_guard;
