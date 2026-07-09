// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

//! Binary layout types for Swiss Ephemeris `.se1` ephemeris files.
//!
//! Low-level internals; exposed for golden tests and advanced use.

// Internal body slot indices (SEI_* — indices into pldat[])
/// Internal: slot index for the Earth-Moon barycenter (and Earth/Sun) in `pldat[]`.
pub const SEI_EMB: i32 = 0;
/// Internal: slot index for the Moon in `pldat[]`.
pub const SEI_MOON: i32 = 1;
/// Internal: slot index for Mercury in `pldat[]`.
pub const SEI_MERCURY: i32 = 2;
/// Internal: slot index for Venus in `pldat[]`.
pub const SEI_VENUS: i32 = 3;
/// Internal: slot index for Mars in `pldat[]`.
pub const SEI_MARS: i32 = 4;
/// Internal: slot index for Jupiter in `pldat[]`.
pub const SEI_JUPITER: i32 = 5;
/// Internal: slot index for Saturn in `pldat[]`.
pub const SEI_SATURN: i32 = 6;
/// Internal: slot index for Uranus in `pldat[]`.
pub const SEI_URANUS: i32 = 7;
/// Internal: slot index for Neptune in `pldat[]`.
pub const SEI_NEPTUNE: i32 = 8;
/// Internal: slot index for Pluto in `pldat[]`.
pub const SEI_PLUTO: i32 = 9;
/// Internal: slot index for the barycentric Sun in `pldat[]`.
pub const SEI_SUNBARY: i32 = 10;
/// Internal: slot index shared by numbered asteroids and planetary moons in `pldat[]`.
pub const SEI_ANYBODY: i32 = 11;
/// Internal: slot index for Chiron in `pldat[]`.
pub const SEI_CHIRON: i32 = 12;
/// Internal: slot index for Pholus in `pldat[]`.
pub const SEI_PHOLUS: i32 = 13;
/// Internal: slot index for Ceres in `pldat[]`.
pub const SEI_CERES: i32 = 14;
/// Internal: slot index for Pallas in `pldat[]`.
pub const SEI_PALLAS: i32 = 15;
/// Internal: slot index for Juno in `pldat[]`.
pub const SEI_JUNO: i32 = 16;
/// Internal: slot index for Vesta in `pldat[]`.
pub const SEI_VESTA: i32 = 17;

// Planet data flags (stored as 1 byte per planet in file)
/// Internal: `iflg` bit indicating coordinates are heliocentric (vs. barycentric).
pub const SEI_FLG_HELIO: u32 = 1;
/// Internal: `iflg` bit indicating Chebyshev coefficients are stored in the orbital-plane frame and need rotation to the reference frame before evaluation.
pub const SEI_FLG_ROTATE: u32 = 2;
/// Internal: `iflg` bit indicating reference ellipse coefficients follow the orbital data in the file.
pub const SEI_FLG_ELLIPSE: u32 = 4;
/// Internal: `iflg` bit indicating a heliocentric Earth record is stored in place of the barycentric Sun.
pub const SEI_FLG_EMBHEL: u32 = 8;

// Public body numbers (SE_* — values stored in file ipl[] arrays)
/// Body-number offset added to a numbered asteroid's MPC number to form its file `ipl[]` entry.
pub const SE_AST_OFFSET: i32 = 10000;
/// Body-number offset added to a planetary moon's encoded id to form its file `ipl[]` entry.
pub const SE_PLMOON_OFFSET: i32 = 9000;

/// Internal: raw 4-byte magic value used to detect a file's byte order ("abc" as `0x616263`).
pub(super) const ENDIAN_TEST: u32 = 0x00616263;

/// The category of body data stored in an `.se1` file, inferred from the filename.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileType {
    /// Sun/Moon/Planets file (`sepl*.se1`, `sepla*.se1`).
    Planet,
    /// Moon file (`semo*.se1`).
    Moon,
    /// Main asteroid file (Chiron, Pholus, Ceres, Pallas, Juno, Vesta).
    MainAsteroid,
    /// Individual numbered asteroid file.
    Asteroid,
    /// Planetary moon file.
    PlanetaryMoon,
}

/// Byte order of an `.se1` file's binary section, detected from the endian test value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteOrder {
    /// File stores multi-byte integers and floats in little-endian order.
    Little,
    /// File stores multi-byte integers and floats in big-endian order.
    Big,
}

impl ByteOrder {
    /// Decode a 2-byte integer according to this byte order.
    pub fn read_i16(self, bytes: [u8; 2]) -> i16 {
        match self {
            Self::Little => i16::from_le_bytes(bytes),
            Self::Big => i16::from_be_bytes(bytes),
        }
    }

    /// Decode a 4-byte integer according to this byte order.
    pub fn read_i32(self, bytes: [u8; 4]) -> i32 {
        match self {
            Self::Little => i32::from_le_bytes(bytes),
            Self::Big => i32::from_be_bytes(bytes),
        }
    }

    /// Decode an 8-byte float according to this byte order.
    pub fn read_f64(self, bytes: [u8; 8]) -> f64 {
        match self {
            Self::Little => f64::from_le_bytes(bytes),
            Self::Big => f64::from_be_bytes(bytes),
        }
    }
}

/// Asteroid metadata parsed from an individual asteroid file's MPC orbital elements line.
#[derive(Debug)]
pub struct AsteroidMeta {
    /// Absolute magnitude H.
    pub h: f64,
    /// Slope parameter G (defaults to 0.15 if absent in the file).
    pub g: f64,
    /// Estimated or recorded diameter in kilometers.
    pub diameter_km: f64,
    /// Asteroid name.
    pub name: String,
}

/// Parsed header of an `.se1` ephemeris file: version, type, time range, and byte order.
#[derive(Debug)]
pub struct FileHeader {
    /// File format version number, parsed from the first text line.
    pub version: i32,
    /// Category of body data stored in the file.
    pub file_type: FileType,
    /// File-level Julian Day range `(tfstart, tfend)` for which the file has data.
    pub time_range: (f64, f64),
    /// JPL DE number the file's data derives from.
    pub denum: i32,
    /// Byte order of the file's binary section.
    pub byte_order: ByteOrder,
    /// Asteroid-specific metadata, present only for individual asteroid files.
    pub asteroid: Option<AsteroidMeta>,
}

/// Per-body metadata and orbital elements parsed from an `.se1` file's per-planet block.
#[derive(Debug)]
pub struct PlanetFileData {
    /// Body number as stored in the file's `ipl[]` array (`SE_*`/`SEI_*` namespace).
    pub body_id: i32,
    /// `SEI_FLG_*` bitfield: helio/bary, rotate, ellipse, embhel.
    pub iflg: u32,
    /// Number of Chebyshev coefficients per segment (polynomial order + 1).
    pub ncoe: usize,
    /// Number of coefficients actually significant for the currently loaded segment (may be less than `ncoe`).
    pub neval: usize,
    /// Normalization factor for the packed Chebyshev coefficients.
    pub rmax: f64,
    /// Number of days covered by one polynomial segment.
    pub dseg: f64,
    /// Earliest Julian Day covered by this body's data.
    pub tfstart: f64,
    /// Latest Julian Day covered by this body's data.
    pub tfend: f64,
    /// File offset of the start of this body's segment index.
    pub lndx0: usize,
    /// Number of segment index entries, computed from `tfstart`, `tfend`, `dseg`.
    pub nndx: usize,
    /// Epoch of the orbital elements used for the orbital-plane rotation.
    pub telem: f64,
    /// Interpolated equinoctial inclination-vector component p.
    pub prot: f64,
    /// Interpolated equinoctial inclination-vector component q.
    pub qrot: f64,
    /// Rate of change of `prot`.
    pub dprot: f64,
    /// Rate of change of `qrot`.
    pub dqrot: f64,
    /// Perihelion longitude, present only when `SEI_FLG_ELLIPSE` is set.
    pub peri: f64,
    /// Rate of change of `peri`, present only when `SEI_FLG_ELLIPSE` is set.
    pub dperi: f64,
    /// Reference ellipse Chebyshev coefficients (2×`ncoe` doubles: X then Y), if `SEI_FLG_ELLIPSE` is set.
    pub refep: Option<Vec<f64>>,
}
