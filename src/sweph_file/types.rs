// Internal body slot indices (SEI_* — indices into pldat[])
pub const SEI_EMB: i32 = 0;
pub const SEI_MOON: i32 = 1;
pub const SEI_MERCURY: i32 = 2;
pub const SEI_VENUS: i32 = 3;
pub const SEI_MARS: i32 = 4;
pub const SEI_JUPITER: i32 = 5;
pub const SEI_SATURN: i32 = 6;
pub const SEI_URANUS: i32 = 7;
pub const SEI_NEPTUNE: i32 = 8;
pub const SEI_PLUTO: i32 = 9;
pub const SEI_SUNBARY: i32 = 10;
pub const SEI_ANYBODY: i32 = 11;
pub const SEI_CHIRON: i32 = 12;
pub const SEI_PHOLUS: i32 = 13;
pub const SEI_CERES: i32 = 14;
pub const SEI_PALLAS: i32 = 15;
pub const SEI_JUNO: i32 = 16;
pub const SEI_VESTA: i32 = 17;

// Planet data flags (stored as 1 byte per planet in file)
pub const SEI_FLG_HELIO: u32 = 1;
pub const SEI_FLG_ROTATE: u32 = 2;
pub const SEI_FLG_ELLIPSE: u32 = 4;
pub const SEI_FLG_EMBHEL: u32 = 8;

// Public body numbers (SE_* — values stored in file ipl[] arrays)
pub const SE_AST_OFFSET: i32 = 10000;
pub const SE_PLMOON_OFFSET: i32 = 9000;

// File endianness magic value
pub(super) const ENDIAN_TEST: u32 = 0x00616263;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileType {
    Planet,
    Moon,
    MainAsteroid,
    Asteroid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteOrder {
    Little,
    Big,
}

impl ByteOrder {
    pub fn read_i16(self, bytes: [u8; 2]) -> i16 {
        match self {
            Self::Little => i16::from_le_bytes(bytes),
            Self::Big => i16::from_be_bytes(bytes),
        }
    }

    pub fn read_i32(self, bytes: [u8; 4]) -> i32 {
        match self {
            Self::Little => i32::from_le_bytes(bytes),
            Self::Big => i32::from_be_bytes(bytes),
        }
    }

    pub fn read_f64(self, bytes: [u8; 8]) -> f64 {
        match self {
            Self::Little => f64::from_le_bytes(bytes),
            Self::Big => f64::from_be_bytes(bytes),
        }
    }
}

#[derive(Debug)]
pub struct FileHeader {
    pub version: i32,
    pub file_type: FileType,
    pub time_range: (f64, f64),
    pub denum: i32,
    pub byte_order: ByteOrder,
}

#[derive(Debug)]
pub struct PlanetFileData {
    pub body_id: i32,
    pub iflg: u32,
    pub ncoe: usize,
    pub neval: usize,
    pub rmax: f64,
    pub dseg: f64,
    pub tfstart: f64,
    pub tfend: f64,
    pub lndx0: usize,
    pub nndx: usize,
    pub telem: f64,
    pub prot: f64,
    pub qrot: f64,
    pub dprot: f64,
    pub dqrot: f64,
    pub peri: f64,
    pub dperi: f64,
    pub refep: Option<Vec<f64>>,
}
