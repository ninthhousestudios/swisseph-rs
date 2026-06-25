use std::ops::{Add, Sub};

use crate::constants;

// ---------------------------------------------------------------------------
// Body ID newtypes — private inner fields enforce range invariants
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FictitiousId(i32);

impl FictitiousId {
    pub fn new(raw_id: i32) -> crate::Result<Self> {
        if (constants::FICT_OFFSET..=constants::FICT_MAX).contains(&raw_id) {
            Ok(Self(raw_id))
        } else {
            Err(crate::Error::InvalidBody(raw_id))
        }
    }

    pub fn raw_id(self) -> i32 {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AsteroidId(i32);

impl AsteroidId {
    pub fn new(mpc_number: i32) -> crate::Result<Self> {
        if mpc_number >= 0 {
            Ok(Self(mpc_number))
        } else {
            Err(crate::Error::InvalidBody(mpc_number))
        }
    }

    pub fn mpc_number(self) -> i32 {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PlanetMoonId(i32);

impl PlanetMoonId {
    pub fn new(encoded: i32) -> crate::Result<Self> {
        if (0..=999).contains(&encoded) {
            Ok(Self(encoded))
        } else {
            Err(crate::Error::InvalidBody(encoded))
        }
    }

    pub fn encoded(self) -> i32 {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CometId(i32);

impl CometId {
    pub fn new(number: i32) -> crate::Result<Self> {
        if (0..=7999).contains(&number) {
            Ok(Self(number))
        } else {
            Err(crate::Error::InvalidBody(number))
        }
    }

    pub fn number(self) -> i32 {
        self.0
    }
}

// ---------------------------------------------------------------------------
// Body
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Body {
    Sun,
    Moon,
    Mercury,
    Venus,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    MeanNode,
    TrueNode,
    MeanApogee,
    OscuApogee,
    Earth,
    Chiron,
    Pholus,
    Ceres,
    Pallas,
    Juno,
    Vesta,
    IntpApogee,
    IntpPerigee,
    Fictitious(FictitiousId),
    Asteroid(AsteroidId),
    PlanetMoon(PlanetMoonId),
    Comet(CometId),
    EclipticNutation,
}

impl Body {
    pub fn fictitious(raw_id: i32) -> crate::Result<Self> {
        Ok(Self::Fictitious(FictitiousId::new(raw_id)?))
    }

    pub fn asteroid(mpc_number: i32) -> crate::Result<Self> {
        Ok(Self::Asteroid(AsteroidId::new(mpc_number)?))
    }

    pub fn planet_moon(encoded: i32) -> crate::Result<Self> {
        Ok(Self::PlanetMoon(PlanetMoonId::new(encoded)?))
    }

    pub fn comet(number: i32) -> crate::Result<Self> {
        Ok(Self::Comet(CometId::new(number)?))
    }

    pub fn to_raw_id(self) -> i32 {
        match self {
            Self::Sun => 0,
            Self::Moon => 1,
            Self::Mercury => 2,
            Self::Venus => 3,
            Self::Mars => 4,
            Self::Jupiter => 5,
            Self::Saturn => 6,
            Self::Uranus => 7,
            Self::Neptune => 8,
            Self::Pluto => 9,
            Self::MeanNode => 10,
            Self::TrueNode => 11,
            Self::MeanApogee => 12,
            Self::OscuApogee => 13,
            Self::Earth => 14,
            Self::Chiron => 15,
            Self::Pholus => 16,
            Self::Ceres => 17,
            Self::Pallas => 18,
            Self::Juno => 19,
            Self::Vesta => 20,
            Self::IntpApogee => 21,
            Self::IntpPerigee => 22,
            Self::Fictitious(id) => id.raw_id(),
            Self::Asteroid(id) => constants::AST_OFFSET + id.mpc_number(),
            Self::PlanetMoon(id) => constants::PLMOON_OFFSET + id.encoded(),
            Self::Comet(id) => constants::COMET_OFFSET + id.number(),
            Self::EclipticNutation => -1,
        }
    }
}

impl TryFrom<i32> for Body {
    type Error = crate::Error;

    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            -1 => Ok(Self::EclipticNutation),
            0 => Ok(Self::Sun),
            1 => Ok(Self::Moon),
            2 => Ok(Self::Mercury),
            3 => Ok(Self::Venus),
            4 => Ok(Self::Mars),
            5 => Ok(Self::Jupiter),
            6 => Ok(Self::Saturn),
            7 => Ok(Self::Uranus),
            8 => Ok(Self::Neptune),
            9 => Ok(Self::Pluto),
            10 => Ok(Self::MeanNode),
            11 => Ok(Self::TrueNode),
            12 => Ok(Self::MeanApogee),
            13 => Ok(Self::OscuApogee),
            14 => Ok(Self::Earth),
            15 => Ok(Self::Chiron),
            16 => Ok(Self::Pholus),
            17 => Ok(Self::Ceres),
            18 => Ok(Self::Pallas),
            19 => Ok(Self::Juno),
            20 => Ok(Self::Vesta),
            21 => Ok(Self::IntpApogee),
            22 => Ok(Self::IntpPerigee),
            40..=999 => Ok(Self::Fictitious(FictitiousId(v))),
            1000..=8999 => Ok(Self::Comet(CometId(v - constants::COMET_OFFSET))),
            9000..=9999 => Ok(Self::PlanetMoon(PlanetMoonId(v - constants::PLMOON_OFFSET))),
            n if n >= constants::AST_OFFSET => {
                Ok(Self::Asteroid(AsteroidId(n - constants::AST_OFFSET)))
            }
            _ => Err(crate::Error::InvalidBody(v)),
        }
    }
}

// ---------------------------------------------------------------------------
// FictitiousBody — named companion for Body::Fictitious (IDs 40-58)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum FictitiousBody {
    Cupido = 40,
    Hades = 41,
    Zeus = 42,
    Kronos = 43,
    Apollon = 44,
    Admetos = 45,
    Vulkanus = 46,
    Poseidon = 47,
    Isis = 48,
    Nibiru = 49,
    Harrington = 50,
    NeptuneLeverrier = 51,
    NeptuneAdams = 52,
    PlutoLowell = 53,
    PlutoPickering = 54,
    Vulcan = 55,
    WhiteMoon = 56,
    Proserpina = 57,
    Waldemath = 58,
}

impl From<FictitiousBody> for Body {
    fn from(f: FictitiousBody) -> Self {
        Body::Fictitious(FictitiousId(f as i32))
    }
}

impl TryFrom<i32> for FictitiousBody {
    type Error = crate::Error;

    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            40 => Ok(Self::Cupido),
            41 => Ok(Self::Hades),
            42 => Ok(Self::Zeus),
            43 => Ok(Self::Kronos),
            44 => Ok(Self::Apollon),
            45 => Ok(Self::Admetos),
            46 => Ok(Self::Vulkanus),
            47 => Ok(Self::Poseidon),
            48 => Ok(Self::Isis),
            49 => Ok(Self::Nibiru),
            50 => Ok(Self::Harrington),
            51 => Ok(Self::NeptuneLeverrier),
            52 => Ok(Self::NeptuneAdams),
            53 => Ok(Self::PlutoLowell),
            54 => Ok(Self::PlutoPickering),
            55 => Ok(Self::Vulcan),
            56 => Ok(Self::WhiteMoon),
            57 => Ok(Self::Proserpina),
            58 => Ok(Self::Waldemath),
            _ => Err(crate::Error::InvalidBody(v)),
        }
    }
}

// ---------------------------------------------------------------------------
// HouseSystem
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HouseSystem {
    Equal,
    Alcabitius,
    Campanus,
    EqualMC,
    Carter,
    Gauquelin,
    Horizon,
    Sunshine,
    SunshineAlt,
    SavardA,
    Koch,
    PullenSD,
    Morinus,
    EqualAries,
    Porphyry,
    Placidus,
    PullenSR,
    Regiomontanus,
    Sripati,
    PolichPage,
    KrusinskiPisaGoelzer,
    Vehlow,
    WholeSign,
    Meridian,
    APC,
}

impl HouseSystem {
    pub fn to_char(self) -> u8 {
        match self {
            Self::Equal => b'A',
            Self::Alcabitius => b'B',
            Self::Campanus => b'C',
            Self::EqualMC => b'D',
            Self::Carter => b'F',
            Self::Gauquelin => b'G',
            Self::Horizon => b'H',
            Self::Sunshine => b'I',
            Self::SunshineAlt => b'i',
            Self::SavardA => b'J',
            Self::Koch => b'K',
            Self::PullenSD => b'L',
            Self::Morinus => b'M',
            Self::EqualAries => b'N',
            Self::Porphyry => b'O',
            Self::Placidus => b'P',
            Self::PullenSR => b'Q',
            Self::Regiomontanus => b'R',
            Self::Sripati => b'S',
            Self::PolichPage => b'T',
            Self::KrusinskiPisaGoelzer => b'U',
            Self::Vehlow => b'V',
            Self::WholeSign => b'W',
            Self::Meridian => b'X',
            Self::APC => b'Y',
        }
    }
}

impl TryFrom<u8> for HouseSystem {
    type Error = crate::Error;

    fn try_from(v: u8) -> std::result::Result<Self, Self::Error> {
        match v {
            b'A' | b'E' => Ok(Self::Equal),
            b'B' => Ok(Self::Alcabitius),
            b'C' => Ok(Self::Campanus),
            b'D' => Ok(Self::EqualMC),
            b'F' => Ok(Self::Carter),
            b'G' => Ok(Self::Gauquelin),
            b'H' => Ok(Self::Horizon),
            b'I' => Ok(Self::Sunshine),
            b'i' => Ok(Self::SunshineAlt),
            b'J' => Ok(Self::SavardA),
            b'K' => Ok(Self::Koch),
            b'L' => Ok(Self::PullenSD),
            b'M' => Ok(Self::Morinus),
            b'N' => Ok(Self::EqualAries),
            b'O' => Ok(Self::Porphyry),
            b'P' => Ok(Self::Placidus),
            b'Q' => Ok(Self::PullenSR),
            b'R' => Ok(Self::Regiomontanus),
            b'S' => Ok(Self::Sripati),
            b'T' => Ok(Self::PolichPage),
            b'U' => Ok(Self::KrusinskiPisaGoelzer),
            b'V' => Ok(Self::Vehlow),
            b'W' => Ok(Self::WholeSign),
            b'X' => Ok(Self::Meridian),
            b'Y' => Ok(Self::APC),
            _ => Err(crate::Error::InvalidHouseSystem(v)),
        }
    }
}

// ---------------------------------------------------------------------------
// CalendarType
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum CalendarType {
    Julian = 0,
    Gregorian = 1,
}

impl TryFrom<i32> for CalendarType {
    type Error = crate::Error;

    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(Self::Julian),
            1 => Ok(Self::Gregorian),
            _ => Err(crate::Error::InvalidCalendarType(v)),
        }
    }
}

// ---------------------------------------------------------------------------
// SiderealMode
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum SiderealMode {
    FaganBradley = 0,
    Lahiri = 1,
    DeLuce = 2,
    Raman = 3,
    Ushashashi = 4,
    Krishnamurti = 5,
    DjwhalKhul = 6,
    Yukteshwar = 7,
    JnBhasin = 8,
    BabylKugler1 = 9,
    BabylKugler2 = 10,
    BabylKugler3 = 11,
    BabylHuber = 12,
    BabylEtpsc = 13,
    Aldebaran15Tau = 14,
    Hipparchos = 15,
    Sassanian = 16,
    GalCent0Sag = 17,
    J2000 = 18,
    J1900 = 19,
    B1950 = 20,
    Suryasiddhanta = 21,
    SuryasiddhantaMsun = 22,
    Aryabhata = 23,
    AryabhataMsun = 24,
    SsRevati = 25,
    SsCitra = 26,
    TrueCitra = 27,
    TrueRevati = 28,
    TruePushya = 29,
    GalCentRgilbrand = 30,
    GalEquIau1958 = 31,
    GalEquTrue = 32,
    GalEquMula = 33,
    GalAlignMardyks = 34,
    TrueMula = 35,
    GalCentMulaWilhelm = 36,
    Aryabhata522 = 37,
    BabylBritton = 38,
    TrueSheoran = 39,
    GalCentCochrane = 40,
    GalEquFiorenza = 41,
    ValensMoon = 42,
    Lahiri1940 = 43,
    LahiriVp285 = 44,
    KrishnamurtiVp291 = 45,
    LahiriIcrc = 46,
    User = 255,
}

impl TryFrom<i32> for SiderealMode {
    type Error = crate::Error;

    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(Self::FaganBradley),
            1 => Ok(Self::Lahiri),
            2 => Ok(Self::DeLuce),
            3 => Ok(Self::Raman),
            4 => Ok(Self::Ushashashi),
            5 => Ok(Self::Krishnamurti),
            6 => Ok(Self::DjwhalKhul),
            7 => Ok(Self::Yukteshwar),
            8 => Ok(Self::JnBhasin),
            9 => Ok(Self::BabylKugler1),
            10 => Ok(Self::BabylKugler2),
            11 => Ok(Self::BabylKugler3),
            12 => Ok(Self::BabylHuber),
            13 => Ok(Self::BabylEtpsc),
            14 => Ok(Self::Aldebaran15Tau),
            15 => Ok(Self::Hipparchos),
            16 => Ok(Self::Sassanian),
            17 => Ok(Self::GalCent0Sag),
            18 => Ok(Self::J2000),
            19 => Ok(Self::J1900),
            20 => Ok(Self::B1950),
            21 => Ok(Self::Suryasiddhanta),
            22 => Ok(Self::SuryasiddhantaMsun),
            23 => Ok(Self::Aryabhata),
            24 => Ok(Self::AryabhataMsun),
            25 => Ok(Self::SsRevati),
            26 => Ok(Self::SsCitra),
            27 => Ok(Self::TrueCitra),
            28 => Ok(Self::TrueRevati),
            29 => Ok(Self::TruePushya),
            30 => Ok(Self::GalCentRgilbrand),
            31 => Ok(Self::GalEquIau1958),
            32 => Ok(Self::GalEquTrue),
            33 => Ok(Self::GalEquMula),
            34 => Ok(Self::GalAlignMardyks),
            35 => Ok(Self::TrueMula),
            36 => Ok(Self::GalCentMulaWilhelm),
            37 => Ok(Self::Aryabhata522),
            38 => Ok(Self::BabylBritton),
            39 => Ok(Self::TrueSheoran),
            40 => Ok(Self::GalCentCochrane),
            41 => Ok(Self::GalEquFiorenza),
            42 => Ok(Self::ValensMoon),
            43 => Ok(Self::Lahiri1940),
            44 => Ok(Self::LahiriVp285),
            45 => Ok(Self::KrishnamurtiVp291),
            46 => Ok(Self::LahiriIcrc),
            255 => Ok(Self::User),
            _ => Err(crate::Error::InvalidSiderealMode(v)),
        }
    }
}

// ---------------------------------------------------------------------------
// EphemerisSource
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EphemerisSource {
    Jpl,
    Swiss,
    Moshier,
}

// ---------------------------------------------------------------------------
// Astronomical model enums
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum PrecessionModel {
    IAU1976 = 1,
    Laskar1986 = 2,
    WillEpsLask = 3,
    Williams1994 = 4,
    Simon1994 = 5,
    IAU2000 = 6,
    Bretagnon2003 = 7,
    IAU2006 = 8,
    Vondrak2011 = 9,
    Owen1990 = 10,
    Newcomb = 11,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum NutationModel {
    IAU1980 = 1,
    IAUCorr1987 = 2,
    IAU2000A = 3,
    IAU2000B = 4,
    Woolard = 5,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum DeltaTModel {
    StephensonMorrison1984 = 1,
    Stephenson1997 = 2,
    StephensonMorrison2004 = 3,
    EspenakMeeus2006 = 4,
    StephensonEtc2016 = 5,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum SiderealTimeModel {
    IAU1976 = 1,
    IAU2006 = 2,
    IersConv2010 = 3,
    Longterm = 4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum BiasModel {
    None = 1,
    IAU2000 = 2,
    IAU2006 = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum JplHorMode {
    LongAgreement = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum JplHoraMode {
    V1 = 1,
    V2 = 2,
    V3 = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AstroModels {
    pub delta_t: DeltaTModel,
    pub prec_longterm: PrecessionModel,
    pub prec_shortterm: PrecessionModel,
    pub nutation: NutationModel,
    pub bias: BiasModel,
    pub jplhor_mode: JplHorMode,
    pub jplhora_mode: JplHoraMode,
    pub sidereal_time: SiderealTimeModel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FrameTransform {
    J2000ToGcrs,
    GcrsToJ2000,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Epsilon {
    pub eps: f64,
    pub sin_eps: f64,
    pub cos_eps: f64,
}

impl Epsilon {
    pub fn new(eps_rad: f64) -> Self {
        Self {
            eps: eps_rad,
            sin_eps: eps_rad.sin(),
            cos_eps: eps_rad.cos(),
        }
    }
}

// ---------------------------------------------------------------------------
// Julian Day newtypes
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct JdTt(pub f64);

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct JdUt1(pub f64);

macro_rules! impl_jd_ops {
    ($T:ty) => {
        impl Add<f64> for $T {
            type Output = Self;
            fn add(self, rhs: f64) -> Self {
                Self(self.0 + rhs)
            }
        }
        impl Sub<f64> for $T {
            type Output = Self;
            fn sub(self, rhs: f64) -> Self {
                Self(self.0 - rhs)
            }
        }
        impl Sub for $T {
            type Output = f64;
            fn sub(self, rhs: Self) -> f64 {
                self.0 - rhs.0
            }
        }
    };
}

impl_jd_ops!(JdTt);
impl_jd_ops!(JdUt1);

// ---------------------------------------------------------------------------
// UTC components
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct UtcComponents {
    pub year: i32,
    pub month: i32,
    pub day: i32,
    pub hour: i32,
    pub minute: i32,
    pub second: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct UtcToJd {
    pub tt: JdTt,
    pub ut1: JdUt1,
}

// ---------------------------------------------------------------------------
// DeltaT trait
// ---------------------------------------------------------------------------

pub trait DeltaT {
    fn delta_t(&self, jd_ut: JdUt1) -> f64;
}

// ---------------------------------------------------------------------------
// DegreeParts — result of swe_split_deg
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DegreeParts {
    pub degrees: i32,
    pub minutes: i32,
    pub seconds: i32,
    pub second_fraction: f64,
    pub sign: i32,
}

// ---------------------------------------------------------------------------
// Astronomical model enums — defaults
// ---------------------------------------------------------------------------

impl Default for AstroModels {
    fn default() -> Self {
        Self {
            delta_t: DeltaTModel::StephensonEtc2016,
            prec_longterm: PrecessionModel::Vondrak2011,
            prec_shortterm: PrecessionModel::Vondrak2011,
            nutation: NutationModel::IAU2000B,
            bias: BiasModel::IAU2006,
            jplhor_mode: JplHorMode::LongAgreement,
            jplhora_mode: JplHoraMode::V3,
            sidereal_time: SiderealTimeModel::Longterm,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fictitious_id_valid_range() {
        assert!(FictitiousId::new(40).is_ok());
        assert!(FictitiousId::new(999).is_ok());
        assert!(FictitiousId::new(500).is_ok());
    }

    #[test]
    fn fictitious_id_rejects_out_of_range() {
        assert!(FictitiousId::new(39).is_err());
        assert!(FictitiousId::new(1000).is_err());
        assert!(FictitiousId::new(0).is_err());
        assert!(FictitiousId::new(-1).is_err());
        assert!(FictitiousId::new(23).is_err());
    }

    #[test]
    fn asteroid_id_rejects_negative() {
        assert!(AsteroidId::new(-1).is_err());
        assert!(AsteroidId::new(-10000).is_err());
        assert!(AsteroidId::new(0).is_ok());
        assert!(AsteroidId::new(1).is_ok());
    }

    #[test]
    fn planet_moon_id_valid_range() {
        assert!(PlanetMoonId::new(0).is_ok());
        assert!(PlanetMoonId::new(999).is_ok());
        assert!(PlanetMoonId::new(1000).is_err());
        assert!(PlanetMoonId::new(-1).is_err());
    }

    #[test]
    fn comet_id_valid_range() {
        assert!(CometId::new(0).is_ok());
        assert!(CometId::new(7999).is_ok());
        assert!(CometId::new(8000).is_err());
        assert!(CometId::new(-1).is_err());
    }

    #[test]
    fn body_constructors_validate() {
        assert!(Body::fictitious(40).is_ok());
        assert!(Body::fictitious(23).is_err());
        assert!(Body::asteroid(0).is_ok());
        assert!(Body::asteroid(-10000).is_err());
        assert!(Body::planet_moon(0).is_ok());
        assert!(Body::planet_moon(1000).is_err());
        assert!(Body::comet(0).is_ok());
        assert!(Body::comet(8000).is_err());
    }

    #[test]
    fn body_no_aliasing_via_constructors() {
        let sun_raw = Body::Sun.to_raw_id();
        for n in [-10000, -1, -100] {
            assert!(Body::asteroid(n).is_err());
        }
        let asteroid_0 = Body::asteroid(0).unwrap();
        assert_ne!(asteroid_0.to_raw_id(), sun_raw);
        assert_eq!(asteroid_0.to_raw_id(), constants::AST_OFFSET);
    }

    #[test]
    fn body_try_from_roundtrip() {
        for raw in [
            40, 58, 500, 999, 1000, 5000, 8999, 9000, 9500, 9999, 10000, 20000,
        ] {
            let body = Body::try_from(raw).unwrap();
            assert_eq!(body.to_raw_id(), raw);
        }
    }

    #[test]
    fn body_try_from_gap_rejected() {
        for raw in [23, 24, 30, 39] {
            assert!(Body::try_from(raw).is_err());
        }
    }
}
