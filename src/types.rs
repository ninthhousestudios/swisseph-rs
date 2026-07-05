//! Core value types: body identifiers, house systems, calendar/sidereal/model enums,
//! Julian Day newtypes, and small result structs shared across the crate's public API.

use std::ops::{Add, Sub};

use crate::constants;

// ---------------------------------------------------------------------------
// Body ID newtypes — private inner fields enforce range invariants
// ---------------------------------------------------------------------------

/// Validated ID for a fictitious (hypothetical) planet, range 40–999.
/// C equivalent: `SE_FICT_OFFSET` + index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FictitiousId(i32);

impl FictitiousId {
    /// Validates and constructs a `FictitiousId` from a raw C body ID (40–999).
    pub fn new(raw_id: i32) -> crate::Result<Self> {
        if (constants::FICT_OFFSET..=constants::FICT_MAX).contains(&raw_id) {
            Ok(Self(raw_id))
        } else {
            Err(crate::Error::InvalidBody(raw_id))
        }
    }

    /// Returns the raw C body ID.
    pub fn raw_id(self) -> i32 {
        self.0
    }
}

/// Validated MPC number for a numbered asteroid (>= 0).
/// C equivalent: `SE_AST_OFFSET` + mpc_number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AsteroidId(i32);

impl AsteroidId {
    /// Validates and constructs an `AsteroidId` from an MPC catalog number (>= 0).
    pub fn new(mpc_number: i32) -> crate::Result<Self> {
        if mpc_number >= 0 {
            Ok(Self(mpc_number))
        } else {
            Err(crate::Error::InvalidBody(mpc_number))
        }
    }

    /// Returns the MPC catalog number.
    pub fn mpc_number(self) -> i32 {
        self.0
    }
}

/// Validated encoded ID for a planetary moon (0–999, where 9n99 = center-of-body).
/// C equivalent: `SE_PLMOON_OFFSET` + encoded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlanetMoonId(i32);

impl PlanetMoonId {
    /// Validates and constructs a `PlanetMoonId` from an encoded planet/moon value (0–999).
    pub fn new(encoded: i32) -> crate::Result<Self> {
        if (0..=999).contains(&encoded) {
            Ok(Self(encoded))
        } else {
            Err(crate::Error::InvalidBody(encoded))
        }
    }

    /// Returns the raw encoded planet/moon value.
    pub fn encoded(self) -> i32 {
        self.0
    }
}

// ---------------------------------------------------------------------------
// Body
// ---------------------------------------------------------------------------

/// A celestial body or computational pseudo-body for ephemeris calculations.
///
/// Fixed variants (Sun through Vesta) correspond to C's `SE_SUN` through `SE_VESTA`.
/// Parameterized variants handle fictitious planets, numbered asteroids, and planetary moons.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Body {
    /// The Sun.
    Sun,
    /// The Moon.
    Moon,
    /// Mercury.
    Mercury,
    /// Venus.
    Venus,
    /// Mars.
    Mars,
    /// Jupiter.
    Jupiter,
    /// Saturn.
    Saturn,
    /// Uranus.
    Uranus,
    /// Neptune.
    Neptune,
    /// Pluto.
    Pluto,
    /// Mean lunar node.
    MeanNode,
    /// True (osculating) lunar node.
    TrueNode,
    /// Mean lunar apogee (mean "Black Moon" / Lilith).
    MeanApogee,
    /// Osculating lunar apogee (true "Black Moon" / Lilith).
    OscuApogee,
    /// The Earth (heliocentric/barycentric calculations).
    Earth,
    /// Asteroid/comet 2060 Chiron.
    Chiron,
    /// Asteroid 5145 Pholus.
    Pholus,
    /// Asteroid 1 Ceres.
    Ceres,
    /// Asteroid 2 Pallas.
    Pallas,
    /// Asteroid 3 Juno.
    Juno,
    /// Asteroid 4 Vesta.
    Vesta,
    /// Interpolated lunar apogee.
    IntpApogee,
    /// Interpolated lunar perigee.
    IntpPerigee,
    /// Named fictitious (hypothetical) planet, keyed by [`FictitiousId`].
    Fictitious(FictitiousId),
    /// Numbered minor planet, keyed by [`AsteroidId`] (MPC catalog number).
    Asteroid(AsteroidId),
    /// Planetary moon, keyed by [`PlanetMoonId`] (encoded planet/moon pair).
    PlanetMoon(PlanetMoonId),
    /// Pseudo-body representing the ecliptic and nutation (not a physical body).
    EclipticNutation,
}

impl Body {
    /// Constructs a `Body::Fictitious` from a raw fictitious-planet ID (40–999).
    pub fn fictitious(raw_id: i32) -> crate::Result<Self> {
        Ok(Self::Fictitious(FictitiousId::new(raw_id)?))
    }

    /// Constructs a `Body::Asteroid` from an MPC catalog number.
    pub fn asteroid(mpc_number: i32) -> crate::Result<Self> {
        Ok(Self::Asteroid(AsteroidId::new(mpc_number)?))
    }

    /// Constructs a `Body::PlanetMoon` from an encoded planet/moon value.
    pub fn planet_moon(encoded: i32) -> crate::Result<Self> {
        Ok(Self::PlanetMoon(PlanetMoonId::new(encoded)?))
    }

    /// Converts this `Body` to the raw C body ID used by `swe_calc` and friends.
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

/// Named fictitious (hypothetical) planets from the built-in Neely catalog (IDs 40–58).
/// Converts to [`Body`] via `From<FictitiousBody>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum FictitiousBody {
    /// Cupido (Hamburg school hypothetical planet).
    Cupido = 40,
    /// Hades (Hamburg school hypothetical planet).
    Hades = 41,
    /// Zeus (Hamburg school hypothetical planet).
    Zeus = 42,
    /// Kronos (Hamburg school hypothetical planet).
    Kronos = 43,
    /// Apollon (Hamburg school hypothetical planet).
    Apollon = 44,
    /// Admetos (Hamburg school hypothetical planet).
    Admetos = 45,
    /// Vulkanus (Hamburg school hypothetical planet).
    Vulkanus = 46,
    /// Poseidon (Hamburg school hypothetical planet).
    Poseidon = 47,
    /// Isis (Sepharial's hypothetical planet).
    Isis = 48,
    /// Nibiru (hypothetical planet).
    Nibiru = 49,
    /// Harrington (hypothetical trans-Neptunian planet).
    Harrington = 50,
    /// Neptune per Le Verrier's original orbital elements.
    NeptuneLeverrier = 51,
    /// Neptune per Adams' original orbital elements.
    NeptuneAdams = 52,
    /// Pluto per Lowell's predicted orbital elements.
    PlutoLowell = 53,
    /// Pluto per Pickering's predicted orbital elements.
    PlutoPickering = 54,
    /// Vulcan (hypothetical intra-Mercurial planet).
    Vulcan = 55,
    /// Selena/White Moon (hypothetical lunar-derived point).
    WhiteMoon = 56,
    /// Proserpina (hypothetical planet).
    Proserpina = 57,
    /// Waldemath's second (hypothetical dark) Moon.
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

/// Astrological house system. Each variant maps to the single-character code used by
/// C's `swe_houses` family (accessible via [`HouseSystem::to_char`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum HouseSystem {
    /// Equal houses from the Ascendant.
    Equal,
    /// Alcabitius house system.
    Alcabitius,
    /// Campanus house system.
    Campanus,
    /// Equal houses from the MC.
    EqualMC,
    /// Carter poli-equatorial house system.
    Carter,
    /// Gauquelin sectors (36 divisions).
    Gauquelin,
    /// Horizon/azimuth-based house system.
    Horizon,
    /// Sunshine house system.
    Sunshine,
    /// Sunshine house system, alternate method.
    SunshineAlt,
    /// Savard-A house system.
    SavardA,
    /// Koch house system.
    Koch,
    /// Pullen SD (sinusoidal delta) house system.
    PullenSD,
    /// Morinus house system.
    Morinus,
    /// Equal houses with house 1 starting at 0° Aries.
    EqualAries,
    /// Porphyry house system.
    Porphyry,
    /// Placidus house system.
    Placidus,
    /// Pullen SR (sinusoidal ratio) house system.
    PullenSR,
    /// Regiomontanus house system.
    Regiomontanus,
    /// Sripati house system.
    Sripati,
    /// Polich/Page (topocentric) house system.
    PolichPage,
    /// Krusinski-Pisa-Goelzer house system.
    KrusinskiPisaGoelzer,
    /// Equal houses, Vehlow variant (cusp 1 at Ascendant - 15°).
    Vehlow,
    /// Whole-sign house system.
    WholeSign,
    /// Axial rotation (Meridian) house system.
    Meridian,
    /// APC (astrological process control) house system.
    APC,
}

impl HouseSystem {
    /// Returns the single-character house system code used by C's `swe_houses` family.
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

    /// Returns the human-readable name of the house system.
    pub fn name(self) -> &'static str {
        match self {
            Self::Equal => "equal",
            Self::Alcabitius => "Alcabitius",
            Self::Campanus => "Campanus",
            Self::EqualMC => "equal (MC)",
            Self::Carter => "Carter poli-equ.",
            Self::Gauquelin => "Gauquelin sectors",
            Self::Horizon => "horizon/azimut",
            Self::Sunshine => "Sunshine",
            Self::SunshineAlt => "Sunshine/alt.",
            Self::SavardA => "Savard-A",
            Self::Koch => "Koch",
            Self::PullenSD => "Pullen SD",
            Self::Morinus => "Morinus",
            Self::EqualAries => "equal/1=Aries",
            Self::Porphyry => "Porphyry",
            Self::Placidus => "Placidus",
            Self::PullenSR => "Pullen SR",
            Self::Regiomontanus => "Regiomontanus",
            Self::Sripati => "Sripati",
            Self::PolichPage => "Polich/Page",
            Self::KrusinskiPisaGoelzer => "Krusinski-Pisa-Goelzer",
            Self::Vehlow => "equal/Vehlow",
            Self::WholeSign => "equal/ whole sign",
            Self::Meridian => "axial rotation system/Meridian houses",
            Self::APC => "APC houses",
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

/// Calendar system for Julian Day conversions. C: `SE_JUL_CAL` / `SE_GREG_CAL`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum CalendarType {
    /// Julian calendar.
    Julian = 0,
    /// Gregorian calendar.
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

/// Sidereal zodiac ayanamsa definition. Each variant defines a fixed reference point
/// that anchors the sidereal zodiac to the sky. C: `SE_SIDM_*` constants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum SiderealMode {
    /// Fagan/Bradley ayanamsa.
    FaganBradley = 0,
    /// Lahiri (Chitrapaksha) ayanamsa.
    Lahiri = 1,
    /// De Luce ayanamsa.
    DeLuce = 2,
    /// Raman ayanamsa.
    Raman = 3,
    /// Usha/Shashi ayanamsa.
    Ushashashi = 4,
    /// Krishnamurti ayanamsa.
    Krishnamurti = 5,
    /// Djwhal Khul ayanamsa.
    DjwhalKhul = 6,
    /// Yukteshwar ayanamsa.
    Yukteshwar = 7,
    /// J.N. Bhasin ayanamsa.
    JnBhasin = 8,
    /// Babylonian/Kugler 1 ayanamsa.
    BabylKugler1 = 9,
    /// Babylonian/Kugler 2 ayanamsa.
    BabylKugler2 = 10,
    /// Babylonian/Kugler 3 ayanamsa.
    BabylKugler3 = 11,
    /// Babylonian/Huber ayanamsa.
    BabylHuber = 12,
    /// Babylonian/Eta Piscium ayanamsa.
    BabylEtpsc = 13,
    /// Babylonian ayanamsa with Aldebaran fixed at 15° Taurus.
    Aldebaran15Tau = 14,
    /// Hipparchos ayanamsa.
    Hipparchos = 15,
    /// Sassanian ayanamsa.
    Sassanian = 16,
    /// Ayanamsa with the Galactic Center fixed at 0° Sagittarius.
    GalCent0Sag = 17,
    /// J2000 ayanamsa (fixed offset from J2000 epoch).
    J2000 = 18,
    /// J1900 ayanamsa (fixed offset from J1900 epoch).
    J1900 = 19,
    /// B1950 ayanamsa (fixed offset from B1950 epoch).
    B1950 = 20,
    /// Suryasiddhanta ayanamsa.
    Suryasiddhanta = 21,
    /// Suryasiddhanta ayanamsa, mean Sun variant.
    SuryasiddhantaMsun = 22,
    /// Aryabhata ayanamsa.
    Aryabhata = 23,
    /// Aryabhata ayanamsa, mean Sun variant.
    AryabhataMsun = 24,
    /// SS (Suryasiddhanta) Revati ayanamsa.
    SsRevati = 25,
    /// SS (Suryasiddhanta) Citra ayanamsa.
    SsCitra = 26,
    /// True Citra ayanamsa.
    TrueCitra = 27,
    /// True Revati ayanamsa.
    TrueRevati = 28,
    /// True Pushya ayanamsa (PVRN Rao).
    TruePushya = 29,
    /// Galactic Center ayanamsa (Gil Brand).
    GalCentRgilbrand = 30,
    /// Galactic Equator ayanamsa (IAU 1958).
    GalEquIau1958 = 31,
    /// Galactic Equator ayanamsa (true).
    GalEquTrue = 32,
    /// Galactic Equator ayanamsa, mid-Mula variant.
    GalEquMula = 33,
    /// Skydram galactic alignment ayanamsa (Mardyks).
    GalAlignMardyks = 34,
    /// True Mula ayanamsa (Chandra Hari).
    TrueMula = 35,
    /// Dhruva/Galactic Center/Mula ayanamsa (Wilhelm).
    GalCentMulaWilhelm = 36,
    /// Aryabhata 522 ayanamsa.
    Aryabhata522 = 37,
    /// Babylonian/Britton ayanamsa.
    BabylBritton = 38,
    /// "Vedic"/Sheoran ayanamsa.
    TrueSheoran = 39,
    /// Cochrane ayanamsa (Galactic Center = 0° Capricorn).
    GalCentCochrane = 40,
    /// Galactic Equator ayanamsa (Fiorenza).
    GalEquFiorenza = 41,
    /// Vettius Valens Moon ayanamsa.
    ValensMoon = 42,
    /// Lahiri 1940 ayanamsa.
    Lahiri1940 = 43,
    /// Lahiri VP285 ayanamsa.
    LahiriVp285 = 44,
    /// Krishnamurti-Senthilathiban ayanamsa (VP291).
    KrishnamurtiVp291 = 45,
    /// Lahiri ICRC ayanamsa.
    LahiriIcrc = 46,
    /// User-defined ayanamsa (custom reference point/date supplied by the caller).
    User = 255,
}

impl SiderealMode {
    /// Returns the human-readable ayanamsa name, or `None` for the user-defined mode.
    pub fn name(self) -> Option<&'static str> {
        match self {
            Self::FaganBradley => Some("Fagan/Bradley"),
            Self::Lahiri => Some("Lahiri"),
            Self::DeLuce => Some("De Luce"),
            Self::Raman => Some("Raman"),
            Self::Ushashashi => Some("Usha/Shashi"),
            Self::Krishnamurti => Some("Krishnamurti"),
            Self::DjwhalKhul => Some("Djwhal Khul"),
            Self::Yukteshwar => Some("Yukteshwar"),
            Self::JnBhasin => Some("J.N. Bhasin"),
            Self::BabylKugler1 => Some("Babylonian/Kugler 1"),
            Self::BabylKugler2 => Some("Babylonian/Kugler 2"),
            Self::BabylKugler3 => Some("Babylonian/Kugler 3"),
            Self::BabylHuber => Some("Babylonian/Huber"),
            Self::BabylEtpsc => Some("Babylonian/Eta Piscium"),
            Self::Aldebaran15Tau => Some("Babylonian/Aldebaran = 15 Tau"),
            Self::Hipparchos => Some("Hipparchos"),
            Self::Sassanian => Some("Sassanian"),
            Self::GalCent0Sag => Some("Galact. Center = 0 Sag"),
            Self::J2000 => Some("J2000"),
            Self::J1900 => Some("J1900"),
            Self::B1950 => Some("B1950"),
            Self::Suryasiddhanta => Some("Suryasiddhanta"),
            Self::SuryasiddhantaMsun => Some("Suryasiddhanta, mean Sun"),
            Self::Aryabhata => Some("Aryabhata"),
            Self::AryabhataMsun => Some("Aryabhata, mean Sun"),
            Self::SsRevati => Some("SS Revati"),
            Self::SsCitra => Some("SS Citra"),
            Self::TrueCitra => Some("True Citra"),
            Self::TrueRevati => Some("True Revati"),
            Self::TruePushya => Some("True Pushya (PVRN Rao)"),
            Self::GalCentRgilbrand => Some("Galactic Center (Gil Brand)"),
            Self::GalEquIau1958 => Some("Galactic Equator (IAU1958)"),
            Self::GalEquTrue => Some("Galactic Equator"),
            Self::GalEquMula => Some("Galactic Equator mid-Mula"),
            Self::GalAlignMardyks => Some("Skydram (Mardyks)"),
            Self::TrueMula => Some("True Mula (Chandra Hari)"),
            Self::GalCentMulaWilhelm => Some("Dhruva/Gal.Center/Mula (Wilhelm)"),
            Self::Aryabhata522 => Some("Aryabhata 522"),
            Self::BabylBritton => Some("Babylonian/Britton"),
            Self::TrueSheoran => Some("\"Vedic\"/Sheoran"),
            Self::GalCentCochrane => Some("Cochrane (Gal.Center = 0 Cap)"),
            Self::GalEquFiorenza => Some("Galactic Equator (Fiorenza)"),
            Self::ValensMoon => Some("Vettius Valens"),
            Self::Lahiri1940 => Some("Lahiri 1940"),
            Self::LahiriVp285 => Some("Lahiri VP285"),
            Self::KrishnamurtiVp291 => Some("Krishnamurti-Senthilathiban"),
            Self::LahiriIcrc => Some("Lahiri ICRC"),
            Self::User => None,
        }
    }
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

/// Ephemeris backend selection. Determines which data source is used for planetary positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EphemerisSource {
    /// JPL Development Ephemeris (DE441 or similar). Highest accuracy, requires a `.eph` file.
    Jpl,
    /// Swiss Ephemeris compressed format (`.se1` files). Near-JPL accuracy, smaller files.
    Swiss,
    /// Moshier semi-analytical series. No files needed; ~1 arcsec accuracy for modern epochs.
    Moshier,
}

// ---------------------------------------------------------------------------
// FileDataKind / FileData
// ---------------------------------------------------------------------------

/// Which category of ephemeris file to query. Mirrors C's `ifno` parameter in
/// `swe_get_current_file_data`. In C, the function reports the file used by the
/// *last* calculation (global state). The stateless Rust equivalent
/// [`Ephemeris::file_data`] takes a `jd` to select which file would serve that
/// epoch instead.
///
/// `Asteroid` and `PlanetMoon` always return `None` in the stateless API because
/// they require knowing which specific body was queried — information that C
/// tracks implicitly via its global "last opened" file slot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum FileDataKind {
    /// Planet file (`sepl*.se1` or the JPL `.eph` file).
    Planet = 0,
    /// Moon file (`semo*.se1`).
    Moon = 1,
    /// Main asteroid file (`seas*.se1`: Chiron, Pholus, Ceres, Pallas, Juno, Vesta).
    MainAsteroid = 2,
    /// Individual numbered asteroid file. Always returns `None` (stateless: no
    /// "last-used" file concept — use the lib API for specific asteroids).
    Asteroid = 3,
    /// Planetary moon file (`sepm*.se1`). Always returns `None` (stateless: no
    /// "last-used" file concept — use the lib API for specific planet moons).
    PlanetMoon = 4,
}

impl TryFrom<i32> for FileDataKind {
    type Error = crate::Error;

    fn try_from(v: i32) -> Result<Self, Self::Error> {
        match v {
            0 => Ok(Self::Planet),
            1 => Ok(Self::Moon),
            2 => Ok(Self::MainAsteroid),
            3 => Ok(Self::Asteroid),
            4 => Ok(Self::PlanetMoon),
            _ => Err(crate::Error::InvalidBody(v)),
        }
    }
}

/// Metadata about an ephemeris file serving a given epoch.
///
/// Returned by [`Ephemeris::file_data`]. This is the stateless equivalent of C's
/// `swe_get_current_file_data` — instead of reporting the file used by the last
/// calculation, the caller provides a Julian Day to select the file.
#[derive(Debug, Clone)]
pub struct FileData {
    /// Filesystem path to the ephemeris file.
    pub path: std::path::PathBuf,
    /// Start of the file's Julian Day coverage.
    pub start_jd: f64,
    /// End of the file's Julian Day coverage.
    pub end_jd: f64,
    /// JPL DE number the file's data derives from.
    pub denum: i32,
}

// ---------------------------------------------------------------------------
// Astronomical model enums
// ---------------------------------------------------------------------------

/// Precession model selection. C: models 1–11 in `swe_set_astro_models` slots 1/2.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum PrecessionModel {
    /// IAU 1976 precession model.
    IAU1976 = 1,
    /// Laskar 1986 precession model.
    Laskar1986 = 2,
    /// Williams 1994 precession model with Laskar-derived obliquity rate.
    WillEpsLask = 3,
    /// Williams 1994 precession model.
    Williams1994 = 4,
    /// Simon 1994 precession model.
    Simon1994 = 5,
    /// IAU 2000 precession model.
    IAU2000 = 6,
    /// Bretagnon 2003 precession model.
    Bretagnon2003 = 7,
    /// IAU 2006 precession model.
    IAU2006 = 8,
    /// Vondrák 2011 long-term precession model.
    Vondrak2011 = 9,
    /// Owen 1990 precession model.
    Owen1990 = 10,
    /// Newcomb precession model.
    Newcomb = 11,
}

/// Nutation model selection. C: `swe_set_astro_models` slot 3.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum NutationModel {
    /// IAU 1980 (Wahr) nutation series, 106 terms.
    IAU1980 = 1,
    /// IAU 1980 with Herring 1987 corrections.
    IAUCorr1987 = 2,
    /// IAU 2000A full nutation model, 1365 terms.
    IAU2000A = 3,
    /// IAU 2000B abridged nutation model, 77 terms. Default.
    IAU2000B = 4,
    /// Woolard 1953 nutation model.
    Woolard = 5,
}

/// Delta T (TT − UT) model selection. C: `swe_set_astro_models` slot 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum DeltaTModel {
    /// Stephenson & Morrison 1984.
    StephensonMorrison1984 = 1,
    /// Stephenson 1997.
    Stephenson1997 = 2,
    /// Stephenson & Morrison 2004.
    StephensonMorrison2004 = 3,
    /// Espenak & Meeus 2006.
    EspenakMeeus2006 = 4,
    /// Stephenson, Morrison & Hohenkerk 2016. Default.
    StephensonEtc2016 = 5,
}

/// Sidereal time (GMST) model selection. C: `swe_set_astro_models` slot 7.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum SiderealTimeModel {
    /// IAU 1976 GMST.
    IAU1976 = 1,
    /// IAU 2006 GMST (Capitaine).
    IAU2006 = 2,
    /// IERS Conventions 2010.
    IersConv2010 = 3,
    /// Long-term model (Vondrák). Default.
    Longterm = 4,
}

/// Frame bias model (GCRS �� J2000 rotation). C: `swe_set_astro_models` slot 4.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum BiasModel {
    /// No frame bias applied.
    None = 1,
    /// IAU 2000 frame bias.
    IAU2000 = 2,
    /// IAU 2006 frame bias. Default.
    IAU2006 = 3,
}

/// JPL Horizons agreement mode. C: `swe_set_astro_models` slot 5.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum JplHorMode {
    /// Long-term agreement with JPL Horizons (dpsi/deps corrections applied).
    LongAgreement = 1,
}

/// JPL Horizons approximate agreement mode variant. C: `swe_set_astro_models` slot 6.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(i32)]
pub enum JplHoraMode {
    /// Version 1 approximation.
    V1 = 1,
    /// Version 2 approximation.
    V2 = 2,
    /// Version 3 approximation. Default.
    V3 = 3,
}

/// Collection of astronomical model overrides. Replaces C's `swe_set_astro_models`
/// (8 comma-separated values in a string). See [`AstroModels::default`] for the recommended
/// modern configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AstroModels {
    /// Delta T (TT − UT) model.
    pub delta_t: DeltaTModel,
    /// Long-term precession model (outside ±CTIES centuries of J2000).
    pub prec_longterm: PrecessionModel,
    /// Short-term precession model (within ±CTIES centuries of J2000).
    pub prec_shortterm: PrecessionModel,
    /// Nutation model.
    pub nutation: NutationModel,
    /// Frame bias model (GCRS ↔ J2000).
    pub bias: BiasModel,
    /// JPL Horizons dpsi/deps mode.
    pub jplhor_mode: JplHorMode,
    /// JPL Horizons approximate agreement variant.
    pub jplhora_mode: JplHoraMode,
    /// Sidereal time (GMST) model.
    pub sidereal_time: SiderealTimeModel,
}

/// Direction of frame-bias rotation (GCRS ↔ J2000).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FrameTransform {
    /// Rotate from J2000 dynamical frame to GCRS (kinematically non-rotating).
    J2000ToGcrs,
    /// Rotate from GCRS to J2000 dynamical frame.
    GcrsToJ2000,
}

/// Direction of precession rotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PrecessionDirection {
    /// Precess from J2000 to the ecliptic/equator of date.
    J2000ToDate,
    /// Precess from the ecliptic/equator of date back to J2000.
    DateToJ2000,
}

/// Mean or true obliquity of the ecliptic with precomputed trig values (radians).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Epsilon {
    /// Obliquity in radians.
    pub eps: f64,
    /// sin(eps).
    pub sin_eps: f64,
    /// cos(eps).
    pub cos_eps: f64,
}

impl Epsilon {
    /// Constructs an `Epsilon` from an obliquity value in radians, precomputing sin/cos.
    pub fn new(eps_rad: f64) -> Self {
        Self {
            eps: eps_rad,
            sin_eps: eps_rad.sin(),
            cos_eps: eps_rad.cos(),
        }
    }
}

/// Nutation angles (radians). `dpsi` = nutation in longitude, `deps` = nutation in obliquity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Nutation {
    /// Nutation in longitude (radians).
    pub dpsi: f64,
    /// Nutation in obliquity (radians).
    pub deps: f64,
}

// ---------------------------------------------------------------------------
// Julian Day newtypes
// ---------------------------------------------------------------------------

/// Julian Day Number on the TT (Terrestrial Time) time scale.
///
/// TT is the uniform time scale used for ephemeris calculations. It differs from UT1 by
/// Delta T (TT = UT1 + ΔT). The newtype prevents accidentally passing a UT value where TT
/// is expected.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct JdTt(
    /// Julian Day number on the TT time scale.
    pub f64,
);

/// Julian Day Number on the UT1 (Universal Time) time scale.
///
/// UT1 tracks the Earth's rotation and is the time scale of civil clocks (approximately).
/// The newtype prevents accidentally passing a TT value where UT is expected.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct JdUt1(
    /// Julian Day number on the UT1 time scale.
    pub f64,
);

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

/// Broken-down UTC date/time for `swe_utc_to_jd` / `swe_jd_to_utc` conversions.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct UtcComponents {
    /// Calendar year.
    pub year: i32,
    /// Calendar month (1–12).
    pub month: i32,
    /// Calendar day of month (1–31).
    pub day: i32,
    /// Hour of day (0–23).
    pub hour: i32,
    /// Minute of hour (0–59).
    pub minute: i32,
    /// Fractional seconds (allows leap-second representation up to 60.999...).
    pub second: f64,
}

/// Result of a UTC → Julian Day conversion, providing both time scales.
#[derive(Debug, Clone, Copy)]
pub struct UtcToJd {
    /// Julian Day on the TT time scale.
    pub tt: JdTt,
    /// Julian Day on the UT1 time scale.
    pub ut1: JdUt1,
}

// ---------------------------------------------------------------------------
// DeltaT trait
// ---------------------------------------------------------------------------

/// Trait for types that can supply a Delta T value (TT − UT1, in days) at a given UT instant.
pub trait DeltaT {
    /// Returns Delta T in days for the given Julian Day (UT1).
    fn delta_t(&self, jd_ut: JdUt1) -> f64;
}

// ---------------------------------------------------------------------------
// DegreeParts — result of swe_split_deg
// ---------------------------------------------------------------------------

/// Decomposed degree value from [`split_degrees`](crate::math::split_degrees).
/// Port of `swe_split_deg` output.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DegreeParts {
    /// Whole degrees (or zodiacal sign index when `ZODIACAL` flag is set).
    pub degrees: i32,
    /// Arc-minutes (0–59).
    pub minutes: i32,
    /// Arc-seconds (0–59).
    pub seconds: i32,
    /// Fractional arc-seconds remainder.
    pub second_fraction: f64,
    /// Sign indicator: 0 = positive, 1 = negative (or zodiacal sign number).
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
    fn body_constructors_validate() {
        assert!(Body::fictitious(40).is_ok());
        assert!(Body::fictitious(23).is_err());
        assert!(Body::asteroid(0).is_ok());
        assert!(Body::asteroid(-10000).is_err());
        assert!(Body::planet_moon(0).is_ok());
        assert!(Body::planet_moon(1000).is_err());
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
        for raw in [40, 58, 500, 999, 9000, 9500, 9999, 10000, 20000] {
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

    #[test]
    fn body_try_from_comet_range_rejected() {
        assert!(Body::try_from(1000).is_err());
        assert!(Body::try_from(5000).is_err());
        assert!(Body::try_from(8999).is_err());
    }
}
