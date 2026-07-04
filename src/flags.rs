//! Bitflag structs for calculation options, eclipse classification, rise/set
//! search, heliacal visibility, and degree-formatting control.

use bitflags::bitflags;

bitflags! {
    /// Calculation control flags passed to [`Ephemeris::calc`](crate::Ephemeris::calc) and
    /// related methods. Selects ephemeris source, reference frame, corrections, and output format.
    /// Corresponds to the C `SEFLG_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct CalcFlags: u32 {
        /// Use JPL Development Ephemeris (DE441). C: `SEFLG_JPLEPH`.
        const JPLEPH        = 1;
        /// Use Swiss Ephemeris (.se1 files). C: `SEFLG_SWIEPH`.
        const SWIEPH        = 2;
        /// Use Moshier analytical ephemeris (no files needed). C: `SEFLG_MOSEPH`.
        const MOSEPH        = 4;
        /// Heliocentric positions. C: `SEFLG_HELCTR`. Forces `NOABERR | NOGDEFL`.
        const HELCTR        = 8;
        /// Geometric position (no light-time correction). C: `SEFLG_TRUEPOS`.
        const TRUEPOS       = 16;
        /// Output in J2000 equatorial frame (no precession to date). C: `SEFLG_J2000`.
        const J2000         = 32;
        /// Suppress nutation. C: `SEFLG_NONUT`.
        const NONUT         = 64;
        /// Compute speed via 3-sample numerical differentiation. C: `SEFLG_SPEED3`.
        /// Auto-set when `SPEED` + `TOPOCTR` + !`NOABERR` (via `plaus_iflag`).
        const SPEED3        = 128;
        /// Compute speed (daily motion). C: `SEFLG_SPEED`. Always use this for speed output.
        const SPEED         = 256;
        /// Suppress gravitational deflection of light. C: `SEFLG_NOGDEFL`.
        const NOGDEFL       = 512;
        /// Suppress annual aberration. C: `SEFLG_NOABERR`.
        const NOABERR       = 1024;
        /// Output in right ascension / declination (equatorial). C: `SEFLG_EQUATORIAL`.
        const EQUATORIAL    = 2048;
        /// Output Cartesian (x,y,z) instead of polar (lon,lat,dist). C: `SEFLG_XYZ`.
        const XYZ           = 4096;
        /// Output in radians instead of degrees. C: `SEFLG_RADIANS`.
        const RADIANS       = 8192;
        /// Barycentric (solar-system barycenter) positions. C: `SEFLG_BARYCTR`.
        /// Swiss/JPL only; Moshier rejects. Forces `NOABERR | NOGDEFL`.
        const BARYCTR       = 16384;
        /// Topocentric positions (requires `EphemerisConfig::topographic`). C: `SEFLG_TOPOCTR`.
        /// Note: bit-aliased as `SEFLG_ORBEL_AA` in orbital-element context.
        const TOPOCTR       = 32768;
        /// Sidereal zodiac (requires `EphemerisConfig::sidereal_mode`). C: `SEFLG_SIDEREAL`.
        const SIDEREAL      = 65536;
        /// ICRS (International Celestial Reference System) frame — no frame bias. C: `SEFLG_ICRS`.
        const ICRS          = 131072;
        /// Use IAU 1980 dpsi/deps. C: `SEFLG_JPLHOR` (name differs for clarity).
        const DPSIDEPS_1980 = 262144;
        /// Approximate JPL Horizons agreement mode. C: `SEFLG_JPLHOR_APPROX`.
        const JPLHOR_APPROX = 524288;
        /// Request center-of-body planetary moon for Jupiter–Pluto. C: `SEFLG_CENTER_BODY`.
        const CENTER_BODY   = 1048576;

        /// Combined `NOABERR | NOGDEFL` — pure geometric (astrometric) position.
        const ASTROMETRIC   = Self::NOABERR.bits() | Self::NOGDEFL.bits();
        /// Default ephemeris source when none is specified: Swiss Ephemeris.
        const DEFAULTEPH    = Self::SWIEPH.bits();
    }
}

bitflags! {
    /// Sidereal mode modifier bits, OR'd into the sidereal mode selection.
    /// Corresponds to C `SE_SIDBIT_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct SiderealBits: u32 {
        /// Project onto the ecliptic of t0 (rigorous geometric sidereal). C: `SE_SIDBIT_ECL_T0`.
        const ECL_T0         = 256;
        /// Project onto the solar-system plane (rigorous geometric sidereal). C: `SE_SIDBIT_SSY_PLANE`.
        const SSY_PLANE      = 512;
        /// Interpret user-supplied t0 as UT, not TT. C: `SE_SIDBIT_USER_UT`.
        const USER_UT        = 1024;
        /// Use the ecliptic of date (standard mode). C: `SE_SIDBIT_ECL_DATE`.
        const ECL_DATE       = 2048;
        /// Do not add the precession offset to the ayanamsa. C: `SE_SIDBIT_NO_PREC_OFFSET`.
        const NO_PREC_OFFSET = 4096;
        /// Use the original precession model for the ayanamsa. C: `SE_SIDBIT_PREC_ORIG`.
        const PREC_ORIG      = 8192;
    }
}

bitflags! {
    /// Method selection for [`Ephemeris::nod_aps`](crate::Ephemeris::nod_aps).
    /// Corresponds to C `SE_NODBIT_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct NodeBits: u32 {
        /// Mean orbital elements (analytic). C: `SE_NODBIT_MEAN`.
        const MEAN     = 1;
        /// Osculating (instantaneous Keplerian) elements. C: `SE_NODBIT_OSCU`.
        const OSCU     = 2;
        /// Osculating, barycentric (for outer planets beyond 6 AU). C: `SE_NODBIT_OSCU_BAR`.
        const OSCU_BAR = 4;
        /// Return the second focal point instead of the aphelion. C: `SE_NODBIT_FOPOINT`.
        const FOPOINT  = 256;
    }
}

bitflags! {
    /// Eclipse/occultation type and visibility flags, used both as input filters (`ifltype`) and
    /// return values. Corresponds to C `SE_ECL_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct EclipseFlags: u32 {
        /// Central eclipse (shadow axis touches Earth surface). C: `SE_ECL_CENTRAL`.
        const CENTRAL           = 1;
        /// Non-central eclipse (shadow axis misses Earth). C: `SE_ECL_NONCENTRAL`.
        const NONCENTRAL        = 2;
        /// Total eclipse. C: `SE_ECL_TOTAL`.
        const TOTAL             = 4;
        /// Annular eclipse. C: `SE_ECL_ANNULAR`.
        const ANNULAR           = 8;
        /// Partial eclipse. C: `SE_ECL_PARTIAL`.
        const PARTIAL           = 16;
        /// Annular-total (hybrid) eclipse. C: `SE_ECL_ANNULAR_TOTAL`.
        const HYBRID            = 32;
        /// Penumbral lunar eclipse. C: `SE_ECL_PENUMBRAL`.
        const PENUMBRAL         = 64;
        /// Eclipse visible from the observer. C: `SE_ECL_VISIBLE`.
        const VISIBLE           = 128;
        /// Maximum of eclipse is visible. C: `SE_ECL_MAX_VISIBLE`.
        const MAX_VISIBLE       = 256;
        /// Beginning of partial phase visible. C: `SE_ECL_PARTBEG_VISIBLE`.
        const PARTBEG_VISIBLE   = 512;
        /// Beginning of totality visible. C: `SE_ECL_TOTBEG_VISIBLE`.
        const TOTBEG_VISIBLE    = 1024;
        /// End of totality visible. C: `SE_ECL_TOTEND_VISIBLE`.
        const TOTEND_VISIBLE    = 2048;
        /// End of partial phase visible. C: `SE_ECL_PARTEND_VISIBLE`.
        const PARTEND_VISIBLE   = 4096;
        /// Beginning of penumbral phase visible. C: `SE_ECL_PENUMBBEG_VISIBLE`.
        const PENUMBBEG_VISIBLE = 8192;
        /// End of penumbral phase visible. C: `SE_ECL_PENUMBEND_VISIBLE`.
        const PENUMBEND_VISIBLE = 16384;
        /// Occultation begins during the day (swephexp.h:329). Numerically identical bit to
        /// [`Self::PENUMBBEG_VISIBLE`] -- the two flag families are mutually exclusive by call
        /// site (lunar-eclipse vs. occultation), not by bit layout; same "shared bit position,
        /// different meaning by context" pattern as [`Self::ONE_TRY`]/`SEFLG_TOPOCTR`.
        const OCC_BEG_DAYLIGHT  = 8192;
        /// Occultation ends during the day (swephexp.h:330). Numerically identical bit to
        /// [`Self::PENUMBEND_VISIBLE`]; see [`Self::OCC_BEG_DAYLIGHT`].
        const OCC_END_DAYLIGHT  = 16384;
        /// Search only one lunation/apparition instead of continuing until a match is found.
        const ONE_TRY           = 32768;

        /// Mask for all solar eclipse type bits.
        const ALLTYPES_SOLAR = Self::CENTRAL.bits()
            | Self::NONCENTRAL.bits()
            | Self::TOTAL.bits()
            | Self::ANNULAR.bits()
            | Self::PARTIAL.bits()
            | Self::HYBRID.bits();
        /// Mask for all lunar eclipse type bits.
        const ALLTYPES_LUNAR = Self::TOTAL.bits()
            | Self::PARTIAL.bits()
            | Self::PENUMBRAL.bits();
    }
}

bitflags! {
    /// Rise/set/transit event selection and algorithm modifiers.
    /// Corresponds to C `SE_CALC_RISE`, `SE_CALC_SET`, `SE_BIT_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct RiseSetFlags: u32 {
        /// Find rising time. C: `SE_CALC_RISE`.
        const RISE              = 1;
        /// Find setting time. C: `SE_CALC_SET`.
        const SET               = 2;
        /// Find upper meridian transit. C: `SE_CALC_MTRANSIT`.
        const MTRANSIT          = 4;
        /// Find lower meridian transit (anti-culmination). C: `SE_CALC_ITRANSIT`.
        const ITRANSIT          = 8;
        /// Use disc center (not upper limb). C: `SE_BIT_DISC_CENTER`.
        const DISC_CENTER       = 256;
        /// Use lower limb of disc. C: `SE_BIT_DISC_BOTTOM`.
        const DISC_BOTTOM       = 8192;
        /// Suppress atmospheric refraction. C: `SE_BIT_NO_REFRACTION`.
        const NO_REFRACTION     = 512;
        /// Civil twilight (Sun 6 deg below horizon). C: `SE_BIT_CIVIL_TWILIGHT`.
        const CIVIL_TWILIGHT    = 1024;
        /// Nautical twilight (Sun 12 deg below horizon). C: `SE_BIT_NAUTIC_TWILIGHT`.
        const NAUTIC_TWILIGHT   = 2048;
        /// Astronomical twilight (Sun 18 deg below horizon). C: `SE_BIT_ASTRO_TWILIGHT`.
        const ASTRO_TWILIGHT    = 4096;
        /// Use a fixed (mean) disc size instead of the apparent one. C: `SE_BIT_FIXED_DISC_SIZE`.
        const FIXED_DISC_SIZE   = 16384;
        /// Force the slow (full-precision) algorithm, bypassing the fast path. C: `SE_BIT_FORCE_SLOW_METHOD`.
        const FORCE_SLOW        = 32768;
        /// Geocentric calculation ignoring ecliptic latitude. C: `SE_BIT_GEOCTR_NO_ECL_LAT`.
        const GEOCTR_NO_ECL_LAT = 128;

        /// Hindu rising: disc center, no refraction, no ecliptic latitude. C: `SE_BIT_HINDU_RISING`.
        const HINDU_RISING = Self::DISC_CENTER.bits()
            | Self::NO_REFRACTION.bits()
            | Self::GEOCTR_NO_ECL_LAT.bits();
    }
}

bitflags! {
    /// Flags for [`split_degrees`](crate::math::split_degrees) output formatting.
    /// Corresponds to C `SE_SPLIT_DEG_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct SplitDegFlags: u32 {
        /// Round to nearest arc-second. C: `SE_SPLIT_DEG_ROUND_SEC`.
        const ROUND_SEC  = 1;
        /// Round to nearest arc-minute. C: `SE_SPLIT_DEG_ROUND_MIN`.
        const ROUND_MIN  = 2;
        /// Round to nearest degree. C: `SE_SPLIT_DEG_ROUND_DEG`.
        const ROUND_DEG  = 4;
        /// Output as zodiacal sign + degree within sign (0-29). C: `SE_SPLIT_DEG_ZODIACAL`.
        const ZODIACAL   = 8;
        /// Preserve the sign for negative input. C: `SE_SPLIT_DEG_KEEP_SIGN`.
        const KEEP_SIGN  = 16;
        /// Preserve degrees > 360 (no normalization). C: `SE_SPLIT_DEG_KEEP_DEG`.
        const KEEP_DEG   = 32;
        /// Divide into 27 nakshatras instead of 12 signs. C: `SE_SPLIT_DEG_NAKSHATRA`.
        const NAKSHATRA  = 1024;
    }
}

bitflags! {
    /// Heliacal visibility calculation flags. Corresponds to C `SE_HELFLAG_*` constants.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct HeliacalFlags: u32 {
        /// Extend the synodic-period retry search up to 1 million days. C: `SE_HELFLAG_LONG_SEARCH`.
        const LONG_SEARCH     = 128;
        /// Use higher-precision iteration. C: `SE_HELFLAG_HIGH_PRECISION`.
        const HIGH_PRECISION  = 256;
        /// Interpret `dobs` as optical instrument parameters. C: `SE_HELFLAG_OPTICAL_PARAMS`.
        const OPTICAL_PARAMS  = 512;
        /// Suppress detailed output fields. C: `SE_HELFLAG_NO_DETAILS`.
        const NO_DETAILS      = 1024;
        /// Reject results more than one synodic period from start. C: `SE_HELFLAG_SEARCH_1_PERIOD`.
        const SEARCH_1_PERIOD = 2048;
        /// Dark-adapted eye (scotopic threshold for vis-limit). C: `SE_HELFLAG_VISLIM_DARK`.
        const VISLIM_DARK     = 4096;
        /// Ignore Moon interference in the limiting-magnitude calculation. C: `SE_HELFLAG_VISLIM_NOMOON`.
        const VISLIM_NOMOON   = 8192;
        /// Force photopic (daylight) vision model. C: `SE_HELFLAG_VISLIM_PHOTOPIC`.
        const VISLIM_PHOTOPIC = 16384;
        /// Force scotopic (night) vision model. C: `SE_HELFLAG_VISLIM_SCOTOPIC`.
        const VISLIM_SCOTOPIC = 32768;
        /// Use arc-visibility (arcus visionis) method instead of vis-limit. C: `SE_HELFLAG_AV`.
        const AV              = 65536;
        /// Arc-vis variant: visibility ratio. C: `SE_HELFLAG_AVKIND_VR`.
        const AVKIND_VR       = 65536;
        /// Arc-vis variant: Ptolemy. C: `SE_HELFLAG_AVKIND_PTO`.
        const AVKIND_PTO      = 131072;
        /// Arc-vis variant: 7th-magnitude fixed depth. C: `SE_HELFLAG_AVKIND_MIN7`.
        const AVKIND_MIN7     = 262144;
        /// Arc-vis variant: 9th-magnitude fixed depth. C: `SE_HELFLAG_AVKIND_MIN9`.
        const AVKIND_MIN9     = 524288;

        /// Mask for all AVKIND bits.
        const AVKIND = Self::AVKIND_VR.bits()
            | Self::AVKIND_PTO.bits()
            | Self::AVKIND_MIN7.bits()
            | Self::AVKIND_MIN9.bits();
    }
}

bitflags! {
    /// Output flags from [`vis_limit_mag`](crate::Ephemeris::vis_limit_mag) indicating the vision
    /// regime used for the limiting-magnitude calculation.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct VisLimFlags: u32 {
        /// Scotopic (dark-adapted rod) vision regime was used.
        const SCOTOPIC = 1;
        /// Mixed (mesopic transition between photopic and scotopic) regime.
        const MIXED    = 2;
    }
}
