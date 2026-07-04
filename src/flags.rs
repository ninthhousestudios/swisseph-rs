use bitflags::bitflags;

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct CalcFlags: u32 {
        const JPLEPH        = 1;
        const SWIEPH        = 2;
        const MOSEPH        = 4;
        const HELCTR        = 8;
        const TRUEPOS       = 16;
        const J2000         = 32;
        const NONUT         = 64;
        const SPEED3        = 128;
        const SPEED         = 256;
        const NOGDEFL       = 512;
        const NOABERR       = 1024;
        const EQUATORIAL    = 2048;
        const XYZ           = 4096;
        const RADIANS       = 8192;
        const BARYCTR       = 16384;
        const TOPOCTR       = 32768;
        const SIDEREAL      = 65536;
        const ICRS          = 131072;
        const DPSIDEPS_1980 = 262144;
        const JPLHOR_APPROX = 524288;
        const CENTER_BODY   = 1048576;

        const ASTROMETRIC   = Self::NOABERR.bits() | Self::NOGDEFL.bits();
        const DEFAULTEPH    = Self::SWIEPH.bits();
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct SiderealBits: u32 {
        const ECL_T0         = 256;
        const SSY_PLANE      = 512;
        const USER_UT        = 1024;
        const ECL_DATE       = 2048;
        const NO_PREC_OFFSET = 4096;
        const PREC_ORIG      = 8192;
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct NodeBits: u32 {
        const MEAN     = 1;
        const OSCU     = 2;
        const OSCU_BAR = 4;
        const FOPOINT  = 256;
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct EclipseFlags: u32 {
        const CENTRAL           = 1;
        const NONCENTRAL        = 2;
        const TOTAL             = 4;
        const ANNULAR           = 8;
        const PARTIAL           = 16;
        const HYBRID            = 32;
        const PENUMBRAL         = 64;
        const VISIBLE           = 128;
        const MAX_VISIBLE       = 256;
        const PARTBEG_VISIBLE   = 512;
        const TOTBEG_VISIBLE    = 1024;
        const TOTEND_VISIBLE    = 2048;
        const PARTEND_VISIBLE   = 4096;
        const PENUMBBEG_VISIBLE = 8192;
        const PENUMBEND_VISIBLE = 16384;
        /// Occultation begins during the day (swephexp.h:329). Numerically identical bit to
        /// [`Self::PENUMBBEG_VISIBLE`] -- the two flag families are mutually exclusive by call
        /// site (lunar-eclipse vs. occultation), not by bit layout; same "shared bit position,
        /// different meaning by context" pattern as [`Self::ONE_TRY`]/`SEFLG_TOPOCTR`.
        const OCC_BEG_DAYLIGHT  = 8192;
        /// Occultation ends during the day (swephexp.h:330). Numerically identical bit to
        /// [`Self::PENUMBEND_VISIBLE`]; see [`Self::OCC_BEG_DAYLIGHT`].
        const OCC_END_DAYLIGHT  = 16384;
        const ONE_TRY           = 32768;

        const ALLTYPES_SOLAR = Self::CENTRAL.bits()
            | Self::NONCENTRAL.bits()
            | Self::TOTAL.bits()
            | Self::ANNULAR.bits()
            | Self::PARTIAL.bits()
            | Self::HYBRID.bits();
        const ALLTYPES_LUNAR = Self::TOTAL.bits()
            | Self::PARTIAL.bits()
            | Self::PENUMBRAL.bits();
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct RiseSetFlags: u32 {
        const RISE              = 1;
        const SET               = 2;
        const MTRANSIT          = 4;
        const ITRANSIT          = 8;
        const DISC_CENTER       = 256;
        const DISC_BOTTOM       = 8192;
        const NO_REFRACTION     = 512;
        const CIVIL_TWILIGHT    = 1024;
        const NAUTIC_TWILIGHT   = 2048;
        const ASTRO_TWILIGHT    = 4096;
        const FIXED_DISC_SIZE   = 16384;
        const FORCE_SLOW        = 32768;
        const GEOCTR_NO_ECL_LAT = 128;

        const HINDU_RISING = Self::DISC_CENTER.bits()
            | Self::NO_REFRACTION.bits()
            | Self::GEOCTR_NO_ECL_LAT.bits();
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct SplitDegFlags: u32 {
        const ROUND_SEC  = 1;
        const ROUND_MIN  = 2;
        const ROUND_DEG  = 4;
        const ZODIACAL   = 8;
        const KEEP_SIGN  = 16;
        const KEEP_DEG   = 32;
        const NAKSHATRA  = 1024;
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct HeliacalFlags: u32 {
        const LONG_SEARCH     = 128;
        const HIGH_PRECISION  = 256;
        const OPTICAL_PARAMS  = 512;
        const NO_DETAILS      = 1024;
        const SEARCH_1_PERIOD = 2048;
        const VISLIM_DARK     = 4096;
        const VISLIM_NOMOON   = 8192;
        const VISLIM_PHOTOPIC = 16384;
        const VISLIM_SCOTOPIC = 32768;
        const AV              = 65536;
        const AVKIND_VR       = 65536;
        const AVKIND_PTO      = 131072;
        const AVKIND_MIN7     = 262144;
        const AVKIND_MIN9     = 524288;

        const AVKIND = Self::AVKIND_VR.bits()
            | Self::AVKIND_PTO.bits()
            | Self::AVKIND_MIN7.bits()
            | Self::AVKIND_MIN9.bits();
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct VisLimFlags: u32 {
        const SCOTOPIC = 1;
        const MIXED    = 2;
    }
}
