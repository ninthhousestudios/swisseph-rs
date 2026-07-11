# Swiss Ephemeris Calculation Library — Catalogue

> **Phase 1: Public API Catalogue**
> Authors: Dieter Koch and Alois Treindl, Astrodienst AG, Switzerland
> License: AGPL or Swiss Ephemeris Professional License
>
> **Phase 2 (not yet done):** For each domain section, expand the `<!-- INTERNALS TODO -->` marker with:
> - All internal `swi_*` functions in that domain — signature, purpose, key logic
> - End-to-end calculation flow traces (e.g., `swe_calc()` → `swecalc()` → `main_planet()` → `sweplan()` → file I/O → Chebyshev interp → corrections → coordinate transform → return)
> - Detailed algorithm explanations: not just "Chebyshev interpolation" but how the coefficients are stored, how segments are selected, how the recurrence works, etc.
> - Cross-domain dependencies (e.g., eclipse code calling rise/set code calling refraction code)

---

## 1. Library Overview

The Swiss Ephemeris is a C library for computing high-precision astronomical positions and events. Originally developed for astrological software, it provides sub-arcsecond accuracy for planetary, lunar, and stellar positions spanning millennia.

### Architecture

```
Application
    │
    ▼
swephexp.h ── Public API (swe_* functions)
    │
    ├── sweph.c ──────── Core dispatcher, fixed stars, ayanamsa, zodiacal crossings
    ├── swecl.c ──────── Eclipses, occultations, rise/set, phenomena, nodes/apsides
    ├── swehouse.c ───── House system calculations (24 systems)
    ├── swehel.c ──────── Heliacal events and visibility
    ├── swedate.c ────── Date/time conversions, Julian Day
    ├── swephlib.c ───── Math, precession, nutation, delta-T, coordinate transforms
    │
    └── Ephemeris backends (selected via SEFLG_* flags):
        ├── swejpl.c ────── JPL DE ephemerides (Chebyshev interpolation of tabulated data)
        ├── sweephe4.c ──── Swiss Ephemeris binary files (compressed Chebyshev coefficients)
        ├── swemplan.c ──── Moshier analytical planetary ephemeris (Fourier series from DE404)
        └── swemmoon.c ──── Moshier analytical lunar ephemeris (ELP2000-85 fit to DE404)
```

### High-Level Calculation Flow

1. **`swe_calc()`** receives a Julian Day, planet ID, and flags
2. Dispatcher selects ephemeris backend: JPL → Swiss Ephemeris → Moshier (cascade fallback)
3. Backend returns heliocentric J2000 ecliptic Cartesian coordinates
4. Corrections applied: light-time, aberration, gravitational deflection
5. Coordinate transformations: precession → nutation → ecliptic/equatorial/horizontal
6. Optional: sidereal projection (ayanamsa subtracted), topocentric parallax

---

## 2. Source File Map

| File | Role | Primary Domain |
|------|------|----------------|
| `sweph.c` | Core dispatcher, fixed stars, ayanamsa, crossings | Planetary Positions, Sidereal, Crossings |
| `swecl.c` | Eclipses, phenomena, rise/set, nodes, orbits, refraction | Eclipses, Phenomena, Rise/Set, Nodes |
| `swephlib.c` | Math utilities, precession, nutation, delta-T, sidereal time | Coordinate Systems, Time Scales, Utilities |
| `swehouse.c` | House cusp calculations (24 systems) | House Systems |
| `swehel.c` | Heliacal events, visibility limit, arcus visionis | Heliacal Events |
| `swedate.c` | Julian Day, calendar conversions, UTC/leap seconds | Time Scales |
| `swejpl.c` | JPL DE ephemeris file reader | Planetary Positions (backend) |
| `sweephe4.c` | Swiss Ephemeris binary file reader | Planetary Positions (backend) |
| `swemplan.c` | Moshier analytical planetary series | Planetary Positions (backend) |
| `swemmoon.c` | Moshier analytical lunar series | Planetary Positions (backend) |
| `swevents.c` | Event detection CLI (unsupported demo) | — |
| `swephexp.h` | Public API declarations, all constants | API |
| `sweph.h` | Internal structures, internal function declarations | Internal |
| `swephlib.h` | Internal math/precession declarations | Internal |
| `swehouse.h` | Internal house calculation declarations | Internal |
| `swedate.h` | Internal date function declarations | Internal |
| `swejpl.h` | Internal JPL reader declarations | Internal |
| `sweephe4.h` | Internal ep4 file reader declarations | Internal |
| `sweodef.h` | Platform defines, type aliases | Internal |
| `swemptab.h` | Moshier planetary series coefficients | Data |
| `swenut2000a.h` | IAU 2000A nutation series (1325 terms) | Data |
| `swevents.h` | Event detection declarations (unsupported) | — |
| `swedll.h` | Windows DLL export macros | Platform |

---

## 3. Calculation Domains

### 3.1 Planetary Positions

**What it computes:** Geocentric, heliocentric, barycentric, topocentric, or planetocentric positions of the Sun, Moon, planets (Mercury–Pluto), asteroids (600,000+), comets, Chiron, Pholus, lunar nodes/apogees, fictitious bodies (Uranian/Hamburg planets, Transpluto, etc.), and fixed stars. Returns ecliptic longitude/latitude/distance plus optional speeds.

#### Public API Functions

**`int32 swe_calc(double tjd, int ipl, int32 iflag, double *xx, char *serr)`**
Core position calculator. Takes Ephemeris Time (TT).
- `tjd`: Julian Day in TT/ET
- `ipl`: Body ID (SE_SUN=0 through SE_PLUTO=9, SE_MEAN_NODE=10, SE_TRUE_NODE=11, SE_CHIRON=15, SE_AST_OFFSET+n for asteroids, SE_COMET_OFFSET+n for comets, SE_FICT_OFFSET+n for fictitious bodies)
- `iflag`: Bit flags controlling ephemeris source, coordinate system, corrections
- `xx[6]`: Output — [longitude, latitude, distance, speed_lon, speed_lat, speed_dist] (degrees, AU, degrees/day)
- Returns: flags indicating ephemeris used; ERR (−1) on failure
- Caches results per (tjd, ipl, iflag) tuple for performance

**`int32 swe_calc_ut(double tjd_ut, int32 ipl, int32 iflag, double *xx, char *serr)`**
Same as `swe_calc()` but takes Universal Time. Internally converts UT→TT via delta-T.

**`int32 swe_calc_pctr(double tjd, int32 ipl, int32 iplctr, int32 iflag, double *xxret, char *serr)`**
Planetocentric calculator — position of `ipl` as seen from `iplctr`. Computes both bodies in barycentric J2000, subtracts, applies light-time iteration, aberration, deflection, precession, and nutation.

**`int32 swe_fixstar(char *star, double tjd, int32 iflag, double *xx, char *serr)`**
Computes fixed star position. `star` is a name ("Sirius"), Bayer designation (",alCMa"), or catalog number. Modified on return to contain full name. Reads from `sefstars.txt` catalog (~9000 stars) or built-in array of major stars. Applies proper motion, precession, parallax from J2000 epoch.

**`int32 swe_fixstar_ut(char *star, double tjd_ut, int32 iflag, double *xx, char *serr)`**
UT wrapper for `swe_fixstar()`.

**`int32 swe_fixstar_mag(char *star, double *mag, char *serr)`**
Returns visual magnitude of a fixed star.

**`int32 swe_fixstar2(char *star, double tjd, int32 iflag, double *xx, char *serr)`**
Improved version of `swe_fixstar()` with struct-based caching. Loads all stars into memory on first call for faster repeated lookups.

**`int32 swe_fixstar2_ut(char *star, double tjd_ut, int32 iflag, double *xx, char *serr)`**
UT wrapper for `swe_fixstar2()`.

**`int32 swe_fixstar2_mag(char *star, double *mag, char *serr)`**
Magnitude lookup for `swe_fixstar2()`.

**`char *swe_get_planet_name(int ipl, char *spname)`**
Returns English name for any body ID. For asteroids, reads from ephemeris file or `seasnam.txt`.

**`char *swe_version(char *s)`**
Returns library version string.

**`char *swe_get_library_path(char *s)`**
Returns absolute path to the loaded shared library.

**`const char *swe_get_current_file_data(int ifno, double *tfstart, double *tfend, int *denum)`**
Returns metadata (filename, time range, DE number) for currently open ephemeris file.

#### Key Flags (SEFLG_*)

| Flag | Value | Effect |
|------|-------|--------|
| `SEFLG_JPLEPH` | 1 | Use JPL ephemeris |
| `SEFLG_SWIEPH` | 2 | Use Swiss Ephemeris files (default) |
| `SEFLG_MOSEPH` | 4 | Use Moshier analytical ephemeris |
| `SEFLG_SPEED` | 256 | Compute high-precision speeds |
| `SEFLG_HELCTR` | 8 | Heliocentric positions |
| `SEFLG_BARYCTR` | 16384 | Barycentric positions |
| `SEFLG_TOPOCTR` | 32768 | Topocentric (requires `swe_set_topo()`) |
| `SEFLG_TRUEPOS` | 16 | Geometric position (no aberration, no light-time) |
| `SEFLG_NOABERR` | 1024 | No annual aberration |
| `SEFLG_NOGDEFL` | 512 | No gravitational deflection |
| `SEFLG_ASTROMETRIC` | 1536 | Astrometric (light-time yes, aberration/deflection no) |
| `SEFLG_J2000` | 32 | J2000 equinox (no precession) |
| `SEFLG_NONUT` | 64 | Mean equinox of date (no nutation) |
| `SEFLG_EQUATORIAL` | 2048 | Equatorial coordinates (RA/Dec) |
| `SEFLG_XYZ` | 4096 | Cartesian coordinates |
| `SEFLG_RADIANS` | 8192 | Radians instead of degrees |
| `SEFLG_SIDEREAL` | 65536 | Sidereal zodiac (requires `swe_set_sid_mode()`) |
| `SEFLG_ICRS` | 131072 | ICRS reference frame (DE406 frame) |
| `SEFLG_JPLHOR` | 262144 | Reproduce JPL Horizons (uses daily dpsi/deps) |
| `SEFLG_JPLHOR_APPROX` | 524288 | Approximate JPL Horizons |
| `SEFLG_CENTER_BODY` | 1048576 | Center-of-body, not system barycenter |

#### Ephemeris Backend Algorithms

**JPL DE Ephemerides** (`swejpl.c`): NASA's Developmental Ephemerides (DE200, DE403–DE406, DE431, DE441) store positions as **Chebyshev polynomial coefficients** over fixed time segments (~32 days). The reader normalizes time to [−1, 1] within each segment, evaluates Chebyshev polynomials via recurrence (T_{n+1} = 2t·T_n − T_{n-1}), and sums weighted coefficients. Velocity, acceleration, and jerk are obtained by differentiating the polynomial. Highest accuracy (~0.001″), requires binary data files (tens to hundreds of MB).

**Swiss Ephemeris Files** (`sweephe4.c`): Compressed binary files storing planetary longitudes in 10-day blocks as first and second differences (integers in centisecond units, ~10:1 compression). Positions are reconstructed by summing differences, then interpolated using **Everett 5th-order interpolation** (Pottenger's method). Fast and compact, accuracy <0.01″ for most bodies.

**Moshier Planetary** (`swemplan.c`): Analytical Fourier series derived by Steve Moshier from JPL DE404. Nine fundamental arguments (planetary mean longitudes from Simon et al. 1994) are combined in harmonic series with polynomial-in-time coefficients: longitude = Σ(C₀ + C₁·T + C₂·T² + ...) · sin/cos(linear combination of arguments). ~1000–2000 terms per body. No data files needed. Accuracy ~0.01″ for inner planets, ~0.1″ for outer, valid ~1800–2100 AD.

**Moshier Lunar** (`swemmoon.c`): Based on Chapront-Touzé & Chapront's **ELP2000-85** lunar theory, adjusted to fit DE404. Sums 1400+ periodic terms using four fundamental arguments (mean anomaly, elongation, distance from node, mean longitude). Each term: coefficient · sin(n₁·M + n₂·D + n₃·F + n₄·L + corrections). Accuracy ~7″ longitude, valid −3000 to +3000 AD.

#### Planet IDs

| ID | Constant | Body |
|----|----------|------|
| 0 | `SE_SUN` | Sun |
| 1 | `SE_MOON` | Moon |
| 2–9 | `SE_MERCURY`–`SE_PLUTO` | Planets |
| 10 | `SE_MEAN_NODE` | Mean lunar node |
| 11 | `SE_TRUE_NODE` | True (osculating) lunar node |
| 12 | `SE_MEAN_APOG` | Mean lunar apogee (Black Moon Lilith) |
| 13 | `SE_OSCU_APOG` | Osculating lunar apogee |
| 14 | `SE_EARTH` | Earth (for heliocentric) |
| 15 | `SE_CHIRON` | Chiron |
| 16 | `SE_PHOLUS` | Pholus |
| 17–20 | `SE_CERES`–`SE_VESTA` | Main-belt asteroids |
| 21 | `SE_INTP_APOG` | Interpolated lunar apogee |
| 22 | `SE_INTP_PERG` | Interpolated lunar perigee |
| 40–58 | `SE_CUPIDO`–`SE_WALDEMATH` | Fictitious / Uranian bodies |
| 10000+n | `SE_AST_OFFSET+n` | Asteroid by MPC number |
| 1000+n | `SE_COMET_OFFSET+n` | Comet |
| 9000+n | `SE_PLMOON_OFFSET+n` | Planetary moons |

<!-- INTERNALS TODO (§3.1): Trace swe_calc() → swecalc() dispatch logic; document swi_pleph(), sweplan(), jplplan() and the cascade fallback chain; detail Chebyshev segment selection and interpolation in swejpl.c; explain Moshier series evaluation loop in swemplan.c; trace fixed star computation from catalog parse through proper motion and precession; document caching/invalidation logic. -->

---

### 3.2 Coordinate Systems & Transformations

**What it computes:** Conversions between ecliptic and equatorial polar coordinates (with optional speed transformation), and between ecliptic/equatorial and horizontal (azimuth/altitude) coordinates including atmospheric refraction.

#### Public API Functions

**`void swe_cotrans(double *xpo, double *xpn, double eps)`**
Transforms between ecliptic and equatorial polar coordinates (position only).
- `xpo[3]`: Input [longitude, latitude, distance] in degrees
- `xpn[3]`: Output [transformed_lon, transformed_lat, distance]
- `eps`: Obliquity of ecliptic in degrees. Sign determines direction: negative = ecliptic→equatorial, positive = equatorial→ecliptic
- Algorithm: Converts to Cartesian, applies rotation matrix around X-axis by angle eps, converts back to polar.

**`void swe_cotrans_sp(double *xpo, double *xpn, double eps)`**
Same as `swe_cotrans()` but also transforms the velocity components (indices 3–5).

**`void swe_azalt(double tjd_ut, int32 calc_flag, double *geopos, double atpress, double attemp, double *xin, double *xaz)`**
Converts ecliptic or equatorial coordinates to horizontal (azimuth/altitude).
- `calc_flag`: `SE_ECL2HOR` (ecliptic→horizontal) or `SE_EQU2HOR` (equatorial→horizontal)
- `geopos[3]`: Observer [longitude, latitude, altitude_meters]
- `atpress`: Atmospheric pressure in mbar (0 = auto-estimate from altitude)
- `attemp`: Temperature in °C
- `xin[2]`: Input [longitude/RA, latitude/Dec]
- `xaz[3]`: Output [azimuth (south=0, west=90), true_altitude, apparent_altitude]
- Algorithm: Computes local sidereal time, hour angle, applies spherical trigonometry, then atmospheric refraction.

**`void swe_azalt_rev(double tjd_ut, int32 calc_flag, double *geopos, double *xin, double *xout)`**
Inverse of `swe_azalt()` — horizontal to ecliptic or equatorial.
- `calc_flag`: `SE_HOR2ECL` or `SE_HOR2EQU`

#### Key Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `SE_ECL2HOR` | 0 | Ecliptic to horizontal |
| `SE_EQU2HOR` | 1 | Equatorial to horizontal |
| `SE_HOR2ECL` | 0 | Horizontal to ecliptic |
| `SE_HOR2EQU` | 1 | Horizontal to equatorial |

<!-- INTERNALS TODO (§3.2): Document swi_coortrf(), swi_cartpol(), swi_polcart() coordinate conversion internals; explain the rotation matrix used in ecliptic↔equatorial; trace swe_azalt() through sidereal time → hour angle → spherical trig → refraction. -->

---

### 3.3 Precession, Nutation & Frame Bias

**What it computes:** The slow drift of Earth's rotational axis (precession), short-period wobbles superimposed on that drift (nutation), and the small constant offset between the dynamical equinox and the ICRS origin (frame bias). These corrections transform between J2000 and the equinox/ecliptic of any date.

#### Precession Models (11)

Selected via `swe_set_astro_models()` with `SE_MODEL_PREC_LONGTERM` and `SE_MODEL_PREC_SHORTTERM`.

| # | Constant | Model | Description |
|---|----------|-------|-------------|
| 1 | `SEMOD_PREC_IAU_1976` | IAU 1976 (Lieske) | Classic rigid-Earth model; polynomial in time from J2000. Standard from 1976–2000. |
| 2 | `SEMOD_PREC_LASKAR_1986` | Laskar 1986 | Long-term analytical theory using secular perturbation theory. Good for ±10,000 years from J2000. Polynomial coefficients from numerical integration of planetary system. |
| 3 | `SEMOD_PREC_WILL_EPS_LASK` | Williams 1994 + Laskar ε | Hybrid: Williams precession rates with Laskar obliquity. Improved accuracy for recent centuries while maintaining long-term validity. |
| 4 | `SEMOD_PREC_WILLIAMS_1994` | Williams 1994 | Updated precession constants from modern VLBI and LLR observations. Better agreement with contemporary data than IAU 1976. |
| 5 | `SEMOD_PREC_SIMON_1994` | Simon 1994 | From "Precession formulae and mean elements for Moon and Planets" (A&A 282, p.675). Extremely long-term analytical expressions derived from planetary perturbation theory. |
| 6 | `SEMOD_PREC_IAU_2000` | IAU 2000 | Lieske 1976 framework with Mathews 2002 corrections from MHB2000 nutation model. Intermediate step between old and new IAU standards. |
| 7 | `SEMOD_PREC_BRETAGNON_2003` | Bretagnon 2003 | From VSOP2000 analytical planetary theory. High-precision for historical and future epochs, accounts for all major planetary perturbations. |
| 8 | `SEMOD_PREC_IAU_2006` | IAU 2006 (Capitaine) | Current IAU standard (P03 precession). Based on SOFA algorithms with improved accuracy for ~1500–2500 AD. Consistent with IAU 2000A nutation. |
| 9 | `SEMOD_PREC_VONDRAK_2011` | Vondrák 2011 **(default)** | Combines Bretagnon analytical terms with modern VLBI data. Best model for all epochs; smooth long-term behavior with high short-term accuracy. |
| 10 | `SEMOD_PREC_OWEN_1990` | Owen 1990 | Used by JPL Horizons before 1799 and after 2202. Bridges ancient and far-future calculations. |
| 11 | `SEMOD_PREC_NEWCOMB` | Newcomb 1895 | Historical model for continuity with 19th-century astronomical tables. |

#### Nutation Models (5)

Selected via `swe_set_astro_models()` with `SE_MODEL_NUT`.

| # | Constant | Model | Description |
|---|----------|-------|-------------|
| 1 | `SEMOD_NUT_IAU_1980` | IAU 1980 (Wahr) | 106 periodic terms modeling forced nutations of a rigid Earth modified for elasticity. Precision ~10 milliarcseconds. |
| 2 | `SEMOD_NUT_IAU_CORR_1987` | IAU 1980 + Herring 1987 | Herring's corrections to IAU 1980 accounting for ocean tide effects and mantle anelasticity. Modest improvement. |
| 3 | `SEMOD_NUT_IAU_2000A` | IAU 2000A | Full-precision model with 1325 periodic terms (678 luni-solar + 687 planetary). Accounts for ocean tides, atmospheric loading, and mantle inelasticity. Precision ~0.1 milliarcseconds. Very computationally expensive. |
| 4 | `SEMOD_NUT_IAU_2000B` | IAU 2000B **(default)** | Truncated IAU 2000A retaining only the ~77 largest-amplitude terms. Precision ~1 milliarcsecond — sufficient for virtually all applications. Much faster than 2000A. |
| 5 | `SEMOD_NUT_WOOLARD` | Woolard 1953 | Historical model; simplest series. Rarely used except for reproducing old calculations. |

#### Frame Bias Models (3)

| # | Constant | Description |
|---|----------|-------------|
| 1 | `SEMOD_BIAS_NONE` | Ignore frame bias |
| 2 | `SEMOD_BIAS_IAU2000` | IAU 2000 bias matrix |
| 3 | `SEMOD_BIAS_IAU2006` | IAU 2006 bias matrix **(default)** |

#### Model Selection Functions

**`void swe_set_astro_models(char *samod, int32 iflag)`**
Selects precession, nutation, delta-T, sidereal time, and frame bias models.
- `samod`: Comma-separated model numbers "D,PL,PS,N,B,J,JA,S" or a version string like "SE2.06" to load historical presets
- `iflag`: Ephemeris flag

**`void swe_get_astro_models(char *samod, char *sdet, int32 iflag)`**
Retrieves current model settings and generates a detailed text description of all active models.

<!-- INTERNALS TODO (§3.3): Document swi_precess(), swi_nutation(), swi_epsiln() internals; trace how each precession model computes psi_A, omega_A, epsilon_A; explain the 1325-term IAU 2000A series evaluation in swenut2000a.h; document frame bias matrix application; explain how models are switched at runtime. -->

---

### 3.4 Time Scales

**What it computes:** Conversions between calendar dates and Julian Day Numbers, between Universal Time (UT) and Terrestrial Time (TT) via Delta-T, sidereal time, equation of time, and local time conversions (LMT↔LAT, UTC↔timezone).

#### Public API Functions

**`double swe_julday(int year, int month, int day, double hour, int gregflag)`**
Converts calendar date/time to Julian Day Number.
- `gregflag`: `SE_GREG_CAL` (Gregorian, post-1582-Oct-15) or `SE_JUL_CAL` (Julian)
- Note: midnight = JD.5, noon = JD.0 (astronomical convention)

**`void swe_revjul(double jd, int gregflag, int *jyear, int *jmon, int *jday, double *jut)`**
Inverse of `swe_julday()`. JD → calendar date + decimal hour.

**`int swe_date_conversion(int y, int m, int d, double utime, char c, double *tjd)`**
Like `swe_julday()` but validates the date. Returns ERR if illegal (e.g., Feb 30).

**`int32 swe_utc_to_jd(int32 iyear, ..., double dsec, int32 gregflag, double *dret, char *serr)`**
Converts UTC to Julian Day, properly handling leap seconds (from 1972 onward).
- `dret[0]`: JD in Terrestrial Time (TT) — for use with `swe_calc()`
- `dret[1]`: JD in UT1 — for use with `swe_calc_ut()`
- Before 1972: input treated as UT1 (UTC did not exist)

**`void swe_jdet_to_utc(double tjd_et, int32 gregflag, int32 *iyear, ..., double *dsec)`**
Converts TT Julian Day back to UTC calendar date. Can output `dsec=60.0` for leap seconds.

**`void swe_jdut1_to_utc(double tjd_ut, int32 gregflag, int32 *iyear, ..., double *dsec)`**
Converts UT1 Julian Day to UTC calendar date.

**`void swe_utc_time_zone(int32 iyear, ..., double dsec, double d_timezone, int32 *iyear_out, ..., double *dsec_out)`**
Applies timezone offset. Handles day-boundary crossings and leap seconds.

**`double swe_deltat(double tjd)`**
Returns Delta-T (TT − UT) in days for a given Julian Day. Uses default ephemeris.

**`double swe_deltat_ex(double tjd, int32 iflag, char *serr)`**
Full-control Delta-T with explicit ephemeris flag. The ephemeris matters because different DE ephemerides assume different lunar tidal acceleration values, which affects the parabolic extrapolation of Delta-T into the past.

**`double swe_sidtime(double tjd_ut)`**
Returns Greenwich Apparent Sidereal Time (GAST) in hours. Automatically computes obliquity and nutation.

**`double swe_sidtime0(double tjd_ut, double eps, double nut)`**
Returns GAST given pre-computed obliquity and nutation (both in degrees).

**`int32 swe_time_equ(double tjd_ut, double *E, char *serr)`**
Equation of time: difference between apparent and mean solar time. Output `E` is in days.

**`int32 swe_lmt_to_lat(double tjd_lmt, double geolon, double *tjd_lat, char *serr)`**
Converts Local Mean Time to Local Apparent Time (adds equation of time).

**`int32 swe_lat_to_lmt(double tjd_lat, double geolon, double *tjd_lmt, char *serr)`**
Converts Local Apparent Time to Local Mean Time (subtracts equation of time, iterative).

**`int swe_day_of_week(double jd)`**
Returns day of week: 0=Monday, 6=Sunday.

#### Delta-T Models (5)

| # | Constant | Range | Method |
|---|----------|-------|--------|
| 1 | `SEMOD_DELTAT_STEPHENSON_MORRISON_1984` | Before 1620 | Polynomial fit to historical eclipse observations. Borkowski formula before 948 AD. |
| 2 | `SEMOD_DELTAT_STEPHENSON_1997` | Before 1600 | Tabulated at 50-year intervals (−500 to 1600) with linear interpolation. Polynomial before −500. |
| 3 | `SEMOD_DELTAT_STEPHENSON_MORRISON_2004` | Before 1600 | Tabulated at 100-year intervals (−1000 to 1600). Parabolic long-term formula before −1000. |
| 4 | `SEMOD_DELTAT_ESPENAK_MEEUS_2006` | Before 1633 | Seven piecewise polynomials with time-varying coefficients for different historical epochs. |
| 5 | `SEMOD_DELTAT_STEPHENSON_ETC_2016` | **(default)** | Cubic spline fit to eclipse/occultation data; polynomial for 2000–2500; parabolic after 2500. Smoothly connected. |

For 1620–present, all models use the same source: tabulated values from the Astronomical Almanac and IERS, interpolated via Bessel's 4th-order method.

All pre-1955 values are adjusted for tidal acceleration differences: `correction = −0.000091 × (tid_acc − tid_acc₀) × (Year − 1955)²`

#### Sidereal Time Models (4)

| # | Constant | Description |
|---|----------|-------------|
| 1 | `SEMOD_SIDT_IAU_1976` | Classic polynomial formula |
| 2 | `SEMOD_SIDT_IAU_2006` | Capitaine 2003 coefficient-based |
| 3 | `SEMOD_SIDT_IERS_CONV_2010` | ERA-based (Earth Rotation Angle) from IERS 2010 |
| 4 | `SEMOD_SIDT_LONGTERM` | **(default)** IERS 2010 for 1850–2050; Simon et al. mean Earth formula outside this range |

<!-- INTERNALS TODO (§3.4): Document Bessel interpolation of tabulated delta-T values; trace how swe_deltat_ex() selects model and applies tidal correction; explain sidereal time polynomial evaluation and the long-term extension; document leap second table and its update mechanism. -->

---

### 3.5 Sidereal/Tropical & Ayanamsa

**What it computes:** The ayanamsa — the angular offset between the tropical zodiac (tied to the vernal equinox, which precesses) and a sidereal zodiac (tied to fixed stars). With 47 predefined ayanamsas plus user-defined, the library can express positions in any sidereal reference frame. Projection modes control whether positions are projected onto the ecliptic of t0, the solar system plane, or the ecliptic of date.

#### Public API Functions

**`void swe_set_sid_mode(int32 sid_mode, double t0, double ayan_t0)`**
Activates sidereal mode and selects which ayanamsa to use.
- `sid_mode`: Predefined ID (0–46) or `SE_SIDM_USER` (255) for custom
- `t0`: Reference epoch for custom ayanamsa (ignored for predefined)
- `ayan_t0`: Ayanamsa value at t0 in degrees (ignored for predefined)
- Can be OR'd with projection bits: `SE_SIDBIT_ECL_T0`, `SE_SIDBIT_SSY_PLANE`, `SE_SIDBIT_ECL_DATE`

**`int32 swe_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr)`**
Returns the ayanamsa at a given time with full ephemeris control.

**`int32 swe_get_ayanamsa_ex_ut(double tjd_ut, int32 iflag, double *daya, char *serr)`**
UT wrapper.

**`double swe_get_ayanamsa(double tjd_et)`**
Simplified ayanamsa (no nutation, no error checking).

**`double swe_get_ayanamsa_ut(double tjd_ut)`**
UT simplified ayanamsa.

**`const char *swe_get_ayanamsa_name(int32 isidmode)`**
Returns the name of a predefined ayanamsa mode.

#### Predefined Ayanamsas (47)

| ID | Constant | Name |
|----|----------|------|
| 0 | `SE_SIDM_FAGAN_BRADLEY` | Fagan/Bradley |
| 1 | `SE_SIDM_LAHIRI` | Lahiri |
| 2 | `SE_SIDM_DELUCE` | De Luce |
| 3 | `SE_SIDM_RAMAN` | Raman |
| 4 | `SE_SIDM_USHASHASHI` | Ushashashi |
| 5 | `SE_SIDM_KRISHNAMURTI` | Krishnamurti |
| 6 | `SE_SIDM_DJWHAL_KHUL` | Djwhal Khul |
| 7 | `SE_SIDM_YUKTESHWAR` | Yukteshwar |
| 8 | `SE_SIDM_JN_BHASIN` | J.N. Bhasin |
| 9–13 | `SE_SIDM_BABYL_KUGLER1`–`SE_SIDM_BABYL_ETPSC` | Babylonian (various) |
| 14 | `SE_SIDM_ALDEBARAN_15TAU` | Aldebaran at 15° Taurus |
| 15 | `SE_SIDM_HIPPARCHOS` | Hipparchos |
| 16 | `SE_SIDM_SASSANIAN` | Sassanian |
| 17 | `SE_SIDM_GALCENT_0SAG` | Galactic Center at 0° Sagittarius |
| 18–20 | `SE_SIDM_J2000`–`SE_SIDM_B1950` | Fixed equinox (J2000, J1900, B1950) |
| 21–24 | `SE_SIDM_SURYASIDDHANTA`–`SE_SIDM_ARYABHATA_MSUN` | Indian classical |
| 25–26 | `SE_SIDM_SS_REVATI`–`SE_SIDM_SS_CITRA` | Suryasiddhanta star-based |
| 27–29 | `SE_SIDM_TRUE_CITRA`–`SE_SIDM_TRUE_PUSHYA` | True star (Spica/Revati/Pushya) |
| 30–34 | `SE_SIDM_GALCENT_RGILBRAND`–`SE_SIDM_GALALIGN_MARDYKS` | Galactic alignment |
| 35 | `SE_SIDM_TRUE_MULA` | True Mula |
| 36 | `SE_SIDM_GALCENT_MULA_WILHELM` | Galactic Center Mula (Wilhelm) |
| 37 | `SE_SIDM_ARYABHATA_522` | Aryabhata 522 |
| 38 | `SE_SIDM_BABYL_BRITTON` | Babylonian (Britton) |
| 39 | `SE_SIDM_TRUE_SHEORAN` | True Sheoran |
| 40 | `SE_SIDM_GALCENT_COCHRANE` | Galactic Center (Cochrane) |
| 41 | `SE_SIDM_GALEQU_FIORENZA` | Galactic Equator (Fiorenza) |
| 42 | `SE_SIDM_VALENS_MOON` | Valens Moon |
| 43 | `SE_SIDM_LAHIRI_1940` | Lahiri 1940 |
| 44 | `SE_SIDM_LAHIRI_VP285` | Lahiri VP285 |
| 45 | `SE_SIDM_KRISHNAMURTI_VP291` | Krishnamurti VP291 |
| 46 | `SE_SIDM_LAHIRI_ICRC` | Lahiri ICRC |
| 255 | `SE_SIDM_USER` | User-defined |

#### Sidereal Projection Bits

| Bit | Constant | Effect |
|-----|----------|--------|
| 256 | `SE_SIDBIT_ECL_T0` | Project onto ecliptic of reference epoch t0 |
| 512 | `SE_SIDBIT_SSY_PLANE` | Project onto solar system invariable plane |
| 1024 | `SE_SIDBIT_USER_UT` | Interpret user-defined t0 as UT (not TT) |
| 2048 | `SE_SIDBIT_ECL_DATE` | Measure ayanamsa on ecliptic of date |
| 4096 | `SE_SIDBIT_NO_PREC_OFFSET` | Don't apply constant precession offset |
| 8192 | `SE_SIDBIT_PREC_ORIG` | Use original precession model of the ayanamsa |

#### Algorithm

Most predefined ayanamsas work by computing a reference fixed star's longitude (e.g., Spica at 180° for True Citra) at J2000, then measuring how far the vernal equinox has precessed from that anchor. "True" ayanamsas use the actual current star position; others use a fixed offset plus precession. User-defined ayanamsas store a reference epoch and value, then apply precession from that epoch to the requested date.

<!-- INTERNALS TODO (§3.5): Document swi_get_ayanamsa_ex() internals; trace how "true star" ayanamsas compute the reference star position and derive precession offset; explain the galactic coordinate transformations; document get_aya_correction() and its precession model handling. -->

---

### 3.6 House Systems

**What it computes:** Astrological house cusps (the division of the local sky into 12 sectors) and key angles (Ascendant, MC, Vertex, etc.) for 24 different house systems. Optionally computes daily speeds of all cusps and angles.

#### Public API Functions

**`int swe_houses(double tjd_ut, double geolat, double geolon, int hsys, double *cusps, double *ascmc)`**
Basic house calculation from date/time and location.
- `cusps[13]`: cusps[1]–cusps[12] are the 12 house cusps (cusps[0] unused). For Gauquelin: cusps[0]–cusps[36] (37 entries).
- `ascmc[10]`: [0]=Ascendant, [1]=MC, [2]=ARMC, [3]=Vertex, [4]=Equatorial Ascendant, [5]=Co-Ascendant (Koch), [6]=Co-Ascendant (Munkasey), [7]=Polar Ascendant

**`int swe_houses_ex(double tjd_ut, int32 iflag, double geolat, double geolon, int hsys, double *cusps, double *ascmc)`**
Extended version with flags for sidereal mode (`SEFLG_SIDEREAL`), nutation suppression (`SEFLG_NONUT`), and ephemeris selection.

**`int swe_houses_ex2(double tjd_ut, int32 iflag, double geolat, double geolon, int hsys, double *cusps, double *ascmc, double *cusp_speed, double *ascmc_speed, char *serr)`**
Full version with optional speed arrays (daily motion of cusps and angles in degrees/day).

**`int swe_houses_armc(double armc, double geolat, double eps, int hsys, double *cusps, double *ascmc)`**
Calculates houses from pre-computed ARMC (Right Ascension of MC in degrees), geographic latitude, and obliquity. Used for composite charts or when ARMC is known.

**`int swe_houses_armc_ex2(double armc, double geolat, double eps, int hsys, double *cusps, double *ascmc, double *cusp_speed, double *ascmc_speed, char *serr)`**
Full ARMC-based calculation with speeds.

**`double swe_house_pos(double armc, double geolat, double eps, int hsys, double *xpin, char *serr)`**
Returns the house position of a planet as a decimal value. E.g., 1.5 = middle of 1st house, 7.0 = cusp of 7th house.
- `xpin[2]`: Planet's ecliptic [longitude, latitude]

**`const char *swe_house_name(int hsys)`**
Returns the name of a house system given its character code.

#### House Systems (24)

| Code | Name | Method |
|------|------|--------|
| P | Placidus **(default)** | Time-based: each house cusp is where a degree on the ecliptic has spent 1/3 of its semi-diurnal or semi-nocturnal arc. Most popular system in Western astrology. |
| K | Koch | Time-based variant: like Placidus but measuring the time for the Ascendant degree's diurnal arc rather than each individual degree. |
| O | Porphyry | Quadrant trisection: divides each quadrant (Asc–MC, MC–Dsc, etc.) into three equal ecliptic arcs. Ancient system. |
| R | Regiomontanus | Space-based: great circles through the north/south celestial pole divide the celestial equator into 30° segments, then project onto the ecliptic. |
| C | Campanus | Space-based: great circles through the north/south points of the horizon divide the prime vertical into 30° segments. |
| A, E | Equal | Each house is exactly 30° of ecliptic longitude from the Ascendant. |
| V | Vehlow Equal | Equal houses but with Ascendant at the middle (not start) of house 1. |
| W | Whole Sign | Each house is an entire zodiac sign; house 1 = the sign containing the Ascendant. |
| D | Equal (MC) | Equal houses measured from MC instead of Ascendant. |
| N | Equal (Aries) | Equal houses starting from 0° Aries. |
| B | Alcabitius | Semi-arc: divides the diurnal and nocturnal semi-arcs of the Ascendant into thirds on the equator, then projects onto ecliptic. |
| M | Morinus | Projects 30° equatorial segments through the ecliptic pole. No latitude dependence. |
| X | Meridian / Axial Rotation | Equatorial house system; projects equator segments onto ecliptic via the ecliptic pole. |
| H | Horizon / Azimuthal | Based on the horizon plane; houses measured from the East point. |
| T | Polich/Page (Topocentric) | Similar to Placidus but uses a topocentric (observer-centered) reference. Proposed by Polich and Page. |
| U | Krusinski-Pisa-Goelzer | Modern system using vertical circles. |
| G | Gauquelin | 36 sectors (not 12 houses); divides diurnal/nocturnal arcs into 18 equal time portions each. Used in statistical astrology research. |
| S | Sripati | Indian Vedic system: midpoints between Porphyry cusps define house boundaries. |
| F | Carter (Poli-Equatorial) | Carter's polar equatorial system. |
| I | Sunshine | Solar-declination-based system; requires Sun's declination. |
| i | Sunshine (alt.) | Alternative Sunshine implementation. |
| J | Savard-A | Savard's variant. |
| L | Pullen SD | Pullen Sinusoidal Delta. |
| Q | Pullen SR | Pullen Sinusoidal Ratio. |
| Y | APC | Astrological Program Center system. |

<!-- INTERNALS TODO (§3.6): Document CalcH() and how each house system algorithm works geometrically (the trigonometric formulas for Placidus semi-arc division, Regiomontanus great circles, Campanus vertical circles, etc.); explain speed computation by finite differences. -->

---

### 3.7 Eclipses & Occultations

**What it computes:** Times, locations, and circumstances of solar eclipses, lunar eclipses, and lunar occultations of planets/stars. Supports global searches (next eclipse anywhere on Earth) and local searches (next eclipse visible from a specific location). Returns contact times, magnitudes, obscuration, and geometric details.

#### Public API Functions

**Solar Eclipses:**

**`int32 swe_sol_eclipse_when_glob(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr)`**
Finds next (or previous) solar eclipse globally.
- `ifltype`: Filter — `SE_ECL_TOTAL`, `SE_ECL_ANNULAR`, `SE_ECL_PARTIAL`, `SE_ECL_CENTRAL`, `SE_ECL_NONCENTRAL`, or combinations
- `tret[10]`: Contact times — [0]=maximum, [2]=first contact, [3]=second, [4]=third, [5]=fourth
- `backward`: 1 = search backward in time
- Algorithm: Steps through lunations using Meeus' eclipse prediction formulas (argument F of lunar latitude), then refines iteratively.

**`int32 swe_sol_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr)`**
Finds next solar eclipse visible from a specific location.
- `geopos[3]`: Observer [longitude, latitude, altitude]
- `attr[]`: Eclipse attributes at that location (magnitude, ratio, obscuration, etc.)

**`int32 swe_sol_eclipse_how(double tjd, int32 ifl, double *geopos, double *attr, char *serr)`**
Computes eclipse circumstances at a specific time and location.
- `attr[]`: [magnitude, diameter_ratio, fraction_obscured, core_shadow_diameter_km, azimuth_sun, true_altitude, apparent_altitude, angular_distance_moon_sun_limbs]

**`int32 swe_sol_eclipse_where(double tjd, int32 ifl, double *geopos, double *attr, char *serr)`**
Returns the geographic location where a solar eclipse is at maximum at a given time.

**Lunar Eclipses:**

**`int32 swe_lun_eclipse_when(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr)`**
Finds next lunar eclipse globally.
- `ifltype`: `SE_ECL_TOTAL`, `SE_ECL_PARTIAL`, `SE_ECL_PENUMBRAL`

**`int32 swe_lun_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr)`**
Finds next lunar eclipse visible from a specific location. Checks Moon rise/set times.
- Returns visibility flags: `SE_ECL_VISIBLE`, `SE_ECL_MAX_VISIBLE`, `SE_ECL_PARTBEG_VISIBLE`, etc.

**`int32 swe_lun_eclipse_how(double tjd_ut, int32 ifl, double *geopos, double *attr, char *serr)`**
Computes lunar eclipse circumstances at a given time.

**Occultations (Moon covering a planet or star):**

**`int32 swe_lun_occult_when_glob(double tjd_start, int32 ipl, char *starname, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr)`**
Finds next occultation of a planet or star by the Moon, globally.

**`int32 swe_lun_occult_when_loc(double tjd_start, int32 ipl, char *starname, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr)`**
Finds next occultation visible from a specific location.

**`int32 swe_lun_occult_where(double tjd, int32 ipl, char *starname, int32 ifl, double *geopos, double *attr, char *serr)`**
Returns geographic location where an occultation is at maximum.

#### Eclipse Type Flags

| Flag | Value | Meaning |
|------|-------|---------|
| `SE_ECL_CENTRAL` | 1 | Central eclipse |
| `SE_ECL_NONCENTRAL` | 2 | Non-central eclipse |
| `SE_ECL_TOTAL` | 4 | Total eclipse |
| `SE_ECL_ANNULAR` | 8 | Annular eclipse |
| `SE_ECL_PARTIAL` | 16 | Partial eclipse |
| `SE_ECL_ANNULAR_TOTAL` / `SE_ECL_HYBRID` | 32 | Hybrid (annular-total) |
| `SE_ECL_PENUMBRAL` | 64 | Penumbral lunar eclipse |
| `SE_ECL_ONE_TRY` | 32768 | Check only the next conjunction (don't search further) |

<!-- INTERNALS TODO (§3.7): Document eclipse_where(), eclipse_how(), lun_eclipse_how() internals; trace the Meeus eclipse prediction algorithm (lunation stepping, F-argument filtering, Besselian elements); explain shadow cone geometry and magnitude/obscuration formulas; document iterative contact-time refinement. -->

---

### 3.8 Rise, Set & Transit

**What it computes:** Times of rising, setting, upper meridian transit, and lower meridian transit for any planet or fixed star, at a given geographic location. Includes atmospheric refraction, disc-size options, and twilight calculations.

#### Public API Functions

**`int32 swe_rise_trans(double tjd_ut, int32 ipl, char *starname, int32 epheflag, int32 rsmi, double *geopos, double atpress, double attemp, double *tret, char *serr)`**
Computes rise/set/transit time.
- `rsmi`: Event type flags (can be OR'd):
  - `SE_CALC_RISE` (1): Rising
  - `SE_CALC_SET` (2): Setting
  - `SE_CALC_MTRANSIT` (4): Upper meridian transit (culmination)
  - `SE_CALC_ITRANSIT` (8): Lower meridian transit (anti-culmination)
- Returns: 0 = found, 1 = circumpolar (never rises), −1 = circumpolar (never sets), ERR = error
- Algorithm: For latitudes <60° (Moon/planets) or <65° (Sun), uses a fast algorithm based on estimated semi-diurnal arc. For higher latitudes, delegates to the slow method.

**`int32 swe_rise_trans_true_hor(double tjd_ut, int32 ipl, char *starname, int32 epheflag, int32 rsmi, double *geopos, double atpress, double attemp, double horhgt, double *tret, char *serr)`**
Same as `swe_rise_trans()` but with a custom horizon height.
- `horhgt`: Horizon height in degrees (negative = below mathematical horizon); −100 = auto-calculate ocean horizon dip from altitude
- Algorithm: Computes positions at 2-hour intervals over the day, finds altitude crossings and culminations via parabolic interpolation, then refines iteratively.

#### Rise/Set Modifier Flags

| Flag | Value | Effect |
|------|-------|--------|
| `SE_BIT_DISC_CENTER` | 256 | Rise/set of disc center (not upper limb) |
| `SE_BIT_DISC_BOTTOM` | 8192 | Rise/set of lower limb |
| `SE_BIT_NO_REFRACTION` | 512 | Ignore atmospheric refraction |
| `SE_BIT_GEOCTR_NO_ECL_LAT` | 128 | Use geocentric position, ignore ecliptic latitude |
| `SE_BIT_CIVIL_TWILIGHT` | 1024 | Sun center at −6° (civil twilight) |
| `SE_BIT_NAUTIC_TWILIGHT` | 2048 | Sun center at −12° (nautical twilight) |
| `SE_BIT_ASTRO_TWILIGHT` | 4096 | Sun center at −18° (astronomical twilight) |
| `SE_BIT_FIXED_DISC_SIZE` | 16384 | Use average disc size, ignore distance variation |
| `SE_BIT_HINDU_RISING` | combined | Hindu rising: disc center, no refraction, geocentric |

<!-- INTERNALS TODO (§3.8): Document the fast vs. slow rise/set algorithms; trace the 2-hour sampling and parabolic interpolation method; explain how circumpolar detection works; document horizon dip calculation for elevated observers. -->

---

### 3.9 Planetary Phenomena

**What it computes:** Phase angle, phase (illuminated fraction), elongation from the Sun, apparent diameter, and visual magnitude for any planet.

#### Public API Functions

**`int32 swe_pheno(double tjd, int32 ipl, int32 iflag, double *attr, char *serr)`**
Computes phenomena at Ephemeris Time.
- `attr[20+]`:
  - [0] = Phase angle (Sun–planet–Earth angle, degrees)
  - [1] = Phase (0 = new/dark, 1 = full/bright): (1 + cos(phase_angle)) / 2
  - [2] = Elongation from Sun (degrees)
  - [3] = Apparent diameter (degrees)
  - [4] = Apparent magnitude
- Algorithm: Computes geocentric and heliocentric positions, derives phase angle from dot product of vectors, magnitude from standard planetary photometric formulas.

**`int32 swe_pheno_ut(double tjd_ut, int32 ipl, int32 iflag, double *attr, char *serr)`**
UT wrapper.

<!-- INTERNALS TODO (§3.9): Document the planetary magnitude formulas (which photometric model is used for each planet); trace the phase angle computation from dot products of position vectors. -->

---

### 3.10 Nodes & Apsides

**What it computes:** Ascending and descending nodes (where a body crosses the ecliptic plane) and apsides (perihelion/aphelion or perigee/apogee) for any planet or asteroid. Can compute mean or osculating values.

#### Public API Functions

**`int32 swe_nod_aps(double tjd_et, int32 ipl, int32 iflag, int32 method, double *xnasc, double *xndsc, double *xperi, double *xaphe, char *serr)`**
Computes nodes and apsides.
- `method`: `SE_NODBIT_MEAN` (mean elements), `SE_NODBIT_OSCU` (osculating), `SE_NODBIT_OSCU_BAR` (osculating about SSB), `SE_NODBIT_FOPOINT` (focal point instead of aphelion)
- Each output array is [6]: [longitude, latitude, distance, speed_lon, speed_lat, speed_dist]
- Algorithm: Mean nodes from analytical orbital element theories (Laskar et al.). Osculating nodes from cross product of position and velocity vectors (the angular momentum vector defines the orbital plane, whose intersection with the ecliptic gives the nodes).

**`int32 swe_nod_aps_ut(double tjd_ut, int32 ipl, int32 iflag, int32 method, double *xnasc, double *xndsc, double *xperi, double *xaphe, char *serr)`**
UT wrapper.

#### Node/Apse Method Flags

| Flag | Value | Meaning |
|------|-------|---------|
| `SE_NODBIT_MEAN` | 1 | Mean nodes/apsides from secular theory |
| `SE_NODBIT_OSCU` | 2 | Osculating (instantaneous) nodes/apsides |
| `SE_NODBIT_OSCU_BAR` | 4 | Osculating about solar system barycenter |
| `SE_NODBIT_FOPOINT` | 256 | Return focal point of orbit instead of aphelion |

<!-- INTERNALS TODO (§3.10): Document how osculating nodes are computed from the angular momentum vector (cross product of r and v); trace mean node computation from Laskar's secular theory; explain the SE_NODBIT_FOPOINT focal point alternative. -->

---

### 3.11 Orbital Elements

**What it computes:** Keplerian osculating orbital elements (semi-major axis, eccentricity, inclination, etc.) and extreme distances for any planet or asteroid.

#### Public API Functions

**`int32 swe_get_orbital_elements(double tjd_et, int32 ipl, int32 iflag, double *dret, char *serr)`**
Returns osculating elements relative to the Sun (or SSB for outer planets with `SEFLG_BARYCTR`).
- `dret[50+]`:
  - [0] = Semi-major axis (AU)
  - [1] = Eccentricity
  - [2] = Inclination (degrees)
  - [3] = Longitude of ascending node
  - [4] = Argument of periapsis
  - [5] = Longitude of periapsis
  - [6] = Mean anomaly
  - [7] = True anomaly
  - [8] = Eccentric anomaly
  - [9] = Mean longitude
  - [10] = Sidereal orbital period (tropical years)
  - [11] = Mean daily motion (degrees/day)
- Not valid for: Sun, Moon, nodes, apsides

**`int32 swe_orbit_max_min_true_distance(double tjd_et, int32 ipl, int32 iflag, double *dmax, double *dmin, double *dtrue, char *serr)`**
Returns maximum, minimum, and current geocentric distances by sampling both orbits in 2° steps and iteratively refining extrema.

<!-- INTERNALS TODO (§3.11): Document the vector-to-elements conversion (how a, e, i, Ω, ω, M are derived from position/velocity vectors); explain Kepler equation solving; trace the orbit-sampling algorithm in swe_orbit_max_min_true_distance(). -->

---

### 3.12 Zodiacal Crossings

**What it computes:** The exact moment when the Sun, Moon, or a planet crosses a specified ecliptic longitude (e.g., 0° Aries for equinoxes). Also: Moon's node crossings (zero-latitude moments).

#### Public API Functions

**`double swe_solcross(double x2cross, double jd_et, int32 flag, char *serr)`**
Next Sun crossing of longitude `x2cross` after `jd_et` (ET). Returns JD of crossing.
- Algorithm: Newton's method — estimates initial time from mean solar speed (360°/365.24 days), then iterates: `t_next = t + (target − current_lon) / speed`, converging to <1 milliarcsecond.

**`double swe_solcross_ut(double x2cross, double jd_ut, int32 flag, char *serr)`**
UT version.

**`double swe_mooncross(double x2cross, double jd_et, int32 flag, char *serr)`**
Next Moon crossing of longitude `x2cross`. Uses mean lunar speed (360°/27.32 days) for initial estimate.

**`double swe_mooncross_ut(double x2cross, double jd_ut, int32 flag, char *serr)`**
UT version.

**`double swe_mooncross_node(double jd_et, int32 flag, double *xlon, double *xlat, char *serr)`**
Next Moon node crossing (latitude = 0). Returns JD and outputs longitude/latitude at crossing.
- Algorithm: Steps forward detecting latitude sign changes, then Newton-refines using latitude and latitude speed.

**`double swe_mooncross_node_ut(double jd_ut, int32 flag, double *xlon, double *xlat, char *serr)`**
UT version.

**`int32 swe_helio_cross(int32 ipl, double x2cross, double jd_et, int32 iflag, int32 dir, double *jd_cross, char *serr)`**
Heliocentric longitude crossing for any planet or asteroid.
- `dir`: ≥0 for forward search, <0 for backward
- Not valid for: Sun, Moon, nodes, apogees

**`int32 swe_helio_cross_ut(int32 ipl, double x2cross, double jd_ut, int32 iflag, int32 dir, double *jd_cross, char *serr)`**
UT version.

<!-- INTERNALS TODO (§3.12): Document the Newton iteration convergence criteria and edge cases (retrograde planets, fast-moving Moon); explain how the initial time estimate is refined. -->

---

### 3.13 Heliacal Events & Visibility

**What it computes:** Heliacal rising and setting dates (when a planet or star first/last becomes visible near the Sun), arcus visionis (minimum Sun depression for visibility), and limiting visual magnitude under given atmospheric and observer conditions. Uses Schaefer's comprehensive atmospheric extinction and human vision model.

#### Public API Functions

**`int32 swe_heliacal_ut(double tjdstart_ut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 iflag, double *dret, char *serr)`**
Searches for the next heliacal event.
- `geopos[3]`: [longitude, latitude, altitude_meters]
- `datm[4]`: [pressure_mbar, temperature_°C, relative_humidity_%, visibility_km]
- `dobs[6]`: [age_years, Snellen_ratio, binocular_flag, magnification, aperture_mm, transmission]
- `TypeEvent`: `SE_HELIACAL_RISING` (1), `SE_HELIACAL_SETTING` (2), `SE_EVENING_FIRST` (3), `SE_MORNING_LAST` (4)
- `dret[0]`: JD of event
- Algorithm: Steps through synodic periods, at each step computing when the Sun is low enough for the object to become visible against twilight (using Schaefer's model for sky brightness, atmospheric extinction, and contrast threshold).

**`int32 swe_heliacal_pheno_ut(double tjd_ut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 helflag, double *darr, char *serr)`**
Detailed visibility phenomena at a specific moment (not a search).
- `darr[28+]`: Comprehensive output including object/Sun altitudes and azimuths, TAV, arcus visionis, DAZ, ARCL, extinction coefficient, first/best/last visibility times, and for the Moon: Yallop crescent parameters (width, q-value, visibility criteria A–F).

**`int32 swe_vis_limit_mag(double tjdut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 helflag, double *dret, char *serr)`**
Calculates the faintest star visible at a given moment (limiting magnitude).
- `dret[8]`: [limiting_mag, obj_alt, obj_az, sun_alt, sun_az, moon_alt, moon_az, obj_mag]
- Algorithm (Schaefer model): Computes background sky brightness from five sources (twilight, airglow, zodiacal light, starlight, moonlight), applies atmospheric extinction, determines contrast threshold based on photopic/scotopic/mesopic vision and observer acuity.

**`int32 swe_heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr)`**
Returns optimal object altitude, arcus visionis, and Sun altitude difference for heliacal visibility.

**`int32 swe_topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr)`**
Calculates the minimum Sun depression (arcus visionis) needed for an object of given magnitude to be visible. Uses bisection search over Sun altitudes from −45° to 0°.

#### Heliacal Event Type Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `SE_HELIACAL_RISING` / `SE_MORNING_FIRST` | 1 | First morning visibility |
| `SE_HELIACAL_SETTING` / `SE_EVENING_LAST` | 2 | Last evening visibility |
| `SE_EVENING_FIRST` | 3 | First evening visibility |
| `SE_MORNING_LAST` | 4 | Last morning visibility |

#### Heliacal Flags

| Flag | Value | Effect |
|------|-------|--------|
| `SE_HELFLAG_LONG_SEARCH` | 128 | Search up to 10 synodic periods |
| `SE_HELFLAG_HIGH_PRECISION` | 256 | Include nutation corrections |
| `SE_HELFLAG_OPTICAL_PARAMS` | 512 | Use optical instrument parameters from dobs |
| `SE_HELFLAG_NO_DETAILS` | 1024 | Skip detailed output |
| `SE_HELFLAG_VISLIM_DARK` | 4096 | Only check astronomical twilight |
| `SE_HELFLAG_VISLIM_NOMOON` | 8192 | Ignore Moon's contribution to sky brightness |

<!-- INTERNALS TODO (§3.13): Document Schaefer's sky brightness model in detail (the five light sources, their altitude dependencies, the extinction integrals); trace VisLimMagn() and TopoArcVisionis() internals; explain the Yallop crescent visibility criteria (A–F grades); document MoonEventJDut() for lunar crescent search. -->

---

### 3.14 Atmospheric & Optical

**What it computes:** Atmospheric refraction corrections (true altitude ↔ apparent altitude), including models for elevated observers, variable lapse rates, and ocean horizon dip.

#### Public API Functions

**`double swe_refrac(double inalt, double atpress, double attemp, int32 calc_flag)`**
Simple refraction correction.
- `calc_flag`: `SE_TRUE_TO_APP` or `SE_APP_TO_TRUE`
- `atpress`: Atmospheric pressure in mbar
- `attemp`: Temperature in °C
- Returns: Corrected altitude
- Algorithm: Saemundsson formula (Meeus) for altitudes >15°; polynomial for −5° to 15°. Pressure/temperature compensation factor applied.

**`double swe_refrac_extended(double inalt, double geoalt, double atpress, double attemp, double lapse_rate, int32 calc_flag, double *dret)`**
Advanced refraction for elevated observers.
- `geoalt`: Observer altitude above sea level (meters)
- `lapse_rate`: Temperature lapse rate (K/m), typically 0.0065
- `dret[4]`: [true_alt, apparent_alt, refraction, dip_of_horizon]
- Handles objects below the geometric horizon as seen from elevated locations; computes ocean horizon dip.

**`void swe_set_lapse_rate(double lapse_rate)`**
Sets the default atmospheric lapse rate for refraction calculations.

#### Refraction Constants

| Constant | Value | Direction |
|----------|-------|-----------|
| `SE_TRUE_TO_APP` | 0 | True altitude → apparent (add refraction) |
| `SE_APP_TO_TRUE` | 1 | Apparent altitude → true (subtract refraction) |

<!-- INTERNALS TODO (§3.14): Document the Saemundsson refraction formula and its polynomial variant for low altitudes; trace the extended refraction model's handling of elevated observers and geometric horizon dip; explain lapse rate effects on refraction. -->

---

### 3.15 Gauquelin Sectors

**What it computes:** The Gauquelin sector position (1–36) of a planet, used in Michel Gauquelin's statistical astrology research. The diurnal and nocturnal arcs are each divided into 18 equal time segments.

#### Public API Functions

**`int32 swe_gauquelin_sector(double t_ut, int32 ipl, char *starname, int32 iflag, int32 imeth, double *geopos, double atpress, double attemp, double *dgsect, char *serr)`**
- `imeth`: Method selector:
  - 0 = Geometric from ecliptic position (uses `swe_house_pos()`)
  - 1 = Geometric without ecliptic latitude
  - 2 = Rise/set-based without refraction
  - 3 = Rise/set-based with refraction
- `dgsect`: Output sector number (1.0–36.0)
- Returns ERR for circumpolar bodies (that never rise or set)

<!-- INTERNALS TODO (§3.15): Document how Gauquelin sectors are computed from rise/set arcs; trace the geometric vs. rise/set methods and their differences. -->

---

### 3.16 Utility & Formatting

**What it computes:** Angle normalization, degree/minute/second splitting, midpoint calculation, centisecond arithmetic (legacy Placalc compatibility), and string formatting for time and coordinates.

#### Public API Functions

**Angle Normalization:**

| Function | Input | Output Range |
|----------|-------|-------------|
| `double swe_degnorm(double x)` | Degrees | [0, 360) |
| `double swe_radnorm(double x)` | Radians | [0, 2π) |
| `double swe_deg_midp(double x1, double x0)` | Two angles (°) | Midpoint on shortest arc |
| `double swe_rad_midp(double x1, double x0)` | Two angles (rad) | Midpoint on shortest arc |

**Degree Splitting:**

**`void swe_split_deg(double ddeg, int32 roundflag, int32 *ideg, int32 *imin, int32 *isec, double *dsecfr, int32 *isgn)`**
Splits decimal degrees into components.
- `roundflag` bits: `SE_SPLIT_DEG_ROUND_SEC` (1), `SE_SPLIT_DEG_ROUND_MIN` (2), `SE_SPLIT_DEG_ROUND_DEG` (4), `SE_SPLIT_DEG_ZODIACAL` (8, output sign 0–11), `SE_SPLIT_DEG_NAKSHATRA` (1024, 27 nakshatras of 13°20′), `SE_SPLIT_DEG_KEEP_SIGN` (16, don't round to next sign), `SE_SPLIT_DEG_KEEP_DEG` (32, don't round to next degree)
- In zodiacal mode: `isgn` = sign number (0–11), `ideg` = degrees within sign (0–29)

**Angular Differences (centisecond legacy):**

| Function | Input | Output Range |
|----------|-------|-------------|
| `centisec swe_csnorm(centisec p)` | Centiseconds | [0, 360°) |
| `centisec swe_difcsn(centisec p1, centisec p2)` | Two angles (cs) | p1−p2 in [0, 360°) |
| `centisec swe_difcs2n(centisec p1, centisec p2)` | Two angles (cs) | p1−p2 in [−180°, 180°) |
| `double swe_difdegn(double p1, double p2)` | Degrees | p1−p2 in [0, 360°) |
| `double swe_difdeg2n(double p1, double p2)` | Degrees | p1−p2 in [−180°, 180°) |
| `double swe_difrad2n(double p1, double p2)` | Radians | p1−p2 in [−π, π) |

**Other Utilities:**

| Function | Description |
|----------|-------------|
| `centisec swe_csroundsec(centisec x)` | Round to nearest second; rounds DOWN at 29°59′30″ to stay within sign |
| `int32 swe_d2l(double x)` | Double to int32 with rounding |

**String Formatting:**

| Function | Description | Example Output |
|----------|-------------|---------------|
| `char *swe_cs2timestr(CSEC t, int sep, AS_BOOL suppressZero, char *a)` | Centiseconds → "HH:MM:SS" | "14:30:00" |
| `char *swe_cs2lonlatstr(CSEC t, char pchar, char mchar, char *s)` | Centiseconds → "DDD°MM'SS" + direction | "122°30'45"E" |
| `char *swe_cs2degstr(CSEC t, char *a)` | Centiseconds → "DD°MM'SS" within sign | "15°30'00" |

<!-- INTERNALS TODO (§3.16): Document centisecond arithmetic internals; trace split_deg_nakshatra(). -->

---

## 4. Constants & Enumerations Reference

### Unit Conversion Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `SE_AUNIT_TO_KM` | 149,597,870.700 | AU to kilometers |
| `SE_AUNIT_TO_LIGHTYEAR` | 1/63,241.077 | AU to light-years |
| `SE_AUNIT_TO_PARSEC` | 1/206,264.806 | AU to parsecs |

### Calendar Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `SE_JUL_CAL` | 0 | Julian calendar |
| `SE_GREG_CAL` | 1 | Gregorian calendar |

### Tidal Acceleration Constants

| Constant | Value (″/cy²) | Ephemeris |
|----------|--------------|-----------|
| `SE_TIDAL_DE200` | −23.8946 | DE200 |
| `SE_TIDAL_DE403` | −25.580 | DE403 |
| `SE_TIDAL_DE404` | −25.580 | DE404 |
| `SE_TIDAL_DE405` | −25.826 | DE405 |
| `SE_TIDAL_DE406` | −25.826 | DE406 |
| `SE_TIDAL_DE421` | −25.85 | DE421 |
| `SE_TIDAL_DE430` | −25.82 | DE430 |
| `SE_TIDAL_DE431` | −25.80 | DE431 (default) |
| `SE_TIDAL_DE441` | −25.936 | DE441 |
| `SE_TIDAL_AUTOMATIC` | 999999 | Auto-select based on ephemeris |

### Ephemeris File Names

| Constant | Value |
|----------|-------|
| `SE_STARFILE` | `"sefstars.txt"` |
| `SE_ASTNAMFILE` | `"seasnam.txt"` |
| `SE_FICTFILE` | `"seorbel.txt"` |
| `SE_FNAME_DFT` | `"de431.eph"` |

### Model Selection Constants

| Constant | Value | Selects |
|----------|-------|---------|
| `SE_MODEL_DELTAT` | 0 | Delta-T model |
| `SE_MODEL_PREC_LONGTERM` | 1 | Long-term precession |
| `SE_MODEL_PREC_SHORTTERM` | 2 | Short-term precession |
| `SE_MODEL_NUT` | 3 | Nutation model |
| `SE_MODEL_BIAS` | 4 | Frame bias model |
| `SE_MODEL_JPLHOR_MODE` | 5 | JPL Horizons method |
| `SE_MODEL_JPLHORA_MODE` | 6 | JPL Horizons approximation |
| `SE_MODEL_SIDT` | 7 | Sidereal time model |

---

## 5. Configuration Functions

### Path & File Configuration

| Function | Effect |
|----------|--------|
| `void swe_set_ephe_path(const char *path)` | Sets directory for ephemeris data files. NULL = use `SE_DATA` env var or default. |
| `void swe_set_jpl_file(const char *fname)` | Sets JPL ephemeris filename (e.g., "de431.eph") |
| `void swe_close(void)` | Closes all ephemeris files and frees memory. Call before program exit or to reset state. |

### Observer Position

| Function | Effect |
|----------|--------|
| `void swe_set_topo(double geolon, double geolat, double geoalt)` | Sets observer position for topocentric calculations. East positive, North positive, meters above sea level. |

### Sidereal Mode

| Function | Effect |
|----------|--------|
| `void swe_set_sid_mode(int32 sid_mode, double t0, double ayan_t0)` | Activates sidereal zodiac with specified ayanamsa |

### Time & Tidal

| Function | Effect |
|----------|--------|
| `void swe_set_tid_acc(double t_acc)` | Sets lunar tidal acceleration for delta-T (″/cy²) |
| `double swe_get_tid_acc(void)` | Gets current tidal acceleration value |
| `void swe_set_delta_t_userdef(double dt)` | Override delta-T with a fixed value (days). `SE_DELTAT_AUTOMATIC` to restore. |

### Atmospheric

| Function | Effect |
|----------|--------|
| `void swe_set_lapse_rate(double lapse_rate)` | Sets default atmospheric lapse rate (K/m) |

### Calculation Models

| Function | Effect |
|----------|--------|
| `void swe_set_astro_models(char *samod, int32 iflag)` | Selects precession, nutation, delta-T, sidereal time, and frame bias models |
| `void swe_get_astro_models(char *samod, char *sdet, int32 iflag)` | Retrieves current model settings and detailed descriptions |

### Nutation

| Function | Effect |
|----------|--------|
| `void swe_set_interpolate_nut(AS_BOOL do_interpolate)` | Enables/disables nutation interpolation (performance vs. accuracy) |
