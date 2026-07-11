# Swiss Ephemeris — Variable Name Glossary

> Quick-reference mapping of recurring short variable names to their physical meanings.
> Organized by domain. Variables appear across multiple source files unless noted.

---

## Position & Velocity Arrays

All are `double[6]`: [x, y, z, vx, vy, vz] in AU and AU/day (equatorial J2000 Cartesian unless noted).

| Variable | Meaning | Context |
|----------|---------|---------|
| `xx` | Working position/velocity vector | General purpose, most functions |
| `xx0` | Original position before light-time correction | `app_pos_etc_plan()` |
| `xxsv` | Saved position for sidereal transforms or iteration | `app_pos_rest()`, `lunar_osc_elem()` |
| `xxsp` | Speed correction from changing light-time delay | `app_pos_etc_plan()` |
| `xp` | Planet position (often a pointer into `pdp->x`) | `sweplan()`, `jplplan()` |
| `xe`, `xxe` | Earth barycentric position | `sweplan()`, `app_pos_etc_moon()` |
| `xs`, `xxs` | Sun (barycentric or heliocentric) position | `sweplan()`, `app_pos_etc_sun()` |
| `xm`, `xxm` | Moon (geocentric or barycentric) position | `sweplan()`, `app_pos_etc_moon()` |
| `xobs` | Observer position (geocentric or topocentric) | `app_pos_etc_plan()` |
| `xobs2` | Observer position at retarded time (t − light_time) | `app_pos_etc_plan()` |
| `xemb` | Earth-Moon barycenter position | `sweplan()` |
| `xxp` | Temporary planet position | `sweplan()`, `jplplan()` |
| `xxctr` | Center body position (planetocentric mode) | `swe_calc_pctr()` |
| `dx` | Displacement vector: planet minus observer | Light-time iteration |
| `xnorm` | Orbital plane normal (angular momentum vector r × v) | `lunar_osc_elem()` |
| `xpos[3][6]` | Three positions at t−dt, t, t+dt | Osculating element speed calc |
| `xreturn[24]` | Final output in 4 coordinate systems (see below) | `pdp->xreturn` |

**`xreturn[24]` layout:**

| Indices | Coordinate System |
|---------|-------------------|
| 0–5 | Ecliptic polar (lon, lat, dist, dlon/dt, dlat/dt, ddist/dt) |
| 6–11 | Ecliptic Cartesian (x, y, z, vx, vy, vz) |
| 12–17 | Equatorial polar (RA, Dec, dist, dRA/dt, dDec/dt, ddist/dt) |
| 18–23 | Equatorial Cartesian |

---

## Eclipse Geometry (swecl.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `rm` | Moon position vector (Cartesian) | AU |
| `rs` | Sun/planet position vector (Cartesian) | AU |
| `rmt`, `rst` | Saved moon/sun position for iteration | AU |
| `e` | Sun-moon unit direction vector (shadow axis) | dimensionless |
| `et` | Sun-moon direction saved for iteration | dimensionless |
| `dm` | Distance of Moon from geocenter | AU |
| `ds` | Distance of Sun from geocenter | AU |
| `dsm` | Distance from Sun to Moon | AU |
| `dsmt` | Distance Sun-Moon (saved copy) | AU |
| `sinf1`, `cosf1` | Half-angle of umbral (core shadow) cone | dimensionless |
| `sinf2`, `cosf2` | Half-angle of penumbral (half-shadow) cone | dimensionless |
| `s0` | Moon's distance from fundamental plane: `−dot(rm, e)` | AU → km |
| `r0` | Shadow axis distance from geocenter: `√(dm² − s0²)` | AU → km |
| `d0` | Umbra diameter on fundamental plane | km |
| `D0` | Penumbra diameter on fundamental plane | km |
| `drad` | Radius of occulting body's disc | AU |
| `rmoon` | Angular or physical radius of Moon | AU or degrees |
| `dmoon` | Diameter of Moon (`2 × RMOON`) | AU |
| `de` | Equatorial radius of Earth | AU |
| `earthobl` | Earth oblateness factor: `1 − EARTH_OBLATENESS` | dimensionless |
| `dctr` | Center-to-center angular separation (Moon–Sun) | degrees |
| `rsun` | Angular radius of Sun | degrees |

**`dcore[]` array:**

| Index | Meaning |
|-------|---------|
| `dcore[0]` | Core shadow width at eclipse point (km) |
| `dcore[1]` | Penumbra width at eclipse point (km) |
| `dcore[2]` | `r0` — shadow axis distance from geocenter (km) |
| `dcore[3]` | `d0` — umbra diameter on fundamental plane (km) |
| `dcore[4]` | `D0` — penumbra diameter on fundamental plane (km) |
| `dcore[5]` | `cosf1` |
| `dcore[6]` | `cosf2` |

**`attr[]` — eclipse/phenomena attributes:**

| Index | Solar Eclipse | Lunar Eclipse | Phenomena |
|-------|--------------|---------------|-----------|
| 0 | Magnitude (fraction of solar diameter covered) | Umbral magnitude | Phase angle (degrees) |
| 1 | Lunar/solar diameter ratio | Penumbral magnitude | Phase (0–1) |
| 2 | Obscuration (fraction of area covered) | — | Elongation (degrees) |
| 3 | Core shadow diameter (km) | — | Apparent diameter (degrees) |
| 4 | Azimuth of Sun (degrees) | — | Apparent magnitude |
| 5 | True altitude of Sun (degrees) | — | — |
| 6 | Apparent altitude of Sun (degrees) | — | — |
| 7 | Angular distance Moon–Sun center (degrees) | — | — |
| 8 | NASA magnitude | — | — |
| 9 | Saros series number | Saros series number | — |
| 10 | Saros member number | Saros member number | — |

**`tret[]` — eclipse contact times (JD UT):**

| Index | Meaning |
|-------|---------|
| 0 | Maximum eclipse |
| 1 | Local apparent noon |
| 2 | Partial phase begin |
| 3 | Partial phase end |
| 4 | Totality/annularity begin |
| 5 | Totality/annularity end |
| 6 | Center line begin |
| 7 | Center line end |

---

## Rise/Set (swecl.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `h[]` | Altitude samples at 2-hour intervals | degrees |
| `xaz` | Azimuth/altitude output from `swe_azalt()` | degrees |
| `rdi` | Refraction + disc-radius correction | degrees |
| `dd` | Physical diameter of body disc | meters or km |
| `sda` | Semi-diurnal arc | degrees |
| `armc` | Right ascension of meridian crossing | degrees |
| `md` | Meridian distance of object | degrees |
| `horhgt` | Horizon height / dip angle | degrees |
| `twohrs` | Two-hour interval: `2.0/24.0` | days |
| `decl` | Declination of object | degrees |

---

## Orbital Elements (swecl.c, sweph.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `sema` | Semi-major axis: `a = 1/(2/r − v²/GM)` | AU |
| `ecce` | Eccentricity: `e = √(1 − p/a)` | dimensionless |
| `incl` | Orbital inclination | degrees |
| `parg` | Argument of perihelion (ω) | degrees |
| `ny` | True anomaly (ν) | radians |
| `uu` | Argument of latitude: `u = ω + ν` | radians |
| `pp` | Semi-latus rectum: `p = h²/GM` | AU |
| `cosE`, `sinE` | Eccentric anomaly trig functions | dimensionless |
| `sinnode`, `cosnode` | Ascending node direction | dimensionless |
| `sinincl`, `cosincl` | Inclination trig functions | dimensionless |
| `sinu`, `cosu` | Argument of latitude trig functions | dimensionless |
| `Gmsm` | GM × total system mass (in AU³/day² units) | AU³/day² |
| `v2` | Velocity squared | (AU/day)² |
| `rxy` | Distance projected onto XY (ecliptic) plane | AU |
| `rxyz` | Total 3D distance | AU |
| `xna` | Ascending node position vector | AU |
| `xnd` | Descending node position vector | AU |
| `xpe` | Perihelion position vector | AU |
| `xap` | Aphelion position vector | AU |

---

## Precession (swephlib.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `pA` | General precession in longitude | radians |
| `Z` | First Euler precession angle | radians |
| `z` | Third Euler precession angle | radians |
| `TH` | Second Euler angle (Θ, theta) | radians |
| `T` | Julian centuries from J2000: `(J − 2451545.0)/36525` | centuries |
| `T2` | T squared | centuries² |
| `R[3]` | Position vector being precessed | AU |
| `pAcof` | Precession polynomial coefficient array | — |
| `nodecof` | Node of moving ecliptic coefficients | — |
| `inclcof` | Inclination of moving ecliptic coefficients | — |

---

## Nutation (swephlib.c)

| Variable | Meaning | Standard Symbol |
|----------|---------|-----------------|
| `MM` | Mean anomaly of Moon | l |
| `MS` | Mean anomaly of Sun | l' |
| `FF` | Moon's argument of latitude | F |
| `DD` | Mean elongation of Moon from Sun | D |
| `OM` | Longitude of Moon's ascending node | Ω |
| `nutlo[0]` | Nutation in longitude | Δψ (radians) |
| `nutlo[1]` | Nutation in obliquity | Δε (radians) |
| `dpsi` | Nutation in longitude | radians |
| `deps` | Nutation in obliquity | radians |
| `ss[5][8]` | Precomputed sines of argument multiples | dimensionless |
| `cc[5][8]` | Precomputed cosines of argument multiples | dimensionless |

---

## Obliquity & Transforms (swephlib.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `eps` | Obliquity of ecliptic (ε) | radians |
| `seps`, `ceps` | sin(ε), cos(ε) — precomputed | dimensionless |
| `xpo` | Input position (polar or Cartesian) | varies |
| `xpn` | Output position (transformed) | varies |
| `l[]` | Polar coords: [lon, lat, radius] | radians, AU |
| `x[]` | Cartesian coords: [x, y, z] | AU |

---

## Aberration & Deflection (sweph.c, swephlib.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `v` | Earth velocity vector (fraction of c) | dimensionless |
| `v2` | v² (velocity squared) | dimensionless |
| `b_1` | Lorentz factor: `√(1 − v²)` | dimensionless |
| `f1` | `(u⃗·v⃗)/|u⃗|` — projection of velocity onto line of sight | dimensionless |
| `f2` | `1 + f1/(1 + b_1)` | dimensionless |
| `ru` | Magnitude of position vector | AU |
| `u` | Planet geocentric unit vector | dimensionless |
| `e` | Earth heliocentric unit vector (in deflection) | dimensionless |
| `q` | Planet heliocentric unit vector (in deflection) | dimensionless |
| `g1` | `2GM/(c²·r_E)` — deflection strength | radians |
| `g2` | `1 + q̂·ê` — geometric factor | dimensionless |
| `uq`, `ue`, `qe` | Dot products for deflection geometry | dimensionless |
| `dt` | Light-time delay | days |
| `dtsave_for_defl` | Saved light-time for deflection step | days |

---

## House Systems (swehouse.c)

| Variable | Meaning | Units |
|----------|---------|-------|
| `th` | Sidereal time (θ, ARMC) | degrees (0–360) |
| `fi` | Geographic latitude (φ) | degrees |
| `ekl` | Obliquity of ecliptic (ε) | degrees |
| `hsy` | House system code | char (A–Z) |
| `hsp` | House system parameters struct | `struct houses` |
| `sine`, `cose` | sin(ε), cos(ε) | dimensionless |
| `tane` | tan(ε) | dimensionless |
| `tanfi` | tan(φ) | dimensionless |
| `sinfi`, `cosfi` | sin(φ), cos(φ) | dimensionless |
| `fh1`, `fh2` | Pole heights for intermediate house cusps | degrees |
| `xh1`, `xh2` | RA offsets from horizon for Campanus | degrees |
| `rectasc` | Right ascension of house cusp being computed | degrees |
| `tant` | tan of intermediate declination angle | dimensionless |
| `f` | Computed pole height for iterative systems | degrees |
| `cuspsv` | Previous iteration's cusp value (convergence check) | degrees |
| `ih` | House index being computed (1–12 or 1–36) | integer |
| `a` | Auxiliary angle: `arcsin(tan(lat) × tan(ε))` | degrees |
| `ac` | Ascendant | degrees |
| `mc` | Midheaven (MC) | degrees |
| `armc` | Right ascension of MC | degrees |
| `acmc` | Angular distance Ascendant − MC | degrees |
| `ad3` | Ascensional difference divided by 3 (Koch system) | degrees |
| `sda` | Semi-diurnal arc | degrees |
| `sd3`, `sn3` | Semi-diurnal/semi-nocturnal arc ÷ 3 (Alcabitius) | degrees |

---

## Heliacal Visibility (swehel.c)

### Sky Brightness (nanoLamberts)

| Variable | Meaning |
|----------|---------|
| `Bday` | Daylight sky brightness |
| `Btwi` | Twilight sky brightness |
| `Bm` | Moonlight contribution |
| `Bn` | Night sky (zodiacal light + starlight) |
| `Bsk`, `Bsky` | Total sky brightness |

### Extinction Coefficients

| Variable | Meaning | Scale Height |
|----------|---------|-------------|
| `kR` | Rayleigh scattering | 8515 m |
| `kt` | Aerosol scattering | 3745 m |
| `kOZ` | Ozone absorption | 20000 m |
| `kW` | Water vapor absorption | 3000 m |
| `kX` | Total extinction at object's altitude | — |

### Geometry & Angles

| Variable | Meaning | Units |
|----------|---------|-------|
| `AltS` | Altitude of Sun | degrees |
| `AltO` | Altitude of object | degrees |
| `AltM` | Altitude of Moon | degrees |
| `AziS` | Azimuth of Sun | degrees |
| `AziO` | Azimuth of object | degrees |
| `AziM` | Azimuth of Moon | degrees |
| `ZendO` | Zenith distance of object: `90 − AltO` | degrees |
| `RS` | Angular distance from object to Sun | degrees |
| `RM` | Angular distance from object to Moon | degrees |
| `TAV` | Topocentric arcus visionis | degrees |
| `DAZ` | Azimuth difference (object − Sun) | degrees |
| `ARCL` | Arc length (longitude separation) | degrees |
| `GeoARCVact` | Geocentric arcus visionis (actual) | degrees |

### Vision & Magnitudes

| Variable | Meaning |
|----------|---------|
| `C1` | Schaefer threshold constant (scotopic: 10⁻⁹·⁸, photopic: 10⁻⁸·³⁵) |
| `C2` | Schaefer background constant (scotopic: 10⁻¹·⁹, photopic: 10⁻⁵·⁹) |
| `Th` | Threshold luminance for detection (nL) |
| `M0` | Reference magnitude constant |
| `MS` | Sun magnitude (−26.74) |
| `MM` | Moon magnitude (phase-dependent) |
| `MagnO` | Object apparent magnitude |
| `Wi` | Lunar crescent width | arcminutes |
| `qYal`, `q` | Yallop visibility criterion | dimensionless |
| `FS` | Forward scattering function | dimensionless |
| `CorrFactor1`, `CorrFactor2` | Observer age/acuity corrections | dimensionless |

---

## Time Variables

| Variable | Meaning | Units |
|----------|---------|-------|
| `tjd` | Julian Day in Terrestrial Time (TT/ET) | days |
| `tjd_ut` | Julian Day in Universal Time (UT1) | days |
| `t` | Working time (often = tjd or tjd − dt) | days |
| `T` | Julian centuries from J2000: `(tjd − 2451545.0)/36525` | centuries |
| `dt` | Time difference (light-time, or iteration step) | days |
| `teval` | Time of last cached evaluation (`pdp->teval`) | days |
| `tseg0`, `tseg1` | Current polynomial segment start/end | days |
| `tfstart`, `tfend` | Ephemeris file time range | days |
| `dseg` | Duration of one polynomial segment | days |
| `telem` | Epoch of orbital elements | days |
| `speed_intv` | Finite-difference interval for speed | days |
| `deltat` | Delta-T (TT − UT) | days or seconds |
| `sidt` | Sidereal time | hours or degrees |

---

## Structure Pointers

| Variable | Type | Points To |
|----------|------|-----------|
| `pdp` | `plan_data *` | `swed.pldat[ipl]` — current planet |
| `pedp` | `plan_data *` | `swed.pldat[SEI_EARTH]` — Earth |
| `psdp` | `plan_data *` | `swed.pldat[SEI_SUNBARY]` — barycentric Sun |
| `pmdp` | `plan_data *` | `swed.pldat[SEI_MOON]` — Moon |
| `pebdp` | `plan_data *` | `swed.pldat[SEI_EMB]` — Earth-Moon barycenter |
| `ndp` | `node_data *` | `swed.nddat[ipl]` — node/apside |
| `fdp` | `file_data *` | `swed.fidat[ifno]` — ephemeris file |
| `oe` | `epsilon *` | `swed.oec` or `swed.oec2000` — obliquity |
| `sip` | `sid_data *` | `swed.sidd` — sidereal config |

---

## Flags & Control

| Variable | Meaning |
|----------|---------|
| `ipl` | Planet/body ID (SE_SUN, SE_MOON, etc.) |
| `ipli` | Internal planet index for ephemeris files |
| `iflag` | Calculation control bitmask (SEFLG_*) |
| `epheflag` | Which ephemeris: `SEFLG_JPLEPH`, `SEFLG_SWIEPH`, `SEFLG_MOSEPH` |
| `ifno` | Ephemeris file number (SEI_FILE_PLANET, etc.) |
| `iephe` | Cached ephemeris type (in `pdp->iephe`) |
| `xflgs` | Cached transformation flags (in `pdp->xflgs`) |
| `retc` | Return code: OK (0), ERR (−1), NOT_AVAILABLE (−2), BEYOND_EPH_LIMITS (−3) |
| `niter` | Number of light-time iterations (0 for Moshier, 1 for JPL/Swiss) |
| `iseg` | Polynomial segment index within ephemeris file |
| `ncoe` | Number of Chebyshev coefficients per component |
| `neval` | Number of coefficients actually evaluated (≤ ncoe) |

---

## Physical Constants

| Variable | Value | Meaning |
|----------|-------|---------|
| `AUNIT` | 1.496×10¹¹ m | Astronomical Unit |
| `CLIGHT` | 2.998×10⁸ m/s | Speed of light |
| `HELGRAVCONST` | 1.327×10²⁰ m³/s² | GM_sun |
| `GEOGCONST` | 3.986×10¹⁴ m³/s² | GM_earth |
| `EARTH_RADIUS` | 6378136.6 m | Equatorial radius |
| `EARTH_OBLATENESS` | 1/298.256 | Flattening |
| `EARTH_MOON_MRAT` | ~81.3 | Earth/Moon mass ratio |
| `J2000` | 2451545.0 | J2000.0 epoch (JD) |
| `STR` | 4.848×10⁻⁶ | Arcseconds to radians |
| `DSUN`, `RSUN` | — | Sun diameter/radius (AU) |
| `DMOON`, `RMOON` | — | Moon diameter/radius (AU) |
| `DEARTH`, `REARTH` | — | Earth diameter/radius (AU) |

---

## Naming Patterns

| Pattern | Meaning | Examples |
|---------|---------|---------|
| `x` + body letter | Position of body | `xs` (Sun), `xm` (Moon), `xe` (Earth) |
| `sin` + name | Precomputed sine | `sinf1`, `sine`, `sinincl` |
| `cos` + name | Precomputed cosine | `cosf1`, `cose`, `cosnode` |
| `tan` + name | Precomputed tangent | `tane`, `tanfi`, `tant` |
| `d` + name | Diameter or delta | `dm` (distance Moon), `dsm` (distance Sun–Moon) |
| `r` + name | Radius or distance | `rmoon`, `rxy`, `rxyz` |
| `k` + letter | Extinction coefficient | `kR` (Rayleigh), `kW` (water) |
| `B` + name | Sky brightness (nL) | `Bday`, `Btwi`, `Bm` |
| `Alt` + letter | Altitude angle | `AltS` (Sun), `AltO` (object) |
| `Azi` + letter | Azimuth angle | `AziS`, `AziM` |
| `p` + `dp` | Planet data pointer | `pdp` (planet), `pedp` (Earth), `psdp` (Sun) |
| `xx` + suffix | Working position variant | `xx0` (original), `xxsv` (saved), `xxsp` (speed correction) |
| `i` + name | Integer index or ID | `ipl` (planet), `iflag` (flags), `ifno` (file number) |
| `n` + name | Count | `niter`, `ncoe`, `neval` |
| `t` + name | Time | `tjd`, `teval`, `tseg0` |
