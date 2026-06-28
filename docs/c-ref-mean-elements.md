# C Reference: Mean Node & Mean Apogee — swemmoon.c / sweph.c

Porting reference for `swi_mean_node`, `swi_mean_apog`, their correction tables, and the
`app_pos_etc_mean` pipeline function. Read this instead of the C source.

---

## Function Map

| C function | Location | Notes |
|---|---|---|
| `swi_mean_node` | swemmoon.c:1493–1534 | Mean lunar node, ecliptic of date |
| `swi_mean_apog` | swemmoon.c:1564–1624 | Mean lunar apogee (Dark Moon / Lilith) |
| `mean_elements` | swemmoon.c:1763–1818 | Moshier mean orbital elements (shared) |
| `corr_mean_node` | swemmoon.c:1470–1486 | Piecewise-linear node correction |
| `corr_mean_apog` | swemmoon.c:1540–1556 | Piecewise-linear apogee correction |
| `mean_node_corr[]` | swemmoon.c:725–763 | 304-entry correction table (degrees) |
| `mean_apsis_corr[]` | swemmoon.c:767–805 | 304-entry correction table (degrees) |
| `app_pos_etc_mean` | sweph.c:4310–4359 | Frame transforms for mean elements |
| `app_pos_rest` | sweph.c:2777–2859 | Nutation, ecliptic, polar, units (shared) |

---

## Constants

| Constant | Value | Source | Meaning |
|---|---|---|---|
| `J2000` | `2451545.0` | sweph.h:67 | J2000.0 Julian day |
| `STR` | `4.8481368110953599359e-6` | sweph.h:663 | Radians per arcsecond |
| `DEGTORAD` | `0.0174532925199433` | sweodef.h:262 | π / 180 |
| `RADTODEG` | `57.2957795130823` | sweodef.h:263 | 180 / π |
| `PI` | `M_PI` | sweph.h:126 | 3.14159… |
| `AUNIT` | `1.49597870700e+11` | sweph.h:273 | AU in metres (DE431 value) |
| `MOON_MEAN_DIST` | `384400000.0` | sweph.h:260 | Mean lunar distance in metres (AA 1996 F2) |
| `MOON_MEAN_INCL` | `5.1453964` | sweph.h:261 | Mean lunar inclination in degrees (AA 1996 D2) |
| `MOON_MEAN_ECC` | `0.054900489` | sweph.h:262 | Mean lunar eccentricity (AA 1996 F2) |
| `MOSHNDEPH_START` | `-3100015.5` | sweph.h:225 | Start of Moshier node ephemeris (JD) |
| `MOSHNDEPH_END` | `8000016.5` | sweph.h:226 | End of Moshier node ephemeris (JD) |
| `MEAN_NODE_SPEED_INTV` | `0.001` | sweph.h:300 | Step size for numerical speed (days) |
| `JPL_DE431_START` | `-3027215.5` | sweph.h:233 | DE431 start (guard for correction tables) |
| `JPL_DE431_END` | `7930192.5` | sweph.h:234 | DE431 end (guard for correction tables) |
| `CORR_MNODE_JD_T0GREG` | `-3063616.5` | swemmoon.c:1466 | Correction table epoch: 1 Jan −13100 Greg |
| `CORR_MAPOG_JD_T0GREG` | `-3063616.5` | swemmoon.c:1536 | Same epoch for apogee table |

---

## 1. `mean_elements()` — Moshier Mean Orbital Elements

**swemmoon.c:1763–1818**

Sets global state variables. Called with `T` and `T2` already set by the caller.

```
Input:  T  (global) = (J - J2000) / 36525.0
        T2 (global) = T * T

fracT = T mod 1  (fractional part only, for high-precision linear term)

M (mean anomaly of sun, arcsec — Laskar):
  M = mods3600(129600000.0 * fracT - 3418.961646 * T + 1287104.76154)
  M += ((((((((
      1.62e-20 * T
    - 1.0390e-17) * T
    - 3.83508e-15) * T
    + 4.237343e-13) * T
    + 8.8555011e-11) * T
    - 4.77258489e-8) * T
    - 1.1297037031e-5) * T
    + 1.4732069041e-4) * T
    - 0.552891801772) * T2

NF (mean distance of moon from ascending node = F, arcsec):
  NF = mods3600(1739232000.0 * fracT + 295263.0983 * T
                - 2.079419901760e-01 * T + 335779.55755)
  NF += ((z[2]*T + z[1])*T + z[0]) * T2

MP (mean anomaly of moon = l, arcsec):
  MP = mods3600(1717200000.0 * fracT + 715923.4728 * T
                - 2.035946368532e-01 * T + 485868.28096)
  MP += ((z[5]*T + z[4])*T + z[3]) * T2

D (mean elongation of moon, arcsec):
  D = mods3600(1601856000.0 * fracT + 1105601.4603 * T
               + 3.962893294503e-01 * T + 1072260.73512)
  D += ((z[8]*T + z[7])*T + z[6]) * T2

SWELP (mean longitude of moon, referred to mean ecliptic/equinox of date, arcsec):
  SWELP = mods3600(1731456000.0 * fracT + 1108372.83264 * T
                   - 6.784914260953e-01 * T + 785939.95571)
  SWELP += ((z[11]*T + z[10])*T + z[9]) * T2
```

The `z[]` coefficients used (DE404 fit, non-MOSH_MOON_200 branch — this is the active branch):

```
z[0]  = -1.312045233711e+01   /* NF (F), t^2 arcsec */
z[1]  = -1.138215912580e-03   /* NF (F), t^3 */
z[2]  = -9.646018347184e-06   /* NF (F), t^4 */
z[3]  =  3.146734198839e+01   /* MP (l), t^2 arcsec */
z[4]  =  4.768357585780e-02   /* MP (l), t^3 */
z[5]  = -3.421689790404e-04   /* MP (l), t^4 */
z[6]  = -6.847070905410e+00   /* D, t^2 arcsec */
z[7]  = -5.834100476561e-03   /* D, t^3 */
z[8]  = -2.905334122698e-04   /* D, t^4 */
z[9]  = -5.663161722088e+00   /* SWELP (L), t^2 arcsec */
z[10] =  5.722859298199e-03   /* SWELP (L), t^3 */
z[11] = -8.466472828815e-05   /* SWELP (L), t^4 */
```

`mods3600(x)` reduces x modulo 1296000 arcsec (one full circle):
```
mods3600(x) = x - 1296000.0 * floor(x / 1296000.0)
```
Result is always in [0, 1296000).

**Note**: The DE200 fit (MOSH_MOON_200 branch, z[] with 71 coefficients) is NOT active.
Only the DE404 fit above (12 z-coefficients) is compiled in.

---

## 2. `corr_mean_node()` — Node Correction Interpolation

**swemmoon.c:1470–1486**

Returns a correction in **degrees**.

```
Input:  J  (Julian day)

J0       = CORR_MNODE_JD_T0GREG  = -3063616.5
dayscty  = 36524.25               (Gregorian century in days)

if J < JPL_DE431_START  →  return 0.0
if J > JPL_DE431_END    →  return 0.0

dJ    = J - J0
i     = floor(dJ / dayscty)        // century index (lower bound)
dfrac = (dJ - i * dayscty) / dayscty

dcor0 = mean_node_corr[i]
dcor1 = mean_node_corr[i + 1]
dcor  = dcor0 + dfrac * (dcor1 - dcor0)   // linear interpolation

return dcor   // in degrees
```

---

## 3. `corr_mean_apog()` — Apogee Correction Interpolation

**swemmoon.c:1540–1556**

Returns a correction in **degrees**. Identical structure to `corr_mean_node` but uses
`mean_apsis_corr[]` and `CORR_MAPOG_JD_T0GREG = -3063616.5`.

```
Input:  J  (Julian day)

J0       = CORR_MAPOG_JD_T0GREG  = -3063616.5
dayscty  = 36524.25

if J < JPL_DE431_START  →  return 0.0
if J > JPL_DE431_END    →  return 0.0

dJ    = J - J0
i     = floor(dJ / dayscty)
dfrac = (dJ - i * dayscty) / dayscty

dcor0 = mean_apsis_corr[i]
dcor1 = mean_apsis_corr[i + 1]
dcor  = dcor0 + dfrac * (dcor1 - dcor0)

return dcor   // in degrees
```

---

## 4. `swi_mean_node()` — Mean Lunar Node

**swemmoon.c:1493–1534**

Returns polar coordinates of the mean lunar node in the ecliptic of date.

```
Input:  J    (Julian day)
        pol  (output array, at least 3 elements)
        serr (error string buffer)

T  = (J - J2000) / 36525.0
T2 = T * T
T3 = T * T2
T4 = T2 * T2

if J < MOSHNDEPH_START or J > MOSHNDEPH_END:
    write error to serr; return ERR

mean_elements()    // sets SWELP, NF (arcsec, in [0, 1296000))

dcor = corr_mean_node(J) * 3600.0   // convert degrees → arcseconds

pol[0] = swi_mod2PI((SWELP - NF - dcor) * STR)   // longitude in radians
pol[1] = 0.0                                       // latitude = 0
pol[2] = MOON_MEAN_DIST / AUNIT                    // distance in AU

return OK
```

`swi_mod2PI(x)` reduces angle to [0, 2π).

**Note**: T3, T4 are set but not used in `swi_mean_node` itself; they are available for
`mean_elements()` if needed (the DE404 z[] corrections only use T2).

---

## 5. `swi_mean_apog()` — Mean Lunar Apogee (Dark Moon / Lilith)

**swemmoon.c:1564–1624**

Returns cartesian ecliptic-of-date position after projecting the apogee onto the ecliptic
through the mean orbital inclination.

```
Input:  J    (Julian day)
        pol  (output array, at least 3 elements)
        serr (error string buffer)

T  = (J - J2000) / 36525.0
T2 = T * T
T3 = T * T2
T4 = T2 * T2

if J < MOSHNDEPH_START or J > MOSHNDEPH_END:
    write error to serr; return ERR

mean_elements()   // sets SWELP, NF, MP (arcsec)

--- Step 1: raw apogee longitude (180° from perigee) ---
pol[0] = swi_mod2PI((SWELP - MP) * STR + PI)   // radians, ecliptic of date
pol[1] = 0
pol[2] = MOON_MEAN_DIST * (1 + MOON_MEAN_ECC) / AUNIT   // apogee distance AU

--- Step 2: apply apogee correction ---
dcor  = corr_mean_apog(J) * DEGTORAD             // degrees → radians
pol[0] = swi_mod2PI(pol[0] - dcor)

--- Step 3: project onto ecliptic plane through inclination ---
// compute corrected node
node  = (SWELP - NF) * STR                      // raw node in radians
dcor  = corr_mean_node(J) * DEGTORAD            // node correction, radians
node  = swi_mod2PI(node - dcor)

// shift apogee argument of latitude (subtract node)
pol[0] = swi_mod2PI(pol[0] - node)

// convert to cartesian, tilt by inclination, convert back
swi_polcart(pol, pol)                            // polar → 3D cartesian
swi_coortrf(pol, pol, -MOON_MEAN_INCL * DEGTORAD)  // rotate by -inclination
swi_cartpol(pol, pol)                            // 3D cartesian → polar

// re-add node to recover ecliptic longitude
pol[0] = swi_mod2PI(pol[0] + node)

return OK
```

**Output format**: `pol[0]` = ecliptic longitude (radians), `pol[1]` = ecliptic latitude
(radians, non-zero after inclination projection), `pol[2]` = distance (AU).

---

## 6. Pipeline: `swecalc()` calling mean elements

**sweph.c:860–928**

Before `app_pos_etc_mean` is called, `swecalc` computes position and numerical speed:

### Mean Node (SE_MEAN_NODE → SEI_MEAN_NODE)

```
if HELCTR or BARYCTR: zero all 24 xreturn entries, return iflag  // not allowed

ndp  = &swed.nddat[SEI_MEAN_NODE]
xp   = ndp->xreturn        // output array [24]
xp2  = ndp->x              // work array  [6]

swi_mean_node(tjd,                       xp2,   serr)  // position  → xp2[0..2]
swi_mean_node(tjd - MEAN_NODE_SPEED_INTV, xp2+3, serr)  // prior pos → xp2[3..5]

// numerical differentiation for longitude speed
xp2[3] = swe_difrad2n(xp2[0], xp2[3]) / MEAN_NODE_SPEED_INTV
xp2[4] = 0   // latitude speed = 0
xp2[5] = 0   // radial speed   = 0

ndp->teval = tjd
ndp->xflgs = -1   // force re-computation in app_pos_etc_mean

app_pos_etc_mean(SEI_MEAN_NODE, iflag, serr)

// post-process: zero latitude/z components to suppress float noise
// (only if NOT sidereal and NOT J2000)
if !(SEFLG_SIDEREAL) and !(SEFLG_J2000):
    ndp->xreturn[1]  = 0.0   // ecl. latitude
    ndp->xreturn[4]  = 0.0   // ecl. lat. speed
    ndp->xreturn[5]  = 0.0   // radial speed
    ndp->xreturn[8]  = 0.0   // z (ecl. cartesian)
    ndp->xreturn[11] = 0.0   // z speed
```

### Mean Apogee (SE_MEAN_APOG → SEI_MEAN_APOG)

```
if HELCTR or BARYCTR: zero all 24 xreturn entries, return iflag  // not allowed

ndp  = &swed.nddat[SEI_MEAN_APOG]
xp   = ndp->xreturn
xp2  = ndp->x

swi_mean_apog(tjd,                       xp2,   serr)  // position
swi_mean_apog(tjd - MEAN_NODE_SPEED_INTV, xp2+3, serr)  // prior pos

// numerical differentiation for lon and lat speeds
for i in 0..1:
    xp2[3+i] = swe_difrad2n(xp2[i], xp2[3+i]) / MEAN_NODE_SPEED_INTV
xp2[5] = 0   // radial speed = 0

ndp->teval = tjd
ndp->xflgs = -1

app_pos_etc_mean(SEI_MEAN_APOG, iflag, serr)

// post-process: radial speed always 0 (apogee distance from mean elements only)
ndp->xreturn[5] = 0.0
```

---

## 7. `app_pos_etc_mean()` — Frame Transforms for Mean Elements

**sweph.c:4310–4359**

Converts `pdp->x[0..5]` (polar ecliptic of date, radians) into the full 24-element
`xreturn` array through frame transforms.

```
Input:
  ipl    = SEI_MEAN_NODE or SEI_MEAN_APOG
  iflag  = caller's flag word
  pdp    = &swed.nddat[ipl]

--- cache check ---
flg1 = iflag  & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ
flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ
if flg1 == flg2:
    pdp->xflgs = iflag
    pdp->iephe = iflag & SEFLG_EPHMASK
    return OK

--- copy input ---
xx[0..5] = pdp->x[0..5]   // polar ecliptic of date (lon, lat, dist, d_lon, d_lat, d_dist)

--- polar → cartesian (ecliptic of date, with speed) ---
swi_polcart_sp(xx, xx)

--- ecliptic cartesian → equatorial cartesian (rotate by −obliquity of date) ---
swi_coortrf2(xx,   xx,   -swed.oec.seps, swed.oec.ceps)
swi_coortrf2(xx+3, xx+3, -swed.oec.seps, swed.oec.ceps)

--- zero speed components if not requested ---
if !(iflag & SEFLG_SPEED):
    xx[3] = xx[4] = xx[5] = 0

--- save J2000 copy for sidereal modes that need it ---
if (SEFLG_SIDEREAL and sidd.sid_mode & SE_SIDBIT_ECL_T0)
   or (sidd.sid_mode & SE_SIDBIT_SSY_PLANE):
    xxsv[0..5] = xx[0..5]
    if pdp->teval != J2000:
        swi_precess(xxsv, pdp->teval, iflag, J_TO_J2000)
        if SEFLG_SPEED:
            swi_precess_speed(xxsv, pdp->teval, iflag, J_TO_J2000)

--- J2000 vs date ---
if SEFLG_J2000:
    // no precession → precess equatorial coords back to J2000
    swi_precess(xx, pdp->teval, iflag, J_TO_J2000)
    if SEFLG_SPEED:
        swi_precess_speed(xx, pdp->teval, iflag, J_TO_J2000)
    oe = &swed.oec2000   // obliquity of J2000
else:
    oe = &swed.oec       // obliquity of date

return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr)
```

**Key difference from planet path**: Mean elements have no ICRS frame bias step and no
light-time, aberration, or gravitational deflection. The input is already in ecliptic of
date (from `swi_mean_node`/`swi_mean_apog`); the function just re-frames it.

---

## 8. `app_pos_rest()` — Shared Nutation/Ecliptic/Polar Stage

**sweph.c:2777–2859**

Called by `app_pos_etc_mean` (and other `app_pos_etc_*` functions). `xx` on entry is
equatorial cartesian; `x2000` is the J2000 equatorial cartesian (for sidereal modes).

```
--- nutation ---
if !(SEFLG_NONUT):
    swi_nutate(xx, iflag, FALSE)   // in-place nutation of equatorial cartesian

--- save equatorial cartesian ---
pdp->xreturn[18..23] = xx[0..5]

--- equatorial → ecliptic cartesian ---
swi_coortrf2(xx,   xx,   oe->seps, oe->ceps)    // rotate by +obliquity
if SEFLG_SPEED:
    swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps)

--- apply nutation rotation (ecliptic nutation) ---
if !(SEFLG_NONUT):
    swi_coortrf2(xx,   xx,   swed.nut.snut, swed.nut.cnut)
    if SEFLG_SPEED:
        swi_coortrf2(xx+3, xx+3, swed.nut.snut, swed.nut.cnut)

--- save ecliptic cartesian ---
pdp->xreturn[6..11] = xx[0..5]

--- sidereal handling (SEFLG_SIDEREAL) ---
if SEFLG_SIDEREAL:
    if sid_mode & SE_SIDBIT_ECL_T0:
        swi_trop_ra2sid_lon(x2000, xreturn+6, xreturn+18, iflag)
    elif sid_mode & SE_SIDBIT_SSY_PLANE:
        swi_trop_ra2sid_lon_sosy(x2000, xreturn+6, iflag)
    else:  // traditional: subtract ayanamsa
        swi_cartpol_sp(xreturn+6, xreturn)
        save xreturn[0..23] to xxsv
        swi_get_ayanamsa_with_speed(pdp->teval, iflag, daya, serr)
        restore xreturn[0..23] from xxsv
        xreturn[0] -= daya[0] * DEGTORAD   // subtract ayanamsa from longitude
        xreturn[3] -= daya[1] * DEGTORAD   // subtract ayanamsa speed
        swi_polcart_sp(xreturn, xreturn+6) // back to cartesian

--- convert to polar ---
swi_cartpol_sp(pdp->xreturn+18, pdp->xreturn+12)  // equatorial cart → polar
swi_cartpol_sp(pdp->xreturn+6,  pdp->xreturn)     // ecliptic  cart → polar

--- radians → degrees ---
for i in 0, 1:
    xreturn[i]    *= RADTODEG   // ecl. lon, lat
    xreturn[i+3]  *= RADTODEG   // ecl. lon speed, lat speed
    xreturn[i+12] *= RADTODEG   // equ. RA, decl
    xreturn[i+15] *= RADTODEG   // equ. RA speed, decl speed
// distances and radial speeds (indices 2, 5, 14, 17) remain in AU and AU/day

--- save ---
pdp->xflgs = iflag
pdp->iephe = iflag & SEFLG_EPHMASK
return OK
```

### `xreturn[24]` Layout

```
[0]  = ecliptic longitude   (degrees)
[1]  = ecliptic latitude    (degrees)
[2]  = distance             (AU)
[3]  = d(longitude)/dt      (degrees/day)
[4]  = d(latitude)/dt       (degrees/day)
[5]  = d(distance)/dt       (AU/day)
[6..8]  = ecliptic cartesian x, y, z
[9..11] = ecliptic cartesian dx, dy, dz (speed)
[12] = right ascension      (degrees)
[13] = declination          (degrees)
[14] = distance             (AU)
[15] = d(RA)/dt             (degrees/day)
[16] = d(decl)/dt           (degrees/day)
[17] = d(distance)/dt       (AU/day)
[18..20] = equatorial cartesian x, y, z
[21..23] = equatorial cartesian dx, dy, dz (speed)
```

---

## 9. `mean_node_corr[]` — Complete Correction Table

**swemmoon.c:725–763**

304 entries, one per Gregorian century, starting at `CORR_MNODE_JD_T0GREG = -3063616.5`
(1 January −13100 Greg.). Values in **degrees**. Step = 36524.25 days = 1 Gregorian century.

Entries covering years 0 to 3000 (indices ~131–160) are zeroed because the corrections
are defined to be negligible in that range. A `#if 0` block with fitted data values for
that period exists in source but is compiled out.

```c
static const double mean_node_corr[] = {
/* index 0: year -13100 */
-2.56,
/* indices 1–130: years -13000 to -100 */
-2.473, -2.392347, -2.316425, -2.239639, -2.167764, -2.095100, -2.024810, -1.957622, -1.890097, -1.826389,
-1.763335, -1.701047, -1.643016, -1.584186, -1.527309, -1.473352, -1.418917, -1.367736, -1.317202, -1.267269,
-1.221121, -1.174218, -1.128862, -1.086214, -1.042998, -1.002491, -0.962635, -0.923176, -0.887191, -0.850403,
-0.814929, -0.782117, -0.748462, -0.717241, -0.686598, -0.656013, -0.628726, -0.600460, -0.573219, -0.548634,
-0.522931, -0.499285, -0.476273, -0.452978, -0.432663, -0.411386, -0.390788, -0.372825, -0.353681, -0.336230,
-0.319520, -0.302343, -0.287794, -0.272262, -0.257166, -0.244534, -0.230635, -0.218126, -0.206365, -0.194000,
-0.183876, -0.172782, -0.161877, -0.153254, -0.143371, -0.134501, -0.126552, -0.117932, -0.111199, -0.103716,
-0.096160, -0.090718, -0.084046, -0.078007, -0.072959, -0.067235, -0.062990, -0.058102, -0.053070, -0.049786,
-0.045381, -0.041317, -0.038165, -0.034501, -0.031871, -0.028844, -0.025701, -0.024018, -0.021427, -0.018881,
-0.017291, -0.015186, -0.013755, -0.012098, -0.010261, -0.009688, -0.008218, -0.006670, -0.005979, -0.004756,
-0.003991, -0.002996, -0.001974, -0.001975, -0.001213, -0.000377, -0.000356, 5.779e-05, 0.000378, 0.000710,
0.001092, 0.000767, 0.000985, 0.001443, 0.001069, 0.001141, 0.001321, 0.001462, 0.001695, 0.001319,
0.001567, 0.001873, 0.001376, 0.001336, 0.001347, 0.001330, 0.001256, 0.000813, 0.000946, 0.001079,
/* indices 131–160: years 0 to 3000 — set to 0 (fitted data compiled out) */
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
/* indices 161–303: years 3100 to 17200 */
-0.000364, -0.000452, -0.001091, -0.001159, -0.001136, -0.001798, -0.002249, -0.002622, -0.002990, -0.003555,
-0.004425, -0.004758, -0.005134, -0.006065, -0.006839, -0.007474, -0.008283, -0.009411, -0.010786, -0.011810,
-0.012989, -0.014825, -0.016426, -0.017922, -0.019774, -0.021881, -0.024194, -0.026190, -0.028440, -0.031285,
-0.033817, -0.036318, -0.039212, -0.042456, -0.045799, -0.048994, -0.052710, -0.056948, -0.061017, -0.065181,
-0.069843, -0.074922, -0.079976, -0.085052, -0.090755, -0.096840, -0.102797, -0.108939, -0.115568, -0.122636,
-0.129593, -0.136683, -0.144641, -0.152825, -0.161044, -0.169758, -0.178916, -0.188712, -0.198401, -0.208312,
-0.219395, -0.230407, -0.241577, -0.253508, -0.265640, -0.278556, -0.291330, -0.304353, -0.318815, -0.332882,
-0.347316, -0.362895, -0.378421, -0.395061, -0.411748, -0.428666, -0.447477, -0.465636, -0.484277, -0.504600,
-0.524405, -0.545533, -0.567020, -0.588404, -0.612099, -0.634965, -0.658262, -0.683866, -0.708526, -0.734719,
-0.761800, -0.788562, -0.818092, -0.846885, -0.876177, -0.908385, -0.939371, -0.972027, -1.006149, -1.039634,
-1.076135, -1.112156, -1.148490, -1.188312, -1.226761, -1.266821, -1.309156, -1.350583, -1.395223, -1.440028,
-1.485047, -1.534104, -1.582023, -1.631506, -1.684031, -1.735687, -1.790421, -1.846039, -1.901951, -1.961872,
-2.021179, -2.081987, -2.146259, -2.210031, -2.276609, -2.344904, -2.413795, -2.486559, -2.559564, -2.634215,
-2.712692, -2.791289, -2.872533, -2.956217, -3.040965, -3.129234, -3.218545, -3.309805, -3.404827, -3.5008,
-3.601, -3.7, -3.8,
};
```

Total entries: **304**. Valid index range for interpolation (accessing `[i]` and `[i+1]`):
max `i` accessed at JPL_DE431_END ≈ 301 → array must be at least 303 entries long.

---

## 10. `mean_apsis_corr[]` — Complete Apogee Correction Table

**swemmoon.c:767–805**

304 entries, same epoch and step as the node table. Values in **degrees**.
Same `#if 0` block structure — the 30 entries covering years 0–3000 are compiled out
and replaced with zeros.

```c
static const double mean_apsis_corr[] = {
/* index 0: year -13100 */
7.525,
/* indices 1–130: years -13000 to -100 */
7.290, 7.057295, 6.830813, 6.611723, 6.396775, 6.189569, 5.985968, 5.788342, 5.597304, 5.410167,
5.229946, 5.053389, 4.882187, 4.716494, 4.553532, 4.396734, 4.243718, 4.094282, 3.950865, 3.810366,
3.674978, 3.543284, 3.414270, 3.290526, 3.168775, 3.050904, 2.937541, 2.826189, 2.719822, 2.616193,
2.515431, 2.419193, 2.323782, 2.232545, 2.143635, 2.056803, 1.974913, 1.893874, 1.816201, 1.741957,
1.668083, 1.598335, 1.529645, 1.463016, 1.399693, 1.336905, 1.278097, 1.220965, 1.165092, 1.113071,
1.060858, 1.011007, 0.963701, 0.916523, 0.872887, 0.829596, 0.788486, 0.750017, 0.711177, 0.675589,
0.640303, 0.605303, 0.573490, 0.541113, 0.511482, 0.483159, 0.455210, 0.430305, 0.404643, 0.380782,
0.358524, 0.335405, 0.315244, 0.295131, 0.275766, 0.259223, 0.241586, 0.225890, 0.210404, 0.194775,
0.181573, 0.167246, 0.154514, 0.143435, 0.131131, 0.121648, 0.111835, 0.102474, 0.094284, 0.085204,
0.078240, 0.070697, 0.063696, 0.058894, 0.052390, 0.047632, 0.043129, 0.037823, 0.034143, 0.029188,
0.025648, 0.021972, 0.018348, 0.017127, 0.013989, 0.011967, 0.011003, 0.007865, 0.007033, 0.005574,
0.004060, 0.003699, 0.002465, 0.002889, 0.002144, 0.001018, 0.001757, -9.67e-05, -0.000734, -0.000392,
-0.001546, -0.000863, -0.001266, -0.000933, -0.000503, -0.001304, 0.000238, -0.000507, -0.000897, 0.000647,
/* indices 131–160: years 0 to 3000 — set to 0 (fitted data compiled out) */
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
/* indices 161–303: years 3100 to 17200 */
0.000514, 0.000683, 0.002228, 0.001974, 0.003485, 0.004280, 0.005409, 0.007468, 0.007938, 0.011012,
0.012525, 0.013757, 0.016757, 0.017932, 0.020780, 0.023416, 0.026386, 0.030428, 0.033512, 0.038789,
0.043126, 0.047778, 0.054175, 0.058891, 0.065878, 0.072345, 0.079668, 0.088238, 0.095307, 0.104873,
0.113533, 0.122336, 0.133205, 0.142922, 0.154871, 0.166488, 0.179234, 0.193928, 0.207262, 0.223089,
0.238736, 0.254907, 0.273232, 0.291085, 0.311046, 0.331025, 0.351955, 0.374422, 0.396341, 0.420772,
0.444867, 0.469984, 0.497448, 0.524717, 0.554752, 0.584581, 0.616272, 0.649744, 0.682947, 0.719405,
0.755834, 0.793780, 0.833875, 0.873893, 0.917340, 0.960429, 1.005471, 1.052384, 1.099317, 1.149508,
1.200130, 1.253038, 1.307672, 1.363480, 1.422592, 1.481900, 1.544111, 1.607982, 1.672954, 1.741025,
1.809727, 1.882038, 1.955243, 2.029956, 2.108428, 2.186805, 2.268697, 2.352071, 2.437370, 2.525903,
2.615415, 2.709082, 2.804198, 2.901704, 3.002606, 3.104412, 3.210406, 3.317733, 3.428386, 3.541634,
3.656634, 3.775988, 3.896306, 4.020480, 4.146814, 4.275356, 4.408257, 4.542282, 4.681174, 4.822524,
4.966424, 5.114948, 5.264973, 5.419906, 5.577056, 5.737688, 5.902347, 6.069138, 6.241065, 6.415155,
6.593317, 6.774853, 6.959322, 7.148845, 7.340334, 7.537156, 7.737358, 7.940882, 8.149932, 8.361576,
8.579150, 8.799591, 9.024378, 9.254584, 9.487362, 9.726535, 9.968784, 10.216089, 10.467716, 10.725293,
10.986, 11.25, 11.52,
};
```

Total entries: **304**.

---

## 11. Differences from Planet Path (`app_pos_etc_plan`)

| Step | Mean elements | Planets |
|---|---|---|
| Light-time correction | None | Yes |
| Aberration | None | Yes |
| Gravitational deflection | None | Yes |
| ICRS frame bias | None | Yes (DE403+) |
| Input position source | `swi_mean_node`/`swi_mean_apog` | Ephemeris file / Moshier |
| Obliquity applied at start | Yes (`swi_coortrf2` in `app_pos_etc_mean`) | Yes (same) |
| Nutation (via `app_pos_rest`) | Yes (same code path) | Yes |
| Sidereal modes | Yes (same code path) | Yes |

---

## 12. Coordinate Helper Summary

These are used directly in the mean element pipeline. No need to port them — they exist
in the Rust codebase already (or should be ported as part of the shared coordinate module).

| Function | Direction | Notes |
|---|---|---|
| `swi_polcart_sp(src, dst)` | polar (lon, lat, dist, d_lon, d_lat, d_dist) → cartesian (x,y,z, dx,dy,dz) | In-place safe (src == dst) |
| `swi_cartpol_sp(src, dst)` | cartesian → polar | In-place safe |
| `swi_polcart(src, dst)` | polar (3 elem) → cartesian (3 elem) | No speed |
| `swi_cartpol(src, dst)` | cartesian (3 elem) → polar (3 elem) | No speed |
| `swi_coortrf2(src, dst, sineps, coseps)` | rotate about x-axis by ±obliquity | Negative angle = ecl→equ, positive = equ→ecl |
| `swi_coortrf(src, dst, eps_rad)` | same but takes angle in radians directly | |
| `swi_nutate(xx, iflag, back)` | apply nutation to equatorial cartesian | `back=FALSE` = forward |
| `swi_mod2PI(x)` | reduce radians to [0, 2π) | |
| `swe_difrad2n(a, b)` | signed angular difference a−b in radians, normalized to (−π, π] | Used for speed calculation |
