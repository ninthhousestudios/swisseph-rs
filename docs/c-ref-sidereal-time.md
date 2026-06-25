# C Reference: Sidereal Time — swephlib.c

Porting reference for Greenwich Apparent Sidereal Time (GAST). Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_sidtime` | swephlib.c:3580–3594 | Yes — public wrapper |
| `swe_sidtime0` | swephlib.c:3464–3556 | Yes — main dispatcher |
| `sidtime_long_term` | swephlib.c:3285–3324 | Yes — long-term model |
| `sidtime_non_polynomial_part` | swephlib.c:3413–3450 | Yes — 33-term Fourier correction |

## Public API Signatures (swephexp.h:956–957)

```c
/* GAST from pre-computed obliquity and nutation.
 * tjd_ut  = Julian Day in Universal Time (UT1)
 * eps     = true obliquity of ecliptic in DEGREES (eps_mean + deps)
 * nut     = nutation in longitude (dpsi) in DEGREES
 * returns: Greenwich Apparent Sidereal Time in decimal HOURS, range [0, 24)
 */
double swe_sidtime0(double tjd_ut, double eps, double nut);

/* GAST with automatic obliquity/nutation computation.
 * tjd_ut  = Julian Day in Universal Time (UT1)
 * returns: Greenwich Apparent Sidereal Time in decimal HOURS, range [0, 24)
 */
double swe_sidtime(double tjd_ut);
```

## Model Constants (swephexp.h:506, 541–545)

```c
#define SE_MODEL_SIDT           7    /* index into astro_models[] */

#define SEMOD_SIDT_IAU_1976         1
#define SEMOD_SIDT_IAU_2006         2
#define SEMOD_SIDT_IERS_CONV_2010   3
#define SEMOD_SIDT_LONGTERM         4
#define SEMOD_SIDT_DEFAULT          SEMOD_SIDT_LONGTERM  /* = 4 */
```

## Long-Term Boundary Constants (swephlib.c:3460–3463)

```c
#define SIDT_LTERM_T0    2396758.5   /* 1 Jan 1850 */
#define SIDT_LTERM_T1    2469807.5   /* 1 Jan 2050 */
#define SIDT_LTERM_OFS0  (0.000378172 / 15.0)   /* offset subtracted before T0, hours */
#define SIDT_LTERM_OFS1  (0.001385646 / 15.0)   /* offset subtracted after  T1, hours */
```

---

## swe_sidtime() — Public Wrapper (swephlib.c:3580–3594)

No model selection logic here. Always:

```
tjd_et = tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL)   // TT
eps    = swi_epsiln(tjd_et, 0) * RADTODEG            // mean obliquity in degrees
nutlo[0..1] = swi_nutation(tjd_et, 0, nutlo)         // dpsi, deps in radians
nutlo[0] *= RADTODEG   // dpsi → degrees
nutlo[1] *= RADTODEG   // deps → degrees
return swe_sidtime0(tjd_ut, eps + nutlo[1], nutlo[0])
//                         ^^^^^^^^^^^^^^^^  ^^^^^^^^
//                         true obliquity    dpsi (degrees)
```

Note: delta-T is called with `iflag = -1`, meaning "use default tidal acceleration for the
current ephemeris".

---

## swe_sidtime0() — Model Dispatcher (swephlib.c:3464–3556)

```
sidt_model = swed.astro_models[SE_MODEL_SIDT]
if sidt_model == 0: sidt_model = SEMOD_SIDT_DEFAULT  // = 4

if sidt_model == SEMOD_SIDT_LONGTERM (4):
    if tjd <= SIDT_LTERM_T0 OR tjd >= SIDT_LTERM_T1:
        gmst = sidtime_long_term(tjd, eps, nut)
        if tjd <= SIDT_LTERM_T0: gmst -= SIDT_LTERM_OFS0
        elif tjd >= SIDT_LTERM_T1: gmst -= SIDT_LTERM_OFS1
        if gmst >= 24: gmst -= 24
        if gmst < 0: gmst += 24
        return gmst           // ← early return, skip EoE below
    // else: fall through to IERS_CONV_2010 path (dates 1850–2050)

compute jd0 (JD at UT midnight) and secs (UT seconds since midnight):
    jd0 = floor(jd)
    secs = tjd - jd0
    if secs < 0.5:
        jd0 -= 0.5; secs += 0.5
    else:
        jd0 += 0.5; secs -= 0.5
    secs *= 86400.0
    tu = (jd0 - J2000) / 36525.0       // UT1 centuries from J2000

if sidt_model in {SEMOD_SIDT_IERS_CONV_2010, SEMOD_SIDT_LONGTERM}:
    → ERA-based model (see below)

elif sidt_model == SEMOD_SIDT_IAU_2006:
    → IAU 2006 polynomial model (see below)

else:   // SEMOD_SIDT_IAU_1976
    → IAU 1976 polynomial model (see below)

// After polynomial path, add Equation of the Equinoxes:
eqeq = 240.0 * nut * cos(eps * DEGTORAD)   // sidereal seconds of time
gmst = gmst + eqeq
gmst = gmst - 86400.0 * floor(gmst / 86400.0)   // modulo 1 sidereal day
gmst /= 3600.0    // seconds → hours
return gmst
```

**Key unit note**: For IAU 1976 and IAU 2006 non-ERA paths, `gmst` is accumulated in **sidereal seconds of time** throughout. For the ERA path, `gmst` is temporarily in **degrees** before converting at line 3510. The long-term path returns directly in **hours**.

**Equation of the Equinoxes derivation**:
- `nut` = dpsi in degrees; `eps` = true obliquity in degrees
- EoE in sidereal arcseconds = dpsi_arcsec × cos(ε) = (dpsi_deg × 3600) × cos(ε)
- EoE in sidereal seconds = dpsi_arcsec × cos(ε) / 15
- Combining: EoE_secs = dpsi_deg × 3600 × cos(ε) / 15 = dpsi_deg × 240 × cos(ε)
- Hence: `eqeq = 240.0 * nut * cos(eps * DEGTORAD)` ✓

---

## ERA-Based Path (IERS Conv 2010 / IAU 2006) (swephlib.c:3500–3510)

Used when `sidt_model` is `SEMOD_SIDT_IERS_CONV_2010` (3) or `SEMOD_SIDT_LONGTERM` (4) in
the 1850–2050 range.

Time variables:
```
jdrel = tjd - J2000                                         // UT days from J2000
tt = (tjd + swe_deltat_ex(tjd, -1, NULL) - J2000) / 36525.0  // TT centuries from J2000
```

Step 1 — Earth Rotation Angle (ERA) in degrees:
```
gmst = swe_degnorm((0.7790572732640 + 1.00273781191135448 * jdrel) * 360.0)
```

Step 2 — GMST correction polynomial (Capitaine et al. 2003), added in degrees:
```
gmst += (0.014506
      + tt * (4612.156534
      + tt * (1.3915817
      + tt * (-0.00000044
      + tt * (-0.000029956
      + tt * -0.0000000368))))) / 3600.0
```
(polynomial result is in arcseconds; /3600 converts to degrees)

Step 3 — Non-polynomial (Fourier) complementary terms in degrees:
```
dadd = sidtime_non_polynomial_part(tt)   // see below; returns degrees
gmst = swe_degnorm(gmst + dadd)
```

Step 4 — Convert to sidereal seconds:
```
gmst = gmst / 15.0 * 3600.0    // degrees → hours → seconds: deg * (1h/15°) * (3600s/1h)
```

Then EoE is added and result is returned as described in the dispatcher.

---

## IAU 2006 Polynomial Path (swephlib.c:3512–3518)

Used when `sidt_model == SEMOD_SIDT_IAU_2006` (2).

```
tt = (jd0 + swe_deltat_ex(jd0, -1, NULL) - J2000) / 36525.0   // TT centuries from J2000
```

GMST at UT midnight (sidereal seconds, Horner form):
```
gmst = (((-0.000000002454 * tt
         - 0.00000199708) * tt
         - 0.0000002926)  * tt
         + 0.092772110)   * tt * tt
      + 307.4771013 * (tt - tu)
      + 8640184.79447825 * tu
      + 24110.5493771
```

Mean sidereal days per solar day at tu (for UT interpolation):
```
msday = 1.0 + ((((-0.000000012270 * tt
                  - 0.00000798832) * tt
                  - 0.0000008778)  * tt
                  + 0.185544220)   * tt
                  + 8640184.79447825) / (86400.0 * 36525.0)
gmst += msday * secs    // secs = UT seconds since midnight
```

Note: `gmst` here is in sidereal seconds. Distinction from IERS path: this uses `jd0` (UT midnight) for TT, not the actual `tjd`.

---

## IAU 1976 Polynomial Path (swephlib.c:3520–3525)

Used when `sidt_model == SEMOD_SIDT_IAU_1976` (1) or as final else.

```
tu = (jd0 - J2000) / 36525.0    // UT1 centuries from J2000
```

GMST at UT midnight (sidereal seconds, Horner form):
```
gmst = ((-6.2e-6 * tu + 9.3104e-2) * tu + 8640184.812866) * tu + 24110.54841
```

Mean sidereal days per solar day (for UT interpolation):
```
msday = 1.0 + ((-1.86e-5 * tu + 0.186208) * tu + 8640184.812866) / (86400.0 * 36525.0)
gmst += msday * secs
```

Note: IAU 1976 does NOT use TT for the polynomial — `tu` is purely UT1-based.

---

## Long-Term Algorithm (swephlib.c:3285–3324)

Used for dates outside [SIDT_LTERM_T0, SIDT_LTERM_T1] (outside 1850–2050) when `SEMOD_SIDT_LONGTERM` is active.

Based on Simon et al. (1994), mean Earth longitude relative to J2000 mean equinox.

### Step-by-step

```
tjd_et  = tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL)        // TT
t       = (tjd_et - J2000) / 365250.0                      // Julian MILLENNIA from J2000
t2 = t*t;  t3 = t*t2
dlt     = AUNIT / CLIGHT / 86400.0                         // light-time: 1 AU in days
```

Mean longitude of Earth (degrees, J2000 ecliptic):
```
dlon = 100.46645683 + (1295977422.83429 * t - 2.04411 * t2 - 0.00523 * t3) / 3600.0
```

Light-time correction (aberration-of-light shift):
```
dlon = swe_degnorm(dlon - dlt * 360.0 / 365.2425)
```

Convert to unit-sphere ecliptic Cartesian (xs[0]=lon_rad, xs[1]=lat_rad=0, xs[2]=r=1):
```
xs[0] = dlon * DEGTORAD;  xs[1] = 0;  xs[2] = 1
swi_polcart(xs, xs)
```

Rotate ecliptic → equatorial J2000:
```
eps_J2000 = swi_epsiln(J2000 + swe_deltat_ex(J2000, -1, NULL), 0) * RADTODEG
swi_coortrf(xs, xs, -eps_J2000 * DEGTORAD)   // negative = ecliptic→equatorial
```

Precess J2000 equatorial → mean equatorial of date:
```
swi_precess(xs, tjd_et, 0, -1)
```

Rotate equatorial of date → ecliptic of date:
```
eps_mean = swi_epsiln(tjd_et, 0) * RADTODEG
swi_nutation(tjd_et, 0, nutlo)               // nutlo[0]=dpsi rad, nutlo[1]=deps rad
eps_true = eps_mean + nutlo[1] * RADTODEG
dpsi_deg = nutlo[0] * RADTODEG
swi_coortrf(xs, xs, eps_mean * DEGTORAD)     // positive = equatorial→ecliptic
swi_cartpol(xs, xs)
xs[0] *= RADTODEG                            // longitude now in degrees
```

Add UT fraction of day (hour angle component):
```
dhour = fmod(tjd_ut - 0.5, 1.0) * 360.0     // degrees [0, 360)
```

Add Equation of the Equinoxes (EoE):
```
if eps == 0 (caller did not supply):
    xs[0] += dpsi_deg * cos(eps_true * DEGTORAD)   // using freshly computed nutation
else:
    xs[0] += nut * cos(eps * DEGTORAD)              // using caller-supplied eps and dpsi
```

Normalize and convert to hours:
```
xs[0] = swe_degnorm(xs[0] + dhour)    // [0, 360)
tsid  = xs[0] / 15.0                  // degrees → hours
```

Boundary offset correction (applied in the dispatcher after calling this function):
- Before 1850: subtract `SIDT_LTERM_OFS0 = 0.000378172 / 15.0` hours
- After  2050: subtract `SIDT_LTERM_OFS1 = 0.001385646 / 15.0` hours

---

## Non-Polynomial Complementary Terms (swephlib.c:3413–3450)

Called only in the ERA-based path. Returns a correction in **degrees**.

Input: `tt` = TT centuries since J2000.

### 14 Fundamental Arguments (delm[0..13], radians)

| Index | Name | Formula |
|---|---|---|
| 0 | l (Moon anomaly) | `swe_radnorm(2.35555598 + 8328.6914269554 * tt)` |
| 1 | l' (Sun anomaly) | `swe_radnorm(6.24006013 + 628.301955 * tt)` |
| 2 | F (Moon lat arg) | `swe_radnorm(1.627905234 + 8433.466158131 * tt)` |
| 3 | D (elongation) | `swe_radnorm(5.198466741 + 7771.3771468121 * tt)` |
| 4 | Om (node) | `swe_radnorm(2.18243920 - 33.757045 * tt)` |
| 5 | L_Me (Mercury) | `swe_radnorm(4.402608842 + 2608.7903141574 * tt)` |
| 6 | L_Ve (Venus) | `swe_radnorm(3.176146697 + 1021.3285546211 * tt)` |
| 7 | L_Ea (Earth) | `swe_radnorm(1.753470314 + 628.3075849991 * tt)` |
| 8 | L_Ma (Mars) | `swe_radnorm(6.203480913 + 334.0612426700 * tt)` |
| 9 | L_Ju (Jupiter) | `swe_radnorm(0.599546497 + 52.9690962641 * tt)` |
| 10 | L_Sa (Saturn) | `swe_radnorm(0.874016757 + 21.3299104960 * tt)` |
| 11 | L_Ur (Uranus) | `swe_radnorm(5.481293871 + 7.4781598567 * tt)` |
| 12 | L_Ne (Neptune) | `swe_radnorm(5.321159000 + 3.8127774000 * tt)` |
| 13 | p_A (precession) | `(0.02438175 + 0.00000538691 * tt) * tt` (quadratic, not normalized) |

Note: `swe_radnorm` normalizes to [0, 2π). The precession argument (index 13) is deliberately not normalized.

### Linear Term Before the Loop

```
dadd = -0.87 * sin(delm[4]) * tt    // dadd in microarcseconds × time, not yet scaled
```

This is applied before the 33-term loop.

### 33-Term Fourier Loop (swephlib.c:3441–3447)

```
for i in 0..33:
    darg = sum(stfarg[i * 14 + j] * delm[j]  for j in 0..14)    // combined angle in radians
    dadd += stcf[i*2] * sin(darg) + stcf[i*2+1] * cos(darg)
```

Unit conversion at end:
```
dadd /= (3600.0 * 1000000.0)    // microarcseconds → degrees
return dadd
```

### Complete stcf Coefficient Table (swephlib.c:3341–3375)

`stcf[SIDTNTERM * 2]`: pairs of (sin_coeff, cos_coeff) in microarcseconds.

```
Row  0:  2640.96   -0.39
Row  1:    63.52   -0.02
Row  2:    11.75    0.01
Row  3:    11.21    0.01
Row  4:    -4.55    0.00
Row  5:     2.02    0.00
Row  6:     1.98    0.00
Row  7:    -1.72    0.00
Row  8:    -1.41   -0.01
Row  9:    -1.26   -0.01
Row 10:    -0.63    0.00
Row 11:    -0.63    0.00
Row 12:     0.46    0.00
Row 13:     0.45    0.00
Row 14:     0.36    0.00
Row 15:    -0.24   -0.12
Row 16:     0.32    0.00
Row 17:     0.28    0.00
Row 18:     0.27    0.00
Row 19:     0.26    0.00
Row 20:    -0.21    0.00
Row 21:     0.19    0.00
Row 22:     0.18    0.00
Row 23:    -0.10    0.05
Row 24:     0.15    0.00
Row 25:    -0.14    0.00
Row 26:     0.14    0.00
Row 27:    -0.14    0.00
Row 28:     0.14    0.00
Row 29:     0.13    0.00
Row 30:    -0.11    0.00
Row 31:     0.11    0.00
Row 32:     0.11    0.00
```

### Complete stfarg Multiplier Table (swephlib.c:3378–3412)

`stfarg[SIDTNTERM * 14]`: integer multipliers for [l, l', F, D, Om, L_Me, L_Ve, L_Ea, L_Ma, L_Ju, L_Sa, L_Ur, L_Ne, p_A].

```
Row  0:   0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  1:   0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  2:   0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  3:   0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  4:   0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  5:   0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  6:   0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  7:   0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  8:   0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row  9:   0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 10:   1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 11:   1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 12:   0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 13:   0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 14:   0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 15:   0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0
Row 16:   0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 17:   0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 18:   1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 19:   1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 20:   0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 21:   0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 22:   0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 23:   0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1
Row 24:   0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 25:   2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 26:   1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 27:   0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 28:   1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 29:   0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 30:   0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 31:   1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0
Row 32:   1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0
```

---

## Model Dispatch Summary

```
swe_sidtime(tjd_ut)
    → computes eps (mean+true obliquity in degrees), dpsi (degrees)
    → calls swe_sidtime0(tjd_ut, eps_true_degrees, dpsi_degrees)

swe_sidtime0(tjd, eps, nut):
    model = swed.astro_models[7]  (default = 4 = LONGTERM)

    LONGTERM outside 1850–2050:
        → sidtime_long_term(tjd, eps, nut)
        apply offset, clamp to [0, 24)
        return immediately (no EoE step)

    LONGTERM inside 1850–2050  OR  IERS_CONV_2010:
        ERA path:
            ERA = swe_degnorm((0.7790572732640 + 1.00273781191135448*(tjd-J2000)) * 360)
            polynomial correction in degrees (arcsec/3600)
            + sidtime_non_polynomial_part(tt_centuries)
            convert degrees → sidereal seconds
        add EoE; normalize; convert to hours

    IAU_2006:
        GMST polynomial using TT at midnight; msday * secs
        add EoE; normalize; convert to hours

    IAU_1976:
        GMST polynomial using UT at midnight; msday * secs
        add EoE; normalize; convert to hours
```

---

## Constants

| Name | Value | Location | Meaning |
|---|---|---|---|
| `J2000` | 2451545.0 | sweph.h:67 | Julian Day of J2000.0 (2000 Jan 1.5 TT) |
| `DEGTORAD` | π/180 | sweodef.h:266 | Degree to radian factor |
| `RADTODEG` | 180/π | sweodef.h:265 | Radian to degree factor |
| `AUNIT` | 1.49597870700e+11 m | sweph.h:273 | 1 AU in metres (DE431) |
| `CLIGHT` | 2.99792458e+8 m/s | sweph.h:274 | Speed of light |
| `dlt` | AUNIT/CLIGHT/86400 | computed | Light-time for 1 AU in days (~0.00578) |
| `SIDT_LTERM_T0` | 2396758.5 | swephlib.c:3460 | 1 Jan 1850 — long-term model boundary |
| `SIDT_LTERM_T1` | 2469807.5 | swephlib.c:3461 | 1 Jan 2050 — long-term model boundary |
| `SIDT_LTERM_OFS0` | 0.000378172/15 | swephlib.c:3462 | Pre-1850 hour offset |
| `SIDT_LTERM_OFS1` | 0.001385646/15 | swephlib.c:3463 | Post-2050 hour offset |
| `SIDTNTERM` | 33 | swephlib.c:3340 | Non-polynomial term count |
| `SIDTNARG` | 14 | swephlib.c:3376 | Fundamental arguments per term |

---

## Unit Flows

### IAU 1976 / IAU 2006 polynomial paths

```
tu         → UT1 centuries
tt         → TT centuries (IAU 2006 only)
gmst_0h    → sidereal SECONDS at UT midnight
secs       → UT seconds since midnight
msday      → dimensionless (≈ 1.00273...)
gmst       → sidereal SECONDS (at tjd_ut)
eqeq       → sidereal SECONDS (= 240 * dpsi_deg * cos(eps))
gmst       → sidereal SECONDS (after EoE)
gmst % 86400 → sidereal SECONDS in [0, 86400)
gmst / 3600  → hours in [0, 24)
```

### ERA path

```
jdrel      → UT days from J2000
tt         → TT centuries from J2000 (uses full tjd, not midnight)
ERA        → degrees [0, 360)
polynomial → degrees (added directly)
dadd       → degrees (from microarcseconds / (3600 * 1e6))
gmst       → degrees [0, 360)
gmst * 240 → sidereal SECONDS (= gmst / 15 * 3600)
eqeq       → sidereal SECONDS
gmst % 86400 → sidereal SECONDS in [0, 86400)
gmst / 3600  → hours in [0, 24)
```

### Long-term path

```
t           → TT Julian MILLENNIA from J2000
dlon        → degrees (mean Earth longitude)
xs[]        → Cartesian (unit sphere)
xs[0] after swi_cartpol → radians
xs[0] * RADTODEG → degrees
dhour       → degrees [0, 360)
xs[0] + dhour → degrees
swe_degnorm(...)  → degrees [0, 360)
/ 15.0      → hours [0, 24)
```

---

## Time Variable Distinctions

| Variable | Type | Epoch | Scale |
|---|---|---|---|
| `tjd` / `tjd_ut` | input | — | UT1 |
| `tjd_et` | computed | — | TT |
| `jd0` | midnight UT | — | UT1 |
| `secs` | time of day | since midnight | UT1 seconds |
| `tu` | centuries | J2000 | UT1 |
| `tt` | centuries | J2000 | TT |
| `t` | millennia | J2000 | TT (long-term only) |
| `jdrel` | days | J2000 | UT1 (ERA calculation) |

The TT conversion always uses `swe_deltat_ex(tjd, -1, NULL)`. Flag `-1` means "use the default tidal acceleration for the current configured ephemeris."

---

## Notes on Calling Conventions

**When `eps == 0` and `nut == 0` are passed to `sidtime_long_term`**: The function detects `eps == 0` and uses freshly computed nutation from `swi_nutation()` for the EoE. If `eps != 0`, the caller-supplied values are used. In the normal flow, `swe_sidtime()` always passes `eps != 0` (it computes eps first), so the long-term function uses caller-supplied values in production.

**`swe_sidtime0` vs. `swe_sidtime`**: The `_sidtime0` form exists so callers that have already computed obliquity and nutation can avoid recomputing them (e.g., house calculations). The `_sidtime` form is the convenience wrapper that computes everything internally.

---

## References

| Source | Used in |
|---|---|
| IAU 1976 (Aoki et al. 1982, A&A 105, 359) | `SEMOD_SIDT_IAU_1976` polynomial |
| Capitaine, Wallace & Chapront 2003, A&A 412, 567–586 | `SEMOD_SIDT_IAU_2006` polynomial and ERA path (p. 582) |
| IERS Conventions 2010, Chapter 5, Table 5.2e | ERA formula and 33 complementary terms |
| Simon et al. 1994, A&A 282, 663–683 | Long-term mean Earth longitude |
| Souchay et al. 1999 | Planetary longitude coefficients in complementary terms |
