# C Reference: Moshier Ephemeris — swemplan.c / swemmoon.c / swemptab.h

Porting reference for the Moshier analytical ephemeris backend. Read this instead of the C source.
Covers planetary series evaluation (`swemplan.c`, `swemptab.h`) and lunar series evaluation (`swemmoon.c`).

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swi_moshplan2` | swemplan.c:134–264 | Yes — core series evaluator |
| `swi_moshplan` | swemplan.c:276–381 | Yes — wrapper, coordinate transform, speed |
| `sscc` (planets) | swemplan.c:387–408 | Yes — sin/cos lookup builder |
| `embofs_mosh` | swemplan.c:416–491 | Yes — EMB → Earth offset |
| `swi_osc_el_plan` | swemplan.c:579–689 | No — fictional planets via osculating elements, separate task |
| `swi_moshmoon2` | swemmoon.c:848–862 | Yes — Moon geometric polar coords |
| `swi_moshmoon` | swemmoon.c:869–934 | Yes — wrapper, to equatorial J2000, speed |
| `mean_elements` | swemmoon.c:1763–1818 | Yes — lunar fundamental arguments |
| `mean_elements_pl` | swemmoon.c:1820–1850 | Yes — planetary longitudes for Moon |
| `moon1` | swemmoon.c:1182–1364 | Yes — T-dependent corrections |
| `moon2` | swemmoon.c:1367–1442 | Yes — additional T^0 planetary terms |
| `moon3` | swemmoon.c:1444–1454 | Yes — main table evaluation, assemble |
| `moon4` | swemmoon.c:1458–1464 | Yes — unit conversion |
| `chewm` | swemmoon.c:1628–1691 | Yes — harmonic table evaluator |
| `sscc` (moon) | swemmoon.c:1696–1714 | Yes — same as planet version, different dims |
| `ecldat_equ2000` | swemmoon.c:1722–1729 | Yes — ecliptic of date → equatorial J2000 |
| `mods3600` (moon) | swemmoon.c:1734–1740 | Yes — arcsec modulo 1296000 |
| `swi_mean_node` | swemmoon.c:1493–1534 | Yes — mean lunar node |
| `swi_mean_apog` | swemmoon.c:1564–1624 | Yes — mean lunar apogee (Lilith) |

**Note**: The C source is compiled without `MOSH_MOON_200`. The DE404 code paths (the `#else` branches) are the ones to port.

---

## Part 1: Planetary Ephemeris (swemplan.c, swemptab.h)

### Constants and Macros

| Name | Value | Location | Notes |
|---|---|---|---|
| `STR` | `4.8481368110953599359e-6` | sweph.h:663 | Radians per arcsecond |
| `TIMESCALE` | `3652500.0` | swemplan.c:67 | Days per 10,000 Julian years |
| `J2000` | `2451545.0` | sweph.h:67 | Julian date of J2000.0 epoch |
| `PLAN_SPEED_INTV` | `0.0001` | sweph.h:299 | Speed interval in days (8.64 s) |
| `MOSHPLEPH_START` | `625000.5` | sweph.h:219 | Ephemeris start JD |
| `MOSHPLEPH_END` | `2818000.5` | sweph.h:220 | Ephemeris end JD |
| `EARTH_MOON_MRAT` | `81.30056` | sweph.h:270 | Earth/Moon mass ratio (DE406) |
| `SEFLG_MOSEPH` | `4` | swephexp.h:188 | Moshier backend flag |
| `mods3600(x)` | `x - 1.296e6 * floor(x/1.296e6)` | swemplan.c:69 | Reduce arcsec to [0, 1296000) |
| `KGAUSS_GEO` | `0.0000298122353216` | swemplan.c:73 | Gaussian for geocentric bodies |

### Planet Index Mapping

The `pnoint2msh[]` array (swemplan.c:84) maps SEI_* constants to the `planets[]` array index:

```
pnoint2msh = {2, 2, 0, 1, 3, 4, 5, 6, 7, 8}
```

| SEI constant | Value | pnoint2msh | Table |
|---|---|---|---|
| SEI_EMB / SEI_EARTH | 0 | 2 | ear404 |
| SEI_MOON | 1 | 2 | (not used via moshplan) |
| SEI_MERCURY | 2 | 0 | mer404 |
| SEI_VENUS | 3 | 1 | ven404 |
| SEI_MARS | 4 | 3 | mar404 |
| SEI_JUPITER | 5 | 4 | jup404 |
| SEI_SATURN | 6 | 5 | sat404 |
| SEI_URANUS | 7 | 6 | ura404 |
| SEI_NEPTUNE | 8 | 7 | nep404 |
| SEI_PLUTO | 9 | 8 | plu404 |

The `planets[]` C array (swemplan.c:116–127) holds: `[mer404, ven404, ear404, mar404, jup404, sat404, ura404, nep404, plu404]`.

### Fundamental Arguments

Nine planetary mean longitudes from Simon et al. (1994). T is in TIMESCALE units (10,000 Julian years from J2000):

```
arg_i = mods3600(freqs[i] * T) + phases[i]    (arcseconds)
```

| i | Body | freqs[i] (arcsec / 10000 yr) | phases[i] (arcsec) |
|---|---|---|---|
| 0 | Mercury | 53810162868.8982 | 252.25090552 × 3600 |
| 1 | Venus | 21066413643.3548 | 181.97980085 × 3600 |
| 2 | Earth | 12959774228.3429 | 100.46645683 × 3600 |
| 3 | Mars | 6890507749.3988 | 355.43299958 × 3600 |
| 4 | Jupiter | 1092566037.7991 | 34.35151874 × 3600 |
| 5 | Saturn | 439960985.5372 | 50.07744430 × 3600 |
| 6 | Uranus | 154248119.3933 | 314.05500511 × 3600 |
| 7 | Neptune | 78655032.0744 | 304.34866548 × 3600 |
| 8 | Pluto | 52272245.1795 | 860492.1546 (directly) |

These are used in `swi_moshplan2` step 1 to build sin/cos lookup tables.

### struct plantbl (sweph.h:698–706)

```c
struct plantbl {
  char max_harmonic[9];    // max harmonic k for each of 9 planetary args
  char max_power_of_t;     // highest polynomial degree in T
  signed char *arg_tbl;    // encoded argument table (see below)
  double *lon_tbl;         // longitude coefficients
  double *lat_tbl;         // latitude coefficients
  double *rad_tbl;         // radius (deviation) coefficients
  double distance;         // mean heliocentric distance in AU
};
```

The `max_harmonic[i]` value for index i means the sscc table for argument i must be built up to harmonic `max_harmonic[i]`. If 0, that argument is not used by this planet.

### Planet Table Instances (swemptab.h)

| Instance | Line | max_harmonic[9] | max_power_of_t | distance (AU) | Total terms |
|---|---|---|---|---|---|
| mer404 | 1073 | {11,14,10,11,4,5,2,0,0} | 6 | 3.8709830979e-1 | 130 |
| ven404 | 1898 | {5,14,13,8,4,5,1,0,0} | 5 | 7.2332982000e-1 | 108 |
| ear404 | 2941 | {1,9,14,17,5,5,2,1,0} | 4 | 1.0 | 135 |
| mar404 | 4500 | {0,5,12,24,9,7,3,2,0} | 5 | 1.5303348827e+0 | 201 |
| jup404 | 5714 | {0,0,1,0,9,16,7,5,0} | 6 | 5.2026032092e+0 | 142 |
| sat404 | 7430 | {0,0,1,0,8,18,9,5,0} | 7 | 9.5575813548e+0 | 215 |
| ura404 | 8793 | {0,0,0,0,5,10,9,12,0} | 6 | 1.9218446061e+1 | 177 |
| nep404 | 9252 | {0,0,0,0,3,8,7,9,0} | 3 | 3.0110386869e+1 | 59 |
| plu404 | 10632 | {0,0,0,0,2,2,9,13,13} | 7 | 3.9540000000e+1 | 173 |

Each instance has three corresponding double arrays: `{body}tabl` (longitude), `{body}tabb` (latitude), `{body}tabr` (radius). For Mercury as a concrete example:

- `mertabl[]` — longitude coefficients, starts swemptab.h:77
- `mertabb[]` — latitude coefficients, starts swemptab.h:364
- `mertabr[]` — radius coefficients, starts swemptab.h:651
- `merargs[]` — arg_tbl, starts swemptab.h:939 (134 signed-char entries + sentinel `-1`)

### arg_tbl Encoding (Critical — Non-Obvious)

The `arg_tbl` is a stream of `signed char` values with three record types. `p` is the read cursor:

**End sentinel:**
```
np = *p++  →  np < 0  →  stop
```

**Polynomial term (np == 0):**
```
np = *p++    // == 0
nt = *p++    // polynomial degree
// Consume (nt+1) doubles from lon_tbl:
//   cu = *pl++;
//   for (ip = 0; ip < nt; ip++) cu = cu * T + *pl++;
//   sl += mods3600(cu);      ← longitude gets mods3600 wrap
// Same pattern for lat_tbl → sb += cu  (no mods3600)
// Same pattern for rad_tbl → sr += cu  (no mods3600)
```

Total doubles consumed per polynomial term: `(nt+1)` each from lon, lat, rad.

**Periodic term (np > 0):**
```
np = *p++    // number of argument pairs (1..n)
for (ip = 0; ip < np; ip++) {
    j = *p++    // harmonic multiplier (signed, nonzero = contributes)
    m = *p++    // planet index, 1-based (subtract 1 to index ss[][])
}
nt = *p++    // highest power of T for amplitude
// Consume 2*(nt+1) doubles from lon_tbl (cos amplitude, sin amplitude, alternating):
//   cu = *pl++;  su = *pl++;
//   for (ip = 0; ip < nt; ip++) { cu = cu * T + *pl++;  su = su * T + *pl++; }
//   sl += cu * cv + su * sv;     ← cv=cos(arg), sv=sin(arg)
// Same pattern for lat_tbl → sb
// Same pattern for rad_tbl → sr
```

Total doubles consumed per periodic term: `2*(nt+1)` each from lon, lat, rad.

**Argument combination** (from np pairs):

```
sv = 0, cv = 0, k1 = 0
for each (j, m) pair:
    if j == 0: skip
    k = |j| - 1
    su = ss[m-1][k]    // sin(|j| * angle_m)
    if j < 0: su = -su
    cu = cc[m-1][k]    // cos(|j| * angle_m)
    if k1 == 0:
        sv = su; cv = cu; k1 = 1
    else:
        t = su*cv + cu*sv
        cv = cu*cv - su*sv
        sv = t
// result: sv = sin(combined), cv = cos(combined)
```

**Example decode** from `merargs[]` (first two records):

```
0, 3,                          // polynomial, degree 3 → lon/lat/rad each consume 4 doubles
3, 1, 1,-10, 3, 11, 4, 0,     // np=3: pairs=(1,1),(-10,3),(11,4); nt=0 → 2 doubles each
//   arg = 1×Mercury - 10×Earth + 11×Mars
```

### sscc() — sin/cos Lookup Builder (swemplan.c:387–408)

Builds `ss[k][0..n-1]` and `cc[k][0..n-1]` where `ss[k][i] = sin((i+1)*arg)`:

```
ss[k][0] = sin(arg);   cc[k][0] = cos(arg)
ss[k][1] = 2*su*cu;    cc[k][1] = cu*cu - su*su
for i = 2 .. n-1:
    s = su*cv + cu*sv
    cv = cu*cv - su*sv
    sv = s
    ss[k][i] = sv;  cc[k][i] = cv
```

Global TLS arrays (swemplan.c:129–130): `ss[9][24]` and `cc[9][24]`. Populated once per `swi_moshplan2` call. Maximum harmonic across all planets is 18 (Saturn index 5).

### swi_moshplan2 (swemplan.c:134–264)

**Signature:** `int swi_moshplan2(double J, int iplm, double *pobj)`

- `J` — Julian day (TDT)
- `iplm` — Moshier planet index (0..8, use `pnoint2msh[ipli]`)
- `pobj[3]` — output: ecliptic heliocentric polar, J2000 ecliptic (rad, rad, AU)

**Algorithm:**

```
T = (J - J2000) / TIMESCALE    // J2000 = 2451545.0, TIMESCALE = 3652500.0

// Step 1: build sin/cos tables for each argument
plan = planets[iplm]
for i = 0..8:
    if plan->max_harmonic[i] > 0:
        sr = (mods3600(freqs[i] * T) + phases[i]) * STR    // arcsec → radians
        sscc(i, sr, plan->max_harmonic[i])

// Step 2: iterate arg_tbl
p = plan->arg_tbl; pl = plan->lon_tbl; pb = plan->lat_tbl; pr = plan->rad_tbl
sl = sb = sr = 0.0
loop:
    np = *p++
    if np < 0: break
    if np == 0:
        nt = *p++
        // evaluate Horner poly of degree nt from pl, pb, pr
        sl += mods3600(lon_poly(T))
        sb += lat_poly(T)
        sr += rad_poly(T)
        continue
    // combine np argument pairs into (sv, cv)
    nt = *p++
    // evaluate Horner poly of degree nt in two streams from pl, pb, pr
    sl += cu*cv + su*sv    // (cu,su) = Horner cos/sin amplitude polynomial
    sb += cu*cv + su*sv
    sr += cu*cv + su*sv

// Step 3: convert outputs
pobj[0] = STR * sl                               // longitude in radians
pobj[1] = STR * sb                               // latitude in radians
pobj[2] = STR * plan->distance * sr + plan->distance   // radius in AU
```

The radius formula: `distance * (1 + STR * sr)` — sr is the normalised deviation from mean distance, in arcseconds, so `STR * sr` is the fractional deviation in radians (pure number).

**Units of sl, sb, sr:** all accumulated in arcseconds. `mods3600` wraps longitude within [0, 1296000).

### swi_moshplan (swemplan.c:276–381)

**Signature:** `int swi_moshplan(double J, int ipli, AS_BOOL do_save, double *xpret, double *xeret, char *serr)`

**Output:** heliocentric equatorial Cartesian J2000, positions in AU, speeds in AU/day, both via 6-element arrays `[x, y, z, vx, vy, vz]`.

**Range check:** `MOSHPLEPH_START - 0.3 ≤ tjd ≤ MOSHPLEPH_END + 0.3` (returns ERR if outside).

**Earth/EMB pipeline:**
```
swi_moshplan2(tjd, ear_msh_idx=2, xe)    // ecliptic heliocentric polar (rad, rad, AU) of EMB
swi_polcart(xe, xe)                       // polar → cartesian
swi_coortrf2(xe, xe, -seps2000, ceps2000) // ecliptic J2000 → equatorial J2000
embofs_mosh(tjd, xe)                      // EMB → geocenter

// Speed: backward difference
swi_moshplan2(tjd - PLAN_SPEED_INTV, 2, x2)
swi_polcart(x2, x2)
swi_coortrf2(x2, x2, -seps2000, ceps2000)
embofs_mosh(tjd - PLAN_SPEED_INTV, x2)
xe[3..5] = (xe[0..2] - x2[0..2]) / PLAN_SPEED_INTV    // AU/day
```

Uses `swed.oec2000.seps` and `.ceps` (sin/cos of J2000 mean obliquity).

**Planet pipeline:** same structure but without `embofs_mosh`. Speed via backward difference at `tjd - PLAN_SPEED_INTV`.

**Caching:** if `tjd == pdp->teval && pdp->iephe == SEFLG_MOSEPH`, returns cached result.

### embofs_mosh (swemplan.c:416–491)

Adjusts EMB position to geocenter by computing a short Moon series and subtracting 1/(EARTH_MOON_MRAT + 1) of the Moon vector.

Epoch: T = (tjd - J1900) / 36525.0 (centuries from J1900 = 2415020.0).

**Moon arguments (all in degrees, then converted to radians):**
- MP (mean anomaly): `swe_degnorm(((1.44e-5*T + 0.009192)*T + 477198.8491)*T + 296.104608)`
- 2D (2× elongation): `2 * swe_degnorm(((1.9e-6*T - 0.001436)*T + 445267.1142)*T + 350.737486)`
- F (mean latitude arg): `swe_degnorm(((-3e-7*T - 0.003211)*T + 483202.0251)*T + 11.250889)`
- 2D−MP: derived from sin/cos product identities

**Moon ecliptic longitude L** (degrees, 6-term approximation):
```
L = ((1.9e-6*T - 0.001133)*T + 481267.8831)*T + 270.434164
M = swe_degnorm((...)*T + 358.475833)   // solar anomaly
L += 6.288750*sin(MP) + 1.274018*sin(2D-MP) + 0.658309*sin(2D)
   + 0.213616*sin(2MP) - 0.185596*sin(M) - 0.114336*sin(2F)
```

**Moon ecliptic latitude B** (degrees, 4 terms), **parallax p** (degrees, 5 terms), then:
```
a = 4.263523e-5 / sin(p)    // distance in AU
xyz = polcart([L_rad, B_rad, a])
xyz = ecliptic_of_date → equatorial_of_date → precess_to_J2000
xe[0..2] -= xyz[0..2] / (EARTH_MOON_MRAT + 1.0)
```

Uses `swed.oec.seps/.ceps` (obliquity of date, not J2000).

---

## Part 2: Coefficient Data Tables (swemptab.h)

### Mercury as Template (Full Detail)

**mer404** (swemptab.h:1073–1081):
```c
static struct plantbl mer404 = {
  { 11, 14, 10, 11,  4,  5,  2,  0,  0,},  // max_harmonic per argument
  6,                                          // max_power_of_t
  merargs,                                    // arg_tbl: signed char[]
  mertabl,                                    // lon_tbl: double[]
  mertabb,                                    // lat_tbl: double[]
  mertabr,                                    // rad_tbl: double[]
  3.8709830979999998e-01,                     // mean distance in AU
};
```

**merargs** (swemptab.h:939–1071): 134 entries total (131 bytes + sentinel `-1`).

Each entry reads as a stream of bytes. Structure is as documented in the arg_tbl encoding section above. The comment at line 1072 says "Total terms = 130, small = 128" meaning 130 rows total (2 are polynomial terms).

**mertabl / mertabb / mertabr** (swemptab.h:77, 364, 651): three separate `double[]` arrays storing longitude, latitude, and radius coefficients respectively. All three arrays are consumed in parallel as `swi_moshplan2` iterates through `merargs`. The layout mirrors the arg_tbl exactly: polynomial terms consume `(nt+1)` doubles, periodic terms consume `2*(nt+1)` doubles.

**First two records decoded:**

Record 1: `0, 3` → polynomial, degree 3.
- mertabl[0..3]: `{35.85255, -163.26379, 53810162857.56026, 908082.18475}` arcsec
- Result: `sl += mods3600(35.85255*T³ - 163.26379*T² + 53810162857.56026*T + 908082.18475)`
- The dominant term `53.8e9 * T` is Mercury's mean longitude rate (arcsec / TIMESCALE).

Record 2: `3, 1, 1,-10, 3, 11, 4, 0` → np=3, pairs: (j=1,m=1), (j=-10,m=3), (j=11,m=4); nt=0.
- Angle: `+1×Mercury - 10×Earth + 11×Mars`
- mertabl[4..5]: `{0.05214, -0.07712}` → `sl += 0.05214*cos(arg) + (-0.07712)*sin(arg)`

### All Planet Tables

| Planet | lon_tbl start | lat_tbl start | rad_tbl start | arg_tbl start |
|---|---|---|---|---|
| Mercury | :77 | :364 | :651 | :939 |
| Venus | :1101 | :1329 | :1557 | venargs after ven404 struct |
| Earth/EMB | :1925 | :2217 | :2509 | earargs |
| Mars | :2968 | :3410 | :3852 | marargs |
| Jupiter | :4529 | :4875 | :5221 | jupargs |
| Saturn | :5743 | :6232 | :6721 | satargs |
| Uranus | :7459 | :7843 | :8227 | uraargs |
| Neptune | :8822 | :8944 | :9066 | nepargs |
| Pluto | :9281 | :9672 | :10063 | pluargs |

All follow the exact same encoding as Mercury. The size of each table is determined by iterating the arg_tbl and counting coefficient consumption, not from any header.

---

## Part 3: Lunar Ephemeris (swemmoon.c)

Based on ELP2000-85 (Chapront-Touzé & Chapront 1988), adjusted by least-squares fit to DE404 over −3000 to +3000. The fit uses 34,247 Lunar positions at 64-day intervals. Maximum discrepancy vs DE404: ~7″ longitude, ~5″ latitude, ~5×10⁻⁸ AU radius.

**Compile-time switch**: `MOSH_MOON_200` is **not defined** in Swiss Ephemeris — use the `#else` branches throughout.

### Constants

| Name | Value | Notes |
|---|---|---|
| `MOSHLUEPH_START` | `625000.5` | Moon ephemeris start JD (with large-range build: `-225000.5`) |
| `MOSHLUEPH_END` | `2818000.5` | Moon ephemeris end JD (with large-range: `3600000.5`) |
| `MOSHNDEPH_START` | `-3100015.5` | Mean node/apog start JD |
| `MOSHNDEPH_END` | `8000016.5` | Mean node/apog end JD |
| `MOON_SPEED_INTV` | `0.00005` | Days (4.32 seconds) for speed central difference |
| `MOON_MEAN_DIST` | `384400000.0` | m, used in mean apogee distance |
| `MOON_MEAN_ECC` | `0.054900489` | Used in mean apogee distance |
| `MOON_MEAN_INCL` | `5.1453964°` | Used in apogee inclination transform |
| `AUNIT` | `1.49597870700e+11` | m per AU (DE431) |
| `STR` | `4.8481368110953599359e-6` | rad/arcsec (shared with planets) |
| `mods3600(x)` | `x - 1296000.0 * floor(x/1296000.0)` | Local function, same as planet macro |

### TLS State Variables (swemmoon.c:811–843)

All declared `static TLS` (thread-local). These are shared across `moon1()` through `moon4()`. In Rust's stateless design these become local variables passed between or captured in a computation context:

| Variable | Type | Role |
|---|---|---|
| `T` | double | Julian centuries from J2000 |
| `T2` | double | T² |
| `T3` | double | T³ (set for node/apogee, zero in main moon) |
| `T4` | double | T⁴ |
| `SWELP` | double | Mean lunar longitude (arcsec) |
| `M` | double | Mean solar anomaly / l′ (arcsec) |
| `MP` | double | Mean lunar anomaly / l (arcsec) |
| `D` | double | Mean elongation (arcsec) |
| `NF` | double | F = mean argument of latitude (arcsec) |
| `Ve` | double | Venus mean longitude (arcsec) |
| `Ea` | double | Earth mean longitude (arcsec) |
| `Ma` | double | Mars mean longitude (arcsec) |
| `Ju` | double | Jupiter mean longitude (arcsec) |
| `Sa` | double | Saturn mean longitude (arcsec) |
| `f` | double | Temp: `18*Ve - 16*Ea` (arcsec) |
| `g` | double | Temp: current argument (radians = STR * arcsec_angle) |
| `cg`, `sg` | double | cos(g), sin(g) |
| `l` | double | Longitude accumulator T^0 term (arcsec) |
| `l1` | double | Longitude T^1 coefficient (units: 10⁻⁵ arcsec, i.e., 0.00001″) |
| `l2` | double | Longitude T^2 coefficient (same units) |
| `l3` | double | Longitude T^3 coefficient (same units) |
| `l4` | double | Longitude T^4 coefficient (same units) |
| `B` | double | Latitude accumulator T^0 term (arcsec) |
| `moonpol[3]` | double[3] | Accumulates [lon, lat, rad] perturbations |
| `ss[5][8]` | double | sin lookup: ss[k][i] = sin((i+1)×arg_k) |
| `cc[5][8]` | double | cos lookup: cc[k][i] = cos((i+1)×arg_k) |

The 5 ss/cc rows correspond to: D (k=0), M (k=1), MP (k=2), NF (k=3), unused (k=4).

### z[] — DE404 Fitting Coefficients (swemmoon.c:284–313)

25 values (DE404 version). Replaces higher-order and planetary secular terms in the mean elements and moon1 perturbations. Unit notes per group:

```c
static const double z[] = {
/* Scaled in arc seconds, time in Julian centuries */
-1.312045233711e+01, /* z[0]:  F,   T^2 coefficient (arcsec/century^2) */
-1.138215912580e-03, /* z[1]:  F,   T^3 */
-9.646018347184e-06, /* z[2]:  F,   T^4 */
 3.146734198839e+01, /* z[3]:  l,   T^2 */
 4.768357585780e-02, /* z[4]:  l,   T^3 */
-3.421689790404e-04, /* z[5]:  l,   T^4 */
-6.847070905410e+00, /* z[6]:  D,   T^2 */
-5.834100476561e-03, /* z[7]:  D,   T^3 */
-2.905334122698e-04, /* z[8]:  D,   T^4 */
-5.663161722088e+00, /* z[9]:  L,   T^2 */
 5.722859298199e-03, /* z[10]: L,   T^3 */
-8.466472828815e-05, /* z[11]: L,   T^4 */
/* Longitude terms in arc seconds × 10^5 (divided by 1e5 in moon3 via *T*1e-5) */
-8.429817796435e+01, /* z[12]: T^2 cos(18V - 16E - l) */
-2.072552484689e+02, /* z[13]: T^2 sin(18V - 16E - l) */
 7.876842214863e+00, /* z[14]: T^2 cos(10V - 3E - l) */
 1.836463749022e+00, /* z[15]: T^2 sin(10V - 3E - l) */
-1.557471855361e+01, /* z[16]: T^2 cos(8V - 13E) */
-2.006969124724e+01, /* z[17]: T^2 sin(8V - 13E) */
 2.152670284757e+01, /* z[18]: T^2 cos(4E - 8M + 3J) */
-6.179946916139e+00, /* z[19]: T^2 sin(4E - 8M + 3J) */
-9.070028191196e-01, /* z[20]: T^2 cos(18V - 16E) */
-1.270848233038e+01, /* z[21]: T^2 sin(18V - 16E) */
-2.145589319058e+00, /* z[22]: T^2 cos(2J - 5S) */
 1.381936399935e+01, /* z[23]: T^2 sin(2J - 5S) */
/* T^3 longitude term */
-1.999840061168e+00, /* z[24]: T^3 sin(l') */
};
```

**z[0..11]** are applied in `mean_elements()` as secular corrections to F, l, D, L respectively. Units: arcseconds/century^n.

**z[12..24]** are applied in `moon1()` as planetary perturbation T-polynomial coefficients. Units: arcseconds × 10⁵ (the `*T*1e-5` in moon3 converts them back to arcseconds).

### mean_elements() (swemmoon.c:1763–1818, DE404 version)

Computes fundamental lunar arguments in arcseconds, stored into TLS variables.

**T** = (J − 2451545.0) / 36525.0 (Julian centuries from J2000).
**fracT** = fmod(T, 1) — used to split the large coefficient to reduce floating-point error.

**M** (mean solar anomaly = l′) — 9th-degree polynomial from Laskar:
```
M = mods3600(129600000.0 * fracT - 3418.961646 * T + 1287104.76154)
  + ((((((((1.62e-20*T - 1.0390e-17)*T - 3.83508e-15)*T
         + 4.237343e-13)*T + 8.8555011e-11)*T
        - 4.77258489e-8)*T - 1.1297037031e-5)*T
       + 1.4732069041e-4)*T - 0.552891801772) * T2
```

**NF** (F = mean argument of latitude, Moon distance from node):
```
NF = mods3600(1739232000.0 * fracT + 295263.0983 * T
              - 2.079419901760e-01 * T + 335779.55755)
NF += ((z[2]*T + z[1])*T + z[0]) * T2
```

**MP** (l = mean lunar anomaly):
```
MP = mods3600(1717200000.0 * fracT + 715923.4728 * T
              - 2.035946368532e-01 * T + 485868.28096)
MP += ((z[5]*T + z[4])*T + z[3]) * T2
```

**D** (mean elongation):
```
D = mods3600(1601856000.0 * fracT + 1105601.4603 * T
             + 3.962893294503e-01 * T + 1072260.73512)
D += ((z[8]*T + z[7])*T + z[6]) * T2
```

**SWELP** (mean lunar longitude = L):
```
SWELP = mods3600(1731456000.0 * fracT + 1108372.83264 * T
                 - 6.784914260953e-01 * T + 785939.95571)
SWELP += ((z[11]*T + z[10])*T + z[9]) * T2
```

All results in arcseconds, in [0, 1296000).

### mean_elements_pl() (swemmoon.c:1820–1850)

Planetary mean longitudes for perturbation terms. High-precision polynomials from Laskar/Bretagnon, in arcseconds:

**Ve** (Venus): `mods3600(210664136.4335482 * T + 655127.283046) + secular_poly_9th_order * T2`

**Ea** (Earth): `mods3600(129597742.26669231 * T + 361679.214649) + secular_poly_9th_order * T2`

**Ma** (Mars): `mods3600(68905077.59284 * T + 1279559.78866) + (-1.043e-5*T + 9.38012e-3) * T2`

**Ju** (Jupiter): `mods3600(10925660.428608 * T + 123665.342120) + (1.543273e-5*T - 3.06037836351e-1) * T2`

**Sa** (Saturn): `mods3600(4399609.65932 * T + 180278.89694) + ((4.475946e-8*T - 6.874806e-5)*T + 7.56161437443e-1) * T2`

The 9th-order secular polynomials for Ve and Ea are Horner-form coefficients in T (see lines 1824–1843 for full coefficients).

### Harmonic Tables

All tables are `const short` arrays. The `chewm()` function reads them.

**Table summary** (DE404 version):

| Name | Size (elements) | Row width | typflg | T power | Content |
|---|---|---|---|---|---|
| `LR[8*NLR]` | 8×118=944 | 8 shorts | 1 | T^0 | Main lon+rad: D,l′,l,F + lon(1″,.0001″) + rad(1km,.0001km) |
| `MB[6*NMB]` | 6×77=462 | 6 shorts | 3 | T^0 | Main lat: D,l′,l,F + lat(1″,.0001″) |
| `LRT[8*NLRT]` | 8×38=304 | 8 shorts | 1 | ×T | Lon+rad: D,l′,l,F + lon(.1″,.00001″) + rad(.1km,.00001km) |
| `BT[5*NBT]` | 5×16=80 | 5 shorts | 4 | ×T | Lat: D,l′,l,F + lat(.00001″) |
| `LRT2[6*NLRT2]` | 6×25=150 | 6 shorts | 2 | ×T² | Lon+rad: D,l′,l,F + lon(.00001″) + rad(.00001km) |
| `BT2[5*NBT2]` | 5×12=60 | 5 shorts | 4 | ×T² | Lat: D,l′,l,F + lat(.00001″) |

`NLR=118`, `NMB=77`, `NLRT=38`, `NBT=16`, `NLRT2=25`, `NBT2=12`.

**LR first row** (`0, 0, 1, 0, 22639, 5858,-20905,-3550`):
- Angle: 0×D + 0×l′ + 1×l + 0×F = l (Moon mean anomaly)
- lon = (10000×22639 + 5858) × sin(l) = 226,395,858 × sin(l)  [×10⁻⁴ arcsec = 22639.5858″]
- rad = (10000×(-20905) + (-3550)) × cos(l)  [×10⁻⁴ km = -20905.355 km]

**MB first row** (`0, 0, 0, 1, 18461, 2387`):
- Angle: F (mean argument of latitude)
- lat = (10000×18461 + 2387) × sin(F) = 184,612,387 × sin(F)  [×10⁻⁴ arcsec = 18461.2387″]

### chewm() — Table Evaluator (swemmoon.c:1628–1691)

**Signature:** `void chewm(const short *pt, int nlines, int nangles, int typflg, double *ans)`

Parameters:
- `pt` — pointer into a harmonic table (advanced through the table)
- `nlines` — number of rows to process
- `nangles` — angular multiplier columns per row (always 4 for Moon: D, M, MP, NF)
- `typflg` — selects coefficient format (see below)
- `ans[3]` — output accumulators: ans[0]=longitude, ans[1]=latitude, ans[2]=radius

**Angle combination** (same as swi_moshplan2, but indexed from ss[m]/cc[m] directly):

```
for each row i of nlines:
    sv = cv = 0; k1 = 0
    for m = 0 .. nangles-1:
        j = *pt++
        if j ≠ 0:
            k = |j| - 1
            su = ss[m][k]; if j < 0: su = -su
            cu = cc[m][k]
            if k1 == 0: sv=su; cv=cu; k1=1
            else: ff=su*cv+cu*sv; cv=cu*cv-su*sv; sv=ff
    // accumulate based on typflg
```

**typflg dispatch:**

```
case 1 (large lon+rad, LR and LRT):
    j=*pt++; k=*pt++
    ans[0] += (10000.0 * j + k) * sv          // sin × (1″,0.0001″ packed)
    j=*pt++; k=*pt++
    if k ≠ 0: ans[2] += (10000.0 * j + k) * cv  // cos × (1km,0.0001km packed); skip if both zero

case 2 (small lon+rad, LRT2):
    j=*pt++; k=*pt++
    ans[0] += j * sv      // sin × 0.00001″
    ans[2] += k * cv      // cos × 0.00001km

case 3 (large lat, MB):
    j=*pt++; k=*pt++
    ans[1] += (10000.0 * j + k) * sv          // sin × (1″,0.0001″ packed)

case 4 (small lat, BT and BT2):
    j=*pt++
    ans[1] += j * sv                            // sin × 0.00001″
```

**Scale note on case 1**: The `if (k)` skip for radius applies when the 0.0001km column is zero, which in practice means the entire radius entry is zero. Do not skip when j≠0 but k=0 — check the actual table; all such rows have j=0 too.

### Unit Flow Through moon1–moon4

Understanding the unit pipeline is critical for correctness:

**moon1() — DE404 version (swemmoon.c:1182–1364)**

Phase A — T² terms:
1. `moonpol[0..2] = 0`
2. `chewm(LRT2, 25, 4, 2, moonpol)` → `moonpol[0]` in 0.00001″, `moonpol[2]` in 0.00001 km
3. `chewm(BT2, 12, 4, 4, moonpol)` → `moonpol[1]` in 0.00001″
4. Planetary perturbation terms: set `l` (arcsec T^0), `l1` (0.00001″/T units), `l2` (0.00001″/T² units) using z[12..24] and hard-coded amplitudes
5. `l2 += moonpol[0]` — absorb LRT2 longitude into l2
6. `moonpol[1] *= T; moonpol[2] *= T` — scale T² terms up by T

Phase B — T¹ terms:
7. `moonpol[0] = 0` (reset)
8. `chewm(BT, 16, 4, 4, moonpol)` → `moonpol[1]` adds T¹ lat in 0.00001″
9. `chewm(LRT, 38, 4, 1, moonpol)` → `moonpol[0]` in 0.0001″ (via 10000-scaling), `moonpol[2]` adds T¹ rad
10. Additional T¹ planetary perturbations: add to `l` and `l1`
11. `l1 += moonpol[0]` — absorb LRT longitude into l1
12. `a = 0.1 * T; moonpol[1] *= a; moonpol[2] *= a`

After moon1, `moonpol[1]` and `moonpol[2]` contain accumulated T-dependent corrections scaled such that multiplying by `1e-4` in moon3 gives arcseconds / km.

**moon2() (swemmoon.c:1367–1442)**

Evaluates ~28 additional T^0 planetary perturbation terms for longitude `l` (arcsec), and ~8 terms for latitude `B` (arcsec). All with hard-coded amplitudes and phase offsets. No table lookup — purely hand-coded sin expressions of combinations of Ve, Ea, Ma, Ju, Sa, SWELP, NF, D, MP.

**moon3() (swemmoon.c:1444–1454)**

Assembly step:
```
moonpol[0] = 0.0     // reset longitude only
chewm(LR, 118, 4, 1, moonpol)   // main lon (×10000): moonpol[0] in 0.0001″
chewm(MB, 77,  4, 3, moonpol)   // main lat (×10000): adds to moonpol[1]

// Polynomial in T for longitude, applied to l:
l += (((l4*T + l3)*T + l2)*T + l1) * T * 1.0e-5   // converts 0.00001″ → arcsec

moonpol[0] = SWELP + l + 1.0e-4 * moonpol[0]     // total longitude in arcsec
moonpol[1] = 1.0e-4 * moonpol[1] + B             // total latitude in arcsec
moonpol[2] = 1.0e-4 * moonpol[2] + 385000.52899  // total radius in km
```

After moon3: moonpol = [lon_arcsec, lat_arcsec, dist_km].

**moon4() (swemmoon.c:1458–1464)**

Unit conversion to Rust-friendly units:
```
moonpol[2] /= (AUNIT / 1000.0)          // km → AU
moonpol[0] = STR * mods3600(moonpol[0]) // arcsec → radians (in [0, 2π))
moonpol[1] = STR * moonpol[1]           // arcsec → radians
B = moonpol[1]                           // save for potential future use
```

Output: `moonpol = [lon_rad, lat_rad, dist_AU]` in ecliptic **of date** (not J2000).

### swi_moshmoon2 (swemmoon.c:848–862)

```
T = (J - J2000) / 36525.0
T2 = T * T
mean_elements()
mean_elements_pl()
moon1(); moon2(); moon3(); moon4()
pol[0..2] = moonpol[0..2]
```

Output: `pol[3]` = ecliptic polar of date [lon_rad, lat_rad, dist_AU]. **Not J2000!**

### swi_moshmoon (swemmoon.c:869–934)

Wrapper that adds coordinate transformation and speed.

**Range check:** `MOSHLUEPH_START - 0.2 ≤ tjd ≤ MOSHLUEPH_END + 0.2`.

**Coordinate transformation** (swemmoon.c:1722–1729, `ecldat_equ2000`):
```
swi_polcart(xpm, xpm)                           // polar → Cartesian
swi_coortrf2(xpm, xpm, -swed.oec.seps, swed.oec.ceps)  // ecliptic of date → equatorial of date
swi_precess(xpm, tjd, 0, J_TO_J2000)            // equatorial of date → J2000
```

Uses `swed.oec` (obliquity of date), **not** `swed.oec2000`. The Moon result is in ecliptic of date, requiring the date obliquity.

**Speed** — second-order central difference (swemmoon.c:915–929):
```
t1 = tjd + MOON_SPEED_INTV
t2 = tjd - MOON_SPEED_INTV
// compute x1 and x2 at t1 and t2 via moshmoon2 + ecldat_equ2000
for i = 0..2:
    b = (x1[i] - x2[i]) / 2
    a = (x1[i] + x2[i]) / 2 - xpm[i]
    xpm[i+3] = (2*a + b) / MOON_SPEED_INTV
```

This is a second-order corrected central difference: `v ≈ f′(t) + f″(t)·dt/2` where the curvature term `a` accounts for the second derivative.

**Caching:** if `tjd == pdp->teval && pdp->iephe == SEFLG_MOSEPH`, skip all computation.

**Output:** equatorial Cartesian J2000, 6-element [x, y, z, vx, vy, vz] in AU and AU/day.

### sscc() — Moon Version (swemmoon.c:1696–1714)

Identical algorithm to the planet version. Operates on `ss[5][8]` / `cc[5][8]`. Called in `moon1()`:
```
sscc(0, STR*D,  6)   // D, harmonics 1..6
sscc(1, STR*M,  4)   // M, harmonics 1..4
sscc(2, STR*MP, 4)   // MP, harmonics 1..4
sscc(3, STR*NF, 4)   // NF, harmonics 1..4
```

The `cc[4][0..7]` and `ss[4][0..7]` slots are unused (zeroed at start of moon1 in DE404 version).

**Bug fix in DE404 version**: moon1() explicitly zeros all ss and cc before calling sscc:
```c
for (i = 0; i < 5; i++)
    for (j = 0; j < 8; j++)
        ss[i][j] = cc[i][j] = 0;
```
This is required because not all 40 entries are initialized by the 4 sscc calls.

### mean_node_corr / mean_apsis_corr Tables

Two large `double[]` tables (swemmoon.c:725–805) providing century-step corrections for the mean lunar node and mean apsis respectively, covering roughly −13100 to +17200 Julian years (in 100-year Gregorian century steps). Used by `swi_mean_node()` and `swi_mean_apog()` to extend node/apogee accuracy outside the standard Moshier range via `corr_mean_node()` / `corr_mean_apog()`.

These corrections are zero for the interval 0–3000 AD (set explicitly in the array). Linear interpolation between adjacent century entries is used.

---

## Coordinate Frame Summary

| Function | Input epoch | Output frame |
|---|---|---|
| `swi_moshplan2` | J (any) | Heliocentric ecliptic J2000 polar (rad, rad, AU) |
| `swi_moshplan` | J (any) | Heliocentric equatorial J2000 Cartesian (AU, AU/day) |
| `swi_moshmoon2` | J (any) | Geocentric **ecliptic of date** polar (rad, rad, AU) |
| `swi_moshmoon` | J (any) | Geocentric equatorial J2000 Cartesian (AU, AU/day) |

## Not Porting

**swi_osc_el_plan**: Computes fictional bodies (Uranian planets, Isis-Transpluto, etc.) via osculating Keplerian elements. Entirely separate from the Moshier series. Port when fictional-body support is added.

**read_elements_file / check_t_terms**: File I/O and T-polynomial parsing for fictional body elements. File-backend concern.

**swi_intp_apsides**: Interpolated true lunar apsides — iterative root-finding that calls moon1–4 repeatedly. Separate task when osculating apogee support is added.

**swi_mean_node / swi_mean_apog**: Mean lunar node and apogee (Lilith). These depend only on `mean_elements()` plus a table lookup correction. Port when node/apogee support is needed.

**embofs_mosh caching**: The original C caches `pedp->teval` and `pedp->iephe` for Earth. In the stateless Rust design, always recompute.
