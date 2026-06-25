# C Reference: Delta-T — swephlib.c

Porting reference for the ΔT (ET − UT) computation. Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_deltat_ex` | swephlib.c:2701–2709 | Thin public wrapper |
| `swe_deltat` | swephlib.c:2712–2716 | Legacy wrapper (guesses ephe flag) |
| `calc_deltat` | swephlib.c:2545–2699 | Main dispatcher — yes |
| `deltat_aa` | swephlib.c:2733–2839 | 1620–present table + future — yes |
| `deltat_longterm_morrison_stephenson` | swephlib.c:2841–2846 | Long-term parabola — yes |
| `deltat_stephenson_morrison_1997_1600` | swephlib.c:2848–2887 | Historical model 1997 — yes |
| `deltat_stephenson_morrison_2004_1600` | swephlib.c:2890–2933 | Historical model 2004 — yes |
| `deltat_stephenson_etc_2016` | swephlib.c:3001–3036 | Historical model 2016 — yes |
| `deltat_espenak_meeus_1620` | swephlib.c:3038–3084 | Historical model EM2006 — yes |
| `adjust_for_tidacc` | swephlib.c:3143–3151 | Tidal correction — yes |
| `bessel` | swephlib.c:2004–2067 | Used by nutation EOP only — skip |
| `init_dt` | swephlib.c:3089–3133 | Reads swe_deltat.txt — skip |
| `swe_set_delta_t_userdef` | swephlib.c:3176–3183 | User override — skip |
| `swe_set_tid_acc` | swephlib.c:3165–3174 | Manual tidal acc — skip |
| `swe_get_tid_acc` | swephlib.c:3154–3157 | Getter — skip |
| `swi_get_tid_acc` | swephlib.c:3198–3240 | Internal tid_acc lookup — yes |
| `swi_set_tid_acc` | swephlib.c:3242–3261 | Internal tid_acc set — yes |

Output unit: **days** (all functions divide by 86400.0 at the end).

## Year Conversion Formulas

Three distinct Y formulas appear — use the right one for each context:

```
J2000 = 2451545.0   // TT noon, 2000 Jan 1 (sweph.h:67)

// Julian year, used in calc_deltat routing and some models
Y     = 2000.0 + (tjd - 2451545.0) / 365.25

// Gregorian year, used in most model functions (Gregorian leap-year spacing)
Ygreg = 2000.0 + (tjd - 2451545.0) / 365.2425

// Table year, used in deltat_aa for table indexing (note: 2451544.5, not J2000!)
Ytab  = 2000.0 + (tjd - 2451544.5) / 365.25
     = Y + 0.5/365.25  ≈ Y + 0.001369
```

The 0.5-day offset in `Ytab` vs `Y` means that table lookups land on Jan 1 midnight rather
than noon. This is a deliberate choice and does NOT affect results in practice beyond the
Bessel interpolation stencil position.

## Main Dispatcher (calc_deltat)

```
calc_deltat(tjd, iflag, *deltat, *serr):
    model = astro_models[SE_MODEL_DELTAT]
    if model == 0: model = SEMOD_DELTAT_DEFAULT  // = 5 = Stephenson2016

    // determine tid_acc from iflag / ephemeris file denum
    tid_acc = swi_get_tid_acc(...)   // see Tidal Acceleration section

    Y     = 2000.0 + (tjd - 2451545.0) / 365.25
    Ygreg = 2000.0 + (tjd - 2451545.0) / 365.2425

    // --- Model routing ---

    // DEFAULT: Stephenson/Morrison/Hohenkerk 2016
    if model == STEPHENSON_ETC_2016 && tjd < 2435108.5:  // before ~1 Jan 1955
        *deltat = deltat_stephenson_etc_2016(tjd, tid_acc)
        // blend 1000-day linear transition into tabulated region
        if tjd >= 2434108.5:   // within 1000 days before 2435108.5
            *deltat += (1.0 - (2435108.5 - tjd) / 1000.0) * 0.6610218 / 86400.0
        return

    // Espenak & Meeus 2006 (used in SE 1.77–2.05.01)
    if model == ESPENAK_MEEUS_2006 && tjd < 2317746.13090277789:  // before ~1633.74
        *deltat = deltat_espenak_meeus_1620(tjd, tid_acc)
        return

    // Stephenson & Morrison 2004 (used in SE 1.72–1.76)
    if model == STEPHENSON_MORRISON_2004 && Y < 1620:
        if Y < 1600:
            *deltat = deltat_stephenson_morrison_2004_1600(tjd, tid_acc)
        else:  // 1600 <= Y < 1620: linear blend from dt2 end to dt start
            B  = 20.0  // TABSTART - TAB2_END
            iy = 26    // (1600 - (-1000)) / 100
            dd = (Y - 1600.0) / B
            ans = dt2[26] + dd * (dt[0] - dt2[26])  // dt2[26]=120, dt[0]=124.0
            ans = adjust_for_tidacc(ans, Ygreg, tid_acc, SE_TIDAL_26, false)
            *deltat = ans / 86400.0
        return

    // Stephenson 1997 (used in SE 1.64–1.71)
    if model == STEPHENSON_1997 && Y < 1620:
        if Y < 1600:
            *deltat = deltat_stephenson_morrison_1997_1600(tjd, tid_acc)
        else:  // 1600 <= Y < 1620: linear blend from dt97 end to dt start
            B  = 20.0  // TABSTART - TAB97_END
            iy = 42    // (1600 - (-500)) / 50
            dd = (Y - 1600.0) / B
            ans = dt97[42] + dd * (dt[0] - dt97[42])  // dt97[42]=110, dt[0]=124.0
            ans = adjust_for_tidacc(ans, Ygreg, tid_acc, SE_TIDAL_26, false)
            *deltat = ans / 86400.0
        return

    // Stephenson & Morrison 1984 (used before SE 1.64)
    if model == STEPHENSON_MORRISON_1984 && Y < 1620:
        if Y >= 948.0:
            B = 0.01 * (Y - 2000.0)
            ans = (23.58 * B + 100.3) * B + 101.6   // seconds
        else:  // before 948 AD
            B = 0.01 * (Y - 2000.0) + 3.75  // = 0.01*(Y - 1625)
            ans = 35.0 * B * B + 40.0        // seconds
        *deltat = ans / 86400.0
        return

    // 1620–present: tabulated AA data + future extrapolation
    if Y >= 1620:
        *deltat = deltat_aa(tjd, tid_acc)
        return

    // (unreachable dead-fall for unknown/legacy model with Y < 1620)
    *deltat = 0.0 / 86400.0
```

**Key JD boundary values:**

| JD | Meaning |
|---|---|
| 2435108.5 | 1 Jan 1955 (switch from 2016 model to IERS table) |
| 2434108.5 | 1000 days before, = ~6 Apr 1952 (start of 2016→table blend) |
| 2317746.130902... | ~1633.74 Gregorian (Espenak-Meeus boundary) |

## Tabulated Delta-T: dt[] Array

`static TLS double dt[TABSIZ_SPACE]` — swephlib.c:2431

- **TABSTART = 1620** (first year)
- **TABEND = 2028** (last hard-coded year)
- **TABSIZ = 2028 − 1620 + 1 = 409** entries
- **TABSIZ_SPACE = 509** (room for 100 extra entries from `swe_deltat.txt` file)
- Units: **seconds** (divide by 86400.0 for days)
- Spacing: **one entry per year** (Jan 1 of each year)

Full table (swephlib.c:2432–2495):

```
/* 1620–1659 */
124.00, 119.00, 115.00, 110.00, 106.00, 102.00, 98.00, 95.00, 91.00, 88.00,
 85.00,  82.00,  79.00,  77.00,  74.00,  72.00, 70.00, 67.00, 65.00, 63.00,
 62.00,  60.00,  58.00,  57.00,  55.00,  54.00, 53.00, 51.00, 50.00, 49.00,
 48.00,  47.00,  46.00,  45.00,  44.00,  43.00, 42.00, 41.00, 40.00, 38.00,
/* 1660–1699 */
 37.00,  36.00,  35.00,  34.00,  33.00,  32.00, 31.00, 30.00, 28.00, 27.00,
 26.00,  25.00,  24.00,  23.00,  22.00,  21.00, 20.00, 19.00, 18.00, 17.00,
 16.00,  15.00,  14.00,  14.00,  13.00,  12.00, 12.00, 11.00, 11.00, 10.00,
 10.00,  10.00,   9.00,   9.00,   9.00,   9.00,  9.00,  9.00,  9.00,  9.00,
/* 1700–1739 */
  9.00,   9.00,   9.00,   9.00,   9.00,   9.00,  9.00,  9.00, 10.00, 10.00,
 10.00,  10.00,  10.00,  10.00,  10.00,  10.00, 10.00, 11.00, 11.00, 11.00,
 11.00,  11.00,  11.00,  11.00,  11.00,  11.00, 11.00, 11.00, 11.00, 11.00,
 11.00,  11.00,  11.00,  11.00,  12.00,  12.00, 12.00, 12.00, 12.00, 12.00,
/* 1740–1779 */
 12.00,  12.00,  12.00,  12.00,  13.00,  13.00, 13.00, 13.00, 13.00, 13.00,
 13.00,  14.00,  14.00,  14.00,  14.00,  14.00, 14.00, 14.00, 15.00, 15.00,
 15.00,  15.00,  15.00,  15.00,  15.00,  16.00, 16.00, 16.00, 16.00, 16.00,
 16.00,  16.00,  16.00,  16.00,  16.00,  17.00, 17.00, 17.00, 17.00, 17.00,
/* 1780–1799 */
 17.00,  17.00,  17.00,  17.00,  17.00,  17.00, 17.00, 17.00, 17.00, 17.00,
 17.00,  17.00,  16.00,  16.00,  16.00,  16.00, 15.00, 15.00, 14.00, 14.00,
/* 1800–1819 */
 13.70,  13.40,  13.10,  12.90,  12.70,  12.60, 12.50, 12.50, 12.50, 12.50,
 12.50,  12.50,  12.50,  12.50,  12.50,  12.50, 12.50, 12.40, 12.30, 12.20,
/* 1820–1859 */
 12.00,  11.70,  11.40,  11.10,  10.60,  10.20,  9.60,  9.10,  8.60,  8.00,
  7.50,   7.00,   6.60,   6.30,   6.00,   5.80,  5.70,  5.60,  5.60,  5.60,
  5.70,   5.80,   5.90,   6.10,   6.20,   6.30,  6.50,  6.60,  6.80,  6.90,
  7.10,   7.20,   7.30,   7.40,   7.50,   7.60,  7.70,  7.70,  7.80,  7.80,
/* 1860–1899 */
  7.88,   7.82,   7.54,   6.97,   6.40,   6.02,  5.41,  4.10,  2.92,  1.82,
  1.61,    .10,  -1.02,  -1.28,  -2.69,  -3.24, -3.64, -4.54, -4.71, -5.11,
 -5.40,  -5.42,  -5.20,  -5.46,  -5.46,  -5.79, -5.63, -5.64, -5.80, -5.66,
 -5.87,  -6.01,  -6.19,  -6.64,  -6.44,  -6.47, -6.09, -5.76, -4.66, -3.74,
/* 1900–1939 */
 -2.72,  -1.54,   -.02,   1.24,   2.64,   3.86,  5.37,  6.14,  7.75,  9.13,
 10.46,  11.53,  13.36,  14.65,  16.01,  17.20, 18.24, 19.06, 20.25, 20.95,
 21.16,  22.25,  22.41,  23.03,  23.49,  23.62, 23.86, 24.49, 24.34, 24.08,
 24.02,  24.00,  23.87,  23.95,  23.86,  23.93, 23.73, 23.92, 23.96, 24.02,
/* 1940–1949 */
 24.33,  24.83,  25.30,  25.70,  26.24,  26.77, 27.28, 27.78, 28.25, 28.71,
/* 1950–1959 */
 29.15,  29.57,  29.97,  30.36,  30.72,  31.07, 31.35, 31.68, 32.18, 32.68,
/* 1960–1969 */
 33.15,  33.59,  34.00,  34.47,  35.03,  35.73, 36.54, 37.43, 38.29, 39.20,
/* 1970–1979 */
 40.18,  41.17,  42.23,  43.37, 44.4841, 45.4761, 46.4567, 47.5214, 48.5344, 49.5862,
/* 1980–1989 */
 50.5387, 51.3808, 52.1668, 52.9565, 53.7882, 54.3427, 54.8713, 55.3222, 55.8197, 56.3000,
/* 1990–1999 */
 56.8553, 57.5653, 58.3092, 59.1218, 59.9845, 60.7854, 61.6287, 62.2951, 62.9659, 63.4673,
/* 2000–2009 */
 63.8285, 64.0908, 64.2998, 64.4734, 64.5736, 64.6876, 64.8452, 65.1464, 65.4574, 65.7768,
/* 2010–2019 */
 66.0699, 66.3246, 66.6030, 66.9069, 67.2810, 67.6439, 68.1024, 68.5927, 68.9676, 69.2202,
/* 2020–2023 */
 69.3612, 69.3593, 69.2945, 69.1833,
/* Extrapolated: 2024–2028 */
 69.10,   69.00,   68.90,   68.80,   68.80,
```

Total hard-coded entries: 409 (1620.0 through 2028.0 inclusive, one per integer year).

## deltat_aa: Tabulated Interpolation + Future Extrapolation

swephlib.c:2733–2839

Handles dates from 1620 to future.

### Year Used for Table Lookup

```
Ytab = 2000.0 + (tjd - 2451544.5) / 365.25   // NOTE: 2451544.5 (midnight), not J2000 (noon)
```

### Bessel 4th-Order Interpolation (1620 to tabend)

Called when `Ytab <= tabend` (last non-zero entry in dt[]).

```
tabsiz = init_dt()         // normally 409 (TABSIZ), may be larger if swe_deltat.txt loaded
tabend = TABSTART + tabsiz - 1   // = 2028 normally

p_int = floor(Ytab)
iy    = (int)(p_int - TABSTART)  // 0-based index into dt[]
p     = Ytab - p_int             // fractional part, 0..1

// Clamp: iy must be in [0, tabsiz-1]
// The caller guarantees Ytab <= tabend, so iy < tabsiz.
// iy can equal tabsiz-1 if Ytab == tabend exactly.

// Zeroth order
ans = dt[iy]

// First order
k = iy + 1
if k >= tabsiz: goto done    // no next entry, return dt[iy]
ans += p * (dt[k] - dt[iy])

// Guard: need iy-1 >= 0 AND iy+2 < tabsiz for second differences
if iy - 1 < 0 || iy + 2 >= tabsiz: goto done

// First differences: d[i] = dt[k+1] - dt[k] for k = iy-2..iy+2
// Out-of-range entries are zeroed (boundary padding)
k = iy - 2
for i = 0..4:
    if k < 0 || k+1 >= tabsiz: d[i] = 0
    else: d[i] = dt[k+1] - dt[k]
    k += 1

// Second differences: 4 values
for i = 0..3: d[i] = d[i+1] - d[i]

B = 0.25 * p * (p - 1.0)
ans += B * (d[1] + d[2])

// Guard for third differences
if iy + 2 >= tabsiz: goto done

// Third differences: 3 values
for i = 0..2: d[i] = d[i+1] - d[i]
B = (2.0 / 3.0) * B              // = (1/6) * p * (p-1)
ans += (p - 0.5) * B * d[1]

// Guard for fourth differences
if iy - 2 < 0 || iy + 3 > tabsiz: goto done

// Fourth differences: 2 values
for i = 0..1: d[i] = d[i+1] - d[i]
B = 0.125 * B * (p + 1.0) * (p - 2.0)
ans += B * (d[0] + d[1])

done:
ans = adjust_for_tidacc(ans, Ytab, tid_acc, SE_TIDAL_26, false)
return ans / 86400.0
```

**CRITICAL boundary conditions:**

| iy value | Order achieved | Reason |
|---|---|---|
| 0 (year 1620) | 1st order only | iy-1 < 0 fails the second-diff guard |
| 1 (year 1621) | 2nd order only | iy-2 < 0 fails the fourth-diff guard |
| 2 (year 1622) | 3rd order | iy-2 == 0 passes third-diff guard; fourth-diff also fails: iy+3 = 5 which is <= tabsiz |
| 2 ≤ iy ≤ tabsiz-4 | 4th order | all guards pass |
| tabsiz-3 | 3rd order | iy+2 == tabsiz-1, iy+3 = tabsiz-2 ≤ tabsiz fails fourth-diff condition |
| tabsiz-2 | 2nd order | iy+2 >= tabsiz fails third-diff guard |
| tabsiz-1 | 1st order | iy+2 >= tabsiz immediately |

Note: the out-of-range zero-padding in the first-differences loop means that even when 2nd
differences are computed, stencil entries that fall before index 0 or at/beyond tabsiz are
treated as zero slope. This is NOT a full symmetric pad — it's asymmetric clamping to zero.

### Future Extrapolation (Ytab > tabend)

Two formula branches based on active model:

**Branch A — Stephenson 2016 model:**
```
B = (Ytab - 2000)
if Ytab < 2500:
    ans = B³ * 121.0/30000000.0 + B² / 1250.0 + B * 521.0/3000.0 + 64.0
    B2  = (tabend - 2000)
    ans2 = B2³ * 121.0/30000000.0 + B2² / 1250.0 + B2 * 521.0/3000.0 + 64.0
else:  // Y >= 2500
    B = 0.01 * (Ytab - 2000)
    ans = B² * 32.5 + 42.5
```

**Branch B — all other models:**
```
B  = 0.01 * (Ytab - 1820)
ans = -20 + 31 * B²
B2 = 0.01 * (tabend - 1820)
ans2 = -20 + 31 * B2²
```

**Slow 100-year transition into tabulated region (both branches, Y < 2500 only):**
```
if Ytab <= tabend + 100:
    ans3 = dt[tabsiz - 1]       // last table entry
    dd   = ans2 - ans3           // gap between formula(tabend) and last table value
    ans  += dd * (Ytab - (tabend + 100)) * 0.01
// At Ytab == tabend:     correction = dd * (-100) * 0.01 = -dd → result = dt[tabsiz-1] ✓
// At Ytab == tabend+100: correction = 0               → result = formula only ✓
```

No tidal adjustment for future extrapolation (only tabulated region is adjusted).

## deltat_longterm_morrison_stephenson

swephlib.c:2841–2846

Parabola used as the deep-past baseline in the 2004 and EM2006 models.

```
Ygreg = 2000.0 + (tjd - 2451545.0) / 365.2425
u = (Ygreg - 1820) / 100.0
return (-20 + 32 * u²)   // seconds (NOT divided by 86400 — caller does it)
```

## Historical Model: StephensonMorrison1984

Handled inline in `calc_deltat` (swephlib.c:2656–2668), no separate function.

```
Y = 2000.0 + (tjd - 2451545.0) / 365.25   // Julian year

if Y >= 948.0:   // 948 AD to 1620
    B   = 0.01 * (Y - 2000.0)
    ans = (23.58 * B + 100.3) * B + 101.6    // seconds
else:            // before 948 AD
    B   = 0.01 * (Y - 2000.0) + 3.75         // = 0.01 * (Y - 1625)
    ans = 35.0 * B² + 40.0                   // seconds

*deltat = ans / 86400.0
```

No tidal correction applied.

## Historical Model: Stephenson1997 (deltat_stephenson_morrison_1997_1600)

swephlib.c:2848–2887

Tables: `dt97[]` — 50-year intervals from −500 to 1600.

```
Y = 2000.0 + (tjd - 2451545.0) / 365.25    // Julian year

// Before -500: Stephenson 1997 p.508 parabola (adjusted)
if Y < -500:
    B   = (Y - 1735) * 0.01
    ans = -20 + 35 * B²
    ans = adjust_for_tidacc(ans, Y, tid_acc, SE_TIDAL_26, false)
    // Blend 100-year transition to match table starting value
    if Y >= -600:   // within 100 years of table start
        ans2 = adjust_for_tidacc(dt97[0], -500.0, tid_acc, SE_TIDAL_26, false)
        B3   = (-500 - 1735) * 0.01
        ans3 = -20 + 35 * B3²
        ans3 = adjust_for_tidacc(ans3, Y, tid_acc, SE_TIDAL_26, false)
        dd   = ans3 - ans2
        B    = (Y - (-600)) * 0.01
        ans  = ans - dd * B

// -500 to 1600: linear interpolation within dt97[]
if -500 <= Y < 1600:
    iy = (int)((floor(Y) - (-500)) / 50.0)
    dd = (Y - (-500 + 50 * iy)) / 50.0
    ans = dt97[iy] + (dt97[iy+1] - dt97[iy]) * dd
    ans = adjust_for_tidacc(ans, Y, tid_acc, SE_TIDAL_26, false)

return ans / 86400.0
```

### dt97[] Table (Stephenson & Morrison 1995)

swephlib.c:2523–2533 — `static const short dt97[43]`

- Range: −500 to 1600, step 50 years
- Units: **seconds**
- 43 entries (index 0 = year −500, index 42 = year 1600)

```
idx  year    sec
 0   -500   16800
 1   -450   16000
 2   -400   15300
 3   -350   14600
 4   -300   14000
 5   -250   13400
 6   -200   12800
 7   -150   12200
 8   -100   11600
 9    -50   11100
10      0   10600
11     50   10100
12    100    9600
13    150    9100
14    200    8600
15    250    8200
16    300    7700
17    350    7200
18    400    6700
19    450    6200
20    500    5700
21    550    5200
22    600    4700
23    650    4300
24    700    3800
25    750    3400
26    800    3000
27    850    2600
28    900    2200
29    950    1900
30   1000    1600
31   1050    1350
32   1100    1100
33   1150     900
34   1200     750
35   1250     600
36   1300     470
37   1350     380
38   1400     300
39   1450     230
40   1500     180
41   1550     140
42   1600     110
```

Note: the C comment says "first value for -550 added from Borkowski" — but the array is indexed
from −500 (TAB97_START). The Borkowski comment refers to the formula used *before* −500.

## Historical Model: StephensonMorrison2004 (deltat_stephenson_morrison_2004_1600)

swephlib.c:2890–2933

```
Y = 2000.0 + (tjd - 2451545.0) / 365.2425    // Gregorian year (note: 365.2425)

// Before -1000: long-term parabola (adjusted)
if Y < -1000:
    ans  = deltat_longterm_morrison_stephenson(tjd)   // returns seconds
    ans  = adjust_for_tidacc(ans, Y, tid_acc, SE_TIDAL_26, false)
    // Blend 100-year transition to table start
    if Y >= -1100:
        ans2 = adjust_for_tidacc(dt2[0], -1000.0, tid_acc, SE_TIDAL_26, false)
        tjd0 = (-1000 - 2000) * 365.2425 + 2451545.0
        ans3 = deltat_longterm_morrison_stephenson(tjd0)
        ans3 = adjust_for_tidacc(ans3, Y, tid_acc, SE_TIDAL_26, false)
        dd   = ans3 - ans2
        B    = (Y - (-1100)) * 0.01
        ans  = ans - dd * B

// -1000 to 1600: linear interpolation within dt2[]
// NOTE: uses Julian year for indexing, not Y (which is Gregorian)!
if -1000 <= Y < 1600:
    Yjul = 2000 + (tjd - 2451557.5) / 365.25   // Julian year (2451557.5 = 2000 Jan 13.0)
    iy   = (int)((floor(Yjul) - (-1000)) / 100.0)
    dd   = (Yjul - (-1000 + 100 * iy)) / 100.0
    ans  = dt2[iy] + (dt2[iy+1] - dt2[iy]) * dd
    ans  = adjust_for_tidacc(ans, Y, tid_acc, SE_TIDAL_26, false)

return ans / 86400.0
```

**Note the mixed-year bug:** outer routing uses `Y` (Gregorian), but table indexing uses `Yjul`
(Julian, with base JD 2451557.5 rather than J2000). Reproduce this exactly.

### dt2[] Table (Morrison & Stephenson 2004)

swephlib.c:2504–2511 — `static const short dt2[27]`

- Range: −1000 to 1600, step 100 years
- Units: **seconds**
- 27 entries (index 0 = year −1000, index 26 = year 1600)

```
idx  year    sec
 0  -1000   25400
 1   -900   23700
 2   -800   22000
 3   -700   21000
 4   -600   19040
 5   -500   17190
 6   -400   15530
 7   -300   14080
 8   -200   12790
 9   -100   11640
10      0   10580
11    100    9600
12    200    8640
13    300    7680
14    400    6700
15    500    5710
16    600    4740
17    700    3810
18    800    2960
19    900    2200
20   1000    1570
21   1100    1090
22   1200     740
23   1300     490
24   1400     320
25   1500     200
26   1600     120
```

## Historical Model: StephensonEtc2016 (deltat_stephenson_etc_2016)

swephlib.c:3001–3036

### Cubic Spline Table (dtcf16)

swephlib.c:2943–3000 — `double dtcf16[54][6]`

54 records. Each record: `[jd_beg, jd_end, c0, c1, c2, c3]`

Evaluation: given `tjd` in `[jd_beg, jd_end)`:
```
t  = (tjd - jd_beg) / (jd_end - jd_beg)   // normalized 0..1
dt = c0 + c1*t + c2*t² + c3*t³             // seconds
```

Full table:

```
idx   jd_beg        jd_end       c0          c1          c2          c3        years
 0  1458085.5   1867156.5   20550.593  -21268.478   11863.418   -4541.129   -720.. 400
 1  1867156.5   2086302.5    6604.404   -5981.266    -505.093    1349.609    400..1000
 2  2086302.5   2268923.5    1467.654   -2452.187    2460.927   -1183.759   1000..1500
 3  2268923.5   2305447.5     292.635    -216.322     -43.614      56.681   1500..1600
 4  2305447.5   2323710.5      89.380     -66.754      31.607     -10.497   1600..1650
 5  2323710.5   2349276.5      43.736     -49.043       0.227      15.811   1650..1720
 6  2349276.5   2378496.5      10.730      -1.321      62.250     -52.946   1720..1800
 7  2378496.5   2382148.5      18.714      -4.457      -1.509       2.507   1800..1810
 8  2382148.5   2385800.5      15.255       0.046       6.012      -4.634   1810..1820
 9  2385800.5   2389453.5      16.679      -1.831      -7.889       3.799   1820..1830
10  2389453.5   2393105.5      10.758      -6.211       3.509      -0.388   1830..1840
11  2393105.5   2396758.5       7.668      -0.357       2.345      -0.338   1840..1850
12  2396758.5   2398584.5       9.317       1.659       0.332      -0.932   1850..1855
13  2398584.5   2400410.5      10.376      -0.472      -2.463       1.596   1855..1860
14  2400410.5   2402237.5       9.038      -0.610       2.325      -2.497   1860..1865
15  2402237.5   2404063.5       8.256      -3.450      -5.166       2.729   1865..1870
16  2404063.5   2405889.5       2.369      -5.596       3.020      -0.919   1870..1875
17  2405889.5   2407715.5      -1.126      -2.312       0.264      -0.037   1875..1880
18  2407715.5   2409542.5      -3.211      -1.894       0.154       0.562   1880..1885
19  2409542.5   2411368.5      -4.388       0.101       1.841      -1.438   1885..1890
20  2411368.5   2413194.5      -3.884      -0.531      -2.473       1.870   1890..1895
21  2413194.5   2415020.5      -5.017       0.134       3.138      -0.232   1895..1900
22  2415020.5   2416846.5      -1.977       5.715       2.443      -1.257   1900..1905
23  2416846.5   2418672.5       4.923       6.828      -1.329       0.720   1905..1910
24  2418672.5   2420498.5      11.142       6.330       0.831      -0.825   1910..1915
25  2420498.5   2422324.5      17.479       5.518      -1.643       0.262   1915..1920
26  2422324.5   2424151.5      21.617       3.020      -0.856       0.008   1920..1925
27  2424151.5   2425977.5      23.789       1.333      -0.831       0.127   1925..1930
28  2425977.5   2427803.5      24.418       0.052      -0.449       0.142   1930..1935
29  2427803.5   2429629.5      24.164      -0.419      -0.022       0.702   1935..1940
30  2429629.5   2431456.5      24.426       1.645       2.086      -1.106   1940..1945
31  2431456.5   2433282.5      27.050       2.499      -1.232       0.614   1945..1950
32  2433282.5   2434378.5      28.932       1.127       0.220      -0.277   1950..1953
33  2434378.5   2435473.5      30.002       0.737      -0.610       0.631   1953..1956
34  2435473.5   2436569.5      30.760       1.409       1.282      -0.799   1956..1959
35  2436569.5   2437665.5      32.652       1.577      -1.115       0.507   1959..1962
36  2437665.5   2438761.5      33.621       0.868       0.406       0.199   1962..1965
37  2438761.5   2439856.5      35.093       2.275       1.002      -0.414   1965..1968
38  2439856.5   2440952.5      37.956       3.035      -0.242       0.202   1968..1971
39  2440952.5   2442048.5      40.951       3.157       0.364      -0.229   1971..1974
40  2442048.5   2443144.5      44.244       3.198      -0.323       0.172   1974..1977
41  2443144.5   2444239.5      47.291       3.069       0.193      -0.192   1977..1980
42  2444239.5   2445335.5      50.361       2.878      -0.384       0.081   1980..1983
43  2445335.5   2446431.5      52.936       2.354      -0.140      -0.166   1983..1986
44  2446431.5   2447527.5      54.984       1.577      -0.637       0.448   1986..1989
45  2447527.5   2448622.5      56.373       1.649       0.709      -0.277   1989..1992
46  2448622.5   2449718.5      58.453       2.235      -0.122       0.111   1992..1995
47  2449718.5   2450814.5      60.677       2.324       0.212      -0.315   1995..1998
48  2450814.5   2451910.5      62.899       1.804      -0.732       0.112   1998..2001
49  2451910.5   2453005.5      64.082       0.675      -0.396       0.193   2001..2004
50  2453005.5   2454101.5      64.555       0.463       0.184      -0.008   2004..2007
51  2454101.5   2455197.5      65.194       0.809       0.161      -0.101   2007..2010
52  2455197.5   2456293.5      66.063       0.828      -0.142       0.168   2010..2013
53  2456293.5   2457388.5      66.917       1.046       0.360      -0.282   2013..2016
```

### Algorithm

```
Ygreg = 2000.0 + (tjd - 2451545.0) / 365.2425

irec = -1
for i = 0..53:
    if tjd < dtcf16[i][0]: break   // past all intervals
    if tjd < dtcf16[i][1]:
        irec = i
        break

if irec >= 0:   // within spline range -720 to 2016
    t  = (tjd - dtcf16[irec][0]) / (dtcf16[irec][1] - dtcf16[irec][0])
    dt = c0 + c1*t + c2*t² + c3*t³

elif Ygreg < -720:   // before -720
    t  = (Ygreg - 1825) / 100.0
    dt = -320 + 32.5 * t²
    dt -= 179.7337208    // continuity offset at 1 Jan -720

else:   // future (after last spline segment, Ygreg >= ~2016)
    t  = (Ygreg - 1825) / 100.0
    dt = -320 + 32.5 * t²
    dt += 269.4790417    // continuity offset at 1 Jan 2016

dt = adjust_for_tidacc(dt, Ygreg, tid_acc, SE_TIDAL_STEPHENSON_2016, true)
return dt / 86400.0
```

**Note:** `adjust_after_1955 = true` here because the Stephenson 2016 spline is based on
occultation data alone, not IERS. So the tidal correction applies to ALL epochs including
post-1955.

This function is only called from `calc_deltat` for `tjd < 2435108.5` (before 1955). The
future branch (`Ygreg >= ~2016`) is unreachable via the normal call path.

## Historical Model: EspenakMeeus2006 (deltat_espenak_meeus_1620)

swephlib.c:3038–3084

7 piecewise polynomials covering −∞ to 2005, all in Gregorian year.

```
Ygreg = 2000.0 + (tjd - 2451545.0) / 365.2425

if Ygreg < -500:
    ans = deltat_longterm_morrison_stephenson(tjd)  // parabola, returns seconds

elif Ygreg < 500:
    u   = Ygreg / 100.0
    ans = (((((0.0090316521 * u + 0.022174192) * u - 0.1798452) * u
            - 5.952053) * u + 33.78311) * u - 1014.41) * u + 10583.6

elif Ygreg < 1600:
    u   = (Ygreg - 1000) / 100.0
    ans = (((((0.0083572073 * u - 0.005050998) * u - 0.8503463) * u
            + 0.319781) * u + 71.23472) * u - 556.01) * u + 1574.2

elif Ygreg < 1700:
    u   = Ygreg - 1600
    ans = 120 - 0.9808 * u - 0.01532 * u² + u³ / 7129.0

elif Ygreg < 1800:
    u   = Ygreg - 1700
    ans = (((-u / 1174000.0 + 0.00013336) * u - 0.0059285) * u + 0.1603) * u + 8.83

elif Ygreg < 1860:
    u   = Ygreg - 1800
    ans = ((((((0.000000000875 * u - 0.0000001699) * u + 0.0000121272) * u
            - 0.00037436) * u + 0.0041116) * u + 0.0068612) * u - 0.332447) * u + 13.72

elif Ygreg < 1900:
    u   = Ygreg - 1860
    ans = ((((u / 233174.0 - 0.0004473624) * u + 0.01680668) * u
           - 0.251754) * u + 0.5737) * u + 7.62

elif Ygreg < 1920:
    u   = Ygreg - 1900
    ans = (((-0.000197 * u + 0.0061966) * u - 0.0598939) * u + 1.494119) * u - 2.79

elif Ygreg < 1941:
    u   = Ygreg - 1920
    ans = 21.20 + 0.84493 * u - 0.076100 * u² + 0.0020936 * u³

elif Ygreg < 1961:
    u   = Ygreg - 1950
    ans = 29.07 + 0.407 * u - u² / 233.0 + u³ / 2547.0

elif Ygreg < 1986:
    u   = Ygreg - 1975
    ans = 45.45 + 1.067 * u - u² / 260.0 - u³ / 718.0

elif Ygreg < 2005:
    u   = Ygreg - 2000
    ans = ((((0.00002373599 * u + 0.000651814) * u + 0.0017275) * u
           - 0.060374) * u + 0.3345) * u + 63.86

// else: Ygreg >= 2005 → ans = 0 (C falls through with ans still 0)
// But this function is only called for tjd < 2317746.13 (≈ 1633.74 Gregorian),
// so Ygreg >= 1700 branches are dead code in practice.

ans = adjust_for_tidacc(ans, Ygreg, tid_acc, SE_TIDAL_26, false)
return ans / 86400.0
```

**Polynomial summary table:**

| Range | Offset variable u | Degree |
|---|---|---|
| −∞ to −500 | long-term parabola | 2 |
| −500 to 500 | u = Ygreg/100 | 6 |
| 500 to 1600 | u = (Ygreg−1000)/100 | 6 |
| 1600 to 1700 | u = Ygreg−1600 | 3 |
| 1700 to 1800 | u = Ygreg−1700 | 4 |
| 1800 to 1860 | u = Ygreg−1800 | 7 |
| 1860 to 1900 | u = Ygreg−1860 | 5 |
| 1900 to 1920 | u = Ygreg−1900 | 4 |
| 1920 to 1941 | u = Ygreg−1920 | 3 |
| 1941 to 1961 | u = Ygreg−1950 | 3 |
| 1961 to 1986 | u = Ygreg−1975 | 3 |
| 1986 to 2005 | u = Ygreg−2000 | 5 |

All ans values are in **seconds**.

## Tidal Acceleration Correction (adjust_for_tidacc)

swephlib.c:3143–3151

Corrects ΔT for the difference between the tidal acceleration embedded in the historical
observations and the one used by the target ephemeris.

Formula (AA page K8):
```
adjust_for_tidacc(ans, Y, tid_acc, tid_acc0, adjust_after_1955):
    if Y < 1955.0 || adjust_after_1955:
        B    = Y - 1955.0
        ans += -0.000091 * (tid_acc - tid_acc0) * B²
    return ans
// Units: ans in seconds in, seconds out
```

Where:
- `tid_acc` = tidal acceleration of the target ephemeris (arcsec/cty²)
- `tid_acc0` = tidal acceleration of the reference data set (arcsec/cty²)
- `B = year - 1955.0` — distance from the reference epoch
- Coefficient `−0.000091` has units: seconds / (arcsec/cty² × yr²)

**`adjust_after_1955` flag:** normally `false` (only correct pre-1955 data). Set `true` for
Stephenson 2016 because that model is based on occultation data throughout and does not use
IERS for any epoch.

### Tidal Acceleration Constants

From swephexp.h:478–493:

| Constant | Value (arcsec/cty²) | Notes |
|---|---|---|
| `SE_TIDAL_DE200` | −23.8946 | JPL DE200 |
| `SE_TIDAL_DE403` | −25.580 | JPL DE403/404 |
| `SE_TIDAL_DE404` | −25.580 | |
| `SE_TIDAL_DE405` | −25.826 | JPL DE405/406 |
| `SE_TIDAL_DE406` | −25.826 | |
| `SE_TIDAL_DE421` | −25.85 | JPL DE421/422 |
| `SE_TIDAL_DE422` | −25.85 | |
| `SE_TIDAL_DE430` | −25.82 | JPL DE430 |
| `SE_TIDAL_DE431` | −25.80 | JPL DE431 |
| `SE_TIDAL_DE441` | −25.936 | JPL DE440/441 |
| `SE_TIDAL_26` | −26.0 | Reference value for historical tables |
| `SE_TIDAL_STEPHENSON_2016` | −25.85 | Reference for 2016 spline |
| `SE_TIDAL_DEFAULT` | = `SE_TIDAL_DE431` = −25.80 | Default ephemeris |

**`tid_acc0` parameter values used in calls:**
- Most historical models: `SE_TIDAL_26` (−26.0) — the AA tables assume ndot = −26
- Stephenson 2016: `SE_TIDAL_STEPHENSON_2016` (−25.85)

### Tidal Acceleration Selection (swi_get_tid_acc)

swephlib.c:3198–3240

```
if manual override set: return swed.tid_acc
if iflag & SEFLG_MOSEPH: return SE_TIDAL_DE404, denum=404
if iflag & SEFLG_JPLEPH && file open: denum from jpldenum
if iflag & SEFLG_SWIEPH && file open: denum from sweph_denum
switch denum:
    200 → DE200,  403/404 → DE403,  405/406 → DE405
    421/422 → DE421,  430 → DE430,  431 → DE431
    440/441 → DE441
    default → SE_TIDAL_DEFAULT (DE431), denum = SE_DE_NUMBER
```

For Moshier ephemeris (SEFLG_MOSEPH): `tid_acc = SE_TIDAL_DE404 = −25.580`.
Default/fallback: `tid_acc = SE_TIDAL_DEFAULT = SE_TIDAL_DE431 = −25.80`.

## Model Enum Values

From swephexp.h:600–606:

| Rust variant | C constant | Value |
|---|---|---|
| `StephensonMorrison1984` | `SEMOD_DELTAT_STEPHENSON_MORRISON_1984` | 1 |
| `Stephenson1997` | `SEMOD_DELTAT_STEPHENSON_1997` | 2 |
| `StephensonMorrison2004` | `SEMOD_DELTAT_STEPHENSON_MORRISON_2004` | 3 |
| `EspenakMeeus2006` | `SEMOD_DELTAT_ESPENAK_MEEUS_2006` | 4 |
| `StephensonEtc2016` | `SEMOD_DELTAT_STEPHENSON_ETC_2016` | 5 |

Default: `SEMOD_DELTAT_DEFAULT = SEMOD_DELTAT_STEPHENSON_ETC_2016 = 5`

## Constants Summary

| Name | Value | Units | Source |
|---|---|---|---|
| J2000 | 2451545.0 | JD | sweph.h:67 |
| J1900 | 2415020.0 | JD | sweph.h:69 |
| TABSTART | 1620 | year | swephlib.c:2426 |
| TABEND | 2028 | year | swephlib.c:2427 |
| TABSIZ | 409 | entries | swephlib.c:2428 |
| TABSIZ_SPACE | 509 | entries | swephlib.c:2430 |
| TAB2_SIZ | 27 | entries | swephlib.c:2497 |
| TAB2_START | −1000 | year | swephlib.c:2498 |
| TAB2_END | 1600 | year | swephlib.c:2499 |
| TAB2_STEP | 100 | years | swephlib.c:2500 |
| TAB97_SIZ | 43 | entries | swephlib.c:2519 |
| TAB97_START | −500 | year | swephlib.c:2520 |
| TAB97_END | 1600 | year | swephlib.c:2521 |
| TAB97_STEP | 50 | years | swephlib.c:2522 |
| NDTCF16 | 54 | spline segments | swephlib.c:2943 |
| JD_2016_boundary | 2435108.5 | JD | swephlib.c:2589 |
| JD_2016_blend_start | 2434108.5 | JD | swephlib.c:2591 |
| JD_EM2006_boundary | 2317746.13090277789 | JD (≈ 1633.74) | swephlib.c:2603 |
| blend_correction | 0.6610218 | seconds | swephlib.c:2592 |
| LTERM_EQUATION_YSTART | 1820 | year | swephlib.c:2501 |
| LTERM_EQUATION_COEFF | 32 | — | swephlib.c:2502 |

## Not Porting

**`init_dt` / `swe_deltat.txt`:** reads additional ΔT values from an external file to extend
the table beyond TABEND. In Rust, the table is static; no file extension. If the table is
sufficient for the use case, skip. If not, provide a different injection mechanism.

**`swe_set_delta_t_userdef`:** user-supplied constant ΔT override (swephlib.c:3176). Skip.
The stateless design means no global override. Expose as an `Option<f64>` on `EphemerisConfig`
if needed.

**`swe_set_tid_acc` / `swe_get_tid_acc`:** global mutable state. `swi_get_tid_acc` logic
(the table lookup by denum) should be ported as a pure function. The manual override path
maps to `EphemerisConfig`.

**Tracing / DEMO paths:** `#if DEMO` printf calls and `#ifdef TRACE` blocks throughout.
Skip entirely.
