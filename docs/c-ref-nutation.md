# C Reference: Nutation — swephlib.c

Porting reference for nutation algorithms. Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swi_nutation` | swephlib.c:2126–2158 | Router only (skip interpolation cache) |
| `calc_nutation` | swephlib.c:2069–2113 | Model dispatch + JPLHOR paths |
| `calc_nutation_iau1980` | swephlib.c:1615–1763 | Yes — 105-term series |
| `calc_nutation_iau2000ab` | swephlib.c:1813–1943 | Yes — shared 2000A/2000B |
| `calc_nutation_woolard` | swephlib.c:1947–2001 | Yes — simple polynomial |
| `nut_matrix` | sweph.c:5073 | Not this task (downstream) |
| `swi_nutate` | sweph.c:3592 | Not this task (downstream) |

## Model Dispatch (calc_nutation)

```
if is_jplhor (DPSIDEPS_1980 flag):
    calc_nutation_iau1980(J, nut_model=IAU1980)
    // add EOP corrections — STUB (needs file I/O)

elif is_jplhor_approx (JPLHOR_APPROX flag + JPLHORA_V3 + J ≤ 2437684.5):
    calc_nutation_iau1980(J, nut_model=IAU1980)
    dpsi += DPSI_IAU1980_TJD0 / 3600.0 * DEGTORAD
    deps += DEPS_IAU1980_TJD0 / 3600.0 * DEGTORAD

elif is_jplhor_approx + JPLHORA_V2 (IAU 2000 path):
    calc_nutation_iau2000ab(J, ...)
    dpsi += -41.7750e-3 / 3600.0 * DEGTORAD   // mas → arcsec → rad
    deps += -6.8192e-3 / 3600.0 * DEGTORAD

else:  // normal dispatch
    match model:
        IAU1980     → calc_nutation_iau1980(J, nut_model=IAU1980)
        IAUCorr1987 → calc_nutation_iau1980(J, nut_model=IAUCorr1987)  // enables Herring rows
        IAU2000A    → calc_nutation_iau2000ab(J, is_2000a=true)
        IAU2000B    → calc_nutation_iau2000ab(J, is_2000a=false)
        Woolard     → calc_nutation_woolard(J)
```

Output: `nutlo[0]` = Δψ (radians), `nutlo[1]` = Δε (radians).

## IAU 1980 Algorithm (calc_nutation_iau1980)

### Delaunay Arguments (Seidelmann 1982 / FK5)

T = (J − 2451545.0) / 36525.0

All computed in arcseconds with this C evaluation form (NOT Horner):
```c
arg = c1*T + c0 + (c3*T + c2)*T*T
```

| Arg | c0 | c1 | c2 | c3 |
|---|---|---|---|---|
| OM (node) | 450160.280 | −6962890.539 | 7.455 | 0.008 |
| MS (Sun anom) | 1287099.804 | 129596581.224 | −0.577 | −0.012 |
| MM (Moon anom) | 485866.733 | 1717915922.633 | 31.310 | 0.064 |
| FF (Moon lat) | 335778.877 | 1739527263.137 | −13.257 | 0.011 |
| DD (elongation) | 1072261.307 | 1602961601.328 | −6.891 | 0.019 |

Conversion: `arg_rad = normalize_degrees(arg_arcsec / 3600.0) * DEGTORAD`

### Sin/Cos Multiple Precomputation

Max multiples: MM→3, MS→2, FF→4, DD→4, OM→2

For each fundamental arg, precompute sin/cos of multiples 1..max:
```
ss[0] = sin(arg),         cc[0] = cos(arg)
ss[1] = 2*ss[0]*cc[0],   cc[1] = cc[0]*cc[0] - ss[0]*ss[0]
ss[n] = ss[0]*cc[n-1] + cc[0]*ss[n-1]
cc[n] = cc[0]*cc[n-1] - ss[0]*ss[n-1]
```

Index i stores sin/cos of (i+1)*arg. To get sin/cos for multiplier m:
- m = 0: sin=0, cos=1 (skip in combination)
- m > 0: ss[m-1], cc[m-1]
- m < 0: −ss[|m|−1], cc[|m|−1]

### Combining Arguments

For a table row with multipliers [mm, ms, ff, dd, om], compute sin/cos of the combined angle via successive angle addition:
```
sin_acc = 0, cos_acc = 1
for each arg k with multiplier m_k ≠ 0:
    s_k = ss_k[|m_k| - 1] * sign(m_k)
    c_k = cc_k[|m_k| - 1]
    new_sin = sin_acc * c_k + cos_acc * s_k
    new_cos = cos_acc * c_k - sin_acc * s_k
    sin_acc = new_sin, cos_acc = new_cos
```

### Hard-Coded Dominant Term

Before the table loop, seed the accumulators with the dominant OM term (NOT in the NT table):
```
dpsi = (-0.01742*T - 17.1996) * sin(OM)    // arcsec
deps = ( 0.00089*T +  9.2025) * cos(OM)    // arcsec
```

### Table Loop (NT array, 112 rows)

Each row: `[MM, MS, FF, DD, OM, LS, LS2, OC, OC2]` as i16.

Standard rows (MM < 100):
```
f = LS * 0.0001 + (if LS2 ≠ 0: 0.00001 * T * LS2, else 0)   // dpsi amplitude (arcsec)
g = OC * 0.0001 + (if OC2 ≠ 0: 0.00001 * T * OC2, else 0)   // deps amplitude (arcsec)
dpsi += f * sin(combined_arg)
deps += g * cos(combined_arg)
```

Herring rows (MM ≥ 100, skip unless model == IAUCorr1987):
- Actual MM multiplier = **0** (C code: `if (j > 100) j = 0; /* p[0] is a flag */`)
- Coefficients are 10× finer: `f *= 0.1; g *= 0.1`
- MM == 102: swap trig — use cos for dpsi, sin for deps
- MM == 101: normal trig (sin for dpsi, cos for deps)

### Final Conversion

```
nutlo[0] = dpsi * DEGTORAD / 3600.0    // arcsec → radians
nutlo[1] = deps * DEGTORAD / 3600.0
```

## IAU 2000A/B Algorithm (calc_nutation_iau2000ab)

### Luni-Solar Delaunay Arguments (Simon et al. 1994)

T = (J − 2451545.0) / 36525.0

Computed in arcseconds, 4th-degree polynomial (Horner form in C):
```c
M = c0 + T*(c1 + T*(c2 + T*(c3 + T*c4)))
```

| Arg | c0 | c1 | c2 | c3 | c4 |
|---|---|---|---|---|---|
| M (Moon anom) | 485868.249036 | 1717915923.2178 | 31.8792 | 0.051635 | −0.00024470 |
| SM (Sun anom) | 1287104.79305 | 129596581.0481 | −0.5532 | 0.000136 | −0.00001149 |
| F (Moon lat) | 335779.526232 | 1739527262.8478 | −12.7512 | −0.001037 | 0.00000417 |
| D (elongation) | 1072260.70369 | 1602961601.2090 | −6.3706 | 0.006593 | −0.00003169 |
| OM (node) | 450160.398036 | −6962890.5431 | 7.4722 | 0.007702 | −0.00005939 |

Conversion: same as IAU 1980 (`normalize_degrees(arcsec / 3600.0) * DEGTORAD`).

### Planetary Arguments (linear, directly in radians)

Only used for IAU 2000A (NPL terms). 14 arguments, 1st-degree polynomial in T:

| Arg | a0 (rad) | a1 (rad/century) |
|---|---|---|
| AL (Moon anom) | 2.35555598 | 8328.6914269554 |
| ALSU (Sun anom) | 6.24006013 | 628.301955 |
| AF (Moon lat) | 1.627905234 | 8433.466158131 |
| AD (elongation) | 5.198466741 | 7771.3771468121 |
| AOM (node) | 2.18243920 | −33.757045 |
| ALME (Mercury) | 4.402608842 | 2608.7903141574 |
| ALVE (Venus) | 3.176146697 | 1021.3285546211 |
| ALEA (Earth) | 1.753470314 | 628.3075849991 |
| ALMA (Mars) | 6.203480913 | 334.0612426700 |
| ALJU (Jupiter) | 0.599546497 | 52.9690962641 |
| ALSA (Saturn) | 0.874016757 | 21.3299104960 |
| ALUR (Uranus) | 5.481293871 | 7.4781598567 |
| ALNE (Neptune) | 5.321159000 | 3.8127774000 |
| APA (precession) | 0.02438175*T + 0.00000538691*T² | — (2nd degree) |

All normalized via `normalize_radians` (mod 2π).

Note: planetary args use deliberately different coefficients from luni-solar args (reproducing MHB2000 code behavior).

### Luni-Solar Series Loop

Tables: NLS (5 multipliers per row) + CLS (6 coefficients per row, i32, in units of 0.1 μas = 1e-7 arcsec).

CLS layout: `[dpsi_sin, dpsi_sin_t, dpsi_cos, deps_cos, deps_cos_t, deps_sin]`

Term count: 678 (2000A) or 77 (2000B). **Loop runs in reverse** (numerical stability).

```
for each row (reverse order):
    darg = normalize_radians(nls[0]*M + nls[1]*SM + nls[2]*F + nls[3]*D + nls[4]*OM)
    dpsi += (cls[0] + cls[1]*T) * sin(darg) + cls[2] * cos(darg)
    deps += (cls[3] + cls[4]*T) * cos(darg) + cls[5] * sin(darg)
```

Accumulator units: 0.1 μas (integer-scaled).

### Planetary Series (2000A only, 687 terms)

Tables: NPL (14 multipliers per row, i16) + ICPL (4 coefficients per row, i16, in 0.1 μas).

ICPL layout: `[dpsi_sin, dpsi_cos, deps_sin, deps_cos]`

```
for each row (reverse order):
    darg = normalize_radians(npl[0]*AL + npl[1]*ALSU + ... + npl[13]*APA)
    dpsi += icpl[0]*sin(darg) + icpl[1]*cos(darg)
    deps += icpl[2]*sin(darg) + icpl[3]*cos(darg)
```

### Unit Conversion (luni-solar + planetary)

```
O1MAS2DEG = 1.0 / 3600.0 / 10_000_000.0   // 0.1 μas → degrees
nutlo[0] = dpsi * O1MAS2DEG    // now in degrees
nutlo[1] = deps * O1MAS2DEG
```

### P03 Precession Correction (unconditional)

Applied after both series, in microarcseconds, using luni-solar OM/F/D:

```
dpsi = -8.1*sin(OM) - 0.6*sin(2*F - 2*D + 2*OM)
     + T*(47.8*sin(OM) + 3.7*sin(2*F - 2*D + 2*OM) + 0.6*sin(2*F + 2*OM) - 0.6*sin(2*OM))
deps = T*(-25.6*cos(OM) - 1.6*cos(2*F - 2*D + 2*OM))

nutlo[0] += dpsi / (3600.0 * 1_000_000.0)    // μas → degrees
nutlo[1] += deps / (3600.0 * 1_000_000.0)
```

### Final Conversion

```
nutlo[0] *= DEGTORAD    // degrees → radians
nutlo[1] *= DEGTORAD
```

## Woolard Algorithm (calc_nutation_woolard)

Labelled "incomplete implementation" in C source.

t = (J − 2415020.0) / 36525.0   (centuries from J1900)

### Five Arguments (degrees)

```
ls = 279.697 + 0.000303*t² + 360*(a - (a as i64) as f64)   where a = 100.0021358*t
ld = 270.434 - 0.001133*t² + 360*(a - (a as i64) as f64)   where a = 1336.855231*t
ms = 358.476 - 0.00015*t²  + 360*(a - (a as i64) as f64)   where a = 99.99736056000026*t
md = 296.105 + 0.009192*t² + 360*(a - (a as i64) as f64)   where a = 13255523.59*t
nm = 259.183 + 0.002078*t² - 360*(a - (a as i64) as f64)   where a = 5.372616667*t
```

Note: `(a as i64) as f64` matches C's `(long)a` truncation toward zero. Do NOT use `a.floor()` — truncation semantics differ for negative values. Results in degrees.

### Nutation Terms

All amplitudes in arcseconds. Intermediate angles:
```
tls = 2*ls * DEGTORAD
nm  = nm * DEGTORAD
tnm = 2*nm                   // note: nm already in radians
ms  = ms * DEGTORAD
tld = 2*ld * DEGTORAD
md  = md * DEGTORAD
```

13 terms for Δψ:
```
dpsi = (-17.2327 - 0.01737*t)*sin(nm) + (-1.2729 - 0.00013*t)*sin(tls)
     + 0.2088*sin(tnm) - 0.2037*sin(tld) + (0.1261 - 0.00031*t)*sin(ms)
     + 0.0675*sin(md) - (0.0497 - 0.00012*t)*sin(tls+ms)
     - 0.0342*sin(tld-nm) - 0.0261*sin(tld+md) + 0.0214*sin(tls-ms)
     - 0.0149*sin(tls-tld+md) + 0.0124*sin(tls-nm) + 0.0114*sin(tld-md)
```

9 terms for Δε:
```
deps = (9.21 + 0.00091*t)*cos(nm) + (0.5522 - 0.00029*t)*cos(tls)
     - 0.0904*cos(tnm) + 0.0884*cos(tld) + 0.0216*cos(tls+ms)
     + 0.0183*cos(tld-nm) + 0.0113*cos(tld+md) - 0.0093*cos(tls-ms)
     - 0.0066*cos(tls-nm)
```

### Final Conversion

```
nutlo[0] = dpsi / 3600.0 * DEGTORAD    // arcsec → radians
nutlo[1] = deps / 3600.0 * DEGTORAD
```

## Constants

| Name | Value | Notes |
|---|---|---|
| DPSI_IAU1980_TJD0 | 0.064284 | arcsec; Horizons constant offset |
| DEPS_IAU1980_TJD0 | 0.006151 | arcsec; Horizons constant offset |
| DPSI_DEPS_IAU1980_TJD0_HORIZONS | 2437684.5 | JD boundary for Horizons path |
| O1MAS2DEG | 1/(3600×10⁷) | 0.1 μas → degrees |
| OFFSET_JPLHORIZONS | −52.3 | mas; used in JPLHOR_APPROX |
| J1900 | 2415020.0 | Woolard epoch |

## Not Porting

**Interpolation cache** (`struct interpol`): 3-point quadratic interpolation at ±1-day intervals. Controlled by `swed.do_interpolate_nut`. Stateless Rust design means always evaluate full series. No port.

**nut_matrix / swi_nutate**: Rotation matrix construction and application. Downstream of the core nutation computation. Separate task.
