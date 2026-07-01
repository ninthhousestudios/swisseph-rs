# C Reference: Lunar Eclipses — swecl.c

Porting reference for the lunar-eclipse subsystem (`swe_lun_eclipse_how`,
`swe_lun_eclipse_when`, `swe_lun_eclipse_when_loc`). Read this instead of the C source.

All line numbers below refer to `swecl.c` unless stated otherwise.

---

## Function Map

| C function | Location | Purpose |
|---|---|---|
| `swe_lun_eclipse_how` | swecl.c:3190–3228 | Public wrapper: geocentric eclipse attributes at `tjd_ut`, plus optional topocentric azimuth/altitude of the Moon at `geopos` |
| `lun_eclipse_how` (static) | swecl.c:3237–3363 | THE CORE: selenocentric shadow-cone geometry, magnitudes, saros lookup |
| `swe_lun_eclipse_when` | swecl.c:3378–3605 | Global search: next/previous lunar eclipse from `tjd_start`, all contact times |
| `swe_lun_eclipse_when_loc` | swecl.c:3633–3728 | Local variant: next/previous lunar eclipse visible (Moon above horizon) from `geopos`, plus moonrise/moonset clipping |
| `find_maximum` (static, shared) | swecl.c:4133–4146 | Parabola-vertex fit through 3 equally-spaced samples — also used by solar-eclipse and other extremum searches |
| `find_zero` (static, shared) | swecl.c:4148–4162 | Parabola zero-crossing fit through 3 equally-spaced samples (returns both roots) — also used by solar-eclipse contact-time searches |
| `swi_dot_prod_unit` | swephlib.c:453–464 | Shared helper: cosine of angle between two vectors (normalizes both internally) |
| `square_sum`, `dot_prod` | sweph.h:308–309 (macros) | `x·x` and `x·y` for 3-vectors |

---

## 1. Constants

All defined near the top of swecl.c (lines 76–86, 114, 305) unless noted.

```c
#define AUNIT        1.49597870700e+11        /* AU in meters, DE431 (sweph.h:273) */
#define DSUN         (1392000000.0 / AUNIT)    /* Sun diameter, in AU (swecl.c:80; alt. value 1391978489.9 under #if 0) */
#define DMOON        (3476300.0 / AUNIT)       /* Moon diameter, in AU */
#define DEARTH       (6378140.0 * 2 / AUNIT)    /* Earth diameter (equatorial), in AU */
#define RSUN         (DSUN / 2)
#define RMOON        (DMOON / 2)
#define REARTH       (DEARTH / 2)

#define SAROS_CYCLE  6585.3213                 /* days; one saros ≈ 18y 11.3d */
#define NSAROS_LUNAR 180                        /* size of saros_data_lunar[] table */

#define SEI_ECL_GEOALT_MIN   (-500.0)           /* m; sweph.h:199 */
#define SEI_ECL_GEOALT_MAX   25000.0            /* m; sweph.h:198 */
```

`RSUN`/`RMOON`/`REARTH` are **radii in AU**, not km — all shadow-cone geometry in
`lun_eclipse_how` is done in AU and converted to degrees only via `asin`/`acos`.

### Eclipse type / visibility bit flags (swephexp.h:307–328)

```c
#define SE_ECL_CENTRAL              1   /* solar only, unused for lunar */
#define SE_ECL_NONCENTRAL            2   /* solar only, unused for lunar */
#define SE_ECL_TOTAL                 4
#define SE_ECL_ANNULAR               8   /* solar only, unused for lunar */
#define SE_ECL_PARTIAL              16
#define SE_ECL_ANNULAR_TOTAL        32   /* solar only, unused for lunar */
#define SE_ECL_PENUMBRAL            64
#define SE_ECL_ALLTYPES_LUNAR       (SE_ECL_TOTAL|SE_ECL_PARTIAL|SE_ECL_PENUMBRAL)
#define SE_ECL_ALLTYPES_SOLAR       (SE_ECL_CENTRAL|SE_ECL_NONCENTRAL|SE_ECL_TOTAL|SE_ECL_ANNULAR|SE_ECL_PARTIAL|SE_ECL_ANNULAR_TOTAL)
#define SE_ECL_VISIBLE               128
#define SE_ECL_MAX_VISIBLE           256
#define SE_ECL_1ST_VISIBLE           512   /* == SE_ECL_PARTBEG_VISIBLE, same bit */
#define SE_ECL_PARTBEG_VISIBLE       512   /* begin of partial phase */
#define SE_ECL_TOTBEG_VISIBLE       1024   /* begin of total phase */
#define SE_ECL_TOTEND_VISIBLE       2048   /* end of total phase */
#define SE_ECL_PARTEND_VISIBLE      4096   /* end of partial phase */
#define SE_ECL_PENUMBBEG_VISIBLE    8192   /* begin of penumbral phase */
#define SE_ECL_PENUMBEND_VISIBLE   16384   /* end of penumbral phase */
```

`retc`/`retflag` for lunar functions is always a subset of
`SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_PENUMBRAL` (`lun_eclipse_how`'s return value), OR'd
with the `*_VISIBLE` bits for `swe_lun_eclipse_when_loc`. `retflag == 0` means "no eclipse".

---

## 2. `lun_eclipse_how` (static core) — swecl.c:3237–3363

```c
static int32 lun_eclipse_how(double tjd_ut, int32 ifl, double *attr, double *dcore, char *serr)
```

Computes the selenocentric shadow-cone geometry and the eclipse magnitudes for one instant
`tjd_ut`. Called directly (not via the public wrapper) by both `swe_lun_eclipse_when` (during
the parabolic contact-time refinement) and `swe_lun_eclipse_how`.

### 2.1 Setup

- Zero `dcore[0..9]` and `attr[0..19]` (line 3255–3258).
- `iflag = SEFLG_SPEED | SEFLG_EQUATORIAL | ifl | SEFLG_XYZ` — equatorial cartesian, with
  speed (speed is requested but not actually used further down; it is a leftover/harmless
  over-request).
- `deltat = swe_deltat_ex(tjd_ut, ifl, serr)`; `tjd = tjd_ut + deltat` — **all subsequent
  positions are computed at TT**, not UT.
- `swe_calc(tjd, SE_MOON, iflag, rm, serr)` → geocentric Moon, cartesian, equatorial, AU.
- `dm = sqrt(square_sum(rm))` — geocentric Moon distance.
- `swe_calc(tjd, SE_SUN, iflag, rs, serr)` → geocentric Sun, cartesian, equatorial, AU.
- `ds = sqrt(square_sum(rs))` — geocentric Sun distance.
- Unit vectors `x1 = rs/ds`, `x2 = rm/dm`; `dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG`
  — geocentric elongation Sun↔Moon (used later only to derive `attr[7]`, distance from
  opposition, as `180 − |dctr|`).

### 2.2 Change of origin: selenocentric frame

```c
rs[i] -= rm[i];   /* selenocentric sun    (sun position minus moon position) */
rm[i] = -rm[i];   /* selenocentric earth  (negate geocentric moon vector) */
```
(lines 3281–3285). After this, `rs` = Sun as seen from the Moon, `rm` = Earth as seen from
the Moon — this is the mirror image of the solar-eclipse geometry (which is built from the
Earth's perspective looking at the Moon's shadow on Earth); here the "shadow" being measured
is the **Earth's** shadow, and the "screen" is the Moon.

```c
e[i] = rm[i] - rs[i];         /* sun→earth vector, in the selenocentric frame */
dsm = sqrt(square_sum(e));     /* distance sun–earth (via selenocentric detour;
                                   numerically ≈ true sun-earth distance)      */
e[i] /= dsm;                    /* unit vector along sun→earth axis            */
```
(lines 3287–3293). `e` is the axis of the earth's shadow cone (antisolar direction, as seen
from the Sun→Earth line).

### 2.3 Umbra/penumbra half-angles

```c
f1 = (RSUN - REARTH) / dsm;   cosf1 = sqrt(1 - f1*f1);   /* umbra cone half-angle sine/cosine */
f2 = (RSUN + REARTH) / dsm;   cosf2 = sqrt(1 - f2*f2);   /* penumbra cone half-angle sine/cosine */
```
(lines 3294–3297). `f1 = sin(umbra half-angle)`, `f2 = sin(penumbra half-angle)` — standard
shadow-cone trig using similar triangles Sun–Earth with radii `RSUN`, `REARTH`.

### 2.4 Position of the Moon relative to the shadow axis

```c
s0 = -dot_prod(rm, e);                    /* distance of earth from the "fundamental plane"
                                              (plane through the moon, perpendicular to e) */
r0 = sqrt(dm * dm - s0 * s0);              /* distance of the shadow axis from the selenocenter
                                              (dm here is |rm|, the selenocentric earth distance,
                                              reused as the hypotenuse) */
```
(lines 3299–3301). This is a right-triangle decomposition: `dm` (= distance Earth↔Moon, since
`rm` after negation still has the same magnitude) is the hypotenuse, `s0` is the leg along the
shadow axis `e`, `r0` is the leg perpendicular to it — i.e. how far off-axis the Moon is from
the shadow center line, projected onto the fundamental plane at the Moon's distance.

**FP hazard**: `dm` here is the *Moon's geocentric distance* computed at §2.1 (`dm = sqrt(square_sum(rm))` before negation — magnitude is unchanged by negation) and is reused directly; it is not recomputed after the frame change. The Rust port must reuse the exact same scalar, not recompute `sqrt(square_sum(new rm))` (same value mathematically, but reuse avoids an extra sqrt and any FP divergence).

### 2.5 Shadow diameters on the fundamental plane, with atmospheric enlargement

```c
/* diameter of core (umbra) shadow on fundamental plane */
/* one 50th is added for effect of atmosphere, AA98, L4 */
d0 = fabs(s0 / dsm * (DSUN - DEARTH) - DEARTH) * (1 + 1.0 / 50.0) / cosf1;
/* diameter of half-shadow (penumbra) on fundamental plane */
D0 = (s0 / dsm * (DSUN + DEARTH) + DEARTH) * (1 + 1.0 / 50.0) / cosf2;
d0 /= cosf1;
D0 /= cosf2;
/* for better agreement with NASA: */
d0 *= 0.99405;
D0 *= 0.98813;
```
(lines 3302–3311). Read this exactly as ordered — it is **not** algebraically simplified in
the C, and the Rust port must preserve the same expression/operation order for bit-fidelity:

1. `d0` starts as the raw cone-similar-triangles umbra diameter at the Moon's distance:
   `|s0/dsm·(DSUN−DEARTH) − DEARTH|`.
2. Multiplied by the atmospheric-enlargement factor `(1 + 1/50)` = `1.02` — this accounts for
   Earth's atmosphere effectively enlarging the shadow (cited source: *Astronomical Almanac*
   1998, section L4).
3. Divided by `cosf1` **once** (`/ cosf1` inline in the same statement).
4. Divided by `cosf1` **again**, a second time, on the next line (`d0 /= cosf1;`) — i.e. the
   umbra diameter ends up divided by `cosf1` **twice** (`cosf1²` total). This is present in the
   C as written; do not "fix" it as a redundancy — replicate exactly, since golden-test
   fidelity is against this exact code path.
5. Symmetric treatment for `D0` (penumbra) with `cosf2` divided twice and the `(DSUN+DEARTH)`
   sign (penumbra cone diverges, so it's a sum, not a difference).
6. Finally both are scaled by empirical fudge factors for agreement with NASA's published
   eclipse circumstances: `d0 *= 0.99405`, `D0 *= 0.98813`. These are unexplained-in-source
   empirical corrections — carry them as opaque literal constants.

`dcore[]` output (only meaningful/populated fields; declared `double dcore[10]`, all others
left 0):

| Index | Meaning |
|---|---|
| `dcore[0]` | `r0` — distance of shadow axis from selenocenter (AU) |
| `dcore[1]` | `d0` — diameter of umbra (core shadow) on the fundamental plane (AU), atmosphere+NASA-corrected |
| `dcore[2]` | `D0` — diameter of penumbra (half-shadow) on the fundamental plane (AU), atmosphere+NASA-corrected |
| `dcore[3]` | `cosf1` — cosine of umbra cone half-angle |
| `dcore[4]` | `cosf2` — cosine of penumbra cone half-angle |
| `dcore[5..9]` | unused, always 0 |

`dcore` is consumed directly (not through `attr`) by `swe_lun_eclipse_when`'s contact-time
refinement (§4.5) — it needs `cosf1`/`cosf2` there too, which `attr` does not carry.

### 2.6 Phase / umbral magnitude

```c
retc = 0;
if (d0/2 >= r0 + rmoon/cosf1) {
    retc = SE_ECL_TOTAL;
    attr[0] = (d0/2 - r0 + rmoon) / dmoon;
} else if (d0/2 >= r0 - rmoon/cosf1) {
    retc = SE_ECL_PARTIAL;
    attr[0] = (d0/2 - r0 + rmoon) / dmoon;
} else if (D0/2 >= r0 - rmoon/cosf2) {
    retc = SE_ECL_PENUMBRAL;
    attr[0] = 0;
} else {
    /* no lunar eclipse at tjd = ... */    /* retc stays 0 */
}
attr[8] = attr[0];
```
(lines 3320–3334; `rmoon = RMOON`, `dmoon = 2*rmoon` set at function entry, line 3252–3253).

- **Total**: the Moon (radius `rmoon`) is entirely inside the umbra circle (radius `d0/2`,
  expanded by `rmoon/cosf1` to allow for the umbra edge being oblique) even at its point of
  furthest offset `r0` from the axis. Umbral magnitude = fraction of the Moon's diameter
  covered by the umbra, `(d0/2 − r0 + rmoon) / dmoon`, and can exceed 1.0 for a deep total
  eclipse.
- **Partial**: the umbra circle overlaps the Moon's disc but doesn't fully cover it — same
  magnitude formula, now `∈ (0, 1)`.
- **Penumbral**: no umbral contact at all, but the penumbra circle (radius `D0/2`) overlaps
  the Moon — `attr[0]` forced to 0 (no umbral magnitude for a penumbral-only eclipse).
- **None**: none of the three circles reach the Moon; `retc` stays 0 and `serr` is set to
  `"no lunar eclipse at tjd = %f"` (informational, not necessarily an error — caller checks
  `retc`, not `serr`, for "no eclipse").
- `attr[8]` mirrors `attr[0]` (umbral magnitude), documented in the public API as a
  convenience duplicate.

### 2.7 Penumbral magnitude (always computed, independent of `retc`)

```c
attr[1] = (D0/2 - r0 + rmoon) / dmoon;
```
(line 3338). Computed unconditionally — even for a total/partial eclipse, `attr[1]` holds the
penumbral magnitude (which is always ≥ umbral magnitude in a real eclipse) — this is *not*
gated behind `retc == SE_ECL_PENUMBRAL`.

### 2.8 Distance from opposition

```c
if (retc != 0)
    attr[7] = 180 - fabs(dctr);
```
(lines 3339–3340). Only set when an eclipse is occurring; `dctr` was the geocentric Sun–Moon
elongation computed in §2.1 — `180 − |dctr|` is how far short of exact opposition (full moon)
the current instant is, in degrees.

### 2.9 Saros series lookup

```c
for (i = 0; i < NSAROS_LUNAR; i++) {
    d = (tjd_ut - saros_data_lunar[i].tstart) / SAROS_CYCLE;
    if (d < 0 && d * SAROS_CYCLE > -2) d = 0.0000001;
    if (d < 0) continue;
    j = (int) d;
    if ((d - j) * SAROS_CYCLE < 2) {           /* within 2 days after a whole-cycle multiple */
        attr[9] = (double) saros_data_lunar[i].series_no;
        attr[10] = (double) j + 1;
        break;
    }
    k = j + 1;
    if ((k - d) * SAROS_CYCLE < 2) {            /* within 2 days before the next multiple */
        attr[9] = (double) saros_data_lunar[i].series_no;
        attr[10] = (double) k + 1;
        break;
    }
}
if (i == NSAROS_LUNAR)
    attr[9] = attr[10] = -99999999;
```
(lines 3342–3361). `saros_data_lunar[NSAROS_LUNAR]` (`NSAROS_LUNAR = 180`, table at
swecl.c:306 onward — parallel to, but a distinct table from, `saros_data_solar[NSAROS_SOLAR]`
at line 116) is `struct saros_data {int series_no; double tstart;}`, one entry per saros
series giving its starting JD. For each series, `d` = number of saros cycles elapsed since
`tstart`; the loop finds the series/member whose predicted eclipse date (an integer multiple
of `SAROS_CYCLE` days after `tstart`) lands within ±2 days of `tjd_ut`, and reports:
`attr[9]` = saros series number, `attr[10]` = member number within that series (1-based,
`j+1` or `k+1`). If no series matches (shouldn't normally happen for valid eclipse dates
within the table's coverage), both are set to the sentinel `-99999999`.

### 2.10 Return value

`lun_eclipse_how` returns `retc` (0, `SE_ECL_TOTAL`, `SE_ECL_PARTIAL`, or `SE_ECL_PENUMBRAL`)
— exactly one of these three bits, never combined, since the `if`/`else if` chain in §2.6 is
mutually exclusive. Returns `ERR` only if the underlying `swe_calc` for Sun or Moon fails.

---

## 3. `swe_lun_eclipse_how` (public wrapper) — swecl.c:3190–3228

```c
int32 CALL_CONV swe_lun_eclipse_how(double tjd_ut, int32 ifl, double *geopos, double *attr, char *serr)
```

1. `geopos` may be `NULL` (geocentric-only query). If non-`NULL`, validate
   `geopos[2]` (altitude, meters) is within `[SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX]` =
   `[-500, 25000]`; else `serr` = `"location for eclipses must be between %.0f and %.0f m above sea"`, return `ERR`.
2. Strip `SEFLG_TOPOCTR` and `SEFLG_JPLHOR`/`SEFLG_JPLHOR_APPROX` from `ifl` (lines 3208–3209)
   — the geocentric core calculation must never run topocentric or JPL-Horizons-matched;
   those refinements are applied only to the separate azimuth/altitude call below.
3. `swi_set_tid_acc(tjd_ut, ifl, 0, serr)` — **STATELESS PORT NOTE:** this sets the C library's
   global tidal-acceleration/ΔT model state (used by `swe_deltat_ex` internally when no
   explicit model is requested). The stateless Rust `Ephemeris` must instead thread the
   ΔT/tidal-acceleration configuration explicitly through `&self` into the `swe_deltat_ex`
   equivalent rather than mutating shared state.
4. `retc = lun_eclipse_how(tjd_ut, ifl, attr, dcore, serr)` — the core geometry (§2). `dcore`
   is a local scratch array here, discarded after the call (not returned to the caller of the
   public API — only `swe_lun_eclipse_when`'s internal refinement loop consumes `dcore`
   directly).
5. If `geopos == NULL`, return `retc` immediately — no topocentric azimuth/altitude.
6. Otherwise, compute topocentric azimuth/altitude of the Moon at `geopos` and `tjd_ut`:
   - `swe_set_topo(geopos[0], geopos[1], geopos[2])` — **STATELESS PORT NOTE:** global
     topocentric-observer state; the Rust port must pass geographic position as an explicit
     parameter into the topocentric position/az-alt calculation instead.
   - `swe_calc_ut(tjd_ut, SE_MOON, ifl | SEFLG_TOPOCTR | SEFLG_EQUATORIAL, lm, serr)` →
     topocentric equatorial Moon position `lm` (apparent, at UT).
   - `swe_azalt(tjd_ut, SE_EQU2HOR, geopos, 0, 10, lm, xaz)` → `xaz[0]` = azimuth,
     `xaz[1]` = true altitude, `xaz[2]` = apparent altitude (the `0, 10` args are
     atpress=0/attemp=10 — i.e. "estimate pressure from `geopos[2]`, use 10°C" per
     `swe_azalt`'s convention).
   - `attr[4] = xaz[0]` (azimuth), `attr[5] = xaz[1]` (true altitude), `attr[6] = xaz[2]`
     (apparent altitude).
   - If `xaz[2] <= 0` (Moon below horizon by apparent altitude), force `retc = 0` — i.e. the
     public wrapper reports "no eclipse" when the Moon is not up, even though the geocentric
     geometry may say an eclipse is in progress. (Note: `attr[0]`/`attr[1]`/magnitudes are
     still left populated from the geocentric calculation; only the return code is zeroed.)

### `attr[]` full index (swecl.c:3172–3189 doc comment; shared by all three public functions)

| Index | Meaning | Populated by |
|---|---|---|
| `attr[0]` | umbral magnitude at `tjd` | `lun_eclipse_how` §2.6 |
| `attr[1]` | penumbral magnitude | `lun_eclipse_how` §2.7 |
| `attr[2]` | *(unused for lunar; reserved — solar eclipse uses this slot for fraction of solar disc covered)* | — |
| `attr[3]` | *(unused for lunar; solar eclipse: ratio of diameters)* | — |
| `attr[4]` | azimuth of Moon at `tjd` | `swe_lun_eclipse_how` §3 step 6 (only if `geopos != NULL`) |
| `attr[5]` | true altitude of Moon above horizon at `tjd` | ditto |
| `attr[6]` | apparent altitude of Moon above horizon at `tjd` | ditto |
| `attr[7]` | distance of Moon from opposition, in degrees | `lun_eclipse_how` §2.8 |
| `attr[8]` | umbral magnitude at `tjd` (duplicate of `attr[0]`) | `lun_eclipse_how` §2.6 |
| `attr[9]` | saros series number | `lun_eclipse_how` §2.9 |
| `attr[10]` | saros series member number | `lun_eclipse_how` §2.9 |
| `attr[11..19]` | unused, left 0 | — |

Caller must declare `attr[20]` minimum (per C doc comment); `attr[2]`/`attr[3]` are vestigial
slots kept for index-parity with the solar-eclipse `attr[]` layout (`swe_sol_eclipse_how`),
not used by the lunar path.

---

## 4. `swe_lun_eclipse_when` — swecl.c:3378–3605

```c
int32 CALL_CONV swe_lun_eclipse_when(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr)
```

Global search for the next (or previous) lunar eclipse of a requested type, starting from
`tjd_start`. No geographic position — purely geocentric.

### 4.1 Setup and `ifltype` normalization

- `ifl &= SEFLG_EPHMASK` (`SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH`) — only the ephemeris-source
  bits of `ifl` survive; sidereal/topocentric/etc. flags are stripped for this search (the
  search only needs ephemeris-source selection).
- `swi_set_tid_acc(tjd_start, ifl, 0, serr)` — **STATELESS PORT NOTE**, same as §3 step 3.
- `iflag = SEFLG_EQUATORIAL | ifl`; `iflagcart = iflag | SEFLG_XYZ`.
- `ifltype`: strip `SE_ECL_CENTRAL|SE_ECL_NONCENTRAL` (solar-only bits, harmless if a caller
  passes a mask copied from solar-eclipse code). If `SE_ECL_ANNULAR|SE_ECL_ANNULAR_TOTAL` bits
  are set, strip them; if that leaves `ifltype == 0`, error out ("annular lunar eclipses don't
  exist") — **this guards against an infinite search loop** for an eclipse type that can never
  occur for the Moon. If `ifltype == 0` after all stripping (i.e. caller passed 0 or only
  solar-only bits), default to `SE_ECL_TOTAL | SE_ECL_PENUMBRAL | SE_ECL_PARTIAL` (search for
  any lunar eclipse type).
- `direction = backward ? -1 : 1`.

### 4.2 Full-moon (synodic-month) stepping via Meeus's lunation number `K`

```c
K = (int)((tjd_start - J2000) / 365.2425 * 12.3685);
K -= direction;
next_try:
```
(lines 3417–3419). `K` is Meeus's lunation number (fractional lunations since the 2000-01-06
"K=0" new moon epoch used by *Astronomical Algorithms* ch. 49); `12.3685` ≈ synodic months per
Julian year. Stepping `K += direction` (line 3433 etc., on every rejection branch) advances or
retreats by exactly one synodic month (~29.53 days) per iteration — full moons only, since
lunar eclipses can only occur at full moon.

### 4.3 Eclipse-possibility pre-filter (latitude argument `F`)

```c
kk = K + 0.5;                    /* full-moon phase offset from Meeus's new-moon-based K */
T = kk / 1236.85; T2=T*T; T3=T2*T; T4=T3*T;
Ff = F = swe_degnorm(160.7108 + 390.67050274*kk - 0.0016341*T2 - 0.00000227*T3 + 0.000000011*T4);
if (Ff > 180) Ff -= 180;
if (Ff > 21 && Ff < 159) {          /* no eclipse possible */
    K += direction;
    goto next_try;
}
```
(lines 3423–3435). `F` is the Moon's argument of latitude (mean elongation from its ascending
node) per Meeus's lunar theory. An eclipse is geometrically impossible unless the full moon
occurs close enough to a node — the coefficients are Meeus, *Astronomical Algorithms*
(German ed.), ch. 54, lunar-eclipse conditions table. `Ff ∈ (21°, 159°)` (i.e. too far from
0°/180°, the node crossings) is rejected outright without any further ephemeris computation —
this is a cheap analytic filter to skip non-eclipse full moons quickly.

### 4.4 Approximate time of maximum eclipse (Meeus periodic series)

Lines 3436–3475: `M` (Sun's mean anomaly), `Mm` (Moon's mean anomaly), `Om` (longitude of
ascending node), `E` (eccentricity correction factor for Earth's orbit), `A1` (Venus
perturbation argument) are all computed via Meeus's polynomial series in `kk`/`T`/`T2`/`T3`
(same structural pattern as the solar-eclipse search — if a `c-ref-eclipse-solar.md` exists,
these are the same formulas/coefficients as its analogous new-moon block; not re-derived
here). `tjd` (approximate JD of geocentric syzygy) accumulates a base polynomial plus a
periodic correction sum of ~17 sine terms in `M`, `Mm`, `F1` (`= F - 0.02665·sin(Om)`, in
degrees pre-conversion), `A1`, `Om`, each with its own amplitude coefficient (0.4075 down to
0.0002) and multiplier of the anomaly/argument — this is Meeus ch. 49's full/new moon time
series, evaluated for the full-moon case (all coefficients here are the ones tabulated for
"Full Moon", distinct from the new-moon coefficient set used by the solar search). All
angles `M, Mm, F, Om` are converted to radians in-place (line 3454–3457) before the
trig evaluation; `A1` converted separately (line 3459). This produces only an **approximate**
JD of geocentric opposition — refined next.

### 4.5 Precise refinement to time of maximum eclipse

```c
dtstart = 0.1;
if (tjd < 2100000 || tjd > 2500000)   /* was tjd < 2000000 until 26-aug-22 */
    dtstart = 5;
dtdiv = 4;
for (j = 0, dt = dtstart; dt > 0.001; j++, dt /= dtdiv) {
    for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {
        /* compute selenocentric-angle proxy dc[i] at t-dt, t, t+dt */
        swe_calc(t, SE_SUN, iflagcart, xs, ...);
        swe_calc(t, SE_MOON, iflagcart, xm, ...);
        xs[m] -= xm[m];  xm[m] = -xm[m];         /* selenocentric sun / earth, as in §2.2 */
        ds = |xs|; dm = |xm|;
        xa = xs/ds; xb = xm/dm;
        dc[i] = acos(swi_dot_prod_unit(xa, xb)) * RADTODEG;
        rearth = asin(REARTH/dm) * RADTODEG;
        rsun = asin(RSUN/ds) * RADTODEG;
        dc[i] -= (rearth + rsun);
    }
    find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
    tjd += dtint + dt;
}
```
(lines 3484–3513). `dc[i]` is a proxy quantity: the selenocentric Sun–Earth angular separation
minus the sum of their angular radii as seen from the Moon — i.e. roughly zero/negative when
the Earth's shadow and the Sun are in maximal alignment (this is the same
"angle minus radii" trick used elsewhere for extremum search, not the exact `r0`/`d0`
geometry — it's cheaper and good enough to locate the parabola vertex). This is minimized
(most negative) at exact opposition/maximum eclipse.

- **Step-size schedule**: starts at `dtstart` (0.1 day normally; **5 days** if `tjd` is
  outside `[2100000, 2500000]` JD — i.e. far from the well-conditioned modern epoch range;
  comment notes this threshold was widened from `2000000` on 26-Aug-2022), divides by
  `dtdiv = 4` each outer iteration, and stops once `dt <= 0.001` days (~86 sec) — i.e. an
  iterative bracket-shrink, not a fixed iteration count; number of iterations depends on
  `dtstart` (with `dtstart=0.1`: `0.1, 0.025, 0.00625, 0.0015625` → 4 iterations before
  `dt ≤ 0.001`; with `dtstart=5`: 6 iterations).
- Each outer iteration: sample `dc` at `t-dt, t, t+dt`, fit a parabola via `find_maximum`
  (§ shared helper below) to get the offset `dtint` of the true vertex from the center sample,
  then `tjd += dtint + dt` — note `+ dt`, **not** `+ dtint` alone: this re-centers `tjd` at the
  *last* sample point (`t` after the loop equals `tjd_prev - dt + 2*dt = tjd_prev + dt`) plus
  the parabola-vertex offset from that point, mirroring `find_maximum`'s convention that
  `dxret` is measured relative to the **third** sample `y2` (see §4.8 below), not the center
  sample `y11`.
- After the loop, three successive ΔT corrections (lines 3514–3516):
  ```c
  tjd2 = tjd - swe_deltat_ex(tjd, ifl, serr);
  tjd2 = tjd - swe_deltat_ex(tjd2, ifl, serr);
  tjd = tjd - swe_deltat_ex(tjd2, ifl, serr);
  ```
  converts the TT-based `tjd` found above back to UT via 3 fixed-point iterations of
  `UT = TT − ΔT(UT)` (ΔT depends on UT, hence the iteration; 3 rounds is empirically enough
  for convergence, not a documented tolerance-based loop).

### 4.6 Confirm eclipse and reject wrong types

```c
if ((retflag = swe_lun_eclipse_how(tjd, ifl, NULL, attr, serr)) == ERR) return retflag;
if (retflag == 0) { K += direction; goto next_try; }
tret[0] = tjd;
if ((backward && tret[0] >= tjd_start - 0.0001) || (!backward && tret[0] <= tjd_start + 0.0001)) {
    K += direction; goto next_try;
}
if (!(ifltype & SE_ECL_PENUMBRAL) && (retflag & SE_ECL_PENUMBRAL)) { K += direction; goto next_try; }
if (!(ifltype & SE_ECL_PARTIAL)   && (retflag & SE_ECL_PARTIAL))   { K += direction; goto next_try; }
if (!(ifltype & SE_ECL_TOTAL)     && (retflag & SE_ECL_TOTAL))     { K += direction; goto next_try; }
```
(lines 3517–3546). Calls the **public** `swe_lun_eclipse_how` (not the static core directly)
with `geopos = NULL` — i.e. purely geocentric confirmation, at the refined `tjd` (UT). If
`retflag == 0`, the candidate full moon isn't actually an eclipse (the §4.3 filter is
necessary-but-not-sufficient) — step to the next lunation and retry. The `tjd_start` boundary
check (with 0.0001-day ≈ 8.6-second tolerance) prevents returning an eclipse at/before the
search start when searching forward (or at/after start when searching backward). Then reject
if the found eclipse's type isn't in the requested `ifltype` mask — note these three checks
are independent, not mutually exclusive rejections (a found eclipse only has exactly one type
bit set per §2.10, so exactly one of the three `if`s can trigger a rejection for any given
found eclipse).

### 4.7 Contact-time computation via `dcore`-based zero search

```c
if (retflag & SE_ECL_PENUMBRAL) o = 0;
else if (retflag & SE_ECL_PARTIAL) o = 1;
else o = 2;
dta = twohr;   /* 2/24 day */
dtb = tenmin;  /* 10/24/60 day, overwritten below */
for (n = 0; n <= o; n++) {
    if (n == 0) { i1 = 6; i2 = 7; }        /* penumbral begin/end */
    else if (n == 1) { i1 = 2; i2 = 3; }    /* partial begin/end */
    else if (n == 2) { i1 = 4; i2 = 5; }    /* total begin/end */
    ...
}
```
(lines 3552–3567). `o` controls how many contact-pairs to compute, based on the eclipse's
actual type: a penumbral-only eclipse (`o=0`) only gets penumbral contacts; a partial eclipse
(`o=1`) gets both penumbral (`n=0`) and partial (`n=1`) contacts; a total eclipse (`o=2`) gets
all three (`n=0,1,2`): penumbral, partial, **and** total contacts. This matches physical
containment: a total eclipse is necessarily also partial and penumbral.

For each `n`, two stages:

**Stage A — coarse bracket** (lines 3568–3582, wrapped in `#if 1`/`#else` with the `#else`
branch — a naive `tjd ± dtb` fallback — dead code, never compiled):
```c
for (i = 0, t = tjd - dta; i <= 2; i += 1, t += dta) {
    lun_eclipse_how(t, ifl, attr, dcore, serr);
    if (n == 0)      dc[i] = dcore[2]/2 + RMOON/dcore[4] - dcore[0];   /* penumbra edge minus axis distance */
    else if (n == 1) dc[i] = dcore[1]/2 + RMOON/dcore[3] - dcore[0];   /* umbra edge (+moon) minus axis distance */
    else if (n == 2) dc[i] = dcore[1]/2 - RMOON/dcore[3] - dcore[0];   /* umbra edge (−moon) minus axis distance */
}
find_zero(dc[0], dc[1], dc[2], dta, &dt1, &dt2);
dtb = (dt1 + dta) / 2;
tret[i1] = tjd + dt1 + dta;
tret[i2] = tjd + dt2 + dta;
```
Samples `dcore` at `tjd - dta, tjd, tjd + dta` (`dta = twohr = 2/24` day, ±2 hours around the
already-refined maximum time `tjd`), builds the "distance from shadow-circle edge" proxy
`dc[i]` — this is exactly the boundary-crossing condition of §2.6's magnitude formulas restated
as "≥ 0 inside, < 0 outside" (sign convention: `dcore[0]` = `r0`, the shadow-axis distance,
subtracted from the relevant circle radius `± RMOON/cosf`):
  - `n=0` (penumbral): circle radius `D0/2` (`dcore[2]/2`), enlarged by `RMOON/cosf2`
    (`dcore[4]`) — first/last penumbral contact (limb of Moon touches edge of penumbra).
  - `n=1` (partial/umbral first contact): circle radius `d0/2` (`dcore[1]/2`), enlarged by
    `RMOON/cosf1` (`dcore[3]`) — limb of Moon touches edge of umbra (partial begin/end).
  - `n=2` (total/umbral full containment): circle radius `d0/2`, **shrunk** by `RMOON/cosf1`
    (subtracted, not added) — far limb of Moon exits/enters the umbra circle entirely
    (totality begin/end).
- `find_zero` fits a parabola through the 3 samples and returns both roots `dt1, dt2` (offsets
  from the **last** sample, `t = tjd + dta`, per the shared helper's convention — see §4.8).
- `tret[i1] = tjd + dt1 + dta` (first/earlier root → begin), `tret[i2] = tjd + dt2 + dta`
  (second/later root → end). `dtb` is set from `dt1` here (`(dt1+dta)/2`) but this value of
  `dtb` is then immediately overwritten inside the refinement loop below on each `m` — it's
  effectively unused output of this branch except within the same statement.

**Stage B — Newton-style bisection refinement** (lines 3587–3602), 3 rounds (`m = 0..2`),
each halving `dt`:
```c
for (m = 0, dt = dtb / 2; m < 3; m++, dt /= 2) {
    for (j = i1; j <= i2; j += (i2 - i1)) {     /* j = i1, then j = i2 (two-point loop) */
        for (i = 0, t = tret[j] - dt; i < 2; i++, t += dt) {
            lun_eclipse_how(t, ifl, attr, dcore, serr);
            dc[i] = <same formula as stage A, per current n>;
        }
        dt1 = dc[1] / ((dc[1] - dc[0]) / dt);    /* secant/linear-interpolation step */
        tret[j] -= dt1;
    }
}
```
For each of the two contact times (`tret[i1]`, `tret[i2]`) independently, samples the same
`dc` proxy at `t - dt` and `t` (only **2** points here, not 3 — this is a secant-method linear
root refinement, not a parabola fit), computes a linear-interpolation correction
`dt1 = dc[1] / ((dc[1]-dc[0])/dt)` (finite-difference slope estimate, then Newton step assuming
the root is at `dc=0`), and subtracts it from `tret[j]`. `dt` starts at `dtb/2` and halves
each of 3 rounds — this progressively tightens the bracket for a fixed 3 iterations (not a
convergence-tolerance loop). **This is a distinct numerical method from `find_zero`'s parabola
fit** — stage A brackets coarsely with a parabola, stage B refines with 3 rounds of 2-point
secant/Newton steps.

### 4.8 Shared parabola helpers — `find_maximum` / `find_zero` (swecl.c:4133–4162)

Both fit a parabola `y = a·x² + b·x + c` through 3 samples `y00, y11, y2` taken at
equally-spaced abscissas `x = 0, 1, 2` (spacing `dx` in real units), using:
```c
c = y11;
b = (y2 - y00) / 2.0;
a = (y2 + y00) / 2.0 - c;
```
(standard 3-point Lagrange/finite-difference parabola coefficients, `x=1` as the reference
point).

- **`find_maximum`**: vertex at `x = -b/(2a)`; `*dxret = (x - 1) * dx` — **offset measured
  relative to the middle sample (`x=1`)**, i.e. `dxret` is added to the time of the *middle*
  sample to get the true extremum time. `*yret = (4ac - b²)/(4a)` — the extremum value itself
  (only computed if `yret != NULL`).
- **`find_zero`**: solves `a·x² + b·x + c = 0` via the quadratic formula; returns `ERR` if the
  discriminant `b² - 4ac < 0` (no real root — parabola doesn't cross zero over the sampled
  range). `*dxret = (x1 - 1) * dx`, `*dxret2 = (x2 - 1) * dx` — **again offsets relative to the
  middle sample** `x=1`, not the third sample. (Re-examine §4.5/§4.7 call sites: `tjd += dtint
  + dt` and `tret[i1] = tjd + dt1 + dta` both add back the *bracket half-width* (`dt`/`dta`) on
  top of the `dxret` offset — this is because the loop that produced the samples starts
  iterating from `t = tjd - dt` (or `- dta`), so by the time the samples are gathered, the
  *local* `tjd`/`t` variable used as the mid-sample reference has NOT been advanced to the
  middle sample's time; the caller's `tjd`/`tjd - dta` is the abscissa of `y00` (first
  sample, `x=0`), and `dxret` is relative to `x=1` (the middle sample) — so the correct
  absolute reconstruction is `first_sample_time + dx + dxret`, i.e. `tjd + dt + dtint`,
  matching the code exactly. Do not re-derive this as "relative to first sample" — it is
  relative to the **middle** sample, and the `+ dt`/`+ dta` term is what converts a
  middle-sample-relative offset back to an absolute time given that the loop variable `tjd`/`t`
  still points at the *first* sample.)

Both are also used by the solar-eclipse search code (same functions, shared, not
lunar-specific) — if porting solar eclipses too, reuse a single Rust `find_maximum`/`find_zero`
(or equivalent parabola-fit) helper rather than duplicating.

### `tret[]` full index (swecl.c:3369–3376 / 3611–3618 doc comments)

| Index | Meaning | Set by |
|---|---|---|
| `tret[0]` | time of maximum eclipse | §4.5/§4.6 (overwritten in `when_loc`, see §5) |
| `tret[1]` | *(unused for lunar eclipses — reserved for index-parity with solar `tret[]`, which uses it for "time of first contact" in some solar variants)* | — |
| `tret[2]` | time of partial phase begin (umbra first touches Moon's limb) | §4.7, `n=1`, `i1=2` |
| `tret[3]` | time of partial phase end | §4.7, `n=1`, `i2=3` |
| `tret[4]` | time of totality begin | §4.7, `n=2`, `i1=4` |
| `tret[5]` | time of totality end | §4.7, `n=2`, `i2=5` |
| `tret[6]` | time of penumbral phase begin | §4.7, `n=0`, `i1=6` |
| `tret[7]` | time of penumbral phase end | §4.7, `n=0`, `i2=7` |
| `tret[8]` | time of moonrise, if it occurs during the eclipse (`swe_lun_eclipse_when_loc` only) | §5 |
| `tret[9]` | time of moonset, if it occurs during the eclipse (`swe_lun_eclipse_when_loc` only) | §5 |

Any `tret[i]` not applicable to the found eclipse's type (e.g. `tret[4]`/`tret[5]` for a
partial-only eclipse) is left at 0 (all zeroed at the top of `next_try:`, line 3421–3422).

Return value: same `retflag` (`SE_ECL_TOTAL`/`SE_ECL_PARTIAL`/`SE_ECL_PENUMBRAL`, single bit)
as `swe_lun_eclipse_how`, or `ERR`.

---

## 5. `swe_lun_eclipse_when_loc` — swecl.c:3633–3728

```c
int32 CALL_CONV swe_lun_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr)
```

Finds the next/previous lunar eclipse (any type) that is at least partly **visible** from
`geopos` (Moon above the horizon during some phase of the eclipse), clipping contact times to
moonrise/moonset as needed.

1. Validate `geopos[2]` range exactly as in §3 step 1 (same error message).
2. Strip `SEFLG_JPLHOR|SEFLG_JPLHOR_APPROX` from `ifl`.
3. **`next_lun_ecl:`** loop label. Call `swe_lun_eclipse_when(tjd_start, ifl, 0, tret, backward, serr)` — `ifltype = 0` → any type (§4.1 default). This gives the full geocentric `tret[]` set (§4's table) with no location filtering yet.
4. **Visibility scan** (lines 3652–3671): for `i` from 7 down to 0 (skip `i==1`, unused slot;
   skip any `tret[i] == 0`, i.e. not-applicable contact), call `swe_lun_eclipse_how(tret[i],
   ifl, geopos, attr, serr)` (the **public** wrapper — this fills `attr[4..6]` az/alt at that
   specific contact time) and check `attr[6] > 0` (apparent altitude of Moon > 0, i.e. above
   horizon accounting for refraction). If visible, OR in `SE_ECL_VISIBLE` plus the
   phase-specific bit:
   `i==0→SE_ECL_MAX_VISIBLE`, `i==2→SE_ECL_PARTBEG_VISIBLE`, `i==3→SE_ECL_PARTEND_VISIBLE`,
   `i==4→SE_ECL_TOTBEG_VISIBLE`, `i==5→SE_ECL_TOTEND_VISIBLE`, `i==6→SE_ECL_PENUMBBEG_VISIBLE`,
   `i==7→SE_ECL_PENUMBEND_VISIBLE`. (Loop runs backward 7→0 but this only affects evaluation
   order, not the OR'd result — order-independent accumulation.)
5. If **no** phase was visible (`!(retflag & SE_ECL_VISIBLE)`): jump `tjd_start` forward (or
   backward) by exactly **25 days** (`tret[0] ± 25`) and `goto next_lun_ecl` — 25 days is
   roughly one synodic-month minus a few days, cheaply skipping past the just-rejected eclipse
   without needing to recompute `K` bookkeeping (the next `swe_lun_eclipse_when` call
   re-derives its own `K` from the new `tjd_start`).
6. **Moonrise/moonset clipping** (lines 3679–3715): `tjd_max = tret[0]` initially.
   - `swe_rise_trans(tret[6] - 0.001, SE_MOON, NULL, ifl, SE_CALC_RISE|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjdr, serr)` — next moonrise (bottom-limb convention) at/after just before penumbral
     begin.
   - Similarly `SE_CALC_SET` for moonset → `tjds`.
   - If `retc < 0` from either call (no rise/set found in the search window — e.g. circumpolar
     Moon at high latitude), skip the clipping block entirely (falls through with `tjd_max`
     unchanged).
   - If moonset `tjds` precedes penumbral begin, **or** moonset falls after moonrise which
     itself is after penumbral end (`tjds < tret[6] || (tjds > tjdr && tjdr > tret[7])`) — i.e.
     degenerate/non-overlapping rise-set ordering relative to the eclipse window — reject this
     candidate entirely: jump `tjd_start` by ±25 days and retry (`goto next_lun_ecl`), same as
     step 5.
   - If moonrise `tjdr` falls strictly inside the penumbral window
     (`tret[6] < tjdr < tret[7]`): the eclipse begins before moonrise, so zero out
     `tret[6]` (penumbral begin no longer applicable — it happened before the Moon rose) and
     any of `tret[2..5]` (partial/total begin markers) that are `< tjdr`; set `tret[8] = tjdr`
     (moonrise time); if `tjdr > tret[0]` (maximum itself occurs after moonrise, i.e. moonrise
     happens mid-eclipse before maximum — wait, actually the condition as written is
     `tjdr > tret[0]`, meaning moonrise is *after* the previously-computed maximum), update
     `tjd_max = tjdr` (the visible "moment of maximum" becomes moonrise itself, since maximum
     already passed before the Moon rose — the visible eclipse effectively "starts already in
     progress" at moonrise, and moonrise stands in for the reportable maximum time).
   - Symmetric handling for moonset `tjds` falling inside the penumbral window: zero
     `tret[7]` and any of `tret[2..5]` that are `> tjds`, set `tret[9] = tjds`, and if
     `tjds < tret[0]` (moonset before original maximum), `tjd_max = tjds`.
7. `tret[0] = tjd_max` (final, possibly rise/set-clipped, maximum/reportable time).
8. Recompute `attr[]` at the final `tjd_max` via `swe_lun_eclipse_how(tjd_max, ifl, geopos,
   attr, serr)` (public wrapper again — refreshes az/alt/magnitudes at the definitive time).
   If this returns 0 (e.g. `tjd_max` got clipped to a moonrise/moonset instant where the
   `swe_lun_eclipse_how`'s own altitude check (§3 step 6) fails), reject and retry with the
   ±25-day jump exactly as steps 5/6.
9. Final return: `retflag |= (retflag2 & SE_ECL_ALLTYPES_LUNAR)` — merges the eclipse-type bit
   (`SE_ECL_TOTAL`/`PARTIAL`/`PENUMBRAL`) from the final `attr`-refresh call into the
   visibility-bits accumulator from step 4, giving a combined return value carrying both the
   eclipse type and which phases were visible.

**FP/logic hazard**: steps 5, the rejection branch inside step 6, and step 8's rejection
branch all use the *same* `tjd_start = tret[0] ± 25; goto next_lun_ecl;` pattern — but note
`tret[0]` at each of these three points may hold a different quantity (the original
geocentric-maximum time from `swe_lun_eclipse_when`, vs. the pre-`tjd_max`-reassignment
`tret[0]` in step 6's rejection, vs... — check which `tret[0]` value is live at each `goto`
site when porting; do not assume they're interchangeable across the three rejection sites).

---

## 6. STATELESS PORT NOTES summary

- `swi_set_tid_acc` (global tidal-acceleration/ΔT model): called in both
  `swe_lun_eclipse_how` and `swe_lun_eclipse_when` before any position calc. Rust: thread the
  equivalent config through `&self`.
- `swe_set_topo` (global topocentric observer lon/lat/alt): called in `swe_lun_eclipse_how`
  step 6 before the topocentric Moon position for az/alt. Rust: pass `geopos` explicitly to
  the position/az-alt call instead of mutating shared state.
- All three public functions internally call `swe_calc`/`swe_calc_ut`/`swe_deltat_ex` which in
  the C library read/write additional global caches (ephemeris file handles, last-computed
  planet positions) — the stateless Rust `Ephemeris` recomputes explicitly on every call; no
  eclipse-specific caching behavior needs to be replicated (unlike the SPEED3/file-boundary
  caveat documented in the top-level `CLAUDE.md`, nothing in this module depends on
  cross-call cache reuse — every value needed is recomputed fresh within each function call in
  the C source itself).

---

## 7. Return-flag conventions recap

- `lun_eclipse_how` / `swe_lun_eclipse_how`: 0 (no eclipse) or exactly one of
  `SE_ECL_TOTAL`/`SE_ECL_PARTIAL`/`SE_ECL_PENUMBRAL`; `ERR` (-1) on ephemeris failure.
- `swe_lun_eclipse_when`: same type bit as above (never 0 on success — the search loop retries
  until it finds a matching eclipse), or `ERR`.
- `swe_lun_eclipse_when_loc`: type bit (from `SE_ECL_ALLTYPES_LUNAR`) OR'd with visibility bits
  (`SE_ECL_VISIBLE` plus zero or more of `*_VISIBLE` per contact), or `ERR`. Never returns a
  "not visible at all" result — the search loop retries with `tjd_start ± 25` days until it
  finds a visible occurrence.
