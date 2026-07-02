# C Reference: Planetary Phenomena — swecl.c

Porting reference for `swe_pheno` / `swe_pheno_ut` (phase angle, phase fraction, elongation,
apparent diameter, apparent magnitude, and — for the Moon — horizontal parallax). Read this
instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_pheno` | swecl.c:3802–4123 | Yes — the whole computation |
| `swe_pheno_ut` | swecl.c:4125–4142 | Yes — UT wrapper, deltaT re-call pattern |
| `swi_dot_prod_unit` | swephlib.c:453–465 | Yes (likely already ported — normalized dot product, clamped to [-1,1]) |
| `mag_elem` table | swecl.c:3773–3801 | Yes — transcribe verbatim (see below) |
| `pla_diam` table | sweph.h:315–333 | Already ported as `PLANETARY_DIAMETERS` in `src/constants.rs` (see Porting Notes) |

## Output Layout (`attr[]`, doc comment at swecl.c:3744–3750)

| Index | Meaning |
|---|---|
| `attr[0]` | phase angle (Earth–planet–Sun), degrees |
| `attr[1]` | phase (illuminated fraction of disc), 0..1 |
| `attr[2]` | elongation of planet (Earth–Sun angle as seen... actually angle Sun–Earth–planet), degrees |
| `attr[3]` | apparent diameter of disc, degrees |
| `attr[4]` | apparent magnitude |
| `attr[5]` | geocentric (or topocentric, if `SEFLG_TOPOCTR`) horizontal parallax — Moon only |
| `attr[6..19]` | unused; C zeroes all of `attr[0..20]` at entry (`for (i = 0; i < 20; i++) attr[i] = 0;`, swecl.c:3817–3818) and never writes indices 6–19 in this function |

Caller must allocate `attr[20]` minimum (C comment: "declare as attr[20] at least!").

## Constants Used

| Name | Value | Location | Notes |
|---|---|---|---|
| `NMAG_ELEM` | `SE_VESTA + 1` = 21 | swecl.c:3759 | size of `mag_elem` table |
| `NDIAM` | `SE_VESTA + 1` = 21 | sweph.h:314 | size of `pla_diam` table (same size, different table) |
| `MAG_MALLAMA_2018` | `1` (compile-time flag, always on in current source) | swecl.c:3760 | selects the Mallama/Hilton 2018 branch; the `#else` branch (old Meeus-only Mercury/Venus/Mars/Jupiter, old Saturn ring formula) is dead code in this checkout — port only the `#if` branch |
| `MAG_MOON_VREIJS` | `1` (compile-time flag, always on) | swecl.c:3761 | selects the two-regime Moon magnitude; `#else` (single Allen formula) is dead code — port only the `#if` branch |
| `EULER` | `2.718281828459` | swecl.c:3758 | **not** `f64::consts::E` (`2.718281828459045...`) — C uses this truncated literal for the Bowell HG phase functions; match it exactly for bit-fidelity |
| `AUNIT` | `1.49597870700e11` m | sweph.h:273 | already in `src/constants.rs` |
| `CLIGHT` | `2.99792458e8` m/s | sweph.h:274 | already in `src/constants.rs` |
| `EARTH_RADIUS` | `6378136.6` m | sweph.h:282 | already in `src/constants.rs` |
| `J2000` | `2451545.0` | sweph.h:67 | already in `src/constants.rs`? (verify; used for Saturn ring `T`) |
| `DEGTORAD` / `RADTODEG` | `PI/180.0` / `180.0/PI` | sweodef.h:262–266 (`M_PI`-based; the hard-coded literals at 262–263 are shadowed/overridden by 265–266 in the actual build) | already in `src/constants.rs` |
| `SE_AST_OFFSET` | `10000` | swephexp.h:128 | asteroid numbering offset |
| `SE_CERES..SE_VESTA` | `17..20` | swephexp.h:118–121 | main-belt asteroids with `mag_elem` rows |
| `SE_CHIRON` | `15` | swephexp.h:116 | boundary: `ipl < SE_CHIRON` is the "planet" `mag_elem` branch; `ipl >= SE_CHIRON` is the Bowell HG branch |
| `SE_PLUTO` | `9` | swephexp.h:110 | |
| body index order 0–20 | Sun,Moon,Mercury,Venus,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto,MeanNode,TrueNode,MeanApog,OscuApog,Earth,Chiron,Pholus,Ceres,Pallas,Juno,Vesta | swephexp.h:101–121 | indexes both `mag_elem` and `pla_diam` |

## `mag_elem[NMAG_ELEM][4]` Table (swecl.c:3773–3801)

Transcribed verbatim. Columns are `[H_or_base, G_or_c1, c2, c3]` — semantics vary by body class
(see algorithm sections below); `99` in column 0 is a sentinel meaning "no simple formula, see
special-case code or ficticious-body branch".

| idx | Body | col0 | col1 | col2 | col3 |
|---|---|---|---|---|---|
| 0 | Sun | −26.86 | 0 | 0 | 0 |
| 1 | Moon | −12.55 | 0 | 0 | 0 |
| 2 | Mercury *(obsolete placeholder, unused — Hilton/Mallama branch used instead; "don't delete this line")* | −0.42 | 3.80 | −2.73 | 2.00 |
| 3 | Venus *(obsolete placeholder, unused)* | −4.40 | 0.09 | 2.39 | −0.65 |
| 4 | Mars | −1.52 | 1.60 | 0 | 0 |
| 5 | Jupiter | −9.40 | 0.5 | 0 | 0 |
| 6 | Saturn | −8.88 | −2.60 | 1.25 | 0.044 |
| 7 | Uranus | −7.19 | 0.0 | 0 | 0 |
| 8 | Neptune | −6.87 | 0.0 | 0 | 0 |
| 9 | Pluto | −1.00 | 0.0 | 0 | 0 |
| 10 | Mean Node | 99 | 0 | 0 | 0 |
| 11 | True Node | 99 | 0 | 0 | 0 |
| 12 | Mean Apogee | 99 | 0 | 0 | 0 |
| 13 | Oscu Apogee | 99 | 0 | 0 | 0 |
| 14 | Earth | 99 | 0 | 0 | 0 |
| 15 | Chiron | 6.5 | 0.15 | 0 | 0 |
| 16 | Pholus | 7.0 | 0.15 | 0 | 0 |
| 17 | Ceres | 3.34 | 0.12 | 0 | 0 |
| 18 | Pallas | 4.13 | 0.11 | 0 | 0 |
| 19 | Juno | 5.33 | 0.32 | 0 | 0 |
| 20 | Vesta | 3.20 | 0.32 | 0 | 0 |

Note rows 2/3 (Mercury/Venus) are dead data — the live code always routes Mercury/Venus through
the dedicated Mallama polynomial branches below, never through `mag_elem[SE_MERCURY]` /
`mag_elem[SE_VENUS]`. Still transcribe them (the comment explicitly says not to delete them —
some other code path or future revert may reference them) but a Rust port only needs them if it
wants byte-for-byte table parity; the numeric values are otherwise inert.

## `pla_diam[]` Table (sweph.h:315–333, meters)

| idx | Body | Diameter (m) |
|---|---|---|
| 0 | Sun | 1,392,000,000.0 |
| 1 | Moon | 3,475,000.0 |
| 2 | Mercury | 2,439,400.0 × 2 = 4,878,800.0 |
| 3 | Venus | 6,051,800.0 × 2 = 12,103,600.0 |
| 4 | Mars | 3,389,500.0 × 2 = 6,779,000.0 |
| 5 | Jupiter | 69,911,000.0 × 2 = 139,822,000.0 |
| 6 | Saturn | 58,232,000.0 × 2 = 116,464,000.0 |
| 7 | Uranus | 25,362,000.0 × 2 = 50,724,000.0 |
| 8 | Neptune | 24,622,000.0 × 2 = 49,244,000.0 |
| 9 | Pluto | 1,188,300.0 × 2 = 2,376,600.0 |
| 10–13 | nodes/apogees | 0 |
| 14 | Earth | 6,371,008.4 × 2 = 12,742,016.8 |
| 15 | Chiron | 271,370.0 |
| 16 | Pholus | 290,000.0 |
| 17 | Ceres | 939,400.0 |
| 18 | Pallas | 545,000.0 |
| 19 | Juno | 246,596.0 |
| 20 | Vesta | 525,400.0 |

**Porting note:** this is almost certainly `src/constants.rs::PLANETARY_DIAMETERS`, already
ported and reused by `riseset.rs::disc_diameter_m` and `eclipse.rs::body_radius_au`. Verify the
values match this table (they should — same C source table) rather than re-deriving; if
`phenomena.rs` needs "diameter in meters for `Body`", prefer extending/reusing
`crate::constants::PLANETARY_DIAMETERS` + the `body.to_raw_id()` pattern from `eclipse.rs:105`
over duplicating a second copy of this table. Constraint: no duplicate logic across files.

## Algorithm: `swe_pheno` (swecl.c:3802–4123)

### 1. Setup and input sanitization (swecl.c:3811–3835)

```c
*serr2 = '\0';
iflag &= ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);
if (ipl == SE_AST_OFFSET + 134340)      // Pluto-as-asteroid → Pluto proper
  ipl = SE_PLUTO;
for (i = 0; i < 20; i++)
  attr[i] = 0;
if (ipl > SE_AST_OFFSET && ipl <= SE_AST_OFFSET + 4)   // Ceres..Vesta given as 10001..10004
  ipl = ipl - SE_AST_OFFSET - 1 + SE_CERES;
```

Then two independently-masked copies of `iflag` are built:

```c
iflag = iflag & (SEFLG_EPHMASK | SEFLG_TRUEPOS | SEFLG_J2000 | SEFLG_NONUT |
                 SEFLG_NOGDEFL | SEFLG_NOABERR | SEFLG_TOPOCTR);
iflagp = iflag & (SEFLG_EPHMASK | SEFLG_TRUEPOS | SEFLG_J2000 | SEFLG_NONUT | SEFLG_NOABERR);
iflagp |= SEFLG_HELCTR;
epheflag = iflag & SEFLG_EPHMASK;
```

Note `iflagp` is derived from the **already-masked** `iflag` (not the raw input), is missing
`SEFLG_NOGDEFL` and `SEFLG_TOPOCTR` relative to `iflag`, and forces `SEFLG_HELCTR`. `iflagp` is
used only for the heliocentric-at-`tjd-dt` positions (§3 below).

### 2. Geocentric planet position, twice (swecl.c:3839–3858)

```c
retflag = swe_calc(tjd, ipl, iflag | SEFLG_XYZ, xx, serr);   // xx: cartesian, for dot products
```
then re-derive `epheflag` from the **returned** flags (`retflag`) and, if the ephemeris that was
actually used differs from what was requested (fallback occurred), patch both `iflag` and
`iflagp`:
```c
epheflag2 = retflag & SEFLG_EPHMASG; // (SEFLG_EPHMASK)
if (epheflag != epheflag2) {
  iflag &= ~epheflag;  iflagp &= ~epheflag;
  iflag |= epheflag2;  iflagp |= epheflag2;
  epheflag = epheflag2;
}
swe_calc(tjd, ipl, iflag, lbr, serr);   // lbr: polar (lon,lat,dist) — same iflag, no XYZ, no SEFLG_XYZ
```
`lbr[2]` is the geocentric distance in AU (polar-coordinate output, third slot).

If `ipl == SE_MOON`, additionally fetch the Sun (cartesian) for the magnitude formula:
```c
swe_calc(tjd, SE_SUN, iflag | SEFLG_XYZ, xxs, serr);
```
(`xxs` is computed but — check — not actually read anywhere later in the function; the Moon
magnitude branch uses `lbr2` (Sun polar distance, from §3) not `xxs`. Confirm during
implementation whether `xxs` is dead; if so, the Rust port can skip this extra `swe_calc` call
for the Moon case, since `xxs` value is discarded. This is a candidate C inefficiency, not a
required side effect — no global state is written by `swe_calc` in this checkout's stateless
successor.)

### 3. Light-time-corrected heliocentric position (swecl.c:3859–3885)

Skipped entirely for `ipl` in `{SE_SUN, SE_EARTH, SE_MEAN_NODE, SE_TRUE_NODE, SE_MEAN_APOG,
SE_OSCU_APOG}` — for these, `attr[0]` and `attr[1]` stay `0` and `dt` stays `0`.

For every other body:
```c
dt = lbr[2] * AUNIT / CLIGHT / 86400.0;     // light time, days (AU·(m/AU)/(m/s)/(s/day))
if (iflag & SEFLG_TRUEPOS)
  dt = 0;
swe_calc(tjd - dt, ipl, iflagp | SEFLG_XYZ, xx2, serr);   // heliocentric cartesian, at tjd-dt
swe_calc(tjd - dt, ipl, iflagp, lbr2, serr);              // heliocentric polar, at tjd-dt
attr[0] = acos(swi_dot_prod_unit(xx, xx2)) * RADTODEG;     // phase angle
attr[1] = (1 + cos(attr[0] * DEGTORAD)) / 2;               // phase (illuminated fraction)
```
`xx` is geocentric-apparent cartesian (from §2, at `tjd`); `xx2` is heliocentric cartesian at
`tjd - dt` (light-time-corrected). The phase angle is the angle between the Earth→planet vector
and the Sun→planet vector, i.e. the standard Sun–planet–Earth phase angle, via the clamped
normalized dot product `swi_dot_prod_unit` (swephlib.c:453–465: `dop = (x·y)/(|x||y|)`, clamped
to `[-1, 1]`).

`lbr2[2]` (heliocentric distance at `tjd-dt`, AU) is reused throughout the magnitude formulas
below as "planet–Sun distance"; `lbr[2]` (geocentric distance at `tjd`, AU) is reused as
"planet–Earth distance". The magnitude formulas' recurring `5 * log10(lbr2[2] * lbr[2])` term is
the standard `5*log10(r*Δ)` brightness-vs-distance term (product form, not `log10(r)+log10(Δ)` —
match the multiplication-then-log grouping literally, not that it's numerically different, but
for FP-order fidelity keep the exact C expression shape).

### 4. Apparent diameter of disc, `attr[3]` (swecl.c:3886–3898)

```c
if (ipl < NDIAM)
  dd = pla_diam[ipl];
else if (ipl > SE_AST_OFFSET)
  dd = swed.ast_diam * 1000;    // km -> m (named-asteroid orbital-element data)
else
  dd = 0;
if (lbr[2] < dd / 2 / AUNIT)
  attr[3] = 180;                 // degenerate: "observer inside the body" — assume on surface
else
  attr[3] = asin(dd / 2 / AUNIT / lbr[2]) * 2 * RADTODEG;
```
Uses `lbr[2]`, the geocentric distance at `tjd` (not light-time corrected — apparent diameter is
evaluated at the apparent geocentric distance).

### 5. Apparent magnitude, `attr[4]` (swecl.c:3899–4068)

Guard: only computed if
```c
ipl > SE_AST_OFFSET || (ipl < NMAG_ELEM && mag_elem[ipl][0] < 99)
```
i.e. named/numbered asteroids beyond the offset, or any of the 21 built-in bodies whose
`mag_elem` row isn't the `99` sentinel (excludes nodes/apogees/Earth).

All branches below are mutually exclusive `if/else if` on `ipl`, evaluated in this order. Only
the live (`MAG_MALLAMA_2018` / `MAG_MOON_VREIJS`) branches are shown; dead `#else` code (old
Saturn ring / old Mercury-Venus-Mars-Jupiter formulas under `#if 0`-equivalent compile flags) is
noted but not to be ported.

#### 5a. Sun (swecl.c:3903–3907)

```c
fac = attr[3] / (asin(pla_diam[SE_SUN] / 2.0 / AUNIT) * 2 * RADTODEG);
fac *= fac;
attr[4] = mag_elem[SE_SUN][0] - 2.5 * log10(fac);
```
i.e. `-26.86 - 2.5*log10((attr[3]/avg_angular_diam)^2)` — brightness scales as inverse-square of
apparent size ratio (Sun's distance varies over the year).

#### 5b. Moon (swecl.c:3908–3929, `MAG_MOON_VREIJS` branch)

```c
double a = attr[0];                    // phase angle, NOT fabs()'d despite the commented-out line above it
if (a <= 147.1385465) {
    /* Allen 1976, Astrophysical Quantities */
    attr[4] = -21.62 + 0.026 * fabs(a) + 0.000000004 * pow(a, 4);
    attr[4] += 5 * log10(lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS);
} else {
    /* Samaha cube-phase-angle formula, VR-adjusted stitch point at 147.1385465° */
    attr[4] = -4.5444 - (2.5 * log10(pow(180 - a, 3)));
    attr[4] += 5 * log10(lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS);
}
```
Note: `pow(a, 4)` uses raw `a` (not `fabs(a)`) even though the linear term uses `fabs(a)`; since
`a` is a phase angle in `[0,180]` from `acos`, this is moot (`a >= 0` always) but transcribe
literally in case of future negative-input paths. The distance term
`lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS` uses `lbr[2]` = Moon geocentric distance (AU),
`lbr2[2]` = Moon heliocentric distance at `tjd-dt` (AU) — i.e. Earth-Moon and Sun-Moon distances,
scaled from AU to Earth-radii via `AUNIT/EARTH_RADIUS`.

The dead `#else` (non-`MAG_MOON_VREIJS`) single-formula variant (swecl.c:3927–3928) is:
```c
attr[4] = -21.62 + 5*log10(lbr[2]*lbr2[2]*AUNIT/EARTH_RADIUS) + 0.026*fabs(attr[0]) + 0.000000004*pow(attr[0],4);
```
Do not port; kept here only so a golden-test regression can be traced back to "which formula".

#### 5c. Mercury (swecl.c:3934–3938, Mallama 2018)

```c
double a = attr[0];
double a2=a*a; a3=a2*a; a4=a3*a; a5=a4*a; a6=a5*a;
attr[4] = -0.613 + a*6.3280E-02 - a2*1.6336E-03 + a3*3.3644E-05
                 - a4*3.4265E-07 + a5*1.6893E-09 - a6*3.0334E-12;
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
```
Coefficients (exact, verbatim): `-0.613`, `6.3280E-02`, `-1.6336E-03`, `3.3644E-05`,
`-3.4265E-07`, `1.6893E-09`, `-3.0334E-12`. Each power computed by repeated multiplication
(`a2=a*a; a3=a2*a; ...`), not `pow()` — match this for FP fidelity (Horner-vs-repeated-multiply
can differ in ULPs).

#### 5d. Venus (swecl.c:3939–3948, Mallama 2018, two regimes)

```c
double a = attr[0];
double a2=a*a; a3=a2*a; a4=a3*a;
if (a <= 163.7)
    attr[4] = -4.384 - a*1.044E-03 + a2*3.687E-04 - a3*2.814E-06 + a4*8.938E-09;
else
    attr[4] = 236.05828 - a*2.81914E+00 + a2*8.39034E-03;
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
if (attr[0] > 179.0)
    sprintf(serr2, "magnitude value for Venus at phase angle i=%.1f is bad; formula is valid only for i < 179.0", attr[0]);
```
Note: the `serr2` warning does **not** abort or change `attr[4]` — it's advisory only, copied
into `serr` at the very end of the function (§8) if non-empty. This is not a fatal error; the
computed (out-of-validity-range) magnitude is still returned.

#### 5e. Mars (swecl.c:3949–3967, Mallama 2018, two regimes; sub-Earth/vernal-equinox terms omitted per C comment)

```c
double a = attr[0]; double a2 = a*a;
if (a <= 50.0)
    attr[4] = -1.601 + a*0.02267 - a2*0.0001302;
else   // irrelevant to earth-centered observation (max Mars phase angle ~45°)
    attr[4] = -0.367 - a*0.02573 + a2*0.0003445;
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
```

#### 5f. Jupiter (swecl.c:3968–3973, Mallama 2018; phase angle never exceeds 12°, single regime)

```c
double a = attr[0]; double a2 = a*a;
attr[4] = -9.395 - a*3.7E-04 + a2*6.16E-04;
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
```

#### 5g. Saturn (swecl.c:3974–3994, Mallama 2018 + Meeus ring geometry)

```c
double a = attr[0]; double sinB2;
T = (tjd - dt - J2000) / 36525.0;
in = (28.075216 - 0.012998*T + 0.000004*T*T) * DEGTORAD;   // ring-plane inclination
om = (169.508470 + 1.394681*T + 0.000412*T*T) * DEGTORAD;  // ascending-node longitude of ring plane
sinB  = sin(in)*cos(lbr[1]*DEGTORAD) *sin(lbr[0]*DEGTORAD - om)  - cos(in)*sin(lbr[1]*DEGTORAD);
sinB2 = sin(in)*cos(lbr2[1]*DEGTORAD)*sin(lbr2[0]*DEGTORAD - om) - cos(in)*sin(lbr2[1]*DEGTORAD);
sinB = fabs(sin((asin(sinB) + asin(sinB2)) / 2.0));
attr[4] = -8.914 - 1.825*sinB + 0.026*a - 0.378*sinB*pow(2.7182818, -2.25*a);
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
```
`sinB` is computed from the Earth-facing geocentric ecliptic lon/lat (`lbr[0]`, `lbr[1]`, at
`tjd`) — "how tilted the rings look from Earth". `sinB2` is the analogous quantity from the
Sun-facing heliocentric lon/lat (`lbr2[0]`, `lbr2[1]`, at `tjd-dt`) — "how tilted the rings look
from the Sun" (i.e. how much of the ring's sunlit face we see). They are averaged via
`asin`→mean→`sin` (angle-averaging, not linear-averaging of the sines), matching Meeus. `T` uses
`tjd - dt` (light-time corrected), not raw `tjd`. `pow(2.7182818, -2.25*a)` uses a **second,
different, shorter** hard-coded Euler's-number literal (`2.7182818`, 7 sig figs) — distinct from
the file-level `EULER` macro (`2.718281828459`, used only in the asteroid HG branch §5j). Do not
conflate the two; use `2.7182818` literally here.

The dead `#else` (pre-Mallama) Saturn branch (swecl.c:4018–4041) computed ring magnitude via
`u1`/`u2` (ring-plane position angles from `atan2`) and `du = swe_degnorm(u1-u2)` folded to
`<=10°`, combined with the `mag_elem[SE_SATURN]` row coefficients — not ported; kept for
traceability only.

#### 5h. Uranus (swecl.c:3995–4005, Mallama 2018, simplified — ignores sub-Earth latitude term)

```c
double a = attr[0]; double a2 = a*a;
double fi_ = 0;  // sub-Earth latitude, ignored (always 0)
attr[4] = -7.110 - 8.4E-04*fi_ + a*6.587E-3 + a2*1.045E-4;
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
attr[4] -= 0.05;    // compensates for fi_ always being 0 (empirical correction, ±0.03m residual per C comment)
```
Note `fi_` is a dead variable always `0` (its term always evaluates to `0`) — the subsequent
unconditional `-= 0.05` is a fixed compensation, not conditional. Both lines still needed for
literal fidelity/ordering even though `- 8.4E-04*fi_` is always exactly `-0`.

#### 5i. Neptune (swecl.c:4006–4016, Mallama 2018, three time regimes)

```c
if (tjd < 2444239.5)
    attr[4] = -6.89;
else if (tjd <= 2451544.5)
    attr[4] = -6.89 - 0.0055 * (tjd - 2444239.5) / 365.25;
else
    attr[4] = -7.00;
attr[4] += 5 * log10(lbr2[2] * lbr[2]);
```
Uses raw `tjd` (not light-time corrected) for the regime boundaries. C comment notes Mallama's
paper actually specifies `0.0054` for the middle-regime slope but JPL Horizons uses `0.0054` too
— yet the code hard-codes `0.0055`, deliberately, to keep the piecewise function continuous at
the `tjd=2444239.5` boundary (using `0.0054` would introduce a discontinuity there). Port
`0.0055` — this is intentional, not a bug to "fix".

#### 5j. `ipl < SE_CHIRON` generic branch (swecl.c:4042–4047)

Any remaining built-in body with index `< SE_CHIRON` (15) that fell through — in practice this is
dead given 5a–5i already cover Sun..Neptune and §3's exclusion list covers nodes/apogees/Earth,
but transcribe for completeness:
```c
attr[4] = 5 * log10(lbr2[2] * lbr[2])
            + mag_elem[ipl][1] * attr[0] / 100.0
            + mag_elem[ipl][2] * attr[0] * attr[0] / 10000.0
            + mag_elem[ipl][3] * attr[0] * attr[0] * attr[0] / 1000000.0
            + mag_elem[ipl][0];
```
(Old-style polynomial-in-phase-angle with mixed-magnitude-unit coefficient scaling: /100, /10000,
/1000000 for linear/quadratic/cubic terms respectively — note `attr[0]*attr[0]` and
`attr[0]*attr[0]*attr[0]` computed by repeated multiplication, not `pow`.)

#### 5k. Bowell HG system — asteroids/Chiron/Pholus (swecl.c:4048–4064)

Condition: `ipl < NMAG_ELEM || ipl > SE_AST_OFFSET` (i.e. `ipl >= SE_CHIRON` and either a
built-in body ≤ Vesta, or a numbered asteroid beyond the offset).
```c
ph1 = pow(EULER, -3.33 * pow(tan(attr[0]*DEGTORAD/2), 0.63));
ph2 = pow(EULER, -1.87 * pow(tan(attr[0]*DEGTORAD/2), 1.22));
if (ipl < NMAG_ELEM) {                       // Chiron/Pholus/Ceres/Pallas/Juno/Vesta
    me[0] = mag_elem[ipl][0];  me[1] = mag_elem[ipl][1];
} else if (ipl == SE_AST_OFFSET + 1566) {    // Icarus: JPL-database H/G override
    me[0] = 16.9;  me[1] = 0.15;
} else {                                     // other numbered asteroids
    me[0] = swed.ast_H;  me[1] = swed.ast_G;
}
attr[4] = 5 * log10(lbr2[2] * lbr[2])
            + me[0]
            - 2.5 * log10((1 - me[1]) * ph1 + me[1] * ph2);
```
`EULER` here is the file-level macro `2.718281828459` (13 sig figs) — do **not** substitute
`std::f64::consts::E`; use the literal to preserve bit-exact parity. `ph1`/`ph2` are the standard
IAU H-G two-term phase functions. `me[0]`/`me[1]` (H, G) come from three sources depending on
`ipl`: the built-in `mag_elem` table (Chiron/Pholus/Ceres/Pallas/Juno/Vesta), a hard-coded
override for asteroid 1566 Icarus, or `swed.ast_H`/`swed.ast_G` (global state populated from the
SE1 asteroid orbital-element file for any other numbered asteroid — see Porting Notes).

#### 5l. Fictitious bodies fallback (swecl.c:4065–4067)

```c
} else {  /* fictitious bodies */
    attr[4] = 0;
}
```
Reached only if none of 5a–5k matched but the outer guard passed — in practice unreachable given
the guard's own condition, but present in C; port as a defensive `0.0` default.

### 6. Elongation, `attr[2]` (swecl.c:4069–4078)

Skipped for `ipl == SE_SUN` or `ipl == SE_EARTH` (stays `0`).
```c
swe_calc(tjd, SE_SUN, iflag | SEFLG_XYZ, xx2, serr);   // re-fetch Sun, cartesian (overwrites xx2 from §3!)
swe_calc(tjd, SE_SUN, iflag, lbr2, serr);              // re-fetch Sun, polar   (overwrites lbr2 from §3!)
attr[2] = acos(swi_dot_prod_unit(xx, xx2)) * RADTODEG;
```
**Important reuse-of-buffer trap**: `xx2` and `lbr2` are the *same* C locals used in §3 for the
light-time-corrected heliocentric planet position — they get clobbered here with the Sun's
geocentric position/distance. This is safe in C only because §5 (magnitude) has already consumed
`lbr2[2]` by this point — the ordering in the function is magnitude (§5) *before* elongation
(§6), specifically so this reuse doesn't corrupt the magnitude calc. **A stateless/Rust port
must use separate named bindings** (not reuse a buffer) since there's no ordering hazard to
preserve — but note the ordering hazard exists in C, meaning if you ever want to double check
against a modified/reordered C build, moving elongation before magnitude would be a bug. `xx` is
the geocentric planet cartesian position from §2 (still valid, not clobbered).

`attr[2]` is described in the header comment as "elongation of planet" — i.e. angular separation
between planet and Sun as seen from Earth (Sun–Earth–planet angle), not the phase angle (which
is Sun–planet–Earth, computed in §3 from the opposite vertex).

### 7. Horizontal parallax, `attr[5]` — Moon only (swecl.c:4079–4119)

```c
if (ipl == SE_MOON) {
    swe_calc(tjd, SE_MOON, epheflag|SEFLG_TRUEPOS|SEFLG_EQUATORIAL|SEFLG_RADIANS, xm, serr);
    sinhp = EARTH_RADIUS / xm[2] / AUNIT;     // xm[2]: true geocentric distance, AU (radians-flag doesn't affect distance units)
    attr[5] = asin(sinhp) / DEGTORAD;          // NOTE: /DEGTORAD not *RADTODEG (equivalent, but literal form differs)
    if (iflag & SEFLG_TOPOCTR) {
        swe_calc(tjd, SE_MOON, epheflag|SEFLG_XYZ|SEFLG_TOPOCTR, xm, serr);   // topocentric cartesian
        swe_calc(tjd, SE_MOON, epheflag|SEFLG_XYZ, xx, serr);                 // geocentric cartesian (overwrites xx from §2/§6!)
        attr[5] = acos(swi_dot_prod_unit(xm, xx)) / DEGTORAD;
    }
}
```
Uses `epheflag` (just the ephemeris-source bits, recomputed/possibly-patched in §2), not the full
`iflag` — deliberately drops `SEFLG_NONUT`/`SEFLG_TRUEPOS` etc. that were in the outer `iflag`.
The geocentric-parallax formula is flagged in C as citing "Expl.Suppl. to the AA 1984, p.400":
`sinhp = EARTH_RADIUS/distance_m` (small-angle parallax), then `asin` for the exact angle. The
topocentric-parallax branch instead computes the *actual* angular displacement between the
geocentric and topocentric apparent directions via the dot-product angle — a direct geometric
measurement rather than the small-angle sine formula, and it clobbers `xx` (which is otherwise
already-consumed by this point — the last read of `xx` was §6's elongation, already complete).

A `#if 0`-disabled block (swecl.c:4096–4117) contains an alternative Expl.Suppl.-1984-p.400
formula using local sidereal time / hour angle / geocentric latitude of the observer, explicitly
noted in the C comment as **not** accounting for the Moon's topocentric distance or the
observer's distance from the geocenter — dead code, not ported.

### 8. Return (swecl.c:4120–4122)

```c
if (*serr2 != '\0' && serr != NULL)
    strcpy(serr, serr2);
return iflag;
```
`iflag` here is the *masked* copy from §1 (possibly patched with `epheflag2` in §2) — this is
the "flags actually used" return convention also seen elsewhere in the C library (fallback
signaling). Any advisory warning built up in `serr2` (currently only the Venus phase-angle
range warning from §5d) is copied into the caller's `serr` buffer, but does **not** change the
return value to `ERR` — it's advisory, not an error.

### Error/boundary handling summary

- No `ipl` value is outright rejected by name in this function; any `swe_calc` failure for any
  of the many internal calls (§2, §3, §5g's helper isn't a separate call, §6, §7) immediately
  `return ERR` with the C-populated `serr`.
- The only body-specific special-casing that changes `ipl` itself (not just which formula
  applies) is: Pluto-as-asteroid (`10000+134340` → `9`), and asteroid-offset Ceres..Vesta
  (`10001..10004` → `17..20`) — both at the very top (§1).
- Venus phase angle > 179° produces a non-fatal advisory `serr` message, not an error return.

## Algorithm: `swe_pheno_ut` (swecl.c:4125–4142)

```c
int32 swe_pheno_ut(double tjd_ut, int32 ipl, int32 iflag, double *attr, char *serr) {
  int32 epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag == 0) {
    epheflag = SEFLG_SWIEPH;
    iflag |= SEFLG_SWIEPH;              // default to SWIEPH if caller specified no ephemeris
  }
  deltat = swe_deltat_ex(tjd_ut, iflag, serr);
  retflag = swe_pheno(tjd_ut + deltat, ipl, iflag, attr, serr);
  if ((retflag & SEFLG_EPHMASK) != epheflag) {
    // ephemeris actually used differs from what was requested (fallback) —
    // deltaT can depend on which ephemeris/model is used, so recompute deltaT
    // with the *actual* flags and re-run swe_pheno with the corrected TT.
    deltat = swe_deltat_ex(tjd_ut, retflag, serr);
    retflag = swe_pheno(tjd_ut + deltat, ipl, iflag, attr, serr);
  }
  return retflag;
}
```
This is the same "conditional re-call after ephemeris-fallback deltaT correction" pattern used
elsewhere in the C library wherever a UT-based wrapper must convert to TT before delegating to
the TT-based core function. The second `swe_pheno` call still passes the *original* `iflag` (not
`retflag`) — only the deltaT input (`tjd_ut + deltat`) changes between the two calls.

## Porting Notes

**Target module**: `src/phenomena.rs` (currently an empty stub per `docs/codebase-map.md:225`).
`swe_calc` → `Ephemeris::calc` (`src/context.rs:193`). `swe_calc_ut`-style deltaT handling has
no direct analog needed here in the stateless port unless a `pheno_ut`-equivalent entry point is
added — but see below re: fallback re-call.

**Global state reads that a stateless port must handle differently:**

1. `swed.ast_diam`, `swed.ast_H`, `swed.ast_G` (sweph.h ~813-815) — populated as a side effect of
   a prior `swe_calc` call for a numbered asteroid (read from the SE1 orbital-element file
   header). In C, `swe_pheno` implicitly depends on `swe_calc(..., ipl, ...)` having *already*
   run (§2) and cached these globals before §4/§5k reads them. A stateless port cannot rely on
   call-order side effects — `Ephemeris::calc` (or whatever asteroid-loading path exists) must
   either return H/G/diameter alongside the position, or `phenomena.rs` must fetch them via an
   explicit, independent lookup rather than assuming a global was just populated. Check whether
   `CalcResult` or the asteroid-loading path already exposes these (search for `ast_diam`/`ast_H`
   in already-ported modules — as of this writing `eclipse.rs:105` explicitly notes this data
   "isn't threaded through a stateless config yet" and stubs asteroid diameter to `0.0`;
   `phenomena.rs` will likely hit the same gap for §4 and §5k's "other asteroids" sub-branch and
   should either share a fix with `eclipse.rs`/`riseset.rs` or explicitly stub/error the same way
   they do, rather than inventing a third stub).
2. `swed.topd.geolon` / `swed.topd.geolat` — read only inside the dead `#if 0` block in §7; not
   needed for the live code path.
3. No other `swed.*` global reads in the live code path of `swe_pheno`/`swe_pheno_ut`.

**Buffer-reuse hazards (C-only, do not replicate in Rust):** §3 and §6 reuse the same C locals
(`xx2`, `lbr2`) for two different physical quantities at two different points in the function,
relying on the fact that §5 (magnitude) fully consumes the §3 values before §6 overwrites them.
§7's topocentric branch reuses `xx` (last legitimately read at end of §6) for yet another
purpose. A Rust port should use distinct, clearly-named bindings for each (helio-planet-at-dt
cartesian/polar, sun-geocentric cartesian/polar, moon-topocentric-vs-geocentric cartesian) —
there is no reason to preserve C's variable reuse, only its *value* semantics at each read site.

**Existing reusable code in this repo:**
- `crate::constants::PLANETARY_DIAMETERS` — likely already the `pla_diam` table; verify values
  against the table above instead of re-adding it.
- `eclipse.rs::body_radius_au` (swecl.c:697-704/1004-1011 pattern) — same
  `PLANETARY_DIAMETERS[raw_id]/2/AUNIT` shape as this function's §4 diameter lookup (minus the
  `swed.ast_diam` fallback, which `eclipse.rs` also stubs to `0.0`). Consider whether §4's
  disc-diameter resolution (in meters, before the `/2/AUNIT` and `asin` steps) can share a helper
  with `eclipse.rs`/`riseset.rs`'s diameter lookups rather than re-deriving a third copy — the
  `constants::PLANETARY_DIAMETERS` + `body.to_raw_id()` pattern from `eclipse.rs:105` and
  `riseset.rs:463` is the established idiom in this codebase.
- `swi_dot_prod_unit` — check if already ported (likely a small free function somewhere shared,
  e.g. `src/math.rs` or similar) before adding a new one; it's used identically here (§3 phase
  angle, §6 elongation, §7 topocentric parallax) and elsewhere in the C source (swehouse.c,
  swecl.c eclipse code).

**Flag-masking fidelity:** §1's two independently-built flag masks (`iflag` vs `iflagp`) are easy
to collapse into one by accident when porting — they differ specifically in `SEFLG_NOGDEFL` and
`SEFLG_TOPOCTR` (present in `iflag`, absent from `iflagp`) and in the forced `SEFLG_HELCTR` (only
on `iflagp`). Keep them as two distinct `CalcFlags` values in the Rust port, derived in the same
order (`iflagp` derived from the masked `iflag`, not from the raw input `iflag` parameter).

**Dead compile-time branches, not to be ported:** the `#else` halves of `MAG_MALLAMA_2018` and
`MAG_MOON_VREIJS` (old single-formula Moon magnitude; old Meeus-only Mercury/Venus/Mars/Jupiter;
old ring-geometry-based Saturn via `u1`/`u2`/`du` and `mag_elem[SE_SATURN]`). Both flags are `1`
unconditionally in this checkout — there is no runtime toggle, so there is nothing to expose in
the Rust API for these; port only the `#if` branches documented in §5.

**Testing surface:** `attr[0..5]` (phase angle, phase, elongation, diameter, magnitude,
Moon-parallax) are all independently checkable against `swetest` golden output per body class:
Sun/Moon/Mercury/Venus/Mars/Jupiter/Saturn/Uranus/Neptune/Pluto each exercise a distinct magnitude
formula branch (§5a–5i), and Chiron/Pholus/Ceres/Pallas/Juno/Vesta exercise the Bowell HG branch
(§5k) via the `mag_elem` table rather than `swed.ast_H/ast_G`. Numbered asteroids beyond Vesta
would additionally exercise the `swed.ast_H/ast_G`/Icarus-override sub-paths of §5k, but those
depend on the still-unresolved asteroid-orbital-element-loading gap noted above.
