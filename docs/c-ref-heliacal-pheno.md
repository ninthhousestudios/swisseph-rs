# C Reference: Heliacal — Pheno & Lunar Crescent (swehel.c part 2)

Porting reference for the second half of the heliacal-visibility machinery: the Moon-crescent
(Yallop) geometry, the small parabola/line-crossing interpolation helpers used to bracket
visibility windows, `DeterTAV` (topocentric arcus-visionis determination for one instant), and
`swe_heliacal_pheno_ut` — the public "phenomena at a given instant" entry point that assembles all
of the above (plus rise/set and magnitude data) into the 28-slot `darr[]` output array.

Sibling doc `docs/c-ref-heliacal-vision.md` covers the internals of `ObjectLoc`, `Magnitude`,
`TopoArcVisionis`, `RiseSet`, `AppAltfromTopoAlt`, `default_heliacal_parameters`, `SunRA`,
`DeterObject`, `kt` and the visual-limiting-magnitude machinery they call into — only the call
**contracts** (signatures, argument meaning/units, return conventions) of those are reproduced here,
as needed to read this file's algorithms.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `WidthMoon` | swehel.c:1715–1720 | Yes |
| `LengthMoon` | swehel.c:1726–1734 | Yes |
| `qYallop` | swehel.c:1741–1745 | Yes |
| `crossing` | swehel.c:1753–1756 | Yes |
| `DeterTAV` | swehel.c:1759–1783 | Yes |
| `x2min` | swehel.c:1791–1797 | Yes |
| `funct2` | swehel.c:1807–1810 | Yes |
| `strcpy_VBsafe` | swehel.c:1812–1824 | No — string-buffer sanitization for the C VB-style fixed API; a Rust port takes an owned `&str`/`String` and has no analogous buffer-safety concern. Document only for traceability (§5). |
| `swe_heliacal_pheno_ut` | swehel.c:1862–2074 | Yes — the whole computation |
| `Sgn` | swehel.c:864–869 | Yes (tiny helper, used once in §4) |
| `mymin` / `mymax` | swehel.c:155–167 | Yes (tiny helpers, used repeatedly in §4) — **not** `f64::min`/`f64::max`: `mymin(a,b)` returns `a` when `a<=b` else `b`; `mymax(a,b)` returns `a` when `a>=b` else `b`. Semantically identical to IEEE min/max for non-NaN inputs, but transcribe as written since no NaNs are expected on this path. |

## Constants Used

| Name | Value | Location | Notes |
|---|---|---|---|
| `AvgRadiusMoon` | `15.541 / 60` deg | swehel.c:112 | "at 2007 CE or BCE" per C comment — average Moon angular radius, degrees. Used only as `LengthMoon`'s fallback diameter (`AvgRadiusMoon * 2`) when `Diamoon == 0`; **`swe_heliacal_pheno_ut` always calls `LengthMoon(WMoon, 0)`**, so this fallback path is always taken in the pheno walkthrough (§4) — the real lunar diameter is never threaded through here. |
| `MaxTryHours` | `4` | swehel.c:84 | main search-loop time bound (hours) around `RiseSetS`, §4 |
| `TimeStepDefault` | `1` | swehel.c:85 | main search-loop step, **minutes** |
| `LocalMinStep` | `8` | swehel.c:86 | local-minimum-confirmation look-ahead step, **minutes**, §4 |
| `TJD_INVALID` | `99999999.0` | swephexp.h:450 | sentinel for "not applicable" JD outputs (`#NA!` in the original VB) |
| `SEI_ECL_GEOALT_MIN` / `_MAX` | `-500.0` / `25000.0` (m) | sweph.h:198–199 | valid observer-altitude bounds, checked at entry |
| `SE_MORNING_FIRST` (`TypeEvent==1`) | `= SE_HELIACAL_RISING = 1` | swephexp.h:424,426 | morning first appearance |
| `SE_EVENING_LAST` (`TypeEvent==2`) | `= SE_HELIACAL_SETTING = 2` | swephexp.h:425,427 | evening last visibility |
| `SE_EVENING_FIRST` (`TypeEvent==3`) | `3` | swephexp.h:428 | evening first appearance (post-conjunction) |
| `SE_MORNING_LAST` (`TypeEvent==4`) | `4` | swephexp.h:429 | morning last visibility (pre-conjunction) |
| `SE_SUN..SE_JUPITER` | `0,1,2,3,4,5` | swephexp.h:101–106 | body-index constants used in the `Planet >= SE_MARS` guard, §4 |
| `AS_MAXCH` | (string buffer size constant, elsewhere in the C headers) | — | size of the local `ObjectName[AS_MAXCH]` buffer in §4; irrelevant to a Rust port using owned strings |

## §1 Moon crescent geometry (`WidthMoon` 1715, `LengthMoon` 1726, `qYallop` 1741)

### `WidthMoon` (swehel.c:1715–1720)

```c
static double WidthMoon(double AltO, double AziO, double AltS, double AziS, double parallax)
{
  /* Yallop 1998, page 3*/
  double GeoAltO = AltO + parallax;
  return 0.27245 * parallax * (1 + sin(GeoAltO * DEGTORAD) * sin(parallax * DEGTORAD))
       * (1 - cos((AltS - GeoAltO) * DEGTORAD) * cos((AziS - AziO) * DEGTORAD));
}
```

- Inputs: `AltO` = topocentric altitude of the Moon (deg), `AziO` = azimuth of the Moon (deg),
  `AltS` = topocentric altitude of Sun (deg), `AziS` = azimuth of Sun (deg), `parallax` = Moon's
  horizontal parallax (deg) — in the caller (§4) this is `ParO = GeoAltO - AltO` (see §4 step 2),
  **not** a value from `swe_pheno`.
- `GeoAltO = AltO + parallax` recovers the geocentric altitude from the topocentric one and the
  parallax (local reconstruction, independent of the `GeoAltO` computed separately in §4 via
  `ObjectLoc(..., Angle=7, ...)`).
- Output: crescent width `W`, degrees. Transcribe the two product terms exactly as grouped
  (`0.27245 * parallax * (...) * (...)`, left-to-right) — do not reassociate.

### `LengthMoon` (swehel.c:1726–1734)

```c
static double LengthMoon(double W, double Diamoon)
{
  double Wi, D;
  if (Diamoon == 0) Diamoon = AvgRadiusMoon * 2;
  Wi = W * 60;
  D = Diamoon * 60;
  /* Crescent length according: http://calendar.ut.ac.ir/Fa/Crescent/Data/Sultan2005.pdf*/
  return (D - 0.3 * (D + Wi) / 2.0 / Wi) / 60.0;
}
```

- `W` deg (crescent width, from `WidthMoon`), `Diamoon` deg (Moon's angular diameter; **always
  passed as `0` by `swe_heliacal_pheno_ut`**, so `Diamoon` is always replaced by
  `AvgRadiusMoon * 2` inside this function on the live call path — see Constants table).
- `Wi`, `D` are in arcminutes (`* 60`). Division order: `0.3 * (D + Wi) / 2.0 / Wi` — left-to-right
  chained division, i.e. `((0.3 * (D + Wi)) / 2.0) / Wi`, not `0.3*(D+Wi)/(2.0*Wi)` — identical
  result mathematically but transcribe the chain form for FP-order fidelity.
- Output: crescent length `LMoon`, degrees (final `/ 60.0` converts back from arcminutes).

### `qYallop` (swehel.c:1741–1745) + grade table (used in §4)

```c
static double qYallop(double W, double GeoARCVact)
{
  double Wi = W * 60;
  return (GeoARCVact - (11.8371 - 6.3226 * Wi + 0.7319 * Wi * Wi - 0.1018 * Wi * Wi * Wi)) / 10;
}
```

- `W` deg (crescent width), `GeoARCVact` deg (geocentric arcus visionis, i.e. `ARCVact` from §4 —
  despite the parameter name saying "Geo", the caller passes the plain `ARCVact` computed in §4,
  which is itself geocentric: `ARCVact = TAVact + ParO`, see §4 step 2).
- `Wi = W * 60` — crescent width in **arcminutes** (the cubic polynomial `11.8371 - 6.3226 Wi +
  0.7319 Wi² - 0.1018 Wi³` is Yallop's empirical fit in arcminutes, matching the coefficients
  quoted in `catalogue-internal.md` §15.5). Each power of `Wi` is built by repeated multiplication
  (`Wi * Wi`, `Wi * Wi * Wi`), not `pow()` — match for FP fidelity.
- Output: dimensionless `q` value, `/10` at the end (Yallop's `q` is scaled by 10 from the raw
  ARCV-minus-polynomial difference in arcminutes-equivalent units).

**Grade table** — assigned only in §4 (swehel.c:1935–1940), immediately after `qYal` is computed,
using **this exact checkout's thresholds** (transcribe these, not the paraphrase in
`catalogue-internal.md` §15.5, which uses different boundary values and a 5-grade A–E scheme; this
C file uses a 6-grade A–F scheme):

```c
if (qYal > 0.216) qCrit = 1;                        /* A */
if (qYal < 0.216 && qYal > -0.014) qCrit = 2;        /* B */
if (qYal < -0.014 && qYal > -0.16) qCrit = 3;        /* C */
if (qYal < -0.16 && qYal > -0.232) qCrit = 4;        /* D */
if (qYal < -0.232 && qYal > -0.293) qCrit = 5;       /* E */
if (qYal < -0.293) qCrit = 6;                        /* F */
```

Notes:
- These are independent `if`s (not `else if`), but the ranges are constructed to be mutually
  exclusive and exhaustive except for the exact boundary values `qYal == 0.216`, `-0.014`, `-0.16`,
  `-0.232`, `-0.293` themselves, at which **no** branch fires and `qCrit` keeps its initial value
  `0` (set at swehel.c:1929, just before the `if (Planet == SE_MOON)` block). A Rust port should
  replicate this open/open interval structure literally (e.g. a `match`/`if` chain with the same
  strict `<`/`>` comparisons) rather than "fixing" the boundaries to be closed on one side — the
  boundary-value behavior (falling through to `0`/no grade) is exact C behavior worth preserving
  test coverage for.
- `qCrit` (and `WMoon`, `qYal`, `LMoon`) are computed **only when `Planet == SE_MOON`**; for every
  other body they are `0` (initialized at swehel.c:1927–1930, just before the Moon check) and stay
  `0` all the way to `darr[16..18]`/`darr[25]` output.

## §2 Interpolation helpers (`crossing` 1753, `x2min` 1791, `funct2` 1807)

These three are generic 2-point-line / 3-point-parabola interpolators, reused for two different
purposes in §4's visibility-search loop: `crossing` finds where a moving quantity first exceeds a
threshold (linear, 2-point); `x2min`/`funct2` find/evaluate the vertex of a parabola through 3
equally-spaced samples (used to refine the "best visibility time" local minimum).

### `crossing` (swehel.c:1753–1756)

```c
/*###################################################################
'A (0,p)
'B (1,q)
'C (0,r)
'D (1,s)
*/
static double crossing(double A, double B, double C, double D)
{
  return (C - A) / ((B - A) - (D - C));
}
```

- Doc-comment convention (from the VB original): two lines, `line 1` through points `A=(0,p)` and
  `B=(1,q)`; `line 2` through points `C=(0,r)` and `D=(1,s)`. Returns the **x**-coordinate (in
  `[0,1]`-normalized units, extrapolable outside that range) where the two lines cross, i.e. the
  solution of `p + (q-p)x = r + (s-r)x`, which rearranges to exactly the expression above (`A=p,
  B=q, C=r, D=s`).
- Call site (§4, swehel.c:2010/2014): `crossing(DeltaAltoud, DeltaAlt, MinTAVoud, MinTAVact)` — line
  1 is "object-minus-sun altitude difference" at the previous sample (`x=0`, value `DeltaAltoud`)
  and current sample (`x=1`, value `DeltaAlt`); line 2 is "visibility threshold" at previous
  (`MinTAVoud`) and current (`MinTAVact`) sample. The returned `crosspoint` is where the
  altitude-difference curve crosses the (also time-varying) visibility threshold, linearly
  interpolated between the previous and current 1-minute samples.
- Converting `crosspoint` (∈ roughly `[0,1]`) back to an absolute JD: `Tc = TimePointer - TimeStep *
  (1 - crosspoint)` — at `crosspoint=1` this gives `Tc = TimePointer` (crossing exactly at the
  current sample); at `crosspoint=0` it gives `Tc = TimePointer - TimeStep` (crossing at the
  previous sample). `TimePointer` and `TimeStep` here are the enclosing loop's current time and
  (signed) step — see §4.

### `x2min` (swehel.c:1791–1797)

```c
/*###################################################################
' A y-value at x=1
' B y-value at x=0
' C y-value at x=-1
' x2min minimum for the quadratic function
*/
static double x2min(double A, double B, double C)
{
  double term = A + C - 2 * B;
  if (term == 0)
    return 0;
  return -(A - C) / 2.0 / term;
}
```

- Three samples of a quadratic at **unit-spaced abscissas** `x = 1` (value `A`), `x = 0` (value
  `B`), `x = -1` (value `C`) — note `B` (not `A`) is the `x=0` sample, i.e. the *middle* of the
  three points chronologically is bound to parameter `B`. Returns the `x`-location of the
  parabola's vertex (min or max, whichever it is) in this same `x∈{-1,0,1}`-centered coordinate
  system. Degenerate case `term == 0` (the 3 points are exactly colinear, no curvature) returns `0`
  (i.e. "vertex at the middle sample") rather than dividing by zero.
- Division order: `-(A - C) / 2.0 / term`, i.e. `(-(A-C) / 2.0) / term` — chained left-to-right;
  transcribe literally, not as `-(A-C)/(2.0*term)`.
- Call site (§4, swehel.c:2004): `x2min(MinTAVact, MinTAVoud, OldestMinTAV)` — `A=MinTAVact` (the
  most recent/current sample, bound to the `x=1` slot), `B=MinTAVoud` (previous sample, `x=0`),
  `C=OldestMinTAV` (two-samples-back, `x=-1`). So despite the docstring's generic "y at x=1/0/-1"
  framing, the *physical* time order is: `C` (oldest) → `B` (middle) → `A` (newest/current), and
  the returned `extrax` is a vertex-x measured with the **previous** sample as origin (`x=0`) and
  one step forward as `x=1`.
- Converting `extrax` to an absolute JD (swehel.c:2005): `TbVR = TimePointer - (1 - extrax) *
  TimeStep` — same conversion shape as `crossing`'s, consistent with `x=1` ↔ `TimePointer`
  (current) and `x=0` ↔ `TimePointer - TimeStep` (previous).

### `funct2` (swehel.c:1807–1810)

```c
/*###################################################################
' A y-value at x=1
' B y-value at x=0
' C y-value at x=-1
' x
' y is y-value of quadratic function
*/
static double funct2(double A, double B, double C, double x)
{
  return (A + C - 2 * B) / 2.0 * x * x + (A - C) / 2.0 * x + B;
}
```

- Evaluates the same 3-point-fit quadratic (same `A`/`B`/`C` convention as `x2min`) at an arbitrary
  `x`. Call site (§4, swehel.c:2006): `funct2(MinTAVact, MinTAVoud, OldestMinTAV, extrax)` — same
  argument binding as the `x2min` call immediately above it, evaluated at the vertex `x=extrax`
  just computed, giving the interpolated minimum arcus-visionis threshold value `MinTAV`.
- Expression grouping: `(A + C - 2*B)/2.0 * x*x + (A - C)/2.0 * x + B` — transcribe the
  division-before-multiplication grouping (`(A+C-2*B)/2.0` computed first, then `* x * x`) exactly.

## §3 `DeterTAV` (swehel.c:1759–1783)

```c
static int32 DeterTAV(double *dobs, double JDNDaysUT, double *dgeo, double *datm,
                       char *ObjectName, int32 helflag, double *dret, char *serr)
{
  double Magn, AltO, AziS, AziO, AziM, AltM;
  double sunra = SunRA(JDNDaysUT, helflag, serr);
  if (Magnitude(JDNDaysUT, dgeo, ObjectName, helflag, &Magn, serr) == ERR)
    return ERR;
  if (ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 0, helflag, &AltO, serr) == ERR)
    return ERR;
  if (ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 1, helflag, &AziO, serr) == ERR)
    return ERR;
  if (strncmp(ObjectName, "moon", 4) == 0) {
    AltM = -90;
    AziM = 0;
  } else {
    if (ObjectLoc(JDNDaysUT, dgeo, datm, "moon", 0, helflag, &AltM, serr) == ERR)
      return ERR;
    if (ObjectLoc(JDNDaysUT, dgeo, datm, "moon", 1, helflag, &AziM, serr) == ERR)
      return ERR;
  }
  if (ObjectLoc(JDNDaysUT, dgeo, datm, "sun", 1, helflag, &AziS, serr) == ERR)
    return ERR;
  if (TopoArcVisionis(Magn, dobs, AltO, AziO, AltM, AziM, JDNDaysUT, AziS, sunra,
                       dgeo[1], dgeo[2], datm, helflag, dret, serr) == ERR)
    return ERR;
  return OK;
}
```

Purpose: for a single instant `JDNDaysUT`, compute the **topocentric arcus visionis threshold**
(the minimum object-above-sun altitude difference at which the object of its computed magnitude
would just become visible) by delegating to `TopoArcVisionis` (see `c-ref-heliacal-vision.md` for
its internals — it runs a bisection search over sun-depression angle against the visual limiting
magnitude model). Writes the threshold value (deg) to `*dret`.

Step by step:
1. `sunra = SunRA(JDNDaysUT, helflag, serr)` — Sun's right ascension at this instant (deg; used
   only to feed `TopoArcVisionis`/extinction downstream, not used directly here).
2. `Magnitude(...)` — apparent magnitude of `ObjectName` at this instant; `ERR` propagates
   immediately.
3. `ObjectLoc(..., Angle=0, ...)` → `AltO`: **topocentric, unrefracted** altitude of the object
   (deg). `ObjectLoc(..., Angle=1, ...)` → `AziO`: topocentric azimuth of the object (deg). (Angle
   codes documented in the call-contract note below.)
4. Moon-relative geometry `AltM`/`AziM`: if `ObjectName` itself **is** "moon" (compared via
   `strncmp(ObjectName, "moon", 4)`), skip the extra lookup and hard-code `AltM = -90`, `AziM = 0`
   — a sentinel meaning "no separate Moon position needed" (the object *is* the Moon, so its
   own alt/az already serve that role downstream in `TopoArcVisionis`). Otherwise, look up the real
   Moon's topocentric altitude/azimuth (`Angle=0`/`Angle=1`) via two more `ObjectLoc` calls — needed
   because `TopoArcVisionis`'s extinction/visibility model accounts for lunar sky-brightening when
   the target is not the Moon itself.
5. `AziS` = Sun's topocentric azimuth (`Angle=1`).
6. Delegate to `TopoArcVisionis(Magn, dobs, AltO, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, dgeo[1]
   (latitude), dgeo[2] (height above sea, m), datm, helflag, dret, serr)`.
7. Any `ERR` from any sub-call short-circuits with an immediate `return ERR` (no partial/garbage
   `*dret` write attempted in that case). On success returns `OK` with the threshold value in
   `*dret`.

Note: this function does **not** re-check `dgeo[2]` altitude bounds itself (that guard lives only
in the public entry points, §4 and `swe_heliacal_ut`) — it is only ever called internally after that
check has already passed.

**Call contracts used above** (see `c-ref-heliacal-vision.md` for internals):
- `SunRA(double JDNDaysUT, int32 helflag, char *serr) -> double` (swehel.c:553) — Sun's RA, degrees.
- `Magnitude(double JDNDaysUT, double *dgeo, char *ObjectName, int32 helflag, double *dmag, char
  *serr) -> int32` (swehel.c:1106) — apparent magnitude of the named object, written to `*dmag`.
- `ObjectLoc(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 Angle, int32
  helflag, double *dret, char *serr) -> int32` (swehel.c:683) — `Angle` selects which quantity is
  written to `*dret`: `0` = topocentric altitude, unrefracted (adds `SEFLG_TOPOCTR`, uses
  `swe_azalt`); `1` = topocentric azimuth (same, then flipped 180° into compass convention); `7` =
  **geocentric** altitude (internally rewritten to the `Angle==0` code path but *without*
  `SEFLG_TOPOCTR`, since the `Angle<5` check that adds that flag runs before the `Angle==7→0`
  rewrite). `ObjectName` is resolved via `DeterObject` (built-in bodies) or `swe_fixstar` fallback
  (unrecognized names, incl. star names) internally.
- `TopoArcVisionis(double Magn, double *dobs, double AltO, double AziO, double AltM, double AziM,
  double JDNDaysUT, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32
  helflag, double *dret, char *serr) -> int32` (swehel.c:1562) — bisection search for the
  sun-depression angle at which the given magnitude becomes visible; writes the arcus-visionis
  threshold (deg) to `*dret`.

## §4 `swe_heliacal_pheno_ut` (swehel.c:1862–2074) — full walkthrough + complete `darr[]` slot table

```c
int32 CALL_CONV swe_heliacal_pheno_ut(double JDNDaysUT, double *dgeo, double *datm, double *dobs,
                                       char *ObjectNameIn, int32 TypeEvent, int32 helflag,
                                       double *darr, char *serr)
```

Parameter units (matching the wider heliacal API, see `c-ref-heliacal-vision.md` for full parameter
docs of `dgeo`/`datm`/`dobs`):
- `dgeo[0..2]` = longitude, latitude (deg), altitude above sea (m).
- `datm[0..3]` = atmospheric pressure (mbar), temperature (°C), relative humidity (%), meteorological
  range VR (km).
- `dobs[0..5]` = observer age (yr), Snellen ratio, binocular flag, optic magnification, optic
  diameter (mm), optic transmission.
- `ObjectNameIn` = object name string (planet/"sun"/"moon"/fixed-star name/numeric asteroid string).
- `TypeEvent` ∈ `{1,2,3,4}` = `SE_MORNING_FIRST`, `SE_EVENING_LAST`, `SE_EVENING_FIRST`,
  `SE_MORNING_LAST` respectively.
- `helflag` = ephemeris + `SE_HELFLAG_*` bits.
- `darr` = output array, caller must allocate ≥ 28 slots (indices `0..27` are written by this
  function on every success path; see slot table below for the two additional slots — `28`,`29` —
  that appear in the **doc-comment header only** and are never actually written).

### Step-by-step

1. **Altitude bound check** (swehel.c:1876–1880): if `dgeo[2] < SEI_ECL_GEOALT_MIN (-500)` or
   `> SEI_ECL_GEOALT_MAX (25000)`, format an error into `serr` and `return ERR` immediately — no
   `darr` slots are written.
2. `swi_set_tid_acc(JDNDaysUT, helflag, 0, serr)` — sets the tidal-acceleration model global (ΔT
   family selection); a stateless Rust port threads this as an explicit parameter/config rather than
   a global (see codebase's existing ΔT handling — this is the same call pattern used by other
   heliacal entry points, e.g. swehel.c:1702, 3392-area).
3. `sunra = SunRA(JDNDaysUT, helflag, serr)` — Sun's RA (deg), used only for the `kt()` extinction
   call in step 8.
4. Name normalization: `strcpy_VBsafe(ObjectName, ObjectNameIn)` (copies only
   alnum/space/`-`/`,` characters, ≤30, into a local buffer — irrelevant to a Rust port using owned
   `String`s; just lowercase/validate as needed) then `tolower_string_star(ObjectName)`
   (lowercases, with fixed-star-specific handling — see `c-ref-heliacal-vision.md`).
5. `default_heliacal_parameters(datm, dgeo, dobs, helflag)` — fills in any zero/unset atmospheric
   and observer defaults in place (ISA pressure/temperature model, default age 36, Snellen 1, etc. —
   see `c-ref-heliacal-vision.md` for the formulas; call contract: mutates `datm`/`dobs` in place,
   reads `dgeo[2]` and `helflag`).
6. `swe_set_topo(dgeo[0], dgeo[1], dgeo[2])` — sets the global topocentric-observer position used
   by subsequent `swe_calc(..., SEFLG_TOPOCTR, ...)` calls inside `ObjectLoc`. **Global state**: a
   stateless Rust port must instead pass topocentric coordinates explicitly into whatever
   `Ephemeris::calc`-equivalent is used downstream (see §5).
7. **Sun and object geometry** (swehel.c:1889–1899), each call short-circuiting on `ERR`:
   - `AziS` = Sun topocentric azimuth (`ObjectLoc(..., "sun", Angle=1, ...)`)
   - `AltS` = Sun topocentric altitude, unrefracted (`Angle=0`)
   - `AziO` = object topocentric azimuth (`Angle=1`)
   - `AltO` = object topocentric altitude, unrefracted (`Angle=0`)
   - `GeoAltO` = object **geocentric** altitude (`Angle=7`)
   - If any of the five returns `ERR`, `return ERR` (no partial `darr` write).
8. Derived instantaneous quantities (swehel.c:1900–1908):
   - `AppAltO = AppAltfromTopoAlt(AltO, datm[1] (temp), datm[0] (pressure), helflag)` — apparent
     (refracted) altitude of the object, via `TopoAltfromAppAlt`'s inverse-Newton iteration (see
     `c-ref-heliacal-vision.md`).
   - `DAZact = AziS - AziO` — azimuth difference (Sun minus object), degrees. **Sign convention**:
     Sun-minus-object, not object-minus-Sun — preserve this order (used again inside `WidthMoon` via
     `AziS - AziO` and in `ARCLact` below).
   - `TAVact = AltO - AltS` — **topocentric** arcus visionis (object altitude minus Sun altitude),
     degrees.
   - `ParO = GeoAltO - AltO` — object's parallax (deg); C comment flags this as "somewhat smaller
     than in Yallop and SkyMap! Needs to be studied" — a known, acknowledged-in-source
     approximation, not a bug to silently "improve" during porting.
   - `Magnitude(JDNDaysUT, dgeo, ObjectName, helflag, &MagnO, serr)` — apparent magnitude; `ERR`
     propagates.
   - `ARCVact = TAVact + ParO` — **geocentric** arcus visionis, degrees (parallax-corrected).
   - `ARCLact = acos(cos(ARCVact * DEGTORAD) * cos(DAZact * DEGTORAD)) / DEGTORAD` — longitude
     difference between object and Sun (great-circle-style combination of the two orthogonal
     angular differences), degrees.
9. **Elongation & illumination** (swehel.c:1909–1918):
   - `Planet = DeterObject(ObjectName)` — returns a body-index constant (`SE_SUN=0`, `SE_MOON=1`,
     …) for recognized planet names or a numeric-asteroid string (→ `atoi(s) + SE_AST_OFFSET`), or
     `-1` for anything else (fixed-star names).
   - If `Planet == -1` (fixed star): `elong = ARCLact` (reuse the longitude-difference already
     computed — stars have no phase), `illum = 100` (assume fully "illuminated"/point source).
   - Else: `iflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)` (ephemeris-selection bits
     only, computed once at function entry, swehel.c:1875) is passed to
     `swe_pheno_ut(JDNDaysUT, Planet, iflag|(SEFLG_TOPOCTR|SEFLG_EQUATORIAL), attr, serr)` (see
     `docs/c-ref-phenomena.md` for `swe_pheno`/`swe_pheno_ut` internals); `ERR` propagates; then
     `elong = attr[2]` (elongation), `illum = attr[1] * 100` (phase fraction → percent). Note this
     branch runs even when `Planet == SE_SUN` (i.e. `ObjectName` is literally "sun") — no special
     case excludes it; `swe_pheno_ut`'s own Sun-handling (elongation stays 0) applies transparently.
10. **Extinction coefficient** (swehel.c:1919): `kact = kt(AltS, sunra, dgeo[1] (lat), dgeo[2]
    (height), datm[1] (temp), datm[2] (RH), datm[3] (VR), 4 /* ExtType */, serr)` — see
    `c-ref-heliacal-vision.md` for `kt`'s internals; `ExtType=4` selects whichever fixed extinction
    model variant that constant denotes there.
11. **Dead code block** (swehel.c:1920–1926): `if ((0)) { darr[26]=kR(...); darr[27]=kW(...);
    darr[28]=kOZ(...); darr[29]=ka(...); darr[30]=darr[26]+darr[27]+darr[28]+darr[29]; }` — compiled
    out unconditionally (`if ((0))`), never executes. **Do not port.** Note this would write
    `darr[26]`/`darr[27]` with extinction-component breakdowns (Rayleigh/water-vapor/ozone/aerosol
    `k` terms + their sum in `darr[30]`, which is out-of-bounds for a 30-element array!) if it were
    ever re-enabled — but since it's dead, `darr[26]`/`darr[27]` are unconditionally overwritten
    later (step 18) with `elong`/`illum` instead, and `darr[28..30]` are simply never touched by
    this function. Flagged here only so a future re-enable of this block isn't mistaken for new
    functionality to port.
12. **Moon-only Yallop crescent block** (swehel.c:1927–1941): `WMoon = qYal = qCrit = LMoon = 0`
    unconditionally first; then only `if (Planet == SE_MOON)`:
    - `WMoon = WidthMoon(AltO, AziO, AltS, AziS, ParO)` (§1)
    - `LMoon = LengthMoon(WMoon, 0)` (§1 — note the literal `0` argument, always triggering
      `LengthMoon`'s `AvgRadiusMoon`-fallback path)
    - `qYal = qYallop(WMoon, ARCVact)` (§1)
    - `qCrit` grade assignment per the six-branch table in §1.
13. **Rise/set of Sun and object** (swehel.c:1942–1959):
    - `RS = 2` by default; `RS = 1` if `TypeEvent == 1` (morning first) or `TypeEvent == 4` (morning
      last) — i.e. `RS=1` means "a rising/morning-oriented search", `RS=2` means
      "setting/evening-oriented".
    - `RiseSet(JDNDaysUT - 4.0/24.0, dgeo, datm, "sun", RS, helflag, 0 /* Rim */, &RiseSetS, serr)`
      — Sun rise/set time near `JDNDaysUT`, searched from 4 hours earlier (a fixed lookback so the
      search always finds the *relevant* rise/set even if `JDNDaysUT` itself is already past it);
      `ERR` propagates. `RS` selects which event (rise vs set) `RiseSet` looks for — see
      `c-ref-heliacal-vision.md`.
    - Same call for `ObjectName` → `RiseSetO`.
    - `TbYallop = TJD_INVALID` initially.
    - If the object's `RiseSet` call returned the special code `-2` ("object does not rise or set"
      at this latitude/date — e.g. circumpolar or never-rises): `Lag = 0`, `noriseO = TRUE`.
    - Else: `Lag = RiseSetO - RiseSetS` (days; object rise/set minus Sun rise/set); if
      `Planet == SE_MOON`: `TbYallop = (RiseSetO * 4 + RiseSetS * 5) / 9.0` — Yallop's empirical
      weighted-average "best observation time" (4:5 weighting toward the Sun's event), only ever
      computed for the Moon.
14. **Early-exit guard** (swehel.c:1960–1967): if `(TypeEvent == 3 || TypeEvent == 4)` **and**
    `(Planet == -1 || Planet >= SE_MARS)` — i.e. an evening-first/morning-last search for a fixed
    star, or for Mars/Jupiter/Saturn/Uranus/Neptune/any numbered asteroid — set `TfirstVR =
    TbVR = TlastVR = TJD_INVALID`, `TvisVR = 0`, `MinTAV = 0`, and `goto output_heliacal_pheno`,
    **skipping the entire visibility-search loop below** (steps 15–17). Only Sun/Moon/Mercury/Venus
    reach the loop for these two `TypeEvent`s; all bodies reach the loop for `TypeEvent` 1/2.
15. **Visibility-window search loop** (swehel.c:1968–2042) — only reached if step 14's guard didn't
    fire. Purpose: scan time in small steps from the Sun's rise/set instant, tracking (a) the
    instantaneous "arcus visionis visibility threshold" `MinTAVact` (via `DeterTAV`, §3, which
    itself calls `TopoArcVisionis`'s bisection search — i.e. this is an outer scan wrapping an inner
    bisection at every sample) and (b) the instantaneous altitude difference `DeltaAlt = AltO2 -
    AltS2`, to find: the best/optimal visibility time `TbVR` (local minimum of the threshold curve,
    parabola-refined), and the first/last crossing times `Tc`/`Ta` where `DeltaAlt` crosses above/
    below the threshold (linearly interpolated).

    Initialization (swehel.c:1970–1978):
    ```
    MinTAVact = 199;  DeltaAlt = 0;  OldestMinTAV = 0;  Ta = 0;  Tc = 0;  TbVR = 0;
    TimeStep = -TimeStepDefault / 24.0 / 60.0;      /* -1 minute, as a fraction of a day */
    if (RS == 2) TimeStep = -TimeStep;               /* flip to +1 minute for RS==2 */
    TimePointer = RiseSetS - TimeStep;
    ```
    So for `RS==1` (rising search) `TimeStep` is negative (stepping **backward** in time from
    `RiseSetS`); for `RS==2` (setting search) it's positive (stepping **forward**). `TimePointer`
    is pre-offset by one step so the first loop iteration's `TimePointer += TimeStep` lands exactly
    on `RiseSetS`.

    Loop body (`do { ... } while (...)`, swehel.c:1979–2017), each iteration:
    ```
    TimePointer += TimeStep;
    OldestMinTAV = MinTAVoud;  MinTAVoud = MinTAVact;  DeltaAltoud = DeltaAlt;   /* shift history */
    AltS2 = ObjectLoc(TimePointer, ..., "sun", Angle=0, ...);      /* topocentric, at TimePointer */
    AltO2 = ObjectLoc(TimePointer, ..., ObjectName, Angle=0, ...);
    DeltaAlt = AltO2 - AltS2;
    DeterTAV(dobs, TimePointer, dgeo, datm, ObjectName, helflag, &MinTAVact, serr);   /* §3 */
    ```
    Any `ObjectLoc`/`DeterTAV` failure → immediate `return ERR`.

    **Local-minimum detection** (swehel.c:1992–2008): if `MinTAVoud < MinTAVact` (the threshold
    curve just started increasing again, i.e. we may have just passed a minimum) **and** `TbVR ==
    0` (haven't found one yet):
    - Look ahead (`TimeCheck = TimePointer + Sgn(TimeStep) * LocalMinStep/24/60`, i.e. one more
      `LocalMinStep` (8) minutes further in the *same* direction the scan is already moving —
      `Sgn(x)` returns `-1` if `x<0` else `1`, i.e. `Sgn(0) = +1`, swehel.c:864–869), clamped against
      `RiseSetO` if the object does rise/set (`TimeCheck = min(TimeCheck, RiseSetO)` if
      `TimeStep>0`, else `max(...)`) so the look-ahead never probes past the object's own
      rise/set.
    - `DeterTAV` at `TimeCheck` → `LocalminCheck`. If `LocalminCheck > MinTAVact` (confirms the
      minimum really is a minimum — the threshold is still higher 8 minutes further out — object
      still above horizon at that check point implicitly via the clamp above): refine via the
      3-point parabola fit (§2): `extrax = x2min(MinTAVact, MinTAVoud, OldestMinTAV)`, `TbVR =
      TimePointer - (1 - extrax) * TimeStep`, `MinTAV = funct2(MinTAVact, MinTAVoud, OldestMinTAV,
      extrax)`.

    **Visibility-start crossing** (swehel.c:2009–2012): if `DeltaAlt > MinTAVact` (altitude
    difference now exceeds the visibility threshold — object visible) **and** `Tc == 0` **and**
    `TbVR == 0` (only look for this crossing before the best-time minimum has been found):
    `crosspoint = crossing(DeltaAltoud, DeltaAlt, MinTAVoud, MinTAVact)` (§2), `Tc = TimePointer -
    TimeStep * (1 - crosspoint)`.

    **Visibility-end crossing** (swehel.c:2013–2016): if `DeltaAlt < MinTAVact` (dropped back below
    threshold) **and** `Ta == 0` **and** `Tc != 0` (only after a start-crossing was already found):
    same `crossing()` call shape, `Ta = TimePointer - TimeStep * (1 - crosspoint)`.

    **Loop termination** (swehel.c:2017):
    ```
    while ( fabs(TimePointer - RiseSetS) <= MaxTryHours/24.0     /* within 4h of Sun rise/set */
            && Ta == 0                                            /* haven't found the end yet */
            && !( TbVR != 0 && (TypeEvent==3 || TypeEvent==4)
                  && ObjectName is not "moon"/"venus"/"mercur[y]" ) );
    ```
    The third clause is an early-out: once the best-visibility time `TbVR` has been found, for
    `TypeEvent` 3/4 (evening-first/morning-last) searches on any object *other than* Moon/Venus/
    Mercury, stop scanning immediately (don't bother looking for `Ta`). Given step 14 already routed
    stars/Mars-and-beyond away from this loop entirely for `TypeEvent` 3/4, in practice this clause
    can only fire here for `ObjectName == "sun"` under `TypeEvent` 3/4 (a degenerate/unlikely
    real-world query, but the logic handles it).
16. **Post-loop assembly** (swehel.c:2018–2042):
    - `RS == 2`: `TfirstVR = Tc`, `TlastVR = Ta`. Else (`RS == 1`): `TfirstVR = Ta`, `TlastVR = Tc`.
    - If both `TfirstVR == 0 && TlastVR == 0` (loop never found *either* crossing — e.g. object was
      already above/below threshold the entire scanned window): fall back to a point estimate around
      `TbVR`: if `RS == 1`, `TfirstVR = TbVR - 0.000001`; else `TlastVR = TbVR + 0.000001` (a
      ~0.0864-second nudge off the best-time instant, used purely so `TfirstVR`/`TlastVR` isn't left
      at the `0` sentinel).
    - If `!noriseO` (object does properly rise/set): clamp against the object's own rise/set time —
      `RS==1`: `TfirstVR = max(TfirstVR, RiseSetO)` (visibility can't start before the object rises);
      `RS==2`: `TlastVR = min(TlastVR, RiseSetO)` (visibility can't outlast the object's set).
    - `TvisVR = TJD_INVALID` initially; if both `TlastVR != 0 && TfirstVR != 0`: `TvisVR = TlastVR -
      TfirstVR` (visibility-window duration, days).
    - Any of `TlastVR`, `TbVR`, `TfirstVR` still `== 0` at this point get remapped to `TJD_INVALID`
      (the `0` sentinel used internally during the search becomes the public "not applicable"
      sentinel on output).
17. `output_heliacal_pheno:` label (swehel.c:2043) — the `goto` target from step 14 also lands here,
    skipping straight to the `darr[]` writes below with whatever values steps 1–13 (but not 15–16)
    established (`MinTAV=0`/`TfirstVR=TbVR=TlastVR=TJD_INVALID`/`TvisVR=0` from step 14 in that
    case).
18. **`darr[]` output writes** (swehel.c:2045–2073) — see the complete slot table below. Function
    returns `OK` (swehel.c:2073) on every path that reaches here (the only `ERR` returns are the
    early ones threaded through steps 1, 7–10, 15).

### Complete `darr[]` slot table

The function's own doc-comment header (swehel.c:1826–1861) enumerates slots `0` through `29`
(30 conceptual slots), but the function body only ever **writes** `darr[0]` through `darr[27]`
(28 slots) — see the "Written?" column below. `darr[28]` and `darr[29]` are documented in the header
comment (as `CVAact [deg]` and `MSk [-]`, apparently planned/aspirational additions) but are dead
labels with no corresponding write in this function; a Rust port's output type should likewise
expose only indices `0..27`, or explicitly mark `28`/`29` as unimplemented/`None` if mirroring the
C array shape for compatibility.

| idx | Name (per C header) | Units | Written? | Computed as (see step above) |
|---|---|---|---|---|
| 0 | `AltO` | deg | yes | topocentric altitude of object, unrefracted (step 7) |
| 1 | `AppAltO` | deg | yes | apparent (refracted) altitude of object (step 8, `AppAltfromTopoAlt`) |
| 2 | `GeoAltO` | deg | yes | geocentric altitude of object (step 7, `Angle=7`) |
| 3 | `AziO` | deg | yes | azimuth of object (step 7, `Angle=1`) |
| 4 | `AltS` | deg | yes | topocentric altitude of Sun (step 7) |
| 5 | `AziS` | deg | yes | azimuth of Sun (step 7) |
| 6 | `TAVact` | deg | yes | actual topocentric arcus visionis = `AltO - AltS` (step 8) |
| 7 | `ARCVact` | deg | yes | actual geocentric arcus visionis = `TAVact + ParO` (step 8) |
| 8 | `DAZact` | deg | yes | azimuth difference `AziS - AziO` (step 8) |
| 9 | `ARCLact` | deg | yes | longitude difference, `acos(cos(ARCVact°)·cos(DAZact°))°` (step 8) |
| 10 | `kact` | (extinction coeff., mag/airmass) | yes | `kt(...)` (step 10) |
| 11 | `MinTAV` | deg | yes | smallest topocentric arcus visionis (parabola-refined local min, step 15; `0` if guard fired at step 14 or no local min found) |
| 12 | `TfistVR` (i.e. `TfirstVR`) | JDN | yes | first time object is visible, per VR search (steps 15–16; `TJD_INVALID` if guard fired) |
| 13 | `TbVR` | JDN | yes | optimum visibility time, per VR search (steps 15–16) |
| 14 | `TlastVR` | JDN | yes | last time object is visible, per VR search (steps 15–16) |
| 15 | `TbYallop` | JDN | yes | Yallop's weighted best-time estimate, Moon only (step 13; `TJD_INVALID` for non-Moon or if object never rises/sets) |
| 16 | `WMoon` | deg | yes | crescent width, Moon only (step 12; `0` otherwise) |
| 17 | `qYal` | (dimensionless) | yes | Yallop q-test value, Moon only (step 12; `0` otherwise) |
| 18 | `qCrit` | (dimensionless, 1–6 / A–F) | yes | Yallop q-test grade, Moon only (step 12; `0` otherwise) |
| 19 | `ParO` | deg | yes | parallax of object = `GeoAltO - AltO` (step 8) |
| 20 | `Magn` | (stellar magnitude) | yes | apparent magnitude of object (step 8, `MagnO`) |
| 21 | `RiseO` | JDN | yes | rise/set time of object (step 13, `RiseSetO`) |
| 22 | `RiseS` | JDN | yes | rise/set time of Sun (step 13, `RiseSetS`) |
| 23 | `Lag` | JDN (i.e. days) | yes | `RiseSetO - RiseSetS` (step 13; `0` if object doesn't rise/set) |
| 24 | `TvisVR` | JDN (days) | yes | visibility duration `TlastVR - TfirstVR` (step 16; `TJD_INVALID` if either endpoint unresolved, `0` if guard fired at step 14) |
| 25 | `LMoon` | deg | yes | crescent length, Moon only (step 12; `0` otherwise) |
| 26 | `CVAact` (header comment name) | deg (header) | yes, but **not** with the header's claimed meaning | actually written with `elong` (step 9) — the header comment's "CVAact" label for slot 26 does not match what the code writes here; trust the code (`darr[26] = elong;`, swehel.c:2071), not the header comment, for this slot |
| 27 | `Illum` (header: `Illum [%] 'new'`) | percent | yes | `illum` (step 9) — `attr[1]*100` for planets, `100` for fixed stars |
| 28 | `CVAact 'new'` (header) | deg (header) | **no** | never written by this function body; header-only, dead label |
| 29 | `MSk` (header) | (dimensionless, header) | **no** | never written by this function body; header-only, dead label |

The header comment literally shows both slot 26 and slot 28 labeled `CVAact [deg]` (with `'new'`
suffixed on 28), which is very likely leftover/stale documentation from an earlier revision of this
function where slot semantics were being reshuffled — the live code's actual slot-26 write is
`elong`, and the plausible original intent for a `CVAact` ("current visibility arc actual"?) slot
was apparently superseded and left undone. Do not attempt to "restore" a `CVAact` computation for
slot 26 or 28 in the Rust port — port exactly what the C code writes (`elong` at 26, `illum` at 27,
nothing at 28/29).

## §5 Porting notes for the stateless Rust port

**Global/mutable state touched by this range, and how a stateless port should handle it:**

1. `swi_set_tid_acc(JDNDaysUT, helflag, 0, serr)` (step 2) — selects the tidal-acceleration/ΔT
   model as global state in C. The existing Rust ΔT handling elsewhere in this codebase already
   threads this as explicit config (see wherever `swi_set_tid_acc`'s Rust analog already lives, e.g.
   near `deltat`/`EphemerisConfig`) — reuse that, don't reintroduce a global.
2. `swe_set_topo(dgeo[0], dgeo[1], dgeo[2])` (step 6) — sets the global topocentric-observer
   position subsequently read by every `swe_calc(..., SEFLG_TOPOCTR, ...)` call inside `ObjectLoc`
   (§3's call contract) for the rest of this function's execution (including every iteration of the
   step-15 search loop, and every nested `DeterTAV`/`TopoArcVisionis` call). A stateless Rust port
   must pass the observer geo-coordinates explicitly into whatever `Ephemeris::calc`-equivalent
   `ObjectLoc`'s Rust counterpart uses, on every call, rather than relying on a prior "set" call.
   This is the same category of issue as `docs/c-ref-phenomena.md`'s notes on `swed.ast_diam`/
   `ast_H`/`ast_G` — check whether the Rust `ObjectLoc`-equivalent (in the sibling
   heliacal-vision port) already threads `dgeo` explicitly per-call before assuming this needs new
   plumbing.
3. No other `swed.*`/global reads or writes occur in this line range beyond the two above and the
   ΔT-family calls buried inside `SunRA`/`ObjectLoc`/`Magnitude`/`RiseSet`/`DeterTAV` (documented in
   `c-ref-heliacal-vision.md`).

**Dead code in this range, not to be ported:**
- swehel.c:1920–1926 — the `if ((0)) { darr[26..30] = kR/kW/kOZ/ka/(sum) }` block (§4 step 11). Note
  it also writes an out-of-bounds `darr[30]` if a 30-element (`0..29`) array were passed — another
  reason this must stay dead.
- The `darr[28]`/`darr[29]` header-comment labels (`CVAact 'new'`, `MSk`) with no corresponding
  code — nothing to port, just don't invent an implementation for them (§4 slot table note).

**Gotchas / easy-to-miss details when porting:**
- `qCrit`'s six-band grading (§1) uses **strict** `<`/`>` on both sides of each band, so the five
  exact boundary values (`0.216, -0.014, -0.16, -0.232, -0.293`) fall through to `qCrit = 0` (no
  grade). Don't silently convert to an exhaustive `match` with closed intervals — preserve the
  gap.
- `LengthMoon` is **always** called with its second argument `0` from `swe_heliacal_pheno_ut` (§4
  step 12), so the "real diameter" parameter path inside `LengthMoon` is dead on this call site —
  the Rust port of the pheno function only ever needs the `AvgRadiusMoon`-fallback branch, though
  `LengthMoon` itself (§1) should still be ported generally (it may be called with a real diameter
  from other, not-yet-ported call sites elsewhere in `swehel.c`).
- `WidthMoon`'s `parallax` parameter is fed `ParO` (`GeoAltO - AltO`, computed in step 8 the normal
  way via two independent `ObjectLoc` calls), **not** anything returned by `swe_pheno`'s
  Moon-parallax output (`attr[5]` in `docs/c-ref-phenomena.md`) — these are two independently
  computed parallax values in the C codebase and are not guaranteed to agree to the last bit; use
  `ParO` here to match this function's behavior, don't "simplify" by reusing a
  `phenomena.rs`-computed parallax even though it's semantically the same quantity.
- The step-15 search loop's early-termination clause (§4 step 15, loop condition) is easy to
  mistranslate: it depends on the **literal object-name string**, checked via `strncmp` against
  `"moon"`/`"venus"`/`"mercur"`-prefixes (note: `"mercury"` in the C `while` condition is checked
  with `strncmp(ObjectName, "mercury", 7)`, 7 chars, matching the full word — distinct from
  `DeterObject`'s own `"mercur"`/6-char prefix match used elsewhere). A Rust port keyed on a
  `Body` enum should reconstruct the equivalent body-identity check (`Body::Moon | Body::Venus |
  Body::Mercury`) rather than re-implementing string prefix matching, but must confirm the
  practical effect is identical (it is, given `DeterObject` already normalized the name upstream).
- `Sgn(0.0) == 1`, not `0` or `-1` (swehel.c:864–869: `if (x<0) return -1; return 1;`). `TimeStep`
  is never actually `0` at the one call site (§4 step 15's local-min look-ahead), so this only
  matters for defensive/edge-case fidelity, not a live bug — but don't reach for a idiomatic
  three-way `signum()` that returns `0` for `0.0` without checking this.
- `RiseSet`'s three-way return convention (`OK` / `ERR` / the special sentinel `-2` meaning "object
  does not rise or set") must be preserved as a tri-state result in the Rust port (e.g. an enum or
  `Result<Option<f64>, Error>` shape) — §4 step 13 branches specifically on the `-2` case
  (`noriseO = TRUE`, `Lag = 0`) distinctly from both success and hard error.
- `swe_pheno_ut`'s `attr` buffer is sized `[30]` locally in this function (`double attr[30]`,
  swehel.c:1869) even though only `attr[0..5]` are ever written by `swe_pheno`/`swe_pheno_ut` (see
  `docs/c-ref-phenomena.md`) — only `attr[1]` and `attr[2]` are actually read back here (step 9).

**Call contracts into part-1 (heliacal-vision) functions used in this range** (signatures only —
see `c-ref-heliacal-vision.md` for algorithm internals):
- `SunRA(double JDNDaysUT, int32 helflag, char *serr) -> double` — swehel.c:553
- `ObjectLoc(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 Angle, int32 helflag, double *dret, char *serr) -> int32` — swehel.c:683; `Angle` codes `0`=topocentric alt (unrefracted), `1`=topocentric az, `7`=geocentric alt (see §3)
- `Magnitude(double JDNDaysUT, double *dgeo, char *ObjectName, int32 helflag, double *dmag, char *serr) -> int32` — swehel.c:1106
- `TopoArcVisionis(double Magn, double *dobs, double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, double *dret, char *serr) -> int32` — swehel.c:1562
- `RiseSet(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 RSEvent, int32 helflag, int32 Rim, double *tret, char *serr) -> int32` — swehel.c:535; tri-state return (`OK`/`ERR`/`-2`)
- `AppAltfromTopoAlt(double TopoAlt, double TempE, double PresE, int32 helflag) -> double` — swehel.c:626
- `default_heliacal_parameters(double *datm, double *dgeo, double *dobs, int helflag)` (void, mutates in place) — swehel.c:1324
- `DeterObject(char *ObjectName) -> int32` — swehel.c:305; returns a body-index constant, or `atoi(name)+SE_AST_OFFSET` for numeric strings, or `-1` for unrecognized/star names
- `kt(double AltS, double sunra, double Lat, double HeightEye, double TempS, double RH, double VR, int32 ExtType, char *serr) -> double` — swehel.c:940
- `swe_pheno_ut` — see `docs/c-ref-phenomena.md` in full
