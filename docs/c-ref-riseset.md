# C Reference: Rise / Set / Meridian Transit — swecl.c

Porting reference for `swe_rise_trans` / `swe_rise_trans_true_hor` and the shared parabola
helpers they use. All line numbers below refer to `swecl.c` unless stated otherwise.

`find_maximum`/`find_zero` are also used throughout the eclipse code (contact-time refinement
in `swe_sol_eclipse_when_loc`, `swe_lun_eclipse_when_loc`, `swe_sol_eclipse_when_glob`,
`swe_lun_eclipse_when`, occultation search, etc. — lines 1298, 1399, 1715, 1856, 2191, 2254,
2309, 2531, 2605, 2662, 3511, 3579). They should be ported once into a shared math module and
reused by both the rise/set module and the (future) eclipse module — do not duplicate.

---

## Function Map

| C function | Location | Purpose |
|---|---|---|
| `find_maximum` | swecl.c:4133–4146 | Parabola-fit extremum (offset + value) through 3 equally-spaced samples. Shared helper. |
| `find_zero` | swecl.c:4148–4162 | Parabola-fit root(s) through 3 equally-spaced samples. Shared helper. |
| `rdi_twilight` | swecl.c:4164–4174 | rsmi twilight bits → target altitude depression (6/12/18°). |
| `get_sun_rad_plus_refr` | swecl.c:4176–4194 | Disc radius (± bottom-limb sign) + refraction, for `rise_set_fast` only. |
| `rise_set_fast` (static) | swecl.c:4203–4325 | Fast algorithm: semi-diurnal-arc estimate + 2–4 Newton iterations. Only for `\|lat\|` ≤ 60/65°, no twilight, no fixstar. |
| `swe_rise_trans` | swecl.c:4355–4383 | Public dispatcher: fast path vs `swe_rise_trans_true_hor(horhgt=0)`. |
| `swe_rise_trans_true_hor` | swecl.c:4387–4686 | Full algorithm: 15-point altitude sampling + culmination insertion + bisection zero-crossing. Handles horizon-dip (`horhgt`). |
| `calc_mer_trans` (static) | swecl.c:4688–4748 | Meridian/anti-meridian transit: 4-iteration RA-vs-armc Newton-like search. |

---

## 1. Shared Parabola Helpers

Both helpers fit a parabola `y(x) = a·x² + b·x + c` through three samples taken at
**equally spaced** abscissas `x ∈ {-1, 0, +1}`, corresponding to values `y00` (at `x=-1`),
`y11` (at `x=0`, the *middle* sample), `y2` (at `x=+1`, the *most recent / rightmost* sample),
with real sample spacing `dx` (i.e. real time = `x·dx` relative to the `y11` sample, or
equivalently `(x-1)·dx` relative to the `y2` sample — see return convention below).

Interpolation coefficients (identical in both functions):
```c
c = y11;
b = (y2 - y00) / 2.0;
a = (y2 + y00) / 2.0 - c;
```
Derivation: solving `a - b + c = y00`, `c = y11`, `a + b + c = y2` for `a, b, c` gives exactly
these three lines. `c` is the value at the middle sample; `b` is half the total rise across the
outer two samples; `a` is the curvature term (residual after removing the linear trend and the
mid value).

### `find_maximum(y00, y11, y2, dx, *dxret, *yret)` — swecl.c:4133–4146

```c
static int find_maximum(double y00, double y11, double y2, double dx,
                        double *dxret, double *yret)
{
  double a, b, c, x, y;
  c = y11;
  b = (y2 - y00) / 2.0;
  a = (y2 + y00) / 2.0 - c;
  x = -b / 2 / a;
  y = (4 * a * c - b * b) / 4 / a;
  *dxret = (x - 1) * dx;
  if (yret != NULL)
    *yret = y;
  return OK;
}
```
Algebra: vertex of `a·x²+b·x+c` is at `dy/dx = 2ax+b = 0 → x = -b/(2a)`; value there is
`y = c - b²/(4a) = (4ac-b²)/(4a)` (both exactly as coded — no branch, so **`a == 0` is a
guaranteed div-by-zero**; callers never hit this because samples always come from real
altitude curves with curvature). No domain/convergence check — this is a single closed-form
evaluation, not an iteration.

**Return convention:** `*dxret = (x - 1) * dx`, i.e. the extremum's time offset is reported
**relative to the `y2` sample** (the third/most-recent one passed in), not relative to `y11`
(the middle one) and not an absolute `x`-space value. Callers exploit this directly: they set a
reference time `t` = the time of the `y2` sample, then do `t_extremum = t + *dxret` (see
`swe_rise_trans_true_hor` culmination refinement below, which literally computes
`tcu = (t - dt); tcu += *dxret + dt;` — algebraically `tcu = t + *dxret`, the `-dt`/`+dt` pair
cancels and exists only for code-historical reasons; replicate the identical arithmetic for FP
fidelity, don't simplify to `t + dxret` directly).

Always returns `OK` (return value is vestigial / unused by any caller — safe to make the Rust
port return the offset+value directly without a `Result`).

### `find_zero(y00, y11, y2, dx, *dxret1, *dxret2)` — swecl.c:4148–4162

```c
static int find_zero(double y00, double y11, double y2, double dx,
                        double *dxret, double *dxret2)
{
  double a, b, c, x1, x2;
  c = y11;
  b = (y2 - y00) / 2.0;
  a = (y2 + y00) / 2.0 - c;
  if (b * b - 4 * a * c < 0)
    return ERR;
  x1 = (-b + sqrt(b * b - 4 * a * c)) / 2 / a;
  x2 = (-b - sqrt(b * b - 4 * a * c)) / 2 / a;
  *dxret = (x1 - 1) * dx;
  *dxret2 = (x2 - 1) * dx;
  return OK;
}
```
Standard quadratic formula on the same `a,b,c`. Returns `ERR` if the discriminant is negative
(parabola never crosses zero within the fitted window — no real root). No `a == 0` guard here
either (div-by-zero if `a` is exactly 0, same as `find_maximum`). Both roots are always
computed and returned (order of `x1`/`x2` is fixed: `x1` uses `+sqrt`, `x2` uses `-sqrt`; the
caller picks whichever lies in-range, or ignores the unused root).

**Return convention:** same as `find_maximum` — both `*dxret`/`*dxret2` are `(x_i - 1) * dx`,
offsets relative to the `y2` (rightmost) sample.

**FP hazard:** `b*b - 4*a*c` is computed twice (once for the `<0` guard, once inside each
`sqrt(...)`) — replicate the duplicate computation rather than caching it in a local, to match
FP rounding exactly if the golden tests are bitwise-exact at this level (they generally aren't,
since callers apply epsilon tolerances to time-of-event results, but keep the structure anyway
for clarity/traceability).

---

## 2. `rsmi` Flag Bits (swephexp.h:335–361)

```c
#define SE_CALC_RISE               1
#define SE_CALC_SET                2
#define SE_CALC_MTRANSIT           4
#define SE_CALC_ITRANSIT           8
#define SE_BIT_DISC_CENTER       256   /* rise/set of disc CENTER instead of upper limb */
#define SE_BIT_DISC_BOTTOM      8192   /* rise/set of disc BOTTOM (lower limb) */
#define SE_BIT_GEOCTR_NO_ECL_LAT 128   /* use geocentric position, ignore ecliptic latitude */
#define SE_BIT_NO_REFRACTION     512   /* ignore atmospheric refraction */
#define SE_BIT_CIVIL_TWILIGHT   1024   /* civil twilight   (target altitude -6°)  */
#define SE_BIT_NAUTIC_TWILIGHT  2048   /* nautical twilight (target altitude -12°) */
#define SE_BIT_ASTRO_TWILIGHT   4096   /* astronomical twilight (target altitude -18°) */
#define SE_BIT_FIXED_DISC_SIZE 16384   /* neglect the effect of distance on disc size
                                        * (uses a fixed nominal distance: Sun=1 AU,
                                        * Moon=0.00257 AU, for the disc-radius calc only) */
#define SE_BIT_FORCE_SLOW_METHOD 32768 /* Astrodienst in-house test flag: forces the full/slow
                                        * algorithm even when the fast path would qualify */
#define SE_BIT_HINDU_RISING  (SE_BIT_DISC_CENTER | SE_BIT_NO_REFRACTION | SE_BIT_GEOCTR_NO_ECL_LAT)
```
`SE_CALC_RISE`/`SE_CALC_SET`/`SE_CALC_MTRANSIT`/`SE_CALC_ITRANSIT` are mutually exclusive event
selectors (low 4 bits); everything else is an OR-able modifier on `SE_CALC_RISE`/`SE_CALC_SET`.
`MTRANSIT` = upper culmination (meridian transit, body crosses the local meridian above the
pole), `ITRANSIT` = lower culmination / "anti-transit" (body crosses the meridian on the far
side, below the pole — historically named after the Latin *infra-*/*imum coeli* transit).

Twilight bits are meaningful only for `ipl == SE_SUN`; setting one forces
`SE_BIT_NO_REFRACTION | SE_BIT_DISC_CENTER` and a negative `horhgt` (see §5.5).

---

## 3. `rise_set_fast` — swecl.c:4203–4325

```c
static int32 rise_set_fast(double tjd_ut, int32 ipl, int32 epheflag, int32 rsmi,
                            double *dgeo, double atpress, double attemp,
                            double *tret, char *serr)
```
Applies to Sun, Moon, planets, nodes/apsides (`ipl` in `SE_SUN..SE_TRUE_NODE`) — never to fixed
stars (dispatcher excludes them, see §4). Comment: doesn't work well above 65°N/S for the Sun,
60°N/S for Moon/planets — hence the dispatcher's latitude gate.

### 3.1 Setup
- `nloop = 2` Newton iterations, except `nloop = 4` if `ipl == SE_MOON` (Moon's faster motion
  needs more refinement).
- `facrise = 1` for rise, `facrise = -1` if `rsmi & SE_CALC_SET`.
- Unless `SE_BIT_GEOCTR_NO_ECL_LAT`: `swe_set_topo(dgeo[0], dgeo[1], dgeo[2])` (STATEFUL — see
  STATELESS PORT NOTE below) and `iflagtopo |= SEFLG_TOPOCTR`. `iflagtopo` always carries
  `SEFLG_EQUATORIAL` initially (right ascension/declination working frame).

### 3.2 Semi-diurnal-arc estimate (`run_rise_again:` label, swecl.c:4232–4272)
1. `swe_calc_ut(tjd_ut, ipl, iflagtopo, xx, serr)` → `decl = xx[1]` (declination at the input
   time — a single snapshot, not iterated on declination drift within this step).
2. Semi-diurnal arc: `sda = -tan(geolat)·tan(decl)` (radians); then:
   - `sda >= 1` → circumpolar-never-rises case: clamp `sda = 10` (not `0` — a literal `0` would
     make `mdrise` degenerate and break the meridian-distance math below; `10°` is an ad-hoc
     placeholder, "to account for refraction").
   - `sda <= -1` → circumpolar-never-sets case: clamp `sda = 180`.
   - else: `sda = acos(sda)` in degrees.
   **This function never returns a circumpolar/`-2` signal** — the clamps above just produce a
   (possibly meaningless) `tr` estimate that the Newton loop then "refines"; only the *full*
   algorithm (`swe_rise_trans_true_hor`) can return `-2`.
3. `armc = swe_degnorm(swe_sidtime(tjd_ut)*15 + dgeo[0])` (local sidereal time in degrees, at
   `tjd_ut` — NOT the Newton-iteration time `tr`; this whole armc/md/mdrise/dmd block is
   evaluated once, at the very start).
4. `md = swe_degnorm(xx[0] - armc)` (meridian distance of the object's RA from local sidereal
   time — hour-angle-like quantity, in degrees).
5. `mdrise = swe_degnorm(sda * facrise)`.
6. `dmd = swe_degnorm(md - mdrise)`; if `dmd > 358`, subtract 360 (guards against "next-day"
   overshoot when the true event is actually in the very near past, within ~2° of `md`).
7. `tr = tjd_ut + dmd/360` — the rough rise/set time estimate (semi-diurnal arc converted from
   degrees of sidereal rotation to a day fraction via the crude `/360`, not the true sidereal
   rate `ARMCS≈360.9856` — acceptable since this is only a Newton starting point).

### 3.3 Newton iteration (swecl.c:4273–4315)
1. Estimate/compute atmospheric pressure if `atpress == 0`:
   `atpress = 1013.25 * (1 - 0.0065·height/288)^5.255` (barometric formula, height in meters).
2. `refr = swe_refrac_extended(0.000001, 0, atpress, attemp, const_lapse_rate, SE_APP_TO_TRUE, xx); refr = xx[1]-xx[0]`
   — refraction at the horizon (apparent altitude ≈ 0), computed **once**, reused every
   iteration (not recomputed at each `tr`).
3. Re-select working flags for the horizontal-conversion path: `SE_BIT_GEOCTR_NO_ECL_LAT` →
   `tohor_flag = SE_ECL2HOR`, `iflagtopo = iflag` (geocentric ecliptic); else
   `tohor_flag = SE_EQU2HOR`, `iflagtopo = iflag | SEFLG_EQUATORIAL | SEFLG_TOPOCTR` (topocentric
   equatorial — "more efficient" per comment) — and **`swe_set_topo` is called again** here.
4. Loop `i = 0..nloop-1`:
   - `swe_calc_ut(tr, ipl, iflagtopo, xx, serr)`.
   - If `SE_BIT_GEOCTR_NO_ECL_LAT`: force `xx[1] = 0` (zero out ecliptic latitude).
   - `rdi = get_sun_rad_plus_refr(ipl, xx[2], rsmi, refr)` — target horizon-altitude offset:
     see §3.4.
   - `swe_azalt(tr, tohor_flag, dgeo, atpress, attemp, xx, xaz)` — altitude at `tr`.
   - `swe_azalt(tr + 0.001, tohor_flag, dgeo, atpress, attemp, xx, xaz2)` — altitude 0.001 day
     (86.4 s) later, **using the same `xx`** (position not recomputed) — i.e. a finite-difference
     slope of the *diurnal* (azimuth/altitude) motion only, holding the object's own
     ecliptic/equatorial position fixed. `dd = xaz2[1] - xaz[1]` (altitude change per 0.001 day).
   - `dalt = xaz[1] + rdi` (signed altitude error: how far above/below the target horizon
     altitude the object currently sits — note `xaz[1]` is *true* altitude, i.e. this compares
     against true altitude, with `rdi` already including refraction, per `get_sun_rad_plus_refr`).
   - `dt = dalt / dd / 1000.0` (Newton step: `dalt / (dd/0.001)` — i.e. `dalt` divided by the
     altitude rate in degrees/day). Clamp `dt` to `[-0.1, 0.1]` days.
   - `tr -= dt`.
   - (Dead code: `if ((0) && ...) nloop++;` — an adaptive-iteration-count branch permanently
     disabled by the `(0) &&`; do not port.)
5. **Retry-next-day guard**: if `tr < tjd_ut0` (the *original* input time) **and** this is not
   already a retry (`!is_second_run`): `tjd_ut += 0.5; is_second_run = TRUE; goto run_rise_again`
   — redo the entire semi-diurnal-arc estimate (§3.2) starting half a day later, then repeat the
   Newton loop once more. At most one retry.
6. `*tret = tr; return OK`.

### 3.4 `get_sun_rad_plus_refr(ipl, dd, rsmi, refr)` — swecl.c:4176–4194
Computes the target-altitude offset (how far below/above true center-of-disc altitude 0° the
horizon-crossing condition should be evaluated), where `dd` is the body's distance (AU):
```c
rdi = 0;
if (SE_BIT_FIXED_DISC_SIZE) { if SE_SUN: dd = 1.0; else if SE_MOON: dd = 0.00257; }
if (!SE_BIT_DISC_CENTER) rdi = asin(pla_diam[ipl] / 2 / AUNIT / dd) * RADTODEG;  // angular radius
if (SE_BIT_DISC_BOTTOM) rdi = -rdi;             // bottom limb: horizon crossing is *earlier*
if (!SE_BIT_NO_REFRACTION) rdi += refr;         // add refraction (positive: raises apparent alt)
return rdi;
```
Semantics: `rdi` is added to the *geometric* altitude to get the *effective* altitude used for
the horizon test — i.e. the body "rises" when `true_altitude + rdi crosses 0` upward, which is
equivalent to comparing `true_altitude` against `-rdi`. Default (upper-limb rise, with
refraction): `rdi = +angular_radius + refr` (target altitude is *below* 0°, since the upper limb
appears at the horizon before the center does).

`pla_diam[]` (sweph.h:315) is a fixed table of body diameters in meters, indexed by `ipl` for
`ipl < NDIAM = SE_VESTA+1`; `AUNIT = 1.49597870700e11` m (sweph.h:273, DE431 value).

**STATELESS PORT NOTE:** `swe_set_topo` mutates `swed.topd` (global observer-position cache),
which `swe_calc_ut(..., SEFLG_TOPOCTR, ...)` reads back internally. The Rust `Ephemeris` is
`&self`/config-only — the port must thread `geopos` explicitly into whatever topocentric
position call replaces `swe_calc_ut`, rather than mutating shared state. Note `rise_set_fast`
calls `swe_set_topo` **twice** (once in setup, once again inside the Newton-loop setup at
swecl.c:4293–4294) with identical arguments — purely a C-state-cache artifact, not a semantic
difference; the Rust port need only pass `geopos` once, structurally.

---

## 4. `swe_rise_trans` — swecl.c:4355–4383 (dispatcher)

```c
int32 CALL_CONV swe_rise_trans(double tjd_ut, int32 ipl, char *starname,
               int32 epheflag, int32 rsmi, double *geopos,
               double atpress, double attemp, double *tret, char *serr)
```
Fast-path eligibility (swecl.c:4372–4378) — **all** of:
1. not a fixed star (`starname` NULL/empty),
2. `rsmi & (SE_CALC_RISE|SE_CALC_SET)` (i.e. not a transit request — those always go through
   `swe_rise_trans_true_hor` → `calc_mer_trans`, see §5.1),
3. `!(rsmi & SE_BIT_FORCE_SLOW_METHOD)`,
4. no twilight bit set (`CIVIL|NAUTIC|ASTRO_TWILIGHT`),
5. `ipl` in `[SE_SUN, SE_TRUE_NODE]` (i.e. Sun..lunar nodes/apsides — the "classic" body range,
   not asteroids beyond `SE_TRUE_NODE`),
6. `|geopos[1]| <= 60`, **or** (`ipl == SE_SUN` and `|geopos[1]| <= 65`).

If eligible → `rise_set_fast(...)`. Otherwise →
`swe_rise_trans_true_hor(tjd_ut, ipl, starname, epheflag, rsmi, geopos, atpress, attemp,
horhgt=0, tret, serr)`.

Return convention (documented in the header comment above `swe_rise_trans`, swecl.c:4327–4353):
`OK` (event found), `ERR` (calculation error, `serr` set), or **`-2`** meaning "the body does not
rise or set" (only ever returned by the full algorithm, §5.6).

---

## 5. `swe_rise_trans_true_hor` — swecl.c:4387–4686 (full algorithm)

```c
int32 CALL_CONV swe_rise_trans_true_hor(double tjd_ut, int32 ipl, char *starname,
               int32 epheflag, int32 rsmi, double *geopos,
               double atpress, double attemp, double horhgt,
               double *tret, char *serr)
```
`horhgt` = height of the local horizon above/below the astronomical (sea-level) horizon, in
degrees; `-100` is a sentinel meaning "use the dip of a sea-level ocean horizon as seen from
`geopos[2]` meters altitude" (§5.5).

### 5.1 Setup (swecl.c:4396–4446)
- Validate `geopos[2]` (altitude) within `[SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX] = [-500, 25000]`
  meters → else `ERR`.
- `horhgt == -100` → `horhgt = 0.0001 + calc_dip(geopos[2], atpress, attemp, const_lapse_rate)`
  (§5.5.1).
- Pluto asteroid-number alias: `ipl == SE_AST_OFFSET+134340` → treat as `SE_PLUTO`.
- `iflag = epheflag & (SEFLG_EPHMASK | SEFLG_NONUT | SEFLG_TRUEPOS)` — **all other input flags
  are dropped** (notably any caller-supplied `SEFLG_TOPOCTR`/`SEFLG_EQUATORIAL`/`SEFLG_SPEED`);
  only the ephemeris selector bits plus `NONUT`/`TRUEPOS` (perf opts) survive.
- `SE_BIT_GEOCTR_NO_ECL_LAT` → `tohor_flag = SE_ECL2HOR` (iflag left geocentric-ecliptic);
  else → `tohor_flag = SE_EQU2HOR`, `iflag |= SEFLG_EQUATORIAL|SEFLG_TOPOCTR`,
  `swe_set_topo(geopos[0], geopos[1], geopos[2])` (STATEFUL, same note as §3.4).
- **Transit dispatch**: `rsmi & (SE_CALC_MTRANSIT|SE_CALC_ITRANSIT)` → return
  `calc_mer_trans(...)` immediately (§6) — the rest of this function (§5.2–5.6) is rise/set only.
- If neither `SE_CALC_RISE` nor `SE_CALC_SET` set, default `rsmi |= SE_CALC_RISE`.
- **Twilight**: if `ipl == SE_SUN` and any twilight bit set:
  `rsmi |= (SE_BIT_NO_REFRACTION | SE_BIT_DISC_CENTER); horhgt = -rdi_twilight(rsmi)` —
  `rdi_twilight` returns `6`/`12`/`18` for civil/nautical/astronomical (first-matching bit wins,
  in that priority order per `rdi_twilight`'s `if` chain at swecl.c:4164–4174 — if multiple
  twilight bits were set simultaneously, `ASTRO_TWILIGHT` (checked last) wins since each `if` is
  unconditional, not `else if`). The target altitude becomes `-6°`/`-12°`/`-18°`, encoded as a
  *negative horhgt* — reusing the same subtraction machinery as a real horizon-dip.

### 5.2 Culmination pre-pass — 15-point sampling (swecl.c:4447–4552)
Rationale (comment, swecl.c:4447–4454): culmination (max/min altitude) points are located first
and inserted into the sample mesh so that near-horizon local extrema (which could otherwise be
missed or misclassified) are captured — meridian-transit times are deliberately *not* used for
this, since in polar regions or short-arc cases the true culmination can deviate substantially
from the meridian crossing.

1. `jmax = 14`; `twohrs = 1/12` (day) `= 2 hours`.
2. Loop `ii = 0..jmax` (15 points), `t = tjd_ut - twohrs + ii·twohrs` — samples span
   `[tjd_ut - 2h, tjd_ut + 26h]` (28 hours total, centered 2h before the input time).
3. Fixed-star position (`swe_fixstar`) is computed **once**, before the loop, at
   `tjd_et = tjd_ut + Δt(tjd_ut)` — for a fixed star, `xc[]` is **never recomputed** inside this
   entire function (all `swe_calc`/`swe_fixstar` re-evaluation is gated by `if (!do_fixstar)`);
   only the diurnal rotation (via `swe_azalt`'s internal `swe_sidtime(t)`) makes the sampled
   altitude vary with `t`. This is a deliberate optimization (proper motion is negligible over
   28 hours) — the Rust port must replicate "compute ecliptic/equatorial position once, vary
   only sidereal time per sample" for stars, not recompute position at every sample.
4. Per sample: compute `xc` (planet, via `swe_calc(t+Δt(t), ...)`) or reuse fixstar `xc`; if
   `SE_BIT_GEOCTR_NO_ECL_LAT`, zero `xc[1]`.
5. Disc diameter `dd` (meters) resolved **once**, at `ii==0`: `0` for fixstar or
   `SE_BIT_DISC_CENTER`; else `pla_diam[ipl]` if `ipl < NDIAM`; else `swed.ast_diam*1000` if
   `ipl > SE_AST_OFFSET`; else `0`.
6. `curdist = xc[2]` (AU), overridden to `1.0` (Sun) / `0.00257` (Moon) if
   `SE_BIT_FIXED_DISC_SIZE`. `rdi = asin(dd/2/AUNIT/curdist) * RADTODEG` (angular radius, degrees).
7. `swe_azalt(t, tohor_flag, geopos, atpress, attemp, xc, xh[ii])` → `xh[ii] = [az, true_alt,
   apparent_alt]`.
8. Disc-limb adjustment: `SE_BIT_DISC_BOTTOM` → `xh[ii][1] -= rdi` (bottom limb); else
   `xh[ii][1] += rdi` (top/upper limb) — applied to the **true** altitude slot.
9. Horizon-height (`horhgt`) + refraction handling — two branches:
   - `SE_BIT_NO_REFRACTION`: `xh[ii][1] -= horhgt; h[ii] = xh[ii][1]` (compare true altitude
     directly against the dip/twilight-adjusted horizon).
   - else: convert the limb-adjusted **true**-altitude point back to equatorial via
     `swe_azalt_rev(t, SE_HOR2EQU, geopos, xh[ii], xc)`, then forward again through
     `swe_azalt(t, SE_EQU2HOR, geopos, atpress, attemp, xc, xh[ii])` — this round-trip
     re-derives the **apparent** (refracted) altitude of the limb-adjusted point (azalt only
     refracts the *center*; the round-trip is how the code gets a refracted altitude for a
     point offset by the disc radius). Then `xh[ii][1] -= horhgt; xh[ii][2] -= horhgt;
     h[ii] = xh[ii][2]` (apparent altitude, dip-adjusted).
10. Culmination detection (only for `ii > 1`, i.e. once 3 points are available):
    `dc = [xh[ii-2][1], xh[ii-1][1], xh[ii][1]]` (the **true**-altitude series, always, even in
    the refraction branch — culmination uses the top/bottom-limb true altitude series
    `xh[*][1]`, not `h[*]`). Local max: `dc[1]>dc[0] && dc[1]>dc[2]`; local min: `dc[1]<dc[0] &&
    dc[1]<dc[2]`.
11. If a culmination is detected: refine via `find_maximum` (§1) with a shrinking window:
    - `dt = twohrs`; `tcu = t - dt`; `find_maximum(dc[0],dc[1],dc[2],dt,&dtint,&dx); tcu += dtint
      + dt` (nets to `tcu = t + dtint`, per §1's return-convention note — `t` here is the
      *current, rightmost* sample time).
    - Then `dt /= 3`, and while `dt > 0.0001` (days, ≈ 8.64 s): resample 3 points at
      `tcu-dt, tcu, tcu+dt` (recomputing `xc`/`swe_azalt`/`horhgt` subtraction identically to
      steps 4–9 but restricted to the no-refraction-style `ah[1] -= horhgt` altitude, regardless
      of the `SE_BIT_NO_REFRACTION` bit — the *culmination refinement* always uses true altitude
      minus `horhgt`, never the refracted round-trip), refit `find_maximum`, update
      `tcu += dtint + dt`, then `dt /= 3` again. With `dt` starting at `twohrs/3 ≈ 0.02778` day
      and dividing by 3 each round, the loop body runs 6 times before `dt` drops below `1e-4`
      day (final `dt ≈ 3.8e-5` day ≈ 3.3 s) — a fixed iteration count since `dt`'s sequence is
      deterministic (not data-dependent).
    - `nculm++; tculm[nculm] = tcu` (at most a handful of culminations found across the 15-point
      pass; `tculm[4]`/`nculm` start at `-1`).

### 5.3 Inserting culminations into the mesh (swecl.c:4556–4610)
For each found culmination `tculm[i]`, find its slot among `tc[1..jmax]` (`tculm[i] < tc[j]`),
shift `tc[j..jmax]`/`h[j..jmax]` up by one, insert `tc[j] = tculm[i]`, recompute a fresh
`xc`/disc-radius/`swe_azalt` (same limb + horhgt + refraction logic as steps 4–9 above) to get
`h[j]`, and `jmax++`. This grows the mesh from 15 points to `15 + (number of culminations found)`
points, each culmination now bracketed by real altitude samples on both sides.

### 5.4 Zero-crossing search — sign change + 20-iteration bisection (swecl.c:4611–4685)
1. Scan `ii = 1..jmax`: skip if `h[ii-1]·h[ii] >= 0` (no sign change — no crossing in this
   interval); skip if the crossing direction doesn't match the requested event
   (`h[ii-1] < h[ii]` is a *rising* crossing — skip unless `SE_CALC_RISE`; `h[ii-1] > h[ii]` is a
   *setting* crossing — skip unless `SE_CALC_SET`).
2. On a matching sign-change interval `[tc[ii-1], tc[ii]]`: **bisection**, not `find_zero`'s
   parabola-root — 20 fixed iterations:
   ```c
   for (i = 0; i < 20; i++) {
     t = (t2[0] + t2[1]) / 2;
     // recompute xc, rdi, swe_azalt, limb + horhgt + refraction-round-trip exactly as §5.2 steps 4-9
     // (aha = resulting altitude, using h[]'s exact formula: xh[1] if NO_REFRACTION else xh[2])
     if (aha * dc[0] <= 0) { dc[1] = aha; t2[1] = t; }
     else                  { dc[0] = aha; t2[0] = t; }
   }
   ```
   Interval halves each iteration; worst case starting width is `twohrs` = 7200 s (culmination
   insertion can only make the bracketing interval *narrower*, never wider). After 20 halvings:
   `7200 / 2^20 ≈ 6.87 ms` final resolution — **not** microsecond-level; treat any claim of
   ~86 µs accuracy for this loop as incorrect (that would require ~30 bisections from a 2-hour
   window, not 20).
3. After the 20 iterations, `t` holds the converged crossing time. Accept only if `t > tjd_ut`
   (must be strictly after the input time — guards against converging to a spurious crossing
   at/before the search start); on acceptance, `*tret = t; return OK`.
4. If no interval yields an accepted crossing across the whole scan: `serr = "rise or set not
   found for planet %d"`; **return `-2`** (circumpolar / does not rise or set — the documented
   `-2` convention from `swe_rise_trans`'s header comment, §4).

### 5.5 Horizon dip / twilight (`horhgt`)
`horhgt` (degrees) is subtracted from every computed altitude before the horizon (zero) test —
a uniform mechanism for three distinct physical effects:
- **True horizon dip** (visible sea horizon below the sensible horizon when observing from
  altitude): sentinel `horhgt == -100` triggers `calc_dip` (§5.5.1), a small **positive** value
  (dip lowers the effective horizon, so the body appears to rise/set later/earlier — subtracting
  a positive `horhgt` from altitude before comparing to 0 is equivalent to comparing against
  `+horhgt`, i.e. requiring the body to clear the dip).
- **Twilight**: `horhgt = -rdi_twilight(rsmi)` — a **negative** value (6/12/18° below the true
  horizon), reusing the same subtraction: comparing `altitude - (-N°) = altitude + N°` against 0
  is equivalent to requiring `altitude >= -N°`, the standard twilight-depression convention.
- **Caller-supplied arbitrary horizon height** (e.g. mountains) via `swe_rise_trans_true_hor`'s
  `horhgt` parameter directly (any sign).

#### 5.5.1 `calc_dip(geoalt, atpress, attemp, lapse_rate)` — swecl.c:3158–3172 (static)
```c
krefr = (0.0342 + lapse_rate) / (0.154 * 0.0238);
d = 1 - 1.8480*krefr*atpress/(273.15+attemp)/(273.15+attemp);
return -180.0/PI * acos(1 / (1 + geoalt/EARTH_RADIUS)) * sqrt(d);
```
Based on A. Thom, *Megalithic Lunar Observations* (1973), metric conversion by V. Reijs (2000).
`lapse_rate` is always `const_lapse_rate` (module-level static, default `SE_LAPSE_RATE =
0.0065` K/m, settable via `swe_set_lapse_rate` — STATELESS PORT NOTE: this is another piece of
global mutable state; the Rust port should take lapse rate as an explicit config/parameter
rather than a module static).

### 5.6 Return-code summary
| Return | Meaning |
|---|---|
| `OK` (0) | Event found; `*tret` set. |
| `ERR` (-1) | Ephemeris/calculation failure (`swe_calc`/`swe_fixstar`/`swe_deltat_ex` error); `serr` set. |
| `-2` | No rise/set/crossing found in the searched window (circumpolar body, or event doesn't occur) — only from §5.4 step 4. `rise_set_fast` (fast path) **never** returns `-2`; it always produces *some* `tr` estimate, even for circumpolar bodies (see §3.2 clamps). |

---

## 6. `calc_mer_trans` — swecl.c:4688–4748 (static; meridian/anti-meridian transit)

```c
static int32 calc_mer_trans(double tjd_ut, int32 ipl, int32 epheflag, int32 rsmi,
               double *geopos, char *starname, double *tret, char *serr)
```
1. `iflag = (epheflag & SEFLG_EPHMASK) | SEFLG_EQUATORIAL | SEFLG_TOPOCTR`.
2. `armc0`: local sidereal time in **hours**, `swe_sidtime(tjd_ut) + geopos[0]/15`, wrapped to
   `[0,24)`, then `*15` → degrees.
3. Initial position `x0` via `swe_fixstar` or `swe_calc` at `tjd_et = tjd_ut + Δt(tjd_ut)`.
4. `x = x0` (RA/dec copied); `t = tjd_ut`; `arxc = armc0`, or `armc0+180` if
   `SE_CALC_ITRANSIT` (lower transit target is the anti-meridian).
5. **4 fixed iterations** (`i = 0..3`):
   - `mdd = swe_degnorm(x[0] - arxc)` (meridian distance of the body's current RA from the
     target meridian, degrees); if `i > 0 && mdd > 180`, `mdd -= 360` (wrap correction — only
     applied on iterations after the first, since the first `mdd` is already normalized to
     `[0,360)` by construction and a >180 value there is expected/correct as "go forward less
     than half a sidereal day").
   - `t += mdd / 361` — Newton-like step converting an RA/armc mismatch (degrees) into a time
     delta, using the constant **361°/day** (a rounded approximation of the sidereal rate
     `ARMCS ≈ 360.9856`; not the precise constant — replicate `361` literally for FP fidelity).
   - Recompute `armc` at the new `t` (same hours→wrap→degrees pattern as step 2), `arxc = armc`
     (or `+180` for ITRANSIT).
   - Recompute `x` via `swe_calc`/`swe_fixstar` at `t + Δt(t)` (skipped for fixed stars — `x`
     stays `x0` for the entire function, same "stars don't move" optimization as §5.2).
6. `*tret = t; return OK`.

No iteration-count/convergence check — always exactly 4 iterations, unconditionally. No `-2`
path (meridian transits always exist, unlike rise/set) — only `OK` or `ERR` (from
`swe_calc`/`swe_fixstar` failure).

---

## 7. Supporting Global-State & Constant Reference

| Symbol | Value / definition | Note |
|---|---|---|
| `AUNIT` | `1.49597870700e11` m (sweph.h:273, DE431) | AU in meters |
| `NDIAM` | `SE_VESTA + 1` (sweph.h:314) | Upper bound of `pla_diam[]` table |
| `pla_diam[]` | sweph.h:315 — body diameters in meters, indexed by `ipl` | e.g. Sun = `1392000000.0` |
| `SEI_ECL_GEOALT_MIN` / `_MAX` | `-500.0` / `25000.0` (sweph.h:198–199) | Valid observer altitude range, meters |
| `SE_LAPSE_RATE` | `0.0065` K/m (sweph.h:306) | Default troposphere lapse rate |
| `SE_ECL2HOR`/`SE_EQU2HOR` | `0`/`1` (swephexp.h:364–365) | `swe_azalt` input-frame selector |
| `SE_HOR2EQU`/(`SE_HOR2ECL`) | `1`/`0` | `swe_azalt_rev` output-frame selector |
| `SE_TRUE_TO_APP`/`SE_APP_TO_TRUE` | `0`/`1` (swephexp.h:370–371) | `swe_refrac`/`swe_refrac_extended` direction |

**STATELESS PORT NOTE (summary):** the entire family relies on two pieces of C global mutable
state that the stateless Rust `Ephemeris` must instead pass explicitly:
1. `swe_set_topo(lon, lat, alt)` → `swed.topd`, read back by `swe_calc`/`swe_calc_ut` whenever
   `SEFLG_TOPOCTR` is set. Called redundantly (same args) at multiple points in both
   `rise_set_fast` and `swe_rise_trans_true_hor` — purely a caching artifact in C, not a
   semantic requirement.
2. `const_lapse_rate` (module-level `static TLS double`, default `SE_LAPSE_RATE`, mutated via
   `swe_set_lapse_rate`, swecl.c:74/2988) — feeds `swe_refrac_extended`/`calc_dip` everywhere in
   this file. Should become an explicit parameter (already part of `EphemerisConfig` if a
   similar knob exists for refraction elsewhere in the Rust port — check `docs/codebase-map.md`
   before adding a duplicate).

`swe_azalt`/`swe_azalt_rev` (swecl.c:2788–2825, 2839–2873) are themselves stateful only through
`const_lapse_rate` and `swe_sidtime`/`swe_calc` (which need Δt) — no additional caches beyond
what's listed above. `xaz = [azimuth (from south, clockwise), true_altitude,
apparent_altitude]` — note the **azimuth convention flip** inside `swe_azalt` (computed
"from east, counterclockwise" internally via `swe_cotrans`, then converted to "from south,
clockwise" for the public `xaz[0]`); replicate both conventions if the Rust port exposes an
internal east-CCW azimuth anywhere (e.g. if `swe_azalt`'s internal helper is factored out and
reused, as the rise/set port likely will do for performance).
