# C Reference: Lunar Occultations — swecl.c

Porting reference for `swe_lun_occult_where`, `swe_lun_occult_when_glob`, and
`swe_lun_occult_when_loc`. All line numbers below refer to `swecl.c`.

**Read `docs/c-ref-eclipse-solar.md` first.** Occultations of a planet or fixed star by
the Moon reuse the *exact same* geometry engine as solar eclipses — the Sun is simply the
`ipl == SE_SUN`, `starname == NULL` special case of a body already generic over
`ipl`/`starname`:

- `eclipse_where()` (swecl.c:640–886) — geographic position + shadow geometry. Already
  parameterized over `ipl`/`starname`; `swe_sol_eclipse_where` and `swe_lun_occult_where`
  both call it directly with no wrapper differences.
- `eclipse_how()` (swecl.c:967–1152) — local circumstances (magnitude, obscuration,
  azimuth/altitude, contact geometry) at a given geographic point. Same: already generic.
- `calc_planet_star()` (swecl.c:888–897) — the one-line dispatcher that is the *entire*
  mechanism for the ipl/starname split: `starname` empty → `swe_calc`; `starname` set →
  `swe_fixstar`. Every occultation function funnels the occulted body through this.
- `find_maximum()` / `find_zero()` (swecl.c:4133–4162) — generic parabola vertex / root
  finder used identically here.

This doc does **not** re-derive any of that. It documents only:
1. The thin `swe_lun_occult_*` wrappers and how they thread `ipl`/`starname` into the
   shared engine (near-zero delta from the solar wrappers).
2. Where the *search strategy* differs — solar eclipse search exploits the fact that the
   Sun's motion is regular (uses a lunation-number (K) closed-form Meeus estimate for the
   next new moon); occultation search cannot do this because `ipl` ranges over Mercury
   through Pluto/asteroids (and fixed stars, which don't move in longitude at all), so it
   falls back to a generic Newton-style bracketing search on Moon–body elongation.
3. Where the point-source (star) vs finite-disc (planet) distinction changes radius
   handling, contact-time semantics, and eclipse-type validity (no annular occultations
   except for the Sun).

---

## Function Map

| C function | Location | Purpose |
|---|---|---|
| `swe_lun_occult_where` | swecl.c:606–630 | Public: geographic position of maximal occultation at given time |
| `eclipse_where` (shared) | swecl.c:640–886 | See `docs/c-ref-eclipse-solar.md` — unmodified, generic over ipl/starname |
| `calc_planet_star` (shared) | swecl.c:888–897 | See solar doc — `swe_calc` vs `swe_fixstar` dispatch |
| `eclipse_how` (shared) | swecl.c:967–1152 | See solar doc — local circumstances at a geo point |
| `swe_lun_occult_when_glob` | swecl.c:1572–1984 | Public: next/previous occultation of ipl/star anywhere on Earth |
| `swe_lun_occult_when_loc` | swecl.c:2071–2098 | Public: thin wrapper — calls `occult_when_loc`, then re-derives core-shadow diameter via `eclipse_where` |
| `occult_when_loc` (static) | swecl.c:2412–2764 | Worker: next occultation of ipl/star visible from a fixed geo location |
| `eclipse_when_loc` (solar sibling) | swecl.c:2100–2410 | For contrast only — see § 3 below; documented fully in solar doc |
| `find_maximum` / `find_zero` (shared) | swecl.c:4133–4162 | See solar doc |

`SE_ECL_ONE_TRY` (swephexp.h:331) `= (32*1024) = SEFLG_TOPOCTR` (swephexp.h:206) — the two
constants are numerically identical by design; a C comment at swecl.c:1544–1548 explains
this is deliberate (`ifl` may already carry `SEFLG_TOPOCTR`, and the code specifically does
*not* fold `SE_ECL_ONE_TRY` into `ifl` to avoid confusing the two). **Port hazard**: do not
represent `SE_ECL_ONE_TRY` and `SEFLG_TOPOCTR` as the same Rust enum/bitflag member — keep
them as textually-separate constants even though the bit position is shared, and never
`ifl |= backward` style casts; keep the `backward: i32` (or bitflag) input strictly
separate from the `flags: EphemerisFlags`-equivalent input.

---

## 1. `swe_lun_occult_where` (swecl.c:606–630)

```c
int32 CALL_CONV swe_lun_occult_where(
    double tjd_ut, int32 ipl, char *starname, int32 ifl,
    double *geopos, double *attr, char *serr)
```

Thin wrapper, essentially identical in shape to `swe_sol_eclipse_where` (see solar doc)
except it threads `ipl`/`starname` instead of hardcoding `SE_SUN`/`NULL`:

1. `if (ipl < 0) ipl = 0;` — defensive clamp (swecl.c:617).
2. `ifl &= SEFLG_EPHMASK;` then `swi_set_tid_acc(tjd_ut, ifl, 0, serr)` — same tidal-acceleration
   global-state setup as the solar path.
   **STATELESS PORT NOTE**: `swi_set_tid_acc` mutates C global state (`swed.tid_acc`) as a
   side effect purely to affect the *subsequent* `swe_deltat_ex` calls inside
   `eclipse_where`/`eclipse_how`. In the Rust port this must instead be an explicit
   `tid_acc` parameter (or derived from `EphemerisConfig`) threaded into the deltaT call —
   do not mutate shared state. This is the same concern already resolved for the solar
   eclipse port; reuse whatever mechanism was adopted there.
3. **Asteroid 134340 aliasing** (swecl.c:620–623): if `ipl == SE_AST_OFFSET + 134340`
   (i.e., asteroid-numbered Pluto), it is silently remapped to `ipl = SE_PLUTO`. This is
   occultation/eclipse-specific special-casing not present in generic `swe_calc` — port it
   as an explicit normalization step at the top of the occultation entry points (present in
   all three: `swe_lun_occult_where`, `swe_lun_occult_when_glob` at swecl.c:1599–1600, and
   `swe_lun_occult_when_loc` at swecl.c:2084–2085).
4. Calls `eclipse_where(tjd_ut, ipl, starname, ifl, geopos, dcore, serr)` (swecl.c:624) —
   fully shared, see solar doc for the shadow-cone / oblateness-iteration algorithm. `ipl`/
   `starname` only affect this function through:
   - the `drad` (occulted-body angular radius in AU) computation at swecl.c:697–704 (see
     § "Radius handling" below — identical logic block appears 4 times in this file:
     swecl.c:697, 1004, 1681, 2489).
   - `calc_planet_star` being invoked instead of a hardcoded `swe_calc(..., SE_SUN, ...)`.
5. Calls `eclipse_how(tjd_ut, ipl, starname, ifl, geopos[0], geopos[1], 0, attr, serr)`
   (swecl.c:626) — geo height passed as `0` (sea level), matching the "where" semantics
   (this is the position of *maximum* occultation, not a specific observer altitude).
6. `attr[3] = dcore[0]` (swecl.c:628) — core "shadow" width in km, overwritten after
   `eclipse_how` already zeroed/filled attr[]. Same convention as solar (see solar doc
   for full `attr[]` index table — unchanged here).
7. Returns `retflag` from `eclipse_where` (the `eclipse_how` result `retflag2` is computed
   only to catch `ERR`; its normal return value is intentionally discarded here — same
   pattern as `swe_sol_eclipse_where`).

### Radius handling: point source (star) vs finite disc (planet)

This exact 4-armed branch recurs verbatim at swecl.c:697–704 (`eclipse_where`), 1004–1011
(`eclipse_how`), 1681–1688 (`swe_lun_occult_when_glob`), and 2489–2496 (`occult_when_loc`).
Port it as a single shared helper (e.g. `occulted_body_radius_au(ipl, starname, ast_diam_km)`)
rather than duplicating it 4×:

```c
if (starname != NULL && *starname != '\0')
    drad = 0;                                    /* fixed star: point source */
else if (ipl < NDIAM)                            /* NDIAM = SE_VESTA + 1, sweph.h:314 */
    drad = pla_diam[ipl] / 2 / AUNIT;             /* named body: half the tabulated diameter */
else if (ipl > SE_AST_OFFSET)
    drad = swed.ast_diam / 2 * 1000 / AUNIT;      /* numbered asteroid: km -> m -> AU, radius from ephemeris file header */
else
    drad = 0;
```

`pla_diam[]` (sweph.h:315–333) is a static table of body diameters **in meters**, indexed
`SE_SUN=0 .. SE_VESTA=17` (`NDIAM = SE_VESTA + 1`). Relevant occultable-body entries: Mercury
`2439400*2`, Venus `6051800*2`, Mars `3389500*2`, Jupiter `69911000*2`, Saturn `58232000*2`,
Uranus `25362000*2`, Neptune `24622000*2`, Pluto `1188300*2`, Chiron `271370`, Pholus
`290000`, Ceres `939400`, Pallas `545000`, Juno `246596`, Vesta `525400` (entries 10–13, the
lunar nodes/apogees, are `0`). `swed.ast_diam` is a **global mutable field populated from the
current asteroid's `.se1`/AST file header** at load time (a stateful side effect of
`swe_calc`/`swi_sweph` for numbered asteroids beyond Vesta) — this is the one place occultation
radius handling depends on ephemeris-file global state beyond tidal acceleration.
**STATELESS PORT NOTE**: the Rust asteroid backend must return the body's radius as part of
its calc result (or from a lookup table / file header read at call time) rather than reading
a cached global; `Ephemeris` has no `ast_diam` cache field, so this must be plumbed through
explicitly, e.g. returned alongside the position from the SE1-file backend or looked up from
the AST file header on each call.

`RMOON = DMOON/2` where `DMOON = 3476300.0 / AUNIT` (swecl.c:82) is a *constant* lunar
diameter (distinct from `pla_diam[SE_MOON] = 3475000.0` used elsewhere — a ~1.3 km / 0.04%
inconsistency baked into the C source between the eclipse-specific `RMOON`/`DMOON` macros
and the generic `pla_diam` table; preserve both constants separately, do not unify them).

A **fixed star has `drad = 0` always** — it is a mathematical point. This has three
downstream consequences documented in §§ 2–3:
- No annular/total occultation distinction is meaningful in the shadow-cone sense used by
  `eclipse_where`'s `SE_ECL_CENTRAL`/`NONCENTRAL`/shadow-width outputs — those still compute
  (the formulas are algebraically well-defined at `drad=0`), but "annular" for `drad=0`
  degenerates to "the star's line grazes exactly through the Moon's disc," which in practice
  the `retc` classification in `eclipse_where` (swecl.c:876–884: `*dcore > 0` → annular,
  else total) still produces one or the other; callers should treat star occultations as
  binary total/no-event since there is no partial "ring" case worth distinguishing from
  total once the occulting body (Moon) is always the finite disc.
- In `swe_lun_occult_when_glob`, annular/annular-total eclipse types are explicitly
  **rejected as an error** for any `ipl != SE_SUN` (see § 2, swecl.c:1614–1623).
- In `occult_when_loc`, contacts 1/4 (penumbra-equivalent, disc-edge-to-disc-edge) collapse
  to being identical to contacts 2/3 (umbra-equivalent) for a point source — see § 3,
  swecl.c:2696–2699.

### `SE_ECL_ONE_TRY` behavior in `_where`

`swe_lun_occult_where` itself takes no `backward`/search-control parameter at all — it
computes local geometry at a single given `tjd_ut`, so `SE_ECL_ONE_TRY` is irrelevant here.
It only matters for the two `_when_*` search functions (§§ 2–3), where it is packed into the
same integer as `backward`.

---

## 2. `swe_lun_occult_when_glob` (swecl.c:1572–1984)

```c
int32 CALL_CONV swe_lun_occult_when_glob(
    double tjd_start, int32 ipl, char *starname, int32 ifl, int32 ifltype,
    double *tret, int32 backward, char *serr)
```

Structurally a near-duplicate of `swe_sol_eclipse_when_glob` (documented in the solar doc) —
same `tret[0..9]` contact-time slots, same `ifltype` bit semantics, same central/non-central/
total/annular/partial classification and rejection-and-retry control flow. Differences below;
everything not called out is identical to the solar version and should reuse that logic
(ideally the same Rust function generic over occulted-body, matching the C's shared helper
functions).

### Setup differences from `swe_sol_eclipse_when_glob`

1. **Asteroid-134340 aliasing** to `SE_PLUTO` (swecl.c:1599–1600) — same as § 1.3.
2. **`one_try` flag** (swecl.c:1593): `int32 one_try = backward & SE_ECL_ONE_TRY;` extracted
   *before* `backward &= 1L` (swecl.c:1605) truncates `backward` down to just its direction
   bit. `one_try` has no equivalent in the solar global search — it lets the caller do a
   single conjunction check without looping until an actual occultation is found (needed
   because non-Sun bodies plus fixed stars can go a long time between occultations — the
   docstring at swecl.c:1541–1543 warns the search "may search successlessly until it
   reaches the end of the ephemeris" without it).
3. **Annular type validity depends on `ipl`** (swecl.c:1614–1628):
   - If `ipl != SE_SUN` and the caller requested `SE_ECL_ANNULAR` or `SE_ECL_ANNULAR_TOTAL`
     explicitly (after masking out `SE_ECL_CENTRAL`/`NONCENTRAL`), this is an **error**:
     `"annular occulation do not exist for object %d %s"`. (Typo "occulation" is in the C
     source; preserve or fix at Rust discretion but note the original spelling for
     traceability.)
   - If `ipl != SE_SUN` and the caller passed `ifltype` containing
     `SE_ECL_ANNULAR|SE_ECL_ANNULAR_TOTAL` alongside other acceptable bits, those bits are
     silently stripped (swecl.c:1622–1623) rather than erroring.
   - Default `ifltype == 0` expansion (swecl.c:1624–1628) omits `SE_ECL_ANNULAR |
     SE_ECL_ANNULAR_TOTAL` entirely unless `ipl == SE_SUN`. **Rust port**: this means the
     default "any occultation type" search space differs by occulted body — encode this as
     a function of `ipl` at the point where you build the default type mask, not a constant.

### Search-stepping algorithm — the actual delta from solar

This is the substantive difference from the solar version, and the reason a "sibling
machinery" doc alone isn't sufficient — solar eclipse search exploits the ~29.53-day
synodic month via a closed-form lunation-number (K) formula (Meeus) to jump directly near
the next new moon. Occultation search must work for **any** body (arbitrary sidereal
period) and for fixed stars (**zero** proper motion in the search timescale), so it uses a
generic bracket-and-refine approach instead:

1. **Rough conjunction in ecliptic longitude** (swecl.c:1640–1666, `next_try:` label):
   - Compute `ls[0]` (occulted body's ecliptic longitude, geocentric, via
     `calc_planet_star(t, ipl, starname, ifl, ls, serr)`) and `lm[0]` (Moon's ecliptic
     longitude via plain `swe_calc`).
   - Guard (swecl.c:1646–1650): if `starname` given and `fabs(ls[1]) > 7`
     (ecliptic latitude of the star, in **degrees**), immediately error — a fixed star this
     far from the ecliptic can never be occulted regardless of lunar parallax/proper
     motion, so the search is aborted rather than run to the end of the ephemeris.
   - `dl = swe_degnorm(ls[0] - lm[0])`; if searching backward (`direction < 0`), subtract
     360 so `dl` is signed consistently with the search direction.
   - **Newton-like linear step**: `while (fabs(dl) > 0.1) { t += dl / 13; ...recompute...; }`
     — `13` here approximates the Moon's mean motion relative to a slow-moving outer body,
     ≈13°/day (the Moon's actual mean motion is ~13.2°/day; the Sun/planets move far slower
     in comparison, so `dl/13` converges to the next Moon–body conjunction in longitude).
     This is **not** a fixed number of iterations — it's a `while` loop to convergence
     within 0.1° of exact conjunction. Recomputes `ls`/`lm` (ecliptic longitude only, no
     cartesian) on every iteration.
   - Sets `tjd = t` after convergence (swecl.c:1666).
2. **Latitude-difference gate** (swecl.c:1667–1677): `drad = fabs(ls[1] - lm[1])` (ecliptic
   latitude difference occulted body vs Moon, in degrees — note `drad` variable name reused,
   unrelated to the later occulted-body angular radius `drad`). If `> 2` degrees, no
   occultation is possible at this conjunction:
   - if `one_try`: return immediately with `retflag=0`, `tret[0] = t + direction` (a
     "suitable date for next try", *not* a valid occultation time — caller is expected to
     re-invoke with this as the new `tjd_start`).
   - else: `t += direction * 20` (jump 20 days in the search direction — roughly 1.5 Moon
     sidereal-month steps, cheap enough to just retry rather than compute the next
     conjunction analytically) and `goto next_try`.
3. **Occulted-body angular radius** `drad` recomputed (swecl.c:1678–1688) — same 4-branch
   logic as § 1's radius handling (note this **overwrites** the `drad` used for the latitude
   gate in step 2 — same variable, different meaning, reused after the gate check).
4. **Refine time of maximum occultation** (swecl.c:1696–1717): parabola-vertex bracketing
   via `find_maximum`, identical numerical scheme to the solar version but with different
   constants:
   - `dtstart = dadd2 = 1` (day), `dtdiv = 3` (solar global search instead starts at
     `dtstart=1` or `5` and divides by `dtdiv=4` — occultation uses a **tighter divisor** of
     3, needing more outer-loop iterations to reach the same `dt > 0.0001` convergence floor,
     compensating for not having as good an initial time estimate as the Meeus K-formula
     gives the solar search).
   - Inner loop evaluates `dc[i] = acos(dot(unit(xs), unit(xm))) * RADTODEG - (rmoon + rsun)`
     at three equally-spaced `t` values (`tjd-dt, tjd, tjd+dt`), where `xs`/`xm` are
     cartesian equatorial vectors (`iflagcart = iflag | SEFLG_XYZ`) of occulted body / Moon,
     `rmoon = asin(RMOON/lm[2])`, `rsun = asin(drad/ls[2])` (apparent angular radii from
     true geocentric distance `lm[2]`/`ls[2]`, in polar-coordinate form from a parallel
     `iflag` call). This `dc[i]` is "angular separation of centers minus sum of apparent
     radii" — its minimum (most negative / most deeply overlapping) is the occultation
     maximum, found via `find_maximum` (same parabola-vertex helper as solar/shared).
   - `tjd += dtint + dt` accumulates the correction each outer iteration; loop continues
     while `dt > 0.0001` (days), dividing by 3 each pass (log₃(1/0.0001) ≈ 9 outer
     iterations from `dt=1`).
5. `tjd -= swe_deltat_ex(tjd, ifl, serr)` (swecl.c:1718) converts the refined **ET** estimate
   back to UT — only **one** deltaT subtraction here (contrast: the solar `_when_loc`
   worker at swecl.c:1301–1303 and `occult_when_loc` at swecl.c:2561–2562 both do a
   **two-step iterated** deltaT correction: `tret[0] = tjd - swe_deltat_ex(tjd,...)` then
   `tret[0] = tjd - swe_deltat_ex(tret[0],...)`. `swe_lun_occult_when_glob` does not
   double-iterate this correction — a minor precision asymmetry vs the `_when_loc` path,
   worth flagging but not "fixing" since it mirrors the C exactly.)
6. From here on (swecl.c:1720 to end, `end_search_global:` at 1968) the control flow is
   **identical in structure** to `swe_sol_eclipse_when_glob`: call `eclipse_where` to get
   type/`dcore`; if no eclipse/occultation, retry (`t = tjd + direction*20`, `goto
   next_try`, honoring `one_try`); re-check against `tjd_start` to avoid returning the
   starting instant; classify against requested `ifltype` bits with retry-on-mismatch for
   each of noncentral/central/annular/partial/total; contact-time refinement for
   `tret[2..7]` via the same `dta=twohr, dtb=tenmin` / `find_zero` two-stage bracket
   (swecl.c:1831–1875, `dcore[]`-threshold formulas unchanged from solar); annular-total
   detection via sign changes of `dcore[0]` at max/tret[4]/tret[5] (swecl.c:1879–1897,
   though as established this branch is unreachable for `ipl != SE_SUN` since annular is
   rejected/stripped earlier); solar-transit / local-noon `tret[1]` computed the same way
   (swecl.c:1920–1967: swecl.c:1929/1931 call `calc_planet_star(tt, ipl, starname, ...)`
   and `swe_calc(tt, SE_MOON, ...)` respectively — i.e. this "local apparent noon" transit
   check is genuinely about when the **occulted body** (not necessarily the Sun) transits
   the meridian relative to the Moon, correctly generalized rather than a copy-paste
   leftover from the solar version).
7. Return value `retflag` — same bit meanings as solar (`SE_ECL_TOTAL/ANNULAR/PARTIAL/
   CENTRAL/NONCENTRAL/ANNULAR_TOTAL`), `tret[0..9]` same slot meanings as documented in the
   solar doc (max, local-noon transit, begin, end, totality-begin, totality-end,
   center-line-begin, center-line-end, unused, unused) — **caveat**: "local apparent noon"
   here means transit of the *occulted body*, and "center line" / "totality" retain their
   shadow-cone meaning from `eclipse_where`, which for a point-source star still computes
   (degenerately) rather than being suppressed.

**FP hazard**: step 4's `find_maximum` inputs use `iflag`/`iflagcart` (equatorial, no
`SEFLG_TOPOCTR`, i.e. **geocentric**) — this is the same flag combination used in the
rough-longitude step 1, but note step 1 uses **ecliptic** (`ifl` alone, no
`SEFLG_EQUATORIAL`) while step 4 uses **equatorial cartesian** (`iflagcart = iflag |
SEFLG_XYZ` where `iflag = SEFLG_EQUATORIAL | ifl`). Do not conflate the two coordinate
frames when porting — the longitude convergence loop and the parabola-fit loop use
different coordinate systems for their respective dot-products/differences.

---

## 3. `swe_lun_occult_when_loc` (swecl.c:2071–2098) + `occult_when_loc` (swecl.c:2412–2764)

```c
int32 CALL_CONV swe_lun_occult_when_loc(
    double tjd_start, int32 ipl, char *starname, int32 ifl,
    double *geopos, double *tret, double *attr, int32 backward, char *serr)
```

### Public wrapper (swecl.c:2071–2098)

1. Validate `geopos[2]` (height) within `[SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX]` =
   `[-500.0, 25000.0]` meters (sweph.h:198–199) — error `"location for occultations must be
   between %.0f and %.0f m above sea"` (contrast solar's wording: "location for eclipses
   must be...") — same bounds, different error string.
2. Asteroid-134340→Pluto aliasing (swecl.c:2084–2085), same as §§ 1, 2.
3. `swi_set_tid_acc` (swecl.c:2087).
4. Calls worker `occult_when_loc(tjd_start, ipl, starname, ifl, geopos, tret, attr,
   backward, serr)` (swecl.c:2088); if result `<= 0`, return immediately (no occultation
   found / error) — **note**: this is `<= 0`, not `== ERR`, meaning a `retflag == 0` success
   with no bits set also short-circuits here (mirrors solar's `swe_sol_eclipse_when_loc`
   exactly, swecl.c:2031).
5. Re-derives core shadow diameter: calls `eclipse_where(tret[0], ipl, starname, ifl,
   geopos2, dcore, serr)` at the found maximum time, ORs `SE_ECL_NONCENTRAL` from that into
   the returned flag, sets `attr[3] = dcore[0]` (swecl.c:2093–2096) — identical pattern to
   `swe_lun_occult_where` step 6 and to `swe_sol_eclipse_when_loc`.

### Worker `occult_when_loc` (swecl.c:2412–2764)

Structurally the occultation counterpart of `eclipse_when_loc` (solar; swecl.c:2100–2410,
documented fully in the solar doc). Key differences:

**Setup**: `iflag = SEFLG_TOPOCTR | ifl` and a derived `iflaggeo = iflag &
~SEFLG_TOPOCTR` (swecl.c:2432–2433) — the solar version doesn't need this split because it
always mixes `SEFLG_EQUATORIAL|SEFLG_TOPOCTR` directly; occultation needs a **geocentric**
variant (`iflaggeo`) specifically for the rough-conjunction bracketing (step 1 below), while
everything downstream of finding the approximate time uses full topocentric
(`iflag`/`iflagcart`). `swe_set_topo(geopos[0], geopos[1], geopos[2])` (swecl.c:2440) is
called once up front (contrast solar's `eclipse_when_loc`, which calls it twice — once at
entry swecl.c:2117 and again at swecl.c:2162 right before the refinement loop; occultation
only needs the one call since it doesn't have the K-formula estimate step to precede it).
**STATELESS PORT NOTE**: `swe_set_topo` sets global observer coordinates
(`swed.topd`) consumed implicitly by later `swe_calc(..., SEFLG_TOPOCTR, ...)` calls. In the
Rust port, `geopos`/observer location is `EphemerisConfig.topographic` (or an explicit
parameter threaded per-call) — every downstream position call in this worker that uses
`SEFLG_TOPOCTR`-equivalent must receive the observer coordinates explicitly rather than
relying on a prior "set" call's side effect.

1. **Rough conjunction search** (swecl.c:2447–2484) — **identical** in structure and
   constants to `swe_lun_occult_when_glob` § 2 step 1: same `dl/13` Newton-style longitude
   convergence loop, same `fabs(ls[1]) > 7` star-latitude rejection, using `iflaggeo`
   (geocentric — the local-visibility topocentric refinement only matters once a rough
   candidate time is known). This block is copy-identical to swecl.c:1640–1666 modulo the
   `iflaggeo` vs `ifl` flag variable name; port as the same shared helper used in § 2.
2. **Latitude-difference gate** (swecl.c:2476–2485) — identical to § 2 step 2 (`drad =
   fabs(ls[1]-lm[1]) > 2` → retry ±20 days or return one-try sentinel).
3. **Occulted-body radius** `drad` (swecl.c:2486–2496) — identical 4-branch logic, same as
   §§ 1–2.
4. **Local-visibility bracket + refine time of maximum** (swecl.c:2497–2537): same
   `find_maximum` parabola scheme as § 2 step 4, but:
   - Uses **topocentric** cartesian/polar (`iflagcart`/`iflag` here include
     `SEFLG_TOPOCTR`), unlike the glob search's geocentric version — this is the "is it
     visible from *this* observer" check baked directly into the bracketing, not deferred
     to a separate `eclipse_how` visibility pass.
   - `dtstart = dadd2 = 1`, `dtdiv = 2` (swecl.c:2498–2499). The loop body contains
     `if (dt < 0.01) dtdiv = 2;` (swecl.c:2503–2504) — a no-op reassignment to the same
     value `dtdiv` already holds, so unlike the solar `eclipse_when_loc` (which starts
     `dtdiv=2` and genuinely switches to `dtdiv=3` once `dt<0.1`, swecl.c:2170–2171),
     `occult_when_loc`'s divisor stays `2` for the entire loop. This looks like a
     copy-paste leftover from editing the solar version (the guard is present but the
     RHS wasn't updated to `3`). **Port faithfully**: use a constant `dtdiv = 2` for this
     loop, matching what the C actually executes, not what the dead conditional suggests
     was intended.
   - **Mid-loop bailout for large latitude separation** (swecl.c:2516–2525): if `dt < 0.1`
     and `fabs(ls[1] - lm[1]) > 2` (topocentric latitude gap too large for occultation) —
     under `one_try` (or if `stop_after_this` already set), just marks `stop_after_this =
     TRUE` and lets the current bracketing pass finish rather than aborting mid-loop;
     otherwise jumps to `t = tjd + direction*20` and `goto next_try` immediately. This
     topocentric re-check has no equivalent step in the geocentric global search (§ 2) —
     it exists because topocentric parallax can push a marginal geocentric near-miss into
     "definitely not occulted from here."
   - `if (stop_after_this)` after the loop (swecl.c:2534–2537): returns `retflag=0` with
     `tret[0] = tjd + direction` as a resume hint, honoring `one_try` semantics.
5. **Reject if not actually occulting** (swecl.c:2551–2560): `dctr > rsplusrm` (separation
   exceeds sum of apparent radii even at best approach) → retry or one-try-return, same as
   § 2.
6. **`tret[0]` via double-iterated deltaT** (swecl.c:2561–2562) — unlike § 2's single
   subtraction, this **does** do the two-step iteration
   (`tret[0]=tjd-Δt(tjd); tret[0]=tjd-Δt(tret[0])`), matching solar `eclipse_when_loc`
   exactly. (Flagged as an asymmetry with § 2 above — worth a golden-test tolerance note if
   `_when_glob` and `_when_loc` results for the same event are ever cross-compared to
   sub-second precision.)
7. **Type classification** `SE_ECL_ANNULAR/TOTAL/PARTIAL` (swecl.c:2574–2580) — same
   `dctr` vs `rsminusrm`/`fabs(rsminusrm)`/`rsplusrm` thresholds as solar. Note: nothing here
   *prevents* `SE_ECL_ANNULAR` from being set for `ipl != SE_SUN` at the `occult_when_loc`
   level (unlike `swe_lun_occult_when_glob`'s explicit rejection in § 2) — this worker
   computes purely from geometry (`rsun` uses the actual `drad` of whatever `ipl` is, which
   is generally smaller than the Moon, so `rsminusrm` is negative and the "annular" branch
   `dctr < rsminusrm` is essentially unreachable in practice for a small occulted body, but
   is not *structurally* guarded the way the glob search guards it). Treat "annular" here as
   "occulted body's disc fits entirely inside the Moon's disc, offset toward one edge" —
   still geometrically meaningful even for e.g. a large asteroid grazing the Moon's limb,
   just not called out with a friendly name.
8. **Contacts 2/3** (2nd/3rd contact = occulted-body disc fully immersed / starts emerging)
   (swecl.c:2582–2641) — identical `find_zero` two-stage refine (`twomin` → `tensec/10^m`)
   to solar, using `drad` (the occulted body) instead of a hardcoded `RSUN`, and the same
   empirical `rmoon *= 0.99916` fudge factor "for better accuracy for 2nd/3rd contacts"
   (swecl.c:2595, 2624 — present verbatim in the solar version too, preserve exactly, do
   not "explain away").
9. **Contacts 1/4 — the point-source branch** (swecl.c:2642–2699): **this is the one
   genuinely occultation-specific structural fork** (not present in the solar path at all,
   since the Sun is never a point source):
   ```c
   if (starname == NULL || *starname == '\0') {
       /* ... full find_zero + tenmin/10^m refine, identical to solar contacts-1/4 code ... */
   } else { /* fixed stars are point sources, contacts 1 and 4 = contacts 2 and 3 */
       tret[1] = tret[2];
       tret[4] = tret[3];
   }
   ```
   For a planet occulting/being-occulted, contacts 1/4 (disc-edge-to-disc-edge, i.e.
   first/last exterior tangency — the "penumbra" analogue) are computed via the same
   `find_zero`+`rsplusrm` bracket as solar (swecl.c:2645–2695, using `drad` instead of
   `RSUN`, no `0.99916` fudge here — that correction is 2nd/3rd-contact-only in both
   paths). For a **fixed star**, there is no disc to have an "edge" — contacts 1/4
   (occulted body's leading/trailing limb vs Moon) are physically identical to contacts 2/3
   (Moon's limb reaching the point), so the code just copies `tret[2]→tret[1]`,
   `tret[3]→tret[4]` rather than re-deriving them. **Port this exactly as a branch on
   `starname.is_some()`**, not as a unification of the math (the underlying `find_zero`
   bracket for a zero-radius source is genuinely degenerate/undefined the way the solar
   contact-1/4 code is written, since it depends on `rsplusrm = rsun + rmoon` with
   `rsun=0` collapsing contacts 1/4 onto 2/3 mathematically anyway — but the C avoids
   computing it a second time and just aliases).
10. **Visibility per contact + `attr[]` fill** (swecl.c:2700–2721): loops `i = 4..0`,
    skipping zero `tret[i]`, calling `eclipse_how(tret[i], ipl, starname, ifl, geopos[0],
    geopos[1], geopos[2], attr, serr)` for each populated contact and OR-ing in
    `SE_ECL_xTH_VISIBLE`/`SE_ECL_MAX_VISIBLE` bits when `attr[6]` (apparent altitude) `> 0`
    at that instant — **exactly** the solar pattern (see solar doc for full `attr[]` index
    table; unchanged here except `ipl`/`starname` threading through `eclipse_how` instead
    of hardcoded `SE_SUN`/`NULL`). Loop order is descending (`i=4` down to `0`) so that
    `attr[]`'s final state (kept across loop iterations, comment "attr for i=0 must be
    kept!!!" swecl.c:2703) reflects contact **0** (maximum), since `i=0` runs last.
11. **Retry if never visible** (swecl.c:2722–2733) — same `goto next_try` /
    `one_try`-sentinel pattern as everywhere else in this function, guarded by `#if 1`.
12. **Occultation-specific daylight flags** (swecl.c:2734–2762) — this block **has no
    solar equivalent at all**:
    - `swe_rise_trans` for the **occulted body** (`ipl`/`starname`, not the Sun) around
      `tret[1]` to find rise/set brackets, filling `tret[5]`/`tret[6]` if the body's
      rise/set falls strictly between contacts 1 and 4 (swecl.c:2734–2743) — mirrors solar's
      analogous rise/set bracketing but for the occulted body rather than the Sun (solar's
      `eclipse_when_loc` uses `SE_SUN` here because the eclipsed body *is* the Sun; for
      occultations the "is the event happening while the target is up" check needs the
      occulted body's own rise/set, hence a separate pair of `swe_rise_trans` calls).
    - **Two additional** `swe_rise_trans(..., SE_SUN, ...)` pairs (swecl.c:2746–2761) check
      whether the **actual Sun** is up at `tret[1]` and at `tret[4]` (contacts 1 and 4) —
      if sunset precedes sunrise in the relevant window, set
      `SE_ECL_OCC_BEG_DAYLIGHT` / `SE_ECL_OCC_END_DAYLIGHT` (swephexp.h:329–330 — numerically
      identical bit values to `SE_ECL_PENUMBBEG_VISIBLE`/`PENUMBEND_VISIBLE` used by lunar
      eclipses; the two flag families are mutually exclusive by call site, not by bit
      layout — same "shared bit position, different meaning by context" pattern as
      `SE_ECL_ONE_TRY`/`SEFLG_TOPOCTR`, flag this the same way in the port). These flags
      tell the caller whether an occultation begins/ends "during the day" (relevant because
      an occultation of e.g. Venus during broad daylight may still be observable, hence the
      docstring note at swecl.c:2054–2055 about telescope/naked-eye visibility).
    - Both blocks use `retc >= 0` gating rather than `!= ERR`: `swe_rise_trans` returns `-2`
      for "circumpolar" bodies, which here is treated as "skip this check" rather than an
      error (contrast solar's `eclipse_when_loc`, swecl.c:2376/2380, which explicitly early
      `return`s on `retc == -2` for the Sun's own rise/set — `occult_when_loc` instead just
      lets `retc >= 0` silently skip populating `tret[5]`/`tret[6]`/the daylight flags when
      circumpolar, without early-returning the whole function).

### `tret[]` slot summary (occultation `_when_loc`)

Same slots as solar `_when_loc` (see solar doc for the authoritative table) with the
addition/reuse below — no new slots introduced beyond what `swe_sol_eclipse_when_loc`
already defines:

| slot | meaning |
|---|---|
| `tret[0]` | time of maximum occultation |
| `tret[1]` | 1st contact (disc-edge tangency in / for star: aliased from `tret[2]`) |
| `tret[2]` | 2nd contact (full immersion begins / for star: same instant as 1st contact) |
| `tret[3]` | 3rd contact (emersion begins / for star: same instant as 4th contact) |
| `tret[4]` | 4th contact (disc-edge tangency out / for star: aliased from `tret[3]`) |
| `tret[5]` | occulted body's rise time, if within [1st,4th] contact window |
| `tret[6]` | occulted body's set time, if within [1st,4th] contact window |
| `tret[7]`, `tret[8]`, `tret[9]` | unused by this worker (zeroed at entry, swecl.c:2441–2442) |

Return-flag bits: `SE_ECL_TOTAL/ANNULAR/PARTIAL` (type), `SE_ECL_VISIBLE` +
`SE_ECL_{MAX,1ST,2ND,3RD,4TH}_VISIBLE` (per-contact visibility from the given geopos),
`SE_ECL_OCC_BEG_DAYLIGHT`/`SE_ECL_OCC_END_DAYLIGHT` (occultation-specific, § step 12),
`SE_ECL_NONCENTRAL` (OR'd in by the public wrapper post-hoc from a fresh `eclipse_where`
call, § "Public wrapper" step 5) — note `SE_ECL_CENTRAL` is *not* similarly OR'd in by the
wrapper (only `NONCENTRAL` is masked in at swecl.c:2095 — same asymmetry exists in the solar
wrapper, swecl.c:2038, so this is not occultation-specific, just noting it's easy to miss
when porting attr/retflag fill-in).
