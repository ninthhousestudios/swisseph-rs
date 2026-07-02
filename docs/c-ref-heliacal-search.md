# C Reference: Heliacal — Event Search (swehel.c part 3)

Porting reference for the top-level heliacal-event search machinery in `swehel.c`: the
oblique-ascension / conjunction bracketing that locates a candidate date, the day-level and
minute-level visibility search that refines it, and the two alternative search strategies
(visibility-limit-based vs. arcus-visionis-based) that `swe_heliacal_ut` dispatches between.
All line numbers refer to `swehel.c` unless stated otherwise.

This doc covers lines 2077–3511 (end of file). It depends on, but does not re-document the
internals of, `swe_vis_limit_mag` (1464), `ObjectLoc` (683), `Magnitude` (1106), `RiseSet` (535,
via `my_rise_trans`/`calc_rise_and_set`), `DeterObject` (305), `default_heliacal_parameters`
(1324), `swe_heliacal_pheno_ut` (1862), `HeliacalAngle` (1636) and `VisLimMagn` — those belong to
sibling "part 1/2" ref docs. Only their call contracts are given here, as needed to understand
the control flow in this file.

---

## Call contracts of dependencies (from earlier in swehel.c)

- **`DeterObject(char *ObjectName) -> int32`** (305–336): lower-cases a copy of the name and
  matches prefixes: `"sun"→SE_SUN`, `"venus"→SE_VENUS`, `"mars"→SE_MARS`, `"mercur"→SE_MERCURY`,
  `"jupiter"→SE_JUPITER`, `"saturn"→SE_SATURN`, `"uranus"→SE_URANUS`, `"neptun"→SE_NEPTUNE`,
  `"moon"→SE_MOON`; else if the string parses as a positive integer, returns
  `SE_AST_OFFSET + atoi(s)` (asteroid number); else returns **`-1`** (meaning "fixed star" — the
  `ObjectName` string itself is the star name, used directly with `swe_fixstar`/`call_swe_fixstar`).
  This `-1` sentinel for "it's a star, not a planet number" is used pervasively throughout this
  file as a branch condition (`if (Planet != -1) … else … call_swe_fixstar(...)`).
- **`call_swe_fixstar(star, tjd, iflag, xx, serr)`** (368–375): thin wrapper around
  `swe_fixstar` that copies the name into a local buffer first (guards against `swe_fixstar`
  overwriting/rewriting the caller's string in place).
- **`my_rise_trans(tjd, ipl, starname, eventtype, helflag, dgeo, datm, tret, serr)`** (508–523):
  the rise/set entry point used throughout this module. If `starname` non-empty, `ipl =
  DeterObject(starname)` first. Then: if `ipl != -1 && |dgeo[1]| < 63°`, use the fast in-module
  `calc_rise_and_set` (not covered here — part 1/2 doc); otherwise (fixed stars, or planets at
  `|lat| >= 63°`) fall through to the public `call_swe_rise_trans` wrapper around
  `swe_rise_trans`/`swe_rise_trans_true_hor` (documented in `docs/c-ref-riseset.md`). Same
  `OK`/`ERR`/`-2` return convention as `swe_rise_trans` (§5.6 of that doc): **`-2` = body does not
  rise/set** (e.g. sun circumpolar at high latitude) — this is threaded through the entire
  heliacal search as a "try a different day" signal, never a hard error.
- **`RiseSet(JDNDaysUT, dgeo, datm, ObjectName, RSEvent, helflag, Rim, tret, serr)`** (535–546):
  `RSEvent` 1=rise/2=set/3=up-transit/4=down-transit; `Rim` 0=disc-center (ORs in
  `SE_BIT_DISC_CENTER`), 1=disc-top(default limb). Resolves `ObjectName` via `DeterObject`, then
  calls `my_rise_trans`.
- **`ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, Angle, helflag, dret, serr) -> OK/ERR`**
  (683–726): computes one scalar angle for `ObjectName` at `JDNDaysUT`. `Angle` code: `0`=topocentric
  apparent altitude (via `swe_azalt`), `1`=azimuth (from-north convention: `xaz[0]+180`, wrapped),
  `2`=topocentric declination, `3`=topocentric RA, `4`=apparent-altitude via the cheaper
  `AppAltfromTopoAlt` Newton approximation instead of a second `swe_azalt` refraction call,
  `5`=geocentric declination, `6`=geocentric RA, `7`=treated identically to `0` (alias, see line
  694 `if (Angle==7) Angle=0`, but importantly `Angle<5` still gates `SEFLG_TOPOCTR` — since 7
  is not `<5`, Angle=7 gets **non-topocentric** position but the **same altitude computation
  branch** as Angle=0; this is a real behavioral distinction from Angle=0, not just an alias —
  used exactly once, `get_acronychal_day` via unused code, not exercised on the vis_lim/arc_vis
  paths documented here). Applies `SEFLG_NONUT|SEFLG_TRUEPOS` unless
  `SE_HELFLAG_HIGH_PRECISION`.
- **`Magnitude(JDNDaysUT, dgeo, ObjectName, helflag, dmag, serr) -> OK/ERR`** (1106–1127): for
  planets, `swe_pheno_ut(...)[4]` (apparent visual magnitude from phenomena calc); for stars,
  `call_swe_fixstar_mag`. Calls `swe_set_topo` internally (STATEFUL — see §8).
- **`swe_vis_limit_mag(tjdut, dgeo, datm, dobs, ObjectName, helflag, dret, serr) -> retval`**
  (1464–1541 area): returns **`-1`**=Error, **`-2`**=object below local horizon (`dret[0] = -100`
  sentinel in that case), **`0`**=OK photopic, bit **`SE_SCOTOPIC_FLAG`(1)** set = OK scotopic,
  bit **`SE_MIXEDOPIC_FLAG`(2)** set = OK but near the photopic/scotopic transition (an
  additional warning bit, not a separate code — check via `retval & SE_MIXEDOPIC_FLAG`, since
  `2` alone never appears without also having `0` or `1` — actually per the header comment the
  encoding is `-1`/`-2` are sentinel returns, otherwise the low bit is `SE_SCOTOPIC_FLAG` and
  another bit is `SE_MIXEDOPIC_FLAG`, i.e. `retval` is a small non-negative flag word when
  `>= 0`). `dret[0]` = visual limiting-magnitude-vs-object delta *is not itself dret[0]* — actual
  layout: `dret[0]` = `VisLimMagn(...)` (the **limiting magnitude** the object would need to
  reach to be visible — NOT a delta), `dret[1..6]` = AltO/AziO/AltS/AziS/AltM/AziM, `dret[7]` =
  object's actual magnitude (from `Magnitude(...)`). **The visibility test used everywhere in
  this file is `darr[0] - darr[7] > 0`** (limiting magnitude minus actual magnitude — since
  smaller/more-negative magnitude is brighter, "actual is brighter than the limit" means
  `darr[7] < darr[0]`, i.e. `darr[0]-darr[7] > 0`). If `AltO < 0` (object below horizon): returns
  `-2` immediately without computing `dret[0]` (only `dret[0] = -100` and `AltO/AziO` etc. are
  left unset/zeroed) — so **do not read `darr[1..7]` when the return is `-2`**.
- **`HeliacalAngle(Magn, dobs, AziO, AltM, AziM, JDNDaysUT, AziS, dgeo, datm, helflag, dangret,
  serr) -> OK/ERR`** (1636): computes the arcus-visionis-method visibility angle;
  `dangret[1]` = required arc of vision (`ArcusVis` in callers), `dangret[2]` = extinction
  coefficient / effective Sun depth used (`ArcusVisPto` in callers, becomes the next iteration's
  `sunsangle` seed — a self-adjusting feedback, see §6).
- **`default_heliacal_parameters(datm, dgeo, dobs, helflag)`** (1324): fills in default
  atmosphere/observer values (pressure/temp from altitude, age=36, SN=1, etc.) for any zero/unset
  input fields — pure value-filling, no reliance on the search logic in this doc.
- **`DeterTAV(dobs, JDNDaysUT, dgeo, datm, ObjectName, helflag, dret, serr) -> OK/ERR`** (1759):
  "Time of Arcus Visionis" — computes, via `HeliacalAngle`, how far the object's actual
  altitude-difference-from-sun (`DeltaAlt = AltO - AziS`... actually `AltO` minus Sun's) is from
  the arc-of-vision threshold at a given instant; `*dret` is the resulting scalar
  (`ArcusVisDelta`-like quantity) that the Moon arc-vis walkthrough (§6) minimizes/roots via
  parabola-fit (`x2min`).
- **`x2min(A, B, C) -> double`** (1791): parabola-vertex-location helper for 3
  **equally-spaced, unit-step** samples `A` (oldest), `B` (middle), `C` (newest) — returns the
  offset (in units of the sample step, presumably in `[0,1]` or similar) of the extremum,
  structurally analogous to `find_maximum` in `docs/c-ref-riseset.md` §1 but with a different
  normalization (callers do `tjd += (1 - x2min(...)) * step`). Not re-derived here; treat as an
  already-available shared helper when porting (check `docs/codebase-map.md` for an existing
  Rust parabola-vertex helper before adding a second one).
- **`Sgn(x) -> int`** (864): standard sign function, `-1/0/+1`.
- **`swe_heliacal_pheno_ut`** (1862): richer per-instant phenomena report; only referenced here
  in **dead** `#if 0` code (line 2377–2384) that would have overwritten `JDNarcvisUT` with
  `darr[13]` from this function — never compiled, skip.

---

## §1 Synodic periods & conjunction table (`get_synodic_period` 2095, `tcon[]` 2566, `find_conjunct_sun` 2579)

### `get_synodic_period(int Planet) -> double` — 2095–2111
Static table of mean synodic periods in days (source: Kelley/Milone/Aveni, *Exploring Ancient
Skies*, p. 43), by `switch(Planet)`:

| Planet | Days |
|---|---|
| `SE_MOON` | 29.530588853 |
| `SE_MERCURY` | 115.8775 |
| `SE_VENUS` | 583.9214 |
| `SE_MARS` | 779.9361 |
| `SE_JUPITER` | 398.8840 |
| `SE_SATURN` | 378.0919 |
| `SE_URANUS` | 369.6560 |
| `SE_NEPTUNE` | 367.4867 |
| `SE_PLUTO` | 366.7207 |
| *default (unmatched, incl. stars and any Planet not listed)* | `366` |

Transcribe these constants exactly (bitwise-exact literals expected).

### `tcon[]` — static const double array, 2566–2577
```c
static const double tcon[] =
{
  0, 0,
  2451550, 2451550,  /* Moon */
  2451604, 2451670,  /* Mercury */
  2451980, 2452280,  /* Venus */
  2451727, 2452074,  /* Mars */
  2451673, 2451877,  /* Jupiter */
  2451675, 2451868,  /* Saturn */
  2451581, 2451768,  /* Uranus */
  2451568, 2451753,  /* Neptune */
};
```
18 doubles = 9 pairs, indexed by **body** (`SE_SUN`=0 .. `SE_NEPTUNE`=8, each body occupying
2 consecutive slots) × 2 (event-type parity, see below). Each pair is a reference Julian Day of
a known conjunction/opposition near J2000, used as the seed epoch for a synodic-period-stepped
search. **There is no `SE_PLUTO` (index 9) row** — the table has only 9 body-pairs (Sun..Neptune).
`get_synodic_period` supports Pluto, but `tcon[]` does not have a slot for it.

### `find_conjunct_sun(tjd_start, ipl, helflag, TypeEvent, *tjd, serr) -> OK/ERR` — 2579–2602
Finds the next Sun–object conjunction (Mercury/Venus, inferior or superior) or
conjunction/opposition (Mars and beyond) after `tjd_start`, via Newton iteration on ecliptic
longitude difference.
1. `daspect = 180` if `ipl >= SE_MARS && TypeEvent >= 3` (i.e. an *opposition*, not conjunction,
   is wanted — outer planets' morning-last/evening-first-type acronychal events pass
   `TypeEvent>=3` here); else `daspect = 0`.
2. **Table index**: `i = (TypeEvent - 1) / 2 + ipl * 2` (integer division — `TypeEvent` 1 or 2 →
   `0`; `TypeEvent` 3 or 4 → `1`). This selects one of the two `tcon[]` slots per body (the two
   slots per body appear to encode two different reference epochs — e.g. for Mercury,
   inferior vs. superior conjunction reference dates — selected by the even/odd-pair split of
   `TypeEvent`).
   **Bounds hazard**: for `ipl == SE_PLUTO` (9), `i` ranges `18` or `19` — **out of bounds** for
   the 18-element `tcon[]` array (valid indices `0..17`). This is a latent C bug (reads
   whatever memory follows the array) that is only reachable if a caller passes a Pluto object
   through this function; in practice `swe_heliacal_ut`'s public API restricts most callers to
   the vis_lim path via `get_asc_obl_with_sun` for outer-planet acronychal events (§5), and
   `find_conjunct_sun` for Pluto's conjunction TypeEvent<=2 case is only reached from
   `heliacal_ut_vis_lim`'s `else` branch when `ipl` is **not** SE_MERCURY/SE_VENUS and
   `TypeEvent > 2` is false is NOT the guard — re-check call site (§5) before deciding whether
   Pluto is actually reachable in the shipped API; **flag this explicitly to the implementer**:
   the Rust port must NOT reproduce an out-of-bounds read — either special-case Pluto (extend
   the table with a documented placeholder, e.g. reuse Neptune's dates, or return an explicit
   error) and note the deviation from C.
3. `tjd0 = tcon[i]`; `dsynperiod = get_synodic_period(ipl)`.
4. `tjdcon = tjd0 + (floor((tjd_start - tjd0) / dsynperiod) + 1) * dsynperiod` — the first
   synodic-period multiple of the reference epoch strictly after `tjd_start`.
5. Newton loop, `ds` initialized to `100` (deg), while `ds > 0.5`:
   - `swe_calc(tjdcon, ipl, epheflag|SEFLG_SPEED, x, serr)`, `swe_calc(tjdcon, SE_SUN,
     epheflag|SEFLG_SPEED, xs, serr)`.
   - `ds = swe_degnorm(x[0] - xs[0] - daspect)`; if `ds > 180`, `ds -= 360` (wrap to
     `(-180,180]`).
   - `tjdcon -= ds / (x[3] - xs[3])` (Newton step using the **speed difference** `x[3]-xs[3]`,
     i.e. ecliptic-longitude-speed of object minus Sun, degrees/day).
   - No iteration cap — an unbounded `while`; relies on quadratic Newton convergence from a
     good-enough seed (`tcon[]` epoch + integer synodic periods). Port as a `loop` with the same
     structure; consider a safety iteration cap for robustness even though C has none (note the
     deviation if added).
6. `*tjd = tjdcon; return OK`. Only `ERR` (propagated from `swe_calc`) or `OK` — no `-2`.

---

## §2 Oblique-ascension machinery (`get_asc_obl` 2452, `get_asc_obl_diff` 2519, `get_asc_obl_with_sun` 2604, `get_asc_obl_acronychal` 2717; dead `_old` variants noted)

### `get_asc_obl(tjd, ipl, star, iflag, dgeo, desc_obl, *daop, serr) -> OK/ERR/-2` — 2452–2484
Computes the **ascensio obliqua** (oblique ascension, or descension) of a body at `tjd` for
observer latitude `dgeo[1]`.
1. Position: `ipl == -1` → `swe_fixstar(star, tjd, epheflag|SEFLG_EQUATORIAL, x, serr)`
   (**note**: calls `swe_fixstar` directly, not the `call_swe_fixstar` string-safety wrapper —
   unlike its `_old` sibling and unlike other callers in this file); else
   `swe_calc(tjd, ipl, epheflag|SEFLG_EQUATORIAL, x, serr)`. `epheflag = iflag &
   (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)` (all other bits of the passed-in `iflag`, e.g. any
   accidental `SEFLG_TOPOCTR`, are dropped before the `swe_calc`/`swe_fixstar` call — always
   geocentric equatorial here).
2. `adp = tan(dgeo[1]·DEGTORAD) · tan(x[1]·DEGTORAD)` (`x[1]` = declination) — the standard
   ascensional-difference formula, `sin(adp_angle) = tan(lat)·tan(decl)`.
3. **Circumpolar guard**: `|adp| > 1` → build a name string (`star` if non-empty, else
   `swe_get_planet_name(ipl,...)`), `serr = "%s is circumpolar, cannot calculate heliacal
   event"`, **return `-2`** (not `ERR`).
4. `adp = asin(adp)/DEGTORAD` (ascensional difference in degrees).
5. `desc_obl == TRUE` → `*daop = x[0] + adp` (descensio obliqua); else `*daop = x[0] - adp`
   (ascensio obliqua). `x[0]` = right ascension. `*daop = swe_degnorm(*daop)`.

**Dead code**: `get_asc_obl_old` (`#if 0`, 2486–2517) is byte-for-byte the same algorithm except
it calls the string-safe `call_swe_fixstar` instead of `swe_fixstar` directly, and its parameter
is `int32 iflag` used identically. Compiled out; **do not port**.

### `get_asc_obl_diff(tjd, ipl, star, iflag, dgeo, desc_obl, is_acronychal, *dsunpl, serr) -> OK/ERR/-2` — 2519–2542
Difference between Sun's and body's oblique ascension (used as the root-finding target for
"Sun and object rise/set together").
1. `aosun` via `get_asc_obl(tjd, SE_SUN, "", iflag, dgeo, desc_obl, &aosun, serr)` — propagate
   `ERR`/`-2` immediately (`retval != OK → return retval`, so a circumpolar `-2` from the Sun
   call is possible and passed straight through; not expected in practice but must be handled).
2. **If `is_acronychal`**: flip `desc_obl` (`TRUE↔FALSE`) before computing the body's
   ascension — i.e. use the Sun's *ascending* horizon-crossing convention paired with the body's
   *opposite* (descending) convention, or vice versa, matching the acronychal (opposite-horizon)
   geometry.
3. `aopl` via `get_asc_obl(tjd, ipl, star, iflag, dgeo, desc_obl [possibly flipped], &aopl,
   serr)` — same error propagation.
4. `*dsunpl = swe_degnorm(aosun - aopl)`.
5. If `is_acronychal`: `*dsunpl = swe_degnorm(*dsunpl - 180)` (shift by half a rotation, since
   acronychal = body on the meridian opposite the Sun's rise/set point).
6. `if (*dsunpl > 180) *dsunpl -= 360` (final wrap to `(-180,180]`, applied unconditionally
   after the acronychal branch too).

**Dead code**: `get_asc_obl_diff_old` (`#if 0`, 2544–2560) — same but no `is_acronychal`
parameter/branch and no final `>180` wrap. Not called; skip.

### `get_asc_obl_with_sun(tjd_start, ipl, star, helflag, evtyp, dperiod, dgeo, *tjdret, serr) -> OK/ERR/-2` — 2604–2677
**This is the live function** used by both the heliacal and acronychal branches of
`heliacal_ut_vis_lim` (§5) to bracket-then-bisect the date on which Sun and object have matching
oblique ascension (cosmic rising/setting together). Two-phase search: coarse forward stepping
to find a sign change, then bisection to `1e-5°` precision.

1. **Mode flags from `evtyp`** (`SE_HELIACAL_RISING`=1/`SE_MORNING_FIRST`, `SE_HELIACAL_SETTING`=2/
   `SE_EVENING_LAST`, `SE_EVENING_FIRST`=3, `SE_MORNING_LAST`=4, `SE_ACRONYCHAL_RISING`=5,
   `SE_ACRONYCHAL_SETTING`=6):
   - `desc_obl = TRUE` if `evtyp == SE_EVENING_LAST(2) || evtyp == SE_EVENING_FIRST(3)`.
   - `retro = TRUE` if `evtyp == SE_MORNING_FIRST(1) || evtyp == SE_EVENING_LAST(2)`.
   - `evtyp == SE_ACRONYCHAL_RISING(5)` → also `desc_obl = TRUE` (overrides/duplicates the
     evening-branch assignment above; net effect: rising(1)→FALSE, setting(2)→TRUE,
     ev-first(3)→TRUE, morning-last(4)→FALSE, acro-rising(5)→TRUE, acro-setting(6)→FALSE).
   - `evtyp == SE_ACRONYCHAL_RISING(5) || evtyp == SE_ACRONYCHAL_SETTING(6)` → `is_acronychal =
     TRUE`; additionally if `ipl != SE_MOON`, `retro = TRUE` (net `retro`: 1→T, 2→T, 3→F(default),
     4→F(default), 5→T(non-moon)/F(moon, unchanged from default F), 6→T(non-moon)/F(moon)).
2. Initial sample: `tjd = tjd_start`; `dsunpl_save = -999999999` (sentinel "no previous sample
   yet"); `get_asc_obl_diff(tjd, ipl, star, epheflag, dgeo, desc_obl, is_acronychal, &dsunpl,
   serr)` — **note**: `epheflag` (ephemeris-selector bits only), not the full `helflag`, is
   passed as the `iflag` argument here — `get_asc_obl` internally re-masks to `epheflag` anyway
   (§2 `get_asc_obl` step 1) so this is inert, but replicate the narrower mask for clarity/fidelity
   if the Rust port threads flags explicitly. Propagate `ERR`/`-2`.
3. **Coarse forward search**, `daystep = 20` days, loop (cap `i <= 5000`, else `ERR` "loop in
   get_asc_obl_with_sun() (1)"):
   ```
   while dsunpl_save == -999999999
      || fabs(dsunpl) + fabs(dsunpl_save) > 180
      || (retro  && !(dsunpl_save < 0 && dsunpl >= 0))
      || (!retro && !(dsunpl_save >= 0 && dsunpl < 0))
   ```
   i.e. keep stepping forward by exactly `+10.0` days each iteration (note: **fixed +10-day
   step, not `daystep`** — `daystep=20` is declared but only used later for the bisection
   bracket width; the coarse-search step is the literal `10.0` at line 2643) until:
   (a) at least one prior sample exists, **and**
   (b) `dsunpl` and `dsunpl_save` are not both near ±180° (`fabs+fabs > 180` guards against a
   spurious "sign change" caused by the ±180° wrap rather than an actual crossing), **and**
   (c) the transition matches the wanted direction: for `retro` events, require the *previous*
   sample negative and the *current* sample non-negative (upward zero-crossing); for
   non-`retro`, require the opposite (downward crossing).
   Each iteration: `dsunpl_save = dsunpl` (save before advancing), `tjd += 10.0`; if `dperiod >
   0 && tjd - tjd_start > dperiod` → **return `-2`** (search-period exceeded, only enforced when
   caller passes a nonzero `dperiod` — both call sites in §5 pass `dperiod = 0`, so this cap is
   presently inert in the shipped vis_lim path, but must still be ported faithfully since it's
   part of the function's public contract). Recompute `dsunpl` via `get_asc_obl_diff`. Propagate
   `ERR`/`-2`.
4. **Bisection**, once the coarse loop exits with a bracketing pair
   `(dsunpl_save at tjd-10, dsunpl at tjd)`:
   - `tjd_start = tjd - daystep` (note: re-anchors using the **20-day** `daystep`, not the 10-day
     coarse-step — i.e. the bracket re-established for bisection is `[tjd-20, tjd]`, twice the
     width of the last coarse step actually taken; this is intentional headroom, not a bug —
     replicate exactly).
   - `daystep /= 2` (→ 10); `tjd = tjd_start + daystep` (midpoint of the 20-day bracket);
     evaluate `dsunpl_test` there.
   - Loop (cap `i <= 5000`, else `ERR` "loop in get_asc_obl_with_sun() (2)") while
     `fabs(dsunpl) > 0.00001`:
     - If `dsunpl_save * dsunpl_test >= 0` (midpoint same sign as the "save" endpoint, i.e. the
       zero is in the *other* half): `dsunpl_save = dsunpl_test; tjd_start = tjd` (move the
       "save" endpoint up to the midpoint).
     - Else: `dsunpl = dsunpl_test` (replace the "current" endpoint value — **note**: `tjd` is
       *not* reassigned in this branch, so the *next* midpoint calculation implicitly keeps the
       old `tjd_start` and just halves `daystep` again from the same `tjd_start` — this is a
       standard bisection where only one of the two bracket values is a live variable
       (`dsunpl`) and the "right edge" `tjd` itself is never separately tracked, only ever
       recomputed as `tjd_start + daystep`).
     - `daystep /= 2.0`; `tjd = tjd_start + daystep`; recompute `dsunpl_test` via
       `get_asc_obl_diff`. Propagate `ERR`/`-2`.
   - Terminates when `fabs(dsunpl) <= 0.00001` (degrees) — **not** a fixed iteration count
     (unlike the riseset bisections); this is a data-dependent convergence loop. Since the
     bracket starts at 20 days and halves each round, ~21 iterations would be needed to reach
     sub-second time resolution, but termination is on the *value* `dsunpl`, not on `daystep`
     magnitude, so the actual iteration count is convergence-dependent (typically fast since
     ascensional-difference varies smoothly, but there's no guaranteed bound beyond the 5000 cap).
5. `*tjdret = tjd; return OK`.

**Dead code**: `get_asc_obl_with_sun_old` (`#if 0`, 2679–2713) — comment "works only for fixed
stars"; older halving-only algorithm (`while (dsunpl < 359.99999) { daystep/=2; ... }`), no
`is_acronychal`/`retro` distinction, no ERR-propagation loop caps, uses a fixed `dsynperiod=367`
day step. Superseded; **do not port**.

**Dead code**: `get_asc_obl_acronychal` (`#if 0`, 2715–2759) — comment "works only for fixed
stars"; a dedicated acronychal-only version using raw `get_asc_obl` calls with `sun_desc`/
`obj_desc` flags fixed by `TypeEvent==4`. Fully superseded by `get_asc_obl_with_sun`'s
`is_acronychal` branch. **Do not port.**

---

## §3 Day-level search (`get_heliacal_day` 2762, `get_acronychal_day` 3043)

### `get_heliacal_day(tjd, dgeo, datm, dobs, ObjectName, helflag, TypeEvent, *thel, serr) -> OK/ERR/-2` — 2762–2921
Starting from a date near conjunction/cosmic-rise (from §1/§2), finds the day and then the
minute at which the object first becomes visible (or, for "last" events, last remains visible)
at sunrise/sunset. This is the **day-stepping search** shared by both Mercury/Venus heliacal
events and (via `heliacal_ut_vis_lim`) all vis_lim-path objects.

1. **Per-`TypeEvent` direction table**:

   | `TypeEvent` | meaning | `is_rise_or_set` | `direct_day` | `direct_time` |
   |---|---|---|---|---|
   | 1 | morning first | `SE_CALC_RISE` | `+1` | `-1` |
   | 2 | evening last | `SE_CALC_SET` | `-1` | `+1` |
   | 3 | evening first | `SE_CALC_SET` | `+1` | `+1` |
   | 4 | morning last | `SE_CALC_RISE` | `-1` | `-1` |

   `is_rise_or_set` selects which Sun event (rise or set) anchors each day's sample;
   `direct_day` is the day-stepping direction (search forward or backward in time from the
   seed); `direct_time` is the minute-level direction used later when refining within a day.

2. **Per-body day-step tuning** (`switch(ipl)`, `ipl = DeterObject(ObjectName)`):

   | Body | `ndays` | pre-adjust to `tjd` | `daystep` | `tfac` | extra |
   |---|---|---|---|---|---|
   | `SE_MOON` | 16 | — | 1 | 1 | |
   | `SE_MERCURY` | 60 | `tjd -= 0·direct_day` (no-op) | 5 | 5 | |
   | `SE_VENUS` | 300 | `tjd -= 30·direct_day` | 5 (15 if `TypeEvent>=3`) | 1 (3 if `TypeEvent>=3`) | |
   | `SE_MARS` | 400 | — | 15 | 5 | |
   | `SE_SATURN` | 300 | — | 20 | 5 | |
   | fixed star (`ipl==-1`) | 300 | — | 15 | 10 (3 if mag<0) | needs `call_swe_fixstar_mag` first (`ERR` on failure); `dmag>2` also sets `daystep=15` (redundant — already 15) |
   | default (Jupiter, Uranus, Neptune, Pluto, asteroids) | 300 | — | 15 | 3 | |

   `tend = tjd + ndays · direct_day`. `retval_old = -2` (init sentinel, "previous day's sun
   rise/set was not found").

3. **Outer day loop**: `for (tday=tjd, i=0; (direct_day>0 && tday<tend) || (direct_day<0 &&
   tday>tend); tday += daystep·direct_day, i++)`:
   - `vdelta = -100` (reset each iteration).
   - **After the first iteration** (`i>0`): `tday -= 0.3·direct_day` — a small backward nudge
     each day (except the very first) to avoid missing an event just before the nominal sample
     point (overlap between consecutive day-samples).
   - `my_rise_trans(tday, SE_SUN, "", is_rise_or_set, helflag, dgeo, datm, &tret, serr)`:
     - `ERR` → propagate.
     - `-2` (sun doesn't rise/set that day, e.g. polar) → `retval_old = -2; continue` (skip to
       next day without evaluating visibility).
   - `swe_vis_limit_mag(tret, dgeo, datm, dobs, ObjectName, helflag, darr, serr)` at the Sun's
     rise/set time `tret`. `ERR` → propagate.
   - **Daystep-shrink-on-first-appearance** (compiled, `#if 1` at 2857): if `retval_old == -2 &&
     retval >= 0 && daystep > 1` (object just transitioned from "not evaluated" to "evaluated
     OK", meaning we've crossed from before-conjunction into after-conjunction visibility
     window): `retval_old = retval; tday -= daystep·direct_day` (roll back one full day-step);
     `daystep = 1` — **except** if `ipl >= SE_MARS || ipl == -1`, in which case `daystep = 5`
     instead of `1` (outer planets/stars: comment notes Mars morning-last periods can be brief
     at high latitude, so 5-day resolution is a deliberate compromise, not 1-day); `continue`
     (retry from the rolled-back day at the finer step, without evaluating visibility this
     iteration). Then unconditionally `retval_old = retval` (this line runs on every iteration
     that reaches it, i.e. even when the shrink branch wasn't taken).
   - If `retval == -2` (object below horizon at Sun's rise/set) → `continue` (try next day).
   - `vdelta = darr[0] - darr[7]` (limiting-mag minus actual mag; see §"call contracts" note on
     `swe_vis_limit_mag`).
   - **Minute-level refinement within the day**: `div = 1440.0` (minutes/day). While `retval !=
     -2 && (vd = darr[0]-darr[7]) < 0` (object not yet/no-longer visible at the current probe
     time): `visible_at_sunsetrise = 0`; step `tret` further into the night/day using an
     **adaptive step size keyed to how far below the visibility threshold `vd` is**:
     - `vd < -1.0` → step `5.0/div · direct_time · tfac`
     - `-1.0 <= vd < -0.5` → step `2.0/div · direct_time · tfac`
     - `-0.5 <= vd < -0.1` → step `1.0/div · direct_time · tfac`
     - else (`-0.1 <= vd < 0`) → step `1.0/div · direct_time` (**no `tfac` multiplier** in the
       finest bracket — replicate this asymmetry exactly).
     Recompute `swe_vis_limit_mag` at the new `tret` each step; `ERR` → propagate. This is an
     **unbounded** `while` (no iteration cap) — relies on the visibility function being
     monotonic enough in the search direction to terminate; port as a `loop`, flag if adding a
     safety cap.
   - **Sunset/sunrise-instant edge nudge** (`visible_at_sunsetrise` flag, set to `1` above only
     if the minute-refinement `while` body never executed, i.e. the object was *already*
     visible right at the Sun's rise/set instant): loop `for i in 0..10`: probe
     `swe_vis_limit_mag(tret + 1.0/div·direct_time, ...)`; if `retval>=0 && darr[0]-darr[7] >
     vd` (visibility margin still improving one minute further out): `vd = darr[0]-darr[7];
     tret += 1.0/div·direct_time` (accept the step). Comment: "vis_limit_mag() has strange
     behaviour" right at sunset/sunrise — this is a documented empirical workaround, not a
     principled correction; port it as-is (fixed 10-iteration greedy hill-climb by 1-minute
     steps).
   - `vdelta = darr[0] - darr[7]` (recomputed after all the above).
   - **Acceptance**: if `vdelta > 0` (object visible): if `(ipl >= SE_MARS || ipl == -1) &&
     daystep > 1` — i.e. outer-planet/star day-search hadn't yet dropped to fine resolution —
     roll back one day-step and set `daystep = 1`, `continue` (redo this day-window at finer
     resolution rather than accepting immediately); else `*thel = tret; return OK`.
4. If the outer loop exhausts `tend` without an accepted return: `serr = "heliacal event does
   not happen"`, **return `-2`**.

### `get_acronychal_day(tjd, dgeo, datm, dobs, ObjectName, helflag, TypeEvent, *thel, serr) -> OK/ERR` — 3043–3105
Refines the acronychal (opposite-Sun) rising/setting day/time, given a seed `tjd` from
`get_asc_obl_with_sun` (called with `is_acronychal=TRUE`).
1. `helflag |= SE_HELFLAG_VISLIM_PHOTOPIC` (force photopic vision model for this search).
2. `TypeEvent == 3 || TypeEvent == 5` → `is_rise_or_set = SE_CALC_RISE`, `direct = -1`; else
   (`TypeEvent` 4 or 6) → `is_rise_or_set = SE_CALC_SET`, `direct = 1`. (Commented-out dead code
   at 3053–3055/3059–3062 shows an abandoned alternative `tret = tjdc ± 3` seed adjustment —
   ignore, not compiled.)
3. `dtret = 999`; loop while `fabs(dtret) > 0.5/1440.0` (i.e. converge to within 0.5 **minutes**
   — note the `#if 0` alternate threshold `0.5` (days) at line 3066 is dead, the live threshold
   is the `#else` branch, minutes):
   - `tjd += 0.7·direct`; if `direct < 0`, additionally `tjd -= 1` (net step `-1.7` days when
     searching backward vs. `+0.7` forward — asymmetric, replicate exactly).
   - `my_rise_trans(tjd, ipl, ObjectName, is_rise_or_set, helflag, dgeo, datm, &tjd, serr)` —
     **note**: `tjd` is both input and output parameter here (passed as both `tjd` arg and
     `&tjd` out-param) — the rise/set time overwrites the seed for the next iteration. `ERR` →
     propagate (**no `-2` handling** — if the object doesn't rise/set, this bubbles up as
     whatever `my_rise_trans` puts in `serr`/return, but the code only checks `== ERR`, not
     `== -2`; a `-2` here falls through and `tjd` is left as `my_rise_trans`'s undefined-in-that-
     case output, a latent gap worth flagging to the implementer as another place to decide the
     Rust port's error-handling stance).
   - `swe_vis_limit_mag(tjd, dgeo, datm, dobs, ObjectName, helflag, darr, serr)`; `ERR` →
     propagate. While `darr[0] < darr[7]` (object **not** visible — limiting mag below actual
     mag): `tjd += 10.0/1440.0 · -direct` (step 10 minutes in the *opposite* of the outer
     direction — walking back toward the rise/set instant until visible), recompute
     `swe_vis_limit_mag`. Unbounded inner loop, no cap.
   - `time_limit_invisible(tjd, ..., helflag|SE_HELFLAG_VISLIM_DARK, direct, &tret_dark, serr)`
     (§4) — visibility-limit time under a "totally dark sky" assumption (Moon/twilight
     suppressed via the `VISLIM_DARK` bit — see swephexp.h flag table, §7 below).
   - `time_limit_invisible(tjd, ..., helflag|SE_HELFLAG_VISLIM_NOMOON, direct, &tret, serr)`
     (§4) — same but only suppressing the Moon's contribution (`VISLIM_NOMOON`), keeping actual
     twilight.
   - `dtret = fabs(tret - tret_dark)` (live `#else` branch at 3092; the `#if 0` branch above it,
     3085–3090, computes an angular separation via `azalt_cart`+dot-product instead — dead,
     don't port).
3. After convergence: `azalt_cart(tret, dgeo, datm, "sun", helflag, darr, serr)` (dead-code
   sibling used only for its `darr[1]` = Sun's altitude output); `*thel = tret`. If Sun's
   altitude `darr[1] < -12°`: `serr = "acronychal rising/setting not available, %f"`
   (**warning only** — still returns `OK`, not an error code, despite the message). Else:
   `serr = "solar altitude, %f"` (**also just informational**, always set on the success path —
   note **`get_acronychal_day` always writes to `serr` on success**, unlike most functions in
   this file which leave `serr` empty on `OK`; the Rust port's error/warning channel should
   surface this as a non-fatal diagnostic, not suppress it).

---

## §4 Visibility timing (`time_optimum_visibility` 2923, `time_limit_invisible` 3000, `get_heliacal_details` 3107)

### `time_optimum_visibility(tjd, dgeo, datm, dobs, ObjectName, helflag, *tret, serr) -> OK/ERR/-2` — 2923–2998
Finds the local-in-time optimum (maximum visibility margin `darr[0]-darr[7]`, i.e. object
brightest relative to the limiting magnitude) near `tjd`, via two independent
"hill-climb-then-shrink-step" searches — one stepping forward in time (`t1`), one backward
(`t2`) — then picks whichever found the larger margin.
1. Seed evaluation at `tjd` itself: `swe_vis_limit_mag(tjd, ...)`; `ERR` → propagate.
   `retval_sv = retval` (saved return-flags word); `phot_scot_opic_sv = retval &
   SE_SCOTOPIC_FLAG` (saved photopic/scotopic bit). `t1 = t2 = tjd`; `vl1 = vl2 = -1`
   (sentinels — any real margin found will exceed the object being below-horizon's effective
   `-1`... actually just an initial "no improvement yet" floor, not derived from the visibility
   scale directly, but functions as a lower bound since valid margins compared are `darr[0]-darr[7]`
   values which start from a real evaluation and only replace `vl1`/`vl2` when strictly greater).
2. **Forward hill-climb** (`t1`): for `i=0, d=100/86400` (100 seconds); `i<3`; `i++, d/=10` (so
   `d` = 100s, 10s, 1s across 3 passes):
   - `t1 += d` (tentative step forward).
   - `t_has_changed = 0`; while `swe_vis_limit_mag(t1-d, ...) >= 0 && darr[0] > darr[7] &&
     darr[0]-darr[7] > vl1`: accept the step-back-by-`d`-from-current position
     (`t1 -= d; vl1 = darr[0]-darr[7]; t_has_changed=1`), save `retval_sv`/
     `phot_scot_opic_sv`. (Net effect: `t1` creeps *backward* in `d`-sized steps as long as
     doing so keeps improving the margin — i.e. despite being called the "forward" search, the
     stepping within each resolution pass is actually backward-seeking from one step ahead; see
     next point.)
   - If no improvement was found this pass (`t_has_changed==0`): **revert** the initial `t1 +=
     d` (`t1 -= d`) — so if the very first probe at `t1-d` (== original `tjd`, after undoing the
     add) wasn scoreasn't better, `t1` ends the pass exactly where it started.
   - `ERR` from the last `swe_vis_limit_mag` call inside the while → propagate.
   This nested "step forward by `d`, then walk backward while improving, else undo" pattern
   at 3 successively finer resolutions (100s→10s→1s) is effectively a coarse-to-fine local
   search **towards decreasing time** bounded below by `tjd` conceptually but not clamped
   explicitly — it can walk arbitrarily far if the margin keeps improving at each `d`, though in
   practice bounded by the visibility function's actual shape.
3. **Backward hill-climb** (`t2`): symmetric, using `t2 -= d` / probing `t2+d` / walking `t2 +=
   d` while improving, i.e. mirror-image of step 2 (searches toward increasing time).
4. `tjd = (vl2 > vl1) ? t2 : t1` (pick whichever direction found the bigger margin — ties go to
   `t1`/forward since the comparison is strict `>`). `*tret = tjd`.
5. **Scotopic/photopic-transition abort**: if the last `retval >= 0`: `phot_scot_opic = retval &
   SE_SCOTOPIC_FLAG`; if it differs from `phot_scot_opic_sv` (the vision mode changed between
   the seed evaluation and the final accepted point) → `printf("hallo -2\n")` (**debug leftover,
   do not port the printf**) and **return `-2`**. Also if `retval_sv & SE_MIXEDOPIC_FLAG` (the
   *saved* — i.e. best-found — point was near the photopic/scotopic boundary) → same `-2` +
   debug printf.
6. Else `return OK`.

### `time_limit_invisible(tjd, dgeo, datm, dobs, ObjectName, helflag, direct, *tret, serr) -> OK/ERR/-2` — 3000–3041
Walks time in direction `direct` (`+1`/`-1`) from `tjd`, extending as far as the object remains
visible, at 3 (4 for the Moon) successively finer resolutions — i.e. finds the boundary of the
visibility window in one direction.
1. `d0 = 100/86400` (100s); if `ObjectName == "moon"` (exact `strcmp`, not prefix): `d0 *= 10`
   (1000s start) and `ncnt = 4` (one extra refinement pass) instead of the default `ncnt = 3`
   (Moon moves faster across the sky/vis-limit boundary, needs both a wider net and finer final
   step).
2. Seed: `swe_vis_limit_mag(tjd + 0·direct, ...)` (i.e. just at `tjd`, the `d*direct` term is `0`
   at this point since `d` isn't initialized until the loop below — effectively `swe_vis_limit_mag(tjd,
   ...)`; `ERR` → propagate. `retval_sv`/`phot_scot_opic_sv` saved.
3. For `i=0, d=d0; i<ncnt; i++, d/=10`: while `swe_vis_limit_mag(tjd + d·direct, ...) >= 0 &&
   darr[0] > darr[7]` (still visible at the probe one `d` further out): accept
   (`tjd += d·direct`), save `retval_sv`/`phot_scot_opic_sv`. No revert-if-no-improvement logic
   here (unlike `time_optimum_visibility`) — this is a straightforward greedy boundary walk,
   not a hill-climb, so there's nothing to undo.
4. `*tret = tjd`; **`*serr = '\0'`** — explicitly clears any warning `swe_vis_limit_mag` may
   have set (specifically: object-below-horizon `-2` warnings are expected/routine at a
   visibility boundary and must be suppressed, per the inline comment "if object disappears at
   setting, retval is -2, but we want it OK, and also suppress the warning").
5. Same scotopic/photopic-transition check as `time_optimum_visibility` step 5 (**no debug
   printf here** — this function's version is silent), returning `-2` on a vision-mode change or
   `SE_MIXEDOPIC_FLAG` on the last accepted (`retval_sv`) point; else `OK`.

### `get_heliacal_details(tday, dgeo, datm, dobs, ObjectName, TypeEvent, helflag, *dret, serr) -> OK/ERR` — 3107–3161
Produces the `dret[0..2]` = [visibility-start, optimum, visibility-end] triple from a day
already known to contain the event (`tday`, from `get_heliacal_day`/`get_acronychal_day`).
1. `dret[1]` (optimum) via `time_optimum_visibility(tday, ..., &dret[1], serr)`. If it returns
   `-2`: treated as **`retval=OK`, `optimum_undefined=TRUE`** (not a hard error — just marks the
   slot as uncertain due to a photopic/scotopic transition, per that function's contract).
2. `direct = 1`, or `-1` if `TypeEvent == 1 || TypeEvent == 4` (morning first / morning last:
   search backward from `tday` for the start-of-visibility boundary, since morning events are
   found working backward from sunrise).
3. `dret[0]` via `time_limit_invisible(tday, ..., direct, &dret[0], serr)`. `-2` → `retval=OK,
   limit_1_undefined=TRUE`.
4. `direct *= -1` (flip); `dret[2]` via `time_limit_invisible(dret[1], ..., direct, &dret[2],
   serr)` — **note: seeded from `dret[1]` (the optimum), not `tday`** — the second boundary
   walk starts from the optimum-visibility instant, walking the opposite direction from the
   first boundary. `-2` → `retval=OK, limit_2_undefined=TRUE`.
5. **Reorder for evening events**: if `TypeEvent == 2 || TypeEvent == 3` (evening last / evening
   first): swap `dret[0]`↔`dret[2]`, and correspondingly swap the `limit_1_undefined`/
   `limit_2_undefined` flags (via an `int i = (int)limit_1_undefined` temp cast through
   `AS_BOOL`). This ensures `dret[0]` is always chronologically first regardless of which
   physical boundary (`time_limit_invisible` direction) it came from.
6. If any of the three flags is set: build a diagnostic `serr` string listing which of `"0,"
   "1," "2,"` are uncertain, suffixed `"] are uncertain due to change between photopic and
   scotopic vision"` — **still returns `OK`** (this is a warning annotation on an otherwise
   successful result, matching the `Warning`-via-`flags_used`/comparable mechanism the Rust port
   uses instead of a separate error type — see project `CLAUDE.md` API pattern).

---

## §5 Event drivers — vis_lim path (`heliacal_ut_vis_lim` 3163, `moon_event_vis_lim` 3249, `MoonEventJDut` 3327)

This is the **default** search strategy (used unless an `SE_HELFLAG_AVKIND_*` bit is set — see
§7 flag table). It combines §1–§4: bracket a candidate date via conjunction/oblique-ascension,
then day-step to visibility, then (optionally) refine start/optimum/end times.

### `heliacal_ut_vis_lim(tjd_start, dgeo, datm, dobs, ObjectName, TypeEventIn, helflag, *dret, serr_ret) -> OK/ERR/-2` — 3163–3246
1. `dret[0..9] = 0`; `*dret = tjd_start` (error-case fallback value); `ipl = DeterObject(ObjectName)`.
2. Seed `tjd = tjd_start - 30` (Mercury) or `tjd_start - 50` (everything else) — pulled back far
   enough that the subsequent forward search can't miss an event that straddles `tjd_start`, at
   the cost of possibly returning an event *before* `tjd_start` (handled by the retry loop in
   `swe_heliacal_ut`, §7).
3. `helflag2 = helflag` (a local copy; the commented-out `&= ~SE_HELFLAG_HIGH_PRECISION` at 3182
   is dead — high precision is *not* stripped for the day-search sub-call, contrary to what an
   older version apparently did).
4. **Branch: heliacal vs. acronychal** — `ipl == SE_MERCURY || ipl == SE_VENUS || TypeEvent <=
   2` selects the **heliacal** branch (true for morning-first/evening-last always, and for
   evening-first/morning-last *only* on Mercury/Venus — outer planets' evening-first/
   morning-last are, per the physics, actually *acronychal* events and go to the `else` branch
   even though `TypeEvent` is nominally 3/4):
   - **Heliacal branch**:
     - `ipl == -1` (star): `get_asc_obl_with_sun(tjd, ipl, ObjectName, helflag, TypeEvent, 0,
       dgeo, &tjd, serr)` (§2) — find the cosmic-rise/set date. Propagate `ERR`/`-2`.
     - else (planet): `find_conjunct_sun(tjd, ipl, helflag, TypeEvent, &tjd, serr)` (§1) — find
       the conjunction date. `ERR` → propagate (note: **`-2` is not a documented return of
       `find_conjunct_sun`**, so only `ERR` is checked here — consistent with §1's analysis
       that it only returns `OK`/`ERR`).
     - `get_heliacal_day(tjd, dgeo, datm, dobs, ObjectName, helflag2, TypeEvent, &tday, serr)`
       (§3). Propagate any non-`OK`.
   - **Acronychal branch** (`else`, i.e. outer planets/stars with `TypeEvent` 3 or 4):
     - Always (the `if (1 || ipl == -1)` at 3206 is dead-condition-simplified — always true):
       `get_asc_obl_with_sun(tjd, ipl, ObjectName, helflag, TypeEvent, 0, dgeo, &tjd, serr)` —
       **same call as the star sub-branch above**, now used for planets too (the commented-out
       `get_asc_obl_acronychal` call and the `find_conjunct_sun` `else` branch below it at
       3211–3215 are dead — unreachable given `1 ||`). Propagate `ERR`/`-2`.
     - `tday = tjd`; `get_acronychal_day(tjd, dgeo, datm, dobs, ObjectName, helflag2, TypeEvent,
       &tday, serr)` (§3). Propagate any non-`OK`.
5. `dret[0] = tday`.
6. **Details refinement** (unless `SE_HELFLAG_NO_DETAILS`):
   - If heliacal-branch conditions (`ipl==SE_MERCURY||ipl==SE_VENUS||TypeEvent<=2`):
     `get_heliacal_details(tday, dgeo, datm, dobs, ObjectName, TypeEvent, helflag2, dret, serr)`
     (§4) — fills `dret[0..2]`. `ERR` → propagate (note: this call's own `-2`/warning path
     already resolves to `OK` internally per §4, so only `ERR` needs handling here).
   - Else: dead code (`else if ((0))`, 3230–3240) — an alternate "walk `*dret` while below
     threshold" refinement, permanently disabled. **Do not port** — acronychal events get no
     `dret[1]`/`dret[2]` refinement (they stay `0`, set at step 1).
7. Return `retval` (the last one set — `OK` unless something jumped via `goto
   swe_heliacal_err`, which itself just copies `serr` and returns whatever `retval` was at the
   jump).

### `moon_event_vis_lim(tjdstart, dgeo, datm, dobs, TypeEvent, helflag, *dret, serr_ret) -> OK/ERR/-2` — 3249–3325
Moon-specific vis_lim search — the Moon has no morning-first/evening-last (its heliacal-type
events are only evening-first/morning-last, `TypeEvent` 3/4, plus event synthesis of the
3-point `dret` differently from planets).
1. `TypeEvent == 1 || TypeEvent == 2` → `ERR` "the moon has no morning first or evening last".
2. `helflag2 = helflag & ~SE_HELFLAG_HIGH_PRECISION` (**high precision IS stripped here**,
   unlike the planet path — inconsistent with `heliacal_ut_vis_lim`'s dead-code comment at
   3182, but this one is live).
3. `tjd = tjdstart - 30`; `find_conjunct_sun(tjd, SE_MOON, helflag, TypeEvent, &tjd, serr)` (§1).
   `ERR` → propagate.
4. `get_heliacal_day(tjd, dgeo, datm, dobs, "moon", helflag2, TypeEvent, &tjd, serr)` (§3, using
   the **stripped** `helflag2`). Non-`OK` → `goto moon_event_err`. `dret[0] = tjd`.
5. `time_optimum_visibility(tjd, ..., helflag [full, not helflag2], &tjd, serr)` (§4, full
   `helflag` here). `ERR` → `goto moon_event_err`. `dret[1] = tjd`.
6. `direct = 1`, or `-1` if `TypeEvent == 4`. `time_limit_invisible(tjd, ..., helflag, direct,
   &tjd, serr)` (§4, full `helflag`). `ERR` → `goto moon_event_err`. `dret[2] = tjd`.
7. `direct *= -1`; `time_limit_invisible(dret[1], ..., helflag, direct, &tjd, serr)` — seeded
   from `dret[1]` like `get_heliacal_details` step 4. `dret[0] = tjd`. `ERR` → `goto
   moon_event_err`.
8. **Sunset/sunrise clamp** (compiled, `#if 1` at 3294): if `TypeEvent == 3` (evening first):
   `my_rise_trans(tjd, SE_SUN, "", SE_CALC_SET, helflag, dgeo, datm, &trise, serr)`; if `trise <
   dret[1]` (sunset earlier than the optimum-visibility instant — i.e. Moon was already visible
   before the Sun even set): `dret[0] = trise` (clamp start-of-visibility to sunset, not the
   earlier computed boundary — comment: "if the moon is visible before sunset, we return sunset
   as start time"; warning suppressed, "it happens too often"). Else (`TypeEvent==4`, morning
   last): `my_rise_trans(dret[1], SE_SUN, "", SE_CALC_RISE, ...)`; if `dret[0] > trise`
   (end-of-visibility later than sunrise): `dret[0] = trise` (clamp — comment: "if the moon is
   visible after sunrise, we return sunrise as end time"; **note this reuses `dret[0]` as the
   "end time" slot for `TypeEvent==4`** — see next step's reorder).
9. **Reorder for `TypeEvent==4`**: swap `dret[0]`↔`dret[2]` (morning-last's boundary walk
   directions produce start/end in the opposite slots compared to evening-first; this swap
   normalizes to the same `[start, optimum, end]` convention as the heliacal path).
10. Return `retval`.

### `MoonEventJDut(JDNDaysUTStart, dgeo, datm, dobs, TypeEvent, helflag, *dret, serr) -> OK/ERR/-2` — 3327–3334
Trivial dispatcher: `avkind = helflag & SE_HELFLAG_AVKIND`; if any AV-kind bit set →
`moon_event_arc_vis` (§6); else → `moon_event_vis_lim` (this section). This is the Moon-specific
analogue of `heliacal_ut` (§7's inner dispatcher).

---

## §6 Event drivers — arc_vis path (`heliacal_ut_arc_vis` 2211, `moon_event_arc_vis` 2114, `HeliacalJDut` 2077)

Selected whenever `helflag & SE_HELFLAG_AVKIND` is nonzero (i.e. one of `AVKIND_VR`,
`AVKIND_PTO`, `AVKIND_MIN7`, `AVKIND_MIN9` is set — see §7 flag table). This is the
"arcus visionis" method: instead of asking "is the object's brightness above the sky's limiting
magnitude" (vis_lim path), it asks "is the object's altitude above the Sun's altitude by at
least the empirically/theoretically required arc" — an older, simpler (pre-`VisLimMagn`) model,
still offered for compatibility/comparison.

### `HeliacalJDut(...)` — 2077–2093 (public-ish VB-compat wrapper, `#endif`-gated at top, i.e.
inside some larger `#if`/`#endif` region started earlier in the file)
Thin shim: packs scalar `Age`/`SN`/`Lat`/`Longitude`/`HeightEye`/`Temperature`/`Pressure`/`RH`/
`VR` into `dgeo[3]`/`datm[4]`/`dobs[6]`, forces `helflag = SE_HELFLAG_HIGH_PRECISION |
SE_HELFLAG_AVKIND_VR`, and calls `swe_heliacal_ut`. Purely a legacy VB-interop convenience
wrapper — **not itself part of the Rust port's public API surface** (the Rust port exposes
`swe_heliacal_ut`'s equivalent directly); mention for completeness only.

### `moon_event_arc_vis(JDNDaysUTStart, dgeo, datm, dobs, TypeEvent, helflag, *dret, serr) -> OK/ERR/-2` — 2114–2209
1. `avkind = helflag & SE_HELFLAG_AVKIND`; default to `SE_HELFLAG_AVKIND_VR` if zero. **Only
   `SE_HELFLAG_AVKIND_VR` is supported for the Moon** — any other avkind combination → `ERR`
   "error: in valid AV kind for the moon" (note: this only checks `avkind !=
   SE_HELFLAG_AVKIND_VR` as a single value, i.e. combining `AVKIND_VR` with e.g. `AVKIND_PTO`
   simultaneously would also be rejected since the combined bitmask wouldn't equal the single
   `AVKIND_VR` value).
2. `TypeEvent == 1 || TypeEvent == 2` → `ERR` "error: the moon has no morning first or evening
   last" (same restriction as the vis_lim Moon path).
3. `iflag = SEFLG_TOPOCTR|SEFLG_EQUATORIAL|epheflag`, `+= SEFLG_NONUT|SEFLG_TRUEPOS` unless
   `SE_HELFLAG_HIGH_PRECISION`.
4. `Daystep = 1`; `TypeEvent == 3` (evening first) → `TypeEvent = 2` (remapped to internal
   "set"-type code); else (`TypeEvent==4`, morning last) → `TypeEvent = 1` (remapped to
   "rise"-type), **and** `Daystep = -Daystep` (search backward for morning-last).
5. **New-moon-date determination**: `JDNDaysUT = JDNDaysUTStart` (`+30` if remapped
   `TypeEvent==1`, i.e. the original evening-first request — start the phase search a synodic
   month later to avoid the immediately-preceding new moon). `swe_pheno_ut(JDNDaysUT, SE_MOON,
   iflag, x, serr)` → `phase2 = x[0]` (phase angle / illuminated-fraction proxy — exact `x[]`
   slot semantics belong to the phenomena ref doc). `goingup = 0`; loop: step `JDNDaysUT +=
   Daystep`, `phase1 = phase2`, recompute `phase2`; `phase2 > phase1` → `goingup = 1`; continue
   `while (goingup==0 || (goingup==1 && phase2>phase1))` — i.e. keep stepping while the phase
   value is still increasing (or hasn't started increasing yet), stop the instant it decreases
   after having increased. This locates the date of **minimum phase** (new moon) by walking
   forward until past the trough. `JDNDaysUT -= Daystep` (back up to the day with smallest
   phase).
6. `JDNDaysUTi = JDNDaysUT` (remember this as the "new moon" anchor for a later 15-day
   bailout bound); `JDNDaysUT -= Daystep` (one step before, to re-enter the next loop's `+=` at
   the anchor); `MinTAVoud = 199` (sentinel large value).
7. **Outer day loop** (do/while): `JDNDaysUT += Daystep`; `RiseSet(JDNDaysUT, dgeo, datm,
   "moon", TypeEvent[remapped], helflag, Rim=0, &tjd_moonevent, serr)` (§ call contracts) — any
   non-`OK` return propagated directly (**including `-2`** — `RiseSet`'s return is passed
   straight through as this function's return, unlike most other call sites in this file which
   distinguish `-2` from `ERR`).
   `tjd_moonevent_start = tjd_moonevent` (anchor for a 120-minute bailout inside the inner loop).
   - **Inner per-minute loop** (do/while): `OldestMinTAV = MinTAVoud; MinTAVoud = MinTAV;
     DeltaAltoud = DeltaAlt` (rotate history); `tjd_moonevent -= (1/60/24)·Sgn(Daystep)`
     (step **one minute** against the day-step's sign — i.e. walk backward in time by a minute
     each inner iteration when `Daystep>0`, forward when `Daystep<0`); `ObjectLoc(...,"sun",
     Angle=0,...)`→`AltS`, `ObjectLoc(...,"moon",Angle=0,...)`→`AltO` (both altitude, `ERR`→
     propagate as `ERR` directly, not via `goto`); `DeltaAlt = AltO - AltS`; `DeterTAV(dobs,
     tjd_moonevent, dgeo, datm, "moon", helflag, &MinTAV, serr)` (arc-of-vision delta at this
     instant); `TimeCheck = tjd_moonevent - (LocalMinStep/60/24)·Sgn(Daystep)` (`LocalMinStep` =
     8, swehel.c:86 — look 8 minutes further in the same walking direction);
     `DeterTAV(dobs, TimeCheck, ..., &LocalminCheck, serr)` (peek-ahead check to distinguish a
     true local minimum from a still-descending trend).
     Continue while `(MinTAV <= MinTAVoud || LocalminCheck < MinTAV) && fabs(tjd_moonevent -
     tjd_moonevent_start) < 120/60/24` (i.e. within 120 minutes of the day's rise/set anchor):
     keep walking as long as the TAV metric is still improving (getting smaller) *or* the
     8-minute-ahead peek suggests it will keep improving (guards against stopping at a local
     wiggle rather than the true minimum).
   - Continue outer loop while `DeltaAltoud < MinTAVoud && fabs(JDNDaysUT - JDNDaysUTi) < 15`
     (still hasn't found a day where the altitude-difference at rise/set already exceeds the
     required arc, and still within 15 days of the new-moon anchor).
8. **Acceptance**: if `fabs(JDNDaysUT - JDNDaysUTi) < 15` (found within the 15-day window):
   `tjd_moonevent += (1 - x2min(MinTAV, MinTAVoud, OldestMinTAV))·Sgn(Daystep)/60/24` — parabola-
   vertex refinement (via the shared `x2min` helper, §"call contracts") across the last three
   per-minute TAV samples, converted to a sub-minute time correction. Else: `ERR` "no date found
   for lunar event".
9. `dret[0] = tjd_moonevent; return OK`.

### `heliacal_ut_arc_vis(JDNDaysUTStart, dgeo, datm, dobs, ObjectName, TypeEventIn, helflag, *dret, serr_ret) -> OK/ERR/-2` — 2211–2450
The general (non-Moon) arcus-visionis search: locate the day the Sun reaches the object's
required arc-of-vision depth below the horizon, at the object's own rise/set time, via a
day-stepping outer search (adaptive step, powers-of-two-ish backoff) plus an optional
`AVKIND_VR` per-minute walkthrough and/or `AVKIND_PTO` symmetric-crossing averaging.

1. `*dret = JDNDaysUTStart` (error fallback); `Planet = DeterObject(ObjectName)`.
2. `Magnitude(JDNDaysUTStart, dgeo, ObjectName, helflag, &objectmagn, serr)` — initial magnitude
   estimate (refined later per-iteration for planets). `ERR` → `goto swe_heliacal_err`.
3. `iflag = SEFLG_TOPOCTR|SEFLG_EQUATORIAL|epheflag`, `+=NONUT|TRUEPOS` unless
   `SE_HELFLAG_HIGH_PRECISION`.
4. **Per-body `DayStep`/`maxlength`** (`switch(Planet)`):

   | Planet | `DayStep` | `maxlength` |
   |---|---|---|
   | `SE_MERCURY` | 1 | 100 |
   | `SE_VENUS` | 64 | 384 |
   | `SE_MARS` | 128 | 640 |
   | `SE_JUPITER` | 64 | 384 |
   | `SE_SATURN` | 64 | 256 |
   | default (Uranus/Neptune/Pluto/stars) | 64 | 256 |

   `eventtype = TypeEvent`. `eventtype==2` (evening last) → `DayStep = -DayStep`. `eventtype==4`
   (morning last) → `eventtype = 1` (remap to morning-first-style rise/set target), `DayStep =
   -DayStep` (double negation with the case-2 branch not applying here since they're
   mutually exclusive `if`s — net: for `eventtype==4`, DayStep flips from its base positive
   value once). `eventtype==3` (evening first) → `eventtype = 2`. Finally
   `eventtype |= SE_BIT_DISC_CENTER` (always search disc-center rise/set, not limb).
5. **Outer adaptive day-stepping search** (nested do/while, 2270–2367):
   - `JDNDaysUTfinal = JDNDaysUTStart + maxlength`; `JDNDaysUT = JDNDaysUTStart - 1`; if
     `DayStep < 0`, swap `JDNDaysUT`↔`JDNDaysUTfinal` (search window direction matches the
     step's sign).
   - `JDNDaysUTstep = JDNDaysUT - DayStep` (pre-position one step before, since the inner loop
     adds `DayStep` before its first real evaluation); `doneoneday = 0`; `ArcusVisDelta = 199`
     (sentinel); `ArcusVisPto = -5.55` (initial extinction-adjusted Sun-depth guess, degrees).
   - **Middle loop** (`do { if (fabs(DayStep)==1) doneoneday=1; <inner loop>; ... } while
     (doneoneday==0 && (JDNDaysUTfinal-JDNDaysUTstep)*Sgn(DayStep) > 0)`) — implements a
     **halving backoff**: start with the body's coarse `DayStep`, and once the inner loop
     (below) finds the sign-change bracket, halve `DayStep` (rounding down via `(int)`) and
     retry from the bracket's *older* endpoint, repeating until `DayStep` reaches `±1` (at which
     point `doneoneday` latches `1` and this becomes the terminal pass).
   - **Inner loop** (per day-step, evaluates the Sun's altitude relative to the object's
     rise/set instant):
     - `JDNDaysUTstepoud = JDNDaysUTstep; ArcusVisDeltaoud = ArcusVisDelta;` (save previous);
       `JDNDaysUTstep += DayStep`.
     - `my_rise_trans(JDNDaysUTstep, SE_SUN, "", eventtype, helflag, dgeo, datm, &tret, serr)` —
       **Sun's** own rise/set (disc-center, per `eventtype`'s `SE_BIT_DISC_CENTER`) at this
       candidate day. `ERR` → `goto swe_heliacal_err` (note: **no `-2` handling** here — if the
       Sun doesn't rise/set that day this would propagate an unexpected `-2` as the function's
       own return, since only `== ERR` is checked; flag as an edge case, mirrors the
       `get_acronychal_day` gap noted in §3).
     - `tjd_tt = tret + swe_deltat_ex(...)`; `swe_calc(tjd_tt, SE_SUN, iflag, x, serr)`;
       `swe_azalt(tret, SE_EQU2HOR, dgeo, Pressure, Temperature, xin, xaz)` → Sun's altitude
       `xaz[1]` at its own rise/set instant (should be ≈0° modulo refraction/disc-center
       settings — used as the reference altitude for the hour-angle formula next).
     - `Trise = HourAngle(xaz[1], x[1], dgeo[1])` (hour angle, in **hours**, at which the Sun
       reaches altitude `xaz[1]` given its declination `x[1]` and observer latitude).
     - `sunsangle = ArcusVisPto` (the running Sun-depth target, self-adjusted each outer
       iteration — see below); overridden to `-7` if `SE_HELFLAG_AVKIND_MIN7`, `-9` if
       `SE_HELFLAG_AVKIND_MIN9` (fixed-depth variants of the arc-vis method, independent of the
       `HeliacalAngle`-computed `ArcusVisPto` feedback).
     - `Theliacal = HourAngle(sunsangle, x[1], dgeo[1])` (hour angle at which the Sun reaches
       the *target* depth).
     - `Tdelta = Theliacal - Trise`; negate if `TypeEvent == 2 || TypeEvent == 3` (evening
       events run the other temporal direction relative to rise/set).
     - `JDNarcvisUT = tret - Tdelta/24` (candidate instant: the object's rise/set time, shifted
       by the hour-angle delta converted to days, landing on the moment the Sun is at the
       target depth).
     - Recompute Sun position/az-alt **at `JDNarcvisUT`** (not `tret`): `AziS = xaz[0]+180`
       (wrapped `<360`), `AltS = xaz[1]`.
     - **Dead code** (`#if 0`, 2323–2333): would have also computed the Moon's `AziM`/`AltM` at
       `JDNarcvisUT` — never compiled; the live code always passes `AltM=-1, AziM=0` to
       `HeliacalAngle` (see below), meaning **Moon-position brightening/interference effects
       are never actually factored into this specific angle calc** despite `HeliacalAngle`'s
       signature accepting Moon az/alt — a real behavioral fact, not just a doc gap.
     - Object/star position at `JDNarcvisUT`: `swe_calc`(planet)/`call_swe_fixstar`(star, `ERR`→
       goto); for planets, also refresh `objectmagn` via `Magnitude(JDNarcvisUT, ...)` (**stars
       do NOT get their magnitude refreshed here** — `objectmagn` stays at the value from step 2,
       computed once at `JDNDaysUTStart`, for the entire star search). `AziO`/`AltO` via
       `swe_azalt`. `DeltaAlt = AltO - AltS`.
     - `HeliacalAngle(objectmagn, dobs, AziO, AltM=-1, AziM=0, JDNarcvisUT, AziS, dgeo, datm,
       helflag, dang, serr)` (commented-out call one line above passes real `AltM`/`AziM` — dead,
       confirms the note above). `ArcusVis = dang[1]`; **`ArcusVisPto = dang[2]`** — this is the
       self-adjusting feedback: next outer-loop iteration's `sunsangle` seed comes from this
       iteration's `HeliacalAngle` output, not a fixed constant (except when overridden by
       `MIN7`/`MIN9`).
     - `ArcusVisDelta = DeltaAlt - ArcusVis` (how far the object's actual altitude-above-Sun
       exceeds the theoretically required arc; negative = not yet visible by this criterion).
     - **Inner-loop continuation**: `while ((ArcusVisDeltaoud > 0 || ArcusVisDelta < 0) &&
       (JDNDaysUTfinal - JDNDaysUTstep)·Sgn(DayStep) > 0)` — keep stepping while either the
       *previous* sample was still "past due" (positive, meaning we overshot on the last step
       — re-check needed) or the *current* sample hasn't yet reached the threshold, and we
       haven't exhausted the search window.
   - **Backoff-on-first-bracket** (2361–2366): once the inner loop exits (bracket found or
     window exhausted), if `doneoneday==0` (haven't yet reached unit day-step) **and** window
     not exhausted: revert to the *previous* delta/step (`ArcusVisDelta = ArcusVisDeltaoud;
     JDNDaysUTstep = JDNDaysUTstepoud`), and halve the step: `DayStep = (int)(fabs(DayStep)/2.0)
     · Sgn(DayStep)` — **integer-truncating halving**, so a `DayStep` of 64 → 32 → 16 → 8 → 4 →
     2 → 1 (six halvings from Venus/Jupiter/Saturn/default's 64; Mars's 128 needs seven), at
     which point the next middle-loop pass's `fabs(DayStep)==1` check sets `doneoneday=1` and
     the following inner-loop pass is the final, finest-resolution one.
6. **Window-exhaustion check** (2369–2375): `d = (JDNDaysUTfinal - JDNDaysUTstep)·Sgn(DayStep)`;
   if `d <= 0 || d >= maxlength` (either ran off the end of the search window, or — oddly —
   still has the *entire* window left, which would indicate the very first inner-loop iteration
   already satisfied the exit condition without narrowing anything, an edge/degenerate case):
   `dret[0] = JDNDaysUTinp` (original input, not the search seed `tjd`); `retval = -2`; `serr =
   "heliacal event not found within maxlength %f\n"`; `goto swe_heliacal_err`.
7. **Dead code** (`#if 0`, 2376–2384): would have called `swe_heliacal_pheno_ut` to overwrite
   `JDNarcvisUT` with `darr[13]` if `AVKIND_VR` — never compiled, the live `AVKIND_VR` handling
   is step 8 below, which is a from-scratch walkthrough, not a call into `swe_heliacal_pheno_ut`.
8. `direct = TimeStepDefault/24/60` (`TimeStepDefault`=1, swehel.c:85 → 1 minute in days),
   negated if `DayStep < 0`.
9. **`AVKIND_VR` per-minute walkthrough** (2387–2418, only if `SE_HELFLAG_AVKIND_VR` set —
   this is the *default* avkind per `swe_heliacal_ut`'s handling, see §7): refine
   `JDNarcvisUT` to the exact minute the Time-of-Arc-of-Vision (`DeterTAV`) metric bottoms out,
   via the same `x2min`-parabola-vertex pattern used in `moon_event_arc_vis` §above:
   - `TimeStep = direct`; `TbVR = 0`; `TimePointer = JDNarcvisUT`.
   - `DeterTAV(..., TimePointer, ..., &OldestMinTAV, serr)`; `TimePointer += TimeStep`;
     `DeterTAV(..., TimePointer, ..., &MinTAVoud, serr)`.
   - If `MinTAVoud > OldestMinTAV` (moved the wrong way — TAV got worse): reset `TimePointer =
     JDNarcvisUT`, negate `TimeStep`, `MinTAVact = OldestMinTAV`. Else: `MinTAVact = MinTAVoud;
     MinTAVoud = OldestMinTAV` (rotate so the search direction that was already correct
     continues).
   - Do/while loop: `TimePointer += TimeStep`; rotate `OldestMinTAV←MinTAVoud←MinTAVact`;
     `DeterTAV(..., TimePointer, ..., &MinTAVact, serr)`; if `MinTAVoud < MinTAVact` (just passed
     the minimum — TAV started increasing again): `extrax = x2min(MinTAVact, MinTAVoud,
     OldestMinTAV)`; `TbVR = TimePointer - (1-extrax)·TimeStep` (sub-step parabola-vertex
     correction). Continue `while (TbVR == 0)` (i.e. loop until the vertex is found —
     **unbounded**, no iteration cap, relies on the metric having exactly one minimum in the
     search direction).
   - `JDNarcvisUT = TbVR`.
10. **`AVKIND_PTO` symmetric-crossing averaging** (2420–2438, only if `SE_HELFLAG_AVKIND_PTO`
    set — **mutually possible in combination with step 9**, since these are independent flag
    bits, though in practice callers likely pick one; the code applies both if both bits are
    set, `PTO` operating on whatever `JDNarcvisUT` step 9 left behind):
    - Do/while: `OudeDatum = JDNarcvisUT`; `JDNarcvisUT -= direct`; recompute object position at
      `JDNarcvisUT`, `Angle = xaz[1]` (object's own altitude, not Sun's, at this instant).
      Continue `while (Angle > 0)` — i.e. step backward by `direct` until the object's altitude
      drops to/below 0° (its own horizon crossing).
    - `JDNarcvisUT = (JDNarcvisUT + OudeDatum) / 2.0` — average the last two samples (one just
      below the horizon, one just above), a crude linear-interpolation-via-averaging estimate of
      the object's own horizon-crossing instant near the arc-vis date.
11. **Sanity bound**: `JDNarcvisUT < -9999999 || > 9999999` → `dret[0] = JDNDaysUT` (note: **not
    `JDNarcvisUT`** — uses the loop variable `JDNDaysUT`, which at this point holds whatever
    `swe_deltat_ex`'s last write left it as, since `JDNDaysUT` isn't otherwise assigned after
    the setup block — this looks like a latent bug/typo in the original C reusing a stale
    variable rather than `JDNDaysUTinp`, flag to the implementer rather than silently
    replicating without comment); `serr = "no heliacal date found"`; `retval = ERR`; `goto
    swe_heliacal_err`.
12. `dret[0] = JDNarcvisUT`. Return `retval` (`OK` unless a `goto` set otherwise).

**Note on `dret[]` for the arc_vis path**: unlike the vis_lim path (§5), `heliacal_ut_arc_vis`
only ever populates `dret[0]` — no optimum/end-of-visibility slots (`dret[1]`/`dret[2]` are
never touched by this function, staying at whatever `swe_heliacal_ut`'s caller-visible
initialization left them — see §7's `dret[]` contract for what that means at the top level;
the caller of `heliacal_ut`/`heliacal_ut_arc_vis` — `swe_heliacal_ut` — does **not**
zero-initialize `dret[]` up front the way `heliacal_ut_vis_lim` does, so `dret[1]`/`dret[2]`
are left as whatever the caller passed in, undefined unless the caller pre-zeroed).

---

## §7 Top-level dispatch (`heliacal_ut` 3336, `swe_heliacal_ut` 3385) — TypeEvent validation, object-class routing, flag preprocessing, `dret[]` contract

### `heliacal_ut(JDNDaysUTStart, dgeo, datm, dobs, ObjectName, TypeEventIn, helflag, *dret, serr_ret) -> OK/ERR/-2` — 3336–3343
Trivial dispatcher, non-Moon: `avkind = helflag & SE_HELFLAG_AVKIND`; nonzero → §6
`heliacal_ut_arc_vis`; zero → §5 `heliacal_ut_vis_lim`.

### `SE_HELFLAG_*` bits relevant to search-strategy selection (swephexp.h:434–449)
```c
#define SE_HELFLAG_LONG_SEARCH        128
#define SE_HELFLAG_HIGH_PRECISION     256
#define SE_HELFLAG_OPTICAL_PARAMS     512
#define SE_HELFLAG_NO_DETAILS        1024
#define SE_HELFLAG_SEARCH_1_PERIOD   (1<<11)  /* 2048 */
#define SE_HELFLAG_VISLIM_DARK       (1<<12)  /* 4096 */
#define SE_HELFLAG_VISLIM_NOMOON     (1<<13)  /* 8192 */
#define SE_HELFLAG_VISLIM_PHOTOPIC   (1<<14)  /* 16384 */
#define SE_HELFLAG_VISLIM_SCOTOPIC   (1<<15)  /* 32768 */
#define SE_HELFLAG_AV                (1<<16)  /* 65536, == AVKIND_VR */
#define SE_HELFLAG_AVKIND_VR          (1<<16)  /* 65536 */
#define SE_HELFLAG_AVKIND_PTO         (1<<17)
#define SE_HELFLAG_AVKIND_MIN7        (1<<18)
#define SE_HELFLAG_AVKIND_MIN9        (1<<19)
#define SE_HELFLAG_AVKIND (SE_HELFLAG_AVKIND_VR|SE_HELFLAG_AVKIND_PTO|SE_HELFLAG_AVKIND_MIN7|SE_HELFLAG_AVKIND_MIN9)
```
**Strategy selection**: `helflag & SE_HELFLAG_AVKIND` nonzero (i.e. **any** of `AVKIND_VR` /
`AVKIND_PTO` / `AVKIND_MIN7` / `AVKIND_MIN9` set) → **arc_vis path** (§6). Zero → **vis_lim
path** (§5, the default/no-flags-set behavior). `SE_HELFLAG_AV` is a legacy alias identical in
value to `AVKIND_VR`.

Within the arc_vis path, `AVKIND_MIN7`/`AVKIND_MIN9` override the self-adjusting `sunsangle`
seed to a fixed `-7°`/`-9°` (§6 step 5 inner loop); `AVKIND_VR` triggers the per-minute
TAV-minimization walkthrough (§6 step 9); `AVKIND_PTO` triggers the symmetric-crossing averaging
(§6 step 10) — these three are independent bits and their effects can combine (VR then PTO
applied in sequence, per the code order) though typical API usage picks one.

**Legacy alias caveat**: `swephexp.h` also defines a `SE_HELIACAL_*` family (455–466) intended
as older-name aliases for the same bits, but `SE_HELIACAL_AVKIND_VR` is defined as `(1<<15)`
(32768) — **one bit position lower than** `SE_HELFLAG_AVKIND_VR` `(1<<16)` — an inconsistency in
the header (the `SE_HELIACAL_VISLIM_SCOTOPIC` bit that exists in the `HELFLAG` family at
`(1<<15)` has no `HELIACAL`-named counterpart, and the `HELIACAL` family's `AVKIND_VR` occupies
that same slot instead). `swehel.c` itself only ever uses the `SE_HELFLAG_*` names. Not
something the Rust port needs to replicate structurally (the Rust API should just expose one
canonical flag set matching `SE_HELFLAG_*`'s semantics), but worth knowing if any caller code
being ported elsewhere uses the `SE_HELIACAL_*` spellings — they are **not** bit-compatible with
`SE_HELFLAG_*` for the AVKIND group.

### `TypeEvent` / `SE_HELIACAL_*` event-type constants (swephexp.h:424–432)
```c
#define SE_HELIACAL_RISING   1   /* = SE_MORNING_FIRST */
#define SE_HELIACAL_SETTING  2   /* = SE_EVENING_LAST */
#define SE_MORNING_FIRST     SE_HELIACAL_RISING    /* 1 */
#define SE_EVENING_LAST      SE_HELIACAL_SETTING   /* 2 */
#define SE_EVENING_FIRST     3
#define SE_MORNING_LAST      4
#define SE_ACRONYCHAL_RISING 5   /* still not implemented [sic — comment stale; it IS handled via arc_vis-path remap, see swe_heliacal_ut below, and vis_lim path via get_asc_obl_with_sun's is_acronychal branch] */
#define SE_ACRONYCHAL_SETTING 6  /* SE_COSMICAL_SETTING alias */
```
The "still not implemented" comments (swephexp.h:430–431) are **stale** — both are handled, just
routed differently depending on `AVKIND`/vis_lim path and object class (see `swe_heliacal_ut`
below).

### `SE_PHOTOPIC_FLAG`/`SE_SCOTOPIC_FLAG`/`SE_MIXEDOPIC_FLAG` (swephexp.h:469–471)
```c
#define SE_PHOTOPIC_FLAG   0
#define SE_SCOTOPIC_FLAG   1
#define SE_MIXEDOPIC_FLAG  2
```
Returned as bits/values in `swe_vis_limit_mag`'s (and by extension `time_optimum_visibility`/
`time_limit_invisible`'s) non-negative return code — see call-contracts section above.

### `swe_heliacal_ut(JDNDaysUTStart, dgeo, datm, dobs, ObjectNameIn, TypeEvent, helflag, *dret, serr_ret) -> OK/ERR/-2` — 3385–3511
Public entry point. Full input-parameter documentation is in the header comment (3345–3384,
transcribed here): `dgeo[3]` = [longitude, latitude, eye-height-m]; `datm[4]` = [pressure hPa,
temperature °C, RH %, VR (meteorological range km, or `1>VR>0` = ktot directly, or `VR=-1` =
compute ktot from other atmospherics)]; `dobs[6]` = [age (default 36), Snellen ratio (default
1), is_binocular, OpticMagn, OpticDia, OpticTrans — the last 4 only apply if
`SE_HELFLAG_OPTICAL_PARAMS` is set]; `TypeEvent` 1–4 (5/6 = acronychal, handled specially, see
below).

1. **Altitude validation**: `dgeo[2]` outside `[SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX]` (same
   `[-500, 25000]` m constants as rise/set, `docs/c-ref-riseset.md` §7) → `ERR`.
2. `swi_set_tid_acc(JDNDaysUTStart, helflag, 0, serr)` — sets the tidal-acceleration model used
   by Δt calculations (STATEFUL — see §8).
3. `MaxCountSynodicPeriod = SE_HELFLAG_LONG_SEARCH ? MAX_COUNT_SYNPER_MAX : MAX_COUNT_SYNPER`
   where `MAX_COUNT_SYNPER = 5` (swehel.c:110, comment "search within 10 synodic periods" — the
   comment is stale/inconsistent with the literal `5`, transcribe the code value `5`, not the
   comment) and `MAX_COUNT_SYNPER_MAX = 1000000` (swehel.c:111, "high, so there is not max
   count" — effectively unbounded). **Dead code** (commented, 3400–3401):
   `SE_HELFLAG_SEARCH_1_PERIOD` would have set `MaxCountSynodicPeriod = 1` directly here — not
   compiled; instead `SEARCH_1_PERIOD`'s effect is applied later (step 8) as a post-hoc
   "reject if the result came from beyond 1.5 periods" check, not a search-space limit.
4. `ObjectName` = lower-cased copy of `ObjectNameIn` (`tolower_string_star`, via
   `strcpy_VBsafe` first — a length-safe copy, since fixed-star name resolution can rewrite the
   buffer to a longer canonical name).
5. `default_heliacal_parameters(datm, dgeo, dobs, helflag)` — fills defaults (see call
   contracts).
6. `swe_set_topo(dgeo[0], dgeo[1], dgeo[2])` — **STATEFUL**, see §8.
7. `Planet = DeterObject(ObjectName)`.
8. **`SE_SUN`** → `ERR` "the sun has no heliacal rising or setting".
9. **`SE_MOON`** branch (3420–3437):
   - `TypeEvent == 1 || 2` → `ERR` (Moon has no morning-first/evening-last).
   - `tjd = tjd0` (= `JDNDaysUTStart`); `MoonEventJDut(tjd, ..., TypeEvent, helflag, dret, serr)`
     (§5/§6 dispatcher).
   - **Retry-forward loop**: `while (retval != -2 && *dret < tjd0) { tjd += 15; ...
     MoonEventJDut(tjd, ...) }` — if a valid (non-`-2`) event was found but it's *before* the
     requested start date (an artifact of the `-30`/`-50`-day backward seed in §5/§6), retry
     15 days later repeatedly until the result is `>= tjd0` or the search itself starts failing
     (`-2`). No iteration cap — relies on eventually converging since each retry moves 15 days
     forward and Moon events recur roughly monthly.
   - Return `retval` (propagating `serr`).
10. **Planets/stars branch** (3438–3511):
    - **Event-type applicability gate** (3441–3454, only when `!(helflag & SE_HELFLAG_AVKIND)`,
      i.e. vis_lim path only): if `Planet == -1 (star) || Planet >= SE_MARS` and
      `TypeEvent == 3 || TypeEvent == 4`: `ERR` — evening-first/morning-last is **not offered**
      via the vis_lim path for outer planets or stars (only Mercury/Venus get those event types
      in the vis_lim path; see §5's branch condition, which sends outer-planet 3/4 to the
      **acronychal** sub-branch instead of erroring — so this validation and that branching
      look at the same condition from two different callers... actually re-examine: this ERR
      fires unconditionally for star/outer-planet + TypeEvent 3/4 regardless of what
      `heliacal_ut_vis_lim` would have done, meaning **this top-level gate is what actually
      prevents that case from ever reaching `heliacal_ut_vis_lim`'s acronychal sub-branch with a
      "3/4"-labeled TypeEvent** — the acronychal sub-branch in §5 is reached only via
      `TypeEvent == SE_ACRONYCHAL_RISING/SETTING` after the remapping in step 11 below, not via
      raw 3/4).
    - **Acronychal TypeEvent remapping for arc_vis path** (3456–3462, only when `helflag &
      SE_HELFLAG_AVKIND`): if `Planet == -1 || Planet >= SE_MARS` and `TypeEvent ==
      SE_ACRONYCHAL_RISING(5)` → `TypeEvent = 3`; `SE_ACRONYCHAL_SETTING(6)` → `TypeEvent = 4`.
      I.e. the arc_vis path (§6, `heliacal_ut_arc_vis`) never sees raw `TypeEvent` 5/6 — it only
      ever receives 1–4, with 3/4 meaning acronychal-for-outer-planets (per §6 step 4's
      `eventtype` remap) or literal evening-first/morning-last for Mercury/Venus.
    - **Acronychal rejection for vis_lim path** (3463–3476, the `else if (1)` — i.e. always
      taken when `AVKIND` is not set): `TypeEvent == SE_ACRONYCHAL_RISING || ==
      SE_ACRONYCHAL_SETTING` → `ERR` "... is not provided for ..." — **raw acronychal
      TypeEvent (5/6) is only ever accepted by the arc_vis path**; the vis_lim path's apparent
      "acronychal branch" in `heliacal_ut_vis_lim` (§5) is reached via `TypeEvent` 3/4 for outer
      planets/stars, which is semantically acronychal but numerically not `5`/`6`. This
      resolves the apparent tension noted in the point above: **`SE_ACRONYCHAL_RISING/SETTING`
      (5/6) themselves are arc_vis-only inputs; the vis_lim path's acronychal handling is
      accessed only through the ordinary 3/4 codes for non-Mercury/Venus bodies.**
    - `dsynperiod = get_synodic_period(Planet)` (§1). `tjdmax = tjd0 + dsynperiod ·
      MaxCountSynodicPeriod`. `tadd = dsynperiod · 0.6`, overridden to `30` (days, not a
      period-fraction) if `Planet == SE_MERCURY`.
    - **Outer synodic-period loop** (3485–3497): `retval = -2` (sentinel, "need another
      period"); `for (tjd=tjd0; tjd<tjdmax && retval==-2; tjd += tadd)`:
      - `heliacal_ut(tjd, dgeo, datm, dobs, ObjectName, TypeEvent, helflag, dret, serr)` (§7
        dispatcher → §5 or §6).
      - **Inner retry-forward loop** (3492–3496): identical pattern to the Moon branch — while
        `retval != -2 && *dret < tjd0`: `tjd += tadd`, retry `heliacal_ut` — pulls a
        too-early result forward by whole `tadd` increments (not the Moon's fixed 15 days;
        Mercury's fixed 30 or others' `0.6·synodic-period`) until it lands at/after `tjd0` or
        the search fails.
      - Outer `for`'s own increment (`tjd += tadd`) then continues *only if* `retval == -2`
        coming out of the inner retry loop (per the `for` condition) — i.e. the outer loop is
        the "try the next synodic period from scratch" fallback when `heliacal_ut` itself
        couldn't find anything in the current period at all (as opposed to finding something
        but it being too early, which the inner loop handles by nudging forward within/across
        periods without abandoning the current period's search state).
11. **Final result classification** (3501–3507):
    - `(helflag & SE_HELFLAG_SEARCH_1_PERIOD) && (retval==-2 || dret[0] > tjd0 + dsynperiod·1.5)`
      → `serr = "no heliacal date found within this synodic period"`, **`retval = -2`** — this
      is `SEARCH_1_PERIOD`'s actual enforcement point (not a search-space limit as the dead code
      at step 3 would have implemented, but a **post-hoc rejection** of any result that took
      more than 1.5 synodic periods to find, treating it as equivalent to "not found").
    - Else `retval == -2` (ran through all `MaxCountSynodicPeriod` periods, nothing found, and
      `SEARCH_1_PERIOD` not overriding): `serr = "no heliacal date found within %d synodic
      periods"`, **`retval = ERR`** (note: **`-2` is promoted to `ERR` here** — the public
      API's final "genuinely not found after exhausting the whole search budget" case is
      reported as a hard error, not the `-2` sentinel used internally throughout the search
      machinery).
12. Return `retval`, with `serr_ret` populated from the local `serr` if non-empty.

### `dret[]` output contract (from the header comment, 3379–3383, and the vis_lim/arc_vis
functions' behavior above)
```
dret[0]: beginning of visibility (Julian day, UT)
dret[1]: optimum visibility (Julian day, UT; 0 if SE_HELFLAG_AV[KIND] any variant set)
dret[2]: end of visibility   (Julian day, UT; 0 if SE_HELFLAG_AV[KIND] any variant set)
```
- **vis_lim path** (§5): all three slots populated when `!SE_HELFLAG_NO_DETAILS` (via
  `get_heliacal_details` for the heliacal branch, or left at `dret[1]=dret[2]=0` — never
  written — for the acronychal branch, since that branch's details-refinement is dead code, see
  §5 step 6). If `SE_HELFLAG_NO_DETAILS` is set, only `dret[0]` is populated (the initial
  `for(i<10) dret[i]=0` zero-fill at the top of `heliacal_ut_vis_lim` means `dret[1]`/`dret[2]`
  are reliably `0` in that case, matching the header comment's "0 if AV" — though the header
  comment's phrasing suggests this is an AVKIND-specific note, in practice it's also true for
  `NO_DETAILS` on the vis_lim path, and always true for the arc_vis path per below).
- **arc_vis path** (§6): only `dret[0]` is ever written by `heliacal_ut_arc_vis`; `dret[1]`/
  `dret[2]` are untouched by that function (not zeroed either — see §6's closing note; whatever
  the caller passed in survives, since `swe_heliacal_ut` itself never zero-fills `dret[]` before
  calling `heliacal_ut`, only `heliacal_ut_vis_lim` does its own internal zero-fill). **A Rust
  port should explicitly zero `dret[1..]` before invoking the arc_vis path**, rather than
  leaving them at caller-supplied garbage, to give deterministic (if not bit-for-bit C-matching)
  output — flag this as an intentional, documented deviation from C's uninitialized-memory
  behavior.
- **Moon path** (`moon_event_vis_lim`, §5 / `moon_event_arc_vis`, §6): `moon_event_vis_lim`
  populates all three slots (start/optimum/end, with the sunset/sunrise clamp and TypeEvent==4
  reorder, §5); `moon_event_arc_vis` populates **only `dret[0]`** (same asymmetry as the
  planet/star arc_vis path).

---

## §8 Porting notes for the stateless Rust port (statics inventory, dead code list, gotchas)

### Global/static state inventory (must become explicit parameters in the Rust port)
1. **`swe_set_topo(lon, lat, alt)`** — called in `swe_heliacal_ut` (step 6 above),
   `rise_set_fast`/`swe_rise_trans_true_hor` (transitively, via `my_rise_trans`/
   `call_swe_rise_trans`, see `docs/c-ref-riseset.md`), `Magnitude`, `swe_vis_limit_mag`. Same
   note as `docs/c-ref-riseset.md` §3.4/§7: thread `dgeo` explicitly instead of mutating a
   shared observer-position cache.
2. **`swi_set_tid_acc(tjd, helflag, 0, serr)`** — sets the tidal-acceleration model used by
   `swe_deltat_ex` calls throughout the search (Δt affects every `tjd_ut → tjd_tt` conversion in
   this file). The Rust `Ephemeris` should take the tidal-acceleration selection as part of its
   config rather than a call-time side effect; verify against `docs/codebase-map.md` whether an
   equivalent already exists for other modules (e.g. eclipse/occultation) before adding a new
   knob.
3. **`SunRA`'s function-local `static TLS double tjdlast/ralast`** (swehel.c:557–558, in the
   part-1/2 region, not re-transcribed here) — a one-entry memoization cache keyed on exact
   `JDNDaysUT` equality. Purely a performance cache (correctness-neutral: recomputing instead of
   hitting the cache gives the identical result, just slower) — the Rust port can simply drop
   this cache (a stateless recompute every call is semantically equivalent; do not port the
   memoization unless profiling of the ported day/minute-stepping loops shows it's needed, since
   `&self` methods have no natural place to stash mutable memo state without introducing
   interior mutability, which the project's stateless-design constraint disallows).
4. **`swi_set_tid_acc`/`swe_set_topo`/`SunRA`'s cache** are the only statics touched by the
   functions in this doc's line range; no other module-level mutable state is read or written
   by `heliacal_ut`/`heliacal_ut_vis_lim`/`heliacal_ut_arc_vis`/`get_heliacal_day`/
   `get_acronychal_day`/`get_asc_obl_with_sun`/`find_conjunct_sun`/`moon_event_*`/
   `time_optimum_visibility`/`time_limit_invisible`/`get_heliacal_details` themselves (they are
   otherwise pure functions of their explicit parameters plus whatever the stateful helpers
   above read).

### Dead code inventory (do NOT port)
- `get_asc_obl_old` (`#if 0`, 2486–2517) — superseded by `get_asc_obl` (2452); only difference
  is `swe_fixstar` vs `call_swe_fixstar` (string-safety wrapper), otherwise byte-identical logic.
- `get_asc_obl_diff_old` (`#if 0`, 2544–2560) — superseded by `get_asc_obl_diff` (2519); lacks
  the `is_acronychal` parameter/branch and the final `>180` wrap.
- `get_asc_obl_with_sun_old` (`#if 0`, 2679–2713) — superseded by `get_asc_obl_with_sun` (2604);
  comment says "works only for fixed stars"; simpler halving-only search, no retro/acronychal
  handling, no loop caps.
- `get_asc_obl_acronychal` (`#if 0`, 2715–2759) — comment says "works only for fixed stars";
  fully superseded by `get_asc_obl_with_sun`'s `is_acronychal` branch.
- `heliacal_ut_arc_vis`'s Moon-position block (`#if 0`, 2323–2333) — would compute `AltM`/`AziM`
  for `HeliacalAngle`; the live code always passes `AltM=-1, AziM=0` instead (a real behavioral
  fact — Moon interference is not modeled in this path — not merely a stripped debug feature).
- `heliacal_ut_arc_vis`'s `swe_heliacal_pheno_ut`-based `JDNarcvisUT` overwrite (`#if 0`,
  2376–2384).
- `heliacal_ut_vis_lim`'s alternate acronychal-details refinement (`else if ((0))`, 3230–3240).
- `get_acronychal_day`'s `azalt_cart`+dot-product `dtret` computation (`#if 0`, 3085–3090) and
  its commented alternate `tret`-seed adjustments (3053–3055, 3059–3062) and outer-loop
  threshold (`#if 0` variant using `0.5` days instead of `0.5/1440` days, 3065–3069 — the `0.5`
  day variant is dead, the minutes variant is live).
- `swe_heliacal_ut`'s commented `SE_HELFLAG_SEARCH_1_PERIOD` direct `MaxCountSynodicPeriod = 1`
  assignment (3400–3401) — superseded by the post-hoc rejection at step 11.
- `HeliacalJDut` (2077–2093) and the dead code around it — legacy VB-interop shim, not part of
  the Rust port's target API (mentioned in §6 for completeness only).
- Miscellaneous debug `printf`/commented `printf` calls throughout (`time_optimum_visibility`'s
  `"hallo -2\n"`, various `/*printf(...)*/` comments in `moon_event_arc_vis` and
  `heliacal_ut_arc_vis`) — never port debug output.

### Other gotchas worth flagging explicitly during implementation review
- **`tcon[]` has no Pluto row** (§1) — a latent C out-of-bounds-read risk if `find_conjunct_sun`
  is ever reached with `ipl == SE_PLUTO`; the Rust port must not reproduce undefined behavior
  here (extend the table with a documented value, or return an explicit error/`Result::Err` for
  Pluto through this code path, and note the deviation).
- **`get_acronychal_day`** (§3) and **`heliacal_ut_arc_vis`'s inner Sun rise/set call** (§6 step
  5) both call `my_rise_trans` without checking for a `-2` (circumpolar/no-rise-or-set) return —
  only `== ERR` is checked. A Rust port should decide explicitly how to handle this (propagate
  `-2` as a distinguished "not found this iteration" the way `get_heliacal_day` does for the
  exact same underlying condition, rather than silently mishandling it) and document the
  decision as a deliberate improvement over the C behavior.
- **`get_acronychal_day`** always writes an informational message to `serr` even on `OK` returns
  (§3) — unlike almost every other function in this file, where `serr` is empty on success. The
  Rust port's warning/diagnostic channel (`CalcResult.flags_used`-style, per project
  `CLAUDE.md`) should surface this, not silently drop it.
- **`heliacal_ut_arc_vis`'s final sanity-bound branch** (§6 step 11) assigns `dret[0] =
  JDNDaysUT` (a stale loop variable) rather than the more obviously-intended `JDNDaysUTinp` —
  likely a latent bug in the original C. Decide and document whether the Rust port replicates
  this exactly (for golden-test fidelity against the reference implementation's actual observed
  output) or fixes it (and if so, under what flag/behavior-versioning scheme the project uses
  for such fixes — check whether `docs/codebase-map.md` or another ref doc has already
  established a convention for "known C bugs we intentionally do not replicate").
- **`heliacal_ut_arc_vis`'s and `heliacal_ut_vis_lim`'s `dret[]` zero-initialization asymmetry**
  (§6/§7 `dret[]` contract) — the arc_vis path never zeros `dret[1..]`. Rust port should
  explicitly zero (or use an `Option`/enum to represent "not computed for this path") rather
  than leave uninitialized memory semantics.
- **Bounds/iteration caps**: most of the fine-grained walkthrough loops in this file
  (`get_asc_obl_with_sun`'s coarse/bisection loops — capped at 5000; `moon_event_arc_vis`'s
  inner minute-loop and `heliacal_ut_arc_vis`'s `AVKIND_VR` walkthrough — **uncapped**;
  `get_heliacal_day`'s minute-refinement `while` — **uncapped**) rely on the underlying physical
  functions being well-behaved (monotonic/single-extremum in the relevant window) to terminate.
  When porting, prefer an explicit conservative iteration cap even where C has none, and flag
  each such addition as a deliberate robustness improvement (not a silent behavior change) per
  the project's `CLAUDE.md` guidance on preserving numerical fidelity — the cap should only ever
  fire as a last-resort safety net, never in normal operation, so it shouldn't affect golden-test
  parity as long as it's set well above any observed real iteration count.
- **Shared parabola-vertex helper**: `x2min` (swehel.c:1791, part-1/2 doc territory) is used by
  both `moon_event_arc_vis` (§6) and `heliacal_ut_arc_vis`'s `AVKIND_VR` walkthrough (§6) — same
  helper, two call sites in this doc's range. Per project `CLAUDE.md` constraints ("reuse
  existing... or extract a shared helper"), port `x2min` once and share it; do not duplicate.
  It is structurally similar to, but distinct in normalization from, `find_maximum`
  (`docs/c-ref-riseset.md` §1) — check whether the Rust port already has a unified
  parabola-vertex utility before adding a second one; if the normalizations genuinely differ
  enough that unifying them would be awkward, that's an acceptable reason to keep them separate
  (per `CLAUDE.md`'s "if the inputs diverge enough... that's fine — use judgment").
