# C Reference: Gauquelin Sector — Rise/Set Branch (`swe_gauquelin_sector`, `imeth` 2–5)

Porting reference for `swe_gauquelin_sector`'s rise/set-based sector methods. The geometric
methods (`imeth` 0/1) are already ported (`Ephemeris::gauquelin_sector_geometric`,
`src/context.rs:616-657`) and are documented in full in `docs/c-ref-houses.md` §10; they are
covered here only briefly, to keep the dispatcher contract complete in one place.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_gauquelin_sector` | swecl.c:6309–6439 | Yes — dispatcher; `imeth` 0/1 already ported, `imeth` 2–5 is this task |
| `swe_rise_trans` | swecl.c:4355–4383 | Already ported: `Ephemeris::rise_trans` (src/context.rs:514–553) |
| `swe_rise_trans_true_hor` | swecl.c:4387–4686 | Already ported: `Ephemeris::rise_trans_true_hor` (src/context.rs:491–506) |
| `swe_deltat_ex` | sweph.c | Already ported: `crate::deltat::calc_deltat` |
| `swe_calc` / `swe_fixstar` | sweph.c / swephlib.c | Already ported: `Ephemeris::calc`, fixstar module |

## 1. Signature & Overall Contract (swecl.c:6309–6324)

```c
int32 swe_gauquelin_sector(
  double t_ut,       /* input time (UT) */
  int32 ipl,         /* planet/moon number; ignored if starname set */
  char *starname,    /* star name, or NULL/empty for a planet/moon */
  int32 iflag,       /* ephemeris flag (SE_SWIEPH/SE_JPLEPH/SE_MOSEPH) + SEFLG_TOPOCTR */
  int32 imeth,       /* 0/1 = Placidus house position (geometric); 2..5 = rise/set-based */
  double *geopos,    /* [geo. longitude, geo. latitude, height above sea] */
  double atpress,    /* mbar; only meaningful for imeth with refraction (3, 5); 0 = auto */
  double attemp,     /* deg C; only meaningful for imeth with refraction (3, 5) */
  double *dgsect,    /* OUT: sector position, 1.0–36.999... */
  char *serr)        /* OUT: error string, may be NULL */
```

Returns `OK` (0) or `ERR` (-1). On error `*dgsect` is left unset by the geometric branch
(the caller's buffer is whatever it was), but the rise/set branch explicitly sets `*dgsect = 0`
on the "rise or set not found" failure path (swecl.c:6434).

**The C doc comment (swecl.c:6289–6308) is stale/incomplete**: it describes only `imeth` 0–3
("use rise and set of body's disc center" / "with refraction") and says nothing about `imeth`
4 and 5. The actual code (swecl.c:6371–6374) supports all four rise/set variants — 4 and 5 are
the "disc edge" (standard visual rise/set) analogs of 2 and 3. A porter reading only the doc
comment would miss two of the four branches entirely.

### `ipl` special-case (swecl.c:6344–6345)

```c
if (ipl == SE_AST_OFFSET + 134340)
  ipl = SE_PLUTO;
```
Asteroid-numbered Pluto (134340) calls are silently redirected to `SE_PLUTO`. Applies to *all*
`imeth` values (this remap happens before the `imeth` branch). Verify this redirect already
exists (or is replicated) wherever `Body`/`ipl` normalization happens upstream of the Rust
dispatcher — if not, it must be added here too.

### `imeth` validation (swecl.c:6337–6341)

```c
if (imeth < 0 || imeth > 5) {
  sprintf(serr, "invalid method: %d", imeth);
  return ERR;
}
```
Simple range check, `0..=5` inclusive. No further validation of `imeth` values within range —
each of the 6 values maps to a defined behavior (see below).

## 2. `imeth` Meaning Table

| `imeth` | Strategy | `SE_BIT_DISC_CENTER` | `SE_BIT_NO_REFRACTION` | Meaning |
|---|---|---|---|---|
| 0 | geometric (`swe_house_pos` 'G') | — | — | Placidus/Gauquelin house position, with ecliptic latitude |
| 1 | geometric | — | — | Same, but ecliptic latitude forced to 0 (projects onto ecliptic) |
| 2 | rise/set | **set** | **set** | Disc-center rise/set, no refraction |
| 3 | rise/set | **set** | not set | Disc-center rise/set, **with** refraction (uses `atpress`/`attemp`) |
| 4 | rise/set | not set | **set** | Disc-edge (standard) rise/set, no refraction |
| 5 | rise/set | not set | not set | Disc-edge (standard) rise/set, **with** refraction (uses `atpress`/`attemp`) |

Derivation (swecl.c:6371–6374):
```c
risemeth = 0;
if (imeth == 2 || imeth == 4) risemeth |= SE_BIT_NO_REFRACTION;
if (imeth == 2 || imeth == 3) risemeth |= SE_BIT_DISC_CENTER;
```
So refraction is *ignored* for `imeth` ∈ {2,4} and *applied* for `imeth` ∈ {3,5}; disc-center
(vs. disc-edge/limb) is used for `imeth` ∈ {2,3} and disc-edge for `imeth` ∈ {4,5}. Only `imeth`
3 and 5 (the refraction-applying variants) actually consult `atpress`/`attemp` — for 2 and 4
those parameters are passed through to `swe_rise_trans` but have no effect since
`SE_BIT_NO_REFRACTION` short-circuits the refraction computation downstream (see
`docs/c-ref-riseset.md` §5.5 / `resolve_atpress`, `azalt.rs:172`).

`SE_BIT_DISC_BOTTOM`, `SE_BIT_GEOCTR_NO_ECL_LAT`, twilight bits, `SE_BIT_FIXED_DISC_SIZE` are
**never** set by this function — none of the 6 `imeth` values touch those flags. Only
`SE_CALC_RISE|risemeth` or `SE_CALC_SET|risemeth` is ever passed to `swe_rise_trans`.

## 3. `atpress`/`attemp` Defaults

`swe_gauquelin_sector` itself does **not** apply any default to `atpress`/`attemp` — it passes
the caller's raw values straight through to every `swe_rise_trans` call (swecl.c:6376, 6395,
6407, 6419, all four calls use the same `atpress, attemp` parameters unmodified). The default
(1013.25 mbar when `atpress == 0`, auto-estimated from `geopos[2]` height via the barometric
formula) is applied **inside** `swe_rise_trans_true_hor`, not here — see `docs/c-ref-riseset.md`
§196–198 and the already-ported `resolve_atpress` (`src/azalt.rs:172`). This matches the C doc
comment's wording ("If imeth=3 and atpress not given (=0), the programm assumes 1013.25 mbar")
— the assumption happens transitively through the rise/set module, not via a check in
`swe_gauquelin_sector` proper. **The Rust port needs no special atpress/attemp handling at the
gauquelin call site** — pass the caller's values through unchanged to `Ephemeris::rise_trans`.

`attemp` has no standalone default in the C source; it defaults implicitly through the same
barometric/refraction formulas (0 °C contributes as `273.15 + attemp` in `calc_dip`/refraction —
see `docs/c-ref-riseset.md` §429).

## 4. `epheflag` Extraction (swecl.c:6333)

```c
int32 epheflag = iflag & SEFLG_EPHMASK;
```
Applied unconditionally at function entry, used for **every** downstream call in *both* branches
(geometric: `swe_deltat_ex`, `swi_epsiln`, `swi_nutation`, `swe_calc`/`swe_fixstar` all receive
the *full* `iflag`, not `epheflag` — only the rise/set branch's `swe_rise_trans` calls use the
masked `epheflag`). This is a subtlety: the geometric branch passes `iflag` (which may include
`SEFLG_TOPOCTR`, `SEFLG_SPEED`, etc.) to `swe_calc`/`swe_fixstar`, while the rise/set branch
passes only the masked ephemeris-source bits to `swe_rise_trans` — topocentric-ness for
rise/set is instead controlled entirely by `geopos` (an observer position always implies
topocentric-like geometry for rise/set) and by the `SE_BIT_GEOCTR_NO_ECL_LAT` bit (never set
here). In Rust terms: `let epheflag = flags & CalcFlags::EPHEMERIS_MASK;` (or equivalent),
passed to `Ephemeris::rise_trans` in place of the full `flags`.

## 5. deltaT / tid_acc Usage

The rise/set branch (`imeth` 2–5) does **not** call `swe_deltat_ex` or `swi_nutation` at all —
those calls exist only in the `imeth` 0/1 geometric branch (swecl.c:6350–6352, already ported).
For `imeth` ≥ 2, all time handling is delegated to `swe_rise_trans`, which internally resolves
its own deltaT/tid_acc via the ephemeris flags it's given (`epheflag`) — see
`docs/c-ref-riseset.md` for that internal resolution. `swe_gauquelin_sector`'s rise/set branch
works entirely in **UT** (`t_ut`, `tret[]`, the interpolation formulas) — there is no `t_et`
variable and no explicit deltaT arithmetic anywhere in this branch.

## 6. Fixed-Star Support (`starname`)

```c
AS_BOOL do_fixstar = (starname != NULL && *starname != '\0');
```
computed once at function entry (swecl.c:6334), used by **both** branches:
- Geometric (`imeth` 0/1): dispatches to `swe_fixstar(starname, t_et, iflag, x0, serr)` vs.
  `swe_calc(t_et, ipl, iflag, x0, serr)` (swecl.c:6356–6362).
- Rise/set (`imeth` 2–5): `starname` (not `do_fixstar`) is passed as-is to every
  `swe_rise_trans` call (swecl.c:6376, 6395, 6407, 6419) — `swe_rise_trans` itself does the
  `starname != NULL && *starname` check internally to decide fixed-star vs. planet dispatch.
  `ipl` is passed alongside `starname` in every call regardless of `do_fixstar`; `swe_rise_trans`
  ignores `ipl` when `starname` is non-empty (same convention as `swe_gauquelin_sector` itself).

In the Rust port this maps directly to `Ephemeris::rise_trans`'s existing `starname:
Option<&str>` parameter — no new logic needed, just thread the same `Option<&str>` through from
the dispatcher's own `starname` argument.

## 7. Rise/Set Search Sequence (swecl.c:6375–6432)

State: `tret[3]` (only `tret[0]` = rise, `tret[1]` = set are used — `tret[2]` unused here),
`rise_found`/`set_found` bools (init `TRUE`), `above_horizon` bool (init `FALSE`).

### 7.1 Step 1 — find the next rising (swecl.c:6376–6393)

```c
retval = swe_rise_trans(t_ut, ipl, starname, epheflag, SE_CALC_RISE|risemeth,
                        geopos, atpress, attemp, &tret[0], serr);
```
- `retval == ERR` (-1): propagate `ERR` immediately (hard failure — bad ephemeris file, invalid
  input, etc.). `*dgsect` is **not** set to 0 on this path (differs from the "not found" path in
  §7.4).
- `retval == -2` (circumpolar — body never rises in the search window): `rise_found = FALSE`.
  This is *not* treated as fatal here; the C comment (swecl.c:6379–6391) explains this is
  deliberately tolerant, with a note that no Gauquelin-sector algorithm exists yet for
  circumpolar bodies (unlike the Placidus/Otto-Ludwig meridian-transit fallback used elsewhere
  for circumpolar house cusps) — search continues to step 2 regardless.
- Otherwise (`retval == 0`, found): `tret[0]` holds the next rise time, `rise_found` stays `TRUE`.

### 7.2 Step 2 — find the next setting (swecl.c:6394–6400)

```c
retval = swe_rise_trans(t_ut, ipl, starname, epheflag, SE_CALC_SET|risemeth,
                        geopos, atpress, attemp, &tret[1], serr);
```
Same `ERR`/`-2` handling as step 1, but sets `set_found = FALSE` on `-2` (does **not** touch
`rise_found`).

### 7.3 Step 3 — bracket determination + one re-search (swecl.c:6401–6425)

This determines whether `t_ut` currently sits **above** or **below** the horizon, and re-derives
whichever of `tret[0]`/`tret[1]` is *not* the immediately-bracketing event, by searching
backward from a point just before the found event.

```c
if (tret[0] < tret[1] && rise_found == TRUE) {
  /* next rise comes before next set → currently BELOW horizon (just set, waiting to rise) */
  above_horizon = FALSE;
  t = t_ut - 1.2;
  if (set_found) t = tret[1] - 1.2;
  set_found = TRUE;
  retval = swe_rise_trans(t, ipl, starname, epheflag, SE_CALC_SET|risemeth,
                          geopos, atpress, attemp, &tret[1], serr);
  if (retval == ERR) return ERR;
  else if (retval == -2) set_found = FALSE;
} else if (tret[0] >= tret[1] && set_found == TRUE) {
  /* next set comes before/at next rise → currently ABOVE horizon */
  above_horizon = TRUE;
  t = t_ut - 1.2;
  if (rise_found) t = tret[0] - 1.2;
  rise_found = TRUE;
  retval = swe_rise_trans(t, ipl, starname, epheflag, SE_CALC_RISE|risemeth,
                          geopos, atpress, attemp, &tret[0], serr);
  if (retval == ERR) return ERR;
  else if (retval == -2) rise_found = FALSE;
}
```

Key details, easy to get wrong in a port:
- The branch condition on the first `if` is `tret[0] < tret[1] && rise_found == TRUE` — note
  `rise_found` is unconditionally reset to `TRUE` *inside* that branch right before the
  re-search (`set_found = TRUE;` — actually it's `set_found`, not `rise_found`, that gets reset
  inside the first branch; the *second* branch resets `rise_found = TRUE`). This means: even if
  the original step-1/2 search failed to find one of the two events (`-2`), the re-search
  re-attempts to find it and optimistically sets the found-flag back to `TRUE` before the retry
  — it is the retry's own `retval` that determines the final found-flag, not this preset.
  A literal Rust port must replicate `set_found = TRUE` (resp. `rise_found = TRUE`) as a
  *pre*-assignment before the retry call, then let the retry's result overwrite it via `-2`
  handling.
- The magic offset **`1.2` (days)** is a fixed backward nudge — not derived from orbital period
  or anything else. It's used only to pick a safe starting epoch for the backward search (so
  the backward `SE_CALC_SET`/`SE_CALC_RISE` search doesn't immediately re-find the same event
  it's trying to move behind). If the *other* event was already found in step 1/2, the search
  point uses that found time minus 1.2 days instead of `t_ut - 1.2`, i.e. it prefers to seed the
  backward search from as close to `t_ut` as safely possible.
- If neither `if` nor `else if` condition is true (i.e. `tret[0] < tret[1]` but
  `rise_found == FALSE`, or `tret[0] >= tret[1]` but `set_found == FALSE`) — no re-search
  happens at all; execution falls through to step 4 with whatever `rise_found`/`set_found`
  state resulted from steps 1–2. This is an easy branch to miss: the re-search is *conditional*
  on the very flag that determines which comparison arm fires, so a body that fails to find
  BOTH rise and set in steps 1-2 never gets a re-search attempt (unsurprising — both are already
  circumpolar-flagged) and falls straight to the "not found" error in step 4. **In the surprising
  case** — `rise_found == FALSE` but the comparison happens to be `tret[0] < tret[1]` — the
  `&& rise_found == TRUE` guard skips the whole re-search block silently (no error raised at this
  point; it's caught by step 4's final check instead).

### 7.4 Step 4 — sector interpolation or final failure (swecl.c:6426–6438)

```c
if (rise_found && set_found) {
  if (above_horizon) {
    *dgsect = (t_ut - tret[0]) / (tret[1] - tret[0]) * 18 + 1;
  } else {
    *dgsect = (t_ut - tret[1]) / (tret[0] - tret[1]) * 18 + 19;
  }
  return OK;
} else {
  *dgsect = 0;
  sprintf(serr, "rise or set not found for planet %d", ipl);
  return ERR;
}
```
- **Above horizon** (sectors 1–18): linear interpolation of `t_ut`'s fractional position between
  the bracketing rise (`tret[0]`) and set (`tret[1]`), scaled to `[0, 18)` and offset by `+1` →
  range `[1, 19)`.
- **Below horizon** (sectors 19–36): linear interpolation between set (`tret[1]`) and rise
  (`tret[0]`), same `*18` scale, offset by `+19` → range `[19, 37)`.
- Note the error message always says `"...for planet %d"` with `ipl` even when `do_fixstar` is
  true and `ipl` was never meaningfully used (cosmetic C wart — `ipl` in that case is whatever
  the caller passed, possibly a garbage/unused value when a star name drove the actual lookup).
  A faithful port can decide whether to preserve this exact wording or produce a
  star-name-aware message via a distinct `Error` variant — either is a legitimate call since the
  string itself is never asserted on in the golden tests (it's a `serr`, not the return-code
  contract).

## 8. Circumpolar / `-2` Contract Summary

`swe_rise_trans` returns:
- `OK` (0): event found, `tret[i]` populated.
- `ERR` (-1): hard failure (propagate as `ERR`/`Err` immediately, no `*dgsect` assignment in the
  geometric-parity sense — Rust: propagate via `?`).
- `-2`: circumpolar in the requested direction — body doesn't rise/set within the search window.
  Treated as soft: the corresponding `found` flag is cleared but the function keeps going. Only
  if **both** rise and set ultimately remain unfound (after the one re-search attempt) does the
  function fail, with `*dgsect = 0` and a formatted `serr`.

In the already-ported Rust `Ephemeris::rise_trans`, C's `-2` is signaled as
`Err(Error::CircumpolarBody)` (src/error.rs:19, 46) rather than a sentinel return value — the
port of `swe_gauquelin_sector`'s rise/set branch must pattern-match on
`Err(Error::CircumpolarBody)` specifically (treat it as the soft "not found" case, matching
`retval == -2`) versus any other `Err(...)` (propagate as a hard error), rather than checking an
integer return code.

## 9. Geometric Branch (`imeth` 0/1) — Summary (swecl.c:6349–6367)

Already ported as `Ephemeris::gauquelin_sector_geometric` (src/context.rs:616–657). Included here
only for dispatcher completeness:

```c
t_et = t_ut + swe_deltat_ex(t_ut, iflag, serr);
eps = swi_epsiln(t_et, iflag) * RADTODEG;
swi_nutation(t_et, iflag, nutlo); nutlo[0..1] *= RADTODEG;
armc = degnorm(swe_sidtime0(t_ut, eps+nutlo[1], nutlo[0]) * 15 + geopos[0]);
x0 = do_fixstar ? swe_fixstar(starname, t_et, iflag, x0, serr)
                : swe_calc(t_et, ipl, iflag, x0, serr);
if (imeth == 1) x0[1] = 0;   /* zero ecliptic latitude */
*dgsect = swe_house_pos(armc, geopos[1], eps + nutlo[1], 'G', x0, NULL);
return OK;
```
Uses the **full** `iflag` (not `epheflag`) throughout — see §4 above for why this differs from
the rise/set branch. See `docs/c-ref-houses.md` §10 for the `swe_house_pos('G', ...)` formula
this delegates to.

## Porting Notes

- `swe_rise_trans` → already ported as `Ephemeris::rise_trans` (`src/context.rs:514–553`).
  C's `-2` circumpolar sentinel → `Err(Error::CircumpolarBody)` (`src/error.rs:19`). C's `ERR`
  (-1) → any other `Err(Error::...)`, propagate with `?`.
- `swe_calc` → already ported as `Ephemeris::calc`.
- The existing geometric branch is `Ephemeris::gauquelin_sector_geometric`
  (`src/context.rs:616–657`); its `imeth ∈ {2,3,4,5}` arm currently returns
  `Err(Error::CError("gauquelin rise/set methods need rise_trans (not ported)".into()))`
  (src/context.rs:628–632) — replace that early-return with the algorithm in §2/§7 above, reusing
  the same function (or splitting into a private helper) rather than duplicating the geometric
  arm. Signature will need `geopos: [f64; 3]` (height is needed for `atpress` auto-resolution
  inside `rise_trans`, unlike the geometric branch which only needs `geolon`/`geolat`) and
  `atpress`/`attemp: f64` parameters threaded straight through to `rise_trans` — see §3 (no
  gauquelin-level defaulting needed, `rise_trans`'s internal `resolve_atpress` handles it).
- `risemeth` derivation (§2) maps directly to `RiseSetFlags::NO_REFRACTION` /
  `RiseSetFlags::DISC_CENTER` (already-defined bits in `crate::flags::RiseSetFlags`, used
  throughout `src/riseset.rs`) OR'd with `RiseSetFlags::RISE`/`RiseSetFlags::SET` per call.
- `epheflag` extraction (§4): use `flags & CalcFlags::EPHEMERIS_MASK` (or whatever the existing
  mask constant/method is named in `src/flags.rs` — check for an existing `.ephemeris_mask()` or
  similar helper before reintroducing the mask bits) for the `epheflag` passed to `rise_trans`,
  but keep the *full* `flags` for the geometric branch's `calc`/fixstar calls (already correct in
  the existing ported code, per its own doc comment at src/context.rs:612–614 about *not*
  matching `swe_houses_ex2`'s deltaT override behavior — same spirit of "read what this call site
  actually does, not what a sibling call site does" applies to the epheflag masking here).
- The `1.2`-day backward-search offset (§7.3) and the "`found` flag pre-set to `TRUE` before
  retry, then possibly re-cleared by `-2`" ordering (§7.3) are the two easiest details to get
  subtly wrong in a literal translation — write them as literal `1.2`-day arithmetic and
  pre-assignment-then-overwrite, don't "clean up" the control flow, since golden tests will
  compare against the exact C sector value at specific epochs including near rise/set boundaries.
- Golden test note: exercise all four rise/set `imeth` values (2–5) plus at least one high
  (or near-polar) latitude case to hit the `-2`/re-search path, since normal mid-latitude bodies
  will almost never trigger the circumpolar branch of §7.3/§8.
