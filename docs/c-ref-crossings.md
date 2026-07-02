# C Reference: Zodiacal/Node/Heliocentric Crossings — sweph.c

Porting reference for the `swe_*cross*` family. Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_solcross` | sweph.c:8321–8343 | Yes — ET variant |
| `swe_solcross_ut` | sweph.c:8355–8377 | Yes — UT wrapper |
| `swe_mooncross` | sweph.c:8389–8411 | Yes — ET variant |
| `swe_mooncross_ut` | sweph.c:8425–8447 | Yes — UT wrapper |
| `swe_mooncross_node` | sweph.c:8456–8486 | Yes — ET variant |
| `swe_mooncross_node_ut` | sweph.c:8493–8523 | Yes — UT wrapper |
| `swe_helio_cross` | sweph.c:8533–8569 | Yes — ET variant |
| `swe_helio_cross_ut` | sweph.c:8579–8615 | Yes — UT wrapper |

Shared helpers used throughout (already ported, see `src/math.rs`):
- `swe_degnorm` (swephlib.c:106–113): `y = fmod(x, 360.0); if fabs(y) < 1e-13 { y = 0.0 }; if y < 0.0 { y += 360.0 }`
- `swe_difdeg2n` (swephlib.c:3819–3824): `dif = swe_degnorm(p1 - p2); if dif >= 180.0 { dif - 360.0 } else { dif }` — signed shortest-arc difference in `(-180, 180]`... actually range is `[-180, 180)` since `dif` from degnorm is `[0,360)` and the `>=180` branch maps to `[-180,0)`, so overall range is `[-180, 180)`.

## Constant

```c
#define CROSS_PRECISION (1 / 3600000.0)   // sweph.c:8308 — one milliarcsecond (in degrees)
```

Used as the convergence threshold in every iteration loop below (`fabs(dist) < CROSS_PRECISION`).

## 1. `swe_solcross` / `swe_solcross_ut` (sweph.c:8321–8377)

Computes the next JD (`jd > jd_et`/`jd > jd_ut`) at which the **Sun's** ecliptic longitude equals `x2cross`.

### Flag handling
- `flag |= SEFLG_SPEED` — forced unconditionally so `x[3]` (longitude speed, deg/day) is available (sweph.c:8329, 8363).
- `ipl = SE_SUN` — hardcoded body.
- Doc comment (sweph.c:8312–8315) states flag covers `SEFLG_HELCTR` (0=geocentric Sun / 1=heliocentric Earth), `SEFLG_TRUEPOS`, `SEFLG_NONUT` — these are simply whatever bits the caller passed through to `swe_calc`/`swe_calc_ut` unmodified. No sidereal mention for solcross (contrast with mooncross_ut, which explicitly documents `SEFLG_SIDEREAL` support) — but nothing in the code strips `SEFLG_SIDEREAL`; it would pass through to `swe_calc` like any other flag.

### Algorithm
1. **Initial calc** at `jd_et`/`jd_ut`: `swe_calc(jd_et, SE_SUN, flag, x, serr)`. If `< 0` (error), **return `jd_et - 1`** (ET) / `jd_ut - 1` (UT) — i.e. one full day before input, guaranteed `< jd_et`.
2. **Mean-speed estimate**: `xlp = 360.0 / 365.24;  /* mean solar speed */` (deg/day, sweph.c:8332/8366).
3. `dist = swe_degnorm(x2cross - x[0])` — unsigned forward distance in `[0, 360)` from current Sun longitude to target.
4. `jd = jd_et + dist / xlp` — linear estimate using mean speed.
5. **Iteration loop** (`for(;;)`, sweph.c:8335–8341 / 8369–8375):
   ```c
   for(;;) {
     if (swe_calc(jd, ipl, flag, x, serr) < 0) return jd_et - 1;
     dist = swe_difdeg2n(x2cross, x[0]);   // signed shortest-arc difference, note argument order (target, current)
     jd += dist / x[3];                     // Newton-like correction using ACTUAL (not mean) speed x[3]
     if (fabs(dist) < CROSS_PRECISION) break;
   }
   return jd;
   ```
   No max-iteration guard — this is an unbounded loop; convergence relies on Newton's method behaving well near a monotonic crossing. Note the break test uses the `dist` computed **before** the `jd` update of that same iteration (i.e. it checks whether the correction just applied was already tiny, not whether the corrected `jd` is exact — one extra `swe_calc` at the converged `jd` is never performed after the last update).

### Error contract
Every `swe_calc`/`swe_calc_ut` failure inside solcross returns `jd_et - 1` (or `jd_ut - 1`), i.e. exactly one day before the *original input* date — **not** relative to the current iterate `jd`. Caller detects error by `return_value < jd_et` (documented at sweph.c:8319: "Errors are indicated by returning a jd < jd_et!").

## 2. `swe_mooncross` / `swe_mooncross_ut` (sweph.c:8389–8447)

Identical structure to solcross, with:
- `ipl = SE_MOON`.
- Mean speed constant: `xlp = 360.0 / 27.32;  /* mean lunar speed */` (sweph.c:8400/8436) — sidereal month approximation, distinct from solcross's tropical-year constant.
- Same `SEFLG_SPEED` forcing, same error convention (`return jd_et - 1` / `jd_ut - 1`), same convergence test.
- `swe_mooncross_ut`'s doc comment (sweph.c:8415–8423) explicitly lists `SEFLG_SIDEREAL` as a recognized bit: "If sidereal is chosen, default mode is Fagan/Bradley. For different ayanamshas, `swe_set_sid_mode()` must be called first." This is purely a documentation note — the sidereal mode is *global state* set via `swe_set_sid_mode` beforehand, not passed through this function; the function itself does nothing special with the bit other than forward it to `swe_calc_ut`.

## 3. `swe_mooncross_node` / `swe_mooncross_node_ut` (sweph.c:8456–8523)

Computes the next JD at which the Moon's **ecliptic latitude** crosses zero (i.e. the Moon crosses its own node — ascending or descending, whichever comes first).

### Flag handling
- `flag |= SEFLG_SPEED` forced (sweph.c:8461/8498) — needed for `x[4]` (latitude speed, deg/day) used in the refinement step.
- `ipl = SE_MOON` hardcoded.
- No explicit bit-masking beyond that; all other flag bits pass through to `swe_calc`/`swe_calc_ut` unchanged.

### Algorithm — two-phase: coarse day-stepping, then Newton refinement

**Phase 1 — bracket the sign change** (sweph.c:8462–8472):
```c
if (swe_calc(jd_et, ipl, flag, x, serr) < 0) return jd_et - 1;
xlat = x[1];              // latitude at the start date
jd = jd_et + 1;
for(;;) {                 // get to sign change
  if (swe_calc(jd, ipl, flag, x, serr) < 0) return jd_et - 1;
  if ((x[1] >= 0 && xlat < 0) || (x[1] < 0 && xlat > 0))
    break;
  jd += 1;                // step size: exactly 1 day
}
```
- Step size is a fixed **1.0 day** — no adaptive stepping.
- Sign-change test: `xlat` (the *previous* sample, initialized to the latitude at `jd_et`) compared against the *current* `x[1]`. Note the asymmetric comparison: `x[1] >= 0 && xlat < 0` (catches ascending crossing, using `>=` so an exact `x[1] == 0` counts as "crossed") vs. `x[1] < 0 && xlat > 0` (descending crossing, strict `<` and `>`, so `xlat == 0` would not trigger — asymmetry is presumably harmless in practice since exact zero latitude at an integer-day sample is measure-zero). **`xlat` is never updated inside this loop** — it stays pinned to the latitude at the original `jd_et`/`jd_ut` for every iteration of the day-stepping loop. This means the loop is really testing "has `x[1]` changed sign relative to the very first sample," not relative to the previous day — for a monotonically increasing/decreasing latitude this is equivalent to catching the first sign flip, but it is a datum a porter must reproduce (not `xlat = x[1]` at the end of the loop body).

**Phase 2 — Newton refinement to zero latitude** (sweph.c:8473–8484):
```c
dist = x[1];               // latitude at the bracketing jd found above
for(;;) {
  jd -= dist / x[4];        // x[4] = latitude speed (deg/day); note MINUS
  if (swe_calc(jd, ipl, flag, x, serr) < 0) return jd_et - 1;
  dist = x[1];
  if (fabs(dist) < CROSS_PRECISION) {
    *xlon = x[0];
    *xla = x[1];
    break;
  }
}
return jd;
```
- Correction direction is `jd -= dist / x[4]` (minus), unlike the longitude-crossing functions which use `jd += dist / x[3]` (plus). This is because here `dist` is the latitude itself (target is 0), and the sign convention of `x[4]` (dlat/dt) relative to `x[1]` (lat) requires subtraction to converge — i.e. it's Newton's method `jd_new = jd - f(jd)/f'(jd)` in its literal form, whereas the longitude functions frame it as `jd_new = jd + (target - current)/speed`, which is algebraically the same Newton step but written with the sign folded into `swe_difdeg2n`'s `(target - current)` order.
- Output params `*xlon`/`*xla` are only written once, inside the converged branch, using the **already-updated** `x[0]`/`x[1]` from the `swe_calc` call at the top of that same iteration (i.e. from the JD that satisfies the precision test, not from a subsequent recompute).
- No max-iteration guard on either phase.

### Error contract
Same as solcross/mooncross: any `swe_calc`/`swe_calc_ut` failure at any point (bracketing phase or refinement phase) returns `jd_et - 1` / `jd_ut - 1` (fixed offset from the *original* input, not the current iterate).

## 4. `swe_helio_cross` / `swe_helio_cross_ut` (sweph.c:8533–8615)

Computes a planet's **heliocentric** longitude crossing of `x2cross`, in either direction (`dir >= 0` → next crossing after `jd_et`; `dir < 0` → previous crossing before `jd_et`). Signature differs from the others: returns `int32` status (`OK`/`ERR`), with the actual JD written through the `jd_cross` out-pointer.

### Flag handling
```c
int flag = iflag | SEFLG_SPEED | SEFLG_HELCTR;   // sweph.c:8537 / 8583
```
- Both `SEFLG_SPEED` and `SEFLG_HELCTR` are **forced on** unconditionally (heliocentric is not optional here — the function is explicitly the heliocentric-crossing API).
- No special handling of `SEFLG_SIDEREAL` — passes through in `iflag` like any other bit.

### Rejected `ipl` values (sweph.c:8538–8547)
```c
if (ipl == SE_SUN
  || ipl == SE_MOON
  || (ipl >= SE_MEAN_NODE && ipl <= SE_OSCU_APOG)
  || (ipl >= SE_INTP_APOG && ipl < SE_NPLANETS)
) {
  char snam[AS_MAXCH];
  swe_get_planet_name(ipl, snam);
  if (serr != NULL) sprintf(serr, "swe_helio_cross: not possible for object %d = %s", ipl, snam);
  return ERR;
}
```
- Rejects: the Sun (no heliocentric "crossing" of itself), the Moon (heliocentric Moon is meaningless/degenerate — Moon's heliocentric motion mirrors Earth's), and the lunar-node/apogee family: everything in the inclusive range `[SE_MEAN_NODE, SE_OSCU_APOG]` (mean/true/osculating node & apogee group) and `[SE_INTP_APOG, SE_NPLANETS)` (interpolated apogee/perigee group, up to but excluding the fictitious-body/asteroid boundary `SE_NPLANETS`).
- Error message is built via `swe_get_planet_name` even though `serr` may be `NULL` (name lookup happens unconditionally; only the `sprintf` is guarded by `serr != NULL`).

### Algorithm
1. Initial calc at `jd_et`/`jd_ut` with the forced `flag`; error → immediately `return ERR` (no `jd_cross` write).
2. **Speed selection** (sweph.c:8550–8552):
   ```c
   xlp = x[3];                  // actual instantaneous longitude speed, NOT a mean-speed constant
   if (ipl == SE_CHIRON)
     xlp = 0.01971;             // use mean speed — special-cased because Chiron's heliocentric
                                 // speed can pass through zero/retrograde near perihelion in ways
                                 // that break a naive initial estimate
   ```
   Note: unlike solcross/mooncross (which always use a hardcoded mean-speed constant for the *initial estimate*), helio_cross uses the just-computed **actual** speed `x[3]` for every body except Chiron, which falls back to a mean constant `0.01971` deg/day.
3. `dist = swe_degnorm(x2cross - x[0])` — same unsigned forward-distance convention as the others.
4. **Direction handling** (sweph.c:8554–8559):
   ```c
   if (dir >= 0) {
     jd = jd_et + dist / xlp;
   } else {
     dist = 360.0 - dist;
     jd = jd_et - dist / xlp;
   }
   ```
   For backward search (`dir < 0`), `dist` is complemented to `360 - dist` (the *backward* forward-distance, i.e. distance if you go the other way around the circle) and then the estimate steps **back** in time (`jd_et - ...`). This assumes `xlp` retains the same sign convention (prograde motion) in both directions; there's no separate handling for a retrograde body at the initial epoch beyond the Chiron mean-speed fallback — a genuinely retrograde `x[3]` at the initial sample would produce a nonsensical initial estimate (this is documented as a known-rough tool: "This should only be used for rough house entry or exit times," sweph.c:8531/8577).
5. **Iteration loop** (sweph.c:8560–8566 / 8606–8612) — identical shape to solcross's loop, using actual speed `x[3]` (not `xlp`) for every correction step, regardless of direction or Chiron special-casing (the Chiron mean-speed constant is *only* used for the initial estimate in step 4, never inside the refinement loop):
   ```c
   for(;;) {
     if (swe_calc(jd, ipl, flag, x, serr) < 0) return ERR;
     dist = swe_difdeg2n(x2cross, x[0]);
     jd += dist / x[3];
     if (fabs(dist) < CROSS_PRECISION) break;
   }
   *jd_cross = jd;
   return OK;
   ```
   No max-iteration guard. No direction-awareness inside the loop itself — once the initial estimate lands on the correct side, ordinary Newton convergence (via `+= dist/x[3]`) is direction-agnostic; there is no explicit check that the converged `jd` is actually on the requested side of `jd_et` (i.e. `dir` only steers the *initial guess*, not the convergence target — a pathological case near a retrograde loop could in principle converge to a crossing on the wrong side, which the C code does not guard against).

### Error contract (differs from the double-returning functions!)
`swe_helio_cross`/`_ut` return an **`int32` status code** (`OK` = 0 conventionally, or `ERR`), *not* a `double` JD with a magnitude-based error sentinel. `jd_cross` (the out-pointer) is **only written on the `OK` path** — on any `ERR` return (rejected `ipl`, or a `swe_calc` failure inside either the initial calc or the loop), `*jd_cross` is left **untouched** (uninitialized from the caller's perspective unless they pre-zeroed it). This is a fundamentally different error contract from the four double-returning functions above.

## Known test vectors

`/home/josh/nhs/soft/astrology/swisseph/setest/suite_10_solcross.c` (TESTSUITE 10) exercises all eight functions (TESTCASE 1–8), but the actual numeric inputs/expected outputs are pulled at runtime via `GET_D`/`GET_I` macros from an external test-data file that is **not present in this repository checkout** (only the `.c` driver exists; no accompanying data file was found under `setest/`). No concrete `(x2cross, jd, iflag, expected_jx)` triples could be extracted.

What the test structure does confirm (useful as invariants for golden tests even without concrete numbers):
1. **TESTCASE 1/2** (`swe_solcross`/`_ut`): after computing `jx`, the test recomputes `swe_calc(jx, SE_SUN, iflag, xx, serr)` and asserts `CHECK_EQUALS_D(xcross, xx[0])` — i.e. the Sun's longitude *at the returned crossing JD* must equal the requested `xcross` value exactly (within the harness's equality tolerance), confirming `CROSS_PRECISION` (1 mas) is tight enough to satisfy the test's equality check.
2. **TESTCASE 3/4** (`swe_mooncross`/`_ut`): identical pattern against `SE_MOON`.
3. **TESTCASE 5/6** (`swe_mooncross_node`/`_ut`): asserts `CHECK_EQUALS_D(xx[1], 0)` — i.e. recomputing the Moon's position at the returned `jx` must give **exactly** latitude 0 (again bounded by the harness tolerance, not literal bit-exact 0.0), confirming the node-crossing refinement converges to true zero latitude, not just "close to zero."
4. **TESTCASE 7/8** (`swe_helio_cross`/`_ut`): only checks `rc` (return status) and `serr`, plus `CHECK_D(jx)` (records the value for regression comparison) — no independent recomputation/re-verification step is performed for the heliocentric case, and `ipl`/`dir` are pulled from the (missing) data file as free parameters, so we don't know which bodies/directions are exercised.
5. `iephe` (ephemeris selector, presumably `SEFLG_JPLEPH`/`SEFLG_SWIEPH`/`SEFLG_MOSEPH`) is set up as `iflag` directly in every test case's `SETUP`, and `swe_set_jpl_file("de431.eph")` is configured globally for the suite — implying at least some cases run against the JPL backend.

Because no data file is available, golden tests for `swisseph-rs`'s `crossings.rs` module should be generated directly against the C library (per `docs/golden-testing.md`) rather than transcribed from this suite.

## Porting notes

- **Global-state reads**: none of these eight functions read `swed` fields directly — the *only* interaction with global state is indirect, through `swe_calc`/`swe_calc_ut`'s own internals (ephemeris file cache, sidereal-mode global set via `swe_set_sid_mode`, nutation cache, etc.), which is out of scope for this doc and already handled by the ported `Ephemeris::calc`/`calc_ut`. `CROSS_PRECISION` is a compile-time constant, not state. No mutable module-level statics are touched by `solcross`/`mooncross`/`mooncross_node`/`helio_cross` themselves.
- **`swe_calc` → `Ephemeris::calc`**: every `swe_calc(jd, ipl, flag, x, serr)` call in this file becomes `self.calc(jd_tt, body, flags)` (src/context.rs:193) returning `Result<CalcResult, Error>`; `swe_calc_ut` → `self.calc_ut(jd_ut, body, flags)` (src/context.rs:236). The array indices used here map onto `CalcResult` fields as: `x[0]` → longitude, `x[1]` → latitude, `x[3]` → longitude speed, `x[4]` → latitude speed (assuming the existing `CalcResult` field layout/order documented in `docs/codebase-map.md`/`src/context.rs` — confirm exact field names when implementing).
- **Forced flags**: `flags::CalcFlags` already defines `SPEED` (bit 256) and `HELCTR` (bit 8) per `src/flags.rs:9,13`. Port `flag |= SEFLG_SPEED` as `flags | CalcFlags::SPEED` (and additionally `| CalcFlags::HELCTR` for the two `helio_cross` functions) — do this once at the top of each function, mirroring the C's single assignment, not per-iteration.
- **Error signaling — the central translation problem**: The four double-returning C functions (`solcross`, `mooncross`, `mooncross_node`, and their `_ut` twins) all signal failure by **returning a JD value less than the original input JD** — specifically, exactly `jd_et - 1` (or `jd_ut - 1`), one day earlier, regardless of *where* in the iteration the underlying `swe_calc` failed. The contract is documented in the C comments verbatim as: *"Errors are indicated by returning a jd < jd_et!"* (sweph.c:8319, 8353, 8387, 8421, 8454, 8491). This is a magnitude-based sentinel with no information about *why* it failed beyond whatever `serr` was populated with by the failing inner `swe_calc` call.
  - For the Rust port: **do not reproduce the `jd - 1` sentinel**. Each function should return `Result<f64, Error>` and propagate the `Error` from the first failing internal `self.calc(...)?` call directly (via `?`), preserving whatever structured `Error` variant that call produced (most likely `Error::BeyondEphemerisLimits` or `Error::EphemerisNotAvailable`, given these functions are typically driven out of the ephemeris's valid date range). This is strictly more informative than the C sentinel and matches the project's `Result<T, Error>` idiom (per `CLAUDE.md` — "Error handling via `Error` enum. No string buffer passing patterns from C.").
  - `swe_helio_cross`/`_ut` already return a structured status (`int32` `OK`/`ERR`) with the real value out-of-band via `jd_cross`. This maps even more directly onto `Result<f64, Error>`: `OK` + `*jd_cross` → `Ok(jd_cross)`; the rejected-`ipl` case → a dedicated `Error` variant (e.g. `Error::InvalidBody(ipl)` or a new variant if the existing enum doesn't cover "valid body but not supported for this operation" — check `src/error.rs` for the closest existing fit, e.g. compare to how other functions reject unsupported body/flag combinations such as `Error::UnsupportedFlags`); any inner `swe_calc` failure → propagate that `Error` via `?`.
- **`mooncross_node`'s frozen `xlat` bug-for-bug detail**: when porting the day-stepping bracket search, do **not** update the comparison variable on each loop iteration — it must stay fixed at the latitude sampled at the *original* input JD for the entire bracketing loop, exactly as the C does (sweph.c:8464, 8501: `xlat = x[1]` is set once, before the loop, and never reassigned inside it). Getting this wrong (e.g. sliding-window sign comparison) would only produce a different result in the edge case of an oscillating/non-monotonic latitude within the first few days, but must still be replicated for bit-for-bit fidelity.
- **Sign convention divergence between longitude-crossing and node-crossing refinement**: longitude crossers use `jd += dist / speed` where `dist = swe_difdeg2n(target, current)`; the node crosser uses `jd -= dist / x[4]` where `dist = x[1]` (current latitude, target implicitly 0). Both are the same Newton's-method step algebraically, but the sign is folded differently — a naive "unify into one helper" refactor must preserve the effective sign in each call site, not just copy one formula. If a shared helper is written (per this repo's `CLAUDE.md` constraint on not duplicating logic), it should probably be `newton_step(target: f64, current: f64, current_speed: f64) -> f64` returning `swe_difdeg2n(target, current) / current_speed` for longitude cases, and the node case can call it as `newton_step(0.0, x[1], x[4])` *if* `swe_difdeg2n(0.0, x[1])` is proven equal to `-x[1]` for the small residuals actually encountered (true for `|x[1]| < 180`, which always holds for ecliptic latitude) — otherwise keep them separate to avoid a forced abstraction (per `CLAUDE.md`: "If the inputs diverge enough that a shared abstraction would be forced or awkward, that's fine — use judgment").
- **No max-iteration guards anywhere**: all five iteration loops (`solcross`, `mooncross`, `mooncross_node`'s two phases, `helio_cross`) are unbounded C `for(;;)` loops relying purely on Newton convergence plus the inner `swe_calc` error path as the only escape hatch. A faithful port could still add a sane iteration cap (e.g. 50) purely as a defensive measure against non-convergence in the Rust port, since C's "infinite loop that never converges" is not a behavior worth reproducing — but this is a judgment call, not a fidelity requirement, since the underlying math (given correct ephemeris data) always converges quadratically in a handful of iterations for these smooth, nearly-linear-speed bodies.
- **UT wrapper pattern**: all four `_ut` variants are byte-for-byte structural copies of their ET counterparts with `swe_calc`→`swe_calc_ut` and `jd_et`→`jd_ut` substituted — no deltaT computation appears explicitly in this file at all (contrast with other UT wrappers elsewhere in the codebase that manually do `jd_et = jd_ut + deltat(jd_ut)`). `swe_calc_ut` itself is responsible for the ET/UT conversion internally. Port each `_ut` function as a thin wrapper delegating to `Ephemeris::calc_ut` exactly as the ET version delegates to `Ephemeris::calc` — do not introduce a manual deltaT calculation in `crossings.rs`.
- **Suggested module**: new `src/crossings.rs`, exposed through `calc.rs`/`Ephemeris` per this repo's architecture rule ("Application modules... go through `calc.rs`, not directly to backends" — `CLAUDE.md`). Six public entry points mirroring the C names: `solcross`/`solcross_ut`, `mooncross`/`mooncross_ut`, `mooncross_node`/`mooncross_node_ut` (returning `Result<(f64, f64, f64), Error>` — jd, xlon, xlat — or a small struct, replacing the two out-pointers), and `helio_cross`/`helio_cross_ut` (returning `Result<f64, Error>`, replacing the `int32` + out-pointer pair).
