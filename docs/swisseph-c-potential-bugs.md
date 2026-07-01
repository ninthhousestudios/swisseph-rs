# Potential bugs in upstream Swiss Ephemeris (C)

Observations from porting. None of these have practical impact at the library's
intended precision, but they're documented for completeness.

## 1. embofs_mosh uses cached obliquity for backward-difference time point

**Location:** `swemplan.c`, `embofs_mosh()` (line ~416)

**What:** `embofs_mosh` reads `swed.oec.seps` / `swed.oec.ceps` — the obliquity
cached from the most recent `swe_calc` call. When `swi_moshplan` computes Earth
velocity by backward difference, it calls `embofs_mosh` twice:

```c
embofs_mosh(tjd, xe);                        // main position
embofs_mosh(tjd - PLAN_SPEED_INTV, x2);     // backward-difference position
```

Both calls read the same `swed.oec` (obliquity at `tjd`). The second call should
use the obliquity at `tjd - PLAN_SPEED_INTV` (0.0001 days earlier).

**Impact:** Over 0.0001 days, obliquity changes by ~6e-16 radians. Applied to the
Moon at ~4.3e-5 AU, position error is ~2.6e-20 AU (~4 picometers). The short Moon
series in `embofs_mosh` is accurate to ~1 arcminute, so this error is 12 orders of
magnitude below the series precision. Zero practical impact.

**Cause:** Global state architecture — `embofs_mosh` is a `static` function that
reads from the `swed` global rather than taking obliquity as a parameter. No
indication this is intentional; no comment suggesting deliberate parameter-holding
for backward-difference stability.

**Our Rust code:** Uses correct per-time-point obliquity. Tests pass at 1e-10
tolerance since the effect is far below that threshold.

## 2. swe_houses vs. swe_houses_ex2 use inconsistent deltaT tidal-acceleration policies

**Location:** `swehouse.c`, `swe_houses()` (line ~139) vs. `swe_houses_ex2()` (line ~220)

**What:** Both functions independently compute `tjde = tjd_ut + swe_deltat_ex(tjd_ut, X, NULL)`
for their own ARMC/obliquity/nutation setup (`swe_houses` does not call `swe_houses_ex2` —
it duplicates the setup and calls `swe_houses_armc_ex2` directly). They differ in what they
pass as `X`:

```c
/* swe_houses, line 139 */
double tjde = tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL);
/* swe_houses_ex2, line 220 */
double tjde = tjd_ut + swe_deltat_ex(tjd_ut, iflag, NULL);
```

`-1` forces `swi_get_tid_acc` down its explicit "default tid_acc" path
(`denum=9999` → switch-default → `SE_TIDAL_DEFAULT`, i.e. DE431's value).
`swe_houses_ex2` instead forwards the caller's raw `iflag`. For a typical house call,
`iflag` only carries `SEFLG_SIDEREAL`/`SEFLG_NONUT`/`SEFLG_RADIANS` — no ephemeris-source
bit — so `swi_get_tid_acc` falls through to the *same* switch-default and the two functions
coincidentally agree. But this is an artifact of the switch-case fallthrough, not a designed
invariant: if a caller passes `SEFLG_MOSEPH` (or `SEFLG_SWIEPH`) in `iflag` to
`swe_houses_ex2` — plausible for a wrapper library that threads the same flags through every
SE call for consistency — `swe_houses_ex2` picks up DE404 (or the SWIEPH file's real denum)
while `swe_houses` on identical `(tjd_ut, geolat, geolon, hsys)` still gives DE431. Two
functions that are meant to be "the same calculation, `_ex2` just exposes more knobs" can
silently disagree.

**Sharper version, same root cause:** within a single `swe_houses_ex2('I'/'i', ...)`
(Sunshine) call, the Sun's declination comes from `swe_calc_ut(tjd_ut, SE_SUN,
SEFLG_SPEED|SEFLG_EQUATORIAL, ...)`, whose deltaT is resolved via `plaus_iflag`/
`swi_guess_ephe_flag` and typically picks up whatever backend is *actually* configured
(real SWIEPH file, real denum). The ARMC's own deltaT, per above, resolves independently
and usually lands on the hardcoded default instead. So the Sun's position and the sidereal
time used to place it relative to the horizon can be computed as of very slightly different
effective clocks within the same call — an internal inconsistency, not just an
API-symmetry gripe.

**Impact:** DeltaT differences between tidal-acceleration models are sub-arcsecond in
ARMC/house-cusp terms even centuries from J2000 (observed ~3e-6° divergence at a
1800-Jan-1 test case during our port — see below). Astronomically inert, same order as
finding #1 and this project's documented stateless-tolerance boundary artifacts
(`CLAUDE.md` `<stateless_tolerance>`).

**Verification status:** Traced from the C source, not empirically run — I have not tested
`swe_houses_ex2` with `SEFLG_MOSEPH` stuffed into `iflag` against a live SWIEPH-file session
to confirm the divergence actually manifests as described. Worth an empirical check before
treating this as confirmed rather than "strongly suggested by the code."

**Our Rust code:** `Ephemeris::houses_ex2` (`src/context.rs`) forces
`tidal_acceleration = Some(TIDAL_DEFAULT)` for its ARMC/eps deltaT computation, matching
C's actual (if seemingly accidental) behavior for the tested cases — this is the correct
fidelity target regardless of whether the C behavior is intentional. Found via golden test
`houses::ut_wrapper` (swisseph-rs/65); see that task's execution_record and the sutra lesson
anchored on `houses_ex2`/`calc_deltat`/`resolve_tidal_acceleration` for the full account.
