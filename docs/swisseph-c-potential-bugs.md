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

## 3. Pluto uses Mercury's mass ratio in get_gmsm / nod_aps (stale ipl_to_elem mapping)

**Location:** `swecl.c:5074` (`ipl_to_elem[]` table); consumed at `swecl.c:5715` and
`swecl.c:5706–5711` (`get_gmsm`, orbital elements) and `swecl.c:5272` (`swe_nod_aps`
osculating branch)

**What:** `ipl_to_elem[15] = {2,0,0,1,3,4,5,6,7,0,0,0,0,0,2}` maps SE_* body ids to rows of
the mean-element tables (`el_node[8][4]` etc., Mercury=0..Neptune=7). Those tables have no
Pluto row — `swe_nod_aps`'s mean branch excludes Pluto — so `ipl_to_elem[SE_PLUTO] = 0`
(Mercury). But `get_gmsm` and the `swe_nod_aps` osculating branch reuse the *same* table to
index the 9-row `plmass[]` array (swecl.c:5063), which *does* have a Pluto row at index 8
(1/136566000). Net effect:

- **Two-body GM** (default): Pluto's `plm = 1/plmass[0] = 1/6023600` (Mercury's Sun/planet
  mass ratio, ≈1.66e-7) instead of `1/plmass[8] ≈ 7.32e-9`. Relative error in `Gmsm` ≈1.6e-7.
- **`SEFLG_ORBEL_AA` branch:** the descending summation loop
  `for (j = ipl; j >= SE_MERCURY; j--) plm += 1/plmass[ipl_to_elem[j]]` double-counts
  Mercury for Pluto (once via the stale mapping at j=SE_PLUTO, once at j=SE_MERCURY) instead
  of adding Pluto's own mass. Same ~1.6e-7 scale.
- **Osculating nodes/apsides** (swecl.c:5272): Pluto's `Gmsm` carries the same wrong mass.

**Impact:** ~1.6e-7 relative in GM → same order in vis-viva semi-major axis; for Pluto at
~39.5 AU that's ~6e-6 AU (≈900 km) in derived element distances. Real but astronomically
negligible, far below any astrological use.

**Cause:** Table reuse across mismatched shapes — a lookup table built for the 8-row mean-
element arrays repurposed for the 9-row mass array without extending the Pluto entry. No
comment suggests intent.

**Our Rust code:** Replicated bit-for-bit (`constants::IPL_TO_ELEM` transcribed verbatim,
PNOC 4/6 — swisseph-rs/85, /87). Golden-test fidelity requires it; see
`docs/c-ref-orbital-elements.md` Constants § quirk note.

## 4. swe_calc_pctr computes light deflection as observed from Earth, regardless of center body

**Location:** `sweph.c:8150–8152` (`swe_calc_pctr` deflection step) vs. `sweph.c:8153–8168`
(aberration step)

**What:** In the planetocentric pipeline, `swi_deflect_light(xx, ...)` takes no observer
parameter — it reads Earth's and the Sun's barycentric position/velocity from the global
cache (`swed.pldat[SEI_EARTH].x`, `swed.pldat[SEI_SUNBARY].x`) and applies the standard
Sun-deflection formula treating the planetocentric vector `xx` as if it were geocentric. So
for "Mars as seen from Jupiter", gravitational light bending is computed as though the
observer were Earth. Annual aberration in the very next step, by contrast, correctly uses
the center body's velocity (`swi_aberr_light(xx, xxctr, ...)`) — the two relativistic
corrections in the same function disagree about who the observer is.

**Impact:** Deflection is ≤1.8 arcsec at the solar limb and typically micro-arcseconds away
from it; the observer-geometry error is a fraction of that. Negligible for any published use
of planetocentric positions, but it is physically wrong for non-Earth centers.

**Cause:** Global-state architecture again — `swi_deflect_light` was written for the
geocentric pipeline and hard-reads the Earth/Sun save slots; `swe_calc_pctr` (a much later
addition) reuses it without threading the center body through.

**Our Rust code:** Replicates the asymmetry deliberately: `calc_pctr` (PNOC 9 —
swisseph-rs/90) calls `corrections::deflect_light` with Earth-observer geometry and
`corrections::aberr_light` with the center body's velocity, per `docs/c-ref-pctr.md` §5–6.

## 5. swe_orbit_max_min_true_distance rough scan covers only half the inner ellipse

**Location:** `swecl.c:6207–6253` (rough grid scan in `swe_orbit_max_min_true_distance`)

**What:** The rough scan samples the outer ellipse's eccentric anomaly at `eano = j*dstep`
(`j=0..181`, `dstep=2` → 0..362°, one harmless duplicate step past 360°), but the inner
ellipse at `eani = (double)i` — `dstep` is never applied — so the inner body is only sampled
over 0..181°. Half the inner ellipse (182°–359°) is never visited in the rough scan. The
asymmetry looks like a straightforward oversight (`dstep` presumably intended for both
loops).

**Impact:** Usually none — the subsequent block-coordinate refinement (swecl.c:6254–6282,
up to 301 alternating passes, 1e-8 AU convergence) walks to a local extremum from whichever
rough candidate was found, and for realistic near-coplanar, low-eccentricity orbit pairs
the true extremes are recoverable from the sampled half. But for geometries whose extremum
lies in the unsampled half of the inner ellipse with a competing local extremum elsewhere,
the refinement can converge to the wrong local extremum. `dmax`/`dmin` would then be
silently wrong in C and (by design) in our port.

**Cause:** Apparent loop-variable oversight; no comment suggests the asymmetry is
intentional.

**Our Rust code:** Ported literally (PNOC 6 — swisseph-rs/87) so the refinement starts from
the same rough candidates C finds; see `docs/c-ref-orbital-elements.md` §8.2's loop-bound
quirk note.
