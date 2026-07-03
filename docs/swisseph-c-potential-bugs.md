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

## 6. lunar_osc_elem node/apogee speed is inconsistent with its own positions (Moshier)

**Location:** `sweph.c:5360–5471` (`lunar_osc_elem` node/apogee central-difference speed),
via the off-center samples' `swi_plan_for_osc_elem` (`sweph.c:5758`) reads of `swed.oec` /
`swed.nut`.

**What:** For `SE_TRUE_NODE` / `SE_OSCU_APOG`, the speed is a central difference of the
node/apogee *position* over `speed_intv` (`(xx[1] - xx[0]) / speed_intv / 2`). The three
position samples are computed at `t = tjd`, `tjd ± speed_intv` and each is rotated into the
ecliptic of date by `swi_plan_for_osc_elem`, which reads the global obliquity/nutation caches
(`swed.oec`, `swed.nut`) with an equality-on-`teps`/`tnut` fast path. The center sample
(`t == tjd`) hits the cache; the off-center samples (`t == tjd ± speed_intv`) miss it and
recompute `calc_epsilon`/`swi_nutation` fresh. The cached and freshly-computed values round
slightly differently, so the off-center node positions carry a ~2.6e-11 AU offset relative to
the same epoch computed via a clean, cache-free path. The net effect: **C's node/apogee speed
does not match a finite difference of C's *own* node/apogee positions.** Verified at
jd=2305447.5 (Moshier): C's internal ecliptic-cartesian node speed is `-3.36864e-6`, but
differencing C's own standalone `swe_calc` node positions at jd±0.1 gives `-3.36877e-6`
(~1.3e-10 discrepancy).

**Impact:** Moshier-backend node/apogee speed differs from a stateless recomputation by up to
~3.6e-6 °/day (≈0.013 arcsec/day) — astronomically negligible. Positions are unaffected
(bit-identical). The Swiss/JPL backends use a 1000× smaller interval (`NODE_CALC_INTV = 1e-4`
vs `NODE_CALC_INTV_MOSH = 0.1`) and are not measurably affected.

**Cause:** Global-state architecture — `swi_plan_for_osc_elem`'s obliquity/nutation cache
fast path makes the center sample take a different arithmetic route than the off-center ones,
so the finite difference mixes cached and freshly-computed frame rotations.

**Our Rust code:** Stateless — obliquity and nutation are recomputed fresh and consistently at
every sample epoch (PNOC 3 — swisseph-rs/84, `calc::lunar_osc_elem` / `calc::plan_for_osc_elem`).
This matches C's node/apogee *positions* bit-for-bit (1e-9) but cannot reproduce the speed
artifact without re-introducing the cache, so the Moshier speed golden tolerance is relaxed to
5e-6 (Swiss stays 1e-7). See `CLAUDE.md <stateless_tolerance>` §3.

## 7. swe_nod_aps mean descending-node distance singularity → ill-conditioned apparent position

**Location:** `swecl.c`, `swe_nod_aps()` mean branch, node-distance block (lines ~5225–5240).

**What:** For a mean node, the heliocentric *distance* of the node is derived from the eccentric
anomaly and then divided by the cosine of the node's argument:

```c
ea  = atan(tan(-parg * DEGTORAD / 2) * sqrt((1-ecce)/(1+ecce))) * 2;
xna[2] = sema * (cos(ea) - ecce) / cos(parg * DEGTORAD);          /* ascending node  */
ea  = atan(tan((180 - parg) * DEGTORAD / 2) * sqrt((1-ecce)/(1+ecce))) * 2;
xnd[2] = sema * (cos(ea) - ecce) / cos((180 - parg) * DEGTORAD);  /* descending node */
```

`parg` is the argument of perihelion measured from the ascending node. When `parg` (ascending) or
`180 − parg` (descending) approaches ±90°, the divisor `cos(...)` approaches zero and the node
"distance" blows up. For the low-inclination outer planets this happens in practice — e.g.
**Jupiter**: `parg ≈ 273.9°`, so `cos((180 − parg)°) = cos(−93.9°) ≈ −0.067`, and the descending
node comes out at a spurious **6.19 AU — larger than Jupiter's 5.45 AU aphelion**. The point is a
mathematical artifact of the tangent-of-half-angle formula near its pole, not a physical location.

**Impact:** The *raw geometry* is still well-defined and reproducible (a stateless port matches C
bit-for-bit under `SEFLG_TRUEPOS` for every body and every point, descending node included). The
problem only surfaces in the **apparent** output: `swe_nod_aps` then applies light deflection and
aberration to this ill-conditioned vector. Deflection alone and aberration alone each remain
bit-exact, but their *combination* amplifies a ~5e-10 floating-point-ordering difference (from the
deflection speed branch) into a visible divergence — up to ~3.5e-4° in longitude and ~1.2e-2°/day
in longitude speed for Jupiter's descending node. Because the seed is FP-ordering, C's own
reference digits there are not reproducible across compilers/optimizers; the *value itself* is
degraded, independent of any port.

**Cause:** The closed-form node-distance formula has an unguarded `1 / cos(θ)` singularity at
θ → 90°. A distance clamp, or computing the node radius from the orbit's semi-latus rectum instead
(`r = a(1−e²)/(1 + e·cos(ν_node))`), would be numerically stable; C does neither.

**Our Rust code:** Ports the formula verbatim (`nodaps::mean_branch`), so it reproduces C's raw
geometry bit-for-bit but inherits the same ill-conditioning in the apparent output. The
`tests/golden/nodaps.rs` tolerances reflect this: `TRUEPOS` geometry is asserted tight (1e-9 / 1e-8),
apparent ascending-node / perihelion / aphelion at 1e-6, and the apparent descending node relaxed
to 1e-3° / 2e-2°/day. (swisseph-rs/85)

## 8. swe_nod_aps osculating branch: node-direction division amplifies backend FP noise (both nodes)

**Location:** `swecl.c`, `swe_nod_aps()` osculating branch, tangent-line node construction (lines
~5300–5309, "A.4.3" in `docs/c-ref-nodaps.md`).

**What:** The ascending-node *direction* is built from a linear extrapolation of the sampled
position back to the ecliptic plane:

```c
fac = xpos[i][2] / xpos[i][5];             /* z / dz — "time" to cross the ecliptic */
xn[i][j] = (xpos[i][j] - fac * xpos[i][j+3]) * sgn;
xs[i][j] = -xn[i][j];                       /* descending = antipode, same fac */
```

`fac = z / ż` amplifies whatever relative floating-point noise is present in `ż` (the sampled
radial velocity) by `fac`'s own magnitude before it ever reaches the node/apsis geometry. Unlike
bug §7 (mean branch), **both** the ascending and descending directions inherit this noise equally
(they're literal antipodes of the same vector) — descending is empirically worse only because the
downstream ellipse-distance rescale (`rn2/ro2`, A.4.4) tends to have a larger ratio there.

**Verified, not guessed:** feeding C's own dumped `xpos[1]` (the ecliptic-of-date sample at exactly
`tjd_et`, read via temporary instrumentation of this function) directly into the Rust port's
per-sample ellipse formula reproduces C's own `uu`/`cosnode`/`sinnode`/`sinincl` to ~12 significant
digits — the formula itself is a faithful, bit-exact port. The measured divergence traces entirely
to the raw position/velocity sample's ~1e-10..1e-11 relative backend noise (Moshier series
evaluation order, or Swiss/JPL Chebyshev interpolation order), which the `fac` division then
amplifies. Observed worst case in the golden battery: Jupiter's `SE_NODBIT_OSCU_BAR` descending
node at ~1.2e-3° position / ~1.9e-2°/day speed.

**Impact:** Same shape as §7 — the raw `SEFLG_TRUEPOS` sample-at-`tjd_et` position is unaffected
(this is a construction *within* the already-sampled position, not a light-effect chain), but the
node/apsis output as a whole (all four points, not just descending) carries more FP-ordering
sensitivity than the mean branch's non-node points. Perihelion/aphelion are far less affected since
they derive from `uu`/`ny`/`sema`/`ecce` directly, picking up node noise only secondhand through
`uu`'s `cosnode`/`sinnode` term.

**Cause:** Same class as §7 — a division (`z/ż` here, `1/cos(θ)` there) with no width/scale guard,
applied to a already-tiny quantity subject to backend FP noise.

**Our Rust code:** Ports the formula verbatim (`nodaps::osculating_branch`). Golden tolerances
(`tests/golden/nodaps.rs::osc_tolerance`) are tiered by point: perihelion/aphelion at 5e-5°/1e-4°
per day, ascending node at 1e-3°/1e-4° per day, descending node at 2e-3°/3e-2° per day. The "oscu"
battery's pre-1900 epoch is additionally nudged off the sepl_18 `.se1` file boundary (see
`osc_epochs` in `tests/c-gen/gen_nodaps.c`) — a separate, already-documented stateless-vs-stateful
artifact (see `CLAUDE.md` `<stateless_tolerance>` §2), not this one. (swisseph-rs/86)

## 9. MOSEPH + asteroid output depends on process call history (stale SEI_SUNBARY)

**Location:** `sweph.c:2332–2343` (`sweph()`'s heliocentric→barycentric conversion, flag-gated
on `SEFLG_JPLEPH || SEFLG_SWIEPH`), `sweph.c:2345–2352` (`pdp->iephe = psdp->iephe` for
asteroid files), `sweph.c:2517–2521` (`app_pos_etc_plan`'s `SEFLG_HELCTR` branch, gated on
`pdp->iephe`), `swemplan.c:276–339` (`swi_moshplan` — never writes `swed.pldat[SEI_SUNBARY]`),
plus `swi_deflect_light`'s global reads of the barycentric Sun.

**What:** The asteroid calc path (`swecalc`'s minor-planet branch, sweph.c:1064–1101) threads
`swed.pldat[SEI_SUNBARY]` into several consumers, but under `SEFLG_MOSEPH` nothing ever
populates that slot — `main_planet(SEI_EARTH, ..., MOSEPH, ...)` calls `swi_moshplan`, which
writes only `SEI_EARTH`. So every consumer reads whatever the *last* SWIEPH/JPLEPH call in the
process happened to leave there (or zeros in a fresh process):

1. `sweph()`'s helio→bary conversion for the asteroid is skipped entirely (the flag gate is
   false in MOSEPH mode), so the asteroid position stays heliocentric. This happens to be
   *self-consistent* in a fresh process: `app_pos_etc_plan`'s observer is `pedp->x`
   (sweph.c:2528–2543), which under MOSEPH is Moshier's heliocentric Earth — heliocentric
   minus heliocentric is a correct geocentric vector.
2. `pdp->iephe` for the asteroid is copied from `psdp->iephe` (`SEI_SUNBARY`'s), not set to
   the ephemeris actually used — under MOSEPH it reports whatever ephemeris last computed the
   barycentric Sun, or garbage-zero in a fresh process.
3. `app_pos_etc_plan`'s `SEFLG_HELCTR` branch subtracts `swed.pldat[SEI_SUNBARY].x` when
   `pdp->iephe` claims SWIEPH/JPLEPH. With the stale `iephe` from (2), a MOSEPH+HELCTR
   asteroid call in a process that previously ran SWIEPH subtracts a stale barycentric Sun
   (evaluated at the *earlier call's epoch*) from a heliocentric position — a genuine
   frame/epoch mixing error.
4. `swi_deflect_light`/`swi_aberr_light` read the barycentric-Sun global for their geometry —
   zero (Sun at origin: consistent with the heliocentric frame) in a fresh process, a
   wrong-epoch Sun after a prior SWIEPH call.

Net: **the same `swe_calc(tjd, SE_CERES, SEFLG_MOSEPH)` call returns (subtly) different
results depending on which calls preceded it in the process.** In a fresh process the zeros
make the whole path a coherent heliocentric-frame computation and the geocentric output is
correct; after a SWIEPH/JPLEPH call at a different epoch, the HELCTR branch and the
deflection/aberration geometry mix frames and epochs.

**Impact:** Geocentric non-HELCTR output: position core unaffected (frame-consistent either
way); only the mas-level deflection/aberration geometry wobbles with stale state.
MOSEPH+HELCTR asteroid output: wrong by the stale barycentric-Sun vector (up to ~0.01 AU)
whenever a SWIEPH/JPLEPH call preceded it. The reported `iephe` is unreliable in all MOSEPH
asteroid cases.

**Cause:** Global-state architecture — `SEI_SUNBARY` is a shared mutable slot with no
invalidation tied to ephemeris mode, and the asteroid save path launders its `iephe` through
it.

**Our Rust code:** models the *fresh-process* semantics, which is the only deterministic
member of C's behavior family: `MoshierEarthProvider` reports heliocentric Earth with an
all-zero `sun_bary`, so the asteroid stays heliocentric against a heliocentric observer and
deflection sees the Sun at the origin of that frame (swisseph-rs/95 decision; dispatch in
swisseph-rs/101). Golden verification uses a dedicated MOSEPH-only generator process so the C
side is guaranteed fresh (`tests/c-gen/gen_asteroid_moseph.c`, swisseph-rs/102). See
`docs/c-ref-asteroid.md` §1.5.

## 10. Fictitious-planet Kepler starting guess: `pow(x, 1/3)` uses C integer division — the cube root is never computed

**Location:** `swemplan.c:638` (`swi_osc_el_plan`, high-eccentricity initial-guess block for
the Kepler solve; the block spans swemplan.c:630–642)

**What:** For fictitious bodies with `ecce > 0.975` near perihelion (after folding the mean
anomaly into a canonical range, `|M2| < 30°`), C refines the Kepler-equation starting value
with a cube-root-based analytic formula before handing off to `swi_kepler`:

```c
alpha = (1 - ecce) / (4 * ecce + 0.5);
beta  = M2 / (8 * ecce + 1);
zeta  = pow(beta + sqrt(beta*beta + alpha*alpha), 1/3);   /* <-- 1/3 is integer division */
sigma = zeta - alpha / 2;
```

`1` and `3` are both `int`, so `1/3` evaluates to `0` at compile time and the call is
`pow(x, 0) = 1.0` for every `x > 0`. `zeta` is unconditionally `1.0`; the intended cube root
is never taken and the entire refinement collapses to a fixed (and mathematically meaningless)
starting value `sigma = 1 - alpha/2 - 0.078·(…)`. Among the 19 named fictitious bodies, only
Nibiru (ecce = 0.981092) can reach this branch at all.

**Impact:** Practically none on output values: the guess only seeds `swi_kepler`
(swephlib.c:4065–4096), whose Newton branch (`ecce >= 0.4`) iterates uncapped to a 1e-12 rad
tolerance and converges to the same root from the degenerate guess — a correct cube-root
guess would merely converge in fewer iterations. But because Newton stops on a step-size
threshold, the *iteration path* determines the final ULPs of `E`, so a "fixed" cube root can
produce bit-level differences in Nibiru positions versus C.

**Cause:** Classic C integer-division literal bug (`1/3` instead of `1.0/3.0`). The
surrounding formula is a standard high-eccentricity starting-value construction, so the cube
root is clearly intended; no comment suggests otherwise. Harmless enough (thanks to the
uncapped solver) that it has evidently never been noticed upstream.

**Our Rust code:** Not yet landed — the fictitious-planets port is tasks swisseph-rs/122–123.
The port must reproduce the bug (use literal `1.0`, not `x.cbrt()`) for bit parity; mandated in
`docs/c-ref-fictitious.md` §4 + Porting Notes and in the swisseph-rs/122 task description.
