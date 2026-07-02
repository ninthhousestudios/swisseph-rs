# C Reference: Planetocentric Positions — `swe_calc_pctr` (sweph.c)

Porting reference for `swe_calc_pctr(tjd, ipl, iplctr, iflag, xxret, serr)`: the position of
body `ipl` as seen from center body `iplctr` (e.g. "Mars as seen from Jupiter"). Read this
instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_calc_pctr` | sweph.c:8042–8283 | Yes — this doc |
| `plaus_iflag` | sweph.c:6066–6148 | Already ported: `crate::calc::plaus_iflag` (calc.rs:29) |
| `swi_aberr_light` | sweph.c:3699–3736 | Already ported: `corrections::aberr_light` |
| `swi_deflect_light` | sweph.c:3743–3891 | Already ported: `corrections::deflect_light` |
| `swi_bias` (frame bias) | swephlib.c | Already ported: `bias::frame_bias` |
| `swi_precess` / `swi_precess_speed` | swephlib.c | Already ported: `precession::precess` |
| `swi_nutate` | sweph.c:3592 | Nutation-rotation tail already ported inside `calc::app_pos_rest` |
| `swi_coortrf2`, `swi_cartpol_sp`, `swi_polcart_sp` | swephlib.c | Already ported: `math::rotate_x_sincos`, `math::cartesian_to_polar_with_speed`, `math::polar_to_cartesian_with_speed` |
| `swi_get_denum` | sweph.c:2407–2443 | Not needed standalone — absorbed into `AstroModels`/`models.bias` construction (see Porting Notes) |
| `swi_get_ayanamsa_with_speed`, `swi_trop_ra2sid_lon*` | sidereal path | Existing sidereal module (`context.rs::apply_sidereal`) — see §9 |

`square_sum(x) = x[0]²+x[1]²+x[2]²`, `dot_prod(x,y) = x·y` (sweph.h:308–309) — plain dot products,
no special evaluation order.

## 0. Input Validation

```c
if (ipl == iplctr) {
    sprintf(serr, "ipl and iplctr (= %d) must not be identical\n", ipl);
    return ERR;
}
```
(sweph.c:8050–8054). This is the **only** iplctr-specific validation in the function. No other
combination of `ipl`/`iplctr` is rejected — e.g. `iplctr = SE_MOON`, `iplctr = SE_SUN`, or asteroid
center bodies are all accepted; `ipl`/`iplctr` validity for the underlying ephemeris lookup is
handled by the ordinary `swe_calc` calls further down (they can return `ERR` and propagate).

`iflag` is then normalized via `plaus_iflag(iflag, ipl, tjd, serr)` (sweph.c:8055) — identical
general-purpose flag normalization used by `swe_calc`/`swe_calc_ut` (topo→turns off helio/bary,
bary turns off helio, helio/bary force NOABERR+NOGDEFL, J2000 forces NONUT, sidereal forces NONUT,
truepos forces NOGDEFL+NOABERR, ephemeris-source bit resolved/defaulted). No pctr-specific
branches inside `plaus_iflag` itself.

## 1. Priming Obliquity/Nutation State (sweph.c:8056–8058)

```c
epheflag = iflag & SEFLG_EPHMASK;
swe_calc(tjd + swe_deltat_ex(tjd, epheflag, serr), SE_ECL_NUT, iflag, xx, serr);
```

This call's *only* purpose is a side effect: it populates the global obliquity/nutation caches
(`swed.oec`, `swed.oec2000`, `swed.nut`) that are read later (§6, §7) — its own `xx` output is
discarded. **The epoch used is `tjd + Δt(tjd)`, not `tjd` and not the light-time-retarded `t`
computed later** — a third, distinct epoch from the other two used in the pipeline. This looks
like `tjd` is expected to already be TT/ET (per `swe_calc_pctr`'s documented contract) and the
extra `+Δt` shift is what the C source actually does — replicate it literally for fidelity, don't
"fix" it as if it were a UT/TT mismatch bug.

If `SEFLG_J2000` is requested, `plaus_iflag` already forced `SEFLG_NONUT`, and the ecliptic
reference (`oe = &swed.oec2000`, §6) is the epoch-independent J2000 mean obliquity constant — so
the odd epoch offset above only actually matters when nutation-of-date is requested.

## 2. Barycentric J2000 Computation of Both Bodies (sweph.c:8059–8072)

```c
iflag &= ~(SEFLG_HELCTR|SEFLG_BARYCTR);
iflag2  = epheflag;
iflag2 |= (SEFLG_BARYCTR|SEFLG_J2000|SEFLG_ICRS|SEFLG_TRUEPOS|SEFLG_EQUATORIAL|SEFLG_XYZ|SEFLG_SPEED);
iflag2 |= (SEFLG_NOABERR|SEFLG_NOGDEFL);
retc = swe_calc(tjd, iplctr, iflag2, xxctr, serr);  /* barycentric, J2000, ICRS, true, equatorial XYZ, speed */
retc = swe_calc(tjd, ipl,    iflag2, xx,    serr);
for (i = 0; i <= 5; i++) xx0[i] = xx[i];            /* xx0 = undisturbed barycentric ipl @ tjd */
```

Note `iflag` (the user's flag, minus HELCTR/BARYCTR) is kept around for later gating; `iflag2` is
a completely separate, hard-forced flag set used **only** for these internal `swe_calc` calls —
it always requests barycentric, J2000, ICRS, **true** position (no light-time/aberration/deflection
inside this inner call — those are done manually below), equatorial, XYZ, with speed, and
explicitly turns off aberration/deflection (which TRUEPOS already implies, but both bits are set
redundantly). Both bodies are fetched at the *unretarded* `tjd`.

If either `swe_calc` call returns `ERR`, `swe_calc_pctr` returns `ERR` immediately (sweph.c:8064–8068).

## 3. Light-Time Iteration, Center Body as Observer (sweph.c:8073–8125)

Gated on `!(iflag & SEFLG_TRUEPOS)`. `niter = 1` (constant — **2 passes total**, `j = 0..=1`).

### 3a. SPEED-only pre-pass: "change of dt" correction seed (sweph.c:8079–8104)

Only runs if `SEFLG_SPEED`. Computes `xxsp` = the daily-motion correction that will later be
subtracted from the planetocentric speed (§5) to account for light-time itself changing over a day:

```c
for (i=0..2) xxsv[i] = xxsp[i] = xx[i] - xx[i+3];        /* "yesterday's" rough position of ipl */
for (j = 0; j <= niter; j++) {
  for (i=0..2) { dx[i] = xxsp[i]; dx[i] -= (xxctr[i] - xxctr[i+3]); }   /* relative to yesterday's iplctr */
  dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
  for (i=0..2) xxsp[i] = xxsv[i] - dt * xx0[i+3];          /* refine "yesterday" apparent pos */
}
for (i=0..2) xxsp[i] = xxsv[i] - xxsp[i];   /* true("yesterday") - apparent("yesterday") */
```

`xxsv`/`xx0` use `ipl`'s **undisturbed barycentric position/velocity from §2** throughout (not
`iplctr`'s), i.e. this is estimating how the *ipl-relative-to-iplctr* light-time changes over one
day using a crude "one day ago" state built from `xx0`'s own velocity, not a true previous-day
ephemeris evaluation.

### 3b. Main dt/apparent-time loop (sweph.c:8105–8117)

```c
for (j = 0; j <= niter; j++) {
  for (i=0..2) { dx[i] = xx[i]; dx[i] -= xxctr[i]; }
  dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
  t = tjd - dt;
  dtsave_for_defl = dt;
  for (i=0..2) xx[i] = xx0[i] - dt * xx0[i+3];   /* rough apparent (light-time-adjusted) ipl position */
}
```

`dx` is computed against `xxctr` — the **center body's** barycentric position (not Earth) —
because `iplctr` is the observer here. `xx0`/`xx0[i+3]` (ipl's own barycentric pos/vel from §2)
are reused every iteration; only `dt`/`t` are refined. After the loop, `dtsave_for_defl` holds the
final `dt` — this exact value, and no other, is what's later passed into deflection (§7).

### 3c. Finalize SPEED change-of-dt correction (sweph.c:8118–8122)

```c
if (iflag & SEFLG_SPEED)
  for (i=0..2) xxsp[i] = xx0[i] - xx[i] - xxsp[i];
```

### 3d. Re-evaluate both bodies at retarded time `t` (sweph.c:8123–8124)

```c
retc = swe_calc(t, iplctr, iflag2, xxctr2, serr);   /* iplctr at retarded time */
retc = swe_calc(t, ipl,    iflag2, xx,     serr);   /* ipl at retarded time — overwrites xx */
```

Same `iflag2` (barycentric/J2000/ICRS/true/equatorial-XYZ/speed) as §2. **This is the last
`swe_calc` call in the function** — as a side effect it repopulates the C global planet-data cache
(`swed.pldat[SEI_EARTH]`, `swed.pldat[SEI_SUNBARY]`) to epoch `t`, which §7 silently depends on.

## 4. Conversion to Planetocentric (sweph.c:8129–8146)

```c
if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR)) {
  for (i=0..5) xx[i] -= xxctr[i];
  if (!(iflag & SEFLG_TRUEPOS) && (iflag & SEFLG_SPEED))
    for (i=3..5) xx[i] -= xxsp[i-3];
}
if (!(iflag & SEFLG_SPEED))
  for (i=3..5) xx[i] = 0;
```

**Quirk:** if the caller's *original* `iflag` (before line 8059 stripped it — but note: the check
here is against the already-stripped local `iflag`, which no longer carries HELCTR/BARYCTR since
line 8059 removed them from `iflag` itself!). Re-reading the C: line 8059 does
`iflag &= ~(SEFLG_HELCTR|SEFLG_BARYCTR)`, so by line 8129 `iflag` **never** has those bits set —
the `!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR)` guard is therefore **always true** in
practice; the subtraction always happens. (The guard reads as dead/defensive code — a porter
should still replicate the always-true branch, but need not implement a "return raw barycentric
ipl" mode gated on these flags, since it's unreachable given line 8059.)

## 5. Gravitational Deflection (sweph.c:8150–8152)

```c
if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOGDEFL))
  swi_deflect_light(xx, dtsave_for_defl, iflag);
```

`xx` at this point is the **ipl-minus-iplctr** vector (planetocentric geometric position), not a
true geocentric vector. **`swi_deflect_light` ignores this distinction entirely**: it does not
take `iplctr`/`xxctr` as a parameter at all. It reads Earth's and the Sun's barycentric
position/velocity from the C global cache (`swed.pldat[SEI_EARTH].x`, `swed.pldat[SEI_SUNBARY].x`
— last written by the §3d call at epoch `t`) and applies the standard Sun-deflection formula
treating `xx` as if it were geocentric. **This means gravitational light deflection in
`swe_calc_pctr` is always computed using the true Earth/Sun geometry, regardless of what
`iplctr` actually is** — e.g. for "Mars as seen from Jupiter", deflection is still computed as
though observed from Earth. This is a faithful description of what the C reference implementation
does, not an approximation to relax in the port — golden-data fidelity requires reproducing it.
See Porting Notes for the stateless equivalent.

The deflection algorithm itself (u/e/q unit vectors, `meff()` solar-limb softening, `HELGRAVCONST`,
speed correction via `DEFL_SPEED_INTV`) is unchanged from the standard pipeline and is already
ported in `corrections::deflect_light` (sweph.c:3743–3891) — no need to re-derive it here.

## 6. Annual Aberration (sweph.c:8153–8168)

```c
if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOABERR)) {
  swi_aberr_light(xx, xxctr, iflag);
  if (iflag & SEFLG_SPEED)
    for (i=3..5) xx[i] += xxctr[i] - xxctr2[i];
}
if (!(iflag & SEFLG_SPEED))
  for (i=3..5) xx[i] = 0;
```

Unlike deflection, aberration correctly uses **`xxctr`** (the center body's velocity, at
unretarded `tjd`) as the observer — `swi_aberr_light(xx, xe, iflag)` only ever reads `xe[3..5]`
(velocity), never `xe[0..2]` (position doesn't matter to the formula). The speed correction term
uses `xxctr - xxctr2` (center body's velocity at `tjd` minus at retarded `t`) — the analogous
"observer velocity changed between emission and reception" correction as the standard pipeline,
but keyed on the center body instead of Earth. `swi_aberr_light` itself is already ported as
`corrections::aberr_light` (sweph.c:3699–3736) with signature `(xx, earth_vel: &[f64;3], has_speed)`
— for pctr, pass `xxctr`'s velocity (not Earth's) as `earth_vel`.

## 7. ICRS → J2000 Frame Bias (sweph.c:8172–8175)

```c
if (!(iflag & SEFLG_ICRS) && swi_get_denum(ipl, epheflag) >= 403)
  swi_bias(xx, t, iflag, FALSE);
```

Uses the **retarded** time `t` (from §3b), not `tjd`. `swi_get_denum` returns 403 unconditionally
for Moshier (`SEFLG_MOSEPH` → hardcoded `return 403`, sweph.c:2410–2411), so for the Moshier-only
port this condition reduces to simply `!ICRS`. Already ported as `bias::frame_bias`; the
`denum >= 403` gate for file-based backends is absorbed into `AstroModels`/`models.bias`
construction elsewhere in the codebase — don't re-derive it inline here.

`xxsv[0..5] = xx` is saved immediately after (sweph.c:8177–8178) as the J2000 coordinates needed
later for sidereal projection (§9).

## 8. Precession (sweph.c:8182–8189)

```c
if (!(iflag & SEFLG_J2000)) {
  swi_precess(xx, tjd, iflag, J2000_TO_J);
  if (iflag & SEFLG_SPEED) swi_precess_speed(xx, tjd, iflag, J2000_TO_J);
  oe = &swed.oec;       /* obliquity of date, epoch tjd+Δt(tjd) from §1 */
} else {
  oe = &swed.oec2000;   /* J2000 mean obliquity, epoch-independent */
}
```

Precession uses the **original** (unretarded) `tjd`, not `t` — asymmetric with §7's use of `t`,
but this asymmetry matches the standard `swe_calc` pipeline too (sweph.c:2758 uses `t` for bias,
2767 uses `pdp->teval` i.e. the reception epoch for precession) — not pctr-specific.

## 9. Nutation, Ecliptic Conversion, Polar Conversion, Degrees (sweph.c:8193–8257)

```c
if (!(iflag & SEFLG_NONUT)) swi_nutate(xx, iflag, FALSE);
xreturn[18..23] = xx;                                    /* equatorial cartesian, save */
swi_coortrf2(xx, xx, oe->seps, oe->ceps);                /* equatorial -> ecliptic (obliquity rotation) */
if (SPEED) swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
if (!(iflag & SEFLG_NONUT)) {
  swi_coortrf2(xx, xx, swed.nut.snut, swed.nut.cnut);    /* nutation-in-obliquity rotation */
  if (SPEED) swi_coortrf2(xx+3, xx+3, swed.nut.snut, swed.nut.cnut);
}
xreturn[6..11] = xx;                                     /* ecliptic cartesian, save */
```

This block is **structurally identical** to the tail of the ordinary `swe_calc` pipeline
(`swi_nutate` → obliquity rotation → nutation rotation → save ecliptic cartesian) — it is exactly
what `calc::app_pos_rest` (calc.rs:251–327) already implements. `swi_nutate` here reads
`swed.nut` populated by §1's priming call (epoch `tjd+Δt(tjd)`), i.e. it does **not** recompute
nutation from `xx`'s own epoch — genuinely global-state-derived, no local recomputation exists in
this function to fall back on.

Sidereal handling (sweph.c:8217–8243, gated on `SEFLG_SIDEREAL`) branches into the same three
sub-algorithms as the standard pipeline (ecl-T0 projection / solar-system-plane projection /
traditional ayanamsha subtraction via `swi_get_ayanamsa_with_speed`) — reuse the existing sidereal
module (`context.rs::apply_sidereal`) rather than reimplementing; no pctr-specific sidereal logic
exists.

Polar conversion + degrees (sweph.c:8247–8257):
```c
swi_cartpol_sp(xreturn+18, xreturn+12);   /* equatorial cartesian -> polar */
swi_cartpol_sp(xreturn+6,  xreturn);      /* ecliptic cartesian -> polar */
for (i=0;i<2;i++) {
  xreturn[i]    *= RADTODEG; xreturn[i+3]  *= RADTODEG;   /* ecliptic lon/lat, lon/lat speed */
  xreturn[i+12] *= RADTODEG; xreturn[i+15] *= RADTODEG;   /* equatorial RA/Dec, RA/Dec speed */
}
```
Distances (`xreturn[2]`, `xreturn[5]`, `xreturn[8]`, `xreturn[11]`, `xreturn[14]`, `xreturn[17]`,
`xreturn[20]`, `xreturn[23]`) are never degree-scaled — matches `app_pos_rest`'s "angles only" step.

## 10. Output Slot Layout and Final Selection (sweph.c:8258–8282)

`xreturn[24]` layout (identical to the standard `xxsv`/`xreturn` convention used elsewhere in
sweph.c):

| Slice | Content |
|---|---|
| `xreturn[0..6]` | Ecliptic polar: lon, lat, dist, lon-speed, lat-speed, dist-speed (degrees) |
| `xreturn[6..12]` | Ecliptic cartesian (equatorial-of-date + obliquity + nutation rotated) |
| `xreturn[12..18]` | Equatorial polar: RA, Dec, dist, RA-speed, Dec-speed, dist-speed (degrees) |
| `xreturn[18..24]` | Equatorial cartesian (post nutation, pre ecliptic rotation) |

Final selection into the caller's `xxret[6]`:
```c
xs = (iflag & SEFLG_EQUATORIAL) ? xreturn+12 : xreturn;
if (iflag & SEFLG_XYZ) xs += 6;             /* switch polar block -> cartesian block */
xxret[0..6] = xs[0..6];
if (!(iflag & SEFLG_SPEED)) xxret[3..6] = 0;
if (iflag & SEFLG_RADIANS) {
  xxret[0..2] *= DEGTORAD;
  if (SPEED) xxret[3..5] *= DEGTORAD;        /* note: index 5 (distance-speed) is NOT converted */
}
```
Distance/dist-speed (`xxret[2]`, `xxret[5]`) are always in AU regardless of `SEFLG_RADIANS` —
only the two angle components and their speeds are unit-converted. Return value is `iflag` (the
normalized flags actually used) on success, `ERR` if any of the internal `swe_calc` calls failed
(checked via the last `retc`, sweph.c:8280–8281 — note only the *last* `swe_calc(t, ipl, ...)`
call's `retc` is actually checked here; the two earlier calls already `return ERR` immediately on
failure at sweph.c:8064–8068, so this final check only matters for the retc set at line 8124).

## Constants

| Name | Value | Notes |
|---|---|---|
| `AUNIT` | 1.49597870700e+11 m | sweph.h:273 (DE431 value) |
| `CLIGHT` | 2.99792458e+8 m/s | sweph.h:274 |
| `HELGRAVCONST` | 1.32712440017987e+20 m³/s² | sweph.h:278 (used inside `deflect_light`) |
| `SUN_RADIUS` | `959.63/3600 * DEGTORAD` rad | sweph.h:281 (used inside `deflect_light`) |
| `DEFL_SPEED_INTV` | 0.0000005 day | sweph.h:304 (used inside `deflect_light`) |
| `niter` | 1 | fixed — 2 total passes in both light-time loops (§3b), no convergence check, no adaptive iteration |
| `iflag2` forced bits | `BARYCTR\|J2000\|ICRS\|TRUEPOS\|EQUATORIAL\|XYZ\|SPEED\|NOABERR\|NOGDEFL` plus resolved `epheflag` | sweph.c:8060–8062 |
| `DEGTORAD`/`RADTODEG` | π/180, 180/π | sweodef.h:262–266 |

## Porting Notes

**Global-state reads and their stateless equivalents:**

1. **`swed.pldat[SEI_EARTH].x` / `swed.pldat[SEI_SUNBARY].x`** (read inside `swi_deflect_light`,
   §5): C reuses whatever barycentric Earth/Sun position was cached as a side effect of the *last*
   `swe_calc` call (§3d, at epoch `t`) — regardless of `iplctr`. To bit-match this in the stateless
   port, `calc_pctr` must **explicitly compute barycentric Earth and Sun at epoch `t`** (a real
   `Ephemeris::calc(t, Body::Earth, BARYCTR|...)` / equivalent Sun-barycentric call, independent of
   `iplctr`) and pass *those* into `deflect_light`'s `earth_helio`/observer argument — **not**
   `xxctr`/`xxctr2`. This mirrors the same known stateless-tolerance pattern already documented in
   the project (`swi_deflect_light` reading `psdp->x` from cache, see CLAUDE.md
   `<stateless_tolerance>` §1) but here it's not just a speed-precision nuance — it changes *which
   body* the deflection geometry is built from. Get this wrong and pctr's deflection will silently
   use the wrong observer for any `iplctr != Earth`.

2. **`swed.oec` / `swed.oec2000` / `swed.nut`** (§1, §8, §9): populated once by the priming
   `swe_calc(tjd+Δt(tjd), SE_ECL_NUT, ...)` call and read later without recomputation. Stateless
   equivalent: compute obliquity(`tjd+Δt(tjd)`) and nutation(`tjd+Δt(tjd)`) once at the top of
   `calc_pctr`, thread the resulting `Epsilon`/`Nutation` values through to the precession/nutation
   tail (§8–§9) exactly as `calc::app_pos_rest` already expects them as parameters — do **not**
   recompute obliquity/nutation at `tjd` or at `t`.

3. **`swi_get_denum(ipl, epheflag) >= 403`** (§7): for Moshier this is always true; for file-backed
   sources it's already captured by `AstroModels`/`models.bias` at `Ephemeris` construction time —
   reuse that, don't add a parallel denum check.

**`CalcFlags::CENTER_BODY` (flags.rs:26, value `1<<20`) is unrelated to this feature.** It ports
C's `SEFLG_CENTER_BODY` (swephexp.h:216), a *different* mechanism used inside `main_planet`/
`app_pos_etc_plan` to select a planet's photocenter vs. system-barycenter (e.g. Pluto vs.
Pluto-system barycenter) for objects with moons (sweph.c:417–441, 2445 `calc_center_body`). It has
no relationship to `iplctr`/planetocentric positions and should not be reused or conflated when
designing `calc_pctr`'s signature — `iplctr` must be a separate `Body` parameter (e.g.
`Ephemeris::calc_pctr(jd_tt, body: Body, center: Body, flags: CalcFlags)`), not a flag bit.

**Reuse, don't reimplement:**
- `crate::calc::plaus_iflag` for the top-level flag normalization (§0).
- `corrections::aberr_light` for §6 (pass `xxctr`'s velocity as the "earth_vel" argument).
- `corrections::deflect_light` for §5 (pass the independently-computed Earth/Sun barycentric
  state per note 1 above, plus the planet's own heliocentric-at-`t` position exactly as
  `calc_planet` already does for the standard pipeline at calc.rs:398–413).
- `calc::app_pos_rest` for §9's nutation/ecliptic-rotation/polar-conversion/degrees tail — the
  logic is identical field-for-field to the standard pipeline's tail; only the `Epsilon`/`Nutation`
  values fed in (from note 2) differ.
- `bias::frame_bias` for §7.
- `context.rs::apply_sidereal` for the sidereal branch (§9) instead of re-deriving
  `swi_trop_ra2sid_lon`/`swi_trop_ra2sid_lon_sosy`/ayanamsha subtraction.

**Not a real convergence loop:** §3b's "light-time iteration" is a fixed 2-pass linear
extrapolation (`niter = 1`), not an iterative solver — no epsilon/convergence check exists to
port; just unroll the loop twice (or keep the `for j in 0..=1` loop verbatim for a single point of
truth against the C).

**Dead branch, implement as unconditional:** §4's `!(iflag & SEFLG_HELCTR) &&
!(iflag & SEFLG_BARYCTR)` guard is unreachable-false given line 8059 strips both bits from `iflag`
immediately beforehand — the planetocentric subtraction always executes. A straightforward port
can therefore simplify this to unconditional subtraction (with a comment noting why), rather than
implementing a dead "raw barycentric passthrough" code path.
