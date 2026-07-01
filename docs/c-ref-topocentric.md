# C Reference: Topocentric Position (SEFLG_TOPOCTR) — sweph.c

Porting reference for the topocentric-observer machinery that `swe_calc`/`swe_calc_ut` use when
`SEFLG_TOPOCTR` is set. The Rust port currently **rejects** `TOPOCTR` outright
(`src/context.rs:178-181`, `unsupported = flags & CalcFlags::TOPOCTR`). This doc is the
prerequisite for lifting that restriction.

All line numbers refer to `sweph.c` unless prefixed otherwise. Line numbers are approximate
(±20 lines vs. the task brief, due to local source drift) — verified against the checked-out
`../swisseph/sweph.c` at doc-writing time.

---

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_set_topo` | sweph.c:7249–7267 | No — replaced by `EphemerisConfig.topographic: Option<TopoPosition>` (already exists, `src/context.rs:12-16, 28`) |
| `swi_force_app_pos_etc` | sweph.c:7269–7280 | No — invalidates C's per-planet/per-node position caches; meaningless in a stateless port (nothing is cached) |
| `swi_get_observer` | sweph.c:7282–7378 | **Yes — the core function.** Computes observer's geocentric offset (position+velocity, [6], AU/AU-day) in the J2000 mean-equatorial frame |
| TOPOCTR branch of `app_pos_etc_plan` | sweph.c:2524–2541, 2697–2708 | Yes — wires `xobs` into planet light-time/parallax/aberration pipeline |
| TOPOCTR branch of `app_pos_etc_sun` | sweph.c:3920–3940 | Yes — same wiring for the Sun (narrower: no light-time-at-t′ observer re-fetch) |
| TOPOCTR branch of `app_pos_etc_moon` | sweph.c:4117–4146, 4186–4200 | Yes — same wiring for the Moon |
| TOPOCTR branch of `swi_deflect_light` | sweph.c:3756–3760 | Read-only reference — shows deflection needs `xobs` too; **do not port the global-cache read**, see §7 |
| SPEED+TOPOCTR+aberration `use_speed3` forcing | sweph.c:407–410 | Yes — a `plaus_iflag`-level flag rule |

---

## 1. The core idea: `xobs` replaces the geocenter

Throughout `app_pos_etc_plan`/`_sun`/`_moon`, there is a local 6-vector `xobs` (position+velocity,
barycentric, AU/AU-day) that represents "the observer". Every place the code would otherwise use
`pedp->x` (Earth's barycentric position, i.e. the geocenter) it uses `xobs` instead:

- **Non-topocentric:** `xobs[i] = pedp->x[i]` for i=0..5 — the observer *is* the geocenter.
- **Topocentric:** `xobs[i] = <topo offset from swi_get_observer>[i] + pedp->x[i]` — the observer
  is the geocenter plus a small (Earth-radius-scale) position+velocity offset.

Everywhere downstream — light-time distance, parallax subtraction (`xx[i] -= xobs[i]`), and the
velocity argument to aberration (`swi_aberr_light(xx, xobs, iflag)`, which reads `xobs[3..6]`) —
is *unchanged code* that simply operates on this generalized `xobs`. This is the whole trick:
topocentric support is "compute a slightly-offset observer vector, then reuse the geocentric
pipeline unmodified."

**Porting implication:** in `src/calc.rs`, `apparent_planet`/`apparent_sun`/`apparent_moon`
currently hardcode `pos.earth_bary` (and `pos.earth_helio`) as the observer. Adding TOPOCTR means:
compute an observer offset vector (§3) once from `EphemerisConfig.topographic` + `jd`, add it to
`pos.earth_bary`/`pos.earth_helio`, and thread that generalized observer through the same
light-time/aberration calls that already exist — not writing new light-time/aberration logic.

---

## 2. `swe_set_topo` / global `topo_data` cache — STATELESS PORT NOTE

```c
/* sweph.c:7249 */
void CALL_CONV swe_set_topo(double geolon, double geolat, double geoalt)
{
  swi_init_swed_if_start();
  if (swed.geopos_is_set == TRUE
    && swed.topd.geolon == geolon && swed.topd.geolat == geolat && swed.topd.geoalt == geoalt) {
    return;
  }
  swed.topd.geolon = geolon;
  swed.topd.geolat = geolat;
  swed.topd.geoalt = geoalt;
  swed.geopos_is_set = TRUE;
  swed.topd.teval = 0;              /* force recompute of xobs next time */
  swi_force_app_pos_etc();          /* invalidate all planet/node position caches */
}
```

`struct topo_data` (sweph.h:758–763):
```c
struct topo_data {
  double geolon, geolat, geoalt;
  double teval;      /* TT epoch at which xobs[] was last computed */
  double tjd_ut;     /* corresponding UT epoch (teval - deltaT) */
  double xobs[6];    /* cached observer offset, J2000 mean equatorial, AU / AU-day */
};
```

**STATELESS PORT NOTE:** The Rust port already stores `(geolon, geolat, geoalt)` statically via
`EphemerisConfig.topographic: Option<TopoPosition>` (`src/context.rs:12-16,28`) — no `swe_set_topo`
call is needed or wanted; it's a constructor field, set once via `EphemerisConfig { topographic:
Some(TopoPosition { longitude, latitude, altitude }), .. }`.

What C caches that Rust must instead **recompute on every call**:
- `swed.topd.teval` / `swed.geopos_is_set == TRUE && teval == pedp->teval` check
  (sweph.c:2526-2527, 3397-3398/4121-4122, 3925-3926) — the "is the cached `xobs` still valid for
  this evaluation time" guard. In the stateless port this guard simply never applies: always call
  the `swi_get_observer` port fresh for the current `jd`.
- `swed.topd.xobs` (sweph.c:7371-7376, `do_save` branch) — the cached observer vector itself.
  Every one of the 3 call sites in §4/§5/§6 that read `swed.topd.xobs` on a cache hit must instead
  just call the observer function again (cheap: no file I/O, pure trig + 2 rotations).
- `swed.topd.tjd_ut` — not read anywhere outside `swi_get_observer` itself; safe to drop entirely.

There is **no** `do_save`/`NO_SAVE` distinction to port: C uses it only to decide whether to write
into `swed.topd`; a pure function has no such side channel, so the Rust observer function should
simply always return the computed vector.

---

## 3. `swi_get_observer` — the core computation

```c
/* sweph.c:7282 */
int swi_get_observer(double tjd, int32 iflag, AS_BOOL do_save, double *xobs, char *serr)
```
- `tjd` — TT (ET) julian day at which to evaluate the observer.
- `iflag` — flags; **critically, every call site in the TOPOCTR pipeline (§4/§5/§6) passes
  `iflag | SEFLG_NONUT`**, forcing the mean-frame branch (see §3.1 below — this is the single
  biggest simplification for the port).
- `do_save` — write result into `swed.topd`; **drop for the stateless port** (§2).
- `xobs` — output `double[6]`: position (AU) + velocity (AU/day), **J2000 mean equatorial frame**,
  i.e. *not yet barycentric* — the caller adds `pedp->x` (or equivalent) afterward.

### 3.1 Why `SEFLG_NONUT` is always forced, and what it skips

`swe_calc`'s three TOPOCTR call sites (§4 `app_pos_etc_plan`, §5 `app_pos_etc_moon`, §6
`app_pos_etc_sun`) all call `swi_get_observer(t, iflag | SEFLG_NONUT, ...)`. Reason (matches the
comment structure of the rest of the pipeline): the observer's position is built directly in the
**mean** equatorial frame of date; nutation is applied exactly once, later, together with the
celestial body itself, inside `app_pos_rest`'s shared `swi_nutate(xx, iflag, FALSE)` call
(sweph.c:2787-ish). Applying nutation twice (once to `xobs` alone, once to the combined result)
would double-count it. Forcing `NONUT` here means the *general* `swi_get_observer` algorithm
(which supports computing a fully-apparent, nutated observer vector, used by non-`swe_calc`
callers such as `swecl.c:5406`'s rise/set/eclipse code, out of scope for this doc) collapses at
these 3 call sites to a much simpler path — steps 4–5 below never execute:

```c
if (iflag & SEFLG_NONUT) {
  nut = 0;                       /* dpsi contribution to sidereal time = 0 */
} else {
  eps += nutlo[1];               /* true obliquity */
  nut = nutlo[0];
}
```
— with `nut == 0`, `swe_sidtime0` computes **mean** sidereal time, and `eps` stays the **mean**
obliquity of date (never gets `+= nutlo[1]`). So the whole "subtract nutation" block at
sweph.c:7356-7363 (guarded by `if (!(iflag & SEFLG_NONUT))`) is dead code on this path.

**Porting implication: the Rust port only needs to implement the NONUT-forced path below.** Do
not port the apparent/nutated branch (steps 4b/5b) unless a future rise/set-style topocentric
function needs it — it is unreachable from `calc()`.

### 3.2 Algorithm (as actually exercised via `calc()`, i.e. with NONUT forced)

Constants (already in `src/constants.rs:12-14`, values match exactly, cite AA = *Astronomical
Almanac*):
```c
#define EARTH_RADIUS      6378136.6                 /* meters, AA 2006 K6 */
#define EARTH_OBLATENESS  (1.0 / 298.25642)          /* flattening f, AA 2006 K6 */
#define EARTH_ROT_SPEED   (7.2921151467e-5 * 86400)  /* rad/day, Explanatory Supplement p.162 */
```

Steps (sweph.c:7300-7377):

1. **Delta-T / UT epoch** (sweph.c:7300-7301):
   ```c
   delt = swe_deltat_ex(tjd, iflag, serr);
   tjd_ut = tjd - delt;
   ```
   Note this literally passes the **TT** `tjd` into the deltaT function (which nominally wants
   UT) — C's own comment (sweph.c:7297-7299) calls this a deliberate, "extremely small" fudge, not
   a bug. Rust: `let tjd_ut = tjd - calc_deltat(tjd, config);` — do not iterate to convergence,
   replicate the one-shot approximation exactly.

2. **Mean obliquity of date** (sweph.c:7302-7307, with `NONUT` forced so the nutation half of the
   cache-check branch is skipped): `eps = swi_epsiln(tjd, iflag)` → Rust: `obliquity(tjd, flags,
   models).eps` (radians). `nut = 0` (sweph.c:7311-7312).

3. **Mean sidereal time, in degrees** (sweph.c:7319-7320):
   ```c
   sidt = swe_sidtime0(tjd_ut, eps * RADTODEG, nut /* == 0 */ * RADTODEG);
   sidt *= 15;   /* hours -> degrees */
   ```
   Rust: `sidereal_time::sidereal_time0(tjd_ut, eps.to_degrees(), 0.0, config) * 15.0` (matches the
   existing call convention at `context.rs:306`).

4. **Geodetic → geocentric flattening correction** (sweph.c:7333-7343). `f = EARTH_OBLATENESS`,
   `re = EARTH_RADIUS` (meters):
   ```c
   cosfi = cos(geolat * DEGTORAD);
   sinfi = sin(geolat * DEGTORAD);
   cc = 1 / sqrt(cosfi*cosfi + (1-f)*(1-f)*sinfi*sinfi);
   ss = (1-f)*(1-f) * cc;
   cosl = cos((geolon + sidt) * DEGTORAD);
   sinl = sin((geolon + sidt) * DEGTORAD);
   h = geoalt;   /* meters, above the reference ellipsoid, NOT mean sea level (see C comment
                    sweph.c:7322-7331: geoid vs. ellipsoid difference is <500m, negligible) */
   xobs[0] = (re*cc + h) * cosfi * cosl;
   xobs[1] = (re*cc + h) * cosfi * sinl;
   xobs[2] = (re*ss + h) * sinfi;
   ```
   `cc`/`ss` are the standard "radius of curvature in the prime vertical" reduction for an oblate
   spheroid; `(geolon + sidt)` folds the site's geographic longitude and Earth's rotation (mean
   sidereal time) into a single angle in the mean-equatorial frame — this is what makes `xobs` an
   **inertial-frame** (not Earth-fixed) vector already at this point, with z = Earth's mean
   rotation axis, x → mean equinox of date. Polar motion is explicitly neglected (comment
   sweph.c:7337-7339, "a few meters").
   FP-fidelity: note `(re*cc + h)` and `(re*ss + h)` — compute `re*cc` (resp. `re*ss`) first,
   *then* add `h`, matching C's left-to-right grouping; do not reorder to `h + re*cc`.

5. **Cartesian → polar, attach diurnal-rotation speed, polar → cartesian** (sweph.c:7348-7352):
   ```c
   swi_cartpol(xobs, xobs);       /* [x,y,z] -> [lon, lat, radius] (radians, radians, meters) */
   xobs[3] = EARTH_ROT_SPEED;     /* rad/day, diurnal rotation: only longitude changes */
   xobs[4] = xobs[5] = 0;         /* no latitude or radial motion */
   swi_polcart_sp(xobs, xobs);    /* -> [x,y,z,vx,vy,vz] meters, meters/day */
   ```
   Rust equivalents already exist: `cartesian_to_polar` (`src/math.rs:204`) and
   `polar_to_cartesian_with_speed` (`src/math.rs:267`) — note the position-only conversion uses
   `cartesian_to_polar`, but the position+speed conversion back uses the **with_speed** variant
   (asymmetric round-trip, matches C's `swi_cartpol` (position-only) then `swi_polcart_sp`
   (position+speed) exactly — do not "simplify" to a single `_sp` round trip, the intermediate
   step legitimately drops any existing velocity, there is none yet at that point anyway).

6. **Convert to AU / AU-day** (sweph.c:7354-7355):
   ```c
   for (i = 0; i <= 5; i++) xobs[i] /= AUNIT;
   ```
   `AUNIT` in meters, already in `src/constants.rs:7`.

7. **Nutation removal** (sweph.c:7357-7363) — **skipped**, since `NONUT` is always forced on this
   path (§3.1). `xobs` is already in the mean-of-date frame at this point.

8. **Precess mean-of-date → J2000** (sweph.c:7365-7367):
   ```c
   swi_precess(xobs, tjd, iflag, J_TO_J2000);          /* position */
   swi_precess_speed(xobs, tjd, iflag, J_TO_J2000);    /* velocity, adds precession-rate term */
   ```
   `J_TO_J2000` (sweph.h:256, `= 1`) means "date → J2000", i.e. exactly
   `PrecessionDirection::DateToJ2000` in `src/precession.rs`/`src/calc.rs:106` — same function,
   same direction enum variant, already implemented; reuse directly. No new precession code
   needed, just call it with `DateToJ2000` instead of the pipeline's usual `J2000ToDate`.

9. **Frame bias skipped** (sweph.c:7368-7369, `/* neglect frame bias (displacement of 45cm) */`)
   — deliberately not applied to `xobs`, even though `frame_bias`/`swi_bias` IS applied to the
   celestial body elsewhere in the pipeline. **Do not call `frame_bias` on the observer offset.**

10. **Cache write, dropped in stateless port** (sweph.c:7370-7376, see §2).

Return `xobs[6]`: position (AU) + velocity (AU/day) offset from the geocenter, in the J2000 mean
equatorial frame — ready for the caller to add to `pedp->x`/`pos.earth_bary` (§1).

---

## 4. TOPOCTR in `app_pos_etc_plan` (planets)

```c
static int app_pos_etc_plan(int ipli, int iplmoon, int32 iflag, char *serr)   /* sweph.c:2465 */
```

**Observer setup** (sweph.c:2524-2541):
```c
if (iflag & SEFLG_TOPOCTR) {
  if (swed.topd.teval != pedp->teval || swed.topd.teval == 0) {
    swi_get_observer(pedp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr);
  } else {
    /* cache hit: copy swed.topd.xobs */
  }
  for (i = 0; i <= 5; i++) xobs[i] = xobs[i] + pedp->x[i];   /* -> barycentric */
} else {
  for (i = 0; i <= 5; i++) xobs[i] = pedp->x[i];             /* observer == geocenter */
}
```
Rust: always call the §3 observer function at `jd` (the current planet evaluation epoch, i.e. the
Earth epoch, not yet light-time-retarded), then `xobs[i] = topo_offset[i] + pos.earth_bary[i]`.

**Light-time loop, aberration part** (sweph.c:2562-2596) — `xobs` (not `pedp->x`) is subtracted
throughout: `dx[i] -= (xobs[i] - xobs[i+3])` (sweph.c:2568, the "speed-influenced-by-dt-change"
pre-pass) and `dx[i] -= xobs[i]` (sweph.c:2586, the main light-time distance). This is identical in
shape to the existing geocentric code in `apparent_planet` (`src/calc.rs:925-947`), just with
`pos.earth_bary` swapped for the generalized `xobs`.

**Geocenter (parallax) subtraction** (sweph.c:2713-2726):
```c
if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR)) {
  for (i = 0; i <= 5; i++) xx[i] -= xobs[i];             /* parallax: geo -> topo shift */
  if (!(iflag & SEFLG_TRUEPOS) && (iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++) xx[i] -= xxsp[i-3];
}
```
This single line, `xx[i] -= xobs[i]`, IS the parallax correction — subtracting the observer's
(rather than the geocenter's) position from the planet's barycentric position directly produces a
topocentric apparent direction, no separate "parallax formula" exists in the C code.

**Observer at retarded time, for the SPEED2 finite-difference term** (sweph.c:2697-2708):
```c
if (iflag & SEFLG_SPEED) {
  if (iflag & SEFLG_TOPOCTR) {
    swi_get_observer(t /* = light-time-retarded epoch */, iflag | SEFLG_NONUT, NO_SAVE, xobs2, serr);
    for (i = 0; i <= 5; i++) xobs2[i] += xearth[i];   /* xearth = Earth bary @ t (retarded) */
  } else {
    for (i = 0; i <= 5; i++) xobs2[i] = xearth[i];
  }
}
```
then, in the aberration block (sweph.c:2748-2751):
```c
if (iflag & SEFLG_SPEED)
  for (i = 3; i <= 5; i++) xx[i] += xobs[i] - xobs2[i];
```
`xobs2` is the observer barycentric vector evaluated **at the retarded epoch** `t = jd - dt`
(not at the current epoch); `xobs[i] - xobs2[i]` for `i=3..5` is a difference of **velocity
components** (index 3..5 of the 6-vector is always velocity) between the current-epoch and
retarded-epoch observer state — i.e. how much the observer's own velocity (Earth's orbital motion
plus the topocentric diurnal-rotation term) changed over the light-time interval — and this delta
is added directly into the reported velocity `xx[3..5]`. Same shape as the existing
`apparent_planet` SPEED2 correction at `src/calc.rs:997-1002` (currently uses
`pos.earth_bary[3..6] - pos_ret.earth_bary[3..6]`) — topocentric port just swaps in the
topo-augmented `xobs`/`xobs2` instead of raw `earth_bary`/`earth_bary_ret`.

---

## 5. TOPOCTR in `app_pos_etc_moon`

```c
static int app_pos_etc_moon(int32 iflag, char *serr)   /* sweph.c:4087 */
```
Structurally identical pattern, slightly different variable layout because the Moon's own
position (`xxm`) is barycentric-relative (Earth+Moon), not Earth-bary + geocentric-Moon:

**Observer setup** (sweph.c:4120-4132):
```c
if (iflag & SEFLG_TOPOCTR) {
  if (swed.topd.teval != pdp->teval || swed.topd.teval == 0)
    swi_get_observer(pdp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr);
  else
    /* cache hit */;
  for (i = 0; i <= 5; i++) xxm[i] -= xobs[i];         /* Moon-relative-to-observer, for light-time dist */
  for (i = 0; i <= 5; i++) xobs[i] += pedp->x[i];     /* -> barycentric observer */
}
```
Note the extra `xxm[i] -= xobs[i]` line (sweph.c:4129-4130) that doesn't appear in the planet/sun
cases — this is because `xxm` here is used purely to get the light-time distance
(`dt = |xxm| * AUNIT/CLIGHT/86400`, sweph.c:4152), computed *before* `xobs` is promoted to
barycentric, i.e. it's `moon_bary - topo_offset`, an intermediate quantity, not reused elsewhere.

**Parallax subtraction** (sweph.c:4205-4206): same pattern, `xx[i] -= xobs[i]` for i=0..5.

**Retarded-time observer for SPEED2 term** (sweph.c:4186-4200, 4210-4222): identical shape to §4 —
`swi_get_observer(t, iflag|SEFLG_NONUT, NO_SAVE, xobs2, ...)`, `xobs2[i] += xe[i]` (Earth bary at
retarded `t`), then `xx[i+3] += xobs[i] - xobs2[i]` inside the aberration block.

**No `swi_deflect_light` call for the Moon** — confirmed by re-reading sweph.c:4087-4246 in full;
the Moon's apparent-position pipeline never calls the deflection function at all (too close to
Earth for the effect to matter). Adding TOPOCTR to the Moon path does **not** need any deflection
hookup — only the light-time/parallax/aberration wiring above.

---

## 6. TOPOCTR in `app_pos_etc_sun`

```c
static int app_pos_etc_sun(int32 iflag, char *serr)   /* sweph.c:3902 */
```

**Observer setup** (sweph.c:3920-3940) — same shape as §4/§5:
```c
if (iflag & SEFLG_TOPOCTR) {
  if (swed.topd.teval != pedp->teval || swed.topd.teval == 0)
    swi_get_observer(pedp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr);
  else
    /* cache hit */;
  for (i = 0; i <= 5; i++) xobs[i] = xobs[i] + pedp->x[i];
} else {
  for (i = 0; i <= 5; i++) xobs[i] = pedp->x[i];
}
```
`xx` (the Sun's position) is then built as `xobs - psdp->x` (heliocentric-Earth-equivalent,
sweph.c:3944-3950) rather than `xobs` being subtracted from a barycentric planet position — same
idea (observer replaces geocenter), different algebraic arrangement because the Sun's own
ephemeris position typically *is* effectively the negative of heliocentric Earth.

**Key asymmetry vs. §4/§5 — no SPEED2 finite-difference observer-at-retarded-time term.**
Re-reading `app_pos_etc_sun` in full (sweph.c:3902-4070): there is no second `swi_get_observer`
call, no `xobs2`, and the aberration block (sweph.c:4045-4048) has no
`xx[i+3] += xobs[i]-xobs2[i]` correction — unlike the plan/moon cases. Only the single `xobs`
(evaluated once, at the Sun's un-retarded evaluation epoch) is used. **Do not add an xobs2 term
when porting the Sun's TOPOCTR path — that would diverge from C.**

Also note: **no deflection call for the Sun either** (a body can't gravitationally deflect its own
light) — same as the Moon, so Sun TOPOCTR support is pure light-time + parallax + aberration, no
deflection hookup.

---

## 7. `swi_deflect_light`'s topocentric dependency — ties to the documented deflection-speed tolerance

`swi_deflect_light` (sweph.c:3743, called only from the **planet** path, §4) reads the observer
offset directly from **global state**, not from a parameter:

```c
/* sweph.c:3756-3760 */
struct plan_data *pedp = &swed.pldat[SEI_EARTH];
for (i = 0; i <= 5; i++) xearth[i] = pedp->x[i];
if (iflag & SEFLG_TOPOCTR)
  for (i = 0; i <= 5; i++) xearth[i] += swed.topd.xobs[i];   /* <-- global cache read */
```
`swed.topd.xobs` here is whatever was last cached by `swi_get_observer(..., DO_SAVE, ...)` at the
**current, un-retarded** evaluation epoch (the same `xobs` computed in §4's observer setup) — so
functionally this is correct (it *is* the same vector `app_pos_etc_plan` computed), it's just
fetched via a second, independent global read rather than being passed as a parameter.

The consequence that matters for the port: `swi_deflect_light`'s **speed** branch
(sweph.c:3819-3887) perturbs time by `dtsp = -DEFL_SPEED_INTV` and reconstructs `e[i] = xearth[i] -
dtsp*xearth[i+3]` — a **linear extrapolation** using `xearth`'s *velocity* (which includes the
observer's cached diurnal-rotation velocity from `swed.topd.xobs[3..6]`), not a fresh
`swi_get_observer` call at `t + dtsp`. This is precisely the mechanism behind this project's
already-documented, already-accepted tolerance:

> **Deflection speed (< 0.06 milliarcseconds)** — "The C `swi_deflect_light` reads the Sun's
> position from a global cache... populated earlier in the same `swe_calc` call. The stateless
> Rust version constructs the deflection geometry from explicitly-passed parameters." (project
> `CLAUDE.md`, `<stateless_tolerance>` §1)

**Porting implication:** when wiring TOPOCTR into `deflect_light` (`src/corrections.rs:214`), pass
`pos.earth_helio + topo_offset` as the `earth_helio` parameter (position **and** velocity) exactly
once, at the current (un-retarded) epoch — do **not** try to re-fetch the observer at the
perturbed `t ± DEFL_SPEED_INTV` epoch to "improve" on C; that would not match C's behavior (C
itself uses the linear extrapolation, not a re-fetch) and would fight the already-accepted
tolerance rather than reproduce it. Four prior debugging sessions already chased this exact
discrepancy to no benefit (see `<stateless_tolerance>` in project CLAUDE.md) — accept the 1e-7°
speed tolerance for TOPOCTR+SPEED the same way it's already accepted for non-topocentric deflection
speed.

---

## 8. SPEED + TOPOCTR + aberration interaction (`use_speed3` forcing)

```c
/* sweph.c:402-410 */
/* high precision speed prevails fast speed */
if ((iflag & SEFLG_SPEED3) && (iflag & SEFLG_SPEED))
  iflag = iflag & ~SEFLG_SPEED3;
if (iflag & SEFLG_SPEED3)
  use_speed3 = TRUE;
/* topocentric with SEFLG_SPEED is not good if aberration is included.
 * in such cases we calculate speed from three positions */
if ((iflag & SEFLG_SPEED) && (iflag & SEFLG_TOPOCTR) && !(iflag & SEFLG_NOABERR))
  use_speed3 = TRUE;
```
Plain-English: normal analytic `SEFLG_SPEED` (the finite-difference-free, direct-from-ephemeris-
derivative approach) is **not trusted** when `TOPOCTR` + aberration are both active — C silently
switches to the 3-point numerical differentiation (`SPEED3`: evaluate positions at `t-dt`, `t`,
`t+dt` and finite-difference) instead, regardless of what the caller asked for. `NOABERR` (or
`TRUEPOS`, which implies it) opts back out.

**Porting implication:** this is a flag-normalization rule, analogous to the existing rules in
`plaus_iflag` (`src/calc.rs:28-66`, e.g. `TOPOCTR` implies removing `HELCTR`/`BARYCTR`,
`SPEED`+`SPEED3` mutual exclusion at line 54-56). Add, in the same function, after the existing
`SPEED`/`SPEED3` exclusion:
```rust
if flags.contains(CalcFlags::SPEED)
    && flags.contains(CalcFlags::TOPOCTR)
    && !flags.contains(CalcFlags::NOABERR)
{
    flags.insert(CalcFlags::SPEED3);
    // Note: unlike the earlier SPEED-prevails-over-SPEED3 rule (line 54-56, which removes
    // SPEED3 when SPEED is set), this rule runs the opposite direction: SPEED3 wins here.
    // C never clears SEFLG_SPEED itself (`use_speed3` is a separate local bool, not a flag
    // mutation) — mirror that: do NOT remove CalcFlags::SPEED, only add SPEED3, and make the
    // SPEED3 dispatch check happen first (context.rs:190) so SPEED3 takes priority when both
    // are set.
}
```
Check against `context.rs:190-192` (`if flags.contains(CalcFlags::SPEED3) { return
self.calc_speed3(...) }`) and `calc_speed3` at `context.rs:608` — the existing SPEED3 dispatch
should already do the right thing once `SPEED3` is set here; no new numerical-differentiation code
should be needed, just correct flag routing before dispatch.

---

## 9. Rust integration points (for quick orientation, not exhaustive)

- `src/context.rs:12-16` — `TopoPosition { longitude, latitude, altitude }`, already defined.
- `src/context.rs:28` — `EphemerisConfig.topographic: Option<TopoPosition>`, already threaded.
- `src/context.rs:178-181` — the `TOPOCTR` rejection to remove once ported.
- `src/constants.rs:7,12-14` — `AUNIT`, `EARTH_RADIUS`, `EARTH_OBLATENESS`, `EARTH_ROT_SPEED`,
  already present with values matching sweph.h exactly.
- `src/math.rs:204,236,267` — `cartesian_to_polar`, `cartesian_to_polar_with_speed`,
  `polar_to_cartesian_with_speed` — reuse for §3 steps 5/6, no new conversion code needed.
  (`rotate_x_sincos` at `math.rs:196` = `swi_coortrf2`, needed only if a future non-`calc()`
  topocentric caller requires the nutation-removal branch skipped per §3.1.)
- `src/sidereal_time.rs:220` — `sidereal_time0(tjd_ut, eps, nut, config)`, reuse for §3 step 3.
- `src/deltat/mod.rs` — `calc_deltat(jd, config)`, reuse for §3 step 1.
- `src/obliquity.rs:30` — `obliquity(jd, flags, models)`, reuse for §3 step 2.
- `src/precession.rs` — `precess`/`PrecessionDirection::DateToJ2000`, reuse for §3 step 8 (this is
  the `J_TO_J2000` direction).
- `src/calc.rs:106` — `precess_speed(..., PrecessionDirection::DateToJ2000)`, reuse for §3 step 8.
- `src/calc.rs:897,1023,1088` — `apparent_planet`, `apparent_sun`, `apparent_moon` — the 3
  functions needing the `xobs` plumbing from §4/§5/§6.
- `src/corrections.rs:121,214` — `aberr_light`, `deflect_light` — existing signatures already take
  the "observer/earth" vector as an explicit parameter (`earth_vel: &[f64;3]`,
  `earth_helio: &[f64;6]`); topocentric support is purely a matter of what gets passed in, no
  signature changes needed.
- `src/calc.rs:28` — `plaus_iflag`, add the `use_speed3` rule from §8.
- `src/context.rs:608` — `calc_speed3`, verify it already handles the flag-forced case correctly
  once §8 is wired in.
