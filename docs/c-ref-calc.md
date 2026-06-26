# C Reference: Calculation Pipeline — sweph.c

Porting reference for the main calculation pipeline (`swe_calc` → `swecalc` → `app_pos_*`).
Read this instead of the C source.

## Function Map

| C function | Location | Notes |
|---|---|---|
| `plaus_iflag` | sweph.c:6066–6148 | Flag validation/normalization |
| `swe_calc` | sweph.c:309–563 | Public entry: cache check, SPEED3 |
| `swecalc` | sweph.c:587–1156 | Body dispatch table |
| `main_planet` | sweph.c:1562–1673 | Ephemeris cascade for main planets |
| `main_planet_bary` | sweph.c:1697–1747 | Barycentric-only variant |
| `sweplan` | sweph.c:1820–1968 | SWISSEPH: planet + Earth + Moon |
| `embofs` | sweph.c:5062–5067 | EMB → Earth conversion |
| `app_pos_etc_plan` | sweph.c:2465–2775 | Planet: light-time, bias, precession |
| `app_pos_etc_sun` | sweph.c:3902–4070 | Sun-specific path |
| `app_pos_etc_moon` | sweph.c:4087–4246 | Moon-specific path |
| `app_pos_rest` | sweph.c:2777–2859 | Shared: nutation, ecl, polar, units |
| `swi_aberr_light` | sweph.c:3699–3736 | Lorentz aberration |
| `swi_deflect_light` | sweph.c:3743–3891 | GR light deflection |
| `meff` | sweph.c:5967–5981 | Effective solar mass factor |
| `swi_mean_node` | swemmoon.c:1493–1534 | Mean lunar node |
| `swi_mean_apog` | swemmoon.c:1564–1624 | Mean lunar apogee (Dark Moon) |
| `mean_elements` | swemmoon.c:1763–1818 | Moshier mean orbital elements |
| `denormalize_positions` | sweph.c:5983–5997 | Unwrap ±180° discontinuities |
| `calc_speed` | sweph.c:5999–6011 | 3-point central-difference speed |

---

## 1. `plaus_iflag()` — Flag Validation/Normalization

**sweph.c:6066–6148**

Called at entry of `swecalc()` (and several other public functions). Normalizes `iflag` by resolving illegal combinations. Returns corrected `iflag`.

### Mutual-exclusion rules (in order applied)

```
1. JPLHOR and JPLHOR_APPROX: if both set, clear JPLHOR_APPROX.
2. TOPOCTR: clear HELCTR and BARYCTR.
3. BARYCTR: clear HELCTR.
4. HELCTR: clear BARYCTR.
5. HELCTR or BARYCTR: force NOABERR | NOGDEFL (aberration and deflection
   are meaningless for helio/bary positions).
6. J2000: force NONUT.
7. SIDEREAL: force NONUT; clear JPLHOR and JPLHOR_APPROX.
8. TRUEPOS: force NOGDEFL | NOABERR.
```

### Default ephemeris selection

```
epheflag = 0
if MOSEPH in iflag  → epheflag = SEFLG_MOSEPH
if SWIEPH in iflag  → epheflag = SEFLG_SWIEPH
if JPLEPH in iflag  → epheflag = SEFLG_JPLEPH
if epheflag == 0    → epheflag = SEFLG_DEFAULTEPH
iflag = (iflag & ~SEFLG_EPHMASK) | epheflag
```

### JPLHOR flag cleanup

```
if !(epheflag & SEFLG_JPLEPH):
    clear JPLHOR and JPLHOR_APPROX

bodies that never use JPLHOR (clear regardless of ephemeris):
    SE_OSCU_APOG, SE_TRUE_NODE, SE_MEAN_APOG, SE_MEAN_NODE,
    SE_INTP_APOG, SE_INTP_PERG, SE_FICT_OFFSET..SE_FICT_MAX

if JPLHOR set but eop_dpsi_loaded ≤ 0:
    → downgrade to JPLHOR_APPROX, set warning in serr

if JPLHOR set (and eop data loaded):
    force SEFLG_ICRS

if JPLHOR_APPROX and jplhora_model == SEMOD_JPLHORA_2:
    force SEFLG_ICRS
```

---

## 2. `swe_calc()` — Public Entry Point

**sweph.c:309–563**

### Pre-processing

```
if ipl == SE_AST_OFFSET + 134340:
    ipl = SE_PLUTO           // Pluto asteroid redirect

if SPEED3 and SPEED both set:
    clear SPEED3             // high-precision speed wins

if SPEED and TOPOCTR and !NOABERR:
    use_speed3 = TRUE        // topocentric with aberration needs 3-call

if XYZ and RADIANS both set:
    clear RADIANS

// Planetary moon / center-of-body number extraction
if SEFLG_CENTER_BODY and ipl <= SE_PLUTO:
    iplmoon = ipl * 100 + 9099
if ipl >= SE_PLMOON_OFFSET and ipl < SE_AST_OFFSET:
    iplmoon = ipl
    ipl = (ipl - 9000) / 100
    iflag |= SEFLG_CENTER_BODY
if SEFLG_CENTER_BODY and ipl <= SE_MARS and iplmoon%100 == 99:
    iplmoon = 0; clear CENTER_BODY
```

### Cache check (save area `sd = swed.savedat[ipl]`)

```
if sd->tsave == tjd && tjd != 0 && ipl == sd->ipl && iplmoon == 0:
    if (sd->iflgsave & ~SEFLG_COORDSYS) == (iflag & ~SEFLG_COORDSYS):
        goto end_swe_calc     // cache hit
```

Coordinate-system flags (EQUATORIAL, XYZ, RADIANS) are excluded from the cache key — the save area stores all coordinate variants.

### SPEED3 mode (3-call central difference)

Body-specific dt values:

| Body | dt (days) | Constant |
|---|---|---|
| SE_MOON | 0.00005 | `MOON_SPEED_INTV` |
| SE_OSCU_APOG, SE_TRUE_NODE | 0.1 | `NODE_CALC_INTV_MOSH` |
| all others | 0.0001 | `PLAN_SPEED_INTV` |

Sequence:
```
swecalc(tjd - dt, ..., x0)
swecalc(tjd + dt, ..., x2)
swecalc(tjd,      ..., sd->xsaves)
denormalize_positions(x0, sd->xsaves, x2)
calc_speed(x0, sd->xsaves, x2, dt)
```

### Output extraction (24-element save area layout)

The save area `sd->xsaves[24]` uses the layout established by `app_pos_rest()`:

| Offset | Content |
|---|---|
| [0..5] | ecliptic polar + speed (lon, lat, dist + rates) |
| [6..11] | ecliptic cartesian + speed |
| [12..17] | equatorial polar + speed (RA, Dec, dist + rates) |
| [18..23] | equatorial cartesian + speed |

Selection:
```
if EQUATORIAL: xs = xsaves + 12
else:          xs = xsaves + 0      (ecliptic)
if XYZ:        xs += 6              (cartesian offset within chosen frame)
if SE_ECL_NUT: copy 4 elements (true_ecl, mean_ecl, dpsi, deps)
else:          copy 3 elements (lon/RA, lat/Dec, dist)
if SPEED3 or SPEED: copy [3..5]
if RADIANS: convert [0..1] to radians; if SE_ECL_NUT convert [0..3]
```

Return value: `iflag` with correct ephemeris bit, with coord-system bits taken from original `iflgsave`.

---

## 3. `swecalc()` — Body Dispatch

**sweph.c:587–1156**

Calls `plaus_iflag()` first. Then dispatches:

### SE_ECL_NUT (lines 656–664)

Does not call any ephemeris function. Reads directly from `swed.oec` and `swed.nut`:
```
x[0] = swed.oec.eps + swed.nut.nutlo[1]    // true ecliptic obliquity (rad)
x[1] = swed.oec.eps                         // mean ecliptic obliquity (rad)
x[2] = swed.nut.nutlo[0]                    // nutation in longitude (rad)
x[3] = swed.nut.nutlo[1]                    // nutation in obliquity (rad)
multiply all by RADTODEG
return iflag
```
These are already computed by `swi_check_ecliptic()` and `swi_check_nutation()` called at the top of `swecalc()`.

### SE_MOON (lines 668–728)

```
ipli = SEI_MOON
switch(epheflag):
    JPLEPH:  jplplan() → if NOT_AVAILABLE try SWIEPH
                       → if BEYOND_EPH_LIMITS try MOSEPH
    SWIEPH:  sweplan(SEI_MOON, SEI_FILE_MOON) → if NOT_AVAILABLE try MOSEPH
    MOSEPH:  swi_moshmoon() + swi_moshplan(SEI_EARTH)
then: app_pos_etc_moon(iflag)
```

### SE_SUN with SEFLG_BARYCTR (lines 733–815)

Special: barycentric sun handled separately from heliocentric sun.
```
JPLEPH:  swi_pleph(J_SUN, J_SBARY) directly
SWIEPH:  sweplan(SEI_EARTH, SEI_FILE_PLANET) — bary sun is a by-product
then: app_pos_etc_sbar(iflag)
Note: pedp->xflgs = -1 to force recomputation if geocentric earth follows
```

### SE_SUN, SE_MERCURY..SE_PLUTO, SE_EARTH (lines 819–854)

```
Guard: heliocentric Sun → return zero vector
Guard: geocentric Earth → return zero vector
ipli = pnoext2int[ipl]
retc = main_planet(tjd, ipli, iplmoon, epheflag, iflag)
iflag = pdp->xflgs
```

`pnoext2int` mapping (line 182):
```
SE_SUN(0)     → SEI_SUN(=SEI_EARTH=0)
SE_MOON(1)    → SEI_MOON(1)
SE_MERCURY(2) → SEI_MERCURY(2)
SE_VENUS(3)   → SEI_VENUS(3)
SE_MARS(4)    → SEI_MARS(4)
SE_JUPITER(5) → SEI_JUPITER(5)
SE_SATURN(6)  → SEI_SATURN(6)
SE_URANUS(7)  → SEI_URANUS(7)
SE_NEPTUNE(8) → SEI_NEPTUNE(8)
SE_PLUTO(9)   → SEI_PLUTO(9)
SE_EARTH(14)  → SEI_EARTH(0)
SE_CHIRON(15) → SEI_CHIRON(15)
SE_PHOLUS(16) → SEI_PHOLUS(16)
SE_CERES(17)  → SEI_CERES(17)
SE_PALLAS(18) → SEI_PALLAS(18)
SE_JUNO(19)   → SEI_JUNO(19)
SE_VESTA(20)  → SEI_VESTA(20)
```

### SE_MEAN_NODE (lines 859–893)

```
ndp = swed.nddat[SEI_MEAN_NODE]
xp2 = ndp->x
swi_mean_node(tjd, xp2)
swi_mean_node(tjd - MEAN_NODE_SPEED_INTV, xp2+3)   // 0.001 day back
xp2[3] = swe_difrad2n(xp2[0], xp2[3]) / MEAN_NODE_SPEED_INTV  // angular speed
xp2[4] = xp2[5] = 0
ndp->teval = tjd; ndp->xflgs = -1
app_pos_etc_mean(SEI_MEAN_NODE, iflag)

// force lat=0: ndp->xreturn[1,4,5,8,11] = 0.0
// (unless SIDEREAL or J2000)
```

### SE_MEAN_APOG (lines 898–927)

```
ndp = swed.nddat[SEI_MEAN_APOG]
swi_mean_apog(tjd, xp2)
swi_mean_apog(tjd - MEAN_NODE_SPEED_INTV, xp2+3)   // 0.001 day back
for i in 0..1:
    xp2[3+i] = swe_difrad2n(xp2[i], xp2[3+i]) / MEAN_NODE_SPEED_INTV
xp2[5] = 0
app_pos_etc_mean(SEI_MEAN_APOG, iflag)
ndp->xreturn[5] = 0.0   // force radial speed = 0
```

`MEAN_NODE_SPEED_INTV = 0.001` days.

---

## 4. `main_planet()` — Ephemeris Cascade for Main Planets

**sweph.c:1562–1673**

Selects ephemeris and handles fallback chain. After position is computed, calls `app_pos_etc_sun()` (for SEI_SUN) or `app_pos_etc_plan()` (for all others).

```
if SEFLG_CENTER_BODY and ipli in MARS..PLUTO:
    sweph(iplmoon, SEI_FILE_ANY_AST)   // center-of-body offset

switch(epheflag):
    JPLEPH:
        jplplan(ipli)
        if NOT_AVAILABLE → switch to SWIEPH, goto sweph_planet
        if BEYOND_EPH_LIMITS and in MOSHPLEPH range → switch to MOSEPH, goto moshier_planet
        if BEYOND_EPH_LIMITS otherwise → ERR
        app_pos_etc_sun or app_pos_etc_plan
        if NOT_AVAILABLE → try SWIEPH
        if BEYOND_EPH_LIMITS → try MOSEPH

    SWIEPH (sweph_planet):
        sweplan(ipli, SEI_FILE_PLANET)
        if NOT_AVAILABLE and in MOSHPLEPH range → switch to MOSEPH
        app_pos_etc_sun or app_pos_etc_plan
        if NOT_AVAILABLE → try MOSEPH

    MOSEPH (moshier_planet):
        swi_moshplan(ipli)
        app_pos_etc_sun or app_pos_etc_plan
```

Range constants: `MOSHPLEPH_START`, `MOSHPLEPH_END` (checked in header).

---

## 5. `main_planet_bary()` — Raw Barycentric

**sweph.c:1697–1747**

Returns raw barycentric positions without `app_pos_*` corrections. Used by `app_pos_etc_plan()` internally for light-time iteration.

```
switch(epheflag):
    JPLEPH:
        jplplan(ipli, do_save, xp, xe, xs)
        if NOT_AVAILABLE → try SWIEPH

    SWIEPH (sweph_planet):
        sweplan(ipli, SEI_FILE_PLANET, do_save, xp, xe, xs, xm)
        if NOT_AVAILABLE and in MOSHPLEPH → try MOSEPH

    MOSEPH (moshier_planet):
        swi_moshplan(ipli, do_save, xp, xe)
        xs[0..5] = 0    // no barycentric sun with Moshier
```

---

## 6. `sweplan()` — SWISSEPH Planet + Earth + Moon Orchestration

**sweph.c:1820–1968**

Computes barycentric positions for planet, Earth, Moon, and Sun (as needed), with caching per sub-body. Returns positions in barycentric equatorial J2000.

### Which sub-bodies to compute

```
do_sunbary = do_save OR ipli==SEI_SUNBARY OR (pdp->iflg & SEI_FLG_HELIO)
             OR xpsret!=NULL OR (iflag & SEFLG_HELCTR)
do_earth   = do_save OR ipli==SEI_EARTH OR xperet!=NULL
do_moon    = do_save OR ipli==SEI_MOON OR ipli==SEI_EARTH OR xperet!=NULL OR xpmret!=NULL
if ipli==SEI_MOON: do_earth=TRUE, do_sunbary=TRUE
```

### Cache check pattern (each sub-body)

```
if tjd == pdp->teval AND pdp->iephe == SEFLG_SWIEPH AND (!speedf2 OR speedf1):
    use cached value
else:
    call sweph()
    update teval, xflgs
```

`speedf1` = cached speed flag, `speedf2` = requested speed flag.

### Moon via Moshier fallback

```
if moon file not found (fidat[SEI_FILE_MOON].fptr == NULL):
    use swi_moshmoon() instead
```

### EMB → Earth conversion (`embofs`)

After loading EMB (Earth-Moon Barycenter), applies:
```c
// sweph.c:5062
for i in 0..2:
    xemb[i] -= xmoon[i] / (EARTH_MOON_MRAT + 1.0)
```

Where `EARTH_MOON_MRAT = 81.30056907419062` (DE431 value, sweph.h:267).

Applied twice if saving (once for position, once for speed):
```
embofs(xpe, xpm)            // position
if xpe==pebdp->x OR SPEED:
    embofs(xpe+3, xpm+3)    // speed
```

### Heliocentric → barycentric planet conversion

```
if pdp->iflg & SEI_FLG_HELIO:
    xp[0..2] += xps[0..2]
    if do_save OR SPEED: xp[3..5] += xps[3..5]
```

### Special planet routing (ipli)

```
SEI_MOON  → xp = xpm
SEI_EARTH → xp = xpe
SEI_SUN   → xp = xps (barycentric sun)
other     → call sweph(ipli, ifno, ...)
```

---

## 7. `app_pos_etc_plan()` — Planet Position Corrections

**sweph.c:2465–2775**

Transforms barycentric J2000 position to apparent position with corrections. Updates `pdp->xreturn[24]`. Called after barycentric position is in `pdp->x`.

### Cache check

```
flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ
flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ
if flg1 == flg2: return OK    // coord-system flags excluded
```

### Observer position

```
if TOPOCTR:
    xobs = swi_get_observer() + pedp->x   // topocentric + barycentric earth
else:
    xobs = pedp->x                         // geocenter
```

### Heliocentric mode

```
if HELCTR and iephe is JPL or SWISSEPH:
    xx[0..5] -= swed.pldat[SEI_SUNBARY].x[0..5]
```

### Light-time iteration

Number of iterations:
- JPL or SWISSEPH: `niter = 1` (one Newton step)
- Moshier or osculating elements: `niter = 0` (linear approximation only — subtraction of `dt * v`)

For speed correction (change of light-delay over time):
```
// Compute apparent position at t-1 day before main loop:
xxsp[0..2] = xx[0..2] - xx[3..5]    // rough position at t-1
for j in 0..niter:
    dx[0..2] = xxsp - xobs(t-1)
    dt = |dx| * AUNIT / CLIGHT / 86400.0
    xxsp[0..2] = xxsv - dt * xx0[3..5]
// then: xxsp = true(t-1) - apparent(t-1)  [correction term]
```

Main light-time loop:
```
for j in 0..niter:
    dx[0..2] = xx - xobs
    dt = |dx| * AUNIT / CLIGHT / 86400.0
    t = pdp->teval - dt              // retarded time
    dtsave_for_defl = dt             // saved for deflection
    xx[0..2] = xx0 - dt * xx0[3..5] // rough apparent position

// For accurate position at t': call ephemeris for body at t'
switch(epheflag):
    JPLEPH:  swi_pleph(t, pnoint2jpl[ipli], J_SBARY, xx)
    SWIEPH:  sweplan(t, ipli, ...) or sweph(t, ipli, ...)
    MOSEPH:  if SPEED: swi_moshplan(t, ...) — only speed taken from this
```

After light-time:
```
if SPEED:
    xx[3..5] -= xxsp[0..2]   // add dt-change correction
```

### Geocentric conversion

```
if !(HELCTR or BARYCTR):
    xx[0..5] -= xobs[0..5]
    if SPEED and !TRUEPOS:
        xx[3..5] -= xxsp[0..2]   // change-of-dt speed correction
if !SPEED: xx[3..5] = 0
```

### Deflection and aberration

```
if !TRUEPOS and !NOGDEFL:
    swi_deflect_light(xx, dtsave_for_defl, iflag)

if !TRUEPOS and !NOABERR:
    swi_aberr_light(xx, xobs, iflag)
    if SPEED:
        xx[3..5] += xobs[0..5 shifted by 3] - xobs2[0..5 shifted by 3]

if !SPEED: xx[3..5] = 0
```

### ICRS → J2000 frame bias

```
if !SEFLG_ICRS and swi_get_denum(ipli, epheflag) >= 403:
    swi_bias(xx, t, iflag, FALSE)    // DE403+: apply frame bias
```

### Precession

```
// Save J2000 for sidereal positions:
xxsv[0..5] = xx[0..5]

if !J2000:
    swi_precess(xx, pdp->teval, iflag, J2000_TO_J)
    if SPEED: swi_precess_speed(xx, pdp->teval, iflag, J2000_TO_J)
    oe = &swed.oec         // obliquity of date
else:
    oe = &swed.oec2000     // J2000 obliquity

return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr)
```

---

## 8. `app_pos_etc_sun()` — Sun-Specific Path

**sweph.c:3902–4070**

Geocentric sun = **negative of Earth's heliocentric position**. The "planet" here is the Earth.

### True heliocentric Earth position

```
if iephe==MOSEPH or BARYCTR:
    xx = xobs              // Moshier: heliocentric
else:
    xx = xobs - psdp->x   // JPL/SWISSEPH: earth - bary_sun = helio_earth
```

### Light-time

For geocentric sun (default):
```
// iterates: at each step, get new sun position at t' = teval - dt
switch(iephe):
    JPLEPH: swi_pleph(t, J_SUN, J_SBARY, xsun)
    SWIEPH: sweph(t, SEI_SUNBARY, SEI_FILE_PLANET, ..., xsun)
    MOSEPH: no new position (subtraction already done)
xx = xearth(t') - xsun(t')    // helio_earth at retarded time
```

For heliocentric or barycentric Earth:
```
niter = 1; iterates earth position backward
```

### Geocentric conversion

```
if !(HELCTR or BARYCTR):
    xx = -xx                // flip sign: helio_earth → geo_sun
```

No gravitational deflection for the Sun (not called).

### Aberration, bias, precession

Same pattern as `app_pos_etc_plan()`:
```
if !TRUEPOS and !NOABERR: swi_aberr_light(xx, xobs, iflag)
if !ICRS and denum>=403: swi_bias(xx, t, iflag, FALSE)
save xxsv = xx (J2000 for sidereal)
if !J2000: swi_precess(xx, pedp->teval, ..., J2000_TO_J)
return app_pos_rest(pedp, iflag, xx, xxsv, oe, serr)
```

---

## 9. `app_pos_etc_moon()` — Moon-Specific Path

**sweph.c:4087–4246**

The Moon position in `pdp->x` is **geocentric** (relative to Earth center), not barycentric.

### Convert to barycentric

```
xx[0..5] = pdp->x[0..5] + pedp->x[0..5]   // moon_geo + earth_bary = moon_bary
xxm[0..5] = pdp->x[0..5]                   // save geocentric for distance calc
```

### Observer position

```
if TOPOCTR:
    xxm -= xobs_topo    // topocentric moon vector
    xobs = xobs_topo + pedp->x
elif BARYCTR:
    xobs = 0
    xxm += pedp->x
elif HELCTR:
    xobs = psdp->x
    xxm += pedp->x - psdp->x
else:
    xobs = pedp->x      // geocenter
```

### Light-time (single step, no iteration)

```
dt = |xxm| * AUNIT / CLIGHT / 86400.0
t = pdp->teval - dt
switch(iephe):
    JPLEPH:  swi_pleph(t, J_MOON, J_EARTH, xx_new)
             swi_pleph(t, J_EARTH, J_SBARY, xe)
             if HELCTR: swi_pleph(t, J_SUN, J_SBARY, xs)
             xx = xx_new + xe             // bary moon at t'
    SWIEPH:  sweplan(t, SEI_MOON, SEI_FILE_MOON, ..., xx, xe, xs, NULL)
             xx += xe
    MOSEPH:  // linear only — "results in an error of a milliarcsec in speed"
             xx[0..2] -= dt * xx[3..5]
             xe[0..2] = pedp->x[0..2] - dt * pedp->x[3..5]
             xe[3..5] = pedp->x[3..5]
             xs = 0
```

No gravitational deflection (Moon too close). The comment notes that Moon's speed comes from central differences already computed by Moshier.

### Aberration

```
if !TRUEPOS and !NOABERR:
    swi_aberr_light(xx - xobs, xobs, iflag)
    if SPEED:
        xx[3..5] += xobs - xobs2    // earth velocity correction
```

### Bias and precession (same as plan)

```
if !ICRS and swi_get_denum(SEI_MOON) >= 403:
    swi_bias(xx, t, iflag, FALSE)
save xxsv = xx
if !J2000: swi_precess(xx, pdp->teval, J2000_TO_J)
return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr)
```

---

## 10. `app_pos_rest()` — Shared Final Pipeline

**sweph.c:2777–2859**

Applies nutation, ecliptic transformation, sidereal adjustment, polar/cartesian conversion, and unit conversion. Fills `pdp->xreturn[24]`.

### xreturn layout (24 doubles)

| Index | Content |
|---|---|
| [0..2] | ecliptic polar: lon (°), lat (°), dist (AU) |
| [3..5] | ecliptic polar speed: dlon/day (°), dlat/day (°), ddist/day (AU) |
| [6..8] | ecliptic cartesian XYZ |
| [9..11] | ecliptic cartesian velocity |
| [12..14] | equatorial polar: RA (°), Dec (°), dist (AU) |
| [15..17] | equatorial polar speed: dRA/day (°), dDec/day (°), ddist/day |
| [18..20] | equatorial cartesian XYZ (also input to this function via `xx`) |
| [21..23] | equatorial cartesian velocity |

### Step 1: Nutation (equatorial cartesian)

```
if !NONUT:
    swi_nutate(xx, iflag, FALSE)
// xx is now equatorial of date (J2000 → equator of date → nutated)
pdp->xreturn[18..23] = xx[0..5]
```

### Step 2: Ecliptic transformation

```
swi_coortrf2(xx, xx, oe->seps, oe->ceps)      // equatorial → ecliptic
if SPEED: swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps)

if !NONUT:
    swi_coortrf2(xx, xx, swed.nut.snut, swed.nut.cnut)
    if SPEED: swi_coortrf2(xx+3, xx+3, ...)

pdp->xreturn[6..11] = xx[0..5]   // ecliptic cartesian
```

### Step 3: Sidereal positions (if SEFLG_SIDEREAL)

Three variants, depending on `swed.sidd.sid_mode`:
```
SE_SIDBIT_ECL_T0:   swi_trop_ra2sid_lon(x2000, xreturn+6, xreturn+18)
SE_SIDBIT_SSY_PLANE: swi_trop_ra2sid_lon_sosy(x2000, xreturn+6)
else (traditional):
    swi_cartpol_sp(xreturn+6, xreturn)
    daya = swi_get_ayanamsa_with_speed(teval, ...)
    xreturn[0] -= daya[0] * DEGTORAD
    xreturn[3] -= daya[1] * DEGTORAD
    swi_polcart_sp(xreturn, xreturn+6)
```

### Step 4: Polar coordinates

```
swi_cartpol_sp(xreturn+18, xreturn+12)   // equatorial: cartesian → polar
swi_cartpol_sp(xreturn+6,  xreturn+0)    // ecliptic:   cartesian → polar
```

### Step 5: Radians → degrees

```
for i in 0..1:
    xreturn[i]    *= RADTODEG    // ecliptic lon, lat
    xreturn[i+3]  *= RADTODEG    // ecliptic lon speed, lat speed
    xreturn[i+12] *= RADTODEG    // equatorial RA, Dec
    xreturn[i+15] *= RADTODEG    // equatorial RA speed, Dec speed
// dist (index 2, 5, 14, 17) stays in AU
// cartesian (xreturn+6, +18) stays in AU
```

### Final save

```
pdp->xflgs = iflag
pdp->iephe = iflag & SEFLG_EPHMASK
```

---

## 11. `swi_aberr_light()` — Lorentz Aberration

**sweph.c:3699–3736**

Applies the exact relativistic (Lorentz) aberration formula. Modifies `xx[0..2]` in place.

### Variables

```
u[0..2] = xx[0..2]           // planet direction vector (barycentric - observer)
ru = |u|                      // distance
v[0..2] = xe[i+3] / (24 * 3600) / CLIGHT * AUNIT   // observer velocity in units of c
v2 = |v|²
b_1 = sqrt(1 - v2)           // reciprocal Lorentz factor denominator
f1 = dot(u, v) / ru           // u·v̂ (projection)
f2 = 1 + f1 / (1 + b_1)
```

### Position formula

```c
xx[i] = (b_1 * xx[i] + f2 * ru * v[i]) / (1.0 + f1)
```

This is the full relativistic formula, not the classical approximation.

### Speed correction (if SEFLG_SPEED)

Uses interval `PLAN_SPEED_INTV = 0.0001` days:
```
u_prev[i] = xxs[i] - intv * xxs[i+3]   // position at t-intv
ru_prev = |u_prev|
f1_prev = dot(u_prev, v) / ru_prev
f2_prev = 1 + f1_prev / (1 + b_1)
xx2[i] = (b_1 * u_prev[i] + f2_prev * ru_prev * v[i]) / (1 + f1_prev)
dx1 = xx[i] - xxs[i]          // aberration correction at t
dx2 = xx2[i] - u_prev[i]      // aberration correction at t-intv
xx[i+3] += (dx1 - dx2) / intv // speed correction
```

Note: `b_1` uses the same value as the main step (observer velocity is assumed constant).

---

## 12. `swi_deflect_light()` — Gravitational Light Deflection

**sweph.c:3743–3891**

Computes GR light bending by the Sun. Modifies `xx[0..2]` in place.

### Vector definitions

```
U = xx[0..2]                           // planet geocentric vector
E = xearth[0..2] - psdp->x[0..2]     // earth heliocentric (JPL/SWI)
    = xearth[0..2]                     // Moshier (already heliocentric)
Q = xx[0..2] + xearth[0..2] - xsun_retarded[0..2]  // planet heliocentric at t-dt

// xsun_retarded = psdp->x - dt * psdp->x+3   (linear extrapolation)

ru = |U|, rq = |Q|, re = |E|
u = U/ru, q = Q/rq, e = E/re    (unit vectors)
uq = dot(u, q)
ue = dot(u, e)
qe = dot(q, e)
```

### Near-sun handling (meff correction)

```
sina = sqrt(1 - ue²)              // sin(angle between sun and planet)
sin_sunr = SUN_RADIUS / re        // SUN_RADIUS = 959.63/3600 * DEGTORAD

if sina < sin_sunr:
    meff_fact = meff(sina / sin_sunr)
else:
    meff_fact = 1.0
```

### Deflection formula (Explanatory Supplement)

```c
g1 = 2.0 * HELGRAVCONST * meff_fact / CLIGHT / CLIGHT / AUNIT / re
g2 = 1.0 + qe

xx2[i] = ru * (u[i] + g1/g2 * (uq * e[i] - ue * q[i]))
```

Constants:
- `HELGRAVCONST = 1.32712440017987e+20` m³/s² (sweph.h:278, G×M_sun, AA 2006 K6)
- `CLIGHT = 2.99792458e+8` m/s
- `AUNIT = 1.49597870700e+11` m

### Speed correction (if SEFLG_SPEED)

Uses `DEFL_SPEED_INTV = 0.0000005` days. Recomputes deflection at `t + dtsp` with perturbed `u`, `e`, `q`, takes finite difference:
```
dtsp = -DEFL_SPEED_INTV
u_pert = xx[0..2] - dtsp * xx[3..5]
e_pert = xearth - dtsp * (xearth+3 - psdp->x+3)    // JPL/SWI
q_pert = u_pert + xearth - xsun_retarded - dtsp * (xearth+3 - xsun+3)
// recompute g1, g2 for perturbed vectors → xx3
dx1 = xx2[i] - xx[i]
dx2 = xx3[i] - u_pert[i] * ru_pert
xx[i+3] += (dx1 - dx2) / dtsp
```

---

## 13. `meff()` — Effective Solar Mass Factor

**sweph.c:5967–5981**

Returns the effective gravitational mass factor for a photon passing at fractional solar radius `r` (0 = center, 1 = surface). Used to smooth light deflection near solar limb.

### Lookup table `eff_arr` (sweph.c:5858–5966)

Format: `{r, m_eff}` pairs, **sorted descending in r** (from 1.000 down to 0.000):

```
{1.000, 1.000000}, {0.990, 0.999979}, {0.980, 0.999940}, {0.970, 0.999881},
{0.960, 0.999811}, {0.950, 0.999724}, {0.940, 0.999622}, {0.930, 0.999497},
{0.920, 0.999354}, {0.910, 0.999192}, {0.900, 0.999000}, {0.890, 0.998786},
{0.880, 0.998535}, {0.870, 0.998242}, {0.860, 0.997919}, {0.850, 0.997571},
{0.840, 0.997198}, {0.830, 0.996792}, {0.820, 0.996316}, {0.810, 0.995791},
{0.800, 0.995226}, {0.790, 0.994625}, {0.780, 0.993991}, {0.770, 0.993326},
{0.760, 0.992598}, {0.750, 0.991770}, {0.740, 0.990873}, {0.730, 0.989919},
{0.720, 0.988912}, {0.710, 0.987856}, {0.700, 0.986755}, {0.690, 0.985610},
{0.680, 0.984398}, {0.670, 0.982986}, {0.660, 0.981437}, {0.650, 0.979779},
{0.640, 0.978024}, {0.630, 0.976182}, {0.620, 0.974256}, {0.610, 0.972253},
{0.600, 0.970174}, {0.590, 0.968024}, {0.580, 0.965594}, {0.570, 0.962797},
{0.560, 0.959758}, {0.550, 0.956515}, {0.540, 0.953088}, {0.530, 0.949495},
{0.520, 0.945741}, {0.510, 0.941838}, {0.500, 0.937790}, {0.490, 0.933563},
{0.480, 0.928668}, {0.470, 0.923288}, {0.460, 0.917527}, {0.450, 0.911432},
{0.440, 0.905035}, {0.430, 0.898353}, {0.420, 0.891022}, {0.410, 0.882940},
{0.400, 0.874312}, {0.390, 0.865206}, {0.380, 0.855423}, {0.370, 0.844619},
{0.360, 0.833074}, {0.350, 0.820876}, {0.340, 0.808031}, {0.330, 0.793962},
{0.320, 0.778931}, {0.310, 0.763021}, {0.300, 0.745815}, {0.290, 0.727557},
{0.280, 0.708234}, {0.270, 0.687583}, {0.260, 0.665741}, {0.250, 0.642597},
{0.240, 0.618252}, {0.230, 0.592586}, {0.220, 0.565747}, {0.210, 0.537697},
{0.200, 0.508554}, {0.190, 0.478420}, {0.180, 0.447322}, {0.170, 0.415454},
{0.160, 0.382892}, {0.150, 0.349955}, {0.140, 0.316691}, {0.130, 0.283565},
{0.120, 0.250431}, {0.110, 0.218327}, {0.100, 0.186794}, {0.090, 0.156287},
{0.080, 0.128421}, {0.070, 0.102237}, {0.060, 0.077393}, {0.050, 0.054833},
{0.040, 0.036361}, {0.030, 0.020953}, {0.020, 0.009645}, {0.010, 0.002767},
{0.000, 0.000000}
```

Table source: sun_model.c with Stix (The Sun, p. 47) mass distribution × 2.

### Algorithm

```
if r <= 0: return 0.0
if r >= 1: return 1.0
// linear search from top: find first entry with r_entry <= r
for i = 0; eff_arr[i].r > r; i++
// eff_arr[i-1].r > r >= eff_arr[i].r
f = (r - eff_arr[i-1].r) / (eff_arr[i].r - eff_arr[i-1].r)
m = eff_arr[i-1].m + f * (eff_arr[i].m - eff_arr[i-1].m)
return m
```

**Note for Rust**: The table is sorted descending. The loop searches forward until `eff_arr[i].r <= r`, then interpolates between `[i-1]` (higher r) and `[i]` (lower r).

---

## 14. `swi_mean_node()` and `swi_mean_apog()`

### Location

`swemmoon.c:1493` and `swemmoon.c:1564`. These call `mean_elements()` internally.

### Validity range check (both functions)

```
MOSHNDEPH_START = -3100015.5 JD   // 15 Aug -13200
MOSHNDEPH_END   =  8000016.5 JD   // 15 Mar 17191
```

Return `ERR` if outside this range.

### `mean_elements()` — Mean Orbital Elements

**swemmoon.c:1763–1818**

Computes module-level globals: `NF`, `MP`, `D`, `SWELP`, `M` (in arcseconds).

```
T = (J - J2000) / 36525.0
fracT = fmod(T, 1)
```

**Mean anomaly of Sun** (line 1768):
```
M = mods3600(129600000.0 * fracT - 3418.961646 * T + 1287104.76154)
  + (polynomial in T from degree T^2 through T^10)   // high-degree secular terms
```

**Mean arguments** (non-MOSH_MOON_200 branch, lines 1793–1810):

```
NF    = mods3600(1739232000.0 * fracT + 295263.0983 * T
                 - 0.2079419901760 * T + 335779.55755)
      + ((z[2]*T + z[1])*T + z[0])*T²

MP    = mods3600(1717200000.0 * fracT + 715923.4728 * T
                 - 0.2035946368532 * T + 485868.28096)
      + ((z[5]*T + z[4])*T + z[3])*T²

D     = mods3600(1601856000.0 * fracT + 1105601.4603 * T
                 + 0.3962893294503 * T + 1072260.73512)
      + ((z[8]*T + z[7])*T + z[6])*T²

SWELP = mods3600(1731456000.0 * fracT + 1108372.83264 * T
                 - 0.6784914260953 * T + 785939.95571)
      + ((z[11]*T + z[10])*T + z[9])*T²
```

Note: the `fracT` split is a numerical precision technique to reduce cancellation error.

**z[] coefficients** (DE404 version, lines 284–313):

| Index | Value | Meaning |
|---|---|---|
| z[0] | -1.312045233711e+01 | F (NF), t² term (arcsec/cty²) |
| z[1] | -1.138215912580e-03 | F (NF), t³ |
| z[2] | -9.646018347184e-06 | F (NF), t⁴ |
| z[3] | +3.146734198839e+01 | l (MP), t² |
| z[4] | +4.768357585780e-02 | l (MP), t³ |
| z[5] | -3.421689790404e-04 | l (MP), t⁴ |
| z[6] | -6.847070905410e+00 | D, t² |
| z[7] | -5.834100476561e-03 | D, t³ |
| z[8] | -2.905334122698e-04 | D, t⁴ |
| z[9] | -5.663161722088e+00 | L (SWELP), t² |
| z[10] | +5.722859298199e-03 | L (SWELP), t³ |
| z[11] | -8.466472828815e-05 | L (SWELP), t⁴ |
| z[12..] | (planetary) | t² cos/sin of 5 planetary combinations |

Constants used:
- `STR = 4.8481368110953599359e-6` rad/arcsec
- `PI = 3.14159...` (from math)
- `MOON_MEAN_DIST = 384400000.0` m (AA 1996 F2)
- `MOON_MEAN_INCL = 5.1453964` degrees (AA 1996 D2)
- `MOON_MEAN_ECC = 0.054900489` (AA 1996 F2)
- `J2000 = 2451545.0`

### `swi_mean_node()` — Mean Lunar Node

**swemmoon.c:1493–1534**

```
mean_elements()                           // compute SWELP, NF, etc.
dcor = corr_mean_node(J) * 3600          // correction in arcseconds

// longitude (radians):
pol[0] = swi_mod2PI((SWELP - NF - dcor) * STR)
pol[1] = 0.0
pol[2] = MOON_MEAN_DIST / AUNIT          // distance in AU
```

The node = `SWELP - NF` is the mean longitude minus the mean argument of latitude.

### `corr_mean_node()` — Node Correction Table

**swemmoon.c:1470–1486**

Piecewise-linear correction table in degrees, indexed by Gregorian centuries since
`CORR_MNODE_JD_T0GREG = -3063616.5` (1 Jan -13100 Greg.):

```
if J < JPL_DE431_START (-3027215.5): return 0
if J > JPL_DE431_END   (7930192.5):  return 0
dJ = J - CORR_MNODE_JD_T0GREG
i = floor(dJ / 36524.25)     // century index
dfrac = (dJ - i * 36524.25) / 36524.25
dcor = mean_node_corr[i] + dfrac * (mean_node_corr[i+1] - mean_node_corr[i])
```

The `mean_node_corr` array (swemmoon.c:725) has one entry per century from -13100 to 17200, initialised to 0 between roughly years 0–3000 (i.e., inside historical DE406 range). The first few values: `-2.56, -2.473, -2.392347, ...` (degrees).

### `swi_mean_apog()` — Mean Lunar Apogee (Dark Moon / Lilith)

**swemmoon.c:1564–1624**

```
mean_elements()

// base apogee longitude: MP is anomaly relative to node direction
pol[0] = swi_mod2PI((SWELP - MP) * STR + PI)   // add π for apogee (opposite to perigee)
pol[1] = 0
pol[2] = MOON_MEAN_DIST * (1 + MOON_MEAN_ECC) / AUNIT  // apogee distance

// apply apogee correction (similar to node correction):
dcor = corr_mean_apog(J) * DEGTORAD
pol[0] = swi_mod2PI(pol[0] - dcor)

// project onto ecliptic (mean orbital plane → ecliptic):
node = (SWELP - NF) * STR
dcor_node = corr_mean_node(J) * DEGTORAD
node = swi_mod2PI(node - dcor_node)
pol[0] = swi_mod2PI(pol[0] - node)      // longitude relative to node
swi_polcart(pol, pol)                    // polar → cartesian
swi_coortrf(pol, pol, -MOON_MEAN_INCL * DEGTORAD)  // rotate by inclination
swi_cartpol(pol, pol)                    // back to polar
pol[0] = swi_mod2PI(pol[0] + node)      // add node back = ecliptic longitude
```

`MOON_MEAN_INCL = 5.1453964°` is the mean orbital inclination to the ecliptic.

---

## 15. Speed Utility Functions

### `denormalize_positions()` (sweph.c:5983–5997)

Ensures longitude continuity across ±180° boundary for SPEED3 mode:
```
for i in {0, 12}:   // ecliptic lon and RA
    if x1[i] - x0[i] < -180:  x0[i] -= 360
    if x1[i] - x0[i] > 180:   x0[i] += 360
    if x1[i] - x2[i] < -180:  x2[i] -= 360
    if x1[i] - x2[i] > 180:   x2[i] += 360
```

### `calc_speed()` (sweph.c:5999–6011)

Computes speeds from 3 positions using a quadratic formula. The formula accounts for non-uniform acceleration:
```
for j in {0, 6, 12}:   // 3 coordinate sets (ecl, ecl+cart, eq, eq+cart skipping)
    for i in 0..2:
        k = j + i
        b = (x2[k] - x0[k]) / 2        // first difference / 2
        a = (x2[k] + x0[k]) / 2 - x1[k]  // second difference
        x1[k+3] = (2*a + b) / dt        // velocity at center point
```

This is correct for a quadratic fit through (t-dt, x0), (t, x1), (t+dt, x2).

---

## 16. ICRS → J2000 Frame Bias (`swi_bias`)

**swephlib.c:2205–2280**

Applied when `!(iflag & SEFLG_ICRS) && swi_get_denum(ipli) >= 403`.

Two bias models:

**IAU 2000** (default, `SEMOD_BIAS_IAU2000`):
```
rb = [
    [+0.9999999999999942, -0.0000000707827974, +0.0000000805621715],
    [+0.0000000707827948, +0.9999999999999969, +0.0000000330604145],
    [-0.0000000805621738, -0.0000000330604088, +0.9999999999999962],
]
```

**IAU 2006** (`SEMOD_BIAS_IAU2006`):
```
rb = [
    [+0.99999999999999412, -0.00000007078368961, +0.00000008056213978],
    [+0.00000007078368695, +0.99999999999999700, +0.00000003306428553],
    [-0.00000008056214212, -0.00000003306427981, +0.99999999999999634],
]
```

The matrix is applied as `x_J2000 = rb × x_ICRS` (row-major, column = source vector).

JPLHOR_APPROX bypass conditions (swephlib.c:2222–2227):
- If `SEMOD_JPLHORA_2`: skip bias entirely.
- If `SEMOD_JPLHORA_3` and `tjd < DPSI_DEPS_IAU1980_TJD0_HORIZONS (2437684.5)`: skip.

---

## 17. Constants Reference

| Constant | Value | Source |
|---|---|---|
| `EARTH_MOON_MRAT` | 81.30056907419062 | sweph.h:267, DE431 |
| `AUNIT` | 1.49597870700e+11 m | sweph.h:273, DE431 |
| `CLIGHT` | 2.99792458e+8 m/s | sweph.h:274, AA 1996 K6 |
| `HELGRAVCONST` | 1.32712440017987e+20 m³/s² | sweph.h:278, AA 2006 K6 |
| `SUN_RADIUS` | 959.63/3600 × DEGTORAD | sweph.h:281, Meeus |
| `MOON_SPEED_INTV` | 0.00005 days (4.32 s) | sweph.h:298 |
| `PLAN_SPEED_INTV` | 0.0001 days (8.64 s) | sweph.h:299 |
| `MEAN_NODE_SPEED_INTV` | 0.001 days | sweph.h:300 |
| `NODE_CALC_INTV` | 0.0001 days | sweph.h:301 |
| `NODE_CALC_INTV_MOSH` | 0.1 days | sweph.h:302 |
| `DEFL_SPEED_INTV` | 0.0000005 days | sweph.h:304 |
| `STR` | 4.8481368110953599359e-6 rad/arcsec | sweph.h:663 |
| `MOON_MEAN_DIST` | 384400000.0 m | sweph.h:260, AA 1996 F2 |
| `MOON_MEAN_INCL` | 5.1453964° | sweph.h:261, AA 1996 D2 |
| `MOON_MEAN_ECC` | 0.054900489 | sweph.h:262, AA 1996 F2 |
| `MOSHNDEPH_START` | -3100015.5 JD | sweph.h:225 |
| `MOSHNDEPH_END` | 8000016.5 JD | sweph.h:226 |
| `JPL_DE431_START` | -3027215.5 JD | sweph.h:233 |
| `JPL_DE431_END` | 7930192.5 JD | sweph.h:234 |
| `CORR_MNODE_JD_T0GREG` | -3063616.5 JD | swemmoon.c:1466 |

---

## 18. Not Porting (Phase 1)

The following are out of scope for the initial Moshier backend port:

- **File I/O backend** (`sweplan` via `sweph()`, `get_new_segment()`, `read_const()`): requires SWISSEPH file reading. Separate task.
- **JPL backend** (`jplplan()`, `jplplan()`, `swi_pleph()`): requires JPL file. Separate task.
- **EOP/JPLHOR corrections** (Earth Orientation Parameters file): file I/O. Stub.
- **Topocentric corrections** (`swi_get_observer()`): downstream of core pipeline.
- **Sidereal positions** (`swi_get_ayanamsa_with_speed()`): downstream.
- **Osculating elements** (`swi_osc_el_plan()`): asteroid-only.
- **Interpolated apsides** (`intp_apsides()`): special lunar feature.
- **Fictitious planets** (`swi_osc_el_plan()`): Uranian hypotheticals.
- **`swi_bias`**: frame rotation matrix — implement as a simple 3×3 matrix multiply. The matrix values above are verbatim sufficient.

For the **Moshier backend** specifically, the relevant path is:

```
swe_calc()
└── swecalc()
    ├── swi_check_ecliptic() + swi_check_nutation()  [already done]
    └── main_planet()
        └── swi_moshplan() / swi_moshmoon()           [already done]
            └── app_pos_etc_plan() / app_pos_etc_sun() / app_pos_etc_moon()
                ├── light-time: 0 iterations (niter=0), linear: xx -= dt * v
                ├── swi_deflect_light()
                ├── swi_aberr_light()
                ├── swi_bias() (if denum >= 403)
                ├── swi_precess()                     [already done]
                └── app_pos_rest()
                    ├── swi_nutate()                  [already done]
                    ├── swi_coortrf2()
                    ├── swi_cartpol_sp()
                    └── degree conversion
```
