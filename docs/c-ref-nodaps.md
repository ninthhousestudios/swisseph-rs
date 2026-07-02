# C Reference: Nodes & Apsides — swecl.c / sweph.c / swemmoon.c

Porting reference for `swe_nod_aps` (mean + osculating nodes/apsides for all bodies) and for
`SE_TRUE_NODE` / `SE_OSCU_APOG` (osculating lunar node & apogee reached through the normal
`swe_calc` pipeline, i.e. `lunar_osc_elem`). Read this instead of the C source.

Mean lunar node/apogee (`SE_MEAN_NODE` / `SE_MEAN_APOG`) are **already ported** — see
`docs/c-ref-mean-elements.md` and `src/calc.rs::calc_mean_node`/`calc_mean_apogee`. This doc
does not re-describe those; it covers the two things not yet ported:

1. `swe_nod_aps` / `swe_nod_aps_ut` — the standalone public API (mean elements for Sun..Neptune,
   Moon; osculating ellipse elements for everything else, or for everything if
   `SE_NODBIT_OSCU`/`SE_NODBIT_OSCU_BAR` requested).
2. `lunar_osc_elem` — the `SE_TRUE_NODE`/`SE_OSCU_APOG` body computation invoked from inside
   `swe_calc`'s main dispatch (`swecalc`, sweph.c:587).

---

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_nod_aps` | swecl.c:5075–5654 | Yes — new module (e.g. `nodaps.rs`) |
| `swe_nod_aps_ut` | swecl.c:5656–5665 | Yes — deltaT wrapper, same pattern as other `_ut` fns |
| `swi_mean_lunar_elements` | swemmoon.c:1742–1761 | Yes — ~20 lines, used only by `swe_nod_aps`'s Moon/MEAN branch (distinct from `swi_mean_node`/`swi_mean_apog` used by `SE_MEAN_NODE`/`SE_MEAN_APOG`) |
| `swi_plan_for_osc_elem` | sweph.c:5758–5856 | Yes — shared by `swe_nod_aps` osculating branch and `lunar_osc_elem` |
| `lunar_osc_elem` | sweph.c:5168–5594 | Yes — `SE_TRUE_NODE`/`SE_OSCU_APOG` computation |
| `intp_apsides` (SE_INTP_APOG/PERG) | sweph.c:5598+ | Not this task |
| dispatch in `swecalc` for `SE_TRUE_NODE`/`SE_OSCU_APOG` | sweph.c:931–967, 1146–1148 | Document only (informs `Ephemeris::calc` wiring) |
| `get_gmsm` / `TEST_ORBEL_AA` block | swecl.c:5667–5700+ | Not ported — test-only feature behind `#ifdef TEST_ORBEL_AA`, not compiled by default |

---

## Constants

| Constant | Value | Header:line | Notes |
|---|---|---|---|
| `GEOGCONST` | `3.98600448e+14` | sweph.h:279 | G·M(earth), m³/s², AA 1996 K6 |
| `HELGRAVCONST` | `1.32712440017987e+20` | sweph.h:278 | G·M(sun), m³/s², AA 2006 K6 (the `1.32712438e+20` AA1996 value at sweph.h:276 is `#if 0`'d out — not active) |
| `EARTH_MOON_MRAT` | `1 / 0.0123000383` = `81.30056903743774...` | sweph.h:265 | AA 2006 K7; the DE431 (`81.30056907419062`) and DE406 (`81.30056`) alternatives at sweph.h:267/270 are `#if 0`'d out |
| `AUNIT` | `1.49597870700e+11` | sweph.h:273 | metres per AU, DE431 |
| `KGAUSS` | `0.01720209895` | sweph.h:280 | Gaussian gravitational constant — **not referenced by this cluster's code**; listed per task spec but unused in `swe_nod_aps`/`lunar_osc_elem` |
| `MOON_MEAN_DIST` | `384400000.0` | sweph.h:260 | metres, AA1996 F2 |
| `MOON_MEAN_INCL` | `5.1453964` | sweph.h:261 | degrees, AA1996 D2 |
| `MOON_MEAN_ECC` | `0.054900489` | sweph.h:262 | AA1996 F2 |
| `NODE_CALC_INTV` | `0.0001` | sweph.h:301 | days; 3-position speed interval for JPL/SWIEPH-precision osculating node/apogee |
| `NODE_CALC_INTV_MOSH` | `0.1` | sweph.h:302 | days; wider interval for Moshier Moon (its own oscillation is coarser — using `NODE_CALC_INTV` would amplify Moshier moon noise) |
| `MEAN_NODE_SPEED_INTV` | `0.001` | sweph.h:300 | days; only used by mean-element path, not this doc's functions |
| `SE_NODBIT_MEAN` | `1` | swephexp.h:291 | mean nodes/apsides |
| `SE_NODBIT_OSCU` | `2` | swephexp.h:292 | osculating nodes/apsides |
| `SE_NODBIT_OSCU_BAR` | `4` | swephexp.h:293 | osculating about barycenter (planets beyond Jupiter only; heliocentric otherwise) |
| `SE_NODBIT_FOPOINT` | `256` | swephexp.h:294 | return 2nd focal point instead of aphelion |
| `square_sum(x)` | macro | sweph.h:308 | `x[0]²+x[1]²+x[2]²` |
| `dot_prod(x,y)` | macro | sweph.h:309 | `x[0]y[0]+x[1]y[1]+x[2]y[2]` |

All of the above already exist as Rust constants in `src/constants.rs` (`AUNIT`, `HELGRAVCONST`,
`GEOGCONST`, `KGAUSS`, `EARTH_MOON_MRAT`, `MOON_MEAN_DIST`, `MOON_MEAN_INCL`, `MOON_MEAN_ECC`) —
reuse them, do not redefine. `NODE_CALC_INTV`/`NODE_CALC_INTV_MOSH`/`SE_NODBIT_*` do not yet
exist in `src/constants.rs` and must be added.

---

## Part A — `swe_nod_aps` (swecl.c:5075–5654)

### A.0 Static element tables (swecl.c:5012–5074)

VSOP87 mean-equinox-of-date elements for Mercury–Neptune. Each row is `[c0, c1, c2, c3]` used
as `c0 + c1·t + c2·t² + c3·t³` where `t = (tjd_et − J2000) / 36525`. Row order: Mercury, Venus,
Earth, Mars, Jupiter, Saturn, Uranus, Neptune (indices 0–7). Earth rows are all-zero (no
sensible ecliptic node/apsis).

```c
el_node[8][4] = {
  { 48.330893,  1.1861890,  0.00017587,  0.000000211},  // Mercury
  { 76.679920,  0.9011190,  0.00040665, -0.000000080},  // Venus
  {  0,         0,          0,           0},            // Earth
  { 49.558093,  0.7720923,  0.00001605,  0.000002325},  // Mars
  {100.464441,  1.0209550,  0.00040117,  0.000000569},  // Jupiter
  {113.665524,  0.8770970, -0.00012067, -0.000002380},  // Saturn
  { 74.005947,  0.5211258,  0.00133982,  0.000018516},  // Uranus
  {131.784057,  1.1022057,  0.00026006, -0.000000636},  // Neptune
};
el_peri[8][4] = {
  { 77.456119,  1.5564775,  0.00029589,  0.000000056},  // Mercury
  {131.563707,  1.4022188, -0.00107337, -0.000005315},  // Venus
  {102.937348,  1.7195269,  0.00045962,  0.000000499},  // Earth
  {336.060234,  1.8410331,  0.00013515,  0.000000318},  // Mars
  { 14.331309,  1.6126668,  0.00103127, -0.000004569},  // Jupiter
  { 93.056787,  1.9637694,  0.00083757,  0.000004899},  // Saturn
  {173.005159,  1.4863784,  0.00021450,  0.000000433},  // Uranus
  { 48.123691,  1.4262677,  0.00037918, -0.000000003},  // Neptune
};
el_incl[8][4] = {
  {  7.004986,  0.0018215, -0.00001809,  0.000000053},  // Mercury
  {  3.394662,  0.0010037, -0.00000088, -0.000000007},  // Venus
  {  0,         0,          0,           0},            // Earth
  {  1.849726, -0.0006010,  0.00001276, -0.000000006},  // Mars
  {  1.303270, -0.0054966,  0.00000465, -0.000000004},  // Jupiter
  {  2.488878, -0.0037363, -0.00001516,  0.000000089},  // Saturn
  {  0.773196,  0.0007744,  0.00003749, -0.000000092},  // Uranus
  {  1.769952, -0.0093082, -0.00000708,  0.000000028},  // Neptune
};
el_ecce[8][4] = {
  {  0.20563175,  0.000020406, -0.0000000284, -0.00000000017},  // Mercury
  {  0.00677188, -0.000047766,  0.0000000975,  0.00000000044},  // Venus
  {  0.01670862, -0.000042037, -0.0000001236,  0.00000000004},  // Earth
  {  0.09340062,  0.000090483, -0.0000000806, -0.00000000035},  // Mars
  {  0.04849485,  0.000163244, -0.0000004719, -0.00000000197},  // Jupiter
  {  0.05550862, -0.000346818, -0.0000006456,  0.00000000338},  // Saturn
  {  0.04629590, -0.000027337,  0.0000000790,  0.00000000025},  // Uranus
  {  0.00898809,  0.000006408, -0.0000000008, -0.00000000005},  // Neptune
};
el_sema[8][4] = {
  {  0.387098310,  0.0,          0.0,          0.0},           // Mercury
  {  0.723329820,  0.0,          0.0,          0.0},           // Venus
  {  1.000001018,  0.0,          0.0,          0.0},           // Earth
  {  1.523679342,  0.0,          0.0,          0.0},           // Mars
  {  5.202603191,  0.0000001913, 0.0,          0.0},           // Jupiter
  {  9.554909596,  0.0000021389, 0.0,          0.0},           // Saturn
  { 19.218446062, -0.0000000372, 0.00000000098, 0.0},          // Uranus
  { 30.110386869, -0.0000001663, 0.00000000069, 0.0},          // Neptune
};
```

`plmass[9]` — Sun-mass/planet-mass ratios, indexed 0..8 = Mercury, Venus, EMB, Mars, Jupiter,
Saturn, Uranus, Neptune, Pluto:

```c
plmass[9] = {
  6023600,        // Mercury
   408523.719,    // Venus
   328900.5,      // Earth and Moon (EMB)
  3098703.59,     // Mars
     1047.348644, // Jupiter
     3497.9018,   // Saturn
    22902.98,     // Uranus
    19412.26,     // Neptune
136566000,        // Pluto
};
```

`ipl_to_elem[15]` — maps a body number `ipl` (indexed directly by the `SE_*` body constant, NOT
a contiguous 0..14 planet count) to a row index into `el_node`/`el_peri`/`el_incl`/`el_ecce`/
`el_sema` (rows 0..7 = Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune) and,
reused, into `plmass` (rows 0..8, same order plus Pluto=8):

```c
ipl_to_elem[15] = {2, 0, 0, 1, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 2};
```

The body-number enum actually is (swephexp.h:101–115): `SE_SUN=0, SE_MOON=1, SE_MERCURY=2,
SE_VENUS=3, SE_MARS=4, SE_JUPITER=5, SE_SATURN=6, SE_URANUS=7, SE_NEPTUNE=8, SE_PLUTO=9,
SE_MEAN_NODE=10, SE_TRUE_NODE=11, SE_MEAN_APOG=12, SE_OSCU_APOG=13, SE_EARTH=14` — note
`SE_EARTH` is **14**, not adjacent to the other planets. Mapping each index to its table row:

| `ipl` (index) | body | → row | table row meaning |
|---|---|---|---|
| 0 | SE_SUN | 2 | Earth's row — the Sun's apparent node/apsis motion mirrors Earth's heliocentric orbit |
| 1 | SE_MOON | 0 | **unused** — Moon never reaches the table lookup; A.3.1 handles it via `swi_mean_lunar_elements` instead |
| 2 | SE_MERCURY | 0 | Mercury |
| 3 | SE_VENUS | 1 | Venus |
| 4 | SE_MARS | 3 | Mars |
| 5 | SE_JUPITER | 4 | Jupiter |
| 6 | SE_SATURN | 5 | Saturn |
| 7 | SE_URANUS | 6 | Uranus |
| 8 | SE_NEPTUNE | 7 | Neptune |
| 9 | SE_PLUTO | 0 | **quirk**: resolves to Mercury's row when used for `plmass[]` in A.4.1's `plm` term (Pluto's *own* row would be `plmass[8]`, but `ipl_to_elem[9]=0`); astronomically inconsequential since `plm` is a `1/mass_ratio` perturbation added to `1` in `Gmsm*(1+plm)` and is negligible regardless of which planet's ratio is picked — port literally, do not "fix" |
| 10–13 | SE_MEAN_NODE..SE_OSCU_APOG | 0 | **unused** — these `ipl` values are rejected outright in A.1, never reach this table |
| 14 | SE_EARTH | 2 | Earth (only relevant when `swe_nod_aps` is called directly with `ipl=SE_EARTH`; A.5.2 additionally zeroes the node/desc-node outputs for Earth since "no nodes for earth" — see the `ipli==SE_EARTH && ij<=1` check) |

### A.1 Rejected `ipl` values (swecl.c:5137–5158)

Before anything else, `swe_nod_aps` rejects and zero-fills all four output vectors, returning
`ERR`, for:
```
ipl == SE_MEAN_NODE || ipl == SE_TRUE_NODE ||
ipl == SE_MEAN_APOG || ipl == SE_OSCU_APOG ||
ipl < 0 ||
(ipl >= SE_NPLANETS && ipl <= SE_AST_OFFSET)
```
i.e. you cannot ask `swe_nod_aps` for the nodes/apsides *of* a node/apogee point itself, nor for
negative or the reserved fictitious-body range. Error text: `"nodes/apsides for planet %5.0f
are not implemented"`.

Special-case remap (swecl.c:5116–5117, *before* the rejection check): asteroid-number Pluto
(`SE_AST_OFFSET + 134340`) is remapped to `SE_PLUTO`.

### A.2 Setup (swecl.c:5075–5160)

```
t = (tjd_et - J2000) / 36525
iflag &= ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX)
xna = xx+0; xnd = xx+6; xpe = xx+12; xap = xx+18   // 4 output slots, 6 doubles each, in one 24-double scratch array
swi_force_app_pos_etc()   // invalidate save-area cache — stateless port: no-op, nothing to invalidate
method %= SE_NODBIT_FOPOINT   // strip the FOPOINT bit (256) out of method for branch dispatch; do_focal_point captured separately
ipli = ipl; if ipl == SE_SUN: ipli = SE_EARTH
if ipl == SE_MOON:
    do_defl = FALSE
    if !(iflag & SEFLG_HELCTR): do_aberr = FALSE
iflg0 = (iflag & (SEFLG_EPHMASK|SEFLG_NONUT)) | SEFLG_SPEED | SEFLG_TRUEPOS
if ipli != SE_MOON: iflg0 |= SEFLG_HELCTR
xx[0..23] = 0
```
`do_aberr = !(iflag & (SEFLG_TRUEPOS | SEFLG_NOABERR))`, `do_defl = !(iflag & SEFLG_TRUEPOS) &&
!(iflag & SEFLG_NOGDEFL)`, `do_focal_point = method & SE_NODBIT_FOPOINT` (tested **before** the
`method %= SE_NODBIT_FOPOINT` line, using the *original* method value).

### A.3 Mean branch (swecl.c:5161–5245)

Condition: `(method == 0 || (method & SE_NODBIT_MEAN)) && ((SE_SUN <= ipl <= SE_NEPTUNE) ||
ipl == SE_EARTH)`. Note this is a *different* condition from A.1's Moon check — Moon *is*
included here (`ipl == SE_MOON` inside `SE_SUN..SE_NEPTUNE`'s numeric range, since `SE_MOON`'s
enum value falls between `SE_SUN` and `SE_NEPTUNE`... actually verify: in C, `SE_SUN=0,
SE_MOON=1, SE_MERCURY=2,...,SE_NEPTUNE=8`, so `SE_MOON` is numerically inside `[SE_SUN,
SE_NEPTUNE]` and matches the range test). Pluto and everything else (asteroids, fictitious)
falls through to the osculating branch (A.4) unconditionally, regardless of `method`.

#### A.3.1 Moon sub-branch
```
swi_mean_lunar_elements(tjd_et, &xna[0], &xna[3], &xpe[0], &xpe[3])  // node lon+speed, peri lon+speed (degrees)
incl = MOON_MEAN_INCL;  vincl = 0
ecce = MOON_MEAN_ECC;   vecce = 0
sema = MOON_MEAN_DIST / AUNIT;  vsema = 0
```

#### A.3.2 Planet sub-branch (all bodies except Moon, in the mean-eligible set)
```
iplx = ipl_to_elem[ipl]
ep = el_incl[iplx]; incl  = ep[0] + ep[1]*t + ep[2]*t*t + ep[3]*t*t*t;  vincl = ep[1] / 36525
ep = el_sema[iplx]; sema  = ep[0] + ep[1]*t + ep[2]*t*t + ep[3]*t*t*t;  vsema = ep[1] / 36525
ep = el_ecce[iplx]; ecce  = ep[0] + ep[1]*t + ep[2]*t*t + ep[3]*t*t*t;  vecce = ep[1] / 36525
ep = el_node[iplx]
xna[0] = ep[0] + ep[1]*t + ep[2]*t*t + ep[3]*t*t*t   // ascending node longitude, degrees
xna[3] = ep[1] / 36525                                // node "speed" — NOTE: this is deg/century-derivative
                                                       // reused as deg/day placeholder; actually see below,
                                                       // it's overwritten via degnorm difference at the end
ep = el_peri[iplx]
xpe[0] = ep[0] + ep[1]*t + ep[2]*t*t + ep[3]*t*t*t   // perihelion longitude, degrees
xpe[3] = ep[1] / 36525
```
`vincl`/`vsema`/`vecce` are **not** speeds in deg/day — they are the linear (`c1`) term divided by
36525, i.e. deg/century → deg/day-ish scaling used only as a finite first-order correction when
computing "distance a fraction of a day later" quantities below (they are consistently paired
with `t` in AU/day-like usage but the C code never multiplies by an actual `dt`; it treats
`sema + vsema` etc. as "value one century later" divided out implicitly through the `xpe[3]`
finite-difference pattern that follows). **Port this literally, do not "fix" the units** — it
reproduces C's exact (slightly informal) linear extrapolation.

#### A.3.3 Shared post-processing (both sub-branches, swecl.c:5195–5245)

```
xnd[0] = swe_degnorm(xna[0] + 180)          // descending node = ascending + 180°
xnd[3] = xna[3]
parg  = xpe[0] = swe_degnorm(xpe[0] - xna[0])                    // arg of perihelion from node
pargx = xpe[3] = swe_degnorm(xpe[0] + xpe[3] - xna[3])            // NOTE: uses the just-updated xpe[0] (parg), not the original
// transform arg-of-perihelion from orbital plane to mean ecliptic of date:
swe_cotrans(xpe, xpe, -incl)          // rotates (lon=xpe[0], lat=0, dist=xpe[2]) by -incl about x-axis-in-lon/lat sense
swe_cotrans(xpe+3, xpe+3, -incl - vincl)   // xpe+3 here is treated as an AUXILIARY POSITION for pargx, NOT a speed vector (comment in C: "xpe+3 is aux. position, not speed!!!")
xpe[0] = swe_degnorm(xpe[0] + xna[0])      // add node back
xpe[3] = swe_degnorm(xpe[3] + xna[0] + xna[3])   // add node + node-speed back to the aux value
xpe[3] = swe_degnorm(xpe[3] - xpe[0])      // finite difference => actual longitude "speed" of perihelion (deg/day, since xna[3] is deg/day)
// heliocentric distance of peri/aphelion:
xpe[2] = sema * (1 - ecce)
xpe[5] = (sema + vsema) * (1 - ecce - vecce) - xpe[2]   // finite-difference "speed" of perihelion distance
xap[0] = swe_degnorm(xpe[0] + 180)   // aphelion longitude = perihelion + 180
xap[1] = -xpe[1]
xap[3] = xpe[3]
xap[4] = -xpe[4]
if do_focal_point:
    xap[2] = sema * ecce * 2
    xap[5] = (sema + vsema) * (ecce + vecce) * 2 - xap[2]
else:
    xap[2] = sema * (1 + ecce)
    xap[5] = (sema + vsema) * (1 + ecce + vecce) - xap[2]
```

**Node/descending-node distances from the osculating ellipse** (this is the key non-obvious
step — node points get a *distance* derived from where the ellipse crosses the node direction,
not simply `sema`):
```
ea  = atan(tan(-parg  * DEGTORAD / 2) * sqrt((1-ecce)/(1+ecce))) * 2         // eccentric anomaly at ascending node
eax = atan(tan(-pargx * DEGTORAD / 2) * sqrt((1-ecce-vecce)/(1+ecce+vecce))) * 2  // same, one day later (for speed)
xna[2] = sema * (cos(ea) - ecce) / cos(parg * DEGTORAD)
xna[5] = (sema+vsema) * (cos(eax) - ecce - vecce) / cos(pargx * DEGTORAD)
xna[5] -= xna[2]      // finite-difference distance speed

ea  = atan(tan((180-parg)  * DEGTORAD / 2) * sqrt((1-ecce)/(1+ecce))) * 2
eax = atan(tan((180-pargx) * DEGTORAD / 2) * sqrt((1-ecce-vecce)/(1+ecce+vecce))) * 2
xnd[2] = sema * (cos(ea) - ecce) / cos((180-parg) * DEGTORAD)
xnd[5] = (sema+vsema) * (cos(eax) - ecce - vecce) / cos((180-pargx) * DEGTORAD)
xnd[5] -= xnd[2]
```
No light-time correction applied here ("speed is extremely small" per comment).

**Degrees → radians → cartesian**, for all four points (`xx` iterated in 4 groups of 6, i.e.
`xna, xnd, xpe, xap`):
```
for each of the 4 six-element blocks xp:
    xp[0] *= DEGTORAD; xp[1] *= DEGTORAD; xp[3] *= DEGTORAD; xp[4] *= DEGTORAD
    swi_polcart_sp(xp, xp)   // polar (lon,lat,dist,dlon,dlat,ddist) -> cartesian, in place
```
This falls through to the **shared output pipeline** (A.5) — it does NOT `return` early.

### A.4 Osculating branch (swecl.c:5249–5400)

Entered when the A.3 condition is false — i.e. `method` requests `SE_NODBIT_OSCU`/
`SE_NODBIT_OSCU_BAR` explicitly, or `ipl` is Pluto/asteroid/fictitious (always osculating
regardless of `method`).

#### A.4.1 Reference heliocentric/barycentric distance and Gmsm

```
if swe_calc(tjd_et, ipli, iflg0, x, serr) == ERR: return ERR   // x[2] = heliocentric distance (AU), used below to pick barycentric-vs-heliocentric and to scale dt
iflJ2000 = (iflag & SEFLG_EPHMASK) | SEFLG_J2000 | SEFLG_EQUATORIAL | SEFLG_XYZ | SEFLG_TRUEPOS | SEFLG_NONUT | SEFLG_SPEED
ellipse_is_bary = FALSE
if ipli != SE_MOON:
    if (method & SE_NODBIT_OSCU_BAR) && x[2] > 6:   // only bodies beyond ~Jupiter's distance (6 AU threshold, hardcoded)
        iflJ2000 |= SEFLG_BARYCTR
        ellipse_is_bary = TRUE
    else:
        iflJ2000 |= SEFLG_HELCTR
```
Note `x[2] > 6` is a **hardcoded AU threshold**, not tied to any planet enum — it's "beyond
Jupiter" in the sense of "further than 6 AU", which happens to sit between Jupiter (~5.2 AU) and
Saturn (~9.5 AU).

```
if ipli == SE_MOON:
    dt = NODE_CALC_INTV                       // 0.0001 days
    dzmin = 1e-15
    Gmsm = GEOGCONST * (1 + 1/EARTH_MOON_MRAT) / AUNIT^3 * 86400^2
else:
    plm = 1/plmass[ipl_to_elem[ipl]] if (SE_MERCURY <= ipli <= SE_PLUTO || ipli == SE_EARTH) else 0
    dt = NODE_CALC_INTV * 10 * x[2]           // scaled by heliocentric distance — farther bodies need a wider dt
    dzmin = 1e-15 * dt / NODE_CALC_INTV
    Gmsm = HELGRAVCONST * (1 + plm) / AUNIT^3 * 86400^2
```
`Gmsm` is GM in units of AU³/day² (the `/AUNIT/AUNIT/AUNIT*86400.0*86400.0` converts from
m³/s² to AU³/day²). Write literally as three separate `/AUNIT` divisions then `*86400.0*86400.0`
if bitwise fidelity to the C evaluation order matters (division is not generally associative
under FP, though here it's benign since it's `((v/A)/A)/A * D * D`).

`dzmin` is a minimum-radial-speed floor used below to avoid division by ~zero.

```
if iflag & SEFLG_SPEED: istart=0; iend=2
else: istart=iend=0; dt=0
```

#### A.4.2 Three-position loop (swecl.c:5286–5299)

```
for i in istart..=iend, t = tjd_et - dt (i=0), tjd_et (i=1 when istart==iend), tjd_et + dt (i=2):
    if istart == iend: t = tjd_et    // single-position (no-speed) case always uses tjd_et exactly
    swe_calc(t, ipli, iflJ2000, xpos[i], serr)   // heliocentric/barycentric J2000 equatorial cartesian, true position, no nutation, with speed
    if ipli == SE_EARTH:
        swe_calc(t, SE_MOON, iflJ2000 & ~(SEFLG_BARYCTR|SEFLG_HELCTR), xposm, serr)
        xpos[i][j] += xposm[j] / (EARTH_MOON_MRAT + 1.0)   for j in 0..=5   // Earth -> EMB correction
    swi_plan_for_osc_elem(iflg0, t, xpos[i])   // rotate J2000 -> ecliptic-of-date frame (see Part B)
```
The loop index arithmetic `t += dt` after each iteration means: i=0 → `tjd_et - dt`, i=1 →
`tjd_et` (since `-dt+dt=tjd_et`), i=2 → `tjd_et + dt`. When `istart==iend==0` the `t = tjd_et`
override inside the loop body makes it a single evaluation at `tjd_et` regardless of `dt`
(which is itself forced to 0 in that case).

#### A.4.3 Node direction (per position, swecl.c:5300–5309)

```
for i in istart..=iend:
    if |xpos[i][5]| < dzmin: xpos[i][5] = dzmin       // floor tiny radial speed to avoid div-by-zero
    fac = xpos[i][2] / xpos[i][5]                      // z / dz — "time to cross ecliptic plane" proxy
    sgn = xpos[i][5] / |xpos[i][5]|                    // sign of radial(z) speed
    for j in 0..=2:
        xn[i][j] = (xpos[i][j] - fac * xpos[i][j+3]) * sgn   // project position back along velocity to z=0 crossing, oriented by sgn
        xs[i][j] = -xn[i][j]                                  // descending node = antipode
```
This gives the **direction** of the node from where the velocity-extrapolated position vector
crosses the ecliptic plane (z=0) — an approximate distance that gets replaced by the true
ellipse-based distance in A.4.4.

#### A.4.4 Osculating ellipse elements (per position, swecl.c:5310–5377)

```
for i in istart..=iend:
    // node longitude direction (from xn, the crossing-point vector)
    rxy = sqrt(xn[i][0]^2 + xn[i][1]^2)
    cosnode = xn[i][0]/rxy;  sinnode = xn[i][1]/rxy

    // inclination, from angular momentum h = r × v
    xnorm = swi_cross_prod(xpos[i], xpos[i]+3)     // xnorm = pos × vel (specific angular momentum direction)
    rxy = xnorm[0]^2 + xnorm[1]^2
    c2  = rxy + xnorm[2]^2                          // |h|^2
    rxyz = sqrt(c2); rxy = sqrt(rxy)
    sinincl = rxy / rxyz
    cosincl = sqrt(1 - sinincl^2)
    if xnorm[2] < 0: cosincl = -cosincl              // retrograde (e.g. 20461 Dioretsa)

    // argument of latitude uu = atan2(sinu, cosu)
    cosu = xpos[i][0]*cosnode + xpos[i][1]*sinnode
    sinu = xpos[i][2] / sinincl
    uu = atan2(sinu, cosu)

    // vis-viva semi-major axis
    rxyz = sqrt(square_sum(xpos[i]))                 // |r|
    v2   = square_sum(xpos[i]+3)                     // |v|^2
    sema = 1 / (2/rxyz - v2/Gmsm)

    // eccentricity from specific angular momentum (semi-latus rectum p = h^2/Gmsm)
    pp = c2 / Gmsm
    ecce = sqrt(1 - pp/sema)

    // eccentric anomaly from r = a(1 - e cosE) and the radial-velocity relation
    cosE = 1/ecce * (1 - rxyz/sema)
    sinE = 1/ecce/sqrt(sema*Gmsm) * dot_prod(xpos[i], xpos[i]+3)

    // true anomaly (half-angle formula)
    ny = 2 * atan( sqrt((1+ecce)/(1-ecce)) * sinE / (1+cosE) )

    // perihelion direction: distance of perihelion from ascending node
    xq[i] = [ swi_mod2PI(uu - ny), 0 /* lat */, sema*(1-ecce) /* dist */ ]
    swi_polcart(xq[i], xq[i])
    swi_coortrf2(xq[i], xq[i], -sinincl, cosincl)     // tilt from orbital plane to ecliptic
    swi_cartpol(xq[i], xq[i])
    xq[i][0] += atan2(sinnode, cosnode)               // add node longitude
    xa[i] = [ swi_mod2PI(xq[i][0] + PI), -xq[i][1],
              do_focal_point ? sema*ecce*2 : sema*(1+ecce) ]   // aphelion (or 2nd focus)
    swi_polcart(xq[i], xq[i])   // xq, xa reconverted to cartesian (both were polar after the atan2 add)
    swi_polcart(xa[i], xa[i])

    // recompute the node's distance from the OSCULATING ELLIPSE (not the tangent-line approx from A.4.3)
    ny  = swi_mod2PI(ny - uu)          // true anomaly AT the ascending node
    ny2 = swi_mod2PI(ny + PI)          // true anomaly AT the descending node
    cosE  = cos(2*atan(tan(ny/2)  / sqrt((1+ecce)/(1-ecce))))
    cosE2 = cos(2*atan(tan(ny2/2) / sqrt((1+ecce)/(1-ecce))))
    rn  = sema*(1 - ecce*cosE)         // true ellipse-based ascending-node distance
    rn2 = sema*(1 - ecce*cosE2)        // descending-node distance
    ro  = sqrt(square_sum(xn[i]))      // old (tangent-line) distance from A.4.3
    ro2 = sqrt(square_sum(xs[i]))
    for j in 0..=2:
        xn[i][j] *= rn/ro              // rescale direction vector to the correct ellipse distance
        xs[i][j] *= rn2/ro2
```

#### A.4.5 Assemble output + speed (swecl.c:5378–5399)

```
for i in 0..=2:
    if iflag & SEFLG_SPEED:
        xpe[i]   = xq[1][i];  xpe[i+3] = (xq[2][i] - xq[0][i]) / dt / 2   // central difference
        xap[i]   = xa[1][i];  xap[i+3] = (xa[2][i] - xa[0][i]) / dt / 2
        xna[i]   = xn[1][i];  xna[i+3] = (xn[2][i] - xn[0][i]) / dt / 2
        xnd[i]   = xs[1][i];  xnd[i+3] = (xs[2][i] - xs[0][i]) / dt / 2
    else:
        xpe[i] = xq[0][i]; xpe[i+3]=0;  xap[i]=xa[0][i]; xap[i+3]=0
        xna[i] = xn[0][i]; xna[i+3]=0;  xnd[i]=xs[0][i]; xnd[i+3]=0
is_true_nodaps = TRUE
```
Note: in the no-speed case (`istart==iend==0`), index `[1]` was never populated by the loop —
all three loop variables end up written only at index 0. But this branch reads `xq[0]` etc., so
that's consistent (`iend=0` means only `i=0` executed in A.4.2–A.4.4, and here the `else` arm
correctly reads index 0, not 1). In the speed case, `xq[1]`/`xa[1]`/`xn[1]`/`xs[1]` are the
values *at* `tjd_et` (the "t+=dt" loop lands i=1 exactly on `tjd_et`), and `[0]`/`[2]` are
`tjd_et∓dt` for the central difference.

### A.5 Shared output-transform pipeline (swecl.c:5401–5652)

Runs for **both** A.3 (mean) and A.4 (osculating) branches — `is_true_nodaps` gates a few
osculating-only steps (nutation-to-equator rotation, and later re-adding it back).

#### A.5.1 Re-establish observer frame (swecl.c:5401–5436)

```
if ipli==SE_MOON && (iflag & (HELCTR|BARYCTR)):
    swi_force_app_pos_etc(); swe_calc(tjd_et, SE_SUN, iflg0, x, serr)   // just to populate save-area sun/earth/nutation — stateless port: recompute directly instead
else:
    swe_calc(tjd_et, ipli, iflg0 | (iflag & SEFLG_TOPOCTR), x, serr)

// observer position xobs:
if iflag & SEFLG_TOPOCTR:
    swi_get_observer(tjd_et, iflag, FALSE, xobs, serr)   // topocentric offset
else:
    xobs = 0

if iflag & (HELCTR|BARYCTR):
    if (iflag & HELCTR) && !(iflag & SEFLG_MOSEPH): xobs = xsun   // xsun = swed.pldat[SEI_SUNBARY].x (global!)
elif ipl == SE_SUN && !(iflag & SEFLG_MOSEPH):
    xobs = xsun
else:
    xobs += xear   // xear = swed.pldat[SEI_EARTH].x (global!)

oe = (iflag & SEFLG_J2000) ? &swed.oec2000 : &swed.oec   // GLOBAL obliquity-of-date/J2000 cache
```

#### A.5.2 Per-point loop (4 points: node, descnode, peri, aphe), swecl.c:5445–5640

```
for ij in 0..=3, xp = one of {xna, xnd, xpe, xap}:
    if ipli == SE_EARTH && ij <= 1:        // no ascending/descending node for Earth itself
        xp[0..5] = 0; continue

    // to equator:
    if is_true_nodaps && !(iflag & SEFLG_NONUT):
        swi_coortrf2(xp, xp, -swed.nut.snut, swed.nut.cnut)      // ecliptic-of-date -> "mean-ecliptic-like" via -nutation
        if SPEED: swi_coortrf2(xp+3, xp+3, -swed.nut.snut, swed.nut.cnut)
    swi_coortrf2(xp,   xp,   -oe->seps, oe->ceps)                 // ecliptic -> equatorial (of date or J2000, per oe)
    swi_coortrf2(xp+3, xp+3, -oe->seps, oe->ceps)

    if is_true_nodaps && !(iflag & SEFLG_NONUT):
        swi_nutate(xp, iflag, TRUE)          // full nutation matrix, back=TRUE (removes nutation, mean->true... actually applies inverse direction)

    // to J2000:
    swi_precess(xp, tjd_et, iflag, J_TO_J2000)
    if SPEED: swi_precess_speed(xp, tjd_et, iflag, J_TO_J2000)

    // to barycenter:
    if ipli == SE_MOON:
        xp += xear   // xear = swed.pldat[SEI_EARTH].x, GLOBAL — geocentric moon-node -> add earth's barycentric position
    else:
        if !(iflag & SEFLG_MOSEPH) && !ellipse_is_bary:
            xp += xsun   // xsun = swed.pldat[SEI_SUNBARY].x, GLOBAL

    // to correct center (subtract observer):
    xp -= xobs
    if ipl == SE_SUN && !(iflag & (HELCTR|BARYCTR)):
        xp = -xp        // geocentric perigee/apogee OF the sun is defined as -(sun's own node/apsis vector), i.e. Earth's apsis mirrored

    // light deflection:
    dt = |xp| * AUNIT / CLIGHT / 86400.0
    if do_defl: swi_deflect_light(xp, dt, iflag)

    // aberration:
    if do_aberr:
        swi_aberr_light(xp, xobs, iflag)
        if SPEED:
            // recompute xobs at tjd_et - dt to get d(xobs)/dt contribution to apparent speed
            swe_calc(tjd_et - dt, ipli, iflg0 | TOPOCTR?, x2, serr)   // side effect: repopulates swed.topd/global caches for t-dt
            xobs2 = (topocentric ? swed.topd.xobs : 0) or xsun or (xear + 0)   // same logic as A.5.1's xobs derivation, evaluated at t-dt
            xp[3..5] += xobs[3..5] - xobs2[3..5]
            swe_calc(tjd_et, SE_SUN, iflg0 | TOPOCTR?, x2, serr)   // restore global save-area state clobbered by the t-dt call above

    // precession back to date (unless J2000 requested):
    x2000 = xp   // save J2000-frame copy for sidereal use
    if !(iflag & SEFLG_J2000):
        swi_precess(xp, tjd_et, iflag, J2000_TO_J)
        if SPEED: swi_precess_speed(xp, tjd_et, iflag, J2000_TO_J)

    // nutation:
    if !(iflag & SEFLG_NONUT): swi_nutate(xp, iflag, FALSE)
    pldat.xreturn[18..23] = xp          // equatorial cartesian, saved

    // to ecliptic:
    swi_coortrf2(xp, xp, oe->seps, oe->ceps)
    if SPEED: swi_coortrf2(xp+3, xp+3, oe->seps, oe->ceps)
    if !(iflag & SEFLG_NONUT):
        swi_coortrf2(xp, xp, swed.nut.snut, swed.nut.cnut)
        if SPEED: swi_coortrf2(xp+3, xp+3, swed.nut.snut, swed.nut.cnut)
    pldat.xreturn[6..11] = xp           // ecliptic cartesian, saved

    // sidereal (SEFLG_SIDEREAL): same three sub-cases as app_pos_rest (ECL_T0 / SSY_PLANE / traditional ayanamsa subtraction) — see docs/c-ref-mean-elements.md §8 for the shared pattern

    // output selection (XYZ/EQUATORIAL combinations short-circuit via `continue`):
    if XYZ && EQUATORIAL: xp = pldat.xreturn[18..23]; continue
    if XYZ:               xp = pldat.xreturn[6..11];  continue
    // else convert both blocks to polar+speed, degrees (unless RADIANS):
    swi_cartpol_sp(pldat.xreturn+18, pldat.xreturn+12)
    swi_cartpol_sp(pldat.xreturn+6,  pldat.xreturn)
    if !(iflag & SEFLG_RADIANS):
        xreturn[0,1,3,4] *= RADTODEG; xreturn[12,13,15,16] *= RADTODEG
    xp = EQUATORIAL ? pldat.xreturn[12..17] : pldat.xreturn[0..5]
```

**Important subtlety**: `pldat` (`struct plan_data pldat;`) is a **local stack variable** in
`swe_nod_aps`, declared fresh at function entry (swecl.c:5098) — its `xreturn` array is not the
global `swed.pldat[...]` save area. So `pldat.xreturn` here is *not* shared/global state; it's
scratch space reused across the 4-point loop. This is different from `lunar_osc_elem` (Part C),
which writes directly into `swed.nddat[SEI_TRUE_NODE]`/`swed.nddat[SEI_OSCU_APOG]`.

#### A.5.3 Final zero-speed cleanup + output copy (swecl.c:5641–5652)

```
for i in 0..=5:
    if i > 2 && !(iflag & SEFLG_SPEED): xna[i]=xnd[i]=xpe[i]=xap[i]=0
    xnasc[i] = xna[i]; xndsc[i] = xnd[i]; xperi[i] = xpe[i]; xaphe[i] = xap[i]   (each only if output ptr non-NULL)
```

### A.6 `swe_nod_aps_ut` (swecl.c:5656–5665)

```
swe_nod_aps_ut(tjd_ut, ipl, iflag, method, xnasc, xndsc, xperi, xaphe, serr):
    return swe_nod_aps(tjd_ut + swe_deltat_ex(tjd_ut, iflag, serr), ipl, iflag, method, ...)
```
Same pattern as every other `_ut` wrapper in the codebase — convert UT→ET via deltaT, delegate.

---

## Part B — `swi_mean_lunar_elements` (swemmoon.c:1742–1761)

Used **only** by `swe_nod_aps`'s Moon/MEAN branch (A.3.1) — distinct from `swi_mean_node`/
`swi_mean_apog` (used by `SE_MEAN_NODE`/`SE_MEAN_APOG`, already ported, `c-ref-mean-elements.md`).
This variant additionally returns **numerical speeds** for both node and perigee longitude in
one call, via a 1/36525-day (1-century... actually re-read: `T -= 1.0/36525` is **1 day**, since
`T` is in Julian centuries and `1/36525` centuries = 1 day) backward finite difference.

```
swi_mean_lunar_elements(tjd, &node, &dnode, &peri, &dperi):
    T = (tjd - J2000) / 36525.0;  T2 = T*T
    mean_elements()                              // sets SWELP, NF, MP (see c-ref-mean-elements.md §1)
    node = swe_degnorm((SWELP - NF) * STR * RADTODEG)     // STR = arcsec-to-radians; result then RADTODEG'd back to degrees
    peri = swe_degnorm((SWELP - MP) * STR * RADTODEG)
    T -= 1.0/36525                                // step back exactly 1 day
    mean_elements()                                // recompute SWELP, NF, MP one day earlier
    dnode = swe_degnorm(node - (SWELP - NF) * STR * RADTODEG)
    dnode -= 360                                   // fold into (-360, 0] so that dnode is a small NEGATIVE-ish daily rate, not a degnorm-wrapped value near 0/360
    dperi = swe_degnorm(peri - (SWELP - MP) * STR * RADTODEG)
    dcor = corr_mean_node(tjd); node = swe_degnorm(node - dcor)
    dcor = corr_mean_apog(tjd); peri = swe_degnorm(peri - dcor)
```
`corr_mean_node`/`corr_mean_apog` are the same piecewise-linear correction tables documented in
`c-ref-mean-elements.md` §2–3, §9–10 (`mean_node_corr[]`/`mean_apsis_corr[]`) — reuse those, do
not re-transcribe.

Note the odd `node`/`peri` (not `dnode`/`dperi`) receive the correction subtraction at the end,
**not** the "one day earlier" values — the corrections are only ever applied to the `tjd`
(current) longitude, not to the backward-stepped one used for the speed estimate. This means
the returned `dnode`/`dperi` speeds are **uncorrected** (a day-to-day difference of the raw,
uncorrected mean longitude), while `node`/`peri` themselves are corrected. This is a **known
minor inconsistency in the C code** — port it exactly as-is (do not "fix" it by also correcting
the day-earlier value), since golden-test parity requires bit-for-bit reproduction of this
asymmetry.

STR (arcsec→radian) and `mean_elements()` itself are documented in `c-ref-mean-elements.md` §1;
`RADTODEG`/`swe_degnorm` are existing helpers (`src/math.rs` / `crate::math`).

---

## Part C — `swi_plan_for_osc_elem` (sweph.c:5758–5856)

Rotates a J2000 equatorial cartesian position+speed vector (position AND speed, 6 doubles) into
the ecliptic-of-date frame needed as input for osculating-ellipse computation. Called once per
position sample, both from `swe_nod_aps`'s osculating branch (A.4.2) and from `lunar_osc_elem`
(Part D).

```
swi_plan_for_osc_elem(iflag, tjd, xx[6]):        // xx modified in place
    // ICRS -> J2000 bias, only if ephemeris is DE403+ and not already ICRS:
    if !(iflag & SEFLG_ICRS) && swi_get_denum(SEI_SUN, iflag) >= 403:
        swi_bias(xx, tjd, iflag, FALSE)

    // precession J2000 -> equator of date (position AND speed, each precessed SEPARATELY —
    // speed is rotated by the SAME epoch's precession matrix as the position, NOT differentiated
    // — the comment says "daily precession 0.137" may not be added", i.e. treat the precession
    // as a pure rotation of the velocity vector, not (d/dt) of the precession matrix itself):
    if !(iflag & SEFLG_J2000):
        swi_precess(xx,   tjd, iflag, J2000_TO_J)
        swi_precess(xx+3, tjd, iflag, J2000_TO_J)     // NOTE: precesses the speed vector via the SAME swi_precess (position) function, not swi_precess_speed
        oe = (tjd == swed.oec.teps) ? &swed.oec
           : (tjd == J2000)         ? &swed.oec2000
           : calc_epsilon(tjd, iflag, &oectmp); oe = &oectmp   // fresh obliquity computation if no cache hit
    else:
        oe = &swed.oec2000

    // nutation (position AND speed rotated by the SAME nutation matrix — no nutation "rate" added):
    if !(iflag & SEFLG_NONUT):
        nutp = (tjd == swed.nut.tnut) ? &swed.nut
             : (tjd == J2000)         ? &swed.nut2000
             : (tjd == swed.nutv.tnut) ? &swed.nutv
             : { swi_nutation(tjd, iflag, nuttmp.nutlo); nuttmp.tnut=tjd;
                 nuttmp.snut=sin(nuttmp.nutlo[1]); nuttmp.cnut=cos(nuttmp.nutlo[1]);
                 nut_matrix(&nuttmp, oe); &nuttmp }
        x[0..2] = xx[0..2] · nutp->matrix   (matrix-vector product, row-major as: x[i] = Σ_k xx[k]*matrix[k][i])
        x[3..5] = xx[3..5] · nutp->matrix   (same matrix, speed is a pure rotation — no nutation-rate term)
        xx = x

    // to ecliptic:
    swi_coortrf2(xx,   xx,   oe->seps, oe->ceps)
    swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps)

    // (SID_TNODE_FROM_ECL_T0 is a compile-time #ifdef; confirmed NOT defined anywhere in the
    // C tree — grep for `SID_TNODE_FROM_ECL_T0` across *.h/*.c finds only its own #ifdef/#endif
    // guards, no #define. Dead code in the default build. Skip it entirely.)
    if !(iflag & SEFLG_NONUT):
        swi_coortrf2(xx, xx, nutp->snut, nutp->cnut)
        swi_coortrf2(xx+3, xx+3, nutp->snut, nutp->cnut)
    return
```

**Global-state reads here that a stateless port must replace with explicit recomputation:**
`swed.oec.teps`/`swed.oec`/`swed.oec2000` (cached obliquity-of-date/J2000 — recompute via
`obliquity(tjd, ...)` every time, no cache-hit shortcut needed), `swed.nut.tnut`/`swed.nut`/
`swed.nut2000`/`swed.nutv` (cached nutation — recompute via `nutation(tjd, ...)` every time),
`swi_get_denum(SEI_SUN, iflag)` (ephemeris DE-number — determines whether frame bias is applied;
port needs an explicit "which JPL DE number is this ephemeris backend" query, or, if the
Moshier/SWIEPH backends in this Rust port never claim DE≥403 bias applicability the same way C
does, confirm against `docs/c-ref-precession.md`/bias-handling doc for how bias is already
gated elsewhere in this codebase before assuming it's always-on or always-off here).

---

## Part D — `lunar_osc_elem` (sweph.c:5168–5594) — `SE_TRUE_NODE` / `SE_OSCU_APOG`

This is invoked from `swecalc`'s dispatch (sweph.c:931–967) when `ipl == SE_TRUE_NODE` or
`ipl == SE_OSCU_APOG`. Both body types are computed **together** in one call (the node is always
needed even when only the apogee was requested, and vice versa — see D.4), each writing into
its own `swed.nddat[SEI_TRUE_NODE]`/`swed.nddat[SEI_OSCU_APOG]` save slot.

### D.0 Caching / early-return (sweph.c:5199–5213)

```
ndp = &swed.nddat[ipl]     // ipl here is the INTERNAL index SEI_TRUE_NODE or SEI_OSCU_APOG
flg1 = iflag & ~EQUATORIAL & ~XYZ;  flg2 = ndp->xflgs & ~EQUATORIAL & ~XYZ
if tjd == ndp->teval && tjd != 0 && flg1==flg2 && (!speedf2 || speedf1):
    ndp->xflgs = iflag; ndp->iephe = iflag & SEFLG_EPHMASK
    return OK     // cache hit — reuse previous xreturn
```
Stateless port: **no cache, always recompute.** This is purely a performance optimization in C
with no numerical effect (same iflag/tjd always yields the same output) — safe to drop entirely.

### D.1 Three lunar positions (sweph.c:5231–5359)

```
epheflag = MOSEPH|SWIEPH|JPLEPH per iflag  (default handled by outer caller)
swed.pldat[SEI_MOON].teval = 0    // force fresh moon computation — GLOBAL STATE RESET, stateless port: n/a, just recompute moon position directly
istart = (iflag & SEFLG_SPEED) ? 0 : 2   // no-speed case only computes index 2 (t = tjd exactly)

for i in istart..=2:
    t = (i==0) ? tjd - speed_intv : (i==1) ? tjd + speed_intv : tjd
    // per-ephemeris moon position, backend-specific:
    switch epheflag:
      JPLEPH:  speed_intv = NODE_CALC_INTV;      xp = jplplan(t, MOON, iflag, ...)
      SWIEPH:  speed_intv = NODE_CALC_INTV;       xp = swemoon(t, iflag|SPEED, ...)
      MOSEPH:  speed_intv = NODE_CALC_INTV_MOSH;  xp = swi_moshmoon(t, ...)   // wider interval — Moshier moon's own short-period terms would otherwise dominate the finite difference
    // light-time correction (JPL/SWIEPH paths only; NOT for TRUEPOS):
    if !(iflag & SEFLG_TRUEPOS) && retc >= OK:
        dt = |xpos[i]| * AUNIT / CLIGHT / 86400.0
        xp = <backend>(t - dt, ...)     // re-evaluate at light-time-corrected epoch (SIMPLE dt*speed subtraction explicitly NOT used — comment: "the error would be greater than the advantage of computation speed")
    // fallback cascade on NOT_AVAILABLE/BEYOND_EPH_LIMITS: JPL->SWIEPH->MOSEPH, retry via `goto three_positions` (restart the WHOLE 3-position loop with the new epheflag, not just the failed sample)
    xpos[i] = swi_plan_for_osc_elem(iflag|SEFLG_SPEED, t, xpos[i])   // Part C — rotates JPL/SWIEPH-frame position (already ecliptic-of-date-ish for these backends — see note below) into final osc-elem-ready frame
```
**Frame note**: unlike `swe_nod_aps`'s osculating branch (which explicitly requests J2000
equatorial cartesian via `iflJ2000` before calling `swi_plan_for_osc_elem`), here the raw
per-backend moon position (`jplplan`/`swemoon`/`swi_moshmoon`) is fed directly into
`swi_plan_for_osc_elem` — these backend functions return **J2000 equatorial cartesian** by
convention in this part of the C codebase (verify against the Moshier/SWIEPH moon backend's
output convention already established in this Rust port, e.g. `src/moshier/moon.rs`, before
assuming — but structurally this is the same J2000-equatorial-cartesian input contract as A.4.2).

`speed_intv` is `NODE_CALC_INTV` (0.0001 d) for JPL/SWIEPH, `NODE_CALC_INTV_MOSH` (0.1 d) for
Moshier — reused as the central-difference `dt` in D.2/D.3.

### D.2 Node with speed (sweph.c:5360–5393)

Same geometry as A.4.3 (tangent-line node direction) applied to the 3 lunar positions:
```
for i in istart..=2:
    if |xpos[i][5]| < 1e-15: xpos[i][5] = 1e-15
    fac = xpos[i][2]/xpos[i][5];  sgn = xpos[i][5]/|xpos[i][5]|
    xx[i][0..2] = (xpos[i][0..2] - fac*xpos[i][3..5]) * sgn
ndnp = &swed.nddat[SEI_TRUE_NODE]
ndnp->x[0..2] = xx[2]                          // position at tjd (index 2, since istart..2 always ends at tjd)
if SPEED:
    b = (xx[1][i]-xx[0][i])/2;  a = (xx[1][i]+xx[0][i])/2 - xx[2][i]
    ndnp->x[i+3] = (2*a+b) / speed_intv          // quadratic-interpolation derivative at t (same 3-point formula as calc_speed / calc_speed_3point in this Rust port — src/calc.rs:662)
else:
    ndnp->x[i+3] = 0
```
**This node distance (`xx[i]`) is later corrected onto the true osculating ellipse in D.3** —
identical pattern to A.4.4's `xn`/`xs` rescale, just computed alongside the apogee this time
(single combined loop, not two separate ones).

### D.3 Apogee with speed + ellipse-corrected node distance (sweph.c:5394–5472)

`Gmsm = GEOGCONST * (1 + 1/EARTH_MOON_MRAT) / AUNIT^3 * 86400^2` — **identical formula** to
A.4.1's Moon case (reuse the same helper/constant, do not re-derive).

```
ndap = &swed.nddat[SEI_OSCU_APOG]
for i in istart..=2:      // per-position ellipse elements — same algebra as A.4.4 but "apogee" (uu-ny+PI) instead of "perihelion" (uu-ny)
    rxy = sqrt(xx[i][0]^2+xx[i][1]^2); cosnode=xx[i][0]/rxy; sinnode=xx[i][1]/rxy
    xnorm = swi_cross_prod(xpos[i], xpos[i]+3)
    rxy = xnorm[0]^2+xnorm[1]^2; c2 = rxy+xnorm[2]^2; rxyz=sqrt(c2); rxy=sqrt(rxy)
    sinincl = rxy/rxyz; cosincl = sqrt(1-sinincl^2)     // NOTE: no retrograde (xnorm[2]<0) sign flip here, unlike A.4.4 — moon's inclination is never retrograde so this case doesn't arise, but PORT LITERALLY (don't add the flip)
    cosu = xpos[i][0]*cosnode + xpos[i][1]*sinnode; sinu = xpos[i][2]/sinincl; uu = atan2(sinu,cosu)
    rxyz = sqrt(square_sum(xpos[i])); v2 = square_sum(xpos[i]+3)
    sema = 1/(2/rxyz - v2/Gmsm)
    pp = c2/Gmsm; ecce = sqrt(1 - pp/sema)
    cosE = 1/ecce*(1-rxyz/sema); sinE = 1/ecce/sqrt(sema*Gmsm)*dot_prod(xpos[i], xpos[i]+3)
    ny = 2*atan(sqrt((1+ecce)/(1-ecce)) * sinE/(1+cosE))
    xxa[i] = [ swi_mod2PI(uu - ny + PI), 0, sema*(1+ecce) ]    // apogee = perihelion + PI, distance = a(1+e) UNCONDITIONALLY (no do_focal_point option here — lunar_osc_elem has no SE_NODBIT_FOPOINT equivalent)
    swi_polcart(xxa[i], xxa[i]); swi_coortrf2(xxa[i], xxa[i], -sinincl, cosincl); swi_cartpol(xxa[i], xxa[i])
    xxa[i][0] += atan2(sinnode, cosnode)
    swi_polcart(xxa[i], xxa[i])
    // ellipse-corrected NODE distance (reusing this position's ecce/sema/uu):
    ny = swi_mod2PI(ny - uu)
    cosE = cos(2*atan(tan(ny/2) / sqrt((1+ecce)/(1-ecce))))
    r0 = sema*(1-ecce*cosE)          // true node distance
    r1 = sqrt(square_sum(xx[i]))    // old tangent-line distance (from D.2)
    xx[i][0..2] *= r0/r1             // rescale in place
// save:
for i in 0..=2:
    ndap->x[i] = xxa[2][i]
    ndap->x[i+3] = SPEED ? (xxa[1][i]-xxa[0][i])/speed_intv/2 : 0    // simple CENTRAL difference here, NOT the quadratic 3-point formula used for the node (D.2) — different speed formulas for apogee vs node!
    ndnp->x[i] = xx[2][i]
    ndnp->x[i+3] = SPEED ? (xx[1][i]-xx[0][i])/speed_intv/2 : 0      // ALSO overwrites node speed with plain central difference here — supersedes the quadratic formula computed in D.2! (D.2's ndnp->x[i+3] assignment is DEAD — this second loop's assignment is what actually survives)
```
**Critical detail a porter would miss**: D.2 computes `ndnp->x[i+3]` using the 3-point quadratic
formula `(2a+b)/speed_intv`, but D.3's final save loop **unconditionally overwrites**
`ndnp->x[i+3]` with the plain central difference `(xx[1][i]-xx[0][i])/speed_intv/2`. The D.2
quadratic-speed computation is therefore dead code for the node (it's recomputed and discarded).
Only `ndap` (apogee) speed uses a single formula (central difference, computed once, in D.3).
**Port only the D.3 formula for node speed** — do not implement D.2's quadratic version at all,
since it's never actually used in the output.

### D.4 Output-frame assembly (sweph.c:5473–5594)

```
for j in 0,1:     // SEI_TRUE_NODE, then SEI_OSCU_APOG
    ndp = swed.nddat[j==0 ? SEI_TRUE_NODE : SEI_OSCU_APOG]
    ndp->xreturn[0..23] = 0
    ndp->xreturn[6..11] = ndp->x[0..5]                       // ecliptic cartesian (already ecliptic-of-date + light-time + no aberration/deflection — see note below)
    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn)             // ecliptic polar
    swi_coortrf2(ndp->xreturn+6, ndp->xreturn+18, -oe->seps, oe->ceps)     // -> equatorial cartesian
    if SPEED: swi_coortrf2(ndp->xreturn+9, ndp->xreturn+21, -oe->seps, oe->ceps)
    // (SID_TNODE_FROM_ECL_T0 branch: compile-time disabled — see Part C note — skip)
    if !(iflag & SEFLG_NONUT):
        swi_coortrf2(ndp->xreturn+18, ndp->xreturn+18, -swed.nut.snut, swed.nut.cnut)
        if SPEED: swi_coortrf2(ndp->xreturn+21, ndp->xreturn+21, -swed.nut.snut, swed.nut.cnut)
    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12)          // equatorial polar
    ndp->xflgs = iflag; ndp->iephe = iflag & SEFLG_EPHMASK

    // SEFLG_SIDEREAL (traditional ayanamsa subtraction) or SEFLG_J2000 re-projection — see comment block at sweph.c:5522-5578; same three-way branch pattern as elsewhere (ECL_T0 rigorous / SSY_PLANE rigorous / traditional), plus a SEPARATE elif branch for plain SEFLG_J2000 (re-precess ndp->xreturn+18 from date to J2000 and rebuild ecliptic/polar from there) — this J2000 elif is NOT present in app_pos_rest (c-ref-mean-elements.md §8), it's specific to lunar_osc_elem because the earlier steps in D.1-D.3 work in ecliptic-OF-DATE throughout, unlike the mean-element/planet paths which branch J2000-vs-date earlier.

    // degrees:
    xreturn[0,1,3,4] *= RADTODEG; xreturn[12,13,15,16] *= RADTODEG
    xreturn[0] = swe_degnorm(xreturn[0]); xreturn[12] = swe_degnorm(xreturn[12])
```

**No light-time/aberration/deflection/ICRS-bias step here** — comment at sweph.c:5473-5477:
"precession and nutation have already been taken into account because the computation is on the
basis of lunar positions that have gone through `swi_plan_for_osc_elem`. light-time is already
contained in lunar positions." I.e. all of that happened per-sample back in D.1 (light-time via
the `t-dt` re-evaluation, precession+nutation via `swi_plan_for_osc_elem`/Part C) — D.4 is purely
a coordinate-representation/output-frame step, not a physical-effects step.

### D.5 Dispatch wiring in `swecalc` (sweph.c:931–967, 1146–1148)

```
} else if ipl == SE_TRUE_NODE:
    if iflag & (HELCTR|BARYCTR): x[0..23]=0; return iflag    // heliocentric/barycentric node is meaningless, rejected
    ndp = &swed.nddat[SEI_TRUE_NODE]; xp = ndp->xreturn
    lunar_osc_elem(tjd, SEI_TRUE_NODE, iflag, serr)
    iflag = ndp->xflgs
    if !(iflag & SIDEREAL) && !(iflag & J2000):
        ndp->xreturn[1]=ndp->xreturn[4]=ndp->xreturn[8]=ndp->xreturn[11] = 0.0   // force exact zero latitude/z (suppress FP noise — node is by definition ON the ecliptic in non-sidereal/non-J2000 output)
} else if ipl == SE_OSCU_APOG:
    if iflag & (HELCTR|BARYCTR): x[0..23]=0; return iflag
    ndp = &swed.nddat[SEI_OSCU_APOG]; xp = ndp->xreturn
    lunar_osc_elem(tjd, SEI_OSCU_APOG, iflag, serr)
    iflag = ndp->xflgs
    // NO latitude-zeroing here — apogee is NOT constrained to the ecliptic (only the node is, by construction)
...
// (falls through to the shared tail, shared by every ipl branch:)
x[0..23] = xp[0..23]     // sweph.c:1146-1147 — final output copy, straight from whichever ndp/pdp->xreturn was populated
return iflag
```
Both branches call `lunar_osc_elem` — since D.1-D.4 compute node AND apogee together in one call
regardless of which was requested, calling it once for `SE_TRUE_NODE` also populates
`swed.nddat[SEI_OSCU_APOG]` as a side effect (and vice versa) in C. **Stateless port**: since
there's no save-area to share, a Rust `calc_true_node`/`calc_oscu_apogee` pair would each run the
full D.1-D.4 computation and simply discard the half they don't need — OR, better, factor D.1-D.3
into one shared helper returning both `(node_xreturn, apogee_xreturn)` and have thin
`Ephemeris::calc` dispatch arms pick the half they want, mirroring the C's "compute once, use
either" structure without needing a cache. This avoids duplicating D.1-D.3's ~150 lines of
ellipse math for two call sites.

---

## Porting Notes — Global State Inventory

A stateless Rust port must replace every one of these with an explicit, request-scoped
computation (no caching across calls):

| C global | Read at | Rust replacement |
|---|---|---|
| `swed.pldat[SEI_SUNBARY].x` (`xsun`) | swecl.c:5099, 5428, 5431, 5485, 5529, 5532 (A.5.1/A.5.2), sweph.c (implicitly via `swe_calc(SE_SUN, BARYCTR...)`) | Explicit `Ephemeris::calc(jd, Body::Sun, BARYCTR\|...)` call to get barycentric sun position at the needed epoch — note A.5.2's aberration-speed correction needs this at BOTH `tjd_et` and `tjd_et - dt` |
| `swed.pldat[SEI_EARTH].x` (`xear`) | swecl.c:5100, 5435, 5479, 5481 (A.5.1/A.5.2), sweph.c D.1 (Earth->EMB correction reads `swed.pldat[SEI_EARTH]` implicitly via nested `swe_calc`) | Explicit heliocentric/barycentric Earth call at the same epoch |
| `swed.oec` / `swed.oec2000` | swecl.c:5438-5441, 5460-5461, 5573 (A.5.1/A.5.2); sweph.c:5764 (`swi_plan_for_osc_elem`), 5190 (`lunar_osc_elem`) | `obliquity(jd, flags, models)` — already exists (`src/obliquity.rs` per other ref docs); no cache needed, just call it at each epoch used |
| `swed.nut` / `swed.nut2000` / `swed.nutv` | swecl.c:5456-5458, 5514-5517, 5576-5579 (A.5.2); sweph.c:5809-5822 (`swi_plan_for_osc_elem`), 5513-5516 (D.4) | `nutation(jd, flags, models)` — already exists; recompute per epoch, no `tnut`-equality cache-hit shortcut needed |
| `swed.topd.xobs` | swecl.c:5521, 5524 (A.5.2 aberration-speed sub-step) | Explicit observer-position helper at `tjd_et - dt` (topocentric parallax) — check if this already exists for other topocentric bodies in this port (e.g. `phenomena.rs`/`houses.rs` observer geometry) before re-deriving |
| `ndp->teval` / `ndp->xflgs` cache-hit checks (D.0, and analogous in A via `swi_force_app_pos_etc`) | sweph.c:5206-5213 | Drop entirely — always recompute; no numerical difference, pure perf optimization in C |
| `swed.pldat[SEI_MOON].teval = 0` forced-reset | sweph.c:5243 | N/A — Rust has no shared moon save-slot to invalidate; just call the moon backend fresh each time |
| `swed.sidd` (sidereal mode/ayanamsa state) | scattered (SEFLG_SIDEREAL branches, A.5.2 and D.4) | Reuse whatever sidereal-mode plumbing already exists for planets (`AstroModels`/`swed.sidd`-equivalent in this Rust port) — this doc doesn't re-derive the sidereal math, see `c-ref-mean-elements.md` §8 for the shared 3-way branch pattern |
| Save-area clobber/restore dance (A.5.2: `swe_calc(tjd_et - dt, ...)` then re-`swe_calc(tjd_et, SE_SUN, ...)` "to restore save area") | swecl.c:5514, 5545 | **Irrelevant in a stateless port** — this entire two-call dance exists purely to repair C's global save-area after a probing call clobbers it. A stateless `Ephemeris::calc` has no shared state to clobber, so just call it once at `tjd_et - dt` for the xobs2 values needed and move on; drop the "restore" call entirely. |
| Cascading ephemeris-fallback `goto three_positions` / `goto sweph_moon` retries on `NOT_AVAILABLE`/`BEYOND_EPH_LIMITS` (D.1) | sweph.c:5251, 5358-5359 | Match this port's existing fallback-flag convention (`CalcResult.flags_used` per project's `<architecture>` conventions) rather than re-deriving C's retry-loop/goto structure — check how `main_planet`'s equivalent fallback is already handled elsewhere in this codebase (e.g. `src/moshier/`, `sweph_file.rs`) for the established pattern |
| `SE_NODBIT_OSCU_BAR`'s `x[2] > 6` hardcoded AU threshold (A.4.1) | swecl.c:5256 | Just a literal constant, not global state — but flag it in code review since "6" isn't named/documented in C either; consider a named constant, e.g. `OSCU_BAR_DISTANCE_THRESHOLD_AU: f64 = 6.0` |

## Porting Notes — Reuse Existing Helpers

The following already exist in this Rust codebase and should be reused rather than
reimplemented:

- `src/math.rs::cross_prod` ↔ `swi_cross_prod` (A.4.4, D.3)
- `src/math.rs::rotate_x_sincos` ↔ `swi_coortrf2` (used pervasively throughout A.5, D.4, Part C)
- `src/math.rs::polar_to_cartesian_with_speed` / `cartesian_to_polar_with_speed` ↔
  `swi_polcart_sp`/`swi_cartpol_sp`
- `src/math.rs::normalize_radians` ↔ `swi_mod2PI`
- `src/calc.rs::calc_speed_3point` (quadratic 3-point derivative) ↔ the `(2a+b)/dt` formula in
  A.4.5/D.2 — **but note D.3's final node-speed save uses a plain central difference, not this
  quadratic formula** (see D.3's critical-detail callout); don't reflexively reuse
  `calc_speed_3point` for the node without checking which of the two formulas actually survives
  to output
- `src/math.rs::diff_radians` — same role as `swe_difrad2n` used elsewhere in the mean-element
  path (not used in this doc's functions, but same pattern if needed for any new speed calc)
- `src/constants.rs` — `AUNIT`, `HELGRAVCONST`, `GEOGCONST`, `EARTH_MOON_MRAT`, `MOON_MEAN_DIST`,
  `MOON_MEAN_INCL`, `MOON_MEAN_ECC` all already defined; `NODE_CALC_INTV` (`0.0001`),
  `NODE_CALC_INTV_MOSH` (`0.1`), and `SE_NODBIT_*` (as a bitflags-style `NodApsMethod` type, not
  raw ints) still need to be added
- `docs/c-ref-mean-elements.md` §1 `mean_elements()`, §2-3 correction-table interpolation, §9-10
  correction tables — reused verbatim by `swi_mean_lunar_elements` (Part B); do not re-transcribe
  the 304-entry tables into this doc's implementation, import/call the existing ported function
- `docs/c-ref-mean-elements.md` §8 `app_pos_rest` sidereal 3-way branch pattern — same shape
  reappears in A.5.2 and D.4, reuse whatever Rust function already implements it for planets if
  one exists, or note it as a needed shared extraction if not
