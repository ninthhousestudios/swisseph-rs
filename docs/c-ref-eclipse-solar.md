# C Reference: Solar Eclipses & Shared Eclipse Helpers — swecl.c

Porting reference for solar-eclipse geometry, local circumstances, and global/local eclipse
search. This is the largest of the eclipse ref docs; it also documents the helpers
(`calc_planet_star`, `eclipse_where`, `eclipse_how`, `find_maximum`, `find_zero`) that are
reused verbatim by the lunar-eclipse and occultation modules — read this doc first before
`c-ref-eclipse-lunar.md` / `c-ref-eclipse-occultation.md` (if/when they exist).

All line numbers refer to `swecl.c` unless stated otherwise.

---

## Function Map

| C function | Location | Purpose |
|---|---|---|
| `calc_planet_star` (static) | swecl.c:888–897 | Shared dispatch: planet via `swe_calc`, star via `swe_fixstar` |
| `eclipse_where` (static) | swecl.c:640–886 | THE CORE: shadow-cone geometry, geographic position of eclipse center, dcore[] |
| `swe_sol_eclipse_where` | swecl.c:565–582 | Public wrapper: `eclipse_where` + `eclipse_how` at the central point |
| `eclipse_how` (static) | swecl.c:967–1152 | Local circumstances (magnitude, obscuration, az/alt, saros) at a given observer |
| `swe_sol_eclipse_how` | swecl.c:922–964 | Public wrapper: `eclipse_how` + redundant az/alt recheck + visibility gate |
| `swe_sol_eclipse_when_glob` | swecl.c:1185–1515 | Global search: next/previous solar eclipse anywhere on Earth |
| `swe_sol_eclipse_when_loc` | swecl.c:2019–2041 | Public wrapper over `eclipse_when_loc` + core-shadow diameter fill |
| `eclipse_when_loc` (static) | swecl.c:2100–2410 | THE CORE: local search — next/previous eclipse visible from a given geo. position |
| `find_maximum` (static) | swecl.c:4133–4146 | Parabolic-interpolation maximum finder (3-point) |
| `find_zero` (static) | swecl.c:4148–4162 | Parabolic-interpolation zero finder (3-point, 2 roots) |

Related but out of scope for this doc (belongs to the occultation/lunar-eclipse ref docs):
`swe_lun_occult_where` (swecl.c:606–630, thin wrapper reusing `eclipse_where`/`eclipse_how`),
`swe_lun_occult_when_glob` (swecl.c:1572+), `occult_when_loc` (swecl.c:2412+),
`swe_lun_occult_when_loc` (swecl.c:2071–2098), `lun_eclipse_how` (lunar eclipses).

---

## 0. Shared Constants

```c
#define DSUN    (1392000000.0 / AUNIT)      // Sun diameter, AU  (swecl.c:80; #if 0 branch uses 1391978489.9)
#define DMOON   (3476300.0 / AUNIT)          // Moon diameter, AU (swecl.c:82)
#define DEARTH  (6378140.0 * 2 / AUNIT)       // Earth equatorial diameter, AU (swecl.c:83)
#define RSUN    (DSUN / 2)
#define RMOON   (DMOON / 2)
#define REARTH  (DEARTH / 2)
```
`AUNIT = 1.49597870700e+11` m (DE431 value, sweph.h:273).

`EARTH_OBLATENESS = 1.0 / 298.25642` (AA 2006 K6 flattening `f`, sweph.h:284).

`SEI_ECL_GEOALT_MIN = -500.0`, `SEI_ECL_GEOALT_MAX = 25000.0` (meters above sea level; valid
range for observer height in eclipse local-circumstance functions, sweph.h:198–199).

`J2000 = 2451545.0` (sweph.h:67). `SAROS_CYCLE = 6585.3213` days, `NSAROS_SOLAR = 181`
(swecl.c:114–115), with a static table `saros_data_solar[]` of `{series_no, tstart}` pairs
(swecl.c:116ff) derived from NASA's eclipse Saros catalogue.

Body physical radius lookup (used generically for solar eclipses AND occultations, since
`eclipse_where`/`eclipse_how` take an arbitrary `ipl`/`starname`):
```c
#define NDIAM (SE_VESTA + 1)
static const double pla_diam[NDIAM] = { 1392000000.0 /*Sun*/, 3475000.0 /*Moon*/, ... };
```
(sweph.h:314–330; note `pla_diam[SE_SUN]` uses `1392000000.0`, matching `DSUN`'s numerator
exactly — so `drad` for `ipl==SE_SUN` reduces to exactly `RSUN`.)

Body radius resolution pattern, used identically in `eclipse_where` (swecl.c:697–704) and
`eclipse_how` (swecl.c:1004–1011):
```c
if (starname != NULL && *starname != '\0')  drad = 0;                              // stars: point source
else if (ipl < NDIAM)                        drad = pla_diam[ipl] / 2 / AUNIT;      // table lookup, meters -> AU
else if (ipl > SE_AST_OFFSET)                drad = swed.ast_diam / 2 * 1000 / AUNIT; // named asteroid, km -> m -> AU
else                                          drad = 0;
```
`swed.ast_diam` is a global set by the asteroid-orbital-element loader (SE1 file read) —
STATELESS PORT NOTE: the Rust port must resolve asteroid diameter from whatever explicit
asteroid-data structure it threads through, not a global.

Shared macros: `square_sum(x) = x0²+x1²+x2²`, `dot_prod(x,y) = x0y0+x1y1+x2y2` (sweph.h:308–309).
`swi_dot_prod_unit(x,y)` (swephlib.c:453–465) computes the dot product of the two vectors
**normalized by their own magnitudes** (i.e. cosine of the angle between them), clamped to
`[-1,1]` before the caller applies `acos`. **FP-fidelity hazard:** several call sites in this
file pre-normalize their vectors (divide by distance) before calling `swi_dot_prod_unit`, which
then re-normalizes internally — this is a redundant (but not incorrect) double division. Port
literally; do not "simplify away" the pre-normalization, as it changes rounding.

`swi_cartpol(x, l)` / `swi_polcart(l, x)` (swephlib.c:314–366): polar `l = [lon_rad, lat_rad, radius]`,
cartesian `x = [x,y,z]`, both in whatever distance unit the input used (AU throughout this file).

Eclipse-type return-flag bitmask (swephexp.h:307–331), used as the return value of nearly every
function in this doc:
```c
SE_ECL_CENTRAL          1        SE_ECL_NONCENTRAL       2
SE_ECL_TOTAL            4        SE_ECL_ANNULAR          8
SE_ECL_PARTIAL          16       SE_ECL_ANNULAR_TOTAL    32   (alias SE_ECL_HYBRID)
SE_ECL_PENUMBRAL        64  (lunar eclipses only)
SE_ECL_VISIBLE          128
SE_ECL_MAX_VISIBLE      256
SE_ECL_1ST_VISIBLE      512      SE_ECL_2ND_VISIBLE      1024
SE_ECL_3RD_VISIBLE      2048     SE_ECL_4TH_VISIBLE      4096
SE_ECL_ONE_TRY          32768    // search-control input flag, ORed into `backward`
```
**Return-value convention**, uniform across this whole module: `ERR` (`-1`) on hard failure
(ephemeris lookup failed); `0` = "no eclipse" (the sentinel meaning "nothing found here/now",
distinct from `ERR`); otherwise a bitwise-OR of the flags above. Callers must test `< 0` for
error, not `!= 0`/`== 0` alone, since `0` is a valid non-error "no eclipse" result.

---

## 1. `calc_planet_star` — shared body/star dispatch (swecl.c:888–897)

```c
static int32 calc_planet_star(double tjd_et, int32 ipl, char *starname, int32 iflag,
                               double *x, char *serr)
{
  int retc = OK;
  if (starname == NULL || *starname == '\0')
    retc = swe_calc(tjd_et, ipl, iflag, x, serr);
  else
    retc = swe_fixstar(starname, tjd_et, iflag, x, serr);
  return retc;
}
```
Trivial but pervasive: every eclipse/occultation function that needs "the position of the
eclipsed/occulted body" (Sun for solar eclipses, any planet/asteroid/star for occultations)
routes through this. `tjd_et` is **always** Ephemeris/Dynamical Time (ET), never UT — callers
convert UT→ET via `swe_deltat_ex` before calling. Note it calls `swe_fixstar` (the "slow" name-
resolving entry point), **not** `swe_fixstar2` — contrary to what a purely-structural reading of
newer swisseph headers might suggest; verified directly at swecl.c:894.

Rust port: this becomes a small helper taking `BodyOrStar` (or equivalent enum) and dispatching
to the corresponding `Ephemeris` method; no change in behavior, just an enum match instead of a
null/empty-string check.

---

## 2. `iflag` construction conventions used throughout this file

Three flag variables recur with a consistent naming pattern:
- **`iflag`** (or `iflag2`): the "polar" flavor — `SEFLG_EQUATORIAL | ifl` (+ `SEFLG_SPEED` where
  speeds are needed), requesting right ascension/declination/distance (**degrees** by default,
  or **radians** if `SEFLG_RADIANS` is additionally ORed in). `ifl` is the caller-supplied
  ephemeris selector (`SEFLG_SWIEPH`/`SEFLG_JPLEPH`/`SEFLG_MOSEPH`, masked to
  `SEFLG_EPHMASK = SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH` at public-API entry points, swecl.c:67).
- **`iflagcart`**: `iflag | SEFLG_XYZ` — same coordinate frame, but cartesian x,y,z (AU).
- **`SEFLG_TOPOCTR`** is added whenever the calculation is for a specific observer location
  (`eclipse_how`, `eclipse_when_loc`) and omitted for geocentric-only work (`eclipse_where`,
  `swe_sol_eclipse_when_glob`).

`eclipse_where` specifically builds (swecl.c:667–669):
```c
iflag  = SEFLG_SPEED | SEFLG_EQUATORIAL | ifl;   // polar, degrees
iflag2 = iflag | SEFLG_RADIANS;                  // polar, radians  (note: NOT reusing the XYZ var name)
iflag  = iflag | SEFLG_XYZ;                       // iflag now becomes cartesian!
```
**Read carefully**: `iflag` is reassigned in place to become the cartesian flavor after
`iflag2` is derived from its pre-XYZ value. This is a real hazard when porting: `iflag2` is
polar+radians, `iflag` ends up cartesian (never combined with `SEFLG_RADIANS`, since cartesian
coordinates have no "degrees vs radians" distinction). Two calls per body then follow — one with
each flavor — because the code wants **both** the cartesian vector (for shadow-cone geometry)
**and** the polar radius/angle in radians (for sidereal-time and distance bookkeeping) from a
single logical position.

`SEFLG_NONUT` (no-nutation) is tested at swecl.c:690 to pick `swe_sidtime0` (mean sidereal time,
explicit obliquity argument) vs `swe_sidtime` (apparent sidereal time) — this only matters when
the caller explicitly suppresses nutation.

---

## 3. `eclipse_where` — shadow geometry core (swecl.c:640–886)

```c
static int32 eclipse_where(double tjd_ut, int32 ipl, char *starname, int32 ifl,
                            double *geopos, double *dcore, char *serr)
```
Computes, for the geocentric instant `tjd_ut` (UT), the shadow-cone geometry of `ipl`
(occulting: always the Moon) with respect to the eclipsed body `ipl`/`starname` (Sun for solar
eclipses), and — if the shadow axis touches the Earth — the geographic longitude/latitude of the
point of greatest eclipse.

**Generality note**: this function is shared between solar eclipses (`ipl=SE_SUN`,
`starname=NULL`) and lunar occultations of any planet/asteroid/star. "Sun" below always means
"the eclipsed/occulted body identified by `ipl`/`starname`."

### 3.1 Inputs & setup (swecl.c:663–704)
- `dcore[0..9]` zeroed.
- `deltat = swe_deltat_ex(tjd_ut, ifl, serr)`; `tjd = tjd_ut + deltat` (ET).
- Moon: cartesian equatorial `rm` (`iflag`, AU) and polar-radians equatorial `lm` (`iflag2`) via
  `swe_calc`.
- Eclipsed body: cartesian `rs` and polar-radians `ls` via `calc_planet_star` (§1).
- `rst`/`rmt` = saved copies of the raw cartesian sun/moon vectors (used later, unmodified by the
  earth-oblateness correction, for the final core-shadow-diameter calc at swecl.c:867–875).
- Sidereal time `sidt` (radians): `swe_sidtime0(...)*15*DEGTORAD` if `SEFLG_NONUT`, else
  `swe_sidtime(tjd_ut)*15*DEGTORAD`.
- `drad` = physical radius of the eclipsed body in AU (§0 lookup pattern; `0` for stars).

### 3.2 Earth-oblateness substitution (swecl.c:705–787, label `iter_where`)
Rather than modeling the Earth as an ellipsoid directly, the code **stretches the z-coordinate**
of the Sun/Moon position vectors so that the Earth can be treated as a **sphere** of radius
`de = 6378140.0 / AUNIT` (equatorial radius) in this rescaled frame:
```c
double earthobl = 1 - EARTH_OBLATENESS;   // first pass: 1 - f
...
swi_polcart(lx, rm);      // rm rebuilt from polar lm (NOT from the cartesian rm computed earlier!)
rm[2] /= earthobl;        // stretch z by 1/(1-f)
dm = sqrt(square_sum(rm));           // distance of moon from geocenter, oblateness-adjusted
swi_polcart(lx, rs);      // same for sun
rs[2] /= earthobl;
```
**FP-fidelity hazard**: `rm`/`rs` cartesian vectors used from here on are **rebuilt from the
polar (`lm`/`ls`) representation** via `swi_polcart`, discarding the directly-computed cartesian
values from `swe_calc(..., iflag, ...)`. This round-trips through degrees→radians→sin/cos and
will not bit-match a naive "just use the cartesian swe_calc output" port.

Sun−Moon geometry:
```c
e[i]  = rm[i] - rs[i];         et[i] = rmt[i] - rst[i];   // moon-minus-sun vector (oblateness-adj. and raw)
dsm   = |e|;                    dsmt  = |et|;
e[i] /= dsm;  et[i] /= dsmt;     // unit vectors
```
Umbra/penumbra half-angles (cone geometry — `sinf1`/`cosf1` = umbra half-angle, tangent cone
converging behind the smaller disc; `sinf2`/`cosf2` = penumbra half-angle, diverging cone):
```c
sinf1 = (drad - rmoon) / dsm;   cosf1 = sqrt(1 - sinf1*sinf1);
sinf2 = (drad + rmoon) / dsm;   cosf2 = sqrt(1 - sinf2*sinf2);
```
where `rmoon = RMOON` (constant, §0) and `dmoon = 2*rmoon`.

Fundamental-plane quantities:
```c
s0 = -dot_prod(rm, e);                     // distance of moon from fundamental plane (along shadow axis)
r0 = sqrt(dm*dm - s0*s0);                   // distance of shadow AXIS from geocenter
d0 = (s0/dsm*(drad*2 - dmoon) - dmoon) / cosf1;   // core-shadow (umbra) diameter on fundamental plane
D0 = (s0/dsm*(drad*2 + dmoon) + dmoon) / cosf2;   // half-shadow (penumbra) diameter on fundamental plane
```
Stored into `dcore[2]=r0`, `dcore[3]=d0`, `dcore[4]=D0`, `dcore[5]=cosf1`, `dcore[6]=cosf2`; then
`dcore[2..4]` (only, **not** `[5]`/`[6]`) are scaled `*= AUNIT/1000.0` (AU → km). **The public
doc comment above `swe_sol_eclipse_where` (swecl.c:635–639) only mentions `dcore[0]`, `[2]`,
`[3]`, `[4]` — `dcore[1]`, `[5]`, `[6]` exist in code but are undocumented in the comment; use
the code, not the comment, as ground truth.**

### 3.3 Central / noncentral / partial classification (swecl.c:764–781)
```c
retc = 0;
if      (de*cosf1 >= r0)                     retc |= SE_ECL_CENTRAL;
else if (r0 <= de*cosf1 + fabs(d0)/2)        retc |= SE_ECL_NONCENTRAL;
else if (r0 <= de*cosf2 + D0/2)              retc |= (SE_ECL_PARTIAL | SE_ECL_NONCENTRAL);
else { /* no eclipse */ retc = 0; no_eclipse = TRUE; geopos[0..1] = 0; *dcore = 0; d = 0; }
```
`de*cosf1` / `de*cosf2` are the Earth's radius projected along the umbra/penumbra cone's tangent
direction. **Central**: shadow axis passes through the Earth (there is a literal center line).
**Noncentral** (non-partial): the umbra/antumbra grazes the Earth off-axis — a full total/annular
eclipse is visible somewhere, but there is no center line. **Partial+Noncentral**: only the
penumbra reaches the Earth — nowhere is the eclipse total/annular. Else: shadow entirely misses
the Earth.

**Important**: the `else` branch does **not** early-return (the `return retc;` is commented out,
swecl.c:780). Execution falls through into §3.4 regardless of `no_eclipse`, so `geopos` is always
populated with the point of *closest approach* of the shadow axis to the Earth's surface even
when no real eclipse occurs there (useful as "if there had been an eclipse, this is where it
would be centered").

### 3.4 Geographic position of eclipse center (swecl.c:782–864)
```c
d = s0*s0 + de*de - dm*dm;   d = (d > 0) ? sqrt(d) : 0;   // fundamental-plane intersection (quadratic root)
s = s0 - d;                                                 // distance moon -> shadow point on earth
xs[i] = rm[i] + s*e[i];                                     // geocentric cartesian shadow point
```
This `d`/`s` pair is literally the positive root of the quadratic for the intersection of the
shadow axis (a ray from the Moon along unit vector `e`) with the sphere of radius `de` — `d` is
the half-chord length, `s0-d` backs up from the moon's fundamental-plane foot-point to the near
intersection with the Earth's surface (or, if `d==0`, to the point of closest approach when the
axis misses the sphere entirely).

`xst` = `xs` with `z *= earthobl` (un-does the earlier z-stretch to get back to true oblate-earth
z), converted to polar via `swi_cartpol`.

**Two-pass oblateness refinement** (swecl.c:842–851, guarded by `if (niter <= 0) { ...; niter++;
goto iter_where; }` — executes the whole §3.2–3.4 block **exactly twice**, not a convergence
loop):
```c
cosfi = cos(xst[1]); sinfi = sin(xst[1]);           // latitude found in pass 1
eobl = EARTH_OBLATENESS;
cc = 1 / sqrt(cosfi*cosfi + (1-eobl)*(1-eobl)*sinfi*sinfi);
ss = (1-eobl)*(1-eobl) * cc;
earthobl = ss;                                       // replaces the crude (1-f) with an ellipsoid-normal factor at this latitude
```
Pass 2 re-runs the entire geometry with this locally-linearized `earthobl`, giving a more exact
ellipsoid correction than the uniform-sphere approximation of pass 1.

Final geographic position (swecl.c:852–864):
```c
swi_cartpol(xs, xs);           // to polar (lon, lat, r), radians
xs[0] -= sidt;                  // subtract Greenwich sidereal time -> geographic longitude
xs[0] *= RADTODEG; xs[1] *= RADTODEG;
xs[0] = swe_degnorm(xs[0]); if (xs[0] > 180) xs[0] -= 360;   // west negative convention
geopos[0] = xs[0];   // central longitude, east positive
geopos[1] = xs[1];   // central latitude, north positive
```

### 3.5 Core-shadow diameter at the point of maximum eclipse (swecl.c:865–875)
Uses the **raw** (non-oblateness-adjusted) saved vectors `rmt` and `xst` (true-z shadow point):
```c
x[i] = rmt[i] - xst[i];   s = |x|;                   // moon -> shadow-point distance (raw)
*dcore    = (s/dsmt*(drad*2 - dmoon) - dmoon) * cosf1;   *dcore    *= AUNIT/1000.0;  // core (umbra) diam, km, SIGNED
dcore[1]  = (s/dsmt*(drad*2 + dmoon) + dmoon) * cosf2;   dcore[1]  *= AUNIT/1000.0;  // penumbra diam, km
```
`*dcore` (`dcore[0]`) carries a **sign**: positive → annular (umbra apex is above the ground,
antumbral shadow only), negative/zero → total (umbra cone reaches past the ground). Final
type refinement (swecl.c:876–884):
```c
if (!(retc & SE_ECL_PARTIAL) && !no_eclipse) {
  if (*dcore > 0) retc |= SE_ECL_ANNULAR;
  else            retc |= SE_ECL_TOTAL;
}
```

### 3.6 `dcore[]` full index summary
| Index | Meaning | Units |
|---|---|---|
| 0 | Core (umbra) shadow diameter at point of max eclipse (signed: >0 annular, <0 total) | km |
| 1 | Penumbra diameter at point of max eclipse | km |
| 2 | `r0`: distance of shadow axis from geocenter | km |
| 3 | `d0`: umbra diameter on fundamental plane | km |
| 4 | `D0`: penumbra diameter on fundamental plane | km |
| 5 | `cosf1`: cosine of umbra half-angle | dimensionless (NOT unit-scaled) |
| 6 | `cosf2`: cosine of penumbra half-angle | dimensionless (NOT unit-scaled) |
| 7–9 | unused, always 0 | — |

### 3.7 Return value
`ERR` on `swe_calc`/`calc_planet_star` failure; else the classification bitmask from §3.3/§3.5
(`0` = no eclipse anywhere on Earth at this instant).

### 3.8 `swe_sol_eclipse_where` — public wrapper (swecl.c:565–582)
```c
int32 CALL_CONV swe_sol_eclipse_where(double tjd_ut, int32 ifl, double *geopos, double *attr, char *serr)
```
1. `ifl &= SEFLG_EPHMASK`; `swi_set_tid_acc(tjd_ut, ifl, 0, serr)` — STATELESS PORT NOTE: mutates
   global tidal-acceleration state used by subsequent internal deltaT/ephemeris calls; the
   stateless Rust port should thread the tidal-acceleration mode explicitly instead.
2. `eclipse_where(tjd_ut, SE_SUN, NULL, ifl, geopos, dcore, serr)` → `retflag`; return early if `ERR`.
3. `eclipse_how(tjd_ut, SE_SUN, NULL, ifl, geopos[0], geopos[1], 0, attr, serr)` (§4) — local
   circumstances **at the just-computed central point**, height fixed at 0 m.
4. `attr[3] = dcore[0]` (overwrite eclipse_how's own `attr[3]`, which is unset/0, with the true
   core-shadow diameter from the geometry pass).
5. Return `retflag` (from step 2 — i.e. the CENTRAL/NONCENTRAL/TOTAL/ANNULAR/PARTIAL classification,
   not anything from `eclipse_how`).

`geopos[2..9]` (northern/southern umbra/penumbra limit lon/lat) are **documented in the public
header comment (swecl.c:527–536) as "not implemented so far"** — always left at whatever the
caller pre-set them to; `eclipse_where` only ever writes `geopos[0]`/`[1]`.

---

## 4. `eclipse_how` — local circumstances (swecl.c:967–1152)

```c
static int32 eclipse_how(double tjd_ut, int32 ipl, char *starname, int32 ifl,
                          double geolon, double geolat, double geohgt,
                          double *attr, char *serr)
```
Computes the eclipse as seen from a specific observer (topocentric).

### 4.1 Setup (swecl.c:986–1011)
```c
attr[0..9] = 0;
geopos = {geolon, geolat, geohgt};
te = tjd_ut + swe_deltat_ex(tjd_ut, ifl, serr);        // ET
swe_set_topo(geolon, geolat, geohgt);                   // STATELESS PORT NOTE: global swed.topd
iflag     = SEFLG_EQUATORIAL | SEFLG_TOPOCTR | ifl;
iflagcart = iflag | SEFLG_XYZ;
ls = calc_planet_star(te, ipl, starname, iflag, ...);      // topocentric polar (deg), eclipsed body
lm = swe_calc(te, SE_MOON, iflag, ...);                    // topocentric polar (deg), Moon
xs = calc_planet_star(te, ipl, starname, iflagcart, ...);  // topocentric cartesian (AU)
xm = swe_calc(te, SE_MOON, iflagcart, ...);
drad = <body radius in AU, §0 lookup, 0 for stars>;
```
`swe_set_topo` here is the **STATELESS PORT NOTE** global-state touch point: C relies on
`swed.topd` being populated before any `SEFLG_TOPOCTR` calc call. The Rust port passes the
observer position explicitly as a parameter into the position calculation instead.

### 4.2 Azimuth/altitude (swecl.c:1012–1029)
```c
swe_azalt(tjd_ut, SE_EQU2HOR, geopos, /*atpress=*/0, /*attemp=*/10, ls, xh);
```
`atpress=0` tells `swe_azalt` to auto-estimate atmospheric pressure from `geopos[2]` (barometric
formula `1013.25 * (1 - 0.0065*h/288)^5.255`); `attemp=10` is a fixed 10 °C. Output convention
(from `swe_azalt`, swecl.c:2788ff): `xh[0]` = azimuth, **measured from south, clockwise via west**;
`xh[1]` = true (geometric) altitude; `xh[2]` = apparent altitude (refraction-corrected via
`swe_refrac_extended`). A dead `#if USE_AZ_NAV` branch (always `0`) computing azimuth
"from north, clockwise via east" exists but is compiled out — ignore it.

### 4.3 Angular radii and center separation (swecl.c:1030–1039)
```c
rmoon = asin(RMOON / lm[2]) * RADTODEG;    // lm[2] = topocentric distance to Moon, AU
rsun  = asin(drad  / ls[2]) * RADTODEG;    // ls[2] = topocentric distance to eclipsed body, AU
rsplusrm  = rsun + rmoon;
rsminusrm = rsun - rmoon;
x1[i] = xs[i]/ls[2];  x2[i] = xm[i]/lm[2];         // manually pre-normalized (see §0 FP hazard)
dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;  // angular separation of centers, topocentric
```

### 4.4 Phase classification (swecl.c:1040–1053)
```c
if      (dctr < rsminusrm)        retc = SE_ECL_ANNULAR;   // only reachable when rsun > rmoon
else if (dctr < fabs(rsminusrm))  retc = SE_ECL_TOTAL;      // reachable when rmoon > rsun (real solar totality)
else if (dctr < rsplusrm)         retc = SE_ECL_PARTIAL;
else                               retc = 0;                 // no eclipse; serr set
```
Uses the **instantaneous** angular radii (which vary with Earth–Moon / Earth–Sun distance), so
the same location can see an annular eclipse near lunar apogee and a total eclipse near lunar
perigee depending purely on `rmoon` vs `rsun` at the time.

### 4.5 Magnitude (`attr[0]`), diameter ratio (`attr[1]`) (swecl.c:1055–1074)
```c
attr[1] = (rsun > 0) ? rmoon/rsun : 0;              // ratio of lunar to solar angular diameter
lsunleft = -dctr + rsun + rmoon;                     // overlap depth (0 at first/last contact)
attr[0] = (rsun > 0) ? lsunleft/rsun/2 : 1;          // fraction of solar diameter covered (IMCCE magnitude)
```
(The literal C computes a `lsun = asin(rsun/2*DEGTORAD)*2` value purely to test its sign as a
proxy for `rsun>0`; algebraically it is monotonic with `rsun` and equivalent to testing
`rsun > 0` directly — port the simpler test, but note the literal expression if bit-exact replay
of this dead computation is ever needed.)

### 4.6 Obscuration (`attr[2]`) — circular-segment lens-area formula (swecl.c:1075–1107)
Fraction of the eclipsed body's **disc area** (not diameter) covered by the Moon:
```c
if (retc == 0 || lsun == 0) {
  attr[2] = 1;                              // sentinel; historically "100", changed to 1 (fraction convention)
} else if (retc == SE_ECL_TOTAL || retc == SE_ECL_ANNULAR) {
  attr[2] = lmoon*lmoon / lsun / lsun;       // one disc entirely inside the other: area ratio
} else {                                     // partial: intersecting-circles lens area
  a = 2*lctr*lmoon;  b = 2*lctr*lsun;
  if (a < 1e-9) {
    attr[2] = lmoon*lmoon / lsun / lsun;     // guard: centers nearly coincident
  } else {
    a = clamp((lctr² + lmoon² - lsun²) / a, -1, 1);   a = acos(a);   // half-angle subtended at Moon's center
    b = clamp((lctr² + lsun²  - lmoon²) / b, -1, 1);   b = acos(b);  // half-angle subtended at Sun's center
    sc1 = (a - sin(a)*cos(a)) * lmoon²/2;    // circular segment area (Moon disc), = (r²/2)(θ - sinθ), θ=2a
    sc2 = (b - sin(b)*cos(b)) * lsun²/2;     // circular segment area (Sun disc),  θ=2b
    attr[2] = (sc1 + sc2) * 2 / PI / lsun²;  // total lens area (2×segments) / full solar-disc area
  }
}
```
where `lsun=rsun`, `lmoon=rmoon`, `lctr=dctr` (all in degrees, but consistently so — the formula
is scale-invariant since it's a ratio of areas). The `*2` factor is required because `sc1`/`sc2`
as coded are each **half** of the standard segment-area formula `(r²/2)(θ - sinθ)` with `θ=2a`
(code computes `(r²/2)(a - sin(a)cos(a)) = (r²/4)(θ - sinθ)`); the final `*2` restores the full
lens area = sum of the two full circular segments.

`attr[7] = dctr` (angular separation of centers, degrees — "elongation").

### 4.7 Visibility threshold (swecl.c:1108–1123)
```c
hmin_appr = -(34.4556 + (1.75 + 0.37) * sqrt(geohgt)) / 60;   // degrees; geohgt in meters
if (xh[1] + rsun + fabs(hmin_appr) >= 0 && retc)
  retc |= SE_ECL_VISIBLE;
attr[4] = xh[0];   // azimuth (south-based, clockwise via west)
attr[5] = xh[1];   // true altitude
attr[6] = xh[2];   // apparent (refracted) altitude
```
`34.4556` arcmin = standard refraction at the horizon (Bennett's formula); `1.75'/√h` = dip of the
horizon for an observer at height `h` meters; `0.37'/√h` = extra refraction between horizon and
observer. The visibility test uses the **true** altitude `xh[1]` (not apparent), padded by the
body's own angular radius and the horizon-depression allowance — i.e. "any part of the disc
could geometrically be above the (refraction- and dip-lowered) horizon."

### 4.8 NASA magnitude and Saros number (swecl.c:1124–1150)
Only computed **if `ipl==SE_SUN` and `starname` empty** (a genuine solar eclipse, not an
occultation reusing this same function):
```c
attr[8] = attr[0];                                    // fraction covered, by default
if (retc & (SE_ECL_TOTAL | SE_ECL_ANNULAR)) attr[8] = attr[1];   // diameter ratio, for total/annular
```
Saros lookup: for each of the 181 entries in `saros_data_solar[]`,
`d = (tjd_ut - series.tstart) / SAROS_CYCLE`; find the integer member `j` (or `j+1`) such that
`tjd_ut` is within **2 days worth of `1/SAROS_CYCLE`** of an exact multiple — practically, within
±2 days of a Saros-cycle boundary for that series. On no match in any series after scanning all
181: `attr[9] = attr[10] = -99999999` (sentinel "no Saros data", e.g. for eclipses outside
[-2955, +2669] roughly). Otherwise `attr[9] = series_no`, `attr[10] = member number` (1-based).

### 4.9 `attr[]` full index summary (solar eclipse / occultation local circumstances)
| Index | Meaning | Units |
|---|---|---|
| 0 | Magnitude: fraction of eclipsed body's diameter covered by Moon (IMCCE convention) | fraction |
| 1 | Ratio of lunar angular diameter to eclipsed body's angular diameter | dimensionless |
| 2 | Obscuration: fraction of eclipsed body's disc **area** covered | fraction |
| 3 | Core (umbra) shadow diameter — filled by the **caller** from `dcore[0]`, not by `eclipse_how` itself | km |
| 4 | Azimuth of eclipsed body, measured from south, clockwise via west | degrees |
| 5 | True altitude of eclipsed body above horizon | degrees |
| 6 | Apparent (refracted) altitude of eclipsed body above horizon | degrees |
| 7 | Angular separation ("elongation") of Moon from eclipsed body's center | degrees |
| 8 | Magnitude per NASA convention (= attr[0] for partial; = attr[1] for total/annular) | fraction / ratio |
| 9 | Saros series number (solar eclipses of the Sun only; `-99999999` if none found) | integer-valued double |
| 10 | Saros series member number, 1-based (solar eclipses of the Sun only) | integer-valued double |

Caller must allocate `attr[20]` minimum (C comment convention) even though only 0–10 are used
here (lunar-eclipse variants use a few more indices, out of scope for this doc).

### 4.10 Return value
`ERR` on ephemeris failure; else `SE_ECL_TOTAL`/`ANNULAR`/`PARTIAL` possibly OR'd with
`SE_ECL_VISIBLE`; `0` if no eclipse visible from this location at this instant.

### 4.11 `swe_sol_eclipse_how` — public wrapper (swecl.c:922–964)
```c
int32 CALL_CONV swe_sol_eclipse_how(double tjd_ut, int32 ifl, double *geopos, double *attr, char *serr)
```
1. Validate `geopos[2]` (height) ∈ `[SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX]` else `ERR` with
   message `"location for eclipses must be between %.0f and %.0f m above sea"`.
2. `ifl &= SEFLG_EPHMASK`; `swi_set_tid_acc(...)`.
3. `retflag = eclipse_how(tjd_ut, SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr)`.
4. `retflag2 = eclipse_where(tjd_ut, SE_SUN, NULL, ifl, geopos2, dcore, serr)` — **into a
   scratch `geopos2[20]`**, purely to obtain `dcore` and the CENTRAL/NONCENTRAL classification
   (this call is geocentric and independent of the observer's `geopos`; it always returns the
   same central-line location for a given `tjd_ut` no matter where `swe_sol_eclipse_how` is
   evaluated from).
5. `if (retflag) retflag |= (retflag2 & (SE_ECL_CENTRAL|SE_ECL_NONCENTRAL))` — only tags
   central/noncentral if the LOCAL observer actually sees an eclipse.
6. `attr[3] = dcore[0]` — **this is the geocentric core-shadow diameter at the moment of maximum
   for the whole Earth, not anything specific to the observer's location or to `tjd_ut`'s local
   circumstances.**
7. **Redundant az/alt recomputation** (a second, independent pass over what `eclipse_how` already
   computed into `attr[4..6]`): `swe_set_topo(...)`; `swe_calc_ut(tjd_ut, SE_SUN, ifl|TOPOCTR|EQUATORIAL, ls, ...)`;
   `swe_azalt(tjd_ut, SE_EQU2HOR, geopos, 0, 10, ls, xaz)`; overwrite `attr[4]=xaz[0]`,
   `attr[5]=xaz[1]`, `attr[6]=xaz[2]`. This uses `swe_calc_ut` (UT + internal deltaT) rather than
   `eclipse_how`'s `calc_planet_star(te, ...)` (pre-converted ET) — numerically equivalent but a
   literal second calc call; port faithfully (do not dedupe away — see step 8's dependency on
   `xaz[2]`).
8. **Visibility gate that can zero out the whole result**:
   ```c
   if (xaz[2] <= 0) retflag = 0;                 // apparent altitude <= 0 -> not visible
   if (retflag == 0) {
     for (i = 0; i <= 3; i++)  attr[i] = 0;        // magnitude/ratio/obscuration/core-diam cleared
     for (i = 8; i <= 10; i++) attr[i] = 0;        // NASA magnitude/saros cleared
     // attr[4..7] (az/alt/elongation) are LEFT populated even when retflag becomes 0
   }
   return retflag;
   ```
   **This means `swe_sol_eclipse_how` can return `0` (no eclipse) purely because the Sun is below
   the horizon, even though `eclipse_how`'s own geometric computation found a real eclipse in
   progress.** The Rust port must replicate this: it is not simply "return `eclipse_how`'s
   result," there is a horizon-visibility override layered on top by the public wrapper.

---

## 5. `swe_sol_eclipse_when_glob` — global eclipse search (swecl.c:1185–1515)

```c
int32 CALL_CONV swe_sol_eclipse_when_glob(double tjd_start, int32 ifl, int32 ifltype,
                                           double *tret, int32 backward, char *serr)
```
Finds the next (or, if `backward`, previous) solar eclipse anywhere on Earth after/before
`tjd_start` (UT), restricted to eclipse types in `ifltype`.

### 5.1 `ifltype` validation & normalization (swecl.c:1210–1226)
- Reject `ifltype == (SE_ECL_PARTIAL|SE_ECL_CENTRAL)` — central partial eclipses cannot exist
  (`ERR`, `"central partial eclipses do not exist"`).
- Reject `ifltype == (SE_ECL_ANNULAR_TOTAL|SE_ECL_NONCENTRAL)` — noncentral hybrids cannot exist
  (`ERR`, message similarly).
- `ifltype == 0` → all types: `TOTAL|ANNULAR|PARTIAL|ANNULAR_TOTAL|NONCENTRAL|CENTRAL`.
- Bare `TOTAL`/`ANNULAR`/`ANNULAR_TOTAL` (no central/noncentral qualifier given) → OR in both
  `NONCENTRAL|CENTRAL` (accept either variant of that type).
- Bare `PARTIAL` → OR in `NONCENTRAL` (partial eclipses are inherently noncentral).

### 5.2 Lunation stepping (Meeus) (swecl.c:1227–1265)
```c
direction = backward ? -1 : 1;
K = (int)((tjd_start - J2000) / 365.2425 * 12.3685) - direction;   // synodic-month index estimate
next_try:
  T = K/1236.85;  T2=T*T; T3=T2*T; T4=T3*T;
  Ff = swe_degnorm(160.7108 + 390.67050274*K - 0.0016341*T2 - 0.00000227*T3 + 0.000000011*T4);
  if (Ff > 180) Ff -= 180;
  if (Ff > 21 && Ff < 159) { K += direction; goto next_try; }    // F-argument filter, see below
  tjd = 2451550.09765 + 29.530588853*K + 0.0001337*T2 - 0.000000150*T3 + 0.00000000073*T4;
  M  = swe_degnorm(2.5534 + 29.10535669*K - 0.0000218*T2 - 0.00000011*T3);     // Sun's mean anomaly
  Mm = swe_degnorm(201.5643 + 385.81693528*K + 0.1017438*T2 + 0.00001239*T3 + 0.000000058*T4); // Moon's mean anomaly
  E  = 1 - 0.002516*T - 0.0000074*T2;                                          // Earth-orbit eccentricity correction
  tjd = tjd - 0.4075*sin(Mm) + 0.1721*E*sin(M);                                // periodic correction
```
`12.3685` = synodic months per Julian year (`365.2425/29.530588853`); `K` is the (signed) synodic
month count since the mean new moon nearest J2000. `Ff` is **Meeus's F argument** (Moon's mean
argument of latitude, i.e. distance from ascending node) evaluated **at the mean conjunction**,
folded into `[0,180)`. **The `[21°, 159°]` exclusion band eliminates lunations where the Moon is
too far from a node for any solar eclipse to be geometrically possible** — this is the primary
efficiency filter, discarding roughly 70% of all new moons without any `swe_calc` call.

`K -= direction` before the loop (so the very first candidate tested is the one nearest
`tjd_start` itself, not one lunation past it); every rejection does `K += direction; goto
next_try;` to advance in the requested search direction.

### 5.3 Iterative refinement to instant of minimum separation (swecl.c:1273–1300)
```c
dtstart = (tjd < 2000000 || tjd > 2500000) ? 5 : 1;   // wider initial window far from "modern" JD range
dtdiv = 4;
for (dt = dtstart; dt > 0.0001; dt /= dtdiv) {
  for (i=0, t=tjd-dt; i<=2; i++, t+=dt) {
    xs,ls = swe_calc(t, SE_SUN, iflagcart/iflag);   xm,lm = swe_calc(t, SE_MOON, ...);
    dc[i] = acos(swi_dot_prod_unit(xs/ls[2], xm/lm[2])) * RADTODEG
            - (rmoon + rsun);       // gap between limbs; negative when overlapping
  }
  find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);   // parabola vertex (most negative = deepest overlap)
  tjd += dtint + dt;
}
```
Note: `t` passed straight to `swe_calc` — the Meeus `tjd` is being treated as an ET/TT instant
throughout this refinement (consistent with Meeus's formula, which is dynamical time), **not**
UT; UT conversion happens once, after convergence (`dt < 0.0001` day ≈ 8.6 s):
```c
tjds = tjd - swe_deltat_ex(tjd, ifl, serr);
tjds = tjd - swe_deltat_ex(tjds, ifl, serr);
tjds = tjd = tjd - swe_deltat_ex(tjds, ifl, serr);      // 3-pass fixed-point ET->UT conversion
```

### 5.4 Confirm & classify (swecl.c:1304–1360)
```c
eclipse_where(tjd, SE_SUN, NULL, ifl, geopos, dcore, serr) -> retflag;
eclipse_how(tjd, SE_SUN, NULL, ifl, geopos[0], geopos[1], 0, attr, serr) -> retflag2;
if (retflag2 == 0) { K += direction; goto next_try; }    // confirm via _how() in case _where() under-detects a tiny eclipse
tret[0] = tjd;
if ((backward && tret[0] >= tjd_start-0.0001) || (!backward && tret[0] <= tjd_start+0.0001))
  { K += direction; goto next_try; }                      // reject candidates not strictly beyond tjd_start
eclipse_where(tjd, ...) -> retflag;                        // re-derive type bits (TOTAL/ANNULAR/PARTIAL/CENTRAL/NONCENTRAL)
if (retflag == 0) { retflag = SE_ECL_PARTIAL|SE_ECL_NONCENTRAL; tret[4]=tret[5]=tjd; dont_times = TRUE; }  // FIXME in upstream C
```
Then a cascade of `ifltype`-bit rejections (advance `K` and retry if the found type isn't
requested): NONCENTRAL, CENTRAL, ANNULAR, PARTIAL, and (provisionally) TOTAL — `ANNULAR_TOTAL` is
only discovered later (§5.6), so a `TOTAL` result is allowed through here if `ANNULAR_TOTAL` was
also requested.

If `dont_times` (the "no eclipse but forced-partial fallback" case above): skip straight to
returning `retflag` without computing any contact times.

### 5.5 Contact-time refinement (swecl.c:1361–1418)
Which contact pairs to compute depends on the found type:
```c
o = (retflag & SE_ECL_PARTIAL) ? 0 : (retflag & SE_ECL_NONCENTRAL) ? 1 : 2;
dta = twohr (2/24 day);  dtb = tenmin/3 (≈3.33 min);
for (n = 0; n <= o; n++) {
  n==0: i1,i2 = 2,3   // eclipse begin/end (anywhere on Earth) — always computed
  n==1: i1,i2 = 4,5   // totality/annularity begin/end — skipped if PARTIAL
  n==2: i1,i2 = 6,7   // center-line begin/end        — skipped if NONCENTRAL
  sample dc[i] at t = tjd-dta, tjd, tjd+dta via eclipse_where(t,...):
    n==0: dc[i] = dcore[4]/2 + de/dcore[5] - dcore[2]     // NOTE: uses cosf1 (dcore[5]), not cosf2 — see hazard below
    n==1: dc[i] = fabs(dcore[3])/2 + de/dcore[6] - dcore[2]
    n==2: dc[i] = de/dcore[6] - dcore[2]
  find_zero(dc[0],dc[1],dc[2],dta,&dt1,&dt2);
  tret[i1] = tjd+dt1+dta;  tret[i2] = tjd+dt2+dta;
  // Newton/secant refinement, 3 passes, dt = dtb, dtb/3, dtb/9:
  for (m=0, dt=dtb; m<3; m++, dt/=3)
    for (j = i1; j <= i2; j += (i2-i1)) {
      sample dc[0] at t=tret[j]-dt, dc[1] at t=tret[j]  (same formula as above, re-eval eclipse_where)
      dt1 = dc[1] / ((dc[1]-dc[0])/dt);   tret[j] -= dt1;
    }
}
```
where `de = 6378.140` km (equatorial radius, **km** here, not AU — different unit convention than
`eclipse_where`'s internal `de`).

**FP-fidelity / literal-quirk hazard**: the `n==0` (eclipse begin/end, driven by the *penumbra*
boundary `D0=dcore[4]`) formula divides by `dcore[5]` (`cosf1`, the **umbra** half-angle cosine),
not `dcore[6]` (`cosf2`, penumbra). Since both half-angles are well under 1°, `cosf1≈cosf2≈1` and
the numerical impact is negligible, but this is what the C literally does — port it exactly as
written (do not "fix" it to use `cosf2`), since the golden tests compare against this exact code
path.

### 5.6 Annular-total (hybrid) detection (swecl.c:1419–1450)
```c
if (retflag & SE_ECL_TOTAL) {
  dc[0] = dcore[0] at tret[0] (max);   dc[1] = dcore[0] at tret[4] (totality begin);   dc[2] = dcore[0] at tret[5] (totality end);
  if (dc[0]*dc[1] < 0 || dc[0]*dc[2] < 0) {
    retflag |= SE_ECL_ANNULAR_TOTAL;   retflag &= ~SE_ECL_TOTAL;
  }
}
```
Recall `dcore[0]` is signed (>0 annular, <0 total, §3.5); a sign change between the maximum and
either edge of totality means the eclipse is annular at one end and total at the other — a hybrid
(annular-total) eclipse. Followed by rejection/retry if `TOTAL`/`ANNULAR_TOTAL` found-but-not-
wanted per `ifltype`.

### 5.7 `tret[1]` — time at "local apparent noon" / RA-conjunction instant (swecl.c:1451–1498)
First checks for a solar-transit sign change: geocentric equatorial `ls[0]`/`lm[0]` (right
ascension) difference, folded to `[-180,180]`, evaluated at `tret[2]` and `tret[3]` (eclipse
begin/end, converted to ET). If no sign change (`dc[0]*dc[1] >= 0`): `tret[1] = 0` (no such
instant within the eclipse window). Otherwise, secant-iterate from `tjds` (the UT max time)
toward the instant where geocentric `RA(Sun) == RA(Moon)` exactly, shrinking `dt` by `/3` each
pass starting from `dt = min(0.1, (tret[3]-tret[2])/4)`, stopping once `dt <= 0.01`.

### 5.8 `tret[]` full index summary (global search)
| Index | Meaning |
|---|---|
| 0 | Time (UT) of maximum eclipse (geocentric minimum Sun–Moon angular separation) |
| 1 | Time when the eclipse's RA-conjunction instant occurs (0 if none within the eclipse window) |
| 2 | Time of eclipse begin (first contact, anywhere on Earth) |
| 3 | Time of eclipse end (last contact, anywhere on Earth) |
| 4 | Time of totality/annularity begin (0 if partial) |
| 5 | Time of totality/annularity end (0 if partial) |
| 6 | Time of center-line begin (0 if noncentral) |
| 7 | Time of center-line end (0 if noncentral) |
| 8 | (annular-total transition to total — **not implemented**, always 0) |
| 9 | (annular-total transition back to annular — **not implemented**, always 0) |

Caller must allocate `tret[10]` minimum.

### 5.9 Return value
`ERR` on ephemeris failure; else the classification bitmask (`CENTRAL`/`NONCENTRAL` combined with
`TOTAL`/`ANNULAR`/`ANNULAR_TOTAL`/`PARTIAL`) of the found eclipse. Never returns `0` (the search
loop retries indefinitely — bounded only by ephemeris range — until a matching eclipse is found).

---

## 6. `swe_sol_eclipse_when_loc` + `eclipse_when_loc` — local eclipse search

### 6.1 `swe_sol_eclipse_when_loc` — public wrapper (swecl.c:2019–2041)
```c
int32 CALL_CONV swe_sol_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos,
                                          double *tret, double *attr, int32 backward, char *serr)
```
1. Validate `geopos[2]` range (same as §4.11 step 1).
2. `ifl &= SEFLG_EPHMASK`; `swi_set_tid_acc(...)`.
3. `retflag = eclipse_when_loc(tjd_start, ifl, geopos, tret, attr, backward, serr)`; return early
   if `retflag <= 0` (covers both `ERR=-1` and the "no eclipse found" `0`, though in practice
   `eclipse_when_loc` retries internally and shouldn't return exactly `0` — see §6.3).
4. `eclipse_where(tret[0], SE_SUN, NULL, ifl, geopos2, dcore, serr)` → merge `SE_ECL_NONCENTRAL`
   bit into `retflag`, and `attr[3] = dcore[0]` (core shadow diameter, same "geocentric, not
   observer-specific" caveat as §4.11 step 6).

### 6.2 `eclipse_when_loc` — worker (swecl.c:2100–2410)

```c
static int32 eclipse_when_loc(double tjd_start, int32 ifl, double *geopos,
                               double *tret, double *attr, int32 backward, char *serr)
```
Finds the next/previous solar eclipse **visible from** `geopos` (topocentric — unlike
`swe_sol_eclipse_when_glob`, which is purely geocentric).

**Note on duplication**: the Meeus lunation-stepping formulas (K estimate, `Ff`/F-argument
filter with the same `[21,159]` threshold, `M`/`Mm`/`E`/periodic-correction `tjd` formula) are
**textually identical** to §5.2 — this function re-derives them independently rather than calling
a shared helper. The Rust port should factor this into one shared function used by both the glob
and loc search entry points.

K/direction handling (swecl.c:2118–2122, spelled with explicit `if`/`else` rather than a
`direction` variable, but equivalent to §5.2's convention):
```c
K = (int)((tjd_start - J2000)/365.2425*12.3685);
if (backward) K++; else K--;                    // pre-loop offset
next_try:
  ... Ff filter: if (backward) K--; else K++; goto next_try;   // in-loop advance (opposite sign of pre-loop offset)
```

`iflag = SEFLG_EQUATORIAL|SEFLG_TOPOCTR|ifl` and `iflagcart = iflag|SEFLG_XYZ` are computed
**once**, at swecl.c:2126–2127, before `swe_set_topo`/`K` are even set up — they are reused
unchanged for the *entire* function (main loop, contacts 2/3, contacts 1/4, everything except the
final `eclipse_how` visibility-scan calls, which build their own topocentric flags internally).

`swe_set_topo(geopos[0], geopos[1], geopos[2])` is called **twice**: once at swecl.c:2128 (before
`K` is even computed, i.e. once per invocation of `eclipse_when_loc`), and again at swecl.c:2173
(immediately after the Meeus periodic-correction `tjd` is finalized, i.e. **inside** the
`next_try:` loop body — every retry re-issues it). Both calls pass identical arguments; idempotent
since `geopos` never changes, but port both call sites faithfully (harmless duplication, not a
bug).

### 6.2.1 Main convergence loop (swecl.c:2174–2204)

```c
dtdiv = 2;
dtstart = 0.5;
if (tjd < 1900000 || tjd > 2500000)   /* because above formula is not good (delta t?) */
  dtstart = 2;
for (dt = dtstart; dt > 0.00001; dt /= dtdiv) {
  if (dt < 0.1)
    dtdiv = 3;
  for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {
    /* this takes some time, but is necessary to avoid missing an eclipse */
    if (swe_calc(t, SE_SUN, iflagcart, xs, serr) == ERR) return ERR;
    if (swe_calc(t, SE_SUN, iflag,     ls, serr) == ERR) return ERR;
    if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR) return ERR;
    if (swe_calc(t, SE_MOON, iflag,     lm, serr) == ERR) return ERR;
    dm = sqrt(square_sum(xm));
    ds = sqrt(square_sum(xs));
    for (k = 0; k < 3; k++) {
      x1[k] = xs[k] / ds;   /*ls[2]*/
      x2[k] = xm[k] / dm;   /*lm[2]*/
    }
    dc[i] = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
  }
  find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
  tjd += dtint + dt;
}
```
`dtstart = 0.5` days normally, `2` days if `tjd` (the Meeus first-guess instant) falls outside
`[1900000, 2500000]` — **a different JD boundary than `swe_sol_eclipse_when_glob`'s
`[2000000,2500000]`** (§5.3); flag this as intentional-but-inconsistent, port per-function, do not
unify. The convergence loop terminates at `dt <= 0.00001` day (~0.86 s) — a tighter threshold than
§5.3's `dt <= 0.0001`. `dtdiv` starts at `2` and switches to `3` **the moment `dt` first drops
below `0.1`** (tested at the *top* of each outer-loop pass, so the very iteration whose `dt < 0.1`
already divides by 3 to produce the *next* `dt`) — a two-stage step-size schedule, unlike §5.3's
constant `/4`.

Sample point stepping is `t = tjd-dt, tjd, tjd+dt` for `i=0,1,2`, and `t` is passed straight to
`swe_calc` — same as §5.3, the Meeus `tjd` is treated as an ET/TT instant throughout this loop
(UT conversion only happens once, after convergence — see below), **not** UT.

**Four `swe_calc` calls per sample** (cartesian+polar × Sun+Moon), but only the cartesian pair
(`xs`, `xm`) is actually used to build `dc[i]` inside the loop — the polar results `ls`/`lm` are
computed and then **discarded** on every iteration of the inner loop (they get overwritten by the
next `i`, and the last-iteration values are never read before being clobbered by the fresh
post-loop calc immediately below). This is dead work inside the loop itself, not a correctness
issue — but do not assume `ls`/`lm` inside the loop carry forward any meaning; only the outer
scope's post-convergence `ls`/`lm` (below) matter.

**FP-fidelity note**: `ds`/`dm` are computed manually via `sqrt(square_sum(xs))`/`sqrt(square_sum(xm))`
rather than reused from `ls[2]`/`lm[2]` (the C source's own inline comments `/*ls[2]*/`/`/*lm[2]*/`
flag this explicitly — the polar distance would be algebraically identical but is not what's coded).
`x1`/`x2` are then pre-normalized by this manually-computed `ds`/`dm` before being passed to
`swi_dot_prod_unit`, which re-normalizes internally — the same redundant-double-division pattern
flagged generically in §0. Port literally.

### 6.2.2 Post-convergence confirmation, ET→UT, rejection (swecl.c:2205–2241)

```c
if (swe_calc(tjd, SE_SUN, iflagcart, xs, serr) == ERR) return ERR;
if (swe_calc(tjd, SE_SUN, iflag,     ls, serr) == ERR) return ERR;
if (swe_calc(tjd, SE_MOON, iflagcart, xm, serr) == ERR) return ERR;
if (swe_calc(tjd, SE_MOON, iflag,     lm, serr) == ERR) return ERR;
dctr = acos(swi_dot_prod_unit(xs, xm)) * RADTODEG;
rmoon = asin(RMOON / lm[2]) * RADTODEG;
rsun  = asin(RSUN  / ls[2]) * RADTODEG;
rsplusrm  = rsun + rmoon;
rsminusrm = rsun - rmoon;
if (dctr > rsplusrm) {
  if (backward) K--; else K++;
  goto next_try;
}
tret[0] = tjd - swe_deltat_ex(tjd, ifl, serr);
tret[0] = tjd - swe_deltat_ex(tret[0], ifl, serr); /* these two lines are an iteration! */
if ((backward && tret[0] >= tjd_start - 0.0001)
  || (!backward && tret[0] <= tjd_start + 0.0001)) {
  if (backward) K--; else K++;
  goto next_try;
}
if (dctr < rsminusrm)
  retflag = SE_ECL_ANNULAR;
else if (dctr < fabs(rsminusrm))
  retflag = SE_ECL_TOTAL;
else if (dctr <= rsplusrm)
  retflag = SE_ECL_PARTIAL;
dctrmin = dctr;
```
This is a **fresh, fifth set** of `swe_calc` calls at the exact converged `tjd` — it does **not**
reuse the loop's last sample (`t = tjd+dt` from the final `i==2` pass), even though that sample is
numerically close. `dctr` here is computed by passing the **raw, non-pre-normalized** `xs`/`xm`
straight into `swi_dot_prod_unit` (unlike the loop's `dc[i]`, which pre-normalizes into `x1`/`x2`
first) — mathematically equivalent (the function normalizes internally regardless) but a distinct
rounding path; port both call shapes exactly as written, do not unify them into one helper that
always pre-normalizes.

This is also the **only** place in the loop/convergence portion where the polar `ls[2]`/`lm[2]`
(topocentric distance) values are actually consumed — via `rmoon`/`rsun`, using **hardcoded
`RSUN`/`RMOON`** constants directly (unlike `eclipse_where`/`eclipse_how`, which resolve body
radius generically via the §0 lookup table) — consistent with this function being Sun/Moon-only
(its sibling `occult_when_loc`, out of scope, generalizes to arbitrary bodies/stars).

ET→UT for `tret[0]` is a **2-pass fixed point** (`tjd - deltaT(tjd)`, then `tjd - deltaT(that)`) —
fewer passes than §5.3's 3-pass version for the global search's `tjd`.

**Rejection/retry points**, both advancing `K` in the search direction and `goto next_try`:
1. `dctr > rsplusrm` — centers too far apart even at closest approach, no eclipse visible from
   this location for this lunation.
2. `tret[0]` not strictly beyond `tjd_start` in the search direction (`±0.0001` day guard,
   textually identical to §5.4's global-search guard).

**Phase classification** (swecl.c:2235–2240) mirrors §4.4's thresholds (`ANNULAR` / `TOTAL` /
`PARTIAL` via `rsminusrm`/`fabs(rsminusrm)`/`rsplusrm`) with **one literal difference**: the
partial-eclipse branch here tests `dctr <= rsplusrm` (`<=`), whereas `eclipse_how`'s §4.4 tests
`dctr < rsplusrm` (strict `<`). Since rejection point 1 above already guarantees `dctr <=
rsplusrm` by this point, the `<=` vs `<` only matters in the boundary case `dctr == rsplusrm`
exactly (measure-zero in practice, but preserve the exact operator when porting — do not silently
"normalize" it to match §4.4).

`dctrmin = dctr` is saved — reused below (in Contacts 2/3 and Contacts 1/4) to avoid resampling
exactly at `tjd` (the already-found minimum) when locating each contact pair.

#### Contacts 2/3 (2nd/3rd contact — umbra/antumbra ingress/egress, i.e. totality or annularity
begin/end **at this location**) (swecl.c:2242–2300)
Skipped (`tret[2]=tret[3]=0`) if `dctr > fabs(rsminusrm)` (only a partial eclipse is visible
here — umbra never reaches this location). Otherwise:
```c
dc[1] = fabs(rsminusrm) - dctrmin;      // REUSED as-is from §6.2.2 — NOT resampled, NOT 0.99916-corrected
for (i = 0, t = tjd - twomin; i <= 2; i += 2, t = tjd + twomin) {
  swe_calc(t, SE_SUN, iflagcart, xs, serr);
  swe_calc(t, SE_MOON, iflagcart, xm, serr);        // cartesian only — no separate polar ls/lm fetch here
  dm = sqrt(square_sum(xm));  ds = sqrt(square_sum(xs));
  rmoon = asin(RMOON / dm) * RADTODEG;
  rmoon *= 0.99916;                     /* gives better accuracy for 2nd/3rd contacts */
  rsun = asin(RSUN / ds) * RADTODEG;
  rsminusrm = rsun - rmoon;
  x1[k] = xs[k]/ds /*ls[2]*/;  x2[k] = xm[k]/dm /*lm[2]*/;
  dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
  dc[i] = fabs(rsminusrm) - dctr;       // i = 0, 2 only — dc[1] set once, above, before this loop
}
find_zero(dc[0], dc[1], dc[2], twomin, &dt1, &dt2);
tret[2] = tjd + dt1 + twomin;
tret[3] = tjd + dt2 + twomin;
for (m = 0, dt = tensec; m < 2; m++, dt /= 10) {          // 2 passes: dt = tensec, tensec/10
  for (j = 2; j <= 3; j++) {
    swe_calc(tret[j], SE_SUN,  iflagcart | SEFLG_SPEED, xs, serr);   // ONE calc call per body per pass
    swe_calc(tret[j], SE_MOON, iflagcart | SEFLG_SPEED, xm, serr);
    for (i = 0; i < 2; i++) {
      if (i == 1) {                                        // i=1: extrapolate BACKWARD by dt using velocity,
        for (k = 0; k < 3; k++) {                           //      NOT a fresh swe_calc — reuses i=0's arrays in place
          xs[k] -= xs[k+3] * dt;
          xm[k] -= xm[k+3] * dt;
        }
      }
      dm = sqrt(square_sum(xm));  ds = sqrt(square_sum(xs));
      rmoon = asin(RMOON / dm) * RADTODEG;
      rmoon *= 0.99916;           /* gives better accuracy for 2nd/3rd contacts */
      rsun = asin(RSUN / ds) * RADTODEG;
      rsminusrm = rsun - rmoon;
      x1[k] = xs[k]/ds /*ls[2]*/;  x2[k] = xm[k]/dm /*lm[2]*/;
      dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
      dc[i] = fabs(rsminusrm) - dctr;
    }
    dt1 = -dc[0] / ((dc[0] - dc[1]) / dt);    // secant step: dc[0] at t=tret[j] (real), dc[1] at t≈tret[j]-dt (extrapolated)
    tret[j] += dt1;
  }
}
tret[2] -= swe_deltat_ex(tret[2], ifl, serr);   // single-pass ET->UT (NOT fixed-point iterated, unlike tret[0])
tret[3] -= swe_deltat_ex(tret[3], ifl, serr);
```
`iflag`/`iflagcart` used here are the same function-top values (§6.2, no `SEFLG_SPEED` in the
initial 3-point sample loop; `SEFLG_SPEED` is added **only** at the two `swe_calc` sites inside
the secant-refinement loop, freshly ORed in at each call site — not hoisted into a shared
variable).

**FP-fidelity / literal-quirk hazard — asymmetric `0.99916` correction**: every place in this
subsection where `rmoon` is freshly computed (`rmoon = asin(RMOON/dm)*RADTODEG`) is immediately
followed by `rmoon *= 0.99916` — in **both** the initial 3-point sample loop (`i=0,2`) **and**
both samples (`i=0,1`) of every secant-refinement pass. The **one** exception is `dc[1]`, which is
not resampled here at all — it is copied from `fabs(rsminusrm) - dctrmin` where `rsminusrm` was
computed back in §6.2.2 using the **uncorrected** `rmoon`. So the parabola/secant fits in this
subsection combine one uncorrected center value with corrected flanking values. This looks like it
could be an oversight in the original C, but the golden test data reflects exactly this — port
literally, do not "fix" it to apply `0.99916` uniformly.

**Secant mechanic**: the `i==1` sample is **not** a second `swe_calc` call. `xs`/`xm` (with their
appended velocity components `xs[3..6]`/`xm[3..6]` from `SEFLG_SPEED`) are fetched once at
`t = tret[j]`, used as-is for `i=0`, then linearly extrapolated backward by `dt` in place
(`xs[k] -= xs[k+3]*dt`) to approximate the position at `t = tret[j] - dt` for `i=1` — cheaper than
a second ephemeris evaluation, at the cost of first-order-only accuracy in the extrapolation (fine
given `dt` shrinks to `tensec/10 ≈ 1 s` by the final pass).

#### Contacts 1/4 (1st/4th contact — penumbra ingress/egress, i.e. visible eclipse begin/end at
this location) (swecl.c:2301–2353)
```c
dc[1] = rsplusrm - dctrmin;             // REUSED from §6.2.2, same pattern as contacts 2/3's dc[1]
for (i = 0, t = tjd - twohr; i <= 2; i += 2, t = tjd + twohr) {
  swe_calc(t, SE_SUN, iflagcart, xs, serr);
  swe_calc(t, SE_MOON, iflagcart, xm, serr);
  dm = sqrt(square_sum(xm));  ds = sqrt(square_sum(xs));
  rmoon = asin(RMOON / dm) * RADTODEG;              // NO 0.99916 correction anywhere in this subsection
  rsun = asin(RSUN / ds) * RADTODEG;
  rsplusrm = rsun + rmoon;
  x1[k] = xs[k]/ds /*ls[2]*/;  x2[k] = xm[k]/dm /*lm[2]*/;
  dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
  dc[i] = rsplusrm - dctr;
}
find_zero(dc[0], dc[1], dc[2], twohr, &dt1, &dt2);
tret[1] = tjd + dt1 + twohr;
tret[4] = tjd + dt2 + twohr;
for (m = 0, dt = tenmin; m < 3; m++, dt /= 10) {         // 3 passes: dt = tenmin, tenmin/10, tenmin/100
  for (j = 1; j <= 4; j += 3) {                          // j = 1, then j = 4 (step 3 skips 2,3)
    swe_calc(tret[j], SE_SUN,  iflagcart | SEFLG_SPEED, xs, serr);
    swe_calc(tret[j], SE_MOON, iflagcart | SEFLG_SPEED, xm, serr);
    for (i = 0; i < 2; i++) {
      if (i == 1) {
        for (k = 0; k < 3; k++) { xs[k] -= xs[k+3]*dt; xm[k] -= xm[k+3]*dt; }
      }
      dm = sqrt(square_sum(xm));  ds = sqrt(square_sum(xs));
      rmoon = asin(RMOON / dm) * RADTODEG;
      rsun = asin(RSUN / ds) * RADTODEG;
      rsplusrm = rsun + rmoon;
      x1[k] = xs[k]/ds /*ls[2]*/;  x2[k] = xm[k]/dm /*lm[2]*/;
      dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
      dc[i] = fabs(rsplusrm) - dctr;      // note: fabs() here even though rsplusrm is never negative
    }
    dt1 = -dc[0] / ((dc[0] - dc[1]) / dt);
    tret[j] += dt1;
  }
}
tret[1] -= swe_deltat_ex(tret[1], ifl, serr);   // single-pass ET->UT, same convention as contacts 2/3
tret[4] -= swe_deltat_ex(tret[4], ifl, serr);
```
Structurally identical mechanics to Contacts 2/3 (reused uncorrected `dc[1]`, one-`swe_calc`-plus-
velocity-extrapolation secant refinement, single-pass ET→UT) but: window `twohr = 2/24` day (not
`twomin`); **no** `0.99916` correction anywhere; `dc[1]/dc[i]` use `rsplusrm` (sum of radii — outer
penumbra boundary) instead of `fabs(rsminusrm)` (umbra boundary); 3 refinement passes instead of 2;
inner loop touches `j=1` and `j=4` (step `+=3`) instead of `j=2,3` (step `+=1`).

#### Visibility scan (swecl.c:2354–2384)
```c
for (i = 4; i >= 0; i--) {          /* attr for i = 0 must be kept !!! */
  if (tret[i] == 0)
    continue;
  if (eclipse_how(tret[i], SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr) == ERR)
    return ERR;
  /*if (retflag2 & SE_ECL_VISIBLE) { could be wrong for 1st/4th contact */
  if (attr[6] > 0) {                /* this is safe, sun above horizon, using app. alt. */
    retflag |= SE_ECL_VISIBLE;
    switch(i) {
    case 0: retflag |= SE_ECL_MAX_VISIBLE; break;
    case 1: retflag |= SE_ECL_1ST_VISIBLE; break;
    case 2: retflag |= SE_ECL_2ND_VISIBLE; break;
    case 3: retflag |= SE_ECL_3RD_VISIBLE; break;
    case 4: retflag |= SE_ECL_4TH_VISIBLE; break;
    default: break;
    }
  }
}
#if 1
if (!(retflag & SE_ECL_VISIBLE)) {
  if (backward) K--; else K++;
  goto next_try;
}
#endif
```
`eclipse_how` is called with the **raw** `ifl` (the ephemeris-selector parameter passed into
`eclipse_when_loc`), **not** `iflag`/`iflagcart` (which already carry `SEFLG_EQUATORIAL|
SEFLG_TOPOCTR`) — correct, since `eclipse_how` (§4.1) ORs in its own `SEFLG_EQUATORIAL|
SEFLG_TOPOCTR` internally; passing the pre-augmented `iflag` here would be redundant (harmless,
since re-OR-ing the same bits is a no-op) but the C code uses plain `ifl`, and every `eclipse_how`
call site in this function (here and in the sunrise/sunset block below) does the same — contrast
with the `swe_rise_trans` calls just below, which do **not** follow this convention (see next
subsection). `attr[6]` is the apparent (refracted) altitude, §4.9 index 6.

**The descending loop order is deliberate** (C comment: `/* attr for i = 0 must be kept !!! */`)
— `attr[]` is a single shared output array overwritten by each `eclipse_how` call; processing
`i=4,3,2,1,0` ensures the **last** write (at `i=0`, the moment of maximum eclipse) is what
remains in `attr[]` when the function returns. Port this ordering exactly. The commented-out
`/*if (retflag2 & SE_ECL_VISIBLE) {...*/` line is dead code preserved from an earlier
implementation attempt — ignore it, but its presence confirms the current `attr[6] > 0` test was a
deliberate simplification over checking `eclipse_how`'s own return flag.

#### Sunrise/sunset interaction (swecl.c:2385–2420)
```c
if ((retc = swe_rise_trans(tret[1] - 0.001, SE_SUN, NULL, iflag,
              SE_CALC_RISE | SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjdr, serr)) == ERR)
  return ERR;
if (retc == -2)                       /* circumpolar sun */
  return retflag;                     // short-circuits BEFORE the SET call; tret[5]/[6] left untouched
if ((retc = swe_rise_trans(tret[1] - 0.001, SE_SUN, NULL, iflag,
              SE_CALC_SET | SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjds, serr)) == ERR)
  return ERR;
if (retc == -2)                       /* circumpolar sun */
  return retflag;                     // independent short-circuit; tret[6] (and [5] if not yet set) untouched
if (tjds < tret[1] || (tjds > tjdr && tjdr > tret[4])) {
  if (backward) K--; else K++;
  goto next_try;                      // whole [1st,4th]-contact window is nighttime -> retry
}
if (tjdr > tret[1] && tjdr < tret[4]) {
  tret[5] = tjdr;
  if (!(retflag & SE_ECL_MAX_VISIBLE)) {
    tret[0] = tjdr;
    if ((retc = eclipse_how(tret[5], SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr)) == ERR)
      return ERR;
    retflag &= ~(SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL);
    retflag |= (retc & (SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL));
  }
}
if (tjds > tret[1] && tjds < tret[4]) {
  tret[6] = tjds;
  if (!(retflag & SE_ECL_MAX_VISIBLE)) {
    tret[0] = tjds;
    if ((retc = eclipse_how(tret[6], SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr)) == ERR)
      return ERR;
    retflag &= ~(SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL);
    retflag |= (retc & (SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL));
  }
}
return retflag;
```
**Literal-quirk hazard**: both `swe_rise_trans` calls pass `iflag` — the function-top
`SEFLG_EQUATORIAL|SEFLG_TOPOCTR|ifl` value (§6.2) — as the ephemeris-flag argument, **not** the
raw `ifl` that every `eclipse_how` call in this function (including the two re-evaluation calls
a few lines below, in this same block) uses. `swe_rise_trans` internally derives its own
topocentric/equatorial handling from the explicit `geopos` argument, so the extra
`SEFLG_EQUATORIAL|SEFLG_TOPOCTR` bits riding along in `iflag` are presumed harmless — but this is
an inconsistency in the original C (the only `swe_rise_trans` call site in this function uses the
augmented flag while every other call site uses the plain one); port the exact bit pattern passed,
do not "normalize" it to `ifl`.

Both calls anchor the search at the same instant, `tret[1] - 0.001` (1st contact, already UT by
this point, backed off by ~86 s) — used for **both** the RISE and the SET search, not
`tret[1]`/`tret[4]` respectively. `atpress`/`atemp`-equivalent trailing args are literal `0, 0`
(auto-estimate / defaults, positionally the 7th/8th parameters of `swe_rise_trans` before the
output pointer).

`retc == -2` (circumpolar — sun does not rise/set within the search range) short-circuits the
**entire function**, returning `retflag` as accumulated so far, **without** setting `tret[5]`/
`tret[6]` (they retain whatever the caller pre-initialized, typically `0`) and without touching
`attr[]` further. This is a distinct non-error control path — the Rust port's equivalent should
be a normal early return (e.g. an `Option`/enum branch, not a `Result::Err`), matching the
"circumpolar body" case documented for `swe_rise_trans` itself, not a failure.

`SE_BIT_DISC_BOTTOM` = rise/set convention using the bottom limb touching the horizon. **Both
the sunrise and sunset blocks can fire** for a single very-long high-latitude event; if so, the
sunset block (evaluated second) wins the final `tret[0]`/`attr` overwrite — a sequential-overwrite
behavior to preserve exactly.

### 6.3 `tret[]` full index summary (local search — DIFFERENT semantics than §5.8!)
| Index | Meaning |
|---|---|
| 0 | Time (UT) of maximum eclipse **as seen from this location** (re-anchored to sunrise/sunset if the true max wasn't visible) |
| 1 | Time of first contact (1st contact — penumbra ingress) |
| 2 | Time of second contact (2nd contact — umbra/antumbra ingress; 0 if only partial visible here) |
| 3 | Time of third contact (3rd contact — umbra/antumbra egress; 0 if only partial visible here) |
| 4 | Time of fourth contact (4th contact — penumbra egress) |
| 5 | Time of sunrise between 1st and 4th contact (0 if none / circumpolar) |
| 6 | Time of sunset between 1st and 4th contact (0 if none / circumpolar) |
| 7–9 | unused, always 0 |

**Contrast with §5.8**: the global search's `tret[2]/[3]` = eclipse begin/end and `[4]/[5]` =
totality begin/end and `[6]/[7]` = center-line begin/end, whereas the local search's
`tret[1]/[4]` = 1st/4th (penumbra) contact and `[2]/[3]` = 2nd/3rd (umbra) contact. **Do not
conflate the two index conventions when designing the Rust struct** — they must be two distinct
named-field structs (e.g. `GlobalEclipseTimes` vs `LocalEclipseTimes`), not one shared `[f64; 10]`
newtype, or the semantic mismatch will be silently ported forward.

### 6.4 Return value
`ERR` on ephemeris failure; else `SE_ECL_TOTAL`/`ANNULAR`/`PARTIAL` OR'd with `SE_ECL_VISIBLE` and
whichever of `MAX/1ST/2ND/3RD/4TH_VISIBLE` applied; `SE_ECL_NONCENTRAL` merged in later by the
`swe_sol_eclipse_when_loc` wrapper (§6.1 step 4). The search loop retries internally until a
visible eclipse is found (bounded by ephemeris range and, for near-polar `geopos`, potentially
very slow — the public doc comment warns of this).

---

## 7. `find_maximum` / `find_zero` — shared parabolic-interpolation helpers (swecl.c:4133–4162)

Both fit a parabola through three equally-spaced samples `(y00, y11, y2)` at `x = -dx, 0, +dx`
(i.e. `y11` is the **center** sample) and solve analytically:
```c
c = y11;
b = (y2 - y00) / 2.0;
a = (y2 + y00) / 2.0 - c;
```
(`a,b,c` are the standard parabola coefficients `y = a*x² + b*x + c` in the local `x∈[-1,1]`
scaled coordinate, i.e. `dx` is factored out.)

```c
static int find_maximum(double y00, double y11, double y2, double dx, double *dxret, double *yret) {
  x = -b / (2*a);                          // vertex location, in units of dx
  y = (4*a*c - b*b) / (4*a);                // vertex value
  *dxret = (x - 1) * dx;                    // offset from the LAST sample (t = t2 = center+dx), not from the center!
  if (yret) *yret = y;
  return OK;
}
```
**Note the `(x - 1) * dx` convention**: callers always sample at `t-dx, t, t+dx` and then advance
`t += dtint + dx` — i.e. `find_maximum` returns the vertex offset measured from the **rightmost**
sample point `t+dx`, not from the center `t`. This is a systematic convention throughout the file
(used identically in `find_zero`); port it exactly (an off-by-`dx` error here silently shifts
every refinement step in the search functions).

```c
static int find_zero(double y00, double y11, double y2, double dx, double *dxret, double *dxret2) {
  if (b*b - 4*a*c < 0) return ERR;                       // no real root (parabola doesn't cross zero)
  x1 = (-b + sqrt(b*b - 4*a*c)) / (2*a);
  x2 = (-b - sqrt(b*b - 4*a*c)) / (2*a);
  *dxret  = (x1 - 1) * dx;
  *dxret2 = (x2 - 1) * dx;
  return OK;
}
```
Standard quadratic formula on the same `a,b,c`; returns **both** roots (used as the two contact
times — e.g. begin and end — from a single 3-point sample of a function that dips through zero
twice). Callers must be prepared for `ERR` (ex: passed to `swe_sol_eclipse_when_glob`/
`eclipse_when_loc` without an explicit check in some call sites — check the actual call sites
during port, since a few do not test the return value and would need an explicit `Result`/`Option`
handling equivalent in Rust that the C silently skips).

---

## 8. Stateless-Port Summary

Global-state touch points in this module that the Rust `Ephemeris` (stateless, `&self`-only)
must instead thread as explicit parameters:

| C global | Set by | Read by | Rust equivalent |
|---|---|---|---|
| `swed.topd` (topocentric observer position) | `swe_set_topo` (eclipse_how §4.1, eclipse_when_loc §6.2, swe_sol_eclipse_how §4.11 step 7) | Any subsequent `SEFLG_TOPOCTR` calc | Pass observer lon/lat/height explicitly into the position-calculation call |
| `swed.tid_acc` (tidal acceleration mode) | `swi_set_tid_acc` (every public entry point: §3.8, §4.11, §5, §6.1) | `swe_deltat_ex` internals | Pass tidal-acceleration mode explicitly via `EphemerisConfig` or a parameter |
| `swed.ast_diam` (asteroid physical diameter) | SE1 asteroid-file loader | §0 body-radius lookup in `eclipse_where`/`eclipse_how` | Thread the resolved asteroid diameter explicitly (already-loaded orbital-element data, not a global) |
| `swed.oec` (mean/true obliquity struct, `&swed.oec` aliased as `oe` in `eclipse_where`, only used for the `SEFLG_NONUT` sidereal-time branch) | nutation/obliquity calc earlier in the same logical pipeline | `eclipse_where` §3.1 sidt calc | Pass the computed obliquity value explicitly (already available from whatever nutation/obliquity call the Rust pipeline made) |

No cached planetary positions are read back from global state in this file (contrast with
`swi_deflect_light`'s use of `psdp->x`, documented in the project's `stateless_tolerance` policy)
— every Sun/Moon position used here comes from a `swe_calc`/`calc_planet_star` call made within
the same function, so this module should port to bit-for-bit (modulo the FP-fidelity hazards
flagged inline above) equivalence more readily than the deflection code did.
