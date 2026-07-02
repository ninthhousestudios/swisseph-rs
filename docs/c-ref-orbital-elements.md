# C Reference: Orbital Elements — swecl.c

Porting reference for osculating (Keplerian) orbital element computation and
max/min/true distance search. Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `get_gmsm` | swecl.c:5687–5742 | Yes — static helper, GM(Sun/Earth-relative) for the body |
| `swe_get_orbital_elements` | swecl.c:5783–5971 | Yes — full state-vector → Kepler-element derivation |
| `osc_get_orbit_constants` | swecl.c:5973–5999 | Yes — static helper, precomputes P/Q rotation vectors |
| `osc_get_ecl_pos` | swecl.c:6001–6015 | Yes — static helper, ecc. anomaly → ecliptic cartesian |
| `get_dist_from_2_vectors` | swecl.c:6017–6024 | Yes — trivial 3-vector distance |
| `osc_iterate_max_dist` | swecl.c:6026–6060 | Yes — static helper, coordinate-descent max search |
| `osc_iterate_min_dist` | swecl.c:6062–6096 | Yes — static helper, coordinate-descent min search |
| `orbit_max_min_true_distance_helio` | swecl.c:6101–6128 | Yes — heliocentric-only branch (Sun/Moon/HELCTR/BARYCTR) |
| `swe_orbit_max_min_true_distance` | swecl.c:6170–6287 | Yes — full geocentric two-ellipse search + heliocentric dispatch |

Not ported (out of scope for this doc): `swe_nod_aps` / `swe_nod_aps_ut` (swecl.c:5075–5665,
5656–5665) — separate nodes/apsides feature, shares the `plmass`/`ipl_to_elem`/`el_*` tables
declared just above this block (swecl.c:5010–5074) but is a distinct C function with its own
mean/osculating dispatch. Only the table *declarations* are read here because `get_gmsm` reuses
`plmass`/`ipl_to_elem`; no sibling ref doc exists yet for `swe_nod_aps` (none found under
`docs/c-ref-*.md` at time of writing) — if/when that module is ported, extract `el_node`,
`el_peri`, `el_incl`, `el_ecce`, `el_sema` (swecl.c:5012–5061) into that doc rather than
duplicating them here.

## Constants

| Name | Value | Location | Notes |
|---|---|---|---|
| `HELGRAVCONST` | 1.32712440017987e+20 | sweph.h:278 | G·M(Sun), m³/s² (AA 2006 K6) |
| `GEOGCONST` | 3.98600448e+14 | sweph.h:279 | G·M(Earth), m³/s² (AA 1996 K6) |
| `AUNIT` | 1.49597870700e+11 | sweph.h:273 | AU in meters (DE431) |
| `KGAUSS` | 0.01720209895 | sweph.h:280 | Gaussian gravitational constant (unused in this file — listed for completeness) |
| `EARTH_MOON_MRAT` | `1 / 0.0123000383` ≈ 81.30056...  | sweph.h:265 | Earth/Moon mass ratio (AA 2006 K7); alternate DE431/DE406 values exist under `#if 0` (sweph.h:267,270) — dead code, ignore |
| `J2000` | 2451545.0 | sweph.h:67 | epoch for `T` |
| `DEGTORAD` / `RADTODEG` | π/180 / 180/π | sweodef.h:262-266 | standard |
| `plmass[9]` | see table below | swecl.c:5063–5073 | Sun-to-planet mass **ratios** (bigger number ⇒ smaller planet mass) |
| `ipl_to_elem[15]` | `{2,0,0,1,3,4,5,6,7,0,0,0,0,0,2}` | swecl.c:5074 | body index → `plmass` row (see quirk below) |
| `SEFLG_ORBEL_AA` | `= SEFLG_TOPOCTR` (32768) | swephexp.h:207 (def), 206 (TOPOCTR) | **bit-aliased** onto TOPOCTR — orbital elements have no topocentric variant, so the bit is repurposed |

### `plmass` table (swecl.c:5063–5073) — Sun/planet mass ratios

| index | value | body |
|---|---|---|
| 0 | 6023600 | Mercury |
| 1 | 408523.719 | Venus |
| 2 | 328900.5 | Earth **and Moon** combined |
| 3 | 3098703.59 | Mars |
| 4 | 1047.348644 | Jupiter |
| 5 | 3497.9018 | Saturn |
| 6 | 22902.98 | Uranus |
| 7 | 19412.26 | Neptune |
| 8 | 136566000 | Pluto |

### `ipl_to_elem` table (swecl.c:5074) — body ID → `plmass` row

Indexed by `SE_*` body id (SE_SUN=0 .. SE_EARTH=14):

| `SE_*` id | 0 (SUN) | 1 (MOON) | 2 (MERCURY) | 3 (VENUS) | 4 (MARS) | 5 (JUPITER) | 6 (SATURN) | 7 (URANUS) | 8 (NEPTUNE) | 9 (PLUTO) | 10-13 (nodes/apsides) | 14 (EARTH) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| value | 2 | 0 | 0 | 1 | 3 | 4 | 5 | 6 | 7 | **0** | 0 | 2 |

**Quirk — `ipl_to_elem[SE_PLUTO] = 0`, not 8.** `ipl_to_elem` was built primarily to index the
8-row `el_node`/`el_peri`/`el_incl`/`el_ecce`/`el_sema` tables used by `swe_nod_aps`
(Mercury=0..Neptune=7, no Pluto row exists there since that function's mean-elements branch
excludes Pluto — see swecl.c:5165-5166 range check `SE_SUN..SE_NEPTUNE`). `get_gmsm` reuses the
*same* table to index the 9-row `plmass` array (which *does* have a Pluto row at index 8), and
inherits the stale mapping. Net effect in `get_gmsm` (swecl.c:5687-5742):
- **Two-body branch** (no `SEFLG_ORBEL_AA`, swecl.c:5714-5717): for Pluto,
  `plm = 1.0 / plmass[ipl_to_elem[SE_PLUTO]] = 1.0 / plmass[0] = 1/6023600 ≈ 1.66e-7`, instead of
  the "correct" `1.0 / plmass[8] = 1/136566000 ≈ 7.32e-9`. Relative error in `Gmsm` for Pluto is
  ≈1.6e-7 (astronomically negligible, but **must be replicated bit-for-bit** for golden-test
  fidelity — do not "fix" it).
- **`SEFLG_ORBEL_AA` branch** (swecl.c:5706-5711): the summation loop
  `for (j = ipl; j >= SE_MERCURY; j--) plm += 1.0/plmass[ipl_to_elem[j]]` double-counts Mercury's
  contribution for Pluto (once at `j = SE_PLUTO` via the stale `0` mapping, once again at
  `j = SE_MERCURY` itself), instead of omitting Pluto's own mass as intended. Same ~1.6e-7-scale
  effect — replicate exactly, do not "fix."

Port: hard-code these two tables as Rust `const` arrays (`PLMASS: [f64; 9]`,
`IPL_TO_ELEM: [usize; 15]` or equivalent `Body → usize` match) exactly as above, including the
Pluto quirk. Do not attempt to derive index 8 for Pluto in the two-body/AA branches.

## 1. `get_gmsm` — GM for the orbit's central body (swecl.c:5687–5742)

Signature: `get_gmsm(tjd_et, ipl, iflag, r, &gmsm, serr) -> i32` where `r` is the body's
heliocentric distance (AU) already computed by the caller.

Sets up (swecl.c:5691-5693):
```c
iflJ2000p = (iflag & (SEFLG_EPHMASK|SEFLG_HELCTR|SEFLG_BARYCTR)) | SEFLG_J2000|SEFLG_TRUEPOS|SEFLG_NONUT;
if (!(iflJ2000p & (SEFLG_HELCTR|SEFLG_BARYCTR)))
    iflJ2000p |= SEFLG_HELCTR;
```
(used only in the asteroid/ORBEL_AA sub-branch to re-query planet positions.)

### Case A — Moon (swecl.c:5694-5695)
```c
Gmsm = GEOGCONST * (1 + 1/EARTH_MOON_MRAT) / AUNIT^3 * 86400^2
```
Geocentric GM scaled by (1 + Moon/Earth mass ratio contribution... actually `1/EARTH_MOON_MRAT`
is Moon-to-Earth fraction added to 1), converted from m³/s² to AU³/day².

### Case B — Mercury..Pluto or Earth (swecl.c:5697-5721)

`plm` (planet-mass-ratio term added to the central mass) depends on `SEFLG_ORBEL_AA`:

- **If `SEFLG_ORBEL_AA` set** (AA "sum masses inside orbit" method):
  - `ipl == SE_EARTH` (swecl.c:5702-5705): `plm = 1/plmass[ipl_to_elem[EARTH]] + 1/plmass[ipl_to_elem[VENUS]] + 1/plmass[ipl_to_elem[MERCURY]]` (explicit 3-term list — Earth+Venus+Mercury, hard-coded, not a loop).
  - otherwise (swecl.c:5706-5711): loop `for (j = ipl; j >= SE_MERCURY; j--) plm += 1/plmass[ipl_to_elem[j]]`, and if `ipl >= SE_MARS` additionally add `1/plmass[ipl_to_elem[SE_EARTH]]` (Earth's mass folded in for objects beyond Mars, since the loop itself stops before reaching Earth's slot — Earth is id 14, outside the `j >= SE_MERCURY` descending loop from `ipl`).
- **Else** (two-body, default): `plm = 1/plmass[ipl_to_elem[ipl]]` (swecl.c:5715) — single term, subject to the Pluto quirk above.

Then (swecl.c:5717): `Gmsm = HELGRAVCONST * (1 + plm) / AUNIT^3 * 86400^2`.

`#ifdef TEST_ORBEL_AA` block (swecl.c:5718-5721, 5667-5686): test-only correction factor table
`Gmsm_factor_AA[9]` — **not compiled by default** (guarded by a build macro not defined in normal
builds). Skip entirely; do not port.

### Case C — asteroid / fictitious body (swecl.c:5722-5738)

- If `SEFLG_ORBEL_AA`: loop `j` from `SE_MERCURY` to `SE_PLUTO` inclusive, calling
  `swe_calc(tjd_et, j, iflJ2000p, x, serr)` for each planet (real `swe_calc`/`Ephemeris::calc`
  invocation, J2000/heliocentric-or-baryctr per `iflJ2000p`), and if `r > x[2]` (the object's
  distance exceeds planet `j`'s heliocentric distance — i.e. planet is interior) add
  `1/plmass[ipl_to_elem[j]]`. Same test then repeated for `SE_EARTH` (its own `swe_calc` call,
  flags with `BARYCTR|HELCTR` stripped per line 5826-pattern — actually uses `iflJ2000p` as-is at
  swecl.c:5732).
- Else: `plm = 0` (pure two-body, no mass correction).
- `Gmsm = HELGRAVCONST * (1 + plm) / AUNIT^3 * 86400^2` (swecl.c:5737), same conversion as Case B.

**Error propagation**: any `swe_calc` failure inside the AA loops returns `ERR` immediately
(swecl.c:5727-5728, 5732-5733).

## 2. `swe_get_orbital_elements` (swecl.c:5783–5971)

Signature: `swe_get_orbital_elements(tjd_et, ipl, iflag, dret[50], serr) -> i32`.

### 2.1 Rejected bodies (swecl.c:5803-5807)

```c
if (ipl <= 0 || ipl == SE_MEAN_NODE || ipl == SE_TRUE_NODE
    || ipl == SE_MEAN_APOG || ipl == SE_OSCU_APOG
    || ipl == SE_INTP_APOG || ipl == SE_INTP_PERG)
  return ERR;  // "object %d not valid"
```
`ipl <= 0` rejects `SE_SUN` (id 0) and any negative id. Nodes/apsides (mean, true, oscu, and the
interpolated apogee/perigee variants — ids 10-13, 21-22) are rejected because osculating elements
of a node/apsis point are not meaningful.

### 2.2 Heliocentric distance probe + center-flag decision (swecl.c:5808-5819)

```c
iflJ2000p = (iflag & SEFLG_EPHMASK)|SEFLG_J2000|SEFLG_TRUEPOS|SEFLG_NONUT|SEFLG_SPEED;
x = swe_calc(tjd_et, ipl, iflJ2000p, ..., serr);  // may ERR
r = x[2];  // heliocentric distance in AU (still default center — see note)
if (ipl != SE_MOON) {
  if ((iflag & SEFLG_BARYCTR) && r > 6)
    iflJ2000 |= SEFLG_BARYCTR;   // barycentric only allowed beyond ~Jupiter's distance
  else
    iflJ2000 |= SEFLG_HELCTR;
}
```
Note: `iflJ2000p` here has **neither** `HELCTR` nor `BARYCTR` set explicitly, so `x[2]` is
whatever `swe_calc`'s default center is (heliocentric by default for planets in this ephemeris
mode). The `r > 6` AU threshold gates `SEFLG_BARYCTR` to "beyond Jupiter" (barycentric only makes
physical sense for the outer planets, per the function's doc comment at swecl.c:5754-5756).

### 2.3 GM and final position query (swecl.c:5820-5830)

```c
get_gmsm(tjd_et, ipl, iflag, r, &Gmsm, serr);   // may ERR
xpos = swe_calc(tjd_et, ipl, iflJ2000, ..., serr);  // iflJ2000 = J2000|XYZ|TRUEPOS|NONUT|SPEED (+ HELCTR/BARYCTR from above)
if (ipl == SE_EARTH) {
  xposm = swe_calc(tjd_et, SE_MOON, iflJ2000 & ~(BARYCTR|HELCTR), ..., serr);  // geocentric Moon
  for j in 0..=5: xpos[j] += xposm[j] / (EARTH_MOON_MRAT + 1.0);   // Earth -> EMB barycenter
}
```
`iflJ2000` requests **cartesian** (`SEFLG_XYZ`), **true/geometric position** (no light-time/
aberration/deflection — `SEFLG_TRUEPOS`), **no nutation** (`SEFLG_NONUT`), **J2000 mean equinox**
(`SEFLG_J2000`), and **speed** (`SEFLG_SPEED`, needed for the velocity components in `xpos[3..6]`
used below). For `SE_EARTH`, the barycentric-corrected position becomes the Earth-Moon
Barycenter (EMB) — orbital elements for "Earth" are actually EMB elements, matching AA convention
(doc comment swecl.c:5824).

### 2.4 Node vector via angular-momentum-perpendicular projection (swecl.c:5831-5840)

```c
fac = xpos[2] / xpos[5];         // r_z / v_z  (NOT literal fraction naming — z-component ratio)
sgn = xpos[5] / fabs(xpos[5]);   // sign of v_z
for j in 0..=2:
  xn[j] = (xpos[j] - fac * xpos[j+3]) * sgn;
  xs[j] = -xn[j];
rxy = sqrt(xn[0]^2 + xn[1]^2)
cosnode = xn[0] / rxy
sinnode = xn[1] / rxy
```
`xn` is the position vector projected onto the orbital plane's node line (intersection with the
xy-plane), scaled/signed so it points toward the ascending node; `xs` is its negation (used later
for the descending-node correction, though the descending node isn't part of `dret[]` here — see
§2.10). This is a first-pass, un-refined node estimate; refined further at §2.9.

### 2.5 Inclination via cross product (swecl.c:5841-5850)

```c
swi_cross_prod(xpos, xpos+3, xnorm);      // xnorm = r × v  (specific angular momentum vector)
rxy = xnorm[0]^2 + xnorm[1]^2
c2  = rxy + xnorm[2]^2                     // |r×v|^2
rxyz = sqrt(c2); rxy = sqrt(rxy)
sinincl = rxy / rxyz
cosincl = sqrt(1 - sinincl^2)
if xnorm[2] < 0: cosincl = -cosincl        // retrograde orbit (e.g. 20461 Dioretsa)
incl = acos(cosincl) * RADTODEG
```
`swi_cross_prod(a,b,x)` (swephlib.c:160-165): `x[0]=a1*b2-a2*b1; x[1]=a2*b0-a0*b2; x[2]=a0*b1-a1*b0`.
`c2` (= `|r×v|²`, the squared specific angular momentum magnitude) is retained for reuse in §2.6.

### 2.6 Argument of latitude, semimajor axis, eccentricity (swecl.c:5851-5864)

```c
cosu = xpos[0]*cosnode + xpos[1]*sinnode
sinu = xpos[2] / sinincl
uu = atan2(sinu, cosu)                      // argument of latitude (radians)
rxyz = sqrt(square_sum(xpos))               // |r|
v2 = square_sum(xpos+3)                     // |v|^2
sema = 1.0 / (2.0/rxyz - v2/Gmsm)           // vis-viva
pp = c2 / Gmsm                              // semi-latus rectum p = h^2/GM
ecce = pp / sema
if ecce > 1: ecce = 1                       // clamp (guards sqrt domain, hyperbolic-orbit safety)
ecce = sqrt(1 - ecce)
```
`square_sum(x)` / `dot_prod(x,y)` are macros at sweph.h:308-309:
`square_sum(x) = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]`, `dot_prod(x,y) = x[0]*y[0]+x[1]*y[1]+x[2]*y[2]`.
Note the reassignment idiom: `ecce` holds `p/a` (= `1 - e²`, NOT the eccentricity) transiently
before the final `sqrt(1 - ecce)` overwrites it with the true eccentricity — preserve this exact
two-step form for FP fidelity (do not algebraically simplify to `sqrt(1 - pp/sema)` in one line
if bitwise parity matters; the intermediate clamp `if (ecce > 1) ecce = 1` operates on `pp/sema`
specifically).

### 2.7 Eccentric and true anomaly (swecl.c:5865-5880)

```c
ecce2 = ecce; if (ecce2 == 0) ecce2 = 0.0000000001;   // avoid div-by-zero for circular orbits
cosE = 1/ecce2 * (1 - rxyz/sema)
sinE = 1/ecce2 / sqrt(sema*Gmsm) * dot_prod(xpos, xpos+3)
eanom = swe_degnorm(atan2(sinE, cosE) * RADTODEG)     // eccentric anomaly, [0,360)

ny = 2 * atan(sqrt((1+ecce)/(1-ecce)) * sinE / (1 + cosE))   // true anomaly (radians), half-angle formula
tanom = swe_degnorm(ny * RADTODEG)
if (eanom > 180 && tanom < 180) tanom += 180
if (eanom < 180 && tanom > 180) tanom -= 180
```
The two `if` correction lines are a quadrant-matching patch: the half-angle `atan` formula for
true anomaly can land in the wrong branch relative to `eanom`; nudging `tanom` by 180° when the
two disagree across the 180° boundary keeps them in the same half-revolution. Commented-out
alternative formulas (`acos(...)`) remain in source at swecl.c:5872, 5880 — dead code, ignore.

### 2.8 Mean anomaly (swecl.c:5881-5882)

```c
manom = swe_degnorm(eanom - ecce * RADTODEG * sin(eanom * DEGTORAD))
```
Kepler's equation `M = E - e·sin(E)`, with `e·sin(E)` computed by first converting `E` to
radians for the `sin`, then scaling the whole `ecce * sin(E_rad)` product by `RADTODEG` (i.e.
`ecce * RADTODEG * sin(eanom * DEGTORAD)`, not `RADTODEG * (ecce * sin(...))` — same value
either way for pure multiplication, so grouping is not FP-sensitive here beyond standard
associativity, but keep the literal factor order for direct-comparison ease against the C.)

### 2.9 Perihelion direction, aphelion, and node/apsis vector refinement (swecl.c:5883-5920)

```c
xq[0] = swi_mod2PI(uu - ny)          // argument of perihelion (radians, from arg-of-latitude minus true anomaly)
parg = xq[0] * RADTODEG
xq[1] = 0                             // latitude = 0 (perihelion lies in orbital plane by definition)
xq[2] = sema * (1 - ecce)             // perihelion distance
swi_polcart(xq, xq)                   // -> cartesian, in orbital-plane-relative frame (lon=arg.peri, lat=0, r=q)
swi_coortrf2(xq, xq, -sinincl, cosincl)   // rotate by inclination (note NEGATED sinincl)
swi_cartpol(xq, xq)                   // back to polar
xq[0] += atan2(sinnode, cosnode)      // add node longitude -> perihelion longitude in ecliptic
xa[0] = swi_mod2PI(xq[0] + PI)        // aphelion longitude = perihelion + 180 deg
xa[1] = -xq[1]                        // aphelion latitude = -perihelion latitude
xa[2] = sema * (1 + ecce)             // aphelion distance
swi_polcart(xq, xq)                   // xq -> cartesian ecliptic perihelion vector
swi_polcart(xa, xa)                   // xa -> cartesian ecliptic aphelion vector
```
`swi_coortrf2(xpo, xpn, sineps, coseps)` (swephlib.c:299-308): rotates about the x-axis:
`xpn[0]=xpo[0]; xpn[1]=xpo[1]*coseps+xpo[2]*sineps; xpn[2]=-xpo[1]*sineps+xpo[2]*coseps`. Passing
`-sinincl` here (rather than `+sinincl`) selects the ecliptic-to-orbital-plane rotation direction
consistent with this coordinate convention — preserve the sign literally, do not "fix" to
`sinincl` by symmetry-reasoning.

A dead `do_focal_point` branch (commented out, swecl.c:5896-5900, would have set
`xa[2] = sema*ecce*2`) — ignore, not compiled.

Then, node/descending-node position vectors are recomputed at higher precision (swecl.c:5903-5920):
```c
ny  = swi_mod2PI(ny - uu)                          // true anomaly AT the node (ascending)
ny2 = swi_mod2PI(ny + PI)                          // true anomaly at descending node
cosE  = cos(2*atan(tan(ny/2)  / sqrt((1+ecce)/(1-ecce))))    // eccentric anomaly at ascending node (via inverse half-angle)
cosE2 = cos(2*atan(tan(ny2/2) / sqrt((1+ecce)/(1-ecce))))    // eccentric anomaly at descending node
rn  = sema * (1 - ecce*cosE)                        // true orbital radius at ascending node
rn2 = sema * (1 - ecce*cosE2)                        // true orbital radius at descending node
ro  = sqrt(square_sum(xn))                          // old (first-pass, §2.4) ascending-node vector magnitude
ro2 = sqrt(square_sum(xs))                          // old descending-node vector magnitude
for j in 0..=2:
  xn[j] *= rn / ro                                   // rescale first-pass node vector to true orbital radius
  xs[j] *= rn2 / ro2
swi_cartpol(xn, xn); swi_cartpol(xq, xq)             // (xq here converts the already-cartesian perihelion vector — re-derives its polar lon for use below)
```
This is a refinement pass: the crude node-line vector from §2.4 (built from position/velocity
z-ratio projection) is rescaled along its own direction to sit exactly on the Kepler ellipse at
the node's true anomaly, giving an exact orbital-plane node distance. Only `xn`'s longitude
(`node`) is actually consumed downstream — `xs`/`xn[2]` (descending node, radius) values are
computed but not exposed via `dret[]` (this function has no descending-node output; that's
`swe_nod_aps`'s job).

### 2.10 Final angle assembly (swecl.c:5921-5925)

```c
node = xn[0] * RADTODEG                    // longitude of ascending node
peri = swe_degnorm(node + parg)            // longitude of perihelion = node + arg.peri
mlon = swe_degnorm(manom + peri)           // mean longitude = mean anomaly + longitude of perihelion
```

### 2.11 Period and daily motion (swecl.c:5926-5948)

```c
csid = sema * sqrt(sema)                    // Kepler's third law: P (years) = a^1.5  (sidereal, in years, for a in AU)
if (ipl == SE_MOON) {
  semam = sema * AUNIT / 383397772.5        // convert Moon's "AU" semimajor axis to units of the Moon's mean distance (383397772.5 mm? — literal constant, treat as a fixed reference distance in meters)
  csid = semam * sqrt(semam)                 // period in sidereal MONTHS (Kepler's law in lunar-distance units)
  csid *= 27.32166 / 365.25636300            // months -> years (27.32166 = sidereal month in days, 365.25636300 = sidereal year in days)
}
dmot = 0.9856076686 / csid                  // mean daily motion, deg/day (0.9856076686 = 360/365.25636 = Earth's mean daily motion)
csid *= 365.25636 / 365.242189              // sidereal period: sidereal years -> "tropical years J2000" units (see literal ratio, not a name change of period type)

T = (tjd_et - J2000) / 365250.0             // Julian millennia from J2000
T2=T*T; T3=T2*T; T4=T3*T; T5=T4*T
pa = (50288.200 + 222.4045*T + 0.2095*T2 - 0.9408*T3 - 0.0090*T4 + 0.0010*T5) / 3600.0 / 365250.0
   // general precession in longitude (Simon et al. 1994), arcsec/millennium -> deg/day

ysid = (1295977422.83429 - 2*2.0441*T - 3*0.00523*T*T) / 3600.0 / 365250.0
ysid = 360.0 / ysid                         // sidereal year rate, deg/day
ytrop = (1296027711.03429 + 2*109.15809*T + 3*0.07207*T2 - 4*0.23530*T3 - 5*0.00180*T4 + 6*0.00020*T5) / 3600.0 / 365250.0
ytrop = 360.0 / ytrop                       // tropical year rate, deg/day

ctro = 360.0 / (dmot + pa) / 365.242189     // tropical period in years (uncorrected)
ctro *= ysid / ytrop                        // corrected to "tropical years J2000" units

csyn = (ipl == SE_EARTH) ? 0
     : 360.0 / (0.9856076686 - dmot)        // synodic period in days
```
The `ysid`/`ytrop` polynomials use the derivative-of-mean-longitude-rate convention: each
coefficient in the source polynomial for mean longitude rate is multiplied by its power's
coefficient (`2*`, `3*`, etc.) — i.e. these are literally *rate* polynomials (already
differentiated), not integrated positions; port the coefficients and multipliers exactly as
written (`2 * 2.0441 * T`, `3 * 0.00523 * T * T`, etc.), do not attempt to re-derive them from a
position polynomial. `csyn` sign convention (per the public dret[13] doc comment,
swecl.c:5777-5778): negative for inner planets (Venus, Mercury) and the Moon — this falls out
naturally from `0.9856076686 - dmot` going negative when `dmot > 0.9856076686` (i.e. the body's
own mean motion exceeds Earth's, meaning it's an inner/faster-orbiting body from Earth's
perspective... note Moon is a special case already excluded from this synodic formula's normal
domain but still uses it, since only `ipl == SE_EARTH` is special-cased to 0).

### 2.12 `dret[]` output — ALL populated slots (swecl.c:5949-5965)

| Slot | Value | C expression |
|---|---|---|
| `dret[0]` | semimajor axis (AU) | `sema` |
| `dret[1]` | eccentricity | `ecce` |
| `dret[2]` | inclination (deg) | `incl` |
| `dret[3]` | longitude of ascending node (deg) | `node` |
| `dret[4]` | argument of perihelion (deg) | `parg` |
| `dret[5]` | longitude of perihelion (deg) | `peri` |
| `dret[6]` | mean anomaly (deg) | `manom` |
| `dret[7]` | true anomaly (deg) | `tanom` |
| `dret[8]` | eccentric anomaly (deg) | `eanom` |
| `dret[9]` | mean longitude (deg) | `mlon` |
| `dret[10]` | sidereal orbital period, tropical years (J2000) | `csid` |
| `dret[11]` | mean daily motion (deg/day) | `dmot` |
| `dret[12]` | tropical period, years | `ctro` |
| `dret[13]` | synodic period, days (negative for inner planets/Moon) | `csyn` |
| `dret[14]` | JD (TT) of perihelion passage | `tjd_et - dret[6]/dmot` i.e. `tjd_et - manom/dmot` |
| `dret[15]` | perihelion distance (AU) | `sema*(1-ecce)` |
| `dret[16]` | aphelion distance (AU) | `sema*(1+ecce)` |

`dret[14]` = `tjd_et - manom / dmot`: walks the mean anomaly back to zero at the body's mean
daily-motion rate, giving the epoch of the most recent perihelion passage. Slots 17+ are **not
populated by this function** (no writes past `dret[16]` in the source — the doc-comment header
at swecl.c:5763-5781 only documents through `dret[16]`, and code confirms no further slots are
touched). Note the `dret[]` doc comment predates slots 13-16 having been added (comment block
originally only described 0-12; slots 13-16 were appended later without updating the numbered
comment structure above them fully — cross-check against actual code, which is authoritative and
matches this table).

### 2.13 Error propagation

Any `swe_calc`/`get_gmsm` failure returns `ERR` immediately (swecl.c:5810-5811, 5820-5821,
5822-5823, 5826-5827) with `serr` populated by the failing call.

## 3. `osc_get_orbit_constants` (swecl.c:5973-5999)

Precomputes a 12-element "PQR" rotation-and-shape-parameter block from `dret[]` (or any
`dp[0..4]` = a, e, i, Ω, ω in AU/dimensionless/degrees), for repeated fast position evaluation at
arbitrary eccentric anomaly (used by the distance-search routines below):

```c
cosnode=cos(node*DEGTORAD); sinnode=sin(node*DEGTORAD)
cosincl=cos(incl*DEGTORAD); sinincl=sin(incl*DEGTORAD)
cosparg=cos(parg*DEGTORAD); sinparg=sin(parg*DEGTORAD)
fac = sqrt((1-ecce)*(1+ecce))     // = sqrt(1 - e^2), grouped as difference-of-squares product (not 1-e*e) — preserve literal form
pqr[0] = cosparg*cosnode - sinparg*cosincl*sinnode
pqr[1] = -sinparg*cosnode - cosparg*cosincl*sinnode
pqr[2] = sinincl*sinnode
pqr[3] = cosparg*sinnode + sinparg*cosincl*cosnode
pqr[4] = -sinparg*sinnode + cosparg*cosincl*cosnode
pqr[5] = -sinincl*cosnode
pqr[6] = sinparg*sinincl
pqr[7] = cosparg*sinincl
pqr[8] = cosincl
pqr[9] = sema
pqr[10] = ecce
pqr[11] = fac
```
`pqr[0..8]` is the classical Gauss P/Q rotation-matrix representation of the orbital-plane basis
vectors in ecliptic coordinates (P = perihelion direction, Q = 90°-advanced direction in the
orbital plane); `pqr[9..11]` cache `sema`, `ecce`, and `fac` for reuse in `osc_get_ecl_pos`.

## 4. `osc_get_ecl_pos` (swecl.c:6001-6015)

Given eccentric anomaly `ean` (degrees) and a `pqr[12]` block from §3, computes the ecliptic
cartesian position on the ellipse:
```c
cose=cos(ean*DEGTORAD); sine=sin(ean*DEGTORAD)
sema=pqr[9]; ecce=pqr[10]; fac=pqr[11]
x0 = sema*(cose - ecce)          // in-plane P-axis coordinate
x1 = sema*fac*sine               // in-plane Q-axis coordinate
xp[0] = pqr[0]*x0 + pqr[1]*x1
xp[1] = pqr[3]*x0 + pqr[4]*x1
xp[2] = pqr[6]*x0 + pqr[7]*x1
```
Pure two-body Kepler-ellipse evaluation — no velocity, no time dependence (elements are already
fixed/osculating; only `ean` varies as the free parameter for the distance searches below).

## 5. `get_dist_from_2_vectors` (swecl.c:6017-6024)

`sqrt((x1_0-x2_0)^2 + (x1_1-x2_1)^2 + (x1_2-x2_2)^2)` — trivial Euclidean distance, no notes.

## 6. `osc_iterate_max_dist` / `osc_iterate_min_dist` (swecl.c:6026-6096)

Coordinate-descent hill-climb over eccentric anomaly `ean` (degrees) to find the local
max/min of `get_dist_from_2_vectors(xb, osc_get_ecl_pos(ean, pqr))`, where `xb` (the other
body's position) is held **fixed** for the duration of this single call — this is a 1-D search
over one ellipse's `ean`, not a joint 2-D search (the caller alternates calls between the two
ellipses — see §8).

Algorithm (identical structure for max/min, only comparison operator and `rmax`/`rmin` naming
differ):
```
dstep_min = high_prec ? 0.000001 : 1     // degrees
ean = 0
xa = osc_get_ecl_pos(ean, pqr)
r = dist(xb, xa)
rmax = r                                  // (rmin for the min version)
dstep = 1                                 // degrees, initial step
while dstep >= dstep_min:
    for i in 0, 1:                        // two directions: +dstep then -dstep
        while r >= rmax:                  // (r <= rmin for min version) — climb while improving
            eansv = ean
            ean += (i==0) ? dstep : -dstep
            xa = osc_get_ecl_pos(ean, pqr)
            r = dist(xb, xa)
            if r > rmax: rmax = r         // (r < rmin for min version)
        ean = eansv                       // step overshot; back off to last-good ean
        r = rmax                          // (rmin)
    ean = eansv
    r = rmax
    dstep /= 10                           // refine: 1, 0.1, 0.01, ... down to dstep_min
*drmax = rmax; *deanopt = eansv
```
This is a per-axis coordinate-descent / golden-ish step-halving(-by-10) search: at each `dstep`
scale, walk in `+` direction until the distance stops increasing (max case) or stops decreasing
(min case), then walk in `-` direction from the same starting point, then shrink `dstep` by 10×
and repeat. `high_prec = TRUE` (used exclusively by `swe_orbit_max_min_true_distance`'s geocentric
path, §8) tightens the final step to 1e-6° instead of 1°; the non-high-precision 1° floor is only
reached in the (unused-here) generic call form — in practice this file only ever calls these
helpers with `high_prec = TRUE` (swecl.c:6263-6264, 6277-6278).

## 7. `orbit_max_min_true_distance_helio` (swecl.c:6101-6128)

Heliocentric-only branch: single-ellipse case (Sun, Moon, or any explicit HELCTR/BARYCTR request).
```c
ipli = (ipl == SE_SUN) ? SE_EARTH : ipl
de = swe_get_orbital_elements(tjd_et, ipli, iflagi, serr)   // may ERR; iflagi = iflag & (EPHMASK|HELCTR|BARYCTR)
dmax = de[16]; dmin = de[15]                                 // aphelion / perihelion distances straight from Kepler elements
pqri = osc_get_orbit_constants(de)
eani = de[8]                                                  // current eccentric anomaly (osculating "now")
xinner = osc_get_ecl_pos(eani, pqri)                          // current heliocentric position from the ellipse model
dtrue = |xinner|                                              // true distance = magnitude of that position vector
```
No iterative search needed here — max/min are read directly off the ellipse's own `a(1∓e)`, and
"true" distance is just the ellipse evaluated at the body's own current eccentric anomaly (which,
being freshly derived from the same `swe_get_orbital_elements` call, exactly reproduces the
body's real instantaneous heliocentric distance).

## 8. `swe_orbit_max_min_true_distance` (swecl.c:6170-6287)

Top-level entry point. Dispatch (swecl.c:6184-6187):
```c
if (ipl == SE_SUN || ipl == SE_MOON || (iflagi & (SEFLG_HELCTR|SEFLG_BARYCTR)))
    return orbit_max_min_true_distance_helio(tjd_et, ipl, iflagi, dmax, dmin, dtrue, serr);
```
where `iflagi = iflag & (SEFLG_EPHMASK | SEFLG_HELCTR | SEFLG_BARYCTR)`. Otherwise, falls through
to the full **geocentric two-ellipse** search (the EMB's orbit vs. the target planet's orbit;
"true distance" = current planet-to-Earth distance, "max/min" = extremal possible planet-to-EMB
distance treating both bodies as free points on their respective *fixed* osculating ellipses).

### 8.1 Setup (swecl.c:6188-6205)

```c
dp = swe_get_orbital_elements(tjd_et, ipl, iflagi, serr)        // may ERR
de = swe_get_orbital_elements(tjd_et, SE_EARTH, iflagi, serr)   // may ERR (EMB elements)
(douter, dinner) = (de[0] > dp[0]) ? (de, dp) : (dp, de)         // whichever has larger semimajor axis is "outer"
pqro = osc_get_orbit_constants(douter)
pqri = osc_get_orbit_constants(dinner)
eano = douter[8]; eani = dinner[8]                                // current eccentric anomalies
xouter = osc_get_ecl_pos(eano, pqro); xinner = osc_get_ecl_pos(eani, pqri)
rtrue = dist(xouter, xinner)                                      // current true separation
```

### 8.2 Rough grid scan (swecl.c:6207-6253)

```c
ncnt = 182; dstep = 2        // degrees
for j in 0..ncnt:             // outer-ellipse eccentric anomaly, 0,2,4,...,362 deg (182 steps of 2 deg)
    eano = j * dstep
    xouter = osc_get_ecl_pos(eano, pqro)
    for i in 0..ncnt:          // inner-ellipse eccentric anomaly, 0,1,2,...,181 deg (182 steps of 1 deg — NOTE: dstep NOT applied to i, i.e. inner loop uses degrees directly, only 182 deg of coverage, not 364!)
        eani = i
        xinner = osc_get_ecl_pos(eani, pqri)
        r = dist(xouter, xinner)
        track (rmax, max_eanisv, max_eanosv, max_xouter, max_xinner) on r > rmax
        track (rmin, min_eanisv, min_eanosv, min_xouter, min_xinner) on r < rmin
```
**Important loop-bound quirk**: the outer loop (`j`, "outer" ellipse) scans `eano = j*dstep` for
`j=0..181`, i.e. `0,2,...,362°` — over-scanning past 360° by one extra 2° step (362°, equivalent
to 2°, a harmless duplicate sample). The inner loop (`i`, "inner" ellipse) sets `eani = (double)i`
directly (**not** `i*dstep`) for `i=0..181`, i.e. only `0°..181°` — **half the ellipse is never
sampled in the rough scan** (182° to 359° for the inner body is skipped). This looks like it could
be a bug (`dstep` presumably intended to apply to both loops symmetrically), but it is what the C
does — **replicate exactly**, since the subsequent per-axis refinement (§8.3) starting from
whatever rough extremum was found in this scan can still converge correctly in practice (the
extremes for realistic orbital geometries tend to occur where both ellipses' near/far points
align, which is captured regardless of this asymmetry) — do not "fix" the inner loop to `i*dstep`
without matching golden-data confirmation that C behaves identically; assume it does not, and
port the loop literally.

Initialization: `max_xouter/max_xinner/min_xouter/min_xinner` zeroed (swecl.c:6220-6225);
`rmax` seeded implicitly at 0 via first-iteration comparisons is NOT explicit — rather, `rmax`
starts uninitialized-to-any-sentinel in this function's locals at declaration (swecl.c:6178:
`rmax = 0`) and `rmin = 100000000` (swecl.c:6178) as a large sentinel, so the very first grid
point always satisfies both `r > rmax` (since real distances are ≪ 1e8 AU... wait rmin sentinel
1e8, rmax sentinel 0 — both correctly bootstrap on first real sample since any real AU-scale `r`
exceeds 0 and is below 1e8).

### 8.3 Refinement (swecl.c:6254-6282)

Two independent refinement loops (max, then min), each alternating `osc_iterate_max_dist` /
`osc_iterate_min_dist` between the two ellipses (§6) until the extremal distance stabilizes:

```c
// maximum:
(eani, eano) = (max_eanisv, max_eanosv); (xinner, xouter) = (max_xinner, max_xouter)
for k in 0..=300:                                    // nitermax = 300 (inclusive bound: 301 iterations max)
    osc_iterate_max_dist(eani, pqri, xinner, xouter, &eani, &rmax, high_prec=TRUE)  // refine inner ellipse's ean, holding xouter fixed; updates xinner in place
    osc_iterate_max_dist(eano, pqro, xouter, xinner, &eano, &rmax, high_prec=TRUE)  // refine outer ellipse's ean, holding (now-updated) xinner fixed
    if k > 0 and |rmax - rmaxsv| < 0.00000001: break   // convergence: 1e-8 AU
    rmaxsv = rmax
// minimum: identical structure with osc_iterate_min_dist / rmin
```
This is a **block coordinate ascent/descent**: alternately hold one ellipse's point fixed and
optimize the other's eccentric anomaly (via the §6 per-axis hill-climb, itself refined to
1e-6° with `high_prec=TRUE`), iterating the two-ellipse alternation up to 300 times or until the
extremal distance changes by less than 1e-8 AU between successive full (inner+outer) passes.
Each call to `osc_iterate_max_dist`/`osc_iterate_min_dist` **mutates** `xinner`/`xouter` in place
(the position arrays passed as both the "current position to compare against" and — via `xa`
inside the helper — the working buffer that ends up holding the found-optimal position), so the
next call sees the updated position immediately (true Gauss-Seidel-style block update, not
Jacobi).

### 8.4 Output (swecl.c:6283-6286)

```c
*dmax = rmax; *dmin = rmin; *dtrue = rtrue
```

## Porting Notes (stateless port)

- **No `swed.*` globals are read anywhere in this file section** (`get_gmsm` through
  `swe_orbit_max_min_true_distance`, swecl.c:5687-6287). Every input is either a function
  parameter or the result of an explicit `swe_calc` call. This module is a clean, fully stateless
  port — no cache-interaction caveats apply (contrast with the deflection-speed /
  SPEED3-at-file-boundary caveats documented in the project's `stateless_tolerance` notes, which
  come from *inside* `swe_calc`'s own state, not from this file).
- Every `swe_calc(tjd_et, ipl, iflag, x, serr)` call in this file maps to
  `Ephemeris::calc(jd_tt, body, flags) -> Result<CalcResult, Error>` (`src/context.rs:193`). Note
  the C flag-construction idiom throughout this file
  (`iflJ2000 = (iflag & SEFLG_EPHMASK) | SEFLG_J2000 | SEFLG_XYZ | SEFLG_TRUEPOS | SEFLG_NONUT | SEFLG_SPEED`,
  swecl.c:5792-5793 etc.) — port as `CalcFlags` bitor construction preserving the ephemeris-source
  bits from the caller's input flags (`flags & (CalcFlags::JPLEPH | CalcFlags::SWIEPH | CalcFlags::MOSEPH)`
  equivalent — check `src/flags.rs` for the current `CalcFlags` bit names; there is presently no
  single `EPHMASK`-equivalent constant defined there, so this port will need one, or explicit
  per-flag masking).
- `SEFLG_ORBEL_AA` is bit-aliased onto `SEFLG_TOPOCTR` (swephexp.h:207, `= SEFLG_TOPOCTR`) because
  orbital elements have no topocentric variant. **This means a caller who mistakenly passes a
  topocentric flag to the orbital-elements API silently gets AA-mass-summing behavior instead of
  an error.** Decide explicitly for the Rust port whether to (a) replicate this bit-reuse (define
  an `ORBEL_AA` flag equal to the `TOPOCTR` bit and accept the same silent-reinterpretation
  behavior — required for exact fidelity if any golden test exercises it) or (b) give `ORBEL_AA`
  its own bit and treat `TOPOCTR` on this API as a hard error (cleaner, but diverges from C
  behavior if a golden test ever passes `SEFLG_TOPOCTR` here). Recommend (a) for fidelity unless
  no golden test covers it, in which case flag clearly in code comments.
- The `plmass`/`ipl_to_elem` Pluto-indexing quirk (see Constants section above) **must be
  replicated exactly**, including the double-counted-Mercury behavior in the `SEFLG_ORBEL_AA`
  summation loop for Pluto. This is a case where "the C is buggy" is not license to fix it in the
  port — golden-data parity requires the same numeric quirk.
- `swe_get_orbital_elements` issues up to 3 real ephemeris calls per invocation (heliocentric-
  distance probe at J2000/TRUEPOS/NONUT/SPEED, the position-with-center-flags call, and — for
  `SE_EARTH` only — a geocentric Moon call to build the EMB). All three must use the *same*
  `tjd_et` and consistent `flags & SEFLG_EPHMASK` passthrough; in the stateless Rust port these
  are independent `Ephemeris::calc` calls with no risk of the C's cache-staleness issues (each
  call is fully self-contained), so no special handling is needed beyond faithfully constructing
  each call's flags per §2.2-2.3.
- `get_gmsm`'s asteroid/`SEFLG_ORBEL_AA` branch (§1, Case C) makes up to 9 additional `swe_calc`
  calls (`Ephemeris::calc`) — Mercury through Pluto plus Earth — purely to compare heliocentric
  distances (`r > x[2]`). This is expensive; if `swe_get_orbital_elements` is called for an
  asteroid with `SEFLG_ORBEL_AA` set, expect ~12 total `Ephemeris::calc` invocations for one
  orbital-elements result. No caching exists in C either (this is recomputed fresh every call) —
  faithful port needs no memoization to match C's behavior, though callers doing bulk asteroid
  queries with `ORBEL_AA` may want to add caching at a higher layer (out of scope for this doc).
- Suggested landing module: `src/orbit.rs` (new). Depends on `crate::context::Ephemeris::calc`,
  `crate::types::{Body, CalcResult}`, `crate::flags::CalcFlags`, and the shared math helpers
  `swi_cross_prod`/`swi_polcart`/`swi_cartpol`/`swi_coortrf2`/`swi_mod2PI`/`square_sum`/`dot_prod`
  — check whether these already exist as ported Rust helpers elsewhere (they are generic
  vector/coordinate utilities used across many C files, e.g. swephlib.c:141-354) before
  reimplementing; if a shared `src/coords.rs`-style module already has equivalents, reuse rather
  than duplicate (per project's shared-logic constraint).
- `Ephemeris::calc` currently returns `CalcResult { data: [f64;6], flags_used: CalcFlags }`
  (`src/context.rs:1563-1566`) with no separate polar/cartesian distinction beyond what the
  request flags encode (`CalcFlags::XYZ` toggles interpretation of `data`) — this matches the C's
  `x[6]` output array convention used throughout this file (`x[0..2]` = position, `x[3..5]` =
  velocity, in whatever polar/cartesian + center + equinox frame the request flags specified).
- FP-sensitive expressions to quote/preserve literally when porting (see inline call-outs above
  for full context): the two-step `ecce = pp/sema; ecce = sqrt(1-ecce)` reassignment (§2.6); the
  `fac = sqrt((1-ecce)*(1+ecce))` grouping in `osc_get_orbit_constants` (§3, not `1-ecce*ecce`);
  the `manom` factor ordering `ecce * RADTODEG * sin(eanom*DEGTORAD)` (§2.8); the rate-polynomial
  coefficients in `ysid`/`ytrop` that are pre-differentiated (`2*2.0441*T`, `3*0.07207*T2`, etc.,
  §2.11) — do not re-derive these from an integrated position polynomial, transcribe the rate
  form directly.
- The inner/outer eccentric-anomaly grid-scan loop-bound asymmetry in
  `swe_orbit_max_min_true_distance` (§8.2: outer loop steps `j*2°` for 182 steps reaching 362°,
  inner loop steps `i*1°` for 182 steps reaching only 181°) is almost certainly an unintentional
  C bug (probably meant `eani = i * dstep` to mirror the outer loop), but **must be ported
  literally** — the subsequent 300-iteration refinement pass (§8.3) is what actually determines
  the returned precision, and changing the rough-scan coverage could shift which local
  extremum the refinement converges to, breaking golden-data parity for edge cases (e.g. bodies
  with two comparably-deep distance minima).
