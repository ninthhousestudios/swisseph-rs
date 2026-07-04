# C Reference: House Systems — swehouse.c

Porting reference for the house-systems subsystem. Read this instead of the C source.
The Rust implementer never needs to open `swehouse.c` or `swehouse.h` directly.
`swe_gauquelin_sector` lives in `swecl.c`, not `swehouse.c`; it is documented here too
since it is the public Gauquelin-sector entry point and structurally belongs with houses.

All line numbers below refer to `swehouse.c` unless prefixed `swecl.c:`.

---

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_houses` | swehouse.c:130–175 | Yes — UT entry, tropical only |
| `swe_houses_ex` | swehouse.c:178–187 | Yes — thin wrapper around `swe_houses_ex2` |
| `swe_houses_ex2` | swehouse.c:207–290 | Yes — primary UT entry, sidereal + speed support |
| `swe_houses_armc` | swehouse.c:590–599 | Yes — thin wrapper around `swe_houses_armc_ex2` |
| `swe_houses_armc_ex2` | swehouse.c:622–774 | Yes — THE DRIVER: ARMC-based entry, calls `CalcH` (+ finite-diff speed pass) |
| `sidereal_houses_ecl_t0` | swehouse.c:318–403 | Yes — sidereal mode "ecliptic of t0" projection |
| `sidereal_houses_ssypl` | swehouse.c:425–532 | Yes — sidereal mode "solar system plane" projection |
| `sidereal_houses_trad` | swehouse.c:535–587 | Yes — sidereal mode "traditional" (subtract ayanamsa) |
| `CalcH` | swehouse.c:892–2050 | Yes — THE CORE: computes `struct houses` for one armc/lat/obliquity/system |
| `Asc1` | swehouse.c:2058–2088 | Yes — quadrant-normalized oblique-ascension trig |
| `Asc2` | swehouse.c:2100–2129 | Yes — raw oblique spherical trig (x ∈ [0,90]) |
| `AscDash` | swehouse.c:2133–2147 | Yes — analytical derivative of `Asc1` w.r.t. armc |
| `swi_armc_to_mc` | swehouse.c:872–888 | Already ported (`src/...` — used by fixed-star ayanamsa `GALCENT_MULA_WILHELM`) |
| `armc_to_mc` (static) | swehouse.c:2149–2166 | Yes — near-duplicate of `swi_armc_to_mc`, used only inside `swe_house_pos` |
| `fix_asc_polar` | swehouse.c:2169–2177 | Yes — keeps Asc on eastern hemisphere near poles |
| `apc_sector` | swehouse.c:782–825 | Yes — APC house system per-cusp formula |
| `swe_house_name` | swehouse.c:827–859 | Yes — house-system char → display name |
| `swe_house_pos` | swehouse.c:2216–2876 | Yes — planet → house-position number, per-system inverse geometry |
| `sunshine_init` | swehouse.c:2878–2904 | Yes — shared ascensional-difference setup for both Sunshine solutions |
| `sunshine_solution_makransky` | swehouse.c:2906–3046 | Yes — Sunshine houses, Makransky method (`hsys='i'`) |
| `sunshine_solution_treindl` | swehouse.c:3048–3143 | Yes — Sunshine houses, Treindl method (`hsys='I'`) |
| `swe_gauquelin_sector` | swecl.c:6298–6428 | Yes — Gauquelin sector lookup; imeth 0/1 geometric, imeth 2–5 rise/set (separate module) |

---

## 1. Header Types & Macros (swehouse.h)

### `struct houses` (swehouse.h:61–84)

```c
struct houses {
    double cusp[37];          // cusp[1..12] normally; cusp[1..36] for Gauquelin ('G')
    double cusp_speed[37];    // daily-motion speed of each cusp, degrees/day
    double ac;                 double ac_speed;
    double mc;                 double mc_speed;
    double armc_speed;         // = ARMCS constant, always (sidereal rate)
    double vertex;              double vertex_speed;
    double equasc;               double equasc_speed;    // "equatorial ascendant"
    double coasc1;                double coasc1_speed;   // co-ascendant (W. Koch)
    double coasc2;                double coasc2_speed;   // co-ascendant (M. Munkasey)
    double polasc;                double polasc_speed;   // polar ascendant (M. Munkasey)
    double sundec;             // declination of Sun, only used by Sunshine houses ('I'/'i')
    AS_BOOL do_speed;          // compute ac/mc/armc/vertex/equasc/coasc1/coasc2/polasc speeds
    AS_BOOL do_hspeed;         // compute cusp_speed[1..12] (subset of do_speed)
    AS_BOOL do_interpol;       // OUTPUT flag: this house system needs finite-diff cusp speeds
    char serr[AS_MAXCH];
};
```

`cusp[0]` and `cusp_speed[0]` are unused (always set to 0 by the caller `swe_houses_armc_ex2`,
swehouse.c:662–663); cusps are 1-indexed, matching the public `ascmc`/`cusp` array convention
used throughout the C library.

`do_speed` vs `do_hspeed`: `do_speed` gates the 8 "special point" speeds (ac, mc, armc, vertex,
equasc, coasc1, coasc2, polasc). `do_hspeed` gates the 12 (or 36) `cusp_speed[]` entries
specifically. `swe_houses_armc_ex2` sets `do_speed = (ascmc_speed != NULL || cusp_speed != NULL)`
and `do_hspeed = (cusp_speed != NULL)` (swehouse.c:644–647) — i.e. `do_hspeed` implies `do_speed`,
but a caller can request `ascmc_speed` without `cusp_speed`, getting only the special-point speeds.

### Constants

```c
#define VERY_SMALL          1E-10                          // generic epsilon
#define VERY_SMALL_PLAC_ITER (1.0 / 360000.0)               // Placidus/Gauquelin Newton convergence, ≈0.01 arcsec
#define MILLIARCSEC          (1.0 / 3600000.0)               // 1 mas in degrees; used as a "nudge into the house" epsilon
#define SOLAR_YEAR            365.24219893
#define ARMCS  ((SOLAR_YEAR+1) / SOLAR_YEAR * 360)            // sidereal rotation rate, ≈360.985647366 °/day (deg/solar day)
```
`ARMCS` is the rate at which the celestial sphere (and hence ARMC) advances per mean solar day —
the standard 360.985647... sidereal constant. It is used directly as `armc_speed` (constant for
every call) and as the multiplier in `AscDash`.

### Degree-wrapped trig macros (swehouse.h:89–98)

```c
#define degtocs(x)  (d2l((x) * DEG))
#define cstodeg(x)  (double)((x) * CS2DEG)
#define sind(x)   sin((x) * DEGTORAD)
#define cosd(x)   cos((x) * DEGTORAD)
#define tand(x)   tan((x) * DEGTORAD)
#define asind(x)  (asin(x) * RADTODEG)
#define acosd(x)  (acos(x) * RADTODEG)
#define atand(x)  (atan(x) * RADTODEG)
#define atan2d(y, x) (atan2(y, x) * RADTODEG)
```
Every angle in `swehouse.c` is in **degrees**, end to end (unlike most of the rest of the C
library, which works in radians internally). The Rust port should keep the same convention for
this module's internal functions, converting to/from radians only at the public-API boundary if
the Rust API differs. `degtocs`/`cstodeg` (centiseconds ↔ degrees) are not used anywhere in
`swehouse.c` and can be ignored.

---

## 2. Core Trig

### `Asc1(x1, f, sine, cose)` — quadrant-normalized oblique ascension (swehouse.c:2058–2088)

Computes where a great circle of "pole height" `f` (e.g. the horizon at geographic latitude `f`,
or a Placidus/Regiomontanus-family auxiliary circle) crosses the ecliptic, given the point `x1`
where that circle crosses the equator (measured as a right-ascension-like angle from the vernal
point). `sine`/`cose` are `sind(ekl)`/`cosd(ekl)` of the ecliptic obliquity, passed in pre-computed.

```c
static double Asc1(double x1, double f, double sine, double cose) {
  x1 = swe_degnorm(x1);
  int n = (int)((x1 / 90) + 1);          // quadrant 1..4 of x1 in [0,360)
  if (fabs(90 - f) < VERY_SMALL) return 180;   // f ≈ +90 (north pole): result pinned to 180°
  if (fabs(90 + f) < VERY_SMALL) return 0;     // f ≈ -90 (south pole): result pinned to 0°
  double ass;
  if (n == 1)      ass = Asc2(x1, f, sine, cose);
  else if (n == 2) ass = 180 - Asc2(180 - x1, -f, sine, cose);
  else if (n == 3) ass = 180 + Asc2(x1 - 180, -f, sine, cose);
  else /* n==4 */  ass = 360 - Asc2(360 - x1, f, sine, cose);
  ass = swe_degnorm(ass);
  // snap to the cardinal points if within VERY_SMALL — guards rounding noise at fi=0,st=0 etc.
  if (fabs(ass - 90) < VERY_SMALL)  ass = 90;
  if (fabs(ass - 180) < VERY_SMALL) ass = 180;
  if (fabs(ass - 270) < VERY_SMALL) ass = 270;
  if (fabs(ass - 360) < VERY_SMALL) ass = 0;
  return ass;
}
```

The quadrant-folding trick reduces every call to `Asc2` to its valid domain `x ∈ [0,90]`, with the
sign of `f` flipped for quadrants 2 and 3 (reflection symmetry of the oblique-ascension triangle).
This is the single workhorse used for nearly every house cusp in the file — always called as
`Asc1(<right-ascension-like angle>, <pole height>, sine, cose)`.

### `Asc2(x, f, sine, cose)` — raw spherical trig, x ∈ [0,90] (swehouse.c:2100–2129)

Solves the spherical triangle formed by the great circle (pole height `f`), the ecliptic
(obliquity `e`), and the equator, via the cotangent four-parts formula (cited in the C comment
as Wikipedia "Spherical trigonometry CT5"): `cot c sin a = cot C sin B + cos a cos B`, with
`a = x`, `B = e` (ecliptic obliquity), `C = 90° + f` (so `cot C = -tan f`).

```c
static double Asc2(double x, double f, double sine, double cose) {
  double ass = -tand(f) * sine + cose * cosd(x);
  if (fabs(ass) < VERY_SMALL) ass = 0;
  double sinx = sind(x);
  if (fabs(sinx) < VERY_SMALL) sinx = 0;
  if (sinx == 0) {
    ass = (ass < 0) ? -VERY_SMALL : VERY_SMALL;     // degenerate: x≈0, avoid atand(0/0)
  } else if (ass == 0) {
    ass = (sinx < 0) ? -90 : 90;                     // degenerate: denominator≈0
  } else {
    ass = atand(sinx / ass);                          // tan c = sin x / (cot c)^-1... i.e. c = atan2-like via atand
  }
  if (ass < 0) ass = 180 + ass;
  return ass;
}
```
Note this is **not** `atan2d(sinx, ass)` — it is `atand(sinx/ass)` followed by a manual
`+180` if negative, which is mathematically equivalent to `atan2d` restricted to the
[0,180) branch (since the two degenerate branches above already special-case the axis
crossings). Replicate the exact branch structure, not a simplified `atan2d` call, to match
FP rounding at the `VERY_SMALL` boundaries bit-for-bit.

### `AscDash(x, f, sine, cose)` — analytical derivative of `Asc1` w.r.t. armc (swehouse.c:2131–2147)

Computes `d(Asc1(x,f,...))/d(armc)`, assuming `x = armc + const` (i.e. `dx/d(armc) = 1`), scaled
by the sidereal rate `ARMCS` to convert "per degree of x" into "per day". Contributed by Graham
Dawson; comment cites `ARMCS ≈ 360.985647366`.

```c
static double AscDash(double x, double f, double sine, double cose) {
  double cosx = cosd(x), sinx = sind(x);
  double sinx2 = sinx * sinx;
  double c = cose * cosx - tand(f) * sine;
  double d = sinx2 + c * c;
  double dudt = (d > VERY_SMALL) ? (cosx * c + cose * sinx2) / d : 0.0;  // on ecliptic axis: speed = 0
  return dudt * ARMCS;
}
```
This is the analytical workhorse for nearly all house-cusp and special-point speeds. It must be
called with **exactly the same `(x, f)` pair** used for the corresponding `Asc1` position call —
the derivative is only valid for that specific (x,f) parametrization (e.g. Koch's `f=fi` constant
across quadrant but `x` varies with `ad3`; Campanus's `f` varies per cusp).

### `swi_armc_to_mc(armc, eps)` — ARMC → MC ecliptic longitude (swehouse.c:872–888)

Already ported (used by fixed-star ayanamsa `GALCENT_MULA_WILHELM`; see `c-ref-fixstar.md`
§"GALCENT_MULA_WILHELM"). Reproduced here because **three near-duplicate implementations of this
formula exist** in `swehouse.c` with subtly different normalization — see FP-fidelity §1 below.

```c
double swi_armc_to_mc(double armc, double eps) {
  double tant, mc;
  if (fabs(armc - 90) > VERY_SMALL && fabs(armc - 270) > VERY_SMALL) {
    tant = tand(armc);
    mc = atand(tant / cosd(eps));
    if (armc > 90 && armc <= 270) mc = swe_degnorm(mc + 180);
    // NOTE: no further swe_degnorm() here — mc can be returned un-normalized (see §FP-fidelity 1)
  } else {
    mc = (fabs(armc - 90) <= VERY_SMALL) ? 90 : 270;
  }
  return mc;
}
```

### `fix_asc_polar(asc, armc, eps, geolat)` (swehouse.c:2169–2177)

Used only inside `swe_house_pos`, to flip the ascendant to the opposite (eastern) hemisphere when
it has fallen onto the wrong side near the poles (mirrors the `acmc < 0` polar-swap logic used
throughout `CalcH`, but expressed via the approximate declination of the MC):

```c
static double fix_asc_polar(double asc, double armc, double eps, double geolat) {
  double demc = atand(sind(armc) * tand(eps));   // approx declination of MC
  if (geolat >= 0 && 90 - geolat + demc < 0) asc = swe_degnorm(asc + 180);
  if (geolat < 0  && -90 - geolat + demc > 0) asc = swe_degnorm(asc + 180);
  return asc;
}
```

---

## 3. Public Entry Points & Call Graph

```
swe_houses(tjd_ut, geolat, geolon, hsys, cusp, ascmc)
  └─ computes armc from tjd_ut+ΔT, mean obliquity, nutation, sidereal time
  └─ if hsys=='I': swe_calc_ut(SE_SUN, SEFLG_SPEED|SEFLG_EQUATORIAL) → ascmc[9] = Sun declination
  └─ swe_houses_armc_ex2(armc, geolat, eps+nutlo[1], hsys, cusp, ascmc, NULL, NULL, NULL)

swe_houses_ex(tjd_ut, iflag, geolat, geolon, hsys, cusp, ascmc)
  └─ swe_houses_ex2(..., cusp_speed=NULL, ascmc_speed=NULL, serr=NULL)

swe_houses_ex2(tjd_ut, iflag, geolat, geolon, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr)
  └─ same armc/eps/nutation setup as swe_houses, but eps via swi_epsiln(tjde, 0) [iflag=0 always]
  └─ if iflag & SEFLG_NONUT: nutlo[0]=nutlo[1]=0
  └─ if hsys=='I': swe_calc_ut(SE_SUN,...) → ascmc[9]; on failure, fall back to hsys='O' (Porphyry)
  └─ if iflag & SEFLG_SIDEREAL:
        dispatch to sidereal_houses_ecl_t0 / sidereal_houses_ssypl / sidereal_houses_trad
        based on swed.sidd.sid_mode bits (SE_SIDBIT_ECL_T0 / SE_SIDBIT_SSY_PLANE / neither)
     else:
        swe_houses_armc_ex2(armc, geolat, eps+nutlo[1], hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr)
  └─ if iflag & SEFLG_RADIANS: convert cusp[1..ito] and ascmc[0..SE_NASCMC) to radians

swe_houses_armc(armc, geolat, eps, hsys, cusp, ascmc)
  └─ swe_houses_armc_ex2(..., cusp_speed=NULL, ascmc_speed=NULL, serr=NULL)

swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr)  ← THE DRIVER
  └─ armc = swe_degnorm(armc)
  └─ h.do_speed = (ascmc_speed != NULL || cusp_speed != NULL)
  └─ h.do_hspeed = (cusp_speed != NULL)
  └─ if hsys=='I': resolve h.sundec from ascmc[9] (== 99 sentinel ⇒ use cached saved_sundec; else
                    ascmc[9] itself, cached into a `static double saved_sundec` for next call)
                    validate -24 <= sundec <= 24, else ERR
  └─ CalcH(armc, geolat, eps, hsys, &h)             ← fills cusps, ac, mc, vertex, etc.
  └─ copy h.cusp[1..ito] → cusp[]; h.ac/mc/armc/vertex/equasc/coasc1/coasc2/polasc → ascmc[0..7]
  └─ ascmc[8] is unused (always 0); ascmc[9] = h.sundec if hsys=='I'
  └─ if h.do_interpol: finite-difference cusp_speed via two extra CalcH(armc∓darmc,...) calls
                        (see §4 Speed Mechanism)
  └─ return retc (OK or ERR, with serr set by CalcH on Placidus/Koch/Gauquelin polar failure)

swe_house_pos(armc, geolat, eps, hsys, xpin, serr)   ← planet→house-position inverse
swe_house_name(hsys)                                  ← char → display string, pure lookup
swe_gauquelin_sector(t_ut, ipl, starname, iflag, imeth, geopos, atpress, attemp, dgsect, serr)
                                                       ← in swecl.c; geometric (imeth 0/1) or
                                                         rise/set-based (imeth 2-5)
```

### `CalcH` top-level structure (swehouse.c:892–2050)

`CalcH(th, fi, ekl, hsy, hsp)` — `th`=ARMC (sidereal time, degrees), `fi`=geographic latitude,
`ekl`=ecliptic obliquity, `hsy`=house-system char. This is **the** function that fills a
`struct houses`. Top-level flow:

1. **Always-run setup** (swehouse.c:942–987), regardless of house system:
   - `cose=cosd(ekl); sine=sind(ekl); tane=tand(ekl)`.
   - Pole clamp: if `|fi|` is within `VERY_SMALL` of 90°, clamp `fi` to `±(90 - VERY_SMALL)`.
   - `tanfi = tand(fi)`.
   - **MC**: closed-form via the `tand(th)/cose` formula (same algebra as `swi_armc_to_mc`, but
     inlined with an unconditional final `swe_degnorm` — see FP-fidelity §1). `mc_speed =
     AscDash(th, 0, sine, cose)` if `do_speed` (note: pole height `f=0` for MC — MC sits on the
     equator's projection, not a latitude-dependent circle).
   - **Ascendant**: `ac = Asc1(th+90, fi, sine, cose)` — the horizon's pole height equals
     geographic latitude; the horizon crosses the equator 90° east of the meridian.
     `ac_speed = AscDash(th+90, fi, sine, cose)` if `do_speed`.
   - If `do_hspeed`: zero `cusp_speed[0..12]` (default, may be overwritten per-system below).
   - `armc_speed = ARMCS` always (constant).
   - `cusp[1] = ac; cusp[10] = mc` (and `cusp_speed[1]=ac_speed; cusp_speed[10]=mc_speed` if
     `do_hspeed`) — pre-seeded; most systems leave these alone, several explicitly re-set them.
   - **Deprecated lowercase handling**: if `hsy > 95` (i.e. lowercase ASCII) and `hsy != 'i'`
     (the one lowercase code that is NOT deprecated — Sunshine/Makransky), set a warning in
     `hsp->serr` and uppercase it: `hsy -= 32`.
2. **`switch (hsy)`** — one case per house system; fills `cusp[1..12]` (or `cusp[1..36]` for `'G'`)
   and (where analytically available) `cusp_speed[1..12]`. See §6 below for every case.
   Several branches `goto porphyry` on polar/iteration failure (label at the `'O'` case,
   swehouse.c:1311).
3. **Post-switch opposite-cusp mirror** (swehouse.c:1985–2000): for every system **except**
   `'G'`, `'Y'`, and `toupper(hsy)=='I'` (i.e. except Gauquelin, APC, and both Sunshine variants —
   these already fill all 12 cusps themselves with independent geometry, not point-opposite
   pairs):
   ```c
   cusp[4]=degnorm(cusp[10]+180); cusp[5]=degnorm(cusp[11]+180); cusp[6]=degnorm(cusp[12]+180);
   cusp[7]=degnorm(cusp[1]+180);  cusp[8]=degnorm(cusp[2]+180);  cusp[9]=degnorm(cusp[3]+180);
   if (do_hspeed && !do_interpol) {
     cusp_speed[4]=cusp_speed[10]; cusp_speed[5]=cusp_speed[11]; cusp_speed[6]=cusp_speed[12];
     cusp_speed[7]=cusp_speed[1];  cusp_speed[8]=cusp_speed[2];  cusp_speed[9]=cusp_speed[3];
   }
   ```
   This mirror runs **even for systems whose `cusp_speed[1..3,10..12]` were never explicitly set
   in the switch** (e.g. `'N'`, `'W'`, `'B'`) — see §4 for the precise (sometimes surprising)
   resulting values.
4. **Special points, always computed** (swehouse.c:2001–2049): vertex, equasc, coasc1, coasc2,
   polasc — see §7 below. These run unconditionally after the switch, for every house system
   including `'G'`/`'Y'`/`'I'`/`'i'`.

---

## 4. THE SPEED MECHANISM

This is the most important structural fact for slicing the Rust port: **speed computation is not
uniform.** There are three distinct mechanisms, and which one applies depends on the house system
and even on which specific cusp.

### 4.1 Special-point speeds (`ascmc_speed[]`) — ALWAYS analytical, ALWAYS computed

`ac_speed`, `mc_speed`, `armc_speed`, `vertex_speed`, `equasc_speed`, `coasc1_speed`,
`coasc2_speed`, `polasc_speed` are computed unconditionally in `CalcH` (gated only by
`hsp->do_speed`, never by house system) via `AscDash` calls with the exact same `(x, f)` argument
pairs used for the corresponding position formula (§7). They are **never** subject to
finite-difference fallback. `armc_speed` is simply the constant `ARMCS`.

### 4.2 Cusp speeds (`cusp_speed[1..12]`) — per-system: analytical | finite-difference | neither

**(a) Analytical via `AscDash`**, computed inline in the `switch` at the same `(x,f)` used for
the position: `'C'` (Campanus), `'H'` (Horizon), `'J'` (Savard-A), `'K'` (Koch, closed-form, no
iteration), `'R'` (Regiomontanus), `'T'` (Polich/Page), `'G'` (Gauquelin, **at the converged
pole height `f` after Newton iteration** — speed is analytical-at-the-converged-point, not
finite-differenced, despite the position itself being iteratively solved), and the default
(`'P'` Placidus, same as Gauquelin — analytical at converged `f`).

**(b) Analytical but NOT via `AscDash`** — `'O'` (Porphyry): cusp speeds are linear
interpolations of the quadrant growth rate, not derivatives of `Asc1`:
```c
q1_speed = ac_speed - mc_speed;                       // rate of growth of quadrant 1
cusp_speed[2]  = ac_speed - q1_speed / 3;
cusp_speed[3]  = ac_speed - q1_speed / 3 * 2;
cusp_speed[11] = ac_speed + q1_speed / 3;
cusp_speed[12] = ac_speed + q1_speed / 3 * 2;
```
(swehouse.c:1326–1333)

**(c) Equal/Vehlow/EqualMC families** — `'A'`/`'E'`/`'V'`: all 12 `cusp_speed[i] = ac_speed`.
`'D'` (Equal-MC): all 12 `cusp_speed[i] = mc_speed`. These are exact (the houses are rigid
30°-spaced rotations of a single moving point, so every cusp has the same angular speed).

**(d) Finite-difference (central difference at the `swe_houses_armc_ex2` driver level)** — the
house system sets `hsp->do_interpol = hsp->do_hspeed;` inside its `switch` case (NOT inside
`CalcH`'s shared code) for: `'I'`/`'i'` (Sunshine, both), `'L'` (Pullen SD), `'Q'` (Pullen SR),
`'S'` (Sripati), `'X'` (Meridian), `'M'` (Morinus), `'F'` (Carter), `'Y'` (APC). For these,
`CalcH` itself does NOT compute `cusp_speed` (leaves it zeroed); instead `swe_houses_armc_ex2`
(swehouse.c:697–723), after seeing `h.do_interpol == TRUE`, makes **two additional full `CalcH`
calls**:
```c
double dt = 1.0 / 86400;          // 1 second, in days
double darmc = dt * ARMCS;        // the armc delta corresponding to that 1 second
hm1 = CalcH(armc - darmc, geolat, eps, hsys, &hm1);   // do_speed=FALSE, do_hspeed=FALSE for these calls
hp1 = CalcH(armc + darmc, geolat, eps, hsys, &hp1);
if (both succeed) {
  // 90°-wrap guard: if the Asc jumped >90° between hp1/hm1 and the center, the interval
  // straddles a degenerate boundary (e.g. polar Asc flip). Fall back to one-sided difference
  // and halve dt.
  if (|difdeg2n(hp1.ac, h.ac)| > 90) { hp1 = h; dt /= 2; }
  else if (|difdeg2n(hm1.ac, h.ac)| > 90) { hm1 = h; dt /= 2; }
  for (i = 1; i <= 12; i++)
    cusp_speed[i] = difdeg2n(hp1.cusp[i], hm1.cusp[i]) / 2 / dt;
}
```
This is a genuine 3-point (or, with the wrap guard, 2-point one-sided) central-difference
derivative of `CalcH`'s own cusp output w.r.t. armc, evaluated **only on `cusp[1..12]`** — never
on `ascmc[]`. Note: this block runs only `if (h.do_interpol)`, and `do_interpol` itself is only
ever set to `TRUE` when `do_hspeed` was already `TRUE` (it's assigned `= hsp->do_hspeed`), so it
is a no-op pass-through when speeds weren't requested at all.

**(e) NEITHER analytical NOR finite-difference — left at zero (or stale pre-switch value)** —
this is a genuine quirk of the C source that golden tests must replicate bit-exactly:
- `'B'` (Alcabitius): the case body never touches `cusp_speed` analytically but DOES set
  `do_interpol = do_hspeed` (swehouse.c:1621), so cusp speeds are computed via the
  driver-level finite-difference path when speeds are requested.
- `'N'` (Equal/1=Aries) and `'W'` (Whole Sign): neither sets `cusp_speed` nor `do_interpol`,
  even though `cusp[1]` (and for `'N'`, every cusp) is reassigned to a value unrelated to `ac`.
  Result: `cusp_speed[1]=ac_speed` (stale — semantically wrong, since whole-sign cusps are a step
  function of `ac`, not `ac` itself), `cusp_speed[10]=mc_speed`, mirrored to `[4]` and `[7]`;
  **all other `cusp_speed[i] = 0`**. Replicate exactly — do not "fix" this in the Rust port.

### 4.3 Summary table

| System | Cusp-speed mechanism |
|---|---|
| P (default), G | Analytical `AscDash` at converged pole height |
| C, H, J, K, R, T | Analytical `AscDash`, closed-form pole height |
| O | Analytical, linear quadrant-rate interpolation |
| A, E, V | All cusps = `ac_speed` |
| D | All cusps = `mc_speed` |
| I, i, L, Q, S, X, M, F, Y | Finite difference (driver-level, 1-second armc step) |
| B, N, W | `cusp_speed[1]`/`[10]`/(mirrored `[4]`/`[7]`) = stale ac/mc speed; rest = 0 |

---

## 5. Per-House-System Algorithms

Each subsection: char code, formula, iteration parameters (if any), polar-circle handling,
line numbers. `th`=armc, `fi`=geolat, `ekl`=obliquity, `sine/cose`=`sind/cosd(ekl)`,
`tane`=`tand(ekl)`, `tanfi`=`tand(fi)`. All angles in degrees.

### A / E — Equal houses (swehouse.c:994–1010)
```
acmc = difdeg2n(ac, mc)
if acmc < 0: ac = degnorm(ac+180); cusp[1] = ac        // polar AC/DC swap
for i=2..12: cusp[i] = degnorm(cusp[1] + (i-1)*30)
cusp_speed[i] = ac_speed for all i (if do_hspeed)
```

### D — Equal, begin at MC (swehouse.c:1011–1027)
```
acmc = difdeg2n(ac, mc); if acmc<0: ac = degnorm(ac+180)   // (cusp[1] NOT reassigned here)
cusp[10] = mc
cusp[11..12] = degnorm(cusp[10] + (i-10)*30)
cusp[1..9]   = degnorm(cusp[10] + (i+2)*30)                 // i.e. cusp[1]=mc+90, cusp[9]=mc+330
cusp_speed[i] = mc_speed for all i (if do_hspeed)
```

### C — Campanus (swehouse.c:1028–1082)
Prime vertical divided into 30° arcs; great circles from north to south celestial pole through
those points intersect the ecliptic. Pole heights via two right-triangle relations:
```
fh1 = asind(sind(fi) / 2)                     // pole height for the 30°/150° points
fh2 = asind(sqrt(3)/2 * sind(fi))              // pole height for the 60°/120° points
cosfi = cosd(fi)
if cosfi == 0: xh1 = xh2 = (fi>0 ? 90 : 270)
else:
  xh1 = atand(sqrt(3) / cosfi)                 // RA offset for 30°/150° points
  xh2 = atand(1/sqrt(3) / cosfi)                // RA offset for 60°/120° points
cusp[11] = Asc1(th+90-xh1, fh1, sine, cose)
cusp[12] = Asc1(th+90-xh2, fh2, sine, cose)
cusp[2]  = Asc1(th+90+xh2, fh2, sine, cose)
cusp[3]  = Asc1(th+90+xh1, fh1, sine, cose)
cusp_speed[11/12/2/3] = AscDash with the same (x,f) pairs (if do_hspeed)
```
**Polar handling** (swehouse.c:1071–1081): if `|fi| >= 90-ekl` (within polar circle) and
`difdeg2n(ac,mc) < 0`: add 180° to `ac`, `mc`, and `cusp[i]` for `i ∈ {1,2,3,10,11,12}` only
(explicitly `continue`s for `4 <= i < 10` — those are filled later by the post-switch mirror, so
shifting them here would be both wrong and redundant).

### H — Horizon / Azimuth (swehouse.c:1083–1155)
Same Campanus-style trisection of the prime vertical, but rotated 180° in `th` and with `fi`
mapped to its co-latitude first (`fi = 90-fi` if `fi>0` else `fi = -90-fi`), then mapped back at
the end. Structurally:
```
fi' = (fi>=0) ? 90-fi : -90-fi          // co-latitude, sign-preserving
clamp |fi'| away from 90 by VERY_SMALL (north/south "equator" degenerate case)
th' = degnorm(th+180)
fh1, fh2, xh1, xh2 — identical formulas to Campanus, using fi'
cusp[11]=Asc1(th'+90-xh1,fh1,..); cusp[12]=Asc1(th'+90-xh2,fh2,..)
cusp[1] =Asc1(th'+90,fi',..)              // note: cusp[1] IS recomputed here (unlike Campanus)
cusp[2] =Asc1(th'+90+xh2,fh2,..); cusp[3]=Asc1(th'+90+xh1,fh1,..)
cusp_speed[11,12,1,2,3] via AscDash, same args
// polar-circle 180° shift exactly as Campanus, on i ∈ {1,2,3,10,11,12}
// then UNCONDITIONALLY: cusp[1..3] += 180°, cusp[11..12] += 180°  (swehouse.c:1141–1144)
// restore fi, th to original values
// final AC/DC sanity check (without MC-shift this time): if difdeg2n(ac,mc)<0: ac += 180
```
The "+180 to cusp[1..3] and [11,12]" step (unconditional, after the polar branch) re-orients the
azimuth-measured cusps into the conventional ecliptic-house ordering. Replicate the exact order
of operations: polar shift happens **before** this unconditional shift, and `fi`/`th` are
restored to their original (non-co-latitude, non-+180) values only at the very end.

### I / i — Sunshine houses (Treindl / Makransky) (swehouse.c:1156–1181)
```
acmc = difdeg2n(ac, mc)
if acmc < 0:
  ac = degnorm(ac+180); cusp[1]=ac
  if !SUNSHINE_KEEP_MC_SOUTH && hsy=='I':       // SUNSHINE_KEEP_MC_SOUTH is a compile-time #define, =0
    mc = degnorm(mc+180); cusp[10]=mc
cusp[4] = degnorm(cusp[10]+180)
cusp[7] = degnorm(cusp[1]+180)
if hsy=='I': retc = sunshine_solution_treindl(th, fi, ekl, hsp)
else:        retc = sunshine_solution_makransky(th, fi, ekl, hsp)
if retc==ERR (Makransky only — Treindl never returns ERR from this call site):
  serr = "within polar circle, switched to Porphyry"; hsy='O'; goto porphyry
do_interpol = do_hspeed       // → finite-difference cusp speeds (§4.2d)
```
See §"Sunshine Houses" below for the two solution algorithms in full.

### J — Savard-A (swehouse.c:1182–1249)
"Albategnius"-style: latitude circles at `2/3·fi` and `fi/3` intersect the prime meridian.
```
sinfi=sind(fi); cosfi=cosd(fi)
if |fi| < VERY_SMALL: xs2 = 1/3; xs1 = 2/3
else: xs2 = asind(sind(fi/3)/sinfi); xs1 = asind(sind(2*fi/3)/sinfi)
if cosfi==0: xh1=xh2=(fi>0?90:270)
else: xh1 = atand(tand(xs1)/cosfi); xh2 = atand(tand(xs2)/cosfi)
fh1 = asind(sind(fi)*sind(90-xs1))
fh2 = asind(sind(fi)*sind(90-xs2))
cusp[12]=Asc1(th+90-xh2,fh2,..); cusp[11]=Asc1(th+90-xh1,fh1,..)
cusp[2] =Asc1(th+90+xh2,fh2,..); cusp[3] =Asc1(th+90+xh1,fh1,..)
cusp_speed[11,12,2,3] via AscDash, same args
// polar shift on i ∈ {1,2,3,10,11,12} exactly as Campanus, when |fi| >= 90-ekl
```

### K — Koch (swehouse.c:1250–1272)
Polar failure: `if |fi| >= 90-ekl: retc=ERR; serr="...switched to Porphyry"; goto porphyry`
(no Newton iteration — Koch fails outright rather than iterating in the polar circle).
Closed-form (no iteration) via:
```
sina = sind(mc) * sine / cosd(fi); clip to [-1,1]
cosa = sqrt(1 - sina^2)                  // always >= 0
c = atand(tanfi / cosa)
ad3 = asind(sind(c) * sina) / 3.0
cusp[11] = Asc1(th+30-2*ad3, fi, sine, cose)
cusp[12] = Asc1(th+60-ad3,   fi, sine, cose)
cusp[2]  = Asc1(th+120+ad3,  fi, sine, cose)
cusp[3]  = Asc1(th+150+2*ad3,fi, sine, cose)
cusp_speed via AscDash, same args (if do_hspeed)
```

### L — Pullen SD "sinusoidal delta" (ex Neo-Porphyry) (swehouse.c:1273–1300)
```
acmc = difdeg2n(ac,mc); if <0: ac+=180; cusp[1]=ac; acmc=difdeg2n(ac,mc)   // recompute after swap
q1 = 180 - acmc
d = (acmc - 90) / 4.0
if acmc <= 30: cusp[11]=cusp[12] = degnorm(mc + acmc/2)     // degenerate: zero-width house 11
else: cusp[11]=degnorm(mc+30+d); cusp[12]=degnorm(mc+60+3*d)
d = (q1 - 90) / 4.0
if q1 <= 30: cusp[2]=cusp[3] = degnorm(ac + q1/2)            // degenerate: zero-width house 2
else: cusp[2]=degnorm(ac+30+d); cusp[3]=degnorm(ac+60+3*d)
do_interpol = do_hspeed
```

### N — Equal, begin at 0° Aries (whole-sign zodiac) (swehouse.c:1301–1309)
```
acmc = difdeg2n(ac,mc); if <0: ac = degnorm(ac+180)   // ac flipped but NOT used for cusps
for i=1..12: cusp[i] = (i-1) * 30.0                     // fixed: 0,30,60,...,330 — ignores ac/mc entirely
```
No `cusp_speed` handling — see §4.2(e).

### O — Porphyry (swehouse.c:1310–1335, label `porphyry:`)
Trisects each quadrant equally (the universal polar-circle fallback target for P/K/G/I).
```
acmc = difdeg2n(ac,mc); if <0: ac=degnorm(ac+180); cusp[1]=ac; acmc=difdeg2n(ac,mc)
cusp[1]=ac; cusp[10]=mc      // re-asserted (may have been clobbered by a failed Gauquelin attempt)
cusp[2] = degnorm(ac + (180-acmc)/3)
cusp[3] = degnorm(ac + (180-acmc)/3*2)
cusp[11]= degnorm(mc + acmc/3)
cusp[12]= degnorm(mc + acmc/3*2)
// cusp_speed: see §4.2(b) — linear quadrant-rate interpolation, NOT AscDash
```

### Q — Pullen SR "sinusoidal ratio" (swehouse.c:1336–1380)
Solves a cubic-derived ratio `r` via a closed-form (Cardano-style) expression, no iteration.
```
third = 1/3; two23 = (2*2)^third  // 2^(2/3)
acmc = difdeg2n(ac,mc); if <0: ac+=180; cusp[1]=ac; acmc=difdeg2n(ac,mc)
q = acmc; if q>90: q = 180-q
if q < 1e-30: x=xr=xr3=0; xr4=180          // degenerate quadrant
else:
  c = (180-q)/q
  csq = c*c
  ccr = (csq - c)^third                     // cuberoot(c²-c)
  cqx = sqrt(2^(2/3)*ccr + 1)
  r1 = 0.5*cqx
  r2 = 0.5*sqrt(-2*(1-2*c)/cqx - two23*ccr + 2)
  r = r1 + r2 - 0.5
  x = q/(2*r+1); xr = r*x; xr3 = xr*r*r; xr4 = xr3*r
if acmc > 90:
  cusp[11]=degnorm(mc+xr3); cusp[12]=degnorm(cusp[11]+xr4)
  cusp[2] =degnorm(ac+xr);  cusp[3] =degnorm(cusp[2]+x)
else:
  cusp[11]=degnorm(mc+xr); cusp[12]=degnorm(cusp[11]+x)
  cusp[2] =degnorm(ac+xr3);cusp[3] =degnorm(cusp[2]+xr4)
do_interpol = do_hspeed
```
Note the `acmc > 90` branch swaps which of `{xr,x}` vs `{xr3,xr4}` apply to the MC-side vs
AC-side houses — replicate exactly, including the multiplication-grouping order (`xr3 = xr*r*r`,
not `xr*r^2` or `r*r*xr` — same value but match the literal expression for FP fidelity).

### R — Regiomontanus (swehouse.c:1381–1409)
```
fh1 = atand(tanfi * 0.5)
fh2 = atand(tanfi * cosd(30))
cusp[11]=Asc1(30+th,fh1,..); cusp[12]=Asc1(60+th,fh2,..)
cusp[2] =Asc1(120+th,fh2,..); cusp[3]=Asc1(150+th,fh1,..)
cusp_speed via AscDash, same args
// polar shift on i ∈ {1,2,3,10,11,12} when |fi| >= 90-ekl, identical pattern to Campanus
```

### S — Sripati (swehouse.c:1410–1431)
Uses Porphyry sector boundaries, then takes the **midpoint** of each sector as the cusp:
```
acmc = difdeg2n(ac,mc); if <0: ac=degnorm(ac+180); acmc=difdeg2n(ac,mc)
q1 = 180 - acmc; s1 = q1/3.0; s4 = acmc/3.0
cusp[1] = degnorm(ac - s4*0.5)
cusp[2] = degnorm(ac + s1*0.5)
cusp[3] = degnorm(ac + s1*1.5)
cusp[10]= degnorm(mc - s1*0.5)
cusp[11]= degnorm(mc + s4*0.5)
cusp[12]= degnorm(mc + s4*1.5)
do_interpol = do_hspeed
```
Note `cusp[1]` and `cusp[10]` are reassigned here (offset from `ac`/`mc` by half a sector), unlike
most systems where `cusp[1]=ac, cusp[10]=mc` exactly.

### T — Polich/Page "topocentric" (swehouse.c:1432–1458)
Structurally identical to Regiomontanus, but with `tanfi/3` and `tanfi*2/3` pole heights instead
of Regiomontanus's `tanfi/2` and `tanfi*cos(30°)`:
```
fh1 = atand(tanfi/3.0); fh2 = atand(tanfi*2.0/3.0)
cusp[11]=Asc1(30+th,fh1,..); cusp[12]=Asc1(60+th,fh2,..)
cusp[2] =Asc1(120+th,fh2,..); cusp[3]=Asc1(150+th,fh1,..)
cusp_speed via AscDash
// polar shift: when |fi| >= 90-ekl: ac+=180, mc+=180, AND for i=1..12: cusp[i]+=180 (ALL 12,
//   not the {1,2,3,10,11,12}-only subset used by C/H/J/R). Cusps 4-9 are not yet meaningfully
//   set at this point (overwritten later by the post-switch mirror), so this is harmless but
//   structurally different from the other quadrant-trisection systems — replicate literally.
```

### V — Vehlow (Equal, ac-15°) (swehouse.c:1459–1473)
```
acmc = difdeg2n(ac,mc); if <0: ac=degnorm(ac+180)
cusp[1] = degnorm(ac - 15)
for i=2..12: cusp[i] = degnorm(cusp[1] + (i-1)*30)
cusp_speed[i] = ac_speed for all i (if do_hspeed)
```

### W — Whole Sign (swehouse.c:1474–1484)
```
acmc = difdeg2n(ac,mc); if <0: ac=degnorm(ac+180); cusp[1]=ac
cusp[1] = ac - fmod(ac, 30)            // snap down to the sign boundary at/below ac
for i=2..12: cusp[i] = degnorm(cusp[1] + (i-1)*30)
```
No `cusp_speed` handling — see §4.2(e). Note `fmod(ac, 30)` (not `swe_degnorm`-wrapped) — `ac`
is already in `[0,360)` from `Asc1`, so `fmod` here is just "distance past the last 30° boundary".

### X — Meridian / axial rotation (swehouse.c:1485–1516)
Ecliptic points whose right ascension is `armc + n·30°`:
```
a = th
for i=1..12:
  j = i+10; if j>12: j -= 12
  a = degnorm(a + 30)
  if |a-90|>VERY_SMALL && |a-270|>VERY_SMALL:
    cusp[j] = atand(tand(a)/cose); if a>90 && a<=270: cusp[j] = degnorm(cusp[j]+180)
  else: cusp[j] = (|a-90|<=VERY_SMALL) ? 90 : 270
  cusp[j] = degnorm(cusp[j])
acmc = difdeg2n(ac,mc); if <0: ac = degnorm(ac+180)
do_interpol = do_hspeed
```
This is literally `swi_armc_to_mc(a, eps)` applied 12 times at `a = th+30, th+60, ..., th+360`,
with results placed into a rotated index `j` (cusp 11 gets the first iteration, cusp 12 the
second, ..., cusp 10 the twelfth — `j = ((i+10-1) mod 12) + 1`).

### M — Morinus (swehouse.c:1517–1540)
Same `armc + n·30` equatorial points as X, but projected onto the ecliptic via a full
`swe_cotrans` (equatorial→ecliptic rotation by `ekl`) instead of the `tand/cose` shortcut:
```
a = th
for i=1..12:
  j = i+10; if j>12: j -= 12
  a = degnorm(a+30)
  x = [a, 0]; swe_cotrans(x, x, ekl)         // equatorial (RA=a, Dec=0) → ecliptic
  cusp[j] = x[0]
acmc = difdeg2n(ac,mc); if <0: ac = degnorm(ac+180)
do_interpol = do_hspeed
```
The `swe_cotrans` call here uses `+ekl` (vs. `-ekl` used elsewhere for ecliptic→equatorial) —
confirm sign convention against the existing `swe_cotrans`/`cotrans` Rust port (`src/math.rs`).

### F — Carter "poli-equatorial" (swehouse.c:1541–1580)
RA of the ascendant is the starting point; cusps are great circles through `(ascendant_RA +
(n-1)·30°)` on the equator and the celestial poles, intersected with the ecliptic.
```
acmc = difdeg2n(ac,mc); if <0: ac=degnorm(ac+180); cusp[1]=ac
x=[ac,0]; swe_cotrans(x,x,-ekl); a = x[0]    // a = RA of ascendant
for i=2..12:
  if i<=3 || i>=10:                           // only houses 2,3,10,11,12 computed this way
    ra = degnorm(a + (i-1)*30)
    if |ra-90|>VERY_SMALL && |ra-270|>VERY_SMALL:
      cusp[i] = atand(tand(ra)/cose); if ra>90&&ra<=270: cusp[i]=degnorm(cusp[i]+180)
    else: cusp[i] = (|ra-90|<=VERY_SMALL) ? 90 : 270
    cusp[i] = degnorm(cusp[i])
  // i in {4..9}: left unset here, filled by post-switch mirror
do_interpol = do_hspeed
```
Note: only `i ∈ {2,3,10,11,12}` are computed in the loop (the `if i<=3||i>=10` guard); `cusp[10]`
was already `mc` from before the switch and is not recomputed. `cusp[1]` is `ac` (already set).

### B — Alcabitius (swehouse.c:1581–1622)
Comment: "created by Alois 17-sep-2000, followed example in Matrix electrical library... This
corresponds to Munkasey's 'The Alcabitius Semiarc House System'."
```
acmc = difdeg2n(ac,mc); if <0: ac=degnorm(ac+180); cusp[1]=ac; acmc=difdeg2n(ac,mc)
dek = asind(sind(ac) * sine)                   // declination of Ascendant
r = -tanfi * tand(dek); clip r to [-1,1]
sda = acosd(r)                                  // semidiurnal arc of Asc, on equator
sna = 180 - sda                                  // seminocturnal arc (complement)
sd3 = sda/3; sn3 = sna/3
cusp[11] = Asc1(degnorm(th + sd3),       0, sine, cose)
cusp[12] = Asc1(degnorm(th + 2*sd3),     0, sine, cose)
cusp[2]  = Asc1(degnorm(th + 180 - 2*sn3), 0, sine, cose)
cusp[3]  = Asc1(degnorm(th + 180 - sn3),   0, sine, cose)
do_interpol — NOT set; cusp_speed left untouched — see §4.2(e)
```
Pole height is always **0** here — the RA offsets are projected onto the ecliptic along the
declination circle, not a latitude-dependent great circle.

### G — 36 Gauquelin sectors (swehouse.c:1623–1730)
Polar failure: `if |fi| >= 90-ekl: retc=ERR; serr="...switched to Porphyry"; hsy='O'; goto porphyry`.
Counted **clockwise**. Two mirrored Newton-iteration loops (4th/2nd quarter, then 1st/3rd
quarter), each producing 8 of the 36 sectors plus their 180°-opposite partners:
```
a = asind(tand(fi) * tane)                       // max declination amplitude factor

// 4th/2nd quarter: ih = 2..9 (sector indices), ih2 = 10-ih
for ih = 2..9:
  ih2 = 10 - ih
  fh1 = atand(sind(a*ih2/9) / tane)                       // initial pole-height guess
  rectasc = degnorm((90/9)*ih2 + th)
  tant = tand(asind(sine * sind(Asc1(rectasc, fh1, sine, cose))))   // seed iteration
  if |tant| < VERY_SMALL:
    cusp[ih] = rectasc; cusp_speed[ih] = armc_speed       // degenerate: cusp coincides with AC/DC axis
  else:
    f = atand(sind(asind(tanfi*tant)*ih2/9) / tant)
    cusp[ih] = Asc1(rectasc, f, sine, cose)
    cuspsv = 0
    for i = 1..niter_max (=100):
      tant = tand(asind(sine * sind(cusp[ih])))
      if |tant| < VERY_SMALL: cusp[ih]=rectasc; cusp_speed[ih]=armc_speed; break
      f = atand(sind(asind(tanfi*tant)*ih2/9) / tant)
      cusp[ih] = Asc1(rectasc, f, sine, cose)
      if i>1 && |difdeg2n(cusp[ih], cuspsv)| < VERY_SMALL_PLAC_ITER: break    // converged
      cuspsv = cusp[ih]
    if i >= niter_max: retc=ERR; hsy='O'; serr="very close to polar circle, switched to Porphyry"; goto porphyry
    cusp_speed[ih] = AscDash(rectasc, f, sine, cose)        // analytical, AT CONVERGED f
  cusp[ih+18] = degnorm(cusp[ih] + 180)                     // opposite sector
  cusp_speed[ih+18] = cusp_speed[ih]

// 1st/3rd quarter: ih = 29..36, ih2 = ih-28 — mirror-image formulas:
//   fh1 = atand(sind(a*ih2/9)/tane)
//   rectasc = degnorm(180 - ih2*90/9 + th)
//   (same Newton loop structure as above)
//   cusp[ih-18] = degnorm(cusp[ih]+180); cusp_speed[ih-18] = cusp_speed[ih]

cusp[1]=ac; cusp[10]=mc; cusp[19]=degnorm(ac+180); cusp[28]=degnorm(mc+180)
cusp_speed[1]=ac_speed; cusp_speed[10]=mc_speed; cusp_speed[19]=ac_speed; cusp_speed[28]=mc_speed
```
**Newton iteration**: `niter_max = 100`, convergence test `|difdeg2n(cusp[ih], cuspsv)| <
VERY_SMALL_PLAC_ITER` (≈0.01 arcsec) starting from the 2nd iteration (`i>1` guard avoids
false-converging on the very first step). Non-convergence after 100 iterations → `ERR`, fallback
to Porphyry (`goto porphyry`, which re-enters the `'O'` case using the *already-corrupted*
`cusp[1]`/`cusp[10]` — hence Porphyry's `cusp[1]=ac; cusp[10]=mc;` re-assertion lines exist
specifically to repair this).

### U — Krusinski-Pisa-Goelzer (swehouse.c:1731–1805)
Written by Bogdan Krusinski (2006). Defines a great circle through the Ascendant and the Zenith
(always 90° apart), divides it into 12 equal 30° arcs, and projects each division point back onto
the ecliptic via meridian (declination) circles. Implemented as a sequence of `swe_cotrans`
rotations (ecliptic↔equatorial↔horizontal↔house-circle), not a closed-form trig formula:
```
acmc = difdeg2n(ac,mc); if <0: ac = degnorm(ac+180)
x = [ac, 0, 1]
swe_cotrans(x,x,-ekl)                  // A1: ecliptic → equatorial
x[0] -= (th - 90)                       // A2: rotate by RA of east point
swe_cotrans(x,x,-(90-fi))              // A3: equatorial → horizontal
krHorizonLon = x[0]                     // save asc longitude on horizon
x[0] -= x[0]                            // A4: rotate to 0 (no-op, x[0]=0)
swe_cotrans(x,x,-90)                   // A5: horizontal → asc-zenith house-circle frame
for i = 0..5:                           // 6 iterations cover cusps 1,2,3,4,5,6 (7-12 via +180)
  x = [30*i, 0]
  swe_cotrans(x,x,90)                   // B1: house-circle → horizontal
  x[0] += krHorizonLon                   // B2: rotate back
  swe_cotrans(x,x,90-fi)                // B3: horizontal → equatorial
  x[0] = degnorm(x[0] + (th-90))         // B4: rotate back → RA of house cusp
  cusp[i+1] = atand(tand(x[0])/cosd(ekl))   // B5: equatorial → ecliptic (MC-style projection)
  if x[0]>90 && x[0]<=270: cusp[i+1] = degnorm(cusp[i+1]+180)
  cusp[i+1] = degnorm(cusp[i+1])
  cusp[i+7] = degnorm(cusp[i+1]+180)
```
No iteration/Newton loop (despite the multi-step rotation chain, it is fully closed-form); no
explicit polar-circle guard in `CalcH` itself (the C comment notes only the exact poles ±90° and
points exactly on the arctic circle where horizon=ecliptic are singular, and that these are
naturally excluded as limits). No `cusp_speed` analytic formula — `'U'` is not listed among the
`do_interpol` systems either, so (like B/N/W) its `cusp_speed[2,3,5,6,8,9,11,12]` end up `0` and
`[1,4,7,10]` carry the stale pre-switch `ac_speed`/`mc_speed`. Double-check this empirically
against golden swetest output before committing the Rust port to this assumption.

### Y — APC houses (swehouse.c:1806–1829, helper `apc_sector` at swehouse.c:782–825)
```
for i=1..12: cusp[i] = apc_sector(i, fi*DEGTORAD, ekl*DEGTORAD, th*DEGTORAD)
cusp[10] = mc                            // overwrites apc_sector's house-10 value with the real MC
                                          // (comment: "MC provided by apc_sector() near lat 90 is not accurate")
cusp[4] = degnorm(mc + 180)
// polar shift on ALL i=1..12 (not a subset) when |fi| >= 90-ekl and difdeg2n(ac,mc)<0
do_interpol = do_hspeed
```
`apc_sector(n, ph, e, az)` (radians in, degrees out) — `ph`=lat, `e`=obliquity, `az`=armc:
```c
if |ph| > 90°-VERY_SMALL (rad): kv = 0; dasc = 0
else:
  kv = atan(tan(ph)*tan(e)*cos(az) / (1 + tan(ph)*tan(e)*sin(az)))   // ascensional difference of Asc
  if |ph| < VERY_SMALL: dasc = ±(90°-VERY_SMALL) (rad, sign of ph)
  else: dasc = atan(sin(kv) / tan(ph))             // declination of the Ascendant
is_below_hor = (n < 8); k = is_below_hor ? n-1 : n-13
a = is_below_hor
    ? kv + az + π/2 + k*(π/2 - kv)/3
    : kv + az + π/2 + k*(π/2 + kv)/3
a = radnorm(a)
dret = atan2(tan(dasc)*tan(ph)*sin(az) + sin(a),
             cos(e)*(tan(dasc)*tan(ph)*cos(az) + cos(a)) + sin(e)*tan(ph)*sin(az-a))
return degnorm(dret * RADTODEG)
```
Note `apc_sector` works in **radians** throughout (unlike everything else in this file) — it is
called with `fi*DEGTORAD, ekl*DEGTORAD, th*DEGTORAD` and internally uses plain `tan`/`atan`/
`atan2`, not the `tand`/`atand` degree macros. The final `atan2` argument order matters for sign
fidelity — replicate literally, do not reorder.

### Default — Placidus (swehouse.c:1830–1983)
Polar failure: `if |fi| >= 90-ekl: retc=ERR; serr="...switched to Porphyry"; goto porphyry`
(no iteration attempted in the polar circle — fails immediately, same as Koch and Gauquelin).
Four independent Newton loops, one per cusp (11, 12, 2, 3), each structurally identical to a
single Gauquelin loop above but with fixed fractional divisors instead of `ih2/9`:

| Cusp | Initial pole-height seed `fh` | RA seed | Iteration divisor |
|---|---|---|---|
| 11 | `fh1 = atand(sind(a/3)/tane)` | `degnorm(30+th)` | `/3` |
| 12 | `fh2 = atand(sind(a*2/3)/tane)` | `degnorm(60+th)` | `/1.5` |
| 2  | `fh2` (same as 12) | `degnorm(120+th)` | `/1.5` |
| 3  | `fh1` (same as 11) | `degnorm(150+th)` | `/3` |

where `a = asind(tand(fi)*tane)` (shared, computed once). Per cusp (example for house 11):
```
rectasc = degnorm(30+th)
tant = tand(asind(sine * sind(Asc1(rectasc, fh1, sine, cose))))
if |tant| < VERY_SMALL: cusp[11]=rectasc; cusp_speed[11]=armc_speed
else:
  f = atand(sind(asind(tanfi*tant)/3) / tant)
  cusp[11] = Asc1(rectasc, f, sine, cose)
  cuspsv = 0
  for i=1..100:
    tant = tand(asind(sine*sind(cusp[11])))
    if |tant|<VERY_SMALL: cusp[11]=rectasc; cusp_speed[11]=armc_speed; break
    f = atand(sind(asind(tanfi*tant)/3) / tant)
    cusp[11] = Asc1(rectasc, f, sine, cose)
    if i>1 && |difdeg2n(cusp[11],cuspsv)| < VERY_SMALL_PLAC_ITER: break
    cuspsv = cusp[11]
  if i>=100: retc=ERR; serr="very close to polar circle, switched to Porphyry"; goto porphyry
  cusp_speed[11] = AscDash(rectasc, f, sine, cose)    // at converged f
```
This is the identical Newton-iteration skeleton used by Gauquelin (§'G' above), just with `ih2/9`
replaced by the fixed divisor `3` (house 11/3 boundary) or `1.5` (house 12/2 boundary). A shared
Rust helper (e.g. `placidus_newton_cusp(rectasc, fh_init, divisor, sine, cose, tane, tanfi, ...)`)
can serve both Placidus and Gauquelin — see `<constraints>` in project CLAUDE.md re: avoiding
duplicated logic.

---

## 6. Sidereal Houses (swe_houses_ex2 dispatch)

These three functions implement the three sidereal projection modes selected by
`swed.sidd.sid_mode` bits (`SE_SIDBIT_ECL_T0`, `SE_SIDBIT_SSY_PLANE`, or neither = traditional).
All compute a tropical house set via `swe_houses_armc_ex2` at an **adjusted** armc/obliquity, then
shift every cusp by a computed ayanamsa-like offset. This logic is genuinely house-system-specific
(not delegated to `sweph.c`'s ayanamsa machinery, except for the "traditional" mode which calls
`swe_get_ayanamsa_ex`, documented in `c-ref-ayanamsa.md`).

### `sidereal_houses_trad` (swehouse.c:535–587) — traditional: subtract ayanamsa
```
ay = swe_get_ayanamsa_ex(tjde, iflag, &ay, NULL)        // see c-ref-ayanamsa.md
ihs = toupper(hsys); ihs2 = (ihs=='W') ? 'E' : ihs        // whole-sign computed as Equal, fixed after
retc = swe_houses_armc_ex2(armc, lat, eps, ihs2, cusp, ascmc, cusp_speed, ascmc_speed, serr)
for i=1..ito: cusp[i] = degnorm(cusp[i] - ay); if ihs=='W': cusp[i] -= fmod(cusp[i], 30)
if ihs=='N': cusp[i] = (i-1)*30 for all i                  // re-fixed to 0,30,60...
for i=0..SE_NASCMC (skip i==2, armc): ascmc[i] = degnorm(ascmc[i] - ay)
```
Note: `ito` is 36 if `hsys=='G'` else 12. The commented-out `- ay - nutl` (replaced by `- ay`
alone) at swehouse.c:570,583 — `nutl` is computed by the caller but unused in the final formula;
this is dead code left in place, do not port a `nutl` parameter usage here.

### `sidereal_houses_ecl_t0` (swehouse.c:318–403) — projection onto ecliptic of t0
Seven-step geometric procedure (full derivation in the C comment block at swehouse.c:292–317):
1. Take the vernal point as a unit vector on the mean ecliptic of `t0` (the ayanamsa epoch),
   convert to the equator of `t0` via `swi_coortrf(x, x, -epst0)` where `epst0 = swi_epsiln(t0,0)`.
2. Precess that point/velocity from `t0` to J2000 to `tjde` (`swi_precess`, both legs).
3. Rotate it onto the **true** equator of `tjde`: rotate by `(eps - nutlo[1])*DEGTORAD`, convert
   to polar, add `nutlo[0]*DEGTORAD` to the longitude (nutation in longitude), convert back to
   Cartesian, then rotate by `-eps*DEGTORAD`.
4. Compute the "auxiliary obliquity" `epsx` as the angle between this point's instantaneous
   orbital plane (`swi_cross_prod(x, x+3, xnorm)`) and the equator: `epsx = asin(rxy/rxyz)*RADTODEG`
   where `rxy = sqrt(xnorm[0]² + xnorm[1]²)`, `rxyz = sqrt(rxy² + xnorm[2]²)`.
5. Compute the "auxiliary vernal point" `xvpx` — the ascending-node-like intersection — via
   `fac = x[2]/x[5]` (with `x[5]` floored to `1e-15` if near zero to avoid division by zero),
   `sgn = sign(x[5])`, `xvpx[j] = (x[j] - fac*x[j+3]) * sgn` for j=0..2.
6. `dvpx = atan2-style longitude of xvpx (via swi_cartpol) * RADTODEG`; `armcx = degnorm(armc -
   dvpx)`.
7. Compute houses at `(armcx, lat, epsx)` via `swe_houses_armc_ex2`.
8. `dvpxe = acosd(swi_dot_prod_unit(x, xvpx))`, sign-flipped if `tjde < sip->t0`; subtract
   `dvpxe + sip->ayan_t0` from every cusp and every `ascmc[i]` (`i != 2`, i.e. skip armc); if
   `hsys=='N'`, re-fix all cusps to `(i-1)*30` afterward.

### `sidereal_houses_ssypl` (swehouse.c:425–532) — projection onto solar system plane
Same nine-step shape as `ecl_t0` (see the C comment at swehouse.c:405–424) but starts from the
**solar-system-plane zero point** instead of the mean-ecliptic-of-t0 vernal point, using the
fixed constants `SSY_PLANE_INCL = 1.578701° (rad)`, `SSY_PLANE_NODE_E2000 = 107.582569° (rad)`,
`SSY_PLANE_NODE = 107.58883388° (rad)` (sweph.h:291–295). It additionally computes a
J2000-referenced ayanamsa term `x00` (the solar-system-plane zero point of `t0`, transformed
through J2000 ecliptic, the SSY plane, and back) and subtracts `dvpxe + ayan_t0 + x00` from every
cusp/ascmc entry. Both this and `ecl_t0` are independently self-contained geometry — neither
delegates to `swe_get_ayanamsa_ex`.

---

## 7. Special Points (`ascmc[3..7]`)

All computed unconditionally at the end of `CalcH` (swehouse.c:2001–2049), for every house
system, using `th`/`fi`/`sine`/`cose` from the original (un-mutated, post-special-case-restore)
call arguments.

### Vertex — `ascmc[3]`
```
f = (fi >= 0) ? 90-fi : -90-fi
vertex = Asc1(th-90, f, sine, cose)
vertex_speed = AscDash(th-90, f, sine, cose)
if |fi| <= ekl:                              // "tropical latitudes" — vertex behaves like polar Asc
  vemc = difdeg2n(vertex, mc)
  if vemc > 0: vertex = degnorm(vertex + 180)
```
The C comment: "with tropical latitudes, the vertex behaves strange, in a similar way as the
ascendant within the polar circle. we keep it always on the western hemisphere." Note: only the
*position* is flipped by this guard — `vertex_speed` is NOT recomputed for the flipped value
(same shift-invariance argument as elsewhere; `AscDash` is periodic in its `x` argument with
period 360°, so a ±180° shift to the position does not require recomputing the speed).

### Equatorial Ascendant — `ascmc[4]` (`equasc`)
This is exactly the MC-formula (`swi_armc_to_mc`-style) but evaluated at `th+90` instead of `th`,
with the pole height implicitly 0 (it is purely an obliquity projection of the RA point 90° east
of the meridian — i.e. the point on the ecliptic with the same right ascension as the East Point):
```
th2 = degnorm(th+90)
if |th2-90|>VERY_SMALL && |th2-270|>VERY_SMALL:
  tant = tand(th2); equasc = atand(tant/cose)
  if th2>90 && th2<=270: equasc = degnorm(equasc+180)
else: equasc = (|th2-90|<=VERY_SMALL) ? 90 : 270
equasc = degnorm(equasc)
equasc_speed = AscDash(th+90, 0, sine, cose)
```

### Co-Ascendant (W. Koch) — `ascmc[5]` (`coasc1`)
```
coasc1 = degnorm(Asc1(th-90, fi, sine, cose) + 180)
coasc1_speed = AscDash(th-90, fi, sine, cose)
```

### Co-Ascendant (M. Munkasey) — `ascmc[6]` (`coasc2`)
```
if fi >= 0: coasc2 = Asc1(th+90, 90-fi, sine, cose); coasc2_speed = AscDash(th+90, 90-fi, sine, cose)
else:        coasc2 = Asc1(th+90, -90-fi, sine, cose); coasc2_speed = AscDash(th+90, -90-fi, sine, cose)
```

### Polar Ascendant (M. Munkasey) — `ascmc[7]` (`polasc`)
```
polasc = Asc1(th-90, fi, sine, cose)
polasc_speed = AscDash(th-90, fi, sine, cose)
```
**Note**: `polasc` uses the *exact same* `(th-90, fi)` argument pair as the pre-`+180` term inside
`coasc1`'s computation — i.e. `polasc == degnorm(coasc1 - 180)` and `polasc_speed ==
coasc1_speed` always. The C source calls `Asc1`/`AscDash` independently for each rather than
reusing the value; since both functions are pure and deterministic, the Rust port may either
replicate the two independent calls (guaranteed bit-identical) or compute `polasc` once and
derive `coasc1` from it — verify bit-exactness against golden data if choosing the latter.

---

## 8. `swe_house_pos` — Planet House Position (swehouse.c:2216–2876)

Inverse problem: given `(armc, geolat, eps, hsys, xpin=[ecl.lon, ecl.lat])`, return a continuous
house position `1.0..13.0` (integer part = house number, fractional part = position within the
house). Per the C doc comment (swehouse.c:2190–2199): **geometrically exact** inverse formulas
exist for `A/E/D/V/W/N/O/B/X/F/M/P/K/C/R/U/T/H/G`; a **simplified** "linear fraction of cusp
span" approximation is used for `Y` (APC), `L` (Pullen SD), `Q` (Pullen SR), `I`/`i` (Sunshine),
`S` (Sripati) — `S` has a dedicated formula in the `'O'/'B'/'S'` case block, `I`/`i`/`Y` share a
dedicated geometric-approximation case block, but `L`/`Q` have **no dedicated case** and fall
through to the generic `default:` bracket-and-interpolate fallback (same fallback also used for
any future/unhandled system code).

### Setup, common to all branches (swehouse.c:2231–2285)
```
hsys = toupper(hsys)
// Pre-check: call swe_houses_armc_ex2(armc,geolat,eps,hsys,...) once; if xpin exactly matches
// a cusp (within MILLIARCSEC) and xpin[1]==0 (zero latitude), return that integer house number
// directly (no further geometry needed). Also extracts dsun for hsys 'I' (from ascmc[9], using
// the ascmc[9]=99 sentinel to read the system's cached/declination value) and for hsys 'Y'
// (declination of the ascendant, via swe_cotrans(ascmc[0],0,1; -eps)).
xeq = cotrans([xpin[0], xpin[1], 1], -eps); ra=xeq[0]; de=xeq[1]
mdd = degnorm(ra - armc); mdn = degnorm(mdd + 180)
if mdd >= 180: mdd -= 360
if mdn >= 180: mdn -= 360
```
`mdd`/`mdn` ("meridian distance, diurnal/nocturnal") are in `(-180, 180]`, the signed RA offset
from the meridian. Every per-system branch below works from `ra, de, mdd, mdn` plus `xpin`.

### Per-system inverse formulas

| System(s) | Method | Key formula | Notes |
|---|---|---|---|
| `N` | trivial | `hpos = xpin[0]/30 + 1` | ignores armc/geolat entirely |
| `A,E,D,V,W` | closed-form | `asc=Asc1(armc+90,geolat,..); mc=armc_to_mc(armc,eps); asc=fix_asc_polar(...); xp0=degnorm(xpin0-asc)`; `V`: `+15`; `W`: `+fmod(asc,30)`; `D`: overridden to `degnorm(xpin0-mc-90)`; `+MILLIARCSEC` nudge; `hpos=xp0/30+1` | shared setup across all 5 |
| `O,S` | closed-form | `xp0=degnorm(xpin0-asc)+nudge`; base `hpos=1` (or `7` if `xp0>=180`, `xp0-=180`); `acmc=difdeg2n(asc,mc)`; `hpos += xp0*3/(180-acmc)` (1st half) or `3+(xp0-180+acmc)*3/acmc` (2nd half); `S` only: `hpos+=0.5`, wrap `>12→1` | `O`=Porphyry inverse, `S`=Sripati = Porphyry inverse shifted half a sector |
| `B` (Alcabitius) | closed-form | `dek=asind(sind(asc)*sine)`; `r=-tanfi*tand(dek)` clip ±1; `sda=acosd(r)`; `sna=180-sda`; 4-branch piecewise-linear map of `mdd` into `[0,360)` using `sda`/`sna` as the two semiarc spans; `hpos=degnorm(hpos-90)/30+1`; wrap `>=13→-=12` | inverse of the semidiurnal-arc trisection |
| `X` | trivial | `hpos = degnorm(mdd-90)/30 + 1` | |
| `F` | closed-form | `x0=Asc1(armc+90,geolat,..)` (fix_asc_polar'd); cotrans to RA; `hpos=degnorm(ra-x0)/30+1` | |
| `M` | closed-form | armc_to_mc-style transform of `xpin[0]` itself, then `hpos=degnorm(result-armc-90)/30+1` | |
| `K` (Koch) | closed-form, circumpolar-aware | `admc=tand(eps)*tand(geolat)*sind(armc)` clip ±1 (flags MC circumpolar); `adp` from `de` vs `geolat` (3-way: circumpolar N / circumpolar S / normal `asind(tand(geolat)*tand(de))`); `samc=90+admc`; east/west branch on sign of `mdd`; `dfac` computed, validity requires `dfac ∈ [0,2]` else **"Koch house position failed in circumpolar area"** (serr set, `hpos=0`); `hpos=xp0/30+1` | can fail (return 0) unlike most other systems |
| `C` (Campanus) | closed-form | `xeq0=degnorm(mdd-90)`; `cotrans(xeq,xp,-geolat)`; `+nudge`; `hpos=xp0/30+1` | |
| `J` (Savard-A) | closed-form + bracket-interpolate | builds 12-entry `hcusp[]` table from `xs1,xs2` (same formula as the `'J'` `CalcH` case); `a=cotrans(degnorm(mdd-90), -geolat)[0]`; bracket `a` between consecutive `hcusp[]` entries (handles "retrograde" cusp ordering via `difdeg2n(hcusp[6],hcusp[1])` sign), linear-interpolate `hpos` within the bracket | same bracket-interpolate pattern as `default:` |
| `U` (Krusinski) | closed-form, multi-rotation | clamps `geolat` away from 0; finds the asc-zenith plane's equator intersection `raaz` and obliquity-to-equator `oblaz` via two 2-step `swe_cotrans` chains; projects both Asc and planet onto that plane; `hpos=xp0/30+1` (plus a declination-circle offset `xp[1]` that is computed but never used) | `xp[1]` computation is vestigial/dead |
| `H` (Horizon) | closed-form | `xeq0=degnorm(mdd-90)`; `cotrans(xeq,xp,90-geolat)` (note: `90-geolat`, not `-geolat` as in Campanus); `+nudge`; `hpos=xp0/30+1` | mirrors Campanus with a different rotation angle |
| `R` (Regiomontanus) | closed-form | exact `mdd≈0→xp0=270`, `mdd≈±180→xp0=90`; else clamp `geolat`,`de` away from ±90; `a=tand(geolat)*tand(de)+cosd(mdd)`; `xp0=degnorm(atand(-a/sind(mdd)))`; `+180` if `mdd<0`; `+nudge`; `hpos=xp0/30+1` | |
| `I,i,Y` | approximate, shared geometric procedure | seed `xp0` via the Regiomontanus-style `a=tand(geolat)*tand(de)+cosd(mdd)` formula; determine `is_above_hor`; compute `harmc`, `darmc`, semi-diurnal-arc `sad`/`san` of `dsun` (Sun declination for `I/i`, Asc declination for `Y`); degenerate circumpolar-`dsun` fast paths; else full spherical-triangle solve for the position-line crossing, `+nudge`; `hpos=xp0/30+1` | full algorithm in §"swe_house_pos Sunshine/APC" below |
| `T` (Polich/Page) | iterative bisection | mirror below-horizon and western-hemisphere points to the canonical quadrant first; binary-search loop (not Newton): `fac` doubles each step (`fac=2,4,8,...`), adjusts `fh`/`ra0` by `±tanfi/fac`/`±90/fac`, recomputes residual `xp[1]` via `cotrans`, stops when `|xp[1]|<1e-6` or `nloop>=1000`; mirror back; `hpos=degnorm(hpos-90)/30+1` | distinct from `CalcH`'s Newton loops — this is a true bisection on the pole-height parameter |
| `P,G` (Placidus/Gauquelin) | closed-form, "Otto Ludwig" circumpolar fallback | circumpolar fast path (`90-|de|<=|geolat|`): `xp0=degnorm(90+mdn/2)` or `degnorm(270+mdd/2)` depending on `sign(de*geolat)`, with a serr warning; else: `sinad=tand(de)*tand(geolat)`, `ad=asind(sinad)`, `is_above_hor` test, `sad=90+ad`, `san=90-ad`, `xp0 = is_above_hor ? (mdd/sad+3)*90 : (mdn/san+1)*90`; `+nudge`; `G`: `xp0=360-xp0` (clockwise), `hpos=xp0/10+1` (10°/sector, not 30°); else `hpos=xp0/30+1` | shared formula for both — `G` just rescales |
| default (anything else, incl. `L`,`Q`) | bracket-and-interpolate | calls `swe_houses_armc_ex2` to get the 12 cusps, then the same bracket/interpolate procedure as `J` above; sets `serr = "swe_house_pos(): using simplified algorithm for system %c"` | |

All branches add `MILLIARCSEC` (1 mas) to the raw fractional position before dividing by 30,
"to make sure that a call with a house cusp position returns a value within the house" (recurring
C comment) — i.e. a planet exactly on cusp N returns house N, not N-1 + 1.0 due to FP rounding.

### `swe_house_pos` Sunshine/APC shared formula (`'I','i','Y'`, swehouse.c:2650–2744)
```
clamp geolat to (-90+MILLIARCSEC, 90-MILLIARCSEC); clamp de away from ±90 by VERY_SMALL
a = tand(geolat)*tand(de) + cosd(mdd)
xp0 = degnorm(atand(-a/sind(mdd))); if mdd<0: xp0+=180; xp0=degnorm(xp0)   // "house pos with hsys='R'"
sinad = tand(de)*tand(geolat); is_above_hor = (sinad + cosd(mdd) >= 0)
harmc = (geolat<0) ? 90+geolat : 90-geolat
darmc = degnorm(xp0-270); is_western_half = darmc>180; if so: darmc = 360-darmc
sinad2 = tand(dsun)*tand(geolat)             // dsun = Sun dec (I/i) or Asc dec (Y)
ad = (sinad2>=1) ? 90 : (sinad2<=-1) ? -90 : asind(sinad2)
sad = 90+ad; san = 90-ad
if sad==0 && is_above_hor: xp0 = 270                    // circumpolar dsun, object above horizon
elif san==0 && !is_above_hor: xp0 = 90                  // circumpolar dsun, object below horizon
else:
  sa = sad
  if !is_above_hor: dsun=-dsun; sa=san; darmc=180-darmc; is_western_half=!is_western_half
  a = acosd(cosd(harmc)*cosd(darmc)); clamp away from 0 by VERY_SMALL
  sinpsi = clip(sind(harmc)/sind(a), -1, 1)
  y = sind(dsun)/sinpsi; y = (y>1)?90-VERY_SMALL:(y<-1)?-(90-VERY_SMALL):asind(y)
  d = acosd(cosd(y)/cosd(dsun)); sign-flip d by sign(dsun) then by sign(geolat)
  darmc += d
  xp0 = is_western_half ? 270 - (darmc/sa)*90 : 270 + (darmc/sa)*90
  if !is_above_hor: xp0 = degnorm(xp0+180)
xp0 += MILLIARCSEC; hpos = xp0/30 + 1
```

---

## 9. Sunshine Houses — Full Solution Algorithms

Both Sunshine variants share `sunshine_init` for the ascensional-difference setup, then diverge.

### `sunshine_init(lat, dec, xh[13])` (swehouse.c:2878–2904)
```c
arg = tand(dec) * tand(lat);
ad = (arg>=1) ? 90-VERY_SMALL : (arg<=-1) ? -90+VERY_SMALL : asind(arg);
nsa = 90 - ad;  dsa = 90 + ad;        // nocturnal / diurnal semiarc of the Sun
xh[2]=-2*nsa/3; xh[3]=-1*nsa/3; xh[5]=1*nsa/3; xh[6]=2*nsa/3;     // night-side offsets (houses 2,3,5,6)
xh[8]=-2*dsa/3; xh[9]=-1*dsa/3; xh[11]=1*dsa/3; xh[12]=2*dsa/3;   // day-side offsets (houses 8,9,11,12)
return (|arg|>=1) ? ERR : OK;          // ERR ⇒ Sun exactly circumpolar at this lat/dec
```
`xh[1,4,7,10]` are left unset (those are the cardinal cusps AC/IC/DC/MC, handled outside this
helper). Houses 2,3,5,6 use the **nocturnal** semiarc; 8,9,11,12 use the **diurnal** semiarc —
matches the convention that houses 1-6 fall below the horizon at sunrise/sunset for a sunshine
clock keyed to the Sun's own diurnal motion (not the local horizon's).

### `sunshine_solution_makransky` (`hsys='i'`, swehouse.c:2906–3046)
For each non-cardinal house `ih` (2,3,5,6,8,9,11,12 — `(ih-1)%3==0` skips 1,4,7,10):
```
md = |xh[ih]|
rah = degnorm(ramc + (ih<=6 ? 180 : 0) + xh[ih]); if lat<0: rah = degnorm(180+rah)
// zenith distance zd of the house circle, measured along the prime vertical:
if md == 90: zd = 90 - atand(sinlat*tandec)
else:
  a = (md<90) ? atand(coslat*tand(md)) : atand(tand(md-90)/coslat)
  b = atand(tanlat*cosd(md))
  c = (ih<=6) ? b+dec : b-dec
  f = atand(sinlat*sind(md)*tand(c))
  zd = a + f
pole = asind(sind(zd)*sinlat)
q = asind(tandec*tand(pole))
w = degnorm(rah + (ih<=3||ih>=11 ? -q : q))
// cu (the cusp) derived from w, pole, ecl via case analysis on w==90, w==270, else general:
//   general case computes m=atand(|tand(pole)/cosd(w)|), then z=m±ecl (sign depends on
//   ih-range and whether w∈(90,270)), then r=atand(|cosd(m)*tand(w)/cosd(z)|), then cu
//   assembled from r with a 4-way quadrant selection on w, with an additional z>90
//   "value will fall away from cancer" quadrant remap (C comment expresses uncertainty
//   about this remap's correctness — replicate verbatim, do not attempt to simplify)
// if lat<0: cu = degnorm(cu+180)
cusp[ih] = cu
```
This is the most structurally intricate single function in the file (swehouse.c:2978–3043) — a
4-to-8-way case split on the quadrant of `w` and the sign of `z-90`. **Port this function nearly
verbatim**, preserving the exact branch order and `r`/`z`/`w` quadrant-test boundaries; do not
attempt to derive a simplified closed form — the C author's own comment expresses uncertainty
about the correctness of the `z>90` remap, meaning byte-for-byte behavioral replication (not
"mathematical intent" replication) is the only safe porting strategy.

### `sunshine_solution_treindl` (`hsys='I'`, swehouse.c:3048–3143)
```
sunshine_init(lat, dec, xh)   // return value ignored here (no ERR short-circuit, unlike Makransky)
mcdec = atand(sind(ramc)*tand(ecl))
mc_under_horizon = |lat - mcdec| > 90
if mc_under_horizon && SUNSHINE_KEEP_MC_SOUTH (compile-time const, =0, so this branch is DEAD CODE):
  for ih=2..12: xh[ih] = -xh[ih]
for each non-cardinal ih (skip (ih-1)%3==0):
  xhs = 2*asind(cosdec * sind(xh[ih]/2))        // great-circle chord length x'
  cosa = tandec * tand(xhs/2); alph = acosd(cosa)
  if ih>7: alpha2 = 180-alph; b = 90-lat+dec       // diurnal side
  else:    alpha2 = alph;     b = 90-lat-dec        // nocturnal side
  cosc = cosd(xhs)*cosd(b) + sind(xhs)*sind(b)*cosd(alpha2); c = acosd(cosc)
  if c < 1e-6: serr = "Sunshine house %d c=%le very small"; retval = ERR     // does NOT early-return
  sinzd = sind(xhs)*sind(alpha2)/sind(c); zd = asind(sinzd)
  rax = atand(coslat*tand(zd))
  pole = asind(sinzd*sinlat)
  if ih<=6: pole = -pole; a = degnorm(rax+ramc+180)
  else:                    a = degnorm(ramc+rax)
  cusp[ih] = Asc1(a, pole, sinecl, cosecl)
if mc_under_horizon && !SUNSHINE_KEEP_MC_SOUTH (i.e. always, given the compile-time const):
  for ih=2..12 (skip cardinals): cusp[ih] = degnorm(cusp[ih]+180)
return retval     // ERR only if any house hit the c<1e-6 degeneracy; cusps are still filled
```
`SUNSHINE_KEEP_MC_SOUTH` is `#define`d to `0` at swehouse.c:870 with the comment "must be 0 or
1" — i.e. it is a compile-time switch the C author left in place but always builds with the
`0` branch. **Port only the `0` (MC-kept-north) behavior**; do not implement the dead `1` branch
unless a future requirement resurrects it. Note also that `sunshine_init`'s `ERR` return is
**ignored** here (no `if (...== ERR) return ERR`, unlike `sunshine_solution_makransky` which does
check it) — Treindl proceeds with the clamped `±(90-VERY_SMALL)` ascensional difference even when
the Sun is exactly circumpolar.

---

## 10. `swe_gauquelin_sector` (swecl.c:6298–6428)

```c
int32 swe_gauquelin_sector(t_ut, ipl, starname, iflag, imeth, geopos, atpress, attemp, dgsect, serr)
```
`imeth` selects one of two structurally distinct strategies:

**`imeth ∈ {0,1}` — geometric, via `swe_house_pos` with `hsys='G'`** (swecl.c:6338–6356):
```
t_et = t_ut + ΔT(t_ut, iflag)
eps = swi_epsiln(t_et, iflag)*RADTODEG; nutlo = swi_nutation(t_et, iflag)*RADTODEG
armc = degnorm(swe_sidtime0(t_ut, eps+nutlo[1], nutlo[0])*15 + geopos[0])
x0 = do_fixstar ? swe_fixstar(starname, t_et, iflag, ..) : swe_calc(t_et, ipl, iflag, ..)
if imeth==1: x0[1] = 0                          // ignore ecliptic latitude
*dgsect = swe_house_pos(armc, geopos[1], eps+nutlo[1], 'G', x0, NULL)
return OK
```
This reuses `swe_house_pos`'s `'P','G'` branch (§8) directly — `imeth=1` zeroes the body's
ecliptic latitude before the call (projects onto the ecliptic plane).

**`imeth ∈ {2,3,4,5}` — from rise/set times** (swecl.c:6357–6428): finds the most recent rising
and the next-or-most-recent setting (or vice versa, depending on which comes first) via
`swe_rise_trans` (`SE_CALC_RISE`/`SE_CALC_SET`, with `SE_BIT_NO_REFRACTION`/`SE_BIT_DISC_CENTER`
flags set per `imeth`), determines whether the body is currently above or below the horizon by
comparing `tret[0]` (rise) vs `tret[1]` (set), then interpolates the sector linearly between the
bracketing rise/set times:
```
if above_horizon: dgsect = (t_ut - t_rise) / (t_set - t_rise) * 18 + 1
else:             dgsect = (t_ut - t_set) / (t_rise - t_set) * 18 + 19
```
(36 sectors over a full diurnal+nocturnal cycle, 18 sectors per half — hence `*18`).

**Dependency note**: `swe_rise_trans` is implemented in `swehel.c`/`swecl.c`'s rise/set module,
which is **not yet ported** in `swisseph-rs` (per `docs/codebase-map.md`, no `rise_set`/`heliacal`
module exists yet — `heliacal.rs` is an empty stub). The Rust port of `swe_gauquelin_sector`
should implement `imeth ∈ {0,1}` first (depends only on `swe_house_pos` + `calc.rs` +
`fixstar2`, all already portable), and stub or defer `imeth ∈ {2,3,4,5}` until rise/set lands.

---

## 11. Sunshine Sun-Declination Input Flow

`struct houses.sundec` (swehouse.h:79) holds the Sun's declination, needed by `CalcH`'s `'I'`/`'i'`
cases (passed through to `sunshine_solution_treindl`/`_makransky` as the `dec` parameter). Flow:

1. **`swe_houses`/`swe_houses_ex2`** (the UT-based, date-aware entry points): when `hsys=='I'`,
   both compute it themselves via `swe_calc_ut(tjd_ut, SE_SUN, SEFLG_SPEED|SEFLG_EQUATORIAL, xp,
   NULL)` and stash `xp[1]` (equatorial latitude = declination) into `ascmc[9]` **before** calling
   `swe_houses_armc_ex2` (swehouse.c:145–155, 260–268). On `swe_calc_ut` failure, `swe_houses`
   falls back to Porphyry (`hsys='O'`) and returns `ERR`; `swe_houses_ex2` does the same
   (swehouse.c:149–153, 263–266).
2. **`swe_houses_armc_ex2`** (the ARMC-based, date-agnostic entry point — no access to a Sun
   ephemeris) reads `h.sundec` from the caller-supplied `ascmc[9]`:
   - If `ascmc[9] == 99` (sentinel meaning "caller didn't supply one"): falls back to a
     `static double saved_sundec` cached from the **previous** call (process-global state!), or
     `0` if no previous call cached one (swehouse.c:648–655).
   - Otherwise: uses `ascmc[9]` directly and updates the static cache for next time.
   - Validates `-24 <= sundec <= 24`, else `ERR` with `serr = "House system I (Sunshine) needs
     valid Sun declination in ascmc[9]"` (swehouse.c:656–659).
3. **`swe_house_pos`** (planet→house inverse) needs Sun declination for `hsys='I'/'i'` too, and
   obtains it the same way: it sets `ascmc[9]=99` before its own internal `swe_houses_armc_ex2`
   call (to force the static-cache fallback), then reads back `ascmc[9]` afterward as `dsun`
   (swehouse.c:2234–2257) — i.e. `swe_house_pos` **relies on the static cache having been
   populated by a prior `swe_houses_armc_ex2`/`swe_houses_ex2` call in the same session**. This
   is fragile, order-dependent global state.

**For the Rust port**: `Ephemeris` is stateless (no mutable cache — see project `CLAUDE.md`
`<architecture>`). The `static double saved_sundec` mechanism must **not** be replicated as
global/thread-local mutable state. Instead:
- The Rust equivalent of `swe_houses_ex2`/`swe_houses` should compute Sun declination internally
  (via the already-ported `calc.rs` Sun pipeline) and pass it explicitly through the call chain —
  no sentinel value needed.
- The Rust equivalent of `swe_houses_armc_ex2` (date-agnostic) should take Sun declination as an
  explicit, **required** parameter (e.g. `Option<f64>` that is `Some` and validated when
  `hsys==Sunshine`, returning an error/`None` cusps when `hsys==Sunshine` and the parameter is
  `None` — never silently falling back to a stale cached value).
- The Rust equivalent of `swe_house_pos` for `hsys==Sunshine`/`APC` must take Sun (or ascendant)
  declination as an explicit parameter rather than relying on a prior call's side effect.

---

## 12. FP-Fidelity Hazards

### 1. Three near-duplicate `armc→mc` implementations with different normalization
- `swi_armc_to_mc` (swehouse.c:872–888, public): final result is **not** unconditionally
  `swe_degnorm`'d — only the `+180` branch is. Can return a value outside `[0,360)` (e.g.
  negative) when `armc ∈ [0,90) ∪ (270,360)` and `atand` returns negative.
- `CalcH`'s inline MC computation (swehouse.c:956–968): identical formula, but with an
  unconditional `hsp->mc = swe_degnorm(hsp->mc);` appended after the if/else — always normalized.
- `armc_to_mc` (static, swehouse.c:2149–2166, used only by `swe_house_pos`): wraps the `atand`
  result in `swe_degnorm` **before** the conditional `+180`, then wraps again after — i.e.
  `mc = swe_degnorm(atand(...)); if (...) mc = swe_degnorm(mc+180);`.

  All three give the same value modulo 360°, but bit-exact golden comparison against C requires
  matching each call site's specific normalization sequence (a `swe_degnorm` on an
  already-in-range value is a no-op bit-for-bit, but on an out-of-range value it changes the
  bit pattern). When porting, implement one shared `armc_to_mc` core, but apply normalization at
  each Rust call site exactly as the corresponding C call site does — do not silently
  "improve" `swi_armc_to_mc` to always normalize, since the fixed-star ayanamsa dispatch
  (`GALCENT_MULA_WILHELM`, see `c-ref-fixstar.md`) immediately re-normalizes its result anyway,
  so this only matters if a future caller relies on `swi_armc_to_mc`'s raw, possibly-unnormalized
  return value.

### 2. `cusp_speed` for B / N / U / W — not analytically meaningful, replicate exactly
See §4.2(e). These systems leave most `cusp_speed[]` entries at `0` (the pre-switch default) or
carrying a stale `ac_speed`/`mc_speed` value unrelated to the actual (reassigned) cusp position.
Do not attempt to derive "correct" speeds for these — golden tests compare against the literal C
output, which is internally inconsistent for these four systems by construction.

### 3. `apc_sector` works in radians; everything else works in degrees
`apc_sector(n, ph, e, az)` is called with `fi*DEGTORAD, ekl*DEGTORAD, th*DEGTORAD` and uses plain
`tan`/`atan`/`atan2`/`sin`/`cos` (not the `tand`/`atand`/`sind`/`cosd` degree-macros used
everywhere else in the file). Converting the whole function to consistently use the degree-macro
style would change which operations round at which step — **port `apc_sector` working in radians
internally**, converting only its inputs (from degree-domain callers) and its output (back to
degrees via `* RADTODEG` inside the final `degnorm`).

### 4. `Asc2`'s degenerate-branch order matters
`Asc2` checks `sinx == 0` (after a `VERY_SMALL` snap-to-zero) **before** checking `ass == 0`, and
both checks happen before the general `atand(sinx/ass)` branch. Near `x=0` or `x=180` (where
`sinx→0`) and simultaneously `ass→0` (which happens near `f=0, ekl=0` or other specific
parameter combinations), the order of these checks determines which `±VERY_SMALL`/`±90`
fallback fires. Replicate the exact `if/elif/else` order, not a mathematically-equivalent
reordering.

### 5. Sripati's `hpos += 0.5` wrap (`swe_house_pos`, swehouse.c:2336–2339)
```c
if (hsys == 'S') {
  hpos += 0.5;
  if (hpos > 12) hpos = 1;
}
```
This is **not** `if (hpos > 13) hpos -= 12` (the usual wrap pattern) — it wraps straight to `1`,
discarding the fractional excess, when `hpos > 12` (not `>= 13`). Since `hpos` is continuous, this
means any `hpos ∈ (12, 12.5]` (which can occur right at the wrap boundary, since `hpos` was
already `<= 13` before the `+0.5`) collapses to exactly `1.0`, not `13.0` or `1.0 + (hpos-12.5)`.
Replicate this exact (slightly lossy) wrap behavior.

### 6. Multiplication grouping in Pullen SR (`'Q'`, swehouse.c:1364–1365)
```c
xr3 = xr * r * r;     // NOT xr * (r*r), evaluated left-to-right: (xr*r)*r
xr4 = xr3 * r;
```
Mathematically `xr*r²` either way, but left-to-right evaluation order affects the last-bit
rounding. Match the literal expression grouping.

### 7. Iteration convergence test uses `swe_difdeg2n`, not raw subtraction
Every Newton loop (Placidus, Gauquelin) tests `|swe_difdeg2n(cusp[ih], cuspsv)| <
VERY_SMALL_PLAC_ITER`, which normalizes the difference into `(-180, 180]` before taking the
absolute value — important near the 0°/360° wrap boundary. A raw `|cusp[ih] - cuspsv|` would give
a near-360° "difference" for a converged value that happens to straddle 0°, falsely failing to
converge (or converging on the wrong iteration count, changing which `f`/`tant` bit pattern is
captured for the subsequent `AscDash` speed call).

### 8. Global/static state to NOT replicate
- `static double saved_sundec` in `swe_houses_armc_ex2` (swehouse.c:636) — see §11.
- Everything else in this module is a pure function of its arguments; no other hidden state.

---

## 13. Constants Reference

| Name | Value | Source | Usage |
|---|---|---|---|
| `VERY_SMALL` | `1e-10` | swehouse.h:87 | generic epsilon (pole snaps, degenerate-axis guards) |
| `VERY_SMALL_PLAC_ITER` | `1/360000 ≈ 2.78e-6°` (0.01″) | swehouse.c:891 | Placidus/Gauquelin Newton convergence test |
| `MILLIARCSEC` | `1/3600000 ≈ 2.78e-7°` | swehouse.c:68 | "nudge into house" epsilon for `swe_house_pos` |
| `SOLAR_YEAR` | `365.24219893` | swehouse.c:69 | tropical year length, days |
| `ARMCS` | `(366.24219893/365.24219893)*360 ≈ 360.985647366` | swehouse.c:70 | sidereal rotation rate, °/day; = `armc_speed` always |
| `niter_max` | `100` | swehouse.c:940 | Placidus/Gauquelin Newton iteration cap |
| `dt` (finite-diff cusp speed) | `1/86400` day (1 second) | swehouse.c:698 | central-difference step for systems I,i,L,Q,S,X,M,F,Y |
| `SSY_PLANE_INCL` | `1.578701°` (rad) | sweph.h:295 | solar-system-plane sidereal projection |
| `SSY_PLANE_NODE_E2000` | `107.582569°` (rad) | sweph.h:291 | solar-system-plane sidereal projection |
| `SSY_PLANE_NODE` | `107.58883388°` (rad) | sweph.h:293 | solar-system-plane sidereal projection |
| `SE_NASCMC` | `8` | swephexp.h:172 | count of `ascmc[]` core entries (0..7); `ascmc[8]` unused, `ascmc[9]` = sundec for `'I'` |
| `SE_SIDBIT_ECL_T0` | `256` | swephexp.h:223 | sidereal mode bit → `sidereal_houses_ecl_t0` |
| `SE_SIDBIT_SSY_PLANE` | `512` | swephexp.h:225 | sidereal mode bit → `sidereal_houses_ssypl` |
| `J_TO_J2000` / `J2000_TO_J` | `1` / `-1` | sweph.h:256–257 | `swi_precess` direction flags |

### House-system character → name (`swe_house_name`, swehouse.c:827–859)

| Char | Name | Char | Name |
|---|---|---|---|
| A | equal | N | equal/1=Aries |
| B | Alcabitius | O | Porphyry |
| C | Campanus | Q | Pullen SR |
| D | equal (MC) | R | Regiomontanus |
| E | equal | S | Sripati |
| F | Carter poli-equ. | T | Polich/Page |
| G | Gauquelin sectors | U | Krusinski-Pisa-Goelzer |
| H | horizon/azimut | V | equal/Vehlow |
| I | Sunshine | W | equal/whole sign |
| i | Sunshine/alt. | X | axial rotation system/Meridian |
| J | Savard-A | Y | APC houses |
| K | Koch | L | Pullen SD |
| M | Morinus | *default* | Placidus |

`swe_house_name` uppercases its input **except** lowercase `'i'` (line 830: `if (h != 'i') h =
toupper(h);`) — it does not apply `CalcH`'s separate "deprecated lowercase → warn + uppercase"
logic; any unrecognized code (after this targeted uppercase) falls to the `default: "Placidus"`
string, matching `CalcH`'s own `switch` default.

---

## References

| Source | Content |
|---|---|
| swehouse.h:61–84 | `struct houses` |
| swehouse.h:87–98 | `VERY_SMALL`, `degtocs`/`cstodeg`, `sind`/`cosd`/`tand`/`asind`/`acosd`/`atand`/`atan2d` macros |
| swehouse.c:68–70 | `MILLIARCSEC`, `SOLAR_YEAR`, `ARMCS` |
| swehouse.c:130–175 | `swe_houses` |
| swehouse.c:178–290 | `swe_houses_ex`, `swe_houses_ex2` |
| swehouse.c:292–532 | sidereal houses comment block + `sidereal_houses_ecl_t0`, `sidereal_houses_ssypl` |
| swehouse.c:535–587 | `sidereal_houses_trad` |
| swehouse.c:590–774 | `swe_houses_armc`, `swe_houses_armc_ex2` (THE DRIVER) |
| swehouse.c:782–825 | `apc_sector` |
| swehouse.c:827–859 | `swe_house_name` |
| swehouse.c:872–888 | `swi_armc_to_mc` |
| swehouse.c:892–2050 | `CalcH` (THE CORE) — all 24 house-system cases + special points |
| swehouse.c:2058–2147 | `Asc1`, `Asc2`, `AscDash` |
| swehouse.c:2149–2177 | `armc_to_mc` (static dup), `fix_asc_polar` |
| swehouse.c:2216–2876 | `swe_house_pos` (planet → house-position inverse, per-system) |
| swehouse.c:2878–2904 | `sunshine_init` |
| swehouse.c:2906–3046 | `sunshine_solution_makransky` |
| swehouse.c:3048–3143 | `sunshine_solution_treindl` |
| swecl.c:6298–6428 | `swe_gauquelin_sector` |
| sweph.h:256–257 | `J_TO_J2000`, `J2000_TO_J` |
| sweph.h:291–295 | `SSY_PLANE_NODE_E2000`, `SSY_PLANE_NODE`, `SSY_PLANE_INCL` |
| sweph.h:765–770 | `struct sid_data` (`sid_mode`, `ayan_t0`, `t0`, `t0_is_UT`) |
| swephexp.h:172 | `SE_NASCMC` |
| swephexp.h:221–235 | `SE_SIDBIT_*` (documented fully in `c-ref-ayanamsa.md`) |
| `c-ref-ayanamsa.md` | `swe_get_ayanamsa_ex` (used by `sidereal_houses_trad`); `swi_armc_to_mc` consumer (`GALCENT_MULA_WILHELM`) |
| `c-ref-fixstar.md` | `swe_fixstar`/`swe_fixstar2` (used by `swe_gauquelin_sector` for `starname != NULL`) |
