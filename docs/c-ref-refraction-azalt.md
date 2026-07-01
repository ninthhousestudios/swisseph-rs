# C Reference: Atmospheric Refraction + Horizontal Coordinates — swecl.c

Porting reference for `swe_azalt`, `swe_azalt_rev`, `swe_refrac`, `swe_refrac_extended`,
`swe_set_lapse_rate`, and their shared static helpers `calc_astronomical_refr` / `calc_dip`.
All are in `swecl.c`. Line numbers below refer to `swecl.c` unless stated otherwise.

---

## Function Map

| C function | Location | Purpose |
|---|---|---|
| `swe_azalt` | swecl.c:2788–2825 | Ecliptic/equatorial → azimuth, true altitude, apparent altitude |
| `swe_azalt_rev` | swecl.c:2839–2873 | Azimuth+true altitude → ecliptic or equatorial coordinates (inverse; no de-refraction) |
| `swe_refrac` | swecl.c:2887–2984 | Simple true↔apparent altitude refraction (Meeus formula; sea-level only, no dip) |
| `swe_set_lapse_rate` | swecl.c:2986–2989 | Sets global `const_lapse_rate` used as `swe_azalt`'s default lapse rate |
| `swe_refrac_extended` | swecl.c:3035–3115 | Elevated-observer true↔apparent altitude, with horizon dip |
| `calc_astronomical_refr` (static) | swecl.c:3124–3148 | Sinclair refraction formula, shared by `swe_refrac_extended` |
| `calc_dip` (static) | swecl.c:3158–3169 | Geometric+refractive dip of the horizon for an elevated observer |

Shared constants (from `sweph.h`, `swephexp.h`):
- `SE_LAPSE_RATE` = `0.0065` (deg K / m) — sweph.h:306. This initializes the static global
  `const_lapse_rate` (swecl.c:74: `static TLS double const_lapse_rate = SE_LAPSE_RATE;`).
- `EARTH_RADIUS` = `6378136.6` meters (AA 2006 K6) — sweph.h:282.

Calc-flag integer constants (swephexp.h:364–371):
```c
#define SE_ECL2HOR   0   /* swe_azalt: input is ecliptic coords */
#define SE_EQU2HOR   1   /* swe_azalt: input is equatorial coords */
#define SE_HOR2ECL   0   /* swe_azalt_rev: output is ecliptic coords */
#define SE_HOR2EQU   1   /* swe_azalt_rev: output is equatorial coords */
#define SE_TRUE_TO_APP  0   /* swe_refrac / swe_refrac_extended direction */
#define SE_APP_TO_TRUE  1
```
Note `SE_ECL2HOR == SE_HOR2ECL == 0` and `SE_EQU2HOR == SE_HOR2EQU == 1` — the same two integers
are reused with different meaning depending on which function they're passed to. Do not conflate
them into a single Rust enum shared between `swe_azalt` and `swe_azalt_rev`; keep two distinct
enums (or a single enum with two clearly-named directions) so the port's type system doesn't
silently allow passing an azalt calc_flag into azalt_rev.

---

## 1. `swe_azalt` — ecliptic/equatorial → azimuth/altitude (swecl.c:2788–2825)

```c
void CALL_CONV swe_azalt(
      double tjd_ut,
      int32  calc_flag,      /* SE_ECL2HOR or SE_EQU2HOR */
      double *geopos,        /* [0]=geolon(east+), [1]=geolat, [2]=height above sea, meters */
      double atpress,        /* hPa; 0 => estimate from geopos[2] */
      double attemp,         /* deg C */
      double *xin,           /* [0]=lon/RA, [1]=lat/dec, degrees (only [0..1] read) */
      double *xaz)           /* out: [0]=azimuth, [1]=true altitude, [2]=apparent altitude */
```

### Algorithm

1. **Sidereal time / ARMC**: `armc = swe_degnorm(swe_sidtime(tjd_ut) * 15 + geopos[0])`.
   `swe_sidtime` returns Greenwich apparent sidereal time in **hours**; `*15` converts to degrees,
   then `+ geopos[0]` (east geographic longitude, degrees) gives the local ARMC. This is exactly
   the same ARMC construction used by the houses module (`swe_houses_ex2`) — reuse the Rust
   equivalent rather than re-deriving it.

2. **Copy input, force homogeneous coord**: `xra[0]=xin[0]; xra[1]=xin[1]; xra[2]=1`.

3. **If `calc_flag == SE_ECL2HOR`** (input is ecliptic longitude/latitude): convert ecliptic →
   equatorial in place via `swe_cotrans(xra, xra, -eps_true)`, where `eps_true` is the **true**
   (nutated) obliquity obtained by `swe_calc(tjd_ut + swe_deltat_ex(tjd_ut,-1,NULL), SE_ECL_NUT, 0, x, NULL); eps_true = x[0]`.
   - `SE_ECL_NUT` calc returns `x[0]` = true obliquity of the ecliptic (mean obliquity + nutation
     in obliquity), `x[1]` = mean obliquity, `x[2]`/`x[3]` = nutation in longitude/obliquity — only
     `x[0]` is used here.
   - `swe_deltat_ex(tjd_ut, -1, NULL)`: ΔT for the given UT with iflag=-1 (auto ephemeris-flag
     selection for ΔT/tidal-acceleration purposes — see the houses module's ΔT note for the same
     pattern).
   - `swe_cotrans(xpo, xpn, eps)` convention (swephlib.c:223–240): rotates about the x-axis by
     angle `eps` (**degrees**), treating `(xpo[0],xpo[1])` as (longitude, latitude) in a polar
     representation; **ecliptic→equatorial requires `eps` negative**, equatorial→ecliptic requires
     `eps` positive (per its own doc comment). `swe_azalt` calls it with `-eps_true` for
     ecl→equ, consistent with that convention.
   - If `calc_flag == SE_EQU2HOR`, this whole step is skipped — `xra` is already equatorial
     (RA, dec) and used as-is.

4. **Local hour angle**: `mdd = swe_degnorm(xra[0] - armc)` (RA minus ARMC = negative hour angle,
   i.e. `mdd` grows as the object moves west... actually `xra[0]` is RA so `mdd` here is used
   directly as an intermediate "meridian distance" angle, not literally HA — see next step).

5. **Rotate to horizontal-like intermediate frame**:
   ```c
   x[0] = swe_degnorm(mdd - 90);
   x[1] = xra[1];              /* declination unchanged */
   x[2] = 1;
   swe_cotrans(x, x, 90 - geolat);   /* "azimuth from east, counterclock" per the C comment */
   ```
   The rotation angle is `90 - geopos[1]` (co-latitude), applied via the same generic
   `swe_cotrans` polar-rotation routine (reused here as an equatorial→horizontal transform, not
   its documented ecliptic/equatorial use — the routine is a generic axis rotation).

6. **Re-orient azimuth to "from south, clockwise via west"**:
   ```c
   x[0] = swe_degnorm(x[0] + 90);
   xaz[0] = 360 - x[0];
   xaz[1] = x[1];    /* true (unrefracted) height */
   ```
   Comment at line 2815: "azimuth from south to west". The final `xaz[0] = 360 - x[0]` flips the
   counterclockwise-from-east convention of step 5 into the public convention: **azimuth measured
   from south, increasing clockwise (i.e. south→west→north→east)**. This is the documented public
   convention for `swe_azalt`'s `xaz[0]`.

7. **Pressure default**: if `atpress == 0`, estimate standard atmosphere pressure at the
   observer's height:
   ```c
   atpress = 1013.25 * pow(1 - 0.0065 * geopos[2] / 288, 5.255);
   ```
   Constants: sea-level pressure `1013.25` hPa, fixed lapse rate `0.0065` K/m (hardcoded here,
   **not** `const_lapse_rate` — this pressure-estimate formula always uses the literal `0.0065`,
   independent of what `swe_set_lapse_rate` was last called with), reference temperature `288` K,
   exponent `5.255` (standard barometric formula `g·M/(R·L)` for the international standard
   atmosphere). `geopos[2]` is the observer's height above sea level in meters.
   FP-fidelity: match the literal expression grouping `1013.25 * pow(1 - 0.0065 * geopos[2] / 288, 5.255)` exactly.

8. **Apparent altitude via refraction**:
   ```c
   xaz[2] = swe_refrac_extended(x[1], geopos[2], atpress, attemp, const_lapse_rate, SE_TRUE_TO_APP, NULL);
   ```
   Uses the **global** `const_lapse_rate` (see §2/§4 below) as the lapse rate — this is the one
   place in the public API where the lapse-rate global is actually consumed.
   STATELESS PORT NOTE: `swe_azalt`'s Rust port must accept `lapse_rate` as an explicit parameter
   (defaulting to `0.0065`), since there is no mutable global in the stateless design. Do not read
   a shared/static default from a `swe_set_lapse_rate`-style setter.
   The commented-out line 2824 (`/* xaz[2] = swe_refrac_extended(xaz[2], ..., SE_APP_TO_TRUE, NULL); */`)
   is dead code — do not port it; it would have been a round-trip sanity check, not part of the
   live computation.

### Summary of xaz[] outputs
- `xaz[0]` — azimuth, degrees, measured from south, positive clockwise via west (south=0, west=90, north=180, east=270).
- `xaz[1]` — true (geometric, unrefracted) altitude, degrees.
- `xaz[2]` — apparent (refracted) altitude, degrees, via `swe_refrac_extended(..., SE_TRUE_TO_APP, ...)`.

---

## 2. `swe_azalt_rev` — azimuth/true-altitude → ecliptic/equatorial (swecl.c:2839–2873)

```c
void CALL_CONV swe_azalt_rev(
      double tjd_ut,
      int32  calc_flag,     /* SE_HOR2ECL or SE_HOR2EQU */
      double *geopos,       /* [0]=geolon, [1]=geolat, [2]=height (unused here) */
      double *xin,          /* [0]=azimuth (from south, clockwise), [1]=TRUE altitude, degrees */
      double *xout)         /* [0]=lon/RA, [1]=lat/dec, degrees */
```

Important: this function inverts **only** the geometric azimuth/altitude → equatorial/ecliptic
transform. It does **not** undo refraction — the input `xin[1]` must already be a **true**
altitude (as documented in the header comment at swecl.c:2830–2838). To go from an apparent
altitude, the caller must first call `swe_refrac(..., SE_APP_TO_TRUE, ...)` (or
`swe_refrac_extended`) themselves.

### Algorithm

1. `armc = swe_degnorm(swe_sidtime(tjd_ut) * 15 + geolon)` — identical construction to
   `swe_azalt` step 1.
2. Copy `xaz[0]=xin[0]` (azimuth), `xaz[1]=xin[1]` (true altitude), `xaz[2]=1`.
3. **Undo the azimuth convention** (inverse of `swe_azalt` step 6): azimuth is from south,
   clockwise; convert to "from east, counterclockwise":
   ```c
   xaz[0] = 360 - xaz[0];
   xaz[0] = swe_degnorm(xaz[0] - 90);
   ```
4. **Inverse horizontal→equatorial rotation**: `dang = geolat - 90` (note: this is the **negative**
   of the forward rotation angle used in `swe_azalt` step 5, which was `90 - geolat` — consistent
   with `swe_cotrans` being applied here as the inverse rotation).
   ```c
   swe_cotrans(xaz, xaz, dang);
   xaz[0] = swe_degnorm(xaz[0] + armc + 90);
   xout[0] = xaz[0];   /* RA */
   xout[1] = xaz[1];   /* dec */
   ```
   At this point `xout` holds **equatorial** (RA, dec) coordinates — this is the final result
   when `calc_flag == SE_HOR2EQU`.
5. **If `calc_flag == SE_HOR2ECL`**: additionally convert equatorial → ecliptic:
   ```c
   swe_calc(tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL), SE_ECL_NUT, 0, x, NULL);
   eps_true = x[0];
   swe_cotrans(xaz, x, eps_true);   /* note: writes into x, reads from xaz; +eps_true (equ→ecl direction) */
   xout[0] = x[0];
   xout[1] = x[1];
   ```
   Same ΔT/`SE_ECL_NUT` pattern as `swe_azalt`. Sign is `+eps_true` here (equatorial→ecliptic),
   opposite of `swe_azalt`'s `-eps_true` (ecliptic→equatorial) — consistent with the `swe_cotrans`
   sign convention documented in swephlib.c:220–222.
   FP-fidelity note: `swe_cotrans(xaz, x, eps_true)` — first argument (input) is `xaz`, second
   (output) is `x`; the call signature is `(xpo, xpn, eps)`, so this reads from `xaz` and writes
   to `x`, not in-place. Preserve this input/output separation; it matters because `xaz[2]` (the
   dummy "1" component) is *not* copied into `x[2]` by this rotation — `swe_cotrans` sets
   `xpn[2] = xpo[2]` internally, so `x[2]` ends up as `xaz[2]==1`, harmless since only `x[0]`/`x[1]`
   are read afterward.

---

## 3. `swe_refrac` — simple sea-level refraction (swecl.c:2887–2984)

```c
double CALL_CONV swe_refrac(double inalt, double atpress, double attemp, int32 calc_flag)
```
`inalt`: input altitude, degrees. `atpress`: hPa. `attemp`: deg C.
`calc_flag`: `SE_TRUE_TO_APP` (0) or `SE_APP_TO_TRUE` (1).

Returns the transformed altitude, degrees. This function assumes a **sea-level observer with an
ideal horizon** — no dip, no elevated-observer geometry (that's what `swe_refrac_extended`
adds). Header comment: "These formulae do not handle the case when the sun is visible below the
geometrical horizon (from a mountain top or an air plane)".

There is a `#if 0`-disabled alternate implementation (S. L. Moshier's, lines 2892–2940, using the
"Almanac for Computers" formula with Newton iteration) that is **dead code**, not compiled — do
not port it. Only the live `#else` branch (Meeus-based, "Meeus, German, p. 114ff.") is active.

### Pressure/temperature scaling factor (computed once, both directions)
```c
double pt_factor = atpress / 1010.0 * 283.0 / (273.0 + attemp);
```
FP-fidelity: exact grouping `(atpress / 1010.0) * 283.0 / (273.0 + attemp)` — left-to-right
division/multiplication as written; do not reorder.

### TRUE → APPARENT (`calc_flag == SE_TRUE_TO_APP`), swecl.c:2944–2965

```c
trualt = inalt;
if (trualt > 15) {
    a = tan((90 - trualt) * DEGTORAD);
    refr = (58.276 * a - 0.0824 * a * a * a);
    refr *= pt_factor / 3600.0;
} else if (trualt > -5) {
    a = trualt + 10.3 / (trualt + 5.11);
    if (a + 1e-10 >= 90)
        refr = 0;
    else
        refr = 1.02 / tan(a * DEGTORAD);
    refr *= pt_factor / 60.0;
} else {
    refr = 0;
}
appalt = trualt;
if (appalt + refr > 0)
    appalt += refr;
return appalt;
```
- **Branch 1** (`trualt > 15`°): coefficients `58.276` and `0.0824`, arcseconds, converted to
  degrees by `/ 3600.0`; cubic-in-`a` refraction term (Meeus's high-altitude polynomial), where
  `a = tan(zenith distance)`.
- **Branch 2** (`-5 < trualt <= 15`): `a = trualt + 10.3/(trualt + 5.11)` (near-singular for
  `trualt` near `-5.11`, guarded only by the domain `trualt > -5`); coefficient `1.02`, arcminutes,
  converted to degrees by `/ 60.0`. Degenerate guard: if `a + 1e-10 >= 90`, `tan(a)` would blow up
  → force `refr = 0` instead.
- **Branch 3** (`trualt <= -5`): `refr = 0` unconditionally (no correction below -5°).
- **Final clamp**: apparent altitude is only nudged by `refr` if doing so keeps it `> 0` — i.e. the
  correction is suppressed (not applied) if `appalt + refr <= 0`. `appalt` (the return value)
  stays at `trualt` unchanged in that case (not clamped to 0).

### APPARENT → TRUE (`calc_flag == SE_APP_TO_TRUE`), swecl.c:2966–2983

```c
appalt = inalt;
a = appalt + 7.31 / (appalt + 4.4);
if (a + 1e-10 >= 90) {
    refr = 0;
} else {
    refr = 1.00 / tan(a * DEGTORAD);
    refr -= 0.06 * sin(14.7 * refr + 13);
}
refr *= pt_factor / 60.0;
trualt = appalt;
if (appalt - refr > 0)
    trualt = appalt - refr;
return trualt;
```
Single unbranched-by-altitude formula (Bennett's formula, refined with a small correction term):
`a = appalt + 7.31/(appalt + 4.4)` (singular near `appalt = -4.4`); base refraction
`refr = cot(a) = 1.00/tan(a·DEGTORAD)`; correction term subtracted: `refr -= 0.06 * sin(14.7*refr + 13)`
(coefficients `0.06`, `14.7`, `13` — the `13` is a phase offset in the same units as `14.7*refr`,
i.e. dimensionless argument to `sin`, degrees since `sin`'s argument here is NOT converted via
DEGTORAD — verify: `sin(14.7 * refr + 13)` — `refr` at this point is in **arcminutes-equivalent
cot units**, and the C `sin()` takes radians; this line does **not** apply `DEGTORAD` — replicate
literally, this is the same non-obviously-unit-converted expression as the original C, a known
quirk of Bennett's empirical formula (the argument to `sin` is an empirical fit variable, not a
geometric angle in a chosen unit — do not "fix" the missing DEGTORAD).
Then `refr *= pt_factor / 60.0` (arcminutes → degrees, same as branch 2 above).
**Final clamp**: symmetric to the true→apparent case — `trualt = appalt - refr` only if that stays
`> 0`; otherwise `trualt` is left at `appalt` (no clamp to 0, just "don't apply the correction if
it would go non-positive").

### Degenerate-domain notes (from the C header comment, swecl.c:2875–2886)
- True→apparent, branch-2 singularity near `trualt ≈ -5.00158` and `≈ 89.89158` (where `a`
  approaches 90° and `tan` blows up) — guarded by the `a + 1e-10 >= 90` check.
- Apparent→true singularity near `inalt ≈ -4.3285` and `≈ 89.9225` — same guard pattern.
- Neither direction handles negative apparent altitude / below-geometric-horizon observers
  (elevated observer) — that's the whole point of `swe_refrac_extended`.

---

## 4. `swe_set_lapse_rate` (swecl.c:2986–2989)

```c
void CALL_CONV swe_set_lapse_rate(double lapse_rate) {
  const_lapse_rate = lapse_rate;
}
```
Sets the file-static global `const_lapse_rate` (declared swecl.c:74:
`static TLS double const_lapse_rate = SE_LAPSE_RATE;`, default `SE_LAPSE_RATE = 0.0065` K/m,
sweph.h:306). `TLS` = thread-local storage qualifier (still process-global-per-thread, i.e.
mutable shared state within a thread's calls).

This global is read in exactly one place in the public API: `swe_azalt`'s call to
`swe_refrac_extended(..., const_lapse_rate, ...)` (swecl.c:2823). It has no effect on
`swe_refrac` (which takes no lapse rate at all) or on direct callers of `swe_refrac_extended`
(which already takes `lapse_rate` as an explicit argument).

STATELESS PORT NOTE: the Rust port has no global to set. `swe_refrac_extended`'s Rust signature
already takes `lapse_rate` as an explicit parameter, so it needs no change. The only port
implication is for a `swe_azalt`-equivalent Rust function: it must take `lapse_rate` as an
explicit parameter (with a `0.0065` default via `EphemerisConfig`/function default, not a
settable global) rather than reading a mutable static. There is no `swe_set_lapse_rate`
equivalent to port as a stateful setter — if a default is wanted, thread it through
`EphemerisConfig` per the project's stateless-config architecture.

---

## 5. `swe_refrac_extended` — elevated-observer refraction with dip (swecl.c:3035–3115)

```c
double CALL_CONV swe_refrac_extended(
      double inalt,       /* altitude of object above geometric horizon, degrees
                              (geometric horizon = plane perpendicular to gravity) */
      double geoalt,       /* observer's altitude above sea level, meters */
      double atpress,      /* hPa; 0 => use standard-atmosphere estimate from geoalt */
      double attemp,       /* deg C */
      double lapse_rate,   /* dT/dh, deg K / m */
      int32  calc_flag,    /* SE_TRUE_TO_APP or SE_APP_TO_TRUE */
      double *dret)        /* out, may be NULL: dret[0]=true alt, [1]=apparent alt,
                                                 [2]=refraction, [3]=dip of horizon */
```
Credit: developed with archaeoastronomer Victor Reijs. Improves on `swe_refrac` by supporting
observer altitude `> 0` (negative apparent heights become meaningful — visible-below-ideal-horizon
case) and by exposing the refraction constant (`lapse_rate`) as a parameter.

Return contract (also documented in the header comment, swecl.c:3013–3033):
- `dret[0]` = true altitude if determinable, else the input value unchanged.
- `dret[1]` = apparent altitude if determinable, else the input value unchanged.
- `dret[2]` = refraction amount, degrees.
- `dret[3]` = dip of the horizon, degrees (always filled, both directions).
- **"The body is above the horizon if `dret[0] != dret[1]`."**

### Setup (always run, regardless of `calc_flag`)
```c
double dip = calc_dip(geoalt, atpress, attemp, lapse_rate);   /* §7 below */
if (inalt > 90)
    inalt = 180 - inalt;     /* fold altitudes >90° back onto [.., 90] domain */
```
Note: `atpress` here is passed through to `calc_dip` **as-is** — if the caller passed `atpress==0`,
`calc_dip` receives `0` too (unlike `swe_azalt`, which pre-resolves the pressure estimate itself
before calling `swe_refrac_extended`; direct callers of `swe_refrac_extended` with `atpress==0`
must be aware `calc_dip`'s internal formula uses `atpress` directly, i.e. `0`, — see §7, there is
no auto-estimate inside `calc_dip` or `swe_refrac_extended` itself; the "pressure-from-altitude
default when atpress==0" behavior belongs to `swe_azalt`'s caller-side logic only, `swe_refrac_extended`
does not itself default a zero pressure to the standard-atmosphere estimate).

### TRUE → APPARENT (`calc_flag == SE_TRUE_TO_APP`), swecl.c:3045–3088

1. **Early exit for very negative altitude**: if `inalt < -10`, refraction is not computed at all:
   ```c
   dret[0]=inalt; dret[1]=inalt; dret[2]=0; dret[3]=dip;
   return inalt;
   ```
2. **Newton iteration to invert `calc_astronomical_refr`** (which is itself parametrized by
   *apparent* altitude — see §6): 5 fixed iterations, no early-exit convergence check:
   ```c
   y = inalt; D = 0.0; yy0 = 0; D0 = D;
   for (i = 0; i < 5; i++) {
     D = calc_astronomical_refr(y, atpress, attemp);
     N = y - yy0;
     yy0 = D - D0 - N;                 /* denominator of derivative */
     if (N != 0.0 && yy0 != 0.0)
       N = y - N * (inalt + D - y) / yy0;   /* Newton step, numerically estimated derivative */
     else
       N = inalt + D;                  /* first pass fallback */
     yy0 = y; D0 = D;
     y = N;
   }
   refr = D;   /* NOTE: refr is D from the LAST loop iteration's calc_astronomical_refr call,
                  evaluated at y from the PREVIOUS iteration — not recomputed at the final y */
   ```
   FP-fidelity: the comment "sic !!! code by Moshier" at swecl.c:3064 flags that this
   numerically-estimated-derivative Newton iteration (reusing consecutive function evaluations as
   a secant-like derivative, rather than an analytic derivative) is deliberate, inherited
   Moshier code — replicate the exact 5-iteration loop and variable reuse pattern, not a
   "cleaned up" fixed-point or bisection substitute. Exactly 5 iterations, unconditional (no
   convergence check/early break).
3. **Below-dip check**: `if (inalt + refr < dip)` → treat as not visible / no refraction:
   ```c
   dret[0]=inalt; dret[1]=inalt; dret[2]=0; dret[3]=dip;
   return inalt;
   ```
4. **Otherwise**, refraction applies:
   ```c
   dret[0] = inalt;
   dret[1] = inalt + refr;
   dret[2] = refr;
   dret[3] = dip;
   return inalt + refr;
   ```

### APPARENT → TRUE (`calc_flag == SE_APP_TO_TRUE`), swecl.c:3089–3114

```c
refr = calc_astronomical_refr(inalt, atpress, attemp);   /* single evaluation, no iteration —
                                                              calc_astronomical_refr's input IS
                                                              apparent altitude already */
trualt = inalt - refr;
if (inalt > dip) {
    dret[0] = trualt;
    dret[1] = inalt;
    dret[2] = refr;
    dret[3] = dip;
} else {
    dret[0] = inalt;
    dret[1] = inalt;
    dret[2] = 0;
    dret[3] = dip;
}
if (inalt >= dip)          /* comment: "bug fix dieter, 4 feb 20" (was `trualt > dip` before) */
    return trualt;
else
    return inalt;
```
FP-fidelity / behavioral note: the `dret[]`-fill condition (`inalt > dip`, strict) and the
**return-value** condition (`inalt >= dip`, inclusive) are subtly different comparisons (`>` vs
`>=`) — this asymmetry is intentional per the inline comment documenting a historical bug fix
(previously the return-value check used `trualt > dip` instead of `inalt >= dip`). Replicate both
comparisons exactly as written, including the boundary-inclusivity difference between the two.

---

## 6. `calc_astronomical_refr` (static, swecl.c:3124–3148)

```c
static double calc_astronomical_refr(double inalt, double atpress, double attemp)
```
`inalt` here is an **apparent** altitude (see the function's own header comment, swecl.c:3117–3123).
Computes refraction in degrees.

There is a `#if 0`-disabled alternate (Bennett 1982 formula, "Journal of Inst. Navigation No. 35,
p.255-259", formula H) — **dead code**, not compiled, do not port.

### Live (`#else`) branch — Sinclair's formula, swecl.c:3136–3147
```c
double r;
if (inalt > 17.904104638432) {          /* chosen so the two branches are C0-continuous at this altitude */
    r = 0.97 / tan(inalt * DEGTORAD);
} else {
    r = (34.46 + 4.23 * inalt + 0.004 * inalt * inalt)
        / (1 + 0.505 * inalt + 0.0845 * inalt * inalt);
}
r = ((atpress - 80) / 930 / (1 + 0.00008 * (r + 39) * (attemp - 10)) * r) / 60.0;
return r;
```
Constants:
- Branch threshold: `17.904104638432`° (comment: "for continuous function, instead of '>15'" —
  i.e. this exact value, not `15`, is where the two sub-formulas cross continuously; do not
  round or approximate it).
- High-altitude branch: coefficient `0.97` (cotangent formula).
- Low-altitude branch: rational polynomial in `inalt` with numerator coefficients
  `34.46, 4.23, 0.004` and denominator coefficients `1, 0.505, 0.0845`.
- Final pressure/temperature scaling (applied identically after either branch):
  `((atpress - 80) / 930 / (1 + 0.00008 * (r + 39) * (attemp - 10)) * r) / 60.0` — converts `r`
  from arcminutes to degrees (`/60.0`) while simultaneously scaling for pressure (baseline `80`
  hPa subtracted, `930` hPa span) and temperature (`0.00008` coefficient, `39` offset added to `r`
  before multiplying by `(attemp - 10)`).
  FP-fidelity: preserve the exact grouping — `(atpress - 80) / 930` computed first, divided by
  `(1 + 0.00008*(r+39)*(attemp-10))`, that whole quotient multiplied by `r`, then the entire
  product divided by `60.0`. Do not reassociate the two divisions.

---

## 7. `calc_dip` (static, swecl.c:3158–3169)

```c
static double calc_dip(double geoalt, double atpress, double attemp, double lapse_rate)
```
Computes the geometric + refractive dip of the horizon (in degrees, returned as a **negative**
number — the horizon appears below the astronomical horizontal plane for `geoalt > 0`) for an
observer at height `geoalt` meters above sea level.

Source: A. Thom, *Megalithic Lunar Observations*, 1973, p.32; metric conversion by V. Reijs, 2000
(archaeocosmology.org/eng/refract.htm#Sea).

```c
double krefr = (0.0342 + lapse_rate) / (0.154 * 0.0238);
double d = 1 - 1.8480 * krefr * atpress / (273.15 + attemp) / (273.15 + attemp);
return -180.0 / PI * acos(1 / (1 + geoalt / EARTH_RADIUS)) * sqrt(d);
```
Constants:
- `krefr` numerator constant `0.0342` added to the caller-supplied `lapse_rate`; denominator
  `0.154 * 0.0238` (product of two literals — compute in that order, or as their product
  `0.0036652`, matching whichever the compiler would; for FP fidelity keep as written:
  `0.154 * 0.0238`).
- `d`: coefficient `1.8480`, and note the **double division** `atpress / (273.15+attemp) / (273.15+attemp)`
  — i.e. `atpress / (273.15+attemp)²` but expressed as two sequential divisions by the same
  denominator, not `atpress / pow(273.15+attemp, 2)`. Replicate the two-division form for FP
  fidelity (order of floating point operations differs subtly from squaring the denominator
  first).
- `273.15` — Celsius→Kelvin conversion, used twice (matches `attemp` in deg C).
- Two commented-out alternate return expressions (lines 3166–3167,
  `-0.03203*sqrt(geoalt)*sqrt(d)` and an `EARTH_RADIUS`-based `acos` without the final
  `sqrt(d)` factor) are dead code — do not port.
- Live return: `-180.0/PI * acos(1/(1 + geoalt/EARTH_RADIUS)) * sqrt(d)` — `PI` is the library's
  standard `M_PI`-equivalent constant (defined elsewhere in the codebase, e.g. `sweph.h`/`swephlib.h`);
  `EARTH_RADIUS = 6378136.6` meters (sweph.h:282, "AA 2006 K6"). The leading minus sign makes the
  dip negative (angle below the horizontal), consistent with "dip of the horizon" convention where
  a horizon-below-observer angle is signed negative.
- **`atpress == 0` behavior**: `calc_dip` does **not** special-case `atpress == 0` — if the caller
  passes `atpress = 0` (e.g. a direct `swe_refrac_extended` call without pre-resolving pressure),
  `d = 1 - 0 = 1`, so `sqrt(d) = 1`, i.e. the refractive contribution to `krefr`'s pressure term
  vanishes and dip reduces to the purely geometric term scaled by `1`. This is different from
  `swe_azalt`, which resolves `atpress` to a standard-atmosphere estimate *before* calling
  `swe_refrac_extended` (see §1 step 7) — so end-to-end, `swe_azalt`'s dip calculation for
  `atpress==0` input does get the full pressure-refraction term, but a bare `swe_refrac_extended`
  call with `atpress==0` does not. Rust port: mirror this exactly — do not silently insert an
  auto-pressure-estimate inside a `calc_dip`/`refrac_extended` port; keep the estimate logic at
  the `azalt`-equivalent call site only, per the C structure.

---

## 8. Cross-references / shared helpers used

- `swe_sidtime(tjd_ut)` — Greenwich apparent sidereal time, **hours** (swephlib.c:3580). Rust
  port should already have this from the houses module (`docs/c-ref-houses.md` references the
  same `armc = swe_degnorm(swe_sidtime(tjd_ut)*15 + geolon)` pattern) — reuse, don't duplicate.
- `swe_deltat_ex(tjd_ut, iflag, serr)` with `iflag = -1` — ΔT with automatic ephemeris-flag
  resolution; used identically in `swe_azalt` and `swe_azalt_rev` to convert `tjd_ut` → `tjd` (TT)
  before the `SE_ECL_NUT` calc.
- `swe_calc(tjd, SE_ECL_NUT, 0, x, NULL)` — returns true obliquity in `x[0]`. Both `swe_azalt`
  (ecl→equ direction) and `swe_azalt_rev` (equ→ecl direction) use this same call.
- `swe_cotrans(xpo, xpn, eps)` (swephlib.c:223–240) — generic polar-coordinate rotation about the
  x-axis by `eps` degrees; ecliptic→equatorial uses `eps` negative, equatorial→ecliptic uses `eps`
  positive (its own doc comment). Also reused generically inside `swe_azalt`/`swe_azalt_rev` for
  the equatorial↔horizontal rotation (by `90 - geolat` / `geolat - 90` respectively) — same
  routine, different physical meaning of the rotation axis in that call site.
- `EARTH_RADIUS` (sweph.h:282) and `SE_LAPSE_RATE` (sweph.h:306) — see Function Map header.

---

## 9. STATELESS PORT NOTES summary

1. `const_lapse_rate` (swecl.c:74, a `static TLS` global defaulting to `SE_LAPSE_RATE = 0.0065`)
   is set by `swe_set_lapse_rate` and read only inside `swe_azalt`. The Rust `azalt`-equivalent
   function must take `lapse_rate: f64` as an explicit parameter (default `0.0065` via config or
   function default), not a mutable global/setter pair.
2. `swe_refrac_extended` in C already takes `lapse_rate` as an explicit argument — no statefulness
   to remove there; port signature 1:1.
3. `calc_dip`'s pressure handling: it does **not** auto-estimate `atpress` when `0` — only
   `swe_azalt`'s call site does that (§1 step 7, §7 last bullet). Keep the auto-estimate logic at
   the same call-site boundary in the Rust port, not inside a `calc_dip`/`refrac_extended` helper.
4. No other global state is touched by this function group (unlike, e.g., `swe_set_topo`'s
   observer-position globals used elsewhere in `swecl.c`) — `swe_azalt`/`swe_azalt_rev` receive
   `geopos` explicitly every call.
