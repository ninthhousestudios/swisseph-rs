# C Reference: Ayanamsa / Sidereal Mode — sweph.c + sweph.h

Porting reference for all ayanamsa (sidereal mode) computation. Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_set_sid_mode` | sweph.c:2861–2928 | Yes — public setter |
| `swi_get_ayanamsa_ex` | sweph.c:3002–3208 | Yes — core internal (no nutation) |
| `swe_get_ayanamsa_ex` | sweph.c:2930–2946 | Yes — public wrapper (adds nutation) |
| `swe_get_ayanamsa_ex_ut` | sweph.c:3226–3244 | Yes — UT-input variant |
| `swe_get_ayanamsa` | sweph.c:3251–3258 | Yes — legacy double return, no nutation |
| `swe_get_ayanamsa_ut` | sweph.c:3260–3266 | Yes — legacy double return, UT input |
| `swi_get_ayanamsa_with_speed` | sweph.c:3210–3224 | Yes — two-point numerical speed |
| `get_aya_correction` | sweph.c:2960–3000 | Yes — static, precession model correction |
| `swi_trop_ra2sid_lon` | sweph.c:3273–3301 | Yes — ECL_T0 projection |
| `swi_trop_ra2sid_lon_sosy` | sweph.c:3308–3356 | Yes — SSY_PLANE projection |

---

## SE_SIDBIT Flag Constants (swephexp.h:221–235)

```c
#define SE_SIDBITS             256   /* mask for the table index (low 8 bits) */
#define SE_SIDBIT_ECL_T0       256   /* project onto ecliptic of t0 */
#define SE_SIDBIT_SSY_PLANE    512   /* project onto solar system equatorial plane */
#define SE_SIDBIT_USER_UT      1024  /* for SE_SIDM_USER: t0 is UT, not TT */
#define SE_SIDBIT_ECL_DATE     2048  /* ayanamsa measured on ecliptic of date */
#define SE_SIDBIT_NO_PREC_OFFSET  4096  /* suppress get_aya_correction() */
#define SE_SIDBIT_PREC_ORIG    8192  /* use ayanamsa's original precession model */
```

**How `sid_mode` is decomposed**: the full `sid_mode` stored in `swed.sidd.sid_mode` carries both
the index (bits 0–7, mask `SE_SIDBITS = 256`) and the projection flags (bits 8+). To get the bare
table index: `sid_mode % SE_SIDBITS` (= `sid_mode & 0xFF`).

---

## SE_SIDM_* Index Constants (swephexp.h:238–286)

```c
#define SE_SIDM_FAGAN_BRADLEY    0
#define SE_SIDM_LAHIRI           1
#define SE_SIDM_DELUCE           2
#define SE_SIDM_RAMAN            3
#define SE_SIDM_USHASHASHI       4
#define SE_SIDM_KRISHNAMURTI     5
#define SE_SIDM_DJWHAL_KHUL      6
#define SE_SIDM_YUKTESHWAR       7
#define SE_SIDM_JN_BHASIN        8
#define SE_SIDM_BABYL_KUGLER1    9
#define SE_SIDM_BABYL_KUGLER2   10
#define SE_SIDM_BABYL_KUGLER3   11
#define SE_SIDM_BABYL_HUBER     12
#define SE_SIDM_BABYL_ETPSC     13
#define SE_SIDM_ALDEBARAN_15TAU 14
#define SE_SIDM_HIPPARCHOS      15
#define SE_SIDM_SASSANIAN       16
#define SE_SIDM_GALCENT_0SAG    17
#define SE_SIDM_J2000           18
#define SE_SIDM_J1900           19
#define SE_SIDM_B1950           20
#define SE_SIDM_SURYASIDDHANTA  21
#define SE_SIDM_SURYASIDDHANTA_MSUN 22
#define SE_SIDM_ARYABHATA       23
#define SE_SIDM_ARYABHATA_MSUN  24
#define SE_SIDM_SS_REVATI       25
#define SE_SIDM_SS_CITRA        26
#define SE_SIDM_TRUE_CITRA      27
#define SE_SIDM_TRUE_REVATI     28
#define SE_SIDM_TRUE_PUSHYA     29
#define SE_SIDM_GALCENT_RGILBRAND 30
#define SE_SIDM_GALEQU_IAU1958  31
#define SE_SIDM_GALEQU_TRUE     32
#define SE_SIDM_GALEQU_MULA     33
#define SE_SIDM_GALALIGN_MARDYKS 34
#define SE_SIDM_TRUE_MULA       35
#define SE_SIDM_GALCENT_MULA_WILHELM 36
#define SE_SIDM_ARYABHATA_522   37
#define SE_SIDM_BABYL_BRITTON   38
#define SE_SIDM_TRUE_SHEORAN    39
#define SE_SIDM_GALCENT_COCHRANE 40
#define SE_SIDM_GALEQU_FIORENZA 41
#define SE_SIDM_VALENS_MOON     42
#define SE_SIDM_LAHIRI_1940     43
#define SE_SIDM_LAHIRI_VP285    44
#define SE_SIDM_KRISHNAMURTI_VP291 45
#define SE_SIDM_LAHIRI_ICRC     46

#define SE_SIDM_USER           255   /* user-defined; t0 is TT unless USER_UT */
#define SE_NSIDM_PREDEF         47   /* count of predefined entries */
```

---

## `struct aya_init` Definition (sweph.h:347–350)

```c
struct aya_init {
    double  t0;          /* epoch (JD), TT or UT depending on t0_is_UT */
    double  ayan_t0;     /* ayanamsa value in DEGREES at t0 */
    AS_BOOL t0_is_UT;    /* TRUE = t0 is UT; FALSE = t0 is TT/ET */
    int     prec_offset; /* precession model the ayanamsa was defined with:
                          *   0  = no correction needed
                          *  -1  = unclear / not applied
                          *   1  = SEMOD_PREC_IAU_1976
                          *  11  = SEMOD_PREC_NEWCOMB */
};
static const struct aya_init ayanamsa[SE_NSIDM_PREDEF];  // sweph.h:351
```

---

## Ayanamsa Data Table (sweph.h:351–596)

All 47 entries, indices 0–46. For indexes where `prec_offset != 0`, the value in
parentheses is the symbolic constant.

Computed numeric values for `ayan_t0` expressions from source:
- Index 1: `23.250182778 - 0.004658035 = 23.245524743`
- Index 3: `360 - 338.98556 = 21.014440`
- Index 4: `360 - 341.33904 = 18.660960`
- Index 5: `360 - 337.636111 = 22.363889`
- Index 6: `360 - 333.0369024 = 26.963098`
- Index 7: `360 - 338.917778 = 21.082222`
- Index 8: `360 - 338.634444 = 21.365556`
- Index 46: `23.25 - 0.00464207 = 23.245358`

Reference epochs as JD:
- `J2000 = 2451545.0` (sweph.h:67)
- `J1900 = 2415020.0` (sweph.h:69)
- `B1950 = 2433282.42345905` (sweph.h:68)

| Idx | Constant | `t0` (JD) | `ayan_t0` (°) | `t0_is_UT` | `prec_offset` | Notes |
|-----|----------|-----------|----------------|------------|---------------|-------|
| 0 | FAGAN_BRADLEY | 2433282.42346 | 24.042044444 | FALSE | 11 (NEWCOMB) | Default |
| 1 | LAHIRI | 2435553.5 | 23.245524743 | FALSE | 1 (IAU_1976) | |
| 2 | DELUCE | 1721057.5 | 0.0 | TRUE | 0 | |
| 3 | RAMAN | 2415020.0 (J1900) | 21.014440 | FALSE | 11 (NEWCOMB) | |
| 4 | USHASHASHI | 2415020.0 (J1900) | 18.660960 | FALSE | -1 | |
| 5 | KRISHNAMURTI | 2415020.0 (J1900) | 22.363889 | FALSE | 11 (NEWCOMB) | |
| 6 | DJWHAL_KHUL | 2415020.0 (J1900) | 26.963098 | FALSE | 0 | |
| 7 | YUKTESHWAR | 2415020.0 (J1900) | 21.082222 | FALSE | -1 | |
| 8 | JN_BHASIN | 2415020.0 (J1900) | 21.365556 | FALSE | -1 | |
| 9 | BABYL_KUGLER1 | 1684532.5 | -5.66667 | TRUE | -1 | |
| 10 | BABYL_KUGLER2 | 1684532.5 | -4.26667 | TRUE | -1 | |
| 11 | BABYL_KUGLER3 | 1684532.5 | -3.41667 | TRUE | -1 | |
| 12 | BABYL_HUBER | 1684532.5 | -4.46667 | TRUE | -1 | |
| 13 | BABYL_ETPSC | 1673941.0 | -5.079167 | TRUE | -1 | |
| 14 | ALDEBARAN_15TAU | 1684532.5 | -4.44138598 | TRUE | 0 | |
| 15 | HIPPARCHOS | 1674484.0 | -9.33333 | TRUE | -1 | |
| 16 | SASSANIAN | 1927135.8747793 | 0.0 | TRUE | -1 | |
| 17 | GALCENT_0SAG | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 18 | J2000 | 2451545.0 (J2000) | 0.0 | FALSE | 0 | ECL_T0 auto-set |
| 19 | J1900 | 2415020.0 (J1900) | 0.0 | FALSE | 0 | ECL_T0 auto-set |
| 20 | B1950 | 2433282.42345905 (B1950) | 0.0 | FALSE | 0 | ECL_T0 auto-set |
| 21 | SURYASIDDHANTA | 1903396.8128654 | 0.0 | TRUE | 0 | |
| 22 | SURYASIDDHANTA_MSUN | 1903396.8128654 | -0.21463395 | TRUE | 0 | |
| 23 | ARYABHATA | 1903396.7895321 | 0.0 | TRUE | 0 | |
| 24 | ARYABHATA_MSUN | 1903396.7895321 | -0.23763238 | TRUE | 0 | |
| 25 | SS_REVATI | 1903396.8128654 | -0.79167046 | TRUE | 0 | |
| 26 | SS_CITRA | 1903396.8128654 | 2.11070444 | TRUE | 0 | |
| 27 | TRUE_CITRA | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 28 | TRUE_REVATI | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 29 | TRUE_PUSHYA | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 30 | GALCENT_RGILBRAND | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 31 | GALEQU_IAU1958 | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 32 | GALEQU_TRUE | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 33 | GALEQU_MULA | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 34 | GALALIGN_MARDYKS | 2451079.734892000 | 30.0 | FALSE | 0 | ECL_T0 auto-set |
| 35 | TRUE_MULA | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 36 | GALCENT_MULA_WILHELM | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 37 | ARYABHATA_522 | 1911797.740782065 | 0.0 | TRUE | 0 | |
| 38 | BABYL_BRITTON | 1721057.5 | -3.2 | TRUE | -1 | |
| 39 | TRUE_SHEORAN | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 40 | GALCENT_COCHRANE | 0.0 | 0.0 | FALSE | 0 | **Fixed-star** override |
| 41 | GALEQU_FIORENZA | 2451544.5 | 25.0 | TRUE | 0 | |
| 42 | VALENS_MOON | 1775845.5 | -2.9422 | TRUE | -1 | |
| 43 | LAHIRI_1940 | 2415020.0 (J1900) | 22.44597222 | FALSE | 11 (NEWCOMB) | |
| 44 | LAHIRI_VP285 | 1825235.2458513028 | 0.0 | FALSE | 0 | |
| 45 | KRISHNAMURTI_VP291 | 1827424.752255678 | 0.0 | FALSE | 0 | |
| 46 | LAHIRI_ICRC | 2435553.5 | 23.245358 | FALSE | 11 (NEWCOMB) | |

**Note on source expressions**: The `ayan_t0` for entries 1, 3–8, 46 are stored in source as
arithmetic expressions (e.g. `360 - 338.98556`). The Rust port must evaluate these to double
precision and store the result. The computed values in the table above are reference — verify by
reading sweph.h:372–594 directly when transcribing to Rust.

**Entries with `t0 == 0.0`**: Indices 17, 27–33, 35, 36, 39, 40. These are all fixed-star
override ayanamsas. The table data is never consulted; the early-return block in
`swi_get_ayanamsa_ex` handles them before any table lookup.

---

## `swe_set_sid_mode()` — Mode Setter (sweph.c:2861–2928)

```c
void swe_set_sid_mode(int32 sid_mode, double t0, double ayan_t0)
```

### Step-by-step

```
sip = &swed.sidd          // struct sid_data pointer
swi_init_swed_if_start()
if sid_mode < 0: sid_mode = 0
sip->sid_mode = sid_mode

// Decompose index (low 8 bits)
idx = sid_mode % SE_SIDBITS   // = sid_mode & 0xFF
```

**Standard equinox auto-flag** (sweph.c:2871–2878):
If `idx` is one of {J2000=18, J1900=19, B1950=20, GALALIGN_MARDYKS=34}:
```
sip->sid_mode = idx | SE_SIDBIT_ECL_T0
```
(Any previously set bits in sid_mode are REPLACED, not ORed in.)

**True-star / galactic clears all projection bits** (sweph.c:2880–2895):
If `idx` is one of {TRUE_CITRA=27, TRUE_REVATI=28, TRUE_PUSHYA=29, TRUE_SHEORAN=39,
TRUE_MULA=35, GALCENT_0SAG=17, GALCENT_COCHRANE=40, GALCENT_RGILBRAND=30,
GALCENT_MULA_WILHELM=36, GALEQU_IAU1958=31, GALEQU_TRUE=32, GALEQU_MULA=33}:
```
sip->sid_mode = idx   // strips ECL_T0, SSY_PLANE, USER_UT bits
```

**Bounds check** (sweph.c:2897–2898):
```
if idx >= SE_NSIDM_PREDEF && idx != SE_SIDM_USER:
    sip->sid_mode = sid_mode = SE_SIDM_FAGAN_BRADLEY  // silent fallback
```

**Mark ayana set** (sweph.c:2899):
```
swed.ayana_is_set = TRUE
```

**Load t0 / ayan_t0** (sweph.c:2900–2909):
```
if sid_mode == SE_SIDM_USER:
    sip->t0      = t0        // from caller parameter
    sip->ayan_t0 = ayan_t0   // from caller parameter
    sip->t0_is_UT = FALSE
    if sip->sid_mode & SE_SIDBIT_USER_UT:
        sip->t0_is_UT = TRUE
else:
    sip->t0      = ayanamsa[idx].t0
    sip->ayan_t0 = ayanamsa[idx].ayan_t0
    sip->t0_is_UT = ayanamsa[idx].t0_is_UT
```

**PREC_ORIG feature** (sweph.c:2912–2926): If `SE_SIDBIT_PREC_ORIG` is set AND
`ayanamsa[idx].prec_offset > 0`:
```
swed.astro_models[SE_MODEL_PREC_LONGTERM]  = ayanamsa[idx].prec_offset
swed.astro_models[SE_MODEL_PREC_SHORTTERM] = ayanamsa[idx].prec_offset
// Also update nutation model to match:
if prec_offset == SEMOD_PREC_NEWCOMB (11):
    swed.astro_models[SE_MODEL_NUT] = SEMOD_NUT_WOOLARD (5)
elif prec_offset == SEMOD_PREC_IAU_1976 (1):
    swed.astro_models[SE_MODEL_NUT] = SEMOD_NUT_IAU_1980 (1)
```

**Force recalculation** (sweph.c:2927):
```
swi_force_app_pos_etc()
```

### `struct sid_data` fields written (sweph.h:767+)

```c
struct sid_data {
    int32   sid_mode;    /* full mode word: index | projection bits */
    double  ayan_t0;     /* ayanamsa value in degrees at t0 */
    double  t0;          /* epoch JD (TT or UT per t0_is_UT) */
    AS_BOOL t0_is_UT;    /* whether t0 is UT */
};
```

---

## `swi_get_ayanamsa_ex()` — Core Ayanamsa (sweph.c:3002–3208)

```c
int32 swi_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr)
```

**Always uses `SEFLG_NONUT` internally.** Nutation is never included here; the public
wrapper `swe_get_ayanamsa_ex()` adds it afterward.

### Entry setup (sweph.c:3002–3028)

```
iflag   = plaus_iflag(iflag, -1, tjd_et, serr)
epheflag = iflag & SEFLG_EPHMASK
otherflag = iflag & ~SEFLG_EPHMASK
*daya = 0.0
iflag = (iflag & SEFLG_EPHMASK) | SEFLG_NONUT   // strips everything except ephe + NONUT
sid_mode = sip->sid_mode % SE_SIDBITS             // bare index
iflag_galequ = iflag | SEFLG_TRUEPOS              // used for galactic pole queries
// iflag_true: passes TRUEPOS / NOABERR / NOGDEFL from otherflag
```

If `!swed.ayana_is_set`: call `swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0)` first.

### Early return: fixed-star modes (sweph.c:3049–3142)

These modes bypass all table/precession logic entirely. Each calls `swe_fixstar()` and
returns immediately. See §"Fixed-Star Ayanamsas" below for details.

### Main computation: Method 1 — default (sweph.c:3143–3173)

Activated when `!(sip->sid_mode & SE_SIDBIT_ECL_DATE)` (the common case).

**Conceptual meaning**: Precess the vernal point of `tjd_et` back to `t0` on the ecliptic
of `t0`, then compute the longitude of that point and add `ayan_t0`.

```
// Start: vernal point at tjd_et, in J2000 equatorial Cartesian
x = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

// Step 1: precess date -> J2000
if tjd_et != J2000:
    swi_precess(x, tjd_et, 0, J_TO_J2000)

// Step 2: precess J2000 -> t0
t0 = sip->t0
if sip->t0_is_UT:
    t0 += swe_deltat_ex(t0, iflag, serr)
swi_precess(x, t0, 0, J2000_TO_J)

// Step 3: rotate equatorial t0 -> ecliptic t0
eps = swi_epsiln(t0, 0)   // mean obliquity at t0, RADIANS
swi_coortrf(x, x, eps)    // positive eps = equatorial -> ecliptic

// Step 4: to polar
swi_cartpol(x, x)         // x[0] = ecliptic longitude in RADIANS

// Step 5: ayanamsa = -(longitude of VP at t0) + initial value
x[0] = -x[0] * RADTODEG + sip->ayan_t0   // result in degrees
```

Then fall through to correction and normalization (§ below).

### Main computation: Method 2 — ECL_DATE (sweph.c:3175–3202)

Activated when `sip->sid_mode & SE_SIDBIT_ECL_DATE`. Programmed 2020-05-15.
Added to SE 2.09. Tracks the zero-point through time on the ecliptic of date.

**Conceptual meaning**: Start from the initial ayanamsa point on the ecliptic of `t0`,
propagate it to the ecliptic of `tjd_et`.

```
// x starts as ecliptic polar: [ayan_t0_rad, 0, 1]
x[0] = swe_degnorm(sip->ayan_t0) * DEGTORAD   // ayanamsa point longitude in radians
x[1] = 0.0; x[2] = 1.0

// Step 1: get t0 (possibly convert UT->TT)
t0 = sip->t0
if sip->t0_is_UT:
    t0 += swe_deltat_ex(t0, iflag, serr)

// Step 2: obliquity at t0
eps = swi_epsiln(t0, 0)   // radians

// Step 3: polar ecliptic -> Cartesian ecliptic -> equatorial at t0
swi_polcart(x, x)         // ecliptic Cartesian at t0
swi_coortrf(x, x, -eps)   // negative = ecliptic -> equatorial

// Step 4: precess equatorial t0 -> J2000
if t0 != J2000:
    swi_precess(x, t0, 0, J_TO_J2000)

// Step 5: precess equatorial J2000 -> date
swi_precess(x, tjd_et, 0, J2000_TO_J)

// Step 6: obliquity at date
eps = swi_epsiln(tjd_et, 0)   // radians

// Step 7: equatorial date -> ecliptic date -> polar
swi_coortrf(x, x, eps)   // positive = equatorial -> ecliptic of date
swi_cartpol(x, x)
x[0] = swe_degnorm(x[0] * RADTODEG)   // ayanamsa in degrees
```

### Correction and normalization (sweph.c:3203–3207, both methods)

```
get_aya_correction(iflag, &corr, serr)
*daya = swe_degnorm(x[0] - corr)
return iflag    // (= SEFLG_EPHMASK | SEFLG_NONUT)
```

---

## `swe_get_ayanamsa_ex()` — Public Wrapper (sweph.c:2930–2946)

```c
int32 swe_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr)
```

Adds nutation to the result from `swi_get_ayanamsa_ex`:

```
retval = swi_get_ayanamsa_ex(tjd_et, iflag, daya, serr)
if !(iflag & SEFLG_NONUT):
    if tjd_et == swed.nut.tnut:
        nutp = &swed.nut           // reuse cached nutation
    else:
        swi_nutation(tjd_et, iflag, nuttmp.nutlo)
        nutp = &nuttmp
    *daya += nutp->nutlo[0] * RADTODEG   // add dpsi in degrees
    retval &= ~SEFLG_NONUT               // remove the internally-forced NONUT bit
return retval
```

Note: `nutp->nutlo[0]` is dpsi in **radians**, so `* RADTODEG` converts to degrees.
The subtracted-then-added nutation delta is exactly the nutation in longitude (dpsi).

---

## `get_aya_correction()` — Precession-Model Correction (sweph.c:2960–3000)

Static helper. Computes a small angular correction (in degrees) to compensate for the fact
that some ayanamsas were originally defined using a precession model different from the
current one (Vondrak 2011 by default).

### Returns 0 (no correction) when (sweph.c:2969–2977)

- `sip->t0 == J2000` (precession starts at J2000, no offset)
- `sip->sid_mode & SE_SIDBIT_NO_PREC_OFFSET` (explicitly disabled)
- `prec_offset == 0` (ayanamsa had no specific model)
- `prec_offset < 0` → forced to 0 (model unclear)
- `prec_model == prec_offset` (already using the matching model)

where:
```
prec_model  = swed.astro_models[SE_MODEL_PREC_LONGTERM]   // current model
prec_offset = ayanamsa[sid_mode].prec_offset               // defined model
```

### Correction computation (sweph.c:2978–2998)

```
t0 = sip->t0
if sip->t0_is_UT:
    t0 += swe_deltat_ex(t0, iflag, serr)

// Vernal point in J2000 equatorial (unit vector)
x = [1.0, 0.0, 0.0]

// Precess t0->J2000 with CURRENT model
swi_precess(x, t0, 0, J_TO_J2000)

// Temporarily switch to the ayanamsa's ORIGINAL model
swed.astro_models[SE_MODEL_PREC_LONGTERM]  = prec_offset
swed.astro_models[SE_MODEL_PREC_SHORTTERM] = prec_offset
// Precess J2000->t0 with ORIGINAL model
swi_precess(x, t0, 0, J2000_TO_J)
// Restore current models
swed.astro_models[SE_MODEL_PREC_LONGTERM]  = prec_model
swed.astro_models[SE_MODEL_PREC_SHORTTERM] = prec_model_short

// To ecliptic of t0
eps = swi_epsiln(t0, 0)   // mean obliquity, radians (no iflag)
swi_coortrf(x, x, eps)    // equatorial -> ecliptic

// To polar
swi_cartpol(x, x)         // x[0] = ecliptic longitude in radians

// Correction in degrees
*corr = x[0] * RADTODEG
if *corr > 350.0:
    *corr -= 360.0    // make it a signed value near 0
```

The correction is **subtracted** in the caller: `*daya = swe_degnorm(x[0] - corr)`.

---

## `swi_get_ayanamsa_with_speed()` — Speed Derivative (sweph.c:3210–3224)

```c
int32 swi_get_ayanamsa_with_speed(double tjd_et, int32 iflag, double *daya, char *serr)
```

Two-point numerical differentiation:

```
tintv = 0.001   // days
t2 = tjd_et - tintv

swi_get_ayanamsa_ex(t2,     iflag, &daya_t2, serr)
swi_get_ayanamsa_ex(tjd_et, iflag, &daya[0], serr)

daya[1] = (daya[0] - daya_t2) / tintv   // degrees/day
```

`daya[0]` = ayanamsa in degrees; `daya[1]` = ayanamsa speed in degrees/day.

---

## Fixed-Star Ayanamsas (sweph.c:3049–3142)

These indices are detected **before** any table or precession logic in `swi_get_ayanamsa_ex`.
They call `swe_fixstar()` to get the current position of a catalog star, then subtract a
constant to define the sidereal zero point.

### Indices, stars, and offsets

| Index | Constant | Catalog key | `iflag_true` or `iflag_galequ` | `*daya` computation |
|---|---|---|---|---|
| 17 | GALCENT_0SAG | `",SgrA*"` | `iflag_true` | `degnorm(x[0] - 240.0)` |
| 27 | TRUE_CITRA | `"Spica"` | `iflag_true` | `degnorm(x[0] - 180.0)` |
| 28 | TRUE_REVATI | `",zePsc"` | `iflag_true` | `degnorm(x[0] - 359.8333333333)` |
| 29 | TRUE_PUSHYA | `",deCnc"` | `iflag_true` | `degnorm(x[0] - 106.0)` |
| 30 | GALCENT_RGILBRAND | `",SgrA*"` | `iflag_true` | `degnorm(x[0] - 210.0 - 90.0*0.3819660113)` |
| 31 | GALEQU_IAU1958 | `",GP1958"` | `iflag_galequ` | `degnorm(x[0] - 150.0)` |
| 32 | GALEQU_TRUE | `",GPol"` | `iflag_galequ` | `degnorm(x[0] - 150.0)` |
| 33 | GALEQU_MULA | `",GPol"` | `iflag_galequ` | `degnorm(x[0] - 150.0 - 6.6666666667)` |
| 35 | TRUE_MULA | `",laSco"` | `iflag_true` | `degnorm(x[0] - 240.0)` |
| 36 | GALCENT_MULA_WILHELM | `",SgrA*"` | `iflag_true | SEFLG_EQUATORIAL` | `swi_armc_to_mc(x[0], eps) - 246.6666666667` (see below) |
| 39 | TRUE_SHEORAN | `",deCnc"` | `iflag_true` | `degnorm(x[0] - 103.49264221625)` |
| 40 | GALCENT_COCHRANE | `",SgrA*"` | `iflag_true` | `degnorm(x[0] - 270.0)` |

`x[0]` is the ecliptic longitude returned by `swe_fixstar()` (in degrees).

### Flag distinction: `iflag_true` vs `iflag_galequ`

- **`iflag_galequ`** `= iflag | SEFLG_TRUEPOS` — used for galactic pole queries (31, 32, 33).
  Always uses true position (no aberration, no light deflection). Required because the galactic
  pole is a geometric direction, not a luminous body.
- **`iflag_true`** — used for everything else. It starts from `iflag` (ephe + NONUT) and
  additionally passes through `SEFLG_TRUEPOS`, `SEFLG_NOABERR`, `SEFLG_NOGDEFL` from the
  caller's `otherflag`.

### GALCENT_MULA_WILHELM special case (index 36, sweph.c:3110–3121)

This index queries the Galactic Centre in **equatorial** coordinates and projects the right
ascension onto the ecliptic via the MC formula:

```
retflag = swe_fixstar(",SgrA*", tjd_et, iflag_true | SEFLG_EQUATORIAL, x)
eps = swi_epsiln(tjd_et, iflag) * RADTODEG   // obliquity in degrees
*daya = swi_armc_to_mc(x[0], eps)            // RA -> MC longitude (see below)
*daya = swe_degnorm(*daya - 246.6666666667)
```

`swi_armc_to_mc(armc, eps)` (swehouse.c:872–888): converts ARMC (RA of MC in degrees) to
MC ecliptic longitude. Formula: `mc = atan(tan(armc) / cos(eps))`, with quadrant adjustment
for armc > 90 && armc <= 270 → add 180°.

### Scope for Rust deferred implementation

All 12 indices listed in the table above require `swe_fixstar()` / catalog access.
**Defer all of them.** The non-fixed-star ayanamsas (indices 0–16, 18–26, 34, 37, 38, 41–46)
are fully computable without any star catalog. Index 34 (GALALIGN_MARDYKS) is table-based
with ECL_T0 despite its name.

---

## Sidereal Projection in the Calc Pipeline

### Where it fires

In `app_pos_rest()` (sweph.c:2777–2837), called from every planet calculation path.
The same pattern appears in `app_pos_etc_sun()` (sweph.c:6616–6642),
`app_pos_etc_moon()`/apside functions, and asteroid paths.

The pipeline saves J2000 equatorial Cartesian coordinates **before** precessing and nutating
to date:
```c
// sweph.c:2760–2762 (example from app_pos_etc_plan)
for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];   // save J2000 equatorial Cartesian + speed
```

Then, after nutation and ecliptic transformation (which populate `pdp->xreturn[6..17]`
with tropical ecliptic and `[18..23]` with equatorial), the sidereal block fires:

### Three dispatch branches (sweph.c:2811–2836)

**Branch 1: `SE_SIDBIT_ECL_T0` set** — project onto ecliptic of epoch t0
```c
swi_trop_ra2sid_lon(x2000, pdp->xreturn+6, pdp->xreturn+18, iflag)
// overwrites xreturn[6..11] (ecliptic Cartesian) and [18..23] (equatorial Cartesian)
```

**Branch 2: `SE_SIDBIT_SSY_PLANE` set** — project onto solar system equatorial plane
```c
swi_trop_ra2sid_lon_sosy(x2000, pdp->xreturn+6, iflag)
// overwrites xreturn[6..11] only
```

**Branch 3: neither** (default, most common) — traditional ayanamsa subtraction
```c
swi_cartpol_sp(pdp->xreturn+6, pdp->xreturn)   // ecliptic Cartesian -> polar [0..5]
// NOTE: swi_get_ayanamsa_ex disturbs cached sun data for TRUE_CHITRA; save/restore:
for i in 0..24: xxsv[i] = pdp->xreturn[i]
swi_get_ayanamsa_with_speed(pdp->teval, iflag, daya, serr)
for i in 0..24: pdp->xreturn[i] = xxsv[i]     // restore
pdp->xreturn[0] -= daya[0] * DEGTORAD          // subtract longitude in RADIANS
pdp->xreturn[3] -= daya[1] * DEGTORAD          // subtract speed in RADIANS
swi_polcart_sp(pdp->xreturn, pdp->xreturn+6)  // back to Cartesian
```

---

## `swi_trop_ra2sid_lon()` — ECL_T0 Projection (sweph.c:3273–3301)

```c
int swi_trop_ra2sid_lon(double *xin, double *xout, double *xoutr, int32 iflag)
```

**Input**: `xin[6]` = J2000 equatorial Cartesian (position + speed).
**Output**: `xout[6]` = ecliptic sidereal Cartesian (relative to ecliptic of t0), `xoutr[6]` = equatorial sidereal.

```
x = xin   // copy

// Step 1: precess J2000 -> t0 (equatorial)
if sip->t0 != J2000:
    swi_precess(x,   sip->t0, 0, J2000_TO_J)   // position
    swi_precess(x+3, sip->t0, 0, J2000_TO_J)   // speed (separate call)
xoutr = x   // equatorial sidereal output (in frame of t0)

// Step 2: equatorial t0 -> ecliptic t0
calc_epsilon(swed.sidd.t0, iflag, &oectmp)   // obliquity at t0
swi_coortrf2(x,   x, oectmp.seps, oectmp.ceps)   // position
if SEFLG_SPEED:
    swi_coortrf2(x+3, x+3, oectmp.seps, oectmp.ceps)  // speed

// Step 3: Cartesian ecliptic -> polar
swi_cartpol_sp(x, x)   // x[0] = ecliptic longitude in RADIANS

// Step 4: subtract ayanamsa initial value and apply correction
get_aya_correction(iflag, &corr, NULL)
x[0] -= sip->ayan_t0 * DEGTORAD        // RADIANS
x[0] = swe_radnorm(x[0] + corr * DEGTORAD)  // normalize in RADIANS

// Step 5: back to Cartesian
swi_polcart_sp(x, xout)
```

**Key unit note**: the subtraction `x[0] -= sip->ayan_t0 * DEGTORAD` operates in
**radians** (unlike the default pipeline which works in degrees then converts).

---

## `swi_trop_ra2sid_lon_sosy()` — SSY_PLANE Projection (sweph.c:3308–3356)

```c
int swi_trop_ra2sid_lon_sosy(double *xin, double *xout, int32 iflag)
```

**Input**: `xin[6]` = J2000 equatorial Cartesian. **Output**: `xout[6]` = Cartesian in
solar system equatorial plane (sidereal).

Constants (sweph.h:291–295):
```c
SSY_PLANE_NODE_E2000 = 107.582569 * DEGTORAD   /* ascending node of SSY plane on ecliptic J2000 */
SSY_PLANE_NODE       = 107.58883388 * DEGTORAD  /* same but at date */
SSY_PLANE_INCL       = 1.578701 * DEGTORAD      /* inclination of SSY plane */
```
The function uses `SSY_PLANE_NODE_E2000` (J2000 value), not `SSY_PLANE_NODE`.

```
oe = swed.oec2000   // J2000 obliquity (precomputed)
plane_node = SSY_PLANE_NODE_E2000
plane_incl = SSY_PLANE_INCL

// === Planet path ===
x = xin

// (a) equatorial J2000 -> ecliptic J2000
swi_coortrf2(x,   x, oe.seps, oe.ceps)
if SEFLG_SPEED: swi_coortrf2(x+3, x+3, oe.seps, oe.ceps)

// (b) ecliptic Cartesian -> polar
swi_cartpol_sp(x, x)

// (c) rotate by -plane_node (longitude shift)
x[0] -= plane_node

// (d) convert back to Cartesian and tilt to SSY plane
swi_polcart_sp(x, x)
swi_coortrf(x,   x,   plane_incl)   // position
swi_coortrf(x+3, x+3, plane_incl)   // speed

// (e) to polar in SSY plane
swi_cartpol_sp(x, x)

// === Zero-point path (vernal point of t0 in SSY plane) ===
x0 = [1.0, 0.0, 0.0]

if sip->t0 != J2000:
    swi_precess(x0, sip->t0, 0, J_TO_J2000)   // precess t0 -> J2000

swi_coortrf2(x0, x0, oe.seps, oe.ceps)   // equatorial J2000 -> ecliptic J2000
swi_cartpol(x0, x0)
x0[0] -= plane_node
swi_polcart(x0, x0)
swi_coortrf(x0, x0, plane_incl)
swi_cartpol(x0, x0)

// === Measure planet relative to zero point ===
x[0] -= x0[0]         // angle difference in RADIANS (polar)
x[0] *= RADTODEG      // now in DEGREES

// Apply ayan_t0 and correction (both in DEGREES)
get_aya_correction(iflag, &corr, NULL)
x[0] -= sip->ayan_t0
x[0] = swe_degnorm(x[0] + corr) * DEGTORAD   // normalize to [0,360) then back to radians

// Back to Cartesian
swi_polcart_sp(x, xout)
```

**Key unit note**: After `x[0] *= RADTODEG`, the subtraction `x[0] -= sip->ayan_t0` and
`swe_degnorm()` work in **degrees**. This is the opposite of `swi_trop_ra2sid_lon` which
stays in radians throughout. The final `* DEGTORAD` converts the result back to radians
before the Cartesian conversion.

---

## Nutation Interaction

`swi_get_ayanamsa_ex()` (internal) always forces `SEFLG_NONUT`. The ayanamsa computed
internally does NOT include nutation.

The public `swe_get_ayanamsa_ex()` adds `dpsi * RADTODEG` (nutation in longitude,
degrees) unless the caller passes `SEFLG_NONUT`.

In the **default calc pipeline** (Branch 3 above), the function called is
`swi_get_ayanamsa_with_speed()` which calls `swi_get_ayanamsa_ex()` — therefore the
ayanamsa subtracted from planetary positions does NOT include nutation. The planet
coordinates at that point already include nutation (via `swi_nutate()`), so subtracting
a non-nutated ayanamsa is correct.

In the **ECL_T0 path** (`swi_trop_ra2sid_lon`), nutation is also absent from the
ayanamsa (the function doesn't compute one). The planet's J2000 coordinates saved
before nutation (`x2000`) are used as input, so nutation is entirely bypassed in
the ecliptic-of-t0 projection.

---

## FP-Fidelity Notes

### 1. Minus-then-add vs subtract in Method 1 (sweph.c:3173)

The source reads:
```c
x[0] = -x[0] * RADTODEG + sip->ayan_t0;
```
NOT:
```c
x[0] = sip->ayan_t0 - x[0] * RADTODEG;   // WRONG order
```

Rust must replicate: `let lon = -x[0] * RADTODEG + sip.ayan_t0;`

IEEE 754 double evaluation: (1) negate `x[0]`, (2) multiply by `RADTODEG`, (3) add
`ayan_t0`. The alternative expression `ayan_t0 - x[0] * RADTODEG` does: (1) multiply
`x[0] * RADTODEG`, (2) subtract. These produce the same mathematical result but different
rounding if the intermediate is not exactly representable. Match the C source expression
exactly.

### 2. ayan_t0 unit context differs between functions

| Function | Where `ayan_t0` is subtracted | Unit at that point |
|---|---|---|
| `swi_get_ayanamsa_ex` | `x[0] = -x[0]*RADTODEG + sip->ayan_t0` | degrees |
| `swi_trop_ra2sid_lon` | `x[0] -= sip->ayan_t0 * DEGTORAD` | radians |
| `swi_trop_ra2sid_lon_sosy` | `x[0] -= sip->ayan_t0` | degrees (after `*RADTODEG`) |

These are three distinct code paths. Do not unify them.

### 3. `get_aya_correction` returns a near-zero signed value

The condition `if (*corr > 350) *corr -= 360` (sweph.c:2997) turns a longitude near 359°
into a small negative correction around -1°. This is intentional: the correction should
be near zero (a few arcseconds typically), not near 360°.

### 4. `swi_epsiln(t0, 0)` vs `swi_epsiln(tjd_et, iflag)`

In Method 1 at line 3168: `eps = swi_epsiln(t0, 0)` — the second argument is 0 (not the
current iflag). This means no special flag handling in the obliquity lookup. The same
applies in `get_aya_correction` at line 2991.

In Method 2 at line 3197: `eps = swi_epsiln(tjd_et, 0)` — also 0 for the date obliquity.

But in `GALCENT_MULA_WILHELM` at line 3116: `eps = swi_epsiln(tjd_et, iflag) * RADTODEG`
— uses the full `iflag`. Match these distinctions.

### 5. `swe_degnorm` vs `swe_radnorm`

- `swi_trop_ra2sid_lon` uses `swe_radnorm` (normalizes radians to [0, 2π))
- `swi_trop_ra2sid_lon_sosy` uses `swe_degnorm` then `* DEGTORAD`
- `swi_get_ayanamsa_ex` uses `swe_degnorm` on the final result

Do not interchange these; unit context differs at each site.

### 6. `ayan_t0` arithmetic expressions in the source

Several entries in the table are written as arithmetic in the C source (e.g.
`360 - 337.636111`, `23.250182778 - 0.004658035`). These are constant-folded at compile
time in C. In Rust, compute and store the folded constant to match the double value.
The values in the table above are provided as reference; independently verify each by
reading sweph.h:372–594 when writing the Rust data.

---

## `swed.sidd` — Runtime State

After `swe_set_sid_mode()`, the global `swed.sidd` holds:

```
swed.sidd.sid_mode   // full mode word (index | projection bits)
swed.sidd.t0         // epoch JD, in TT (or UT if t0_is_UT)
swed.sidd.ayan_t0    // ayanamsa at t0, degrees
swed.sidd.t0_is_UT   // whether t0 is UT
swed.ayana_is_set    // TRUE after first swe_set_sid_mode call
```

In the stateless Rust port, all of these are fields of `EphemerisConfig` / `SidData`
passed explicitly. No global state.

---

## Constants Summary

| Name | Value | Location | Meaning |
|---|---|---|---|
| `SE_SIDBITS` | 256 | swephexp.h:221 | Mask for the 8-bit ayanamsa index |
| `SE_SIDM_USER` | 255 | swephexp.h:286 | User-defined ayanamsa |
| `SE_NSIDM_PREDEF` | 47 | swephexp.h:288 | Number of predefined ayanamsas |
| `J2000` | 2451545.0 | sweph.h:67 | Julian Day of J2000.0 |
| `J1900` | 2415020.0 | sweph.h:69 | Julian Day of J1900.0 |
| `B1950` | 2433282.42345905 | sweph.h:68 | Julian Day of B1950 |
| `SEMOD_PREC_IAU_1976` | 1 | swephexp.h:511 | Precession model IAU 1976 |
| `SEMOD_PREC_NEWCOMB` | 11 | swephexp.h:521 | Newcomb precession model |
| `SEMOD_PREC_DEFAULT` | 9 (Vondrak 2011) | swephexp.h:522 | Current default model |
| `SEMOD_NUT_WOOLARD` | 5 | swephexp.h:536 | Woolard nutation model |
| `SEMOD_NUT_IAU_1980` | 1 | swephexp.h:531 | IAU 1980 nutation model |
| `SSY_PLANE_NODE_E2000` | 107.582569° in radians | sweph.h:291 | SSY plane node on ecliptic J2000 |
| `SSY_PLANE_INCL` | 1.578701° in radians | sweph.h:295 | SSY plane inclination |
| `DEGTORAD` | π/180 | sweodef.h:266 | Degree → radian |
| `RADTODEG` | 180/π | sweodef.h:265 | Radian → degree |

---

## Ayanamsas Deferred (Fixed-Star Dependent)

The following 12 indices require `swe_fixstar()` / star catalog access and are **not
implementable without a fixed-star subsystem**. Defer these:

| Index | Constant |
|---|---|
| 17 | SE_SIDM_GALCENT_0SAG |
| 27 | SE_SIDM_TRUE_CITRA |
| 28 | SE_SIDM_TRUE_REVATI |
| 29 | SE_SIDM_TRUE_PUSHYA |
| 30 | SE_SIDM_GALCENT_RGILBRAND |
| 31 | SE_SIDM_GALEQU_IAU1958 |
| 32 | SE_SIDM_GALEQU_TRUE |
| 33 | SE_SIDM_GALEQU_MULA |
| 35 | SE_SIDM_TRUE_MULA |
| 36 | SE_SIDM_GALCENT_MULA_WILHELM |
| 39 | SE_SIDM_TRUE_SHEORAN |
| 40 | SE_SIDM_GALCENT_COCHRANE |

All other 35 indices (0–16, 18–26, 34, 37, 38, 41–46) are computable from precession
and table data alone.

---

## References

| Source | Used in |
|---|---|
| sweph.h:336–596 | `struct aya_init` definition and full ayanamsa table |
| swephexp.h:221–288 | SE_SIDBIT_* and SE_SIDM_* constants |
| sweph.c:2861–2928 | `swe_set_sid_mode` |
| sweph.c:2930–2946 | `swe_get_ayanamsa_ex` (public, adds nutation) |
| sweph.c:2960–3000 | `get_aya_correction` |
| sweph.c:3002–3208 | `swi_get_ayanamsa_ex` (core, no nutation) |
| sweph.c:3210–3224 | `swi_get_ayanamsa_with_speed` |
| sweph.c:3273–3301 | `swi_trop_ra2sid_lon` (ECL_T0 projection) |
| sweph.c:3308–3356 | `swi_trop_ra2sid_lon_sosy` (SSY_PLANE projection) |
| sweph.c:2777–2836 | `app_pos_rest` — sidereal dispatch in planet pipeline |
| sweph.c:6616–6642 | Same pattern in `app_pos_etc_sun` |
| swehouse.c:872–888 | `swi_armc_to_mc` (used by GALCENT_MULA_WILHELM) |
| sweph.h:291–295 | SSY_PLANE_* constants |
