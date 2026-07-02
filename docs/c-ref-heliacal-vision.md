# C Reference: Heliacal — Vision & Sky Model (swehel.c part 1)

Porting reference for the vision/sky-brightness/atmospheric-extinction core of the heliacal
module: everything from the file-header constants (line 88) through `swe_heliacal_angle` (line
1695-1705). All line numbers refer to `swehel.c` unless stated otherwise. Lines 1715 onward
(`WidthMoon`, Yallop crescent criteria, `swe_heliacal_ut`/`swe_heliacal_pheno_ut` search drivers)
are covered by a sibling doc — do not duplicate.

This module has no direct C analogue split — everything here is `static` (module-private) except
the three `CALL_CONV` public entry points `swe_vis_limit_mag`, `swe_topo_arcus_visionis`,
`swe_heliacal_angle`. The public functions in this range are explicitly marked in swephexp.h as
**"secret, for Victor Reijs'"** (undocumented API, swephexp.h:680-682) — only `swe_vis_limit_mag`
is a fully public/documented entry point.

---

## Function Map

| C function | Location | Purpose |
|---|---|---|
| `mymin`/`mymax` | 155–167 | trivial min/max |
| `Tanh` | 170–173 | manual tanh via exp |
| `CVA` | 180–195 | contrast threshold angle (arcmin acuity) |
| `PupilDia` | 202–206 | pupil diameter vs. age/brightness |
| `OpticFactor` | 224–301 | optics/eye correction factors for vis. threshold |
| `DeterObject` | 305–336 | object-name string → planet number |
| `call_swe_calc` | 338–364 | **dead code** (`#if 0`), memoizing wrapper, never compiled |
| `call_swe_fixstar` | 368–375 | star-name-safe wrapper around `swe_fixstar` |
| `call_swe_fixstar_mag` | 379–394 | star-name-safe + memoized wrapper around `swe_fixstar_mag` |
| `call_swe_rise_trans` | 398–406 | star-name-safe wrapper around `swe_rise_trans` |
| `calc_rise_and_set` | 416–506 | fast Newton-iteration rise/set (planets/Sun/Moon, `\|lat\|`<63°) |
| `my_rise_trans` | 508–523 | dispatch: fast path vs. full `swe_rise_trans` |
| `RiseSet` | 535–546 | public-ish wrapper: object-name → `my_rise_trans` |
| `SunRA` | 553–583 | Sun's right ascension (memoized; **crude-approximation branch always active**, see §11) |
| `Kelvin` | 589–593 | °C → K |
| `TopoAltfromAppAlt` | 601–616 | apparent altitude → topocentric altitude (fast refraction) |
| `AppAltfromTopoAlt` | 626–654 | topocentric altitude → apparent altitude (Newton inversion of the above) |
| `HourAngle` | 662–672 | hour angle from altitude/declination/latitude |
| `ObjectLoc` | 683–726 | generic object position/angle accessor (7 angle modes) |
| `azalt_cart` | 737–771 | az/alt + cartesian apparent-altitude vector |
| `DistanceAngle` | 780–793 | haversine angular distance |
| `kW` | 801–807 | water-vapor extinction coefficient |
| `kOZ` | 815–841 | ozone extinction coefficient (memoized) |
| `kR` | 848–862 | Rayleigh extinction coefficient |
| `Sgn` | 864–869 | sign (±1, never 0) |
| `ka` | 881–927 | aerosol extinction coefficient (memoized) |
| `kt` | 940–957 | total/selectable extinction coefficient dispatcher |
| `Airmass` | 964–972 | Kasten-style airmass (reachable only if `staticAirmass≠0`, see §11) |
| `Xext` | 980–983 | optical-path length through an exponential atmosphere layer (non-thin) |
| `Xlay` | 991–996 | optical-path length through a thin/high layer (ozone) |
| `TempEfromTempS` | 1005–1008 | station → eye-height temperature via lapse rate |
| `PresEfromPresS` | 1016–1019 | station → eye-height pressure (barometric formula) |
| `Deltam` | 1033–1059 | total atmospheric extinction magnitude Δm (memoized) |
| `Bn` | 1073–1096 | night-sky background brightness (zodiacal+starlight, solar-cycle modulated) |
| `Magnitude` | 1106–1127 | object visual magnitude via `swe_pheno_ut`/`swe_fixstar_mag` |
| `fast_magnitude` | 1129–1152 | **dead code** (`#if 0`), memoizing wrapper, never compiled |
| `MoonsBrightness` | 1159–1164 | Moon's magnitude vs. distance/phase |
| `MoonPhase` | 1172–1181 | Moon phase angle from alt/az of Moon and Sun |
| `Bm` | 1186–1213 | moonlight sky brightness contribution |
| `Btwi` | 1218–1234 | twilight sky brightness contribution |
| `Bday` | 1246–1261 | daylight sky brightness contribution |
| `Bcity` | 1268–1274 | light-pollution brightness passthrough |
| `Bsky` | 1279–1307 | total sky brightness (dispatches day/twilight/night + adds moon/city/starlight) |
| `default_heliacal_parameters` | 1324–1361 | fills in defaults for `datm`/`dgeo`/`dobs` |
| `VisLimMagn` | 1382–1443 | limiting visual magnitude (contrast-threshold model) |
| `tolower_string_star` | 1446–1452 | lowercase a star name, preserving Bayer designation after `,` |
| `swe_vis_limit_mag` | 1464–1541 | **public**: limiting magnitude + object/Sun/Moon geometry |
| `TopoArcVisionis` | 1562–1599 | bisection: Sun depression at which object becomes visible |
| `swe_topo_arcus_visionis` | 1601–1610 | public (secret) wrapper around `TopoArcVisionis` |
| `HeliacalAngle` | 1636–1693 | 2-D bisection: optimal object-altitude / arcus-visionis pair |
| `swe_heliacal_angle` | 1695–1705 | public (secret) wrapper around `HeliacalAngle` |

---

## §1 Constants & flags

### 1.1 File-header constants (swehel.c:76–154)

```c
#define PLSV   0            /* if Planet, Lunar and Stellar Visibility formula is needed PLSV=1 */
#define criticalangle   0.0 /* [deg] */
#define BNIGHT   1479.0     /* [nL] */
#define BNIGHT_FACTOR   1.0
#define PI   M_PI
#define Min2Deg   (1.0 / 60.0)
#define SWEHEL_DEBUG  0
#define DONE  1
#define MaxTryHours   4
#define TimeStepDefault 1
#define LocalMinStep 8
```
`PLSV` gates a dead branch inside `HeliacalAngle` (§10 — always 0, so the branch never executes).
`criticalangle` (0.0°) is used only inside that dead `PLSV==1` branch. `DONE` is entirely unused
anywhere in the 3511-line file (dead). `MaxTryHours`, `TimeStepDefault`, `LocalMinStep` are used
only past line 1714 (sibling doc's territory — heliacal-event search loop), not in this range.

```c
/* time constants */
#define Y2D   365.25          /* [Day] */
#define D2Y   (1 / Y2D)       /* [Year] */
#define D2H   24.0            /* [Hour] */
#define H2S   3600.0          /* [sec] */
#define D2S   (D2H * H2S)     /* [sec] */
#define S2H   (1.0 / H2S)     /* [Hour] */
#define JC2D   36525.0        /* [Day] */
#define M2S   60.0            /* [sec] */
```
**Entirely dead**: every one of `Y2D, D2Y, D2H, H2S, D2S, S2H, JC2D, M2S` is referenced nowhere in
`swehel.c` outside its own `#define` (verified by grep across the whole 3511-line file — each
symbol's only occurrence is its own definition, or another dead constant's definition that chains
to it). Do not port; a Rust implementer can skip this block entirely.

```c
/* Determines which algorithms are used */
#define REFR_SINCLAIR    0
#define REFR_BENNETTH    1
#define FormAstroRefrac   REFR_SINCLAIR  /* for Astronomical refraction can be "bennetth" or "sinclair" */
#define GravitySource   2   /* 0=RGO, 1=Wikipedia, 2=Exp. Suppl. 1992, 3=van der Werf */
#define REarthSource   1    /* 0=RGO (constant), 1=WGS84 method */
```
`FormAstroRefrac`/`REFR_SINCLAIR`/`REFR_BENNETTH` are used only past line 1714 (sibling doc). Note
this fixes the algorithm selection at compile time — there is no runtime branch. `GravitySource`
and `REarthSource` are **entirely dead** (never referenced anywhere in the file) — they document
which physical-constant source was chosen for constants that are hardcoded directly (`Ra`, `Rb`
below), not selected via these defines. Skip both in the port.

```c
#define StartYear   1820                    /* [year] */
#define Average   1.80546834626888          /* [msec/cy] */
#define Periodicy   1443.67123144531        /* [year] */
#define Amplitude   3.75606495492684        /* [msec] */
#define phase   0                           /* [deg] */
#define MAX_COUNT_SYNPER           5        /* search within 10 synodic periods */
#define MAX_COUNT_SYNPER_MAX 1000000        /* high, so there is not max count */
#define AvgRadiusMoon  (15.541 / 60)        /* [Deg] at 2007 CE or BCE */
```
All of these (`StartYear`, `Average`, `Periodicy`, `Amplitude`, `phase`, `MAX_COUNT_SYNPER`,
`MAX_COUNT_SYNPER_MAX`) are used only past line 1714 — day-length/ΔT drift model and heliacal
search-loop iteration caps (sibling doc's territory). `AvgRadiusMoon` is likewise first used at
line 1729 (`WidthMoon`, sibling doc). None are referenced in the functions covered by this doc.

```c
/* WGS84 ellipsoid constants — http://w3sli.wcape.gov.za/Surveys/Mapping/wgs84.htm */
#define Ra   6378136.6     /* [m] */
#define Rb   6356752.314   /* [m] */
```
`Ra` (equatorial radius) **is** used in this range: `MoonsBrightness` (§6) divides distance by
`Ra/1000` (km). `Rb` is unused anywhere in the file (dead — the polar radius is defined but never
consulted; flattening/oblateness is not modeled anywhere in this module).

```c
/* choices in Schaefer's model */
#define nL2erg   (1.02E-15)
#define erg2nL   (1 / nL2erg)              /* erg2nL to nLambert */
#define MoonDistance   384410.4978         /* [km] */
#define scaleHwater   3000.0     /* [m] Ricchiazzi [1997] 8200 Schaefer [2000] */
#define scaleHrayleigh   8515.0  /* [m] Su [2003] 8200 Schaefer [2000] */
#define scaleHaerosol   3745.0   /* [m] Su [2003] 1500 Schaefer [2000] */
#define scaleHozone   20000.0    /* [m] Schaefer [2000] */
#define astr2tau   0.921034037197618   /* LN(10 ^ 0.4) */
#define tau2astr   1 / astr2tau
```
All of these are live and used throughout §4/§6. **FP hazard**: `tau2astr` is defined as
`1 / astr2tau` *without parentheses* — because `/` binds normally this evaluates as
`(1) / (astr2tau)` regardless (no ambiguity here since it's a single division), so textually
replicate as `1.0 / astr2tau`, computed once as a constant (matches C's compile-time constant
folding — no repeated re-division at each use site needed, but bitwise it's the same value either
way).

```c
/* meteorological constants */
#define C2K   273.15         /* [K] */
#define DELTA   18.36
#define TempNulDiff   0.000001
#define PressRef   1000      /* [mbar] */
#define MD   28.964           /* [kg] Mol weight of dry air van der Werf */
#define MW   18.016           /* [kg] Mol weight of water vapor */
#define GCR   8314.472        /* [L/kmol/K] van der Werf */
#define LapseSA   0.0065      /* [K/m] standard atmosphere */
#define LapseDA   0.0098      /* [K/m] dry adiabatic */
```
`C2K` is live (used by `Kelvin`, §3). **`DELTA`, `TempNulDiff`, `PressRef`, `MD`, `MW`, `GCR`,
`LapseDA` are entirely dead** — none is referenced anywhere in the 3511-line file. `LapseSA` is
live (used in `Deltam`/`Bn` to convert station temperature to eye-height temperature via
`TempEfromTempS`, §3/§4). Do not port the dead ones.

```c
/* lowest apparent altitude to provide */
#define LowestAppAlt   -3.5   /* [Deg] */

/* optimization delta */
#define epsilon   0.001
/* for Airmass usage */
#define staticAirmass   0     /* use staticAirmass=1 instead depending on difference k's */

/* optic stuff */
#define GOpticMag   1         /* telescope magnification */
#define GOpticTrans   0.8     /* telescope transmission */
#define GBinocular   1        /* 1-binocular 0=monocular */
#define GOpticDia   50        /* telescope diameter [mm] */
```
`LowestAppAlt` is live (clamp bound in `TopoAltfromAppAlt`/`AppAltfromTopoAlt`, §3). `epsilon` is
live (bisection tolerance in `TopoArcVisionis`, §9). `staticAirmass` is live as a **compile-time
constant** consulted via `if (staticAirmass == 0)` in `Deltam` (§4) — since it is always `0`, the
`else` branch (which calls `Airmass()`) is unreachable dead code under the shipped build; `Airmass`
itself remains a real (but effectively unreachable) function. **`GOpticMag`, `GOpticTrans`,
`GBinocular`, `GOpticDia` are entirely dead** — never referenced anywhere in the file; the actual
optics defaults are hardcoded directly inside `default_heliacal_parameters` (§7, `dobs[2]=1`,
`dobs[3]=1`) rather than sourced from these constants.

### 1.2 `SE_HELFLAG_*` and related flags (swephexp.h:423–471)

| Constant | Value | Meaning |
|---|---|---|
| `SE_HELIACAL_RISING` | 1 | = `SE_MORNING_FIRST` |
| `SE_HELIACAL_SETTING` | 2 | = `SE_EVENING_LAST` |
| `SE_EVENING_FIRST` | 3 | |
| `SE_MORNING_LAST` | 4 | |
| `SE_ACRONYCHAL_RISING` | 5 | not implemented |
| `SE_ACRONYCHAL_SETTING` | 6 | = `SE_COSMICAL_SETTING`; not implemented |
| `SE_HELFLAG_LONG_SEARCH` | 128 | |
| `SE_HELFLAG_HIGH_PRECISION` | 256 | gates ephemeris-based (vs. crude) computations throughout this module |
| `SE_HELFLAG_OPTICAL_PARAMS` | 512 | caller is supplying `dobs[2..5]` (optics) explicitly |
| `SE_HELFLAG_NO_DETAILS` | 1024 | |
| `SE_HELFLAG_SEARCH_1_PERIOD` | 2048 (`1<<11`) | |
| `SE_HELFLAG_VISLIM_DARK` | 4096 (`1<<12`) | ignore Sun/twilight brightness entirely (§8) |
| `SE_HELFLAG_VISLIM_NOMOON` | 8192 (`1<<13`) | ignore Moon brightness contribution |
| `SE_HELFLAG_VISLIM_PHOTOPIC` | 16384 (`1<<14`) | force photopic (cone) vision model; undocumented/test-only |
| `SE_HELFLAG_VISLIM_SCOTOPIC` | 32768 (`1<<15`) | force scotopic (rod) vision model; undocumented/test-only |
| `SE_HELFLAG_AV` / `SE_HELFLAG_AVKIND_VR` | 65536 (`1<<16`) | |
| `SE_HELFLAG_AVKIND_PTO` | `1<<17` | |
| `SE_HELFLAG_AVKIND_MIN7` | `1<<18` | |
| `SE_HELFLAG_AVKIND_MIN9` | `1<<19` | |
| `SE_HELFLAG_AVKIND` | OR of the 4 above | |
| `SE_PHOTOPIC_FLAG` | 0 | return-value bit meaning for `swe_vis_limit_mag` |
| `SE_SCOTOPIC_FLAG` | 1 | ditto |
| `SE_MIXEDOPIC_FLAG` | 2 | ditto (OR-ed in, "near the scotopic/photopic transition") |

Note: `swephexp.h:453-467` also defines a **dead, `#if 0`-wrapped** parallel set of
`SE_HELIACAL_LONG_SEARCH`/`SE_HELIACAL_HIGH_PRECISION`/etc. constants ("unused and redundant" per
the header comment) — these never compile; only the `SE_HELFLAG_*` names above are real. Do not
port the `SE_HELIACAL_*` flag-bit aliases.

Also relevant: `SE_ECL2HOR = 0`, `SE_EQU2HOR = 1` (swephexp.h:364-365) — `swe_azalt` input-frame
selector, used throughout §3/§7.

`SIMULATE_VICTORVB` (swephexp.h:451, `#define SIMULATE_VICTORVB 1`) is **always defined** (there
is no build configuration in this codebase that undefines it) — see §11 for the significant
consequences this has on which code paths in `SunRA`, `ka`, and `default_heliacal_parameters`
actually execute.

### 1.3 Input array layouts (documented only by convention/usage, not by header comment)

**`dgeo[3]`** (observer geographic position):
- `dgeo[0]` — geographic longitude, degrees (east positive)
- `dgeo[1]` — geographic latitude, degrees (north positive)
- `dgeo[2]` — eye height above sea level, meters

**`datm[4]`** (atmospheric conditions):
- `datm[0]` — atmospheric pressure, mbar (`<= 0` triggers ISA-model default, §7 `default_heliacal_parameters`)
- `datm[1]` — temperature, °C (station/sea-level value; `0` triggers a lapse-rate-based default when `datm[0] <= 0`)
- `datm[2]` — relative humidity, % (`0` triggers default `40`; range clamped internally in some paths, see §11)
- `datm[3]` — either the total broadband extinction coefficient (VR<1 usage) or Meteorological Range in km (VR≥1 usage) — dual-use field, dispatched inside `ka()` (§4) by comparing against `1`

**`dobs[6]`** (observer/optics definition):
- `dobs[0]` — observer age, years (default 36 via `default_heliacal_parameters`)
- `dobs[1]` — Snellen ratio / visual acuity (default 1)
- `dobs[2]` — `Binocular` flag (1 = binocular, 0 = monocular)
- `dobs[3]` — `OpticMag` telescope magnification (1 = "use eye", i.e. naked-eye mode)
- `dobs[4]` — `OpticDia` telescope/instrument diameter, mm
- `dobs[5]` — `OpticTrans` telescope transmission factor (0–1)

`dobs[2..5]` are zeroed by `default_heliacal_parameters` unless `SE_HELFLAG_OPTICAL_PARAMS` is
set (§7); when `dobs[3] == 0` after that, the function further forces `dobs[2]=1, dobs[3]=1`
(binocular=1/eye-mode) as the naked-eye default.

### 1.4 `default_heliacal_parameters` (swehel.c:1324) — see full treatment in §7.

---

## §2 Object determination & ephemeris call wrappers

### `DeterObject(char *ObjectName)` — swehel.c:305–336

```c
static int32 DeterObject(char *ObjectName)
{
  char s[AS_MAXCH];
  char *sp;
  int32 ipl;
  strcpy(s, ObjectName);
  for (sp = s; *sp != '\0'; sp++)
    *sp = tolower(*sp);
  if (strncmp(s, "sun", 3) == 0) return SE_SUN;
  if (strncmp(s, "venus", 5) == 0) return SE_VENUS;
  if (strncmp(s, "mars", 4) == 0) return SE_MARS;
  if (strncmp(s, "mercur", 6) == 0) return SE_MERCURY;
  if (strncmp(s, "jupiter", 7) == 0) return SE_JUPITER;
  if (strncmp(s, "saturn", 6) == 0) return SE_SATURN;
  if (strncmp(s, "uranus", 6) == 0) return SE_URANUS;
  if (strncmp(s, "neptun", 6) == 0) return SE_NEPTUNE;
  if (strncmp(s, "moon", 4) == 0) return SE_MOON;
  if ((ipl = atoi(s)) > 0) { ipl += SE_AST_OFFSET; return ipl; }
  return -1;
}
```
Lower-cases a private copy, then matches **prefixes** in the fixed order above (note: "mercur"
matches "mercury"/"mercurius" etc.; "neptun" matches "neptune"/"neptunus"). If none match, tries
`atoi(s)` — a leading-numeric string is treated as an asteroid catalog number and offset by
`SE_AST_OFFSET`; if `atoi` returns 0 (non-numeric or literal "0"), returns `-1` meaning "not a
built-in planet — treat `ObjectName` as a fixed-star name instead." Pure string logic, no
ephemeris/global state. Rust port: a `match`/prefix-check returning an `Option<Body>`-like enum,
falling through to fixed-star lookup on `None`.

### `call_swe_fixstar(char *star, double tjd, int32 iflag, double *xx, char *serr)` — swehel.c:368–375
Pure delegation to `swe_fixstar`, copying `star` into a local buffer first because `swe_fixstar`
may overwrite/normalize the string in place (comment: "avoids problems with star name string that
may be overwritten"). No caching, no numeric transformation. Rust port: irrelevant — Rust owns
`&str`/`String` and has no aliasing/mutation hazard here; just call the stateless fixed-star
lookup directly.

### `call_swe_fixstar_mag(char *star, double *mag, char *serr)` — swehel.c:379–394
```c
static int32 call_swe_fixstar_mag(char *star, double *mag, char *serr)
{
  int32 retval;
  char star2[AS_MAXCH];
  static TLS double dmag;
  static TLS char star_save[AS_MAXCH];
  if (strcmp(star, star_save) == 0) { *mag = dmag; return OK; }
  strcpy(star_save, star);
  strcpy(star2, star);
  retval = swe_fixstar_mag(star2, &dmag, serr);
  *mag = dmag;
  return retval;
}
```
**Cross-call static state**: `star_save`/`dmag` persist across calls (thread-local via `TLS`).
Pure memoization keyed on exact string equality of the star name (not epoch-dependent — fixed-star
magnitude is treated as constant, which is astronomically correct: `swe_fixstar_mag` doesn't take
a time argument either). Cache hit skips the `swe_fixstar_mag` call entirely. **This is safe to
drop in a stateless Rust port** — recomputing `swe_fixstar_mag(star)` every call gives bitwise
identical results (no hidden epoch/state dependency), just costs a redundant catalog lookup.

### `call_swe_rise_trans(double tjd, int32 ipl, char *star, int32 helflag, int32 eventtype, double *dgeo, double atpress, double attemp, double *tret, char *serr)` — swehel.c:398–406
```c
static int32 call_swe_rise_trans(double tjd, int32 ipl, char *star, int32 helflag, int32 eventtype, double *dgeo, double atpress, double attemp, double *tret, char *serr)
{
  int32 retval;
  int32 iflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  char star2[AS_MAXCH];
  strcpy(star2, star);
  retval = swe_rise_trans(tjd, ipl, star2, iflag, eventtype, dgeo, atpress, attemp, tret, serr);
  return retval;
}
```
`helflag → iflag` conversion: **only the ephemeris-selector bits** (`SEFLG_JPLEPH|SEFLG_SWIEPH|
SEFLG_MOSEPH`) survive; every other `helflag` bit (precision, vislim mode, etc.) is masked out
before calling `swe_rise_trans`. This narrow mask pattern (`helflag & SEFLG_EPHMASK`-equivalent,
though written out explicitly as the 3 individual bits rather than a named `SEFLG_EPHMASK`
constant) recurs identically in `calc_rise_and_set`, `SunRA`, `ObjectLoc`, `azalt_cart`,
`Magnitude` (see each below) — it is the module's standard "extract just the ephemeris backend
choice" idiom. Rust port: this maps directly to selecting which `Ephemeris` backend/config to use;
no other flag needs to cross this boundary.

### `calc_rise_and_set(double tjd_start, int32 ipl, double *dgeo, double *datm, int32 eventflag, int32 helflag, double *trise, char *serr)` — swehel.c:416–506
Fast Newton-style rise/set specifically for the heliacal module (distinct from, but structurally
similar to, `swecl.c`'s `rise_set_fast` — see `docs/c-ref-riseset.md` §3; **do not merge/dedupe
these two into one Rust function** unless a review of both confirms the underlying math and
disc-radius handling truly coincide — they diverge in disc-radius source values, iteration count,
and day-anchoring logic, see below).

1. `iflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)`. If `!(helflag &
   SE_HELFLAG_HIGH_PRECISION)`: `iflag |= SEFLG_NONUT|SEFLG_TRUEPOS` (perf shortcut: skip nutation,
   use "true" i.e. non-aberration/light-time-adjusted position — an approximation used throughout
   this module whenever high precision isn't requested).
2. `tjdnoon = (int)tjd0 - dgeo[0]/15.0/24.0` (local-noon estimate: truncate `tjd0` to integer JD,
   then shift by longitude/15/24 — converting longitude degrees to a day-fraction via the
   15°/hour sidereal-ish rate).
3. `iflag |= SEFLG_EQUATORIAL`; compute Sun (`xs`) and object (`xx`) at `tjd0` via `swe_calc_ut`.
   `tjdnoon -= swe_degnorm(xs[0]-xx[0])/360.0 + 0` (adjust the noon estimate by the RA difference
   between Sun and object, converted to a day-fraction at the crude 360°/day rate — the `+ 0` is
   a no-op, presumably vestigial).
4. `swe_azalt(tjd0, SE_EQU2HOR, dgeo, datm[0], datm[1], xx, xaz)` — is the object currently above
   (`xaz[2] > 0`, apparent altitude) or below the horizon?
5. **Day-anchoring while-loops** (lines 441–457): for rise events, if currently above horizon,
   push `tjdnoon` into `(tjd0+0.5, tjd0+1.5)`; if below, into `(tjd0, tjd0+1)`. For set events (the
   `else` branch), mirrored with different half-open interval edges (`(tjd0-0.5,tjd0+0.5]` above
   horizon vs. `(tjd0-1,tjd0]` below horizon — note the asymmetric strict/non-strict comparisons:
   `>` vs `<` on the lower/upper bound, replicate exactly). These are `while` loops (not `if`),
   incrementing/decrementing by whole days until inside the target window — always terminates in
   ≤2 iterations for any reasonable `tjd0`/`tjdnoon` separation.
6. Recompute `xx` at `tjdnoon` (position "at local noon", used only to get a representative
   declination `xx[1]` for the semi-diurnal-arc formula — not iterated further at this stage).
7. Disc radius `rdi`: `0` if `SE_BIT_DISC_CENTER` set; else `asin(696000000.0 / 1.49597870691e+11
   / xx[2]) / DEGTORAD` for `SE_SUN` (Sun radius 696,000 km, hardcoded — **note**: uses
   `1.49597870691e11` here, a *different* AU value than the riseset module's `AUNIT =
   1.49597870700e11`, DE431 — this is an intentional/historical inconsistency in the C source;
   replicate the literal `1.49597870691e+11` for FP fidelity in this function specifically, do
   not substitute the riseset module's `AUNIT`), or `asin(1737000.0 / 1.49597870691e+11 / xx[2]) /
   DEGTORAD` for `SE_MOON` (Moon radius 1737 km); `0` for any other body (no disc-radius
   correction for planets/stars in this fast path).
8. `rh = -(34.5/60.0 + rdi)` — target true-altitude at rise/set: a **fixed** refraction constant
   of 34.5 arcminutes (not looked up from `datm`/pressure/temperature at all — this differs from
   the riseset module's dynamic `swe_refrac_extended` call) plus the disc radius, negated (the
   body's center must be `rh` degrees *below* the true horizon for the limb to touch it).
9. `sda = acos(-tan(dgeo[1]*DEGTORAD) * tan(xx[1]*DEGTORAD)) * RADTODEG` — semi-diurnal arc,
   degrees (no clamping for circumpolar cases here — if the `acos` argument is out of `[-1,1]`,
   this yields `NaN`, propagating through to a `NaN` result; **no circumpolar guard**, unlike
   `swecl.c`'s `rise_set_fast`, §3.2 of the riseset doc).
10. `tjdrise = tjdnoon - sda/360.0` (rise) or `tjdnoon + sda/360.0` (set) — rough estimate.
11. **Refinement loop**, `iflag = epheflag|SEFLG_SPEED|SEFLG_EQUATORIAL`; `SEFLG_TOPOCTR` added
    only for `SE_MOON`; `SEFLG_NONUT|SEFLG_TRUEPOS` added again unless high-precision. Exactly
    **2 iterations** (`for (i=0;i<2;i++)`, fixed, no convergence check):
    - `swe_calc_ut(tjdrise, ipl, iflag, xx, serr)` (with speed).
    - `swe_azalt(tjdrise, SE_EQU2HOR, dgeo, datm[0], datm[1], xx, xaz)`.
    - `dfac = 1/365.25` (one Julian year in days — reused here purely as a small time-step
      constant, **not** as a calendar year; a slightly-obscure choice of step size ≈ 0.00274 days
      ≈ 3.9 minutes).
    - `xx[0] -= xx[3]*dfac; xx[1] -= xx[4]*dfac` — back-propagate RA/decl by one `dfac` step using
      the speed components `xx[3]`(dRA/day)/`xx[4]`(ddecl/day) obtained from `SEFLG_SPEED`
      (linear extrapolation, not a second `swe_calc_ut` call — cheaper).
    - `swe_azalt(tjdrise - dfac, SE_EQU2HOR, dgeo, datm[0], datm[1], xx, xaz2)` — altitude at the
      back-stepped time, using the linearly-extrapolated position (not a fresh ephemeris call).
    - `tjdrise -= (xaz[1]-rh) / (xaz[1]-xaz2[1]) * dfac` — secant-style update: `(xaz[1]-rh)` is
      the altitude error against the target `rh`; `(xaz[1]-xaz2[1])/dfac` is the finite-difference
      altitude rate; dividing error by rate gives the time correction (no explicit clamp on the
      step size here, unlike the riseset module's `dt` clamp to `[-0.1,0.1]`).
12. `*trise = tjdrise; return retc` (`retc` is always `OK` on this path — the two early-`return
    ERR` are only from `swe_calc_ut` failures on the Sun/object calls).

### `my_rise_trans(double tjd, int32 ipl, char *starname, int32 eventtype, int32 helflag, double *dgeo, double *datm, double *tret, char *serr)` — swehel.c:508–523
Dispatcher: if `starname` non-empty, resolve `ipl = DeterObject(starname)` first (so a numeric or
recognized-planet-name "star" string still routes through the planet path). Then: **if `ipl != -1`
and `|dgeo[1]| < 63`** (latitude threshold, degrees) → `calc_rise_and_set` (fast path, §2 above);
**else** → `call_swe_rise_trans` (full/slow path, §2 above) — this covers both fixed stars
(`ipl==-1`) and high-latitude planets. Note the threshold is `63°`, distinct from `swecl.c`'s
`rise_set_fast` gate of `60°`/`65°` (riseset doc §4) — **do not unify the two latitude constants**,
they are independently-tuned per module.

### `RiseSet(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 RSEvent, int32 helflag, int32 Rim, double *tret, char *serr)` — swehel.c:535–546
Thin wrapper: `Rim==0` → OR in `SE_BIT_DISC_CENTER` into the event flags (rise/set of disc center,
not limb). Resolves `ObjectName` via `DeterObject`; if recognized, calls `my_rise_trans` with
`starname=""`; else calls `my_rise_trans` with `ipl=-1` and the literal `ObjectName` as star name.
Pure dispatch, no computation of its own.

---

## §3 Meteorological & coordinate helpers

### `mymin`/`mymax` — swehel.c:155–167
Trivial: `a<=b ? a : b` / `a>=b ? a : b`. Note the `<=`/`>=` (not `<`/`>`) — on ties, returns `a`
in both cases; irrelevant for floats except NaN propagation semantics (if `a` is NaN, `a<=b` is
false, so `mymin` returns `b`; if `b` is NaN, `a<=b` is false too, so `mymin` also returns `b`
— i.e. `mymin(NaN, x) = x` always, `mymin(x, NaN) = NaN` always; **not commutative under NaN**,
replicate exactly with a matching non-commutative Rust helper, not `f64::min` which has different
NaN handling).

### `Tanh(double x)` — swehel.c:170–173
```c
return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
```
Manual hyperbolic tangent via two `exp` calls (not libm `tanh`) — for FP fidelity, replicate this
exact expression (4 `exp` calls total per invocation: note `exp(x)` and `exp(-x)` are each
computed twice, once in numerator once in denominator — do **not** hoist into locals, replicate
literally) rather than calling Rust's `f64::tanh`.

### `Kelvin(double Temp)` — swehel.c:589–593
`Temp + C2K` (`C2K = 273.15`). Trivial.

### `TopoAltfromAppAlt(double AppAlt, double TempE, double PresE)` — swehel.c:601–616
```c
double R = 0, retalt = 0;
if (AppAlt >= LowestAppAlt) {                 // LowestAppAlt = -3.5
  if (AppAlt > 17.904104638432)
    R = 0.97 / tan(AppAlt * DEGTORAD);
  else
    R = (34.46 + 4.23*AppAlt + 0.004*AppAlt*AppAlt) / (1 + 0.505*AppAlt + 0.0845*AppAlt*AppAlt);
  R = (PresE - 80) / 930 / (1 + 0.00008*(R+39)*(TempE-10)) * R;
  retalt = AppAlt - R * Min2Deg;              // Min2Deg = 1/60
} else {
  retalt = AppAlt;
}
return retalt;
```
Refraction model (Sinclair/Bennett-style rational approximation): apparent → topocentric (true)
altitude, degrees. `R` is refraction in **arcminutes** before the final `* Min2Deg` conversion to
degrees. Below `LowestAppAlt = -3.5°`, no refraction correction is applied at all (returns
`AppAlt` unchanged) — this is the "lowest apparent altitude" clamp mentioned in the file-header
constants. The `17.904104638432°` breakpoint switches between a simple cotangent formula (high
altitude) and a rational polynomial fit (low altitude, where simple refraction models break down
near the horizon). Pressure/temperature correction factor `(PresE-80)/930/(1+0.00008*(R+39)*
(TempE-10))` is applied multiplicatively to the geometric `R` in both branches identically.

### `AppAltfromTopoAlt(double TopoAlt, double TempE, double PresE, int32 helflag)` — swehel.c:626–654
Iterative Newton-style inversion of `TopoAltfromAppAlt` (comment: "call this instead of
`swe_azalt()`, because it is faster (lower precision is required)"). `nloop = 2` normally, `5` if
`SE_HELFLAG_HIGH_PRECISION`.
```c
newAppAlt = TopoAlt; newTopoAlt = 0.0; oudAppAlt = newAppAlt; oudTopoAlt = newTopoAlt;
for (i = 0; i <= nloop; i++) {                 // note: <=, so nloop+1 iterations (3 or 6)
  newTopoAlt = newAppAlt - TopoAltfromAppAlt(newAppAlt, TempE, PresE);
  verschil = newAppAlt - oudAppAlt;
  oudAppAlt = newTopoAlt - oudTopoAlt - verschil;
  if ((verschil != 0) && (oudAppAlt != 0))
    verschil = newAppAlt - verschil * (TopoAlt + newTopoAlt - newAppAlt) / oudAppAlt;
  else
    verschil = TopoAlt + newTopoAlt;
  oudAppAlt = newAppAlt; oudTopoAlt = newTopoAlt; newAppAlt = verschil;
}
retalt = TopoAlt + newTopoAlt;
if (retalt < LowestAppAlt) retalt = TopoAlt;
return retalt;
```
**Iteration count is `nloop+1`, not `nloop`** (loop condition `i <= nloop`) — 3 iterations
normally, 6 under high precision; replicate the off-by-one exactly. This is a secant-like scheme
reusing the previous two iterates (`oudAppAlt`/`oudTopoAlt` hold prior-iteration state — genuine
per-call local state, not a global/static; safe to port directly as a local loop). Final clamp:
if the converged result is below `LowestAppAlt`, fall back to returning `TopoAlt` unchanged
(same rationale as the forward function's cutoff).

### `HourAngle(double TopoAlt, double TopoDecl, double Lat)` — swehel.c:662–672
```c
double ha = (sin(Alti) - sin(Lati)*sin(decli)) / cos(Lati) / cos(decli);
if (ha < -1) ha = -1;
if (ha > 1) ha = 1;
return acos(ha) / DEGTORAD / 15.0;
```
Standard spherical-astronomy hour-angle formula; `Alti`/`decli`/`Lati` are `TopoAlt`/`TopoDecl`/
`Lat` converted to radians via `* DEGTORAD` first. Clamped to `[-1,1]` before `acos` (guards
float rounding pushing the argument minutely outside domain). Result is in **hours** (`acos(...)`
in degrees `/15`, not radians — the `/DEGTORAD` converts the `acos` radian result to degrees
first, then `/15.0` converts degrees to hours at the sidereal 15°/hour rate).

### `DistanceAngle(double LatA, double LongA, double LatB, double LongB)` — swehel.c:780–793
Haversine great-circle angular distance, **radians in and out**:
```c
double dlon = LongB - LongA, dlat = LatB - LatA;
double sindlat2 = sin(dlat/2), sindlon2 = sin(dlon/2);
double corde = sindlat2*sindlat2 + cos(LatA)*cos(LatB)*sindlon2*sindlon2;
if (corde > 1) corde = 1;
return 2 * asin(sqrt(corde));
```
Clamp `corde` to `≤1` before `sqrt`/`asin` (float-rounding guard, standard haversine hygiene). All
callers in this file pass degrees-to-radians-converted args (`Alt*DEGTORAD`, `Azi*DEGTORAD`) and
divide the radian result `/DEGTORAD` back to degrees at the call site — `DistanceAngle` itself is
unit-agnostic radians-in/radians-out.

### `TempEfromTempS(double TempS, double HeightEye, double Lapse)` — swehel.c:1005–1008
`TempS - Lapse * HeightEye`. Trivial linear lapse-rate correction; `Lapse` is always passed as
`LapseSA = 0.0065` K/m by every caller in this file.

### `PresEfromPresS(double TempS, double Press, double HeightEye)` — swehel.c:1016–1019
```c
return Press * exp(-9.80665 * 0.0289644 / (Kelvin(TempS) + 3.25 * HeightEye / 1000) / 8.31441 * HeightEye);
```
Barometric-formula pressure correction from station pressure to eye height. Coefficients:
`9.80665` (standard gravity, m/s²), `0.0289644` (molar mass of dry air, kg/mol — note this is a
*different* value than the file-header `MD = 28.964` [kg, presumably intended as g/mol] constant,
and `MD` is dead/unused anyway, §1), `8.31441` (universal gas constant, J/(mol·K) — again distinct
from the dead `GCR = 8314.472` header constant, different units/unused). `Kelvin(TempS) + 3.25 *
HeightEye/1000` approximates the mean temperature of the air column between station and eye
height (station temp plus a fixed lapse adjustment for half the column, roughly). Evaluate the
expression exactly in this left-to-right grouping: `exp(-9.80665 * 0.0289644 / (denom) / 8.31441 *
HeightEye)` — i.e. `((-9.80665 * 0.0289644) / denom) / 8.31441) * HeightEye` as the exponent,
standard left-to-right `*`/`/` associativity.

---

## §5 Optics & vision helpers

(Presented before §4 in file order but grouped per the doc's required section numbering — §4
Atmospheric extinction depends conceptually on nothing here, and §5 depends on nothing in §4
either; order in this doc follows the requested §-numbering, not file order.)

### `CVA(double B, double SN, int32 helflag)` — swehel.c:180–195
"Contrast visual acuity" / critical angle, degrees. Citation: Schaefer, *Astronomy and the limits
of vision*, Archaeoastronomy, 1993.
```c
AS_BOOL is_scotopic = FALSE;
if (B < 1394) is_scotopic = TRUE;                     // NOT BNIGHT (1479) — see note below
if (helflag & SE_HELFLAG_VISLIM_PHOTOPIC) is_scotopic = FALSE;
if (helflag & SE_HELFLAG_VISLIM_SCOTOPIC) is_scotopic = TRUE;
if (is_scotopic)
  return mymin(900, 380 / SN * pow(10, (0.3 * pow(B, (-0.29))))) / 60.0 / 60.0;
else
  return (40.0 / SN) * pow(10, (8.28 * pow(B, (-0.29)))) / 60.0 / 60.0;
```
**Threshold quirk**: the comment explicitly says `//if (B < BNIGHT)` was replaced by the literal
`1394` "to make the function continuous" — `BNIGHT` is `1479.0`, a *different* value from the
`1645`/`1645` thresholds used in `OpticFactor` and `VisLimMagn` (§5/§8) — **all three
scotopic/photopic brightness thresholds in this file (`1394` here, `1645` in `OpticFactor`, `1645`
in `VisLimMagn`) are distinct hand-tuned literals, not all equal to `BNIGHT`; do not consolidate
them into a single named constant** — replicate each literal exactly at its own call site.
`SE_HELFLAG_VISLIM_PHOTOPIC`/`SCOTOPIC` force-override the brightness-based default (undocumented
test flags, §1.2). Scotopic branch clamps the arcsecond-scale result to `≤900` before the final
`/3600` (arcsec→deg) conversion; photopic branch has no such clamp.

### `PupilDia(double Age, double B)` — swehel.c:202–206
```c
return (0.534 - 0.00211*Age - (0.236 - 0.00127*Age) * Tanh(0.4*log(B)/log(10) - 2.2)) * 10;
```
Pupil diameter, mm, from age (years) and background brightness `B` (nL) — Garstang [2000]
age-dependency model. `log(B)/log(10)` is a manual `log10(B)` (not `log10()` directly — replicate
the two-`log`-calls-and-divide form for FP fidelity, do not substitute `f64::log10`).

### `OpticFactor(double Bback, double kX, double *dobs, double JDNDaysUT, char *ObjectName, int TypeFactor, int helflag)` — swehel.c:224–301
Computes one of two composite correction factors (intensity-factor when `TypeFactor==0`,
background-factor when `TypeFactor==1`) folding in optics, atmosphere, color, and pupil-size
effects. Inputs unpacked from `dobs`: `Age=dobs[0]`, `SN=dobs[1]` (clamped to `≥1e-8` internally
as `SNi`), `Binocular=dobs[2]`, `OpticMag=dobs[3]`, `OpticDia=dobs[4]`, `OpticTrans=dobs[5]`.
`JDNDaysUT` is accepted but **explicitly unused** (`JDNDaysUT += 0.0; /* currently not used,
statement prevents compiler warning */`) — drop this parameter in the Rust port signature.

1. `Pst = PupilDia(23, Bback)` — pupil diameter at a **fixed reference age of 23** (Garstang's
   standard), used as the baseline against which the actual-age pupil (`PupilDia(Age, Bback)`,
   used later in `Fp`) is compared.
2. `if (OpticMag == 1) { OpticTrans = 1; OpticDia = Pst; }` — "using eye" mode: override
   caller-supplied transmission/diameter with ideal-eye values. **Dead code**: an `#if 0`-wrapped
   block immediately below (lines 243–250) duplicates this same logic for `OpticMag == 0`
   ("undefined") with an added comment "is done in default_heliacal_parameters()" — confirming
   that path was moved to `default_heliacal_parameters` (§7) and this local copy was disabled, not
   deleted. Do not port the `#if 0` block.
3. `CIb = 0.7` (background color index, "from Ben Sugerman"), `CIi = 0.5` (object color index for
   white light, "should be function of ObjectName" per comment but is **hardcoded constant**
   regardless of `ObjectName` — the `if (strcmp(ObjectName,"moon")==0) { ; }` branch is a literal
   no-op, `ObjectSize` stays `0` always). `ObjectSize = 0` unconditionally.
4. `Fb = (Binocular == 0) ? 1.41 : 1` (binocular gain factor).
5. Scotopic/photopic branch — **third distinct threshold value**: `if (Bback < 1645)
   is_scotopic = TRUE` (not `1394` from `CVA`, not `BNIGHT=1479` — see the consolidation warning
   under `CVA` above), again override-able by `SE_HELFLAG_VISLIM_PHOTOPIC`/`SCOTOPIC`.
   - Scotopic: `Fe = pow(10, 0.48*kX)`; `Fsc = mymin(1, (1 - (Pst/124.4)^4) / (1 -
     (OpticDia/OpticMag/124.4)^4))`; `Fci = pow(10, -0.4*(1 - CIi/2.0))`; `Fcb = pow(10,
     -0.4*(1 - CIb/2.0))`.
   - Photopic: `Fe = pow(10, 0.4*kX)`; `Fsc = mymin(1, (OpticDia/OpticMag/Pst)^2 * (1 -
     exp(-(Pst/6.2)^2)) / (1 - exp(-(OpticDia/OpticMag/6.2)^2)))`; `Fci = 1`; `Fcb = 1`.
6. `Ft = 1/OpticTrans`; `Fp = mymax(1, (Pst/(OpticMag*PupilDia(Age,Bback)))^2)`; `Fa =
   (Pst/OpticDia)^2`; `Fr = (1 + 0.03*(OpticMag*ObjectSize/CVA(Bback,SNi,helflag))^2) / SNi^2`
   (note `ObjectSize` is always `0`, so the `0.03*(...)^2` term is always `0`, making `Fr =
   1/SNi^2` in practice — but replicate the full formula structure for fidelity/traceability, in
   case `ObjectSize` is ever wired up); `Fm = OpticMag^2`.
7. Debug-only `#if SWEHEL_DEBUG` block (dead under normal build, `SWEHEL_DEBUG=0`) — skip.
8. Return: `TypeFactor==0` → `Fb*Fe*Ft*Fp*Fa*Fr*Fsc*Fci` (intensity factor, includes `Fr`/`Fci`,
   excludes `Fm`/`Fcb`); else → `Fb*Ft*Fp*Fa*Fm*Fsc*Fcb` (background factor, includes `Fm`/`Fcb`,
   excludes `Fr`/`Fci`). **Multiplication order matters for FP bit-fidelity** — replicate the
   exact left-to-right order shown (C evaluates `*` left-to-right for equal-precedence operators).

No global/static state in `OpticFactor` itself — pure function of its arguments (it calls
`PupilDia` and `CVA`, both also pure).

---

## §4 Atmospheric extinction

Four exponential-atmosphere-layer extinction coefficients, combined via `kt`/`Deltam` into a
total magnitude-of-extinction `Δm` at a given true altitude and Sun altitude/lat/eye-height.

### `kW(double HeightEye, double TempS, double RH)` — swehel.c:801–807
```c
double WT = 0.031;
WT *= 0.94 * (RH/100.0) * exp(TempS/15) * exp(-1*HeightEye/scaleHwater);
return WT;
```
Water-vapor extinction coefficient (Schaefer, Archaeoastronomy XV, 2000, p.128). `scaleHwater =
3000.0` m. Pure function, no state.

### `kOZ(double AltS, double sunra, double Lat)` — swehel.c:815–841
```c
static TLS double koz_last, alts_last, sunra_last;
if (AltS == alts_last && sunra == sunra_last) return koz_last;
alts_last = AltS; sunra_last = sunra;
OZ = 0.031;
LT = Lat * DEGTORAD;
kOZret = OZ * (3.0 + 0.4*(LT*cos(sunra*DEGTORAD) - cos(3*LT))) / 3.0;
altslim = -AltS - 12; if (altslim < 0) altslim = 0;
CHANGEKO = (100 - 11.6 * mymin(6, altslim)) / 100;
koz_last = kOZret * CHANGEKO;
return koz_last;
```
**Memoization caveat**: cache key is `(AltS, sunra)` only — **`Lat` is not part of the cache key**,
even though `kOZret` depends on `Lat`. If `Lat` changes between calls while `AltS`/`sunra` happen
to repeat, the C code returns a **stale value computed with the previous `Lat`** — this is a
genuine (if obscure/unlikely-to-matter-in-practice, since `Lat` is normally constant across a
whole computation) state-dependent bug/quirk in the original C. **Pure memoization would require
`Lat` in the key**; since the Rust port is stateless, simply recompute every time (equivalent to
"always miss the cache", which produces the *mathematically correct* result — the Rust output can
legitimately differ from a pathological C call sequence that varies `Lat` while `AltS`/`sunra`
repeat, but note this in case golden-test parity requires deliberately replicating the C caching
bug for a specific test sequence). Under `if ((0))` at line 833: dead debug-print block, skip. Note
`CHANGEKO`'s "day/night vision" transition uses altitude of the **Sun** (`AltS`) depressed below
`-12°` (astronomical-twilight-adjacent threshold), scaling ozone absorption from 100% down to a
floor of `100 - 11.6*6 = 30.4%` at `AltS ≤ -18°`.

### `kR(double AltS, double HeightEye)` — swehel.c:848–862
```c
double val = -AltS - 12; if (val < 0) val = 0; if (val > 6) val = 6;
CHANGEK = (1 - 0.166667 * val);
LAMBDA = 0.55 + (CHANGEK - 1) * 0.04;
return 0.1066 * exp(-1*HeightEye/scaleHrayleigh) * pow(LAMBDA/0.55, -4);
```
Rayleigh extinction (Schaefer p.128). `scaleHrayleigh = 8515.0` m. Same day/night `AltS`-based
`[0,6]` clamp pattern as `kOZ`'s `CHANGEKO`, but expressed with explicit `if` clamps here instead
of `mymin` (functionally identical to `mymin(6, mymax(-AltS-12, 0))`, just written out — note the
commented-out alternate form directly above in the source, `/*CHANGEK = (1 - 0.166667 *
Min(6, Max(-AltS - 12, 0)));*/`, confirming equivalence). No caching/state here (unlike `kOZ`/`ka`).
`0.166667` ≈ `1/6` (literal, not `1.0/6.0` — replicate the literal decimal for FP fidelity, it is
*not* bit-identical to `1.0/6.0`).

### `Sgn(double x)` — swehel.c:864–869
`x < 0 ? -1 : 1` — **note: `Sgn(0) = 1`, not `0`** (no zero case; only two return values, `int`).

### `ka(double AltS, double sunra, double Lat, double HeightEye, double TempS, double RH, double VR, char *serr)` — swehel.c:881–927
Aerosol extinction coefficient. **Memoized** on `(AltS, sunra)` only (same non-`Lat`-keyed caveat
as `kOZ` — `Lat`, `HeightEye`, `TempS`, `RH`, `VR` are all NOT part of the cache key despite the
formula depending on several of them; same "recompute always" resolution recommended for the
Rust port, see `kOZ` note above).
```c
SL = Sgn(Lat);
CHANGEKA = (1 - 0.166667 * mymin(6, mymax(-AltS - 12, 0)));
LAMBDA = 0.55 + (CHANGEKA - 1) * 0.04;
if (VR != 0) {
  if (VR >= 1) {
    // "prevailing visibility" / meteorological range interpretation, VR in km
    BetaVr = 3.912 / VR;
    Betaa = BetaVr - (kW(HeightEye,TempS,RH)/scaleHwater + kR(AltS,HeightEye)/scaleHrayleigh) * 1000 * astr2tau;
    kaact = Betaa * scaleHaerosol / 1000 * tau2astr;
    if (kaact < 0) { /* serr = warning: "Meteorological range is too long..." */ }
  } else {
    // VR interpreted directly as a (0,1) broadband extinction coefficient ktot
    kaact = VR - kW(HeightEye,TempS,RH) - kR(AltS,HeightEye) - kOZ(AltS,sunra,Lat);
    if (kaact < 0) { /* serr = warning: "atmospheric coefficient (ktot) is too low..." */ }
  }
} else {
  // From Schaefer, Archaeoastronomy XV, 2000, page 128 — humidity-based aerosol model
#ifdef SIMULATE_VICTORVB    // ALWAYS ACTIVE — SIMULATE_VICTORVB is always defined, see §11
  if (RH <= 0.00000001) RH = 0.00000001;
  if (RH >= 99.99999999) RH = 99.99999999;
#endif
  kaact = 0.1 * exp(-1*HeightEye/scaleHaerosol) * pow(1 - 0.32/log(RH/100.0), 1.33)
          * (1 + 0.33*SL*sin(sunra*DEGTORAD));
  kaact = kaact * pow(LAMBDA/0.55, -1.3);
}
ka_last = kaact;
return kaact;
```
The `VR` (datm[3]) dual-use dispatch: `VR==0` → humidity-based Schaefer model (default, when
caller doesn't supply an extinction/visibility override); `0 < VR < 1` → `VR` is a direct
broadband `ktot` coefficient, decompose by subtracting the other three components; `VR >= 1` →
`VR` is a Meteorological Range in km, converted via the Koschmieder relation (`3.912/VR` — note
this is the **standard** Koschmieder constant `ln(50)≈3.912` for "meteorological range" defined at
2% contrast threshold; comment cites `MOR=2.995/ke` from an ICAO doc as an *alternative*
definition not used here — the `2.995` value does not appear in the code, only in the comment as
context). **The `SIMULATE_VICTORVB` RH-clamp block is always compiled and always executes** in the
`VR==0` branch (§11) — clamping `RH` into `(1e-8, 99.999999998)` to avoid `log(0)`/`log(negative-
after-100-subtraction)` singularities in `pow(1 - 0.32/log(RH/100.0), 1.33)`. Warnings
(`kaact < 0`) are advisory only (`serr` set but **the function still returns the negative value**
— callers like `kt` clamp `kaact` to `≥0` afterward, see below).

### `kt(double AltS, double sunra, double Lat, double HeightEye, double TempS, double RH, double VR, int32 ExtType, char *serr)` — swehel.c:940–957
Dispatcher selecting which subset of `{kR, kW, kOZ, ka}` to sum, by `ExtType`:
`0`→`ka` only, `1`→`kW` only, `2`→`kR` only, `3`→`kOZ` only, `4`→all four (`ktot`). After
selection, **`if (kaact < 0) kaact = 0`** — clamps only the aerosol term to non-negative (matches
the "warning, not hard error" semantics of `ka`'s internal negative-value paths) before summing.
Return `kWact + kRact + kOZact + kaact` (all four terms present in the sum regardless of
`ExtType`, but non-selected terms are left at their `= 0` initializers — so effectively only the
selected subset contributes).

### `Airmass(double AppAltO, double Press)` — swehel.c:964–972
```c
zend = (90 - AppAltO) * DEGTORAD; if (zend > PI/2) zend = PI/2;
airm = 1 / (cos(zend) + 0.025 * exp(-11*cos(zend)));
return Press / 1013 * airm;
```
Kasten-style airmass approximation with pressure scaling (`1013` mbar reference, note: not
`1013.25` — a different rounded literal than `PressRef`, which is itself dead/unused, §1). Zenith
distance clamped to `≤90°` (`π/2` rad) before the `cos`. **Reachability**: only called from
`Deltam`'s `else` branch (`staticAirmass != 0`), which is unreachable since `staticAirmass` is the
compile-time constant `0` (§1) — i.e. under the shipped build, `Airmass` is *never actually
invoked*. Still port the function (it's a real, correct implementation, just dead under current
constants — porting it costs nothing and preserves parity if `staticAirmass` is ever flipped), but
do not expect any golden-test coverage to exercise it via `Deltam`.

### `Xext(double scaleH, double zend, double Press)` — swehel.c:980–983
```c
return Press / 1013.0 / (cos(zend) + 0.01 * sqrt(scaleH/1000.0) * exp(-30.0/sqrt(scaleH/1000.0) * cos(zend)));
```
Optical path length (airmass-like) through an exponential atmosphere layer of scale height
`scaleH` (meters), at zenith distance `zend` (radians) and station pressure `Press` (mbar).
`sqrt(scaleH/1000.0)` appears twice — replicate as two separate `sqrt` calls (matches C source
literally; do not hoist to a shared local for FP-fidelity purposes, though the value is
identical either way since `sqrt` is deterministic — the concern here is purely
readability/traceability parity with the C, not rounding, since caching a pure function's result
doesn't change its value).

### `Xlay(double scaleH, double zend, double Press)` — swehel.c:991–996
```c
double a = sin(zend) / (1.0 + (scaleH / Ra));   // Ra = 6378136.6 m, WGS84 equatorial radius
return Press / 1013.0 / sqrt(1.0 - a * a);
```
Optical path through a thin/high layer (used for ozone, `scaleHozone=20000`) accounting for
Earth's curvature via `Ra`. A commented-out equivalent one-liner directly above (`/*return Press /
1013.0 / sqrt(1.0 - pow(sin(zend) / (1.0 + (scaleH / Ra)), 2));*/`) confirms `a*a` is exactly
`pow(...,2)` refactored for perf — same value, replicate the live (uncommented) `a*a` form.

### `Deltam(double AltO, double AltS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)` — swehel.c:1033–1059
Total extinction magnitude Δm for an object at true altitude `AltO`, given Sun altitude `AltS`.
```c
PresE = PresEfromPresS(datm[1], datm[0], HeightEye);
TempE = TempEfromTempS(datm[1], HeightEye, LapseSA);
AppAltO = AppAltfromTopoAlt(AltO, TempE, PresE, helflag);
static TLS double alts_last, alto_last, sunra_last, deltam_last;
if (AltS==alts_last && AltO==alto_last && sunra==sunra_last) return deltam_last;
alts_last=AltS; alto_last=AltO; sunra_last=sunra;
if (staticAirmass == 0) {                              // ALWAYS true (compile-time constant)
  zend = (90 - AppAltO) * DEGTORAD; if (zend > PI/2) zend = PI/2;
  xR  = Xext(scaleHrayleigh, zend, datm[0]);
  XW  = Xext(scaleHwater,    zend, datm[0]);
  Xa  = Xext(scaleHaerosol,  zend, datm[0]);
  XOZ = Xlay(scaleHozone,    zend, datm[0]);
  deltam = kR(AltS,HeightEye)*xR
         + kt(AltS,sunra,Lat,HeightEye,datm[1],datm[2],datm[3],0,serr)*Xa
         + kOZ(AltS,sunra,Lat)*XOZ
         + kW(HeightEye,datm[1],datm[2])*XW;
} else {                                                // unreachable, see Airmass note above
  deltam = kt(AltS,sunra,Lat,HeightEye,datm[1],datm[2],datm[3],4,serr) * Airmass(AppAltO, datm[0]);
}
deltam_last = deltam;
return deltam;
```
**Memoization caveat, again**: cache key `(AltS, AltO, sunra)` excludes `Lat`, `HeightEye`,
`datm[*]`, `helflag` — same non-total-key issue as `kOZ`/`ka`. Recommended Rust-port resolution:
recompute always (drop the cache; mathematically correct, matches C only when the excluded
parameters are in fact held constant across the calling sequence, which is the normal/intended
usage pattern in this module — see §11 for the general policy on these caches).
`kt(...,0,...)` in the live branch selects `ka` only (`ExtType=0`, aerosol) for the `Xa` term's
coefficient — i.e. `kt` here is used purely as "give me `ka` with clamping," not as a multi-term
sum (since `ExtType=0` zeros the other three internally before summing, per `kt`'s definition
above).

---

## §6 Sky brightness model

All five functions below return **nanoLamberts (nL)**; internally most work in `erg`-equivalent
units (Schaefer's original units) and convert via `erg2nL = 1/nL2erg = 1/1.02e-15` at the very end
(`mymax(result, 0) * erg2nL`).

### `Bn(double AltO, double JDNDayUT, double AltS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)` — swehel.c:1073–1096
Night-sky brightness (zodiacal light + starlight), modulated by an 11.1-year solar-activity cycle.
```c
PresE = PresEfromPresS(datm[1], datm[0], HeightEye);
TempE = TempEfromTempS(datm[1], HeightEye, LapseSA);
AppAltO = AppAltfromTopoAlt(AltO, TempE, PresE, helflag);
if (AppAltO < 10) AppAltO = 10;                        // floor: Bn constant below 10° altitude
zend = (90 - AppAltO) * DEGTORAD;
swe_revjul(JDNDayUT, SE_GREG_CAL, &iyar, &imon, &iday, &dut);
YearB = iyar; MonthB = imon; DayB = iday;
B0 = 0.0000000000001;                                  // 1e-13
Bna = B0 * (1 + 0.3 * cos(6.283 * (YearB + ((DayB-1)/30.4 + MonthB-1)/12 - 1990.33) / 11.1));
kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
Bnb = Bna * (0.4 + 0.6/sqrt(1 - 0.96*pow(sin(zend),2))) * pow(10, -0.4*kX);
return mymax(Bnb, 0) * erg2nL;
```
`6.283` ≈ `2π` (literal, **not** `2*PI` or `2*M_PI` — replicate the rounded literal exactly for FP
fidelity). `1990.33` is the reference epoch for the solar-cycle phase. `(DayB-1)/30.4` (integer
`DayB`/`MonthB` promoted to double via the `double YearB,MonthB,DayB` locals) approximates
day-of-year fraction using a fixed 30.4-day month, matching the same `30.4` divisor used in
`SunRA`'s crude formula (§7) and `PupilDia`-adjacent code elsewhere — a recurring "average month
length" constant in this module, always `30.4`, never `30.44` or a true value. `swe_revjul` is the
**only ephemeris/calendar call** in `Bn` (used purely to extract Y/M/D from the UT Julian day for
the solar-cycle phase — no position/ephemeris data needed).

### `MoonsBrightness(double dist, double phasemoon)` — swehel.c:1159–1164
```c
double log10 = 2.302585092994;                          // this is actually ln(10), used as a divisor — see below
return -21.62 + 5*log(dist/(Ra/1000)) / log10 + 0.026*fabs(phasemoon) + 0.000000004*pow(phasemoon,4);
```
**Naming trap in the C source**: the local variable is named `log10` but holds `ln(10) =
2.302585092994` (used as `log(x)/log10` = manual `log10(x)`, same pattern as `PupilDia`) — it is
NOT `f64::log10` or `libm log10()`. Replicate the manual `ln(x)/ln(10)` computation; do not
substitute a direct `log10` call (same value mathematically, but keep the pattern explicit per
this module's general style and to match any golden-test FP tolerance expectations). `Ra/1000` =
Earth radius in km (`6378.1366`), used as a reference distance for the `5*log10(dist/Ra_km)`
distance-brightness term (standard `-5 log10(d)` inverse-square-ish scaling for a reflected-light
body). `0.000000004` = `4e-9`; `phasemoon` in degrees (magnitude phase-angle correction terms).

### `MoonPhase(double AltM, double AziM, double AltS, double AziS)` — swehel.c:1172–1181
```c
double MoonAvgPar = 0.95;                               // average Moon parallax, degrees
return 180 - acos(cos(AziSi-AziMi-MoonAvgPar*DEGTORAD) * cos(AltMi+MoonAvgPar*DEGTORAD) * cos(AltSi)
                  + sin(AltSi)*sin(AltMi+MoonAvgPar*DEGTORAD)) / DEGTORAD;
```
A commented-out simpler variant immediately above (without the `- MoonAvgPar*DEGTORAD` term inside
the first `cos`) shows this `AziSi-AziMi-MoonAvgPar*DEGTORAD` adjustment was a deliberate
correction added later — **use the live (uncommented) formula**, not the commented one. Computes
the Sun-Moon phase angle (degrees) from topocentric alt/az of both bodies, folding in a fixed
average parallax correction (`0.95°`) rather than a per-epoch computed lunar parallax.

### `Bm(double AltO, double AziO, double AltM, double AziM, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)` — swehel.c:1186–1213
Moonlight contribution to sky brightness at object position `(AltO,AziO)`.
```c
double M0 = -11.05;
double lunar_radius = 0.25 * DEGTORAD;
object_is_moon = (AltO==AltM && AziO==AziM);
Bm = 0;
if (AltM > -0.26 && !object_is_moon) {                  // Moon must be above horizon (with a small margin) and object != Moon
  RM = DistanceAngle(AltO*DEGTORAD, AziO*DEGTORAD, AltM*DEGTORAD, AziM*DEGTORAD) / DEGTORAD;
  if (RM <= lunar_radius) RM = lunar_radius;             // avoid singularity when object is inside/at the Moon's disc
  kXM = Deltam(AltM, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
  kX  = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
  C3 = pow(10, -0.4*kXM);
  FM = 62000000.0/RM/RM + pow(10, 6.15 - RM/40) + pow(10,5.36)*(1.06 + pow(cos(RM*DEGTORAD),2));
  Bm = FM*C3 + 440000*(1 - C3);
  phasemoon = MoonPhase(AltM,AziM,AltS,AziS);
  MM = MoonsBrightness(MoonDistance, phasemoon);         // MoonDistance = 384410.4978 km, fixed (not epoch distance)
  Bm = Bm * pow(10, -0.4*(MM - M0 + 43.27));
  Bm = Bm * (1 - pow(10, -0.4*kX));
}
Bm = mymax(Bm,0) * erg2nL;
return Bm;
```
`AltM > -0.26` threshold (not `0`) — allows the Moon to contribute light slightly below the
geometric horizon (accounts for the Moon's own angular radius/refraction near the horizon,
roughly). Note `MoonDistance` is the file-header **constant** `384410.4978` km — the Moon's actual
current geocentric distance (which would vary with epoch/`swe_calc`) is never queried here; this
is a fixed nominal Earth-Moon distance used only for the `MoonsBrightness` magnitude-vs-distance
term, not a live ephemeris value. `62000000.0` (`FM`'s first term numerator) is the **same
constant** used in `Bday`'s `FS` (§ below) — the Schaefer daylight/moonlight scattering formula is
structurally identical between `Bm` and `Bday`, differing only in which body's magnitude feeds the
final `10^(-0.4*(M-M0+43.27))` scaling term (`MM` here vs. the Sun's `MS=-26.74` in `Bday`).

### `Btwi(double AltO, double AziO, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)` — swehel.c:1218–1234
```c
M0 = -11.05; MS = -26.74;
PresE = PresEfromPresS(datm[1],datm[0],HeightEye);
TempE = TempEfromTempS(datm[1],HeightEye,LapseSA);
AppAltO = AppAltfromTopoAlt(AltO, TempE, PresE, helflag);
ZendO = 90 - AppAltO;
RS = DistanceAngle(AltO*DEGTORAD,AziO*DEGTORAD,AltS*DEGTORAD,AziS*DEGTORAD) / DEGTORAD;
kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
k  = kt(AltS, sunra, Lat, HeightEye, datm[1], datm[2], datm[3], 4, serr);   // ExtType=4: full ktot
Btwi = pow(10, -0.4*(MS - M0 + 32.5 - AltS - (ZendO/(360*k))));
Btwi = Btwi * (100/RS) * (1 - pow(10,-0.4*kX));
Btwi = mymax(Btwi,0) * erg2nL;
return Btwi;
```
Twilight sky brightness (Schaefer p.129). Uses the **total** extinction coefficient `k`
(`ExtType=4`) unlike `Bm`/`Bday` which use `Deltam` (the altitude-integrated version); this `k` is
the bare per-airmass coefficient sum, used directly in the `ZendO/(360*k)` twilight-decay term
(zenith distance in **degrees** here, not radians — `ZendO = 90 - AppAltO` with no `DEGTORAD`
conversion, since it's divided by the dimensionless `360*k`, an empirical twilight-glow decay
constant, not a trig argument).

### `Bday(double AltO, double AziO, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)` — swehel.c:1246–1261
```c
M0 = -11.05; MS = -26.74;
RS = DistanceAngle(AltO*DEGTORAD,AziO*DEGTORAD,AltS*DEGTORAD,AziS*DEGTORAD) / DEGTORAD;
kXS = Deltam(AltS, AltS, sunra, Lat, HeightEye, datm, helflag, serr);   // NOTE: object-altitude arg = AltS (Sun's own altitude), not AltO
kX  = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
C4 = pow(10, -0.4*kXS);
FS = 62000000.0/RS/RS + pow(10, 6.15 - RS/40) + pow(10,5.36)*(1.06 + pow(cos(RS*DEGTORAD),2));
Bday = FS*C4 + 440000.0*(1 - C4);
Bday = Bday * pow(10, -0.4*(MS - M0 + 43.27));
Bday = Bday * (1 - pow(10, -0.4*kX));
Bday = mymax(Bday,0) * erg2nL;
return Bday;
```
**Important asymmetry vs. `Bm`**: `kXS = Deltam(AltS, AltS, ...)` passes the Sun's own altitude as
*both* the "object altitude" and "sun altitude" arguments to `Deltam` — i.e. it computes the
extinction *along the Sun's own line of sight* (how much the Sun itself is extinguished/scattered
overhead), as opposed to `Bm`'s `kXM = Deltam(AltM, AltS, ...)` which correctly separates Moon
altitude from Sun altitude as two different arguments. This is intentional (per the accompanying
old-BASIC-source comment block at lines 1236–1245, `2310 C4=10.0^(-.4*K(I)*XS)` — `K(I)` there is
the Sun's own extinction), not a copy-paste bug — replicate exactly (`Deltam(AltS, AltS, ...)`, not
`Deltam(AltO, AltS, ...)`, for the `kXS` term specifically).

### `Bcity(double Value, double Press)` — swehel.c:1268–1274
```c
Press += 0.0;  /* unused; statement prevents compiler warning */
return mymax(Value, 0);
```
Trivial passthrough/clamp — light-pollution brightness is currently always supplied as a literal
`0` by every call site in this line range (`Bsky` calls `Bcity(0, datm[0])`, §below) — i.e. light
pollution is **not modeled** at all in this range of the file; `Bcity` exists as a hook for a
future/external light-pollution value but is fed a constant zero here. `Press` parameter is
unused (drop from Rust signature, or keep for API-shape parity with a comment).

### `Bsky(double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)` — swehel.c:1279–1307
```c
Bsky = 0;
if (AltS < -3) {
  Bsky += Btwi(...);
} else if (AltS > 4) {
  Bsky += Bday(...);
} else {
  Bsky += mymin(Bday(...), Btwi(...));                  // transition zone: take the smaller of the two models
}
if (Bsky < 200000000.0)                                 // 2e8 nL: "if max. Bm [1E7] <5% of Bsky don't add Bm"
  Bsky += Bm(...);
if (AltS <= 0)
  Bsky += Bcity(0, datm[0]);                            // always adds 0 (Bcity clamp of literal 0), see above
if (Bsky < 5000)                                         // "if max. Bn [250] <5% of Bsky don't add Bn"
  Bsky = Bsky + Bn(...);
return Bsky;
```
Note the C comment thresholds (`200000000.0`/`5000`) reference "5% of Bm's max [1E7]" and "5% of
Bn's max [250]" respectively — the actual coded thresholds (`2e8`, `5000`) are **20×** the cited
"max" values (`1e7 * 20 = 2e8`; `250 * 20 = 5000`), i.e. the optimization skips adding the smaller
brightness source once the dominant sky brightness is already ≥20× that source's own maximum
possible contribution (so the smaller term could contribute at most ~5%, hence the comment) —
replicate the literal thresholds `200000000.0` and `5000`, not a recomputed "20×" formula (do not
derive these from `1e7`/`250` symbolically; they are independent literals in the source, chosen
to *not* need updating symbolically if the cited values were revised elsewhere). Dead `if ((0))`
debug-print block at line 1289 — skip. Both `Bday`/`Btwi` are called **twice each** in the
transition-zone (`else`) branch (`AltS` between -3° and 4°) — once inside the `mymin(...)` call's
two arguments — no result caching between the two calls within `Bsky` itself, though `Deltam`'s
internal memoization (§4) will short-circuit the *actual extinction* sub-computation on the second
call if `(AltO,AltS,sunra)` are unchanged (as they are here, called twice with identical args).

---

## §7 Object location & magnitude

### `SunRA(double JDNDaysUT, int32 helflag, char *serr)` — swehel.c:553–583
```c
static TLS double tjdlast, ralast;
if (JDNDaysUT == tjdlast) return ralast;                 // memoized on tjd only (helflag not part of key!)
#ifndef SIMULATE_VICTORVB                                 // DEAD — SIMULATE_VICTORVB always defined, see §11
  ... swe_calc(SE_SUN, iflag=epheflag|SEFLG_EQUATORIAL|SEFLG_NONUT|SEFLG_TRUEPOS) at tjd_tt ...
  ralast = x[0]; tjdlast = JDNDaysUT; return ralast;
#endif
swe_revjul(JDNDaysUT, SE_GREG_CAL, &iyar, &imon, &iday, &dut);
tjdlast = JDNDaysUT;
ralast = swe_degnorm((imon + (iday-1)/30.4 - 3.69) * 30);
return ralast;
```
**§11-flagged finding**: because `SIMULATE_VICTORVB` is always `#define`d (swephexp.h:451), the
`#ifndef SIMULATE_VICTORVB` block (the "real" ephemeris-based Sun-RA computation via `swe_calc`)
**never compiles** — `SunRA` **always** uses the crude calendar-only approximation:
`swe_degnorm((month + (day-1)/30.4 - 3.69) * 30)` — a 12-step, 30°-per-month linear model of the
Sun's right ascension based purely on calendar month/day (not even the year), completely ignoring
`helflag`'s `SE_HELFLAG_HIGH_PRECISION` bit (the `if (1) { ... }` guard immediately inside the dead
block, `helflag & SE_HELFLAG_HIGH_PRECISION` commented out in favor of a literal `if(1)`, is itself
moot since the whole block is unreachable). **A Rust port must replicate the crude formula as the
actual behavior**, not the "intended" ephemeris-based one — porting the `swe_calc`-based branch
instead would silently diverge from the shipped C library's real output. `swe_revjul` is used
"because it seems much faster than calling swe_revjul()" per the comment (sic — the comment
appears to be leftover from an earlier refactor and no longer makes literal sense; the real
rationale is presumably that this crude method avoids a full ephemeris call for a value used
"1000s of times" per the comment). Memoization key is `JDNDaysUT` only — `helflag` is deliberately
irrelevant here since the live path never consults it anyway.

### `ObjectLoc(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 Angle, int32 helflag, double *dret, char *serr)` — swehel.c:683–726
Generic accessor returning one of 7 angle types for a named object:
`Angle`: `0`=TopoAlt, `1`=Azi, `2`=Topo Declination, `3`=Topo Rectascension, `4`=AppAlt,
`5`=Geo Declination, `6`=Geo Rectascension, `7`=(alias for `0`, remapped at entry: `if (Angle==7)
Angle=0`).
```c
iflag = SEFLG_EQUATORIAL | (helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH));
if (!(helflag & SE_HELFLAG_HIGH_PRECISION)) iflag |= SEFLG_NONUT|SEFLG_TRUEPOS;
if (Angle < 5) iflag |= SEFLG_TOPOCTR;                   // topocentric for Angle 0-4, geocentric for 5/6
tjd_tt = JDNDaysUT + swe_deltat_ex(JDNDaysUT, epheflag, serr);
Planet = DeterObject(ObjectName);
x = (Planet != -1) ? swe_calc(tjd_tt, Planet, iflag, x, serr) : call_swe_fixstar(ObjectName, tjd_tt, iflag, x, serr);
if (Angle==2 || Angle==5) *dret = x[1];                  // declination directly from equatorial x[]
else if (Angle==3 || Angle==6) *dret = x[0];             // RA directly from equatorial x[]
else {
  xin[0]=x[0]; xin[1]=x[1];
  swe_azalt(JDNDaysUT, SE_EQU2HOR, dgeo, datm[0], datm[1], xin, xaz);
  if (Angle==0) *dret = xaz[1];                          // true altitude
  if (Angle==4) *dret = AppAltfromTopoAlt(xaz[1], datm[0], datm[1], helflag);  // NOTE arg order below
  if (Angle==1) { xaz[0]+=180; if (xaz[0]>=360) xaz[0]-=360; *dret = xaz[0]; } // azimuth, flipped 180°
}
return OK;
```
**Argument-order gotcha**: the call `AppAltfromTopoAlt(xaz[1], datm[0], datm[1], helflag)` passes
`datm[0]` (pressure) into the function's `TempE` parameter slot and `datm[1]` (temperature) into
`PresE` — **this looks like a swapped-argument bug relative to `AppAltfromTopoAlt`'s own signature
`(TopoAlt, TempE, PresE, helflag)`** (compare `Deltam`'s correct usage:
`AppAltfromTopoAlt(AltO, TempE, PresE, helflag)` computed from proper locals). **Verify against the
live C behavior, not the "intended" signature, when porting `ObjectLoc`'s `Angle==4` branch** —
replicate the swapped `datm[0]`/`datm[1]` order exactly as coded here (this is what
`swe_heliacal_pheno_ut`, past line 1714, actually receives when it asks for Angle-4/AppAlt via this
function) even though it differs from every other caller of `AppAltfromTopoAlt` in the file.
`Angle==1` (azimuth) flips the `swe_azalt` "from south" convention to "from north" by adding 180°
(wrapped to `[0,360)`) — a convention difference from `azalt_cart` below, which returns the raw
(un-flipped) `swe_azalt` azimuth. `iflag` mask: same narrow ephemeris-bit extraction pattern as
§2's `call_swe_rise_trans`.

### `azalt_cart(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 helflag, double *dret, char *serr)` — swehel.c:737–771
Returns 6 values in `dret[0..5]`: `[0]`=azimuth (raw `swe_azalt` convention, not flipped), `[1]`=
true altitude, `[2]`=apparent altitude, `[3..5]`=Cartesian unit vector of the **apparent**-altitude
direction (`xaz[1]=xaz[2]` (swap in apparent altitude), `xaz[2]=1` (unit radius), then
`swi_polcart(xaz,xaz)` converts az/alt/1 spherical → Cartesian in place). Always forces
`iflag |= SEFLG_TOPOCTR` unconditionally (unlike `ObjectLoc`, which conditions topocentric-ness on
`Angle<5`) — `azalt_cart` is always topocentric. Otherwise structurally identical setup to
`ObjectLoc` (`epheflag` mask, `SEFLG_NONUT|TRUEPOS` unless high-precision, `swe_deltat_ex` for
`tjd_tt`, `DeterObject`/`call_swe_fixstar` dispatch).

### `Magnitude(double JDNDaysUT, double *dgeo, char *ObjectName, int32 helflag, double *dmag, char *serr)` — swehel.c:1106–1127
```c
*dmag = -99.0;
Planet = DeterObject(ObjectName);
iflag = SEFLG_TOPOCTR | SEFLG_EQUATORIAL | (helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH));
if (!(helflag & SE_HELFLAG_HIGH_PRECISION)) iflag |= SEFLG_NONUT|SEFLG_TRUEPOS;
if (Planet != -1) {
  swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);               // STATEFUL — see §11
  if (swe_pheno_ut(JDNDaysUT, Planet, iflag, x, serr) == ERR) return ERR;
  *dmag = x[4];                                          // swe_pheno_ut's magnitude output slot
} else {
  if (call_swe_fixstar_mag(ObjectName, dmag, serr) == ERR) return ERR;
}
return OK;
```
**`swe_set_topo` call**: this is a genuine C global-state mutation (`swed.topd`) that
`swe_pheno_ut` reads back internally when `SEFLG_TOPOCTR` is set. **Stateless Rust port note**:
must thread `dgeo` explicitly into whatever phenomena/magnitude computation replaces
`swe_pheno_ut` here, rather than mutating shared state — same pattern as documented in
`docs/c-ref-riseset.md` §3.4's STATELESS PORT NOTE for `swe_set_topo`. `*dmag` default `-99.0`
(sentinel, overwritten on success) is never actually returned to the caller on the error paths
(function returns `ERR` before use) — effectively dead initialization outside of successful
completion, but harmless/idiomatic defensive-init; no need to replicate the literal `-99.0`
sentinel unless a caller elsewhere in the file inspects `*dmag` after an `ERR` return (grep confirms
no caller in this range does).

### `fast_magnitude` — swehel.c:1129–1152 — **dead code** (`#if 0`)
Memoizing wrapper around `Magnitude`, keyed on `(ipl, helflag, tjd within 5/1440 day)`, with the
same `ipl > SE_MOON → ipli = 2` planet-index-collapsing trick as the also-dead `call_swe_calc`
(§2). Never compiles; skip entirely — do not port.

---

## §8 `VisLimMagn` (1382) & `swe_vis_limit_mag` (1464)

### `VisLimMagn(double *dobs, double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, int32 *scotopic_flag, char *serr)` — swehel.c:1382–1443
```c
Bsk = Bsky(AltO,AziO,AltM,AziM,JDNDaysUT,AltS,AziS,sunra,Lat,HeightEye,datm,helflag,serr);
kX  = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
CorrFactor1 = OpticFactor(Bsk, kX, dobs, JDNDaysUT, "", 1, helflag);   // TypeFactor=1: background factor
CorrFactor2 = OpticFactor(Bsk, kX, dobs, JDNDaysUT, "", 0, helflag);   // TypeFactor=0: intensity factor
is_scotopic = (Bsk < 1645);                              // THIRD distinct threshold literal, see CVA note (§5)
if (helflag & SE_HELFLAG_VISLIM_PHOTOPIC) is_scotopic = FALSE;
if (helflag & SE_HELFLAG_VISLIM_SCOTOPIC) is_scotopic = TRUE;
if (is_scotopic) {
  C1 = 1.5848931924611e-10;   /* == pow(10,-9.8), precomputed literal */
  C2 = 0.012589254117942;     /* == pow(10,-1.9) */
  if (scotopic_flag) *scotopic_flag = 1;                 // SE_SCOTOPIC_FLAG = 1
} else {
  C1 = 4.4668359215096e-9;    /* == pow(10,-8.35) */
  C2 = 1.2589254117942e-6;    /* == pow(10,-5.9) */
  if (scotopic_flag) *scotopic_flag = 0;                 // SE_PHOTOPIC_FLAG = 0
}
if (scotopic_flag) {
  if (BNIGHT*BNIGHT_FACTOR > Bsk && BNIGHT/BNIGHT_FACTOR < Bsk)   // BNIGHT=1479, BNIGHT_FACTOR=1.0
    *scotopic_flag |= 2;                                  // SE_MIXEDOPIC_FLAG bit, OR-ed in
}
Bsk = Bsk * CorrFactor1;                                  // NOTE: multiply, despite a commented-out "/CorrFactor1" alternative directly above
Th = C1 * pow(1 + sqrt(C2 * Bsk), 2) * CorrFactor2;
return -16.57 - 2.5 * (log(Th) / log10);                  // log10 here = manual ln(10)=2.302585092994 local, NOT libm log10()
```
**`log10` naming trap again** (same pattern as `MoonsBrightness`, §6): the local `double log10 =
2.302585092994;` is `ln(10)`, used as a divisor for a manual `log10(Th)` — do not call
`f64::log10` directly, replicate `Th.ln() / 2.302585092994_f64.ln()` (or the literal
`2.302585092994`) form. `C1`/`C2` are given both as decimal literals **and** as `pow(10,...)`
comments confirming the precomputed values — replicate the **decimal literal**, not a fresh
`pow(10.0, -9.8)` call (matches C's compile-time-equivalent precision; C wrote out the literal
specifically, presumably to avoid a runtime `pow` call — bit-identical either way in practice for
these particular exponents, but replicate the literal for direct traceability). Dead `#if 0` block
at lines 1437–1441 (an alternate return incorporating `SN`/Snellen-ratio directly into the
magnitude formula) — skip, the live `return` statement (line 1442) is the one actually compiled.
`*scotopic_flag` semantics: base value `0`(photopic)/`1`(scotopic) per `SE_PHOTOPIC_FLAG`/
`SE_SCOTOPIC_FLAG`, with bit `2` (`SE_MIXEDOPIC_FLAG`) OR-ed in when `Bsk` falls strictly inside
`(BNIGHT/BNIGHT_FACTOR, BNIGHT*BNIGHT_FACTOR) = (1479, 1479)` — since `BNIGHT_FACTOR=1.0` this
open interval is **empty** (`1479 > Bsk` and `Bsk > 1479` can never both hold), so **the
`SE_MIXEDOPIC_FLAG` bit can never actually be set under the current `BNIGHT_FACTOR=1.0` constant**
— it is live code but permanently unreachable given the current constant value (would become
reachable if `BNIGHT_FACTOR` were set to e.g. `1.1`, widening the interval to `(1344.5,1626.9)`).
Port the logic faithfully (don't hardcode "flag 2 never sets") since it depends on a named,
independently-adjustable constant, not a structural dead branch.

### `swe_vis_limit_mag(double tjdut, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 helflag, double *dret, char *serr)` — swehel.c:1464–1541
```c
int32 retval = OK, scotopic_flag = 0;
for (i=0;i<7;i++) dret[i]=0;
tolower_string_star(ObjectName);                          // mutates caller's ObjectName buffer in place!
if (DeterObject(ObjectName) == SE_SUN) { serr = "..."; return ERR; }
swi_set_tid_acc(tjdut, helflag, 0, serr);                 // sets tidal-acceleration model, STATEFUL global
sunra = SunRA(tjdut, helflag, serr);
default_heliacal_parameters(datm, dgeo, dobs, helflag);   // mutates datm/dobs in place — fills defaults
swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);                  // STATEFUL, see §11
if (ObjectLoc(tjdut, dgeo, datm, ObjectName, 0, helflag, &AltO, serr) == ERR) return ERR;
if (AltO < 0) { serr="object is below local horizon"; *dret = -100; return -2; }
if (ObjectLoc(..., 1, ..., &AziO, ...) == ERR) return ERR;
if (helflag & SE_HELFLAG_VISLIM_DARK) { AltS=-90; AziS=0; }
else { AltS = ObjectLoc(...,"sun",0,...); AziS = ObjectLoc(...,"sun",1,...); }
if (starts-with "moon" || VISLIM_DARK || VISLIM_NOMOON) { AltM=-90; AziM=0; }
else { AltM = ObjectLoc(...,"moon",0,...); AziM = ObjectLoc(...,"moon",1,...); }
dret[0] = VisLimMagn(dobs, AltO,AziO,AltM,AziM, tjdut, AltS,AziS, sunra, dgeo[1], dgeo[2], datm, helflag, &scotopic_flag, serr);
dret[1]=AltO; dret[2]=AziO; dret[3]=AltS; dret[4]=AziS; dret[5]=AltM; dret[6]=AziM;
if (Magnitude(tjdut, dgeo, ObjectName, helflag, &(dret[7]), serr) == ERR) return ERR;
retval = scotopic_flag;
return retval;
```
**Return-value semantics** (from the function's own doc comment, lines 1454–1463):
`-1`=Error, `-2`=Object below horizon, `0`=OK photopic, `|1`=OK scotopic, `|2`=OK near
photopic/scotopic limit (i.e. return value is `scotopic_flag`, which itself is `0|1` OR-able with
`2`, per `VisLimMagn`'s semantics above — so possible successful returns are `0, 1, 2, 3`
(`3 = 1|2`, scotopic+mixed), never negative except the two explicit error/below-horizon sentinels).

**`dret[0..7]` layout**:
| Index | Meaning |
|---|---|
| `dret[0]` | limiting visual magnitude (from `VisLimMagn`) |
| `dret[1]` | object's true altitude, degrees |
| `dret[2]` | object's azimuth, degrees |
| `dret[3]` | Sun's altitude, degrees (`-90` if `VISLIM_DARK`) |
| `dret[4]` | Sun's azimuth, degrees (`0` if `VISLIM_DARK`) |
| `dret[5]` | Moon's altitude, degrees (`-90` if object is Moon, or `VISLIM_DARK`/`VISLIM_NOMOON`) |
| `dret[6]` | Moon's azimuth, degrees (`0` under the same conditions) |
| `dret[7]` | object's actual visual magnitude (from `Magnitude`) |

**Global state touched**: `swi_set_tid_acc` (tidal-acceleration/ΔT model selection — affects
subsequent `swe_deltat_ex` calls globally) and `swe_set_topo` (observer position cache) are both
mutated here before the rest of the computation reads them back indirectly via `ObjectLoc`/
`Magnitude`/`swe_pheno_ut`. **Stateless Rust port**: both must become explicit parameters/config
threaded through the equivalent call chain rather than global mutation — `swi_set_tid_acc`'s
effect should already be covered by whatever ΔT/tidal-acceleration config exists in
`EphemerisConfig` (check `docs/codebase-map.md` before adding a duplicate knob); `swe_set_topo`
follows the same pattern documented in `docs/c-ref-riseset.md`'s STATELESS PORT NOTE. Also note
`tolower_string_star(ObjectName)` **mutates the caller's string buffer in place** in C (in-place
lowercase up to the first comma) — in Rust this is just "lowercase the input `&str` into an owned
`String` before use," no aliasing concern.

---

## §9 `TopoArcVisionis` (1562) & `swe_topo_arcus_visionis` (1601)

### `TopoArcVisionis(double Magn, double *dobs, double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, double *dret, char *serr)` — swehel.c:1562–1599
Bisection search for the **arcus visionis** (Sun-depression-below-object angle, `Xl`/`xR` in
degrees) at which an object of known magnitude `Magn` becomes exactly visible.
```c
Xl = 45; xR = 0;                                          // search bracket: Sun 45° below object .. Sun at object's altitude
Yl = Magn - VisLimMagn(dobs, AltO,AziO,AltM,AziM, JDNDaysUT, AltO-Xl, AziS, sunra, Lat, HeightEye, datm, helflag, NULL, serr);
Yr = Magn - VisLimMagn(dobs, AltO,AziO,AltM,AziM, JDNDaysUT, AltO-xR, AziS, sunra, Lat, HeightEye, datm, helflag, NULL, serr);
if (Yl * Yr <= 0) {                                       // sign change present → root exists in bracket
  while (fabs(xR - Xl) > epsilon) {                       // epsilon = 0.001
    Xm = (xR + Xl) / 2.0;
    AltSi = AltO - Xm; AziSi = AziS;
    Ym = Magn - VisLimMagn(dobs, AltO,AziO,AltM,AziM, JDNDaysUT, AltSi, AziSi, sunra, Lat, HeightEye, datm, helflag, NULL, serr);
    if (Yl * Ym > 0) { Xl = Xm; Yl = Ym; }                 // same sign as left endpoint → discard left half
    else             { xR = Xm; Yr = Ym; }                 // discard right half
  }
  Xm = (xR + Xl) / 2.0;
} else {
  Xm = 99;                                                 // no sign change in [0,45]° bracket → sentinel "never visible in this range"
}
if (Xm < AltO) Xm = AltO;                                  // clamp: arc can't be less than the object's own altitude
*dret = Xm;
return OK;
```
Note the **Sun's position is parameterized only by altitude** here (`AltO - Xm` as the Sun's
altitude, `AziS` held fixed at the caller-supplied value throughout the whole bisection) — azimuth
is never varied during the search, only depression angle. `scotopic_flag` pointer is passed as
`NULL` throughout (this caller doesn't need the vision-mode classification, only the numeric
magnitude comparison). Sentinel `Xm = 99` (degrees) signals "no visibility crossing found in
`[0°,45°]` Sun-depression" — always `≥ AltO` after the final clamp, so downstream callers (e.g.
`HeliacalAngle`, §10) can treat any `Xm` this large as "not visible in this configuration" without
a separate error code. No global state, no ephemeris calls — pure numerical bisection over
`VisLimMagn` (which itself is pure given its already-resolved inputs).

### `swe_topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr)` — swehel.c:1601–1610
```c
sunra = SunRA(tjdut, helflag, serr);
if (serr != NULL && *serr != '\0') return ERR;             // NOTE: SunRA never actually sets *serr in the live code path, see §7 — this check is effectively dead
default_heliacal_parameters(datm, dgeo, dobs, helflag);
return TopoArcVisionis(mag, dobs, alt_obj, azi_obj, alt_moon, azi_moon, tjdut, azi_sun, sunra, dgeo[1], dgeo[2], datm, helflag, dret, serr);
```
Thin wrapper: resolves `sunra` and default parameters, then delegates. Note the caller supplies
`alt_obj`/`azi_obj`/`azi_sun`/`azi_moon`/`alt_moon` directly (no ephemeris lookup inside this
function at all — geometry is entirely caller-supplied, unlike `swe_vis_limit_mag` which computes
object/Sun/Moon positions itself via `ObjectLoc`). `swi_set_tid_acc` is **not** called here (unlike
`swe_vis_limit_mag`/`swe_heliacal_angle`) — this function relies on whatever tidal-acceleration
model is already globally active from a prior call; **stateless Rust port**: since this function
takes no `tjd`-dependent ephemeris positions as input beyond `sunra` (itself using the crude,
ΔT-independent calendar formula per §7), the tidal-acceleration config barely matters here anyway
— thread it through for API consistency but it has no numerical effect on this function's own
computation path.

---

## §10 `HeliacalAngle` (1636) & `swe_heliacal_angle` (1695)

### `HeliacalAngle(double Magn, double *dobs, double AziO, double AltM, double AziM, double JDNDaysUT, double AziS, double *dgeo, double *datm, int32 helflag, double *dangret, char *serr)` — swehel.c:1636–1693
Finds the **optimal** object-altitude/arcus-visionis pair (2-D search: for each candidate object
altitude `x`, `TopoArcVisionis` gives the arcus visionis `Arc`; minimize `Arc` over `x`) — this is
the "critical/most favorable altitude for first visibility" search.
```c
sunra = SunRA(JDNDaysUT, helflag, serr);
Lat = dgeo[1]; HeightEye = dgeo[2];
if (PLSV == 1) {                                          // DEAD — PLSV is the file-header constant 0, never 1
  dangret[0] = criticalangle;                              // = 0.0
  dangret[1] = criticalangle + Magn*2.492 + 13.447;
  dangret[2] = -(Magn*2.492 + 13.447);
  return OK;
}
// coarse scan, integer altitudes 2..20 degrees inclusive
minx=2; maxx=20; xmin=0; ymin=10000;
for (x = minx; x <= maxx; x++) {                           // x stepped as a double, but effectively integers 2,3,...,20 (19 samples)
  TopoArcVisionis(Magn, dobs, x, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, Lat, HeightEye, datm, helflag, &Arc, serr);
  if (Arc < ymin) { ymin = Arc; xmin = x; }
}
Xl = xmin - 1; xR = xmin + 1;                              // bracket the coarse minimum by ±1°
TopoArcVisionis(..., xR, ..., &Yr, ...);
TopoArcVisionis(..., Xl, ..., &Yl, ...);
// golden-section-like bisection refining the MINIMUM (not a zero-crossing) of Arc(x)
while (fabs(xR - Xl) > 0.1) {                              // coarser tolerance than TopoArcVisionis's own 0.001
  Xm = (xR + Xl) / 2.0;
  DELTAx = 0.025;
  xmd = Xm + DELTAx;
  TopoArcVisionis(..., Xm,  ..., &Ym,  ...);
  TopoArcVisionis(..., xmd, ..., &ymd, ...);
  if (Ym >= ymd) { Xl = Xm; Yl = Ym; }                      // Arc increasing rightward at Xm → minimum is to the right → discard left
  else           { xR = Xm; Yr = Ym; }                      // discard right
}
Xm = (xR + Xl) / 2.0;
Ym = (Yr + Yl) / 2.0;                                       // NOTE: averages the last iteration's Yl/Yr, not a fresh TopoArcVisionis(Xm) evaluation
dangret[1] = Ym;                                            // arcus visionis at the optimum
dangret[2] = Xm - Ym;                                       // Sun's altitude at the optimum (object_alt - arcus_visionis)
dangret[0] = Xm;                                            // object's altitude at the optimum
return OK;
```
**Coarse-scan loop variable is `double x` incremented by integer steps** (`for (x=minx; x<=maxx;
x++)` with `double minx=2, maxx=20` — the `x++` on a `double` still increments by exactly `1.0`,
so this reliably samples `x = 2.0, 3.0, ..., 20.0`, 19 evaluations; replicate as an integer loop
`2..=20` cast to `f64` in Rust, bit-identical). **`Ym` in the final `dangret` output is the
*averaged* `(Yr+Yl)/2`, not a re-evaluation of `TopoArcVisionis` at the converged `Xm`** — replicate
this averaging exactly (do not "simplify" by calling `TopoArcVisionis` once more at the final
`Xm`, which would generally give a slightly different value than the averaged `Yl`/`Yr` from the
last bisection iteration before convergence). This is a **minimum-finding bisection** (comparing
`Ym` vs. `ymd` — a one-sided finite-difference slope check — not a sign-change/root bisection like
`TopoArcVisionis`'s own search) — it assumes `Arc(x)` is unimodal (single minimum) over `[Xl,xR]`,
which holds physically (arcus visionis has one favorable-altitude minimum) given the coarse
1°-grid pre-scan already bracketed it. `dangret[2] = Xm - Ym` (object altitude minus arcus visionis
= implied Sun altitude at the optimum) is derived arithmetically, not from a fresh `AltS` lookup.

### `swe_heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr)` — swehel.c:1695–1705
```c
if (dgeo[2] < SEI_ECL_GEOALT_MIN || dgeo[2] > SEI_ECL_GEOALT_MAX) {   // [-500, 25000] meters, sweph.h
  serr = "location for heliacal events must be between %.0f and %.0f m above sea";
  return ERR;
}
swi_set_tid_acc(tjdut, helflag, 0, serr);
default_heliacal_parameters(datm, dgeo, dobs, helflag);
return HeliacalAngle(mag, dobs, azi_obj, alt_moon, azi_moon, tjdut, azi_sun, dgeo, datm, helflag, dret, serr);
```
Altitude-range validation (`SEI_ECL_GEOALT_MIN`/`MAX = -500/25000` m, same constants as
`docs/c-ref-riseset.md` §5.1's `swe_rise_trans_true_hor` bound) happens **before** any other work.
`swi_set_tid_acc` is called here (unlike `swe_topo_arcus_visionis`, §9) — consistent with
`swe_vis_limit_mag`'s pattern (§8) of setting the tidal-acceleration model whenever the function
internally calls `SunRA`/other ΔT-dependent paths (even though, per §7, `SunRA`'s live path
ignores ΔT entirely — the call is defensive/for-consistency rather than load-bearing given the
current `SunRA` implementation). No `swe_set_topo` call here (unlike `swe_vis_limit_mag`) —
`HeliacalAngle`/`TopoArcVisionis` take all geometry as caller-supplied azimuths/altitudes, no
internal position lookups requiring topocentric ephemeris state.

**`dangret`/`dret[0..2]` layout** (same 3-element output for both `HeliacalAngle` and
`swe_heliacal_angle`, just forwarded): `dret[0]` = object's altitude at the optimum (degrees),
`dret[1]` = arcus visionis at the optimum (degrees), `dret[2]` = implied Sun altitude at the
optimum (degrees, `= dret[0] - dret[1]`).

---

## §11 Porting notes for the stateless Rust port

### 11.1 The `SIMULATE_VICTORVB` finding (the single most important gotcha in this range)
`swephexp.h:451` unconditionally `#define`s `SIMULATE_VICTORVB 1` — there is no build path in this
codebase where it is undefined. This flips the sense of every `#ifdef`/`#ifndef` on it:
- **`SunRA` (swehel.c:563, `#ifndef SIMULATE_VICTORVB`)**: dead. The "real" `swe_calc`-based
  solar-RA computation never runs; **the shipped library always uses the crude
  `swe_degnorm((month + (day-1)/30.4 - 3.69) * 30)` calendar approximation** for `sunra` — a value
  threaded into `kOZ`, `ka`, `Bn`'s solar-cycle phase indirectly (no — `Bn` uses `JDNDaysUT`
  directly for its own phase, not `sunra`), `Bday`/`Btwi`/`Bm` (via `Deltam`→`kOZ`/`ka`, and
  `Btwi`'s own `kt` call), and `HeliacalAngle`/`TopoArcVisionis` transitively. **Port the crude
  formula as the real, correct behavior** — do not "improve" it to a real ephemeris call, or every
  downstream extinction/brightness golden-test value will diverge from the real C library.
- **`ka` (swehel.c:918, `#ifdef SIMULATE_VICTORVB`)**: live. RH is clamped to `(1e-8,
  99.999999998)` before the humidity-based aerosol formula whenever `VR==0` (the default,
  no-explicit-visibility-override case).
- **`default_heliacal_parameters` (swehel.c:1339, `#ifndef SIMULATE_VICTORVB`)**: dead. When the
  caller supplies an explicit `datm[0] > 0` (pressure), **no RH clamping happens at
  `default_heliacal_parameters` time** — the only RH clamp that actually executes anywhere in this
  call chain is the one inside `ka()` above, and only on the `VR==0` path.

### 11.2 Memoization caches: safe to drop, with one caveat
Every `static TLS` cache in this range (`call_swe_fixstar_mag`'s star-name/magnitude cache, `kOZ`'s
`(AltS,sunra)`→result cache, `ka`'s `(AltS,sunra)`→result cache, `Deltam`'s `(AltS,AltO,sunra)`→
result cache, `SunRA`'s `tjd`→result cache) is **pure memoization under the module's normal calling
convention** (where `Lat`/`HeightEye`/`datm`/`helflag` are held constant across a single logical
computation) — recomputing on every call in the stateless Rust port gives identical numeric
results in all realistic usage. The one caveat: `kOZ`/`ka`/`Deltam`'s cache keys omit some
parameters the formulas actually depend on (`Lat` for `kOZ`/`ka`; `Lat`,`HeightEye`,`datm`,
`helflag` for `Deltam`) — if a golden-test harness ever calls the C functions with those
"non-keyed" parameters varying while the keyed ones repeat, the C output would be **stale/wrong**
relative to a fresh computation; the Rust port (always recomputing) would then legitimately differ
from that pathological C sequence. This should never arise from `swe_vis_limit_mag`/
`swe_topo_arcus_visionis`/`swe_heliacal_angle`'s own call patterns (each holds `Lat`/`HeightEye`/
`datm`/`helflag` fixed throughout a single call), so it's not expected to surface in practice, but
document it if a golden-test mismatch on `Deltam`/`kOZ`/`ka` ever appears.

### 11.3 Global state genuinely threaded through (not just memoization)
- **`swe_set_topo(dgeo[0], dgeo[1], dgeo[2])`**: called in `Magnitude` (swehel.c:1118) and
  `swe_vis_limit_mag` (swehel.c:1481). Mutates `swed.topd`, read back by `swe_pheno_ut`/
  `swe_calc`-family functions internally whenever `SEFLG_TOPOCTR` is set. **Must become an explicit
  parameter** in the Rust port's equivalent phenomena/position calls — same pattern as
  `docs/c-ref-riseset.md`'s STATELESS PORT NOTE.
- **`swi_set_tid_acc(tjdut, helflag, 0, serr)`**: called in `swe_vis_limit_mag` (1478) and
  `swe_heliacal_angle` (1702), but notably **not** in `swe_topo_arcus_visionis` (§9) — sets the
  tidal-acceleration/ΔT model globally, consumed by `swe_deltat_ex` inside `ObjectLoc`/
  `azalt_cart`/`Magnitude`. Should map to whatever ΔT/tidal-acceleration config already exists in
  `EphemerisConfig` — check `docs/codebase-map.md` before adding a duplicate knob. Given §7's
  finding that `SunRA` ignores ΔT entirely on its live path, this setting's *only* live effect in
  this file is on the `swe_deltat_ex`-based `tjd_tt` conversion inside `ObjectLoc`/`azalt_cart`/
  `Magnitude` — i.e. it affects the ephemeris epoch used for Sun/Moon/object/planet positions, not
  `sunra` itself.
- **`swe_calc_ut`/`swe_calc`/`swe_fixstar`/`swe_pheno_ut`/`swe_azalt`/`swe_deltat_ex`**: all route
  through the existing stateless `Ephemeris` API in the Rust port (per this repo's `calc.rs`
  architecture) — the functions in this doc that call them (`calc_rise_and_set`, `SunRA`'s dead
  branch, `ObjectLoc`, `azalt_cart`, `Magnitude`) are the boundary where this module's pure math
  hands off to the ephemeris backend.

### 11.4 Pure-math vs. ephemeris-calling classification
**Pure math (no ephemeris calls)**: `mymin`, `mymax`, `Tanh`, `CVA`, `PupilDia`, `OpticFactor`,
`DeterObject`, `Kelvin`, `TopoAltfromAppAlt`, `AppAltfromTopoAlt`, `HourAngle`, `DistanceAngle`,
`kW`, `kOZ`, `kR`, `Sgn`, `ka`, `kt`, `Airmass`, `Xext`, `Xlay`, `TempEfromTempS`, `PresEfromPresS`,
`Deltam`, `Bn` (except `swe_revjul` for calendar decomposition — pure calendar math, not ephemeris
position lookup), `MoonsBrightness`, `MoonPhase`, `Bm`, `Btwi`, `Bday`, `Bcity`, `Bsky`,
`default_heliacal_parameters`, `VisLimMagn`, `tolower_string_star`, `TopoArcVisionis`,
`HeliacalAngle`. These form the bulk of the vision/sky-brightness model and should port as plain
functions with no `Ephemeris` dependency at all.

**Calls the ephemeris backend** (via `swe_calc`/`swe_calc_ut`/`swe_fixstar`/`swe_fixstar_mag`/
`swe_pheno_ut`/`swe_azalt`/`swe_rise_trans`/`swe_deltat_ex`): `call_swe_fixstar`,
`call_swe_fixstar_mag`, `call_swe_rise_trans`, `calc_rise_and_set`, `my_rise_trans`, `RiseSet`,
`SunRA` (dead branch only — live branch is pure `swe_revjul` calendar math), `ObjectLoc`,
`azalt_cart`, `Magnitude`, `swe_vis_limit_mag` (via `ObjectLoc`/`Magnitude`/`SunRA`),
`swe_heliacal_angle`/`swe_topo_arcus_visionis` (via `SunRA` only — dead branch — plus whatever
`default_heliacal_parameters` needs, which is itself pure). These route through the existing
stateless `Ephemeris` API rather than being reimplemented locally.

### 11.5 Dead code inventory (do not port)
- `call_swe_calc` (swehel.c:338–364, `#if 0`) — memoizing `swe_calc` wrapper.
- `fast_magnitude` (swehel.c:1129–1152, `#if 0`) — memoizing `Magnitude` wrapper.
- `OpticFactor`'s `#if 0` block (swehel.c:243–250) — superseded `OpticMag==0` handling, moved to
  `default_heliacal_parameters`.
- `VisLimMagn`'s `#if 0` block (swehel.c:1437–1441) — alternate Snellen-ratio-aware return.
- `HeliacalAngle`'s `PLSV==1` branch (swehel.c:1643–1648) — `PLSV` is the constant `0`, never `1`.
- `Airmass`/`Deltam`'s `else` branch (swehel.c:1054–1056) — `staticAirmass` is the constant `0`,
  so `Deltam` always takes the `if` branch; `Airmass` itself becomes unreachable in practice (but
  is a complete, correct function — port it anyway for parity if `staticAirmass` is ever flipped).
- Debug-only `#if SWEHEL_DEBUG` blocks (`SWEHEL_DEBUG=0`) in `OpticFactor` (284–296) and
  `VisLimMagn` (1429–1435) — `fprintf(stderr,...)` diagnostics only.
- `if ((0)) { ... }` blocks (not preprocessor `#if 0`, but a runtime `if` on a literal `0`) in
  `kOZ` (833–838), `Bsky` (1289–1294), `VisLimMagn` (1392–1397) — one-shot debug `printf`s guarded
  by a `static int a` latch; dead under any compiler that constant-folds `if(0)`, and even if not
  folded, never observably affects return values.
- The parallel `#if 0`-wrapped `SE_HELIACAL_LONG_SEARCH`/etc. flag-bit aliases (swephexp.h:453–467)
  — redundant with the live `SE_HELFLAG_*` names.
- All eight time constants `Y2D,D2Y,D2H,H2S,D2S,S2H,JC2D,M2S` and meteorological constants
  `DELTA,TempNulDiff,PressRef,MD,MW,GCR,LapseDA` and optics constants `GOpticMag,GOpticTrans,
  GBinocular,GOpticDia` and misc `GravitySource,REarthSource,DONE,Rb` — entirely unreferenced
  anywhere in the 3511-line file (§1).

### 11.6 Notable non-obvious quirks to preserve exactly
1. Three **different, non-consolidatable** scotopic/photopic brightness thresholds coexist:
   `1394` (`CVA`), `1645` (`OpticFactor`, `VisLimMagn`), and `BNIGHT=1479`/`BNIGHT_FACTOR=1.0`
   (the mixed-flag interval in `VisLimMagn`, which is currently always empty — see §8).
2. `ObjectLoc`'s `Angle==4` (AppAlt) branch calls `AppAltfromTopoAlt(xaz[1], datm[0], datm[1],
   helflag)` with `datm[0]`/`datm[1]` (pressure/temperature) apparently swapped relative to the
   function's own `(TopoAlt, TempE, PresE, helflag)` signature and every other call site in the
   file — replicate exactly, this is what the shipped library actually computes (§7).
3. `Bday`'s `kXS = Deltam(AltS, AltS, ...)` deliberately passes the Sun's own altitude as both
   `Deltam` arguments (extinction along the Sun's own sightline) — not a bug, matches the
   accompanying legacy BASIC-source comment (§6).
4. `MoonsBrightness`/`VisLimMagn` both declare a local named `log10` that actually holds `ln(10) =
   2.302585092994`, used as `log(x)/log10` (manual `log10`) — never call `f64::log10` as a
   substitute; replicate the two-`ln`-calls-divided form.
5. `calc_rise_and_set` uses the AU literal `1.49597870691e+11` for its own Sun/Moon disc-radius
   calc — a *different* AU value than the riseset module's DE431 `AUNIT = 1.49597870700e11` — do
   not unify these two AU constants across modules.
6. The `30.4`-day "average month" constant recurs in both `SunRA` and `Bn` — always exactly `30.4`,
   never a more precise value; keep as a shared-by-convention (not shared-by-code) literal in each
   function, since the two functions aren't otherwise related.
7. `my_rise_trans`'s latitude gate is `63°`; `swecl.c`'s `rise_set_fast` gate (riseset doc §4) is
   `60°`/`65°` — two independently-tuned thresholds in two different modules, do not merge.
