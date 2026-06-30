# C Reference: Fixed Stars — sweph.c

Porting reference for the fixed-star subsystem. Read this instead of the C source.
The Rust implementer never needs to open `sweph.c`, `sweph.h`, or `sefstars.txt` directly.

---

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_fixstar2` | sweph.c:6818–6876 | Yes — primary public entry |
| `swe_fixstar2_ut` | sweph.c:6878–6898 | Yes — UT-input variant |
| `swe_fixstar2_mag` | sweph.c:6911–6944 | Yes — magnitude lookup |
| `swe_fixstar` | sweph.c:7896–7955 | Legacy (file-scan approach); same math via `swi_fixstar_calc_from_record` |
| `swe_fixstar_ut` | sweph.c:7957–7977 | Legacy UT variant |
| `swe_fixstar_mag` | sweph.c:7990–8030 | Legacy magnitude lookup |
| `load_all_fixed_stars` | sweph.c:6324–6395 | Yes — catalog loader (called once) |
| `fixstar_format_search_name` | sweph.c:6154–6174 | Yes — normalizes input name |
| `fixstar_cut_string` | sweph.c:6211–6306 | Yes — parses one CSV record into struct |
| `save_star_in_struct` | sweph.c:6178–6190 | Yes — appends to in-memory array |
| `fixedstar_name_compare` | sweph.c:6193–6198 | Yes — qsort comparator (strcmp on skey) |
| `fstar_node_compare` | sweph.c:6201–6206 | Yes — bsearch comparator (strcmp on skey) |
| `search_star_in_list` | sweph.c:6674–6748 | Yes — looks up star by name/number/bayer |
| `get_builtin_star` | sweph.c:6750–6803 | Yes — hardcoded fallback for 8 ayanamsa stars |
| `fixstar_calc_from_struct` | sweph.c:6407–6669 | Yes — THE CORE: catalog → xx[6] |

---

## sefstars.txt File Format (sweph.c:6208–6306, ephe/sefstars.txt)

### Column specification

Each data line is a comma-separated record with **at least 14 fields** (fields 15–16 are
optional and ignored). Fields are numbered 0–15:

| Field | Index | Content | Units |
|---|---|---|---|
| traditional name | 0 | Traditional star name, or same as Bayer if none | string (≤ 40 chars) |
| Bayer/Flamsteed | 1 | Bayer or Flamsteed designation | string (≤ 40 chars) |
| epoch | 2 | Coordinate epoch: `"1950"`, `"2000"`, or `"ICRS"` | string parsed via `atof` |
| RA hours | 3 | Right ascension hours | integer/float |
| RA minutes | 4 | Right ascension minutes | integer/float |
| RA seconds | 5 | Right ascension seconds | float |
| Dec degrees | 6 | Declination degrees (sign from leading `-`) | signed integer/float |
| Dec minutes | 7 | Declination arcminutes | float |
| Dec seconds | 8 | Declination arcseconds | float |
| RA proper motion | 9 | Proper motion in RA, in 0.001 arcsec/yr × cos(dec) | float |
| Dec proper motion | 10 | Proper motion in Dec, in 0.001 arcsec/yr | float |
| radial velocity | 11 | Radial velocity | km/s |
| parallax | 12 | Annual parallax | 0.001 arcsec (mas) |
| magnitude | 13 | Visual magnitude V | float (may be negative or 999.99) |
| DM zone | 14 | Durchmusterung zone (ignored since SE 2.07) | integer |
| DM number | 15 | Durchmusterung number (ignored since SE 2.07) | integer |

**Comment lines**: Lines starting with `#` are skipped entirely (sweph.c:6350).

**Blank-line handling**: Lines starting with `\n`, `\r`, or `\0` are skipped in `load_all_fixed_stars`
(sweph.c:6351–6353). NOTE: the legacy `swe_fixstar` file-scan path (swi_fixstar_load_record) does NOT
skip blank lines; a blank line would trigger "data corrupted" because `swi_cutstr` finds no comma. The
Rust port uses the `load_all_fixed_stars` approach and MUST skip blank lines silently.

**Field count**: If fewer than 14 fields are found after splitting, the record is rejected with
`serr = "data of star 'name,bayer' incomplete"` (sweph.c:6225–6234).

**Field 2 — epoch parsing**: `epoch = atof(cpos[2])` (sweph.c:6248). Since `atof("ICRS") == 0.0`,
`atof("2000") == 2000.0`, `atof("1950") == 1950.0`, the epoch is a float discriminator:
- `0.0` → ICRS (already ICRS J2000, no FK conversion needed)
- `2000.0` → FK5 J2000 (needs FK5→ICRS conversion)
- `1950.0` → FK4 B1950 (needs FK4→FK5, precess B1950→J2000, then FK5→ICRS)

**Field 6 — declination sign**: Sign is read from the raw string `sde_d = cpos[6]` using
`strchr(sde_d, '-')`. This means `"-00"` gives negative declination even if `atof` would
return 0.0 (sweph.c:6267–6271):
```c
if (strchr(sde_d, '-') == NULL) {
    de = de_s / 3600.0 + de_m / 60.0 + de_d;
} else {
    de = -de_s / 3600.0 - de_m / 60.0 + de_d;
}
```
All three components of negative declination are negated; `de_d` already carries the sign.

**Parallax sign fix**: `if (parall < 0) parall = -parall` (sweph.c:6260). The field is always
treated as positive (handles historical bug in old files like "old Rasalgheti").

### Verbatim sample lines

```
# comment line — skipped entirely
alfTau,alfTau,ICRS,04,35,55.23907,+16,30,33.4885,63.45,-188.94,54.398,48.94,0.86
Aldebaran,alfTau,ICRS,04,35,55.23907,+16,30,33.4885,63.45,-188.94,54.398,48.94,0.86
alfSco,alfSco,ICRS,16,29,24.45970,-26,25,55.2094,-12.11,-23.30,-3.50,5.89,0.91
Antares,alfSco,ICRS,16,29,24.45970,-26,25,55.2094,-12.11,-23.30,-3.50,5.89,0.91
alfCMa,alfCMa,ICRS,06,45,08.91728,-16,42,58.0171,-546.01,-1223.07,-5.50,379.21,-1.46
Sirius,alfCMa,ICRS,06,45,08.91728,-16,42,58.0171,-546.01,-1223.07,-5.50,379.21,-1.46
SgrA*,SgrA*,ICRS,17,45,40.03599,-29,00,28.1699,-2.755718425,-5.547,0.0,0.125,999.99,0,0
Galactic Center,SgrA*,ICRS,17,45,40.03599,-29,00,28.1699,-2.755718425,-5.547,0.0,0.125,999.99,0,0
GA,GA,2000,16,15,02.836,-60,53,22.54,0.000,0.00,0.0,0.0000159,999.99,0,0
alfHer,alfHer,ICRS,17,14,38.85818,+14,23,25.2262,-7.32,36.07,-32.09,9.07,3.06, 14, 3207
```

Note: `alfHer` and similar lines have trailing DM fields with spaces around the comma. The
`swi_cutstr` function handles any whitespace around commas correctly.

### Unit conversions inside `fixstar_cut_string` (sweph.c:6264–6304)

```
// RA and Dec → degrees
ra = (ra_s / 3600.0 + ra_m / 60.0 + ra_h) * 15.0    // hours→degrees
de = de_s/3600.0 + de_m/60.0 + de_d                  // arcmin/sec → degrees (with sign)

// For new file format (sefstars.txt):
// file value in 0.001 arcsec/yr = mas/yr; convert to degrees/century:
ra_pm = ra_pm / 10.0 / 3600.0     // (0.001 arcsec/yr × 100 yr/century) / 3600 arcsec/° = deg/century
de_pm = de_pm / 10.0 / 3600.0     // same
parall /= 1000.0                   // 0.001 arcsec → arcsec

// For old file format (fixstars.cat):
ra_pm = ra_pm * 15 / 3600.0       // seconds-of-time/yr → degrees/yr (NOT century)
de_pm = de_pm / 3600.0            // arcsec/yr → degrees/yr (NOT century)

// Parallax → degrees
if (parall > 1):
    parall = (1 / parall / 3600.0)   // arcsec parallax → degrees (large value = nearby star)
else:
    parall /= 3600.0                  // arcsec → degrees

// Radial velocity → AU/century
radv *= KM_S_TO_AU_CTY   // = 21.095 (km/s → AU/century)

// All angular quantities → radians
ra *= DEGTORAD
de *= DEGTORAD
ra_pm *= DEGTORAD
de_pm *= DEGTORAD
ra_pm /= cos(de)          // catalog stores RA proper motion × cos(dec); divide to get true spherical rate
parall *= DEGTORAD
```

The `ra_pm /= cos(de)` step (sweph.c:6295) converts from "RA motion projected on the great
circle" (which is what catalogs store) to the pure spherical RA rate. This is the inverse of
the catalog convention `μα × cos(δ)`.

---

## `struct fixed_star` (sweph.h:773–779)

```c
#define SWI_STAR_LENGTH 40   // sweph.h:772

struct fixed_star {
    char   skey[SWI_STAR_LENGTH + 2];      // search key (42 bytes); comma-prefixed for Bayer entries
    char   starname[SWI_STAR_LENGTH + 1];  // traditional name (41 bytes), or "" if none
    char   starbayer[SWI_STAR_LENGTH + 1]; // Bayer/Flamsteed designation (41 bytes)
    char   starno[10];                     // NOT used in new path; always "" in load_all_fixed_stars
    double epoch;   // 0.0=ICRS, 2000.0=FK5 J2000, 1950.0=FK4 B1950; from atof(epoch_field)
    double ra;      // right ascension in RADIANS (J2000 equatorial)
    double de;      // declination in RADIANS
    double ramot;   // RA proper motion in RADIANS/century (already /cos(de) applied)
    double demot;   // Dec proper motion in RADIANS/century
    double radvel;  // radial velocity in AU/century
    double parall;  // annual parallax in RADIANS (positive; 0 if not known)
    double mag;     // visual magnitude V (999.99 = no magnitude)
};
```

The `starno` field exists in the struct definition but is never populated by `load_all_fixed_stars`
(only by the old file-scan path). In the new path, sequential-number lookup works by direct
array index, not by `skey`.

In the `swe_data` global (sweph.h:843–846):
```c
AS_BOOL n_fixstars_real;      // count of distinct Bayer records (unique stars)
AS_BOOL n_fixstars_named;     // count of records with traditional-name keys
AS_BOOL n_fixstars_records;   // total records (n_fixstars_real + n_fixstars_named)
struct fixed_star *fixed_stars; // heap-allocated array, sorted by skey
```

(Note: the types are declared `AS_BOOL` in C, which is `int` — these are actually integer counts
stored in a boolean-typed field. This is a C bug; the Rust port should use `usize` or `u32`.)

---

## Catalog Loading — `load_all_fixed_stars` (sweph.c:6324–6395)

### Call semantics

`swe_fixstar2` calls `load_all_fixed_stars(serr)` on every invocation (sweph.c:6836).
The function is idempotent: if `swed.n_fixstars_records > 0`, it returns `-2` immediately
(sweph.c:6333–6335). In the Rust stateless design, load once and store in `Ephemeris`.

### File search

```
Try SE_STARFILE ("sefstars.txt") via swi_fopen(SEI_FILE_FIXSTAR=4, ..., ephepath):
    if not found:
        set is_old_starfile = TRUE
        try SE_STARFILE_OLD ("fixstars.cat"):
            if not found: set is_old_starfile = FALSE; return ERR
```

`swi_fopen` searches `ephepath` (a colon-separated path list on Unix, semicolon on Windows).
The Rust port reuses the same path-search logic already in the ephemeris backend.

### Parsing loop (sweph.c:6348–6393)

```
nstars = 0      // unique Bayer records
nrecs = 0       // total array entries
nnamed = 0      // traditional-name records
last_starbayer = ""

for each line in file:
    skip if starts with '#', '\n', '\r', '\0'
    parse line → fstdata via fixstar_cut_string()
    if fixstar_cut_string returns ERR: return ERR

    // Traditional-name record (if field 0 is non-empty)
    if fstdata.starname != "":
        nrecs++; nnamed++
        fstdata.skey = lowercase(starname) with spaces removed
        append copy to fixed_stars array

    // Bayer-designation record
    // Dedup: skip if same Bayer as previous star (consecutive duplicates only)
    if fstdata.starbayer == last_starbayer:
        continue
    nstars++; nrecs++
    fstdata.skey = "," + starbayer with spaces removed   // comma prefix
    last_starbayer = fstdata.starbayer
    append copy to fixed_stars array

// Sort the entire array by skey (strcmp order)
qsort(fixed_stars, nrecs, sizeof(struct fixed_star), fixedstar_name_compare)
```

**Array layout after sort** (sorted lexicographically by `skey`):

The sorted array is a flat array of all records. Because comma `','` (ASCII 44) sorts BEFORE
letters, Bayer entries (prefix `,`) come first, then traditional-name entries (lowercase). This
means searching Bayer vs. traditional-name requires either knowing the range or scanning all.

In practice, `search_star_in_list` uses the known counts:
- Bayer entries: `fixed_stars[0..n_fixstars_real]`
- Named entries: `fixed_stars[n_fixstars_real..n_fixstars_real+n_fixstars_named]`

File size: `sefstars.txt` has approximately 2869 non-comment data lines → roughly 1000–1500
unique Bayer entries and 1000–1500 named entries (many stars have 3–5 name variants).

### Old star file (`fixstars.cat`)

When `is_old_starfile == TRUE`, the proper-motion and parallax units differ (see §File Format
above). The Rust port should check the filename found and set a flag. If the old file is not
present and `sefstars.txt` is not found, return an error.

---

## Search Logic — `search_star_in_list` (sweph.c:6674–6748)

### Entry normalization: `fixstar_format_search_name` (sweph.c:6154–6174)

Applied to the caller's `star` string before any search:
1. Copy to `sstar` (max `SWI_STAR_LENGTH` chars)
2. Remove all space characters (in-place)
3. For the part BEFORE the first comma: convert to lowercase
4. Part after first comma (Bayer designation): LEFT UNCHANGED (case-sensitive!)
5. If result is empty: return ERR "swe_fixstar(): star name empty"

Example: `"Alpha  Virginis,alfVir"` → `"alphavirginis,alfVir"`

### Three search modes (sweph.c:6682–6747)

**Mode 1: Sequential number** — first char of `sstar` is a digit

```
star_nr = atoi(sstar)
if star_nr > n_fixstars_real: return ERR "sequential fixed star number N is not available"
*stardata = fixed_stars[star_nr - 1]    // 1-indexed; uses Bayer-record entries (first n_fixstars_real)
return OK
```

The first `n_fixstars_real` entries in the sorted array are the Bayer-key entries. After sort
they begin with `,alfCMa`, `,alfLeo`, etc. Entry `[0]` is the alphabetically-first Bayer key.

**Mode 2: Bayer designation** — `sstar` starts with `,` OR `sstar` contains `,` after the trad. name

When the raw input starts with `,` (e.g. `",alfVir"`, `",SgrA*"`):
- `is_bayer = TRUE`
- `searchkey = sstar` (already comma-prefixed)
- Binary-search `fixed_stars[0..n_fixstars_real]` (the Bayer-key range)

When the input is `"trad_name,bayer"` (both parts present):
- After `fixstar_format_search_name`, if `sstar` contains a comma, `search_star_in_list`
  strips the trad-name part: `swi_strcpy(sstar, strchr(sstar, ','))` → `sstar = ",bayer_part"`
- `is_bayer = TRUE`, then binary-search the Bayer range

**Mode 3: Traditional name** — no leading digit, no comma in sstar

```
searchkey = sstar    // lowercase, space-stripped
binary-search fixed_stars[n_fixstars_real .. n_fixstars_real+n_fixstars_named]
```

**Wildcard mode (`%`-suffix)** — `sstar` ends with `%` (sweph.c:6701–6721):

Linear scan of the named-name range (not binary search). The `%` is stripped; the remaining
prefix must match via `strncmp`. Only traditional-name range is searched. Must have `%` as
the LAST character; otherwise returns ERR "invalid search string".

### Builtin-star fallback — `get_builtin_star` (sweph.c:6750–6803)

Called in `swe_fixstar2` BEFORE `search_star_in_list`. If the file is unavailable, or for the
named ayanamsa stars, these hardcoded records are used. The catalog file does NOT need to be
present for these stars. Conditions (sweph.c:6754–6801):

| Test on `star` (original input, not normalized) | Bayer key | Hardcoded record |
|---|---|---|
| starts with "spica" or "Spica" | `"spica"` | `"Spica,alVir,ICRS,13,25,11.57937,-11,09,40.7501,-42.35,-30.67,1,13.06,0.97,-10,3672"` |
| contains `",zePsc"` or starts with "revati"/"Revati" | `"revati"` | `"Revati,zePsc,ICRS,01,13,43.88735,+07,34,31.2745,145,-55.69,15,18.76,5.187,06,174"` |
| contains `",deCnc"` or starts with "pushya"/"Pushya" | `"pushya"` | `"Pushya,deCnc,ICRS,08,44,41.09921,+18,09,15.5034,-17.67,-229.26,17.14,24.98,3.94,18,2027"` |
| contains `",deCnc"` (second check for Sheoran) | `"pushya"` | same as above |
| contains `",laSco"` or starts with "mula"/"Mula" | `"mula"` | `"Mula,laSco,ICRS,17,33,36.52012,-37,06,13.7648,-8.53,-30.8,-3,5.71,1.62,-37,11673"` |
| contains `",SgrA*"` | `",SgrA*"` | `"Gal. Center,SgrA*,2000,17,45,40.03599,-29,00,28.1699,-2.755718425,-5.547,0.0,0.125,999.99,0,0"` |
| contains `",GP1958"` | `",GP1958"` | `"Gal. Pole IAU1958,GP1958,1950,12,49,0.0,27,24,0.0,0.0,0.0,0.0,0.0,0.0,0,0"` |
| contains `",GPol"` | `",GPol"` | `"Gal. Pole,GPol,ICRS,12,51,36.7151981,27,06,11.193172,0.0,0.0,0.0,0.0,0.0,0,0"` |

Note: The `",deCnc"` check appears twice (for Pushya and Sheoran); both return the same record.
The `",GPol"` check also appears twice (for GALEQU_TRUE and GALEQU_MULA); both return `GPol`.

Note: The Galactic Centre builtin has **epoch `"2000"`** (not `"ICRS"`), meaning FK5→ICRS
conversion will be applied even though it is close to the file's ICRS data (sweph.c:6783).
This is a subtle inconsistency in the C code — the ICRS position of SgrA* from SIMBAD is
stored with epoch `"2000"` in the builtin, triggering FK5→ICRS. Replicate exactly.

### Cache (per-call)

`swe_fixstar2` caches the last-used star in `static TLS struct fixed_star last_stardata` and
`static TLS char slast_starname[]` (sweph.c:6826–6827). If the normalized search name matches
`slast_starname`, the struct is reused without searching. In the stateless Rust port, this
cache is not needed (callers are stateless); omit it or replace with a local variable.

---

## Position Computation — `fixstar_calc_from_struct` (sweph.c:6407–6669)

This is the core algorithm. All field values are already converted to radians/AU/century
by `fixstar_cut_string`. Line numbers below are from `fixstar_calc_from_struct`.

### Setup (sweph.c:6419–6453)

```
iflgsave = iflag              // save for output (speed zeroing, ephe mask)
iflag |= SEFLG_SPEED          // always compute speed internally
iflag = plaus_iflag(iflag, -1, tjd, serr)
epheflag = iflag & SEFLG_EPHMASK

// Reset ephemeris files if epheflag changed
if swed.last_epheflag != epheflag:
    free_planets(), close_jpl_file(), close_sweph_files()
    swed.last_epheflag = epheflag

// Default sidereal mode if needed
if (iflag & SEFLG_SIDEREAL) && !swed.ayana_is_set:
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0)

// Precompute ecliptic obliquity (oec2000, oec) and nutation (swed.nut)
swi_check_ecliptic(tjd, iflag)
swi_check_nutation(tjd, iflag)
```

### Step 1: Epoch offset (sweph.c:6459–6463)

```
if epoch == 1950:
    t = tjd - B1950    // B1950 = 2433282.42345905 JD
else:                   // epoch 2000.0 OR ICRS (epoch 0.0)
    t = tjd - J2000    // J2000 = 2451545.0 JD
// t is in DAYS since epoch
```

### Step 2: Build initial polar+speed vector (sweph.c:6464–6477)

```
x[0] = ra          // RA in radians (epoch J2000 or B1950)
x[1] = de          // Dec in radians
x[2] = rdist       // distance in AU
x[3] = ramot / 36525.0   // RA proper motion in rad/day (ramot is rad/century)
x[4] = demot / 36525.0   // Dec proper motion in rad/day
x[5] = radvel / 36525.0  // radial velocity in AU/day (radvel is AU/century)

// Distance:
if parall == 0:
    rdist = 1000000000   // 1e9 AU (effectively at infinity)
else:
    rdist = (1.0 / (parall * RADTODEG * 3600)) * PARSEC_TO_AUNIT
    // parall is in radians; × RADTODEG × 3600 converts to arcsec;
    // 1/arcsec = parallax in arcsec = distance in parsecs;
    // × PARSEC_TO_AUNIT (206264.8062471) = AU
```

`PARSEC_TO_AUNIT = 206264.8062471` (= 648000/π per IAU 2016 B2 resolution, sweph.h:288).

### Step 3: Polar → Cartesian with space-motion (sweph.c:6479)

```
swi_polcart_sp(x, x)
// x[0..2] = Cartesian position at epoch (in AU)
// x[3..5] = Cartesian space-motion vector (AU/day)
// The velocity components are the full 3-D space velocity, not just proper motion.
```

`swi_polcart_sp` converts from polar (lon, lat, r, dlon, dlat, dr) to Cartesian (X, Y, Z, dX, dY, dZ),
correctly combining angular and radial motion components.

### Step 4: FK4/FK5/ICRS frame corrections (sweph.c:6483–6496)

```
if epoch == 1950:
    swi_FK4_FK5(x, B1950)           // FK4 → FK5 rotation at B1950
    swi_precess(x,   B1950, 0, J_TO_J2000)   // precess position B1950 → J2000
    swi_precess(x+3, B1950, 0, J_TO_J2000)   // precess velocity B1950 → J2000

// For epoch != 0 (FK4 or FK5 data, not ICRS):
if epoch != 0:
    swi_icrs2fk5(x, iflag, TRUE)     // TRUE=backward: FK5 → ICRS
    if swi_get_denum(SEI_SUN, iflag) >= 403:
        swi_bias(x, J2000, SEFLG_SPEED, FALSE)  // ICRS → J2000 frame bias at J2000
// For epoch == 0 (ICRS): skip both steps; frame bias applied later (step 9)
```

After this step, `x[0..5]` is in the equatorial J2000 ICRS frame (or just after frame bias for non-ICRS data).

### Step 5: Earth and Sun positions (sweph.c:6501–6508)

```
dt = PLAN_SPEED_INTV * 0.1   // = 0.00001 days (about 0.864 seconds)

// Earth/Sun barycentric positions for parallax, light deflection, aberration:
if not (SEFLG_BARYCTR) and not (SEFLG_HELCTR and SEFLG_MOSEPH):
    main_planet_bary(tjd - dt, SEI_EARTH, epheflag, iflag, NO_SAVE,
                     xearth_dt, xearth_dt, xsun_dt, NULL, serr)
    main_planet_bary(tjd,      SEI_EARTH, epheflag, iflag, DO_SAVE,
                     xearth,   xearth,   xsun,   NULL, serr)
```

Both `xearth` and `xsun` are 6-element arrays (position + velocity), barycentric equatorial J2000.

### Step 6: Observer position (sweph.c:6513–6529)

```
if SEFLG_TOPOCTR:
    xobs_dt = swi_get_observer(tjd - dt, iflag | SEFLG_NONUT, NO_SAVE)
    xobs    = swi_get_observer(tjd,      iflag | SEFLG_NONUT, NO_SAVE)
    xobs[i] += xearth[i]        // topocentric → barycentric
    xobs_dt[i] += xearth_dt[i]
else if not SEFLG_BARYCTR and not (SEFLG_HELCTR+SEFLG_MOSEPH):
    xobs[i] = xearth[i]         // geocenter (barycentric)
    xobs_dt[i] = xearth_dt[i]
// else: barycentric → xobs stays unset (xpo = NULL in step 7)
```

### Step 7: Apply proper motion and parallax (sweph.c:6534–6557)

Select the parallax reference point (`xpo`):
```
if SEFLG_HELCTR and SEFLG_MOSEPH: xpo = NULL           // heliocentric Moshier: no parallax
elif SEFLG_HELCTR:                 xpo = xsun / xsun_dt  // heliocentric: parallax from Sun
elif SEFLG_BARYCTR:                xpo = NULL           // barycentric: no parallax
else:                              xpo = xobs / xobs_dt  // geocentric: parallax from observer
```

Apply motion over time `t` (in days from catalog epoch):
```
if xpo == NULL:
    for i in 0..3:
        x[i] += t * x[i+3]        // position only: add proper motion × elapsed time
else:
    for i in 0..3:
        x[i] += t * x[i+3]        // add proper motion
        x[i] -= xpo[i]            // subtract observer position (parallax)
        x[i+3] -= xpo[i+3]        // subtract observer velocity
```

After this step, `x[0..2]` is the geocentric direction vector to the star (in AU), and
`x[3..5]` is the rate of change of that direction.

### Step 8: Gravitational deflection (sweph.c:6561–6563)

```
if not SEFLG_TRUEPOS and not SEFLG_NOGDEFL:
    swi_deflect_light(x, 0, iflag & SEFLG_SPEED)
    // Note: second argument is dt=0 (not DEFL_SPEED_INTV as for planets).
    // Speed computed from position differences at dt=0 internally.
```

### Step 9: Annual aberration (sweph.c:6568–6569)

```
if not SEFLG_TRUEPOS and not SEFLG_NOABERR:
    swi_aberr_light_ex(x, xpo, xpo_dt, dt, iflag & SEFLG_SPEED)
    // Uses PLAN_SPEED_INTV * 0.1 = 0.00001 days for speed approximation.
    // Comment in C: "speed is incorrect !!!" — the speed approximation is inexact here.
```

### Step 10: ICRS → J2000 frame bias (sweph.c:6571–6573)

```
if not SEFLG_ICRS and (denum >= 403 or SEFLG_BARYCTR):
    swi_bias(x, tjd, iflag, FALSE)
    // Applied ONLY for ICRS catalog data (epoch == 0) or all data with DE < 403.
    // For epoch != 0 with DE >= 403, frame bias was already applied in Step 4.
```

### Step 11: Save J2000 equatorial Cartesian (sweph.c:6574–6576)

```
for i in 0..6:
    xxsv[i] = x[i]   // saved for sidereal transform below
```

### Step 12: Precession J2000 → equinox of date (sweph.c:6581–6587)

```
if not SEFLG_J2000:
    swi_precess(x,   tjd, iflag, J2000_TO_J)
    if SEFLG_SPEED:
        swi_precess_speed(x, tjd, iflag, J2000_TO_J)
    oe = &swed.oec        // obliquity at date
else:
    oe = &swed.oec2000    // obliquity at J2000
```

### Step 13: Nutation (sweph.c:6591–6592)

```
if not SEFLG_NONUT:
    swi_nutate(x, iflag, FALSE)   // applies nutation to equatorial Cartesian of date
```

### Step 14: Equatorial → ecliptic (sweph.c:6602–6611)

```
if not SEFLG_EQUATORIAL:
    swi_coortrf2(x,   x, oe->seps, oe->ceps)    // rotate by obliquity (position)
    if SEFLG_SPEED:
        swi_coortrf2(x+3, x+3, oe->seps, oe->ceps)  // (velocity)
    if not SEFLG_NONUT:
        swi_coortrf2(x,   x, swed.nut.snut, swed.nut.cnut)   // nutation in ecliptic
        if SEFLG_SPEED:
            swi_coortrf2(x+3, x+3, swed.nut.snut, swed.nut.cnut)
```

After this step, `x` is in ecliptic Cartesian coordinates of date.

### Step 15: Sidereal transform (sweph.c:6616–6642)

```
if SEFLG_SIDEREAL:
    if sidd.sid_mode & SE_SIDBIT_ECL_T0:
        // Rigorous: project onto ecliptic of t0
        swi_trop_ra2sid_lon(xxsv, x, xxsv, iflag)
        if SEFLG_EQUATORIAL:
            x[i] = xxsv[i]   // use equatorial sidereal output
    elif sidd.sid_mode & SE_SIDBIT_SSY_PLANE:
        // Project onto solar system equatorial plane
        swi_trop_ra2sid_lon_sosy(xxsv, x, iflag)
        if SEFLG_EQUATORIAL:
            x[i] = xxsv[i]
    else:
        // Traditional: subtract ayanamsa longitude
        swi_cartpol_sp(x, x)                // ecliptic Cartesian → polar
        swi_get_ayanamsa_with_speed(tjd, iflag, daya, serr)   // daya[0]=lon, daya[1]=speed (deg)
        x[0] -= daya[0] * DEGTORAD          // subtract ayanamsa (radians)
        x[3] -= daya[1] * DEGTORAD          // subtract ayanamsa speed
        swi_polcart_sp(x, x)                // back to Cartesian
```

The three sidereal paths are identical to the planet pipeline (see `c-ref-ayanamsa.md`
§"Sidereal Projection in the Calc Pipeline").

### Step 16: Cartesian → polar (sweph.c:6647–6648)

```
if not SEFLG_XYZ:
    swi_cartpol_sp(x, x)
    // x[0] = longitude (rad), x[1] = latitude (rad), x[2] = distance (AU)
    // x[3] = lon speed (rad/day), x[4] = lat speed (rad/day), x[5] = dist speed (AU/day)
```

### Step 17: Radians → degrees (sweph.c:6652–6657)

```
if not SEFLG_RADIANS and not SEFLG_XYZ:
    for i in [0, 1]:
        x[i] *= RADTODEG    // position angles: lon, lat
        x[i+3] *= RADTODEG  // speed components: dlon/dt, dlat/dt
    // x[2] (distance AU) and x[5] (dist speed) stay in AU / AU/day
```

### Step 18: Copy to output and zero speeds if not requested (sweph.c:6658–6667)

```
for i in 0..6:
    xx[i] = x[i]

if not (iflgsave & SEFLG_SPEED):
    xx[3] = xx[4] = xx[5] = 0.0   // zero speeds if caller didn't request them

if (iflgsave & SEFLG_EPHMASK) == 0:
    iflag &= ~SEFLG_DEFAULTEPH    // don't return chosen ephe if none was requested
iflag &= ~SEFLG_SPEED             // remove internally-forced speed bit from return value
return iflag
```

### Output `xx[6]` contents

Depends on flags. Default (ecliptic, polar, degrees):
- `xx[0]` = ecliptic longitude in degrees
- `xx[1]` = ecliptic latitude in degrees
- `xx[2]` = distance in AU
- `xx[3]` = ecliptic longitude speed in degrees/day (0 if `SEFLG_SPEED` not set)
- `xx[4]` = ecliptic latitude speed in degrees/day
- `xx[5]` = distance speed in AU/day

With `SEFLG_EQUATORIAL`: RA/Dec instead of lon/lat.
With `SEFLG_XYZ`: Cartesian XYZ in AU (speeds in AU/day).
With `SEFLG_RADIANS`: angles in radians instead of degrees.

---

## Public Entry Points

### `swe_fixstar2` (sweph.c:6818–6876)

```
swe_fixstar2(star, tjd, iflag, xx, serr):
    load_all_fixed_stars(serr)   // idempotent; loads on first call
    fixstar_format_search_name(star, sstar, serr)   // normalize
    if sstar == slast_starname: stardata = last_stardata; goto found   // cache hit
    if get_builtin_star(star, sstar, srecord):       // builtin fallback
        fixstar_cut_string(srecord, star, &stardata, serr)
        goto found
    search_star_in_list(sstar, &stardata, serr)      // bsearch/sequential
    found:
        slast_starname = sstar; last_stardata = stardata
        fixstar_calc_from_struct(&stardata, tjd, iflag, star, xx, serr)
    // star is modified in-place: set to "tradname,bayerdesig"
    // return: iflag on success, ERR on failure, xx zeroed on error
```

Note: `star` (the input buffer) is modified: on return it contains `"tradname,bayer"` or just
the matching record's fields. The caller must provide a buffer of at least `2 * SE_MAX_STNAME`
(= 512) bytes per the API spec (swephexp.h:299–303).

### `swe_fixstar2_ut` (sweph.c:6878–6898)

```
swe_fixstar2_ut(star, tjd_ut, iflag, xx, serr):
    iflag = plaus_iflag(iflag, -1, tjd_ut, serr)
    epheflag = iflag & SEFLG_EPHMASK
    if epheflag == 0: epheflag = SEFLG_SWIEPH; iflag |= SEFLG_SWIEPH
    deltat = swe_deltat_ex(tjd_ut, iflag, serr)
    retflag = swe_fixstar2(star, tjd_ut + deltat, iflag, xx, serr)
    // If ephemeris used differs from requested, recompute delta-T with correct ephe:
    if retflag != ERR and (retflag & SEFLG_EPHMASK) != epheflag:
        deltat = swe_deltat_ex(tjd_ut, retflag, NULL)
        retflag = swe_fixstar2(star, tjd_ut + deltat, iflag, xx, NULL)
    return retflag
```

### `swe_fixstar2_mag` (sweph.c:6911–6944)

No position calculation. Same lookup path as `swe_fixstar2` (builtin fallback → list search),
but only extracts `stardata.mag`. Returns `OK` or `ERR`. On success, also updates `star` in-place
with `"tradname,bayer"`.

```
swe_fixstar2_mag(star, mag, serr):
    load_all_fixed_stars(serr)
    fixstar_format_search_name(star, sstar, serr)
    if sstar == slast_starname: stardata = last_stardata; goto found
    search_star_in_list(sstar, &stardata, serr)   // no builtin-star check here!
    found:
        last_stardata = stardata
        *mag = stardata.mag
        star = "tradname,bayer"
        return OK
```

Note: `swe_fixstar2_mag` does NOT call `get_builtin_star` (unlike `swe_fixstar2`). If the
catalog file is not available, builtin stars will not be found via this function.

### Legacy `swe_fixstar` / `swe_fixstar_ut` / `swe_fixstar_mag` (sweph.c:7896–8030)

These are the older API variants. They do NOT call `load_all_fixed_stars`; instead they
scan the file line-by-line on every call (via `swi_fixstar_load_record`) and call
`swi_fixstar_calc_from_record` (which re-parses the record every time). The math is identical
to `fixstar_calc_from_struct`. Port `swe_fixstar2` and friends; the legacy functions can be
thin wrappers that call `swe_fixstar2` for the Rust port, or can be omitted if not needed.

---

## `swe_fixstar2_mag` — Magnitude Lookup

See §Public Entry Points above. Return value is `OK` (0) on success, `ERR` (-1) on failure.
A magnitude of `999.99` signals "no magnitude available" in the catalog (see SgrA* in the file).

---

## Fixed-Star Ayanamsa Integration (sweph.c:3002–3142)

The 12 fixed-star ayanamsas in `swi_get_ayanamsa_ex` each call `swe_fixstar()` (the legacy
API). The Rust port can call `swe_fixstar2` (or its equivalent) instead — the results are
identical. See `c-ref-ayanamsa.md` §"Fixed-Star Ayanamsas" for the complete table.

### Flag construction (sweph.c:3007–3028)

After `plaus_iflag`:
```
iflag = (iflag & SEFLG_EPHMASK) | SEFLG_NONUT   // ephe + nonut only
epheflag = iflag & SEFLG_EPHMASK
otherflag = original_iflag & ~SEFLG_EPHMASK

iflag_galequ = iflag | SEFLG_TRUEPOS     // for galactic pole queries (always true position)

iflag_true = iflag                        // starts from stripped iflag
if otherflag & SEFLG_TRUEPOS:  iflag_true |= SEFLG_TRUEPOS
if otherflag & SEFLG_NOABERR:  iflag_true |= SEFLG_NOABERR
if otherflag & SEFLG_NOGDEFL:  iflag_true |= SEFLG_NOGDEFL
```

### Per-mode dispatch (sweph.c:3049–3142)

| Index | Constant | star key | iflag arg | `*daya` formula |
|---|---|---|---|---|
| 27 | TRUE_CITRA | `"Spica"` | `iflag_true` | `degnorm(x[0] - 180)` |
| 28 | TRUE_REVATI | `",zePsc"` | `iflag_true` | `degnorm(x[0] - 359.8333333333)` |
| 29 | TRUE_PUSHYA | `",deCnc"` | `iflag_true` | `degnorm(x[0] - 106)` |
| 39 | TRUE_SHEORAN | `",deCnc"` | `iflag_true` | `degnorm(x[0] - 103.49264221625)` |
| 35 | TRUE_MULA | `",laSco"` | `iflag_true` | `degnorm(x[0] - 240)` |
| 17 | GALCENT_0SAG | `",SgrA*"` | `iflag_true` | `degnorm(x[0] - 240.0)` |
| 40 | GALCENT_COCHRANE | `",SgrA*"` | `iflag_true` | `degnorm(x[0] - 270.0)` |
| 30 | GALCENT_RGILBRAND | `",SgrA*"` | `iflag_true` | `degnorm(x[0] - 210.0 - 90.0*0.3819660113)` |
| 36 | GALCENT_MULA_WILHELM | `",SgrA*"` | `iflag_true \| SEFLG_EQUATORIAL` | `degnorm(swi_armc_to_mc(x[0], eps) - 246.6666666667)` |
| 31 | GALEQU_IAU1958 | `",GP1958"` | `iflag_galequ` | `degnorm(x[0] - 150)` |
| 32 | GALEQU_TRUE | `",GPol"` | `iflag_galequ` | `degnorm(x[0] - 150)` |
| 33 | GALEQU_MULA | `",GPol"` | `iflag_galequ` | `degnorm(x[0] - 150 - 6.6666666667)` |

`x[0]` is the ecliptic longitude in **degrees** returned by `swe_fixstar()`.

All 12 dispatch branches `return (retflag & SEFLG_EPHMASK)` immediately after computing `*daya`.

### GALCENT_MULA_WILHELM — RA-to-MC projection (sweph.c:3110–3121)

```
swe_fixstar(",SgrA*", tjd_et, iflag_true | SEFLG_EQUATORIAL, x, serr)
eps = swi_epsiln(tjd_et, iflag) * RADTODEG    // obliquity of date, in DEGREES
*daya = swi_armc_to_mc(x[0], eps)             // x[0] = RA in degrees
*daya = swe_degnorm(*daya - 246.6666666667)
```

`swi_armc_to_mc(armc, eps)` (swehouse.c:872–888): converts ARMC (RA of Midheaven, degrees) to
MC ecliptic longitude. Formula:
```
mc = atan2(tan(armc_rad), cos(eps_rad)) in degrees
if armc > 90 && armc <= 270: mc += 180   // quadrant correction
mc = degnorm(mc)
```
This uses `eps` in degrees, converted to radians internally. Note `swi_epsiln(tjd_et, iflag)` uses
the full `iflag` (not `iflag=0`) here — this is the ONLY ayanamsa case where `iflag` is passed to
`swi_epsiln` (see FP Fidelity §4 in `c-ref-ayanamsa.md`).

---

## Constants and Magic Numbers

| Name | Value | Source | Usage |
|---|---|---|---|
| `SWI_STAR_LENGTH` | 40 | sweph.h:772 | Max chars for starname/starbayer fields |
| `SE_MAX_STNAME` | 256 | swephexp.h:299 | Max input/output star name buffer (caller must provide 2× this) |
| `SE_STARFILE` | `"sefstars.txt"` | swephexp.h:387 | Primary catalog filename |
| `SE_STARFILE_OLD` | `"fixstars.cat"` | swephexp.h:386 | Legacy fallback catalog |
| `SEI_FILE_FIXSTAR` | 4 | sweph.h:177 | File slot index in `swed.fidat[]` |
| `J2000` | `2451545.0` | sweph.h:67 | Julian Day of J2000.0 |
| `B1950` | `2433282.42345905` | sweph.h:68 | Julian Day of B1950.0 |
| `PARSEC_TO_AUNIT` | `206264.8062471` | sweph.h:288 | 648000/π per IAU B2 2016 |
| `KM_S_TO_AU_CTY` | `21.095` | sweph.h:297 | km/s → AU/century |
| `PLAN_SPEED_INTV` | `0.0001` | sweph.h:299 | days; ×0.1 = 0.00001 for fixstar dt |
| `AUNIT` | `1.49597870700e+11` | sweph.h:273 | Meters per AU (DE431) |
| `DEGTORAD` | `π/180` | sweodef.h | Degree → radian |
| `RADTODEG` | `180/π` | sweodef.h | Radian → degree |
| `rdist (no parallax)` | `1000000000` AU | sweph.c:6468 | 1e9 AU = effectively infinite distance |

The proper-motion storage format in `struct fixed_star`:
- `ramot`, `demot`: RADIANS per Julian century (already divided by cos(dec) for RA)
- `radvel`: AU per Julian century
- `parall`: RADIANS (positive; 0 means unknown parallax)
- `epoch`: float; `0.0` = ICRS, `1950.0` = FK4 B1950, `2000.0` = FK5 J2000

---

## FP Fidelity Notes

### 1. dt for Earth/Sun positions

`dt = PLAN_SPEED_INTV * 0.1 = 0.0001 * 0.1 = 0.00001` days (sweph.c:6416). This is 10× smaller
than the planet pipeline's `PLAN_SPEED_INTV`. The speed of a fixed star is dominated by proper
motion, which is tiny, so a smaller `dt` suffices. Match exactly.

### 2. `swi_deflect_light(x, 0, ...)` — dt=0

Fixed stars use `dt = 0` as the second argument to `swi_deflect_light` (sweph.c:6562), unlike
planets which pass `dtsave_for_defl`. This means the deflection uses the Sun position from the
same epoch, with no retardation correction. Do not pass any non-zero dt here.

### 3. Speed comment "speed is incorrect"

The C source at line 6567 comments "speed is incorrect !!!" for aberration speed. This refers to
the fact that `swi_aberr_light_ex` uses `dt = PLAN_SPEED_INTV * 0.1` to numerically differentiate
the aberrated position, but the star's proper motion makes the two positions non-symmetrically
displaced. The error is small. Replicate the C behavior; do not attempt to fix it.

### 4. Proper motion over days (not centuries)

The `t * x[i+3]` application in step 7 uses `t` in DAYS and `x[i+3]` in AU/day. The struct
fields `ramot/demot` are in rad/century; they are divided by `36525.0` to get rad/day and stored in
`x[3..5]` before `swi_polcart_sp`. The final `x[i+3]` after polcart is in AU/day. Multiplying
by `t` (days) gives the displacement in AU. This is correct.

### 5. `ra_pm /= cos(de)` — applied once, at parse time

The catalog stores RA proper motion as the great-circle rate `μα × cos(δ)` (this is the
Hipparcos/ICRS convention). The code in `fixstar_cut_string` divides once by `cos(de)` to
recover the pure spherical RA rate `μα`. When the position is then converted to Cartesian via
`swi_polcart_sp`, the proper motion is treated as pure angular motion in (RA, Dec) space.
Do not apply this conversion again in the computation stage.

### 6. Epoch `0.0` (ICRS) treats all two decimal-point epochs as J2000

The C code's `else` branch in the epoch check (`epoch != 1950` → use J2000) means that
`epoch = 2000.0` and `epoch = 0.0` both use `t = tjd - J2000`. This is correct: ICRS is
defined with origin at J2000.0, and the epoch only differs in the FK4/FK5 correction step.

### 7. Sequential number lookup uses pre-sort array position

The sorted array has Bayer-key entries first (indices 0 to `n_fixstars_real - 1`) and named
entries after. Sequential-number lookup (`star_nr`) accesses `fixed_stars[star_nr - 1]`, which
is the (star_nr)th entry in the SORTED Bayer range. The sequential numbers therefore do NOT
correspond to file order after sorting — they are positions in the Bayer-sorted half.

### 8. Builtin `",SgrA*"` uses epoch `"2000"` not `"ICRS"`

The hardcoded SgrA* record is `"Gal. Center,SgrA*,2000,..."` (sweph.c:6783). The `atof("2000")`
gives epoch `2000.0`, triggering FK5→ICRS conversion via `swi_icrs2fk5` even though the SIMBAD
coordinates are already ICRS. The sefstars.txt version correctly uses `ICRS`. This means the
builtin SgrA* and the file version will give slightly different results. Replicate the C
behavior for the builtin exactly (keep epoch 2000.0).

### 9. `swi_polcart_sp` interprets x[1] as latitude, not declination

After setting `x[0]=ra, x[1]=de, x[2]=rdist`, `swi_polcart_sp` is called treating x[0] as
longitude and x[1] as latitude in spherical coordinates. RA and declination are used as if
they were ecliptic lon/lat. This is correct because the function is a generic (lon, lat, r) →
(X, Y, Z) converter and doesn't care about the semantics of the input frame.

---

## References

| Source | Content |
|---|---|
| sweph.c:6154–6174 | `fixstar_format_search_name` — input normalization |
| sweph.c:6178–6190 | `save_star_in_struct` — append to array |
| sweph.c:6193–6206 | `fixedstar_name_compare`, `fstar_node_compare` — sort/search comparators |
| sweph.c:6211–6306 | `fixstar_cut_string` — CSV parsing and unit conversion |
| sweph.c:6324–6395 | `load_all_fixed_stars` — catalog loader and array builder |
| sweph.c:6407–6669 | `fixstar_calc_from_struct` — position computation (THE CORE) |
| sweph.c:6674–6748 | `search_star_in_list` — lookup by name/number/bayer |
| sweph.c:6750–6803 | `get_builtin_star` — hardcoded fallback records |
| sweph.c:6818–6876 | `swe_fixstar2` — primary public entry |
| sweph.c:6878–6898 | `swe_fixstar2_ut` — UT-input variant |
| sweph.c:6911–6944 | `swe_fixstar2_mag` — magnitude lookup |
| sweph.c:3002–3142 | `swi_get_ayanamsa_ex` — fixed-star ayanamsa dispatch |
| sweph.h:773–779 | `struct fixed_star` definition |
| sweph.h:67–68 | `J2000`, `B1950` JD constants |
| sweph.h:257, 273, 288, 297, 299 | `AUNIT`, `PARSEC_TO_AUNIT`, `KM_S_TO_AU_CTY`, `PLAN_SPEED_INTV` |
| sweph.h:772 | `SWI_STAR_LENGTH` |
| sweph.h:843–846 | `n_fixstars_real`, `n_fixstars_named`, `n_fixstars_records` in `swe_data` |
| swephexp.h:299, 386–387 | `SE_MAX_STNAME`, `SE_STARFILE`, `SE_STARFILE_OLD` |
| ephe/sefstars.txt | Actual catalog file (4245 lines, ~2869 data records) |

---

## Frame Transforms: swi_icrs2fk5 and swi_FK4_FK5 (swephlib.c)

This section covers the two frame-transform functions called from `fixstar_calc_from_struct`
Step 4. Both are in `swephlib.c`. Neither is the same as `swi_bias` (`src/bias.rs`):

- `swi_bias` — ICRS↔GCRS **frame bias** (IAU 2000/2006 tiny rotation, already ported to Rust)
- `swi_icrs2fk5` — GCRS/ICRS↔**FK5** rotation (different, larger matrix, documented here)
- `swi_FK4_FK5` — FK4↔**FK5** RA coordinate correction (not a matrix rotation, documented here)

---

### `swi_icrs2fk5` (swephlib.c:2292–2333)

#### Signature

```c
void swi_icrs2fk5(double *x, int32 iflag, AS_BOOL backward)
```

Comment above the function: `/* GCRS to FK5 */`.

#### Rotation matrix

One 3×3 matrix `rb[3][3]` is hardcoded (swephlib.c:2302–2310). Indices are `rb[row][col]`:

```c
rb[0][0] = +0.9999999999999928;
rb[0][1] = +0.0000001110223287;
rb[0][2] = +0.0000000441180557;
rb[1][0] = -0.0000001110223330;
rb[1][1] = +0.9999999999999891;
rb[1][2] = +0.0000000964779176;
rb[2][0] = -0.0000000441180450;
rb[2][1] = -0.0000000964779225;
rb[2][2] = +0.9999999999999943;
```

| Constant | Row | Col | Value |
|---|---|---|---|
| `rb[0][0]` | 0 | 0 | `+0.9999999999999928` |
| `rb[0][1]` | 0 | 1 | `+0.0000001110223287` |
| `rb[0][2]` | 0 | 2 | `+0.0000000441180557` |
| `rb[1][0]` | 1 | 0 | `-0.0000001110223330` |
| `rb[1][1]` | 1 | 1 | `+0.9999999999999891` |
| `rb[1][2]` | 1 | 2 | `+0.0000000964779176` |
| `rb[2][0]` | 2 | 0 | `-0.0000000441180450` |
| `rb[2][1]` | 2 | 1 | `-0.0000000964779225` |
| `rb[2][2]` | 2 | 2 | `+0.9999999999999943` |

#### Computation

The matrix is orthogonal. `rb` itself rotates FK5→ICRS; `rb^T` (transpose) rotates ICRS→FK5.
The C uses two distinct access patterns:

**Forward (backward=FALSE) — ICRS → FK5** (swephlib.c:2322–2330):

```c
for (i = 0; i <= 2; i++) {
    xx[i] = x[0] * rb[0][i] + x[1] * rb[1][i] + x[2] * rb[2][i];
    if (iflag & SEFLG_SPEED)
        xx[i+3] = x[3] * rb[0][i] + x[4] * rb[1][i] + x[5] * rb[2][i];
}
```

This is `xx = rb^T · x` (column-major access = transpose multiplication).

**Backward (backward=TRUE) — FK5 → ICRS** (swephlib.c:2312–2320):

```c
for (i = 0; i <= 2; i++) {
    xx[i] = x[0] * rb[i][0] + x[1] * rb[i][1] + x[2] * rb[i][2];
    if (iflag & SEFLG_SPEED)
        xx[i+3] = x[3] * rb[i][0] + x[4] * rb[i][1] + x[5] * rb[i][2];
}
```

This is `xx = rb · x` (row-major access = standard matrix multiplication).

**Write-back** (swephlib.c:2332):

```c
for (i = 0; i <= 5; i++) x[i] = xx[i];
```

All six components are **always** written back, regardless of `SEFLG_SPEED`. If speed was not
requested (`!(iflag & SEFLG_SPEED)`), then `xx[3..5]` were never populated from the stack and
contain garbage. In `fixstar_calc_from_struct`, `iflag |= SEFLG_SPEED` is forced at the start
(Step 0), so this path never triggers in practice; do not try to reproduce the uninitialized
behaviour, but be aware of it. `swi_bias` does NOT have this issue — it conditionally copies
speed components.

#### Units

`x[0..2]` — equatorial Cartesian position in AU. `x[3..5]` — AU/day. The same matrix is applied
to both position and velocity (correct for a pure rotation with no time-varying component).

#### Call in `fixstar_calc_from_struct`

Step 4 (sweph.c:6488): `swi_icrs2fk5(x, iflag, TRUE)`

Called only when `epoch != 0.0` (FK5 J2000 or FK4 B1950 after precession). `backward=TRUE` →
applies `rb` → FK5→ICRS direction. The immediately following `swi_bias(x, J2000, SEFLG_SPEED,
FALSE)` then applies the tiny ICRS↔GCRS frame bias when `denum >= 403`.

For ICRS-epoch catalog data (`epoch == 0.0`), `swi_icrs2fk5` is **skipped**; `swi_bias` is
applied later at Step 10 instead.

---

### `swi_FK4_FK5` (swephlib.c:4098–4112) and `swi_FK5_FK4` (swephlib.c:4114–4123)

#### Signatures

```c
void swi_FK4_FK5(double *xp, double tjd)
void swi_FK5_FK4(double *xp, double tjd)
```

`swi_FK5_FK4` is the exact inverse; both are documented together here.

#### Algorithm — `swi_FK4_FK5`

The correction is applied in **spherical polar space**, not Cartesian. The function converts,
adjusts RA only, then converts back. Reference: Explanatory Supplement to the Astronomical
Almanac, p. 167f.

```
// 1. Guard: zero position → return unchanged
if xp[0] == 0 && xp[1] == 0 && xp[2] == 0: return

// 2. Detect whether speed is genuine
correct_speed = (xp[3] != 0)

// 3. Cartesian → polar+speed   (swephlib.c:4106)
swi_cartpol_sp(xp, xp)
// After: xp[0] = RA (radians), xp[1] = Dec (radians), xp[2] = radius (AU)
//        xp[3] = dRA/dt (rad/day), xp[4] = dDec/dt (rad/day), xp[5] = dr/dt (AU/day)

// 4. Apply RA correction   (swephlib.c:4108)
xp[0] += (0.035 + 0.085 * (tjd - B1950) / 36524.2198782) / 3600.0 * 15.0 * DEGTORAD

// 5. Apply RA speed correction (only when original speed was non-zero)   (swephlib.c:4109–4110)
if correct_speed:
    xp[3] += (0.085 / 36524.2198782) / 3600.0 * 15.0 * DEGTORAD

// 6. Polar+speed → Cartesian   (swephlib.c:4111)
swi_polcart_sp(xp, xp)
```

#### RA correction formula — units breakdown

The raw addend before DEGTORAD is `(0.035 + 0.085 * (tjd - B1950) / 36524.2198782) / 3600 * 15`.

| Term | Value | Unit |
|---|---|---|
| `0.035` | `0.035` | arcseconds of time (time-independent offset) |
| `0.085` | `0.085` | arcseconds of time per tropical century |
| `36524.2198782` | `36524.2198782` | days per tropical Julian century |
| `/ 3600` | — | arcseconds of time → degrees of time (hours) |
| `* 15` | — | degrees of time → degrees of arc |
| `* DEGTORAD` | `M_PI / 180.0` | degrees → radians |

Combined: adds `(0.035 + 0.085 * t_cty)` arcsec-time of RA, where `t_cty = (tjd - B1950) / 36524.2198782`
is the number of tropical Julian centuries since B1950.

Speed addend: `(0.085 / 36524.2198782) / 3600 * 15 * DEGTORAD` — the derivative of the
time-varying term with respect to Julian day (constant), in rad/day.

#### Constants table

| Symbol | Value | Source |
|---|---|---|
| `B1950` | `2433282.42345905` JD | sweph.h:68 ("1950 January 0.923") |
| `36524.2198782` | days/tropical-century | swephlib.c:4108 (comment: "trop. centuries") |
| constant RA offset | `0.035` arcsec-time | swephlib.c:4108 (Expl.Suppl. p. 167f) |
| RA rate | `0.085` arcsec-time/century | swephlib.c:4108 |
| `DEGTORAD` | `M_PI / 180.0` | sweodef.h:266 |

#### Algorithm — `swi_FK5_FK4` (inverse)

```
// 1. Guard
if xp[0] == 0 && xp[1] == 0 && xp[2] == 0: return

// 2. Cartesian → polar+speed
swi_cartpol_sp(xp, xp)

// 3. Subtract RA correction   (swephlib.c:4120)
xp[0] -= (0.035 + 0.085 * (tjd - B1950) / 36524.2198782) / 3600.0 * 15.0 * DEGTORAD

// 4. Subtract RA speed (ALWAYS — no correct_speed guard)   (swephlib.c:4121)
xp[3] -= (0.085 / 36524.2198782) / 3600.0 * 15.0 * DEGTORAD

// 5. Polar+speed → Cartesian
swi_polcart_sp(xp, xp)
```

Asymmetry: `swi_FK4_FK5` skips the speed correction when `xp[3] == 0`; `swi_FK5_FK4` always
subtracts it (no guard at swephlib.c:4121). Replicate this exactly.

#### E-terms of aberration

`swi_FK4_FK5` does **not** remove E-terms of aberration. It applies only the RA rotation
correction from the Explanatory Supplement. E-terms of aberration are a separate effect
(see Woolard & Clemence; FK4→FK5 sometimes includes E-term removal in other implementations).
The Swiss Ephemeris simplified FK4→FK5 path omits E-terms. Do not add them.

#### Helper functions

Both functions call `swi_cartpol_sp` (swephlib.c:362–413) and `swi_polcart_sp` (swephlib.c:420–451).
These convert between Cartesian `(X, Y, Z, dX/dt, dY/dt, dZ/dt)` and spherical polar
`(lon, lat, r, dlon/dt, dlat/dt, dr/dt)` where angles are in radians and distances in AU.
The angular speed components in the polar representation are the instantaneous rates
of change in radians per day (not per Julian century). The same helper functions appear
in Step 3 (`swi_polcart_sp`) and Step 16 (`swi_cartpol_sp`) of `fixstar_calc_from_struct`.

#### Call in `fixstar_calc_from_struct`

Step 4 (sweph.c:6484): `swi_FK4_FK5(x, B1950)`

Called only when `epoch == 1950.0`. At the call site, `x[0..5]` is the Cartesian
position+velocity vector at B1950 in the FK4 frame. After the call it is in FK5.
The next two calls precess the result to J2000 and then apply `swi_icrs2fk5`:

```
swi_FK4_FK5(x, B1950)                             // FK4 B1950 → FK5 B1950
swi_precess(x,   B1950, 0, J_TO_J2000)            // precess position  FK5 B1950 → J2000
swi_precess(x+3, B1950, 0, J_TO_J2000)            // precess velocity  FK5 B1950 → J2000
swi_icrs2fk5(x, iflag, TRUE)                      // FK5 J2000 → ICRS
if denum >= 403:
    swi_bias(x, J2000, SEFLG_SPEED, FALSE)        // ICRS → GCRS frame bias
```

`swi_FK5_FK4` is not called in `fixstar_calc_from_struct` but is the public inverse.

---

### Relationship to `swi_bias` / `frame_bias` (src/bias.rs)

| Function | Frame pair | Method | Matrix / coefficients |
|---|---|---|---|
| `swi_bias` | ICRS ↔ GCRS | 3×3 rotation | IAU 2000 or IAU 2006 matrix (two variants), swephlib.c:2229–2250 |
| `swi_icrs2fk5` | GCRS/ICRS ↔ FK5 | 3×3 rotation | `rb[3][3]` above, swephlib.c:2302–2310 |
| `swi_FK4_FK5` | FK4 ↔ FK5 | RA-only shift in polar space | 0.035 + 0.085×t_cty arcsec-time |

`frame_bias` in `src/bias.rs` is the Rust port of `swi_bias` only. The Rust porter must
implement `icrs2fk5` and `fk4_fk5` separately from `frame_bias`.

---

### FP Fidelity Notes for Frame Transforms

1. **`swi_icrs2fk5` multiply order**: the inner product is written as
   `x[0]*rb[j][i] + x[1]*rb[j][i] + ...` (left-to-right, no fused multiply-add). Rust's
   default IEEE 754 evaluation matches this as long as LLVM does not reorder across
   operator precedence boundaries.

2. **`swi_FK4_FK5` RA expression**: `(0.035 + 0.085 * (tjd - B1950) / 36524.2198782) / 3600 * 15 * DEGTORAD`.
   The division by `36524.2198782` is on the second term only; `0.035` is NOT divided by
   `36524.2198782`. Parentheses in the C force this: `(0.035 + <time-varying part>)`.
   The outer `/3600 * 15 * DEGTORAD` is then applied to the whole sum. Match exactly —
   especially the `/ 3600 * 15` ordering (left-to-right sequential, not a single `/ 240`).

3. **`swi_FK4_FK5` zero-position guard**: checks `xp[0] == 0 && xp[1] == 0 && xp[2] == 0`
   using exact equality (not epsilon). Replicate this — it is a guard against an all-zero
   Cartesian input, not a near-zero check.

4. **`swi_FK4_FK5` correct_speed gate**: checks `xp[3] == 0` exactly. A star with genuine
   zero x-velocity would have its speed correction suppressed. In practice this does not occur
   for real catalog entries (which always have non-zero proper motion after `swi_polcart_sp`).

---

### References

| Source | Content |
|---|---|
| swephlib.c:2292–2333 | `swi_icrs2fk5` — GCRS/ICRS↔FK5 rotation |
| swephlib.c:4098–4112 | `swi_FK4_FK5` — FK4→FK5 RA correction |
| swephlib.c:4114–4123 | `swi_FK5_FK4` — FK5→FK4 RA correction (inverse) |
| swephlib.c:362–413 | `swi_cartpol_sp` — Cartesian→polar+speed |
| swephlib.c:420–451 | `swi_polcart_sp` — polar+speed→Cartesian |
| swephlib.c:2205–2289 | `swi_bias` — ICRS↔GCRS frame bias (already in src/bias.rs) |
| sweph.h:68 | `B1950 = 2433282.42345905` |
| sweodef.h:266 | `DEGTORAD = M_PI / 180.0` |
