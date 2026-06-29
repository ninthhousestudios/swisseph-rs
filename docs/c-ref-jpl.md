# C Reference: JPL Binary Ephemeris Reader — swejpl.c / swejpl.h

Porting reference for the JPL DE binary ephemeris reader (`.bsp`-equivalent flat binary
format). Read this instead of the C source.

---

## Constants

### JPL Body Indices (swejpl.h:68–83)

These are the `J_*` indices used throughout `swejpl.c` as array slots in `pv[]` and
as `ntarg`/`ncent` arguments to `swi_pleph`. They are **not** the same as `SEI_*` or
`SE_*` constants.

| Constant | Value | Notes |
|---|---|---|
| `J_MERCURY` | 0 | |
| `J_VENUS` | 1 | |
| `J_EARTH` | 2 | In the file this is EMB (Earth-Moon Barycenter) |
| `J_MARS` | 3 | |
| `J_JUPITER` | 4 | |
| `J_SATURN` | 5 | |
| `J_URANUS` | 6 | |
| `J_NEPTUNE` | 7 | |
| `J_PLUTO` | 8 | |
| `J_MOON` | 9 | Geocentric Moon in the file |
| `J_SUN` | 10 | Barycentric Sun (SSB–Sun vector, same as pvsun) |
| `J_SBARY` | 11 | Solar System Barycenter — zero origin |
| `J_EMB` | 12 | Alias: pv[6×J_EMB] = pv[6×J_EARTH] after assembly |
| `J_NUT` | 13 | Nutations (longitude, obliquity) — 2 components |
| `J_LIB` | 14 | Librations — 3 components |

`pv[]` array size: 13 bodies × 6 doubles = 78 (`js->pv[78]`), swejpl.c:110.

### `pnoint2jpl` — SEI to JPL Index Map (sweph.h:311, sweph.c:180)

Maps the SwissEph internal body ID (`SEI_*`) to the JPL index for `swi_pleph` calls:

```c
#define PNOINT2JPL {J_EARTH, J_MOON, J_MERCURY, J_VENUS, J_MARS, \
                    J_JUPITER, J_SATURN, J_URANUS, J_NEPTUNE, J_PLUTO, J_SUN}
```

| SEI index | SEI name | JPL index |
|---|---|---|
| 0 | SEI_EMB / SEI_EARTH | J_EARTH (= EMB in file) |
| 1 | SEI_MOON | J_MOON |
| 2 | SEI_MERCURY | J_MERCURY |
| 3 | SEI_VENUS | J_VENUS |
| 4 | SEI_MARS | J_MARS |
| 5 | SEI_JUPITER | J_JUPITER |
| 6 | SEI_SATURN | J_SATURN |
| 7 | SEI_URANUS | J_URANUS |
| 8 | SEI_NEPTUNE | J_NEPTUNE |
| 9 | SEI_PLUTO | J_PLUTO |
| 10 | SEI_SUNBARY | J_SUN |

---

## Data Structures

### `struct jpl_save` (swejpl.c:97–111)

Single global instance, heap-allocated on first open, freed on close. Thread-local
(`TLS`) in threaded builds.

| Field | C Type | Purpose |
|---|---|---|
| `jplfname` | `char *` | Heap-allocated ephemeris filename |
| `jplfpath` | `char *` | Heap-allocated search path |
| `jplfptr` | `FILE *` | Open file handle; `NULL` when closed |
| `do_reorder` | `short` | 1 = byte-swap all reads; 0 = no swap |
| `eh_cval[400]` | `double[400]` | Named constant values from record 1 |
| `eh_ss[3]` | `double[3]` | `[start_jd, end_jd, segment_days]` |
| `eh_au` | `double` | AU in km |
| `eh_emrat` | `double` | Earth/Moon mass ratio |
| `eh_denum` | `int32` | DE number (e.g. 405) |
| `eh_ncon` | `int32` | Number of named constants |
| `eh_ipt[39]` | `int32[39]` | Body coefficient table: 13 bodies × 3 entries |
| `ch_cnam[2400]` | `char[2400]` | 400 constant names, 6 chars each |
| `pv[78]` | `double[78]` | 13 body state vectors (6 doubles each) |
| `pvsun[6]` | `double[6]` | Barycentric Sun state (populated on every call) |
| `buf[1500]` | `double[1500]` | Coefficient buffer for the currently-read record |
| `pc[18]` | `double[18]` | Chebyshev position polynomials T_0..T_{ncf-1} |
| `vc[18]` | `double[18]` | Chebyshev velocity polynomials |
| `ac[18]` | `double[18]` | Chebyshev acceleration polynomials |
| `jc[18]` | `double[18]` | Chebyshev jerk polynomials |
| `do_km` | `short` | Output units: 0 = AU + AU/day (default); 1 = km + km/sec |

`js` is a file-level `static TLS struct jpl_save *` (swejpl.c:113).

---

## File Format

### Overview

A JPL DE binary file has no ASCII header. It is a flat binary file consisting of
fixed-size records of `irecsz` bytes each, in native floating-point format (either
big-endian or little-endian; detected by value plausibility).

```
Record 0:   file header (title, constant names, ipt[] table, etc.)
Record 1:   constant values (cval[])
Record 2:   first data segment
Record 3:   second data segment
...
Record N+1: last data segment
```

`irecsz` = `ksize × 4` bytes (ksize is in 32-bit word units).
`ncoeffs` = `ksize / 2` (number of `double` values per record).

### Known ksize / record sizes

| DE versions | ksize (32-bit words) | irecsz (bytes) | ncoeffs (doubles) |
|---|---|---|---|
| 403, 405, 410, 413, 414, 418, 421 | 2036 | 8144 | 1018 |
| 404, 406 | 1456 | 5824 | 728 |
| 200 | 1652 | 6608 | 826 |
| 102 | 1652 (padded from 1546) | 6608 | 826 |

`ksize` is computed by `fsizer()` from the `ipt[]` table rather than hard-coded
(swejpl.c:275–292). DE102 files are padded to match DE200 size (swejpl.c:292–293).

### Record 0 — File Header

Read sequentially by `fsizer()` (swejpl.c:189–328) and again in full by `state()`
(swejpl.c:668–730). The `fread` calls happen against `js->jplfptr` starting from
offset 0.

```
Offset  Size (bytes)  Type       Field / Notes
------  ------------  --------   -----
0       252           char[252]  Title block: 3 × 84-byte ASCII strings
                                 (e.g. "JPL Planetary Ephemeris DE404/LE404 ...")
252     2400          char[2400] Constant names: 400 × 6-char null-padded strings
                                 stored in js->ch_cnam
2652    24            double[3]  ss[0]=start_jd, ss[1]=end_jd, ss[2]=seg_days
2676    4             int32      ncon = number of named constants
2680    8             double     au = AU in km
2688    8             double     emrat = Earth/Moon mass ratio
2696    144           int32[36]  ipt[0..35] — body coefficient table, first 12 bodies
2840    4             int32      numde = DE version number
2844    12            int32[3]   lpt[0..2] — libration table → copied to ipt[36..38]
```

After the above, `fsizer` rewinds; `state()` also reads record 1 at `offset = 1 × irecsz`:

```
Offset 1×irecsz   3200   double[400]   cval[0..399] = constant values
```

### Record 0 `ipt[]` Layout — Body Coefficient Table

`eh_ipt[39]` holds three entries per body. For body `i` (0-indexed):

| Field | Array index | Meaning |
|---|---|---|
| offset | `ipt[i*3+0]` | **1-based** index of first coefficient in `buf[]` |
| ncf | `ipt[i*3+1]` | Number of Chebyshev coefficients per component |
| na | `ipt[i*3+2]` | Number of sub-intervals per segment |

Body ordering in ipt[]:

| i | Body | JPL name | ncm (components) |
|---|---|---|---|
| 0 | 0 | Mercury | 3 |
| 1 | 1 | Venus | 3 |
| 2 | 2 | Earth-Moon Barycenter | 3 |
| 3 | 3 | Mars | 3 |
| 4 | 4 | Jupiter | 3 |
| 5 | 5 | Saturn | 3 |
| 6 | 6 | Uranus | 3 |
| 7 | 7 | Neptune | 3 |
| 8 | 8 | Pluto | 3 |
| 9 | 9 | Moon (geocentric) | 3 |
| 10 | 10 | SSBary Sun | 3 |
| 11 | 11 | Nutations | **2** (longitude, obliquity only) |
| 12 | 12 | Librations | 3 |

Libration table is stored separately as `lpt[3]` in record 0 and copied to `ipt[36..38]`
after reading (swejpl.c:269–271 in `fsizer`, swejpl.c:729–730 in `state`).

### DE403 Example Values (from swejpl.c:123–165 comment)

```
ipt[ 0, 1, 2] = [  3, 14, 4]  Mercury: starts buf[2],  14 coeff, 4 sub-intervals → 14×4×3=168
ipt[ 3, 4, 5] = [171, 10, 2]  Venus:   starts buf[170], 10 coeff, 2 sub-intervals → 10×2×3=60
ipt[ 6, 7, 8] = [231, 13, 2]  Earth:   starts buf[230], 13 coeff, 2 sub-intervals → 13×2×3=78
ipt[ 9,10,11] = [309, 11, 1]  Mars:    starts buf[308], 11 coeff, 1 sub-interval  → 11×1×3=33
ipt[12,13,14] = [342,  8, 1]  Jupiter: starts buf[341],  8 coeff, 1 sub-interval  → 8×1×3=24
ipt[15,16,17] = [366,  7, 1]  Saturn:  starts buf[365],  7 coeff, 1 sub-interval  → 7×1×3=21
ipt[18,19,20] = [387,  6, 1]  Uranus:  starts buf[386],  6 coeff, 1 sub-interval  → 6×1×3=18
ipt[21,22,23] = [405,  6, 1]  Neptune: starts buf[404],  6 coeff, 1 sub-interval  → 6×1×3=18
ipt[24,25,26] = [423,  6, 1]  Pluto:   starts buf[422],  6 coeff, 1 sub-interval  → 6×1×3=18
ipt[27,28,29] = [441, 13, 8]  Moon:    starts buf[440], 13 coeff, 8 sub-intervals → 13×8×3=312
ipt[30,31,32] = [753, 11, 2]  Sun:     starts buf[752], 11 coeff, 2 sub-intervals → 11×2×3=66
ipt[33,34,35] = [819, 10, 4]  Nut:     starts buf[818], 10 coeff, 4 sub-intervals → 10×4×2=80
ipt[36,37,38] = [899, 10, 4]  Lib:     starts buf[898], 10 coeff, 4 sub-intervals → 10×4×3=120
                                        last element at buf[1017]
```

Segment size for DE403: 32 days.

### Record 2+ — Data Records

Each data record:
```
buf[0]:  start JD of this segment
buf[1]:  end JD of this segment
buf[2..ncoeffs-1]: Chebyshev coefficients for all bodies
```

For body `i`, sub-interval `k` (0-based), component `c`, coefficient `j`:
```
buf_index = (ipt[i*3] - 1) + j + (c + k * ncm) * ncf
```

where `ipt[i*3]` is 1-based, so subtract 1 for 0-based indexing.

---

## Byte Order / Endianness

### Detection Strategy (swejpl.c:217–226)

Unlike SE1 files (which use a magic number), JPL files use **value plausibility**:

```c
// After fread of ss[3] without any reorder:
if (ss[2] < 1.0 || ss[2] > 200.0)   // segment length implausible
    js->do_reorder = TRUE;
else
    js->do_reorder = FALSE;
```

There is no separate "big-endian" / "little-endian" flag. `do_reorder` is a `short`
(0 or 1), not a bitfield. Validation bounds: `ss[0]` in [-5583942, 9025909],
`ss[1]` likewise, `ss[2]` in [1, 200] (swejpl.c:228–236).

### `reorder()` — In-Place Byte Reversal (swejpl.c:895–908)

```c
static void reorder(char *x, int size, int number) {
    char s[8];
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < size; j++)
            s[j] = x[size - j - 1];
        for (int j = 0; j < size; j++)
            x[j] = s[j];
        x += size;
    }
}
```

Reverses bytes of `number` items, each `size` bytes wide, in place.
Called immediately after each `fread` when `js->do_reorder` is set.
No size mismatch (file element size = memory element size, always).

Pattern throughout the code (swejpl.c:240–268 and similar):
```c
fread(&field, sizeof(T), count, js->jplfptr);
if (js->do_reorder)
    reorder((char *)&field, sizeof(T), count);
```

For the coefficient buffer in `state()` (swejpl.c:806–814), each double is read and
conditionally reordered individually in a loop.

---

## Record Selection and Time Normalisation (`state()`)

Source: swejpl.c:783–815.

### Epoch Decomposition (swejpl.c:783–797)

The input epoch `et` is split into an integer part (`et_mn`, always a half-integer JD)
and fractional part (`et_fr`):

```c
s = et - 0.5;
et_mn = floor(s);       // integer part of (et - 0.5)
et_fr = s - et_mn;      // fractional days since previous midnight
et_mn += 0.5;           // restore: et_mn is now the preceding midnight (half-integer JD)
```

This decomposes `et` as `et = et_mn + et_fr` where `et_mn` is a half-integer (noon-based
Julian Day) and `0 ≤ et_fr < 1`.

### Record Number (swejpl.c:794–797)

```c
nr = (int32)((et_mn - js->eh_ss[0]) / js->eh_ss[2]) + 2;
if (et_mn == js->eh_ss[1]) --nr;   // clamp: use last record for end epoch
```

`+ 2` because records 0 and 1 are the header and constants; data starts at record 2.
The end-of-range clamp prevents reading past the last segment.

### Normalised Time Within Segment (swejpl.c:797)

```c
t = (et_mn - ((nr - 2) * js->eh_ss[2] + js->eh_ss[0]) + et_fr) / js->eh_ss[2];
```

Simplified: `t = (et - seg_start) / seg_length`, where `seg_start = (nr-2)*ss[2] + ss[0]`.
Result `t ∈ [0, 1)`.

### Record Read (swejpl.c:799–815)

```c
if (nr != nrl) {           // nrl: last-read record number (static, cached)
    nrl = nr;
    fseeko(jplfptr, (off_t64)(nr * irecsz), SEEK_SET);
    for (k = 1; k <= ncoeffs; ++k) {
        fread(&buf[k-1], sizeof(double), 1, jplfptr);
        if (js->do_reorder) reorder(&buf[k-1], sizeof(double), 1);
    }
}
```

`nrl` is a `static TLS int32` initialized to 0. The read is skipped if the same record
was read on the previous call (single-segment cache).

---

## `interp()` — Chebyshev Evaluation (swejpl.c:472–591)

### Signature (swejpl.c:472–473)

```c
static int interp(double *buf, double t, double intv, int32 ncfin,
                  int32 ncmin, int32 nain, int32 ifl, double *pv)
```

| Parameter | Meaning |
|---|---|
| `buf` | Pointer into `js->buf` at first coefficient for this body (0-based: `&buf[ipt[i*3]-1]`) |
| `t` | Normalised segment time `∈ [0, 1)` |
| `intv` | Segment length: `ss[2]` (days) if `do_km=FALSE`; `ss[2]*86400` (seconds) if `do_km=TRUE` |
| `ncfin` | `ncf` — coefficients per component |
| `ncmin` | `ncm` — number of components (3 or 2 for nutations) |
| `nain` | `na` — sub-intervals per segment |
| `ifl` | 1=position only, 2=pos+vel, 3=pos+vel+acc, 4=pos+vel+acc+jerk |
| `pv` | Output: `pv[0..ncm-1]` positions, `pv[ncm..2*ncm-1]` velocities, etc. |

### Sub-Interval Selection (swejpl.c:497–504)

```c
// dt1 = floor(t), which is 0 for t ∈ [0,1) or 1 for t=1 (end-clamped)
dt1 = floor(t);           // handles negative t safely
temp = na * t;
ni = (int)(temp - dt1);   // sub-interval index, 0-based
// Normalised Chebyshev time within sub-interval, ∈ [-1, 1]:
tc = (fmod(temp, 1.0) + dt1) * 2.0 - 1.0;
```

`ni` selects which set of `ncf` coefficients (out of `na` sets) to use.
`tc` maps the sub-interval to the Chebyshev domain [-1, 1].

### Position Polynomial Recurrence (swejpl.c:511–533)

`pc[]` stores Chebyshev polynomials T_0(tc), T_1(tc), ..., T_{ncf-1}(tc).
Lazily extended: only recomputed when `tc` changes from previous call.

```c
// Initial values (set in swi_open_jpl_file):
pc[0] = 1.0;   // T_0(x) = 1
pc[1] = tc;    // T_1(x) = x  (updated whenever tc changes)
twot = tc * 2.0;

// Build up to ncf terms:
for (i = np; i < ncf; ++i)
    pc[i] = twot * pc[i-1] - pc[i-2];   // Chebyshev recurrence
np = ncf;
```

Position evaluation for component `i`:

```c
pv[i] = 0;
for (j = ncf-1; j >= 0; --j)
    pv[i] += pc[j] * buf[j + (i + ni*ncm)*ncf];
```

### Velocity Polynomial Recurrence (swejpl.c:540–553)

`vc[]` stores derivative polynomials (not the same as `pc[]`). The vc[] recurrence is:

```c
// Initial values (set in swi_open_jpl_file):
vc[1] = 1.0;

// Update at tc change:
vc[2] = twot + twot;   // = 4*tc (derivative of T_2 scaled)

// Build up:
for (i = nv; i < ncf; ++i)
    vc[i] = twot * vc[i-1] + 2.0*pc[i-1] - vc[i-2];
```

This is the forward Chebyshev derivative recurrence:
`U_i(x) = T'_i(x) / i` scaled such that `vc[j] = j * T'_j(tc) / (j)` — the exact
form that, when dotted with coefficients, gives `d/dtc(sum c_j T_j)`.

Velocity scaling factor:
```c
bma = (na + na) / intv;   // = 2*na / intv
```

This accounts for the chain rule: `d/dt_phys = (d/dtc) * (d_tc/d_t) * (d_t/d_t_phys)`.
Since `tc = 2*na*t - 1` (within sub-interval), `d_tc/d_t_phys = 2*na/intv`.

```c
pv[i+ncm] = 0;
for (j = ncf-1; j >= 1; --j)
    pv[i+ncm] += vc[j] * buf[j + (i + ni*ncm)*ncf];
pv[i+ncm] *= bma;
```

### Acceleration and Jerk (swejpl.c:556–590)

Similarly built with `ac[]` and `jc[]`:
```c
bma2 = bma * bma;
ac[3] = pc[1] * 24.0;   // initial value
for (i = nac; i < ncf; ++i)
    ac[i] = twot * ac[i-1] + 4.0*vc[i-1] - ac[i-2];

bma3 = bma * bma2;
jc[4] = pc[1] * 192.0;
for (i = njk; i < ncf; ++i)
    jc[i] = twot * jc[i-1] + 6.0*ac[i-1] - jc[i-2];
```

Scale factors: `bma2` for acceleration, `bma3` for jerk.

### Comparison with SE1 Chebyshev (`chebyshev_eval` / `chebyshev_deriv` in src/math.rs)

Both implement Chebyshev series evaluation over `[-1, 1]`, but in different styles:

| Aspect | JPL `interp()` | SE1 `swi_echeb` / `chebyshev_eval` |
|---|---|---|
| Position | Forward recurrence, dot product | Backward Clenshaw recurrence |
| Derivative | Forward `vc[]` recurrence, then dot | Backward Clenshaw on modified coefficients |
| Sub-intervals | `na` sub-intervals; `ni` selects which | Always 1 sub-interval per segment |
| Velocity scaling | `2*na/intv` | `2/dseg` |
| Input domain | `t ∈ [0,1)` (whole segment), mapped to `tc` | `t ∈ [-1,1]` (already mapped) |

For the Rust port, `interp()` needs its own implementation. The existing
`chebyshev_eval` / `chebyshev_deriv` in `src/math.rs` use the SE1 convention (takes
pre-normalised `t ∈ [-1,1]`). The JPL version adds sub-interval selection on top of
the same Chebyshev math. The polynomial basis is identical (`T_j` of the first kind);
only the organisation differs.

---

## `state()` — File Reading and Body Interpolation (swejpl.c:652–853)

### First-Call Initialisation (swejpl.c:668–780)

Triggered when `js->jplfptr == NULL`. Steps:

1. Call `fsizer()` to open the file, detect byte order, compute `ksize`
2. Re-read record 0 in full: title (252 bytes), cnam (2400 bytes), ss[3], ncon, au, emrat, ipt[0..35], numde, lpt[3] → ipt[36..38]
3. Read record 1 at `fseeko(fp, 1*irecsz, SEEK_SET)`: cval[400]
4. Validate file length: compute expected length from nseg and coefficient counts, compare against actual file size (swejpl.c:732–762)
5. Verify first/last segment start/end dates match header ss[0..1] (swejpl.c:765–779)
6. Set `nrl = 0` (no record cached yet)

If `list == NULL` on entry, return after initialisation (used by `read_const_jpl`).

### Interpolation Loop (swejpl.c:823–851)

After reading the correct record into `buf[]`:

```c
// Unit factor: if do_km=FALSE, aufac = 1/au (convert km-based file to AU)
aufac = js->do_km ? 1.0 : 1.0 / js->eh_au;
intv  = js->do_km ? js->eh_ss[2] * 86400.0 : js->eh_ss[2];

// Always compute SSBary Sun first (needed as heliocentric origin):
interp(&buf[ipt[30]-1], t, intv, ipt[31], 3, ipt[32], 2, pvsun);
for (i = 0; i < 6; ++i) pvsun[i] *= aufac;

// Planets + Moon (bodies 0..9):
for (i = 0; i < 10; ++i) {
    if (list[i] > 0) {
        interp(&buf[ipt[i*3]-1], t, intv, ipt[i*3+1], 3, ipt[i*3+2], list[i], &pv[i*6]);
        for (j = 0; j < 6; ++j) {
            if (i < 9 && !do_bary)
                pv[j + i*6] = pv[j + i*6]*aufac - pvsun[j];  // heliocentric
            else
                pv[j + i*6] *= aufac;                          // barycentric
        }
    }
}

// Nutations (2 components):
if (list[10] > 0 && ipt[34] > 0)
    interp(&buf[ipt[33]-1], t, intv, ipt[34], 2, ipt[35], list[10], nut);

// Librations (3 components):
if (list[11] > 0 && ipt[37] > 0)
    interp(&buf[ipt[36]-1], t, intv, ipt[37], 3, ipt[38], list[1], &pv[60]);
    // Note: list[1] in source is a C bug (should be list[11]), but librations
    // are rarely requested in practice.
```

Critical: `pv[6*J_EARTH..]` from `state()` is the **EMB** (Earth-Moon Barycenter),
not Earth. `pv[6*J_MOON..]` is the geocentric Moon. Conversion to barycentric Earth
and barycentric Moon happens in `swi_pleph()` after `state()` returns.

---

## `swi_pleph()` — Body Position Assembly (swejpl.c:362–449)

### Signature (swejpl.c:362)

```c
int swi_pleph(double et, int ntarg, int ncent, double *rrd, char *serr)
```

Returns position+velocity of `ntarg` relative to `ncent` in `rrd[6]`.

### list[] Setup (swejpl.c:374–416)

The `list[]` array tells `state()` which bodies to interpolate. Special dependencies:

```c
// Body needs Earth-Moon pair to derive each other:
if (ntarg == J_MOON)  list[J_EARTH] = 2;
if (ntarg == J_EARTH) list[J_MOON]  = 2;
if (ntarg == J_EMB)   list[J_EARTH] = 2;
// Same logic for ncent
```

Nutations and librations bypass this path (swejpl.c:375–398).

State is called with `do_bary = TRUE` for all regular bodies (swejpl.c:416):
```c
retc = state(et, list, TRUE, pv, pvsun, rrd, serr);
```

### Post-State Assembly (swejpl.c:418–447)

After `state()` returns, `pv[]` holds barycentric EMB (at J_EARTH slot) and
geocentric Moon (at J_MOON slot):

```c
// Sun slot: fill from pvsun
if (ntarg == J_SUN || ncent == J_SUN)
    for (i = 0; i < 6; ++i) pv[i + 6*J_SUN] = pvsun[i];

// SBARY slot: zero (it's the origin)
if (ntarg == J_SBARY || ncent == J_SBARY)
    for (i = 0; i < 6; ++i) pv[i + 6*J_SBARY] = 0.0;

// EMB slot: copy from Earth slot
if (ntarg == J_EMB || ncent == J_EMB)
    for (i = 0; i < 6; ++i) pv[i + 6*J_EMB] = pv[i + 6*J_EARTH];

// Earth/Moon conversion (swejpl.c:431–445):
if ((ntarg==J_EARTH && ncent==J_MOON) || (ntarg==J_MOON && ncent==J_EARTH)) {
    // Earth-Moon relative: zero the Earth slot; Moon is already geocentric.
    for (i = 0; i < 6; ++i) pv[i + 6*J_EARTH] = 0.0;
} else {
    // Convert EMB → barycentric Earth:
    if (list[J_EARTH] == 2)
        for (i = 0; i < 6; ++i)
            pv[i + 6*J_EARTH] -= pv[i + 6*J_MOON] / (js->eh_emrat + 1.0);
    // Convert geocentric Moon → barycentric Moon:
    if (list[J_MOON] == 2)
        for (i = 0; i < 6; ++i)
            pv[i + 6*J_MOON] += pv[i + 6*J_EARTH];
}

// Final: relative vector
for (i = 0; i < 6; ++i)
    rrd[i] = pv[i + ntarg*6] - pv[i + ncent*6];
```

### Earth Derivation Formula

```
Earth_bary = EMB_bary - Moon_geo / (emrat + 1)
Moon_bary  = Moon_geo + Earth_bary
```

where `emrat = eh_emrat` (e.g. 81.30056 for DE405).

---

## `fsizer()` — Record Size Computation (swejpl.c:189–328)

Called once from `state()` on first open. Steps:

1. Open file via `swi_fopen(SEI_FILE_PLANET, ...)` (swejpl.c:200)
2. `fread` 252 bytes title + 2400 bytes constant names
3. `fread` `ss[3]`; detect `do_reorder` by plausibility
4. `fread` ncon, au, emrat, `ipt[0..35]`, numde, `lpt[0..2]`; reorder each if needed
5. `rewind(jplfptr)`
6. Compute ksize from ipt[] (swejpl.c:275–291):
   ```c
   // Find body with highest starting offset:
   kmx = 0; khi = 0;
   for (i = 0; i < 13; i++) {
       if (ipt[i*3] > kmx) { kmx = ipt[i*3]; khi = i+1; }
   }
   nd = (khi == 12) ? 2 : 3;  // nutations have 2 components
   ksize = (ipt[khi*3-3] + nd * ipt[khi*3-2] * ipt[khi*3-1] - 1) * 2;
   ```
7. Special-case DE102: if `ksize == 1546`, set `ksize = 1652`
8. Validate: `1000 <= ksize <= 5000`
9. Return ksize

The caller (`state()`) computes `irecsz = 4 * ksize` and `ncoeffs = ksize / 2`.

---

## `read_const_jpl()` — Public Entry (swejpl.c:859–893)

```c
static int read_const_jpl(double *ss, char *serr) {
    int retc = state(0.0, NULL, FALSE, NULL, NULL, NULL, serr);
    if (retc != OK) return retc;
    for (int i = 0; i < 3; i++) ss[i] = js->eh_ss[i];
    return OK;
}
```

Triggers first-call initialisation in `state()` (by passing `list=NULL`), then copies
`ss[0..2]` (start JD, end JD, segment days) out to the caller.

---

## `swi_open_jpl_file()` and `swi_close_jpl_file()` (swejpl.c:924–958, 910–922)

### Open (swejpl.c:924–952)

```c
int swi_open_jpl_file(double *ss, char *fname, char *fpath, char *serr) {
    if (js != NULL && js->jplfptr != NULL) return OK;  // already open
    js = CALLOC(1, sizeof(struct jpl_save));
    js->jplfname = MALLOC(strlen(fname)+1); strcpy(js->jplfname, fname);
    js->jplfpath = MALLOC(strlen(fpath)+1); strcpy(js->jplfpath, fpath);
    retc = read_const_jpl(ss, serr);
    if (retc != OK) { swi_close_jpl_file(); return retc; }
    // Initialise polynomial arrays (lazy recurrence seeds):
    js->pc[0] = 1;   js->pc[1] = 2;  // pc[1]=2 is NOT tc; it's a sentinel so tc!=pc[1] on first call
    js->vc[1] = 1;
    js->ac[2] = 4;
    js->jc[3] = 24;
    return OK;
}
```

The initial `pc[1] = 2` is a dummy value (not a valid tc, which is in [-1,1]). It
ensures that the first real call to `interp()` always enters the `if (tc != pc[1])`
branch and initialises `twot`, `np`, `nv`, `nac`, `njk` correctly.

### Close (swejpl.c:910–922)

Closes file, frees `jplfname`, `jplfpath`, the `jpl_save` struct itself, sets `js = NULL`.

---

## Integration Contract — `jplplan()` in `sweph.c`

Source: sweph.c:1989–2106.

### What `jplplan()` Does

`jplplan()` is the only caller of `swi_pleph()` in `sweph.c`. It bridges the internal
SwissEph body ID space (`SEI_*`) to the JPL body space (`J_*`).

**Signature:**
```c
static int jplplan(double tjd, int ipli, int32 iflag, AS_BOOL do_save,
                   double *xpret, double *xperet, double *xpsret, char *serr)
```

| Parameter | Purpose |
|---|---|
| `ipli` | SEI body index (e.g. SEI_MERCURY=2, SEI_MOON=1) |
| `do_save` | Save result to `swed.pldat[ipli].x[]` and set `teval`, `xflgs`, `iephe` |
| `xpret` | Output `double[6]` for the requested body (may be NULL) |
| `xperet` | Output `double[6]` for Earth (always fetched; may be NULL) |
| `xpsret` | Output `double[6]` for barycentric Sun (always fetched; may be NULL) |

**Side effects:** `jplplan()` always fetches barycentric Earth and barycentric Sun as
a side effect, saving them in `swed.pldat[SEI_EARTH].x[]` and
`swed.pldat[SEI_SUNBARY].x[]` respectively (when `do_save=TRUE`).

### Body Mapping (sweph.c:2089)

```c
retc = swi_pleph(tjd, pnoint2jpl[ipli], ictr, xp, serr);
```

Center (`ictr`) selection:
- Moon: `ictr = J_EARTH` → geocentric Moon position
- All others: `ictr = J_SBARY` → barycentric position

Earth and Sun calls:
```c
swi_pleph(tjd, J_EARTH, J_SBARY, xpe, serr);   // barycentric Earth
swi_pleph(tjd, J_SUN,   J_SBARY, xps, serr);   // barycentric Sun
```

### Output Frame and Units

`jplplan()` places raw `swi_pleph()` output directly into `swed.pldat[ipli].x[6]`
with no frame conversion. The output is:

- **Frame**: barycentric equatorial J2000 / ICRF (for DE403+)
  (sweph.c:1677: "In barycentric equatorial position of the J2000 equinox")
- **Units**: AU (position) and AU/day (velocity); `do_km` is always FALSE in this context
- **Moon exception**: geocentric equatorial J2000, not barycentric

The comment in `state()` (swejpl.c:630–631) says "earth mean equator and equinox of
epoch" — this is the original Fortran heritage text. In practice, for DE403+, the actual
frame is ICRF (≈ J2000 equatorial). The downstream code handles the FK5/ICRF distinction
explicitly.

### Post-Processing Pipeline (after `jplplan` returns)

`main_planet()` (sweph.c:1562–1673) calls `jplplan()` then passes control to
`app_pos_etc_plan()` (sweph.c:2465+) or `app_pos_etc_sun()` for the remaining pipeline.
The pipeline applies to the raw barycentric equatorial J2000 positions in `pdp->x[]`:

1. **Light-time iteration**: iterative re-evaluation at `t - light_time`
2. **Frame bias** (DE403+ only, sweph.c:6488–6495): `swi_icrs2fk5()` backward (to ICRF),
   then `swi_bias()` for the ~22 mas J2000 → GCRS frame rotation
3. **Aberration**: stellar aberration correction
4. **Nutation**: IAU nutation in longitude and obliquity
5. **Precession**: from J2000 to date (ecliptic or equatorial)

The raw `pdp->x[]` (barycentric equatorial J2000) is what the Rust JPL reader must
produce. Everything downstream is shared with the SE1 pipeline.

---

## `swi_IERS_FK5()` (swejpl.h:102)

Declared in `swejpl.h` but not defined in the current codebase — it is a legacy stub
that was never implemented. Not used anywhere in `sweph.c`. Do not port.

---

## File Validation in `state()` (swejpl.c:732–779)

Two checks after first-call setup:

**1. File length** (swejpl.c:732–762):
```
nseg = (ss[1] - ss[0]) / ss[2]         // number of segments
expected_bytes = 0
for i in 0..12:
    k = if i==11 then 2 else 3          // nutations have 2 components
    expected_bytes += ipt[i*3+1] * ipt[i*3+2] * k * nseg  // coefficients
expected_bytes += 2 * nseg              // buf[0] and buf[1] (start/end JD) per segment
expected_bytes *= 8                     // doubles to bytes
expected_bytes += 2 * ksize * 4        // header + constants records
```
Accepts `flen == expected_bytes` or `flen == expected_bytes + ksize*4` (one extra record
at end, which some files have).

**2. Start/end date cross-check** (swejpl.c:765–779):
Reads `buf[0..1]` from the first data record (record 2) and from the last data record,
verifies they match `ss[0]` and `ss[1]` respectively.

---

## Summary: Key Differences from SE1 Backend

| Aspect | JPL backend | SE1 backend |
|---|---|---|
| Endian detection | Value plausibility of `ss[2]` | Magic number 0x616263 |
| Byte swap | `reorder()`: always reverses | `do_fread()`: handles size mismatches |
| Segment cache | One record cached by index (`nrl`) | Segment cached by JD range (`tseg0/tseg1`) |
| Coefficient storage | Raw doubles in `buf[]` per body/component/sub-interval | Packed variable-precision integers |
| Sub-intervals | Up to 8 per segment (e.g. Moon) | Always 1 per segment |
| Chebyshev style | Forward recurrence `pc[]`, dot product | Backward Clenshaw |
| Velocity scaling | `2*na/intv` applied in `interp()` | `2/dseg` applied in `sweph()` caller |
| Body assembly | Earth = EMB - Moon/(emrat+1) inside `swi_pleph` | Direct from segment, EMBHEL special case |
| Output frame | Barycentric equatorial J2000/ICRF | Barycentric ecliptic J2000 (planets) / equatorial J2000 (Moon) |
| Mutable global state | `static TLS struct jpl_save *js` | `swed.pldat[]`, `swed.fidat[]` in `swed` struct |
