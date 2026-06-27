# C Reference: SE1 File Format — sweph.c / sweph.h

Porting reference for `.se1` ephemeris file header parsing and metadata extraction.
Read this instead of the C source.

---

## Constants

### Internal Body IDs (sweph.h:131–152)

These are the `SEI_*` values used internally by the C library as indices into
`swed.pldat[]`. They are **not** the same as the public `SE_*` body numbers,
though many overlap in value.

| Constant | Value | Notes |
|---|---|---|
| `SEI_EPSILON` | -2 | Obliquity record |
| `SEI_NUTATION` | -1 | Nutation record |
| `SEI_EMB` | 0 | Earth-Moon barycenter |
| `SEI_EARTH` | 0 | Same slot as EMB/SUN |
| `SEI_SUN` | 0 | Same slot as EMB/EARTH |
| `SEI_MOON` | 1 | |
| `SEI_MERCURY` | 2 | |
| `SEI_VENUS` | 3 | |
| `SEI_MARS` | 4 | |
| `SEI_JUPITER` | 5 | |
| `SEI_SATURN` | 6 | |
| `SEI_URANUS` | 7 | |
| `SEI_NEPTUNE` | 8 | |
| `SEI_PLUTO` | 9 | |
| `SEI_SUNBARY` | 10 | Barycentric Sun |
| `SEI_ANYBODY` | 11 | Any asteroid |
| `SEI_CHIRON` | 12 | |
| `SEI_PHOLUS` | 13 | |
| `SEI_CERES` | 14 | |
| `SEI_PALLAS` | 15 | |
| `SEI_JUNO` | 16 | |
| `SEI_VESTA` | 17 | |
| `SEI_NPLANETS` | 18 | Size of `swed.pldat[]` array |

Node/apogee body IDs (stored in a **separate** array `swed.nddat[]`, not `swed.pldat[]`):

| Constant | Value |
|---|---|
| `SEI_MEAN_NODE` | 0 |
| `SEI_TRUE_NODE` | 1 |
| `SEI_MEAN_APOG` | 2 |
| `SEI_OSCU_APOG` | 3 |
| `SEI_INTP_APOG` | 4 |
| `SEI_INTP_PERG` | 5 |
| `SEI_NNODE_ETC` | 6 |

### Public Body Numbers (swephexp.h:101–128)

These are the `SE_*` values used in the public API (`swe_calc`). They partially
overlap with `SEI_*` values but are a distinct namespace.

| Constant | Value |
|---|---|
| `SE_SUN` | 0 |
| `SE_MOON` | 1 |
| `SE_MERCURY` | 2 |
| `SE_VENUS` | 3 |
| `SE_MARS` | 4 |
| `SE_JUPITER` | 5 |
| `SE_SATURN` | 6 |
| `SE_URANUS` | 7 |
| `SE_NEPTUNE` | 8 |
| `SE_PLUTO` | 9 |
| `SE_MEAN_NODE` | 10 |
| `SE_TRUE_NODE` | 11 |
| `SE_MEAN_APOG` | 12 |
| `SE_OSCU_APOG` | 13 |
| `SE_EARTH` | 14 |
| `SE_CHIRON` | 15 |
| `SE_PHOLUS` | 16 |
| `SE_CERES` | 17 |
| `SE_PALLAS` | 18 |
| `SE_JUNO` | 19 |
| `SE_VESTA` | 20 |
| `SE_INTP_APOG` | 21 |
| `SE_INTP_PERG` | 22 |
| `SE_NPLANETS` | 23 |
| `SE_PLMOON_OFFSET` | 9000 |
| `SE_AST_OFFSET` | 10000 |

### File Type Constants (sweph.h:173–178)

Used as `ifno` parameter — the index into `swed.fidat[]`.

| Constant | Value | Contents |
|---|---|---|
| `SEI_FILE_PLANET` | 0 | Sun/Moon/Planets (seplm*.se1, sepla*.se1) |
| `SEI_FILE_MOON` | 1 | Moon (semo*.se1) |
| `SEI_FILE_MAIN_AST` | 2 | Main asteroids Chiron, Pholus, Ceres, etc. |
| `SEI_FILE_ANY_AST` | 3 | Individual numbered asteroid |
| `SEI_FILE_FIXSTAR` | 4 | Fixed stars (different format, not `.se1`) |
| `SEI_FILE_PLMOON` | 5 | Planetary moons |

`SEI_NEPHFILES = 7` (sweph.h:194) — size of `swed.fidat[]` array (some slots reserved).

### Byte Order / Endianness Constants (sweph.h:183–187)

| Constant | Value | Meaning |
|---|---|---|
| `SEI_FILE_TEST_ENDIAN` | `0x616263L` | Magic value = "abc" bytes 0x61 0x62 0x63 |
| `SEI_FILE_BIGENDIAN` | 0 | File is big-endian |
| `SEI_FILE_NOREORD` | 0 | No byte reorder needed |
| `SEI_FILE_LITENDIAN` | 1 | File is little-endian |
| `SEI_FILE_REORD` | 2 | Byte reorder needed (file endian ≠ host endian) |

`fdp->iflg` is a bitfield: `freord | fendian`. Bits:
- bit 0 (`SEI_FILE_LITENDIAN`): 1 = file is little-endian, 0 = file is big-endian
- bit 1 (`SEI_FILE_REORD`): 1 = must reorder bytes when reading

### Other File Constants (sweph.h:189–195)

| Constant | Value | Meaning |
|---|---|---|
| `SEI_FILE_NMAXPLAN` | 50 | Max planets per file (size of `fdp->ipl[]`) |
| `SEI_FILE_EFPOSBEGIN` | 500 | (Reference offset, not used in read_const) |
| `SEI_CURR_FPOS` | -1 | Sentinel: use current file position (don't seek) |
| `SE_FILE_SUFFIX` | `"se1"` | File extension |

### Plan Data Flags (sweph.h:165–171)

Stored in `pdp->iflg` per planet; read as 1 byte from file.

| Constant | Value | Meaning |
|---|---|---|
| `SEI_FLG_HELIO` | 1 | Coordinates are heliocentric (vs. barycentric) |
| `SEI_FLG_ROTATE` | 2 | Chebyshev coefficients are in orbital-plane frame |
| `SEI_FLG_ELLIPSE` | 4 | Reference ellipse coefficients follow in file |
| `SEI_FLG_EMBHEL` | 8 | Heliocentric Earth given instead of barycentric Sun |

---

## Data Structures

### `struct file_data` (sweph.h:708–720)

One instance per open `.se1` file, stored in `swed.fidat[ifno]`.

| Field | C Type | Size | Purpose |
|---|---|---|---|
| `fnam` | `char[AS_MAXCH]` | variable | Full path to ephemeris file |
| `fversion` | `int` | 4 | Version number parsed from header line 1 |
| `astnam` | `char[50]` | 50 | Asteroid name (only for `SEI_FILE_ANY_AST`) |
| `sweph_denum` | `int32` | 4 | DE number of JPL ephemeris this derives from |
| `fptr` | `FILE *` | ptr | Open file handle |
| `tfstart` | `double` | 8 | File usable from this Julian Day |
| `tfend` | `double` | 8 | File usable through this Julian Day |
| `iflg` | `int32` | 4 | Byte order flags: `freord | fendian` |
| `npl` | `short` | 2 | Number of planets stored in file |
| `ipl` | `int[50]` | 200 | Array of planet body numbers (SE_* or SEI_*) |

`AS_MAXCH` is 256 in the C library.

### `struct plan_data` (sweph.h:610–656)

One instance per body, stored in `swed.pldat[SEI_*]`. Fields fall into three groups.

**Group 1: File metadata** — read once when file is opened via `read_const()`:

| Field | C Type | Source | Purpose |
|---|---|---|---|
| `ibdy` | `int` | assigned | Internal body number (`ipli` from file's `ipl[]`) |
| `iflg` | `int32` | file (1 byte) | `SEI_FLG_*` flags: helio/bary, rotate, ellipse |
| `ncoe` | `int` | file (1 byte) | Number of Chebyshev coefficients per segment (= order + 1) |
| `lndx0` | `int32` | file (4 bytes) | File offset of start of this planet's segment index |
| `nndx` | `int32` | computed | Number of index entries: `(tfend - tfstart + 0.1) / dseg` |
| `tfstart` | `double` | file (8 bytes) | Earliest JD covered (for this planet) |
| `tfend` | `double` | file (8 bytes) | Latest JD covered (for this planet) |
| `dseg` | `double` | file (8 bytes) | Days covered by one polynomial segment |
| `telem` | `double` | file (8 bytes) | Epoch of orbital elements |
| `prot` | `double` | file (8 bytes) | Orbital element: inclination vector component p |
| `dprot` | `double` | file (8 bytes) | Rate of change of `prot` |
| `qrot` | `double` | file (8 bytes) | Orbital element: inclination vector component q |
| `dqrot` | `double` | file (8 bytes) | Rate of change of `qrot` |
| `rmax` | `double` | computed | Normalisation factor for Chebyshev coefficients |
| `peri` | `double` | file (8 bytes) | Perihelion longitude (only if `SEI_FLG_ELLIPSE`) |
| `dperi` | `double` | file (8 bytes) | Rate of change of `peri` (only if `SEI_FLG_ELLIPSE`) |
| `refep` | `double *` | heap | Pointer to reference ellipse coefficients (2×ncoe doubles); `NULL` if no ellipse |

`rmax` derivation: read as `int32 lng` from file, then:
- Normal bodies: `pdp->rmax = lng / 1000.0`
- Planetary moon center-of-body: `pdp->rmax = lng / 1000000.0`
  (condition: `ipli >= SE_PLMOON_OFFSET` and `(ipli % 100) == 99` or `(ipli - 9000) / 100 == SE_MARS`)

`nndx` is not stored in the file; it is computed from `tfstart`, `tfend`, `dseg`.

**Group 2: Segment cache** — updated each time a new segment is loaded:

| Field | C Type | Purpose |
|---|---|---|
| `tseg0` | `double` | Start JD of cached segment |
| `tseg1` | `double` | End JD of cached segment |
| `segp` | `double *` | Pointer to unpacked Chebyshev coefficients (3×ncoe doubles for x,y,z) |
| `neval` | `int` | Coefficients to evaluate (may be < ncoe) |

**Group 3: Most recent evaluation result:**

| Field | C Type | Purpose |
|---|---|---|
| `teval` | `double` | JD for which `x[]` was computed |
| `iephe` | `int32` | Which ephemeris computed this |
| `x[6]` | `double` | Position + velocity, equatorial J2000 |
| `xflgs` | `int32` | Flags used in last computation |
| `xreturn[24]` | `double` | Results in 4 coordinate systems (6 doubles each) |

`xreturn` layout:
- `[0..5]`: ecliptic polar
- `[6..11]`: ecliptic cartesian
- `[12..17]`: equatorial polar
- `[18..23]`: equatorial cartesian

---

## File Format: Binary Layout

### Header (text section — always big-endian/ASCII)

The first portion of every `.se1` file is text with CRLF (`\r\n`) line endings.
Number of text lines depends on file type.

```
Line 1:  version string            e.g. "SE_EPHE_VERSION 2.00\r\n"
Line 2:  expected filename         e.g. "sepl_18.se1\r\n"
Line 3:  copyright line            e.g. "(c) Astrodienst AG...\r\n"
Line 4:  [asteroid file only]      MPC orbital elements line
```

The version is extracted from line 1 by scanning for the first run of digits.

### Binary Section

Immediately following the last text line (no padding):

```
Offset  Size  Type    Field
------  ----  ------  -----
+0      4     uint32  Endian test: value 0x616263 (big-endian) or
                      0x636261 (little-endian). Raw fread, no swap.
+4      4     int32   File length in bytes (sanity check)
+8      4     int32   JPL DE number (sweph_denum)
+12     8     double  tfstart  (file time range start, JD)
+20     8     double  tfend    (file time range end, JD)
+28     2     int16   nplan    (number of planets in file)
                      if nplan > 256: extended mode (nbytes_ipl=4, actual nplan = nplan % 256)
+30     2×N   int16[] ipl[nplan]  planet body numbers (2 bytes each, normal mode)
   OR   4×N   int32[] ipl[nplan]  planet body numbers (4 bytes each, extended mode)
```

After the planet ID array (if `SEI_FILE_ANY_AST`):
```
+?      30    char[30]  Asteroid name field (30 bytes, null-padded)
                        May be skipped/overwritten if name taken from elements line
```

Then for all files:
```
+?      4     uint32   CRC32 of all bytes from file start through here (exclusive)
+?      40    double×5 Physical constants: clight, aunit, helgravconst, ratme, sunradius
```

### Per-Planet Metadata Block

Repeated `npl` times (one per planet in file). For each planet:

```
Size  C type   Field         Notes
----  -------  -----         -----
4     int32    lndx0         File position of this planet's segment index
1     byte     iflg          SEI_FLG_* flags (stored in int32 with sign extension)
1     byte     ncoe          Number of Chebyshev coefficients
4     int32    rmax_raw      rmax * 1000 (or * 1000000 for some planetary moons)
80    double×10 orbital_data  [tfstart, tfend, dseg, telem, prot, dprot, qrot, dqrot, peri, dperi]
```

If `iflg & SEI_FLG_ELLIPSE`, immediately following:
```
ncoe×2×8  double[]  refep    Reference ellipse Chebyshev coefficients
                              First ncoe doubles = X component
                              Second ncoe doubles = Y component
```

All binary fields are subject to byte-order conversion via `do_fread()`.

---

## read_const() Algorithm (sweph.c:4510–4888)

Step-by-step walkthrough of header parsing. Called once per file, immediately
after the file is successfully opened.

### Step 1 — Version string (sweph.c:4535–4549)

```c
sp = fgets(s, AS_MAXCH, fp);   // read first text line (must contain \r\n)
fdp->fversion = atoi(sp);      // scan forward past non-digits, parse first integer
```

Error if `fgets` returns NULL or line lacks `\r\n`.

### Step 2 — Filename validation (sweph.c:4553–4582)

```c
sp = fgets(s, AS_MAXCH, fp);  // second line: stored expected filename
```

Extracts basename of `fdp->fnam` (after last DIR separator), converts both to
lowercase, compares. Error if mismatch.

### Step 3 — Copyright line (sweph.c:4586–4590)

```c
sp = fgets(s, AS_MAXCH, fp);  // third line: content discarded
```

Error only if `fgets` returns NULL or line lacks `\r\n`.

### Step 4 — Asteroid orbital elements (sweph.c:4594–4622, only if SEI_FILE_ANY_AST)

```c
sp = fgets(s, AS_MAXCH * 2, fp);  // fourth line: MPC elements record
```

Parses:
- `swed.astelem = s` (full line saved)
- `swed.ast_H = atof(s + 35 + i)` (absolute magnitude)
- `swed.ast_G = atof(s + 42 + i)` (slope parameter; defaults to 0.15 if 0)
- `swed.ast_diam = atof(s[51+i .. 51+i+7])` (diameter km; estimated from H if 0)

where `i` = number of leading characters before the asteroid name starts.

### Step 5 — Byte order detection (sweph.c:4626–4653)

```c
fread(&testendian, 4, 1, fp);          // raw read, no swap
if (testendian == SEI_FILE_TEST_ENDIAN) {
    freord = SEI_FILE_NOREORD;         // host and file agree
} else {
    freord = SEI_FILE_REORD;           // must swap bytes
    // byte-reverse testendian, verify == SEI_FILE_TEST_ENDIAN
}
// detect file endianness:
c = (char *)&testendian;
c2 = SEI_FILE_TEST_ENDIAN / 16777216L; // = 0x61 = 'a' (MSB of big-endian)
if (*c == c2)
    fendian = SEI_FILE_BIGENDIAN;      // first byte is MSB → big-endian file
else
    fendian = SEI_FILE_LITENDIAN;
fdp->iflg = (int32)freord | fendian;
```

The value `16777216 = 0x01000000 = 2^24`, so `0x616263 / 0x1000000 = 0x61`.

### Step 6 — File length check (sweph.c:4657–4670)

```c
do_fread(&lng, 4, 1, 4, fp, SEI_CURR_FPOS, ...);  // int32 from file
fpos = ftell(fp);                                   // save position for next read
fseek(fp, 0, SEEK_END);
flen = ftell(fp);
assert(lng == flen);                                // damage check 'h'
```

### Step 7 — DE number (sweph.c:4674–4677)

```c
do_fread(&fdp->sweph_denum, 4, 1, 4, fp, fpos, ...);  // seek to saved fpos
```

### Step 8 — File time range (sweph.c:4681–4688)

```c
do_fread(&fdp->tfstart, 8, 1, 8, fp, SEI_CURR_FPOS, ...);
do_fread(&fdp->tfend,   8, 1, 8, fp, SEI_CURR_FPOS, ...);
```

### Step 9 — Planet count (sweph.c:4692–4703)

```c
do_fread(&nplan, 2, 1, 2, fp, SEI_CURR_FPOS, ...);
if (nplan > 256) {
    nbytes_ipl = 4;       // extended format: 4-byte planet IDs
    nplan %= 256;
}
assert(1 <= nplan <= 20);
fdp->npl = nplan;
```

### Step 10 — Planet ID array (sweph.c:4705–4708)

```c
do_fread(fdp->ipl, nbytes_ipl, nplan, sizeof(int), fp, SEI_CURR_FPOS, ...);
```

Reads `nplan` integers, each `nbytes_ipl` bytes wide in the file, sign-extended
into `sizeof(int)` in memory. `sizeof(int)` on typical 64-bit Linux = 4 bytes.

### Step 11 — Asteroid name (sweph.c:4712–4753, only if SEI_FILE_ANY_AST)

Parses MPC number from `sastnam` (the name extracted from orbital elements line).
- If asteroid number matches `fdp->ipl[0] - SE_AST_OFFSET` or `fdp->ipl[0]` directly
  (planetary moon case): name is taken from elements record, and 30 bytes are read
  from file and discarded.
- Otherwise (older format): 30 bytes are read directly into `fdp->astnam`.

Name is then right-trimmed and double-space terminated.

### Step 12 — CRC32 check (sweph.c:4757–4780)

```c
fpos = ftell(fp);                              // position just before CRC field
do_fread(&ulng, 4, 1, 4, fp, SEI_CURR_FPOS, ...);  // read stored CRC
fseek(fp, 0, SEEK_SET);
fread(s, fpos, 1, fp);                         // re-read header bytes
assert(swi_crc32(s, fpos) == ulng);            // verify
fseek(fp, fpos + 4, SEEK_SET);                 // resume after CRC field
```

The CRC covers all bytes from start of file up to (not including) the CRC field itself.

### Step 13 — Physical constants (sweph.c:4786–4794)

```c
do_fread(doubles, 8, 5, 8, fp, SEI_CURR_FPOS, ...);  // 5 doubles
swed.gcdat.clight       = doubles[0];
swed.gcdat.aunit        = doubles[1];
swed.gcdat.helgravconst = doubles[2];
swed.gcdat.ratme        = doubles[3];
swed.gcdat.sunradius    = doubles[4];
```

These constants are written to the global state but are not currently used in
calculations.

### Step 14 — Per-planet loop (sweph.c:4798–4872)

For each planet `kpl` from 0 to `npl-1`:

```c
ipli = fdp->ipl[kpl];
if (ipli >= SE_AST_OFFSET)       pdp = &swed.pldat[SEI_ANYBODY];  // numbered asteroid
else if (ipli >= SE_PLMOON_OFFSET) pdp = &swed.pldat[SEI_ANYBODY];  // planetary moon
else                               pdp = &swed.pldat[ipli];          // normal planet
pdp->ibdy = ipli;
```

**a. lndx0** (4 bytes, int32):
```c
do_fread(&pdp->lndx0, 4, 1, 4, fp, SEI_CURR_FPOS, ...);
```

**b. iflg** (1 byte from file → stored as int32):
```c
do_fread(&pdp->iflg, 1, 1, sizeof(int32), fp, SEI_CURR_FPOS, ...);
```

**c. ncoe** (1 byte from file → stored as int):
```c
do_fread(&pdp->ncoe, 1, 1, sizeof(int), fp, SEI_CURR_FPOS, ...);
```

**d. rmax** (4 bytes → int32, then divided):
```c
do_fread(&lng, 4, 1, 4, fp, SEI_CURR_FPOS, ...);
pdp->rmax = lng / 1000.0;
// For planetary moon center-of-body records:
if (ipli >= SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET)
    if ((ipli % 100) == 99 || (ipli - 9000) / 100 == SE_MARS)
        pdp->rmax = lng / 1000000.0;
```

**e. 10 orbital doubles** (80 bytes):
```c
do_fread(doubles, 8, 10, 8, fp, SEI_CURR_FPOS, ...);
pdp->tfstart = doubles[0];
pdp->tfend   = doubles[1];
pdp->dseg    = doubles[2];
// nndx computed (not read):
pdp->nndx    = (int32)((doubles[1] - doubles[0] + 0.1) / doubles[2]);
pdp->telem   = doubles[3];
pdp->prot    = doubles[4];
pdp->dprot   = doubles[5];
pdp->qrot    = doubles[6];
pdp->dqrot   = doubles[7];
pdp->peri    = doubles[8];
pdp->dperi   = doubles[9];
```

**f. Reference ellipse** (only if `pdp->iflg & SEI_FLG_ELLIPSE`):
```c
pdp->refep = malloc(ncoe * 2 * 8);  // 2 × ncoe doubles
do_fread(pdp->refep, 8, 2 * ncoe, 8, fp, SEI_CURR_FPOS, ...);
// Layout: refep[0..ncoe-1] = X ellipse coefficients
//         refep[ncoe..2*ncoe-1] = Y ellipse coefficients
```

---

## Byte Order Handling

### `do_fread()` (sweph.c:4904–4955)

Signature:
```c
static int do_fread(void *trg, int size, int count, int corrsize,
                    FILE *fp, int32 fpos, int freord, int fendian, int ifno, char *serr)
```

- `size`: bytes per item **in the file**
- `count`: number of items
- `corrsize`: bytes per item **in memory** (may differ, e.g. 1-byte char → 4-byte int32)
- `fpos`: if `>= 0`, seek to this position first; if `== SEI_CURR_FPOS (-1)`, no seek

**Fast path** (no reorder, size == corrsize):
```c
fread(targ, totsize, 1, fp);
```

**Slow path** (reorder and/or size conversion):
1. Read all data into temp buffer `space[]`
2. For each item `i`, for each byte `j` from size-1 down to 0:
   - If `freord`: source byte `k = size - j - 1` (reversed)
   - Else: source byte `k = j` (same order)
   - If `size != corrsize`: adjust target position for big/little endian alignment
   - Copy: `targ[i*corrsize + k] = space[i*size + j]`

When size < corrsize (e.g. 1-byte to 4-byte):
- Zero-initialise the target first (`memset`)
- Alignment: if file is big-endian without reorder, or little-endian with reorder,
  shift source bytes toward high end: `k += corrsize - size`

---

## Segment Index and Segment Data (get_new_segment, sweph.c:4367–4503)

The segment index is a compact array of 3-byte (24-bit) file offsets.
To find segment `n` for a given JD:

```
n = (int32)((tjd - pdp->tfstart) / pdp->dseg)
index_entry_pos = pdp->lndx0 + n * 3
do_fread(&fpos, 3, 1, 4, fp, index_entry_pos, ...)  // read 3 bytes → int32
fseek(fp, fpos, SEEK_SET)
```

### Packed Chebyshev Coefficient Format

For each of 3 coordinates (x, y, z), coefficients are packed with variable precision.

**Header** (2 or 4 bytes):
```
Byte 0, bit 7 = 0:  4 size fields in 2 bytes
  nsize[0] = byte0 >> 4
  nsize[1] = byte0 & 0x0f
  nsize[2] = byte1 >> 4
  nsize[3] = byte1 & 0x0f

Byte 0, bit 7 = 1:  6 size fields in 4 bytes
  nsize[0] = byte1 >> 4
  nsize[1] = byte1 & 0x0f
  nsize[2] = byte2 >> 4
  nsize[3] = byte2 & 0x0f
  nsize[4] = byte3 >> 4
  nsize[5] = byte3 & 0x0f
```

`nsize[i]` = number of coefficients packed at precision level `i`.

**Precision levels** (i = 0..5):
| i | Bytes per coeff | Scale factor |
|---|---|---|
| 0 | 4 | `1e-9 * rmax / 2` |
| 1 | 3 | `1e-9 * rmax / 2` |
| 2 | 2 | `1e-9 * rmax / 2` |
| 3 | 1 | `1e-9 * rmax / 2` |
| 4 | half-byte (4 bits) | `rmax / 2 / 1e9` |
| 5 | quarter-byte (2 bits) | `rmax / 2 / 1e9` |

**Sign encoding** (levels 0–3): integers are sign-magnitude encoded.
Odd value → negative: `coeff = -((v + 1) / 2) / 1e9 * rmax / 2`
Even value → positive: `coeff = (v / 2) / 1e9 * rmax / 2`

**Half-byte (i=4)**: read `(nsize[4] + 1) / 2` bytes; two 4-bit values per byte.
Each nibble: if `(v & 16)` → negative, similar sign-magnitude.

**Quarter-byte (i=5)**: read `(nsize[5] + 3) / 4` bytes; four 2-bit values per byte.
Each pair: if `(v & 64)` → negative.

The 3×ncoe doubles are laid out as: all X coefficients, then all Y, then all Z.

---

## File Open / Selection (sweph.c:2180–2231)

Files are opened lazily on first use. The caller (`sweph_calc` or similar):
1. Calls `swi_gen_filename(tjd, ipli, fname)` to compute filename from JD + body
2. Calls `swi_fopen(ifno, fname, swed.ephepath, serr)` → returns `FILE *`
3. Calls `read_const(ifno, serr)` to parse the header

For numbered asteroids (`ipli > SE_AST_OFFSET`), if the normal file fails to open,
the code retries with an `s` inserted before the `.se1` suffix (short file variant).
For planetary moons (`ipli > SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET`), the code
retries without the subdirectory prefix (e.g. `sat/` for Saturn moons).

`swed.fidat[ifno].fptr == NULL` signals that no file is open for that slot.
