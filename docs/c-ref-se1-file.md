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

---

## sweph() — Segment Evaluation Pipeline (sweph.c:2125–2358)

`sweph()` is the core evaluation routine. It takes a Julian Day and internal body ID,
ensures the correct `.se1` file is open, loads the right Chebyshev segment, rotates
the packed coefficients to a usable frame, evaluates the polynomial, and returns
Cartesian position + velocity.

### Signature (sweph.c:2125)

```c
static int sweph(double tjd, int ipli, int ifno, int32 iflag,
                 double *xsunb, AS_BOOL do_save, double *xpret, char *serr)
```

| Parameter | Purpose |
|---|---|
| `tjd` | Julian Day to evaluate |
| `ipli` | Internal body ID (`SEI_*`) or asteroid/moon offset value |
| `ifno` | File slot index (`SEI_FILE_*`) |
| `iflag` | Computation flags (`SEFLG_SPEED`, `SEFLG_JPLEPH`, etc.) |
| `xsunb` | Barycentric Sun position (AU); used to convert heliocentric asteroids to barycentric; may be NULL |
| `do_save` | If true, write result into `pdp->x[]` and update `pdp->teval` |
| `xpret` | Output array [6] for position+velocity; may be NULL |
| `serr` | Error string buffer; may be NULL |

Returns `OK`, `ERR`, or `NOT_AVAILABLE`.

### Step 1 — Body ID Aliasing (sweph.c:2137–2142)

Asteroids and planetary moons share one `plan_data` slot (`SEI_ANYBODY`):

```c
ipl = ipli;
if (ipli > SE_AST_OFFSET)   ipl = SEI_ANYBODY;
if (ipli > SE_PLMOON_OFFSET) ipl = SEI_ANYBODY;
pdp = &swed.pldat[ipl];
```

`ipl` is used as the slot index throughout; `ipli` is kept for file lookup.

### Step 2 — Cache Hit Check (sweph.c:2150–2160)

```c
speedf1 = pdp->xflgs & SEFLG_SPEED;
speedf2 = iflag & SEFLG_SPEED;
if (tjd == pdp->teval
    && pdp->iephe == SEFLG_SWIEPH
    && (!speedf2 || speedf1)
    && ipl < SEI_ANYBODY) {
  for (i = 0; i <= 5; i++) xpret[i] = pdp->x[i];
  return OK;
}
```

Skipped for asteroids/planetary moons (`ipl == SEI_ANYBODY`): they always recompute.

### Step 3 — File Management (sweph.c:2164–2231)

If the open file's time range no longer covers `tjd`, or the body ID changed (asteroid),
close the file and free cached segment/refep:

```c
if (fdp->fptr != NULL) {
  if (tjd < fdp->tfstart || tjd > fdp->tfend
      || (ipl == SEI_ANYBODY && ipli != pdp->ibdy)) {
    fclose(fdp->fptr);  fdp->fptr = NULL;
    free(pdp->refep);   pdp->refep = NULL;
    free(pdp->segp);    pdp->segp = NULL;
  }
}
```

If no file is open: call `swi_gen_filename(tjd, ipli, fname)` to derive filename from
JD and body; open via `swi_fopen(ifno, fname, swed.ephepath, serr)`. On open failure:
- Planetary moons: retry without subdirectory prefix (e.g. `sat/`)
- Numbered asteroids: retry with `s` inserted before `.se1`; if still failing, try
  without the `astN/` subdirectory as well
- Otherwise: return `NOT_AVAILABLE`

After successful open: call `read_const(ifno, serr)` to parse the header.

### Step 4 — Segment Loading (sweph.c:2272–2283)

Load a new segment if the cache is empty or `tjd` is outside the cached range:

```c
if (pdp->segp == NULL || tjd < pdp->tseg0 || tjd > pdp->tseg1) {
  retc = get_new_segment(tjd, ipl, ifno, serr);
  if (retc != OK) return retc;
  if (pdp->iflg & SEI_FLG_ROTATE) {
    rot_back(ipl);
  } else {
    pdp->neval = pdp->ncoe;
  }
}
```

After `get_new_segment()`, if `SEI_FLG_ROTATE` is set: call `rot_back()` to transform
the Chebyshev coefficients from orbital-plane frame to the reference frame and set
`pdp->neval` (the last significant coefficient index). Otherwise set `neval = ncoe`
(use all coefficients).

### Step 5 — Time Normalisation (sweph.c:2285–2287)

Map `tjd` to the Chebyshev domain [-1, 1]:

```
t = (tjd - pdp->tseg0) / pdp->dseg   // ∈ [0, 1]
t = t * 2 - 1                          // ∈ [-1, 1]
```

No explicit boundary clamping is applied; the caller ensures `tjd` is in range.

### Step 6 — Chebyshev Evaluation (sweph.c:2294–2302)

```c
need_speed = (do_save || (iflag & SEFLG_SPEED));
for (i = 0; i <= 2; i++) {
  xp[i]   = swi_echeb(t, pdp->segp + i*pdp->ncoe, pdp->neval);
  if (need_speed)
    xp[i+3] = swi_edcheb(t, pdp->segp + i*pdp->ncoe, pdp->neval) / pdp->dseg * 2;
  else
    xp[i+3] = 0;
}
```

Velocity scaling: `swi_edcheb` returns the derivative with respect to the normalised
time `t`. Because `dt/d(tjd) = 2 / dseg`, the velocity in AU/day is:

```
xp[i+3] = (d/dt)[Chebyshev] * (2 / dseg)
```

`segp` layout after `get_new_segment()` + optional `rot_back()`:
- `segp[0 .. ncoe-1]`         — X (or right-ascension component) coefficients
- `segp[ncoe .. 2*ncoe-1]`    — Y coefficients
- `segp[2*ncoe .. 3*ncoe-1]`  — Z coefficients

### Step 7 — EMBHEL Special Case (sweph.c:2312–2331)

Current `.se1` files do not have a direct barycentric Sun record. Instead they store
heliocentric Earth. When `ipl == SEI_SUNBARY && (pdp->iflg & SEI_FLG_EMBHEL)`:

```c
// Force re-evaluation of EMB (don't use cached Earth value)
tsv = pedp->teval;
pedp->teval = 0;
retc = sweph(tjd, SEI_EMB, ifno, iflag | SEFLG_SPEED, NULL, NO_SAVE, xemb, serr);
pedp->teval = tsv;
// barycentric Sun = barycentric EMB - heliocentric Earth
for (i = 0; i <= 2; i++) xp[i] = xemb[i] - xp[i];
if (need_speed)
  for (i = 3; i <= 5; i++) xp[i] = xemb[i] - xp[i];
```

This is a recursive call to `sweph()` for `SEI_EMB` within the same function.

### Step 8 — Asteroid Heliocentric → Barycentric (sweph.c:2334–2343)

Asteroid positions in `.se1` files are heliocentric. When using JPL or SWISSEPH
ephemeris and `xsunb != NULL`:

```c
if (ipl >= SEI_ANYBODY) {
  for (i = 0; i <= 2; i++) xp[i] += xsunb[i];
  if (need_speed)
    for (i = 3; i <= 5; i++) xp[i] += xsunb[i];
}
```

### Step 9 — Save and Return (sweph.c:2345–2358)

```c
if (do_save) {
  pdp->teval = tjd;
  pdp->xflgs = -1;   // invalidate coordinate-system cache
  pdp->iephe = SEFLG_SWIEPH;  // (or psdp->iephe for asteroid files)
}
if (xpret != NULL)
  for (i = 0; i <= 5; i++) xpret[i] = xp[i];
return OK;
```

Output coordinate frame: for `SEI_FLG_ROTATE` bodies, `rot_back()` has already
converted to barycentric/heliocentric rectangular coordinates. For the Moon the frame
is equatorial J2000; for planets it is ecliptic J2000 (further transforms to date and
geocentric happen in the caller `app_pos_etc_plan`).

---

## get_new_segment() — Coefficient Unpacking (sweph.c:4367–4503)

Reads one Chebyshev segment from the open `.se1` file into `pdp->segp`.

### Signature (sweph.c:4367)

```c
static int get_new_segment(double tjd, int ipli, int ifno, char *serr)
```

Returns `OK` or `ERR`. On error: closes `fdp->fptr`, frees all planet data, returns `ERR`.

### Step 1 — Segment Number and Boundaries (sweph.c:4383–4387)

```c
iseg = (int32)((tjd - pdp->tfstart) / pdp->dseg);
pdp->tseg0 = pdp->tfstart + iseg * pdp->dseg;
pdp->tseg1 = pdp->tseg0 + pdp->dseg;
```

No bounds check on `iseg`; the caller ensures `tjd` is within `[pdp->tfstart, pdp->tfend]`.

### Step 2 — Locate Segment Data in File (sweph.c:4389–4393)

The segment index is a packed array of 3-byte (24-bit) absolute file offsets:

```c
fpos = pdp->lndx0 + iseg * 3;          // position of index entry
do_fread(&fpos, 3, 1, 4, fp, fpos, ...); // read 3 bytes → int32
fseek(fp, fpos, SEEK_SET);              // seek to coefficient data
```

### Step 3 — Allocate and Zero segp (sweph.c:4395–4397)

```c
if (pdp->segp == NULL)
  pdp->segp = malloc(pdp->ncoe * 3 * 8);   // 3 coords × ncoe × sizeof(double)
memset(pdp->segp, 0, pdp->ncoe * 3 * 8);
```

The zero-fill ensures trailing coefficients (for unpacked slots not present in the file)
default to 0.0.

### Step 4 — Per-Coordinate Header Parsing (sweph.c:4399–4425)

For each coordinate `icoord` in `{0=X, 1=Y, 2=Z}`:

```c
idbl = icoord * pdp->ncoe;   // write offset into segp[]
do_fread(c, 1, 2, 1, fp, SEI_CURR_FPOS, ...);  // read 2 header bytes

if (c[0] & 128) {   // bit 7 set → 6 precision levels, read 2 more bytes
  do_fread(c+2, 1, 2, 1, fp, SEI_CURR_FPOS, ...);
  nsizes = 6;
  nsize[0] = c[1] >> 4;   nsize[1] = c[1] & 0x0f;
  nsize[2] = c[2] >> 4;   nsize[3] = c[2] & 0x0f;
  nsize[4] = c[3] >> 4;   nsize[5] = c[3] & 0x0f;
} else {            // bit 7 clear → 4 precision levels, 2-byte header
  nsizes = 4;
  nsize[0] = c[0] >> 4;   nsize[1] = c[0] & 0x0f;
  nsize[2] = c[1] >> 4;   nsize[3] = c[1] & 0x0f;
}
nco = nsize[0] + ... + nsize[nsizes-1];
assert(nco <= pdp->ncoe);   // sanity check
```

`nsize[i]` = number of coefficients packed at precision level `i`.
`nco` = total coefficients written for this coordinate (may be less than `ncoe`; the
rest remain 0.0 from the memset).

### Step 5 — Unpacking Loop (sweph.c:4440–4494)

For each precision level `i`:

#### Levels 0–3: sign-magnitude integers (sweph.c:4443–4454)

```c
j = 4 - i;          // bytes per coefficient in file: 4, 3, 2, 1
k = nsize[i];       // number of coefficients at this level
do_fread(longs, j, k, 4, fp, SEI_CURR_FPOS, ...);
// Reads k items of j bytes each, sign-extended to uint32
for (m = 0; m < k; m++, idbl++) {
  if (longs[m] & 1)          // odd → negative
    segp[idbl] = -(((longs[m]+1) / 2) / 1e9 * rmax / 2);
  else                       // even → positive
    segp[idbl] = (longs[m] / 2) / 1e9 * rmax / 2;
}
```

Sign encoding: the integer is the sign-magnitude value doubled. Bit 0 is the sign bit.
Magnitude = `floor(v / 2)` for even, `floor((v+1) / 2)` for odd. Then scale by
`rmax / 2 / 1e9`.

#### Level 4: half-byte (4-bit) packing (sweph.c:4455–4473)

```c
k = (nsize[4] + 1) / 2;   // bytes to read (2 nibbles per byte)
do_fread(longs, 1, k, 4, fp, SEI_CURR_FPOS, ...);
// Two nibbles per byte: o=16 for high nibble, o=1 for low nibble
for (m = 0, j = 0; m < k && j < nsize[4]; m++) {
  for (n = 0, o = 16; n < 2 && j < nsize[4]; n++, j++, idbl++, longs[m] %= o, o /= 16) {
    if (longs[m] & o)
      segp[idbl] = -(((longs[m]+o) / o / 2) * rmax / 2 / 1e9);
    else
      segp[idbl] =  ((longs[m]    / o / 2) * rmax / 2 / 1e9);
  }
}
```

Each byte holds two 4-bit values. The sign bit is the `o`-bit of `longs[m]` before
the modulo strips the already-processed nibble. Magnitude is `floor((v & (o-1)) / 2)`
where `v` is the current byte value after previous nibbles are stripped.

#### Level 5: quarter-byte (2-bit) packing (sweph.c:4474–4493)

```c
k = (nsize[5] + 3) / 4;   // bytes to read (4 pairs per byte)
do_fread(longs, 1, k, 4, fp, SEI_CURR_FPOS, ...);
// Four 2-bit values per byte: o=64, 16, 4, 1
for (m = 0, j = 0; m < k && j < nsize[5]; m++) {
  for (n = 0, o = 64; n < 4 && j < nsize[5]; n++, j++, idbl++, longs[m] %= o, o /= 4) {
    if (longs[m] & o)
      segp[idbl] = -(((longs[m]+o) / o / 2) * rmax / 2 / 1e9);
    else
      segp[idbl] =  ((longs[m]    / o / 2) * rmax / 2 / 1e9);
  }
}
```

Same sign-magnitude pattern as level 4, but with 2-bit granularity. Four values per
byte, using `o = 64, 16, 4, 1` in successive iterations.

### Scaling Summary

All levels use the same final formula:

```
coefficient = magnitude_integer * (rmax / 2) / 1e9
```

where `rmax = pdp->rmax` (read from file; see `rmax` derivation in Per-Planet Metadata
section above). This maps the integer mantissa to AU (or radians for angular quantities).

### Result Layout in segp[]

After the three-coordinate loop:

```
segp[0    .. ncoe-1]      X Chebyshev coefficients (degrees, low to high)
segp[ncoe .. 2*ncoe-1]    Y Chebyshev coefficients
segp[2*ncoe .. 3*ncoe-1]  Z Chebyshev coefficients
```

Trailing coefficients for each coordinate (beyond `nco` for that coordinate) remain
0.0 from the memset.

---

## rot_back() — Orbital Plane to Reference Frame Rotation (sweph.c:4963–5054)

`rot_back()` is called immediately after `get_new_segment()` when `pdp->iflg & SEI_FLG_ROTATE`.
It transforms the Chebyshev coefficient arrays in `pdp->segp` **in place** from the
orbital-plane frame to the ecliptic (or equatorial, for the Moon) J2000 frame, and
optionally adds the reference ellipse.

### Constants (sweph.c:4976–4977)

```c
double seps2000 = 0.39777715572793088;   // sin(eps2000) — obliquity of J2000 ecliptic
double ceps2000 = 0.91748206215761929;   // cos(eps2000)
// eps2000 ≈ 0.409092804 radians ≈ 23°26'21"
```

### Step 1 — Segment Midpoint and Time Difference (sweph.c:4980–4984)

The orbital elements are interpolated at the **midpoint** of the segment, not at `tjd`:

```c
t = pdp->tseg0 + pdp->dseg / 2;
tdiff = (t - pdp->telem) / 365250.0;   // Julian millennia from element epoch
```

### Step 2 — Interpolated Equinoctal Elements (sweph.c:4985–4994)

**Moon** (ipli == SEI_MOON):

```c
dn = pdp->prot + tdiff * pdp->dprot;
dn -= floor(dn / TWOPI) * TWOPI;         // reduce mod 2π
qav = (pdp->qrot + tdiff * pdp->dqrot) * cos(dn);
pav = (pdp->qrot + tdiff * pdp->dqrot) * sin(dn);
```

**All other bodies**:

```c
qav = pdp->qrot + tdiff * pdp->dqrot;
pav = pdp->prot + pdp->dprot * tdiff;
```

Here `pav` and `qav` are the interpolated equinoctal inclination variables.

### Step 3 — Reference Ellipse Addition (sweph.c:5001–5013)

Only executed when `pdp->iflg & SEI_FLG_ELLIPSE`:

```c
omtild = pdp->peri + tdiff * pdp->dperi;
omtild -= floor(omtild / TWOPI) * TWOPI;   // reduce mod 2π
com = cos(omtild);
som = sin(omtild);
for (i = 0; i < nco; i++) {
  x[i][0] = chcfx[i] + com * refepx[i] - som * refepy[i];
  x[i][1] = chcfy[i] + com * refepy[i] + som * refepx[i];
  // x[i][2] unchanged (chcfz[i])
}
```

`refepx = pdp->refep` (first `ncoe` doubles), `refepy = pdp->refep + ncoe`.
`omtild` is the perihelion longitude interpolated from `peri` and `dperi`.
The rotation by `omtild` mixes the two ellipse coefficient arrays before adding
them to the Chebyshev coefficients. Only X and Y are modified; Z is left as-is.

Without `SEI_FLG_ELLIPSE`, `x[i][0..2]` is simply a copy of `chcfx/chcfy/chcfz`.

### Step 4 — Equinoctal Frame Construction (sweph.c:5019–5032)

From `pav` (p) and `qav` (q), three orthonormal basis vectors are computed:

```c
cosih2 = 1.0 / (1.0 + q*q + p*p);

// Orbit pole (normal to orbital plane):
uiz[0] =  2.0 * p * cosih2;
uiz[1] = -2.0 * q * cosih2;
uiz[2] =  (1.0 - q*q - p*p) * cosih2;

// Origin of longitudes (in-plane, reference direction):
uix[0] =  (1.0 + q*q - p*p) * cosih2;
uix[1] =  2.0 * q * p * cosih2;
uix[2] = -2.0 * p * cosih2;

// In-plane vector perpendicular to uix:
uiy[0] =  2.0 * q * p * cosih2;
uiy[1] =  (1.0 - q*q + p*p) * cosih2;
uiy[2] =  2.0 * q * cosih2;
```

These are the standard equinoctal element basis vectors. For planets, `p` and `q`
are referenced to the ecliptic J2000 frame; for the Moon, they are referenced to
the Moon's ecliptic orbital frame.

### Step 5 — Rotation and neval Update (sweph.c:5034–5048)

For each coefficient index `i` from 0 to `nco-1`:

```c
xrot = x[i][0]*uix[0] + x[i][1]*uiy[0] + x[i][2]*uiz[0];
yrot = x[i][0]*uix[1] + x[i][1]*uiy[1] + x[i][2]*uiz[1];
zrot = x[i][0]*uix[2] + x[i][1]*uiy[2] + x[i][2]*uiz[2];

if (fabs(xrot) + fabs(yrot) + fabs(zrot) >= 1e-14)
  pdp->neval = i;    // track last significant coefficient index (0-based)
```

The rotation matrix is `[uix | uiy | uiz]ᵀ`: each output component is the dot
product of the input vector `x[i]` with one of the basis vectors.

**Moon only** — additional ecliptic → equatorial J2000 rotation (sweph.c:5043–5047):

```c
if (ipli == SEI_MOON) {
  x[i][1] = ceps2000 * yrot - seps2000 * zrot;
  x[i][2] = seps2000 * yrot + ceps2000 * zrot;
  // x[i][0] = xrot (unchanged, x-axis is the same)
}
```

This is a rotation about the X axis by the obliquity of the ecliptic (eps2000):
```
| 1      0          0      |   | xrot |
| 0   cos(ε)   -sin(ε)    | × | yrot |
| 0   sin(ε)    cos(ε)    |   | zrot |
```

Non-Moon bodies: `x[i][0..2] = xrot, yrot, zrot` (no further rotation).

### Step 6 — Write Back (sweph.c:5049–5053)

```c
for (i = 0; i < nco; i++) {
  chcfx[i] = x[i][0];
  chcfy[i] = x[i][1];
  chcfz[i] = x[i][2];
}
```

`chcfx/chcfy/chcfz` are pointers into `pdp->segp` (same layout as after `get_new_segment()`),
so the rotation is performed in-place.

### neval semantics

`pdp->neval` after `rot_back()` holds the 0-based index of the last Chebyshev
coefficient with magnitude ≥ 1e-14 (in L1 norm of all three components after
rotation). `swi_echeb` and `swi_edcheb` are called with `neval` as `ncf`, which
evaluates indices 0..neval-1 — consistently truncating one term below the tracked
threshold, which is negligible for the precision goal.

In the non-ROTATE path: `pdp->neval = pdp->ncoe` (evaluate all coefficients).

---

## Chebyshev Evaluation — swi_echeb / swi_edcheb (swephlib.c:171–213)

Both functions use Clenshaw's recurrence (Roger Broucke's ACM algorithm 446, 1973).

### swi_echeb — Position (swephlib.c:171–185)

```c
double swi_echeb(double x, double *coef, int ncf)
{
  double x2 = x * 2.0;
  double br = 0., brp2 = 0., brpp = 0.;
  for (int j = ncf - 1; j >= 0; j--) {
    brp2 = brpp;
    brpp = br;
    br = x2 * brpp - brp2 + coef[j];
  }
  return (br - brp2) * 0.5;
}
```

Evaluates `coef[0..ncf-1]` at `x ∈ [-1,1]`. This is the standard backward Clenshaw
recurrence for `Σ cⱼ Tⱼ(x)` where `Tⱼ` are Chebyshev polynomials of the first kind.

### swi_edcheb — Derivative (swephlib.c:190–213)

```c
double swi_edcheb(double x, double *coef, int ncf)
{
  double x2 = x * 2.0;
  double xjp2 = 0., xjpl = 0., bjp2 = 0., bjpl = 0., bf = 0., bj = 0.;
  for (int j = ncf - 1; j >= 1; j--) {
    double dj = (double)(j + j);
    double xj = coef[j] * dj + xjp2;
    bj = x2 * bjpl - bjp2 + xj;
    bf = bjp2;
    bjp2 = bjpl;
    bjpl = bj;
    xjp2 = xjpl;
    xjpl = xj;
  }
  return (bj - bf) * 0.5;
}
```

Returns `d/dx [Σ cⱼ Tⱼ(x)]`. The result is the derivative with respect to the
normalised time `t ∈ [-1,1]`, not with respect to JD. The caller in `sweph()` applies
the chain rule:

```
velocity [AU/day] = swi_edcheb(t, coef, neval) * (2 / dseg)
```
