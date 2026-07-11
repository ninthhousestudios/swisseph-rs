# Swiss Ephemeris — Internal Architecture & Algorithms

> **Phase 2: Internal Implementation Catalogue**
> Companion to `catalogue-public.md` (Phase 1: Public API).
> Covers all `swi_*` internal functions, end-to-end calculation flow traces,
> algorithm details, data structures, and cross-domain dependencies.

---

## 1. Global State & Data Structures

### 1.1 The `swed` Global (struct swe_data)

All library state lives in a single thread-local global:

```c
extern TLS struct swe_data swed;
```

Key fields:

| Field | Type | Purpose |
|-------|------|---------|
| `pldat[SEI_NPLANETS]` | `plan_data[18]` | Cached planetary ephemeris data (one per internal body) |
| `nddat[SEI_NNODE_ETC]` | `node_data[6]` | Cached node/apside data |
| `savedat[SE_NPLANETS+1]` | `save_positions[]` | Final cached results keyed by (tjd, ipl, iflag) |
| `fidat[SEI_NEPHFILES]` | `file_data[7]` | Open ephemeris file descriptors |
| `oec`, `oec2000` | `struct epsilon` | Obliquity of ecliptic (date and J2000) |
| `nut`, `nut2000`, `nutv` | `struct nut` | Nutation (date, J2000, speed) |
| `topd` | `struct topo_data` | Observer geographic position |
| `sidd` | `struct sid_data` | Sidereal mode configuration |
| `interpol` | `struct interpol` | Nutation interpolation cache (3-point quadratic) |
| `*fixed_stars` | `fixed_star[]` | Dynamically loaded star catalog |
| `*dpsi`, `*deps` | `double[]` | IERS Earth Orientation Parameters (36525 entries) |
| `astro_models[8]` | `int32[]` | Selected precession/nutation/delta-T/sidereal-time/bias models |
| `tid_acc` | `double` | Lunar tidal acceleration (″/cy²) |
| `ephepath[]` | `char[256]` | Ephemeris data file path |
| `jplfnam[]` | `char[256]` | JPL file name |
| `last_epheflag` | `int32` | Cached ephemeris flag for invalidation |

### 1.2 struct plan_data — Planetary Ephemeris Cache

One per internal body. Stores both file metadata and the last evaluation result:

**File metadata (from ephemeris file header):**

| Field | Purpose |
|-------|---------|
| `ibdy` | Internal body number |
| `iflg` | Bit flags: `SEI_FLG_HELIO` (heliocentric), `SEI_FLG_ROTATE` (orbital plane coords), `SEI_FLG_ELLIPSE` (reference ellipse) |
| `ncoe` | Number of Chebyshev polynomial coefficients |
| `lndx0`, `nndx` | File position and count of index entries |
| `tfstart`, `tfend` | Julian Day range covered |
| `dseg` | Days per polynomial segment |
| `rmax` | Normalization factor for Chebyshev coefficients |
| `telem` | Epoch of orbital elements |
| `prot`, `qrot`, `dprot`, `dqrot` | Orbital plane rotation parameters and derivatives |
| `peri`, `dperi` | Perihelion parameters |
| `*refep` | Reference ellipse Chebyshev coefficients (2 × ncoe) |

**Current segment data:**

| Field | Purpose |
|-------|---------|
| `tseg0`, `tseg1` | Start/end JD of currently loaded polynomial segment |
| `*segp` | Unpacked Chebyshev coefficients (3 × ncoe doubles) |
| `neval` | Coefficients actually used in evaluation (≤ ncoe) |

**Last evaluation results:**

| Field | Purpose |
|-------|---------|
| `teval` | JD of last calculation (cache key) |
| `iephe` | Which ephemeris was used |
| `x[6]` | Position and velocity in equatorial J2000 |
| `xflgs` | Flags applied |
| `xreturn[24]` | Final output in 4 coordinate systems: ecliptic polar [0–5], ecliptic Cartesian [6–11], equatorial polar [12–17], equatorial Cartesian [18–23] |

### 1.3 Other Internal Structures

**struct epsilon** — Obliquity of ecliptic:
- `teps` (JD), `eps` (radians), `seps`/`ceps` (precomputed sin/cos)

**struct nut** — Nutation:
- `tnut` (JD), `nutlo[2]` (Δψ, Δε in radians), `snut`/`cnut` (sin/cos of Δε), `matrix[3][3]` (rotation matrix)

**struct interpol** — Nutation interpolation cache:
- Three-point quadratic: `tjd_nut0`/`tjd_nut2` bracket times, `nut_dpsi0`/`1`/`2`, `nut_deps0`/`1`/`2`

**struct topo_data** — Observer position:
- `geolon`, `geolat`, `geoalt`, `teval`, `tjd_ut`, `xobs[6]` (Cartesian observer in ecliptic coords)

**struct sid_data** — Sidereal mode:
- `sid_mode`, `ayan_t0`, `t0`, `t0_is_UT`

**struct file_data** — Open ephemeris file:
- `fnam[]`, `fversion`, `astnam[]`, `sweph_denum`, `*fptr`, `tfstart`, `tfend`, `iflg` (byte order), `npl`, `ipl[50]`

**struct fixed_star** — Star catalog record:
- `skey[]`, `starname[]`, `starbayer[]`, `starno[]`, `epoch`, `ra`, `de`, `ramot`, `demot`, `radvel`, `parall`, `mag`

**struct gen_const** — Physical constants:
- `clight`, `aunit`, `helgravconst`, `ratme`, `sunradius`

### 1.4 Key Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `AUNIT` | 1.49597870700×10¹¹ m | Astronomical Unit (DE431) |
| `CLIGHT` | 2.99792458×10⁸ m/s | Speed of light |
| `HELGRAVCONST` | 1.32712440018×10²⁰ m³/s² | GM_sun |
| `GEOGCONST` | 3.98600448×10¹⁴ m³/s² | GM_earth |
| `EARTH_RADIUS` | 6378136.6 m | Equatorial radius |
| `EARTH_OBLATENESS` | 1/298.25642 | Flattening |
| `MOON_SPEED_INTV` | 0.00005 days | Speed finite-difference interval (Moon) |
| `PLAN_SPEED_INTV` | 0.0001 days | Speed finite-difference interval (planets) |
| `NODE_CALC_INTV` | 0.0001 days | Node speed interval |
| `NODE_CALC_INTV_MOSH` | 0.1 days | Node interval for Moshier (nodes oscillate wildly) |
| `DEFL_SPEED_INTV` | 0.0000005 days | Deflection speed interval |
| `NUT_SPEED_INTV` | 0.0001 days | Nutation speed interval |

### 1.5 Platform Defines (sweodef.h)

| Define | Purpose |
|--------|---------|
| `TLS` | Thread-local storage: `__thread` (GCC), `__declspec(thread)` (MSVC), or empty |
| `int32` | `int` on 64-bit, `long` on 32-bit |
| `centisec` | `int32` — angles/times in centiseconds |
| `AS_MAXCH` | 256 — max string length |
| `DEGTORAD` | π/180 |
| `RADTODEG` | 180/π |
| `DIR_GLUE` | `"/"` (Unix), `"\\"` (Windows) |

---

## 2. Core Calculation Flow

### 2.1 swe_calc() → swecalc() Dispatch

**Entry point: `swe_calc()` (sweph.c:309)**

1. Validates and normalizes flags via `plaus_iflag()`
2. **Cache check**: if `savedat[ipl].tsave == tjd && ipl == sd->ipl && flags match` → return cached `xreturn[]`
3. **Speed mode selection**:
   - Default (SEFLG_SPEED): single call to `swecalc()` requesting speed
   - SEFLG_SPEED3: three calls at t−dt, t, t+dt with body-dependent interval; speed from central difference
4. Calls `swecalc()` for the actual computation
5. Copies result into `savedat[ipl]` for future cache hits

**Internal dispatcher: `swecalc()` (sweph.c:587)**

Validates flags, then dispatches by body number:

| Body ID | Branch | Functions Called |
|---------|--------|-----------------|
| SE_ECL_NUT (−1) | Direct | `swi_nutation()`, `swi_epsiln()` → returns obliquity + nutation |
| SE_MOON (1) | Moon | `jplplan()` → `sweplan()` → `swi_moshmoon()`, then `app_pos_etc_moon()` |
| SE_SUN (0) | Sun | `sweplan()` for barycentric sun, then `app_pos_etc_sun()` |
| Mercury–Pluto, Earth | Planets | `main_planet()` → backend + `app_pos_etc_plan()` |
| SE_MEAN_NODE (10) | Mean node | `swi_mean_node()` at two times, speed from difference |
| SE_MEAN_APOG (12) | Mean apogee | `swi_mean_apog()` similarly |
| SE_TRUE_NODE (11) | True node | `lunar_osc_elem()` with Newton iteration |
| SE_OSCU_APOG (13) | Osc. apogee | `lunar_osc_elem()` |
| SE_INTP_APOG (21) | Interp. apogee | `intp_apsides()` |
| 10000+ | Asteroids | `sweph()` with cascade fallback |
| SE_FICT_OFFSET+ | Fictitious | `swi_osc_el_plan()` from orbital elements |

### 2.2 Ephemeris Cascade Fallback

When the preferred ephemeris fails, the library automatically falls back:

```
Try preferred ephemeris (JPL, Swiss, or Moshier per flags):
  if NOT_AVAILABLE (-2):  → try next: JPL → Swiss → Moshier
  if BEYOND_EPH_LIMITS (-3) and date within Moshier range:  → switch to Moshier
  else:  → return error
```

The actual ephemeris used is recorded in the return flags so the caller can detect fallback.

### 2.3 main_planet() and sweplan()

**`main_planet()` (sweph.c:1562)**: Routes major planets through ephemeris selection. Calls `jplplan()` → `sweplan()` → `swi_moshplan()` with cascade. Then calls `app_pos_etc_plan()` for corrections.

**`main_planet_bary()` (sweph.c:1697)**: Returns raw barycentric positions without corrections. Used when multiple bodies are needed at the same time.

**`sweplan()` (sweph.c:1820)**: Orchestrates planet + moon + barycentric sun calculation. Handles caching for each body independently. Key conversions:
- Heliocentric → barycentric: adds Sun position
- EMB → Earth: `earth = emb − moon/(EARTH_MOON_MRAT + 1)`

**`jplplan()` (sweph.c:1989)**: Calls JPL kernel `swi_pleph()`. Returns barycentric equatorial J2000. Computes Earth and Sun as by-products when needed.

### 2.4 Apparent Position Pipeline

After raw coordinates are obtained from a backend, corrections are applied in sequence by `app_pos_etc_plan()` (sweph.c:2465) and `app_pos_rest()` (sweph.c:2777):

```
Raw barycentric equatorial J2000 position
  │
  ├─ 1. Light-time iteration  (Newton: t' = t − |Δx|/c, 1–2 iterations)
  ├─ 2. Aberration             (special-relativistic Lorentz transform)
  ├─ 3. Gravitational deflection (GR light bending by Sun)
  ├─ 4. Precession             (J2000 → equinox of date)
  ├─ 5. Nutation               (mean → true equinox of date)
  ├─ 6. Ecliptic transformation (equatorial → ecliptic via obliquity rotation)
  ├─ 7. Sidereal transformation (subtract ayanamsa, if SEFLG_SIDEREAL)
  ├─ 8. Polar conversion       (Cartesian → longitude, latitude, distance)
  └─ 9. Unit conversion        (radians → degrees, unless SEFLG_RADIANS)
```

Each correction step also adjusts the velocity components, typically via finite differences at a body-appropriate interval.

**The four output coordinate slots in `xreturn[24]`:**
- `[0–5]`: Ecliptic polar (lon, lat, dist, speed_lon, speed_lat, speed_dist)
- `[6–11]`: Ecliptic Cartesian (x, y, z, vx, vy, vz)
- `[12–17]`: Equatorial polar (RA, Dec, dist, ...)
- `[18–23]`: Equatorial Cartesian

The caller's `SEFLG_EQUATORIAL` / `SEFLG_XYZ` flags select which slot is copied to the output `xx[6]`.

---

## 3. Relativistic Corrections

### 3.1 Light-Time Iteration

**Location**: `app_pos_etc_plan()` (sweph.c:2545–2596)

Newton iteration to find retarded position:

```
For niter iterations (1 for JPL/Swiss, 0 for Moshier):
  1. dx[i] = planet[i] − observer[i]
  2. dt = |dx| × AUNIT / CLIGHT / 86400        (light-time in days)
  3. Apparent position: xx[i] = xx0[i] − dt × xx0[i+3]
     (true position minus velocity × light-time)
  4. Recompute dx with corrected position; iterate
```

Speed correction accounts for the changing light-delay:
- Three positions at t−1, t, t+1 yield d(dt)/dt
- Apparent speed = true speed × (1 − d(dt)/dt)

### 3.2 Annual Aberration: swi_aberr_light()

**Location**: sweph.c:3699

Special-relativistic aberration (Lorentz transformation):

```
v[i] = Earth_velocity[i] / (24 × 3600 × CLIGHT × AUNIT)    // in units of c
v² = v[0]² + v[1]² + v[2]²
β = √(1 − v²)                                                // Lorentz factor

f₁ = (u⃗ · v⃗) / |u⃗|
f₂ = 1 + f₁/(1 + β)

x_apparent[i] = (β × x[i] + f₂ × |u⃗| × v[i]) / (1 + f₁)
```

Speed correction uses finite differences at `PLAN_SPEED_INTV`.

### 3.3 Gravitational Deflection: swi_deflect_light()

**Location**: sweph.c:3743

GR light bending by the Sun:

```
U⃗ = planet geocentric       E⃗ = Earth heliocentric       Q⃗ = planet heliocentric
û, ê, q̂ = unit vectors      rE = |E⃗|

g₁ = 2 × GM_sun × m_eff / (c² × AUNIT × rE)
g₂ = 1 + (q̂ · ê)

x_deflected[i] = |U⃗| × (û[i] + (g₁/g₂) × (û·q̂ × ê[i] − û·ê × q̂[i]))
```

The **effective mass factor** `meff()` smooths the singularity near the solar disc by treating the Sun as an extended body. When the planet is near or behind the Sun, `meff` reduces the deflection to prevent divergence.

Speed correction uses `DEFL_SPEED_INTV` (0.0000005 days ≈ 43 ms).

---

## 4. Ephemeris Backends — Internal Algorithms

### 4.1 JPL DE Ephemerides (swejpl.c)

#### File Structure

JPL DE files contain fixed-size records of Chebyshev polynomial coefficients:

| Record | Contents |
|--------|----------|
| 0–1 | Headers: title, constant names, epoch range (`ss[0]`, `ss[1]`, `ss[2]`), AU, Earth-Moon mass ratio, body index table (`ipt[39]`), DE version number |
| 2+ | Data records: Chebyshev coefficients for all bodies over one segment (typically 32 days) |

**Body index table `ipt[i×3]`**: For each of 13 bodies (Mercury–Pluto, Moon, Sun, nutations, librations):
- `ipt[i×3+0]`: Starting buffer index for coefficients
- `ipt[i×3+1]`: Number of Chebyshev coefficients per component (`ncf`)
- `ipt[i×3+2]`: Number of sub-intervals within the segment (`na`)

Total coefficients per body = ncf × na × ncomponents (3 for position, 2 for nutations).

#### Segment Selection

```
et = Julian Ephemeris Date
nr = floor((et − ss[0]) / ss[2]) + 2       // record number (1-indexed, skip 2 header records)
t  = (et − (nr−2)×ss[2]) / ss[2]           // normalized time in [0, 1]
```

File seek: `fseek(fp, nr × record_size, SEEK_SET)` then read entire record into buffer.

#### Sub-Interval Selection

Each segment may be divided into `na` sub-intervals:

```
temp = na × t
ni   = (int)(temp)                           // which sub-interval [0, na)
tc   = 2 × frac(temp) − 1                   // Chebyshev parameter in [−1, 1]
```

#### Chebyshev Evaluation: interp()

**Position** — Clenshaw recurrence for T_n(tc):

```c
pc[0] = 1.0;  pc[1] = tc;
for (i = 2; i < ncf; i++)
    pc[i] = 2×tc × pc[i−1] − pc[i−2];

position[comp] = Σ pc[j] × coeff[j + (comp + ni×ncm)×ncf]
```

**Velocity** — derivative polynomials via modified Chebyshev:

```c
vc[0] = 0;  vc[1] = 1;  vc[2] = 4×tc;
for (i = 3; i < ncf; i++)
    vc[i] = 2×tc × vc[i−1] + 2×pc[i−1] − vc[i−2];

velocity[comp] = Σ vc[j] × coeff[...] × (2×na / intv)
```

**Acceleration** and **jerk** use further derivative recurrences (`ac[]`, `jc[]`), scaled by `(2na/intv)²` and `(2na/intv)³` respectively. The `ifl` parameter controls how many derivatives to compute (1=pos, 2=pos+vel, 3=+accel, 4=+jerk).

### 4.2 Swiss Ephemeris Files (sweephe4.c)

#### EP4 File Format

Files named `sep4_<N>` store daily ephemerides using **first and second differences** for compression (~10:1):

```c
struct ep4 {
    short j_10000, j_rest;           // Julian day encoding
    short ecl0m, ecl0s;             // True ecliptic day 0 (minutes, 0.01″)
    short ecld1[9];                  // First differences of ecliptic, days 1–9
    short nuts[10];                  // Nutation, days 0–9
    struct elon elo[14];             // Longitudes for 14 bodies
};

struct elon {
    short p0m, p0s;                  // Day 0 longitude (minutes, 0.01″)
    short pd1m, pd1s;                // First difference day 1 (minutes, 0.01″)
    short pd2[8];                    // Second differences days 2–9
};
```

Each file covers 10,000 days. Each record covers 10 days (NDB = 10).

#### Position Reconstruction from Differences

```
p₀ = p0m × 6000 + p0s               (full value in centiseconds)
d₁ = pd1m × 6000 + pd1s             (first difference)
p₁ = p₀ + d₁
For i = 2..9:
    dᵢ = dᵢ₋₁ + pd2[i−2] × scale   (scale = 10 for Moon/Mercury, 1 otherwise)
    pᵢ = pᵢ₋₁ + dᵢ
```

Wrapping handled at ±180° boundaries.

#### Everett 5th-Order Interpolation (Pottenger's Method)

The `inpolq()` function interpolates between tabulated daily values using 6 points (n−2 through n+3):

Given fractional day offset p ∈ [0, 1] and complementary q = 1−p:

**3rd-order Everett formula:**

```
result = q×x[n] + q₃×δ²₀ + p×x[n+1] + p₃×δ²₊₁

where:
    q₃ = (q+1)×q×(q−1)/6
    p₃ = (p+1)×p×(p−1)/6
    δ²₀ = x[n+1] − 2×x[n] + x[n−1]       (2nd central difference)
    δ²₊₁ = x[n+2] − 2×x[n+1] + x[n]
```

**5th-order extension** adds 4th-difference terms:

```
result += p₅×δ⁴₊₂ + q₅×δ⁴₊₁

where:
    p₅ = (p+2)×p₃×(p−2)/20
    q₅ = (q+2)×q₃×(q−2)/20
```

**Velocity** from derivative of interpolation formula:

```
deriv = d₀ + u×δ²₊₁ − u₀×δ²₀

where:
    u  = (3p² − 1)/6
    u₀ = (3q² − 1)/6
    d₀ = x[n+1] − x[n]              (first difference)
```

**Interpolation orders** per body: Planets use order 5; mean node, Chiron, Lilith use order 3.

#### Buffer Management

A static buffer `lon[14][20]` holds 20 days of data. Buffer reload triggers when the index falls below 2 or exceeds 16, with partial reload shifting the upper half down to avoid full re-reads.

### 4.3 Moshier Planetary Ephemeris (swemplan.c)

#### Fundamental Arguments (Simon et al. 1994)

Nine planetary mean longitudes computed as linear functions of time:

```c
freq[i] × T + phase[i]    (arcseconds, T in units of 10,000 Julian years)
```

| Index | Body | Frequency (″/10⁴ yr) | Phase (″) |
|-------|------|----------------------|-----------|
| 0 | Mercury | 53810162868.90 | 908127.69 |
| 1 | Venus | 21066413643.35 | 655127.28 |
| 2 | Earth | 12959774228.34 | 361679.24 |
| 3 | Mars | 6890507749.40 | 1279559.79 |
| 4 | Jupiter | 1092566037.80 | 123665.34 |
| 5 | Saturn | 439960985.54 | 180278.90 |
| 6 | Uranus | 154248119.39 | 1130597.56 |
| 7 | Neptune | 78655032.07 | 1095459.30 |
| 8 | Moon (anomaly) | 52272245.18 | 860492.15 |

#### Series Evaluation (swi_moshplan2)

**Step 1 — Build sin/cos lookup tables:**

For each of the 9 arguments, compute sin(k×L) and cos(k×L) for k = 0..max_harmonic using the recurrence:

```
sin((k+1)L) = sin(L)×cos(kL) + cos(L)×sin(kL)
cos((k+1)L) = cos(L)×cos(kL) − sin(L)×sin(kL)
```

**Step 2 — Sum harmonic series:**

The argument table (`arg_tbl`) encodes which fundamental arguments enter each term and with what multiplier. For each term:

1. Combine arguments: `W = Σ mᵢ × Lᵢ` using product-of-sums trig identity:
   ```
   sin(A+B) = sin(A)cos(B) + cos(A)sin(B)
   ```
2. Evaluate polynomial-in-T amplitude: `C(T) = c₀ + c₁T + c₂T² + ...`
3. Accumulate: `longitude += C_cos × cos(W) + C_sin × sin(W)`

**Step 3 — Convert:**

```
longitude = STR × sl        (STR = arcsec → radians)
latitude  = STR × sb
distance  = STR × plan->distance × sr + plan->distance
```

Each planet has ~1000–2000 terms. The `plantbl` struct holds pointers to coefficient tables in `swemptab.h`.

### 4.4 Moshier Lunar Ephemeris (swemmoon.c)

Based on Chapront-Touzé & Chapront's **ELP2000-85**, adjusted to fit DE404.

#### Four Fundamental Arguments

Computed in `mean_elements()`:

| Argument | Symbol | Description |
|----------|--------|-------------|
| l | MP | Mean anomaly of Moon |
| l' | M | Mean anomaly of Sun |
| D | D | Mean elongation |
| F | NF | Mean distance from ascending node |
| Lp | SWELP | Mean longitude of Moon |

Each is a high-precision polynomial in T (centuries from J2000) plus secular correction terms from 24 fitting coefficients (`z[]` array) that adjust for best agreement with DE404 from −3000 to +3000.

Additional planetary longitudes (Ve, Ea, Ma, Ju, Sa) enter perturbation terms.

#### Harmonic Term Tables

| Table | Terms | Content |
|-------|-------|---------|
| LR | 118 | Main longitude and radius (format: D, l', l, F, lon_coeff, rad_coeff) |
| MB | 77 | Main latitude |
| LRT | 38 | Linear-in-T longitude/radius corrections |
| BT | 16 | Linear-in-T latitude corrections |
| LRT2 | 25 | Quadratic-in-T longitude/radius corrections |
| BT2 | 12 | Quadratic-in-T latitude corrections |

Each term contributes: `coefficient × sin(n₁M + n₂D + n₃F + n₄Lp + corrections)`.

The computation is split across four functions (`moon1()` through `moon4()`) that accumulate into static thread-local variables `l` (longitude), `B` (latitude), and `moonpol[3]`.

**Speed**: computed by central differences at ±MOON_SPEED_INTV (0.00005 days ≈ 4.3 seconds).

### 4.5 Swiss Ephemeris Chebyshev Files (sweph.c)

#### Segment Loading: get_new_segment()

**Location**: sweph.c:4367

Segment selection:

```
iseg = (int)((tjd − tfstart) / dseg)
tseg0 = tfstart + iseg × dseg
tseg1 = tseg0 + dseg
fpos = lndx0 + iseg × 3               (3-byte file offset per segment)
```

**Coefficient unpacking** — four packing modes based on header byte:

| nsize | Format | Bytes per coefficient |
|-------|--------|----------------------|
| 0–3 | Integer packing | 1, 2, 3, or 4 bytes |
| 4 | Half-byte | 4-bit nibbles (2 per byte) |
| 5 | Quarter-byte | 2-bit values (4 per byte) |

Each coefficient is unpacked as: `coeff = ±((packed/2) / 10⁹ × rmax/2)`

Three coordinate sets (X, Y, Z) are packed separately within each segment.

#### Polynomial Evaluation: sweph()

**Location**: sweph.c:2285

```c
t = 2 × (tjd − tseg0) / dseg − 1         // normalize to [−1, 1]
position[i] = swi_echeb(t, segp + i×ncoe, neval)
velocity[i] = swi_edcheb(t, segp + i×ncoe, neval) / dseg × 2
```

Where `swi_echeb()` and `swi_edcheb()` use Clenshaw recurrence for Chebyshev evaluation and its derivative.

---

## 5. Precession, Nutation & Frame Bias Internals

### 5.1 Precession: swi_precess()

**Location**: swephlib.c:1373

**Router function** that selects the precession algorithm based on the configured model and time range:

| Time Range | Models Available | Function |
|------------|-----------------|----------|
| ±75–200 centuries | IAU 1976, IAU 2000, IAU 2006, Bretagnon 2003, Newcomb | `precess_1()` |
| Long-term | Laskar 1986, Simon 1994, Williams 1994 | `precess_2()` |
| Long-term | Vondrák 2011 (default), Owen 1990 | `precess_3()` |

#### precess_1() — Euler Angle Method

Computes three Euler angles (Z, z, Θ) as polynomials in T = (J − J2000)/36525:

```
For IAU 2006 (P03):
  Z  = (−3.173×10⁻⁷T − 5.971×10⁻⁶T + 0.01801828T + 0.2988499T + 2306.083227)T + 2.650545) / 3600°
  z  = (−2.904×10⁻⁷T − 2.8596×10⁻⁵T + 0.01826837T + 1.0927348T + 2306.077181)T − 2.650545) / 3600°
  Θ  = (−1.1274×10⁻⁷T − 7.089×10⁻⁶T − 0.04182264T − 0.4294934T + 2004.191903)T / 3600°
```

Applies rotation matrix: R_z(Z) × R_y(−Θ) × R_z(z)

#### precess_2() — Ecliptic-Based Long-Term

1. Computes obliquity at date
2. Rotates equatorial → ecliptic
3. Computes precession in longitude (pA), node, and inclination of moving ecliptic from polynomial expansions up to T¹⁰
4. Applies three sequential rotations on the ecliptic plane
5. Rotates back to equatorial

Uses 11-term polynomial coefficients from Laskar 1986, Williams 1994, or Simon 1994.

#### precess_3() — Vondrák 2011 (Default)

Builds precession matrix from equator and ecliptic pole positions:

1. `pre_pequ()`: equator pole (x, y) = polynomial + 14 periodic terms with periods 203–4043 years
2. `pre_pecl()`: ecliptic pole (x, y) = polynomial + 10 periodic terms
3. Equinox direction = cross product: equinox = peq × pecl (normalized)
4. Matrix columns = [equinox, perpendicular, equator_pole]

### 5.2 Nutation: swi_nutation()

**Location**: swephlib.c:2126

**Features**:
- Optional quadratic interpolation from 3 cached points (max error ~3 milliarcseconds)
- Delegates to `calc_nutation()` for full computation

#### IAU 1980 (Wahr): calc_nutation_iau1980()

Five fundamental Delaunay arguments:

```
MM = mean anomaly of Moon       (1717915922.633×T + 485866.733)/3600°
MS = mean anomaly of Sun        (129596581.224×T + 1287099.804)/3600°
FF = Moon's argument of latitude (1739527263.137×T + 335778.877)/3600°
DD = Moon's mean elongation      (1602961601.328×T + 1072261.307)/3600°
OM = Moon's ascending node       (−6962890.539×T + 450160.280)/3600°
```

Pre-computes sin/cos for multiples of each argument (up to 4×). Evaluates 105 terms:

```
Δψ = Σ (a + b×T) × sin(i×MM + j×MS + k×FF + l×DD + m×OM)
Δε = Σ (c + d×T) × cos(i×MM + j×MS + k×FF + l×DD + m×OM)
```

Optional Herring 1987 corrections add 7 additional terms.

#### IAU 2000A/2000B: calc_nutation_iau2000ab()

**Two-stage process:**

**Stage 1 — Luni-solar nutation** (678 terms for 2000A, 488 for 2000B):
- Uses Simon et al. 1994 Delaunay arguments with higher-order polynomial coefficients
- Series stored in `swenut2000a.h`

**Stage 2 — Planetary nutation** (687 terms, only for 2000A):
- Uses 14 fundamental arguments including planetary mean longitudes and accumulated general precession
- Adds very small corrections from planetary perturbations

**P03 adjustment**: Small corrections applied for consistency with IAU 2006 precession.

#### JPL Horizons Mode

When `SEFLG_JPLHOR` is set:
1. Base: IAU 1980 nutation series
2. Plus: Bessel interpolation of IERS Earth Orientation Parameters (`dpsi`, `deps` from `eop_1962_today.txt`)
3. Before 1962 or after EOP data: uses IAU 2000B with approximate corrections via `swi_approx_jplhor()`

### 5.3 Obliquity: swi_epsiln()

**Location**: swephlib.c:887

Computes mean obliquity of the ecliptic ε. Model-dependent polynomials in T:

| Model | Formula |
|-------|---------|
| IAU 1976 | ε = (1.813×10⁻³T − 5.9×10⁻⁴T − 46.8150)T + 84381.448″ |
| IAU 2006 | ε = (−4.34×10⁻⁸T⁵ − ... − 46.836769T + 84381.406″ |
| Laskar 1986 | 10th-degree polynomial in T/10 (decamillennial) |
| Vondrák 2011 | Via `swi_ldp_peps()`: polynomial + 10 periodic terms |

Returns ε in radians with precomputed sin(ε), cos(ε).

### 5.4 Frame Bias: swi_bias()

**Location**: swephlib.c:2205

Applies GCRS ↔ J2000 (FK5) offset matrix:

```
IAU 2000 bias matrix rb:
  [+0.9999999999999942  −7.078×10⁻⁸  +8.056×10⁻⁸]
  [+7.078×10⁻⁸          +0.9999999999999969  +3.306×10⁻⁸]
  [−8.056×10⁻⁸          −3.306×10⁻⁸  +0.9999999999999962]
```

Applies as matrix multiplication; `backward` flag transposes for inverse transformation.

---

## 6. Coordinate Transformations

### 6.1 Core Transform Functions

**`swi_coortrf(xpo, xpn, eps)`** (swephlib.c:279) — Rotation about x-axis by angle ε:

```
xpn[0] = xpo[0]
xpn[1] = xpo[1]×cos(ε) + xpo[2]×sin(ε)
xpn[2] = −xpo[1]×sin(ε) + xpo[2]×cos(ε)
```

Used for ecliptic ↔ equatorial: negative ε = ecliptic→equatorial, positive = equatorial→ecliptic.

**`swi_coortrf2(xpo, xpn, sineps, coseps)`** (swephlib.c:299) — Same but accepts precomputed sin/cos for efficiency.

**`swi_cartpol(x, l)`** (swephlib.c:314) — Cartesian → polar:

```
r   = √(x² + y² + z²)
rxy = √(x² + y²)
lon = atan2(y, x)        ∈ [0, 2π)
lat = atan(z / rxy)      ∈ [−π/2, π/2]
```

**`swi_polcart(l, x)`** (swephlib.c:343) — Polar → Cartesian:

```
x = r × cos(lat) × cos(lon)
y = r × cos(lat) × sin(lon)
z = r × sin(lat)
```

**`swi_cartpol_sp(x, l)`** (swephlib.c:362) — Cartesian → polar with velocity transformation via Jacobian.

**`swi_polcart_sp(l, x)`** (swephlib.c:420) — Inverse with velocity.

### 6.2 Coordinate Transform Sequence in app_pos_rest()

**Location**: sweph.c:2777

Applied after all physical corrections:

1. **Precession**: `swi_precess()` on equatorial J2000 Cartesian
2. **Nutation**: multiply by `nut.matrix[3][3]` to get true equinox of date
3. **Ecliptic transform**: `swi_coortrf2()` with sin/cos of true obliquity (ε + Δε)
4. **Sidereal** (if `SEFLG_SIDEREAL`):
   - Standard: subtract ayanamsa from ecliptic longitude
   - `SE_SIDBIT_ECL_T0`: project onto ecliptic of reference epoch via `swi_trop_ra2sid_lon()`
   - `SE_SIDBIT_SSY_PLANE`: project onto solar system invariable plane via `swi_trop_ra2sid_lon_sosy()`
5. **Polar conversion**: `swi_cartpol()` or `swi_cartpol_sp()`
6. **Degree conversion**: radians → degrees (unless `SEFLG_RADIANS`)

### 6.3 Horizontal Coordinates: swe_azalt()

Converts ecliptic or equatorial → horizontal (azimuth/altitude):

1. If ecliptic input: first convert to equatorial via `swe_cotrans()`
2. Compute local sidereal time → hour angle = LST − RA
3. Spherical trigonometry for altitude/azimuth:
   ```
   sin(alt) = sin(dec)×sin(lat) + cos(dec)×cos(lat)×cos(HA)
   ```
4. Apply atmospheric refraction via `swe_refrac_extended()`

---

## 7. Time Scale Internals

### 7.1 Julian Day: swe_julday()

**Location**: swedate.c:159

```
u = year;  if (month < 3) u -= 1
u₀ = u + 4712
u₁ = month + 1;  if (u₁ < 4) u₁ += 12

JD = floor(u₀ × 365.25) + floor(30.6 × u₁ + 0.000001) + day + hour/24 − 63.5

if (Gregorian):
    u₂ = floor(|u|/100) − floor(|u|/400)
    JD -= u₂ + 2
```

### 7.2 Delta-T: swe_deltat_ex()

**Location**: swephlib.c:2772

**Dispatcher logic**:

| Year Range | Method |
|------------|--------|
| Before model start (e.g., −720 for default) | Long-term parabola |
| Model-specific historical range | Model's formula or spline fit |
| 1620–present | Tabulated values with Bessel interpolation |
| After table end | Polynomial extrapolation with smooth transition |

#### Bessel Interpolation of Tabulated Delta-T

**Location**: swephlib.c:2004

4th-order Bessel (Stirling) formula for equidistant data:

```
Given array v[n] at 1-year spacing, fractional index t:

p = floor(t);  frac = t − p
ans = v[p]
ans += frac × (v[p+1] − v[p])                         // 1st order

Construct difference table:
d¹[i] = v[i+1] − v[i]
d²[i] = d¹[i+1] − d¹[i]
d³[i], d⁴[i] similarly

B = frac×(frac−1)/4
ans += B × (d²[1] + d²[2])                            // 2nd order

B = 2B/3
ans += (frac−0.5) × B × d³[1]                         // 3rd order

B = B × (frac+1) × (frac−2) / 8
ans += B × (d⁴[0] + d⁴[1])                            // 4th order
```

#### Tidal Acceleration Adjustment

Different DE ephemerides assume different lunar tidal acceleration values. Delta-T is adjusted:

```
correction = −0.000091 × (tid_acc − tid_acc_ref) × (Year − 1955)²
```

### 7.3 Sidereal Time: swe_sidtime0()

**Location**: swephlib.c:3464

Returns Greenwich Apparent Sidereal Time (hours):

1. Compute Earth Rotation Angle (ERA) from UT1
2. Add equation of equinoxes (33-term Fourier expansion evaluating planetary arguments)
3. For long-term model: uses Simon et al. mean Earth formula outside 1850–2050

The equation of equinoxes accounts for nutation's effect on the relationship between UT and sidereal time. Its 33 terms include arguments from lunar/solar Delaunay elements plus planetary longitudes.

### 7.4 Leap Second Handling

**Location**: swedate.c:273

Static table of 27 leap second dates (1972–2016) plus optional extension from `seleapsec.txt`.

`swe_utc_to_jd()`:
- Before 1972: input treated as UT1
- From 1972: validates against leap second table
- Allows `dsec ≥ 60.0` only on leap second dates
- TT = UTC + accumulated_leap_seconds + 32.184 seconds
- Returns both TT and UT1 Julian Days

---

## 8. Sidereal/Ayanamsa Internals

### 8.1 swi_get_ayanamsa_ex()

**Location**: sweph.c:3002

**Two computation methods** selected by `SE_SIDBIT_ECL_DATE`:

**Method 1 (traditional, default):**

```
1. Start with vernal point at J2000: x⃗ = (1, 0, 0)
2. Precess to observation epoch:     x⃗ → equatorial of date
3. Precess to reference epoch t₀:    x⃗ → equatorial of t₀
4. Convert to ecliptic of t₀:        apply ε(t₀)
5. Convert to polar:                  get longitude
6. ayanamsa = −longitude + initial_value
```

Measures how far the vernal point has moved on the ecliptic of t₀.

**Method 2 (ecliptic of date):**

```
1. Start with ayanamsa direction at t₀: x⃗ = (ayan_t₀, 0, 1)
2. Convert to equatorial of t₀:         remove ε(t₀)
3. Precess to J2000:                     x⃗ → J2000
4. Precess to observation epoch:         x⃗ → date
5. Apply ε(date):                        back to ecliptic
6. Convert to polar:                     get updated ayanamsa
```

More consistent: ayanamsa measured on ecliptic of date.

### 8.2 True Star Ayanamsas

For "true" ayanamsas (SE_SIDM_TRUE_CITRA, etc.):

1. Compute reference star position with `swe_fixstar()` at the observation epoch
2. Subtract the fixed ecliptic degree assigned to that star:
   - True Citra: Spica at exactly 180° (0° Libra)
   - True Revati: ζ Psc at 359°50′
   - True Pushya: δ Cnc at 106°
   - True Mula: λ Sco at 240°
3. The difference = current ayanamsa

### 8.3 Galactic Ayanamsas

For galactic-center-based ayanamsas:
- Uses galactic pole position (fixed star catalog)
- Requires `SEFLG_TRUEPOS` (no aberration/deflection) for pole accuracy
- Galactic equator measured from intersection of galactic pole great circle with ecliptic

### 8.4 Ayanamsa Speed

**Location**: sweph.c:3210

Two-point numerical derivative:

```
speed = (ayanamsa(t) − ayanamsa(t − 0.001)) / 0.001
```

Captures both precession and nutation contributions.

---

## 9. Eclipse Computation Internals

### 9.1 Solar Eclipse Search: swe_sol_eclipse_when_glob()

**Location**: swecl.c:1185

**Meeus Algorithm — Lunation Stepping:**

1. **Lunation number estimation**: K = (tjd − J2000) / 365.2425 × 12.3685
2. **F-argument filtering** (eliminates ~70% of lunations):
   ```
   F = 160.7108 + 390.67050274×K − 0.0016341×T² − ...
   F = F mod 180°
   if (F > 21° and F < 159°): skip — no eclipse possible
   ```
3. **Approximate maximum time** (Meeus formulas):
   ```
   tjd = 2451550.09765 + 29.530588853×K + ...
   M = Sun's mean anomaly;  Mm = Moon's mean anomaly;  E = eccentricity factor
   tjd -= 0.4075×sin(Mm) + 0.1721×E×sin(M)
   ```
4. **Iterative refinement**: quadratic interpolation via `find_maximum()` with progressively finer dt (starting at hours, converging to ~8 seconds)

### 9.2 eclipse_where() — Shadow Geometry

**Location**: swecl.c:640

**Shadow cone calculations:**

```
e⃗ = (moon − sun) / |moon − sun|        (shadow axis direction)
sinf₁ = (R_sun − R_moon) / d_sm        (umbra half-angle)
sinf₂ = (R_sun + R_moon) / d_sm        (penumbra half-angle)

s₀ = −(r⃗_moon · e⃗)                    (Moon distance from fundamental plane)
r₀ = √(d_m² − s₀²)                     (shadow axis distance from geocenter)

d₀ = |s₀/d_sm × (D_sun − D_moon) − D_moon| / cosf₁    (umbra diameter on fundamental plane)
D₀ = (s₀/d_sm × (D_sun + D_moon) + D_moon) / cosf₂    (penumbra diameter)
```

**Eclipse type detection:**
- Central: if Earth_radius × cosf₁ ≥ r₀
- Noncentral: if r₀ ≤ Earth_radius × cosf₁ + |d₀|/2
- Partial: if r₀ ≤ Earth_radius × cosf₂ + D₀/2

**Earth oblateness**: handled by scaling the z-coordinate: `z /= (1 − EARTH_OBLATENESS)`, with one refinement iteration.

### 9.3 eclipse_how() — Local Eclipse Circumstances

**Location**: swecl.c:967

Computes at a specific observer location using topocentric positions:

**Magnitude**: `attr[0] = (−dctr + rsun + rmoon) / (2 × rsun)`

Where dctr = center-to-center angular separation.

**Obscuration**: For partial eclipses, uses circular segment intersection:
```
sc₁ = area of Moon cap beyond Sun center
sc₂ = area of Sun cap beyond Moon center
attr[2] = (sc₁ + sc₂) × 2 / (π × rsun²)
```

For total: `attr[2] = 1.0`. For annular: `attr[2] = (rmoon/rsun)²`.

### 9.4 Lunar Eclipse: lun_eclipse_how()

**Location**: swecl.c:3237

Computes shadow cone from selenocentric perspective (Sun and Earth as seen from Moon):

```
f₁ = (R_sun − R_earth) / d_sm           (umbra half-angle)
f₂ = (R_sun + R_earth) / d_sm           (penumbra half-angle)
d₀ = umbra diameter × (1 + 1/50) × 0.99405     (atmospheric correction, NASA agreement)
D₀ = penumbra diameter × (1 + 1/50) × 0.98813
```

The factors 1/50 and 0.99405/0.98813 account for Earth's atmospheric enlargement of the shadow.

### 9.5 Contact Time Refinement

Three phases of contact detection:

| Phase | Condition | Physical Meaning |
|-------|-----------|------------------|
| Penumbra | D₀/2 + R_earth/cosf₂ − r₀ | First/last contact (P1/P4) |
| Umbra | |d₀|/2 + R_earth/cosf₁ − r₀ | Start/end of totality (U1/U4) |
| Central | R_earth/cosf₁ − r₀ | Center line begin/end |

Each uses `find_zero()` (quadratic root finding) with 2-hour initial sampling, then 3 refinement iterations narrowing to ~10 seconds.

---

## 10. Rise, Set & Transit Internals

### 10.1 Fast Algorithm: rise_set_fast()

**Location**: swecl.c:4203

Used when: Sun/Moon/planets at |lat| ≤ 60° (65° for Sun), no twilight.

1. **Semi-diurnal arc**: `sda = arccos(−tan(lat) × tan(dec))`
2. **Initial estimate**: from meridian distance and mean sidereal rate
3. **Iterative refinement** (2–4 iterations):
   ```
   Compute topocentric altitude at trial time
   dd = altitude rate of change (from 1.44-minute finite difference)
   dt = (altitude + disc_radius + refraction) / dd
   trial_time −= dt
   ```

### 10.2 Full Algorithm: swe_rise_trans_true_hor()

**Location**: swecl.c:4387

1. **Sample 14–15 points** at 2-hour intervals over ~28 hours
2. At each point, compute apparent altitude including:
   - Atmospheric refraction (`swe_refrac_extended()`)
   - Disc radius adjustment (upper limb, center, or lower limb)
   - Horizon height correction (geometric dip, custom horizon, twilight angles)
3. **Detect zero crossings**: sign changes in altitude array
4. **Binary search refinement**: 20 iterations converge to ~86 μs accuracy
5. **Culmination detection**: parabolic interpolation of altitude maxima/minima

### 10.3 Meridian Transit

**Location**: swecl.c:4688

```
For 4 iterations:
    armc = local sidereal time + geolon/15
    mdd  = RA − armc (mod 360°)
    t   += mdd / 361                        (361° not 360° to account for Earth's orbital motion)
    Recompute RA at new t
```

### 10.4 Circumpolar Detection

Objects that never rise or set are detected when:
- `sda = arccos(−tan(lat) × tan(dec))` has no solution (argument > 1 or < −1)
- Returns −2 (always above horizon) or +1 (always below)

---

## 11. Atmospheric Refraction

### 11.1 swe_refrac() / swe_refrac_extended()

**Location**: swecl.c:3035

**Saemundsson formula** (altitudes > 17.9°):

```
r = 0.97 / tan(alt)                        (arcminutes)
```

**Low-altitude polynomial** (−5° to 17.9°):

```
r = (34.46 + 4.23×alt + 0.004×alt²) / (1 + 0.505×alt + 0.0845×alt²)
```

**Atmospheric correction**: `r × (P − 80)/930 / (1 + 0.00008×(r+39)×(T−10)) / 60`

### 11.2 Horizon Dip for Elevated Observers

**calc_dip()** uses Thom's megalithic formula:

```
k_refr = (0.0342 + lapse_rate) / (0.154 × 0.0238)
d = 1 − 1.8480 × k_refr × P / (273.15 + T)²
dip = −arccos(1 / (1 + h/R_earth)) × √d
```

Where h = observer altitude, lapse_rate = temperature gradient (typically 0.0065 K/m).

---

## 12. Planetary Phenomena Internals

### 12.1 Phase Angle and Phase

**Location**: swecl.c:3791

```
phase_angle = arccos(geocentric_unit⃗ · heliocentric_unit⃗)
phase = (1 + cos(phase_angle)) / 2
elongation = arccos(planet_unit⃗ · sun_unit⃗)
apparent_diameter = 2 × arcsin(physical_diameter / (2 × AUNIT × distance))
```

### 12.2 Magnitude Formulas

| Body | Model |
|------|-------|
| Sun | Base −26.86, adjusted by (apparent/mean diameter)² |
| Moon | −21.62 + 0.026×|α| + 4×10⁻⁹×α⁴ (α ≤ 147°); transition to −4.5444 − 2.5×log₁₀((180−α)³) |
| Mercury, Venus | Mallama 2018: high-order polynomials in phase angle |
| Mars | Two-regime polynomial (α ≤ 50° and α > 50°) |
| Jupiter | Polynomial with phase angle |
| Saturn | Includes ring inclination parameter B |
| Uranus, Neptune | Quadratic in phase angle |
| Asteroids | Bowell HG system: ph₁ = e^(−3.33×tan(α/2)^0.63), ph₂ = e^(−1.87×tan(α/2)^1.22); V = H + 5×log₁₀(r×Δ) − 2.5×log₁₀((1−G)×ph₁ + G×ph₂) |

---

## 13. Nodes & Apsides Internals

### 13.1 Mean Elements

**Location**: swecl.c:5154

Static arrays `el_node[]`, `el_peri[]`, `el_incl[]`, `el_ecce[]`, `el_sema[]` contain 4-term polynomial coefficients for each planet's mean orbital elements.

```
element(T) = a₀ + a₁T + a₂T² + a₃T³      (T in centuries from J2000)
```

Perihelion direction is transformed from the orbital plane to the ecliptic:
```
ω_ecliptic = ω_orbit rotated by inclination around node direction
```

Node distance computed from the osculating ellipse using eccentric anomaly at the node crossing.

### 13.2 Osculating (True) Nodes and Apsides

**Location**: swecl.c:5238

**Angular momentum vector method:**

1. Compute 3 positions at t−dt, t, t+dt (dt scales with distance for outer bodies)
2. **Angular momentum**: h⃗ = r⃗ × v⃗ (cross product)
3. **Inclination**: sin(i) = |h_xy| / |h⃗|
4. **Node direction**: intersection of orbital plane (⊥ h⃗) with ecliptic (z = 0)
5. **Eccentricity and semi-major axis** from vis-viva:
   ```
   a = 1 / (2/r − v²/GM)
   p = h²/GM                     (semi-latus rectum)
   e = √(1 − p/a)
   ```
6. **True anomaly** from eccentric anomaly:
   ```
   cos(E) = (1 − r/a) / e
   sin(E) = (r⃗ · v⃗) / (e × √(a × GM))
   ν = 2 × arctan(√((1+e)/(1−e)) × sin(E) / (1 + cos(E)))
   ```
7. **Perihelion argument**: ω = u − ν (argument of latitude minus true anomaly)
8. **Perihelion distance**: a × (1 − e)

Node distances are recomputed using the osculating ellipse rather than raw intersection distances (which can be unrealistically large).

Speeds from central differences: `v = (pos(t+dt) − pos(t−dt)) / (2×dt)`

### 13.3 Lunar Osculating Elements

**Location**: sweph.c:5168 (`lunar_osc_elem()`)

Special handling for the Moon:
- Computes 3 positions with light-time correction
- Uses `GEOGCONST × (1 + 1/EARTH_MOON_MRAT)` for GM
- Moshier backend requires larger dt (NODE_CALC_INTV_MOSH = 0.1 days) because osculating lunar nodes "oscillate wildly"
- Node vector computed by projecting Moon position onto ecliptic (finding where tangent crosses z = 0)

---

## 14. House System Internals

### 14.1 Core Trigonometric Functions

**Location**: swehouse.c

**Asc1(x, f, sine, cose)** — quadrant-aware wrapper that normalizes the input rectascension to [0°, 90°] and delegates to Asc2(), then maps the result back to the correct quadrant.

**Asc2(x, f, sine, cose)** — core oblique spherical trigonometry:

```
cot(λ) = (−tan(f) × sin(ε) + cos(ε) × cos(x)) / sin(x)
λ = arctan(sin(x) / numerator)
```

This solves the spherical triangle formed by the equator, ecliptic, and a great circle with pole height f.

**AscDash(x, f, sine, cose)** — derivative dλ/dt for cusp speed:

```
c = cos(ε)×cos(x) − tan(f)×sin(ε)
d = sin²(x) + c²
dλ/dt = (cos(x)×c + cos(ε)×sin²(x)) / d × 360.985647°/day
```

### 14.2 Ascendant and MC

**MC**: `tan(MC) = tan(ARMC) / cos(ε)` — intersection of meridian with ecliptic.

**Ascendant**: `Asc1(ARMC + 90°, latitude, sin(ε), cos(ε))` — intersection of horizon with ecliptic.

### 14.3 House System Algorithms

**Placidus** (semi-arc division, swehouse.c:1830):
- Computes auxiliary angle: a = arcsin(tan(lat) × tan(ε))
- Divides semi-arc into thirds: fh₁ = arctan(sin(a/3) / tan(ε)), fh₂ = arctan(sin(2a/3) / tan(ε))
- **Iterative Newton-Raphson** (max 100 iterations) to place each cusp on its correct declination circle
- Fails in polar regions (|lat| ≥ 90° − ε)

**Koch** (swehouse.c:1250):
- Simpler single pole height instead of iterative refinement
- Interpolates between meridian (pole = 0) and horizon (pole = latitude)
- `ad₃ = arcsin(sin(c) × sin_a) / 3`, then cusps at ARMC ± offsets

**Regiomontanus** (swehouse.c:1381):
- Fixed pole heights: `fh₁ = arctan(tan(lat) × 0.5)`, `fh₂ = arctan(tan(lat) × cos(30°))`
- No iteration needed

**Campanus** (swehouse.c:1028):
- Pole heights from prime vertical: `fh₁ = arcsin(sin(lat)/2)`, `fh₂ = arcsin(√3/2 × sin(lat))`
- Rectascension offsets from horizon: `xh₁ = arctan(√3 / cos(lat))`

**Topocentric / Polich-Page** (swehouse.c:1432):
- `fh₁ = arctan(tan(lat)/3)`, `fh₂ = arctan(2×tan(lat)/3)`
- Same RA intervals as Regiomontanus but different pole heights

**Alcabitius** (swehouse.c:1581):
- Computes Ascendant's declination, then semi-diurnal/semi-nocturnal arcs
- Divides arcs into thirds, projects onto ecliptic with pole height = 0

**Equal** (swehouse.c:994): `cusp[i] = Asc + (i−1) × 30°`

**Whole Sign** (swehouse.c:1474): `cusp[1] = Asc − (Asc mod 30°)`, then 30° increments

**Porphyry** (swehouse.c:1310): Trisects each quadrant arc (Asc–MC, MC–Dsc, etc.)

**Gauquelin 36 Sectors** (swehouse.c:1623):
- Iterative refinement placing 18 sectors in each semi-arc
- Uses same Placidus-like Newton iteration for each of 36 division points

**Meridian/Axial** (swehouse.c:1485): `λ = arctan(tan(ARMC + i×30°) / cos(ε))` — ecliptic longitude at equatorial intervals

### 14.4 Special Points

| Point | Formula |
|-------|---------|
| Vertex | `Asc1(ARMC − 90°, 90° − lat, sin(ε), cos(ε))` |
| Equatorial Ascendant | `arctan(tan(ARMC + 90°) / cos(ε))` |
| Co-Ascendant (Koch) | `Asc1(ARMC − 90°, lat) + 180°` |
| Co-Ascendant (Munkasey) | `Asc1(ARMC + 90°, 90° − lat)` |
| Polar Ascendant | `Asc1(ARMC − 90°, lat)` |

### 14.5 Speed Computation

For most house systems: analytical derivative via `AscDash()`.

For complex systems (Sunshine, APC, Pullen): finite differences `(cusp(t+dt) − cusp(t−dt)) / (2×dt)`.

---

## 15. Heliacal Visibility Internals

### 15.1 Schaefer's Sky Brightness Model

**Location**: swehel.c

Five light sources contribute to twilight sky brightness (in nanoLamberts):

**1. Daylight (Bday)** — Sun altitude > 4°:

```
F_S = 62000000/RS² + 10^(6.15 − RS/40) + 10^5.36 × (1.06 + cos²(RS))
B_day = F_S × 10^(−0.4×kX_sun) + 440000 × (1 − 10^(−0.4×kX_sun))
B_day × 10^(−0.4×(M_sun − M₀ + 43.27)) × (1 − 10^(−0.4×kX_obj))
```

**2. Twilight (Btwi)** — Sun altitude < −3°:

```
B_twi = 10^(−0.4×(M_sun − M₀ + 32.5 − alt_sun − ZD/(360×k)))
B_twi × (100/RS) × (1 − 10^(−0.4×kX_obj))
```

**3. Moonlight (Bm)** — same scattering formula as daylight with Moon's magnitude.

**4. Zodiacal + starlight (Bn)** — small, modulated by 11.1-year solar cycle.

**5. Light pollution (Bcity)** — configurable.

**Total**: In transition zone (Sun −3° to 4°): `min(B_day, B_twi)`. Moonlight added when B > 200M nL. Zodiacal added when B > 5000 nL.

### 15.2 Atmospheric Extinction

Four exponential atmosphere layers integrated along the line of sight:

| Layer | Scale Height | Coefficient |
|-------|-------------|-------------|
| Rayleigh scattering | 8515 m | k_R |
| Aerosol scattering | 3745 m | k_t (humidity-dependent) |
| Ozone absorption | 20000 m (thin layer) | k_OZ |
| Water vapor | 3000 m | k_W (humidity-dependent) |

Total extinction at zenith distance z: `kX = k_R × X_R(z) + k_t × X_t(z) + k_OZ × X_OZ(z) + k_W × X_W(z)`

where X functions integrate optical depth through each layer's exponential profile.

### 15.3 VisLimMagn() — Limiting Visual Magnitude

**Location**: swehel.c:1382

```
B_sky = total sky brightness (from 5 sources above)
B_sky × CorrFactor₁                            (age/acuity correction)

Vision mode:
    if B_sky < 1645 nL:  scotopic (rod cells)
        C₁ = 10^(−9.8),  C₂ = 10^(−1.9)
    else: photopic (cone cells)
        C₁ = 10^(−8.35), C₂ = 10^(−5.9)

Threshold = C₁ × (1 + √(C₂ × B_sky))² × CorrFactor₂

Limiting magnitude = −16.57 − 2.5 × log₁₀(Threshold)
```

### 15.4 TopoArcVisionis() — Bisection Search

**Location**: swehel.c:1562

Finds the Sun depression angle at which an object of given magnitude becomes visible:

```
bracket: X_L = 45° (deep twilight), X_R = 0° (sunset)
while |X_R − X_L| > ε:
    X_M = (X_R + X_L) / 2
    Y_M = object_mag − VisLimMagn(alt_obj − X_M)
    if sign(Y_M) == sign(Y_L): X_L = X_M else: X_R = X_M
```

### 15.5 Yallop Lunar Crescent Criteria

**Location**: swehel.c:1741

```
q = (ARCV − (11.8371 − 6.3226W' + 0.7319W'² − 0.1018W'³)) / 10
```

Where W' = crescent width in arcminutes, ARCV = geocentric arcus visionis.

| q Range | Grade | Visibility |
|---------|-------|-----------|
| q < −0.14 | A | Impossible |
| −0.14 < q < 0.11 | B | Very difficult (rare) |
| 0.11 < q < 0.34 | C | Possible under good conditions |
| 0.34 < q < 1.15 | D | Easy |
| q > 1.15 | E | Cannot miss |

Crescent width computed from:
```
W = 0.27245 × π × (1 + sin(alt)×sin(π)) × (1 − cos(Δalt)×cos(Δaz))
```

where π = lunar parallax (also approximately the Moon's angular semi-diameter).

---

## 16. Zodiacal Crossings Internals

### 16.1 Newton Iteration (swe_solcross, swe_mooncross)

**Location**: sweph.c

```
Initial estimate: t₀ = tjd + (target_lon − current_lon) / mean_speed

Iterate:
    Compute position and speed at t
    Δlon = swe_degnorm(target − current + 180°) − 180°    (shortest-arc difference)
    dt = Δlon / speed
    t += dt
Until |dt| < 1 milliarcsecond equivalent
```

**Mean speeds used for initial estimate:**
- Sun: 360°/365.24 days
- Moon: 360°/27.32 days

**Edge cases**: Retrograde planets may require the iteration to handle sign changes in speed. The `swe_degnorm(...) − 180°` wrapping ensures convergence from either direction.

### 16.2 Moon Node Crossing (swe_mooncross_node)

Steps forward detecting latitude sign changes, then Newton-refines using latitude and latitude speed:

```
dt = −latitude / latitude_speed
t += dt
```

---

## 17. Orbital Elements Internals

### 17.1 swe_get_orbital_elements()

Derives Keplerian elements from position/velocity state vectors:

```
h⃗ = r⃗ × v⃗                              (angular momentum)
GM = HELGRAVCONST                          (or + planet mass contribution)

a = 1 / (2/r − v²/GM)                     (semi-major axis, vis-viva)
p = |h⃗|² / GM                             (semi-latus rectum)
e = √(1 − p/a)                            (eccentricity)
i = arccos(h_z / |h⃗|)                     (inclination)
Ω = arctan2(h_x, −h_y)                    (longitude of ascending node)
u = arctan2(z/sin(i), x×cos(Ω)+y×sin(Ω)) (argument of latitude)
ν = from vis-viva and radial velocity      (true anomaly)
ω = u − ν                                 (argument of periapsis)
E = 2×arctan(√((1−e)/(1+e)) × tan(ν/2))  (eccentric anomaly)
M = E − e×sin(E)                          (mean anomaly, Kepler's equation)
```

### 17.2 swe_orbit_max_min_true_distance()

Samples both the planet's orbit and Earth's orbit in 2° steps, computes geocentric distance at each point, then iteratively refines the maximum and minimum using quadratic interpolation.

---

## 18. Planetocentric Coordinates

### 18.1 swe_calc_pctr()

**Location**: sweph.c:8042

Computes position of body `ipl` as seen from center body `iplctr`:

1. Compute both bodies in barycentric J2000
2. `pos_relative = pos_planet − pos_center`
3. Light-time iteration using center body as observer (not Earth)
4. Apply same aberration/deflection/precession/nutation pipeline
5. Returns position in center body's coordinate frame

---

## 19. Gauquelin Sector Internals

### 19.1 Geometric Method (imeth 0–1)

```
armc = swe_sidtime0(...) × 15 + geolon
sector = swe_house_pos(armc, geolat, ε, 'G', planet_xyz)
```

Uses the Gauquelin house system ('G') which divides the diurnal arc into 18 sectors.

### 19.2 Rise/Set Method (imeth 2–3)

1. Find rise and set times bracketing the current moment
2. Determine if object is above or below horizon
3. Linear interpolation:
   - Above horizon: `sector = (t − t_rise) / (t_set − t_rise) × 18 + 1`
   - Below horizon: `sector = (t − t_set) / (t_rise − t_set) × 18 + 19`

imeth 2 excludes refraction; imeth 3 includes it.

---

## 20. Fixed Star Computation Internals

### 20.1 Catalog Loading

**`swe_fixstar2()`** loads the entire star catalog into memory on first call via `load_all_fixed_stars()`. Each star record (`struct fixed_star`) stores: search key, traditional name, Bayer designation, catalog number, epoch, RA, Dec, proper motion, radial velocity, parallax, magnitude.

### 20.2 Position Computation

For each star at observation epoch:

1. **Epoch transformation** (if catalog epoch ≠ J2000):
   ```
   RA_J2000 = RA_catalog + RA_proper_motion × (J2000 − epoch)
   Dec_J2000 = Dec_catalog + Dec_proper_motion × (J2000 − epoch)
   ```
2. **Precession**: J2000 → observation date via `swi_precess()`
3. **Proper motion velocity vector** also precessed
4. **Annual aberration** and **gravitational deflection**: same pipeline as planets
5. **Sidereal transformation**: if `SEFLG_SIDEREAL`, subtract ayanamsa

### 20.3 Built-in Stars

A small array of major stars is compiled into the library for use when the catalog file is unavailable. These are searched before the file-based catalog.

---

## 21. Cross-Domain Dependencies

```
swe_sol_eclipse_when_loc()
  └─ swe_rise_trans()             (eclipse visibility from location needs rise/set)
       └─ swe_calc_ut()           (planet positions for altitude)
            └─ swe_deltat_ex()    (UT→TT conversion)
                 └─ tidal acceleration model depends on ephemeris backend

swe_heliacal_ut()
  ├─ swe_rise_trans()             (find sunrise/sunset)
  ├─ swe_calc_ut()                (planet positions)
  ├─ swe_azalt()                  (horizontal coordinates)
  │    └─ swe_sidtime()           (sidereal time)
  │         └─ swi_nutation()     (equation of equinoxes)
  ├─ swe_pheno_ut()               (planet magnitude)
  └─ swe_refrac_extended()        (atmospheric refraction)

swe_houses_ex2()
  ├─ swe_sidtime()                (ARMC from sidereal time)
  ├─ swi_nutation()               (for NONUT suppression)
  └─ swe_calc_ut(SE_SUN)          (Sun declination for Sunshine houses)

swe_gauquelin_sector()
  ├─ swe_house_pos()              (geometric method)
  └─ swe_rise_trans()             (rise/set method)

swe_lun_eclipse_when_loc()
  └─ swe_rise_trans(SE_MOON)      (Moon visibility from location)
```
