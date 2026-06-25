# Codebase Map

Quick reference for implementation agents. Prevents 60k-token exploration sweeps.

## Module Layout

```
src/
├── lib.rs              — pub mod declarations (lines 1–21), re-exports (lines 23–30)
├── types.rs            — all domain types, enums, newtypes (778 lines)
├── constants.rs        — physical constants, epochs, unit conversions (198 lines)
├── flags.rs            — bitflags! structs: CalcFlags, SiderealBits, etc. (146 lines)
├── error.rs            — Error enum
├── context.rs          — Ephemeris, EphemerisConfig, CalcResult
├── math.rs             — pure math functions: normalize, chebyshev, cartpol, cotrans, poly_eval
├── date.rs             — Julian Day ↔ calendar conversion, delta-T, UTC
├── obliquity.rs        — swi_epsiln port: all 11 obliquity models
├── bias.rs             — swi_bias port: GCRS↔J2000 frame rotation
├── precession.rs       — swi_precess port: 3 algorithm families, 11 models, JPLHOR paths
├── deltat/
│   ├── mod.rs          — calc_deltat dispatcher, 5 historical models, Bessel interpolation, future extrapolation, tidal correction
│   └── data.rs         — static tables: DT (409 entries 1620–2028), DT97 (43), DT2 (27), DTCF16 (54×6 spline)
├── nutation/
│   ├── mod.rs          — router + 5 algorithms: IAU 1980, Herring 1987, IAU 2000A/B, Woolard
│   └── data.rs         — generated nutation term tables (IAU 2000A, 2000B, 1980)
├── calc.rs             — EMPTY stub
├── moshier/mod.rs      — EMPTY stub
├── jpl.rs              — EMPTY stub
├── sweph_file.rs       — EMPTY stub
├── houses.rs           — EMPTY stub
├── eclipse.rs          — EMPTY stub
├── ayanamsa.rs         — EMPTY stub
├── heliacal.rs         — EMPTY stub
├── phenomena.rs        — EMPTY stub
└── stars.rs            — EMPTY stub

tests/
├── golden/
│   ├── main.rs         — test harness: golden_data_path(), assert_f64_exact(), assert_f64_eps()
│   ├── math.rs         — golden tests for math module
│   ├── date.rs         — golden tests for date module
│   ├── obliquity_bias.rs — golden tests for obliquity + bias
│   ├── precession.rs  — golden tests for precession (374 cases)
│   ├── nutation.rs    — golden tests for nutation (80 cases + router tests)
│   └── deltat.rs      — golden tests for delta-T (217 cases: 5 models × 43 epochs)
├── golden-data/
│   ├── math.json       — C-generated reference data for math
│   ├── date.json       — C-generated reference data for date
│   ├── obliquity_bias.json — C-generated reference data for obliquity/bias
│   ├── precession.json — C-generated reference data for precession
│   ├── nutation.json   — C-generated reference data for nutation
│   └── deltat.json     — C-generated reference data for delta-T
└── c-gen/
    ├── gen_obliquity_bias.c — C harness to regenerate obliquity_bias.json
    ├── gen_precession.c — C harness to regenerate precession.json
    ├── gen_nutation.c  — C harness to regenerate nutation.json
    └── gen_deltat.c    — C harness to regenerate deltat.json
```

## Key Types in types.rs

### Astronomical model enums (lines 517–597)

| Type | Lines | Variants | repr |
|---|---|---|---|
| `PrecessionModel` | 522–534 | IAU1976=1..Newcomb=11 (11 variants) | i32 |
| `NutationModel` | 538–544 | IAU1980=1..Woolard=5 | i32 |
| `DeltaTModel` | 548–554 | 5 variants | i32 |
| `SiderealTimeModel` | 558–563 | IAU1976=1..Longterm=4 | i32 |
| `BiasModel` | 567–571 | None=1, IAU2000=2, IAU2006=3 | i32 |
| `JplHorMode` | 575–577 | LongAgreement=1 | i32 |
| `JplHoraMode` | 581–585 | V1=1, V2=2, V3=3 | i32 |
| `AstroModels` | 588–597 | 8 fields: delta_t, prec_longterm, prec_shortterm, nutation, bias, jplhor_mode, jplhora_mode, sidereal_time |
| `Default` impl | 680–693 | longterm=shortterm=Vondrak2011, nutation=IAU2000B, bias=IAU2006, jplhora=V3, sidereal=Longterm |

### Frame/obliquity types (after AstroModels, before JdTt)

| Type | Lines | Notes |
|---|---|---|
| `FrameTransform` | ~599 | J2000ToGcrs, GcrsToJ2000 |
| `PrecessionDirection` | ~603 | J2000ToDate, DateToJ2000 |
| `Epsilon` | ~607 | eps, sin_eps, cos_eps + `Epsilon::new(eps_rad)` |
| `Nutation` | ~628 | dpsi, deps (radians) |

### Julian Day newtypes

| Type | Line | Notes |
|---|---|---|
| `JdTt` | ~624 | newtype `(pub f64)` |
| `JdUt1` | ~627 | newtype `(pub f64)` |

### Other key types

| Type | Lines | Notes |
|---|---|---|
| `Body` | 82–111 | enum: Sun..Vesta + Fictitious/Asteroid/PlanetMoon/Comet |
| `HouseSystem` | 272–298 | 22 house system variants |
| `CalendarType` | 373–376 | Julian, Gregorian |
| `SiderealMode` | 396–445 | 42 sidereal mode variants |
| `EphemerisSource` | 510–514 | Jpl, Swisseph, Moshier |
| `UtcComponents` | 640–647 | year, month, day, hour, min, sec |
| `DeltaT` | 659–661 | trait |
| `DegreeParts` | 668–674 | degrees, minutes, seconds, second_fraction, sign |

## Flags (src/flags.rs)

Six `bitflags!` structs. Most relevant:

**CalcFlags (u32)** — lines 3–31:
- JPLEPH=1, SWIEPH=2, MOSEPH=4, HELCTR=8, TRUEPOS=16, J2000=32
- NONUT=64, SPEED3=128, SPEED=256, NOGDEFL=512, NOABERR=1024
- EQUATORIAL=2048, XYZ=4096, RADIANS=8192, BARYCTR=16384
- TOPOCTR=32768, SIDEREAL=65536, ICRS=131072
- **DPSIDEPS_1980=262144** (C: SEFLG_JPLHOR)
- **JPLHOR_APPROX=524288**
- CENTER_BODY=1048576

## Constants (src/constants.rs)

Key constants for quick reference:

| Name | Value | Line |
|---|---|---|
| DEGTORAD | π/180 | 29 |
| STR | 4.8481368e-6 (arcsec→rad) | 38 |
| TWOPI | 2π | 34 |
| J2000 | 2451545.0 | 59 |
| B1950 | 2433282.42345905 | 60 |
| J1900 | 2415020.0 | 61 |
| B1850 | 2396758.2035810 | 62 |
| PREC_IAU_1976_CTIES | 2.0 | ~64 |
| PREC_IAU_2000_CTIES | 2.0 | ~65 |
| PREC_IAU_2006_CTIES | 75.0 | ~66 |
| DPSI_DEPS_IAU1980_TJD0_HORIZONS | 2437684.5 | ~70 |

## Math Functions (src/math.rs)

All `pub fn`. Key functions and their line ranges:

| Function | Lines | Signature |
|---|---|---|
| normalize_degrees | 11–20 | (f64) → f64 |
| normalize_radians | 22–31 | (f64) → f64 |
| mod_2pi | 33–39 | (f64) → f64 |
| diff_degrees_norm | 45–47 | (f64, f64) → f64 |
| diff_degrees | 49–52 | (f64, f64) → f64 |
| diff_radians | 54–57 | (f64, f64) → f64 |
| midpoint_degrees | 63–66 | (f64, f64) → f64 |
| midpoint_radians | 68–70 | (f64, f64) → f64 |
| csnorm | 76–93 | (i32) → i32 |
| difcsn | 95–97 | (i32, i32) → i32 |
| difcs2n | 99–102 | (i32, i32) → i32 |
| d2l | 108–114 | (f64) → i32 |
| chebyshev_eval | 120–131 | (f64, &[f64]) → f64 |
| chebyshev_deriv | 133–156 | (f64, &[f64]) → f64 |
| rotate_x | 162–166 | ([f64;3], f64) → [f64;3] |
| rotate_x_sincos | 168–174 | ([f64;3], f64, f64) → [f64;3] |
| cartesian_to_polar | 176–193 | ([f64;3]) → [f64;3] |
| polar_to_cartesian | 195–202 | ([f64;3]) → [f64;3] |
| cartesian_to_polar_with_speed | 208–237 | ([f64;6]) → [f64;6] |
| polar_to_cartesian_with_speed | 239–261 | ([f64;6]) → [f64;6] |
| cotrans | 267–274 | ([f64;3], f64) → [f64;3] |
| cotrans_with_speed | 276–290 | ([f64;6], f64) → [f64;6] |
| split_degrees | 304–364 | (f64, SplitDegFlags) → DegreeParts |
| poly_eval | 366–368 | (&[f64], f64) → f64 — Horner's method |
| OWEN_T0S | 374 | [f64; 5] — Owen interval boundaries |
| owen_t0_icof | 376 | (f64) → (f64, usize) — Owen interval + index |
| owen_chebyshev_basis | 390 | (f64) → (usize, [f64; 10]) — shared by obliquity + precession |
| **unit tests** | 410+ | |

## Golden Test Pattern

1. JSON data in `tests/golden-data/<name>.json` — top-level object, keys are test groups, values are arrays of case structs
2. Rust test file in `tests/golden/<name>.rs`:
   - `#[derive(Deserialize)]` case structs matching JSON fields
   - Top-level struct aggregating all case vectors
   - `fn load() -> TopStruct` using `super::golden_data_path()`
   - `#[test]` per group calling `super::assert_f64_exact()` or `super::assert_f64_eps()`
3. Add `mod <name>;` to `tests/golden/main.rs`
4. C generator in `tests/c-gen/` compiled against `../swisseph/libswe.a`

## C Source Reference

- C repo: `../swisseph/`
- C catalogues: `../swisseph/claude/catalogue-{public,internal}.md`
- Ephemeris data: `../swisseph/ephe/`
- Pre-built library: `../swisseph/libswe.a`
- Key headers: `swephexp.h` (public API, model enums), `swephlib.h` (internal functions, PREC constants), `sweph.h` (internal constants)
- `swe_set_astro_models(char *samod, ...)` — comma-separated string: "delta_t,prec_long,prec_short,nutation,bias,jplhor_mode,jplhora_mode,sidt"

## Floating-Point Fidelity Notes

- C polynomial models compute `result * DEGTORAD / 3600` (two runtime ops). Rust must match: `poly_eval(...) * DEGTORAD / 3600.0` — NOT `* STR`.
- C `eps *= DEGTORAD/3600.0` (Laskar, Vondrák) folds to single multiply. Rust: `* (DEGTORAD / 3600.0)` with parens.
- Owen 1990 returns degrees (not arcsec): multiply by `DEGTORAD`, no `/3600`.
- Vondrák `swi_ldp_peps` returns radians directly via `* AS2R` = `* (DEGTORAD / 3600.0)`.

## Insertion Points for New Modules

| What | Where | After |
|---|---|---|
| New model enums | src/types.rs | ~585 (after JplHoraMode) |
| New shared types | src/types.rs | after AstroModels block |
| New AstroModels field | src/types.rs | line 596 + update Default at 681 |
| New constants | src/constants.rs | after existing epoch block |
| New `pub mod` | src/lib.rs | alphabetical in lines 7–21 |
| New re-export | src/lib.rs | inside `pub use types::{...}` block |
| New golden test mod | tests/golden/main.rs | after existing mod declarations |
