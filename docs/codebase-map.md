# Codebase Map

Quick reference for implementation agents. Prevents 60k-token exploration sweeps.

## Module Layout

```
src/
├── lib.rs              — pub mod declarations (lines 1–21), re-exports (lines 23–30)
├── types.rs            — all domain types, enums, newtypes (778 lines)
├── constants.rs        — physical constants, epochs, unit conversions (199 lines)
├── corrections.rs      — relativistic corrections: meff (lookup), aberr_light (Lorentz), deflect_light (GR bending)
├── flags.rs            — bitflags! structs: CalcFlags, SiderealBits, etc. (146 lines)
├── error.rs            — Error enum
├── context.rs          — Ephemeris (calc, calc_ut, calc_inner, calc_speed3, extract_for_body, fixstar2, fixstar2_ut, fixstar2_mag, calc_fixstar, houses/houses_ex/houses_ex2 — UT-based house wrappers: ARMC+eps+nutation setup, Sun-declination resolution for Sunshine, traditional-sidereal dispatch, RADIANS conversion, gauquelin_sector_geometric — swe_gauquelin_sector imeth 0/1 port: own ARMC/eps/nutation setup using the caller's flags directly (NOT houses_ex2's forced TIDAL_DEFAULT — swe_gauquelin_sector's C deltaT call genuinely carries the caller's ephemeris-source iflag), calc() the body, house_pos with hsys=Gauquelin; imeth 2-5 (rise/set) return Err, not yet ported), EphemerisConfig, CalcResult; stars: StarCatalog field on Ephemeris
├── math.rs             — pure math functions: normalize, chebyshev, cartpol, cotrans, poly_eval
├── date.rs             — Julian Day ↔ calendar conversion, delta-T, UTC
├── obliquity.rs        — swi_epsiln port: all 11 obliquity models
├── bias.rs             — swi_bias port: GCRS↔J2000 frame rotation; icrs2fk5 (RB matrix), fk4_fk5 (B1950 RA correction)
├── precession.rs       — swi_precess port: 3 algorithm families, 11 models, JPLHOR paths
├── deltat/
│   ├── mod.rs          — calc_deltat dispatcher, 5 historical models, Bessel interpolation, future extrapolation, tidal correction
│   └── data.rs         — static tables: DT (409 entries 1620–2028), DT97 (43), DT2 (27), DTCF16 (54×6 spline)
├── nutation/
│   ├── mod.rs          — router + 5 algorithms: IAU 1980, Herring 1987, IAU 2000A/B, Woolard
│   └── data.rs         — generated nutation term tables (IAU 2000A, 2000B, 1980)
├── sidereal_time.rs    — swe_sidtime0/swe_sidtime port: 4 GMST models, 33-term EoE, long-term model
│                          (sidtime_long_term, pre-1850/post-2050 dates) — its 2 internal deltaT
│                          calls force TIDAL_DEFAULT via a local deltat_config override, matching
│                          C's swe_deltat_ex(tjd, -1, NULL) sentinel (swephlib.c:3291,3301), which
│                          always resolves SE_TIDAL_DEFAULT regardless of the actually-configured
│                          ephemeris backend — the same deltaT/tid_acc-inconsistency pattern as
│                          houses_ex2's; discovered via a gauquelin_sector golden-test mismatch at
│                          1600/1800 AD epochs (swisseph-rs/66)
├── calc.rs             — calc_planet, calc_sun, calc_moon, calc_mean_node, calc_mean_apogee, calc_ecl_nut, extract_output, extract_ecl_nut, plaus_iflag, speed3_interval, denormalize_positions, calc_speed_3point: light-time, retarded velocity, aberration, deflection pipeline + mean element pipeline + SPEED3 helpers
├── moshier/
│   ├── mod.rs          — PlantTbl struct, PLANETS array re-export, element-count tests
│   ├── backend.rs      — compute() public API, compute_pipeline() for calc.rs, embofs_mosh, planet/earth velocity helpers, Body dispatch
│   ├── moon.rs         — moshmoon2() lunar series evaluator: MeanElements, mean_elements(), chewm(), moon1–4, mean_node(), mean_apogee(), correction interpolation
│   ├── moon_tables.rs  — generated const arrays: LR/MB/LRT/BT/LRT2/BT2 + z[25] + MEAN_NODE_CORR[304] + MEAN_APSIS_CORR[304]
│   ├── planets.rs      — moshplan2() series evaluator, sscc() harmonic recurrence, fundamental argument constants
│   └── tables.rs       — generated const arrays: 9 planet tables (do not hand-edit, see scripts/gen_moshier_tables.py)
├── jpl/
│   ├── mod.rs          — JplFile (mmap + JplHeader), JplFile::open, byte_order/header/bytes accessors. Re-exports ByteOrder, JplHeader. J_* body index constants. pub fn jpl_pleph (body assembly entry point).
│   ├── header.rs       — ByteOrder enum + Reader cursor, detect_byte_order (plausibility of ss[2]), parse_header (record-0 offsets), compute_ksize (ipt[] algorithm), validate_file_length, JplHeader struct
│   └── interp.rs       — read_record (mmap→Vec<f64>), interp (JPL forward-recurrence Chebyshev eval + sub-interval selection), state (record selection + body interpolation loop)
├── sweph_file/
│   ├── mod.rs          — SwissEphFile (mmap-based .se1 reader), body_file_id(Body → ipl value), evaluate_body re-export
│   ├── types.rs        — FileHeader, PlanetFileData, FileType, ByteOrder, SEI_*/SE_* body constants, SEI_FLG_* flags
│   ├── parse.rs        — binary format parser: Reader cursor, detect_byte_order, parse_file (header + per-planet metadata)
│   ├── segment.rs      — Chebyshev coefficient unpacking from mmap'd bytes: 6 packing modes (4/3/2/1-byte, nibble, quarter-byte)
│   └── evaluate.rs     — rot_back (orbital-plane→ecliptic/equatorial transform), evaluate_body (public API: file + body_id + jd → [x,y,z,vx,vy,vz])
├── houses.rs           — AscMc, HouseResult (public types); houses_armc driver (swe_houses_armc_ex2 port);
│                          calc_h (CalcH core; systems implemented: A/D/N/V/W equal-family,
│                          O/S/X/M/F quadrant-arithmetic, R/C/T/H/J great-circle/pole-height,
│                          P/K/G Newton-iteration (Placidus/Koch/Gauquelin-36),
│                          U/Y/L/Q closed-form misc (Krusinski cotrans chain, APC apc_sector,
│                          Pullen SD/SR), I/i Sunshine (Treindl + Makransky, sunshine_init,
│                          finite-diff cusp speeds, Porphyry fallback on Makransky ERR) —
│                          B still stubbed Err);
│                          Asc1/Asc2/AscDash core trig, fix_asc_polar, mc_like (shared MC/equasc
│                          projection), polar_shift_subset (shared C/H/J/R polar-circle 180° flip),
│                          NewtonCusp/placidus_newton_cusp (shared P/G Newton-iteration skeleton),
│                          apc_sector (radians-domain helper for Y);
│                          sidereal_houses_trad (traditional sidereal: houses_armc at tropical
│                          armc/eps then subtract ayanamsa from cusps/ascmc except armc; W routed
│                          through Equal + re-fixed to 30° multiples; N re-fixed to (i-1)*30
│                          unconditionally) — ayanamsa is passed in, kept pure;
│                          sidereal_houses_ecl_t0/sidereal_houses_ssypl (geometric sidereal
│                          projections onto the ecliptic-of-t0/solar-system-plane: build a vernal-
│                          point-like moving vector, precess/nutate it to tjde's true equator,
│                          derive an auxiliary obliquity+vernal-point from its orbital-plane normal
│                          via cross_prod/dot_prod_unit, compute houses at the resulting
│                          armcx/epsx, then subtract dvpxe+ayan_t0(+x00 for ssypl) from every
│                          cusp/ascmc except armc) — t0/ayan_t0 passed in raw (unlike
│                          ayanamsa::resolve_t0's callers, C's own sidereal_houses_ecl_t0/ssypl
│                          never apply the t0_is_UT deltaT adjustment); shared helpers
│                          rotate_to_true_equator, sidereal_houses_geom_core, apply_sidereal_shift;
│                          house_pos (swe_house_pos port: planet ecl.[lon,lat] → continuous house
│                          position 1.0..13.0, all 25 HouseSystem variants incl. Alcabitius which
│                          calc_h doesn't implement yet — house_pos's Alcabitius branch is
│                          self-contained; Koch returns Err on circumpolar dfac-out-of-range;
│                          Sunshine requires Some(sundec) in range or Err) — helpers:
│                          armc_to_mc_house_pos (degnorm-before-and-after-+180 variant, distinct
│                          from crate::math::armc_to_mc and mc_like — c-ref-houses.md §12.1),
│                          mc_transform_raw (Morinus's un-normalized tand/cose transform),
│                          bracket_interpolate_12 (shared J/default fallback), koch_house_pos,
│                          topocentric_house_pos, sunshine_apc_house_pos (shared I/i/Y formula)
├── eclipse.rs          — EMPTY stub
├── ayanamsa.rs         — EMPTY stub
├── heliacal.rs         — EMPTY stub
├── phenomena.rs        — EMPTY stub
└── stars.rs            — StarCatalog, Star, load_catalog, builtin_star (8 ayanamsa ref stars), search, parse

tests/
├── golden/
│   ├── main.rs         — test harness: golden_data_path(), assert_f64_exact(), assert_f64_eps()
│   ├── calc.rs        — golden tests for calc pipeline (1176 cases: 14 bodies × 7 epochs × 12 flag combos incl. SPEED3, no_speed)
│   ├── corrections.rs — golden tests for corrections (30 meff + 40 aberr + 15 pipeline)
│   ├── math.rs         — golden tests for math module
│   ├── date.rs         — golden tests for date module
│   ├── obliquity_bias.rs — golden tests for obliquity + bias
│   ├── precession.rs  — golden tests for precession (374 cases)
│   ├── nutation.rs    — golden tests for nutation (80 cases + router tests)
│   ├── deltat.rs      — golden tests for delta-T (217 cases: 5 models × 43 epochs)
│   ├── sidereal_time.rs — golden tests for sidereal time (128 cases: 4 models × 32 epochs)
│   ├── mean_elements.rs — golden tests for mean node, mean apogee, ECL_NUT (165 cases: 3 bodies × 11 epochs × 5 flag combos)
│   ├── moshier_backend.rs — golden tests for backend::compute (110 cases: 10 bodies × 11 epochs + Earth zero-check)
│   ├── moshier_moon.rs — golden tests for moshmoon2 (11 cases: Moon at 11 epochs)
│   ├── moshier_planet.rs — golden tests for moshplan2 (81 cases: 9 planets × 9 epochs)
│   ├── se1_header.rs  — golden tests for SE1 file parsing (11 planet metadata fields, byte-order detection on 84 files)
│   ├── sweph_eval.rs  — golden tests for evaluate_body (80 cases: 10 bodies × 8 epochs, bitwise-exact positions + velocities)
│   ├── jpl_pleph.rs   — golden tests for jpl_pleph (84 cases: 11 bodies × 7 epochs barycentric + 7 geocentric Moon, 1e-9 eps)
│   └── houses.rs      — golden tests for houses_armc (angles_special: 30 cases system-independent special points;
│                         equal_family: 150 cases, 5 systems A/D/N/V/W × 30 battery cases, bitwise-exact;
│                         quad_arith: 150 cases, 5 systems O/S/X/M/F × 30 battery cases, eps 1e-9/1e-7;
│                         great_circle: 150 cases, 5 systems R/C/T/H/J × 30 battery cases, eps 1e-9;
│                         iterative: 84 cases, 2 systems P/K × 6 armc × 7 geolat (incl. ±78 polar) ×
│                         1 eps, eps 1e-9 cusps/1e-7 speeds; gauquelin36: 42 cases, G × 6 armc × 7 geolat
│                         × 1 eps, cusps[1..36], same eps; closed_form_misc: 120 cases, 4 systems
│                         U/Y/L/Q × 30 battery cases, eps 1e-9 cusps/1e-7 speeds (U speeds eps 0 —
│                         stale pre-switch values, asserted exactly per c-ref-houses.md §4.2e);
│                         sunshine: 76 cases — 60: 2 systems I/i × 6 armc × 5 geolat (1 sundec per
│                         case, rotated through {-23,-10,0,10,23}), eps 1e-9 cusps (1e-8 for
│                         Makransky 'i')/1e-7 speeds (driver-level finite-diff); + 16: 2 systems ×
│                         2 armc × geolat {70,-70} × sundec {23,-23} (all four combinations satisfy
│                         |tand(geolat)·tand(sundec)|≥1, triggering Makransky's circumpolar ERR →
│                         Porphyry fallback; Treindl never short-circuits on it, included at the
│                         same combos for contrast), eps 1e-9 cusps/speeds for the 'i' fallback
│                         subset (fill_porphyry is closed-form, bitwise-exact elsewhere);
│                         sunshine_requires_sundec: negative test, Sunshine + sundec=None
│                         returns Err; ut_wrapper: 42 cases — 36: 6 (tjd_ut,geolat,geolon) triples
│                         × 6 systems P/K/C/R/W/I + 6: 1 triple × 6 systems with SEFLG_NONUT,
│                         eps 1e-7 cusps/ascmc, 1e-6 speeds via Ephemeris::houses_ex2 (compounds
│                         deltaT/obliquity/nutation/sidtime — looser than the pure-armc tests);
│                         sidereal_trad: 9 cases, systems P/W/E × 3 triples, SEFLG_SIDEREAL +
│                         Lahiri, same tolerances; sidereal_geom: 18 cases, systems P/C/W x 3
│                         triples x 2 sid_modes (Lahiri|SE_SIDBIT_ECL_T0, Lahiri|SE_SIDBIT_SSY_PLANE),
│                         same tolerances, exercises sidereal_houses_ecl_t0/ssypl; house_pos: 150
│                         cases, all 25 house-system chars
│                         × 2 armc/geolat/eps triples (1 temperate, 1 polar — armc=105/geolat=67
│                         chosen so 2/3 xpin succeed and 1/3 hits Koch's genuine hpos==0 circumpolar
│                         failure) × 3 xpin, eps 1e-7 hpos, "err" cases (1: Koch) assert Err instead
│                         (err is driven by hpos==0.0, not serr non-empty — several systems set an
│                         informational serr on a valid nonzero hpos); gauquelin_sector: 36 cases,
│                         6 (tjd_ut,ipl,imeth) combos × 3 bodies (Sun/Moon/Mars) × 2 imeth, eps 1e-6
│                         dgsect via Ephemeris::gauquelin_sector_geometric — this test caught a
│                         pre-existing sidtime_long_term deltaT/tid_acc bug at 1600/1800 AD epochs,
│                         fixed in sidereal_time.rs (see that file's codebase-map entry))
├── golden-data/
│   ├── calc.json       — C-generated reference data for calc pipeline (swe_calc full pipeline)
│   ├── corrections.json — C-generated reference data for corrections (meff, aberr_light, pipeline)
│   ├── math.json       — C-generated reference data for math
│   ├── date.json       — C-generated reference data for date
│   ├── obliquity_bias.json — C-generated reference data for obliquity/bias
│   ├── precession.json — C-generated reference data for precession
│   ├── nutation.json   — C-generated reference data for nutation
│   ├── deltat.json     — C-generated reference data for delta-T
│   ├── sidereal_time.json — C-generated reference data for sidereal time
│   ├── mean_elements.json — C-generated reference data for mean node, mean apogee, ECL_NUT
│   ├── moshier_backend.json — C-generated reference data for backend::compute (swe_calc with ICRS)
│   ├── moshier_moon.json — C-generated reference data for moshmoon2
│   ├── moshier_planet.json — C-generated reference data for moshplan2
│   ├── se1_header.json — C-generated reference data for SE1 file headers (sepl_18, semo_18)
│   ├── sweph_eval.json — C-generated reference data for evaluate_body (raw Chebyshev eval + rot_back + ecl→equ rotation)
│   ├── jpl_pleph.json  — C-generated reference data for jpl_pleph (84 cases via swi_pleph against de441.eph)
│   ├── fixstar.json    — C-generated reference data for swe_fixstar2 (196 position cases + 4 mag cases, 7 stars × 4 epochs × 7 flags)
│   └── houses.json     — C-generated reference data for swe_houses_armc_ex2 (battery: 6 armc × 5 geolat × 1 eps, reused across all houses sub-tasks; iterative/gauquelin36 keys add a 7th/8th polar geolat (±78) to exercise the Placidus/Koch/Gauquelin Porphyry fallback; closed_form_misc key reuses the standard 5-geolat battery for U/Y/L/Q; sunshine key reuses the standard 6 armc × 5 geolat battery for I/i, crossed with a rotated (not full cross-product) Sun-declination set {-23,-10,0,10,23}, plus a dedicated circumpolar-Sun sub-battery (geolat {70,-70} × sundec {23,-23}) to exercise Makransky's ERR→Porphyry fallback; ut_wrapper key: swe_houses_ex2 (UT-based) over 6 (tjd_ut,geolat,geolon) triples × 6 systems, + a SEFLG_NONUT variant at 1 triple; sidereal_trad key: swe_houses_ex2 with SEFLG_SIDEREAL + swe_set_sid_mode(SE_SIDM_LAHIRI) over 3 triples × 3 systems P/W/E; house_pos key: swe_house_pos over all 25 house-system chars × 2 (armc,geolat,eps) triples × 3 xpin, "err" field is hpos==0.0 (Koch's real failure sentinel), NOT serr-non-empty (P/G/J/L/Q/default set an informational serr on valid results) — the static sundec cache 'I'/'i' need is primed via a preceding swe_houses_armc_ex2(ascmc[9]=sundec) call; gauquelin_sector key: swe_gauquelin_sector imeth∈{0,1} over 6 ut_triples × 3 bodies (Sun/Moon/Mars))
└── c-gen/
    ├── gen_calc.c      — C harness to regenerate calc.json (full swe_calc pipeline, 14 bodies × 7 epochs × 12 flags, ECL_NUT cleanup)
    ├── gen_mean_elements.c — C harness to regenerate mean_elements.json (mean node, mean apogee, ECL_NUT)
    ├── gen_corrections.c — C harness to regenerate corrections.json (meff copied from sweph.c, swi_aberr_light direct, pipeline via swe_calc)
    ├── gen_obliquity_bias.c — C harness to regenerate obliquity_bias.json
    ├── gen_precession.c — C harness to regenerate precession.json
    ├── gen_nutation.c  — C harness to regenerate nutation.json
    ├── gen_deltat.c    — C harness to regenerate deltat.json
    ├── gen_sidereal_time.c — C harness to regenerate sidereal_time.json
    ├── gen_moshier_backend.c — C harness to regenerate moshier_backend.json (swe_calc with all corrections disabled + ICRS)
    ├── gen_moshier_moon.c — C harness to regenerate moshier_moon.json
    ├── gen_moshier_planet.c — C harness to regenerate moshier_planet.json
    ├── gen_sweph_eval.c — C harness to regenerate sweph_eval.json (raw SE1 Chebyshev eval via swed.pldat internals)
    ├── gen_se1_header.c — standalone binary parser, dumps header + planet metadata as JSON
    ├── gen_jpl_pleph.c  — C harness to regenerate jpl_pleph.json (swi_pleph direct calls against de441.eph)
    ├── gen_fixstar.c    — C harness to regenerate fixstar.json (swe_fixstar2: 7 stars × 4 epochs × 7 flags + 4 mag cases)
    └── gen_houses.c     — C harness to regenerate houses.json (swe_houses_armc_ex2: angles_special,
                            equal_family, quad_arith, great_circle, iterative, gauquelin36,
                            closed_form_misc, sunshine — sunshine key sets ascmc[9]=sundec before
                            calling, per c-ref-houses.md §11; ut_wrapper/sidereal_trad keys use
                            swe_houses_ex2 (UT-based) instead, over a 6-triple
                            (tjd_ut,geolat,geolon) battery; house_pos key: swe_house_pos over all
                            25 house-system chars, primes the 'I'/'i' static sundec cache via a
                            swe_houses_armc_ex2(ascmc[9]=sundec) call immediately before each
                            swe_house_pos call; gauquelin_sector key: swe_gauquelin_sector
                            imeth∈{0,1} reusing the ut_wrapper triples × 3 bodies)
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
| mods3600 | 41–43 | (f64) → f64 — arcsec modulo 1296000 |
| diff_degrees_norm | 49–51 | (f64, f64) → f64 |
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
| cross_prod | 167–173 | ([f64;3], [f64;3]) → [f64;3] — swi_cross_prod |
| dot_prod_unit | 177–185 | ([f64;3], [f64;3]) → f64 — swi_dot_prod_unit, clamped to [-1,1] |
| rotate_x | 190–194 | ([f64;3], f64) → [f64;3] |
| rotate_x_sincos | 196–202 | ([f64;3], f64, f64) → [f64;3] |
| cartesian_to_polar | 204–221 | ([f64;3]) → [f64;3] |
| polar_to_cartesian | 223–230 | ([f64;3]) → [f64;3] |
| cartesian_to_polar_with_speed | 236–265 | ([f64;6]) → [f64;6] |
| polar_to_cartesian_with_speed | 267–289 | ([f64;6]) → [f64;6] |
| cotrans | 295–302 | ([f64;3], f64) → [f64;3] |
| cotrans_with_speed | 304–318 | ([f64;6], f64) → [f64;6] |
| split_degrees | 332–392 | (f64, SplitDegFlags) → DegreeParts |
| poly_eval | 394–396 | (&[f64], f64) → f64 — Horner's method |
| OWEN_T0S | 402 | [f64; 5] — Owen interval boundaries |
| owen_t0_icof | 404 | (f64) → (f64, usize) — Owen interval + index |
| owen_chebyshev_basis | 418 | (f64) → (usize, [f64; 10]) — shared by obliquity + precession |
| **unit tests** | 438+ | |

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
- **`+=` vs `= x +`**: C's `L = L + a + b + c` accumulates left-to-right with L in each step. Rust's `l += a + b + c` evaluates `a + b + c` first, then adds to l. When L is large (~481k) and corrections are small (~6), the different accumulation order produces ULP-level rounding differences that propagate through backward-difference velocity (÷1e-4) and deflection speed (÷5e-7). Always use `l = l + ...` to match C's evaluation order.
- **Multiplication order matters**: `2.0 * x * DEGTORAD` ≠ `2.0 * DEGTORAD * x` due to FP non-associativity. Match C's grouping exactly.

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
