# Golden Data Testing Harness

Differential tests that assert our Rust output matches the C Swiss Ephemeris library bitwise-exactly (or within epsilon for iterative functions).

## How it works

1. A C program (`tools/golden_gen.c`) links against `../swisseph/libswe.a` and calls C functions with curated edge-case inputs
2. It outputs JSON golden data files to `tests/golden-data/`
3. Rust integration tests (`tests/golden/`) deserialize the JSON and compare against our implementations

## Structure

```
tools/
  golden_gen.c      # C generator — one gen_<module>() function per module
  Makefile           # builds golden_gen, runs it to produce JSON

tests/
  golden/
    main.rs          # test crate root — shared assertion helpers
    date.rs          # date module golden tests
  golden-data/
    date.json        # generated fixture data (checked into git)
```

## Adding a new module

1. **C side**: add a `gen_<module>()` function in `golden_gen.c` and a new `argv[1]` case in `main()`. Add a Makefile target for the new JSON file.
2. **Generate**: `cd tools && make generate`
3. **Rust side**: add `tests/golden/<module>.rs` with typed `#[derive(Deserialize)]` case structs, a `load()` helper, and one `#[test]` per function. Add `mod <module>;` to `tests/golden/main.rs`.

## Assertion helpers (in `tests/golden/main.rs`)

- `assert_f64_exact(label, expected, actual)` — bitwise comparison; reports hex bits on failure. Use for pure-math functions (julday, revjul).
- `assert_f64_eps(label, expected, actual, eps)` — epsilon comparison; reports diff magnitude. Use for iterative/convergent functions.

## Regenerating golden data

Requires the C library built at `../swisseph/`:

```sh
cd ../swisseph && make libswe.a   # if not already built
cd ../swisseph.rs/tools && make generate
```

Golden data files are checked into git so CI doesn't need the C toolchain. Regenerate when:
- Adding new test vectors or a new module
- The C library version changes

## Current coverage

| Module | Functions | Cases | Comparison |
|--------|-----------|-------|------------|
| date   | julday, revjul, date_conversion, day_of_week | 292 | bitwise-exact |

## Deferred

- `utc_to_jd` / `jdet_to_utc` — these depend on a `DeltaT` provider. Golden testing deferred until our delta-T implementation matches the C library's `swe_deltat`. The JSON structure and Rust test stubs will be added at that point.

## Edge case coverage strategy

Cases are curated to hit every code branch rather than brute-force sweeps:
- Negative years, year 0, far past/future (-5000 to +5000)
- Gregorian/Julian calendar switch (Oct 1582)
- Century correction for negative years (the `u/100 == floor(u/100)` branch in julday)
- Leap year boundaries (2000, 1900, 1600, -4)
- Month < 3 / >= 3 paths
- revjul Gregorian correction boundary (JD 1830691.5)
- Invalid dates for date_conversion (Feb 30, Gregorian gap, month 13)
