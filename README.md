# swisseph

Pure-Rust, stateless port of the [Swiss Ephemeris](https://www.astro.com/swisseph/)
astronomical calculation library.

Computes planetary and lunar positions, house cusps, eclipses, occultations, rise/set
times, heliacal events, nodes/apsides, orbital elements, longitude crossings, and fixed
stars. Bit-compatible with upstream C **version 2.10.03** — the golden-test oracle is
a vendored 2.10.03 checkout; upstream v3 is in development and may change quirks this
port deliberately replicates (see `docs/swisseph-c-potential-bugs.md`).

## Installation

```toml
[dependencies]
swisseph = "0.1"
```

### Feature flags

| Feature | Default | Effect |
|---------|---------|--------|
| `swisseph-files` | yes | Enables the Swiss Ephemeris `.se1` file reader |
| `jpl` | yes | Enables the JPL DE ephemeris reader |
| `serde` | no | Derives `Serialize`/`Deserialize` on all public types |
| `star-tools` | no | Builds the `make-swe-stars` binary (catalog maintenance) |

Disable both `swisseph-files` and `jpl` for a zero-file-IO, pure-Moshier build:

```toml
swisseph = { version = "0.1", default-features = false }
```

The Moshier backend is always available and requires no data files. For higher precision
(especially for the Moon and outer planets beyond ~3000 BCE / 3000 CE), you need Swiss
Ephemeris data files (`.se1`) or a JPL DE file (e.g. `de441.eph`). These are available
from [Astrodienst's download page](https://www.astro.com/ftp/swisseph/ephe/).

## Quick start

The Moshier backend is self-contained — no data files needed:

```rust
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig};

let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();
let jd_ut = 2451545.0; // J2000.0
let result = eph.calc_ut(jd_ut, Body::Sun, CalcFlags::SPEED).unwrap();
let longitude = result.data[0];
assert!((280.0..=281.0).contains(&longitude));
```

With Swiss Ephemeris files for higher precision:

```rust
use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig};
use swisseph::types::EphemerisSource;
use std::path::PathBuf;

let config = EphemerisConfig {
    ephemeris_source: EphemerisSource::Swiss,
    ephe_path: Some(PathBuf::from("/path/to/ephe")),
    ..Default::default()
};
let eph = Ephemeris::new(config).unwrap();
let result = eph.calc_ut(2451545.0, Body::Moon, CalcFlags::SPEED).unwrap();
```

## Design: stateless vs the C original

This is the architectural reason this port exists. The C Swiss Ephemeris works; it is
numerically excellent. But its internal design — a process-global mutable struct — makes
it painful to use in modern concurrent applications.

### How C works

The C library stores all state in one process-global `swed` struct:

- **Configuration via mutation**: `swe_set_ephe_path`, `swe_set_topo`, `swe_set_sid_mode`,
  `swe_set_tid_acc` — each mutates global state that subsequent `swe_calc` calls read.
- **Cached file handles and positions**: open `.se1` files, last-computed planetary
  positions, and nutation values persist across calls.
- **Consequences**: not thread-safe, one configuration per process, and bindings to
  higher-level languages must serialize all calls or maintain one copy of the shared
  library per thread/isolate. The
  [swisseph.dart](https://github.com/nickvdyck/swisseph.dart) package, for example,
  loads a separate copy of the native library per Dart isolate — the direct motivation
  for this stateless port.

### How this port works

- `EphemerisConfig` is plain data, immutable after construction.
- `Ephemeris` holds only read-only config + read-only mmap'd file handles.
- Every method takes `&self`, never `&mut self`.
- The calculation pipeline is pure: inputs → math → output, no side effects.
- Per-call configuration variance (ephemeris source, topographic position) goes through
  explicit parameters or `CalcFlags` bits, not set-then-call.
- `Ephemeris` is `Send + Sync` — one instance shared across any number of threads with
  zero synchronization.

### What statelessness costs

Three documented precision boundaries arise from the absence of C's global caches.
All are astronomically negligible:

1. **Deflection speed** (< 0.06 milliarcseconds): C reads a cached Sun position for
   light-deflection geometry; stateless Rust constructs it from explicit parameters.
2. **SPEED3 at `.se1` file boundaries**: C's file cache means the three evaluations
   can use different files; stateless Rust independently selects per evaluation.
3. **Moshier osculating node/apogee speed** (< 4e-6 deg/day): C's global
   obliquity/nutation cache rounds differently than clean per-epoch recomputation.

Details in `CLAUDE.md` § "Stateless vs Stateful: Known Precision Boundaries".

## Numerical fidelity

Every module has golden differential tests against the C oracle — thousands of cases
covering all backends (Moshier, Swiss, JPL), flag combinations, and edge cases. Pure
math is bitwise-exact; iterative searches use documented epsilons. See
`docs/golden-testing.md`.

## API mapping

The most-used C functions and their Rust equivalents:

| C function | Rust method |
|-----------|-------------|
| `swe_calc_ut` | `Ephemeris::calc_ut` |
| `swe_calc` | `Ephemeris::calc` |
| `swe_houses_ex2` | `Ephemeris::houses_ex2` |
| `swe_fixstar2_ut` | `Ephemeris::fixstar2_ut` |
| `swe_rise_trans` | `Ephemeris::rise_trans` |
| `swe_azalt` | `Ephemeris::azalt` |
| `swe_sol_eclipse_when_glob` | `Ephemeris::sol_eclipse_when_glob` |
| `swe_lun_eclipse_when` | `Ephemeris::lun_eclipse_when` |
| `swe_heliacal_ut` | `Ephemeris::heliacal_ut` |
| `swe_pheno_ut` | `Ephemeris::pheno_ut` |
| `swe_nod_aps_ut` | `Ephemeris::nod_aps_ut` |
| `swe_get_orbital_elements` | `Ephemeris::get_orbital_elements` |
| `swe_julday` | `date::julday` |
| `swe_revjul` | `date::revjul` |
| `swe_sidtime` | `sidereal_time::sidereal_time` |
| `swe_house_pos` | `houses::house_pos` |
| `swe_split_deg` | `math::split_degrees` |
| `swe_cotrans` | `math::cotrans` |
| `swe_gauquelin_sector` | `Ephemeris::gauquelin_sector` |

Every method carries a `#[doc(alias = "swe_...")]` attribute, so
[docs.rs](https://docs.rs) search works with C function names.

## Not ported

| Feature | Reason |
|---------|--------|
| EP4 compressed ephemeris reader | Not planned — standard `.se1` files cover all use cases |
| `swe_set_interpolate_nut` | Inherently a stateful cache optimization — intentionally unsupported |
| `swe_get_current_file_data` | Stateful / not applicable to this design |
| `swe_get_library_path` | Stateful / not applicable |
| `swe_fixstar` (v1 API) | Deprecated upstream; v2 (`fixstar2`) is ported |

## Examples

Four runnable examples in `examples/` (all Moshier-only, no data files needed):

```sh
cargo run --example natal           # 10 planets + Placidus cusps
cargo run --example rise_set        # sunrise/sunset/moonrise/moonset
cargo run --example eclipse_search  # next 3 solar + lunar eclipses
cargo run --example sidereal        # tropical vs Lahiri side by side
```

## Performance

| Workload | Rust | C | Ratio |
|----------|------|---|-------|
| Full chart, Swiss files (10 planets) | 102 µs | 40 µs | 2.6× |
| Placidus houses | 5.4 µs | 9.6 µs | **0.56×** |
| Solar eclipse search | 2.84 ms | 2.61 ms | 1.09× |

Thread scaling (one shared `&Ephemeris`, Moshier, 10k charts):
1 thread → 5,912 charts/sec, 4 threads → 20,120 charts/sec (**3.4× scaling**),
8 threads → 24,691 charts/sec. Zero synchronization — `Ephemeris` is `Send + Sync`.

The Moshier backend's single-call overhead is higher than C (C caches Earth/nutation
globally across calls), but scales linearly with threads where C cannot.
See [`docs/benchmarks.md`](docs/benchmarks.md) for methodology and analysis.

## License

This crate is a derivative work of the Swiss Ephemeris.

**Original C source:**
Copyright (C) 1997 - 2021 Astrodienst AG, Switzerland.
Authors: Dieter Koch and Alois Treindl.
The original is dual-licensed under AGPL-3.0-or-later or the Swiss Ephemeris
Professional License.

**This Rust port** is distributed under the
[GNU Affero General Public License v3.0 or later](LICENSE).
