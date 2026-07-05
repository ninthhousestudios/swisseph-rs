# swisseph-ffi

C ABI wrapper for the `swisseph` Rust ephemeris library. Produces a shared
library (`.so`/`.dylib`/`.dll`) and a static library (`.a`) with a generated
C header (`include/swisseph.h`).

## Building

```sh
cargo build -p swisseph-ffi --release
```

Artifacts:

| Platform | Shared library | Static library |
|---|---|---|
| Linux | `target/release/libswisseph_ffi.so` | `target/release/libswisseph_ffi.a` |
| macOS | `target/release/libswisseph_ffi.dylib` | `target/release/libswisseph_ffi.a` |
| Windows | `target/release/swisseph_ffi.dll` | `target/release/swisseph_ffi.lib` |

The generated header is `swisseph-ffi/include/swisseph.h`. Copy it alongside
the library for C/Dart consumers.

## Linking example

```sh
cc my_app.c -I swisseph-ffi/include -L target/release -lswisseph_ffi -lm -lpthread -ldl -o my_app
LD_LIBRARY_PATH=target/release ./my_app
```

## Handle model

All computation goes through an opaque `SweEphemeris*` handle:

```c
SweConfig config;
swisseph_config_default(&config);
// config.ephemeris_source = 1;  /* Swiss */
// config.ephe_path = "/path/to/ephe";

SweEphemeris *eph = NULL;
char err[256];
int ret = swisseph_new(&config, &eph, err, sizeof(err));
if (ret != 0) { /* handle error */ }

double xx[6];
int flags_used;
ret = swisseph_calc_ut(eph, 2451545.0, 0 /* Sun */, 256 /* SPEED */,
                       NULL, NULL, xx, &flags_used, err, sizeof(err));

swisseph_free(eph);  /* NULL-safe */
```

The handle is **immutable after construction** -- there are no `set_*` mutators.
Per-call variance (topographic position, sidereal mode) is passed as explicit
parameters to each function call.

## Isolate/thread sharing (Dart pattern)

Because `Ephemeris` is `Send + Sync` with no interior mutability, a single
handle can be shared across threads (or Dart isolates) with zero synchronization:

```c
// Create once
SweEphemeris *eph = /* swisseph_new(...) */;
uintptr_t handle_int = (uintptr_t)eph;

// In each thread/isolate: recover the pointer from the integer
SweEphemeris *eph = (SweEphemeris*)(uintptr_t)handle_int;
swisseph_calc_ut(eph, ...);

// Free once, from any single thread, after all others are done
swisseph_free(eph);
```

This is tested by `tests/concurrency.rs`: 8 threads calling `swisseph_calc_ut`
through a `usize`-round-tripped pointer produce results bitwise-identical to
serial computation.

## Return-code conventions

Three conventions are used, matching the C Swiss Ephemeris where possible:

| Convention | Functions | Success | Failure |
|---|---|---|---|
| 0 = OK | Most functions (`calc_ut`, `houses`, `new`, ...) | `0` | Negative `SweErrorCode` |
| Positive flags | Eclipse/occultation functions (`sol_eclipse_where`, ...) | Positive `EclipseFlags` bitmask | Negative `SweErrorCode`; `0` = no eclipse |
| Vision flags | `vis_limit_mag` | `0`=photopic, `1`=scotopic, `2`=mixed | `-2` = below horizon; other negative = error |

All functions that can fail accept an optional `err_buf`/`err_cap` pair for the
error message. Pass `NULL`/`0` to suppress messages. Error messages are
NUL-terminated and silently truncated at a UTF-8 boundary when the buffer is
too small.

## Thread safety

The handle is safe to use from multiple threads concurrently without any
locking. This is guaranteed by the Rust type system (`Ephemeris: Send + Sync`)
and verified at runtime by `tests/concurrency.rs` (bitwise determinism) and
`tests/send_sync.rs` (compile-time trait assertion + runtime determinism).

## License

AGPL-3.0-or-later, same as the underlying `swisseph` crate and the original
C Swiss Ephemeris library.
