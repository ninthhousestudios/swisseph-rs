# Benchmarks

## Methodology

- **Hardware**: AMD Ryzen 5 7535U (6 cores / 12 threads), laptop
- **Rust**: rustc 1.96.1, `cargo bench` (release profile, criterion 0.5)
- **C**: GCC 16.1.1 (`-O2`), linked against vendored `libswe.a`
- **Workloads**: identical inputs (J2000.0, same bodies/coordinates)
- **Exclusions**: first-call file-open costs excluded in both (warm-up phase)

## Rust vs C — single-thread

| Workload | Rust | C | Ratio |
|----------|------|---|-------|
| Full chart, Moshier (10 planets) | 150 µs | 3.4 µs | **44×** |
| Full chart, Swiss files (10 planets) | 102 µs | 40 µs | 2.6× |
| Placidus houses (Zurich) | 5.4 µs | 9.6 µs | **0.56×** |
| Solar eclipse search (next from J2000) | 2.84 ms | 2.61 ms | 1.09× |

### Why Moshier is 44× slower

This is the expected cost of stateless architecture.  C's global `swed` struct caches
the Earth position, nutation, and obliquity across all 10 `swe_calc` calls within a
chart.  The Rust version recomputes them from scratch on every call — by design, since
`Ephemeris` is `Send + Sync` with no mutable state.

The 44× gap is specific to the "compute 10 bodies at the same epoch" pattern. For
single-body calls or different epochs, the caching advantage disappears.  And as the
thread-scaling section shows, the architectural payoff arrives the moment you need
concurrency.

For the Swiss Ephemeris file backend, the gap narrows to 2.6× because file I/O
(Chebyshev coefficient reads from mmap'd `.se1` files) dominates, and both
implementations do roughly the same work per body.

Houses and eclipse search show Rust at parity or faster — these are compute-heavy
paths where Rust's optimizer excels and C's caching provides no advantage.

## Micro-benchmarks (hot kernels)

| Kernel | Time |
|--------|------|
| `nutation::nutation` (IAU 2000B) | 1.6 µs |
| `moshier::moon::moshmoon2` (lunar series) | 2.3 µs |

These are the profile targets. The Moon series evaluation alone accounts for ~1.5% of
a full Moshier chart, and nutation ~1.1%.  The dominant cost in a full Moshier chart
is the repeated Earth computation (needed for every body's aberration/deflection
correction) — a per-call cache for the current epoch's Earth position would eliminate
most of the 44× gap without sacrificing thread safety.

## Thread scaling (Moshier, 10,000 epochs × 10 bodies)

| Threads | Wall time | Speedup | Charts/sec |
|---------|-----------|---------|------------|
| 1 | 1.69 s | 1.0× | 5,912 |
| 2 | 881 ms | 1.92× | 11,351 |
| 4 | 497 ms | 3.40× | 20,120 |
| 8 | 405 ms | 4.17× | 24,691 |

Near-linear scaling through the physical core count (6 cores). The 8-thread result
shows diminishing returns from hyperthreading. No locks, no contention — one shared
`&Ephemeris` across all threads with zero synchronization.

**The architectural payoff**: a C consumer needing thread safety must either:
- serialize all `swe_calc` calls behind a mutex (1× throughput), or
- maintain one copy of the library + its global state per thread (6× memory, complex lifecycle)

The Rust port achieves 4× throughput with one shared instance, zero synchronization,
and `Arc<Ephemeris>` costing nothing beyond the reference count at clone time.

## Honesty notes

1. **C timings exclude locking**: a real multi-threaded C consumer would add mutex
   overhead to every call.  The single-thread C numbers are a best-case that no
   concurrent C program can achieve without per-thread library copies.

2. **First-call costs excluded**: both Rust and C timings use a warm-up phase. The
   first `calc_ut` with Swiss files is slower (mmap fault / file header parse).

3. **Profile hint for Moshier gap**: `cargo flamegraph --bench ephemeris -- --bench
   "calc_moshier"` will show the repeated Earth/nutation computation. An epoch-scoped
   computation cache (interior `Cell`-based or a separate `EpochCache` struct passed
   per chart) is the obvious optimization — a separate task.
