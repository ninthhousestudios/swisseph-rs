## Architecture Invariants

### Stateless Design
- `Ephemeris` holds only read-only configuration. No mutable cache, no internal state.
- All methods on `Ephemeris` take `&self`, never `&mut self`.
- The calculation pipeline is pure: inputs -> math -> output, no side effects.

### API Patterns
- All C integer-to-enum conversions use `TryFrom`. Never transmute or cast blindly.
- Error handling via `Error` enum. No string buffer passing patterns from C.
- Warning/fallback signaling via `CalcResult.flags_used` -- compare requested vs used flags. No separate `Warning` type.
- Construction via `EphemerisConfig` struct with `Default`. No builder pattern.

### Module Structure
- Module tree mirrors C file structure. New functionality goes in the corresponding module.
- Three ephemeris backends (`moshier/`, `jpl.rs`, `sweph_file.rs`) are independent -- no cross-dependencies. Enforced by sutra `forbidden_dep` constraints.
- Application modules (`houses`, `eclipse`, `phenomena`, etc.) go through `calc.rs`, not directly to backends.

### Process
- Moshier backend first (self-contained, no file I/O). File backends after core pipeline.
- Each arc phase: study the C source -> implement in Rust -> verify against C `swetest` golden data.
- Numerical fidelity is non-negotiable. Every phase includes golden-data comparison tests.

### Testing
- Golden data harness: `docs/golden-testing.md` — how to add tests for new modules, regenerate data, assertion patterns.
- Every module must have golden-data differential tests against the C library. Bitwise-exact for pure math, epsilon for iterative functions.

### C Source Reference
- Original C repo: `../swisseph/`
- C catalogues: `../swisseph/claude/catalogue-{public,internal}.md`
- Ephemeris data: `../swisseph/ephe/`
- C test suite: `../swisseph/setest/`

### Codebase Map
- `docs/codebase-map.md` — module layout, key types with line numbers, golden test patterns, insertion points. Consult before exploring; update after landing new modules.
