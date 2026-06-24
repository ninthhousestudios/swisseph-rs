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

### C Source Reference
- Original C repo: `../swisseph/`
- C catalogues: `../swisseph/claude/catalogue-{public,internal}.md`
- Ephemeris data: `../swisseph/ephe/`
- C test suite: `../swisseph/setest/`
