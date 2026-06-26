<c_source_budget>
HARD CONSTRAINT: Do not read the C source (`../swisseph/`) directly during planning or implementation.
For each module being ported, a C reference doc exists at `docs/c-ref-{module}.md` containing
algorithm details, coefficients, loop structures, unit conversions, and line numbers — everything
needed to implement without re-reading the original.

If no C ref doc exists for the module you're working on, launch a single dedicated agent
(model: sonnet) to read the C source and produce the ref doc. The agent writes
`docs/c-ref-{module}.md`; the planning/implementation agent never touches the C files. This
isolates the 50-80k tokens of C reading into a disposable context. Brief the agent to:
- The ref doc agent must NOT spawn subagents — it does all reading and writing itself
  in one context (catalogue lookup, C source reading, doc writing)
- Extract: algorithm structure, coefficients/tables, loop logic, unit conversions,
  boundary handling, line numbers for traceability
- Follow the format of an existing sibling ref doc (point it at e.g. `docs/c-ref-nutation.md`)

The C catalogues (`../swisseph/claude/catalogue-{public,internal}.md`) are still fine for
high-level function discovery. The prohibition is on reading the actual `.c` / `.h` implementation
files during plan or implement phases.

**Exception — FP fidelity debugging**: When golden tests fail and the cause is floating-point
evaluation order (not algorithm bugs), read the C source directly. Ref docs capture algorithm
structure but not character-level expression order (`+=` vs `= x +`, multiplication grouping),
which is exactly what matters for FP fidelity. See `docs/golden-testing.md` § "Debugging FP
fidelity failures" for the instrument-first protocol.
</c_source_budget>

<codebase_map>
`docs/codebase-map.md` — module layout, key types with line numbers, golden test patterns,
insertion points. Read this file BEFORE launching any exploration agents or broad codebase
searches. It exists to prevent expensive exploration sweeps. Agents that skip it and re-discover
what the map already documents are wasting tokens. Update after landing new modules.
</codebase_map>

<architecture>
## Stateless Design
- `Ephemeris` holds only read-only configuration. No mutable cache, no internal state.
- All methods on `Ephemeris` take `&self`, never `&mut self`.
- The calculation pipeline is pure: inputs → math → output, no side effects.

## API Patterns
- All C integer-to-enum conversions use `TryFrom`. Never transmute or cast blindly.
- Error handling via `Error` enum. No string buffer passing patterns from C.
- Warning/fallback signaling via `CalcResult.flags_used` — compare requested vs used flags. No separate `Warning` type.
- Construction via `EphemerisConfig` struct with `Default`. No builder pattern.

## Module Structure
- Module tree mirrors C file structure. New functionality goes in the corresponding module.
- Three ephemeris backends (`moshier/`, `jpl.rs`, `sweph_file.rs`) are independent — no cross-dependencies. Enforced by sutra `forbidden_dep` constraints.
- Application modules (`houses`, `eclipse`, `phenomena`, etc.) go through `calc.rs`, not directly to backends.
</architecture>

<process>
- Moshier backend first (self-contained, no file I/O). File backends after core pipeline.
- Each arc phase: study the C source → implement in Rust → verify against C `swetest` golden data.
- Numerical fidelity is non-negotiable. Every phase includes golden-data comparison tests.
</process>

<testing>
- Golden data harness: `docs/golden-testing.md` — how to add tests for new modules, regenerate data, assertion patterns.
- Every module must have golden-data differential tests against the C library. Bitwise-exact for pure math, epsilon for iterative functions.
</testing>

<c_source_reference>
- Original C repo: `../swisseph/`
- C catalogues: `../swisseph/claude/catalogue-{public,internal}.md`
- C reference docs: `docs/c-ref-*.md` (per-module porting references — use these, not the raw C)
- Ephemeris data: `../swisseph/ephe/`
- C test suite: `../swisseph/setest/`
</c_source_reference>
