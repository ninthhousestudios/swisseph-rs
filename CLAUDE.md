<c_source_budget>
## When to use what

| Phase | Source | Why |
|---|---|---|
| **Planning / implementation** | `docs/c-ref-{module}.md` | Algorithm structure, coefficients, loop logic — everything needed to write the Rust port without re-reading the C |
| **High-level function discovery** | `../swisseph/claude/catalogue-{public,internal}.md` | Finding which C functions exist and what they do |
| **Golden test failures — algorithm bugs** | C source directly | Ref docs capture structure but omit details like global state interactions, implicit zeroing, which struct fields are actually populated at a given call site. When the test failure is a logic/algorithm bug (errors > ~1e-8), read the C source — the ref doc's abstraction hides exactly the details you need. |
| **Golden test failures — FP fidelity** | C source directly | Ref docs don't capture expression evaluation order (`+=` vs `= x +`, multiplication grouping). See `docs/golden-testing.md` § "Debugging FP fidelity failures". |

## Ref doc creation

Do not read the C source (`../swisseph/`) directly during planning or implementation.
For each module being ported, a C reference doc exists at `docs/c-ref-{module}.md`.

If no C ref doc exists for the module you're working on, launch a single dedicated agent
(model: sonnet) to read the C source and produce the ref doc. The agent writes
`docs/c-ref-{module}.md`; the planning/implementation agent never touches the C files. This
isolates the 50-80k tokens of C reading into a disposable context. Brief the agent to:
- The ref doc agent must NOT spawn subagents — it does all reading and writing itself
  in one context (catalogue lookup, C source reading, doc writing)
- Extract: algorithm structure, coefficients/tables, loop logic, unit conversions,
  boundary handling, line numbers for traceability
- Follow the format of an existing sibling ref doc (point it at e.g. `docs/c-ref-nutation.md`)

## The key lesson

Ref docs are a **planning and implementation** tool. They save you from re-reading 50k tokens
of C to understand an algorithm's structure. But they are an abstraction, and abstractions hide
details. When debugging golden test failures, the hidden details (global state, implicit zeros,
expression order, which fields are populated when) are often the exact cause. Don't reason from
the ref doc about what the C "should" do — read the C and see what it actually does.
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

<stateless_tolerance>
## Stateless vs Stateful: Known Precision Boundaries

The C Swiss Ephemeris uses global mutable state (cached planetary positions, open file handles)
that subtly affects intermediate computations. Our stateless Rust port produces slightly different
results in two specific areas. These are NOT bugs — they are inherent consequences of stateless
architecture and are astronomically negligible.

### 1. Deflection speed (< 0.06 milliarcseconds)
The C `swi_deflect_light` reads the Sun's position from a global cache (`psdp->x`) populated
earlier in the same `swe_calc` call. The stateless Rust version constructs the deflection
geometry from explicitly-passed parameters. The resulting speed values differ by up to ~1e-7
degrees/day (0.06 mas) — roughly 500,000x below any practical astrological significance.

Golden test tolerance: 1e-7 degrees for speed components, 1e-10 for positions.

### 2. SPEED3 at SE1 file boundaries
SPEED3 evaluates three positions (t-dt, t, t+dt) within a single `swe_calc` call. C's stateful
file caching means the three evaluations can use different .se1 files (the first evaluation opens
a file that stays cached for subsequent ones). Stateless Rust independently selects files for each
evaluation. At file boundaries (e.g. jd=2378496.5 = sepl_18's tfstart), this produces large
SPEED3 divergence — up to 0.5 degrees for Moon. This is a C state artifact, not an error.

Golden test tolerance: relaxed at file-boundary epochs for SPEED3.

### Do not chase sub-milliarcsecond deflection speed differences
Four debugging sessions were spent on swisseph-rs/41 trying to exactly match C's deflection
speed output. The root cause is that C's global state produces a slightly different geometric
construction than stateless Rust. The difference is astronomically meaningless. Accept the
1e-7 tolerance and move on. (See task swisseph-rs/41 execution_record for full post-mortem.)
</stateless_tolerance>

<c_source_reference>
- Original C repo: `../swisseph/`
- C catalogues: `../swisseph/claude/catalogue-{public,internal}.md`
- C reference docs: `docs/c-ref-*.md` (per-module porting references — use these, not the raw C)
- Ephemeris data: `../swisseph/ephe/`
- C test suite: `../swisseph/setest/`
</c_source_reference>
