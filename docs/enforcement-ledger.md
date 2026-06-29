# Enforcement Ledger

Seeded from PRD (yojana:swisseph-rs) on 2026-06-24.

Every architectural claim is routed to exactly one enforcement mechanism.
A claim with no row is a routing bug.

## Claims

| # | Claim | Source | Bucket | Mechanism | Status |
|---|---|---|---|---|---|
| 1 | Eliminate global state (explicit context struct) | project desc | (c) | CLAUDE.md: stateless design | live |
| 2 | Ephemeris has no mutable fields; methods take &self | phase-1-considerations.md | (c) | CLAUDE.md: stateless design | live |
| 3 | Completely stateless library for thread safety | phase-1-considerations.md | (c) | CLAUDE.md: stateless design | live |
| 4 | Ephemeris must be Send + Sync | phase-1-considerations.md | (d) | compile-time assert | test (add with Phase 1) |
| 5 | Continuous testing against C swetest for numerical fidelity | project desc | (d) | golden-data tests per phase | test (per-phase) |
| 6 | TryFrom for all C integer-to-enum conversions | phase-1-considerations.md | (c) | CLAUDE.md: API patterns | live |
| 7 | Error enum replaces C serr string buffers | phase-1-considerations.md | (c) | CLAUDE.md: API patterns | live |
| 8 | Warning handling via flags_used, no Warning type | decision swisseph-rs/1 | (c) | CLAUDE.md: API discipline | live |
| 9 | Moshier backend first, file backends later | project handoff | (c) | CLAUDE.md: process | live |
| 10 | Each arc phase: study, implement, verify | project handoff | (c) | CLAUDE.md: process | live |
| 11 | Backends (moshier, jpl, sweph_file) independent | module tree design | (a) | sutra: forbidden_dep x6 | live — sweph↔moshier bound; jpl globs pre-pointed to `src/jpl/**`, inert until swisseph-rs/43 (re-pointed checkpoint:swisseph-rs/8) |
| 12 | App modules go through calc, not backends | module tree design | (b) | sutra: forbidden_dep app->backend | deferred: first app module impl (Phase 8+) — re-keyed checkpoint:swisseph-rs/8 (calc dispatcher exists since Phase 5, but app modules are empty stubs; binding now = inert) |
| 13 | App modules independent of each other | module tree design | (b) | sutra: forbidden_dep app<->app | deferred: Phase 8+ |
| 14 | No dependency cycles in module graph | skill default | (b) | sutra: no_cycles scope src/ | live — bound checkpoint:swisseph-rs/8 (2026-06-29), 35 files, 0 violations |
| 15 | Module tree mirrors C file structure | scaffold task | (c) | CLAUDE.md: module structure | live |
| 16 | EphemerisConfig with Default, no builder | decision swisseph-rs/1 | (c) | CLAUDE.md: API discipline | live |

## Maintenance

- **Guard from first commit**: blocking constraints are live the moment code exists.
- **Per-task review**: `sutra_review` runs per task; bucket (c) items are review checklist material.
- **First review checkpoint**: run vidhi-sutra-tend to fire due triggers, run initial convention triage, set fan-in guardrails, constrain module interiors.
- **New track PRDs**: re-run vidhi-sutra-seed additively (append rows, never regenerate).

## Tend passes

### 2026-06-29 — checkpoint:swisseph-rs/8 (JPL backend planning)

- **Drift repaired**: the four `backend-isolation:*→sweph` / `sweph→*` globs targeting `src/sweph_file.rs` had gone silently inert (`dead_constraint`) when sweph_file became a directory module at an earlier phase — silent unenforcement. Re-pointed to `src/sweph_file/**`; `moshier↔sweph` now bind.
- **Pre-point for swisseph-rs/43**: jpl is converting `src/jpl.rs` → `src/jpl/` in that task. The four jpl-touching globs were re-pointed to `src/jpl/**` ahead of the conversion. They are intentionally inert now (empty stub, nothing to guard) and bind automatically once 43 lands. A combined `src/jpl{.rs,/**}` brace was attempted first but globset rejects it (recursive wildcard must be a lone path component).
- **Trigger fired**: row 14 `no_cycles` (deferred since "first phase with real imports", Phase 1 — six phases overdue). Bound at `scope = src/`, blocking; 35 files, 0 violations. Caveat: sutra's graph captures `use`-based edges only, so the pre-existing `context`↔`deltat` coupling (via `crate::`-qualified calls) is not seen and is out of scope.
- **Trigger re-keyed, not fired**: row 12 app→backend. Its "Phase 5 calc dispatcher exists" trigger is technically met, but houses/eclipse/phenomena/heliacal/stars/ayanamsa are empty stubs (zero import edges) — binding the 18 rules now would only add inert `dead_constraint` noise. Re-keyed to "first app module implementation (Phase 8+)".
- **Not added**: `max_fan_in` — analysis tier disabled; manual fan-in read shows only foundational modules (types/constants/math/flags/error) high, all legitimate infrastructure; `calc.rs` is the intended dispatcher (high fan-out, low fan-in). No emergent non-infra hub.
- **Conventions**: none triaged — no FCA/conventions data this pass.
