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
| 11 | Backends (moshier, jpl, sweph_file) independent | module tree design | (a) | sutra: forbidden_dep x6 | live (dead_constraint until backends have imports) |
| 12 | App modules go through calc, not backends | module tree design | (b) | sutra: forbidden_dep app->backend | deferred: Phase 5 (calc dispatcher) |
| 13 | App modules independent of each other | module tree design | (b) | sutra: forbidden_dep app<->app | deferred: Phase 8+ |
| 14 | No dependency cycles in module graph | skill default | (b) | sutra: no_cycles scope src/ | deferred: first phase with real imports |
| 15 | Module tree mirrors C file structure | scaffold task | (c) | CLAUDE.md: module structure | live |
| 16 | EphemerisConfig with Default, no builder | decision swisseph-rs/1 | (c) | CLAUDE.md: API discipline | live |

## Maintenance

- **Guard from first commit**: blocking constraints are live the moment code exists.
- **Per-task review**: `sutra_review` runs per task; bucket (c) items are review checklist material.
- **First review checkpoint**: run vidhi-sutra-tend to fire due triggers, run initial convention triage, set fan-in guardrails, constrain module interiors.
- **New track PRDs**: re-run vidhi-sutra-seed additively (append rows, never regenerate).
