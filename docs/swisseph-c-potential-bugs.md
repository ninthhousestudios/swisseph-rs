# Potential bugs in upstream Swiss Ephemeris (C)

Observations from porting. None of these have practical impact at the library's
intended precision, but they're documented for completeness.

## 1. embofs_mosh uses cached obliquity for backward-difference time point

**Location:** `swemplan.c`, `embofs_mosh()` (line ~416)

**What:** `embofs_mosh` reads `swed.oec.seps` / `swed.oec.ceps` — the obliquity
cached from the most recent `swe_calc` call. When `swi_moshplan` computes Earth
velocity by backward difference, it calls `embofs_mosh` twice:

```c
embofs_mosh(tjd, xe);                        // main position
embofs_mosh(tjd - PLAN_SPEED_INTV, x2);     // backward-difference position
```

Both calls read the same `swed.oec` (obliquity at `tjd`). The second call should
use the obliquity at `tjd - PLAN_SPEED_INTV` (0.0001 days earlier).

**Impact:** Over 0.0001 days, obliquity changes by ~6e-16 radians. Applied to the
Moon at ~4.3e-5 AU, position error is ~2.6e-20 AU (~4 picometers). The short Moon
series in `embofs_mosh` is accurate to ~1 arcminute, so this error is 12 orders of
magnitude below the series precision. Zero practical impact.

**Cause:** Global state architecture — `embofs_mosh` is a `static` function that
reads from the `swed` global rather than taking obliquity as a parameter. No
indication this is intentional; no comment suggesting deliberate parameter-holding
for backward-difference stability.

**Our Rust code:** Uses correct per-time-point obliquity. Tests pass at 1e-10
tolerance since the effect is far below that threshold.
