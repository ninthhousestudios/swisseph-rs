# C Reference: Asteroid Calculation Path — sweph.c

Porting reference for the asteroid-specific parts of `swe_calc()`: dispatch, filename
generation, MPC orbital-element parsing, and downstream H/G/diameter consumption.

**Prerequisite reading** (not duplicated here):
- `docs/c-ref-se1-file.md` — `.se1` binary format, `read_const()` steps 1–3/5–14,
  `do_fread`, `get_new_segment()`, `rot_back()`, `sweph()`'s segment-evaluation pipeline
  (including its asteroid heliocentric→barycentric step and `SEI_ANYBODY` slot aliasing),
  `swi_echeb`/`swi_edcheb`.
- `docs/c-ref-calc.md` — `swe_calc()`, `swecalc()` main-planet branches, `main_planet()`
  ephemeris cascade, `sweplan()`, `app_pos_etc_plan()` (generic light-time/aberration/
  deflection/precession pipeline shared by planets and asteroids).

This doc covers only what those two omit: the `swecalc()` minor-planet branch, asteroid
file naming, the MPC elements-line column arithmetic, and where `ast_H`/`ast_G`/`ast_diam`
end up being read.

---

## 1. `swecalc()` — Minor Planet Branch (sweph.c:1016–1101)

### 1.1 Branch condition and internal-ID mapping (sweph.c:1019–1042)

```c
} else if (ipl == SE_CHIRON || ipl == SE_PHOLUS || ipl == SE_CERES
    || ipl == SE_PALLAS || ipl == SE_JUNO || ipl == SE_VESTA
    || ipl > SE_PLMOON_OFFSET
    || ipl > SE_AST_OFFSET  // obsolete after previous condition
    ) {
  /* internal planet number */
  if (ipl < SE_NPLANETS) {
    ipli = pnoext2int[ipl];
  } else if (ipl <= SE_AST_OFFSET + MPC_VESTA && ipl > SE_AST_OFFSET) {
    ipli = SEI_CERES + ipl - SE_AST_OFFSET - 1;
    ipl  = SE_CERES  + ipl - SE_AST_OFFSET - 1;
  } else {                    /* any asteroid except*/
    ipli = SEI_ANYBODY;
  }
  if (ipli == SEI_ANYBODY) {
    ipli_ast = ipl;
  } else {
    ipli_ast = ipli;
  }
  pdp = &swed.pldat[ipli];
  xp = pdp->xreturn;
```

Three cases, resolved in order:

| Public `ipl` | `ipli` (internal slot) | `ipli_ast` (passed to `sweph()`) | Notes |
|---|---|---|---|
| `SE_CHIRON..SE_VESTA` (15–20) | `pnoext2int[ipl]` = `SEI_CHIRON..SEI_VESTA` (12–17) | = `ipli` | `pnoext2int[]` table, sweph.c:182 |
| `SE_AST_OFFSET+1 .. SE_AST_OFFSET+4` (10001–10004) | `SEI_CERES + (ipl-SE_AST_OFFSET-1)` | = `ipli` | **`ipl` itself is rewritten** to `SE_CERES + (ipl-SE_AST_OFFSET-1)` (17–20) — Ceres/Pallas/Juno/Vesta given by MPC number canonicalize to the same internal slot and same rewritten public `ipl` as the named constants. This is the "Ceres-as-10001" case. |
| any other `ipl > SE_AST_OFFSET` (numbered asteroid) | `SEI_ANYBODY` (11) | = `ipl` (i.e. `SE_AST_OFFSET + MPC#`, unmodified) | `pdp` aliases the shared `SEI_ANYBODY` slot — see c-ref-se1-file.md § `sweph()` Step 1 for the file-level aliasing this implies. |

`pdp->xreturn` is the eventual output buffer regardless of branch.

**Note on the Pluto-as-asteroid remap**: `ipl == SE_AST_OFFSET + 134340` is intercepted
*earlier*, in `swe_calc()` itself (sweph.c:365, `if (ipl == SE_AST_OFFSET + 134340) ipl = SE_PLUTO;`,
run before `swecalc()` is even entered), so it never reaches this branch — it takes the
ordinary main-planet path (`SE_PLUTO`, full JPL/SWIEPH/MOSEPH cascade via `main_planet()`),
not the asteroid path described here. `swe_get_planet_name()` applies the identical remap
independently at sweph.c:6958.

### 1.2 File selection: `ifno` (sweph.c:1044–1050)

```c
if (ipli_ast > SE_AST_OFFSET) {
  ifno = SEI_FILE_ANY_AST;
} else if (ipli_ast > SE_PLMOON_OFFSET) {
  ifno = SEI_FILE_ANY_AST;    // unreachable given ipli_ast's possible values here; legacy dup
} else {
  ifno = SEI_FILE_MAIN_AST;
}
```

Net effect: `ipli_ast` is either `SEI_CHIRON..SEI_VESTA` (12–17, always `≤ SE_PLMOON_OFFSET`)
→ `SEI_FILE_MAIN_AST` (`seas*.se1`), or `SE_AST_OFFSET + MPC#` → `SEI_FILE_ANY_AST`
(`astN/seNNNNN.se1`). The middle branch is dead code for this call site (kept from when
planetary moons shared this code path more directly); see `app_pos_etc_plan()` (§1.6 below)
where the equivalent three-way split (`IS_PLANET` / `IS_MAIN_ASTEROID` / `IS_ANY_BODY`) is
live and non-dead.

### 1.3 Chiron/Pholus time-range guards (sweph.c:1051–1063)

```c
if (ipli == SEI_CHIRON && (tjd < CHIRON_START || tjd > CHIRON_END)) {
  if (serr != NULL)
    sprintf(serr, "Chiron's ephemeris is restricted to JD %8.1f - JD %8.1f",
            CHIRON_START, CHIRON_END);
  return ERR;
}
if (ipli == SEI_PHOLUS && (tjd < PHOLUS_START || tjd > PHOLUS_END)) {
  if (serr != NULL)
    sprintf(serr, "Pholus's ephemeris is restricted to JD %8.1f - JD %8.1f",
            PHOLUS_START, PHOLUS_END);
  return ERR;
}
```

Hard `ERR` (not `NOT_AVAILABLE`, no fallback) outside range. Exact constants — see §5.
No equivalent guard exists for Ceres/Pallas/Juno/Vesta or for numbered asteroids; their
range is enforced only by the `.se1` file's own `tfstart`/`tfend` inside `sweph()`
(→ `NOT_AVAILABLE`, see c-ref-se1-file.md § `sweph()` Step 3/§ File Open).

### 1.4 Earth/Sun prerequisite + asteroid evaluation (sweph.c:1064–1101)

```c
do_asteroid:
  /* earth and sun are also needed */
  retc = main_planet(tjd, SEI_EARTH, 0, epheflag, iflag, serr);
  if (retc == ERR)
    goto return_error;
  /* iflag (ephemeris bit) has possibly changed in main_planet() */
  iflag = swed.pldat[SEI_EARTH].xflgs;
  /* asteroid */
  if (serr != NULL) { strcpy(serr2, serr); *serr = '\0'; }
  retc = sweph(tjd, ipli_ast, ifno, iflag, psdp->x, DO_SAVE, NULL, serr);
  if (retc == ERR || retc == NOT_AVAILABLE)
    goto return_error;
  retc = app_pos_etc_plan(ipli_ast, 0, iflag, serr);
  if (retc == ERR)
    goto return_error;
  /* app_pos_etc_plan() might have failed, if t(light-time)
   * is beyond ephemeris range. in this case redo with Moshier */
  if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS) {
    if (epheflag != SEFLG_MOSEPH) {
      iflag = (iflag & ~SEFLG_EPHMASK) | SEFLG_MOSEPH;
      epheflag = SEFLG_MOSEPH;
      if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
        strcat(serr, "\nusing Moshier eph.; ");
      goto do_asteroid;
    } else
      goto return_error;
  }
```

`psdp = &swed.pldat[SEI_SUNBARY]` (barycentric Sun save slot, set up as a local at the top
of `swecalc()`, sweph.c:595).

**This is the single most important structural fact about asteroid calc: the asteroid
position ITSELF is always fetched via `sweph(..., ifno=SEI_FILE_MAIN_AST|SEI_FILE_ANY_AST, ...)`
— i.e. always read from the `.se1` file — regardless of `epheflag` (`SEFLG_MOSEPH`/`SEFLG_SWIEPH`/
`SEFLG_JPLEPH`). There is no Moshier (or JPL) asteroid model; `SEFLG_MOSEPH` for an asteroid body
does not mean "compute the asteroid with Moshier," it means "compute Earth (and hence the
geocentric conversion) with Moshier, but still read the asteroid itself from the `.se1` file."**
This directly contradicts the naive assumption that `SEFLG_MOSEPH` avoids all file I/O for
asteroids — it does not, and if the `.se1` file is unavailable in MOSEPH mode, asteroid calc
still fails (see below).

`epheflag` only governs how `main_planet(SEI_EARTH, ...)` computes Earth (full JPL→SWIEPH→
MOSEPH cascade documented in c-ref-calc.md §4), and — critically — whether `sweph()`'s
heliocentric→barycentric conversion step actually fires for the asteroid (see §1.5).

**Fallback cascade — what actually retries, and what doesn't:**
- If `sweph()` itself returns `NOT_AVAILABLE` or `ERR` (asteroid `.se1` file missing, or `tjd`
  outside the file's `tfstart..tfend`) → straight to `return_error`. **No MOSEPH retry** — the
  `do_asteroid:` retry loop is *not* entered in this case, because the `goto do_asteroid` is
  only reached from the `app_pos_etc_plan()` failure branch below, not from the `sweph()` call.
  A missing/out-of-range asteroid file is therefore a hard failure even if the ephemeris
  directory otherwise supports Moshier fallback for planets.
- If `app_pos_etc_plan()` returns `NOT_AVAILABLE`/`BEYOND_EPH_LIMITS` (light-time iteration
  needed a *second* `sweph()`/`sweplan()` evaluation at `t - dt` that fell outside file range —
  see c-ref-calc.md §7 and §1.6 below) — **and** `epheflag != SEFLG_MOSEPH` yet — then `iflag`/
  `epheflag` are forced to `SEFLG_MOSEPH` and control jumps back to `do_asteroid:`. This re-runs
  `main_planet(SEI_EARTH, ..., SEFLG_MOSEPH, ...)` (Earth via Moshier) and **re-runs `sweph()`
  for the asteroid again** (still reading the `.se1` file — the retry does not change how the
  asteroid itself is fetched, only how Earth is fetched and — per §1.5 — whether the
  barycentric-Sun correction is applied).
- If already `epheflag == SEFLG_MOSEPH` when `app_pos_etc_plan()` fails this way → hard error,
  no further fallback.

**JPLEPH mode**: `main_planet(SEI_EARTH, ..., SEFLG_JPLEPH, ...)` computes Earth (and, as a
side effect of `main_planet`'s internal `jplplan()`/`sweplan()` call for Earth with `DO_SAVE`,
also populates barycentric Sun — see §1.5) from the JPL ephemeris; the asteroid is still read
from `.se1` via the same unconditional `sweph()` call. So "JPLEPH + asteroid" in practice means
Earth/Sun from the JPL DE file, asteroid from Swiss Ephemeris `.se1` file — a genuine mixed-source
result, not an error and not a fallback; this is normal, intended behavior (the Swiss Ephemeris
`.se1` asteroid files are themselves JPL-derived orbit fits, so mixing is numerically consistent
to the file's fit tolerance).

### 1.5 The Sun-barycentric-vector interaction: `SEFLG_MOSEPH` + asteroid is subtly wrong in C

`psdp->x` (barycentric Sun, `swed.pldat[SEI_SUNBARY].x`) is passed as `xsunb` into `sweph()`,
which (per c-ref-se1-file.md § `sweph()` Step 8, sweph.c:2332–2343) does:

```c
if (xsunb != NULL && ((iflag & SEFLG_JPLEPH) || (iflag & SEFLG_SWIEPH))) {
  if (ipl >= SEI_ANYBODY) {
    for (i = 0; i <= 2; i++) xp[i] += xsunb[i];
    if (need_speed) for (i = 3; i <= 5; i++) xp[i] += xsunb[i];
  }
}
```

Note the *exact* condition is `xsunb != NULL && ((iflag & SEFLG_JPLEPH) || (iflag & SEFLG_SWIEPH))`
— **not just `xsunb != NULL`** (c-ref-se1-file.md's summary of this step glosses over the flag
check; the flag check is load-bearing here). When `epheflag == SEFLG_MOSEPH`:

1. `psdp` (`SEI_SUNBARY`) is **never populated** in MOSEPH mode: `main_planet()`'s MOSEPH branch
   calls `swi_moshplan(tjd, SEI_EARTH, DO_SAVE, NULL, NULL, serr)` (swemplan.c:276), which
   writes only `swed.pldat[SEI_EARTH]` — it never touches `swed.pldat[SEI_SUNBARY]`. (Contrast
   with the SWIEPH path, where `sweplan()`'s `do_sunbary` flag is forced true whenever
   `do_save` is true — see c-ref-calc.md §6 — so calling `main_planet(SEI_EARTH, ...)` with
   SWIEPH *does* refresh `SEI_SUNBARY` as a side effect.) So in MOSEPH mode `psdp->x` is
   **stale**: either zero-initialized (never computed this process) or left over from
   whatever the *last* SWIEPH/JPLEPH call happened to compute.
2. Even setting staleness aside: the `iflag` passed into the asteroid's `sweph()` call carries
   the `SEFLG_MOSEPH` bit (not `SWIEPH`/`JPLEPH`), so the flag check above is false regardless
   of what `psdp->x` contains — **the heliocentric→barycentric conversion is skipped entirely**.
   The asteroid position returned from `sweph()` remains purely heliocentric (as stored in the
   `.se1` file). This is NOT a position error for plain geocentric output, though (an earlier
   draft of this section overstated it): `app_pos_etc_plan()`'s observer is `pedp->x`
   (sweph.c:2528–2543, verified directly), and under MOSEPH `swi_moshplan` fills
   `pldat[SEI_EARTH]` with *heliocentric* Earth — heliocentric asteroid minus heliocentric
   Earth is a frame-consistent geocentric vector. The genuine defects are confined to the
   `SEI_SUNBARY` consumers: (a) the `SEFLG_HELCTR` branch (sweph.c:2517–2521) subtracts a
   stale barycentric Sun when the stale `pdp->iephe` (point below) claims SWIEPH/JPLEPH —
   wrong by up to ~0.01 AU after a prior SWIEPH/JPLEPH call in the same process; (b)
   deflection/aberration geometry reads the zero-or-stale barycentric-Sun global (zero =
   Sun-at-origin, consistent with the heliocentric frame; stale = wrong-epoch Sun, mas-level
   wobble). Full write-up: `docs/swisseph-c-potential-bugs.md` §9. Net: MOSEPH+asteroid
   output is deterministic and essentially correct in a *fresh process*, and call-history-
   dependent otherwise.

This is a real quirk/bug in the C library's stateful design, not a documentation
simplification — verified directly against sweph.c:2335 and swemplan.c:276-339. It is exactly
the class of issue this project's stateless architecture sidesteps by construction (see
`CLAUDE.md` § Stateless Tolerance), but a porter aiming for bit-exact golden-test fidelity
against C's `SEFLG_MOSEPH`-asteroid output needs to know this discrepancy exists in the
reference implementation itself.

`pdp->iephe` for an asteroid (sweph.c:2345–2352, inside `sweph()`'s save step) is set to
`psdp->iephe` (**not** unconditionally `SEFLG_SWIEPH`, unlike planet/moon files where
`ifno == SEI_FILE_PLANET || SEI_FILE_MOON` forces `pdp->iephe = SEFLG_SWIEPH`) — so an
asteroid's recorded "which ephemeris computed this" flag is inherited from whatever
`SEI_SUNBARY` last recorded, which (per point 1 above) may itself be stale/wrong in MOSEPH mode.

### 1.6 Asteroid-specific behavior inside `app_pos_etc_plan()` (sweph.c:2465–2775)

Covered generically in c-ref-calc.md §7. Asteroid-specific deltas, all in `app_pos_etc_plan()`:

**File/slot classification** (sweph.c:2480–2497) — a three-way split, parallel to but not
identical in shape to §1.2's `ifno` selection:

```c
if (ipli > SE_PLMOON_OFFSET || ipli > SE_AST_OFFSET) { // 2nd condition obsolete
  ifno = SEI_FILE_ANY_AST;  ibody = IS_ANY_BODY;        pdp = &swed.pldat[SEI_ANYBODY];
} else if (ipli == SEI_CHIRON || ipli == SEI_PHOLUS || ipli == SEI_CERES
        || ipli == SEI_PALLAS || ipli == SEI_JUNO || ipli == SEI_VESTA) {
  ifno = SEI_FILE_MAIN_AST; ibody = IS_MAIN_ASTEROID;   pdp = &swed.pldat[ipli];
} else {
  ifno = SEI_FILE_PLANET;   ibody = IS_PLANET;          pdp = &swed.pldat[ipli];
}
```

Note main asteroids (`IS_MAIN_ASTEROID`) get **their own dedicated `pldat[]` slots** (12–17),
*not* the shared `SEI_ANYBODY` slot — only numbered asteroids (`IS_ANY_BODY`) alias
`SEI_ANYBODY`. This matches c-ref-se1-file.md's `sweph()` Step 1 aliasing rule (`ipli >
SE_AST_OFFSET → SEI_ANYBODY`), which only triggers for `ipli > SE_AST_OFFSET`, i.e. numbered
asteroids — `SEI_CHIRON..SEI_VESTA` (12–17) never satisfy that.

**Light-time re-evaluation** (sweph.c:2604–2690, the loop that recomputes position at
`t - dt` for accurate light-time): dispatches on `epheflag` exactly as `main_planet()`'s main
cascade does, but with an asteroid-specific inner call in every branch:

| `epheflag` | Planet (`ibody == IS_PLANET`) | Asteroid (`ibody != IS_PLANET`) |
|---|---|---|
| `SEFLG_JPLEPH` | `swi_pleph(t, pnoint2jpl[ipli], J_SBARY, xx)` | `swi_pleph(t, J_SUN, J_SBARY, xsun)` then `sweph(t, ipli, ifno, iflag, xsun, NO_SAVE, xx)` — asteroid **always** re-read from `.se1`, Sun from JPL |
| `SEFLG_SWIEPH` | `sweplan(t, ipli, ifno, ..., NO_SAVE, xx, xearth, xsun, NULL)` | `sweplan(t, SEI_EARTH, SEI_FILE_PLANET, ..., NO_SAVE, xearth, NULL, xsun, NULL)` then `sweph(t, ipli, ifno, iflag, xsun, NO_SAVE, xx)` |
| `SEFLG_MOSEPH` (only if `SEFLG_SPEED` requested) | `swi_moshplan(t, ipli, NO_SAVE, xxsv, xearth)` | `sweph(t, ipli, ifno, iflag, NULL, NO_SAVE, xxsv)` (note: **`xsunb=NULL` here** — no barycentric conversion is even attempted, consistent with §1.5) then `swi_moshplan(t, SEI_EARTH, NO_SAVE, xearth, xearth)`; only `xx[3..5]` (speed) is kept from this recomputation, position uses the earlier `dt`-shifted estimate |

So: **light-time/aberration/deflection machinery itself is identical between planets and
asteroids** (same `app_pos_etc_plan()` code, same formulas — c-ref-calc.md §7/§11/§12 apply
unchanged); the only difference is *how the position at the shifted evaluation time `t-dt` is
re-fetched* — asteroids always go back through `sweph()`+`.se1`, planets go through
`sweplan()`/`jplplan()`/`swi_moshplan()` as appropriate. Aberration and deflection are applied
identically to both.

### 1.7 `swe_get_planet_name()` for asteroids (sweph.c:7046–7078, brief)

```c
if (ipl > SE_PLMOON_OFFSET || ipl > SE_AST_OFFSET) { // 2nd condition obsolete
  if (ipl == swed.fidat[SEI_FILE_ANY_AST].ipl[0]) {
    strcpy(s, swed.fidat[SEI_FILE_ANY_AST].astnam);      // already cached from a prior file read
  } else {
    retc = sweph(J2000, ipl, SEI_FILE_ANY_AST, 0, NULL, NO_SAVE, xp, NULL);
    if (retc != ERR && retc != NOT_AVAILABLE)
      strcpy(s, swed.fidat[SEI_FILE_ANY_AST].astnam);    // opened the file purely to read its header
    else
      sprintf(s, "%d: not found (asteroid)", ipl - SE_AST_OFFSET);
  }
  /* '?' or provisional-designation names get a secondary lookup in seasnam.txt */
  ...
}
```

If the name isn't already cached (`fdp->ipl[0]` matches), it triggers a full `sweph()` call
at epoch J2000 with `NO_SAVE` — this opens/reads the `.se1` file (running `read_const()`,
which populates `fdp->astnam`) purely as a side effect of getting the name string; the
computed position (`xp`) is discarded. If the name is `'?'` or starts with a digit
(provisional designation placeholder), a secondary lookup scans `seasnam.txt` for a
user-updatable name mapping keyed by MPC catalog number. Out of scope beyond this summary —
no position-calculation logic here.

---

## 2. `swi_gen_filename()` — Filename Generation (swephlib.c:3610–3691)

```c
void swi_gen_filename(double tjd, int ipli, char *fname)
```

### 2.1 Body → prefix dispatch (swephlib.c:3618–3653)

```c
switch(ipli) {
  case SEI_MOON:
    strcpy(fname, "semo"); break;
  case SEI_EMB: case SEI_MERCURY: case SEI_VENUS: case SEI_MARS: case SEI_JUPITER:
  case SEI_SATURN: case SEI_URANUS: case SEI_NEPTUNE: case SEI_PLUTO: case SEI_SUNBARY:
    strcpy(fname, "sepl"); break;
  case SEI_CERES: case SEI_PALLAS: case SEI_JUNO: case SEI_VESTA:
  case SEI_CHIRON: case SEI_PHOLUS:
    strcpy(fname, "seas"); break;
  default:  /* asteroid or planetary moon */
    if (ipli > SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET) {
      sprintf(fname, "sat%ssepm%d.%s", DIR_GLUE, ipli, SE_FILE_SUFFIX);
    } else {
      sform = "ast%d%sse%05d.%s";
      if (ipli - SE_AST_OFFSET > 99999)
        sform = "ast%d%ss%06d.%s";
      sprintf(fname, sform, (ipli - SE_AST_OFFSET) / 1000, DIR_GLUE,
              ipli - SE_AST_OFFSET, SE_FILE_SUFFIX);
    }
    return;  /* asteroids or planetary moons: only one file 3000 bc - 3000 ad */
}
```

Planetary-moon naming (`sat/sepm{N}.se1`) — out of scope, one line only, not covered further.

### 2.2 Main asteroids (`seas`) and planets (`sepl`) — century arithmetic (swephlib.c:3654–3691)

Falls through from the `sepl`/`seas`/`semo` cases only (the `default:` branch `return`s early).

```c
if (tjd >= 2305447.5) { gregflag = TRUE;  swe_revjul(tjd, gregflag, &jyear, &jmon, &jday, &jut); }
else                  { gregflag = FALSE; swe_revjul(tjd, gregflag, &jyear, &jmon, &jday, &jut); }

if (jyear < 0) sgn = -1; else sgn = 1;
icty = jyear / 100;
if (sgn < 0 && jyear % 100 != 0)
  icty -= 1;
while (icty % ncties != 0)   // ncties = (int) NCTIES = 6
  icty--;

if (icty < 0) strcat(fname, "m"); else strcat(fname, "_");
icty = abs(icty);
sprintf(fname + strlen(fname), "%02d.%s", icty, SE_FILE_SUFFIX);
```

`NCTIES = 6.0` (sweph.h:249, "number of centuries per eph. file") — files span **6 centuries
= 600 years**, and `icty` (the century index, `jyear/100`) is rounded *down* to the nearest
multiple of 6 (via the `while (icty % ncties != 0) icty--` loop) to find the file's start
century. Examples: JD for year 1850 → `icty = 18` → already a multiple of 6 → file `sepl_18.se1`
/ `seas_18.se1` (covers 1800–2400). Year 1700 → `icty = 17` → decremented to 12 →
`sepl_12.se1` (covers 1200–1800). Negative years get `icty -= 1` first when not an exact
century boundary (so BC year -50, i.e. `jyear=-50`, `icty = -50/100 = 0` in C integer
division, then `icty -= 1` → `-1`, then rounded down to a multiple of 6 → `-6`) — prefixed
with `"m"` (minus/BC) instead of `"_"` (AD), then `%02d` of `abs(icty)`, e.g. `seplm06.se1`.

This is identical logic and identical prefix table for both `sepl` (main planets, out of
scope, covered structurally by c-ref-calc.md) and `seas` (Chiron/Pholus/Ceres/Pallas/Juno/
Vesta) — main asteroid files use exactly the same 600-year century-bucketing scheme as the
main planet files, just with `seas` instead of `sepl` and no per-planet variation (all six
main-asteroid bodies share one `seas*.se1` file per century-bucket, distinguished internally
by `fdp->ipl[]` entries — see §6 below and c-ref-se1-file.md's per-planet loop).

### 2.3 Numbered asteroids: subdirectory + filename formula

```
"ast%d%sse%05d.%s"   with args: (ipli - SE_AST_OFFSET) / 1000, DIR_GLUE, ipli - SE_AST_OFFSET, "se1"
```

e.g. MPC# 433 (Eros) → `ast0/se00433.se1` (subdir `ast{n/1000}`, filename `se{n:05d}.se1`).
For MPC numbers > 99999, the format switches to 6-digit zero-padding with an `s` already
baked into the format string (`"ast%d%ss%06d.%s"`) — this is a *different* mechanism from the
short-file-variant `s` insertion described next; it's unconditional for very high MPC numbers.
`DIR_GLUE` is `"/"` on POSIX, `"\\"` on Windows (sweodef.h:304/319).

### 2.4 File-open retry order (sweph.c:2179–2231, referenced generically in
c-ref-se1-file.md § File Open — this is the exact sequence for **numbered asteroids**
specifically, verified against the source)

The retry loop lives in `sweph()`, not in `swi_gen_filename()` itself — `swi_gen_filename()`
only produces the *first* candidate name; `sweph()` mutates the string in place and retries:

```c
again:
  fdp->fptr = swi_fopen(ifno, s, swed.ephepath, serr);
  if (fdp->fptr == NULL) {
    if (ipli > SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET) {
      /* planetary moon: strip "sat/" once, retry */
    } else if (ipli > SE_AST_OFFSET) {
      char *spp = strchr(s, '.');
      if (spp > s && *(spp-1) != 's') {          /* no 's' before '.' yet? */
        sprintf(spp, "s.%s", SE_FILE_SUFFIX);     /* insert 's': seNNNNN.se1 -> seNNNNNs.se1 */
        goto again;
      }
      spp--;                                       /* spp now points at the 's' */
      swi_strcpy(spp, spp + 1);                     /* remove it: back to seNNNNN.se1 */
      if (subdirlen > 0 && strncmp(s, subdirnam, (size_t) subdirlen) == 0) {
        swi_strcpy(s, s + subdirlen + 1);           /* strip "astN/" prefix, retry */
        goto again;
      }
    }
    return(NOT_AVAILABLE);
  }
```

Tracing through this state machine for a numbered asteroid gives exactly **four** attempts,
in this order, before giving up:

1. `astN/seNNNNN.se1` (subdir + long name — `swi_gen_filename()`'s original output)
2. `astN/seNNNNNs.se1` (subdir + short/`s`-suffixed name — first retry: insert `s`)
3. `seNNNNN.se1` (no subdir, long name — second retry: the `s` insertion is detected on
   re-entry to the `else if` block since `*(spp-1) == 's'` now, so instead the `s` is
   *removed* and the subdir prefix is stripped)
4. `seNNNNNs.se1` (no subdir, short name — third retry: back at the top of the `else if`
   block with no subdir and no `s`, so `s` is inserted again)

After attempt 4 fails, `*(spp-1) == 's'` again, so the code strips the `s` (back to
`seNNNNN.se1`, no subdir) and checks the subdir-prefix condition — which is now false (`s` no
longer starts with `subdirnam`) — so **no fifth attempt** is made; `NOT_AVAILABLE` is returned.

**Main asteroids (Chiron/Pholus/Ceres/Pallas/Juno/Vesta) have NO retry logic at all.** Their
`ipli` (`SEI_CHIRON..SEI_VESTA`, 12–17) satisfies neither `ipli > SE_PLMOON_OFFSET && ipli <
SE_AST_OFFSET` nor `ipli > SE_AST_OFFSET`, so a failed `swi_fopen()` for `seas_NN.se1` returns
`NOT_AVAILABLE` immediately, on the first attempt.

---

## 3. MPC Elements Line Parsing — Exact Column Arithmetic (`read_const()`, sweph.c:4594–4622,
4712–4753)

c-ref-se1-file.md § Step 4 gives the loose summary (`atof(s + 35 + i)`); this section gives
the precise byte offsets and the diameter-estimation formula.

### 3.1 Locals (sweph.c:4514, 4518)

```c
char sastnam[41];
int lastnam = 19;
```

### 3.2 Step 4 — H, G, diameter (sweph.c:4594–4622, only if `ifno == SEI_FILE_ANY_AST`)

```c
if (ifno == SEI_FILE_ANY_AST) {
  sp = fgets(s, AS_MAXCH * 2, fp);          /* fourth header line: MPC-format elements record */
  if (sp == NULL || strstr(sp, "\r\n") == NULL) { smsg = "d"; goto file_damage; }

  /* find where the name starts: skip leading spaces, skip the MPC-number digits, skip one space */
  while (*sp == ' ') sp++;
  while (isdigit((int) *sp)) sp++;
  sp++;
  i = (int) (sp - s);                       /* i = byte offset (from start of line) where name begins */

  strncpy(sastnam, s, lastnam + i);          /* copy [leading spaces][MPC#][space][name, up to 19 chars] */
  *(sastnam + lastnam + i) = '\0';

  strcpy(swed.astelem, s);                   /* full raw line, kept for swe_plan_pheno() */

  swed.ast_H = atof(s + 35 + i);             /* absolute magnitude */
  swed.ast_G = atof(s + 42 + i);             /* slope parameter */
  if (swed.ast_G == 0) swed.ast_G = 0.15;    /* default slope parameter */

  strncpy(s2, s + 51 + i, 7);                /* diameter field: 7 chars, column 51+i */
  *(s2 + 7) = '\0';
  swed.ast_diam = atof(s2);
  if (swed.ast_diam == 0) {
    /* estimate the diameter from magnitude; assume albedo = 0.15 */
    swed.ast_diam = 1329 / sqrt(0.15) * pow(10, -0.2 * swed.ast_H);
  }
}
```

`i` is **not** a fixed constant — it is computed per-line as the number of leading bytes
(spaces + MPC-number digits + one separating space) before the asteroid name text starts.
The column offsets `35`, `42`, `51` for H/G/diameter are then measured **relative to `i`**,
i.e. relative to where the name field begins, not relative to the start of the line — because
the fixed-width MPC/astorb.dat-style columns for H/G/diameter are defined relative to the name
field, and the name field's own start position varies with the MPC-number's digit count.

**Diameter-from-magnitude formula** (used whenever the file's diameter field is `0`/blank):

```
D (km) = 1329 / sqrt(0.15) * 10^(-0.2 * H)
```

This is the standard albedo-based asteroid diameter estimate (`D = 1329/√albedo · 10^(-0.2H)`
km) with albedo fixed at 0.15 — the same numeric constant as the default slope parameter `G`,
but conceptually unrelated (physical albedo assumption vs. photometric phase-slope parameter);
they happen to share the value 0.15 by convention, not by code sharing.

**Worked example, verified against a real file** (`ephe/ast0/se00433s.se1`, line 4 of the
text header):

```
000433 Eros               L.H. Wasserman  10.38  0.15                    4   0 ...
```

Skip leading spaces (none), skip digits (`000433`, 6 chars), skip one space → `i = 7`
(name "Eros" starts at byte 7). Then `ast_H = atof(s + 42)` = `10.38`, `ast_G = atof(s + 49)`
= `0.15`, diameter field `s[58..65]` is blank → `atof` = 0 → H-based estimate fires:
`1329/sqrt(0.15) * 10^(-0.2*10.38)` ≈ 28.8 km. The same file's binary header has
`nplan = 1`, `ipl[0] = 10433` (= `SE_AST_OFFSET + 433`), and (short variant) time range
`tfstart/tfend = 2268922.5/2488522.5` (≈1500–2100 AD) — short (`s`-suffix) files cover
roughly 1500–2100 AD, so golden-test epochs for them must stay inside that window.

### 3.3 Step 11 — asteroid name extraction and trimming (sweph.c:4712–4753, only if
`ifno == SEI_FILE_ANY_AST`)

```c
char sastno[12];
int j = 4;                                   /* old astorb.dat: 4-digit MPC# field */
while (sastnam[j] != ' ' && j < 10)           /* new astorb.dat: 5-digit MPC# field */
  j++;
strncpy(sastno, sastnam, j);
sastno[j] = '\0';
i = (int) atol(sastno);                       /* re-parse the MPC number from sastnam */

if (i == fdp->ipl[0] - SE_AST_OFFSET || i == fdp->ipl[0] /* planetary moon */) {
  /* current-format elements record: name comes from the record we already parsed */
  strncpy(fdp->astnam, sastnam + j + 1, lastnam);
  fdp->astnam[lastnam] = '\0';
  fread((void *) s, 30, 1, fp);               /* consume+discard the old 30-byte name field */
} else {
  /* older file format: name comes from the dedicated 30-byte field instead */
  fread((void *) fdp->astnam, 30, 1, fp);
}

/* right-trim trailing spaces, then double-space-terminate */
i = (int) (strlen(fdp->astnam) - 1);
if (i < 0) i = 0;
sp = fdp->astnam + i;
while (*sp == ' ') sp--;
sp[1] = '\0';
if ((sp = strstr(fdp->astnam, "  ")) != NULL)
  *sp = '\0';
```

This cross-checks the MPC number re-parsed from the elements-line name field against
`fdp->ipl[0] - SE_AST_OFFSET` (or `fdp->ipl[0]` directly, for the planetary-moon case sharing
this code path) to decide which of two name sources to trust — a Bowell-database-era
elements record (name embedded in the line already parsed in §3.2) vs. an older format where
the name lives in a fixed 30-byte field read separately at this point in the header. Either
way, 30 bytes are consumed from the file stream here (either discarded, or read directly into
`fdp->astnam`), keeping the file cursor aligned for the CRC/physical-constants reads that
follow (c-ref-se1-file.md § Step 12).

---

## 4. Downstream Consumption of `ast_H` / `ast_G` / `ast_diam`

### 4.1 Magnitude (swecl.c, Bowell HG system)

Fully covered by `docs/c-ref-phenomena.md` §5k (`swe_pheno`, swecl.c:4048–4064) — cite that
section rather than duplicating. Summary of the split it documents: main asteroids
(`SE_CHIRON..SE_VESTA`, `ipl < NMAG_ELEM`) use the hardcoded `mag_elem[ipl]` table row (H, G
built into the table, swecl.c:3773–3801); numbered asteroids beyond Vesta
(`ipl >= NMAG_ELEM`) use `swed.ast_H`/`swed.ast_G` from the last-read `.se1` file — **except**
asteroid 1566 Icarus, which gets a hardcoded override (`me[0]=16.9, me[1]=0.15`, "elements
from JPL database", swecl.c:4053–4055) instead of reading `swed.ast_H/G`. **Main asteroids
never populate or consult `swed.ast_H`/`ast_G` at all** — their magnitude data lives entirely
in the static `mag_elem` table, sourced independently of anything in the `seas*.se1` file
(which, per §3.2, only has an MPC-elements header line — and hence populates
`ast_H`/`ast_G`/`ast_diam` — when `ifno == SEI_FILE_ANY_AST`; `SEI_FILE_MAIN_AST` files carry
no such header line at all).

### 4.2 Body radius / disc diameter (`ast_diam`)

Same pattern repeated verbatim at six call sites in `swecl.c`, all gated on
`ipl > SE_AST_OFFSET` (i.e. numbered asteroids only — never main asteroids, never fictitious
bodies):

```c
if (ipl < NDIAM)
  dd = pla_diam[ipl];                    /* or drad = pla_diam[ipl]/2/AUNIT, depending on call site */
else if (ipl > SE_AST_OFFSET)
  dd = swed.ast_diam * 1000;             /* km -> m */
else
  dd = 0;
```

Call sites (function, line, variable name/scaling used at that site):

| Function | Line | Variable | Expression |
|---|---|---|---|
| `eclipse_where()` | swecl.c:702 | `drad` | `swed.ast_diam / 2 * 1000 / AUNIT` (km → m → AU, already halved to radius) |
| `eclipse_how()` | swecl.c:1020 | `drad` | same |
| `swe_lun_occult_when_glob()` | swecl.c:1697 | `drad` | same |
| `occult_when_loc()` | swecl.c:2505 | `drad` | same |
| `swe_pheno()` | swecl.c:3892 | `dd` | `swed.ast_diam * 1000` (km → m, **not** halved — `attr[3]` computation halves it later via `dd/2/AUNIT`) — see c-ref-phenomena.md §4 |
| `swe_rise_trans_true_hor()` | swecl.c:4489 | `dd` | `swed.ast_diam * 1000` (km → m, not halved, same downstream pattern as `swe_pheno`) |

**Fallback when `ast_diam` is 0** (uninitialized global, or a numbered-asteroid `.se1` file
whose elements line genuinely had a blank/zero diameter field that wasn't itself
back-filled by the H-based estimate — cannot happen after §3.2's fallback runs, but *can*
happen if `ast_diam` was never populated at all, i.e. no `SEI_FILE_ANY_AST` file has been
read yet in this process): **no error** — the body is simply treated as a point source
(`drad`/`dd` = 0 → `asin(0/lbr[2])` → apparent disc diameter 0, eclipse/occultation geometry
degenerates to point-body case). This is a silent, non-erroring fallback, not a warning.

### 4.3 When `ast_H`/`ast_G`/`ast_diam` are populated — and the staleness hazard

These three globals (plus `swed.astelem`, the full raw elements line) are populated **only**
inside `read_const()`'s `ifno == SEI_FILE_ANY_AST` branch (sweph.c:4594, guarded exactly as
shown in §3.2) — i.e. **only** when a numbered-asteroid `.se1` file is opened. They are:

- **Never** touched when reading a `SEI_FILE_MAIN_AST` (`seas*.se1`) file — Chiron/Pholus/
  Ceres/Pallas/Juno/Vesta computations never write these globals, and (per §4.1) never read
  them either.
- **Never** touched when reading `SEI_FILE_PLANET`/`SEI_FILE_MOON`/`SEI_FILE_PLMOON` files.
- Initialized to `0.0`/`""` at process start (sweph.c:114–117, the static initializer for
  `swed`).
- **Left in place** (never cleared) after the file that populated them is closed or a
  different numbered asteroid's file is opened for a different body — i.e. these are
  **process-global, not per-file or per-body** state. Calling `swe_calc()` for asteroid A,
  then for main-planet Mars, then calling `swe_pheno()` for asteroid A again, still works
  correctly only because nothing overwrote `ast_H`/`ast_G`/`ast_diam` in between — but calling
  `swe_calc()` for a *different* numbered asteroid B and then reading `swed.ast_diam`
  expecting asteroid A's value would silently return B's value instead. C code that interleaves
  `swe_calc()` calls for two different numbered asteroids before reading these globals gets
  the wrong metadata for one of them with no error indication.

**Porting implication** (already flagged generically in `docs/c-ref-phenomena.md` § Porting
Notes point 1, `eclipse.rs:105` and `riseset.rs:463` already stub `ast_diam` to `0.0`): a
stateless Rust port must key H/G/diameter metadata to the specific asteroid/file being
evaluated (e.g. return it alongside the position in whatever result type the asteroid
calc path produces), not to a shared mutable global — there is no way to reproduce C's
"whichever numbered-asteroid file was opened most recently" semantics safely, and no reason
to try; the correct behavior is simply "the H/G/diameter that came from *this* body's file."

---

## 5. Constants Table

| Constant | Value | Source |
|---|---|---|
| `SE_AST_OFFSET` | `10000` | swephexp.h:128 (see c-ref-se1-file.md) |
| `SE_PLMOON_OFFSET` | `9000` | swephexp.h |
| `MPC_CERES` | `1` | sweph.h:72 |
| `MPC_PALLAS` | `2` | sweph.h:73 |
| `MPC_JUNO` | `3` | sweph.h:74 |
| `MPC_VESTA` | `4` | sweph.h:75 |
| `MPC_CHIRON` | `2060` | sweph.h:76 |
| `MPC_PHOLUS` | `5145` | sweph.h:77 |
| `CHIRON_START` | `1967601.5` (1 Jan 675) | sweph.h:207 (superseded old limit `1958470.5` / 1 Jan 650, commented out at sweph.h:206) |
| `CHIRON_END` | `3419437.5` (1 Jan 4650) | sweph.h:208 |
| `PHOLUS_START` | `640648.5` (1 Jan −2958 Julian) | sweph.h:216 (superseded old limit `314845.5` / 1 Jan −3850, commented out at sweph.h:215) |
| `PHOLUS_END` | `4390617.5` (1 Jan 7309) | sweph.h:217 |
| `NCTIES` | `6.0` | sweph.h:249 — centuries per `sepl`/`seas` file (600-year span) |
| `SE_FILE_SUFFIX` | `"se1"` | sweph.h:192 |
| `SEI_ANYBODY` | `11` | sweph.h (see c-ref-se1-file.md) |
| `SEI_CHIRON..SEI_VESTA` | `12..17` | sweph.h (see c-ref-se1-file.md) |
| `NMAG_ELEM` | `SE_VESTA + 1 = 21` | swecl.c:3759 (see c-ref-phenomena.md) |
| Diameter-from-H albedo assumption | `0.15` | sweph.c:4619 (same literal as the default `ast_G`, coincidental) |
| SE_AST_OFFSET+1566 (Icarus) magnitude override | `H=16.9, G=0.15` | swecl.c:4053–4055 |

---

## 6. Porting Notes

- **Asteroids never use Moshier or JPL position models.** `SEFLG_MOSEPH`/`SEFLG_JPLEPH`/
  `SEFLG_SWIEPH` for an asteroid body only select how *Earth* (and, indirectly, the
  barycentric Sun) is computed; the asteroid's own position always comes from `sweph()` +
  the `.se1` file (`SEI_FILE_MAIN_AST` or `SEI_FILE_ANY_AST`). A missing/out-of-range asteroid
  file is a hard `NOT_AVAILABLE`/`ERR` with **no fallback to Moshier**, even if
  `SEFLG_MOSEPH` was requested and Earth successfully falls back. Do not model this as "three
  interchangeable asteroid backends" the way main planets have three backends — there is
  exactly one asteroid position source.

- **`SEFLG_MOSEPH` + asteroid has a genuine C-side correctness quirk** (§1.5): the
  heliocentric→barycentric conversion for the asteroid position is skipped whenever the
  computation is in MOSEPH mode (both because `SEI_SUNBARY` is never freshly computed by
  `swi_moshplan()`, and independently because `sweph()`'s conversion guard explicitly checks
  for the `SWIEPH`/`JPLEPH` flag bits). If golden-test parity against C's `SEFLG_MOSEPH`
  asteroid output is required, this needs to be replicated bit-for-bit (or, more likely,
  called out as a known-divergence in the golden harness the way `stateless_tolerance`
  documents the deflection-speed and SPEED3 divergences) — a "correct" stateless
  implementation that always converts to barycentric would deliberately *not* match C here.
  Get an explicit decision before implementing rather than assuming either "replicate the bug"
  or "fix it" is obviously right.

- **`seas*.se1` (main-asteroid) files store six bodies per file**, distinguished by
  `fdp->ipl[]` entries. **Hex-dump-verified against a real `ephe/seas_18.se1`** (binary
  header at offset 116, little-endian, nplan=6, tfstart/tfend = 2378496.5/2597996.5):
  `ipl[] = [12, 13, 14, 15, 16, 17]` — the **internal `SEI_*` numbers**
  (`SEI_CHIRON=12, SEI_PHOLUS=13, SEI_CERES=14, SEI_PALLAS=15, SEI_JUNO=16, SEI_VESTA=17`),
  in that order, NOT the public `SE_*` numbers (15–20). `read_const()`'s per-planet loop
  routes each directly to its own slot via the `else` branch (`pdp = &swed.pldat[ipli]`
  with `ipli` = 12..17), consistent with §1.6's observation that `app_pos_etc_plan()`
  never routes main asteroids through `SEI_ANYBODY`. A Rust port's file-level body-id
  lookup for main asteroids must therefore use 12..17 (matching `sweph()`'s `ipli_ast`,
  which for main asteroids IS the `SEI_*` value — §1.1), the same convention as the
  existing `body_file_id()` returning SEI-style ids for planet files.

- **Numbered-asteroid `ipl[0]`** is `SE_AST_OFFSET + MPC#` (per §1.1/§3.3's cross-check
  against `fdp->ipl[0] - SE_AST_OFFSET`), and each numbered-asteroid `.se1` file contains
  exactly one body (`fdp->npl == 1`), unlike `sepl`/`seas` files which pack several bodies
  per century-file.

- **Output frame**: asteroid `.se1` files use `SEI_FLG_ROTATE` (equinoctal-element Chebyshev
  packing, per c-ref-se1-file.md § `rot_back()`), just like the main planets — `rot_back()`'s
  non-Moon branch rotates into the ecliptic J2000 frame (no further Moon-specific equatorial
  rotation applies, since `ipli == SEI_MOON` is false for every asteroid). So `sweph()`
  returns asteroids in **heliocentric ecliptic-rectangular J2000** coordinates before the
  barycentric-conversion step (§1.5) — same frame convention as heliocentric main planets
  (`SEI_FLG_HELIO` semantics, c-ref-se1-file.md § Plan Data Flags), consistent with the
  "asteroids are heliocentric" comment directly in the C source (sweph.c:2333).
  `app_pos_etc_plan()` then applies the same barycentric/geocentric/topocentric conversions,
  light-time, aberration, deflection, frame bias, and precession as it does for planets
  (c-ref-calc.md §7) — no asteroid-specific branching exists past the position-fetch stage
  except the light-time re-evaluation dispatch table in §1.6.

- **`ifno` selection has two near-duplicate implementations** that must be kept in sync
  conceptually: `swecalc()`'s dispatch (§1.2, `ipli_ast`-based, with one dead branch) and
  `app_pos_etc_plan()`'s dispatch (§1.6, `ipli`-based, three live branches with `ibody`
  classification). Both reduce to the same two-way split (main asteroid → `SEI_FILE_MAIN_AST`,
  numbered asteroid → `SEI_FILE_ANY_AST`) but arrive at it via different variable names and
  different (partially dead) conditionals. A Rust port should implement this classification
  **once** (e.g. an `AsteroidKind { Main(SeiId), Numbered(u32) }` enum or similar) and reuse it
  at both call sites, rather than porting each C dispatch site's dead-branch structure
  literally.

- **Pluto-as-asteroid (`SE_AST_OFFSET + 134340`) remap happens twice, independently, in two
  different public entry points** (`swe_calc()` at sweph.c:365, `swe_get_planet_name()` at
  sweph.c:6958) — both must be ported if either is, and a shared helper is warranted per this
  project's `constraints` guidance against duplicating logic.
