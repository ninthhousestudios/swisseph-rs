# C Reference: Planetary Moons & Center-of-Body (SEFLG_CENTER_BODY) — sweph.c

Porting reference for planetary-moon bodies (`SE_PLMOON_OFFSET` = 9000, `ipl` 9000–9999,
`sat/sepm*.se1` files) and the `SEFLG_CENTER_BODY` flag, both introduced in SE 2.10 (checkout
verified at `SE_VERSION "2.10.03"`, sweph.h:65).

**Prerequisite reading** (not duplicated here):
- `docs/c-ref-asteroid.md` — the numbered-asteroid `.se1` path. Planetary moons **reuse** most
  of this machinery (`SEI_FILE_ANY_AST`, `SEI_ANYBODY` slot aliasing, `sweph()`'s generic file-open/
  segment/Chebyshev-decode pipeline). This doc calls out every point of reuse and every point of
  divergence explicitly — see § Porting Notes for the summary table.
- `docs/c-ref-se1-file.md` — `.se1` binary format, `read_const()`, `get_new_segment()`,
  `rot_back()`, `do_fread()`, `swi_echeb`/`swi_edcheb`.
- `docs/c-ref-calc.md` — `swe_calc()`, `swecalc()` main-planet branches, `main_planet()`
  ephemeris cascade, `sweplan()`, `app_pos_etc_plan()`'s generic light-time/aberration/
  deflection/precession pipeline (shared unchanged by plmoon bodies — see §6).

Bodies covered (`../swisseph/ephe/sat/plmolist.txt`): 9401/9402 Mars moons (Phobos/Deimos),
9501–9504 Jupiter moons (Io/Europa/Ganymede/Callisto), 9599 Jupiter COB, 9601–9608 Saturn moons,
9699 Saturn COB, 9701–9705 Uranus moons, 9799 Uranus COB, 9801/9802/9808 Neptune moons (Triton/
Nereid/Proteus), 9899 Neptune COB, 9901–9905 Pluto moons (Charon/Nix/Hydra/Kerberos/Styx), 9999
Pluto COB. Numbering scheme: `9pmm` where `p` = parent-planet `SE_*` number (`SE_MARS=4` …
`SE_PLUTO=9`) and `mm` = moon index 01.. within that planet, or `99` reserved for "center of
body" (COB).

---

## 1. Constants

| Constant | Value | Source |
|---|---|---|
| `SE_PLMOON_OFFSET` | `9000` | swephexp.h:127 |
| `SE_AST_OFFSET` | `10000` | swephexp.h:128 (plmoon range is `[9000, 10000)`, strictly below asteroids) |
| `SEFLG_CENTER_BODY` | `1024*1024` (bit 20) | swephexp.h:216, "calculate position of center of body (COB)" |
| `SEFLG_TEST_PLMOON` | `2*1024*1024 \| SEFLG_J2000 \| SEFLG_ICRS \| SEFLG_HELCTR \| SEFLG_TRUEPOS` | swephexp.h:218, "test raw data in files sepm9*" — debug-only flag bundle, see §7.4 |
| `SEI_FILE_ANY_AST` | `3` | sweph.h:176 — the file-cache slot plmoon reads actually use |
| `SEI_FILE_PLMOON` | `5` | sweph.h:178 — **defined but never referenced anywhere else in the codebase** (verified: no other hit in any `.c`/`.h`). Dead/vestigial constant; do not treat its existence as evidence that plmoon files get a dedicated file-cache slot — they don't (§3). |
| `SEI_ANYBODY` | `11` | sweph.h:146 — the shared `pldat[]` slot for both numbered asteroids and plmoon bodies |
| `SEI_MARS`..`SEI_PLUTO` | `4`..`9` | sweph.h:139–144, numerically identical to `SE_MARS`..`SE_PLUTO` (swephexp.h:105–110) — the two constant families interchange freely in the plmoon guards below |
| `SE_NPLANETS` | `23` | swephexp.h:125 — used by `swe_nod_aps`'s reject range (§8) |

---

## 2. `swe_calc()` Front-Door Translation (sweph.c:416–437)

This is the **only** place the public `ipl` (9000–9999, or a main-planet `ipl` combined with
`SEFLG_CENTER_BODY`) gets normalized into an internal `(ipl=parent planet, iplmoon=9pmm)` pair.
Nothing downstream re-derives this from scratch — `iplmoon` is threaded as an explicit extra
parameter through `swecalc()` → `main_planet()` → `app_pos_etc_plan()`.

```c
// sweph.c:416-437
/* planetary center of body or planetary moon: either planet is called
 * with SEFLG_CENTER_BODY or center of body with ipl = 9n99 is called.
 * we want to handle both cases the same way. */
// planet is called with SE_PLUTO etc. and SEFLG_CENTER_BODY:
// get number of center of body
if ((iflag & SEFLG_CENTER_BODY) && ipl <= SE_PLUTO && (iflag & SEFLG_TEST_PLMOON) != SEFLG_TEST_PLMOON) {
  iplmoon = ipl * 100 + 9099; // planetary center of body
}
// planet center of body or planetary moon is called using 9... number:
// moon number and planet number
if (ipl >= SE_PLMOON_OFFSET && ipl < SE_AST_OFFSET && (iflag & SEFLG_TEST_PLMOON) != SEFLG_TEST_PLMOON) {
  iplmoon = ipl; // planetary center of body or planetary moon
  ipl = (int) ((ipl - 9000) / 100);
  iflag |= SEFLG_CENTER_BODY;
}
// with Mercury to Mars, we do not have center of body different from barycenter
if ((iflag & SEFLG_CENTER_BODY) && ipl <= SE_MARS && (iplmoon % 100) == 99) {
  iplmoon = 0;
  iflag &= ~SEFLG_CENTER_BODY;
}
if ((iflag & SEFLG_CENTER_BODY) || iplmoon > 0)
  swi_force_app_pos_etc();
```

Two entry paths converge on the same `(ipl, iplmoon)` representation:

1. **Planet ipl + `SEFLG_CENTER_BODY`** (e.g. `swe_calc(t, SE_JUPITER, SEFLG_CENTER_BODY, ...)`):
   `iplmoon = ipl*100 + 9099` (line 422) — always synthesizes the **COB** moon-number (`...99`)
   for that planet, never a specific named moon. There is no way to request "Io's offset,
   applied to a plain-Jupiter call with CENTER_BODY" — CENTER_BODY on a bare planet ipl always
   means COB specifically.
2. **Direct 9pmm `ipl`** (e.g. `swe_calc(t, 9503, 0, ...)` for Ganymede, or `swe_calc(t, 9599, 0,
   ...)` for Jupiter/COB): `iplmoon = ipl` unchanged (retains the full moon/COB number), and
   `ipl` is **overwritten in place** to the parent planet number via integer division:
   `ipl = (ipl - 9000) / 100`. `SEFLG_CENTER_BODY` is set **unconditionally** for *any* value in
   `[9000, 10000)` — real moons (`...01`..`...08`) as well as COB (`...99`). The flag name is
   generic "center body machinery is active", not "this is specifically a COB query" — see
   §5's `calc_center_body()` which just does `xx += xcom` regardless of which moon `xcom` holds.
3. **Guard for Mercury–Mars COB** (line 432–435): if the *derived* parent is Mercury/Venus/Earth/
   Mars (`ipl <= SE_MARS`, note `SE_MARS=4` so this is Sun(0)/Moon(1)/Mercury(2)/Venus(3)/Mars(4))
   **and** the moon-suffix is specifically `99` (COB), the whole thing is cancelled:
   `iplmoon=0`, `SEFLG_CENTER_BODY` cleared. Rationale (comment): "we do not have center of body
   different from barycenter" for these bodies — no COB offset file exists for Sun/Moon/Mercury/
   Venus, and although real Mars moons (9401/9402 Phobos/Deimos) **do** exist and are *not*
   caught by this guard (suffix `01`/`02` ≠ `99`), a bare `9499` (hypothetical "Mars COB") would be
   silently discarded to a plain Mars calc if it existed. In practice no `94xx` COB file is
   shipped, so this guard is defensive rather than load-bearing against real data.
4. **`SEFLG_TEST_PLMOON` bypass**: if this exact bundle (`SEFLG_TEST_PLMOON` — not just any
   subset of its bits, note the `!=` compares the *masked* `iflag` against the *full* constant
   value) is present, **both** translation branches above are skipped entirely — `ipl` stays as
   the raw 9pmm value and `iplmoon` stays 0. This routes the 9pmm value into the **ordinary
   minor-planet branch** of `swecalc()` (§3.1's `ipl > SE_PLMOON_OFFSET` condition), i.e. treats
   the sepm file as a self-contained heliocentric body exactly like a numbered asteroid, for
   testing/inspecting the raw file contents without the planet-relative-offset semantics. This
   is a debug/test-only path (`swetest -tpm`, swetest.c:930–931) — not part of the normal public
   API contract, but a Rust port should preserve it if bit-exact swetest parity is a golden-test
   goal, since it changes which code path (§3.1 vs §3.2/§4) actually executes for the same `ipl`.

**Save-area interaction**: `sd->tsave == tjd && ... && iplmoon == 0` (sweph.c:456) — the
fast-path cache hit requires `iplmoon == 0`. Any plmoon/COB call (`iplmoon != 0`) **always**
forces a fresh `swecalc()` invocation, never reusing the per-`ipl` save slot. This matters for
a stateless Rust port only as a behavioral note (the Rust port has no save-area at all), but it
means the C reference output for repeated identical plmoon queries in the same process is not
"free" the way a repeated plain-planet query is — every plmoon call re-walks the full pipeline.

---

## 3. Dispatch Through `main_planet()` / `app_pos_etc_plan()` / `calc_center_body()`

### 3.1 `main_planet()` (sweph.c:1562–1573): fetch the moon/COB offset once, before the planet cascade

```c
static int main_planet(double tjd, int ipli, int iplmoon, int32 epheflag, int32 iflag, char *serr)
{
  int retc;
  if ((iflag & SEFLG_CENTER_BODY)
    && ipli >= SE_MARS && ipli <= SE_PLUTO) {
    /* jupiter center of body, relative to jupiter barycenter */
    retc = sweph(tjd, iplmoon, SEI_FILE_ANY_AST, iflag, NULL, DO_SAVE, NULL, serr);
    if (retc == ERR || retc == NOT_AVAILABLE)
      return ERR;
  }
  switch(epheflag) { /* ... ordinary JPL/SWIEPH/MOSEPH planet cascade, unmodified ... */ }
}
```

`ipli` here is the **parent planet's** `SEI_*` number (Mars=4 .. Pluto=9 — note the guard
excludes Sun/Moon/Mercury/Venus, consistent with §2's Mercury–Mars-COB carve-out, but *does*
include Mars itself, so real Mars moons 9401/9402 correctly reach this branch even though
Mars-COB does not exist).

This `sweph()` call reads the sepm file **directly at full epoch `tjd`**, `do_save=DO_SAVE`
(→ result cached in `swed.pldat[SEI_ANYBODY].x`, since `ipli=iplmoon > SE_PLMOON_OFFSET` aliases
to `SEI_ANYBODY` inside `sweph()` itself — see §4), and crucially **`xsunb=NULL`** — no
heliocentric→barycentric conversion is attempted for the moon/COB offset (see §5 for why this is
correct, not an oversight). Only after this does the ordinary JPL/SWIEPH/MOSEPH cascade for the
**parent planet's own barycentric position** run, completely unchanged from a plain planet call.

**Failure propagation**: if the moon/COB file is missing or `tjd` is outside its `tfstart..tfend`,
`main_planet()` returns `ERR` immediately — **before even attempting the parent planet's own
position**. A missing/out-of-range plmoon file is a hard failure with no fallback to "just
return the parent planet instead", exactly parallel to the asteroid doc's "no Moshier fallback for
missing asteroid files" (c-ref-asteroid.md §1.4).

### 3.2 `calc_center_body()` (sweph.c:2445–2455): the actual `xx += xcom` addition

```c
static int calc_center_body(int32 ipli, int32 iflag, double *xx, double *xcom, char *serr)
{
  if (!(iflag & SEFLG_CENTER_BODY))
    return OK;
  if (ipli < SEI_MARS || ipli > SEI_PLUTO)
    return OK;
  for (i = 0; i <= 5; i++)
    xx[i] += xcom[i];
  return OK;
}
```

Trivial vector addition (position **and** speed, `i<=5`), gated only on the flag and the same
Mars–Pluto parent range. It is called from `app_pos_etc_plan()` at **two points** with two
different sources for `xcom`:

- **Before** the light-time loop (sweph.c:2512): `calc_center_body(ipli, iflag, xx,
  swed.pldat[SEI_ANYBODY].x, serr)` — uses the offset **already cached** by `main_planet()`'s
  `sweph()` call at full-precision `tjd` (§3.1). This gives the "true position" (pre-light-time)
  moon/COB barycentric vector, `xx0` in the surrounding code.
- **After** the light-time re-evaluation at `t - dt` (sweph.c:2605–2691): a **fresh** `sweph(t,
  iplmoon, SEI_FILE_ANY_AST, iflag, NULL, NO_SAVE, xcom, serr)` call re-reads the offset at the
  light-time-shifted epoch, then `calc_center_body(ipli, iflag, xx, xcom, serr)` adds *that*
  offset to the freshly-fetched planet position at `t-dt`. This is the same "re-fetch at t-dt"
  pattern the asteroid doc documents for asteroid positions (c-ref-asteroid.md §1.6) — plmoon
  bodies get an equivalent second file read for light-time accuracy, structurally parallel to
  (but a separate, simpler call than) the asteroid epheflag-dispatch table, since the moon-offset
  re-fetch is unconditionally via `sweph()` regardless of `epheflag` (whereas the *planet*
  position re-fetch below it in the same function still dispatches JPL/SWIEPH/MOSEPH per usual).

**Ordering relative to HELCTR**: both `calc_center_body()` calls happen **before** the
`SEFLG_HELCTR` barycentric-Sun subtraction that follows them in the same function (sweph.c:2516–
2520 and 2692–2696). So the moon offset is added to the **barycentric** planet position, and
*then* (if HELCTR requested) barycentric Sun is subtracted from the combined (planet+moon)
vector — algebraically equivalent to "heliocentric planet + moon offset" since the offset itself
is translation-invariant between barycentric and heliocentric frames. This is why `xcom` need
not itself be barycentric-corrected (§5).

### 3.3 `app_pos_etc_plan()` file/slot classification is unaffected by plmoon

The three-way `ifno`/`ibody` classification documented in c-ref-asteroid.md §1.6
(`ipli > SE_PLMOON_OFFSET || ipli > SE_AST_OFFSET` → `SEI_FILE_ANY_AST`/`IS_ANY_BODY`) governs
the classification of `ipli`, which for plmoon calls is always the **parent planet's** `SEI_*`
number (4–9), never the moon number — so plmoon calls always fall into the **plain
`IS_PLANET`** branch of that classification (`ifno = SEI_FILE_PLANET`), exactly like an ordinary
planet call. The moon/COB-specific file I/O is entirely confined to the two `calc_center_body`-
adjacent `sweph(..., iplmoon, SEI_FILE_ANY_AST, ...)` calls described in §3.1/§3.2 — it does not
touch the main `ifno`/`ibody` dispatch that governs the *planet's own* position fetch at all.

---

## 4. `sweph()` — Slot Aliasing (sweph.c:2137–2141)

```c
static int sweph(double tjd, int ipli, int ifno, int32 iflag, double *xsunb, AS_BOOL do_save, double *xpret, char *serr)
{
  ...
  ipl = ipli;
  if (ipli > SE_AST_OFFSET)
    ipl = SEI_ANYBODY;
  if (ipli > SE_PLMOON_OFFSET)
    ipl = SEI_ANYBODY;
  pdp = &swed.pldat[ipl];
  ...
```

Both numbered-asteroid `ipli` (`> SE_AST_OFFSET`) **and** plmoon `ipli` (`> SE_PLMOON_OFFSET`,
which includes everything in `[9001, 9999]` since that's `> 9000`) alias to the same
`swed.pldat[SEI_ANYBODY]` (slot 11) storage — **identical aliasing to numbered asteroids**, not
a separate slot. Consequence, also identical to the asteroid case: only **one** plmoon body (or
numbered asteroid — they share the slot with each other too) can be "current" at a time; calling
for a different 9pmm value (or a numbered asteroid) invalidates the cache
(`ipl == SEI_ANYBODY && ipli != pdp->ibdy` at sweph.c:2168 forces the file to close and reopen).

**File-cache slot** (`ifno`, the *file*-level index into `swed.fidat[]`, distinct from the
*body*-level `pldat[]` slot above) is **also** shared: every plmoon call uses `ifno =
SEI_FILE_ANY_AST` (index 3) — the identical file-cache slot numbered asteroids use, **not**
the unused `SEI_FILE_PLMOON` (index 5) constant. So a plmoon read and a numbered-asteroid read
in the same process contend for the same open-file-handle slot too, on top of sharing the
`pldat[]` position-cache slot. A Rust port's file-cache design should treat "plmoon" and
"numbered asteroid" as the same cache-eviction class if it chooses to model C's single-slot
caching behavior at all (the stateless Rust architecture likely sidesteps this entirely per
`CLAUDE.md`'s Stateless Design principle, but it's worth knowing the C behavior if a golden test
interleaves plmoon and asteroid calls and expects to see file-reopen side effects/timing).

---

## 5. File Content Semantics — Planetocentric, Not Heliocentric

This is the most important divergence from the asteroid `.se1` convention and is **verified
directly against real shipped files**, not inferred from comments.

### 5.1 Hex/struct-level verification

Parsed `ephe/sat/sepm9599.se1` (Jupiter/COB), `sepm9401.se1` (Phobos/Mars), `sepm9699.se1`
(Saturn/COB), and `sepm9501.se1` (Io/Jupiter) directly (Python struct-unpack replicating
`read_const()`'s exact byte layout — file length matched exactly for all four, confirming the
parse is byte-accurate):

| File | `ipl[0]` | `tfstart`/`tfend` | `dseg` (days/segment) | `ncoe` | `rmax_raw` | `iflg` |
|---|---|---|---|---|---|---|
| sepm9599.se1 (Jupiter/COB) | 9599 | 2378491.5 / 2524599.5 (~1800–2200) | 4.0 | 39 | 10 | `0x08` |
| sepm9401.se1 (Phobos/Mars) | 9401 | 2415015.5 / 2469082.5 (~1900–2048) | 1.0 | 39 | 10000 | `0x08` |
| sepm9699.se1 (Saturn/COB) | 9699 | 2378491.5 / 2524599.5 | 4.0 | 40 | 10 | `0x08` |
| sepm9501.se1 (Io/Jupiter) | 9501 | 2378491.5 / 2524599.5 | 4.0 | 40 | 10 | `0x08` |

`iflg = 0x08` for all four → decoded against `SEI_FLG_*` (sweph.h:165–168, `HELIO=1, ROTATE=2,
ELLIPSE=4, EMBHEL=8`): **`HELIO` bit clear, `ROTATE` bit clear, `ELLIPSE` bit clear, `EMBHEL` bit
set**.

### 5.2 What each bit implies

- **`SEI_FLG_ROTATE` clear**: `sweph()`'s segment-load step (sweph.c:2277–2284) only calls
  `rot_back()` (the equinoctal-elements-to-rectangular rotation asteroid/planet files rely on)
  `if (pdp->iflg & SEI_FLG_ROTATE)`. Since this bit is clear for plmoon files, `rot_back()` is
  **never invoked** — `pdp->neval = pdp->ncoe` directly (sweph.c:2282) and the raw Chebyshev
  coefficients decoded by `swi_echeb`/`swi_edcheb` (sweph.c:2295–2302) are used **as-is** as
  rectangular X/Y/Z (and their time-derivatives for speed). **This directly contradicts the
  asteroid pattern** — c-ref-asteroid.md's Porting Notes state "asteroid `.se1` files use
  `SEI_FLG_ROTATE` ... just like the main planets" (equinoctal-element packing). Plmoon files do
  **not**: they store direct rectangular coefficients, no orbital-plane rotation step at all.
- **`SEI_FLG_HELIO` clear**: per the field's own comment (c-ref-se1-file.md:132, "true if helio,
  false if bary"), this nominally marks the data as *not* heliocentric-in-the-file-storage sense.
  However, this bit is **never actually read** in the code path plmoon calls take — the
  helio→barycentric conversion in `sweph()` (sweph.c:2335, `if (xsunb != NULL && ((iflag &
  SEFLG_JPLEPH) || (iflag & SEFLG_SWIEPH)))`) is gated on the **caller-supplied `xsunb` pointer**,
  not on `pdp->iflg & SEI_FLG_HELIO`. Both plmoon-reading call sites (main_planet.c:1570 and
  app_pos_etc_plan sweph.c:2609) pass **`xsunb = NULL` explicitly** — so this conversion branch
  is unconditionally skipped for plmoon reads regardless of what the `HELIO` bit says. The
  `HELIO` bit in plmoon files' `iflg` byte is therefore **inert/vestigial** for the actual
  computation (it may just reflect how Astrodienst's file-generation tooling reused the general
  per-planet header writer without a plmoon-specific mode).
- **`SEI_FLG_EMBHEL` set**: per sweph.c:2312 (`if (ipl == SEI_SUNBARY && (pdp->iflg &
  SEI_FLG_EMBHEL))`), this bit is only ever consulted when `ipl == SEI_SUNBARY` (10). Since
  plmoon reads always alias `ipl` to `SEI_ANYBODY` (11, §4), this branch can never fire for
  plmoon data regardless of the bit's value — **also inert** for plmoon computation. Likely
  vestigial for the same reason as `HELIO` above.

**Conclusion (the "breaks the naive assumption" finding)**: plmoon `.se1` Chebyshev coefficients
decode directly to **rectangular offset vectors, in whatever frame the parent-planet's own
barycentric-equatorial-J2000 position is stored in** (since `calc_center_body()` simply adds
them, §3.2) — i.e. **planetocentric rectangular offsets** (moon/COB relative to parent planet,
or COB relative to planet barycenter), not heliocentric orbital-element-packed asteroid-style
data. No sun-vector addition, no orbital-plane rotation. This is a fundamentally simpler
decode than the asteroid path, not a directory-relabeled copy of it.

### 5.3 `rmax` scaling — an exact, narrower condition than "planetary moons get finer scale"

```c
// sweph.c:4830-4835
pdp->rmax = lng / 1000.0;
// planet's center of body, e.g. 9599 for Jupiter or Mars moons
if (ipli >= SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET) {
  if ((ipli % 100) == 99 || (ipli - 9000) / 100 == SE_MARS)
    pdp->rmax = lng / 1000000.0;
}
```

The finer (`/1000000.0`, i.e. 1000× smaller unit) scale applies **only** to (a) any COB entry
(`...99` suffix, any parent planet) **or** (b) **Mars's own moons specifically** (Phobos/Deimos) —
verified numerically: `sepm9599.se1`/`sepm9699.se1` (both COB) and `sepm9401.se1` (Mars moon,
non-COB) all hit the `/1e6` branch (`rmax` = 1e-5, 1e-5, 0.01 respectively), while `sepm9501.se1`
(Io, a Jupiter moon, non-COB) uses the ordinary `/1000.0` branch (`rmax` = 0.01, same numeric
value here but via the coarser divisor on a 10× larger `rmax_raw`... no — raw was `10` for both
sepm9599 and sepm9501; the divisors differ, `10/1e6=1e-5` vs `10/1e3=0.01`). Moons of
Jupiter/Saturn/Uranus/Neptune/Pluto (all *non*-COB, i.e. real named moons other than Mars's) use
the standard `/1000.0` divisor identical to asteroid/planet `rmax` derivation
(c-ref-se1-file.md:186–188 already documents this rule, but summarizes it loosely as "planetary
moon center-of-body" — the precise condition additionally includes Mars's real moons, not just
COB entries, which this section's C-verified reading makes explicit).

### 5.4 Time segmentation is per-file, not fixed

`dseg` (days per Chebyshev segment) differs by body: 4.0 days for Jupiter/Saturn moons and all
COB files sampled, but **1.0 day** for Phobos (Mars's inner moon orbits in ~7.66 hours — needs
much finer sampling than the 4-day segments adequate for slower-moving outer-planet moons).
Likewise `tfstart`/`tfend` differ per file — see §9.

### 5.5 Elements-line parsing collateral damage (H/G/diam garbage) — genuinely new finding

`read_const()`'s MPC-orbital-elements-line parser (c-ref-asteroid.md §3.2, `swed.ast_H`/
`ast_G`/`ast_diam`) is gated **only** on `ifno == SEI_FILE_ANY_AST` (sweph.c:4594) — with
**no additional check** that the body being read is actually a numbered asteroid rather than a
plmoon body. Since plmoon reads always use `ifno = SEI_FILE_ANY_AST` (§4), this parser **runs
unconditionally on every plmoon file open too** — but the plmoon file's 4th header line is *not*
an MPC elements record; it's free-form text (`"009599 Jupiter/COB  Ephemeris / WWW_USER Fri Dec
18 13:37:26 2020 Pasadena, USA ... "`). Reading the exact fixed byte offsets (`s+35+i`, `s+42+i`,
`s+51+i..+58+i` where `i=7` for a 6-digit-number + 1-space prefix) against this differently-
shaped line yields **garbage, not zero**: verified directly against the `sepm9599.se1` line —
`ast_H = atof("ri Dec ")` = `0.0` (stops at non-numeric `'r'`), `ast_G = atof("18 13:37:")` =
`18.0` (!, a nonsensical slope-parameter value, and critically **not 0** so the `if (ast_G==0)
ast_G=0.15` default-fallback at sweph.c:4597 does **not** fire), `ast_diam = atof("26 2020")` =
`26.0` (also nonzero, so the H-based diameter-estimate fallback at sweph.c:4599–4602 does not
fire either). These globals genuinely get overwritten with meaningless values parsed out of a
date/time string every time a plmoon file is opened.

**Why this doesn't visibly break anything**: every downstream consumer of `ast_diam`
(c-ref-asteroid.md §4.2's six `swecl.c` call sites) gates its use on `ipl > SE_AST_OFFSET`
(10000) — which **excludes** every plmoon `ipl` (9000–9999) by construction. Likewise the
`ast_H`/`ast_G` magnitude consumer (c-ref-asteroid.md §4.1) is gated on `ipl >= NMAG_ELEM` (21)
which plmoon `ipl` values also satisfy... **but** magnitude computation for plmoon bodies is
moot in practice since no call site in the codebase invokes `swe_pheno`/eclipse magnitude logic
keyed on a *plmoon* `ipl` in a way that would read these particular globals (magnitude/diameter
consumers key on the `ipl` passed to `swe_pheno` etc. directly, and `swe_pheno` itself only
special-cases `ipl > SE_AST_OFFSET`/`ipl < NMAG_ELEM`, treating any plmoon `ipl` as neither — see
§8). So this parsing bug is **currently inert for the plmoon body's own output**, but it *does*
pollute the process-global `ast_H`/`ast_G`/`ast_diam`/`astelem` state for whatever numbered
asteroid is queried *next* in the same process — exactly the staleness hazard
c-ref-asteroid.md §4.3 already documents, just with a previously-undocumented garbage-data
source (plmoon file opens) feeding into it. A stateless Rust port that (per that section's
existing porting implication) keys H/G/diameter metadata to the specific body being evaluated
rather than a shared global sidesteps this bug by construction — but should not attempt to
port the parsing logic to plmoon files at all, since there is no meaningful data to parse there.

**Interesting asymmetry**: the very next parsing step in the same function — asteroid-name
extraction (c-ref-asteroid.md §3.3, sweph.c:4712–4753) — **is** explicitly plmoon-aware: its
cross-check is `if (i == fdp->ipl[0] - SE_AST_OFFSET || i == fdp->ipl[0] /* planetary moon */)`
(sweph.c:4723–4724, comment present in the source itself). For a plmoon file, the MPC-number-
shaped prefix parsed from the (differently-formatted, but still numeric-prefixed) line is
`i = fdp->ipl[0]` directly (e.g. `9599`), which matches the second arm of this OR — so plmoon
**names** ("Jupiter/COB", "Phobos/Mars") are correctly extracted into `fdp->astnam`, while the
H/G/diameter fields three statements earlier in the same `if (ifno == SEI_FILE_ANY_AST)` block
are not given equivalent treatment. One code path was updated for plmoon awareness when this
feature was added in 2.10; the other, functionally adjacent, was not.

---

## 6. Center-of-Body (`9n99`) Bodies

COB entries exist for Jupiter (9599), Saturn (9699), Uranus (9799), Neptune (9899), and Pluto
(9999) — **not** for Sun/Moon/Mercury/Venus/Mars/Earth (§2's Mercury–Mars guard; also no `94xx`,
`90xx`..`93xx` COB files are shipped in `ephe/sat/`). A COB position represents the offset between
the *photocenter/body-center* and the *system barycenter* used for the planet's own orbital
integration (relevant because outer planets' "position" in the main `sepl*.se1`/JPL ephemerides
is nominally the system barycenter, not the visible disc center — for Jupiter with its four large
Galilean moons this offset is non-negligible at high precision).

Structurally, COB is **not a special case** anywhere in the calc pipeline past §2's `iplmoon =
ipl*100 + 9099` synthesis — it flows through `main_planet()`/`calc_center_body()`/`sweph()`
exactly like a real named moon; the `...99` suffix only matters for: (a) the Mercury–Mars
cancellation guard (§2), (b) the `rmax` fine-scale condition (§5.3, though real Mars moons also
qualify), and (c) error-message text generation (sweph.c:2245–2249 distinguishes `"plan. COB No.
%d"` vs `"plan. moon No. %d"` purely by testing whether the resolved filename contains the
substring `"99."` — a string-matching heuristic on the generated filename, not a semantic
`iplmoon % 100 == 99` check, so it would misfire if a hypothetical future moon file were ever
named ending in literal digits `99` before the extension for unrelated reasons; harmless in
practice given the fixed `sepm{iplmoon}.se1` naming scheme guarantees the substring only appears
for genuine `...99` COB filenames).

Time ranges for COB files sampled: 1800–2200 (same 400-year span as the outer-planet moons of
that same parent) — see §9.

---

## 7. `SEFLG_CENTER_BODY` — Complete Semantics

### 7.1 Accepted `ipl` values

- Any main-planet `ipl` in `{SE_JUPITER..SE_PLUTO}` (5–9) combined with the flag → synthesizes
  COB (§2 case 1). `{SE_SUN..SE_MARS}` (0–4) combined with the flag → flag is silently cleared,
  **no error**, behaves as a plain call (§2 case 3).
- Any direct plmoon `ipl` (9000–9999) → flag is **force-set** internally regardless of whether
  the caller passed it explicitly (§2 case 2) — passing `SEFLG_CENTER_BODY` explicitly alongside
  a `9xyz` ipl is redundant but harmless (the bit is already going to be set).
- Any other `ipl` (Sun/Moon-family constants below `SE_PLMOON_OFFSET`, mean node, fictitious
  bodies, numbered asteroids) → the flag is simply **never inspected** by anything in `swecalc()`
  outside the main-planet branch (`calc_center_body()`'s own guard, `ipli < SEI_MARS || ipli >
  SEI_PLUTO`, independently protects against misuse even if the flag somehow reached that
  function for a non-planet body) — effectively a silent no-op, not an error.

### 7.2 Interaction with other flags

- **`SEFLG_HELCTR`/`SEFLG_BARYCTR`**: orthogonal — CENTER_BODY's `xx += xcom` addition (§3.2)
  happens *before* the HELCTR/BARYCTR frame decision in `app_pos_etc_plan()`, so COB/moon
  positions are available in heliocentric, barycentric, or (the default) geocentric/topocentric
  output exactly as any planet would be.
- **`SEFLG_TOPOCTR`**: unaffected — topocentric observer-position logic in `app_pos_etc_plan()`
  operates on `pedp->x`/`xobs` independent of the CENTER_BODY branch.
- **`SEFLG_J2000`**: orthogonal, ordinary precession-frame flag, applied to the final combined
  (planet+offset) vector like any other body.
- **`SEFLG_SPEED`/`SEFLG_SPEED3`**: speed is included in the `xx[i] += xcom[i]` addition
  (`i<=5`, §3.2) — the moon/COB's own orbital speed relative to its parent is added to the
  parent's speed, giving the combined body's absolute speed. No special-casing needed; ordinary
  vector addition of two velocity vectors is exactly correct here (unlike position, which needs
  no correction either since it's a simple translation).
- **`SEFLG_TEST_PLMOON`**: mutually exclusive by construction with normal CENTER_BODY translation
  — see §2 case 4; when this bundle is present, plain CENTER_BODY-style handling never engages
  at all (the 9pmm value is instead routed as if it were a self-contained heliocentric body via
  the ordinary minor-planet branch, c-ref-asteroid.md §1.1's condition `ipl > SE_PLMOON_OFFSET`).

### 7.3 Behavior for unsupported bodies

No error is ever raised for "CENTER_BODY requested but not supported" — it degrades silently to
a plain calculation (§2 case 3, Mercury–Mars) or is simply ignored (§7.1, non-planet bodies).
There is no `flags_used`-style C mechanism to signal this to the caller beyond the flag's
presence/absence in the returned `iflag` (§7.4) — a Rust port choosing to surface this via
`CalcResult.flags_used` (per this project's architecture) would need to explicitly detect the
Mercury–Mars-COB cancellation case itself, since C gives no other signal.

### 7.4 Survives into `retflag`? **Yes.**

`swe_calc()`'s final return value is built (sweph.c:544–548) as `sd->iflgsave & ~SEFLG_COORDSYS`
combined with the *original* caller's coordinate-system bits (`SEFLG_EQUATORIAL|SEFLG_XYZ|
SEFLG_RADIANS`). `sd->iflgsave` is exactly the `iflag` value `swecalc()` returned, which (absent
any code that strips `SEFLG_CENTER_BODY` — none was found in `swecalc()`, `main_planet()`,
or `app_pos_etc_plan()` past the initial §2 translation) still carries the bit that was forced
on in `swe_calc()` itself. So a caller reading the returned `iflag` from `swe_calc()` **can**
detect "CENTER_BODY was applied" by testing the bit in the return value, including for the
"9pmm was passed in directly" case where the caller never set the bit explicitly — an
implicit-but-real bit of API feedback.

---

## 8. Other Entry Points — Accept/Reject Summary (brief, no deep dive)

| Function | Accepts plmoon `ipl` (9000–9999)? | Accepts `SEFLG_CENTER_BODY` on a planet ipl? | Notes |
|---|---|---|---|
| `swe_pheno`/`swe_pheno_ut` (swecl.c:3802) | Yes, transparently — delegates position to `swe_calc(ipl, iflag\|SEFLG_XYZ, ...)` internally (swecl.c:3839), which applies the full §2 translation. | Yes, same mechanism. | Magnitude/diameter tables (`mag_elem[]`, `ast_diam`) have **no plmoon-specific entries**: `ipl` (9xxx) fails both `ipl < NMAG_ELEM` (21) and `ipl > SE_AST_OFFSET` (10000) gates in the diameter/magnitude code (swecl.c:3891, 3902, 4048), so `dd`/magnitude computation for a plmoon body degenerates to the point-source/no-table-entry fallback — position is correct, photometric attributes are not meaningful. |
| `swe_nod_aps`/`swe_nod_aps_ut` (swecl.c:5075) | **No — hard `ERR`.** Explicit reject: `(ipl >= SE_NPLANETS && ipl <= SE_AST_OFFSET)` (swecl.c:5138–5141, `SE_NPLANETS=23`) covers the entire plmoon range and returns `ERR` with `serr = "nodes/apsides for planet %5.0f are not implemented"`. | Not applicable — flag is never inspected in this function (no `PLMOON`/`CENTER_BODY` reference anywhere in `swecl.c` outside the functions in this table); a CENTER_BODY bit set on an ordinary planet `ipl` passed to `swe_nod_aps` is silently carried but has no effect on the node/apsides computation itself. | Clean, unambiguous, well-defined rejection — no silent wrong-answer risk here. |
| `swe_get_orbital_elements` (swecl.c:5783) | **Yes, silently — and produces physically dubious output.** Only rejects `ipl<=0` or mean-node/apogee/interpolated-apsides constants (swecl.c:5803); a plmoon `ipl` passes through to internal `swe_calc(tjd, ipl, ...)` calls (which apply §2's translation, so the fetched position is the moon/COB's own heliocentric-ish position) and to `get_gmsm()` (swecl.c:5687), whose mass lookup only special-cases `SE_MOON` and `{SE_MERCURY..SE_PLUTO}`/`SE_EARTH` (swecl.c:5694–5717); any other `ipl` (including plmoon) falls into the generic "asteroid or fictitious object" branch (swecl.c:5722+) which assumes **zero perturbing mass** and a **pure heliocentric two-body** model. Net effect: orbital elements are computed for the moon/COB's heliocentric position as if it directly orbited the Sun with `GM = HELGRAVCONST` — essentially the *parent planet's* heliocentric orbital elements (the moon-relative offset is astronomically negligible at this level), not "moon around planet" elements. Not an error, but not the elements a caller might expect from the name of the function applied to e.g. Io. | N/A (this function does not take an ephemeris flag path that inspects CENTER_BODY at all). | Worth flagging for a Rust port only as "if implemented, match this quirk or explicitly reject" — no evidence C intends this as a real feature for plmoon bodies. |
| `swe_fixstar*`, eclipse/occultation (`swe_lun_occult_when_glob`, `eclipse_where`, `eclipse_how`, `swe_rise_trans_true_hor`) | Not applicable to fixstars. Eclipse/occultation functions accept any `ipl` transparently via internal `swe_calc()` calls (so plmoon positions **are** computed correctly), but their disc-diameter logic uses the same `ipl > SE_AST_OFFSET` gate documented in c-ref-asteroid.md §4.2 — plmoon `ipl` always fails it, so these bodies are always treated as **point sources** (`drad`/`dd` = 0) in eclipse/occultation geometry, silently, same fallback behavior asteroids get when `ast_diam` is unset. | Not specially handled; ordinary flag pass-through. | No crash, no error — just always-point-source geometry for any plmoon body passed to these functions. |

---

## 9. Time Ranges

Sampled directly from file headers (§5.1 table); general pattern:

- **Jupiter/Saturn/Uranus/Neptune moons and all COB files**: `tfstart=2378491.5` (~1 Jan 1800),
  `tfend=2524599.5` (~1 Jan 2200) — a uniform 400-year span across every file of this class
  sampled (9599, 9699, 9501 all matched exactly).
- **Mars moons (Phobos/Deimos)**: narrower — `tfstart=2415015.5` (~1900), `tfend=2469082.5`
  (~2048), roughly 148 years — presumably reflecting a shorter high-precision numerical
  integration arc for the fast-orbiting, closely-perturbed Martian satellite system (also
  consistent with the finer `dseg=1.0`-day segmentation, §5.4, vs 4.0 days elsewhere).
- Pluto's moons were not independently sampled in this pass but are expected to follow the
  general-outer-planet 400-year pattern (all Pluto-system files are unusually small on disk,
  e.g. `sepm9902.se1`/`sepm9905.se1`, consistent with a small system needing fewer distinct
  perturbation terms rather than a different time span — not independently verified).

**Out-of-range behavior**: no plmoon-specific hard guard exists (unlike Chiron/Pholus's explicit
`CHIRON_START`/`CHIRON_END` checks, c-ref-asteroid.md §1.3) — range enforcement is purely the
generic `sweph()` file `tfstart`/`tfend` check (sweph.c:2235, `NOT_AVAILABLE`, with an
error-message branch at sweph.c:2245–2249 that specifically distinguishes "plan. COB No. %d" vs
"plan. moon No. %d" text, §6). Since `main_planet()` fetches the moon/COB offset **before**
attempting the parent planet's own position (§3.1), an out-of-range plmoon query fails with
`NOT_AVAILABLE`/`ERR` even if the *parent planet's* position would have been computable via
JPL/Moshier for that same epoch (main planet ephemerides span millennia; plmoon files span only
~150–400 years) — i.e. requesting Jupiter/COB for 1000 AD fails outright, even though plain
`swe_calc(1000AD, SE_JUPITER, 0, ...)` would succeed via any backend.

---

## 10. Apparent-Position Pipeline — Deltas vs. Standard Planet/Asteroid Path

Per §3.3, `app_pos_etc_plan()`'s generic light-time/aberration/deflection/precession machinery
(c-ref-calc.md §7/§11/§12) is applied **completely unchanged** to the combined (planet+moon-
offset) vector — there is no plmoon-specific branching anywhere in that machinery beyond the two
`calc_center_body()` call sites (§3.2) that inject the offset before the standard pipeline runs
on the result. Concretely:

- **Light-time**: computed from the combined vector's distance to the observer, exactly as for
  a plain planet — i.e. light-time delay reflects the *parent planet's* distance (the
  planetocentric offset is many orders of magnitude too small to meaningfully change light-time
  iteration convergence), but the *position itself* used downstream is the moon/COB's true
  combined position at the light-time-corrected epoch (via the second `calc_center_body()` call
  using a freshly-fetched offset at `t-dt`, §3.2).
- **Aberration/light deflection**: applied to the final combined vector by the same code that
  handles any planet — no moon-specific geometry construction.
- **Speed**: combined via simple vector addition (§7.2) at both the "before light-time" and
  "after light-time re-fetch" stages — no special SPEED3/central-difference handling beyond
  what ordinary planets already get.

**No new fidelity concerns are introduced by this project's stateless architecture** for plmoon
bodies beyond what already applies to the parent planet (per `CLAUDE.md`'s Stateless Tolerance
section) — the moon/COB offset itself is a pure, stateless Chebyshev evaluation (§5.2, no
rotation/sun-vector steps to diverge on), so any C-vs-Rust divergence here should be no larger
than the parent planet's own known divergence sources (deflection speed, SPEED3 boundaries).

---

## 11. Porting Notes

1. **Plmoon is not "asteroid with a different directory."** The filename lives in a different
   subdirectory (`sat/` vs `astN/`) and reuses the `SEI_FILE_ANY_AST`/`SEI_ANYBODY` cache slots
   (identical to numbered asteroids, §4), but the **file content semantics are fundamentally
   different**: planetocentric rectangular offsets with no rotation/sun-addition step (§5),
   vs. asteroids' heliocentric orbital-element-packed (`SEI_FLG_ROTATE`) data requiring a
   sun-vector addition. A Rust port must implement a genuinely separate decode path for plmoon
   bodies, not parameterize the asteroid decode with a different directory string.

2. **The result is computed by addition, not by direct file read.** `swe_calc(9599, ...)` never
   directly returns what's in `sepm9599.se1` — it returns `(Jupiter's own barycentric position,
   via the ordinary JPL/SWIEPH/MOSEPH planet cascade) + (offset read from sepm9599.se1)`. Both
   halves must be computed; the plmoon file alone is insufficient. (Exception: `SEFLG_TEST_PLMOON`
   bypasses this, §2 case 4 — a debug-only path, not the normal contract.)

3. **`SEFLG_CENTER_BODY` and direct-9pmm-ipl are two names for the same underlying computation.**
   `swe_calc(t, SE_JUPITER, SEFLG_CENTER_BODY, ...)` and `swe_calc(t, 9599, 0, ...)` compute
   identical output (both resolve to `ipl=SE_JUPITER, iplmoon=9599` internally, §2). A Rust API
   should probably normalize both call shapes to one internal representation early, exactly as C
   does, rather than maintaining two independent code paths.

4. **`rmax` fine-scale condition is `COB OR Mars-moon`, not `all plmoon`** (§5.3) — get this
   exact condition right; it's easy to over-generalize from the C source comment
   ("e.g. 9599 for Jupiter or Mars moons") into "all plmoon files use the fine scale," which is
   false for Io/Europa/Ganymede/Callisto/Titan/etc.

5. **Do not port the H/G/diameter elements-line parser to plmoon files at all** (§5.5) — it's a
   genuine C-side bug (parses garbage from an unrelated text field) that happens to be inert only
   because of an unrelated downstream gate (`ipl > SE_AST_OFFSET`). A stateless Rust port with
   no shared globals to pollute has no reason to replicate this; simply never populate
   H/G/diameter metadata for plmoon bodies (there is no legitimate source for it in these files).

6. **Mercury/Venus/Sun/Moon/Earth silently ignore `SEFLG_CENTER_BODY`** (§2 case 3, §7.1) — no
   error. If a Rust port wants to surface "this flag had no effect" via `flags_used` (per this
   project's architecture pattern), it needs to detect this case explicitly; C gives no signal
   beyond the bit's absence in the returned flags.

7. **`swe_nod_aps` hard-rejects all plmoon `ipl`; `swe_get_orbital_elements` silently accepts
   but produces parent-planet-equivalent (not moon-around-planet) elements** (§8) — decide
   explicitly whether the Rust port's equivalent functions should replicate the reject (clean)
   or the silent-degenerate-accept (matches C bit-for-bit but is arguably not useful) for
   `get_orbital_elements`; this needs a product decision, not an assumption either way.

8. **Time ranges are much narrower than main-planet ephemerides and are enforced with no
   fallback** (§9) — a plmoon query can fail even when the parent planet succeeds. A Rust port's
   error type should distinguish "parent planet unavailable" from "moon/COB offset unavailable
   for otherwise-valid planet epoch," since C conflates them into one `ERR`/`NOT_AVAILABLE` but a
   Rust `Error` enum could usefully preserve which half failed.

9. **`SEI_FILE_PLMOON` (sweph.h:178) is a dead constant** — don't be misled into modeling a
   dedicated file-cache slot for plmoon files; they use `SEI_FILE_ANY_AST`, sharing cache
   contention with numbered asteroids (§4).

10. **Filename retry is single-step (strip `sat/`), not the asteroid's four-attempt cascade**
    (§12 table) — get the retry-order difference right if porting file-resolution fallback logic.

### IDENTICAL vs. DIFFERS summary (plmoon vs. numbered-asteroid pattern)

| Aspect | Numbered Asteroid | Planetary Moon / COB |
|---|---|---|
| `pldat[]` slot | `SEI_ANYBODY` (shared across all numbered asteroids) | **IDENTICAL** — same `SEI_ANYBODY` slot, shared with asteroids too |
| File-cache slot (`ifno`) | `SEI_FILE_ANY_AST` | **IDENTICAL** — same slot, `SEI_FILE_PLMOON` const unused |
| Filename directory | `astN/` (subdir by MPC#÷1000) | **DIFFERS** — `sat/` flat, no subdirectory tiering |
| Filename pattern | `se{MPC#:05d}.se1` (or `s{MPC#:06d}` for MPC#>99999) | **DIFFERS** — `sepm{9pmm}.se1`, no digit-count branching |
| File-open retry cascade | 4 attempts: subdir+long, subdir+short(`s`), no-subdir+long, no-subdir+short (c-ref-asteroid.md §2.4) | **DIFFERS** — 1 retry only: strip `sat/` prefix, retry flat filename; no `s`-suffix short-file variant exists for plmoon |
| Coefficient frame | Heliocentric ecliptic-rectangular (equinoctal-element packed) | **DIFFERS** — planetocentric rectangular offset, no orbital-element packing |
| `SEI_FLG_ROTATE` | Set — `rot_back()` runs | **DIFFERS** — clear; `rot_back()` never runs, raw Chebyshev = final rectangular coords |
| Helio→barycentric conversion (`xsunb` addition in `sweph()`) | Applied (JPL/SWIEPH modes) | **DIFFERS** — never applied; both call sites pass `xsunb=NULL` (offset added directly to parent's barycentric position instead, §3.2/§5.2) |
| `rmax` scale | Always `lng/1000.0` | **DIFFERS** — `lng/1000000.0` for COB or Mars-moon entries specifically (§5.3), else same `/1000.0` |
| 4th-header-line format | Genuine MPC orbital-elements record | **DIFFERS in content, IDENTICAL in code path** — same parser runs, but the line isn't MPC-shaped, producing garbage H/G/diam (§5.5); name-extraction cross-check *is* plmoon-aware and works correctly |
| Position computed via | Direct file read (converted to barycentric) | **DIFFERS** — file read gives an *offset*, added to a separately-computed parent-planet position (§3.1/§3.2) |
| Light-time/aberration/deflection pipeline | `app_pos_etc_plan()`, standard | **IDENTICAL** — same function, same formulas, applied to the combined vector (§10) |
| Public entry-point ipl range | `SE_AST_OFFSET + MPC#` (>10000) | **DIFFERS** — `SE_PLMOON_OFFSET..SE_AST_OFFSET-1` (9000–9999), or a main-planet ipl + `SEFLG_CENTER_BODY` |
| Time-range hard guard | None (Chiron/Pholus only, not numbered asteroids) | **IDENTICAL** — none; generic file `tfstart`/`tfend` only (§9) |
| Missing-file fallback | Hard `NOT_AVAILABLE`/`ERR`, no Moshier retry | **IDENTICAL** in kind — hard failure, no fallback (§3.1) |
