# C Reference: Small Public Utilities

Porting reference for the remaining small public-API utility functions not yet ported: equation
of time / LMT↔LAT conversion, body/house/ayanamsa name lookups, the delta-T user-override, and
the centisecond string-formatting family. Read this instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swe_time_equ` | sweph.c:7387–7413 | Yes |
| `swe_lmt_to_lat` | sweph.c:7415–7423 | Yes |
| `swe_lat_to_lmt` | sweph.c:7425–7436 | Yes |
| `swe_get_planet_name` | sweph.c:6946–7125 | Yes |
| `swe_house_name` | swehouse.c:827–859 | Yes |
| `swe_get_ayanamsa_name` | sweph.c:7127–7133 (table sweph.c:130–179) | Yes |
| `swe_set_delta_t_userdef` | swephlib.c:3176–3184 | Yes |
| `swe_cs2timestr` | swephlib.c:3864–3886 | Yes |
| `swe_cs2lonlatstr` | swephlib.c:3888–3916 | Yes |
| `swe_cs2degstr` | swephlib.c:3918–3929 | Yes |
| `swe_csroundsec` | swephlib.c:3836–3843 | Yes |

Supporting/related pieces read for context (not separately ported, or already ported elsewhere):
- `swi_get_fict_name` — swemplan.c:513–520 (fictitious-body name resolution, called from
  `swe_get_planet_name`)
- `read_elements_file` — swemplan.c:694–908 (reads `seorbel.txt`; supplies both the fictitious
  planet's orbital elements *and* its display name)
- `plan_fict_nam[SE_NFICT_ELEM]` — swemplan.c:506–511 (built-in fallback names, 15 entries)
- Asteroid-name-from-.se1-file extraction (`sweph()`'s `read_const` MPC-header parsing) —
  **already ported** in `src/sweph_file/parse.rs` (`extract_asteroid_name`, `AsteroidMeta.name`,
  documented in `docs/c-ref-asteroid.md` §3.3). This doc only covers the piece that still needs
  porting: the `seasnam.txt` fallback lookup and the public dispatcher that ties body-ID →
  name together.

## 1. `swe_time_equ` (sweph.c:7387–7413)

Computes the **equation of time**: `E = LAT − LMT` (Local Apparent Time minus Local Mean Time),
returned in **days**. Input `tjd_ut` is UT.

```c
int32 CALL_CONV swe_time_equ(double tjd_ut, double *E, char *serr)
{
  int32 retval;
  double t, dt, x[6];
  double sidt = swe_sidtime(tjd_ut);
  int32 iflag = SEFLG_EQUATORIAL;
  iflag = plaus_iflag(iflag, -1, tjd_ut, serr);
  if (swi_init_swed_if_start() == 1 && !(iflag & SEFLG_MOSEPH) && serr != NULL) {
    strcpy(serr, "Please call swe_set_ephe_path() or swe_set_jplfile() before calling "
                 "swe_time_equ(), swe_lmt_to_lat() or swe_lat_to_lmt()");
  }
  if (swed.jpl_file_is_open)
    iflag |= SEFLG_JPLEPH;
  t = tjd_ut + 0.5;
  dt = t - floor(t);
  sidt -= dt * 24;
  sidt *= 15;
  if ((retval = swe_calc_ut(tjd_ut, SE_SUN, iflag, x, serr)) == ERR) {
    *E = 0;
    return ERR;
  }
  dt = swe_degnorm(sidt - x[0] - 180);
  if (dt > 180)
    dt -= 360;
  dt *= 4;
  *E = dt / 1440.0;
  return OK;
}
```

### Algorithm, step by step
1. **Sidereal time**: `sidt = swe_sidtime(tjd_ut)` — Greenwich Apparent Sidereal Time in
   **hours** (already includes nutation/obliquity internally; see `src/sidereal_time.rs`).
2. **Flags**: start from `SEFLG_EQUATORIAL` (note: this flag is set but never actually used —
   the code proceeds to read `x[0]`, the *ecliptic* longitude, from `swe_calc_ut`'s output;
   `SEFLG_EQUATORIAL` has no effect on `x[0]`'s meaning here since it isn't OR'd into anything
   that changes the array layout in a way this function reads — it is dead/vestigial in this
   call). Then `plaus_iflag(iflag, -1, tjd_ut, serr)` sanitizes/defaults the flag (`ipl = -1`
   means "no specific body" for the plausibility check). If a JPL file happens to already be
   open globally (`swed.jpl_file_is_open`), `SEFLG_JPLEPH` is OR'd in — this is a **global-state
   read**, not something the caller controls via arguments.
3. **Local sidereal time at Greenwich, in degrees, aligned to the current UT fractional day**:
   ```c
   t = tjd_ut + 0.5;        // shift to start-of-day epoch convention
   dt = t - floor(t);       // fractional part of the day
   sidt -= dt * 24;         // subtract elapsed sidereal hours corresponding to elapsed UT hours
   sidt *= 15;              // hours -> degrees (15 deg/hour)
   ```
   This produces `sidt` as GAST **at 0h UT of the current day**, expressed in degrees — i.e. it
   removes the diurnal rotation component already accumulated since midnight, leaving the
   "sidereal time at 0h" reference angle. This is the trick that lets the subsequent longitude
   comparison directly yield the equation of time without an explicit RA-of-Sun computation.
4. **Sun's ecliptic longitude**: `swe_calc_ut(tjd_ut, SE_SUN, iflag, x, serr)`. On error, `*E = 0`
   and return `ERR` immediately (note: `*E` is explicitly zeroed on the error path, unlike many
   other functions in this file that leave output params untouched on error).
5. **Core formula**:
   ```c
   dt = swe_degnorm(sidt - x[0] - 180);
   if (dt > 180) dt -= 360;   // fold into [-180, 180]
   dt *= 4;                  // degrees -> minutes of time (1 deg = 4 min, since 360deg = 1440 min)
   *E = dt / 1440.0;          // minutes -> days (1440 min/day)
   ```
   `sidt - x[0] - 180` is (sidereal-time-at-0h minus Sun's mean-frame longitude, offset by 180°)
   — algebraically equivalent to the standard "hour angle of the mean Sun minus hour angle of
   the true Sun" formulation of the equation of time, folded through the degnorm/fold-to-signed
   idiom used throughout this codebase (compare `swe_difdeg2n`, already ported in `src/math.rs`,
   though this is not literally `difdeg2n` — it's `degnorm(...)` then a manual `> 180` fold, not
   `>= 180`; **note the boundary difference**: `swe_difdeg2n` uses `>= 180.0` while this inline
   code uses `> 180`, so at exactly `dt == 180.0` the two would disagree by choosing opposite
   signs of the same magnitude — negligible in practice for this application but worth
   preserving literally rather than reusing `difdeg2n` for bit-for-bit fidelity).
6. Return `OK` on success.

### Global-state reads
- `swed.jpl_file_is_open` — whether a JPL ephemeris file happens to already be open, silently
  changing which ephemeris backend computes the Sun's position. **A stateless port cannot
  reproduce this** — the Rust port must take the ephemeris source explicitly from
  `EphemerisConfig`/the flags the caller passes, and should NOT silently upgrade to JPL based on
  incidental prior state. This is a deliberate, documented divergence (see Porting notes).
- `swi_init_swed_if_start()` — lazy global-state initializer; only affects a warning message when
  `serr` is requested and no ephemeris path/JPL file has been configured yet. Not relevant to a
  stateless port (the equivalent condition — no ephemeris configured — should surface as a
  regular `Error`, not a warning string, per this project's error-handling convention).

### Error contract
Returns `ERR`/`OK` (int32); `*E` is zeroed on error. Port as `Result<f64, Error>`, propagating
whatever `Error` the inner `calc_ut(tjd_ut, Body::Sun, flags)` produces.

## 2. `swe_lmt_to_lat` / `swe_lat_to_lmt` (sweph.c:7415–7436)

Convert between **Local Mean Time** (civil clock time, tied to mean solar motion) and **Local
Apparent Time** (sundial time, tied to the true Sun's position). Both take `tjd_*` as a Julian
Day already shifted into the *local* time zone implied by `geolon` (i.e. these are not UTC
conversions — the caller is expected to have already applied a fixed zone offset, and this
function applies the *variable*, date-dependent equation-of-time correction on top).

`geolon`: **east positive** (standard Swiss Ephemeris convention — matches `swe_set_topo`).
`tjd_lmt0 = tjd_lmt - geolon / 360.0` converts a **local** mean time to the **Greenwich**-referenced
UT instant needed by `swe_time_equ` (dividing degrees of longitude by 360 converts to a fraction
of a day; east-positive longitude means local time runs *ahead* of Greenwich, hence subtraction
here — local civil noon at a place east of Greenwich corresponds to an *earlier* Greenwich UT).

### `swe_lmt_to_lat` — no iteration
```c
int32 CALL_CONV swe_lmt_to_lat(double tjd_lmt, double geolon, double *tjd_lat, char *serr)
{
  int32 retval;
  double E, tjd_lmt0;
  tjd_lmt0 = tjd_lmt - geolon / 360.0;
  retval = swe_time_equ(tjd_lmt0, &E, serr);
  *tjd_lat = tjd_lmt + E;
  return retval;
}
```
Single call to `swe_time_equ` at the (converted-to-Greenwich) input instant; `E` (LAT − LMT) is
added directly: `tjd_lat = tjd_lmt + E`. No iteration is needed here because the equation of time
is evaluated at the *known* LMT/UT instant — going from mean time to apparent time is a direct
lookup.

### `swe_lat_to_lmt` — 2 extra iterations (3 evaluations total)
```c
int32 CALL_CONV swe_lat_to_lmt(double tjd_lat, double geolon, double *tjd_lmt, char *serr)
{
  int32 retval;
  double E, tjd_lmt0;
  tjd_lmt0 = tjd_lat - geolon / 360.0;
  retval = swe_time_equ(tjd_lmt0, &E, serr);
  /* iteration */
  retval = swe_time_equ(tjd_lmt0 - E, &E, serr);
  retval = swe_time_equ(tjd_lmt0 - E, &E, serr);
  *tjd_lmt = tjd_lat - E;
  return retval;
}
```
Going the other direction requires iteration because the equation of time `E` itself is a
function of the (unknown) mean-time instant, not the apparent-time instant the caller supplied.
The C code performs a **fixed 3 total evaluations, unconditionally** (no convergence check, no
early exit, no tolerance test):
1. `E₀ = time_equ(tjd_lmt0)` — first guess, evaluated at the apparent-time instant itself
   (treating LAT ≈ LMT for the first pass).
2. `E₁ = time_equ(tjd_lmt0 - E₀)` — refine using the previous `E` as a correction to the
   evaluation instant.
3. `E₂ = time_equ(tjd_lmt0 - E₁)` — one more refinement pass.
4. `*tjd_lmt = tjd_lat - E₂` (note: subtracts the **last-computed** `E`, i.e. `E₂`, from the
   *original* `tjd_lat`, not from `tjd_lmt0`).

This is fixed-point iteration (not Newton's method) with exactly 2 refinement steps beyond the
initial guess — **hardcoded, not a convergence loop**. A faithful port must call the equivalent
of `time_equ` exactly 3 times in this exact sequence, not "iterate until converged" — the number
3 is load-bearing for bit-for-bit fidelity (the equation of time changes by at most ~30
seconds/day near its extrema, so 2 refinements over a same-day evaluation converge far below
double precision in practice, but the C code doesn't check — it just stops after 3).

### Error contract
`retval` is simply the last `swe_time_equ` return code — errors from earlier iterations are
**silently discarded** if a later iteration happens to return `OK` (though in practice, if the
first call fails, the ephemeris/flag problem will persist through the following calls with the
same behavior, so this is mostly academic). Port as `Result<f64, Error>` propagating the *last*
call's `Result` via `?`, matching this discard-earlier-errors behavior — or, more defensibly,
propagate the *first* error immediately via `?` after each call (short-circuiting) since in
practice a failing `time_equ` at any iteration means the ephemeris source itself is broken for
every subsequent call too. Either is defensible; prefer short-circuiting on first error since
it's simpler and the C behavior only "recovers" in the vanishingly unlikely case where
intermediate ephemeris state changes between three back-to-back calls at nearly the same epoch.

## 3. `swe_get_planet_name` (sweph.c:6946–7125)

Resolves any body ID to a display name. Full switch/dispatch:

### Hardcoded names (verbatim, from sweph.h:79–96 — confirmed by direct read)
```
SE_NAME_SUN        = "Sun"
SE_NAME_MOON       = "Moon"
SE_NAME_MERCURY    = "Mercury"
SE_NAME_VENUS      = "Venus"
SE_NAME_MARS       = "Mars"
SE_NAME_JUPITER    = "Jupiter"
SE_NAME_SATURN     = "Saturn"
SE_NAME_URANUS     = "Uranus"
SE_NAME_NEPTUNE    = "Neptune"
SE_NAME_PLUTO      = "Pluto"
SE_NAME_MEAN_NODE  = "mean Node"
SE_NAME_TRUE_NODE  = "true Node"
SE_NAME_MEAN_APOG  = "mean Apogee"
SE_NAME_OSCU_APOG  = "osc. Apogee"
SE_NAME_INTP_APOG  = "intp. Apogee"
SE_NAME_INTP_PERG  = "intp. Perigee"
SE_NAME_EARTH      = "Earth"
SE_NAME_CHIRON     = "Chiron"
SE_NAME_PHOLUS     = "Pholus"
SE_NAME_CERES      = "Ceres"
SE_NAME_PALLAS     = "Pallas"
SE_NAME_JUNO       = "Juno"
SE_NAME_VESTA      = "Vesta"
```
Note the lowercase-first-word style on the node/apogee entries ("mean Node", "true Node", "osc.
Apogee", etc.) and the abbreviations "osc."/"intp." — copy these exactly, they are easy to
mistype from memory.

### Dispatch logic, in order
1. **Pluto alias** (sweph.c:6958): `if ipl == SE_AST_OFFSET + 134340` (i.e. asteroid #134340,
   Pluto's MPC number) → treat as `SE_PLUTO`. This mirrors `normalize_asteroid_aliases`
   (`src/context.rs`, swisseph-rs/101) already ported for `calc_inner` — **reuse that same
   normalization helper here** rather than re-implementing the alias check, per this project's
   no-duplicate-logic constraint.
2. **Name cache** (sweph.c:6960–6963, a **global-state read**): `swed.i_saved_planet_name` /
   `swed.saved_planet_name` — a single-entry memoization cache (last-body-name-computed). If
   `ipl` matches the cached ID, the cached string is returned directly, **skipping all further
   logic including the asteroid-file/seasnam.txt read**. This is pure C-side performance
   optimization with **no observable effect on output** (the cache always holds the last value
   computed by this exact same function, so hitting it produces byte-identical output to
   recomputing) — **do not port the cache**; a stateless Rust function that always recomputes is
   behaviorally equivalent and simpler. Note this as a deliberate non-port, not an oversight.
3. **Fixed bodies**: Sun through Vesta (including the MPC-number aliases for
   Chiron/Pholus/Ceres/Pallas/Juno/Vesta, e.g. `SE_AST_OFFSET + MPC_CHIRON` == `SE_AST_OFFSET +
   2060` also maps to `"Chiron"`) — direct string constant, switch statement, sweph.c:6965–7039.
4. **Fictitious planets** (sweph.c:7041–7045): `if SE_FICT_OFFSET <= ipl <= SE_FICT_MAX` (40 to
   999) → `swi_get_fict_name(ipl - SE_FICT_OFFSET, s)` (see §3a below).
5. **Asteroids / planetary moons** (sweph.c:7047–7109): `if ipl > SE_PLMOON_OFFSET (9000) || ipl
   > SE_AST_OFFSET (10000)` (the second condition is redundant/subsumed by the first per the C
   comment "obsolete", since anything `> SE_AST_OFFSET` is already `> SE_PLMOON_OFFSET`) — see
   §3b below.
6. **Fallback** (sweph.c:7110–7113): any other `ipl` (i.e. an ID in `1 <= ipl <= SE_PLMOON_OFFSET`
   not matched by any case above — practically this only reaches planetary-moon-range IDs below
   9000 that aren't real bodies) → `sprintf(s, "%d", ipl)`, i.e. the **numeric ID itself**,
   stringified, becomes the "name".
7. **Post-processing** (sweph.c:7120–7123): if the resulting string is under 80 chars, update the
   name cache (see point 2 — not worth porting).

### 3a. Fictitious body names — `swi_get_fict_name` (swemplan.c:513–520)

```c
char *swi_get_fict_name(int32 ipl, char *snam)
{
  if (read_elements_file(ipl, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                          snam, NULL, NULL) == ERR)
    strcpy(snam, "name not found");
  return snam;
}
```
`ipl` here is **already offset-adjusted** (`raw_ipl - SE_FICT_OFFSET`, i.e. 0-based index into
the fictitious-body table). Delegates entirely to `read_elements_file`, requesting **only** the
`pname` output (all the orbital-element out-pointers are `NULL`, meaning `read_elements_file`
skips populating them — see its `if (x != NULL)` guards throughout, swemplan.c:694–908). On
failure, name is hardcoded to the literal string `"name not found"`.

`read_elements_file` (swemplan.c:694–908) — name resolution branches on whether `seorbel.txt`
(`SE_FICTFILE`) is present at all in the ephemeris path (`swi_fopen(-1, SE_FICTFILE, ...)`):

- **File absent** (swemplan.c:707–733): fall back to the **built-in table**
  `plan_fict_nam[SE_NFICT_ELEM]` (15 entries, swemplan.c:506–511, verbatim):
  ```
  "Cupido", "Hades", "Zeus", "Kronos", "Apollon", "Admetos", "Vulkanus", "Poseidon",
  "Isis-Transpluto", "Nibiru", "Harrington", "Leverrier", "Adams", "Lowell", "Pickering"
  ```
  (indices 0–14, corresponding to raw `ipl` 40–54, i.e. `SE_FICT_OFFSET..SE_FICT_OFFSET+14`).
  If `ipl >= SE_NFICT_ELEM` (15) and the file is absent, **error** ("no elements for fictitious
  body no ...") — there is no built-in data for indices ≥15 (e.g. Vulcan=55, White
  Moon/Selena=56, Proserpina=57, Waldemath=58 in this project's `FictitiousBody` enum, which has
  more variants than the C built-in table covers).
- **File present** (swemplan.c:737–908): parses `seorbel.txt` line by line (comment lines start
  with `#`, blank lines skipped), counting non-comment lines (`iplan++`) until `iplan == ipl`;
  the **9th comma-separated field** (`cpos[8]`) on the matching line, trimmed, is the name
  (swemplan.c:882–888). A trailing 10th field (`cpos[9]`) may contain the literal substring
  `"geo"` (case-insensitive) setting `FICT_GEO` — irrelevant to name resolution but shares the
  same line-parsing pass.

  **Critical fidelity detail**: this project's `ephe/seorbel.txt` **exists** (verified:
  `/home/josh/nhs/soft/astrology/swisseph-rs/ephe/seorbel.txt`), which means in this codebase's
  actual runtime behavior the file-present branch is **always** taken — the built-in
  `plan_fict_nam` table is **dead code** whenever a standard `ephe/` directory ships alongside
  the library. This matters because the file's entries can literally differ from the built-in
  array (e.g. seorbel.txt line 7 names index 6 `"Vulcanus"` — same as the built-in table — but a
  user-edited `seorbel.txt` could rename any of them, and the file takes priority unconditionally
  whenever it is found at all, not just for indices beyond the built-in table's 15 entries). A
  faithful Rust port must attempt to read `seorbel.txt` **first**, for every fictitious-body
  index including 0–14, and only fall back to the hardcoded 15-name table if the file genuinely
  cannot be opened — mirroring `swi_fopen`'s search failing entirely (not "field missing" —
  ipl-not-found *within* an opened file is a separate, later error, sweph.c/swemplan.c:899–903,
  "elements for planet ... not found", which does **not** fall back to the built-in table; only
  total file-absence triggers the fallback).

### 3b. Asteroid names — .se1 file MPC header, then `seasnam.txt` fallback (sweph.c:7047–7109)

1. **Already-cached name check** (sweph.c:7049–7050): `swed.fidat[SEI_FILE_ANY_AST].ipl[0] ==
   ipl` → reuse `swed.fidat[SEI_FILE_ANY_AST].astnam` directly (another global-state
   memoization, this time keyed on "which single-asteroid file is currently open" rather than a
   dedicated name cache — **do not port**; this is again purely a perf shortcut with no
   observable effect since it always mirrors what a fresh file read would produce).
2. **Otherwise, trigger an ephemeris read** (sweph.c:7053): `sweph(J2000, ipl, SEI_FILE_ANY_AST,
   0, NULL, NO_SAVE, xp, NULL)` — this opens (or reuses an already-open) `.se1` asteroid file
   purely as a side effect of computing a throwaway position at J2000, and as a side effect of
   *that*, `read_const` (sweph.c:4505–4772, specifically the MPC-header parse at
   sweph.c:4593–4622 and 4709–4753) populates `fdp->astnam` from the file's 4th text line (the
   MPC elements line) or, for older-format files, from a raw 30-byte name field. **This exact
   extraction logic is already ported** — see `src/sweph_file/parse.rs`
   (`extract_asteroid_name`, feeding `AsteroidMeta.name`), documented in full in
   `docs/c-ref-asteroid.md` §3.3. The Rust port of `swe_get_planet_name`'s asteroid branch should
   call into the existing asteroid-file-open path and read `AsteroidMeta.name` — **do not
   reimplement the MPC-header parsing**.
3. **If the file read fails entirely** (sweph.c:7056–7062): synthesize a fallback string —
   `"%d: not found (asteroid)"` (ipl − `SE_AST_OFFSET`) or `"%d: not found (planetary moon)"`
   (raw ipl), depending on which ID range was requested.
4. **`seasnam.txt` override** (sweph.c:7064–7109) — **this is the piece that still needs
   porting** on top of the already-ported file-name extraction. Triggered only when:
   - `ipl > SE_AST_OFFSET` (true asteroids only, not planetary moons), **and**
   - the name just extracted from the `.se1` file starts with `'?'` (old-format placeholder) OR
     its *second character* is a digit (`isdigit(s[1])`) — heuristic for "this looks like a bare
     provisional designation rather than a real name" (e.g. a name literally being a designation
     string like `"1989 UR"` would have `s[1]` = a digit if `s[0]` happens to be a single digit
     too, or more generally provisional designations follow patterns where an early character is
     numeric).

   File format (`seasnam.txt`, `SE_ASTNAMFILE`): plain text, `#`-prefixed comment lines allowed
   anywhere, at least two whitespace/bracket-delimited columns: **(1)** asteroid catalog number
   (optionally wrapped in `(`, `[`, or `{`), **(2)** the name, running to the next `#` or line
   end, right-trimmed. Verified against this project's `ephe/seasnam.txt`:
   ```
   000022
   000048
   000024
   000001  Ceres
   000002  Pallas
   000003  Juno
   000004  Vesta
   000005  Astraea
   ...
   ```
   (leading zero-padded 6-digit numbers; some entries, e.g. `000022`/`000048`/`000024`, have
   **no name** — a lookup for those numbers should leave the previously-extracted `.se1`-file
   name unmodified, matching the C loop: it keeps scanning `while(ipli != iplf && ...)`, and only
   overwrites `s` with `strcpy(s, sp)` when a **non-empty** trimmed name is found after the
   number — if `*sp == '\0'` after skipping whitespace it `continue`s without touching `s`, per
   sweph.c:7102–7103, but note the loop condition `ipli != iplf` means it **keeps reading past** a
   matching-but-nameless line looking for a *later* line with the same number and a name — this
   is unusual and easy to get wrong: the loop does **not** stop at the first `iplf == ipli` match
   if that line has no name field after the number).
   Parsing detail: after matching the catalog number (`atoi`), the pointer advances past the
   number via `strpbrk(sp, " \t")` (first whitespace after the digits), skips further
   whitespace, then everything up to the next `#`/`\r`/`\n` (or end of string) is the name,
   right-trimmed (`swi_right_trim`). If no whitespace exists after the number at all (`sp ==
   NULL`), the line is skipped (`continue`) — "there is no name" on that line.

### Global-state summary for `swe_get_planet_name`
| Global read | Purpose | Port as |
|---|---|---|
| `swed.i_saved_planet_name` / `swed.saved_planet_name` | last-computed-name cache | do not port (no observable effect) |
| `swed.fidat[SEI_FILE_ANY_AST].ipl[0]` / `.astnam` | "is the right asteroid file already open" cache | do not port; always resolve via the (already-ported) file-open + parse path |
| `swed.ephepath` | search directory for `seorbel.txt`/`seasnam.txt` | thread explicitly from `EphemerisConfig.ephe_path` |

## 4. `swe_house_name` (swehouse.c:827–859)

```c
const char *CALL_CONV swe_house_name(int hsys)
{
  int h = hsys;
  if (h != 'i') h = toupper(h);
  switch (h) {
  case 'A': return "equal";
  case 'B': return "Alcabitius";
  case 'C': return "Campanus";
  case 'D': return "equal (MC)";
  case 'E': return "equal";
  case 'F': return "Carter poli-equ.";
  case 'G': return "Gauquelin sectors";
  case 'H': return "horizon/azimut";
  case 'I': return "Sunshine";
  case 'i': return "Sunshine/alt.";
  case 'J': return "Savard-A";
  case 'K': return "Koch";
  case 'L': return "Pullen SD";
  case 'M': return "Morinus";
  case 'N': return "equal/1=Aries";
  case 'O': return "Porphyry";
  case 'Q': return "Pullen SR";
  case 'R': return "Regiomontanus";
  case 'S': return "Sripati";
  case 'T': return "Polich/Page";
  case 'U': return "Krusinski-Pisa-Goelzer";
  case 'V': return "equal/Vehlow";
  case 'W': return "equal/ whole sign";
  case 'X': return "axial rotation system/Meridian houses";
  case 'Y': return "APC houses";
  default: return "Placidus";
  }
}
```
Names are all **lowercase-first-word style** except proper nouns (contrast with the catalogue
doc's Title-Case summaries — the actual return strings are exactly as shown above, e.g.
`"equal"` not `"Equal"`, `"horizon/azimut"` not `"Horizon/Azimuthal"` — **use these literal
strings**, not the catalogue's paraphrased names).

Case handling: `'i'` (lowercase) is a **distinct** house system from `'I'` (Sunshine vs. Sunshine
alt.) — the one case where case matters. Every other input is uppercased via `toupper` before
matching, so `'p'`, `'p'` or even garbage/unmatched codes **all fall through to the `default`
case, returning `"Placidus"`** — there is no error path; any unrecognized `hsys` silently reports
as Placidus.

### Porting note
This project's `HouseSystem` enum (`src/types.rs:272–298`) already has a closed set of variants
and a `to_char()` method (`src/types.rs:301+`) mapping variant → C character code. Because
`HouseSystem` is validated at construction (via whatever `TryFrom<u8>`/parse path already exists
— confirm in `src/types.rs`), the C function's "any garbage char defaults to Placidus" fallback
path is **unreachable** in the Rust port: every `HouseSystem` value is already one of the 24 known
systems. Implement as `impl HouseSystem { pub fn name(self) -> &'static str { ... } }` — a total
function over the closed enum, no `Option`/`Result` needed, and the `default => "Placidus"` arm
simply never triggers because there's no invalid variant to hit it with. Note `'i'`/`'I'`
distinctness is already preserved since `HouseSystem::Sunshine`/`HouseSystem::SunshineAlt` are
separate variants.

## 5. `swe_get_ayanamsa_name` (sweph.c:7127–7133, table sweph.c:130–179)

```c
const char *CALL_CONV swe_get_ayanamsa_name(int32 isidmode)
{
  isidmode %= SE_SIDBITS;         // SE_SIDBITS = 256
  if (isidmode < SE_NSIDM_PREDEF) // SE_NSIDM_PREDEF = 47
    return ayanamsa_name[isidmode];
  return NULL;
}
```

`isidmode %= SE_SIDBITS` strips off any of the sidereal **projection bits**
(`SE_SIDBIT_ECL_T0`=256, `SE_SIDBIT_SSY_PLANE`=512, `SE_SIDBIT_USER_UT`=1024,
`SE_SIDBIT_ECL_DATE`=2048, `SE_SIDBIT_NO_PREC_OFFSET`=4096, `SE_SIDBIT_PREC_ORIG`=8192) that a
raw `sid_mode` value passed to `swe_set_sid_mode` may have OR'd in — this function accepts the
same combined value and needs to recover just the base mode ID (0–46, or 255 for
`SE_SIDM_USER`). **Critical edge case**: `SE_SIDM_USER` = 255 is a **valid** sidereal mode
(`swed.sidd.sid_mode` can legitimately be 255), but `255 % 256 == 255`, and `255 < 47` is
**false** — so `swe_get_ayanamsa_name(SE_SIDM_USER)` returns **`NULL`**, not a
`"User-defined"` string. There is no name-table entry for the user-defined mode at all.

### Name table — verbatim, `ayanamsa_name[47]` (sweph.c:130–178)
```
 0  Fagan/Bradley
 1  Lahiri
 2  De Luce
 3  Raman
 4  Usha/Shashi
 5  Krishnamurti
 6  Djwhal Khul
 7  Yukteshwar
 8  J.N. Bhasin
 9  Babylonian/Kugler 1
10  Babylonian/Kugler 2
11  Babylonian/Kugler 3
12  Babylonian/Huber
13  Babylonian/Eta Piscium
14  Babylonian/Aldebaran = 15 Tau
15  Hipparchos
16  Sassanian
17  Galact. Center = 0 Sag
18  J2000
19  J1900
20  B1950
21  Suryasiddhanta
22  Suryasiddhanta, mean Sun
23  Aryabhata
24  Aryabhata, mean Sun
25  SS Revati
26  SS Citra
27  True Citra
28  True Revati
29  True Pushya (PVRN Rao)
30  Galactic Center (Gil Brand)
31  Galactic Equator (IAU1958)
32  Galactic Equator
33  Galactic Equator mid-Mula
34  Skydram (Mardyks)
35  True Mula (Chandra Hari)
36  Dhruva/Gal.Center/Mula (Wilhelm)
37  Aryabhata 522
38  Babylonian/Britton
39  "Vedic"/Sheoran
40  Cochrane (Gal.Center = 0 Cap)
41  Galactic Equator (Fiorenza)
42  Vettius Valens
43  Lahiri 1940
44  Lahiri VP285
45  Krishnamurti-Senthilathiban
46  Lahiri ICRC
```
Note several strings **differ from the catalogue doc's paraphrased names**, e.g. index 4 is
`"Usha/Shashi"` (not "Ushashashi"), index 34 is `"Skydram (Mardyks)"` (not "Galactic alignment
(Mardyks)"), index 39 is the literal `"\"Vedic\"/Sheoran"` (embedded double-quotes around
"Vedic"), index 42 is `"Vettius Valens"` (not "Valens Moon"), index 45 is
`"Krishnamurti-Senthilathiban"` (not "Krishnamurti VP291"). **Use these exact strings** — they
are the real return values, independent of the more descriptive constant names
(`SE_SIDM_*`) used elsewhere.

### Porting note
This project's `SiderealMode` enum (`src/types.rs:396–445`) already has all 48 variants
(0–46 plus `User = 255`) with a `TryFrom<i32>` impl. Implement `impl SiderealMode { pub fn
name(self) -> Option<&'static str> }`, returning `None` **specifically** for
`SiderealMode::User` (mirroring the C `NULL`), and `Some(&'static str)` for every other variant.
Since the Rust enum is already a closed, validated type, the modulo-256 bit-stripping step is
moot (no projection bits can leak into a `SiderealMode` value) — that logic only existed in C to
cope with the flat `int32` calling convention where projection bits and mode ID share one word.

## 6. `swe_set_delta_t_userdef` (swephlib.c:3176–3184)

```c
void CALL_CONV swe_set_delta_t_userdef(double dt)
{
  if (dt == SE_DELTAT_AUTOMATIC) {
    swed.delta_t_userdef_is_set = FALSE;
  } else {
    swed.delta_t_userdef_is_set = TRUE;
    swed.delta_t_userdef = dt;
  }
}
```
`SE_DELTAT_AUTOMATIC` = `-1E-10` (swephexp.h:497) — a sentinel value chosen because it's an
implausible real deltaT (deltaT is always much larger in magnitude, or positive) — passing it
**clears** the override (reverts to computed deltaT). Any other `dt` value (interpreted as
**days**, not seconds — consistent with every other deltaT value in this library) is stored
verbatim and used **unconditionally** from then on.

### Where the override is consumed — the single choke point
```c
double CALL_CONV swe_deltat_ex(double tjd, int32 iflag, char *serr)
{
  double deltat;
  if (swed.delta_t_userdef_is_set)
    return swed.delta_t_userdef;
  if (serr != NULL)
    *serr = '\0';
  calc_deltat(tjd, iflag, &deltat, serr);
  return deltat;
}
```
(swephlib.c:2701–2710). This is the **only** place the override is checked. Critically:
- The check happens **before** `calc_deltat` is even called — `iflag` (ephemeris source),
  `serr`, and `tjd` are **all ignored** when the override is active. The returned value does not
  depend on which epoch was requested, which ephemeris backend is in use, or anything else.
- **No interaction with tidal acceleration whatsoever**: the override completely bypasses the
  tidal-acceleration-dependent parabolic/tabulated deltaT computation
  (`calc_deltat`/`deltat_aa`, which is where `swi_get_tid_acc`'s ephemeris-specific tidal
  acceleration value would normally factor in). Setting a user-defined deltaT makes
  `swe_set_tid_acc` (and the ephemeris-dependent tidal defaults) **entirely moot** for any code
  path that goes through `swe_deltat_ex` — there is no blending, no partial application, nothing.
- **Every UT↔ET conversion in the library goes through `swe_deltat_ex`** (confirmed by grep:
  `swe_calc_ut`, `swe_houses`/`swe_houses_ex`/`swe_houses_ex2` (swehouse.c:139,220),
  `swe_utc_to_jd`/`swe_jdet_to_utc`/`swe_jdut1_to_utc` family (swedate.c:404,431,464–466,
  496–498,553–554,585), all eclipse/occultation search functions (swecl.c), heliacal-visibility
  functions (swehel.c), and the ayanamsa functions (sweph.c:3264)). There is **no separate
  bypass** — every one of these honors the override identically, because they all funnel through
  this one function. A stateless Rust port that threads `EphemerisConfig` everywhere gets this
  "honored everywhere" property for free, **provided** the check is placed at the single
  `calc_deltat` entry point rather than duplicated at each call site.

### Porting note
`src/deltat/mod.rs:372` already has `pub fn calc_deltat(tjd: f64, config: &EphemerisConfig) ->
f64` as the single dispatcher every call site in this codebase goes through (mirroring C's
`swe_deltat_ex` choke point exactly). The natural port is:
1. Add `pub delta_t_userdef: Option<f64>` to `EphemerisConfig` (`src/config.rs`) — `None` = C's
   `SE_DELTAT_AUTOMATIC`/`delta_t_userdef_is_set == FALSE`; `Some(dt)` = an active override. This
   mirrors the existing `tidal_acceleration: Option<f64>` field's sentinel-free pattern exactly
   (`src/config.rs:30`) — **do not introduce a magic-number sentinel** in the Rust API; the
   `Option` already expresses "automatic vs. overridden" without one.
2. At the very top of `calc_deltat`, short-circuit: `if let Some(dt) = config.delta_t_userdef {
   return dt; }` before any tidal-acceleration/model-selection logic runs — this reproduces the
   "ignores tjd/iflag/tid_acc entirely" behavior exactly, and since every UT↔ET conversion in this
   codebase already calls `calc_deltat(tjd, config)`, the override is automatically honored
   everywhere with a single-line change, matching the "one choke point" property of the C code.
3. There is no separate "setter" function needed the way C has `swe_set_delta_t_userdef` mutating
   global state — in the stateless design, the caller simply sets the field on the
   `EphemerisConfig` they pass to `Ephemeris::new`/wherever configs are threaded, per this
   project's `EphemerisConfig` + `Default` construction pattern (`CLAUDE.md`: "Construction via
   `EphemerisConfig` struct with `Default`. No builder pattern.").

## 7. Centisecond formatting family (swephlib.c:3836–3929, plus 3785–3831 for the related
   normalize/diff helpers)

### Unit: what a "centisecond" (`CSEC`/`centisec`, `typedef int32`) means
```c
#define DEG      360000           /* degree expressed in centiseconds */  (sweodef.h:272)
#define DEG30    (30 * DEG)       /* = 10,800,000 */
#define DEG180   (180 * DEG)      /* = 64,800,000 */
#define DEG360   (360 * DEG)      /* = 129,600,000 */
```
`DEG = 360000` units per "large unit" — this is **1/100 of an arcsecond** when the large unit is
a *degree* (3600 arcsec/degree × 100 = 360000), and **equivalently 1/100 of a time-second** when
the large unit is treated as an *hour* (3600 sec/hour × 100 = 360000) — the same numeric scale
factor serves both interpretations because both "degree" and "hour" divide into 3600
sub-units. There is no type-level distinction between "angle-centiseconds" and "time-centiseconds"
in C — `centisec` is just an `int32`, and which interpretation applies is determined entirely by
**which formatting function is called** and what the caller put into it:
- `swe_cs2timestr` treats its input as **time-of-day** centiseconds, wrapping at 24 hours
  (`% (24*3600)` **seconds**, after first dividing the raw centisecond value by 100).
- `swe_cs2degstr`/`swe_cs2lonlatstr` treat input as **arc** centiseconds (degrees), with
  `swe_csnorm`/`DEG`-family constants defining the degree-based wraparound.

A porter should **not** create a single `Centiseconds` newtype conflating both semantics without
a clear doc comment — consider two distinct wrapper types (e.g. `ArcCentiseconds`,
`TimeCentiseconds`) or a single type with the unit made explicit at each call site, since mixing
them up produces silently wrong results (same underlying `i32` representation, incompatible
meaning).

### `swe_csnorm` (swephlib.c:3785–3792) — normalize to `[0, DEG360)`
```c
centisec CALL_CONV swe_csnorm(centisec p)
{
  if (p < 0)
    do { p += DEG360; } while (p < 0);
  else if (p >= DEG360)
    do { p -= DEG360; } while (p >= DEG360);
  return (p);
}
```
Loop-based (not `%`) normalization — for `int32` centiseconds this only matters for inputs many
multiples of `DEG360` away from `[0, DEG360)`; a Rust port can use `rem_euclid(DEG360)` for the
same result more efficiently, since the loop is doing exactly what Euclidean modulo does for
integer types (no FP rounding concerns unlike the `double`-based `swe_degnorm`, which is already
ported in `src/math.rs` and has its own `fabs(y) < 1e-13 -> 0` snap-to-zero guard that does
**not** apply here since `centisec` is an exact integer type).

### `swe_csroundsec` (swephlib.c:3836–3843) — round to nearest arcsecond, with a sign-boundary guard
```c
centisec CALL_CONV swe_csroundsec(centisec x)
{
  centisec t;
  t = (x + 50) / 100 * 100L;         /* round to nearest 100 (= nearest arcsecond) */
  if (t > x && t % DEG30 == 0)       /* rounded UP across a 30-degree (zodiac sign) boundary */
    t = x / 100 * 100L;             /* round DOWN instead, to the last whole arcsecond in-sign */
  return (t);
}
```
Rounds to the nearest whole arcsecond (nearest multiple of 100 centiseconds) using integer
truncating division with a `+50` bias (standard round-half-up for positive values; **for
negative `x`, integer division truncation direction matters** — C's `/` truncates toward zero,
so `(x + 50) / 100` for negative `x` does not behave like a symmetric round-half-up; a Rust port
must replicate C's truncating-toward-zero integer division exactly, not use `div_euclid` or
float-based rounding, to stay bit-for-bit faithful for negative centisecond values).

The guard clause is the interesting part: if rounding **up** pushed the value onto an exact
30-degree boundary (`t % DEG30 == 0`, i.e. `t` landed exactly on a zodiac sign cusp, since
`DEG30` = one zodiac sign's width in centiseconds) **and** the rounding direction was actually
upward (`t > x`), the function instead rounds **down** to the last whole arcsecond *before* that
boundary — i.e. `29°59'59.6"` rounds to `29°59'60"` = `30°00'00"` normally, but this function
special-cases exactly that scenario to stay at `29°59'59"` instead, keeping the value within the
original zodiac sign rather than spilling into the next one. This only triggers when the rounded
result lands **exactly** on a multiple of 30 degrees — an extremely narrow window (within 0.5
arcsec of a sign boundary) but a real, documented behavior (catalogue doc calls it out too: "Round
to nearest second; rounds DOWN at 29°59′30″ to stay within sign").

### `swe_cs2timestr` (swephlib.c:3864–3886) — `"HH:MM:SS"` (or with a custom separator)
```c
char *CALL_CONV swe_cs2timestr(CSEC t, int sep, AS_BOOL suppressZero, char *a)
{
  centisec h, m, s;
  strcpy(a, "        ");            /* 8 spaces, pre-fills the 8-char buffer incl. NUL slot */
  a[2] = a[5] = sep;                 /* separator chars at fixed positions */
  t = ((t + 50) / 100) % (24L * 3600L); /* round to whole seconds, then wrap at 24h */
  s = t % 60L;
  m = (t / 60) % 60L;
  h = t / 3600 % 100L;               /* NOTE: allows up to 99 "hours" (no true 24h clamp on h) */
  if (s == 0 && suppressZero)
    a[5] = '\0';                     /* cuts the string short right after MM, dropping ":SS" entirely */
  else {
    a[6] = (char)(s / 10 + '0');
    a[7] = (char)(s % 10 + '0');
  }
  a[0] = (char)(h / 10 + '0');
  a[1] = (char)(h % 10 + '0');
  a[3] = (char)(m / 10 + '0');
  a[4] = (char)(m % 10 + '0');
  return (a);
}
```
- **Rounds** to the nearest whole second (`+50`, truncating divide by 100), **then** wraps modulo
  86400 seconds (24h) — so the *rounding* happens before the *wraparound*; a value that rounds up
  to exactly 24:00:00 wraps to `00:00:00`, not `24:00:00`.
- `sep` is a **raw `int`** written directly into `char` array slots 2 and 5 — typically `':'`
  (0x3A), but passing `0` produces a NUL byte at both separator positions, which (combined with
  `strcpy`'s original all-space buffer) would produce a string that *looks* like `"HH\0MM\0SS"` in
  memory but reads as just `"HH"` through any C-string API — **an easy foot-gun if `sep` is ever
  0**; a Rust port should take a `char`/`u8` separator argument, not silently allow embedding NUL.
- `suppressZero` **only ever suppresses the seconds field** (never hours or minutes) — the C
  comment says exactly this ("does not suppress zeros in hours or minutes"). When active and
  `s == 0`, the string is truncated right after the minutes field (position 5 becomes the
  terminator, discarding the separator that would have preceded seconds too) — output is
  `"HH:MM"` (5 chars + implicit-from-buffer trailing spaces that are irrelevant since the NUL cuts
  the string there).
- `h` can range 0–99 (`% 100L`) — this function is **not** exclusively for 24-hour clock times;
  it's reused elsewhere in the library for arbitrary HH:MM:SS-formatted durations/values that can
  exceed 24 in the hours field (though the `% (24L*3600L)` wrap on the raw `t` value *does* impose
  a 24-hour wraparound on the total value before splitting into h/m/s — so in practice `h` will
  never exceed 23 through this code path; the `%100L` on `h` is defensive/vestigial given the
  prior wrap already bounds `t` to `< 86400`, making `h < 24` always).

### `swe_cs2lonlatstr` (swephlib.c:3888–3916) — `"DDD°MM'SS\""` + direction letter, up to 999°, no wraparound
```c
char *CALL_CONV swe_cs2lonlatstr(CSEC t, char pchar, char mchar, char *sp)
{
  char a[10];
  char *aa;
  centisec h, m, s;
  strcpy(a, "      '  ");           /* mask: dddEmm'ss" -- note embedded literal apostrophe at [6] */
  if (t < 0) pchar = mchar;         /* negative input -> use the "minus" direction letter instead */
  t = (ABS4(t) + 50) / 100;         /* round to whole arcseconds; ABS4 = labs (long abs) */
  s = t % 60L;
  m = t / 60 % 60L;
  h = t / 3600 % 1000L;             /* up to 999 degrees -- NO modulo-360 wraparound at all */
  if (s == 0)
    a[6] = '\0';                    /* cut off seconds (always, unconditionally -- no suppressZero param here) */
  else {
    a[7] = (char)(s / 10 + '0');
    a[8] = (char)(s % 10 + '0');
  }
  a[3] = pchar;                     /* direction letter goes at position 3, between degrees and minutes */
  if (h > 99) a[0] = (char)(h / 100 + '0');
  if (h > 9)  a[1] = (char)(h % 100 / 10 + '0');
  a[2] = (char)(h % 10 + '0');
  a[4] = (char)(m / 10 + '0');
  a[5] = (char)(m % 10 + '0');
  aa = a;
  while (*aa == ' ') aa++;          /* skip leading blank degree-digit positions */
  strcpy(sp, aa);
  return (sp);
}
```
- Takes the **sign of `t`** to pick between `pchar` (positive direction letter, e.g. `'E'` for
  east longitude or `'N'` for north latitude) and `mchar` (negative direction letter, `'W'`/`'S'`)
  — the direction letter itself is then always printed, and the numeric part is the **absolute
  value** (`ABS4`).
- **Always** drops the seconds field when `s == 0` — unlike `cs2timestr`, there's no
  `suppressZero` parameter; seconds suppression on exact-zero is unconditional here.
- **Leading-zero suppression on the degrees field only** happens through the "skip leading
  spaces" scan (`while (*aa == ' ') aa++`) combined with the buffer being pre-filled with spaces
  and only positions `[0]`/`[1]` conditionally written when `h > 99`/`h > 9` — so e.g. `h = 5`
  leaves `a[0]=a[1]=' '` (untouched from the initial `strcpy`), and only `a[2]='5'` is written;
  the leading-space-skip then makes the direction-letter position `a[3]` (`pchar`/`mchar`) the
  first character the final string exposes if degrees is a single digit... **wait**: re-examine
  — `aa` skips leading spaces starting from `a[0]`; if `h < 10`, `a[0]` and `a[1]` are both still
  `' '` (never written), so `aa` advances past both, landing on `a[2]` (the single ones-digit of
  degrees) as the first real character — meaning single-digit degrees are NOT zero-padded and the
  direction letter correctly follows immediately after the digit(s), e.g. `"5°30'00"E"` not
  `"005°30'00\"E"`. Degrees up to 999 are supported (`h > 99` writes the hundreds digit, `h > 9`
  writes the tens digit), consistent with this being used for absolute ecliptic longitude-like
  values that can exceed 360 in some contexts (or simply never gets modulo-wrapped — caller is
  responsible for normalizing to `[0,360)` beforehand if that's desired; this function will
  happily print `450°` as `"450°..."` verbatim).
- Output buffer is only 10 bytes (`char a[10]`) — `"XXX°MM'SS\"X\0"` — wait, actually there's no
  degree symbol character in this format at all (contrast with `cs2degstr` which uses
  `ODEGREE_STRING`) — the raw layout is `d d d P m m ' s s` (10 slots: indices 0-2 degrees, 3
  direction char, 4-5 minutes, 6 apostrophe-or-NUL, 7-8 seconds) with **no closing quote/degree
  symbol appended** — the catalogue doc's example `"122°30'45"E"` appears to show a degree symbol
  and closing double-quote that **do not actually appear in the C output**; the real format is
  closer to `122E30'45` (direction letter immediately after degrees, apostrophe before seconds,
  no unit symbols) — verify this precisely against the actual C string layout
  (`"      '  "` is exactly 9 characters + NUL = the 10-byte `a[10]`) when implementing; do not
  trust the catalogue doc's example string for exact punctuation, trust the `strcpy(a, "...")`
  mask literal instead.

### `swe_cs2degstr` (swephlib.c:3918–3929) — `"DD°MM'SS"` within a single sign, **truncates** (does not round)
```c
char *CALL_CONV swe_cs2degstr(CSEC t, char *a)
{
  centisec h, m, s;
  t = t / 100 % (30L * 3600L);      /* TRUNCATE to whole arcseconds, then wrap into [0, 30 deg) */
  s = t % 60L;
  m = t / 60 % 60L;
  h = t / 3600 % 100L;              /* defensive mask; already < 30 from the prior wrap */
  sprintf(a, "%2d%s%02d'%02d", h, ODEGREE_STRING, m, s);
  return (a);
}
```
- **No `+50` rounding bias** — this is the one formatter in the family that **truncates** rather
  than rounds (contrast explicitly with `cs2timestr`/`cs2lonlatstr`, both of which add 50 before
  the integer divide). Must not add a rounding step when porting this one.
- Wraps into `[0, 30*3600)` **arcseconds**, i.e. always displays degrees-within-sign (0–29°), not
  the full 0–359° range — intended for "position within zodiac sign" display, hence the name
  `cs2degstr` (as opposed to `cs2lonlatstr`'s full-range longitude/latitude use).
  `ODEGREE_STRING` = `"°"` (UTF-8 degree sign, sweodef.h:246).
- Output format via `sprintf("%2d%s%02d'%02d", ...)`: degrees is `%2d` (space-padded to 2 chars,
  **not** zero-padded — `" 5°30'00"` for `h=5`, not `"05°30'00"`), minutes/seconds are `%02d`
  (zero-padded to 2 digits each). No leading/trailing sign-direction letter (unlike
  `cs2lonlatstr`) since this is purely "offset within sign", not an absolute geographic
  coordinate.

### Related normalize/diff helpers, for completeness (already summarized in the catalogue; verbatim here since this doc is the authoritative source for the centisecond family)
```c
centisec swe_difcsn(centisec p1, centisec p2)  { return swe_csnorm(p1 - p2); }         // [0, 360deg)
centisec swe_difcs2n(centisec p1, centisec p2) { d = swe_csnorm(p1-p2); return d>=DEG180 ? d-DEG360 : d; } // [-180,180)deg
double   swe_difdegn(double p1, double p2)     { return swe_degnorm(p1 - p2); }        // degree analog of difcsn
double   swe_difdeg2n(double p1, double p2)    { d = swe_degnorm(p1-p2); return d>=180.0 ? d-360.0 : d; }  // already ported, src/math.rs
double   swe_difrad2n(double p1, double p2)    { d = swe_radnorm(p1-p2); return d>=PI ? d-TWOPI : d; }     // radian analog
int32    swe_d2l(double x) { return x>=0 ? (int32)(x+0.5) : -(int32)(0.5-x); }         // round-half-away-from-zero, no overflow check
```
`swe_difdeg2n` is already ported (`src/math.rs`, per `docs/c-ref-crossings.md`'s Function Map).
`swe_difcsn`/`swe_difcs2n` are the integer-centisecond analogs and are new for this doc's scope
if the centisecond family is ported at all — but note **none of the callers surveyed in this
codebase's C reference docs so far actually need `centisec`-typed arithmetic** (all internal math
in this project uses `f64` degrees throughout, per `CLAUDE.md`'s architecture notes) — these only
matter if/when the formatting functions themselves are ported, as inputs to them.

## Porting notes

- **`swe_get_planet_name` is a dispatcher over already-ported pieces**: the `.se1`-file asteroid
  name extraction (`extract_asteroid_name`/`AsteroidMeta.name`, `src/sweph_file/parse.rs`) is
  done. What's missing is (a) the top-level name-by-`Body`-variant switch for the fixed bodies,
  (b) the `seorbel.txt`-based fictitious-name lookup (§3a — note the **file-always-wins**
  behavior since this repo ships `ephe/seorbel.txt`, making the 15-entry built-in table
  effectively dead code here), and (c) the `seasnam.txt` fallback-name override for asteroids
  whose `.se1`-embedded name looks like a bare provisional designation (§3b). None of (a)/(b)/(c)
  should attempt to reproduce the two C-side memoization caches (`swed.i_saved_planet_name`,
  `swed.fidat[...].astnam`-as-cache) — both are pure performance shortcuts with zero effect on
  the returned string, confirmed by tracing what each cache is keyed on and what it's populated
  from.
- **`swi_get_fict_name`/`read_elements_file` share one file-parsing pass with the *orbital
  elements* computation** (`swi_osc_el_plan`, swemplan.c:579+ — not covered by this doc, but
  visible in the code read above) — both the name and the 7 numeric orbital elements come from
  the *same* line of `seorbel.txt`/the same built-in table row. If a future ref doc covers
  fictitious-planet position computation, the file-parsing logic should be a **single shared
  helper** returning both the elements and the name together (per this project's constraint
  against duplicating logic) rather than parsing the file twice (once for elements, once for the
  name) the way the C code's separate `swi_get_fict_name` call effectively does when a caller
  wants *just* the name.
- **`swe_house_name`/`swe_get_ayanamsa_name` are simple total/partial functions over already-
  ported closed enums** (`HouseSystem`, `SiderealMode`) — no file I/O, no global state, pure
  string tables. `house_name` becomes a total function (`HouseSystem` has no invalid-code case to
  fall back from); `ayanamsa_name` becomes `Option`-returning specifically because
  `SiderealMode::User` has no table entry (mirroring C's `NULL` for `SE_SIDM_USER`) — this is the
  one genuine "partial function" edge case in an otherwise total mapping, worth a unit test
  (`SiderealMode::User.name() == None`).
- **`swe_set_delta_t_userdef` needs no mutator function in the stateless design** — just a new
  `Option<f64>` field on `EphemerisConfig` (mirroring the existing `tidal_acceleration` field's
  pattern exactly) plus a one-line short-circuit at the top of the existing `calc_deltat`
  dispatcher (`src/deltat/mod.rs:372`). Because every UT↔ET conversion in this codebase already
  funnels through `calc_deltat`, this single change makes the override "honored everywhere"
  automatically — matching the C library's single-choke-point design without needing to hunt down
  every call site individually. **Do not** invent a sentinel float value (`SE_DELTAT_AUTOMATIC` =
  `-1e-10`) in the Rust API; `Option<f64>` already expresses the automatic/overridden distinction
  cleanly and the C sentinel exists only because C has no `Option` type.
- **Centisecond formatting is the one area with a genuine FP/integer-fidelity hazard**: three of
  the four formatters round (`+50` bias before truncating divide) and one truncates
  (`cs2degstr`) — get this backwards and every value within half an arcsecond of a boundary will
  differ. `swe_csroundsec`'s sign-boundary special case (round-up-across-a-30°-cusp gets pulled
  back down) is a narrow-but-real behavior that a naive "round to nearest 100" port would miss
  entirely. `swe_cs2lonlatstr`'s exact output layout (no degree symbol, no closing quote, despite
  what the catalogue doc's example string suggests) should be re-verified against the literal
  `strcpy(a, "      '  ")` mask and field-index writes, not against prose descriptions, when
  implementing — write a golden test that captures the raw byte layout before trusting any
  paraphrase (including this doc's own prose above).
- **Two distinct centisecond semantics share one C type** (`typedef int32 centisec`): time-of-day
  centiseconds (wraps at 24h, consumed by `cs2timestr`) vs. arc-degree centiseconds (wraps at
  360° or 30°, consumed by `csnorm`/`cs2degstr`/`cs2lonlatstr`). A Rust port should make this
  distinction visible in the type system (two newtypes, or explicit unit documentation at each
  call site) rather than reusing a single bare `i32`/`f64` the way C does — this project's
  existing convention of using distinct newtypes for distinct ID spaces (`FictitiousId`,
  `AsteroidId`, `PlanetMoonId` in `src/types.rs`) suggests the same treatment would fit here if
  these formatters are ported with real call sites (as opposed to being pure leaf
  utility/display functions with no internal callers, in which case a bare `f64` degrees /
  `f64` days input — converted to centiseconds only at the function boundary — may be simpler
  and more consistent with "all internal math is `f64` degrees" per `CLAUDE.md`).
- **`swe_time_equ`'s `swed.jpl_file_is_open` read cannot be reproduced statelessly** — a Rust
  port must take the ephemeris source explicitly (from `EphemerisConfig`/the caller's flags) and
  must NOT silently prefer JPL just because some *other*, unrelated call earlier in the process
  happened to open a JPL file. This is a deliberate behavioral divergence from C, not an
  oversight — document it at the call site the way other stateless-vs-stateful divergences are
  documented per `CLAUDE.md`'s `<stateless_tolerance>` section (though this one is a source-
  selection difference, not a numerical-precision one, so it may not need golden-test tolerance
  relaxation — it needs the caller to pass the right flag instead).
- **`swe_lat_to_lmt`'s fixed 3-call iteration is a magic number, not a convergence loop** — port
  literally as 3 sequential calls, not a `while |delta| > eps` loop, for bit-for-bit fidelity
  (the C code has no tolerance check at all).
