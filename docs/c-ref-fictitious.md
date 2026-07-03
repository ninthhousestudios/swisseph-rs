# C Reference: Fictitious Planets (Uranian/Hamburg bodies, Transpluto, etc.) — swemplan.c / sweph.c

Porting reference for the code path behind `swe_calc()`/`swe_calc_ut()` for bodies
`SE_FICT_OFFSET` (40) .. `SE_FICT_OFFSET + SE_NFICT_ELEM - 1` and beyond, up to `SE_FICT_MAX`
(999): Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon, Isis-Transpluto,
Nibiru, Harrington, Neptune-Leverrier, Neptune-Adams, Pluto-Lowell, Pluto-Pickering, Vulcan,
White Moon/Selena, Proserpina, Waldemath — plus any user-added rows in `seorbel.txt`. Read this
instead of the C source.

## Function Map

| C function | Location | Port? |
|---|---|---|
| `swi_osc_el_plan` | swemplan.c:579–690 | Yes — Kepler-element → J2000-equatorial-barycentric position/velocity |
| `read_elements_file` | swemplan.c:694–913 (static) | Yes — seorbel.txt parser + built-in table fallback |
| `check_t_terms` | swemplan.c:916–967 (static) | Yes — polynomial-in-T expression evaluator |
| `swi_get_fict_name` | swemplan.c:513–520 | Yes — thin wrapper, name lookup only |
| `swi_kepler` | swephlib.c:4065–4096 | Yes — Kepler equation solver (shared with other callers; **not yet ported**, grep confirms no existing Rust equivalent) |
| dispatch in `swecalc()`-equivalent | sweph.c:1106–1136 | Yes — `ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX` branch |
| `app_pos_etc_plan_osc` | sweph.c:3365–3547 | Yes — apparent-position pipeline specific to osculating-element bodies |
| `app_pos_rest` | sweph.c:2777–2859 | Shared tail (nutation, ecliptic transform, sidereal, polar/degrees) — likely already ported as part of the existing `apparent_planet`/`app_pos_rest`-equivalent in `calc.rs`; confirm before duplicating |

## Constants (swephexp.h / sweph.h)

```c
#define SE_FICT_OFFSET      40      // swephexp.h:131 — first fictitious-body ipl
#define SE_FICT_OFFSET_1    39      // swephexp.h:132 — used only in the seorbel.txt comment formula
#define SE_FICT_MAX         999     // swephexp.h:133 — highest legal ipl (960 max user-added rows per file header comment)
#define SE_NFICT_ELEM        15     // swephexp.h:134 — count of BUILT-IN table rows (Cupido..Pickering)
#define SE_NALL_NAT_POINTS  (SE_NPLANETS + SE_NFICT_ELEM)  // swephexp.h:138
```

`SE_CUPIDO`=40 … `SE_PLUTO_PICKERING`=54 (rows 0–14 of `plan_oscu_elem[]`), `SE_VULCAN`=55,
`SE_WHITE_MOON`=56, `SE_PROSERPINA`=57, `SE_WALDEMATH`=58 (swephexp.h:141–160). **Bodies 55–58
have NO built-in table row** — `SE_NFICT_ELEM` is 15, not 19. They exist only as rows 16–19 of
`seorbel.txt` (`ephe/seorbel.txt`, see §3). If `seorbel.txt` is missing, `ipl >= SE_NFICT_ELEM`
(i.e. ipl-40 ≥ 15, meaning Vulcan/White-Moon/Proserpina/Waldemath) hits the `ERR` branch in
`read_elements_file` (swemplan.c:709–713) — **these four bodies are unusable without the file**,
unlike the first 15 which have a built-in fallback.

Other constants used: `J1900`=2415020.0, `J2000`=2451545.0, `B1950`=2433282.42345905
(sweph.h:67–69); `DEGTORAD`=π/180, `RADTODEG`=180/π (sweodef.h:265–266); `AUNIT`=1.49597870700e11 m
(DE431 value, sweph.h:273); `CLIGHT`=2.99792458e8 m/s (sweph.h:274); `KGAUSS`=0.01720209895
(Gaussian gravitational constant, sweph.h:280); `SUN_EARTH_MRAT`=332946.050895 (Sun/Earth-only
mass ratio, AA 2006 K7, sweph.h:264 — note a commented-out Sun/(Earth+Moon) variant at sweph.h:263
is dead code, not used).

```c
#define FICT_GEO 1                                  // swemplan.c:71 — bit in fict_ifl
#define KGAUSS_GEO 0.0000298122353216                // swemplan.c:72 — "Earth only" Gaussian const for geocentric fictitious bodies
/* #define KGAUSS_GEO 0.00002999502129737 */          // swemplan.c:73 — commented-out Earth+Moon variant, dead code
```

## §1. The built-in table `plan_oscu_elem[SE_NFICT_ELEM][8]` (swemplan.c:522–571)

Columns: **epoch** (JD or `J1900` sentinel value), **equinox** (JD), **mean anomaly at epoch**
(deg), **semi-major axis** (AU), **eccentricity**, **argument of perihelion** (deg), **ascending
node** (deg), **inclination** (deg). All angle columns are plain degrees in the table; conversion
to radians happens in `read_elements_file`, not here.

Two variants exist behind `#ifdef SE_NEELY` (swemplan.c:504 — **this macro is always defined**,
so the `#else` branch, rows without Neely's revisions, is **dead code** in any standard build;
transcribe both since a port might expose a config toggle, but the live values are the Neely set):

### Neely-revised Uranian ("Witte/Sieggrün") elements — LIVE (swemplan.c:524–531)

| # | ipl | Name | Epoch | Equinox | Mean anomaly | Semi-axis (AU) | Eccentricity | Arg. peri. | Asc. node | Inclination |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 40 | Cupido | J1900 | J1900 | 163.7409 | 40.99837 | 0.00460 | 171.4333 | 129.8325 | 1.0833 |
| 1 | 41 | Hades | J1900 | J1900 | 27.6496 | 50.66744 | 0.00245 | 148.1796 | 161.3339 | 1.0500 |
| 2 | 42 | Zeus | J1900 | J1900 | 165.1232 | 59.21436 | 0.00120 | 299.0440 | 0.0000 | 0.0000 |
| 3 | 43 | Kronos | J1900 | J1900 | 169.0193 | 64.81960 | 0.00305 | 208.8801 | 0.0000 | 0.0000 |
| 4 | 44 | Apollon | J1900 | J1900 | 138.0533 | 70.29949 | 0.00000 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 45 | Admetos | J1900 | J1900 | 351.3350 | 73.62765 | 0.00000 | 0.0000 | 0.0000 | 0.0000 |
| 6 | 46 | Vulkanus | J1900 | J1900 | 55.8983 | 77.25568 | 0.00000 | 0.0000 | 0.0000 | 0.0000 |
| 7 | 47 | Poseidon | J1900 | J1900 | 165.5163 | 83.66907 | 0.00000 | 0.0000 | 0.0000 | 0.0000 |

Note row 3 (Kronos): the built-in table's semi-axis is `64.81960`, but `ephe/seorbel.txt` line 31
has `64.81690` for the same body — **the file and the built-in table disagree** in the 4th decimal
(`64.81960` vs `64.81690`). Since `seorbel.txt` exists in every real deployment (`ephe/` directory
ships with it), the file's value wins for any installation that has it; the built-in table is only
reached when the file is absent entirely. Transcribe the built-in value in the Rust constant table
exactly as printed above — do not "fix" it to match the file.

### Non-Neely elements — DEAD CODE, `#else` branch (swemplan.c:533–540, never compiled with `SE_NEELY` defined)

| # | ipl | Name | Epoch | Equinox | Mean anomaly | Semi-axis (AU) | Ecc. | Arg. peri. | Node | Incl. |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 40 | Cupido | J1900 | J1900 | 104.5959 | 40.99837 | 0 | 0 | 0 | 0 |
| 1 | 41 | Hades | J1900 | J1900 | 337.4517 | 50.667443 | 0 | 0 | 0 | 0 |
| 2 | 42 | Zeus | J1900 | J1900 | 104.0904 | 59.214362 | 0 | 0 | 0 | 0 |
| 3 | 43 | Kronos | J1900 | J1900 | 17.7346 | 64.816896 | 0 | 0 | 0 | 0 |
| 4 | 44 | Apollon | J1900 | J1900 | 138.0354 | 70.361652 | 0 | 0 | 0 | 0 |
| 5 | 45 | Admetos | J1900 | J1900 | -8.678 | 73.736476 | 0 | 0 | 0 | 0 |
| 6 | 46 | Vulkanus | J1900 | J1900 | 55.9826 | 77.445895 | 0 | 0 | 0 | 0 |
| 7 | 47 | Poseidon | J1900 | J1900 | 165.3595 | 83.493733 | 0 | 0 | 0 | 0 |

**Porting decision: only transcribe the Neely (live) set as the Rust built-in constant table.**
The non-Neely set is unreachable in the reference C build and cannot be golden-tested against it;
including it would be dead weight. Flag this explicitly in code comments so a future reader doesn't
wonder why only 8 rows of "the C table" were ported.

### Remaining 7 built-in rows — same regardless of `SE_NEELY` (swemplan.c:542–560)

| # | ipl | Name | Epoch (JD) | Equinox (JD) | Mean anomaly | Semi-axis (AU) | Ecc. | Arg. peri. | Node | Incl. |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 48 | Isis-Transpluto | 2368547.66 | 2431456.5 | 0.0 | 77.775 | 0.3 | 0.7 | 0 | 0 |
| 9 | 49 | Nibiru | 1856113.380954 | 1856113.380954 | 0.0 | 234.8921 | 0.981092 | 103.966 | -44.567 | 158.708 |
| 10 | 50 | Harrington | 2374696.5 | J2000 (2451545.0) | 0.0 | 101.2 | 0.411 | 208.5 | 275.4 | 32.4 |
| 11 | 51 | Leverrier (Neptune) | 2395662.5 | 2395662.5 | 34.05 | 36.15 | 0.10761 | 284.75 | 0 | 0 |
| 12 | 52 | Adams (Neptune) | 2395662.5 | 2395662.5 | 24.28 | 37.25 | 0.12062 | 299.11 | 0 | 0 |
| 13 | 53 | Lowell (Pluto) | 2425977.5 | 2425977.5 | 281 | 43.0 | 0.202 | 204.9 | 0 | 0 |
| 14 | 54 | Pickering (Pluto) | 2425977.5 | 2425977.5 | 48.95 | 55.1 | 0.31 | 280.1 | 100 | 15 |

Row 9 (Nibiru) has **eccentricity 0.981092 > 0.975** — this triggers the high-eccentricity
initial-guess branch in `swi_osc_el_plan` (§4, and see the FP-fidelity hazard in Porting Notes).

Two additional commented-out rows exist at swemplan.c:561–570 inside `#if 0` (a JPL-fit Ceres
initial-elements experiment and a Chiron Bowell-database entry) — **not compiled, not part of any
live table, do not port**.

`plan_fict_nam[SE_NFICT_ELEM]` (swemplan.c:506–511) — the 15 built-in names, index-parallel to
`plan_oscu_elem`: `"Cupido", "Hades", "Zeus", "Kronos", "Apollon", "Admetos", "Vulkanus",
"Poseidon", "Isis-Transpluto", "Nibiru", "Harrington", "Leverrier", "Adams", "Lowell",
"Pickering"`. Note these differ cosmetically from the Rust `FictitiousBody` enum names (e.g. C
`"Leverrier"` vs Rust `NeptuneLeverrier`) — this table is used only for the human-readable name
returned by `swe_get_planet_name`/`swi_get_fict_name`, not for dispatch.

## §2. `swi_get_fict_name` (swemplan.c:513–520)

```c
char *swi_get_fict_name(int32 ipl, char *snam) {
  if (read_elements_file(ipl, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                          snam, NULL, NULL) == ERR)
    strcpy(snam, "name not found");
  return snam;
}
```
Calls `read_elements_file` with every numeric out-param `NULL` except `pname` — i.e. it reuses
the exact same file-search/fallback logic as position computation purely to extract column 9
(the name). `tjd` is passed as `0`, which is harmless because `tjd` only matters for T-term
evaluation of the *numeric* columns (none of which are requested here).

## §3. `seorbel.txt` format and `read_elements_file` (swemplan.c:694–913)

### File location & fallback

`swi_fopen(-1, SE_FICTFILE, swed.ephepath, serr)` (swemplan.c:707) — searched in the configured
ephemeris path(s); `SE_FICTFILE` = `"seorbel.txt"` (per catalogue-public.md:951). The `-1` file-slot
argument means **the file handle is never cached** — it is reopened and closed (line 908/911) on
every single call. If the file cannot be opened at all, C falls back to the built-in table (§1) —
but only for `ipl < SE_NFICT_ELEM` (swemplan.c:709); for `ipl >= 15` (Vulcan and beyond) this is a
hard error: `"error no elements for fictitious body no %7.0f"` (swemplan.c:711).

If the file *is* found but the given row/`ipl` isn't present, or a line is malformed, this is
always an error too (`elem_found` stays `FALSE` → `"elements for planet %7.0f not found"`,
swemplan.c:902–906) — **there is no partial fallback to the built-in table once the file is
successfully opened**, even for `ipl < 15`. Opening the file at all commits you to the file's
data exclusively for that call.

### Line syntax

Real content lines are comma-separated with **exactly 9 required fields** (a 10th optional field
for the `geo` flag):

```
epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name[, geo]
```

Parsing (swemplan.c:739–763):
1. Read line via `fgets`; strip leading spaces/tabs (`sp` advanced past ` `/`\t`, then
   `swi_strcpy(s, sp)` — collapses leading whitespace only, not internal).
2. Skip a line entirely if, after that trim, the first character is `#`, `\r`, `\n`, or `\0`
   (blank/comment lines).
3. Truncate at an embedded `#` (`strchr(s, '#')` → `*sp = '\0'`) — **inline trailing comments
   are supported** (e.g. `... Cupido   # 1` in the real file, swemplan.c:753–754).
4. `swi_cutstr(s, ",", cpos, 20)` splits on commas into up to 20 fields (`cpos[0..ncpos-1]`).
5. If `ncpos < 9` → error `"nine elements required"` (swemplan.c:757–761) for **that row** — note
   this check happens for every non-blank/non-comment line regardless of which `ipl` is being
   searched for, so a malformed row anywhere in the file before the target row aborts the whole
   lookup.
6. Row index `iplan` increments once per successfully-parsed line (`iplan++` at line 763) — **this
   is a zero-based sequential counter over content lines, not tied to any numeric ID printed in the
   file**. `ipl = SE_FICT_OFFSET_1 + number_of_elements_set` per the file's own header comment
   (swemplan.c / seorbel.txt:24–25), i.e. row *N* (1-based, "# N" trailing comment is purely
   documentation) maps to `ipl = 39 + N`. If `iplan != ipl` the row is skipped (`continue`,
   swemplan.c:764–765) — parsing continues to completion (or first error) even after a match is
   NOT yet found; once `iplan == ipl`, `elem_found = TRUE` and the row's 9 columns are parsed,
   then `break` (line 900) — **no further lines are read** once the target row is located, even if
   there are formatting errors later in the file.

### Column 1 — epoch (swemplan.c:768–786)

```c
sp = cpos[0];
for (i = 0; i < 5; i++) sp[i] = tolower(sp[i]);   // lower-cases first 5 chars UNCONDITIONALLY
if (strncmp(sp, "j2000", 5) == OK) *tjd0 = J2000;
else if (strncmp(sp, "b1950", 5) == OK) *tjd0 = B1950;
else if (strncmp(sp, "j1900", 5) == OK) *tjd0 = J1900;
else if (*sp == 'j' || *sp == 'b') { /* error: invalid epoch */ }
else *tjd0 = atof(sp);
tt = tjd - *tjd0;   // used by ALL subsequent check_t_terms() calls on this row
```
**Hazard**: `sp[i] = tolower(sp[i])` for `i` in `0..5` runs even if the field is shorter than 5
characters or is a bare number like `"2368547.66"` — this mutates in place up to 5 bytes past the
field start, which is safe here only because `cpos[]` entries point into the shared line buffer
`s` (there is always at least a comma or more numeric text following in-bounds). A Rust port
should NOT reproduce blind 5-byte lowercasing on a raw slice; instead branch on
`field.trim().len() >= 5 && field[..5].eq_ignore_ascii_case(...)` or similar, checking length
first, and only fall through to `atof`-equivalent parsing for numeric epochs.

`tt` (swemplan.c:785) is computed **once, right after the epoch column**, using the file's raw
`tjd0` (before any later override — see mean-anomaly handling below) and is reused unchanged for
every other column's T-term evaluation on that row (`check_t_terms(tt, ...)` calls at lines 812,
831, 841, 851, 863, 875). Evaluation order matters: don't recompute `tt` per-column.

### Column 2 — equinox (swemplan.c:788–809)

Same `j2000`/`b1950`/`j1900` sentinels as epoch, **plus** a 4th: `"jdate"` → `*tequ = tjd` (the
date being queried — "equinox of date"). This is how Vulcan/White-Moon/Proserpina/Waldemath
declare their elements in the equinox-of-date frame (all four use `JDATE` in `seorbel.txt`).
Falls through to `atof(sp)` for a literal JD equinox otherwise.

### Columns 3–8 — mano/sema/ecce/parg/node/incl, all via `check_t_terms` (§3.1) (swemplan.c:810–884)

Each column: `retc = check_t_terms(tt, cpos[N], &out)`; angle columns (mano, parg, node, incl)
are `swe_degnorm`'d then `*= DEGTORAD`; `sema`/`ecce` get range checks (`sema <= 0` → error,
`ecce >= 1 || ecce < 0` → error, "no parabolic or hyperbolic orbits allowed"). `retc == ERR` from
`check_t_terms` (malformed expression) is always an error for that column, with a
column-specific message (`"mean anomaly value invalid"`, `"semi-axis value invalid"`, etc.).

**Special case — mean anomaly only** (swemplan.c:820–826):
```c
/* if mean anomaly has t terms (which happens with fictitious
 * planet Vulcan), we set
 * epoch = tjd, so that no motion will be added anymore
 * equinox = tjd */
if (retc == 1) {
  *tjd0 = tjd;
}
```
`check_t_terms` returns `1` (not `ERR`) when the expression contains any `+`/`-` (i.e. has
additional polynomial terms beyond a bare constant — see §3.1). When that's true for the **mean
anomaly** column specifically, `*tjd0` (the epoch, already used to compute `tt` above) is
overwritten to the query date `tjd`. Effect: back in `swi_osc_el_plan` (§4), the daily-motion term
`(tjd - tjd0) * dmot` becomes exactly zero, because the T-polynomial already encodes the fully
secularly-evolved mean anomaly at `tjd` — adding `dmot * (tjd - tjd0)` on top would double-count
the motion. **The code comment's claim that "equinox = tjd" is also set here is misleading — this
block only touches `tjd0` (epoch), never `tequ` (equinox).** The equinox is separately forced to
`tjd` via the *file's own* `JDATE` sentinel in column 2 (all four T-term bodies — Vulcan, White
Moon, Proserpina*, Waldemath — pair a `JDATE` equinox with polynomial columns in `seorbel.txt`;
*Proserpina's mean anomaly is actually a bare constant, `170.73`, with no T terms, so this override
never fires for Proserpina — only Vulcan/White-Moon/Waldemath have literal `+ ... * T` in column 3).
This override applies **only to the mano column**; `sema`/`ecce`/`parg`/`node`/`incl` T-terms (used
by Waldemath's `parg`/`node`, which do have polynomial expressions) are evaluated once against the
pre-override `tt` and used as static per-call values — there is no analogous "extra motion" term
for those elements in `swi_osc_el_plan` to suppress.

### Column 9 — name (swemplan.c:886–892)

Trimmed of leading whitespace and right-trimmed (`swi_right_trim`); copied verbatim including
embedded text like `"Leverrier (Neptune)"` or `"Selena/White Moon"`.

### Column 10 (optional) — `geo` flag (swemplan.c:893–899)

```c
if (fict_ifl != NULL && ncpos > 9) {
  for (sp = cpos[9]; *sp != '\0'; sp++) *sp = tolower(*sp);
  if (strstr(cpos[9], "geo") != NULL) *fict_ifl |= FICT_GEO;
}
```
Whole-field lower-cased in place, then substring-matched for `"geo"` anywhere in the field (so
e.g. `" geo"`, `"geocentric"`, `"geo "` would all match — it's a substring test, not equality).
In the real file, only **White Moon/Selena** (line 57) and **Waldemath** (line 68) carry `, geo`
as a trailing 10th field; all 15 built-in bodies and the file's Vulcan/Proserpina rows are
heliocentric-elements (`fict_ifl` stays 0). `FICT_GEO` changes both the Kepler daily-motion
constant and the Gaussian constant (§4), and which barycentric anchor (`xearth` vs `xsun`) the
final position is added to.

### §3.1 `check_t_terms` — polynomial-in-T expression parser (swemplan.c:916–967)

```c
tt[0] = t / 36525;   // Julian centuries since epoch — t = tjd - tjd0 (or tjd - tjd (=0) if overridden)
tt[1] = tt[0];        // T^1 (duplicate of tt[0], NOT T^0 — see hazard below)
tt[2] = tt[1] * tt[1]; // T^2
tt[3] = tt[2] * tt[1]; // T^3
tt[4] = tt[3] * tt[1]; // T^4
retc = strpbrk(sinp, "+-") != NULL ? 1 : 0;  // "has additional terms" flag, computed up front
```

State machine (single left-to-right pass, `fac` = running product for the *current* term,
`*doutp` = accumulated sum, `z` = term counter used only to suppress adding `fac` before any term
has been read):

```c
sp = sinp; *doutp = 0; fac = 1; z = 0;
while (1) {
  skip spaces/tabs;
  if (*sp is '+' or '-' or '\0') {
    if (z > 0) *doutp += fac;        // commit the term that was just accumulated in `fac`
    isgn = (*sp == '-') ? -1 : 1;
    fac = 1 * isgn;                   // reset accumulator, carrying the sign into the NEXT term
    if (*sp == '\0') return retc;     // done
    sp++;                             // consume the sign character
  } else {
    skip '*'/space/tab (implicit-multiplication separator);
    if next non-skipped char is 't'/'T' {
      sp++;
      if next char is '+'/'-'/end-of-term-implied: fac *= tt[0];     // bare "T" -> power 1
      else: i = atoi(sp); if (0 <= i <= 4) fac *= tt[i];             // "T2","T3","T4" -> that power (see hazard)
    } else {
      fac *= atof(sp);                 // a numeric literal (coefficient or bare constant)
    }
    skip digits/'.' of the number just consumed (only relevant in the numeric-literal branch);
  }
  z++;
}
```

Worked example — Vulcan's mean anomaly, `"252.8987988 + 707550.7341 * T"`:
1. Token `252.8987988` (number) → `fac = 252.8987988`.
2. Token `+` → commit: `doutp = 252.8987988`; `fac` reset to `1` (sign `+`).
3. Token `707550.7341` (number) → `fac = 707550.7341`.
4. Token `*` then `T` (bare, no trailing digit) → `fac *= tt[0]` → `fac = 707550.7341 * tt[0]`.
5. End of string → commit: `doutp = 252.8987988 + 707550.7341 * (t/36525)`. `return retc=1`
   (since a `+` was present).

Result: `mano(t) = 252.8987988° + 707550.7341° · T` where `T = (tjd − tjd0)/36525` — a standard
linear secular-rate expression in Julian centuries.

## §4. `swi_osc_el_plan` — elements → J2000 barycentric-equatorial state vector (swemplan.c:579–690)

Signature: `swi_osc_el_plan(double tjd, double *xp, int ipl, int ipli, double *xearth, double
*xsun, char *serr)`. `ipl` here is **already offset-adjusted** (caller passes `ipl - SE_FICT_OFFSET`,
sweph.c:1117/3479) — i.e. it's the 0-based row index into `plan_oscu_elem`/`seorbel.txt`, matching
`read_elements_file`'s `ipl` parameter directly. `xearth`/`xsun` are **already-computed
barycentric** Earth and Sun state vectors (6-vectors, position+velocity), supplied by the caller
(read from global `swed.pldat[SEI_EARTH].x` / `swed.pldat[SEI_SUNBARY].x` at the two call sites —
see Porting Notes on global state).

### 4.1 Fetch elements (swemplan.c:594–597)
`read_elements_file(ipl, tjd, &tjd0, &tequ, &mano, &sema, &ecce, &parg, &node, &incl, NULL,
&fict_ifl, serr)` — per §3. Failure propagates as `ERR` immediately.

### 4.2 Daily motion & Gaussian constant (swemplan.c:598–600, 647–650)
```c
dmot = 0.9856076686 * DEGTORAD / sema / sqrt(sema);   /* daily motion, deg/day -> rad/day */
if (fict_ifl & FICT_GEO)
  dmot /= sqrt(SUN_EARTH_MRAT);
...
K = (fict_ifl & FICT_GEO) ? KGAUSS_GEO / sqrt(sema) : KGAUSS / sqrt(sema);
```
Evaluation order: `dmot` is `((0.9856076686 * DEGTORAD) / sema) / sqrt(sema)` — two sequential
divisions, not a single `/ pow(sema, 1.5)` — preserve the exact operation order for FP fidelity.
`0.9856076686` deg/day is the mean daily motion of a 1-AU-period body (derived from the Gaussian
constant); dividing further by `sqrt(SUN_EARTH_MRAT)` rescales it for elements defined around
Earth's mass instead of the Sun's (only for `FICT_GEO` bodies — White Moon, Waldemath).

### 4.3 Gaussian P/Q/R rotation vectors (swemplan.c:601–616)

Standard Keplerian orbital-plane → ecliptic rotation, built from `cos`/`sin` of node/inclination/
argument-of-perihelion:
```c
cosnode = cos(node); sinnode = sin(node);
cosincl = cos(incl); sinincl = sin(incl);
cosparg = cos(parg); sinparg = sin(parg);
pqr[0] = cosparg*cosnode - sinparg*cosincl*sinnode;
pqr[1] = -sinparg*cosnode - cosparg*cosincl*sinnode;
pqr[2] = sinincl*sinnode;
pqr[3] = cosparg*sinnode + sinparg*cosincl*cosnode;
pqr[4] = -sinparg*sinnode + cosparg*cosincl*cosnode;
pqr[5] = -sinincl*cosnode;
pqr[6] = sinparg*sinincl;
pqr[7] = cosparg*sinincl;
pqr[8] = cosincl;
```
(`pqr[0..2]` = P vector components x/y/z... actually laid out as two rows of 3 used below as
`[Px,Py,Pz, Qx,Qy,Qz, Rx,Ry,Rz]` — note only `pqr[0],pqr[1],pqr[3],pqr[4],pqr[6],pqr[7]` are
actually *used* in the position/velocity transform at §4.5 below; `pqr[2]`, `pqr[5]`, `pqr[8]`
are computed but never read in this function — dead values, harmless but worth noting so a port
doesn't feel obligated to "use" them.)

### 4.4 Kepler equation (swemplan.c:617–645)

```c
E = M = swi_mod2PI(mano + (tjd - tjd0) * dmot);   /* mean anomaly of date, radians */
```
(For a T-term mano with the epoch override from §3, `tjd - tjd0 == 0` exactly, so `M = mano`
verbatim, already secularly evolved.)

**High-eccentricity initial-guess refinement** (only when `ecce > 0.975`, e.g. Nibiru's 0.981092):
```c
if (ecce > 0.975) {
  M2 = M * RADTODEG;
  if (M2 > 150 && M2 < 210) { M2 -= 180; M_180_or_0 = 180; } else M_180_or_0 = 0;
  if (M2 > 330) M2 -= 360;
  if (M2 < 0) { M2 = -M2; Msgn = -1; } else Msgn = 1;
  if (M2 < 30) {
    M2 *= DEGTORAD;
    alpha = (1 - ecce) / (4 * ecce + 0.5);
    beta = M2 / (8 * ecce + 1);
    zeta = pow(beta + sqrt(beta*beta + alpha*alpha), 1/3);   /* SEE HAZARD BELOW */
    sigma = zeta - alpha / 2;
    sigma = sigma - 0.078 * sigma*sigma*sigma*sigma*sigma / (1 + ecce);
    E = Msgn * (M2 + ecce * (3*sigma - 4*sigma*sigma*sigma)) + M_180_or_0;
  }
}
E = swi_kepler(E, M, ecce);
```
This is a cube-root-based analytic starting-value formula (valid only very near perihelion,
`|M2| < 30°` after folding into a canonical range) intended to give `swi_kepler` a good initial
guess for high-eccentricity orbits where the naive `E0 = M` guess converges slowly. **`1/3` is
integer division in C, evaluating to `0`** — so `pow(x, 1/3)` is actually `pow(x, 0) = 1.0` for
any `x`, i.e. `zeta` is unconditionally `1.0`, not the intended cube root. This is a real bug in
the C source, silently neutering the entire refinement block into "always start from `zeta=1`."
See Porting Notes for why this must be reproduced bit-for-bit rather than "fixed."

`swi_kepler` (swephlib.c:4065–4096) — Newton/fixed-point solver, convergence tolerance `1e-12`
radians, **no iteration cap**:
```c
double swi_kepler(double E, double M, double ecce) {
  double dE = 1, E0, x;
  if (ecce < 0.4) {
    while (dE > 1e-12) {
      E0 = E;
      E = M + ecce * sin(E0);          /* simple fixed-point iteration */
      dE = fabs(E - E0);
    }
  } else {
    while (dE > 1e-12) {
      E0 = E;
      x = (M + ecce*sin(E0) - E0) / (1 - ecce*cos(E0));   /* Newton step */
      dE = fabs(x);
      if (dE < 1e-2) {
        E = E0 + x;                     /* skip swi_mod2PI for small steps (gcc optimizer workaround, see comment) */
      } else {
        E = swi_mod2PI(E0 + x);
        dE = fabs(E - E0);
      }
    }
  }
  return E;
}
```
Branch is chosen by **eccentricity**, not by which of the 19 bodies it is — every fictitious body
with `ecce < 0.4` (most of them; several have `ecce == 0` exactly, converging in 1 iteration since
`dE = fabs(M + 0 - E0)`... for `ecce=0`, `E = M` immediately, `dE = |M - E0|` — first iteration
already exact if `E0` started at `M`, which it does per line 618, so **zero-eccentricity bodies
converge in exactly one loop pass**) uses the simple fixed-point form; `ecce >= 0.4` (only Nibiru
at 0.981092 among the 19 named bodies — none of the others exceed roughly 0.31) uses Newton's
method.

### 4.5 Position/velocity in orbital plane, then rotate to ecliptic (swemplan.c:646–665)

```c
cose = cos(E); sine = sin(E);
fac = sqrt((1-ecce)*(1+ecce));       /* = sqrt(1-ecce^2), written as a product not 1-e^2 directly */
rho = 1 - ecce*cose;
x[0] = sema*(cose - ecce);
x[1] = sema*fac*sine;
x[3] = -K*sine/rho;
x[4] = K*fac*cose/rho;
/* ecliptic (still equinox `tequ`) */
xp[0] = pqr[0]*x[0] + pqr[1]*x[1];
xp[1] = pqr[3]*x[0] + pqr[4]*x[1];
xp[2] = pqr[6]*x[0] + pqr[7]*x[1];
xp[3] = pqr[0]*x[3] + pqr[1]*x[4];
xp[4] = pqr[3]*x[3] + pqr[4]*x[4];
xp[5] = pqr[6]*x[3] + pqr[7]*x[4];
```
`x[2]`/`x[5]` (out-of-plane components) are never set — implicitly zero, consistent with pure
two-body Kepler motion confined to the orbital plane.

### 4.6 Equator + precession + barycentric shift (swemplan.c:666–689)

```c
eps = swi_epsiln(tequ, 0);                    /* mean obliquity at the elements' own equinox */
swi_coortrf(xp, xp, -eps);                    /* ecliptic -> equatorial (note: NEGATIVE eps) */
swi_coortrf(xp+3, xp+3, -eps);
if (tequ != J2000) {
  swi_precess(xp, tequ, 0, J_TO_J2000);       /* precess position: tequ -> J2000 */
  swi_precess(xp+3, tequ, 0, J_TO_J2000);     /* precess velocity likewise */
}
if (fict_ifl & FICT_GEO) {
  for (i = 0; i <= 5; i++) xp[i] += xearth[i];   /* geocentric elements -> add barycentric Earth */
} else {
  for (i = 0; i <= 5; i++) xp[i] += xsun[i];     /* heliocentric elements -> add barycentric Sun */
}
if (pdp->x == xp) {                            /* only true when xp IS the cached plan_data slot */
  pdp->teval = tjd;                            /* global-state bookkeeping, see Porting Notes */
  pdp->iephe = pedp->iephe;
}
return OK;
```
Note `swi_precess(xp+3, ...)` precesses the *velocity* using the **same plain `swi_precess`** used
for position (not a dedicated `precess_speed` variant) — because at this stage velocity is still a
simple Cartesian rotation of the position-frame, not yet subject to the light-time-dependent speed
corrections that come later in `app_pos_etc_plan_osc` (§5).

## §5. Dispatch & apparent-position pipeline (sweph.c)

### 5.1 Dispatch (sweph.c:1106–1136)

```c
} else if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX) {
  ipli = SEI_ANYBODY;
  pdp = &swed.pldat[ipli];
  xp = pdp->xreturn;
do_fict_plan:
  retc = main_planet(tjd, SEI_EARTH, 0, epheflag, iflag, serr);   /* populates swed.pldat[SEI_EARTH] and SEI_SUNBARY as a side effect */
  iflag = swed.pldat[SEI_EARTH].xflgs;
  if (swi_osc_el_plan(tjd, pdp->x, ipl-SE_FICT_OFFSET, ipli, pedp->x, psdp->x, serr) != OK)
    goto return_error;
  if (retc == ERR)
    goto return_error;
  retc = app_pos_etc_plan_osc(ipl, ipli, iflag, serr);
  if (retc == ERR)
    goto return_error;
  if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS) {
    if (epheflag != SEFLG_MOSEPH) {
      iflag = (iflag & ~SEFLG_EPHMASK) | SEFLG_MOSEPH;
      epheflag = SEFLG_MOSEPH;
      strcat(serr, "\nusing Moshier eph.; ");
      goto do_fict_plan;                        /* retry entirely with Moshier Earth/Sun */
    } else
      goto return_error;
  }
}
```
Every fictitious body routes through the **single shared internal body slot** `SEI_ANYBODY`
(`ipli = SEI_ANYBODY` unconditionally, sweph.c:1108) — this is the same slot asteroids use, so
there is no per-fictitious-body persistent cache entry; two different fictitious bodies computed
back-to-back will each overwrite the shared `SEI_ANYBODY` `plan_data` struct. `main_planet` for
Earth is always computed first, regardless of what `epheflag` (Moshier/Swiss/JPL) is configured,
because the fictitious-body Kepler math needs a real Earth+Sun barycentric position to shift into
(§4.6) even though the fictitious body's *own* motion is pure analytic Kepler, no ephemeris file
of its own. On light-time-range failure (`BEYOND_EPH_LIMITS`/`NOT_AVAILABLE` from the *Earth*
lookup, propagated through `app_pos_etc_plan_osc`'s re-evaluation at `t-dt` in §5.2), the whole
`do_fict_plan` block retries once with forced `SEFLG_MOSEPH` — same fallback pattern as the
asteroid path just above it (sweph.c:1086–1095).

Invalid-`ipl` fallthrough (sweph.c:1140–1144): any `ipl` not matching any branch (including
`ipl > SE_FICT_MAX`) → `"illegal planet number %d."`.

Name lookup for error messages / `swe_get_planet_name` (sweph.c:7041–7043):
```c
if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX) {
  swi_get_fict_name(ipl - SE_FICT_OFFSET, s);
}
```

### 5.2 `app_pos_etc_plan_osc` (sweph.c:3365–3547) — the apparent-position pipeline

This is a **separate, dedicated function** from the general `app_pos_etc_plan` (sweph.c:2465,
used for real planets/asteroids/Moon) — it is NOT a code-path variant of the general one; the two
have diverged in a few specific ways documented below. `app_pos_etc_plan_osc` is called **only**
from the fictitious-planet dispatch (§5.1); nothing else invokes it.

**Observer position** (sweph.c:3396–3422): identical structure to the general path — topocentric
adds `swi_get_observer` offset to barycentric Earth; barycentric request zeroes the observer;
heliocentric request uses barycentric Sun (`psdp->x`, zero under Moshier since Moshier has no
independent Sun-barycentric position); otherwise plain barycentric Earth (`pedp->x`).

**Light-time** (sweph.c:3426–3493) — skipped entirely under `SEFLG_TRUEPOS`. Otherwise:
```c
niter = 1;   /* HARDCODED — always 2 light-time passes (j=0,1), regardless of ephemeris backend */
```
This is the key divergence from `app_pos_etc_plan`, whose niter is backend-dependent
(`niter=1` for JPL/SWIEPH, `niter=0` for Moshier/osculating — sweph.c:2547–2551, though that
comment's mention of "or planet from osculating elements" is stale/misleading since osculating
bodies never actually reach `app_pos_etc_plan`, only `app_pos_etc_plan_osc`, which always uses
`niter=1` unconditionally). **A Rust port must give fictitious bodies exactly 2 light-time
iterations always**, not the backend-conditional count used elsewhere.

Speed pre-pass (`SEFLG_SPEED`, sweph.c:3428–3454): estimates `dt` from a "true minus rough-apparent"
position delta using `pdp->x[i+3]` (the body's own velocity) to extrapolate a rough apparent
position, iterated `niter+1` times (i.e. twice). Main pass (sweph.c:3456–3471): recomputes `dt`
from the actual light-time distance and produces `xx[i] = pdp->x[i] - dt*pdp->x[i+3]` — **a linear
extrapolation using the already-computed instantaneous velocity, not a re-evaluation of the Kepler
orbit at `t-dt`** for the *position* itself.

**Speed refinement — re-derive the whole Kepler solution at `t-dt`** (sweph.c:3472–3492, only
under `SEFLG_SPEED`):
```c
t = pdp->teval - dt;
retc = main_planet_bary(t, SEI_EARTH, epheflag, iflag, NO_SAVE, xearth, xearth, xsun, xmoon, serr);
if (swi_osc_el_plan(t, xx, ipl-SE_FICT_OFFSET, ipli, xearth, xsun, serr) != OK)
  return ERR;
if (retc != OK) return retc;
```
Unlike the position computation (which uses a linear velocity extrapolation), speed accuracy
requires **re-running the full analytic Kepler pipeline** (§4) at `t = teval - dt`, with a freshly
computed Earth/Sun barycentric state at that earlier epoch (`main_planet_bary`, `NO_SAVE` — doesn't
touch the cached `swed.pldat[SEI_EARTH]`). The resulting `xx` (from this second `swi_osc_el_plan`
call, at `t-dt`) plus the topocentric offset at `t-dt` (`xobs2`) feeds into the "part of daily
motion resulting from change of dt" correction later applied to `xx[3..5]` (sweph.c:3505–3507) —
same aberration-of-motion technique documented in the crossings/other planet ref docs' shared infra.

**Geocentric conversion, deflection, aberration, precession, and the `app_pos_rest` tail**
(sweph.c:3494–3546) are **byte-for-byte identical in structure** to the general planet apparent-
position pipeline already covered by this codebase's existing `apparent_planet`
(`calc.rs`/`context.rs`, per `docs/codebase-map.md`): `swi_deflect_light`, `swi_aberr_light`,
`swi_precess`/`swi_precess_speed`, then `app_pos_rest` (nutation, ecliptic transform, sidereal,
polar conversion, degrees). **Reuse that existing pipeline machinery rather than re-deriving it —
the only genuinely new logic for fictitious bodies is §4 (`swi_osc_el_plan`) and the file-parsing
of §3; the apparent-position tail is infrastructure this codebase already has.**

### Error strings surfaced by the fictitious-body path

- `"error no elements for fictitious body no %7.0f"` — file absent, `ipl >= SE_NFICT_ELEM` (swemplan.c:711).
- `"error in file %s, line %7.0f: nine elements required"` — malformed row, any row scanned before/at the match (swemplan.c:757–761).
- `"... invalid epoch"` / `"... invalid equinox"` — unrecognized epoch/equinox sentinel string starting with `j`/`b` that isn't `j1900`/`j2000`/`b1950`(/`jdate` for equinox) (swemplan.c:778–782, 802–806).
- `"... mean anomaly value invalid"` / `"... semi-axis value invalid"` / `"... eccentricity invalid (no parabolic or hyperbolic orbits allowed)"` / `"... perihelion argument value invalid"` / `"... node value invalid"` / `"... inclination value invalid"` — per-column `check_t_terms` failure or (for sema/ecce) out-of-range value (swemplan.c:814–883).
- `"... elements for planet %7.0f not found"` — file opened successfully but target row never reached (swemplan.c:902–906).
- `"illegal planet number %d."` — `ipl` outside every recognized range, including `> SE_FICT_MAX` (sweph.c:1142).

### Boundary conditions

No date-range restriction analogous to `CHIRON_START`/`CHIRON_END` or `MOSHLUEPH_START`/`_END`
exists anywhere in this path — Kepler propagation from `tjd0` is unbounded in time (it will happily
extrapolate centuries or millennia away; the resulting position becomes astronomically meaningless
far from the fitted epoch, but the C code enforces no cutoff and neither should the Rust port).
**Vulcan has no special validity-range check in the C source** — despite being a hypothetical
intra-Mercurial body, it is propagated with the same unconditional analytic Kepler formula as
every other fictitious body; its only "special" treatment is the T-term epoch-override mechanic
in §3 (shared with White Moon and Waldemath), not a date guard.

## Porting notes

- **Global-state reads/writes a stateless port must eliminate**:
  - `swi_osc_el_plan` reads `swed.pldat[SEI_EARTH]` and `swed.pldat[SEI_SUNBARY]` *indirectly*
    only insofar as its caller passes `pedp->x`/`psdp->x` as the `xearth`/`xsun` parameters — the
    function itself is otherwise pure given those two 6-vectors plus `tjd`/`ipl`. Port the
    signature as-is: `fn osc_el_plan(tjd, body_row, xearth: [f64;6], xsun: [f64;6]) -> Result<[f64;6], Error>`,
    with the caller (the fictitious-body branch of `calc.rs`) responsible for first computing
    Earth's barycentric state via the existing Moshier/Swiss/JPL Earth provider (same one
    `apparent_planet` already uses) and passing it in explicitly.
  - `swi_osc_el_plan`'s tail (`if (pdp->x == xp) { pdp->teval = tjd; pdp->iephe = pedp->iephe; }`,
    swemplan.c:685–688) is pure C-side caching bookkeeping — irrelevant to a stateless port, do not
    replicate.
  - `read_elements_file` reopens+closes `seorbel.txt` on every call (no caching, swemplan.c:707,
    "-1, because file information is not saved, file is always closed"). A Rust port is free to
    (and should) cache the parsed catalog once per `Ephemeris` construction or lazily-and-cached,
    matching the existing `stars.rs`/`load_catalog` pattern (`src/stars.rs:298`) — this is a
    behavior-preserving optimization since re-parsing produces identical results every time (the
    file's own mtime isn't checked either way).
  - `main_planet`/`main_planet_bary` for Earth mutate `swed.pldat[SEI_EARTH]`/`[SEI_SUNBARY]` as a
    side effect in C (`DO_SAVE` vs `NO_SAVE` distinguishes the two call sites, sweph.c:1113 vs
    sweph.c:3478) — in the stateless Rust design this is simply "call the Earth-position function
    and use its return value," no distinct save/no-save modes needed.
  - All fictitious bodies share the single `SEI_ANYBODY` internal slot in C (no per-body cache) —
    not relevant to a stateless port (there is no cache to alias in the first place), but explains
    why the C code has no special-casing for "is this the same fictitious body as last time."

- **The `pow(x, 1/3)` integer-division bug (swemplan.c:638) — reproduce, do not fix.** The
  high-eccentricity refinement branch (only reachable for `ecce > 0.975`, i.e. only Nibiru among
  the 19 named bodies, and only when the mean anomaly of date falls within the narrow
  post-folding window `|M2| < 30°` near perihelion) computes `zeta = pow(beta + sqrt(beta*beta +
  alpha*alpha), 1/3)`. Because `1` and `3` are both `int` in C, `1/3` evaluates to integer `0` at
  compile time, making this call `pow(x, 0.0) == 1.0` for any `x > 0` — the cube root is never
  actually taken. The subsequent `sigma`/`E` formula is then evaluated using this constant
  `zeta=1.0`, producing a specific (not mathematically "correct" cube-root-based) initial guess
  for `E` that is then handed to `swi_kepler` for full Newton refinement to `1e-12` rad. Because
  `swi_kepler`'s convergence criterion is an absolute tolerance rather than a fixed iteration
  count, a *mathematically correct* cube-root starting guess would very likely converge to the
  same final `E` to well within double precision — but "very likely" is not "certain" for a
  bitwise-exact golden test, and the C behavior is unambiguous and cheap to replicate exactly.
  **Port `pow(x, 0.0)` (i.e. just use the literal `1.0`) rather than `x.cbrt()` or `x.powf(1.0/3.0)`**,
  and comment why, citing this doc and swemplan.c:638.
- **`check_t_terms`'s `"T0"` quirk**: `tt[0]` and `tt[1]` are both `T^1` (only `tt[2..4]` are true
  powers `T^2..T^4`); a bare `atoi(sp)` result of `0` after a `T` (e.g. a hypothetical `"T0"` token)
  indexes `tt[0]`, which is `T^1`, **not** `T^0 = 1`. No entry in the shipped `seorbel.txt` actually
  writes `T0`, so this is latent, but a Rust parser should model the table as `[T^1, T^1, T^2, T^3,
  T^4]` (or equivalently special-case index 0 and 1 to the same value) to stay bug-compatible if a
  user-supplied custom `seorbel.txt` ever used it.
- **`check_t_terms` output-accumulation timing**: the running product `fac` for a term is only
  folded into `*doutp` when the parser *next* encounters a `+`, `-`, or end-of-string — i.e. the
  final (or only) term's contribution is added on loop exit, not incrementally. A straightforward
  Rust re-implementation (e.g. `sinp.split_inclusive(['+','-'])`-based tokenization with explicit
  per-term evaluation) is fine and arguably clearer, but must reduce to the exact same worked
  trace as the §3.1 example for every row in `seorbel.txt` — validate against all 19 rows (plus any
  additional custom rows a golden test wants to exercise), not just the simple non-T rows.
- **Row-index-based lookup, not ID-based**: `read_elements_file`'s `ipl` parameter is matched
  against a **sequential content-line counter** (`iplan`, incremented once per successfully parsed
  line, §3), not against any value written in the file. A Rust catalog loader should build a
  `Vec<FictElements>` indexed 0-based by line-parse-order and look up by `ipl` (already
  `ipl - SE_FICT_OFFSET`-adjusted by the caller) as a plain vector index, mirroring the C
  behavior — the `"# N"` trailing comments in the shipped file are documentation only.
- **`FICT_GEO` (the `geo` 10th-field flag) changes three things together**, not just the additive
  anchor: (1) `dmot` divided by `sqrt(SUN_EARTH_MRAT)`, (2) `K` uses `KGAUSS_GEO` instead of
  `KGAUSS`, (3) the final barycentric shift adds `xearth` instead of `xsun`. All three must flip
  together on the same flag — there is no way to get e.g. the geocentric daily-motion scaling with
  a heliocentric anchor in the C model.
- **Which built-in rows to embed**: only transcribe the **Neely-revised** 8-row Uranian block
  (§1, "LIVE") plus the 7 non-Neely-variant rows (Isis-Transpluto..Pickering) as the Rust built-in
  constant table — `SE_NEELY` is unconditionally defined in every real build, making the `#else`
  block dead code with no golden-test signal to validate against. Do not port it.
- **Built-in table vs. shipped `seorbel.txt` numeric mismatch (Kronos semi-axis, §1)**: the
  built-in fallback table and the ephemeris-directory file disagree by 0.0027 AU for Kronos. Since
  `ephe/seorbel.txt` ships with every real deployment, golden tests run against a real ephemeris
  path will exercise the *file* value, not the *built-in* value — make sure a built-in-table-only
  golden test (simulating "no seorbel.txt present") is exercised separately if bit-fidelity of the
  fallback table itself matters, since the two code paths are not cross-validated against each
  other by any existing C test.
- **`niter=1` hardcoded for the light-time loop is a genuine, deliberate divergence** from the
  ephemeris-backend-conditional `niter` used by the general planet/asteroid path — do not
  parameterize this by backend for fictitious bodies; it is always exactly 2 passes (§5.2).
- **Reuse existing infrastructure, do not re-derive**: `swi_deflect_light`, `swi_aberr_light`,
  `swi_precess`/`swi_precess_speed`, and the `app_pos_rest` tail (nutation, ecliptic transform,
  sidereal, polar/degree conversion) are structurally identical to what `calc.rs`'s existing
  `apparent_planet` pipeline already implements (per `docs/codebase-map.md`'s description of
  "light-time, retarded velocity, aberration, deflection pipeline"). The fictitious-body port's
  genuinely new surface area is: (a) the seorbel.txt parser + built-in table (§3, sibling to
  `stars.rs`'s `parse_catalog`), (b) `swi_kepler` (§4.4, a new Kepler-equation solver — grep
  confirms no existing equivalent in `src/`), and (c) `swi_osc_el_plan` itself (§4) as the
  position/velocity source feeding into the existing apparent-position machinery, structured
  analogously to how `orbit.rs` already separates "derive elements from state" from "the shared
  calc/flags/error plumbing." `swi_osc_el_plan`'s Gaussian-vector rotation (§4.3) and Kepler
  position/velocity construction (§4.5) are net-new math distinct from `orbit.rs`'s (inverse
  direction: state→elements, not elements→state) — do not attempt to force a shared abstraction
  between them; per this repo's `CLAUDE.md` constraint, diverging inputs (forward vs. inverse
  Kepler problems) justify separate implementations.
- **Bodies 55–58 (Vulcan, White Moon, Proserpina, Waldemath) require `seorbel.txt`** — there is no
  built-in-table fallback for them (`SE_NFICT_ELEM == 15`). A Rust `Ephemeris` configured without
  an ephemeris path containing `seorbel.txt` should surface `Error::FileNotFound` (or the existing
  fictitious-specific analogue) for these four IDs specifically, while bodies 40–54 should still
  resolve via the embedded built-in table.
