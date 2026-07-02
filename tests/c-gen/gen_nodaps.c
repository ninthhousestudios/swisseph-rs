/*
 * Generates golden reference data for swe_nod_aps: mean (swisseph-rs/85,
 * PNOC 4) and osculating (swisseph-rs/86, PNOC 5) branches.
 *
 * "mean" (SE_NODBIT_MEAN):
 *   Bodies: SE_SUN, SE_MOON, SE_MERCURY..SE_NEPTUNE, SE_EARTH (10 bodies —
 *           the mean-eligible set: Sun..Neptune numeric range + Earth).
 *   Flags:  MOSEPH|SPEED, MOSEPH|SPEED|EQUATORIAL, MOSEPH,
 *           MOSEPH|SPEED|TRUEPOS, MOSEPH|SPEED|EQUATORIAL|TRUEPOS.
 *           The two TRUEPOS combos give the pure geometry (no light
 *           deflection/aberration) — bit-exact vs the Rust port, so the
 *           golden test asserts them tightly; the apparent (light-effect)
 *           combos are asserted loosely for the descending node only (see
 *           the test's tolerance note).
 *   Epochs: 4 of the gen_calc.c epochs incl. one pre-1900 (1800-Jan-1).
 *   => 10 bodies * 4 epochs * 5 flags = 200 cases.
 *
 * "oscu" (SE_NODBIT_OSCU):
 *   Bodies: SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN,
 *           SE_URANUS, SE_NEPTUNE, SE_PLUTO (9 bodies).
 *   Flags:  MOSEPH|SPEED, SWIEPH|SPEED (2 — exercises both `nodaps_observer`/
 *           `nodaps_osc_body_j2000` backend arms).
 *   Epochs: all 4 (same `epochs[]` as "mean").
 *   => 9 bodies * 4 epochs * 2 flags = 72 cases.
 *
 * "oscu_bar" (SE_NODBIT_OSCU_BAR): only meaningful for a backend with a real
 *   solar-system barycenter distinct from the Sun (SWIEPH) — Moshier has no
 *   such frame (Rust rejects OSCU_BAR there, matching calc_inner's BARYCTR
 *   gate), so this battery is SWIEPH-only.
 *   Bodies: SE_JUPITER, SE_SATURN, SE_PLUTO (beyond the 6 AU threshold),
 *           SE_MERCURY (inside it) — covers both sides of the threshold.
 *   Flags:  SWIEPH|SPEED.
 *   Epochs: epochs[0..1] (J2000, 2024-Jan-1).
 *   => 4 bodies * 2 epochs * 1 flag = 8 cases.
 *
 * "fopoint" (SE_NODBIT_OSCU | SE_NODBIT_FOPOINT):
 *   Bodies: SE_MOON, SE_MARS, SE_JUPITER.
 *   Flags:  MOSEPH|SPEED.
 *   Epochs: epochs[0..1] (J2000, 2024-Jan-1).
 *   => 3 bodies * 2 epochs * 1 flag = 6 cases.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_nodaps tests/c-gen/gen_nodaps.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_nodaps > tests/golden-data/nodaps.json
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swephexp.h"

static int bodies[] = {SE_SUN,    SE_MOON,   SE_MERCURY, SE_VENUS,  SE_MARS,
                       SE_JUPITER, SE_SATURN, SE_URANUS,  SE_NEPTUNE, SE_EARTH};
static const char *body_names[] = {"Sun",     "Moon",   "Mercury", "Venus",   "Mars",
                                   "Jupiter", "Saturn", "Uranus",  "Neptune", "Earth"};
#define NBODIES 10

static double epochs[] = {
    2451545.0, /* J2000.0 */
    2460310.5, /* 2024-Jan-1 */
    2433282.5, /* 1950-Jan-1 */
    2378496.5, /* 1800-Jan-1 (pre-1900) */
};
#define NEPOCHS 4

/* Same as `epochs[]` but with the pre-1900 epoch nudged off the sepl_18 .se1
 * file's `tfstart` boundary (2378496.5) — the "oscu" battery's 3-position
 * SWIEPH sampling (t-dt/t/t+dt) hits the documented stateful-vs-stateless
 * file-selection divergence at an exact boundary (CLAUDE.md
 * <stateless_tolerance> §2; same fix as tests/golden/pheno.rs). "mean" keeps
 * the un-nudged epoch since its closed-form polynomials never touch file
 * selection. */
static double osc_epochs[] = {
    2451545.0, /* J2000.0 */
    2460310.5, /* 2024-Jan-1 */
    2433282.5, /* 1950-Jan-1 */
    2378500.5, /* 1800-Jan-5 (pre-1900, off sepl_18's tfstart) */
};

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    {SEFLG_MOSEPH | SEFLG_SPEED, "SPEED"},
    {SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL, "EQUATORIAL"},
    {SEFLG_MOSEPH, "no_speed"},
    {SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_TRUEPOS, "TRUEPOS"},
    {SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_TRUEPOS, "TRUEPOS_EQU"},
};
#define NFLAGS 5

/* "oscu" battery — Moon + Mercury..Neptune (SE_NODBIT_OSCU forces the
 * osculating branch even though these are also mean-eligible) plus Pluto
 * (always osculating, regardless of method). Skips Sun/Earth, which have no
 * two-body ellipse of their own in this API. */
static int osc_bodies[] = {SE_MOON,   SE_MERCURY, SE_VENUS,  SE_MARS,   SE_JUPITER,
                           SE_SATURN, SE_URANUS,  SE_NEPTUNE, SE_PLUTO};
static const char *osc_body_names[] = {"Moon",    "Mercury", "Venus",  "Mars",   "Jupiter",
                                       "Saturn",  "Uranus",  "Neptune", "Pluto"};
#define NOSCBODIES 9

static struct flag_combo osc_flag_combos[] = {
    {SEFLG_MOSEPH | SEFLG_SPEED, "MOSEPH_SPEED"},
    {SEFLG_SWIEPH | SEFLG_SPEED, "SWIEPH_SPEED"},
};
#define NOSCFLAGS 2

/* "oscu_bar" battery — Jupiter/Saturn/Pluto (beyond the 6 AU threshold) and
 * Mercury (inside it), SWIEPH only (Moshier has no real barycenter). */
static int bar_bodies[] = {SE_JUPITER, SE_SATURN, SE_PLUTO, SE_MERCURY};
static const char *bar_body_names[] = {"Jupiter", "Saturn", "Pluto", "Mercury"};
#define NBARBODIES 4

/* "fopoint" battery — 2nd focal point instead of aphelion. */
static int fop_bodies[] = {SE_MOON, SE_MARS, SE_JUPITER};
static const char *fop_body_names[] = {"Moon", "Mars", "Jupiter"};
#define NFOPBODIES 3

#define NOSCEPOCHS 2 /* epochs[0..1]: J2000, 2024-Jan-1 */

static int first = 1;

static void emit(int body, const char *body_name, double jd, int flags,
                 const char *flag_name, int rc, const double *asc,
                 const double *desc, const double *peri, const double *aphe) {
    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, \"body_name\": \"%s\", "
           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", \"retflag\": %d, "
           "\"asc\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
           "\"desc\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
           "\"peri\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
           "\"aphe\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
           body, body_name, jd, flags, flag_name, rc,
           asc[0], asc[1], asc[2], asc[3], asc[4], asc[5],
           desc[0], desc[1], desc[2], desc[3], desc[4], desc[5],
           peri[0], peri[1], peri[2], peri[3], peri[4], peri[5],
           aphe[0], aphe[1], aphe[2], aphe[3], aphe[4], aphe[5]);
}

/* Shared swe_nod_aps call + emit, reused by every battery below. */
static void run_case(double jd, int body, const char *body_name, int flags,
                     const char *flag_name, int method) {
    char serr[256];
    double asc[6], desc[6], peri[6], aphe[6];
    memset(asc, 0, sizeof(asc));
    memset(desc, 0, sizeof(desc));
    memset(peri, 0, sizeof(peri));
    memset(aphe, 0, sizeof(aphe));
    int rc = swe_nod_aps(jd, body, flags, method, asc, desc, peri, aphe, serr);
    if (rc < 0) {
        fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x method=%d\n",
                serr, body, jd, flags, method);
        exit(1);
    }
    emit(body, body_name, jd, flags, flag_name, rc, asc, desc, peri, aphe);
}

int main(void) {
    swe_set_ephe_path("../swisseph/ephe");
    printf("{\n\"mean\": [\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                run_case(epochs[ie], bodies[ib], body_names[ib],
                         flag_combos[ifl].flag, flag_combos[ifl].name,
                         SE_NODBIT_MEAN);
            }
        }
    }
    printf("\n],\n\"oscu\": [\n");
    first = 1;
    for (int ib = 0; ib < NOSCBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NOSCFLAGS; ifl++) {
                run_case(osc_epochs[ie], osc_bodies[ib], osc_body_names[ib],
                         osc_flag_combos[ifl].flag, osc_flag_combos[ifl].name,
                         SE_NODBIT_OSCU);
            }
        }
    }
    printf("\n],\n\"oscu_bar\": [\n");
    first = 1;
    for (int ib = 0; ib < NBARBODIES; ib++) {
        for (int ie = 0; ie < NOSCEPOCHS; ie++) {
            run_case(epochs[ie], bar_bodies[ib], bar_body_names[ib],
                     SEFLG_SWIEPH | SEFLG_SPEED, "SWIEPH_SPEED",
                     SE_NODBIT_OSCU_BAR);
        }
    }
    printf("\n],\n\"fopoint\": [\n");
    first = 1;
    for (int ib = 0; ib < NFOPBODIES; ib++) {
        for (int ie = 0; ie < NOSCEPOCHS; ie++) {
            run_case(epochs[ie], fop_bodies[ib], fop_body_names[ib],
                     SEFLG_MOSEPH | SEFLG_SPEED, "MOSEPH_SPEED",
                     SE_NODBIT_OSCU | SE_NODBIT_FOPOINT);
        }
    }
    printf("\n]\n}\n");
    return 0;
}
