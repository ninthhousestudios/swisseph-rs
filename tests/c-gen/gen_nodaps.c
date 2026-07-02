/*
 * Generates golden reference data for the mean nodes & apsides
 * (swe_nod_aps with method = SE_NODBIT_MEAN), swisseph-rs/85 (PNOC 4).
 *
 * Bodies: SE_SUN, SE_MOON, SE_MERCURY..SE_NEPTUNE, SE_EARTH (10 bodies — the
 *         mean-eligible set: Sun..Neptune numeric range + Earth).
 * Method: SE_NODBIT_MEAN.
 * Flags:  MOSEPH|SPEED, MOSEPH|SPEED|EQUATORIAL, MOSEPH,
 *         MOSEPH|SPEED|TRUEPOS, MOSEPH|SPEED|EQUATORIAL|TRUEPOS.
 *         The two TRUEPOS combos give the pure geometry (no light deflection /
 *         aberration) — bit-exact vs the Rust port, so the golden test asserts
 *         them tightly; the apparent (light-effect) combos are asserted loosely
 *         for the descending node only (see the test's tolerance note).
 * Epochs: 4 of the gen_calc.c epochs incl. one pre-1900 (1800-Jan-1).
 * => 10 bodies * 4 epochs * 5 flags = 200 cases, each recording all four
 *    node/apsis state vectors (asc, desc, peri, aphe) [0..6) + retflag.
 *
 * Emitted under the top-level key "mean"; PNOC 5 will add "oscu"/"oscu_bar"/
 * "fopoint" keys to the same file.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_nodaps tests/c-gen/gen_nodaps.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_nodaps > tests/golden-data/nodaps.json
 */

#include <stdio.h>
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

int main(void) {
    char serr[256];
    swe_set_ephe_path("../swisseph/ephe");
    printf("{\n\"mean\": [\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = flag_combos[ifl].flag;
                double asc[6], desc[6], peri[6], aphe[6];
                memset(asc, 0, sizeof(asc));
                memset(desc, 0, sizeof(desc));
                memset(peri, 0, sizeof(peri));
                memset(aphe, 0, sizeof(aphe));
                int rc = swe_nod_aps(epochs[ie], bodies[ib], flags, SE_NODBIT_MEAN,
                                     asc, desc, peri, aphe, serr);
                if (rc < 0) {
                    fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flags);
                    return 1;
                }
                emit(bodies[ib], body_names[ib], epochs[ie], flags,
                     flag_combos[ifl].name, rc, asc, desc, peri, aphe);
            }
        }
    }
    printf("\n]\n}\n");
    return 0;
}
