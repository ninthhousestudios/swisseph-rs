/*
 * Generates golden reference data for asteroid calc pipeline.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_asteroid gen_asteroid.c \
 *      ../../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_asteroid > tests/golden-data/asteroid.json
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "swephexp.h"

static int bodies[] = {
    SE_CHIRON, SE_PHOLUS, SE_CERES, SE_PALLAS, SE_JUNO, SE_VESTA
};
#define NBODIES 6

static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2460310.5,   /* 2024-Jan-1 */
    2433282.5,   /* 1950-Jan-1 */
    2415020.5,   /* 1900-Jan-1 */
    2305447.5,   /* 1600-Jan-1 */
};
#define NEPOCHS 5

struct flag_combo {
    int flag;
    const char *name;
};

/* NOTE: no MOSEPH cases here -- MOSEPH + asteroid output depends on process
 * call history (see docs/c-ref notes / commit 2c2bbcf), so it is unsafe to
 * generate deterministic stateless golden data for it. */
static struct flag_combo flag_combos[] = {
    { SEFLG_SWIEPH,                                           "SWIEPH" },
    { SEFLG_SWIEPH | SEFLG_SPEED,                              "SWIEPH_SPEED" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL,           "SWIEPH_SPEED_EQUATORIAL" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_XYZ,                  "SWIEPH_SPEED_XYZ" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000 | SEFLG_NONUT,  "SWIEPH_SPEED_J2000_NONUT" },
    { SEFLG_SWIEPH | SEFLG_TRUEPOS,                            "SWIEPH_TRUEPOS" },
    { SEFLG_SWIEPH | SEFLG_NONUT,                              "SWIEPH_NONUT" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_TOPOCTR,              "SWIEPH_SPEED_TOPOCTR" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR,               "SWIEPH_SPEED_HELCTR" },
    { SEFLG_JPLEPH | SEFLG_SPEED,                              "JPLEPH_SPEED" },
};
#define NFLAGS 10

int main(void) {
    char serr[256];
    char pname[AS_MAXCH];
    int first = 1;

    swe_set_ephe_path("ephe");
    swe_set_topo(8.55, 47.37, 500.0);

    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = flag_combos[ifl].flag;
                double xx[6];
                memset(xx, 0, sizeof(xx));

                /* Reset C library state before each call so file caching does
                 * not carry over between test cases. This gives deterministic,
                 * stateless golden data that matches our stateless Rust
                 * implementation. */
                swe_close();
                swe_set_ephe_path("ephe");
                swe_set_topo(8.55, 47.37, 500.0);

                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr,
                            "swe_calc error: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flags);
                    return 1;
                }

                swe_get_planet_name(bodies[ib], pname);

                if (!first) printf(",\n");
                first = 0;
                printf("  {\"body\": %d, \"body_name\": \"%s\", "
                       "\"jd\": %.17g, \"flags\": %d, \"flag_name\": \"%s\", "
                       "\"retflag\": %d, "
                       "\"output\": [%.17g, %.17g, %.17g, %.17g, %.17g, %.17g]}",
                       bodies[ib], pname,
                       epochs[ie], flags, flag_combos[ifl].name,
                       rc,
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
            }
        }
    }
    printf("\n]\n");
    swe_close();
    return 0;
}
