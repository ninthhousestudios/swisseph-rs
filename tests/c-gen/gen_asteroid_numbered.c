/*
 * Generates golden reference data for numbered asteroid calc pipeline.
 *
 * Tests Body::Asteroid(n) with MPC numbers: 433 (Eros), 7066 (Nessus),
 * 136199 (Eris, >99999 s%06d naming), 2060 (Chiron-as-numbered,
 * SEI_FILE_ANY_AST path distinct from SE_CHIRON via seas).
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_asteroid_numbered gen_asteroid_numbered.c \
 *      ../../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_asteroid_numbered > tests/golden-data/asteroid_numbered.json
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "swephexp.h"

static int bodies[] = {
    SE_AST_OFFSET + 433,     /* Eros */
    SE_AST_OFFSET + 7066,    /* Nessus */
    SE_AST_OFFSET + 136199,  /* Eris (>99999, exercises s%06d naming) */
    SE_AST_OFFSET + 2060,    /* Chiron-as-numbered (SEI_FILE_ANY_AST, not SE_CHIRON) */
};
#define NBODIES 4

static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2460310.5,   /* 2024-Jan-1 */
    2433282.5,   /* 1950-Jan-1 */
};
#define NEPOCHS 3

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_SWIEPH | SEFLG_SPEED,                              "SWIEPH_SPEED" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL,           "SWIEPH_SPEED_EQUATORIAL" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_XYZ,                  "SWIEPH_SPEED_XYZ" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000 | SEFLG_NONUT,  "SWIEPH_SPEED_J2000_NONUT" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR,               "SWIEPH_SPEED_HELCTR" },
    { SEFLG_JPLEPH | SEFLG_SPEED,                              "JPLEPH_SPEED" },
};
#define NFLAGS 6

int main(void) {
    char serr[256];
    int first = 1;

    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = flag_combos[ifl].flag;
                double xx[6];
                memset(xx, 0, sizeof(xx));

                swe_close();
                swe_set_ephe_path("../../ephe");

                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr,
                            "swe_calc error: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flags);
                    return 1;
                }

                if (!first) printf(",\n");
                first = 0;
                printf("  {\"body\": %d, "
                       "\"jd\": %.17g, \"flags\": %d, \"flag_name\": \"%s\", "
                       "\"retflag\": %d, "
                       "\"output\": [%.17g, %.17g, %.17g, %.17g, %.17g, %.17g]}",
                       bodies[ib],
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
