/*
 * Generates golden reference data for MOSEPH-mode asteroid calc.
 *
 * CRITICAL: This must be a SEPARATE binary that issues ONLY SEFLG_MOSEPH calls.
 * C's MOSEPH asteroid output depends on the process-global swed.pldat[SEI_SUNBARY],
 * which stays zero only while the process has NEVER executed a SWIEPH/JPLEPH calc.
 * Do NOT merge this generator into any other generator binary — doing so would
 * contaminate the global state and produce non-reproducible golden data.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_asteroid_moseph gen_asteroid_moseph.c \
 *      ../../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_asteroid_moseph > tests/golden-data/asteroid_moseph.json
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "swephexp.h"

static int bodies[] = {
    SE_CERES,               /* main asteroid via seas file */
    SE_VESTA,               /* main asteroid via seas file */
    SE_AST_OFFSET + 433,    /* numbered asteroid (Eros) via ast0/ file */
};
#define NBODIES 3

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
    { SEFLG_MOSEPH,                        "MOSEPH" },
    { SEFLG_MOSEPH | SEFLG_SPEED,          "MOSEPH_SPEED" },
    { SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_XYZ, "MOSEPH_SPEED_XYZ" },
};
#define NFLAGS 3

int main(void) {
    char serr[256];
    int first = 1;

    swe_set_ephe_path("../../ephe");

    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = flag_combos[ifl].flag;
                double xx[6];
                memset(xx, 0, sizeof(xx));

                /* Do NOT call swe_close() between cases -- we need the
                 * process-global state to remain in its virgin MOSEPH-only
                 * condition throughout. The only state that matters is
                 * swed.pldat[SEI_SUNBARY] staying zero. */

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
