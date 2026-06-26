/*
 * Generates golden reference data for the full calc pipeline (swe_calc).
 * Tests multiple bodies, epochs, and flag combinations through the
 * Moshier backend.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_calc gen_calc.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_calc > ../golden-data/calc.json
 */

#include <stdio.h>
#include "swephexp.h"

static int bodies[] = {
    SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN
};
static const char *body_names[] = {
    "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"
};
#define NBODIES 7

static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2460310.5,    /* 2024-Jan-1 */
    2433282.5,    /* 1950-Jan-1 */
    2488069.5,    /* 2100-Jan-1 */
    2378496.5,    /* 1800-Jan-1 */
};
#define NEPOCHS 5

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { 0,             "default" },
    { SEFLG_TRUEPOS, "TRUEPOS" },
    { SEFLG_NOABERR, "NOABERR" },
    { SEFLG_NOGDEFL, "NOGDEFL" },
    { SEFLG_J2000,   "J2000" },
    { SEFLG_NONUT,   "NONUT" },
    { SEFLG_EQUATORIAL, "EQUATORIAL" },
    { SEFLG_XYZ,     "XYZ" },
    { SEFLG_RADIANS, "RADIANS" },
    { SEFLG_ICRS,    "ICRS" },
};
#define NFLAGS 10

int main(void) {
    char serr[256];
    swe_set_ephe_path(NULL);
    int first = 1;
    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = SEFLG_MOSEPH | SEFLG_SPEED | flag_combos[ifl].flag;
                double xx[6];
                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "error: %s body=%d jd=%.1f flags=%d\n",
                            serr, bodies[ib], epochs[ie], flags);
                    return 1;
                }
                if (!first) printf(",\n");
                first = 0;
                printf("  {\"body\": %d, \"body_name\": \"%s\", "
                       "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
                       "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                       bodies[ib], body_names[ib],
                       epochs[ie], flags, flag_combos[ifl].name,
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
            }
        }
    }
    printf("\n]\n");
    return 0;
}
