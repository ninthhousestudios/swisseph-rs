/*
 * Generates golden reference data for the SwissEph (.se1) calc pipeline.
 * Same structure as gen_calc.c but uses SEFLG_SWIEPH instead of SEFLG_MOSEPH.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_calc_sweph gen_calc_sweph.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_calc_sweph > ../golden-data/calc_sweph.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static int bodies[] = {
    SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER,
    SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO
};
static const char *body_names[] = {
    "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter",
    "Saturn", "Uranus", "Neptune", "Pluto"
};
#define NBODIES 10

static double epochs[] = {
    2451545.0,                        /* J2000.0 */
    2460310.5,                        /* 2024-Jan-1 */
    2433282.5,                        /* 1950-Jan-1 */
    2488069.5,                        /* 2100-Jan-1 */
    2378496.5,                        /* 1800-Jan-1 */
    2451545.0 + 0.5,                  /* J2000 + half day */
    2451545.0 + 27.3,                 /* ~1 lunar month */
    2451545.0 + 365.25,               /* +1yr */
};
#define NEPOCHS 8

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_SPEED,                           "default" },
    { SEFLG_SPEED | SEFLG_J2000,             "J2000" },
    { SEFLG_SPEED | SEFLG_NONUT,             "NONUT" },
    { SEFLG_SPEED | SEFLG_EQUATORIAL,        "EQUATORIAL" },
    { SEFLG_SPEED | SEFLG_XYZ,              "XYZ" },
    { SEFLG_SPEED | SEFLG_RADIANS,           "RADIANS" },
    { SEFLG_SPEED | SEFLG_TRUEPOS,           "TRUEPOS" },
    { SEFLG_SPEED | SEFLG_NOABERR,           "NOABERR" },
    { SEFLG_SPEED | SEFLG_NOGDEFL,           "NOGDEFL" },
    { SEFLG_SPEED | SEFLG_ICRS,              "ICRS" },
    { SEFLG_SPEED3,                          "SPEED3" },
    { 0,                                     "no_speed" },
};
#define NFLAGS 12

int main(void) {
    char serr[256];
    int first = 1;
    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = SEFLG_SWIEPH | flag_combos[ifl].flag;
                double xx[6];
                memset(xx, 0, sizeof(xx));
                /* Reset C library state before each call so file caching does not
                 * carry over between test cases. This gives deterministic, stateless
                 * golden data that matches our stateless Rust implementation. */
                swe_close();
                swe_set_ephe_path("../../ephe");
                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "skipping: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flags);
                    continue;
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
    swe_close();
    return 0;
}
