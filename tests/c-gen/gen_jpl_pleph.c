/*
 * Generates golden reference data for jpl_pleph (swi_pleph) interpolation.
 * Forces the JPL file open via swe_calc, then calls swi_pleph directly for a
 * grid of bodies and epochs, emitting JSON consumed by tests/golden/jpl_pleph.rs.
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../../swisseph -o gen_jpl_pleph gen_jpl_pleph.c \
 *      ../../../swisseph/libswe.a -lm
 * Run (from tests/c-gen/):
 *   ./gen_jpl_pleph > ../golden-data/jpl_pleph.json
 */

#include <stdio.h>
#include <math.h>
#include "swephexp.h"
#include "swejpl.h"

/* J_* body index constants (swejpl.h:68-83) */
#define J_MERCURY 0
#define J_VENUS   1
#define J_EARTH   2
#define J_MARS    3
#define J_JUPITER 4
#define J_SATURN  5
#define J_URANUS  6
#define J_NEPTUNE 7
#define J_PLUTO   8
#define J_MOON    9
#define J_SUN    10
#define J_SBARY  11
#define J_EMB    12

/* Reuse the 7 epochs from calc.json */
static double epochs[] = {
    2159345.5,   /* ~1200 CE */
    2305447.5,   /* ~1600 CE */
    2378496.5,   /* 1800-Jan-1 */
    2433282.5,   /* ~B1950 */
    2451545.0,   /* J2000 */
    2460310.5,   /* 2024-Jan-1 */
    2488069.5,   /* 2100-Jan-1 */
};
#define NEPOCHS 7

int main(void) {
    double xx[6];
    double rrd[6];
    char serr[256];
    int first = 1;

    swe_set_ephe_path("../../ephe");
    swe_set_jpl_file("de441.eph");

    /* Force-open the JPL file */
    if (swe_calc(epochs[0], SE_SUN, SEFLG_JPLEPH | SEFLG_SPEED, xx, serr) < 0) {
        fprintf(stderr, "swe_calc failed: %s\n", serr);
        return 1;
    }

    printf("{\"cases\": [\n");

    for (int ie = 0; ie < NEPOCHS; ie++) {
        double jd = epochs[ie];

        /* Barycentric cases: all bodies vs J_SBARY */
        int ntargs[] = {J_MERCURY, J_VENUS, J_EARTH, J_MARS, J_JUPITER,
                        J_SATURN, J_URANUS, J_NEPTUNE, J_PLUTO, J_MOON, J_SUN};
        int nntargs = (int)(sizeof(ntargs) / sizeof(ntargs[0]));

        for (int it = 0; it < nntargs; it++) {
            int ntarg = ntargs[it];
            int ret = swi_pleph(jd, ntarg, J_SBARY, rrd, serr);
            if (ret != 0) {
                fprintf(stderr, "swi_pleph failed ntarg=%d jd=%.1f: %s\n",
                        ntarg, jd, serr);
                continue;
            }
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"ntarg\": %d, \"ncent\": %d, \"jd\": %.20e,"
                   " \"rrd\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                   ntarg, J_SBARY, jd,
                   rrd[0], rrd[1], rrd[2], rrd[3], rrd[4], rrd[5]);
        }

        /* Geocentric Moon */
        int ret = swi_pleph(jd, J_MOON, J_EARTH, rrd, serr);
        if (ret != 0) {
            fprintf(stderr, "swi_pleph failed Moon/Earth jd=%.1f: %s\n", jd, serr);
        } else {
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"ntarg\": %d, \"ncent\": %d, \"jd\": %.20e,"
                   " \"rrd\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                   J_MOON, J_EARTH, jd,
                   rrd[0], rrd[1], rrd[2], rrd[3], rrd[4], rrd[5]);
        }
    }

    printf("\n]}\n");
    swe_close();
    return 0;
}
