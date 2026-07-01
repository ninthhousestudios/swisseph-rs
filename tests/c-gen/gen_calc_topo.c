/*
 * Generates golden reference data for topocentric swe_calc (SEFLG_TOPOCTR).
 * Same body/epoch/flag shape as gen_calc.c but iterates over several
 * observer positions via swe_set_topo.
 *
 * Compile:
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_calc_topo tests/c-gen/gen_calc_topo.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_calc_topo > tests/golden-data/calc_topo.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

struct observer {
    double lon, lat, alt;
};

static struct observer observers[] = {
    { 8.55, 47.37, 500.0 },     /* Zurich */
    { -74.0, 40.7, 10.0 },      /* New York */
    { 139.7, 35.7, 40.0 },      /* Tokyo */
};
#define NOBSERVERS 3

static int bodies[] = { SE_SUN, SE_MOON, SE_MERCURY, SE_MARS, SE_VENUS };
static const char *body_names[] = { "Sun", "Moon", "Mercury", "Mars", "Venus" };
#define NBODIES 5

static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2378496.5,    /* 1800-Jan-1 — sepl_18 tfstart, SPEED3 file-boundary case */
    2469807.5,    /* 2050-Jan-1 */
};
#define NEPOCHS 3

int main(void) {
    char serr[256];
    swe_set_ephe_path(NULL);
    int flags = SEFLG_MOSEPH | SEFLG_TOPOCTR | SEFLG_EQUATORIAL | SEFLG_SPEED;
    int first = 1;
    printf("[\n");
    for (int io = 0; io < NOBSERVERS; io++) {
        swe_set_topo(observers[io].lon, observers[io].lat, observers[io].alt);
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS; ie++) {
                double xx[6];
                memset(xx, 0, sizeof(xx));
                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "skipping: %s body=%d jd=%.1f\n",
                            serr, bodies[ib], epochs[ie]);
                    continue;
                }
                if (!first) printf(",\n");
                first = 0;
                printf("  {\"lon\": %.20e, \"lat\": %.20e, \"alt\": %.20e, "
                       "\"body\": %d, \"body_name\": \"%s\", \"jd\": %.20e, "
                       "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                       observers[io].lon, observers[io].lat, observers[io].alt,
                       bodies[ib], body_names[ib], epochs[ie],
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
            }
        }
    }
    printf("\n]\n");
    return 0;
}
