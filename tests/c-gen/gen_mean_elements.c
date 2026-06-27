/*
 * Generates golden reference data for mean node, mean apogee, and ECL_NUT.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_mean_elements gen_mean_elements.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_mean_elements > ../golden-data/mean_elements.json
 */

#include <stdio.h>
#include "swephexp.h"

static int bodies[] = {
    SE_MEAN_NODE, SE_MEAN_APOG, SE_ECL_NUT
};
static const char *body_names[] = {
    "MeanNode", "MeanApogee", "EclipticNutation"
};
#define NBODIES 3

static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2460310.5,    /* 2024-Jan-1 */
    2433282.5,    /* 1950-Jan-1 */
    2488069.5,    /* 2100-Jan-1 */
    2378496.5,    /* 1800-Jan-1 */
    2305447.5,    /* 1600-Jan-1 */
    2159345.5,    /* 1200-Jan-1 */
    2013243.5,    /* 800-Jan-1 */
    1720693.5,    /* 0 CE */
    1538187.5,    /* -500 */
    990557.5,     /* -2000 */
};
#define NEPOCHS 11

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { 0,                "default" },
    { SEFLG_J2000,      "J2000" },
    { SEFLG_NONUT,      "NONUT" },
    { SEFLG_EQUATORIAL, "EQUATORIAL" },
    { SEFLG_XYZ,        "XYZ" },
};
#define NFLAGS 5

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
