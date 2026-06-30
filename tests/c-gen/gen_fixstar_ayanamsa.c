/*
 * Generates golden reference data for the 12 fixed-star ayanamsa modes.
 * Tests swe_get_ayanamsa_ex for each mode, plus swe_calc with SEFLG_SIDEREAL
 * for a couple of indices to validate the apply_sidereal path.
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../../swisseph -o gen_fixstar_ayanamsa gen_fixstar_ayanamsa.c \
 *      ../../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_fixstar_ayanamsa > ../golden-data/fixstar_ayanamsa.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static int fixed_star_indices[] = {17, 27, 28, 29, 30, 31, 32, 33, 35, 36, 39, 40};
#define NINDICES 12

static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2415020.0,   /* J1900.0 */
    2469807.5,   /* 2050 Jan 1 */
    2086302.5,   /* ~1000 AD */
};
#define NEPOCHS 4

/* Indices to also test via swe_calc with SEFLG_SIDEREAL */
static int calc_indices[] = {17, 27};
#define NCALC_INDICES 2

int main(void) {
    double daya;
    double xx[6];
    char serr[512];
    int retflag;
    int i, j;
    int first_outer;

    swe_set_ephe_path("../../../swisseph/ephe");

    printf("{\n");

    /* === ayanamsa cases === */
    printf("  \"ayanamsa\": [\n");
    first_outer = 1;
    for (i = 0; i < NINDICES; i++) {
        int idx = fixed_star_indices[i];
        swe_set_sid_mode(idx, 0.0, 0.0);
        for (j = 0; j < NEPOCHS; j++) {
            double tjd = epochs[j];
            daya = 0.0;
            memset(serr, 0, sizeof(serr));
            retflag = swe_get_ayanamsa_ex(tjd, SEFLG_MOSEPH, &daya, serr);
            if (!first_outer) printf(",\n");
            first_outer = 0;
            printf("    {\"idx\": %d, \"tjd\": %.1f, \"daya\": %.17g, \"retflag\": %d}",
                   idx, tjd, daya, retflag);
        }
    }
    printf("\n  ],\n");

    /* === calc (SIDEREAL) cases — validate apply_sidereal path === */
    printf("  \"calc\": [\n");
    first_outer = 1;
    for (i = 0; i < NCALC_INDICES; i++) {
        int idx = calc_indices[i];
        swe_set_sid_mode(idx, 0.0, 0.0);
        for (j = 0; j < NEPOCHS; j++) {
            double tjd = epochs[j];
            memset(xx, 0, sizeof(xx));
            memset(serr, 0, sizeof(serr));
            swe_calc(tjd, SE_SUN, SEFLG_MOSEPH | SEFLG_SIDEREAL, xx, serr);
            if (!first_outer) printf(",\n");
            first_outer = 0;
            printf("    {\"idx\": %d, \"tjd\": %.1f, \"lon\": %.17g}",
                   idx, tjd, xx[0]);
        }
    }
    printf("\n  ]\n");

    printf("}\n");
    swe_close();
    return 0;
}
