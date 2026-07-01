/*
 * Generates golden reference data for the eclipse module: swe_sol_eclipse_where
 * (RSE 5, swisseph-rs/72). Later RSE tasks (6-12) add more keys to this same file.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_eclipse tests/c-gen/gen_eclipse.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_eclipse > tests/golden-data/eclipse.json
 */

#include <stdio.h>
#include "swephexp.h"

/* Instants of maximum eclipse (UT) for known central solar eclipses, plus a plain-noon epoch
 * that is nowhere near a solar conjunction (exercises the "no eclipse anywhere" / retval==0
 * path with a well-behaved geometry, unlike an arbitrary near-miss date). */
static double sol_where_tjd_uts[] = {
    2451401.9604166667,   /* 1999-08-11 11:03 UT, total */
    2459375.9458333333,   /* 2021-06-10 10:42 UT, annular */
    2460586.28125,        /* 2024-10-02 18:45 UT, annular */
    2451545.0,            /* 2000-01-01 12:00 UT, no eclipse */
};
static int sol_where_nonut[] = { 0, 0, 1, 0 };
#define N_SOL_WHERE (sizeof(sol_where_tjd_uts) / sizeof(sol_where_tjd_uts[0]))

int main(void) {
    int first;
    swe_set_ephe_path(NULL);

    printf("{\n");

    /* === sol_where === */
    printf("  \"sol_where\": [\n");
    first = 1;
    for (size_t i = 0; i < N_SOL_WHERE; i++) {
        double tjd_ut = sol_where_tjd_uts[i];
        int32 ifl = SEFLG_MOSEPH | (sol_where_nonut[i] ? SEFLG_NONUT : 0);
        double geopos[10] = { 0 };
        double attr[20] = { 0 };
        char serr[256] = { 0 };
        int32 retval = swe_sol_eclipse_where(tjd_ut, ifl, geopos, attr, serr);
        if (!first) printf(",\n");
        first = 0;
        printf("    {\"tjd_ut\": %.17g, \"nonut\": %s, \"retval\": %d, \"geopos\": [",
               tjd_ut, sol_where_nonut[i] ? "true" : "false", retval);
        for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", geopos[k]);
        printf("], \"attr\": [");
        for (int k = 0; k < 11; k++) printf("%s%.20e", k ? ", " : "", attr[k]);
        printf("]}");
    }
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
