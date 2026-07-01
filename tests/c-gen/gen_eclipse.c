/*
 * Generates golden reference data for the eclipse module: swe_sol_eclipse_where
 * (RSE 5, swisseph-rs/72), swe_sol_eclipse_how (RSE 6, swisseph-rs/73),
 * swe_sol_eclipse_when_glob (RSE 7, swisseph-rs/74), and swe_sol_eclipse_when_loc
 * (RSE 8, swisseph-rs/75). Later RSE tasks (9-12) add more keys to this same file.
 *
 * The "dcore" field in each sol_where case comes from swi_test_eclipse_where_dcore, a
 * non-static test-only hook added to ../swisseph/swecl.c (right after calc_planet_star) that
 * forwards the full dcore[0..9] array out of the static eclipse_where(). C's public
 * swe_sol_eclipse_where only exposes dcore[0] (via attr[3]); the rest of EclipseWhere's fields
 * (penumbra_diameter_km, shadow_axis_distance_km, both fundamental-plane diameters, both cone
 * half-angle cosines) have no other C oracle. Requires libswe.a to be rebuilt after that patch
 * (`cc -g -Wall -fPIC -c swecl.c -o swecl.o && ar r libswe.a swecl.o` from ../swisseph).
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_eclipse tests/c-gen/gen_eclipse.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_eclipse > tests/golden-data/eclipse.json
 */

#include <stdio.h>
#include "swephexp.h"

/* swi_test_eclipse_where_dcore: see ../swisseph/swecl.c patch note above. Not declared in any
 * public header. */
extern int32 swi_test_eclipse_where_dcore(double tjd_ut, int32 ipl, char *starname, int32 ifl,
                                           double *geopos, double *dcore, char *serr);

/* Local to swecl.c (not in any public header); mirrored here to mask ifl the same way
 * swe_sol_eclipse_where does before its own eclipse_where call (swecl.c:67,574). */
#define SEFLG_EPHMASK (SEFLG_JPLEPH | SEFLG_SWIEPH | SEFLG_MOSEPH)

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

/* Observers for sol_how: one near-central (close to the 1999/2021/2024 eclipse tracks), one
 * off-track -- crossed with the same sol_where epochs (incl. the no-eclipse epoch, to exercise
 * the "no eclipse visible here" clearing path). */
static double sol_how_geopos[][3] = {
    { 8.55, 47.37, 500.0 },
    { -100.0, 40.0, 0.0 },
};
#define N_SOL_HOW_GEOPOS (sizeof(sol_how_geopos) / sizeof(sol_how_geopos[0]))

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

        /* swe_sol_eclipse_where masks ifl &= SEFLG_EPHMASK before calling eclipse_where
         * (swecl.c:574) -- replicate that here so the dcore we capture matches what
         * swe_sol_eclipse_where's own eclipse_where call actually saw. */
        double geopos2[10] = { 0 };
        double dcore[10] = { 0 };
        char serr2[256] = { 0 };
        swi_test_eclipse_where_dcore(tjd_ut, SE_SUN, NULL, ifl & SEFLG_EPHMASK, geopos2, dcore, serr2);

        if (!first) printf(",\n");
        first = 0;
        printf("    {\"tjd_ut\": %.17g, \"nonut\": %s, \"retval\": %d, \"geopos\": [",
               tjd_ut, sol_where_nonut[i] ? "true" : "false", retval);
        for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", geopos[k]);
        printf("], \"attr\": [");
        for (int k = 0; k < 11; k++) printf("%s%.20e", k ? ", " : "", attr[k]);
        printf("], \"dcore\": [");
        for (int k = 0; k < 7; k++) printf("%s%.20e", k ? ", " : "", dcore[k]);
        printf("]}");
    }
    printf("\n  ],\n");

    /* === sol_how === */
    printf("  \"sol_how\": [\n");
    first = 1;
    for (size_t i = 0; i < N_SOL_WHERE; i++) {
        double tjd_ut = sol_where_tjd_uts[i];
        int32 ifl = SEFLG_MOSEPH;
        for (size_t g = 0; g < N_SOL_HOW_GEOPOS; g++) {
            double geopos[3] = { sol_how_geopos[g][0], sol_how_geopos[g][1], sol_how_geopos[g][2] };
            double attr[20] = { 0 };
            char serr[256] = { 0 };
            int32 retval = swe_sol_eclipse_how(tjd_ut, ifl, geopos, attr, serr);

            if (!first) printf(",\n");
            first = 0;
            printf("    {\"tjd_ut\": %.17g, \"geopos\": [%.17g, %.17g, %.17g], \"retval\": %d, \"attr\": [",
                   tjd_ut, geopos[0], geopos[1], geopos[2], retval);
            for (int k = 0; k < 11; k++) printf("%s%.20e", k ? ", " : "", attr[k]);
            printf("]}");
        }
    }
    printf("\n  ],\n");

    /* === sol_when_glob === */
    static double glob_tjd_starts[] = { 2451545.0, 2459000.5 };
    #define N_GLOB_START (sizeof(glob_tjd_starts) / sizeof(glob_tjd_starts[0]))
    printf("  \"sol_when_glob\": [\n");
    first = 1;
    for (size_t i = 0; i < N_GLOB_START; i++) {
        for (int backward = 0; backward <= 1; backward++) {
            double tjd_start = glob_tjd_starts[i];
            int32 ifltype = 0;
            double tret[10] = { 0 };
            char serr[256] = { 0 };
            int32 retval = swe_sol_eclipse_when_glob(tjd_start, SEFLG_MOSEPH, ifltype, tret,
                                                       backward, serr);

            if (!first) printf(",\n");
            first = 0;
            printf("    {\"tjd_start\": %.17g, \"backward\": %s, \"retval\": %d, \"tret\": [",
                   tjd_start, backward ? "true" : "false", retval);
            for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", tret[k]);
            printf("]}");
        }
    }
    printf("\n  ],\n");

    /* === sol_when_loc === */
    static double loc_geopos[][3] = {
        { 8.55, 47.37, 500.0 },     /* near-central for the 1999/2021/2024 tracks (sol_where set) */
        { -71.0, -33.0, 500.0 },    /* Chile: near the 2019/2020 total/annular tracks */
    };
    #define N_LOC_GEOPOS (sizeof(loc_geopos) / sizeof(loc_geopos[0]))
    static double loc_tjd_starts[] = { 2451545.0, 2458800.5 };
    #define N_LOC_START (sizeof(loc_tjd_starts) / sizeof(loc_tjd_starts[0]))
    printf("  \"sol_when_loc\": [\n");
    first = 1;
    for (size_t g = 0; g < N_LOC_GEOPOS; g++) {
        for (size_t i = 0; i < N_LOC_START; i++) {
            for (int backward = 0; backward <= 1; backward++) {
                double geopos[3] = { loc_geopos[g][0], loc_geopos[g][1], loc_geopos[g][2] };
                double tjd_start = loc_tjd_starts[i];
                double tret[10] = { 0 };
                double attr[20] = { 0 };
                char serr[256] = { 0 };
                int32 retval = swe_sol_eclipse_when_loc(tjd_start, SEFLG_MOSEPH, geopos, tret,
                                                          attr, backward, serr);

                if (!first) printf(",\n");
                first = 0;
                printf("    {\"geopos\": [%.17g, %.17g, %.17g], \"tjd_start\": %.17g, "
                       "\"backward\": %s, \"retval\": %d, \"tret\": [",
                       geopos[0], geopos[1], geopos[2], tjd_start, backward ? "true" : "false",
                       retval);
                for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", tret[k]);
                printf("], \"attr\": [");
                for (int k = 0; k < 11; k++) printf("%s%.20e", k ? ", " : "", attr[k]);
                printf("]}");
            }
        }
    }
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
