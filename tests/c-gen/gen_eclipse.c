/*
 * Generates golden reference data for the eclipse module: swe_sol_eclipse_where
 * (RSE 5, swisseph-rs/72), swe_sol_eclipse_how (RSE 6, swisseph-rs/73),
 * swe_sol_eclipse_when_glob (RSE 7, swisseph-rs/74), swe_sol_eclipse_when_loc
 * (RSE 8, swisseph-rs/75), swe_lun_eclipse_how (RSE 9, swisseph-rs/76),
 * swe_lun_eclipse_when + swe_lun_eclipse_when_loc (RSE 10, swisseph-rs/77),
 * swe_lun_occult_where + swe_lun_occult_when_glob (RSE 11, swisseph-rs/78), and
 * swe_lun_occult_when_loc (RSE 12, swisseph-rs/79).
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
#include <string.h>
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

/* Occulted bodies for occ_where/occ_when_glob: two planets (Venus, Mars), one fixed star
 * (Aldebaran, via starname -- ipl is ignored by calc_planet_star/body_radius_au once starname is
 * non-empty, so any placeholder ipl works; -1 mirrors C's own "ipl<0 -> clamp to SE_SUN"
 * convention for star-only calls), and numbered-asteroid Pluto (SE_AST_OFFSET + 134340), which
 * all three swe_lun_occult_* entry points alias to SE_PLUTO before computing -- exercises that
 * aliasing end-to-end (swisseph-rs/92; with SEFLG_MOSEPH a bare asteroid would otherwise have no
 * ephemeris, so a wrong/absent alias fails loudly). */
struct occ_body { int32 ipl; const char *starname; };
static struct occ_body occ_bodies[] = {
    { SE_VENUS, NULL },
    { SE_MARS, NULL },
    { -1, "Aldebaran" },
    { SE_AST_OFFSET + 134340, NULL },
};
#define N_OCC_BODIES (sizeof(occ_bodies) / sizeof(occ_bodies[0]))

int main(void) {
    int first;
    /* Explicit path (rather than NULL's compiled-in default) so the occ_where/occ_when_glob
     * Aldebaran cases can find sefstars.txt (swisseph-rs/78) -- everything else in this file only
     * uses SEFLG_MOSEPH, which needs no ephemeris files. */
    swe_set_ephe_path("../../ephe");

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
    printf("\n  ],\n");

    /* === lun_how === */
    static double lun_how_tjd_uts[] = {
        2451919.347916667,   /* 2001-01-09 20:21 UT, total */
        2459360.971527778,   /* 2021-05-26 11:19 UT, total */
        2460571.613888889,   /* 2024-09-18 02:44 UT, partial */
    };
    #define N_LUN_HOW (sizeof(lun_how_tjd_uts) / sizeof(lun_how_tjd_uts[0]))
    static double lun_how_geopos[][3] = {
        { 8.55, 47.37, 500.0 },
    };
    #define N_LUN_HOW_GEOPOS (sizeof(lun_how_geopos) / sizeof(lun_how_geopos[0]))
    printf("  \"lun_how\": [\n");
    first = 1;
    for (size_t i = 0; i < N_LUN_HOW; i++) {
        for (size_t g = 0; g < N_LUN_HOW_GEOPOS; g++) {
            double tjd_ut = lun_how_tjd_uts[i];
            double geopos[3] = { lun_how_geopos[g][0], lun_how_geopos[g][1], lun_how_geopos[g][2] };
            double attr[20] = { 0 };
            char serr[256] = { 0 };
            int32 retval = swe_lun_eclipse_how(tjd_ut, SEFLG_MOSEPH, geopos, attr, serr);

            if (!first) printf(",\n");
            first = 0;
            printf("    {\"tjd_ut\": %.17g, \"geopos\": [%.17g, %.17g, %.17g], \"retval\": %d, \"attr\": [",
                   tjd_ut, geopos[0], geopos[1], geopos[2], retval);
            for (int k = 0; k < 11; k++) printf("%s%.20e", k ? ", " : "", attr[k]);
            printf("]}");
        }
    }
    printf("\n  ],\n");

    /* === lun_when === */
    static double lun_when_starts[] = { 2451545.0, 2459000.5 };
    #define N_LUN_WHEN_START (sizeof(lun_when_starts) / sizeof(lun_when_starts[0]))
    printf("  \"lun_when\": [\n");
    first = 1;
    for (size_t i = 0; i < N_LUN_WHEN_START; i++) {
        for (int backward = 0; backward <= 1; backward++) {
            double tjd_start = lun_when_starts[i];
            int32 ifltype = 0;
            double tret[10] = { 0 };
            char serr[256] = { 0 };
            int32 retval = swe_lun_eclipse_when(tjd_start, SEFLG_MOSEPH, ifltype, tret,
                                                  backward, serr);

            if (!first) printf(",\n");
            first = 0;
            printf("    {\"tjd_start\": %.17g, \"backward\": %s, \"retval\": %d, \"tret\": [",
                   tjd_start, backward ? "true" : "false", retval);
            for (int k = 0; k < 8; k++) printf("%s%.20e", k ? ", " : "", tret[k]);
            printf("]}");
        }
    }
    printf("\n  ],\n");

    /* === lun_when_loc === */
    static double lun_when_loc_geopos[][3] = {
        { 8.55, 47.37, 500.0 },
    };
    #define N_LUN_WHEN_LOC_GEOPOS (sizeof(lun_when_loc_geopos) / sizeof(lun_when_loc_geopos[0]))
    printf("  \"lun_when_loc\": [\n");
    first = 1;
    for (size_t g = 0; g < N_LUN_WHEN_LOC_GEOPOS; g++) {
        for (size_t i = 0; i < N_LUN_WHEN_START; i++) {
            for (int backward = 0; backward <= 1; backward++) {
                double geopos[3] = {
                    lun_when_loc_geopos[g][0], lun_when_loc_geopos[g][1], lun_when_loc_geopos[g][2]
                };
                double tjd_start = lun_when_starts[i];
                double tret[10] = { 0 };
                double attr[20] = { 0 };
                char serr[256] = { 0 };
                int32 retval = swe_lun_eclipse_when_loc(tjd_start, SEFLG_MOSEPH, geopos, tret,
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
    printf("\n  ],\n");

    /* === occ_where === */
    static double occ_where_tjd_ut = 2458800.5;
    printf("  \"occ_where\": [\n");
    first = 1;
    for (size_t i = 0; i < N_OCC_BODIES; i++) {
        double tjd_ut = occ_where_tjd_ut;
        int32 ipl = occ_bodies[i].ipl;
        /* swe_fixstar (reached via calc_planet_star when starname is set) strcpy()s into this
         * buffer in place (e.g. resolving traditional names) -- must be a writable buffer, not a
         * string-literal pointer, or it segfaults. */
        char starbuf[AS_MAXCH] = { 0 };
        if (occ_bodies[i].starname != NULL)
            strcpy(starbuf, occ_bodies[i].starname);
        char *starname = occ_bodies[i].starname != NULL ? starbuf : NULL;
        int32 ifl = SEFLG_MOSEPH;
        double geopos[10] = { 0 };
        double attr[20] = { 0 };
        char serr[256] = { 0 };
        int32 retval = swe_lun_occult_where(tjd_ut, ipl, starname, ifl, geopos, attr, serr);

        /* swe_lun_occult_where masks ifl &= SEFLG_EPHMASK AND aliases numbered-asteroid Pluto
         * (134340) to SE_PLUTO (swecl.c:619-623) before calling eclipse_where -- replicate BOTH
         * here so the dcore we capture (via the same test-only hook used for sol_where) matches
         * what swe_lun_occult_where's own internal eclipse_where call actually saw. Without the
         * alias, eclipse_where would take the `ipl > SE_AST_OFFSET` branch and read the
         * unpopulated swed.ast_diam (== 0 under MOSEPH), yielding a point-source drad=0 shadow
         * that disagrees with the aliased-to-Pluto geometry the geopos[] above already reflects. */
        int32 ipl_dcore = (ipl == SE_AST_OFFSET + 134340) ? SE_PLUTO : ipl;
        double geopos2[10] = { 0 };
        double dcore[10] = { 0 };
        char serr2[256] = { 0 };
        swi_test_eclipse_where_dcore(tjd_ut, ipl_dcore, starname, ifl & SEFLG_EPHMASK, geopos2,
                                      dcore, serr2);

        if (!first) printf(",\n");
        first = 0;
        printf("    {\"tjd_ut\": %.17g, \"ipl\": %d, \"starname\": %s, \"retval\": %d, \"geopos\": [",
               tjd_ut, ipl, starname ? "\"Aldebaran\"" : "null", retval);
        for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", geopos[k]);
        printf("], \"dcore\": [");
        for (int k = 0; k < 7; k++) printf("%s%.20e", k ? ", " : "", dcore[k]);
        printf("]}");
    }
    printf("\n  ],\n");

    /* === occ_when_glob === */
    static double occ_when_glob_tjd_start = 2451545.0;
    printf("  \"occ_when_glob\": [\n");
    first = 1;
    for (size_t i = 0; i < N_OCC_BODIES; i++) {
        for (int backward = 0; backward <= 1; backward++) {
            double tjd_start = occ_when_glob_tjd_start;
            int32 ipl = occ_bodies[i].ipl;
            char starbuf[AS_MAXCH] = { 0 };
            if (occ_bodies[i].starname != NULL)
                strcpy(starbuf, occ_bodies[i].starname);
            char *starname = occ_bodies[i].starname != NULL ? starbuf : NULL;
            int32 ifltype = 0;
            double tret[10] = { 0 };
            char serr[256] = { 0 };
            int32 retval = swe_lun_occult_when_glob(tjd_start, ipl, starname, SEFLG_MOSEPH,
                                                      ifltype, tret, backward, serr);

            if (!first) printf(",\n");
            first = 0;
            printf("    {\"tjd_start\": %.17g, \"ipl\": %d, \"starname\": %s, \"backward\": %s, "
                   "\"retval\": %d, \"tret\": [",
                   tjd_start, ipl, starname ? "\"Aldebaran\"" : "null",
                   backward ? "true" : "false", retval);
            for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", tret[k]);
            printf("]}");
        }
    }
    printf("\n  ],\n");

    /* === occ_when_loc === */
    /* Venus (planet, finite disc) + Aldebaran (star, point source) -- same two occulted-body
     * shapes as occ_where/occ_when_glob's first and third entries, skipping Mars to keep the
     * battery small (RSE 12 spec). */
    struct occ_body occ_when_loc_bodies[] = { occ_bodies[0], occ_bodies[2] };
    #define N_OCC_WHEN_LOC_BODIES \
        (sizeof(occ_when_loc_bodies) / sizeof(occ_when_loc_bodies[0]))
    static double occ_when_loc_geopos[3] = { 8.55, 47.37, 500.0 };
    static double occ_when_loc_tjd_start = 2451545.0;
    printf("  \"occ_when_loc\": [\n");
    first = 1;
    for (size_t i = 0; i < N_OCC_WHEN_LOC_BODIES; i++) {
        for (int backward = 0; backward <= 1; backward++) {
            double geopos[3] = {
                occ_when_loc_geopos[0], occ_when_loc_geopos[1], occ_when_loc_geopos[2]
            };
            double tjd_start = occ_when_loc_tjd_start;
            int32 ipl = occ_when_loc_bodies[i].ipl;
            char starbuf[AS_MAXCH] = { 0 };
            if (occ_when_loc_bodies[i].starname != NULL)
                strcpy(starbuf, occ_when_loc_bodies[i].starname);
            char *starname = occ_when_loc_bodies[i].starname != NULL ? starbuf : NULL;
            double tret[10] = { 0 };
            double attr[20] = { 0 };
            char serr[256] = { 0 };
            int32 retval = swe_lun_occult_when_loc(tjd_start, ipl, starname, SEFLG_MOSEPH,
                                                     geopos, tret, attr, backward, serr);

            if (!first) printf(",\n");
            first = 0;
            printf("    {\"geopos\": [%.17g, %.17g, %.17g], \"tjd_start\": %.17g, \"ipl\": %d, "
                   "\"starname\": %s, \"backward\": %s, \"retval\": %d, \"tret\": [",
                   geopos[0], geopos[1], geopos[2], tjd_start, ipl,
                   starname ? "\"Aldebaran\"" : "null", backward ? "true" : "false", retval);
            for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", tret[k]);
            printf("], \"attr\": [");
            for (int k = 0; k < 11; k++) printf("%s%.20e", k ? ", " : "", attr[k]);
            printf("]}");
        }
    }
    printf("\n  ],\n");

    /* === occ_when_glob_ifltype === */
    /* Exercises the ifltype filter/validation logic in swe_lun_occult_when_glob that the
     * ifltype=0 occ_when_glob battery never touches (swisseph-rs/92): the type-match retry that
     * skips non-matching occultations (TOTAL vs PARTIAL), and the two hard-error returns
     * (annular-for-planet, central-partial). retval < 0 (ERR) marks an expected-error case; tret
     * is then meaningless (left zeroed). All forward from J2000. */
    struct { int32 ipl; const char *starname; int32 ifltype; } occ_ifltype_cases[] = {
        { SE_VENUS, NULL, SE_ECL_TOTAL },                    /* matches 1st Venus occ (total) */
        { SE_VENUS, NULL, SE_ECL_PARTIAL },                  /* retries past totals to a partial */
        { SE_VENUS, NULL, SE_ECL_ANNULAR },                  /* ERR: annular impossible for planet */
        { SE_VENUS, NULL, SE_ECL_PARTIAL | SE_ECL_CENTRAL }, /* ERR: central-partial impossible */
        { -1, "Aldebaran", SE_ECL_TOTAL },                   /* star point-source is always total */
    };
    #define N_OCC_IFLTYPE_CASES (sizeof(occ_ifltype_cases) / sizeof(occ_ifltype_cases[0]))
    static double occ_ifltype_tjd_start = 2451545.0;
    printf("  \"occ_when_glob_ifltype\": [\n");
    first = 1;
    for (size_t i = 0; i < N_OCC_IFLTYPE_CASES; i++) {
        double tjd_start = occ_ifltype_tjd_start;
        int32 ipl = occ_ifltype_cases[i].ipl;
        char starbuf[AS_MAXCH] = { 0 };
        if (occ_ifltype_cases[i].starname != NULL)
            strcpy(starbuf, occ_ifltype_cases[i].starname);
        char *starname = occ_ifltype_cases[i].starname != NULL ? starbuf : NULL;
        int32 ifltype = occ_ifltype_cases[i].ifltype;
        double tret[10] = { 0 };
        char serr[256] = { 0 };
        int32 retval = swe_lun_occult_when_glob(tjd_start, ipl, starname, SEFLG_MOSEPH,
                                                  ifltype, tret, 0, serr);

        if (!first) printf(",\n");
        first = 0;
        printf("    {\"tjd_start\": %.17g, \"ipl\": %d, \"starname\": %s, \"ifltype\": %d, "
               "\"retval\": %d, \"tret\": [",
               tjd_start, ipl, starname ? "\"Aldebaran\"" : "null", ifltype, retval);
        for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", tret[k]);
        printf("]}");
    }
    printf("\n  ],\n");

    /* === occ_where_asteroid === */
    /* Eros (433) via SWIEPH, exercising body_radius_au's asteroid-metadata branch. A dummy
     * swe_calc pre-populates C's swed.ast_diam global (which body_radius_au reads) — without it,
     * the first eclipse_where call for this asteroid would see drad=0 (file not yet opened).
     * The stateless Rust port reads from the already-loaded file and doesn't have this ordering
     * dependency. */
    printf("  \"occ_where_asteroid\": [\n");
    first = 1;
    {
        double x[6];
        char serr[256] = { 0 };
        swe_calc(2458800.5, SE_AST_OFFSET + 433, SEFLG_SWIEPH, x, serr);

        double tjd_ut = 2458800.5;
        int32 ipl = SE_AST_OFFSET + 433;
        int32 ifl = SEFLG_SWIEPH;
        double geopos[10] = { 0 };
        double attr[20] = { 0 };
        swi_test_eclipse_where_dcore(tjd_ut, ipl, NULL, ifl & SEFLG_EPHMASK, geopos, attr, serr);

        double geopos2[10] = { 0 };
        double attr2[20] = { 0 };
        int32 retval = swe_lun_occult_where(tjd_ut, ipl, NULL, ifl, geopos2, attr2, serr);

        double dcore[10] = { 0 };
        swi_test_eclipse_where_dcore(tjd_ut, ipl, NULL, ifl & SEFLG_EPHMASK, geopos, dcore, serr);

        printf("    {\"tjd_ut\": %.17g, \"ipl\": %d, \"starname\": null, \"retval\": %d, \"geopos\": [",
               tjd_ut, ipl, retval);
        for (int k = 0; k < 10; k++) printf("%s%.20e", k ? ", " : "", geopos2[k]);
        printf("], \"dcore\": [");
        for (int k = 0; k < 7; k++) printf("%s%.20e", k ? ", " : "", dcore[k]);
        printf("]}");
    }
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
