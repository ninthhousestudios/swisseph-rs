/*
 * Generates golden reference data for sidereal calc via swe_calc with
 * SEFLG_SIDEREAL (default-branch modes only — no ECL_T0 / SSY_PLANE).
 *
 * Compile:
 *   cc -Wall -I../../../swisseph -o gen_ayanamsa_calc gen_ayanamsa_calc.c \
 *      ../../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_ayanamsa_calc > ../golden-data/ayanamsa_calc.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

/* Default-branch sidereal indices (no ECL_T0 / SSY_PLANE) */
static int indices[] = {0, 1, 3, 5, 26, 43};
#define NINDICES 6

static int bodies[] = {SE_SUN, SE_MOON, SE_MARS, SE_JUPITER, SE_MEAN_NODE};
static const char *body_names[] = {"Sun", "Moon", "Mars", "Jupiter", "MeanNode"};
#define NBODIES 5

static double epochs[] = {
    2415020.0,   /* J1900.0 */
    2433282.5,   /* ~B1950  */
    2451545.0,   /* J2000.0 */
    2455197.5,   /* 2010-Jan-1 */
    2458849.5,   /* 2020-Jan-1 */
};
#define NEPOCHS 5

int main(void) {
    char serr[256];
    double xx[6];
    int first_case = 1;

    printf("{\n");
    printf("  \"cases\": [\n");

    for (int ii = 0; ii < NINDICES; ii++) {
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS; ie++) {
                swe_close();
                swe_set_ephe_path(NULL);
                swe_set_sid_mode(indices[ii], 0, 0);

                int flags = SEFLG_MOSEPH | SEFLG_SIDEREAL | SEFLG_SPEED;
                memset(xx, 0, sizeof(xx));
                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "Error: idx=%d body=%d jd=%.1f: %s\n",
                            indices[ii], bodies[ib], epochs[ie], serr);
                    continue;
                }

                if (!first_case) printf(",\n");
                first_case = 0;

                printf("    {\"index\":%d,\"body\":\"%s\",\"tjd\":%.1f,"
                       "\"lon\":%.17g,\"lat\":%.17g,\"dist\":%.17g,"
                       "\"lon_speed\":%.17g}",
                       indices[ii], body_names[ib], epochs[ie],
                       xx[0], xx[1], xx[2], xx[3]);
            }
        }
    }

    printf("\n  ],\n");

    /* Equatorial group: Sun at J2000 under Lahiri (index 1).
     * Confirms equatorial output is tropical (same as calling without SIDEREAL). */
    printf("  \"equ\": [\n");

    int equ_first = 1;
    double equ_epochs[] = {2451545.0, 2451900.0};
    int n_equ = 2;

    for (int ie = 0; ie < n_equ; ie++) {
        /* Sidereal + equatorial */
        swe_close();
        swe_set_ephe_path(NULL);
        swe_set_sid_mode(1, 0, 0); /* Lahiri */

        int flags = SEFLG_MOSEPH | SEFLG_SIDEREAL | SEFLG_SPEED | SEFLG_EQUATORIAL;
        memset(xx, 0, sizeof(xx));
        int rc = swe_calc(equ_epochs[ie], SE_SUN, flags, xx, serr);
        if (rc < 0) { fprintf(stderr, "equ error: %s\n", serr); continue; }

        double sid_ra = xx[0], sid_dec = xx[1], sid_dist = xx[2];

        /* Tropical + equatorial (no SIDEREAL) */
        swe_close();
        swe_set_ephe_path(NULL);
        flags = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL;
        memset(xx, 0, sizeof(xx));
        rc = swe_calc(equ_epochs[ie], SE_SUN, flags, xx, serr);
        if (rc < 0) { fprintf(stderr, "trop equ error: %s\n", serr); continue; }

        double trop_ra = xx[0], trop_dec = xx[1];

        if (!equ_first) printf(",\n");
        equ_first = 0;

        printf("    {\"tjd\":%.1f,\"sid_ra\":%.17g,\"trop_ra\":%.17g,"
               "\"sid_dec\":%.17g,\"trop_dec\":%.17g,\"dist\":%.17g}",
               equ_epochs[ie], sid_ra, trop_ra, sid_dec, trop_dec, sid_dist);
    }

    printf("\n  ]\n}\n");
    return 0;
}
