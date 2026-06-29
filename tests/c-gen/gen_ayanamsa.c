/*
 * gen_ayanamsa.c — generate golden data for ayanamsa tests
 *
 * Compile:
 *   cc -Wall -I../../../swisseph -o gen_ayanamsa gen_ayanamsa.c \
 *       ../../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_ayanamsa > ../golden-data/ayanamsa.json
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "swephexp.h"

/* Declare the internal function not in the public header */
extern int32 swi_get_ayanamsa_with_speed(double tjd_et, int32 iflag,
                                          double *daya, char *serr);

/* Non-deferred ayanamsa indices (fixed-star indices excluded) */
static int indices[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    /* 17 deferred */
    18, 19, 20, 21, 22, 23, 24, 25, 26,
    /* 27-33 deferred */
    34,
    /* 35-36 deferred */
    37, 38,
    /* 39-40 deferred */
    41, 42, 43, 44, 45, 46
};
static int n_indices = sizeof(indices) / sizeof(indices[0]);

/* Test epochs (TT/ET) */
static double epochs[] = {
    1356175.5,   /* ~-1100 CE */
    1757583.5,   /* ~0 CE    */
    2415020.0,   /* J1900    */
    2433282.5,   /* B1950    */
    2451545.0,   /* J2000    */
    2455197.5,   /* 2010     */
    2458849.5,   /* 2020     */
    2488070.0,   /* 2100     */
};
static int n_epochs = sizeof(epochs) / sizeof(epochs[0]);

int main(void) {
    char serr[256];
    double daya, dd[2];
    int first_case = 1;

    printf("{\n");
    printf("  \"cases\": [\n");

    for (int ii = 0; ii < n_indices; ii++) {
        int idx = indices[ii];
        for (int ej = 0; ej < n_epochs; ej++) {
            double tjd = epochs[ej];

            /* Reset state for each measurement */
            swe_close();
            swe_set_ephe_path(NULL);
            swe_set_sid_mode(idx, 0, 0);

            /* with_nut: default (includes nutation) */
            daya = 0.0;
            memset(serr, 0, sizeof(serr));
            swe_get_ayanamsa_ex(tjd, SEFLG_MOSEPH, &daya, serr);
            double with_nut = daya;

            /* no_nut */
            swe_close();
            swe_set_ephe_path(NULL);
            swe_set_sid_mode(idx, 0, 0);
            daya = 0.0;
            memset(serr, 0, sizeof(serr));
            swe_get_ayanamsa_ex(tjd, SEFLG_MOSEPH | SEFLG_NONUT, &daya, serr);
            double no_nut = daya;

            /* speed */
            swe_close();
            swe_set_ephe_path(NULL);
            swe_set_sid_mode(idx, 0, 0);
            dd[0] = 0.0; dd[1] = 0.0;
            memset(serr, 0, sizeof(serr));
            swi_get_ayanamsa_with_speed(tjd, SEFLG_MOSEPH, dd, serr);
            double speed = dd[1];

            if (!first_case) printf(",\n");
            first_case = 0;

            printf("    {\"index\":%d,\"tjd\":%.10f,"
                   "\"with_nut\":%.15g,"
                   "\"no_nut\":%.15g,"
                   "\"speed\":%.15g}",
                   idx, tjd, with_nut, no_nut, speed);
        }
    }

    printf("\n  ],\n");

    /* USER mode case */
    printf("  \"user\": [\n");
    swe_close();
    swe_set_ephe_path(NULL);
    swe_set_sid_mode(SE_SIDM_USER, 2433282.5, 24.0);
    daya = 0.0;
    memset(serr, 0, sizeof(serr));
    swe_get_ayanamsa_ex(2451545.0, SEFLG_MOSEPH | SEFLG_NONUT, &daya, serr);
    printf("    {\"t0\":2433282.5,\"ayan_t0\":24.0,\"tjd\":2451545.0,\"no_nut\":%.15g}\n",
           daya);
    printf("  ]\n");
    printf("}\n");

    swe_close();
    return 0;
}
