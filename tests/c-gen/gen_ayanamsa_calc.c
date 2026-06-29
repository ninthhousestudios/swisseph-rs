/*
 * Generates golden reference data for sidereal calc via swe_calc with
 * SEFLG_SIDEREAL (default-branch modes, ECL_T0, and SSY_PLANE).
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

/* ECL_T0 / SSY_PLANE bodies and epochs */
static int ecl_bodies[] = {SE_SUN, SE_MOON, SE_MARS};
static const char *ecl_body_names[] = {"Sun", "Moon", "Mars"};
#define NECL_BODIES 3

static double ecl_epochs[] = {2415020.0, 2451545.0, 2458849.5};
#define NECL_EPOCHS 3

/* ECL_T0 ayanamsa indices (auto-set ECL_T0 bit) */
static int ecl_t0_indices[] = {
    SE_SIDM_J2000,           /* 18 */
    SE_SIDM_J1900,           /* 19 */
    SE_SIDM_B1950,           /* 20 */
    SE_SIDM_GALALIGN_MARDYKS /* 34 */
};
#define NECL_T0_INDICES 4

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

    printf("\n  ],\n");

    /* ECL_T0 group: modes 18/19/20/34 x bodies {Sun,Moon,Mars} x epochs */
    printf("  \"ecl_t0\": [\n");

    int ecl_first = 1;
    for (int ii = 0; ii < NECL_T0_INDICES; ii++) {
        for (int ib = 0; ib < NECL_BODIES; ib++) {
            for (int ie = 0; ie < NECL_EPOCHS; ie++) {
                swe_close();
                swe_set_ephe_path(NULL);
                swe_set_sid_mode(ecl_t0_indices[ii], 0, 0);

                int flags = SEFLG_MOSEPH | SEFLG_SIDEREAL | SEFLG_SPEED;
                memset(xx, 0, sizeof(xx));
                int rc = swe_calc(ecl_epochs[ie], ecl_bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "ECL_T0 error: idx=%d body=%d jd=%.1f: %s\n",
                            ecl_t0_indices[ii], ecl_bodies[ib], ecl_epochs[ie], serr);
                    continue;
                }

                if (!ecl_first) printf(",\n");
                ecl_first = 0;

                printf("    {\"index\":%d,\"body\":\"%s\",\"tjd\":%.1f,"
                       "\"lon\":%.17g,\"lat\":%.17g,\"dist\":%.17g,"
                       "\"lon_speed\":%.17g}",
                       ecl_t0_indices[ii], ecl_body_names[ib], ecl_epochs[ie],
                       xx[0], xx[1], xx[2], xx[3]);
            }
        }
    }

    printf("\n  ],\n");

    /* USER-ECL_T0 group: SE_SIDM_USER|SE_SIDBIT_ECL_T0, t0=J2000, ayan=25.0 */
    printf("  \"user_ecl_t0\": [\n");

    int user_first = 1;
    for (int ib = 0; ib < NECL_BODIES; ib++) {
        for (int ie = 0; ie < NECL_EPOCHS; ie++) {
            swe_close();
            swe_set_ephe_path(NULL);
            swe_set_sid_mode(SE_SIDM_USER | SE_SIDBIT_ECL_T0, 2451545.0, 25.0);

            int flags = SEFLG_MOSEPH | SEFLG_SIDEREAL | SEFLG_SPEED;
            memset(xx, 0, sizeof(xx));
            int rc = swe_calc(ecl_epochs[ie], ecl_bodies[ib], flags, xx, serr);
            if (rc < 0) {
                fprintf(stderr, "USER-ECL_T0 error: body=%d jd=%.1f: %s\n",
                        ecl_bodies[ib], ecl_epochs[ie], serr);
                continue;
            }

            if (!user_first) printf(",\n");
            user_first = 0;

            printf("    {\"body\":\"%s\",\"tjd\":%.1f,"
                   "\"lon\":%.17g,\"lat\":%.17g,\"dist\":%.17g,"
                   "\"lon_speed\":%.17g}",
                   ecl_body_names[ib], ecl_epochs[ie],
                   xx[0], xx[1], xx[2], xx[3]);
        }
    }

    printf("\n  ],\n");

    /* SSY_PLANE group: SE_SIDM_LAHIRI|SE_SIDBIT_SSY_PLANE, Sun at J2000 */
    printf("  \"ssy\": [\n");

    swe_close();
    swe_set_ephe_path(NULL);
    swe_set_sid_mode(SE_SIDM_LAHIRI | SE_SIDBIT_SSY_PLANE, 0, 0);

    int flags = SEFLG_MOSEPH | SEFLG_SIDEREAL | SEFLG_SPEED;
    memset(xx, 0, sizeof(xx));
    int rc = swe_calc(2451545.0, SE_SUN, flags, xx, serr);
    if (rc < 0) {
        fprintf(stderr, "SSY_PLANE error: %s\n", serr);
    } else {
        printf("    {\"body\":\"Sun\",\"tjd\":2451545.0,"
               "\"lon\":%.17g,\"lat\":%.17g,\"dist\":%.17g,"
               "\"lon_speed\":%.17g}",
               xx[0], xx[1], xx[2], xx[3]);
    }

    printf("\n  ],\n");

    /* SPEED3 group: SEFLG_SPEED3 (no SEFLG_SPEED) triggers C's use_speed3,
     * which calls swecalc 3x with SEFLG_SIDEREAL and differences the projected
     * positions. Covers a default-branch mode (Lahiri) and an ECL_T0 mode
     * (J2000) so the projected-then-differenced speed is exercised on both. */
    printf("  \"speed3\": [\n");

    int s3_indices[] = {SE_SIDM_LAHIRI, SE_SIDM_J2000}; /* 1 (default), 18 (ECL_T0) */
    int s3_bodies[] = {SE_SUN, SE_MOON, SE_MARS};
    const char *s3_body_names[] = {"Sun", "Moon", "Mars"};
    double s3_epochs[] = {2451545.0, 2458849.5};
    int s3_first = 1;

    for (int ii = 0; ii < 2; ii++) {
        for (int ib = 0; ib < 3; ib++) {
            for (int ie = 0; ie < 2; ie++) {
                swe_close();
                swe_set_ephe_path(NULL);
                swe_set_sid_mode(s3_indices[ii], 0, 0);

                int f3 = SEFLG_MOSEPH | SEFLG_SIDEREAL | SEFLG_SPEED3;
                memset(xx, 0, sizeof(xx));
                int rc3 = swe_calc(s3_epochs[ie], s3_bodies[ib], f3, xx, serr);
                if (rc3 < 0) {
                    fprintf(stderr, "SPEED3 error: idx=%d body=%d jd=%.1f: %s\n",
                            s3_indices[ii], s3_bodies[ib], s3_epochs[ie], serr);
                    continue;
                }

                if (!s3_first) printf(",\n");
                s3_first = 0;

                printf("    {\"index\":%d,\"body\":\"%s\",\"tjd\":%.1f,"
                       "\"lon\":%.17g,\"lat\":%.17g,\"dist\":%.17g,"
                       "\"lon_speed\":%.17g}",
                       s3_indices[ii], s3_body_names[ib], s3_epochs[ie],
                       xx[0], xx[1], xx[2], xx[3]);
            }
        }
    }

    printf("\n  ]\n}\n");
    return 0;
}
