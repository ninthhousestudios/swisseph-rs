/*
 * Generates golden reference data for swe_houses_armc_ex2 (house systems).
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_houses tests/c-gen/gen_houses.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_houses > tests/golden-data/houses.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

/* Battery reused across all houses sub-tasks. */
static double armcs[] = { 0.0, 47.5, 123.4, 215.7, 290.1, 350.0 };
#define N_ARMC (sizeof(armcs) / sizeof(armcs[0]))

static double geolats[] = { 51.5, -33.87, 0.0, 64.0, -64.0 };
#define N_GEOLAT (sizeof(geolats) / sizeof(geolats[0]))

static double epss[] = { 23.4392911 };
#define N_EPS (sizeof(epss) / sizeof(epss[0]))

static char equal_family_systems[] = { 'A', 'D', 'N', 'V', 'W' };
#define N_EQUAL_FAMILY (sizeof(equal_family_systems) / sizeof(equal_family_systems[0]))

static char quad_arith_systems[] = { 'O', 'S', 'X', 'M', 'F' };
#define N_QUAD_ARITH (sizeof(quad_arith_systems) / sizeof(quad_arith_systems[0]))

static char great_circle_systems[] = { 'R', 'C', 'T', 'H', 'J' };
#define N_GREAT_CIRCLE (sizeof(great_circle_systems) / sizeof(great_circle_systems[0]))

static char iterative_systems[] = { 'P', 'K' };
#define N_ITERATIVE (sizeof(iterative_systems) / sizeof(iterative_systems[0]))

static char closed_form_misc_systems[] = { 'U', 'Y', 'L', 'Q' };
#define N_CLOSED_FORM_MISC (sizeof(closed_form_misc_systems) / sizeof(closed_form_misc_systems[0]))

static char sunshine_systems[] = { 'I', 'i' };
#define N_SUNSHINE (sizeof(sunshine_systems) / sizeof(sunshine_systems[0]))

/* Sun declinations spanning the year, one assigned per (armc, geolat) case (rotated) rather
 * than full cross-product, to bound case count to N_SUNSHINE * N_ARMC * N_GEOLAT = 60. */
static double sundecs[] = { -23.0, -10.0, 0.0, 10.0, 23.0 };
#define N_SUNDEC (sizeof(sundecs) / sizeof(sundecs[0]))

/* Polar-circle geolats, added to this task's battery only, to exercise the Placidus/Koch/
 * Gauquelin Porphyry fallback (|fi| >= 90-eps, eps=23.4392911 => cutoff ~66.56 deg). */
static double polar_geolats[] = { 51.5, -33.87, 0.0, 64.0, -64.0, 78.0, -78.0 };
#define N_POLAR_GEOLAT (sizeof(polar_geolats) / sizeof(polar_geolats[0]))

/* --- Houses 7/9 (Ephemeris UT wrappers + traditional sidereal) --- */

struct ut_triple { double tjd_ut, geolat, geolon; };

static struct ut_triple ut_triples[] = {
    { 2451545.0,  51.5,    -0.13  },  /* J2000.0, London */
    { 2460310.5, -33.87,  151.21  },  /* 2024-Jan-1, Sydney */
    { 2433282.5,   0.0,     0.0   },  /* 1950-Jan-1, equator */
    { 2488069.5,  40.71,  -74.01  },  /* 2100-Jan-1, New York */
    { 2378496.5,  64.0,    10.0   },  /* 1800-Jan-1, near-polar N */
    { 2305447.5, -64.0,   151.0   },  /* 1600-Jan-1, near-polar S */
};
#define N_UT_TRIPLE (sizeof(ut_triples) / sizeof(ut_triples[0]))

static char ut_wrapper_systems[] = { 'P', 'K', 'C', 'R', 'W', 'I' };
#define N_UT_WRAPPER_SYS (sizeof(ut_wrapper_systems) / sizeof(ut_wrapper_systems[0]))

static char sidereal_trad_systems[] = { 'P', 'W', 'E' };
#define N_SIDEREAL_TRAD_SYS (sizeof(sidereal_trad_systems) / sizeof(sidereal_trad_systems[0]))

int main(void) {
    int ia, ig, ie, is;
    int first;
    double cusp[40], cusp_speed[40], ascmc[10], ascmc_speed[10];
    char serr[256];

    printf("{\n");

    /* --- angles_special: the 8 ascmc special points + speeds, system-independent --- */
    printf("  \"angles_special\": [\n");
    first = 1;
    for (ia = 0; ia < N_ARMC; ia++) {
        for (ig = 0; ig < N_GEOLAT; ig++) {
            for (ie = 0; ie < N_EPS; ie++) {
                double armc = armcs[ia];
                double geolat = geolats[ig];
                double eps = epss[ie];
                int retc;

                memset(cusp, 0, sizeof(cusp));
                memset(cusp_speed, 0, sizeof(cusp_speed));
                memset(ascmc, 0, sizeof(ascmc));
                memset(ascmc_speed, 0, sizeof(ascmc_speed));
                serr[0] = '\0';

                retc = swe_houses_armc_ex2(armc, geolat, eps, 'P', cusp, ascmc,
                                            cusp_speed, ascmc_speed, serr);
                (void)retc;

                if (!first) printf(",\n");
                first = 0;
                printf("    {\"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                       "\"ascmc\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
                       "\"ascmc_speed\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                       armc, geolat, eps,
                       ascmc[0], ascmc[1], ascmc[2], ascmc[3], ascmc[4], ascmc[5], ascmc[6], ascmc[7],
                       ascmc_speed[0], ascmc_speed[1], ascmc_speed[2], ascmc_speed[3],
                       ascmc_speed[4], ascmc_speed[5], ascmc_speed[6], ascmc_speed[7]);
            }
        }
    }
    printf("\n  ],\n");

    /* --- equal_family: cusps[1..12] + speeds for A/D/N/V/W --- */
    printf("  \"equal_family\": [\n");
    first = 1;
    for (is = 0; is < N_EQUAL_FAMILY; is++) {
        char hsys = equal_family_systems[is];
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ig = 0; ig < N_GEOLAT; ig++) {
                for (ie = 0; ie < N_EPS; ie++) {
                    double armc = armcs[ia];
                    double geolat = geolats[ig];
                    double eps = epss[ie];
                    int retc, i;

                    memset(cusp, 0, sizeof(cusp));
                    memset(cusp_speed, 0, sizeof(cusp_speed));
                    memset(ascmc, 0, sizeof(ascmc));
                    memset(ascmc_speed, 0, sizeof(ascmc_speed));
                    serr[0] = '\0';

                    retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                                cusp_speed, ascmc_speed, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"cusps\": [", hsys, armc, geolat, eps);
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                    }
                    printf("], \"cusp_speed\": [");
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                    }
                    printf("]}");
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- quad_arith: cusps[1..12] + speeds for O/S/X/M/F --- */
    printf("  \"quad_arith\": [\n");
    first = 1;
    for (is = 0; is < N_QUAD_ARITH; is++) {
        char hsys = quad_arith_systems[is];
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ig = 0; ig < N_GEOLAT; ig++) {
                for (ie = 0; ie < N_EPS; ie++) {
                    double armc = armcs[ia];
                    double geolat = geolats[ig];
                    double eps = epss[ie];
                    int retc, i;

                    memset(cusp, 0, sizeof(cusp));
                    memset(cusp_speed, 0, sizeof(cusp_speed));
                    memset(ascmc, 0, sizeof(ascmc));
                    memset(ascmc_speed, 0, sizeof(ascmc_speed));
                    serr[0] = '\0';

                    retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                                cusp_speed, ascmc_speed, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"cusps\": [", hsys, armc, geolat, eps);
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                    }
                    printf("], \"cusp_speed\": [");
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                    }
                    printf("]}");
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- great_circle: cusps[1..12] + speeds for R/C/T/H/J --- */
    printf("  \"great_circle\": [\n");
    first = 1;
    for (is = 0; is < N_GREAT_CIRCLE; is++) {
        char hsys = great_circle_systems[is];
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ig = 0; ig < N_GEOLAT; ig++) {
                for (ie = 0; ie < N_EPS; ie++) {
                    double armc = armcs[ia];
                    double geolat = geolats[ig];
                    double eps = epss[ie];
                    int retc, i;

                    memset(cusp, 0, sizeof(cusp));
                    memset(cusp_speed, 0, sizeof(cusp_speed));
                    memset(ascmc, 0, sizeof(ascmc));
                    memset(ascmc_speed, 0, sizeof(ascmc_speed));
                    serr[0] = '\0';

                    retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                                cusp_speed, ascmc_speed, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"cusps\": [", hsys, armc, geolat, eps);
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                    }
                    printf("], \"cusp_speed\": [");
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                    }
                    printf("]}");
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- iterative: cusps[1..12] + speeds for P/K, incl. polar-circle geolats --- */
    printf("  \"iterative\": [\n");
    first = 1;
    for (is = 0; is < N_ITERATIVE; is++) {
        char hsys = iterative_systems[is];
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ig = 0; ig < N_POLAR_GEOLAT; ig++) {
                for (ie = 0; ie < N_EPS; ie++) {
                    double armc = armcs[ia];
                    double geolat = polar_geolats[ig];
                    double eps = epss[ie];
                    int retc, i;

                    memset(cusp, 0, sizeof(cusp));
                    memset(cusp_speed, 0, sizeof(cusp_speed));
                    memset(ascmc, 0, sizeof(ascmc));
                    memset(ascmc_speed, 0, sizeof(ascmc_speed));
                    serr[0] = '\0';

                    retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                                cusp_speed, ascmc_speed, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"cusps\": [", hsys, armc, geolat, eps);
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                    }
                    printf("], \"cusp_speed\": [");
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                    }
                    printf("]}");
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- gauquelin36: cusps[1..36] + speeds for G, incl. polar-circle geolats --- */
    printf("  \"gauquelin36\": [\n");
    first = 1;
    for (ia = 0; ia < N_ARMC; ia++) {
        for (ig = 0; ig < N_POLAR_GEOLAT; ig++) {
            for (ie = 0; ie < N_EPS; ie++) {
                double armc = armcs[ia];
                double geolat = polar_geolats[ig];
                double eps = epss[ie];
                int retc, i;

                memset(cusp, 0, sizeof(cusp));
                memset(cusp_speed, 0, sizeof(cusp_speed));
                memset(ascmc, 0, sizeof(ascmc));
                memset(ascmc_speed, 0, sizeof(ascmc_speed));
                serr[0] = '\0';

                retc = swe_houses_armc_ex2(armc, geolat, eps, 'G', cusp, ascmc,
                                            cusp_speed, ascmc_speed, serr);
                (void)retc;

                if (!first) printf(",\n");
                first = 0;
                printf("    {\"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                       "\"cusps\": [", armc, geolat, eps);
                for (i = 1; i <= 36; i++) {
                    printf("%.20e%s", cusp[i], (i < 36) ? ", " : "");
                }
                printf("], \"cusp_speed\": [");
                for (i = 1; i <= 36; i++) {
                    printf("%.20e%s", cusp_speed[i], (i < 36) ? ", " : "");
                }
                printf("]}");
            }
        }
    }
    printf("\n  ],\n");

    /* --- closed_form_misc: cusps[1..12] + speeds for U/Y/L/Q --- */
    printf("  \"closed_form_misc\": [\n");
    first = 1;
    for (is = 0; is < N_CLOSED_FORM_MISC; is++) {
        char hsys = closed_form_misc_systems[is];
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ig = 0; ig < N_GEOLAT; ig++) {
                for (ie = 0; ie < N_EPS; ie++) {
                    double armc = armcs[ia];
                    double geolat = geolats[ig];
                    double eps = epss[ie];
                    int retc, i;

                    memset(cusp, 0, sizeof(cusp));
                    memset(cusp_speed, 0, sizeof(cusp_speed));
                    memset(ascmc, 0, sizeof(ascmc));
                    memset(ascmc_speed, 0, sizeof(ascmc_speed));
                    serr[0] = '\0';

                    retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                                cusp_speed, ascmc_speed, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"cusps\": [", hsys, armc, geolat, eps);
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                    }
                    printf("], \"cusp_speed\": [");
                    for (i = 1; i <= 12; i++) {
                        printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                    }
                    printf("]}");
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- sunshine: cusps[1..12] + speeds for I/i, across a battery of Sun declinations.
     * ascmc[9] must be set BEFORE calling swe_houses_armc_ex2 -- that is how the ARMC-based
     * entry point receives Sun declination (c-ref-houses.md S11). */
    printf("  \"sunshine\": [\n");
    first = 1;
    for (is = 0; is < N_SUNSHINE; is++) {
        char hsys = sunshine_systems[is];
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ig = 0; ig < N_GEOLAT; ig++) {
                double armc = armcs[ia];
                double geolat = geolats[ig];
                double eps = epss[0];
                double sundec = sundecs[(ia * N_GEOLAT + ig) % N_SUNDEC];
                int retc, i;

                memset(cusp, 0, sizeof(cusp));
                memset(cusp_speed, 0, sizeof(cusp_speed));
                memset(ascmc, 0, sizeof(ascmc));
                memset(ascmc_speed, 0, sizeof(ascmc_speed));
                serr[0] = '\0';
                ascmc[9] = sundec;

                retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                            cusp_speed, ascmc_speed, serr);
                (void)retc;

                if (!first) printf(",\n");
                first = 0;
                printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                       "\"sundec\": %.20e, \"cusps\": [", hsys, armc, geolat, eps, sundec);
                for (i = 1; i <= 12; i++) {
                    printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                }
                printf("], \"cusp_speed\": [");
                for (i = 1; i <= 12; i++) {
                    printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                }
                printf("]}");
            }
        }
    }
    /* Circumpolar-Sun combinations (|tand(geolat)*tand(sundec)| >= 1): exercises Makransky's
     * sunshine_init ERR -> Porphyry fallback (swehouse.c:1175-1180, c-ref-houses.md S5 "I/i").
     * Treindl is included at the same combinations for contrast -- it never short-circuits on
     * this condition (sunshine_init's ERR is ignored), so 'I' should compute normally here. */
    {
        double polar_lats[] = { 70.0, -70.0 };
        double polar_decs[] = { 23.0, -23.0 };
        double polar_armcs[] = { 0.0, 215.7 };
        int ip, id, ipa;
        for (is = 0; is < N_SUNSHINE; is++) {
            char hsys = sunshine_systems[is];
            for (ipa = 0; ipa < 2; ipa++) {
                for (ip = 0; ip < 2; ip++) {
                    for (id = 0; id < 2; id++) {
                        double armc = polar_armcs[ipa];
                        double geolat = polar_lats[ip];
                        double eps = epss[0];
                        double sundec = polar_decs[id];
                        int retc, i;

                        memset(cusp, 0, sizeof(cusp));
                        memset(cusp_speed, 0, sizeof(cusp_speed));
                        memset(ascmc, 0, sizeof(ascmc));
                        memset(ascmc_speed, 0, sizeof(ascmc_speed));
                        serr[0] = '\0';
                        ascmc[9] = sundec;

                        retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc,
                                                    cusp_speed, ascmc_speed, serr);
                        (void)retc;

                        printf(",\n");
                        printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                               "\"sundec\": %.20e, \"cusps\": [", hsys, armc, geolat, eps, sundec);
                        for (i = 1; i <= 12; i++) {
                            printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                        }
                        printf("], \"cusp_speed\": [");
                        for (i = 1; i <= 12; i++) {
                            printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                        }
                        printf("]}");
                    }
                }
            }
        }
    }
    printf("\n  ]\n");

    /* --- ut_wrapper: swe_houses_ex2 (UT-based, date-aware entry point). Exercises
     * Ephemeris::houses_ex2's own ARMC/obliquity/nutation/sidtime setup, plus the
     * self-computed Sun declination for hsys='I'. One triple gets a NONUT variant. */
    printf(",\n  \"ut_wrapper\": [\n");
    first = 1;
    {
        int it, is2;
        for (it = 0; it < N_UT_TRIPLE; it++) {
            for (is2 = 0; is2 < N_UT_WRAPPER_SYS; is2++) {
                double tjd_ut = ut_triples[it].tjd_ut;
                double geolat = ut_triples[it].geolat;
                double geolon = ut_triples[it].geolon;
                char hsys = ut_wrapper_systems[is2];
                int retc, i;

                memset(cusp, 0, sizeof(cusp));
                memset(cusp_speed, 0, sizeof(cusp_speed));
                memset(ascmc, 0, sizeof(ascmc));
                memset(ascmc_speed, 0, sizeof(ascmc_speed));
                serr[0] = '\0';

                retc = swe_houses_ex2(tjd_ut, 0, geolat, geolon, hsys, cusp, ascmc,
                                       cusp_speed, ascmc_speed, serr);
                (void)retc;

                if (!first) printf(",\n");
                first = 0;
                printf("    {\"tjd_ut\": %.20e, \"geolat\": %.20e, \"geolon\": %.20e, "
                       "\"hsys\": \"%c\", \"nonut\": false, \"cusps\": [",
                       tjd_ut, geolat, geolon, hsys);
                for (i = 1; i <= 12; i++) {
                    printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                }
                printf("], \"cusp_speed\": [");
                for (i = 1; i <= 12; i++) {
                    printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                }
                printf("], \"ascmc\": [");
                for (i = 0; i < 8; i++) {
                    printf("%.20e%s", ascmc[i], (i < 7) ? ", " : "");
                }
                printf("], \"ascmc_speed\": [");
                for (i = 0; i < 8; i++) {
                    printf("%.20e%s", ascmc_speed[i], (i < 7) ? ", " : "");
                }
                printf("]}");
            }
        }

        /* NONUT variant: one triple (index 0), all systems in the subset. */
        for (is2 = 0; is2 < N_UT_WRAPPER_SYS; is2++) {
            double tjd_ut = ut_triples[0].tjd_ut;
            double geolat = ut_triples[0].geolat;
            double geolon = ut_triples[0].geolon;
            char hsys = ut_wrapper_systems[is2];
            int retc, i;

            memset(cusp, 0, sizeof(cusp));
            memset(cusp_speed, 0, sizeof(cusp_speed));
            memset(ascmc, 0, sizeof(ascmc));
            memset(ascmc_speed, 0, sizeof(ascmc_speed));
            serr[0] = '\0';

            retc = swe_houses_ex2(tjd_ut, SEFLG_NONUT, geolat, geolon, hsys, cusp, ascmc,
                                   cusp_speed, ascmc_speed, serr);
            (void)retc;

            printf(",\n");
            printf("    {\"tjd_ut\": %.20e, \"geolat\": %.20e, \"geolon\": %.20e, "
                   "\"hsys\": \"%c\", \"nonut\": true, \"cusps\": [",
                   tjd_ut, geolat, geolon, hsys);
            for (i = 1; i <= 12; i++) {
                printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
            }
            printf("], \"cusp_speed\": [");
            for (i = 1; i <= 12; i++) {
                printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
            }
            printf("], \"ascmc\": [");
            for (i = 0; i < 8; i++) {
                printf("%.20e%s", ascmc[i], (i < 7) ? ", " : "");
            }
            printf("], \"ascmc_speed\": [");
            for (i = 0; i < 8; i++) {
                printf("%.20e%s", ascmc_speed[i], (i < 7) ? ", " : "");
            }
            printf("]}");
        }
    }
    printf("\n  ]\n");

    /* --- sidereal_trad: swe_houses_ex2 with SEFLG_SIDEREAL, traditional (non ECL_T0/SSY_PLANE)
     * mode, Lahiri ayanamsa. Reuses the first 3 ut_triples. */
    swe_set_sid_mode(SE_SIDM_LAHIRI, 0, 0);
    printf(",\n  \"sidereal_trad\": [\n");
    first = 1;
    {
        int it, is2;
        for (it = 0; it < 3; it++) {
            for (is2 = 0; is2 < N_SIDEREAL_TRAD_SYS; is2++) {
                double tjd_ut = ut_triples[it].tjd_ut;
                double geolat = ut_triples[it].geolat;
                double geolon = ut_triples[it].geolon;
                char hsys = sidereal_trad_systems[is2];
                int retc, i;

                memset(cusp, 0, sizeof(cusp));
                memset(cusp_speed, 0, sizeof(cusp_speed));
                memset(ascmc, 0, sizeof(ascmc));
                memset(ascmc_speed, 0, sizeof(ascmc_speed));
                serr[0] = '\0';

                retc = swe_houses_ex2(tjd_ut, SEFLG_SIDEREAL, geolat, geolon, hsys, cusp, ascmc,
                                       cusp_speed, ascmc_speed, serr);
                (void)retc;

                if (!first) printf(",\n");
                first = 0;
                printf("    {\"tjd_ut\": %.20e, \"geolat\": %.20e, \"geolon\": %.20e, "
                       "\"hsys\": \"%c\", \"sid_mode\": %d, \"cusps\": [",
                       tjd_ut, geolat, geolon, hsys, SE_SIDM_LAHIRI);
                for (i = 1; i <= 12; i++) {
                    printf("%.20e%s", cusp[i], (i < 12) ? ", " : "");
                }
                printf("], \"cusp_speed\": [");
                for (i = 1; i <= 12; i++) {
                    printf("%.20e%s", cusp_speed[i], (i < 12) ? ", " : "");
                }
                printf("], \"ascmc\": [");
                for (i = 0; i < 8; i++) {
                    printf("%.20e%s", ascmc[i], (i < 7) ? ", " : "");
                }
                printf("], \"ascmc_speed\": [");
                for (i = 0; i < 8; i++) {
                    printf("%.20e%s", ascmc_speed[i], (i < 7) ? ", " : "");
                }
                printf("]}");
            }
        }
    }
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
