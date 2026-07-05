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

static char quad_arith_systems[] = { 'O', 'S', 'X', 'M', 'F', 'B' };
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

/* --- Houses 9/9 (ECL_T0 / SSY_PLANE geometric sidereal projections) --- */

static char sidereal_geom_systems[] = { 'P', 'C', 'W' };
#define N_SIDEREAL_GEOM_SYS (sizeof(sidereal_geom_systems) / sizeof(sidereal_geom_systems[0]))

static int sidereal_geom_modes[] = {
    SE_SIDM_LAHIRI | SE_SIDBIT_ECL_T0,
    SE_SIDM_LAHIRI | SE_SIDBIT_SSY_PLANE,
};
#define N_SIDEREAL_GEOM_MODES (sizeof(sidereal_geom_modes) / sizeof(sidereal_geom_modes[0]))

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
    /* Alcabitius near-polar supplement: exercise the r-clamp boundary (|r| → 1). */
    {
        double polar_lats[] = { 78.0, -78.0 };
        int np = sizeof(polar_lats) / sizeof(polar_lats[0]);
        int ip;
        for (ia = 0; ia < N_ARMC; ia++) {
            for (ip = 0; ip < np; ip++) {
                for (ie = 0; ie < N_EPS; ie++) {
                    double armc = armcs[ia];
                    double geolat = polar_lats[ip];
                    double eps = epss[ie];
                    int retc, i;

                    memset(cusp, 0, sizeof(cusp));
                    memset(cusp_speed, 0, sizeof(cusp_speed));
                    memset(ascmc, 0, sizeof(ascmc));
                    memset(ascmc_speed, 0, sizeof(ascmc_speed));
                    serr[0] = '\0';

                    retc = swe_houses_armc_ex2(armc, geolat, eps, 'B', cusp, ascmc,
                                                cusp_speed, ascmc_speed, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"B\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"cusps\": [", armc, geolat, eps);
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
    printf("\n  ],\n");

    /* --- sidereal_geom: swe_houses_ex2 with SEFLG_SIDEREAL, ECL_T0/SSY_PLANE geometric
     * projection modes (swehouse.c:318-532), Lahiri t0/ayan_t0. Reuses the first 3 ut_triples. */
    printf("  \"sidereal_geom\": [\n");
    first = 1;
    {
        int it, is2, im;
        for (im = 0; im < N_SIDEREAL_GEOM_MODES; im++) {
            swe_set_sid_mode(sidereal_geom_modes[im], 0, 0);
            for (it = 0; it < 3; it++) {
                for (is2 = 0; is2 < N_SIDEREAL_GEOM_SYS; is2++) {
                    double tjd_ut = ut_triples[it].tjd_ut;
                    double geolat = ut_triples[it].geolat;
                    double geolon = ut_triples[it].geolon;
                    char hsys = sidereal_geom_systems[is2];
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
                           tjd_ut, geolat, geolon, hsys, sidereal_geom_modes[im]);
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
    }
    printf("\n  ],\n");

    /* --- house_pos: swe_house_pos (planet -> continuous house position) across all house
     * system chars, a couple of (armc, geolat, eps) triples (one temperate, one polar -- chosen
     * to exercise Koch's genuine circumpolar-failure branch, confirmed by a targeted scan: at
     * armc=105/geolat=67 two of the three xpin below succeed and one hits Koch's hpos=0
     * sentinel), and a few planet positions. The static sundec cache used internally by 'I'/'i'
     * (c-ref-houses.md S11) is primed via a swe_houses_armc_ex2 call with ascmc[9]=sundec set,
     * immediately before each swe_house_pos call -- harmless for non-Sunshine systems.
     *
     * "err" is driven by hpos==0.0, NOT by serr being non-empty: several systems (P/G's "Otto
     * Ludwig" circumpolar note, J/L/Q/default's "using simplified algorithm" note) set an
     * informational serr on a perfectly valid, non-zero hpos. Koch's hpos==0.0 is the one
     * genuine failure sentinel in this function -- no successful branch produces exactly 0.0. */
    printf("  \"house_pos\": [\n");
    first = 1;
    {
        static char hp_systems[] = {
            'A', 'B', 'C', 'D', 'F', 'G', 'H', 'I', 'i', 'J', 'K', 'L', 'M', 'N',
            'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y'
        };
        int n_hp_systems = sizeof(hp_systems) / sizeof(hp_systems[0]);
        struct hp_triple { double armc, geolat, eps; };
        struct hp_triple hp_triples[] = {
            { 47.5,  51.5, 23.4392911 },
            { 105.0, 67.0, 23.4392911 },
        };
        int n_hp_triples = sizeof(hp_triples) / sizeof(hp_triples[0]);
        double hp_xpins[][2] = { { 10.0, 0.0 }, { 123.4, 2.5 }, { 280.0, -1.0 } };
        int n_hp_xpins = sizeof(hp_xpins) / sizeof(hp_xpins[0]);
        double sundec = 10.0;
        int is3, it3, ix;

        for (is3 = 0; is3 < n_hp_systems; is3++) {
            char hsys = hp_systems[is3];
            for (it3 = 0; it3 < n_hp_triples; it3++) {
                double armc = hp_triples[it3].armc;
                double geolat = hp_triples[it3].geolat;
                double eps = hp_triples[it3].eps;
                for (ix = 0; ix < n_hp_xpins; ix++) {
                    double xpin[6];
                    double hpos;
                    int has_err;

                    xpin[0] = hp_xpins[ix][0];
                    xpin[1] = hp_xpins[ix][1];

                    serr[0] = '\0';
                    memset(cusp, 0, sizeof(cusp));
                    memset(ascmc, 0, sizeof(ascmc));
                    ascmc[9] = sundec;
                    swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc, NULL, NULL, serr);

                    serr[0] = '\0';
                    hpos = swe_house_pos(armc, geolat, eps, hsys, xpin, serr);
                    has_err = (hpos == 0.0);

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"hsys\": \"%c\", \"armc\": %.20e, \"geolat\": %.20e, \"eps\": %.20e, "
                           "\"xpin\": [%.20e, %.20e], \"sundec\": %.20e, \"hpos\": %.20e, "
                           "\"err\": %s}",
                           hsys, armc, geolat, eps, xpin[0], xpin[1], sundec, hpos,
                           has_err ? "true" : "false");
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- gauquelin_riseset: swe_gauquelin_sector, imeth in {0,1,2,3,4,5} via the full dispatcher.
     * Planet cases: 6 ut_triples x 3 bodies x 4 imeth(2-5) = 72 cases.
     * Star cases: 2 ut_triples x Aldebaran x 3 imeth(0,1,2) = 6 cases. */
    printf("  \"gauquelin_riseset\": [\n");
    first = 1;
    {
        int gq_bodies[] = { SE_SUN, SE_MOON, SE_MARS };
        int n_gq_bodies = sizeof(gq_bodies) / sizeof(gq_bodies[0]);
        int it5, ib5, im5;
        int gq_imeths[] = { 2, 3, 4, 5 };
        int n_gq_imeths = sizeof(gq_imeths) / sizeof(gq_imeths[0]);

        for (it5 = 0; it5 < N_UT_TRIPLE; it5++) {
            for (ib5 = 0; ib5 < n_gq_bodies; ib5++) {
                for (im5 = 0; im5 < n_gq_imeths; im5++) {
                    double t_ut = ut_triples[it5].tjd_ut;
                    double geopos[3];
                    double dgsect = 0;
                    int32 retc;

                    geopos[0] = ut_triples[it5].geolon;
                    geopos[1] = ut_triples[it5].geolat;
                    geopos[2] = 0.0;
                    serr[0] = '\0';

                    retc = swe_gauquelin_sector(t_ut, gq_bodies[ib5], NULL, SEFLG_MOSEPH,
                                                 gq_imeths[im5], geopos, 0, 0, &dgsect, serr);

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"tjd_ut\": %.20e, \"ipl\": %d, \"imeth\": %d, "
                           "\"geolon\": %.20e, \"geolat\": %.20e, "
                           "\"dgsect\": %.20e, \"retval\": %d}",
                           t_ut, gq_bodies[ib5], gq_imeths[im5], geopos[0], geopos[1],
                           dgsect, retc);
                }
            }
        }

        /* Fixed-star cases: Aldebaran at 2 mid-latitude triples x imeth {0, 1, 2}.
         * swe_fixstar strcpy()s into its starname argument, so use a local char[].
         * Set ephe_path here (not at top) so sefstars.txt is findable without
         * contaminating earlier Sunshine (I/i) cases that rely on Moshier fallback. */
        swe_set_ephe_path("../swisseph/ephe");
        {
            int star_triples[] = { 0, 1 };
            int star_imeths[] = { 0, 1, 2 };
            int n_star_triples = sizeof(star_triples) / sizeof(star_triples[0]);
            int n_star_imeths = sizeof(star_imeths) / sizeof(star_imeths[0]);
            int ist, ism;

            for (ist = 0; ist < n_star_triples; ist++) {
                for (ism = 0; ism < n_star_imeths; ism++) {
                    int idx = star_triples[ist];
                    double t_ut = ut_triples[idx].tjd_ut;
                    double geopos[3];
                    double dgsect = 0;
                    int32 retc;
                    char starname[256];

                    strcpy(starname, "Aldebaran");
                    geopos[0] = ut_triples[idx].geolon;
                    geopos[1] = ut_triples[idx].geolat;
                    geopos[2] = 0.0;
                    serr[0] = '\0';

                    retc = swe_gauquelin_sector(t_ut, SE_SUN, starname, SEFLG_MOSEPH,
                                                 star_imeths[ism], geopos, 0, 0, &dgsect, serr);

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"tjd_ut\": %.20e, \"ipl\": %d, \"imeth\": %d, "
                           "\"starname\": \"Aldebaran\", "
                           "\"geolon\": %.20e, \"geolat\": %.20e, "
                           "\"dgsect\": %.20e, \"retval\": %d}",
                           t_ut, SE_SUN, star_imeths[ism], geopos[0], geopos[1],
                           dgsect, retc);
                }
            }
        }
    }
    printf("\n  ],\n");

    /* --- gauquelin_sector: swe_gauquelin_sector, imeth in {0,1} (geometric, via house_pos 'G'). */
    printf("  \"gauquelin_sector\": [\n");
    first = 1;
    {
        int gq_bodies[] = { SE_SUN, SE_MOON, SE_MARS };
        int n_gq_bodies = sizeof(gq_bodies) / sizeof(gq_bodies[0]);
        int it4, ib4, im4;

        for (it4 = 0; it4 < N_UT_TRIPLE; it4++) {
            for (ib4 = 0; ib4 < n_gq_bodies; ib4++) {
                for (im4 = 0; im4 < 2; im4++) {
                    double t_ut = ut_triples[it4].tjd_ut;
                    double geopos[3];
                    double dgsect;
                    int32 retc;

                    geopos[0] = ut_triples[it4].geolon;
                    geopos[1] = ut_triples[it4].geolat;
                    geopos[2] = 0.0;
                    serr[0] = '\0';

                    retc = swe_gauquelin_sector(t_ut, gq_bodies[ib4], NULL, SEFLG_MOSEPH,
                                                 im4, geopos, 0, 0, &dgsect, serr);
                    (void)retc;

                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"tjd_ut\": %.20e, \"ipl\": %d, \"imeth\": %d, "
                           "\"geolon\": %.20e, \"geolat\": %.20e, \"dgsect\": %.20e}",
                           t_ut, gq_bodies[ib4], im4, geopos[0], geopos[1], dgsect);
                }
            }
        }
    }
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
