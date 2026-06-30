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
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
