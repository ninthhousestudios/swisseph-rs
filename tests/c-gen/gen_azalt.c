/*
 * Generates golden reference data for atmospheric refraction and horizontal-coordinate
 * transforms: swe_refrac, swe_refrac_extended, swe_azalt, swe_azalt_rev.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_azalt tests/c-gen/gen_azalt.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_azalt > tests/golden-data/azalt.json
 */

#include <stdio.h>
#include "swephexp.h"

static double inalts[] = { -1.0, 0.0, 5.0, 17.9, 18.0, 45.0, 89.0 };
#define N_INALT (sizeof(inalts) / sizeof(inalts[0]))

static double atpresses[] = { 1013.25, 0.0 };
#define N_ATPRESS (sizeof(atpresses) / sizeof(atpresses[0]))

static double geoalts[] = { 0.0, 3000.0 };
#define N_GEOALT (sizeof(geoalts) / sizeof(geoalts[0]))

static int dirs[] = { SE_TRUE_TO_APP, SE_APP_TO_TRUE };
static const char *dir_names[] = { "TrueToApp", "AppToTrue" };
#define N_DIR 2

static double azalt_tjd_uts[] = { 2451545.0, 2459000.5, -1000000.0 };
/* -1000000.0 = far-past epoch (matches obliquity_bias.json's "Far past" epoch):
 * exercises the tidal-acceleration correction in the ARMC's sidereal-time
 * computation (adjust_for_tidacc scales with (year-1955)^2 -- see
 * azalt_armc_eps's forced-TIDAL_DEFAULT deltaT). At this epoch, an ARMC computed
 * with Moshier's default tid_acc (TIDAL_DE404) instead of the C-mandated
 * TIDAL_DEFAULT (TIDAL_DE431) diverges by ~4.5e-7 degrees -- enough to trip the
 * azalt/azalt_rev golden tests' 1e-7 tolerance. */
#define N_AZALT_TJD (sizeof(azalt_tjd_uts) / sizeof(azalt_tjd_uts[0]))

struct geopos_t { double lon, lat, height; };
static struct geopos_t geoposs[] = {
    { 8.55, 47.37, 500.0 },     /* Zurich */
    { -74.0, 40.7, 10.0 },      /* NYC */
};
#define N_GEOPOS (sizeof(geoposs) / sizeof(geoposs[0]))

static int azalt_dirs[] = { SE_ECL2HOR, SE_EQU2HOR };
static const char *azalt_dir_names[] = { "EclToHor", "EquToHor" };
#define N_AZALT_DIR 2

static int hor_dirs[] = { SE_HOR2ECL, SE_HOR2EQU };
static const char *hor_dir_names[] = { "HorToEcl", "HorToEqu" };
#define N_HOR_DIR 2

int main(void) {
    int first;
    swe_set_ephe_path(NULL);

    printf("{\n");

    /* === refrac === */
    printf("  \"refrac\": [\n");
    first = 1;
    for (size_t ii = 0; ii < N_INALT; ii++) {
        for (size_t ip = 0; ip < N_ATPRESS; ip++) {
            for (int id = 0; id < N_DIR; id++) {
                double inalt = inalts[ii];
                double atpress = atpresses[ip];
                double attemp = 15.0;
                double out = swe_refrac(inalt, atpress, attemp, dirs[id]);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"inalt\": %.17g, \"atpress\": %.17g, \"attemp\": %.17g, "
                       "\"dir\": \"%s\", \"out\": %.20e}",
                       inalt, atpress, attemp, dir_names[id], out);
            }
        }
    }
    printf("\n  ],\n");

    /* === refrac_ext === */
    printf("  \"refrac_ext\": [\n");
    first = 1;
    for (size_t ii = 0; ii < N_INALT; ii++) {
        for (size_t ip = 0; ip < N_ATPRESS; ip++) {
            for (size_t ig = 0; ig < N_GEOALT; ig++) {
                for (int id = 0; id < N_DIR; id++) {
                    double inalt = inalts[ii];
                    double atpress = atpresses[ip];
                    double attemp = 15.0;
                    double geoalt = geoalts[ig];
                    double lapse_rate = 0.0065;
                    double dret[4] = { 0, 0, 0, 0 };
                    double out = swe_refrac_extended(inalt, geoalt, atpress, attemp,
                                                       lapse_rate, dirs[id], dret);
                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"inalt\": %.17g, \"geoalt\": %.17g, \"atpress\": %.17g, "
                           "\"attemp\": %.17g, \"lapse_rate\": %.17g, \"dir\": \"%s\", "
                           "\"out\": %.20e, \"dret\": [%.20e, %.20e, %.20e, %.20e]}",
                           inalt, geoalt, atpress, attemp, lapse_rate, dir_names[id], out,
                           dret[0], dret[1], dret[2], dret[3]);
                }
            }
        }
    }
    printf("\n  ],\n");

    /* === azalt === */
    printf("  \"azalt\": [\n");
    first = 1;
    for (size_t it = 0; it < N_AZALT_TJD; it++) {
        for (size_t ig = 0; ig < N_GEOPOS; ig++) {
            for (int id = 0; id < N_AZALT_DIR; id++) {
                double tjd_ut = azalt_tjd_uts[it];
                double geopos[3] = { geoposs[ig].lon, geoposs[ig].lat, geoposs[ig].height };
                double xin[2] = { 120.0, 5.0 };
                double atpress = 0.0;
                double attemp = 15.0;
                double xaz[3] = { 0, 0, 0 };
                swe_azalt(tjd_ut, azalt_dirs[id], geopos, atpress, attemp, xin, xaz);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"tjd_ut\": %.17g, \"geopos\": [%.17g, %.17g, %.17g], "
                       "\"dir\": \"%s\", \"xin\": [%.17g, %.17g], "
                       "\"xaz\": [%.20e, %.20e, %.20e]}",
                       tjd_ut, geopos[0], geopos[1], geopos[2], azalt_dir_names[id],
                       xin[0], xin[1], xaz[0], xaz[1], xaz[2]);
            }
        }
    }
    printf("\n  ],\n");

    /* === azalt_rev === */
    printf("  \"azalt_rev\": [\n");
    first = 1;
    for (size_t it = 0; it < N_AZALT_TJD; it++) {
        for (size_t ig = 0; ig < N_GEOPOS; ig++) {
            for (int id = 0; id < N_HOR_DIR; id++) {
                double tjd_ut = azalt_tjd_uts[it];
                double geopos[3] = { geoposs[ig].lon, geoposs[ig].lat, geoposs[ig].height };
                double xin[2] = { 120.0, 30.0 };
                double xout[3] = { 0, 0, 0 };
                swe_azalt_rev(tjd_ut, hor_dirs[id], geopos, xin, xout);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"tjd_ut\": %.17g, \"geopos\": [%.17g, %.17g, %.17g], "
                       "\"dir\": \"%s\", \"xin\": [%.17g, %.17g], "
                       "\"xout\": [%.20e, %.20e]}",
                       tjd_ut, geopos[0], geopos[1], geopos[2], hor_dir_names[id],
                       xin[0], xin[1], xout[0], xout[1]);
            }
        }
    }
    printf("\n  ]\n");

    printf("}\n");
    return 0;
}
