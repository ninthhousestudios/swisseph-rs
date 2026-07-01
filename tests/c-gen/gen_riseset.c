/*
 * Generates golden reference data for the rise/set/meridian-transit full algorithm:
 * swe_rise_trans_true_hor (SE_BIT_FORCE_SLOW_METHOD is included in rsmi for clarity, though it
 * has no effect on this function -- it only matters to the swe_rise_trans dispatcher, RSE 4).
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_riseset tests/c-gen/gen_riseset.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_riseset > tests/golden-data/riseset.json
 */

#include <stdio.h>
#include "swephexp.h"

struct geopos_t { double lon, lat, height; const char *name; };
static struct geopos_t geoposs[] = {
    { 8.55, 47.37, 500.0, "Zurich" },
    { 0.0, 0.0, 0.0, "Null Island" },
    { 18.95, 69.65, 10.0, "Tromso" },
};
#define N_GEOPOS (sizeof(geoposs) / sizeof(geoposs[0]))

static int bodies[] = { SE_SUN, SE_MOON };
static const char *body_names[] = { "Sun", "Moon" };
#define N_BODY (sizeof(bodies) / sizeof(bodies[0]))

static double tjd_uts[] = { 2451545.0, 2459000.5 };
#define N_TJD (sizeof(tjd_uts) / sizeof(tjd_uts[0]))

static int rsmis[] = { SE_CALC_RISE, SE_CALC_SET, SE_CALC_MTRANSIT };
static const char *rsmi_names[] = { "RISE", "SET", "MTRANSIT" };
#define N_RSMI (sizeof(rsmis) / sizeof(rsmis[0]))

int main(void) {
    int first;
    swe_set_ephe_path(NULL);

    printf("{\n");
    printf("  \"full\": [\n");
    first = 1;
    for (size_t ig = 0; ig < N_GEOPOS; ig++) {
        for (size_t ib = 0; ib < N_BODY; ib++) {
            for (size_t it = 0; it < N_TJD; it++) {
                for (size_t ir = 0; ir < N_RSMI; ir++) {
                    double geopos[3] = { geoposs[ig].lon, geoposs[ig].lat, geoposs[ig].height };
                    double tjd_ut = tjd_uts[it];
                    int32 ipl = bodies[ib];
                    int32 rsmi = rsmis[ir] | SE_BIT_FORCE_SLOW_METHOD;
                    double tret[10] = { 0 };
                    char serr[256] = { 0 };
                    int32 retval = swe_rise_trans_true_hor(
                        tjd_ut, ipl, NULL, SEFLG_MOSEPH, rsmi, geopos,
                        1013.25, 15.0, 0.0, tret, serr);
                    if (!first) printf(",\n");
                    first = 0;
                    printf("    {\"geopos\": [%.17g, %.17g, %.17g], \"geopos_name\": \"%s\", "
                           "\"body\": \"%s\", \"tjd_ut\": %.17g, \"rsmi\": \"%s\", "
                           "\"retval\": %d, \"tret0\": %.20e}",
                           geopos[0], geopos[1], geopos[2], geoposs[ig].name,
                           body_names[ib], tjd_ut, rsmi_names[ir], (int)retval, tret[0]);
                }
            }
        }
    }
    printf("\n  ]\n");
    printf("}\n");
    return 0;
}
