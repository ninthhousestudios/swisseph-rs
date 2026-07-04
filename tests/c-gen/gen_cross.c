/*
 * Generates golden reference data for the crossings module (PNOC 7):
 *   swe_solcross / swe_solcross_ut        ("solcross" key, 18 cases)
 *   swe_mooncross / swe_mooncross_ut      ("mooncross" key, 18 cases)
 *   swe_mooncross_node / _ut              ("mooncross_node" key, 6 cases)
 *   swe_helio_cross / _ut                 ("helio_cross" key, 24 cases)
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_cross tests/c-gen/gen_cross.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_cross > tests/golden-data/crossings.json
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swephexp.h"

static int first_item = 1;

static void comma(void) {
    if (!first_item) printf(",");
    first_item = 0;
    printf("\n");
}

int main(void) {
    char serr[256];
    double x2cross_arr[] = { 0.0, 180.0, 359.5 };
    double jd_starts[] = { 2451500.0, 2440000.0, 2460600.0 };
    int n_x2cross = 3;
    int n_starts = 3;

    swe_set_ephe_path("../swisseph/ephe");

    printf("{\n");

    /* ---- solcross ---- */
    printf("\"solcross\": [");
    first_item = 1;
    for (int ix = 0; ix < n_x2cross; ix++) {
        for (int ij = 0; ij < n_starts; ij++) {
            double x2cross = x2cross_arr[ix];
            double jd_start = jd_starts[ij];

            /* ET variant */
            {
                double jd_result = swe_solcross(x2cross, jd_start, SEFLG_MOSEPH, serr);
                int ok = (jd_result >= jd_start) ? 1 : 0;
                comma();
                printf("  {\"x2cross\":%.1f,\"jd_start\":%.1f,\"variant\":\"et\","
                       "\"jd_result\":%.17g,\"ok\":%d}",
                       x2cross, jd_start, jd_result, ok);
            }
            /* UT variant */
            {
                double jd_result = swe_solcross_ut(x2cross, jd_start, SEFLG_MOSEPH, serr);
                int ok = (jd_result >= jd_start) ? 1 : 0;
                comma();
                printf("  {\"x2cross\":%.1f,\"jd_start\":%.1f,\"variant\":\"ut\","
                       "\"jd_result\":%.17g,\"ok\":%d}",
                       x2cross, jd_start, jd_result, ok);
            }
        }
    }
    printf("\n],\n");

    /* ---- mooncross ---- */
    printf("\"mooncross\": [");
    first_item = 1;
    for (int ix = 0; ix < n_x2cross; ix++) {
        for (int ij = 0; ij < n_starts; ij++) {
            double x2cross = x2cross_arr[ix];
            double jd_start = jd_starts[ij];

            /* ET */
            {
                double jd_result = swe_mooncross(x2cross, jd_start, SEFLG_MOSEPH, serr);
                int ok = (jd_result >= jd_start) ? 1 : 0;
                comma();
                printf("  {\"x2cross\":%.1f,\"jd_start\":%.1f,\"variant\":\"et\","
                       "\"jd_result\":%.17g,\"ok\":%d}",
                       x2cross, jd_start, jd_result, ok);
            }
            /* UT */
            {
                double jd_result = swe_mooncross_ut(x2cross, jd_start, SEFLG_MOSEPH, serr);
                int ok = (jd_result >= jd_start) ? 1 : 0;
                comma();
                printf("  {\"x2cross\":%.1f,\"jd_start\":%.1f,\"variant\":\"ut\","
                       "\"jd_result\":%.17g,\"ok\":%d}",
                       x2cross, jd_start, jd_result, ok);
            }
        }
    }
    printf("\n],\n");

    /* ---- mooncross_node ---- */
    printf("\"mooncross_node\": [");
    first_item = 1;
    double node_jds[] = { 2451545.0, 2440000.0, 2460600.0 };
    int n_node_jds = 3;
    for (int ij = 0; ij < n_node_jds; ij++) {
        double jd_start = node_jds[ij];

        /* ET */
        {
            double xlon = 0, xlat = 0;
            double jd_result = swe_mooncross_node(jd_start, SEFLG_MOSEPH, &xlon, &xlat, serr);
            int ok = (jd_result >= jd_start) ? 1 : 0;
            comma();
            printf("  {\"jd_start\":%.1f,\"variant\":\"et\","
                   "\"jd_result\":%.17g,\"xlon\":%.17g,\"xlat\":%.17g,\"ok\":%d}",
                   jd_start, jd_result, xlon, xlat, ok);
        }
        /* UT */
        {
            double xlon = 0, xlat = 0;
            double jd_result = swe_mooncross_node_ut(jd_start, SEFLG_MOSEPH, &xlon, &xlat, serr);
            int ok = (jd_result >= jd_start) ? 1 : 0;
            comma();
            printf("  {\"jd_start\":%.1f,\"variant\":\"ut\","
                   "\"jd_result\":%.17g,\"xlon\":%.17g,\"xlat\":%.17g,\"ok\":%d}",
                   jd_start, jd_result, xlon, xlat, ok);
        }
    }
    printf("\n],\n");

    /* ---- helio_cross ---- */
    printf("\"helio_cross\": [");
    first_item = 1;
    int helio_ipls[] = { SE_MERCURY, SE_MARS, SE_JUPITER, SE_EARTH };
    int n_helio_ipls = 4;
    double helio_x2cross[] = { 0.0, 120.5 };
    int n_helio_x2cross = 2;
    int dirs[] = { 1, -1 };
    int n_dirs = 2;
    double helio_jd_start = 2451545.0;

    for (int ip = 0; ip < n_helio_ipls; ip++) {
        for (int ix = 0; ix < n_helio_x2cross; ix++) {
            for (int id = 0; id < n_dirs; id++) {
                int ipl = helio_ipls[ip];
                double x2cross = helio_x2cross[ix];
                int dir = dirs[id];

                /* ET */
                {
                    double jd_cross = 0;
                    int rc = swe_helio_cross(ipl, x2cross, helio_jd_start,
                                             SEFLG_MOSEPH, dir, &jd_cross, serr);
                    comma();
                    printf("  {\"ipl\":%d,\"x2cross\":%.1f,\"jd_start\":%.1f,"
                           "\"dir\":%d,\"variant\":\"et\","
                           "\"jd_result\":%.17g,\"rc\":%d}",
                           ipl, x2cross, helio_jd_start, dir, jd_cross, rc);
                }
                /* UT */
                {
                    double jd_cross = 0;
                    int rc = swe_helio_cross_ut(ipl, x2cross, helio_jd_start,
                                                SEFLG_MOSEPH, dir, &jd_cross, serr);
                    comma();
                    printf("  {\"ipl\":%d,\"x2cross\":%.1f,\"jd_start\":%.1f,"
                           "\"dir\":%d,\"variant\":\"ut\","
                           "\"jd_result\":%.17g,\"rc\":%d}",
                           ipl, x2cross, helio_jd_start, dir, jd_cross, rc);
                }
            }
        }
    }
    printf("\n]\n");

    printf("}\n");
    swe_close();
    return 0;
}
