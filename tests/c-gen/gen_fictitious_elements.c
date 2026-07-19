/*
 * Generates golden reference data for the fictitious-planet element layer:
 * resolved elements (after check_t_terms) and raw swi_osc_el_plan output.
 *
 * Because read_elements_file and check_t_terms are `static` in swemplan.c,
 * we include the source directly (same pattern as gen_heliacal_internals.c).
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../../swisseph -o gen_fictitious_elements \
 *      gen_fictitious_elements.c ../../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_fictitious_elements > ../golden-data/fictitious_elements.json
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

/* Include swemplan.c to get access to static functions. */
#include "../../../swisseph/swemplan.c"

static int g_first = 1;

static void comma(void) {
    if (!g_first) printf(",\n");
    g_first = 0;
}

int main(void) {
    char serr[256];
    swe_set_ephe_path("../../ephe");

    double epochs[] = {2451545.0, 2415020.0, 2460600.5, 2378506.5, 2500000.5, 1720010.5};
    int n_epochs = (int)(sizeof(epochs) / sizeof(epochs[0]));

    printf("[\n");

    for (int row = 0; row <= 18; row++) {
        for (int ie = 0; ie < n_epochs; ie++) {
            double tjd = epochs[ie];

            /* ── Element resolution ─────────────────────────────────── */
            double tjd0, tequ, mano, sema, ecce, parg, node_val, incl;
            char pname[256];
            int32 fict_ifl = 0;

            memset(pname, 0, sizeof(pname));
            int retc = read_elements_file(row, tjd, &tjd0, &tequ, &mano, &sema, &ecce,
                                          &parg, &node_val, &incl, pname, &fict_ifl, serr);
            if (retc == ERR) continue;

            /* ── Prime Earth/Sun barycentric ─────────────────────────── */
            double x_dummy[6];
            swe_calc(tjd, SE_SUN, SEFLG_SPEED, x_dummy, serr);

            double xearth[6], xsun[6];
            for (int i = 0; i < 6; i++) {
                xearth[i] = swed.pldat[SEI_EARTH].x[i];
                xsun[i]   = swed.pldat[SEI_SUNBARY].x[i];
            }

            /* ── swi_osc_el_plan ─────────────────────────────────────── */
            double xp[6] = {0, 0, 0, 0, 0, 0};
            retc = swi_osc_el_plan(tjd, xp, row, SEI_ANYBODY,
                                   xearth, xsun, serr);
            if (retc != OK) continue;

            comma();
            printf("{\"row\":%d,\"ipl\":%d,\"tjd\":%.1f,\"name\":\"%s\",\"is_geo\":%d,\n",
                   row, row + SE_FICT_OFFSET, tjd, pname, (fict_ifl & FICT_GEO) ? 1 : 0);

            /* Resolved elements — angles in radians as returned by read_elements_file */
            printf(" \"elem\":{\"tjd0\":%.20e,\"tequ\":%.20e,\"mano\":%.20e,\"sema\":%.20e,"
                   "\"ecce\":%.20e,\"parg\":%.20e,\"node\":%.20e,\"incl\":%.20e},\n",
                   tjd0, tequ, mano, sema, ecce, parg, node_val, incl);

            /* Earth/Sun barycentric state vectors used by osc_el_plan */
            printf(" \"xearth\":[%.20e,%.20e,%.20e,%.20e,%.20e,%.20e],\n",
                   xearth[0], xearth[1], xearth[2], xearth[3], xearth[4], xearth[5]);
            printf(" \"xsun\":[%.20e,%.20e,%.20e,%.20e,%.20e,%.20e],\n",
                   xsun[0], xsun[1], xsun[2], xsun[3], xsun[4], xsun[5]);

            /* Raw heliocentric output of swi_osc_el_plan */
            printf(" \"xp\":[%.20e,%.20e,%.20e,%.20e,%.20e,%.20e]}",
                   xp[0], xp[1], xp[2], xp[3], xp[4], xp[5]);
        }
    }

    printf("\n]\n");
    swe_close();
    return 0;
}
