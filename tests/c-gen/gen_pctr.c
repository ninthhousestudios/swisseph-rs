/*
 * Generates golden reference data for swe_calc_pctr (planet-centric
 * calculation: position of ipl as seen from iplctr).
 *
 * Test matrix: 6 (ipl, iplctr) pairs x 3 epochs x 5 flag combos = 90 cases.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_pctr tests/c-gen/gen_pctr.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_pctr > tests/golden-data/pctr.json
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

struct body_pair {
    int ipl;
    int iplctr;
};

static struct body_pair pairs[] = {
    { SE_SUN,     SE_MARS },
    { SE_EARTH,   SE_MARS },
    { SE_JUPITER, SE_MARS },
    { SE_MOON,    SE_VENUS },
    { SE_SATURN,  SE_JUPITER },
    { SE_EARTH,   SE_MOON },
};
#define NPAIRS 6

static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2460600.0,   /* modern */
    2415020.5,   /* 1900-Jan-0.5 */
};
#define NEPOCHS 3

static int flag_combos[] = {
    SEFLG_MOSEPH | SEFLG_SPEED,
    SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL,
    SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_XYZ,
    SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_J2000,
    SEFLG_SWIEPH | SEFLG_SPEED,
};
#define NFLAGS 5

int main(void) {
    char serr[256];
    int count = 0;

    swe_set_ephe_path("../swisseph/ephe");

    printf("{\n");
    printf("\"pctr\": [");
    first_item = 1;

    for (int ip = 0; ip < NPAIRS; ip++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int ipl = pairs[ip].ipl;
                int iplctr = pairs[ip].iplctr;
                double tjd = epochs[ie];
                int iflag = flag_combos[ifl];
                double xx[6];
                memset(xx, 0, sizeof(xx));
                memset(serr, 0, sizeof(serr));

                int rc = swe_calc_pctr(tjd, ipl, iplctr, iflag, xx, serr);
                int ok = (rc == ERR) ? 0 : 1;

                comma();
                printf("  {\"ipl\":%d,\"iplctr\":%d,\"tjd\":%.17g,"
                       "\"iflag\":%d,\"retflag\":%d,"
                       "\"xx\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],"
                       "\"ok\":%d}",
                       ipl, iplctr, tjd, iflag, rc,
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5],
                       ok);
                count++;
            }
        }
    }

    printf("\n]\n");
    printf("}\n");

    fprintf(stderr, "gen_pctr: emitted %d cases\n", count);

    swe_close();
    return 0;
}
