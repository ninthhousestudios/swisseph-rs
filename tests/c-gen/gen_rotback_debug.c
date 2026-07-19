#include <stdio.h>
#include <math.h>
#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"

static double epochs[] = {2451545.0, 2460310.5, 2433282.5, 2488069.5};
static int bodies[] = {SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS,
                        SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO};

int main(void) {
    double xx[6];
    char serr[256];
    double jd = 2378496.5;

    swe_set_ephe_path("../../ephe");

    /* Process earlier epochs (like gen_calc_sweph order) */
    for (int ie = 0; ie < 4; ie++) {
        for (int ib = 0; ib < 10; ib++) {
            swe_calc(epochs[ie], bodies[ib], SEFLG_SWIEPH | SEFLG_SPEED, xx, serr);
        }
    }

    /* Now process jd=2378496.5 and track which file each body uses */
    for (int ib = 0; ib < 10; ib++) {
        /* Track the file BEFORE the computation */
        const char *file_before = swed.fidat[SEI_FILE_PLANET].fptr ?
            swed.fidat[SEI_FILE_PLANET].fnam + strlen("../../../swisseph/ephe/") : "NULL";

        swe_calc(jd, bodies[ib], SEFLG_SWIEPH | SEFLG_SPEED, xx, serr);

        const char *file_after = swed.fidat[SEI_FILE_PLANET].fptr ?
            swed.fidat[SEI_FILE_PLANET].fnam + strlen("../../../swisseph/ephe/") : "NULL";

        printf("body%d (%d->%s): lon=%.15e lat=%.15e dist=%.15e\n",
               ib, 0, file_after, xx[0], xx[1], xx[2]);
    }

    /* Also check what file each planet's pldat came from */
    printf("\nPldat sources:\n");
    for (int i = 0; i < 10; i++) {
        printf("  pldat[%d].tfstart=%.5f\n", i, swed.pldat[i].tfstart);
    }

    swe_close();
    return 0;
}
