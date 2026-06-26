#include <stdio.h>
#include <math.h>
#include "swephexp.h"
#include "sweph.h"

/*
 * Calls swi_moshmoon2() directly for the Moon at 11 test epochs.
 * pol[3] is filled with [lon, lat, dist] in radians/AU.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_moshier_moon gen_moshier_moon.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_moshier_moon > ../golden-data/moshier_moon.json
 */

static double epochs[] = {
    2451545.0 - 3000.0 * 365.25,
    2451545.0 - 1000.0 * 365.25,
    2451545.0 - 500.0 * 365.25,
    2451545.0 - 100.0 * 365.25,
    2451545.0,
    2451545.0 + 100.0 * 365.25,
    2451545.0 + 500.0 * 365.25,
    2451545.0 + 1000.0 * 365.25,
    2451545.0 + 999.0 * 365.25,
    2451545.0 + 0.5,
    2451545.0 + 27.3,
};
#define NEPOCHS 11

int main(void) {
    printf("{\"moon\": [\n");
    for (int ie = 0; ie < NEPOCHS; ie++) {
        double pol[3];
        swi_moshmoon2(epochs[ie], pol);
        printf("  {\"jd\": %.20e, \"lon\": %.20e, \"lat\": %.20e, \"dist\": %.20e}",
               epochs[ie], pol[0], pol[1], pol[2]);
        if (ie < NEPOCHS - 1) printf(",");
        printf("\n");
    }
    printf("]}\n");
    return 0;
}
