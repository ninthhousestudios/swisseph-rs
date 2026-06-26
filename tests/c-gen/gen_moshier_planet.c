#include <stdio.h>
#include <math.h>
#include "swephexp.h"
#include "sweph.h"

/*
 * Calls swi_moshplan2() directly for all 9 planets at 9 test epochs.
 * Planet index maps to planets[] array: 0=mer, 1=ven, 2=ear, 3=mar,
 * 4=jup, 5=sat, 6=ura, 7=nep, 8=plu.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_moshier_planet gen_moshier_planet.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_moshier_planet > ../golden-data/moshier_planet.json
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
};
#define NEPOCHS 9

static const char *planet_names[] = {
    "mercury", "venus", "earth", "mars", "jupiter",
    "saturn", "uranus", "neptune", "pluto"
};
#define NPLANETS 9

int main(void) {
    printf("{\n");
    for (int ip = 0; ip < NPLANETS; ip++) {
        printf("  \"%s\": [\n", planet_names[ip]);
        for (int ie = 0; ie < NEPOCHS; ie++) {
            double pobj[3];
            swi_moshplan2(epochs[ie], ip, pobj);
            printf("    {\"jd\": %.20e, \"lon\": %.20e, \"lat\": %.20e, \"dist\": %.20e}",
                   epochs[ie], pobj[0], pobj[1], pobj[2]);
            if (ie < NEPOCHS - 1) printf(",");
            printf("\n");
        }
        printf("  ]");
        if (ip < NPLANETS - 1) printf(",");
        printf("\n");
    }
    printf("}\n");
    return 0;
}
