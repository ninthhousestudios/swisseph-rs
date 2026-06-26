/*
 * Generates golden reference data for the Moshier backend wrapper layer.
 * Calls swe_calc with flags that disable all corrections, requesting
 * equatorial Cartesian J2000 output — matching our backend::compute().
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_moshier_backend gen_moshier_backend.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_moshier_backend > ../golden-data/moshier_backend.json
 */

#include <stdio.h>
#include "swephexp.h"

#define FLAGS (SEFLG_MOSEPH | SEFLG_J2000 | SEFLG_NONUT | SEFLG_TRUEPOS \
             | SEFLG_EQUATORIAL | SEFLG_SPEED | SEFLG_NOGDEFL \
             | SEFLG_NOABERR | SEFLG_XYZ | SEFLG_ICRS)

static int bodies[] = {
    SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS,
    SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO
};
static const char *body_names[] = {
    "sun", "moon", "mercury", "venus", "mars",
    "jupiter", "saturn", "uranus", "neptune", "pluto"
};
#define NBODIES 10

static double epochs[] = {
    2451545.0,                        /* J2000 */
    2451545.0 + 100.0 * 365.25,      /* +100yr */
    2451545.0 - 100.0 * 365.25,      /* -100yr */
    2451545.0 + 500.0 * 365.25,      /* +500yr */
    2451545.0 - 500.0 * 365.25,      /* -500yr */
    2451545.0 + 1000.0 * 365.25,     /* +1000yr */
    2451545.0 - 1000.0 * 365.25,     /* -1000yr */
    2451545.0 + 0.5,                  /* half day offset */
    2451545.0 + 27.3,                 /* ~1 lunar month */
    2451545.0 + 365.25,               /* +1yr */
    2451545.0 - 365.25,               /* -1yr */
};
#define NEPOCHS 11

int main(void) {
    char serr[256];
    swe_set_ephe_path(NULL);
    printf("{\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        printf("  \"%s\": [\n", body_names[ib]);
        for (int ie = 0; ie < NEPOCHS; ie++) {
            double xx[6];
            int rc = swe_calc(epochs[ie], bodies[ib], FLAGS, xx, serr);
            if (rc < 0) {
                fprintf(stderr, "error %s body=%d jd=%.1f\n",
                        serr, bodies[ib], epochs[ie]);
                return 1;
            }
            printf("    {\"jd\": %.20e, \"x\": %.20e, \"y\": %.20e, "
                   "\"z\": %.20e, \"vx\": %.20e, \"vy\": %.20e, "
                   "\"vz\": %.20e}",
                   epochs[ie], xx[0], xx[1], xx[2],
                   xx[3], xx[4], xx[5]);
            if (ie < NEPOCHS - 1) printf(",");
            printf("\n");
        }
        printf("  ]");
        if (ib < NBODIES - 1) printf(",");
        printf("\n");
    }
    printf("}\n");
    return 0;
}
