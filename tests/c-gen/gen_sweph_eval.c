/*
 * Generates golden reference data for SE1 segment evaluation.
 * Forces segment loading via swe_calc, then evaluates Chebyshev polynomials
 * directly from the internal segp array to get raw SE1 backend positions.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_sweph_eval gen_sweph_eval.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_sweph_eval > ../golden-data/sweph_eval.json
 */

#include <stdio.h>
#include <math.h>
#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"

/* Internal body IDs matching the SE1 file body_id values */
static int body_ids[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
/* SE_* API body to call swe_calc with, to trigger loading each internal body */
static int se_bodies[] = {SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS,
                          SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO};
#define NBODIES 10

static double epochs[] = {
    2451545.0,                        /* J2000 */
    2451545.0 + 0.5,                  /* half day */
    2451545.0 + 27.3,                 /* ~1 lunar month */
    2451545.0 + 365.25,               /* +1yr */
    2451545.0 - 365.25,               /* -1yr */
    2451545.0 + 100.0 * 365.25,       /* +100yr */
    2451545.0 - 100.0 * 365.25,       /* -100yr */
    2451545.0 + 500.0 * 365.25,       /* +500yr */
    2378496.5,                        /* 1800-Jan-1 */
};
#define NEPOCHS 9

int main(void) {
    double xx[6];
    char serr[256];
    int first = 1;

    printf("{\"cases\": [\n");

    for (int ie = 0; ie < NEPOCHS; ie++) {
        double jd = epochs[ie];
        for (int ib = 0; ib < NBODIES; ib++) {
            int body_id = body_ids[ib];
            int se_body = se_bodies[ib];
            /* Reset C library state before each call so file caching does not
             * carry over between test cases. This gives deterministic, stateless
             * golden data matching our stateless Rust implementation. */
            swe_close();
            swe_set_ephe_path("../../ephe");
            /* Force file loading and segment caching */
            int ret = swe_calc(jd, se_body, SEFLG_SWIEPH | SEFLG_SPEED, xx, serr);
            if (ret < 0) {
                fprintf(stderr, "swe_calc failed for body %d at JD %.1f: %s\n",
                        se_body, jd, serr);
                continue;
            }

            /* Access internal plan_data for this body */
            struct plan_data *pdp = &swed.pldat[body_id];
            if (pdp->segp == NULL) {
                fprintf(stderr, "segp is NULL for body %d at JD %.1f\n", body_id, jd);
                continue;
            }

            /* Evaluate Chebyshev from raw coefficients */
            double t = (jd - pdp->tseg0) / pdp->dseg * 2.0 - 1.0;
            double pos[6] = {0};
            for (int k = 0; k < 3; k++) {
                pos[k] = swi_echeb(t, pdp->segp + k * pdp->ncoe, pdp->neval);
                pos[k+3] = swi_edcheb(t, pdp->segp + k * pdp->ncoe, pdp->neval)
                           / pdp->dseg * 2.0;
            }

            if (!first) printf(",\n");
            first = 0;
            printf("  {\"body_id\": %d, \"jd\": %.20e, "
                   "\"x\": %.20e, \"y\": %.20e, \"z\": %.20e, "
                   "\"vx\": %.20e, \"vy\": %.20e, \"vz\": %.20e, "
                   "\"neval\": %d}",
                   body_id, jd,
                   pos[0], pos[1], pos[2],
                   pos[3], pos[4], pos[5],
                   pdp->neval);
        }
    }

    printf("\n]}\n");
    swe_close();
    return 0;
}
