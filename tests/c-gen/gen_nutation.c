#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "swephexp.h"
#include "swephlib.h"
#include "sweph.h"

/* Model indices for swe_set_astro_models string:
 * 0=delta_t, 1=prec_longterm, 2=prec_shortterm, 3=nutation,
 * 4=bias, 5=jplhor_mode, 6=jplhora_mode, 7=sidt
 */

static void set_nutation_model(int nut_model) {
    char buf[64];
    swe_close();
    swe_set_ephe_path(NULL);
    snprintf(buf, sizeof(buf), "0,0,0,%d,0,0,0,0", nut_model);
    swe_set_astro_models(buf, 0);
    swed.do_interpolate_nut = FALSE;
}

static double test_epochs[] = {
    2451545.0,       /* J2000 */
    2415020.0,       /* J1900 */
    2488070.0,       /* J2000 + 1 century */
    2414995.0,       /* ~J2000 - 1 century */
    2524595.0,       /* J2000 + 2 centuries */
    2378495.0,       /* J2000 - 2 centuries */
    2451545.5,       /* J2000 + 0.5 day */
    2451544.5,       /* J2000 - 0.5 day */
    2460000.0,       /* recent epoch */
    2440000.0,       /* ~1968 */
    2430000.0,       /* ~1941 */
    2500000.0,       /* ~2132 */
    2550000.0,       /* ~2269 */
    2400000.0,       /* ~1858 */
    2437684.5,       /* HORIZONS TJD0 */
    2451545.25,      /* J2000 + 0.25 day */
};
static int n_epochs = sizeof(test_epochs) / sizeof(test_epochs[0]);

struct model_info {
    const char *name;
    int id;
};

static struct model_info nut_models[] = {
    {"IAU1980", 1},
    {"IAUCorr1987", 2},
    {"IAU2000A", 3},
    {"IAU2000B", 4},
    {"Woolard", 5},
};
static int n_nut_models = sizeof(nut_models) / sizeof(nut_models[0]);

int main(void) {
    int i, m;
    int first = 1;
    double nutlo[2];

    printf("{\n");
    printf("  \"nutation\": [\n");

    /* All 5 models × all epochs, flags = 0 (normal dispatch) */
    for (m = 0; m < n_nut_models; m++) {
        set_nutation_model(nut_models[m].id);
        for (i = 0; i < n_epochs; i++) {
            double jd = test_epochs[i];
            swi_nutation(jd, 0, nutlo);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"model\": \"%s\", \"jd\": %.6f, \"flags\": 0, "
                   "\"dpsi\": %.20e, \"deps\": %.20e}",
                   nut_models[m].name, jd,
                   nutlo[0], nutlo[1]);
        }
    }

    /* JPLHOR and JPLHOR_APPROX paths skipped — they crash without EOP data files.
     * The JPLHOR path is IAU1980 + EOP corrections (stub in Rust).
     * The JPLHOR_APPROX paths are model output + constant offsets,
     * verified by Rust-only unit tests. */

    printf("\n  ]\n");
    printf("}\n");

    return 0;
}
