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

static void set_sidt_model(int sidt_model) {
    char buf[64];
    swe_close();
    swe_set_ephe_path(NULL);
    snprintf(buf, sizeof(buf), "0,0,0,0,0,0,0,%d", sidt_model);
    swe_set_astro_models(buf, 0);
}

struct model_info {
    const char *name;
    int id;
};

static struct model_info sidt_models[] = {
    {"IAU1976", 1},
    {"IAU2006", 2},
    {"IersConv2010", 3},
    {"Longterm", 4},
};
static int n_sidt_models = sizeof(sidt_models) / sizeof(sidt_models[0]);

/* Test epochs spanning all time ranges */
static double test_epochs[] = {
    /* Deep past — long-term model territory */
    625000.5,           /* ~-3100 */
    990575.5,           /* ~-2100 */
    1356175.5,          /* ~-1100 */
    1757583.5,          /* ~0 CE */
    2086302.5,          /* ~1000 */
    2268923.5,          /* ~1500 */

    /* Near long-term boundary (1850) */
    2396758.0,          /* just before 1 Jan 1850 */
    2396758.5,          /* exactly SIDT_LTERM_T0 */
    2396759.0,          /* just after */

    /* Historical range */
    2415020.0,          /* J1900 */
    2433282.5,          /* B1950 */

    /* Modern era */
    2451545.0,          /* J2000.0 */
    2451545.5,          /* J2000 + 0.5 day */
    2451727.25,         /* ~2000.5 */
    2453371.5,          /* ~2005.0 */
    2455197.5,          /* ~2010.0 */
    2457023.5,          /* ~2015.0 */
    2458849.5,          /* ~2020.0 */
    2460000.0,          /* ~2023.13 */
    2460676.5,          /* ~2025.0 */

    /* Near long-term boundary (2050) */
    2469807.0,          /* just before 1 Jan 2050 */
    2469807.5,          /* exactly SIDT_LTERM_T1 */
    2469808.0,          /* just after */

    /* Far future — long-term model territory */
    2488070.0,          /* ~2100 */
    2524595.0,          /* ~2200 */
    2634370.0,          /* ~2500 */
    2816787.5,          /* ~3000 */

    /* Fractional days to test midnight splitting */
    2451545.25,         /* J2000 + 6h */
    2451545.75,         /* J2000 + 18h */
    2460000.3,          /* fractional */
    2460000.7,          /* fractional */
    2460000.9999,       /* near midnight */
};
static int n_epochs = sizeof(test_epochs) / sizeof(test_epochs[0]);

int main(void) {
    int i, m;
    int first = 1;

    printf("{\n");
    printf("  \"sidtime\": [\n");

    for (m = 0; m < n_sidt_models; m++) {
        set_sidt_model(sidt_models[m].id);
        for (i = 0; i < n_epochs; i++) {
            double tjd = test_epochs[i];
            double result = swe_sidtime(tjd);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"model\": \"%s\", \"tjd_ut\": %.20e, \"expected\": %.20e}",
                   sidt_models[m].name, tjd, result);
        }
    }

    printf("\n  ]\n");
    printf("}\n");

    return 0;
}
