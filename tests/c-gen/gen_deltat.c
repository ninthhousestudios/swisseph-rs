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

static void set_deltat_model(int dt_model) {
    char buf[64];
    swe_close();
    swe_set_ephe_path(NULL);
    snprintf(buf, sizeof(buf), "%d,0,0,0,0,0,0,0", dt_model);
    swe_set_astro_models(buf, 0);
}

struct model_info {
    const char *name;
    int id;
};

static struct model_info dt_models[] = {
    {"StephensonMorrison1984", 1},
    {"Stephenson1997", 2},
    {"StephensonMorrison2004", 3},
    {"EspenakMeeus2006", 4},
    {"StephensonEtc2016", 5},
};
static int n_dt_models = sizeof(dt_models) / sizeof(dt_models[0]);

/* Test epochs spanning all year ranges */
static double test_epochs[] = {
    /* Deep past (before -1000) */
    625000.5,           /* ~-3100 (Moshier start) */
    990575.5,           /* ~-2100 */
    1356175.5,          /* ~-1100 */

    /* Historical range (-1000 to -500) */
    1356200.5,          /* just past -1100 */
    1500000.5,          /* ~-707 */
    1538395.5,          /* ~-602 */
    1541700.5,          /* ~-593 */
    1574945.5,          /* ~-502 */

    /* Historical range (-500 to 1600) */
    1575000.5,          /* ~-500 */
    1757583.5,          /* ~0 */
    2086302.5,          /* ~1000 */
    2159557.5,          /* ~1200.5 */
    2268923.5,          /* ~1500 */
    2299160.5,          /* Gregorian reform */

    /* Blend zone 1600-1620 */
    2305447.5,          /* ~1600 */
    2310000.5,          /* ~1612 */
    2312447.5,          /* ~1619 */

    /* Tabulated range: boundaries + interior */
    2312814.5,          /* 1620.0 (TABSTART) */
    2312815.5,          /* 1620 + 1 day */
    2313179.5,          /* ~1621.0 */
    2313544.5,          /* ~1622.0 */
    2313909.5,          /* ~1623.0 */

    /* Tabulated range: interior fractional years */
    2415020.0,          /* ~1899.998 (J1900) */
    2433282.5,          /* ~1950.0 (B1950) */
    2451545.0,          /* 2000.0 (J2000) */
    2451545.5,          /* J2000 + 0.5 day */
    2451727.25,         /* ~2000.5 */
    2460000.0,          /* ~2023.13 */

    /* Near tabulated end */
    2461788.5,          /* ~2028.0 (near TABEND) */

    /* Stephenson 2016 blend zone around JD 2435108.5 */
    2434108.5,          /* blend start */
    2434600.5,          /* mid-blend */
    2435100.5,          /* near blend end */
    2435108.5,          /* blend end boundary */
    2435110.5,          /* just past blend */

    /* Future extrapolation */
    2470000.5,          /* ~2050 */
    2488070.0,          /* ~2100 */
    2524595.0,          /* ~2200 */
    2561120.0,          /* ~2300 */
    2634370.0,          /* ~2500 */
    2670895.0,          /* ~2600 */

    /* Espenak-Meeus boundary */
    2317746.13090277789, /* exact EM2006 boundary */
    2317740.5,          /* just before */
    2317750.5,          /* just after */
};
static int n_epochs = sizeof(test_epochs) / sizeof(test_epochs[0]);

int main(void) {
    int i, m;
    int first = 1;
    char serr[256];
    double dt;

    printf("{\n");
    printf("  \"deltat\": [\n");

    for (m = 0; m < n_dt_models; m++) {
        set_deltat_model(dt_models[m].id);
        for (i = 0; i < n_epochs; i++) {
            double tjd = test_epochs[i];
            dt = swe_deltat_ex(tjd, SEFLG_MOSEPH, serr);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"model\": \"%s\", \"tjd\": %.20e, \"expected\": %.20e}",
                   dt_models[m].name, tjd, dt);
        }
    }

    printf("\n  ]\n");
    printf("}\n");

    return 0;
}
