#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "swephexp.h"
#include "swephlib.h"

#define J2000 2451545.0

/* Model indices (SE_MODEL_* slots):
 * 0=delta_t, 1=prec_longterm, 2=prec_shortterm, 3=nutation,
 * 4=bias, 5=jplhor_mode, 6=jplhora_mode, 7=sidt
 */

static void set_prec_model(int longterm, int shortterm) {
    char buf[64];
    snprintf(buf, sizeof(buf), "0,%d,%d,0,0,0,0,0", longterm, shortterm);
    swe_set_astro_models(buf, 0);
}

/* Test epochs */
static double test_epochs[] = {
    2415020.0,       /* J1900 */
    2488070.0,       /* J2000 + 1 century */
    2414995.0,       /* J2000 - 1 century (approx) */
    2524595.0,       /* J2000 + 2 centuries */
    2378495.0,       /* J2000 - 2 centuries */
    2817045.0,       /* J2000 + 10 centuries */
    2086045.0,       /* J2000 - 10 centuries */
    4278045.0,       /* J2000 + 50 centuries */
    625045.0,        /* J2000 - 50 centuries */
    990544.5,        /* Owen boundary */
    3912544.5,       /* Owen boundary */
    2437684.5,       /* HORIZONS TJD0 */
    2378131.5,       /* JPLHOR range start */
    2525323.5,       /* JPLHOR range end */
    2396758.0,       /* Newcomb epoch */
};
static int n_epochs = sizeof(test_epochs) / sizeof(test_epochs[0]);

/* Precession model names and IDs */
struct model_info {
    const char *name;
    int id;
};

static struct model_info prec_models[] = {
    {"IAU1976", 1},
    {"Laskar1986", 2},
    {"WillEpsLask", 3},
    {"Williams1994", 4},
    {"Simon1994", 5},
    {"IAU2000", 6},
    {"Bretagnon2003", 7},
    {"IAU2006", 8},
    {"Vondrak2011", 9},
    {"Owen1990", 10},
    {"Newcomb", 11},
};
static int n_prec_models = sizeof(prec_models) / sizeof(prec_models[0]);

int main(void) {
    int i, m, d;
    int first = 1;
    double input[3] = {1.0, 0.1, 0.05};
    /* directions: C uses -1 = J2000->Date, +1 = Date->J2000 */
    int directions[] = {-1, 1};
    const char *dir_names[] = {"J2000ToDate", "DateToJ2000"};

    printf("{\n");
    printf("  \"precession\": [\n");

    /* All 11 models × both directions × all epochs */
    for (m = 0; m < n_prec_models; m++) {
        set_prec_model(prec_models[m].id, prec_models[m].id);
        for (d = 0; d < 2; d++) {
            for (i = 0; i < n_epochs; i++) {
                double jd = test_epochs[i];
                double x[3] = {input[0], input[1], input[2]};
                swi_precess(x, jd, 0, directions[d]);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"model\": \"%s\", \"direction\": \"%s\", "
                       "\"jd\": %.1f, \"flags\": 0, "
                       "\"input\": [%.20e, %.20e, %.20e], "
                       "\"output\": [%.20e, %.20e, %.20e]}",
                       prec_models[m].name, dir_names[d], jd,
                       input[0], input[1], input[2],
                       x[0], x[1], x[2]);
            }
        }
    }

    /* JPLHOR (DPSIDEPS_1980 = 262144) — in-range uses IAU1976, out-of-range uses Owen */
    {
        double jplhor_epochs[] = {
            2451545.0,   /* J2000 — in range */
            2415020.0,   /* J1900 — in range */
            2488070.0,   /* J2000 + 1 century — in range */
            2000000.0,   /* Far past — out of range, uses Owen */
            2600000.0,   /* Far future — out of range, uses Owen */
        };
        int n_jplhor = sizeof(jplhor_epochs) / sizeof(jplhor_epochs[0]);
        set_prec_model(9, 9);  /* Vondrak default (doesn't matter, JPLHOR overrides) */
        for (d = 0; d < 2; d++) {
            for (i = 0; i < n_jplhor; i++) {
                double jd = jplhor_epochs[i];
                double x[3] = {input[0], input[1], input[2]};
                swi_precess(x, jd, 262144, directions[d]);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"model\": \"JPLHOR\", \"direction\": \"%s\", "
                       "\"jd\": %.1f, \"flags\": 262144, "
                       "\"input\": [%.20e, %.20e, %.20e], "
                       "\"output\": [%.20e, %.20e, %.20e]}",
                       dir_names[d], jd,
                       input[0], input[1], input[2],
                       x[0], x[1], x[2]);
            }
        }
    }

    /* JPLHOR_APPROX (524288) with Vondrak default */
    {
        double jplhor_approx_epochs[] = {
            2451545.0,   /* J2000 */
            2437684.5,   /* HORIZONS TJD0 — at boundary */
            2415020.0,   /* J1900 — below boundary, is_jplhor=true */
            2460000.0,   /* Recent — above boundary, not is_jplhor */
            2488070.0,   /* J2000 + 1 century — above boundary */
            2000000.0,   /* Far past — below boundary, Owen path */
        };
        int n_approx = sizeof(jplhor_approx_epochs) / sizeof(jplhor_approx_epochs[0]);
        set_prec_model(9, 9);
        for (d = 0; d < 2; d++) {
            for (i = 0; i < n_approx; i++) {
                double jd = jplhor_approx_epochs[i];
                double x[3] = {input[0], input[1], input[2]};
                swi_precess(x, jd, 524288, directions[d]);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"model\": \"JPLHOR_APPROX\", \"direction\": \"%s\", "
                       "\"jd\": %.1f, \"flags\": 524288, "
                       "\"input\": [%.20e, %.20e, %.20e], "
                       "\"output\": [%.20e, %.20e, %.20e]}",
                       dir_names[d], jd,
                       input[0], input[1], input[2],
                       x[0], x[1], x[2]);
            }
        }
    }

    /* Roundtrip tests: precess J2000->Date then Date->J2000 should recover original */
    {
        double rt_epochs[] = {2415020.0, 2488070.0, 2817045.0, 990544.5};
        int rt_models[] = {1, 2, 8, 9, 10};
        const char *rt_names[] = {"IAU1976", "Laskar1986", "IAU2006", "Vondrak2011", "Owen1990"};
        int n_rt_epochs = sizeof(rt_epochs) / sizeof(rt_epochs[0]);
        int n_rt_models = sizeof(rt_models) / sizeof(rt_models[0]);
        int rm, ri;
        for (rm = 0; rm < n_rt_models; rm++) {
            set_prec_model(rt_models[rm], rt_models[rm]);
            for (ri = 0; ri < n_rt_epochs; ri++) {
                double jd = rt_epochs[ri];
                double x[3] = {input[0], input[1], input[2]};
                swi_precess(x, jd, 0, -1);   /* J2000 -> Date */
                swi_precess(x, jd, 0, 1);    /* Date -> J2000 */
                printf(",\n    {\"model\": \"roundtrip_%s\", \"direction\": \"roundtrip\", "
                       "\"jd\": %.1f, \"flags\": 0, "
                       "\"input\": [%.20e, %.20e, %.20e], "
                       "\"output\": [%.20e, %.20e, %.20e]}",
                       rt_names[rm], jd,
                       input[0], input[1], input[2],
                       x[0], x[1], x[2]);
            }
        }
    }

    printf("\n  ]\n");
    printf("}\n");

    return 0;
}
