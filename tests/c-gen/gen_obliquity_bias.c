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
    /* delta_t=default, prec_long, prec_short, rest=default */
    snprintf(buf, sizeof(buf), "0,%d,%d,0,0,0,0,0", longterm, shortterm);
    swe_set_astro_models(buf, 0);
}

static void set_bias_model(int bias_model) {
    char buf[64];
    snprintf(buf, sizeof(buf), "0,0,0,0,%d,0,0,0", bias_model);
    swe_set_astro_models(buf, 0);
}

/* Test epochs */
static double test_epochs[] = {
    2451545.0,       /* J2000.0 */
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
    -1000000.0,      /* Far past */
    5000000.0,       /* Far future */
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
    int i, m;

    printf("{\n");

    /* === OBLIQUITY TESTS === */
    printf("  \"obliquity\": [\n");
    int first = 1;
    for (m = 0; m < n_prec_models; m++) {
        /* Set both longterm and shortterm to the same model
         * so the model is used regardless of time range */
        set_prec_model(prec_models[m].id, prec_models[m].id);
        for (i = 0; i < n_epochs; i++) {
            double jd = test_epochs[i];
            double eps = swi_epsiln(jd, 0);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"model\": \"%s\", \"jd\": %.1f, \"eps\": %.20e}",
                   prec_models[m].name, jd, eps);
        }
    }

    /* Test shortterm/longterm fallback:
     * shortterm=IAU2006 (75 century range), longterm=Vondrak2011
     * At T=76 centuries, should fall through to Vondrak */
    set_prec_model(9, 8);  /* long=Vondrak, short=IAU2006 */
    {
        double jd_in_range = J2000 + 50.0 * 36525.0;   /* 50 centuries, in IAU2006 range */
        double jd_out_range = J2000 + 76.0 * 36525.0;  /* 76 centuries, outside IAU2006 range */
        double eps;
        eps = swi_epsiln(jd_in_range, 0);
        printf(",\n    {\"model\": \"shortterm_IAU2006\", \"jd\": %.1f, \"eps\": %.20e}", jd_in_range, eps);
        eps = swi_epsiln(jd_out_range, 0);
        printf(",\n    {\"model\": \"fallback_Vondrak\", \"jd\": %.1f, \"eps\": %.20e}", jd_out_range, eps);
    }

    /* JPLHOR_APPROX path with Vondrak + dcor_eps_jpl correction */
    set_prec_model(9, 9);  /* both Vondrak */
    {
        double eps;
        /* SEFLG_JPLHOR_APPROX = 512*1024 = 524288 */
        /* jplhora_mode defaults to V3 */
        eps = swi_epsiln(J2000, 524288);
        printf(",\n    {\"model\": \"Vondrak_JPLHOR_APPROX\", \"jd\": %.1f, \"eps\": %.20e}", J2000, eps);
        eps = swi_epsiln(2437846.5, 524288);
        printf(",\n    {\"model\": \"Vondrak_JPLHOR_APPROX\", \"jd\": %.1f, \"eps\": %.20e}", 2437846.5, eps);
        eps = swi_epsiln(2460000.0, 524288);
        printf(",\n    {\"model\": \"Vondrak_JPLHOR_APPROX\", \"jd\": %.1f, \"eps\": %.20e}", 2460000.0, eps);
    }

    /* JPLHOR (DPSIDEPS_1980) path */
    {
        double eps;
        /* SEFLG_JPLHOR = SEFLG_DPSIDEPS_1980 = 256*1024 = 262144 */
        eps = swi_epsiln(J2000, 262144);  /* in 1799-2202 range: uses IAU1976 */
        printf(",\n    {\"model\": \"JPLHOR_IAU1976\", \"jd\": %.1f, \"eps\": %.20e}", J2000, eps);
        eps = swi_epsiln(2000000.0, 262144);  /* outside range: uses Owen */
        printf(",\n    {\"model\": \"JPLHOR_Owen\", \"jd\": %.1f, \"eps\": %.20e}", 2000000.0, eps);
    }

    printf("\n  ],\n");

    /* === BIAS TESTS === */
    printf("  \"bias\": [\n");
    first = 1;

    /* Test vectors */
    double test_vecs[][6] = {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 3.0, 0.1, 0.2, 0.3},
    };
    int n_vecs = sizeof(test_vecs) / sizeof(test_vecs[0]);

    int bias_models[] = {2, 3};  /* IAU2000=2, IAU2006=3 */
    const char *bias_names[] = {"IAU2000", "IAU2006"};
    int directions[] = {0, 1};  /* 0=forward(GCRS->J2000), 1=backward(J2000->GCRS) */
    const char *dir_names[] = {"GcrsToJ2000", "J2000ToGcrs"};
    int bm, d, v;

    for (bm = 0; bm < 2; bm++) {
        set_bias_model(bias_models[bm]);
        for (d = 0; d < 2; d++) {
            for (v = 0; v < n_vecs; v++) {
                double x[6];
                for (i = 0; i < 6; i++) x[i] = test_vecs[v][i];
                /* iflag = SEFLG_SPEED (256) to include velocity */
                swi_bias(x, J2000, 256, directions[d]);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"bias_model\": \"%s\", \"direction\": \"%s\", "
                       "\"jd\": %.1f, \"flags\": %d, "
                       "\"input\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
                       "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                       bias_names[bm], dir_names[d], J2000, 256,
                       test_vecs[v][0], test_vecs[v][1], test_vecs[v][2],
                       test_vecs[v][3], test_vecs[v][4], test_vecs[v][5],
                       x[0], x[1], x[2], x[3], x[4], x[5]);
            }
        }
    }

    /* BiasModel::None should be identity */
    set_bias_model(1);  /* NONE=1 */
    {
        double x[6] = {1.0, 2.0, 3.0, 0.1, 0.2, 0.3};
        swi_bias(x, J2000, 256, 0);
        printf(",\n    {\"bias_model\": \"None\", \"direction\": \"GcrsToJ2000\", "
               "\"jd\": %.1f, \"flags\": %d, "
               "\"input\": [1.0, 2.0, 3.0, 0.1, 0.2, 0.3], "
               "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
               J2000, 256,
               x[0], x[1], x[2], x[3], x[4], x[5]);
    }

    /* === JPLHOR_APPROX BIAS TESTS === */
    /* SEFLG_JPLHOR_APPROX = 524288, SEFLG_SPEED = 256 */
    /* jplhora_mode defaults to V3 (SEMOD_JPLHORA_3) */
    {
        double jplhor_dates[] = {
            2451545.0,   /* J2000 - mid table */
            2437846.5,   /* Correction table start */
            2447000.0,   /* Mid table */
            2430000.0,   /* Before table (clamps to first entry) */
            2460000.0,   /* After table (clamps to last entry) */
        };
        int n_dates = sizeof(jplhor_dates) / sizeof(jplhor_dates[0]);
        double jplhor_vec[6] = {1.0, 2.0, 3.0, 0.1, 0.2, 0.3};
        int flags_combos[] = {524288 | 256, 524288};  /* with SPEED, without SPEED */
        int n_flags = 2;
        int fi, di, ji;

        set_bias_model(3);  /* IAU2006 */
        /* jplhora_mode = 0 means default (V3) */
        for (fi = 0; fi < n_flags; fi++) {
            for (di = 0; di < 2; di++) {
                for (ji = 0; ji < n_dates; ji++) {
                    double x[6];
                    for (i = 0; i < 6; i++) x[i] = jplhor_vec[i];
                    swi_bias(x, jplhor_dates[ji], flags_combos[fi], directions[di]);
                    printf(",\n    {\"bias_model\": \"IAU2006\", \"direction\": \"%s\", "
                           "\"jd\": %.1f, \"flags\": %d, "
                           "\"input\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
                           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                           dir_names[di], jplhor_dates[ji], flags_combos[fi],
                           jplhor_vec[0], jplhor_vec[1], jplhor_vec[2],
                           jplhor_vec[3], jplhor_vec[4], jplhor_vec[5],
                           x[0], x[1], x[2], x[3], x[4], x[5]);
                }
            }
        }
    }

    printf("\n  ]\n");
    printf("}\n");

    return 0;
}
