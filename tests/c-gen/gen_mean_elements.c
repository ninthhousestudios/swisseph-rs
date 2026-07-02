/*
 * Generates golden reference data for mean node, mean apogee, and ECL_NUT.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_mean_elements gen_mean_elements.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_mean_elements > ../golden-data/mean_elements.json
 */

#include <stdio.h>
#include "swephexp.h"

static int bodies[] = {
    SE_MEAN_NODE, SE_MEAN_APOG, SE_ECL_NUT
};
static const char *body_names[] = {
    "MeanNode", "MeanApogee", "EclipticNutation"
};
#define NBODIES 3

static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2460310.5,    /* 2024-Jan-1 */
    2433282.5,    /* 1950-Jan-1 */
    2488069.5,    /* 2100-Jan-1 */
    2378496.5,    /* 1800-Jan-1 */
    2305447.5,    /* 1600-Jan-1 */
    2159345.5,    /* 1200-Jan-1 */
    2013243.5,    /* 800-Jan-1 */
    1720693.5,    /* 0 CE */
    1538187.5,    /* -500 */
    990557.5,     /* -2000 */
};
#define NEPOCHS 11

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { 0,                "default" },
    { SEFLG_J2000,      "J2000" },
    { SEFLG_NONUT,      "NONUT" },
    { SEFLG_EQUATORIAL, "EQUATORIAL" },
    { SEFLG_XYZ,        "XYZ" },
};
#define NFLAGS 5

/* Sidereal cases apply to the two lunar elements only (not ECL_NUT). */
static int sid_bodies[] = { SE_MEAN_NODE, SE_MEAN_APOG };
static const char *sid_body_names[] = { "MeanNode", "MeanApogee" };
#define NSIDBODIES 2

struct sid_combo {
    int sid_mode;
    const char *name;
};
static struct sid_combo sid_combos[] = {
    { SE_SIDM_LAHIRI,                         "LAHIRI" },     /* traditional subtraction */
    { SE_SIDM_LAHIRI | SE_SIDBIT_ECL_T0,      "ECL_T0" },     /* rigorous, ecliptic of t0 */
    { SE_SIDM_LAHIRI | SE_SIDBIT_SSY_PLANE,   "SSY_PLANE" },  /* rigorous, solar-system plane */
};
#define NSID 3

static int first = 1;

static void emit(int body, const char *body_name, double jd, int flags,
                 const char *flag_name, int sid_mode, const double *xx) {
    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, \"body_name\": \"%s\", "
           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
           "\"sid_mode\": %d, "
           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
           body, body_name, jd, flags, flag_name, sid_mode,
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
}

int main(void) {
    char serr[256];
    swe_set_ephe_path(NULL);
    printf("[\n");
    /* Tropical cases (sid_mode = 0). */
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = SEFLG_MOSEPH | SEFLG_SPEED | flag_combos[ifl].flag;
                double xx[6];
                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "error: %s body=%d jd=%.1f flags=%d\n",
                            serr, bodies[ib], epochs[ie], flags);
                    return 1;
                }
                emit(bodies[ib], body_names[ib], epochs[ie], flags,
                     flag_combos[ifl].name, 0, xx);
            }
        }
    }
    /* Sidereal cases: SEFLG_SIDEREAL | SEFLG_SPEED, one per sid_mode. */
    for (int is = 0; is < NSID; is++) {
        swe_set_sid_mode(sid_combos[is].sid_mode, 0, 0);
        for (int ib = 0; ib < NSIDBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS; ie++) {
                int flags = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_SIDEREAL;
                double xx[6];
                int rc = swe_calc(epochs[ie], sid_bodies[ib], flags, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "error: %s body=%d jd=%.1f sid=%s\n",
                            serr, sid_bodies[ib], epochs[ie], sid_combos[is].name);
                    return 1;
                }
                emit(sid_bodies[ib], sid_body_names[ib], epochs[ie], flags,
                     sid_combos[is].name, sid_combos[is].sid_mode, xx);
            }
        }
    }
    printf("\n]\n");
    return 0;
}
