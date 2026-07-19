/*
 * Generates golden reference data for the osculating lunar node (SE_TRUE_NODE)
 * and osculating apogee / "true Lilith" (SE_OSCU_APOG) through swe_calc
 * (lunar_osc_elem + swi_plan_for_osc_elem), swisseph-rs/84.
 *
 * Bodies: SE_TRUE_NODE (11), SE_OSCU_APOG (13).
 * Backends: SEFLG_MOSEPH, SEFLG_SWIEPH.
 * Flags:   SPEED, SPEED|EQUATORIAL, SPEED|XYZ, SPEED|NONUT, SPEED|J2000, 0.
 * Epochs:  the 7 gen_calc.c epochs.
 * => 2 bodies * 2 backends * 6 flags * 7 epochs = 168 tropical cases.
 *
 * Sidereal (SEFLG_SIDEREAL|SEFLG_SPEED) across three sid_modes — Lahiri
 * (traditional ayanamsa subtraction), Lahiri|ECL_T0 and Lahiri|SSY_PLANE (the
 * two "rigorous" projections that read the J2000 equatorial vector) — over the
 * same bodies/backends/epochs: 3 * 2 * 2 * 7 = 84 sidereal cases (swisseph-rs/84
 * review follow-up). 252 cases total.
 *
 * Records xx[0..6) + the swe_calc return flag + sid_mode (0 = tropical).
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_truenode tests/c-gen/gen_truenode.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_truenode > tests/golden-data/truenode.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static int bodies[] = {SE_TRUE_NODE, SE_OSCU_APOG};
static const char *body_names[] = {"TrueNode", "OscuApogee"};
#define NBODIES 2

static double epochs[] = {
    2451545.0, /* J2000.0 */
    2460310.5, /* 2024-Jan-1 */
    2433282.5, /* 1950-Jan-1 */
    2488069.5, /* 2100-Jan-1 */
    2378496.5, /* 1800-Jan-1 */
    2305447.5, /* 1600-Jan-1 */
    2159345.5, /* 1200-Jan-1 */
};
#define NEPOCHS 7

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    {SEFLG_SPEED, "SPEED"},
    {SEFLG_SPEED | SEFLG_EQUATORIAL, "EQUATORIAL"},
    {SEFLG_SPEED | SEFLG_XYZ, "XYZ"},
    {SEFLG_SPEED | SEFLG_NONUT, "NONUT"},
    {SEFLG_SPEED | SEFLG_J2000, "J2000"},
    {0, "no_speed"},
};
#define NFLAGS 6

struct eph_combo {
    int flag;
    const char *name;
};

static struct eph_combo eph_combos[] = {
    {SEFLG_MOSEPH, "MOSEPH"},
    {SEFLG_SWIEPH, "SWIEPH"},
};
#define NEPH 2

/* Sidereal mode dimension: sid_mode passed to swe_set_sid_mode, plus a label. */
struct sid_combo {
    int sid_mode;
    const char *name;
};

static struct sid_combo sid_combos[] = {
    {SE_SIDM_LAHIRI, "LAHIRI"},                          /* traditional subtraction */
    {SE_SIDM_LAHIRI | SE_SIDBIT_ECL_T0, "ECL_T0"},       /* rigorous, ecliptic of t0 */
    {SE_SIDM_LAHIRI | SE_SIDBIT_SSY_PLANE, "SSY_PLANE"}, /* rigorous, solar-system plane */
};
#define NSID 3

static int first = 1;

static void emit(int body, const char *body_name, double jd, int flags,
                 const char *flag_name, const char *eph_name, int sid_mode,
                 int rc, const double *xx) {
    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, \"body_name\": \"%s\", "
           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
           "\"eph_name\": \"%s\", \"sid_mode\": %d, \"retflag\": %d, "
           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
           body, body_name, jd, flags, flag_name, eph_name, sid_mode, rc,
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
}

int main(void) {
    char serr[256];
    swe_set_ephe_path("../../ephe");
    printf("[\n");
    /* Tropical cases (sid_mode = 0). */
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int iep = 0; iep < NEPH; iep++) {
            for (int ie = 0; ie < NEPOCHS; ie++) {
                for (int ifl = 0; ifl < NFLAGS; ifl++) {
                    int flags = eph_combos[iep].flag | flag_combos[ifl].flag;
                    double xx[6];
                    memset(xx, 0, sizeof(xx));
                    int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                    if (rc < 0) {
                        fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x\n",
                                serr, bodies[ib], epochs[ie], flags);
                        return 1;
                    }
                    emit(bodies[ib], body_names[ib], epochs[ie], flags,
                         flag_combos[ifl].name, eph_combos[iep].name, 0, rc, xx);
                }
            }
        }
    }
    /* Sidereal cases: SEFLG_SIDEREAL | SEFLG_SPEED, one per sid_mode. */
    for (int is = 0; is < NSID; is++) {
        swe_set_sid_mode(sid_combos[is].sid_mode, 0, 0);
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int iep = 0; iep < NEPH; iep++) {
                for (int ie = 0; ie < NEPOCHS; ie++) {
                    int flags = eph_combos[iep].flag | SEFLG_SIDEREAL | SEFLG_SPEED;
                    double xx[6];
                    memset(xx, 0, sizeof(xx));
                    int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                    if (rc < 0) {
                        fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x sid=%s\n",
                                serr, bodies[ib], epochs[ie], flags, sid_combos[is].name);
                        return 1;
                    }
                    emit(bodies[ib], body_names[ib], epochs[ie], flags,
                         sid_combos[is].name, eph_combos[iep].name,
                         sid_combos[is].sid_mode, rc, xx);
                }
            }
        }
    }
    printf("\n]\n");
    return 0;
}
