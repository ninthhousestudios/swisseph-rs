/*
 * Generates golden reference data for the osculating lunar node (SE_TRUE_NODE)
 * and osculating apogee / "true Lilith" (SE_OSCU_APOG) through swe_calc
 * (lunar_osc_elem + swi_plan_for_osc_elem), swisseph-rs/84.
 *
 * Bodies: SE_TRUE_NODE (11), SE_OSCU_APOG (13).
 * Backends: SEFLG_MOSEPH, SEFLG_SWIEPH.
 * Flags:   SPEED, SPEED|EQUATORIAL, SPEED|XYZ, SPEED|NONUT, SPEED|J2000, 0.
 * Epochs:  the 7 gen_calc.c epochs.
 * => 2 bodies * 2 backends * 6 flags * 7 epochs = 168 cases.
 *
 * Records xx[0..6) + the swe_calc return flag.
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

int main(void) {
    char serr[256];
    swe_set_ephe_path("../swisseph/ephe");
    int first = 1;
    printf("[\n");
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
                    if (!first) printf(",\n");
                    first = 0;
                    printf("  {\"body\": %d, \"body_name\": \"%s\", "
                           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
                           "\"eph_name\": \"%s\", \"retflag\": %d, "
                           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                           bodies[ib], body_names[ib], epochs[ie], flags,
                           flag_combos[ifl].name, eph_combos[iep].name, rc,
                           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
                }
            }
        }
    }
    printf("\n]\n");
    return 0;
}
