/*
 * Generates golden reference data for fictitious planets (swe_calc, ipl 40–58).
 * Battery: 12 bodies × 4 epochs × 6 flag combos ≈ 288 cases (minus C errors).
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_fictitious gen_fictitious.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_fictitious > ../golden-data/fictitious.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static int bodies[] = {
    SE_CUPIDO, SE_HADES, SE_ZEUS, SE_KRONOS,
    SE_APOLLON, SE_ADMETOS, SE_VULKANUS, SE_POSEIDON,
    SE_ISIS, SE_NIBIRU,
    SE_VULCAN, SE_WHITE_MOON,
};
static const char *body_names[] = {
    "Cupido", "Hades", "Zeus", "Kronos",
    "Apollon", "Admetos", "Vulkanus", "Poseidon",
    "Isis", "Nibiru",
    "Vulcan", "WhiteMoon",
};
#define NBODIES 12

static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2415020.0,    /* J1900 */
    2460600.5,    /* 2024-Oct */
    2305447.5,    /* 1600-Jan-1 */
};
#define NEPOCHS 4

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_MOSEPH | SEFLG_SPEED,                          "moseph_speed" },
    { SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL,      "moseph_equatorial" },
    { SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_XYZ,             "moseph_xyz" },
    { SEFLG_MOSEPH,                                        "moseph_nospeed" },
    { SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_HELCTR,          "moseph_helctr" },
    { SEFLG_SWIEPH | SEFLG_SPEED,                         "swieph_speed" },
};
#define NFLAGS 6

int main(void) {
    char serr[256];
    swe_set_ephe_path("../../ephe");
    int first = 1;
    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                double xx[6];
                memset(xx, 0, sizeof(xx));
                memset(serr, 0, sizeof(serr));
                int rc = swe_calc(epochs[ie], bodies[ib], flag_combos[ifl].flag, xx, serr);
                if (!first) printf(",\n");
                first = 0;
                if (rc < 0) {
                    printf("  {\"body\": %d, \"body_name\": \"%s\", "
                           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
                           "\"error\": \"%s\"}",
                           bodies[ib], body_names[ib],
                           epochs[ie], flag_combos[ifl].flag, flag_combos[ifl].name,
                           serr);
                    fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flag_combos[ifl].flag);
                } else {
                    printf("  {\"body\": %d, \"body_name\": \"%s\", "
                           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
                           "\"retflag\": %d, "
                           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                           bodies[ib], body_names[ib],
                           epochs[ie], flag_combos[ifl].flag, flag_combos[ifl].name,
                           rc,
                           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
                }
            }
        }
    }
    printf("\n]\n");
    swe_close();
    return 0;
}
