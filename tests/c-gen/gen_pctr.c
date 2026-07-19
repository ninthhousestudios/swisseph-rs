/*
 * Generates golden reference data for swe_calc_pctr (planet-centric
 * calculation: position of ipl as seen from iplctr).
 *
 * Matrix (swisseph-rs/98):
 *   - 6 (ipl, iplctr) pairs x 3 epochs x 10 non-sidereal flag combos
 *     (1 MOSEPH reject-path + 9 SWIEPH success across every output flag:
 *      ecliptic/equatorial x polar/xyz, J2000 precession-skip, TRUEPOS,
 *      NOABERR, NOGDEFL). = 180 cases.
 *   - A sidereal block (SWIEPH|SIDEREAL): traditional Lahiri, ECL_T0 (J2000),
 *     SSY_PLANE, and USER|ECL_T0, over 2 pairs x 2 epochs. = 16 cases.
 *
 * Each case records sid_mode (-1 = tropical) plus sid_t0/sid_ayan so the Rust
 * side can reproduce swe_set_sid_mode exactly.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_pctr tests/c-gen/gen_pctr.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_pctr > tests/golden-data/pctr.json
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swephexp.h"

static int first_item = 1;

static void comma(void) {
    if (!first_item) printf(",");
    first_item = 0;
    printf("\n");
}

struct body_pair {
    int ipl;
    int iplctr;
};

static struct body_pair pairs[] = {
    { SE_SUN,     SE_MARS },
    { SE_EARTH,   SE_MARS },
    { SE_JUPITER, SE_MARS },
    { SE_MOON,    SE_VENUS },
    { SE_SATURN,  SE_JUPITER },
    { SE_EARTH,   SE_MOON },
};
#define NPAIRS 6

static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2460600.0,   /* modern */
    2415020.5,   /* 1900-Jan-0.5 */
};
#define NEPOCHS 3

/* Non-sidereal flag combos. First is the MOSEPH reject-path regression; the
 * rest are SWIEPH successes covering every output-frame branch. */
static int flag_combos[] = {
    SEFLG_MOSEPH | SEFLG_SPEED,
    SEFLG_SWIEPH | SEFLG_SPEED,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_XYZ,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_XYZ,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000 | SEFLG_EQUATORIAL | SEFLG_XYZ,
    SEFLG_SWIEPH | SEFLG_TRUEPOS,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NOABERR,
    SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NOGDEFL,
};
#define NFLAGS 10

static int count = 0;

static void emit_case(int ipl, int iplctr, double tjd, int iflag,
                      int sid_mode, double sid_t0, double sid_ayan) {
    char serr[256];
    double xx[6];
    memset(xx, 0, sizeof(xx));
    memset(serr, 0, sizeof(serr));

    if (sid_mode >= 0)
        swe_set_sid_mode(sid_mode, sid_t0, sid_ayan);

    int rc = swe_calc_pctr(tjd, ipl, iplctr, iflag, xx, serr);
    int ok = (rc == ERR) ? 0 : 1;

    comma();
    printf("  {\"ipl\":%d,\"iplctr\":%d,\"tjd\":%.17g,"
           "\"iflag\":%d,\"retflag\":%d,"
           "\"sid_mode\":%d,\"sid_t0\":%.17g,\"sid_ayan\":%.17g,"
           "\"xx\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],"
           "\"ok\":%d}",
           ipl, iplctr, tjd, iflag, rc,
           sid_mode, sid_t0, sid_ayan,
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5],
           ok);
    count++;
}

/* Sidereal sub-cases: (sid_mode with SE_SIDBIT bits, t0, ayan). */
struct sid_case {
    int sid_mode;
    double t0;
    double ayan;
};

static struct sid_case sid_cases[] = {
    { SE_SIDM_LAHIRI, 0, 0 },                            /* traditional ayanamsha */
    { SE_SIDM_J2000, 0, 0 },                             /* ECL_T0 (auto-set for 18) */
    { SE_SIDM_LAHIRI | SE_SIDBIT_SSY_PLANE, 0, 0 },      /* solar-system-plane */
    { SE_SIDM_USER | SE_SIDBIT_ECL_T0, 2451545.0, 25.0 },/* user ECL_T0 */
};
#define NSID 4

static struct body_pair sid_pairs[] = {
    { SE_SUN,    SE_MARS },
    { SE_SATURN, SE_JUPITER },
};
#define NSIDPAIRS 2

static double sid_epochs[] = { 2451545.0, 2460600.0 };
#define NSIDEPOCHS 2

int main(void) {
    swe_set_ephe_path("../../ephe");

    printf("{\n");
    printf("\"pctr\": [");
    first_item = 1;

    for (int ip = 0; ip < NPAIRS; ip++)
        for (int ie = 0; ie < NEPOCHS; ie++)
            for (int ifl = 0; ifl < NFLAGS; ifl++)
                emit_case(pairs[ip].ipl, pairs[ip].iplctr, epochs[ie],
                          flag_combos[ifl], -1, 0, 0);

    for (int ip = 0; ip < NSIDPAIRS; ip++)
        for (int ie = 0; ie < NSIDEPOCHS; ie++)
            for (int is = 0; is < NSID; is++)
                emit_case(sid_pairs[ip].ipl, sid_pairs[ip].iplctr, sid_epochs[ie],
                          SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED,
                          sid_cases[is].sid_mode, sid_cases[is].t0, sid_cases[is].ayan);

    printf("\n]\n");
    printf("}\n");

    fprintf(stderr, "gen_pctr: emitted %d cases\n", count);

    swe_close();
    return 0;
}
