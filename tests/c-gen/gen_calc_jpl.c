/*
 * Generates golden reference data for the JPL DE ephemeris calc pipeline.
 * Uses swe_calc with SEFLG_JPLEPH against de441.eph stored in swisseph-rs/ephe/.
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../../swisseph -o gen_calc_jpl gen_calc_jpl.c \
 *      ../../../swisseph/libswe.a -lm
 * Run (from tests/c-gen/):
 *   ./gen_calc_jpl > ../golden-data/calc_jpl.json
 */

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "swephexp.h"

/* Redirect stdout to /dev/null around swe_calc to suppress any debug
 * printf()s compiled into libswe.a, then restore it. */
static int suppress_stdout(void) {
    fflush(stdout);
    int saved = dup(STDOUT_FILENO);
    int null_fd = open("/dev/null", O_WRONLY);
    dup2(null_fd, STDOUT_FILENO);
    close(null_fd);
    return saved;
}

static void restore_stdout(int saved) {
    fflush(stdout);
    dup2(saved, STDOUT_FILENO);
    close(saved);
}

static int bodies[] = {
    SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER,
    SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO
};
static const char *body_names[] = {
    "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter",
    "Saturn", "Uranus", "Neptune", "Pluto"
};
#define NBODIES 10

static double epochs[] = {
    2159345.5,   /* ~1200 CE */
    2305447.5,   /* ~1600 CE */
    2378496.5,   /* 1800-Jan-1 */
    2433282.5,   /* ~B1950 */
    2451545.0,   /* J2000 */
    2460310.5,   /* 2024-Jan-1 */
    2488069.5,   /* 2100-Jan-1 */
};
#define NEPOCHS 7

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_SPEED,                    "default" },
    { SEFLG_SPEED | SEFLG_J2000,      "J2000" },
    { SEFLG_SPEED | SEFLG_NONUT,      "NONUT" },
    { SEFLG_SPEED | SEFLG_EQUATORIAL, "EQUATORIAL" },
    { SEFLG_SPEED | SEFLG_XYZ,        "XYZ" },
    { SEFLG_SPEED | SEFLG_RADIANS,    "RADIANS" },
    { SEFLG_SPEED | SEFLG_TRUEPOS,    "TRUEPOS" },
    { SEFLG_SPEED | SEFLG_NOABERR,    "NOABERR" },
    { SEFLG_SPEED | SEFLG_NOGDEFL,    "NOGDEFL" },
    { SEFLG_SPEED | SEFLG_ICRS,       "ICRS" },
    { SEFLG_SPEED3,                   "SPEED3" },
    { 0,                              "no_speed" },
};
#define NFLAGS 12

int main(void) {
    char serr[256];
    int first = 1;
    printf("[\n");
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = SEFLG_JPLEPH | flag_combos[ifl].flag;
                double xx[6];
                memset(xx, 0, sizeof(xx));
                /* Reset C library state before each call so file caching does
                 * not carry over between test cases. */
                swe_close();
                swe_set_ephe_path("../../ephe");
                swe_set_jpl_file("de441.eph");
                int saved = suppress_stdout();
                int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                restore_stdout(saved);
                if (rc < 0) {
                    fprintf(stderr, "skipping: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flags);
                    continue;
                }
                if (!first) printf(",\n");
                first = 0;
                printf("  {\"body\": %d, \"body_name\": \"%s\", "
                       "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", "
                       "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                       bodies[ib], body_names[ib],
                       epochs[ie], flags, flag_combos[ifl].name,
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
            }
        }
    }
    printf("\n]\n");
    swe_close();
    return 0;
}
