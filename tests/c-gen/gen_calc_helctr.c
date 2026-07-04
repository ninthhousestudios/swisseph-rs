/*
 * Generates golden reference data for the heliocentric calc pipeline
 * (swe_calc with SEFLG_HELCTR) across all three backends in one file.
 *
 * SEFLG_HELCTR was ported for phenomena (swisseph-rs/83) but only ever
 * exercised transitively via the pheno battery, which never sets SEFLG_SPEED
 * on its heliocentric calls and never touches the JPL Moon. This generator
 * gives HELCTR first-class coverage: Sun..Pluto + Moon, the polar and XYZ
 * frames, with/without J2000/EQUATORIAL, and with/without SPEED, over the
 * Moshier / Swiss / JPL backends.
 *
 * Each record carries a "backend" string ("moshier"/"sweph"/"jpl") so the
 * Rust side can route to the matching Ephemeris. JPL rows are still emitted
 * even if de441.eph is absent here (they simply won't be — swe_calc returns
 * an error and the row is skipped); the Rust test independently skips the
 * jpl rows when the file is missing.
 *
 * Epochs avoid the exact sepl_18 tfstart (2378496.5): at that .se1 file
 * boundary C's stateful file caching diverges from a stateless port (see the
 * <stateless_tolerance> note in CLAUDE.md and sutra lesson 019f2337 §4).
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../../swisseph -o gen_calc_helctr gen_calc_helctr.c \
 *      ../../../swisseph/libswe.a -lm
 * Run (from tests/c-gen/):
 *   ./gen_calc_helctr > ../golden-data/calc_helctr.json
 */

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "swephexp.h"

/* Redirect stdout to /dev/null around swe_calc to suppress any debug
 * printf()s compiled into libswe.a (the JPL reader emits some), then
 * restore it. */
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

/* Sun..Pluto + Moon + Earth (SE_SUN=0 .. SE_PLUTO=9; Moon=1; Earth=14). */
static int bodies[] = {
    SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER,
    SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO, SE_EARTH
};
static const char *body_names[] = {
    "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter",
    "Saturn", "Uranus", "Neptune", "Pluto", "Earth"
};
#define NBODIES 11

static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2460310.5,    /* 2024-Jan-1 */
    2433282.5,    /* 1950-Jan-1 */
    2378500.5,    /* 1800-Jan-5 (off the sepl_18 file boundary) */
    2488069.5,    /* 2100-Jan-1 */
};
#define NEPOCHS 5

struct flag_combo {
    int flag;
    const char *name;
};

/* Every combo carries SEFLG_HELCTR or SEFLG_BARYCTR. plaus_iflag forces
 * NOABERR|NOGDEFL for both, so these differ only in output frame, center
 * flag, and whether speed is requested. BARYCTR works for planets + Moon +
 * Earth on Swiss/JPL (Moshier rejects BARYCTR and is skipped via rc<0).
 * Sun BARYCTR skipped — C's app_pos_etc_sun doesn't handle it cleanly. */
static struct flag_combo flag_combos[] = {
    { SEFLG_HELCTR | SEFLG_SPEED,                                   "polar" },
    { SEFLG_HELCTR,                                                 "polar_nospeed" },
    { SEFLG_HELCTR | SEFLG_SPEED | SEFLG_XYZ,                       "xyz" },
    { SEFLG_HELCTR | SEFLG_XYZ,                                     "xyz_nospeed" },
    { SEFLG_HELCTR | SEFLG_SPEED | SEFLG_J2000,                     "j2000" },
    { SEFLG_HELCTR | SEFLG_SPEED | SEFLG_EQUATORIAL,               "equatorial" },
    { SEFLG_HELCTR | SEFLG_SPEED | SEFLG_XYZ | SEFLG_J2000,        "xyz_j2000" },
    { SEFLG_HELCTR | SEFLG_SPEED | SEFLG_XYZ | SEFLG_EQUATORIAL,   "xyz_equatorial" },
    { SEFLG_BARYCTR | SEFLG_SPEED,                                  "bary_polar" },
    { SEFLG_BARYCTR | SEFLG_SPEED | SEFLG_XYZ,                      "bary_xyz" },
    { SEFLG_BARYCTR | SEFLG_SPEED | SEFLG_J2000,                    "bary_j2000" },
    { SEFLG_BARYCTR | SEFLG_SPEED | SEFLG_XYZ | SEFLG_J2000,       "bary_xyz_j2000" },
};
#define NFLAGS 12

struct backend {
    int flag;
    const char *name;
    int is_jpl;
};

static struct backend backends[] = {
    { SEFLG_MOSEPH, "moshier", 0 },
    { SEFLG_SWIEPH, "sweph",   0 },
    { SEFLG_JPLEPH, "jpl",     1 },
};
#define NBACKENDS 3

int main(void) {
    char serr[256];
    int first = 1;
    printf("[\n");
    for (int be = 0; be < NBACKENDS; be++) {
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS; ie++) {
                for (int ifl = 0; ifl < NFLAGS; ifl++) {
                    int flags = backends[be].flag | flag_combos[ifl].flag;
                    /* Skip BARYCTR for the Sun — apparent_sun's frame
                     * construction doesn't handle BARYCTR Sun. */
                    if ((flag_combos[ifl].flag & SEFLG_BARYCTR) && bodies[ib] == SE_SUN)
                        continue;
                    double xx[6];
                    memset(xx, 0, sizeof(xx));
                    /* Reset C library state before each call so file caching
                     * does not carry over between cases — deterministic,
                     * stateless golden data matching the Rust port. */
                    swe_close();
                    if (backends[be].is_jpl) {
                        swe_set_ephe_path("../../ephe");
                        swe_set_jpl_file("de441.eph");
                    } else {
                        swe_set_ephe_path("../../../swisseph/ephe");
                    }
                    int saved = suppress_stdout();
                    int rc = swe_calc(epochs[ie], bodies[ib], flags, xx, serr);
                    restore_stdout(saved);
                    if (rc < 0) {
                        fprintf(stderr,
                                "skipping: %s backend=%s body=%d jd=%.1f flags=0x%x\n",
                                serr, backends[be].name, bodies[ib],
                                epochs[ie], flags);
                        continue;
                    }
                    if (!first) printf(",\n");
                    first = 0;
                    printf("  {\"backend\": \"%s\", \"body\": %d, "
                           "\"body_name\": \"%s\", \"jd\": %.20e, "
                           "\"flags\": %d, \"flag_name\": \"%s\", "
                           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                           backends[be].name, bodies[ib], body_names[ib],
                           epochs[ie], flags, flag_combos[ifl].name,
                           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
                }
            }
        }
    }
    printf("\n]\n");
    swe_close();
    return 0;
}
