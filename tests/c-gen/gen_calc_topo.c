/*
 * Generates golden reference data for topocentric swe_calc (SEFLG_TOPOCTR).
 *
 * Three sub-matrices:
 *   - moshier: full observer x body x epoch x flag-shape matrix (no file I/O,
 *     exercises the Moshier-specific topo_offset wiring in calc_planet/
 *     calc_sun/calc_moon).
 *   - sweph / jpl: reduced observer set (still all bodies, both flag shapes,
 *     including the SPEED3 file-boundary epoch) exercising the generic
 *     apparent_planet/apparent_sun/apparent_moon xobs wiring shared by the
 *     SWIEPH and JPL backends — previously untested topocentrically
 *     (see swisseph-rs/68 Codex review, swisseph-rs/80).
 *
 * flag shapes: "speed" (TOPOCTR|EQUATORIAL|SPEED, forces SPEED3 since
 * !NOABERR) and "speed_noaberr" (adds NOABERR, which skips the SPEED3-forcing
 * rule and exercises the analytic-speed non-SPEED3 path).
 *
 * Compile:
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_calc_topo tests/c-gen/gen_calc_topo.c \
 *      -L../swisseph -lswe -lm
 * Run (from tests/c-gen/, so relative ephe path resolves):
 *   ./gen_calc_topo > ../golden-data/calc_topo.json
 */

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "swephexp.h"

struct observer {
    double lon, lat, alt;
};

static struct observer observers_full[] = {
    { 8.55, 47.37, 500.0 },     /* Zurich */
    { -74.0, 40.7, 10.0 },      /* New York */
    { 139.7, 35.7, 40.0 },      /* Tokyo */
};
#define NOBS_FULL 3

static struct observer observers_reduced[] = {
    { 8.55, 47.37, 500.0 },     /* Zurich */
    { 139.7, 35.7, 40.0 },      /* Tokyo */
};
#define NOBS_REDUCED 2

static int bodies[] = { SE_SUN, SE_MOON, SE_MERCURY, SE_MARS, SE_VENUS };
static const char *body_names[] = { "Sun", "Moon", "Mercury", "Mars", "Venus" };
#define NBODIES 5

static double epochs_full[] = {
    2451545.0,    /* J2000.0 */
    2378496.5,    /* 1800-Jan-1 — sepl_18 tfstart, SPEED3 file-boundary case */
    2469807.5,    /* 2050-Jan-1 */
};
#define NEPOCHS_FULL 3

static double epochs_reduced[] = {
    2451545.0,    /* J2000.0 */
    2378496.5,    /* 1800-Jan-1 — sepl_18 tfstart, SPEED3 file-boundary case */
};
#define NEPOCHS_REDUCED 2

struct flag_shape {
    int flag;
    const char *name;
};

static struct flag_shape flag_shapes[] = {
    { SEFLG_TOPOCTR | SEFLG_EQUATORIAL | SEFLG_SPEED,             "speed" },
    { SEFLG_TOPOCTR | SEFLG_EQUATORIAL | SEFLG_SPEED | SEFLG_NOABERR, "speed_noaberr" },
};
#define NSHAPES 2

/* Redirect stdout to /dev/null around swe_calc to suppress any debug
 * printf()s compiled into libswe.a, then restore it. Used for the JPL
 * sub-matrix, matching gen_calc_jpl.c. */
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

static void emit(int *first, struct observer *o, int body, const char *body_name,
                  double jd, int flags, const char *flag_name, const char *ephem_name,
                  double *xx) {
    if (!*first) printf(",\n");
    *first = 0;
    printf("  {\"lon\": %.20e, \"lat\": %.20e, \"alt\": %.20e, "
           "\"body\": %d, \"body_name\": \"%s\", \"jd\": %.20e, "
           "\"flags\": %d, \"flag_name\": \"%s\", \"ephemeris\": \"%s\", "
           "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
           o->lon, o->lat, o->alt, body, body_name, jd, flags, flag_name, ephem_name,
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
}

int main(void) {
    char serr[256];
    int first = 1;
    printf("[\n");

    /* --- Moshier: full matrix, no file I/O. --- */
    swe_set_ephe_path(NULL);
    for (int io = 0; io < NOBS_FULL; io++) {
        swe_set_topo(observers_full[io].lon, observers_full[io].lat, observers_full[io].alt);
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS_FULL; ie++) {
                for (int is = 0; is < NSHAPES; is++) {
                    int flags = SEFLG_MOSEPH | flag_shapes[is].flag;
                    double xx[6];
                    memset(xx, 0, sizeof(xx));
                    int rc = swe_calc(epochs_full[ie], bodies[ib], flags, xx, serr);
                    if (rc < 0) {
                        fprintf(stderr, "skipping moshier: %s body=%d jd=%.1f flags=0x%x\n",
                                serr, bodies[ib], epochs_full[ie], flags);
                        continue;
                    }
                    emit(&first, &observers_full[io], bodies[ib], body_names[ib],
                         epochs_full[ie], flags, flag_shapes[is].name, "moshier", xx);
                }
            }
        }
    }

    /* --- SWIEPH: reduced observers, all bodies/epochs/flag-shapes. --- */
    for (int io = 0; io < NOBS_REDUCED; io++) {
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS_REDUCED; ie++) {
                for (int is = 0; is < NSHAPES; is++) {
                    /* Reset C library state before each call so file caching does not
                     * carry over between test cases (matches gen_calc_sweph.c). */
                    swe_close();
                    swe_set_ephe_path("../../../swisseph/ephe");
                    swe_set_topo(observers_reduced[io].lon, observers_reduced[io].lat,
                                 observers_reduced[io].alt);
                    int flags = SEFLG_SWIEPH | flag_shapes[is].flag;
                    double xx[6];
                    memset(xx, 0, sizeof(xx));
                    int rc = swe_calc(epochs_reduced[ie], bodies[ib], flags, xx, serr);
                    if (rc < 0) {
                        fprintf(stderr, "skipping sweph: %s body=%d jd=%.1f flags=0x%x\n",
                                serr, bodies[ib], epochs_reduced[ie], flags);
                        continue;
                    }
                    emit(&first, &observers_reduced[io], bodies[ib], body_names[ib],
                         epochs_reduced[ie], flags, flag_shapes[is].name, "sweph", xx);
                }
            }
        }
    }

    /* --- JPL: reduced observers, all bodies/epochs/flag-shapes. --- */
    for (int io = 0; io < NOBS_REDUCED; io++) {
        for (int ib = 0; ib < NBODIES; ib++) {
            for (int ie = 0; ie < NEPOCHS_REDUCED; ie++) {
                for (int is = 0; is < NSHAPES; is++) {
                    swe_close();
                    swe_set_ephe_path("../../../swisseph/ephe");
                    swe_set_jpl_file("de441.eph");
                    swe_set_topo(observers_reduced[io].lon, observers_reduced[io].lat,
                                 observers_reduced[io].alt);
                    int flags = SEFLG_JPLEPH | flag_shapes[is].flag;
                    double xx[6];
                    memset(xx, 0, sizeof(xx));
                    int saved = suppress_stdout();
                    int rc = swe_calc(epochs_reduced[ie], bodies[ib], flags, xx, serr);
                    restore_stdout(saved);
                    if (rc < 0) {
                        fprintf(stderr, "skipping jpl: %s body=%d jd=%.1f flags=0x%x\n",
                                serr, bodies[ib], epochs_reduced[ie], flags);
                        continue;
                    }
                    emit(&first, &observers_reduced[io], bodies[ib], body_names[ib],
                         epochs_reduced[ie], flags, flag_shapes[is].name, "jpl", xx);
                }
            }
        }
    }

    printf("\n]\n");
    swe_close();
    return 0;
}
