/*
 * Generates golden reference data for planetary moon calc pipeline.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_plmoon gen_plmoon.c \
 *      ../../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_plmoon > tests/golden-data/plmoon.json
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "swephexp.h"

#ifndef SEFLG_CENTER_BODY
#define SEFLG_CENTER_BODY (1024*1024)
#endif

/* Full 10-flag matrix bodies: one representative moon per planet family,
 * plus all five center-of-body (9n99) pseudo-bodies. */
static int bodies_full[] = {
    9401, /* Phobos */
    9501, /* Io */
    9606, /* Titan */
    9705, /* Miranda */
    9801, /* Triton */
    9901, /* Charon */
    9599, /* Jupiter COB */
    9699, /* Saturn COB */
    9799, /* Uranus COB */
    9899, /* Neptune COB */
    9999, /* Pluto COB */
};
#define NBODIES_FULL 11

/* Reduced matrix bodies: remaining 21 moons, SWIEPH|SPEED only. */
static int bodies_reduced[] = {
    9402, /* Deimos */
    9502, 9503, 9504,                   /* Europa, Ganymede, Callisto */
    9601, 9602, 9603, 9604, 9605,       /* Mimas, Enceladus, Tethys, Dione, Rhea */
    9607, 9608,                          /* Hyperion, Iapetus */
    9701, 9702, 9703, 9704,             /* Ariel, Umbriel, Titania, Oberon */
    9802, 9808,                          /* Nereid, Proteus */
    9902, 9903, 9904, 9905,             /* Nix, Hydra, Kerberos, Styx */
};
#define NBODIES_REDUCED 21

static double epochs_nonmars[] = {
    2451545.0,   /* J2000.0 */
    2460310.5,   /* 2024-Jan-1 */
    2433282.5,   /* 1950-Jan-1 */
    2416846.5,   /* 1905-Jan-1 -- within all sepm9*.se1 coverage (earliest lower bound 2415029.5) */
    2465424.5,   /* 2038-Jan-1 -- within all sepm9*.se1 coverage (earliest upper bound 2469081.5) */
};

static double epochs_mars[] = {
    2451545.0,   /* J2000.0 */
    2460310.5,   /* 2024-Jan-1 */
    2433282.5,   /* 1950-Jan-1 */
    2422324.5,   /* 1920-Jan-1 */
    2466162.5,   /* 2040-Jan-1 */
};
#define NEPOCHS 5

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_SWIEPH,                                           "SWIEPH" },
    { SEFLG_SWIEPH | SEFLG_SPEED,                              "SWIEPH_SPEED" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL,           "SWIEPH_SPEED_EQUATORIAL" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_XYZ,                  "SWIEPH_SPEED_XYZ" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000 | SEFLG_NONUT,  "SWIEPH_SPEED_J2000_NONUT" },
    { SEFLG_SWIEPH | SEFLG_TRUEPOS,                            "SWIEPH_TRUEPOS" },
    { SEFLG_SWIEPH | SEFLG_NONUT,                              "SWIEPH_NONUT" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_TOPOCTR,              "SWIEPH_SPEED_TOPOCTR" },
    { SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR,               "SWIEPH_SPEED_HELCTR" },
    { SEFLG_JPLEPH | SEFLG_SPEED,                              "JPLEPH_SPEED" },
};
#define NFLAGS 10

/* Center-of-body equivalence: planet ipl + CENTER_BODY should equal the
 * corresponding 9n99 pseudo-body's plain SWIEPH|SPEED output. */
static int cob_ipl[] = { SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO };
static int cob_body[] = { 9599, 9699, 9799, 9899, 9999 };
#define NCOB 5

static double cob_epochs[] = { 2451545.0, 2460310.5, 2433282.5 };
#define NCOB_EPOCHS 3

/* Cancellation rows: main-planet CENTER_BODY should reduce to a plain-planet
 * call (no observable center-of-body offset for the planet itself). */
static int cancel_ipl[] = { SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS };
#define NCANCEL 5

static int first = 1;

static void reset_state(void) {
    /* Reset C library state before each call so file caching does not carry
     * over between test cases. This gives deterministic, stateless golden
     * data that matches our stateless Rust implementation. */
    swe_close();
    swe_set_ephe_path("ephe");
    swe_set_topo(8.55, 47.37, 500.0);
}

/* Returns 1 if row was emitted, 0 if skipped (out-of-range). */
static int emit_row(int body, double jd, int flags, const char *flag_name) {
    char serr[256];
    char pname[AS_MAXCH];
    double xx[6];
    memset(xx, 0, sizeof(xx));

    reset_state();

    int rc = swe_calc(jd, body, flags, xx, serr);
    if (rc < 0) {
        /* Out-of-range for this body/epoch — skip gracefully. */
        fprintf(stderr,
                "skipping: body=%d jd=%.1f flags=0x%x: %s\n",
                body, jd, flags, serr);
        return 0;
    }

    swe_get_planet_name(body, pname);

    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, \"body_name\": \"%s\", "
           "\"jd\": %.17g, \"flags\": %d, \"flag_name\": \"%s\", "
           "\"retflag\": %d, "
           "\"output\": [%.17g, %.17g, %.17g, %.17g, %.17g, %.17g]}",
           body, pname,
           jd, flags, flag_name,
           rc,
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
    return 1;
}

int main(void) {
    swe_set_ephe_path("ephe");
    swe_set_topo(8.55, 47.37, 500.0);

    printf("[\n");

    /* Full 10-flag matrix x 5 epochs x 11 bodies. */
    for (int ib = 0; ib < NBODIES_FULL; ib++) {
        int body = bodies_full[ib];
        double *epochs = (body == 9401) ? epochs_mars : epochs_nonmars;
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                emit_row(body, epochs[ie], flag_combos[ifl].flag,
                         flag_combos[ifl].name);
            }
        }
    }

    /* Reduced matrix (SWIEPH|SPEED only) x 5 epochs x 21 remaining moons. */
    for (int ib = 0; ib < NBODIES_REDUCED; ib++) {
        int body = bodies_reduced[ib];
        double *epochs = (body == 9402) ? epochs_mars : epochs_nonmars;
        for (int ie = 0; ie < NEPOCHS; ie++) {
            emit_row(body, epochs[ie], SEFLG_SWIEPH | SEFLG_SPEED,
                     "SWIEPH_SPEED");
        }
    }

    /* Center-of-body equivalence rows. */
    for (int ic = 0; ic < NCOB; ic++) {
        for (int ie = 0; ie < NCOB_EPOCHS; ie++) {
            emit_row(cob_ipl[ic], cob_epochs[ie],
                     SEFLG_CENTER_BODY | SEFLG_SWIEPH | SEFLG_SPEED,
                     "CENTER_BODY_planet");
            emit_row(cob_body[ic], cob_epochs[ie],
                     SEFLG_SWIEPH | SEFLG_SPEED,
                     "COB_equiv");
        }
    }

    /* Cancellation rows: main planets with CENTER_BODY at J2000. */
    for (int ic = 0; ic < NCANCEL; ic++) {
        emit_row(cancel_ipl[ic], 2451545.0,
                 SEFLG_CENTER_BODY | SEFLG_SWIEPH | SEFLG_SPEED,
                 "CENTER_BODY_cancel");
    }

    printf("\n]\n");
    swe_close();
    return 0;
}
