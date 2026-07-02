/*
 * Generates golden reference data for the orbital-elements module (PNOC 6):
 *   swe_get_orbital_elements  ("elements" key)
 *   swe_orbit_max_min_true_distance  ("maxmin" key)
 *
 * "elements": records dret[0..17) + retflag.
 *
 *   Battery A (default / HELCTR / ORBEL_AA on Moshier):
 *     Bodies: SE_MERCURY..SE_PLUTO (8) + SE_EARTH (9 bodies).
 *     Flags:  MOSEPH|SPEED, MOSEPH|SPEED|HELCTR,
 *             MOSEPH|SPEED|TOPOCTR  (TOPOCTR is bit-aliased to SEFLG_ORBEL_AA —
 *             "sum masses inside the orbit"; exercises get_gmsm's AA branch).
 *     Epochs: 4 (incl. pre-1900).  => 9 * 4 * 3 = 108 cases.
 *
 *   Battery B (BARYCTR gate, inside-6-AU fallback to HELCTR):
 *     Bodies: SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_EARTH — every body
 *             whose heliocentric distance is < 6 AU, so the `if (BARYCTR && r>6)`
 *             gate falls through to HELCTR and produces a valid result. (Beyond
 *             6 AU, MOSEPH|BARYCTR errors: "barycentric Moshier positions are
 *             not supported" — that path is unsupported by both C-Moshier and
 *             this codebase's calc pipeline, so it is not exercised.)
 *     Flags:  MOSEPH|SPEED|BARYCTR.
 *     Epochs: epochs[0..1].  => 5 * 2 = 10 cases.
 *
 *   Battery C (Swiss backend):
 *     Bodies: SE_MERCURY, SE_JUPITER, SE_PLUTO.
 *     Flags:  SWIEPH|SPEED, SWIEPH|SPEED|HELCTR.
 *     Epochs: epochs[0..1] (J2000, 2024 — off any .se1 file boundary).
 *             => 3 * 2 * 2 = 12 cases.
 *
 *   Total elements = 130 cases.
 *
 * "maxmin": records dmax/dmin/dtrue + retflag.
 *   Bodies: SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_PLUTO.
 *   Flags:  MOSEPH (geocentric two-ellipse search),
 *           MOSEPH|HELCTR (single-ellipse helio branch).
 *   Epochs: epochs[0..2] (J2000, 2024, 1950).  => 5 * 3 * 2 = 30 cases.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_orbit tests/c-gen/gen_orbit.c \
 *      ../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_orbit > tests/golden-data/orbit.json
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swephexp.h"

static double epochs[] = {
    2451545.0, /* J2000.0 */
    2460310.5, /* 2024-Jan-1 */
    2433282.5, /* 1950-Jan-1 */
    2378496.5, /* 1800-Jan-1 (pre-1900) */
};

struct flag_combo {
    int flag;
    const char *name;
};

static int first = 1;

static void emit_elements(int body, const char *body_name, double jd, int flags,
                          const char *flag_name, int rc, const double *dret) {
    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, \"body_name\": \"%s\", "
           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", \"retflag\": %d, "
           "\"dret\": [",
           body, body_name, jd, flags, flag_name, rc);
    for (int i = 0; i < 17; i++) {
        printf("%.20e%s", dret[i], i < 16 ? ", " : "");
    }
    printf("]}");
}

static void run_elements(double jd, int body, const char *body_name, int flags,
                         const char *flag_name) {
    char serr[256];
    double dret[50];
    memset(dret, 0, sizeof(dret));
    int rc = swe_get_orbital_elements(jd, body, flags, dret, serr);
    if (rc < 0) {
        fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x\n", serr, body, jd, flags);
        exit(1);
    }
    emit_elements(body, body_name, jd, flags, flag_name, rc, dret);
}

static void emit_maxmin(int body, const char *body_name, double jd, int flags,
                        const char *flag_name, int rc, double dmax, double dmin,
                        double dtrue) {
    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, \"body_name\": \"%s\", "
           "\"jd\": %.20e, \"flags\": %d, \"flag_name\": \"%s\", \"retflag\": %d, "
           "\"dmax\": %.20e, \"dmin\": %.20e, \"dtrue\": %.20e}",
           body, body_name, jd, flags, flag_name, rc, dmax, dmin, dtrue);
}

static void run_maxmin(double jd, int body, const char *body_name, int flags,
                       const char *flag_name) {
    char serr[256];
    double dmax = 0, dmin = 0, dtrue = 0;
    int rc = swe_orbit_max_min_true_distance(jd, body, flags, &dmax, &dmin, &dtrue, serr);
    if (rc < 0) {
        fprintf(stderr, "error: %s body=%d jd=%.1f flags=0x%x\n", serr, body, jd, flags);
        exit(1);
    }
    emit_maxmin(body, body_name, jd, flags, flag_name, rc, dmax, dmin, dtrue);
}

int main(void) {
    swe_set_ephe_path("../swisseph/ephe");

    /* --- elements --- */
    /* Battery A: the full Mercury..Pluto range (SE_MERCURY..SE_PLUTO) + Earth. */
    static int bodyA[] = {SE_MERCURY, SE_VENUS,   SE_MARS,   SE_JUPITER, SE_SATURN,
                          SE_URANUS,  SE_NEPTUNE, SE_PLUTO,  SE_EARTH};
    static const char *bodyA_names[] = {"Mercury", "Venus",   "Mars",    "Jupiter", "Saturn",
                                        "Uranus",  "Neptune", "Pluto",   "Earth"};
    static struct flag_combo flagsA[] = {
        {SEFLG_MOSEPH | SEFLG_SPEED, "MOSEPH_SPEED"},
        {SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_HELCTR, "MOSEPH_HELCTR"},
        {SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_TOPOCTR, "MOSEPH_ORBEL_AA"},
    };
    int nbodyA = sizeof(bodyA) / sizeof(bodyA[0]);
    int nflagsA = sizeof(flagsA) / sizeof(flagsA[0]);
    int nepochs = sizeof(epochs) / sizeof(epochs[0]);

    printf("{\n\"elements\": [\n");
    for (int ib = 0; ib < nbodyA; ib++)
        for (int ie = 0; ie < nepochs; ie++)
            for (int ifl = 0; ifl < nflagsA; ifl++)
                run_elements(epochs[ie], bodyA[ib], bodyA_names[ib],
                             flagsA[ifl].flag, flagsA[ifl].name);

    /* Battery B: BARYCTR gate, inside-6-AU bodies only. */
    static int bodyB[] = {SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_EARTH};
    static const char *bodyB_names[] = {"Mercury", "Venus", "Mars", "Jupiter", "Earth"};
    int nbodyB = sizeof(bodyB) / sizeof(bodyB[0]);
    for (int ib = 0; ib < nbodyB; ib++)
        for (int ie = 0; ie < 2; ie++)
            run_elements(epochs[ie], bodyB[ib], bodyB_names[ib],
                         SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_BARYCTR, "MOSEPH_BARYCTR");

    /* Battery C: Swiss backend. */
    static int bodyC[] = {SE_MERCURY, SE_JUPITER, SE_PLUTO};
    static const char *bodyC_names[] = {"Mercury", "Jupiter", "Pluto"};
    static struct flag_combo flagsC[] = {
        {SEFLG_SWIEPH | SEFLG_SPEED, "SWIEPH_SPEED"},
        {SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR, "SWIEPH_HELCTR"},
    };
    int nbodyC = sizeof(bodyC) / sizeof(bodyC[0]);
    int nflagsC = sizeof(flagsC) / sizeof(flagsC[0]);
    for (int ib = 0; ib < nbodyC; ib++)
        for (int ie = 0; ie < 2; ie++)
            for (int ifl = 0; ifl < nflagsC; ifl++)
                run_elements(epochs[ie], bodyC[ib], bodyC_names[ib],
                             flagsC[ifl].flag, flagsC[ifl].name);

    /* --- maxmin --- */
    printf("\n],\n\"maxmin\": [\n");
    first = 1;
    static int bodyM[] = {SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_PLUTO};
    static const char *bodyM_names[] = {"Mercury", "Venus", "Mars", "Jupiter", "Pluto"};
    static struct flag_combo flagsM[] = {
        {SEFLG_MOSEPH, "MOSEPH"},
        {SEFLG_MOSEPH | SEFLG_HELCTR, "MOSEPH_HELCTR"},
    };
    int nbodyM = sizeof(bodyM) / sizeof(bodyM[0]);
    int nflagsM = sizeof(flagsM) / sizeof(flagsM[0]);
    for (int ib = 0; ib < nbodyM; ib++)
        for (int ie = 0; ie < 3; ie++)
            for (int ifl = 0; ifl < nflagsM; ifl++)
                run_maxmin(epochs[ie], bodyM[ib], bodyM_names[ib],
                           flagsM[ifl].flag, flagsM[ifl].name);

    printf("\n]\n}\n");
    return 0;
}
