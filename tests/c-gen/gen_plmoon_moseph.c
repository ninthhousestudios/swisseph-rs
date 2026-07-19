/*
 * Generates golden reference data for MOSEPH-mode planetary moon calc,
 * plus the §2 normalization quirk rows.
 *
 * CRITICAL: This must be a SEPARATE binary that issues ONLY SEFLG_MOSEPH calls.
 * C's MOSEPH plmoon path reads the parent planet via Moshier (compute_pipeline)
 * which depends on process-global swed.pldat[SEI_SUNBARY] staying zero — any
 * prior SWIEPH/JPLEPH call would contaminate it. Same hazard as asteroids
 * (gen_asteroid_moseph.c, swisseph-rs/101 decision 1).
 *
 * The sepm*.se1 file read (moon offset) is inherent to C's plmoon path even
 * under MOSEPH — purity means no SWIEPH *planet* calls, not no file I/O.
 *
 * Compile (from repo root):
 *   cc -Wall -I../../swisseph -o tests/c-gen/gen_plmoon_moseph \
 *      tests/c-gen/gen_plmoon_moseph.c ../../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_plmoon_moseph > tests/golden-data/plmoon_moseph.json
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "swephexp.h"

#ifndef SEFLG_CENTER_BODY
#define SEFLG_CENTER_BODY (1024*1024)
#endif

/* Bodies: one per planet family + two COBs. */
static int bodies[] = {
    9401,  /* Phobos (Mars) */
    9501,  /* Io (Jupiter) */
    9599,  /* Jupiter COB */
    9699,  /* Saturn COB */
    9901,  /* Charon (Pluto) */
};
#define NBODIES 5

/* Epochs inside ALL ranges: Phobos is the tightest (2415015.5–2469082.5). */
static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2460310.5,   /* 2024-Jan-1 */
    2433282.5,   /* 1950-Jan-1 */
};
#define NEPOCHS 3

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_MOSEPH,                                "MOSEPH" },
    { SEFLG_MOSEPH | SEFLG_SPEED,                  "MOSEPH_SPEED" },
    { SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL, "MOSEPH_SPEED_EQUATORIAL" },
};
#define NFLAGS 3

static int first = 1;

static void emit(int body, double jd, int flags, const char *flag_name,
                 int rc, const double *xx) {
    if (!first) printf(",\n");
    first = 0;
    printf("  {\"body\": %d, "
           "\"jd\": %.17g, \"flags\": %d, \"flag_name\": \"%s\", "
           "\"retflag\": %d, "
           "\"output\": [%.17g, %.17g, %.17g, %.17g, %.17g, %.17g]}",
           body, jd, flags, flag_name, rc,
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
}

int main(void) {
    char serr[256];
    double xx[6];

    swe_set_ephe_path("../../ephe");

    printf("{\"moseph\": [\n");

    /* MOSEPH calc battery: 5 bodies × 3 epochs × 3 flags = 45 cases. */
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                memset(xx, 0, sizeof(xx));
                int rc = swe_calc(epochs[ie], bodies[ib],
                                  flag_combos[ifl].flag, xx, serr);
                if (rc < 0) {
                    fprintf(stderr, "ERROR: body=%d jd=%.1f flags=0x%x: %s\n",
                            bodies[ib], epochs[ie], flag_combos[ifl].flag, serr);
                    return 1;
                }
                emit(bodies[ib], epochs[ie], flag_combos[ifl].flag,
                     flag_combos[ifl].name, rc, xx);
            }
        }
    }

    printf("\n],\n\"quirks\": [\n");
    first = 1;

    /* §2 Quirk cases: verify normalization behavior for edge-case ipls.
     *
     * 9099 (Sun COB): clause (iii) cancels → plain Sun, CENTER_BODY cleared.
     * 9201 (Mercury moon #01): parent < Mars, no file opened → plain Mercury,
     *   CENTER_BODY survives (inert flag).
     * 9499 (Mars COB): clause (iii) cancels → plain Mars, CENTER_BODY cleared.
     *
     * Each quirk row includes both the raw-ipl call AND the expected
     * plain-planet call for bitwise comparison. */
    struct { int ipl; int plain; const char *desc; } quirks[] = {
        { 9099, SE_SUN,     "9099_sun_cob_cancel" },
        { 9201, SE_MERCURY, "9201_mercury_inert" },
        { 9499, SE_MARS,    "9499_mars_cob_cancel" },
    };
    int nquirks = 3;
    double quirk_jd = 2451545.0; /* J2000 */
    int quirk_flags = SEFLG_MOSEPH | SEFLG_SPEED;

    for (int iq = 0; iq < nquirks; iq++) {
        /* The quirk ipl. */
        memset(xx, 0, sizeof(xx));
        int rc = swe_calc(quirk_jd, quirks[iq].ipl, quirk_flags, xx, serr);
        if (rc < 0) {
            fprintf(stderr, "ERROR quirk: ipl=%d: %s\n", quirks[iq].ipl, serr);
            return 1;
        }
        emit(quirks[iq].ipl, quirk_jd, quirk_flags, quirks[iq].desc, rc, xx);

        /* The plain planet for comparison. */
        double xx_plain[6];
        memset(xx_plain, 0, sizeof(xx_plain));
        int rc_plain = swe_calc(quirk_jd, quirks[iq].plain, quirk_flags, xx_plain, serr);
        if (rc_plain < 0) {
            fprintf(stderr, "ERROR plain: ipl=%d: %s\n", quirks[iq].plain, serr);
            return 1;
        }
        char plain_name[64];
        snprintf(plain_name, sizeof(plain_name), "%s_plain", quirks[iq].desc);
        emit(quirks[iq].plain, quirk_jd, quirk_flags, plain_name, rc_plain, xx_plain);
    }

    printf("\n],\n\"pheno\": [\n");
    first = 1;

    /* Pheno cases for plmoon bodies: 9501, 9599 at J2000 with MOSEPH. */
    int pheno_bodies[] = { 9501, 9599 };
    int npheno = 2;
    int pheno_flags = SEFLG_MOSEPH;

    for (int ip = 0; ip < npheno; ip++) {
        double attr[20];
        memset(attr, 0, sizeof(attr));
        int rc = swe_pheno(quirk_jd, pheno_bodies[ip], pheno_flags, attr, serr);
        if (rc < 0) {
            fprintf(stderr, "ERROR pheno: ipl=%d: %s\n", pheno_bodies[ip], serr);
            return 1;
        }
        if (!first) printf(",\n");
        first = 0;
        printf("  {\"ipl\": %d, \"jd\": %.17g, \"flags\": %d, \"retflag\": %d, "
               "\"attr\": [%.17g, %.17g, %.17g, %.17g, %.17g, %.17g]}",
               pheno_bodies[ip], quirk_jd, pheno_flags, rc,
               attr[0], attr[1], attr[2], attr[3], attr[4], attr[5]);
    }

    printf("\n],\n\"orbit\": [\n");
    first = 1;

    /* Orbital elements for 9501 at J2000 with MOSEPH|SPEED.
     * C produces NaN for some elements (degenerate near-parent heliocentric
     * orbit) — output null for NaN to keep JSON valid. */
    {
        double dret[50];
        memset(dret, 0, sizeof(dret));
        int orbit_flags = SEFLG_MOSEPH | SEFLG_SPEED;
        int rc = swe_get_orbital_elements(quirk_jd, 9501, orbit_flags, dret, serr);
        if (rc < 0) {
            fprintf(stderr, "ERROR orbit: ipl=9501: %s\n", serr);
            return 1;
        }
        printf("  {\"ipl\": 9501, \"jd\": %.17g, \"flags\": %d, \"retflag\": %d, "
               "\"dret\": [", quirk_jd, orbit_flags, rc);
        for (int i = 0; i < 17; i++) {
            if (isnan(dret[i]) || isinf(dret[i]))
                printf("null");
            else
                printf("%.17g", dret[i]);
            if (i < 16) printf(", ");
        }
        printf("]}");
    }

    printf("\n]}\n");
    swe_close();
    return 0;
}
