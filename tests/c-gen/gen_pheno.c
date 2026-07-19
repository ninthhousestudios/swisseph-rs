/*
 * Generates golden reference data for the phenomena module: swe_pheno
 * (phase angle, illuminated fraction, elongation, apparent diameter,
 * apparent magnitude, Moon horizontal parallax), swisseph-rs/83 + /103.
 *
 * Bodies: Sun..Pluto (10) exercise magnitude branches 5a-5j. Main asteroids
 * Chiron/Pholus/Ceres/Pallas/Juno/Vesta exercise §5k via mag_elem (SWIEPH
 * only — MOSEPH asteroid pheno would read C's global sun_bary cache that an
 * earlier SWIEPH call populated, making the result process-state-dependent;
 * c-ref-asteroid.md §1.5). Numbered asteroids 433/7066/136199 exercise §5k
 * via the SE1 orbital-element H/G metadata.
 *
 * Compile (from repo root):
 *   cc -O2 -I../swisseph -o tests/c-gen/gen_pheno tests/c-gen/gen_pheno.c \
 *      -L../swisseph -lswe -lm
 * Run:
 *   ./tests/c-gen/gen_pheno > tests/golden-data/pheno.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static int bodies[] = {
    SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS,
    SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO
};
static const char *body_names[] = {
    "Sun", "Moon", "Mercury", "Venus", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"
};
#define NBODIES 10

/* Main asteroids (use MAG_ELEM, SWIEPH only). */
static int main_ast_bodies[] = {
    SE_CHIRON, SE_PHOLUS, SE_CERES, SE_PALLAS, SE_JUNO, SE_VESTA
};
static const char *main_ast_names[] = {
    "Chiron", "Pholus", "Ceres", "Pallas", "Juno", "Vesta"
};
#define NMAIN_AST 6

/* Numbered asteroids (use SE1 orbital-element H/G, SWIEPH only). */
static int numbered_ast_bodies[] = {
    SE_AST_OFFSET + 433, SE_AST_OFFSET + 7066, SE_AST_OFFSET + 136199
};
static const char *numbered_ast_names[] = {
    "Eros_433", "Nessus_7066", "Eris_136199"
};
#define NNUMBERED_AST 3

/* Includes J2000 and a pre-1900 epoch. NB: the pre-1900 epoch is nudged 4 days off gen_calc.c's
 * exact 1800-Jan-1 (2378496.5), which is sepl_18's tfstart: at that exact file boundary C's
 * swe_pheno reuses the .se1 file its own internal planet call cached, so its elongation differs
 * from a clean "calc-the-Sun-fresh" elongation by ~5e-8 -- a documented stateless-vs-stateful
 * file-boundary artifact (CLAUDE.md "Stateless vs Stateful"), not a phenomena bug. Sampling a few
 * days into the file removes the ambiguity so the whole battery holds 1e-9. */
static double epochs[] = {
    2451545.0,    /* J2000.0 */
    2460310.5,    /* 2024-Jan-1 */
    2433282.5,    /* 1950-Jan-1 */
    2378500.5,    /* 1800-Jan-5 (pre-1900, off the sepl_18 file boundary) */
};
#define NEPOCHS 4

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_MOSEPH,                 "MOSEPH" },
    { SEFLG_MOSEPH | SEFLG_TRUEPOS, "MOSEPH_TRUEPOS" },
    { SEFLG_SWIEPH,                 "SWIEPH" },
};
#define NFLAGS 3

static void emit_case(int *first, double tjd, int ipl, const char *name,
                       int flags, const char *flagname, double *attr, int rc) {
    if (!*first) printf(",\n");
    *first = 0;
    printf("  {\"tjd_et\": %.20e, \"ipl\": %d, \"body_name\": \"%s\", "
           "\"iflag\": %d, \"flag_name\": \"%s\", \"retflag\": %d, "
           "\"attr\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
           tjd, ipl, name, flags, flagname, rc,
           attr[0], attr[1], attr[2], attr[3], attr[4], attr[5]);
}

int main(void) {
    char serr[256];
    swe_set_ephe_path("../../ephe");
    int first = 1;
    printf("[\n");

    /* Original planet battery: Sun..Pluto × 4 epochs × 3 flags = 120 cases. */
    for (int ib = 0; ib < NBODIES; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            for (int ifl = 0; ifl < NFLAGS; ifl++) {
                int flags = flag_combos[ifl].flag;
                double attr[20];
                memset(attr, 0, sizeof(attr));
                int rc = swe_pheno(epochs[ie], bodies[ib], flags, attr, serr);
                if (rc < 0) {
                    fprintf(stderr, "skipping: %s body=%d jd=%.1f flags=0x%x\n",
                            serr, bodies[ib], epochs[ie], flags);
                    continue;
                }
                emit_case(&first, epochs[ie], bodies[ib], body_names[ib],
                          flags, flag_combos[ifl].name, attr, rc);
            }
        }
    }

    /* Main asteroids: SWIEPH only (MOSEPH asteroid pheno depends on C's global
     * sun_bary cache state — c-ref-asteroid.md §1.5). 6 × 4 = 24 cases. */
    for (int ib = 0; ib < NMAIN_AST; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            int flags = SEFLG_SWIEPH;
            double attr[20];
            memset(attr, 0, sizeof(attr));
            int rc = swe_pheno(epochs[ie], main_ast_bodies[ib], flags, attr, serr);
            if (rc < 0) {
                fprintf(stderr, "skipping: %s body=%d jd=%.1f flags=0x%x\n",
                        serr, main_ast_bodies[ib], epochs[ie], flags);
                continue;
            }
            emit_case(&first, epochs[ie], main_ast_bodies[ib], main_ast_names[ib],
                      flags, "SWIEPH", attr, rc);
        }
    }

    /* Numbered asteroids: SWIEPH only. 3 × 4 = 12 cases. */
    for (int ib = 0; ib < NNUMBERED_AST; ib++) {
        for (int ie = 0; ie < NEPOCHS; ie++) {
            int flags = SEFLG_SWIEPH;
            double attr[20];
            memset(attr, 0, sizeof(attr));
            int rc = swe_pheno(epochs[ie], numbered_ast_bodies[ib], flags, attr, serr);
            if (rc < 0) {
                fprintf(stderr, "skipping: %s body=%d jd=%.1f flags=0x%x\n",
                        serr, numbered_ast_bodies[ib], epochs[ie], flags);
                continue;
            }
            emit_case(&first, epochs[ie], numbered_ast_bodies[ib], numbered_ast_names[ib],
                      flags, "SWIEPH", attr, rc);
        }
    }

    printf("\n]\n");
    return 0;
}
