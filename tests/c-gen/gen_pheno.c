/*
 * Generates golden reference data for the phenomena module: swe_pheno
 * (phase angle, illuminated fraction, elongation, apparent diameter,
 * apparent magnitude, Moon horizontal parallax), swisseph-rs/83.
 *
 * Bodies: Sun..Pluto (10). These exercise every magnitude branch of swe_pheno
 * EXCEPT the Bowell H-G asteroid branch (§5k): Sun(5a) Moon(5b) Mercury(5c)
 * Venus(5d) Mars(5e) Jupiter(5f) Saturn(5g) Uranus(5h) Neptune(5i) Pluto(5j
 * generic). Ceres/Pallas/Juno/Vesta (which would exercise §5k via mag_elem)
 * are intentionally OMITTED: no ephemeris backend in the Rust port can compute
 * an asteroid position yet (Moshier has no asteroids; the SE1 sweph backend maps
 * Ceres.. -> None; JPL de441 has no main-belt asteroids). The §5k branch is
 * implemented but stays golden-uncovered until asteroid SE1 position support
 * lands -- the same underlying gap as eclipse.rs::body_radius_au's asteroid stub.
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

int main(void) {
    char serr[256];
    swe_set_ephe_path("../swisseph/ephe");
    int first = 1;
    printf("[\n");
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
                if (!first) printf(",\n");
                first = 0;
                printf("  {\"tjd_et\": %.20e, \"ipl\": %d, \"body_name\": \"%s\", "
                       "\"iflag\": %d, \"flag_name\": \"%s\", \"retflag\": %d, "
                       "\"attr\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                       epochs[ie], bodies[ib], body_names[ib],
                       flags, flag_combos[ifl].name, rc,
                       attr[0], attr[1], attr[2], attr[3], attr[4], attr[5]);
            }
        }
    }
    printf("\n]\n");
    return 0;
}
