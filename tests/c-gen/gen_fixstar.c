/*
 * Generates golden reference data for swe_fixstar2 (Moshier backend).
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../../swisseph -o gen_fixstar gen_fixstar.c \
 *      ../../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_fixstar > ../golden-data/fixstar.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static const char *stars[] = {
    "Sirius",
    "Spica",
    "Aldebaran",
    "Regulus",
    ",SgrA*",
    ",GPol",
    ",GP1958",
};
#define NSTARS 7

static double epochs[] = {
    2451545.0,   /* J2000.0 */
    2415020.0,   /* J1900.0 */
    2469807.5,   /* 2050 Jan 1 */
    2086302.5,   /* ~1000 AD */
};
#define NEPOCHS 4

struct flag_combo {
    int flag;
    const char *name;
};

static struct flag_combo flag_combos[] = {
    { SEFLG_MOSEPH,                           "base" },
    { SEFLG_MOSEPH | SEFLG_SPEED,             "SPEED" },
    { SEFLG_MOSEPH | SEFLG_EQUATORIAL,        "EQUATORIAL" },
    { SEFLG_MOSEPH | SEFLG_TRUEPOS,           "TRUEPOS" },
    { SEFLG_MOSEPH | SEFLG_J2000,             "J2000" },
    { SEFLG_MOSEPH | SEFLG_NONUT,             "NONUT" },
    { SEFLG_MOSEPH | SEFLG_XYZ,              "XYZ" },
};
#define NFLAGS 7

int main(void) {
    double xx[6];
    char serr[256];
    char star_buf[512];
    int retflag;
    int i, j, k;
    int first;

    swe_set_ephe_path("../../../swisseph/ephe");

    printf("{\n");

    /* ---- position cases ---- */
    printf("  \"fixstar\": [\n");
    first = 1;
    for (i = 0; i < NSTARS; i++) {
        for (j = 0; j < NEPOCHS; j++) {
            for (k = 0; k < NFLAGS; k++) {
                strncpy(star_buf, stars[i], sizeof(star_buf) - 1);
                star_buf[sizeof(star_buf) - 1] = '\0';

                retflag = swe_fixstar2(star_buf, epochs[j], flag_combos[k].flag,
                                       xx, serr);
                if (!first) printf(",\n");
                first = 0;

                printf("    {\"star\": \"%s\", \"tjd\": %.1f, \"iflag\": %d,"
                       " \"flag_name\": \"%s\","
                       " \"xx\": [%.17g, %.17g, %.17g, %.17g, %.17g, %.17g],"
                       " \"retflag\": %d,"
                       " \"star_out\": \"%s\"}",
                       stars[i], epochs[j], flag_combos[k].flag,
                       flag_combos[k].name,
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5],
                       retflag, star_buf);
            }
        }
    }
    printf("\n  ],\n");

    /* ---- magnitude cases (catalog only — no builtins via fixstar2_mag) ---- */
    printf("  \"mag\": [\n");
    {
        const char *mag_stars[] = { "Sirius", "Aldebaran", "Regulus", "Spica" };
        int nmag = 4;
        first = 1;
        for (i = 0; i < nmag; i++) {
            double mag = 0.0;
            strncpy(star_buf, mag_stars[i], sizeof(star_buf) - 1);
            star_buf[sizeof(star_buf) - 1] = '\0';
            retflag = swe_fixstar2_mag(star_buf, &mag, serr);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"star\": \"%s\", \"mag\": %.17g, \"retflag\": %d,"
                   " \"star_out\": \"%s\"}",
                   mag_stars[i], mag, retflag, star_buf);
        }
        printf("\n");
    }
    printf("  ]\n");

    printf("}\n");
    swe_close();
    return 0;
}
