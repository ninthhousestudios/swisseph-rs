/*
 * Generates golden reference data for utility functions (swisseph-rs/124):
 * swe_time_equ, swe_lmt_to_lat, swe_lat_to_lmt, swe_get_planet_name,
 * swe_house_name, swe_get_ayanamsa_name, swe_cs2timestr, swe_cs2lonlatstr,
 * swe_cs2degstr, swe_csroundsec.
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_utilities gen_utilities.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_utilities > ../golden-data/utilities.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

static void print_str(const char *key, const char *val) {
    printf("\"%s\": ", key);
    if (val == NULL) {
        printf("null");
    } else {
        printf("\"");
        for (const char *p = val; *p; p++) {
            if (*p == '"') printf("\\\"");
            else if (*p == '\\') printf("\\\\");
            else putchar(*p);
        }
        printf("\"");
    }
}

int main(void) {
    char serr[256];
    swe_set_ephe_path("../../ephe");
    int first = 1;
    printf("{\n");

    /* ---- time_equ ---- */
    {
        double epochs[] = {
            2451545.0,  /* J2000 */
            2415020.0,  /* J1900 */
            2460600.5,  /* 2024-Oct */
            2305447.5,  /* 1600-Jan-1 */
            2488069.5,  /* ~2100-Jan-1 */
            2378496.5,  /* ~1800-Jan-1 */
            2436204.5,  /* 1958-Jan-1 */
            2440587.5,  /* 1970-Jan-1 */
            2444239.5,  /* 1980-Jan-1 */
            2448622.5,  /* 1992-Jan-1 */
            2453371.5,  /* 2005-Jan-1 */
            2459580.5,  /* 2022-Jan-1 */
        };
        int n = sizeof(epochs)/sizeof(epochs[0]);
        printf("  \"time_equ\": [\n");
        for (int i = 0; i < n; i++) {
            double E;
            memset(serr, 0, sizeof(serr));
            int rc = swe_time_equ(epochs[i], &E, serr);
            if (i > 0) printf(",\n");
            printf("    {\"jd_ut\": %.20e, \"E\": %.20e, \"rc\": %d}", epochs[i], E, rc);
        }
        printf("\n  ],\n");
    }

    /* ---- lmt_to_lat / lat_to_lmt ---- */
    {
        struct { double jd; double geolon; } pairs[] = {
            {2451545.0,   8.55},
            {2451545.0, -71.0},
            {2460600.5, 139.75},
            {2415020.0, -122.4},
            {2305447.5,  0.0},
            {2488069.5, 180.0},
        };
        int n = sizeof(pairs)/sizeof(pairs[0]);
        printf("  \"lmt_to_lat\": [\n");
        for (int i = 0; i < n; i++) {
            double tjd_lat;
            memset(serr, 0, sizeof(serr));
            int rc = swe_lmt_to_lat(pairs[i].jd, pairs[i].geolon, &tjd_lat, serr);
            if (i > 0) printf(",\n");
            printf("    {\"jd_lmt\": %.20e, \"geolon\": %.20e, \"tjd_lat\": %.20e, \"rc\": %d}",
                   pairs[i].jd, pairs[i].geolon, tjd_lat, rc);
        }
        printf("\n  ],\n");

        printf("  \"lat_to_lmt\": [\n");
        for (int i = 0; i < n; i++) {
            double tjd_lmt;
            memset(serr, 0, sizeof(serr));
            int rc = swe_lat_to_lmt(pairs[i].jd, pairs[i].geolon, &tjd_lmt, serr);
            if (i > 0) printf(",\n");
            printf("    {\"jd_lat\": %.20e, \"geolon\": %.20e, \"tjd_lmt\": %.20e, \"rc\": %d}",
                   pairs[i].jd, pairs[i].geolon, tjd_lmt, rc);
        }
        printf("\n  ],\n");
    }

    /* ---- get_planet_name ---- */
    {
        int ipls[] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,       /* Sun..Pluto */
            10, 11, 12, 13, 14, 15,               /* MeanNode..IntpPerigee */
            16,                                    /* Earth (SE_ECL_NUT actually, but let's use 16) */
            17, 18, 19, 20, 21, 22,               /* Ceres..IntpPerigee */
            40, 41, 42, 43, 44, 45, 46, 47,       /* Fictitious: Cupido..Poseidon */
            48, 49, 50, 51, 52, 53, 54,           /* Fictitious: Isis..Pickering */
            55, 56, 57, 58,                        /* Vulcan..Waldemath */
            10000+433, 10000+5, 10000+7066,       /* Asteroids */
        };
        int n = sizeof(ipls)/sizeof(ipls[0]);
        printf("  \"planet_name\": [\n");
        for (int i = 0; i < n; i++) {
            char name[256];
            swe_get_planet_name(ipls[i], name);
            if (i > 0) printf(",\n");
            printf("    {\"ipl\": %d, ", ipls[i]);
            print_str("name", name);
            printf("}");
        }
        printf("\n  ],\n");
    }

    /* ---- house_name ---- */
    {
        char systems[] = "ABCDEFGHIiJKLMNOPQRSTUVWXY";
        int n = strlen(systems);
        printf("  \"house_name\": [\n");
        for (int i = 0; i < n; i++) {
            const char *name = swe_house_name(systems[i]);
            if (i > 0) printf(",\n");
            printf("    {\"hsys\": \"%c\", ", systems[i]);
            print_str("name", name);
            printf("}");
        }
        printf("\n  ],\n");
    }

    /* ---- ayanamsa_name ---- */
    {
        printf("  \"ayanamsa_name\": [\n");
        for (int i = 0; i <= 46; i++) {
            const char *name = swe_get_ayanamsa_name(i);
            if (i > 0) printf(",\n");
            printf("    {\"sidm\": %d, ", i);
            print_str("name", name);
            printf("}");
        }
        /* SE_SIDM_USER = 255 */
        printf(",\n    {\"sidm\": 255, ");
        print_str("name", swe_get_ayanamsa_name(255));
        printf("}");
        printf("\n  ],\n");
    }

    /* ---- centisecond formatters ---- */
    {
        int32 cs_values[] = {
            0,
            100,           /* 1 arcsecond */
            149,
            150,
            6000,          /* 1 arcminute */
            360000,        /* 1 degree */
            4500000,       /* 12h30m = 12*360000 + 30*6000 */
            4500000 + 4500,/* 12h30m45s */
            10800000 - 50, /* just below 30-degree boundary */
            10800000,      /* exact 30-degree boundary */
            10800000 + 50, /* just above 30-degree boundary */
            129600000 - 50,/* just below 360-degree */
            129600000,     /* exact 360-degree */
            1,             /* sub-second */
            99,
            43956789,      /* arbitrary value: 122d05m47.89s */
        };
        int n = sizeof(cs_values)/sizeof(cs_values[0]);
        printf("  \"csroundsec\": [\n");
        for (int i = 0; i < n; i++) {
            if (i > 0) printf(",\n");
            printf("    {\"input\": %d, \"output\": %d}", cs_values[i], swe_csroundsec(cs_values[i]));
        }
        printf("\n  ],\n");

        printf("  \"cs2timestr\": [\n");
        first = 1;
        for (int i = 0; i < n; i++) {
            for (int suppress = 0; suppress <= 1; suppress++) {
                char buf[80];
                swe_cs2timestr(cs_values[i], ':', suppress, buf);
                if (!first) printf(",\n");
                first = 0;
                printf("    {\"input\": %d, \"suppress_zero\": %s, ", cs_values[i], suppress ? "true" : "false");
                print_str("output", buf);
                printf("}");
            }
        }
        printf("\n  ],\n");

        printf("  \"cs2lonlatstr\": [\n");
        first = 1;
        for (int i = 0; i < n; i++) {
            char buf[80];
            swe_cs2lonlatstr(cs_values[i], 'E', 'W', buf);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"input\": %d, \"pchar\": \"E\", \"mchar\": \"W\", ", cs_values[i]);
            print_str("output", buf);
            printf("}");
        }
        printf("\n  ],\n");

        printf("  \"cs2degstr\": [\n");
        first = 1;
        for (int i = 0; i < n; i++) {
            char buf[80];
            swe_cs2degstr(cs_values[i], buf);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\"input\": %d, ", cs_values[i]);
            print_str("output", buf);
            printf("}");
        }
        printf("\n  ]\n");
    }

    printf("}\n");
    swe_close();
    return 0;
}
