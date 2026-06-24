/*
 * Golden data generator for swisseph.rs differential tests.
 *
 * Links against the C libswe.a and generates JSON golden data by calling
 * C Swiss Ephemeris functions with curated edge-case inputs.
 *
 * Usage: ./golden_gen <module>
 *   module = "date"  (more modules added as Rust implementation progresses)
 */

#include "swephexp.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* JSON output helpers — track whether we need a leading comma */
static int first_item;

static void arr_start(void) { first_item = 1; }

static void arr_sep(void) {
    if (!first_item) printf(",");
    first_item = 0;
}

/* --------------------------------------------------------------------------
 * julday cases
 * -------------------------------------------------------------------------- */

static void emit_julday(int y, int m, int d, double h, int g) {
    double jd = swe_julday(y, m, d, h, g);
    arr_sep();
    printf("{\"y\":%d,\"m\":%d,\"d\":%d,\"h\":%.20g,\"g\":%d,\"jd\":%.20g}",
           y, m, d, h, g, jd);
}

static void gen_julday(void) {
    printf("\"julday\":[");
    arr_start();

    /* --- Branch: month < 3 (Jan, Feb) vs >= 3, both calendars --- */
    int months[] = {1, 2, 3, 6, 10, 12};
    int nmonths = sizeof(months) / sizeof(months[0]);
    for (int mi = 0; mi < nmonths; mi++) {
        for (int g = 0; g <= 1; g++) {
            emit_julday(2000, months[mi], 15, 12.0, g);
        }
    }

    /* --- Branch: negative years, year 0, positive years --- */
    int years[] = {-5000, -4713, -1000, -500, -400, -200, -100, -4, -1,
                   0, 1, 4, 100, 200, 400, 500, 1000, 1582, 1583, 1900,
                   1970, 2000, 2024, 3000, 5000};
    int nyears = sizeof(years) / sizeof(years[0]);
    for (int yi = 0; yi < nyears; yi++) {
        for (int g = 0; g <= 1; g++) {
            emit_julday(years[yi], 1, 1, 0.0, g);
            emit_julday(years[yi], 6, 15, 12.0, g);
        }
    }

    /* --- Branch: century correction for negative years ---
     * Triggers: u < 0 && u/100 == floor(u/100) && u/400 != floor(u/400)
     * i.e., negative century years not divisible by 400 */
    int century_years[] = {-100, -200, -300, -500, -700, -900, -1100};
    int ncentury = sizeof(century_years) / sizeof(century_years[0]);
    for (int ci = 0; ci < ncentury; ci++) {
        emit_julday(century_years[ci], 3, 1, 0.0, 1);  /* Gregorian only */
    }
    /* Negative century years divisible by 400 (should NOT trigger correction) */
    int century400[] = {-400, -800, -1200, -2000};
    int ncent400 = sizeof(century400) / sizeof(century400[0]);
    for (int ci = 0; ci < ncent400; ci++) {
        emit_julday(century400[ci], 3, 1, 0.0, 1);
    }

    /* --- Gregorian/Julian switch boundary --- */
    emit_julday(1582, 10, 4, 12.0, 0);   /* last Julian day */
    emit_julday(1582, 10, 4, 12.0, 1);   /* same date, Gregorian */
    emit_julday(1582, 10, 15, 12.0, 0);  /* Julian */
    emit_julday(1582, 10, 15, 12.0, 1);  /* first Gregorian day */
    emit_julday(1582, 10, 10, 12.0, 0);  /* in the gap, Julian */
    emit_julday(1582, 10, 10, 12.0, 1);  /* in the gap, Gregorian */

    /* --- Leap year dates --- */
    int leap_years[] = {2000, 1900, 1600, 4, -4, 2024, 1996, 400, -400};
    int nleap = sizeof(leap_years) / sizeof(leap_years[0]);
    for (int li = 0; li < nleap; li++) {
        emit_julday(leap_years[li], 2, 28, 12.0, 1);
        emit_julday(leap_years[li], 2, 29, 12.0, 1);
        emit_julday(leap_years[li], 3, 1, 12.0, 1);
    }

    /* --- C setest's own vectors --- */
    emit_julday(19, 2, 1964, 22.5, 0);
    emit_julday(19, 2, 1964, 22.5, 1);
    emit_julday(19, 2, 1965, 22.5, 0);
    emit_julday(19, 2, 1965, 22.5, 1);

    /* --- Various hours --- */
    double hours[] = {0.0, 6.0, 12.0, 18.0, 22.5, 23.999999};
    int nhours = sizeof(hours) / sizeof(hours[0]);
    for (int hi = 0; hi < nhours; hi++) {
        emit_julday(2000, 1, 1, hours[hi], 1);
    }

    /* --- J2000 and Unix epoch --- */
    emit_julday(2000, 1, 1, 12.0, 1);   /* J2000.0 */
    emit_julday(1970, 1, 1, 0.0, 1);    /* Unix epoch */
    emit_julday(1858, 11, 17, 0.0, 1);  /* MJD epoch */

    /* --- Year 0 months sweep --- */
    for (int m = 1; m <= 12; m++) {
        emit_julday(0, m, 1, 0.0, 1);
    }

    printf("]");
}

/* --------------------------------------------------------------------------
 * revjul cases
 * -------------------------------------------------------------------------- */

static void emit_revjul(double jd, int g) {
    int y, m, d;
    double h;
    swe_revjul(jd, g, &y, &m, &d, &h);
    arr_sep();
    printf("{\"jd\":%.20g,\"g\":%d,\"y\":%d,\"m\":%d,\"d\":%d,\"h\":%.20g}",
           jd, g, y, m, d, h);
}

static void gen_revjul(void) {
    printf("\"revjul\":[");
    arr_start();

    /* --- Key JD values, both calendars --- */
    double jds[] = {
        0.0,                /* JD epoch */
        1.0,
        100000.0,
        500000.0,
        1000000.0,
        1721425.5,          /* Jan 1, 1 CE Julian */
        1830691.5,          /* revjul Gregorian correction boundary */
        1830692.5,          /* just past boundary */
        2000000.0,
        2299160.5,          /* Oct 15, 1582 Gregorian (switch day) */
        2299159.5,          /* Oct 4, 1582 Julian (last Julian day) */
        2299161.5,          /* Oct 16, 1582 */
        2341524.0,          /* C setest vector */
        2415020.0,          /* Jan 0.5, 1900 */
        2440587.5,          /* Unix epoch: Jan 1, 1970 */
        2451545.0,          /* J2000.0 */
        2451545.5,          /* J2000.0 + 0.5 */
        2500000.0,
    };
    int njds = sizeof(jds) / sizeof(jds[0]);
    for (int i = 0; i < njds; i++) {
        for (int g = 0; g <= 1; g++) {
            emit_revjul(jds[i], g);
        }
    }

    /* --- Negative JD values --- */
    double neg_jds[] = {-1.0, -100.0, -10000.0, -100000.0, -1000000.0};
    int nneg = sizeof(neg_jds) / sizeof(neg_jds[0]);
    for (int i = 0; i < nneg; i++) {
        for (int g = 0; g <= 1; g++) {
            emit_revjul(neg_jds[i], g);
        }
    }

    /* --- Fractional hours (mid-day, midnight, etc.) --- */
    emit_revjul(2451544.5, 1);   /* J2000 midnight */
    emit_revjul(2451545.25, 1);  /* J2000 + 6h */
    emit_revjul(2451545.75, 1);  /* J2000 + 18h */

    printf("]");
}

/* --------------------------------------------------------------------------
 * date_conversion cases
 * -------------------------------------------------------------------------- */

static void emit_date_conv(int y, int m, int d, double h, int g) {
    double tjd = 0;
    char cal = g ? 'g' : 'j';
    int rc = swe_date_conversion(y, m, d, h, cal, &tjd);
    int valid = (rc == 0);
    arr_sep();
    if (valid) {
        printf("{\"y\":%d,\"m\":%d,\"d\":%d,\"h\":%.20g,\"g\":%d,\"valid\":true,\"jd\":%.20g}",
               y, m, d, h, g, tjd);
    } else {
        printf("{\"y\":%d,\"m\":%d,\"d\":%d,\"h\":%.20g,\"g\":%d,\"valid\":false,\"jd\":null}",
               y, m, d, h, g);
    }
}

static void gen_date_conversion(void) {
    printf("\"date_conversion\":[");
    arr_start();

    /* --- Valid dates --- */
    emit_date_conv(2000, 1, 1, 0.0, 1);
    emit_date_conv(2000, 1, 1, 12.0, 1);
    emit_date_conv(2000, 2, 29, 0.0, 1);  /* leap year */
    emit_date_conv(1900, 2, 28, 0.0, 1);  /* non-leap */
    emit_date_conv(1582, 10, 15, 0.0, 1); /* first Gregorian */
    emit_date_conv(1582, 10, 4, 0.0, 0);  /* last Julian */
    emit_date_conv(-1, 1, 1, 0.0, 0);
    emit_date_conv(0, 1, 1, 0.0, 0);
    emit_date_conv(0, 1, 1, 0.0, 1);
    emit_date_conv(-5000, 6, 15, 0.0, 0);
    emit_date_conv(5000, 6, 15, 0.0, 1);

    /* --- Invalid dates --- */
    emit_date_conv(2000, 2, 30, 0.0, 1);  /* Feb 30 */
    emit_date_conv(1900, 2, 29, 0.0, 1);  /* Feb 29 non-leap */
    emit_date_conv(2000, 13, 1, 0.0, 1);  /* month 13 */
    emit_date_conv(2000, 0, 1, 0.0, 1);   /* month 0 */
    emit_date_conv(2000, 1, 0, 0.0, 1);   /* day 0 */
    emit_date_conv(2000, 1, 32, 0.0, 1);  /* day 32 */
    emit_date_conv(2000, 4, 31, 0.0, 1);  /* Apr 31 */
    emit_date_conv(1582, 10, 10, 0.0, 1); /* in Gregorian gap */
    emit_date_conv(1582, 10, 5, 0.0, 1);  /* in Gregorian gap */
    emit_date_conv(1582, 10, 14, 0.0, 1); /* last day in gap */

    printf("]");
}

/* --------------------------------------------------------------------------
 * day_of_week cases
 * -------------------------------------------------------------------------- */

static void emit_dow(double jd) {
    int dow = swe_day_of_week(jd);
    arr_sep();
    printf("{\"jd\":%.20g,\"dow\":%d}", jd, dow);
}

static void gen_day_of_week(void) {
    printf("\"day_of_week\":[");
    arr_start();

    /* Known reference dates:
     * 2000-01-01 (Sat=6), 2000-01-02 (Sun=0), 2000-01-03 (Mon=1)
     * 1970-01-01 (Thu=4), 2024-06-24 (Mon=1)
     * Full week from J2000 */
    double ref_jds[] = {
        2451544.5,  /* 2000-01-01 00:00 */
        2451545.5,  /* 2000-01-02 00:00 */
        2451546.5,  /* 2000-01-03 00:00 */
        2451547.5,  /* 2000-01-04 */
        2451548.5,  /* 2000-01-05 */
        2451549.5,  /* 2000-01-06 */
        2451550.5,  /* 2000-01-07 */
        2440587.5,  /* 1970-01-01 */
        2460486.5,  /* 2024-06-24 */
        2299160.5,  /* 1582-10-15 (Gregorian switch, Friday) */
        0.0,
        1.0,
        2451545.0,  /* J2000.0 noon */
    };
    int nref = sizeof(ref_jds) / sizeof(ref_jds[0]);
    for (int i = 0; i < nref; i++) {
        emit_dow(ref_jds[i]);
    }

    /* Sequence of 7 consecutive days from various epochs */
    double epochs[] = {-1000000.0, 0.0, 1000000.0, 2451545.0};
    int nepochs = sizeof(epochs) / sizeof(epochs[0]);
    for (int ei = 0; ei < nepochs; ei++) {
        for (int d = 0; d < 7; d++) {
            emit_dow(epochs[ei] + d);
        }
    }

    printf("]");
}

/* --------------------------------------------------------------------------
 * main
 * -------------------------------------------------------------------------- */

static void gen_date(void) {
    printf("{");
    gen_julday();
    printf(",");
    gen_revjul();
    printf(",");
    gen_date_conversion();
    printf(",");
    gen_day_of_week();
    printf("}\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <module>\n  module: date\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "date") == 0) {
        gen_date();
    } else {
        fprintf(stderr, "Unknown module: %s\n", argv[1]);
        return 1;
    }
    return 0;
}
