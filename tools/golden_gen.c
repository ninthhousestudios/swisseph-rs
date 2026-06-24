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
#include "swephlib.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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

/* ==========================================================================
 * math module
 * ========================================================================== */

/* --------------------------------------------------------------------------
 * degnorm / radnorm
 * -------------------------------------------------------------------------- */

static void emit_degnorm(double x) {
    arr_sep();
    printf("{\"input\":%.20g,\"output\":%.20g}", x, swe_degnorm(x));
}

static void emit_radnorm(double x) {
    arr_sep();
    printf("{\"input\":%.20g,\"output\":%.20g}", x, swe_radnorm(x));
}

static void gen_degnorm(void) {
    printf("\"degnorm\":[");
    arr_start();
    double vals[] = {0, 360, -360, 720, -720, 90, -90, 180, -180, 270,
                     1e-14, -1e-14, 1e-13, -1e-13, 1e-12,
                     359.9999999999999, 360.0000000000001,
                     1e12, -1e12, 1e15, -1e15,
                     0.1, 45.5, 123.456789, 359.999, 360.001};
    int n = sizeof(vals) / sizeof(vals[0]);
    for (int i = 0; i < n; i++) emit_degnorm(vals[i]);
    printf("]");
}

static void gen_radnorm(void) {
    printf("\"radnorm\":[");
    arr_start();
    double pi = M_PI;
    double twopi = 2.0 * pi;
    double vals[] = {0, twopi, -twopi, 2*twopi, -2*twopi, pi, -pi,
                     pi/2, -pi/2, 3*pi/2,
                     1e-14, -1e-14, 1e-13, -1e-13,
                     twopi - 1e-15, twopi + 1e-15,
                     1e12, -1e12, 0.1, 1.0, 2.5};
    int n = sizeof(vals) / sizeof(vals[0]);
    for (int i = 0; i < n; i++) emit_radnorm(vals[i]);
    printf("]");
}

/* --------------------------------------------------------------------------
 * diff / midpoint
 * -------------------------------------------------------------------------- */

static void emit_difdeg2n(double p1, double p2) {
    arr_sep();
    printf("{\"p1\":%.20g,\"p2\":%.20g,\"output\":%.20g}", p1, p2, swe_difdeg2n(p1, p2));
}

static void emit_difdegn(double p1, double p2) {
    arr_sep();
    printf("{\"p1\":%.20g,\"p2\":%.20g,\"output\":%.20g}", p1, p2, swe_difdegn(p1, p2));
}

static void emit_difrad2n(double p1, double p2) {
    arr_sep();
    printf("{\"p1\":%.20g,\"p2\":%.20g,\"output\":%.20g}", p1, p2, swe_difrad2n(p1, p2));
}

static void emit_midp_deg(double x1, double x0) {
    arr_sep();
    printf("{\"x1\":%.20g,\"x0\":%.20g,\"output\":%.20g}", x1, x0, swe_deg_midp(x1, x0));
}

static void emit_midp_rad(double x1, double x0) {
    arr_sep();
    printf("{\"x1\":%.20g,\"x0\":%.20g,\"output\":%.20g}", x1, x0, swe_rad_midp(x1, x0));
}

static void gen_diff_midpoint(void) {
    double pairs[][2] = {
        {0, 180}, {180, 0}, {350, 10}, {10, 350}, {0, 0}, {180, 180},
        {90, 270}, {270, 90}, {1, 359}, {359, 1}, {0.001, 359.999},
        {45, 315}, {315, 45}, {179.999, 180.001}, {720, 360}
    };
    int np = sizeof(pairs) / sizeof(pairs[0]);

    printf("\"difdeg2n\":[");
    arr_start();
    for (int i = 0; i < np; i++) emit_difdeg2n(pairs[i][0], pairs[i][1]);
    printf("]");

    printf(",\"difdegn\":[");
    arr_start();
    for (int i = 0; i < np; i++) emit_difdegn(pairs[i][0], pairs[i][1]);
    printf("]");

    printf(",\"difrad2n\":[");
    arr_start();
    double pi = M_PI;
    double rpairs[][2] = {
        {0, pi}, {pi, 0}, {6.0, 0.2}, {0.2, 6.0}, {0, 0},
        {pi/2, 3*pi/2}, {1.0, 5.0}, {5.0, 1.0}
    };
    int nrp = sizeof(rpairs) / sizeof(rpairs[0]);
    for (int i = 0; i < nrp; i++) emit_difrad2n(rpairs[i][0], rpairs[i][1]);
    printf("]");

    printf(",\"midp_deg\":[");
    arr_start();
    for (int i = 0; i < np; i++) emit_midp_deg(pairs[i][0], pairs[i][1]);
    printf("]");

    printf(",\"midp_rad\":[");
    arr_start();
    for (int i = 0; i < nrp; i++) emit_midp_rad(rpairs[i][0], rpairs[i][1]);
    printf("]");
}

/* --------------------------------------------------------------------------
 * centisecond functions
 * -------------------------------------------------------------------------- */

static void emit_csnorm(int32 p) {
    arr_sep();
    printf("{\"input\":%d,\"output\":%d}", p, swe_csnorm(p));
}

static void emit_difcsn(int32 p1, int32 p2) {
    arr_sep();
    printf("{\"p1\":%d,\"p2\":%d,\"output\":%d}", p1, p2, swe_difcsn(p1, p2));
}

static void emit_difcs2n(int32 p1, int32 p2) {
    arr_sep();
    printf("{\"p1\":%d,\"p2\":%d,\"output\":%d}", p1, p2, swe_difcs2n(p1, p2));
}

static void gen_centisec(void) {
    int32 csvals[] = {0, 1, -1, 129600000, -129600000, 64800000, -64800000,
                      129599999, 129600001, -129599999, 259200000, -259200000,
                      100000, -100000, 360000, 10800000};
    int ncs = sizeof(csvals) / sizeof(csvals[0]);

    printf("\"csnorm\":[");
    arr_start();
    for (int i = 0; i < ncs; i++) emit_csnorm(csvals[i]);
    printf("]");

    int32 cspairs[][2] = {
        {0, 64800000}, {64800000, 0}, {129599999, 1}, {1, 129599999},
        {0, 0}, {100000, 50000}, {50000, 100000},
        {-100000, 100000}, {100000, -100000}
    };
    int ncp = sizeof(cspairs) / sizeof(cspairs[0]);

    printf(",\"difcsn\":[");
    arr_start();
    for (int i = 0; i < ncp; i++) emit_difcsn(cspairs[i][0], cspairs[i][1]);
    printf("]");

    printf(",\"difcs2n\":[");
    arr_start();
    for (int i = 0; i < ncp; i++) emit_difcs2n(cspairs[i][0], cspairs[i][1]);
    printf("]");
}

/* --------------------------------------------------------------------------
 * d2l
 * -------------------------------------------------------------------------- */

static void emit_d2l(double x) {
    arr_sep();
    printf("{\"input\":%.20g,\"output\":%d}", x, swe_d2l(x));
}

static void gen_d2l(void) {
    printf("\"d2l\":[");
    arr_start();
    double vals[] = {0.0, 0.5, -0.5, 1.0, -1.0, 1.5, -1.5, 0.4999999, -0.4999999,
                     0.50000001, -0.50000001, 100.7, -100.7, 1e6, -1e6};
    int n = sizeof(vals) / sizeof(vals[0]);
    for (int i = 0; i < n; i++) emit_d2l(vals[i]);
    printf("]");
}

/* --------------------------------------------------------------------------
 * Chebyshev
 * -------------------------------------------------------------------------- */

static void emit_echeb(double t, double *coef, int ncf, const char *label) {
    arr_sep();
    printf("{\"t\":%.20g,\"label\":\"%s\",\"ncf\":%d,\"coef\":[", t, label, ncf);
    for (int i = 0; i < ncf; i++) {
        if (i > 0) printf(",");
        printf("%.20g", coef[i]);
    }
    printf("],\"value\":%.20g,\"deriv\":%.20g}",
           swi_echeb(t, coef, ncf), swi_edcheb(t, coef, ncf));
}

static void gen_chebyshev(void) {
    printf("\"chebyshev\":[");
    arr_start();

    /* T0 = 1: Broucke coeffs [2.0] (halved convention) */
    double c_t0[] = {2.0};
    double ts[] = {0.0, 0.5, -0.5, 1.0, -1.0, 0.3, -0.7};
    int nt = sizeof(ts) / sizeof(ts[0]);
    for (int i = 0; i < nt; i++) emit_echeb(ts[i], c_t0, 1, "T0");

    /* T1 = t: coeffs [0, 1] */
    double c_t1[] = {0.0, 1.0};
    for (int i = 0; i < nt; i++) emit_echeb(ts[i], c_t1, 2, "T1");

    /* T2 = 2t^2-1: coeffs [0, 0, 1] */
    double c_t2[] = {0.0, 0.0, 1.0};
    for (int i = 0; i < nt; i++) emit_echeb(ts[i], c_t2, 3, "T2");

    /* Realistic: 5 coefficients */
    double c5[] = {1.5, -0.3, 0.7, -0.1, 0.05};
    for (int i = 0; i < nt; i++) emit_echeb(ts[i], c5, 5, "poly5");

    /* 8 coefficients */
    double c8[] = {2.1, 0.5, -0.8, 0.3, -0.15, 0.07, -0.03, 0.01};
    for (int i = 0; i < nt; i++) emit_echeb(ts[i], c8, 8, "poly8");

    printf("]");
}

/* --------------------------------------------------------------------------
 * Coordinate transforms
 * -------------------------------------------------------------------------- */

static void emit_cartpol(double x, double y, double z) {
    double inp[3] = {x, y, z};
    double out[3];
    swi_cartpol(inp, out);
    arr_sep();
    printf("{\"x\":%.20g,\"y\":%.20g,\"z\":%.20g,"
           "\"lon\":%.20g,\"lat\":%.20g,\"dist\":%.20g}",
           x, y, z, out[0], out[1], out[2]);
}

static void emit_polcart(double lon, double lat, double dist) {
    double inp[3] = {lon, lat, dist};
    double out[3];
    swi_polcart(inp, out);
    arr_sep();
    printf("{\"lon\":%.20g,\"lat\":%.20g,\"dist\":%.20g,"
           "\"x\":%.20g,\"y\":%.20g,\"z\":%.20g}",
           lon, lat, dist, out[0], out[1], out[2]);
}

static void emit_cartpol_sp(double x0, double x1, double x2,
                             double x3, double x4, double x5) {
    double inp[6] = {x0, x1, x2, x3, x4, x5};
    double out[6];
    swi_cartpol_sp(inp, out);
    arr_sep();
    printf("{\"xi\":[%.20g,%.20g,%.20g,%.20g,%.20g,%.20g],"
           "\"xo\":[%.20g,%.20g,%.20g,%.20g,%.20g,%.20g]}",
           x0, x1, x2, x3, x4, x5,
           out[0], out[1], out[2], out[3], out[4], out[5]);
}

static void emit_polcart_sp(double l0, double l1, double l2,
                             double l3, double l4, double l5) {
    double inp[6] = {l0, l1, l2, l3, l4, l5};
    double out[6];
    swi_polcart_sp(inp, out);
    arr_sep();
    printf("{\"li\":[%.20g,%.20g,%.20g,%.20g,%.20g,%.20g],"
           "\"xo\":[%.20g,%.20g,%.20g,%.20g,%.20g,%.20g]}",
           l0, l1, l2, l3, l4, l5,
           out[0], out[1], out[2], out[3], out[4], out[5]);
}

static void emit_coortrf(double x0, double x1, double x2, double eps) {
    double inp[3] = {x0, x1, x2};
    double out[3];
    swi_coortrf(inp, out, eps);
    arr_sep();
    printf("{\"xi\":[%.20g,%.20g,%.20g],\"eps\":%.20g,"
           "\"xo\":[%.20g,%.20g,%.20g]}",
           x0, x1, x2, eps, out[0], out[1], out[2]);
}

static void emit_cotrans(double lon, double lat, double dist, double eps) {
    double inp[3] = {lon, lat, dist};
    double out[3];
    swe_cotrans(inp, out, eps);
    arr_sep();
    printf("{\"xi\":[%.20g,%.20g,%.20g],\"eps\":%.20g,"
           "\"xo\":[%.20g,%.20g,%.20g]}",
           lon, lat, dist, eps, out[0], out[1], out[2]);
}

static void emit_cotrans_sp(double x0, double x1, double x2,
                             double x3, double x4, double x5, double eps) {
    double inp[6] = {x0, x1, x2, x3, x4, x5};
    double out[6];
    swe_cotrans_sp(inp, out, eps);
    arr_sep();
    printf("{\"xi\":[%.20g,%.20g,%.20g,%.20g,%.20g,%.20g],\"eps\":%.20g,"
           "\"xo\":[%.20g,%.20g,%.20g,%.20g,%.20g,%.20g]}",
           x0, x1, x2, x3, x4, x5, eps,
           out[0], out[1], out[2], out[3], out[4], out[5]);
}

static void gen_coords(void) {
    double pi = M_PI;

    /* --- cartpol --- */
    printf("\"cartpol\":[");
    arr_start();
    emit_cartpol(1, 0, 0);
    emit_cartpol(0, 1, 0);
    emit_cartpol(0, 0, 1);
    emit_cartpol(0, 0, -1);
    emit_cartpol(0, 0, 0);
    emit_cartpol(1, 1, 1);
    emit_cartpol(-1, -1, -1);
    emit_cartpol(1, 2, 3);
    emit_cartpol(-3, 2, -1);
    emit_cartpol(0.001, 0, 0.001);
    emit_cartpol(1e-15, 1e-15, 1e-15);
    emit_cartpol(1e10, 1e10, 1e10);
    printf("]");

    /* --- polcart --- */
    printf(",\"polcart\":[");
    arr_start();
    emit_polcart(0, 0, 1);
    emit_polcart(pi/2, 0, 1);
    emit_polcart(pi, 0, 1);
    emit_polcart(0, pi/2, 1);
    emit_polcart(0, -pi/2, 1);
    emit_polcart(pi/4, pi/4, 2);
    emit_polcart(1.0, 0.5, 3.0);
    emit_polcart(5.0, -0.3, 0.5);
    emit_polcart(0, 0, 0);
    printf("]");

    /* --- coortrf --- */
    printf(",\"coortrf\":[");
    arr_start();
    double obliquity = 23.4393 * pi / 180.0;
    emit_coortrf(1, 0, 0, 0);
    emit_coortrf(0, 1, 0, 0);
    emit_coortrf(0, 0, 1, 0);
    emit_coortrf(1, 2, 3, obliquity);
    emit_coortrf(1, 2, 3, -obliquity);
    emit_coortrf(0.5, -0.3, 0.7, pi/6);
    emit_coortrf(0.5, -0.3, 0.7, -pi/6);
    emit_coortrf(1, 1, 1, pi/4);
    emit_coortrf(1, 1, 1, pi/2);
    printf("]");

    /* --- cartpol_sp --- */
    printf(",\"cartpol_sp\":[");
    arr_start();
    emit_cartpol_sp(1, 2, 3, 0.1, 0.2, 0.3);
    emit_cartpol_sp(0, 0, 0, 1, 0, 0);
    emit_cartpol_sp(0, 0, 0, 0, 1, 0);
    emit_cartpol_sp(0, 0, 0, 0, 0, 1);
    emit_cartpol_sp(0, 0, 0, 0, 0, 0);
    emit_cartpol_sp(1, 0, 0, 0, 0, 0);
    emit_cartpol_sp(0, 0, 1, 0, 0, 0);
    emit_cartpol_sp(1, 1, 0, 0.5, -0.5, 0.1);
    emit_cartpol_sp(-1, 2, -3, 0.01, -0.02, 0.03);
    emit_cartpol_sp(1e-5, 1e-5, 1e-5, 1e-6, 1e-6, 1e-6);
    printf("]");

    /* --- polcart_sp --- */
    printf(",\"polcart_sp\":[");
    arr_start();
    emit_polcart_sp(1.0, 0.5, 2.0, 0.1, 0.05, 0.3);
    emit_polcart_sp(0, 0, 1, 0, 0, 0);
    emit_polcart_sp(pi, 0, 1, 0, 0, 0);
    emit_polcart_sp(0, 0, 0, 0.1, 0.2, 0.3);
    emit_polcart_sp(pi/4, pi/6, 3.0, 0.01, -0.02, 0.5);
    emit_polcart_sp(5.0, -0.3, 0.5, -0.1, 0.05, -0.2);
    printf("]");

    /* --- cotrans --- */
    printf(",\"cotrans\":[");
    arr_start();
    double eps_vals[] = {23.4393, -23.4393, 0, 45, -45, 90};
    int neps = sizeof(eps_vals) / sizeof(eps_vals[0]);
    double coords[][3] = {
        {0, 0, 1}, {90, 0, 1}, {180, 0, 1}, {270, 0, 1},
        {0, 45, 1}, {0, -45, 1}, {0, 90, 1}, {0, -90, 1},
        {45, 30, 2.5}, {123.456, -23.789, 1.5},
        {350, 10, 0.5}, {359.999, 0.001, 1}
    };
    int nc = sizeof(coords) / sizeof(coords[0]);
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < neps; j++)
            emit_cotrans(coords[i][0], coords[i][1], coords[i][2], eps_vals[j]);
    printf("]");

    /* --- cotrans_sp --- */
    printf(",\"cotrans_sp\":[");
    arr_start();
    emit_cotrans_sp(45, 30, 2.5, 1.0, -0.5, 0.3, 23.4393);
    emit_cotrans_sp(45, 30, 2.5, 1.0, -0.5, 0.3, -23.4393);
    emit_cotrans_sp(0, 0, 1, 0, 0, 0, 23.4393);
    emit_cotrans_sp(180, 0, 1, 0.1, 0.2, 0.3, 23.4393);
    emit_cotrans_sp(0, 90, 1, 1.0, 0.0, 0.0, 23.4393);
    emit_cotrans_sp(123.456, -23.789, 1.5, 0.5, -0.3, 0.1, 23.4393);
    emit_cotrans_sp(350, 10, 0.5, -0.1, 0.05, -0.02, -23.4393);
    printf("]");
}

/* --------------------------------------------------------------------------
 * split_deg
 * -------------------------------------------------------------------------- */

static void emit_split_deg(double ddeg, int32 flags) {
    int32 ideg, imin, isec, isgn;
    double dsecfr;
    swe_split_deg(ddeg, flags, &ideg, &imin, &isec, &dsecfr, &isgn);
    arr_sep();
    printf("{\"ddeg\":%.20g,\"flags\":%d,"
           "\"deg\":%d,\"min\":%d,\"sec\":%d,\"secfr\":%.20g,\"sign\":%d}",
           ddeg, flags, ideg, imin, isec, dsecfr, isgn);
}

static void gen_split_deg(void) {
    printf("\"split_deg\":[");
    arr_start();

    double degs[] = {0, 123.456789, 359.9999999, 29.9999, 30.0, 30.0001,
                     -45.5, -123.456789, 0.000001, 89.99999999, 360.0};
    int nd = sizeof(degs) / sizeof(degs[0]);

    /* SE_SPLIT_DEG flags (without NAKSHATRA) */
    int32 flag_combos[] = {
        0,      /* no flags */
        1,      /* ROUND_SEC */
        2,      /* ROUND_MIN */
        4,      /* ROUND_DEG */
        8,      /* ZODIACAL */
        8|16,   /* ZODIACAL|KEEP_SIGN */
        8|32,   /* ZODIACAL|KEEP_DEG */
        1|8,    /* ROUND_SEC|ZODIACAL */
        2|8,    /* ROUND_MIN|ZODIACAL */
        4|8,    /* ROUND_DEG|ZODIACAL */
        1|8|16, /* ROUND_SEC|ZODIACAL|KEEP_SIGN */
        4|32,   /* ROUND_DEG|KEEP_DEG */
    };
    int nf = sizeof(flag_combos) / sizeof(flag_combos[0]);

    for (int i = 0; i < nd; i++)
        for (int j = 0; j < nf; j++)
            emit_split_deg(degs[i], flag_combos[j]);

    printf("]");
}

/* --------------------------------------------------------------------------
 * math module entry
 * -------------------------------------------------------------------------- */

static void gen_math(void) {
    printf("{");
    gen_degnorm();
    printf(",");
    gen_radnorm();
    printf(",");
    gen_diff_midpoint();
    printf(",");
    gen_centisec();
    printf(",");
    gen_d2l();
    printf(",");
    gen_chebyshev();
    printf(",");
    gen_coords();
    printf(",");
    gen_split_deg();
    printf("}\n");
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
        fprintf(stderr, "Usage: %s <module>\n  module: date | math\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "date") == 0) {
        gen_date();
    } else if (strcmp(argv[1], "math") == 0) {
        gen_math();
    } else {
        fprintf(stderr, "Unknown module: %s\n", argv[1]);
        return 1;
    }
    return 0;
}
