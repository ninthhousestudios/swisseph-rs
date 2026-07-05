/*
 * C timing harness for swisseph-rs benchmark comparison.
 *
 * Same four core workloads as benches/ephemeris.rs:
 *   1. calc_moshier: 10 planets at J2000, Moshier, with speed
 *   2. calc_swiss:   10 planets at J2000, Swiss Ephemeris files, with speed
 *   3. houses_placidus: Placidus houses at Zurich, J2000
 *   4. eclipse_search: next solar eclipse from J2000
 *
 * Note: C numbers are single-thread only — libswe is not thread-safe.
 * That asymmetry IS the point of the thread-scaling capstone in Rust.
 *
 * Build: make -C tools bench_c
 * Run:   tools/bench_c
 */

#include "swephexp.h"
#include <stdio.h>
#include <time.h>

#define J2000 2451545.0

static double elapsed_ns(struct timespec *start, struct timespec *end) {
    return (double)(end->tv_sec - start->tv_sec) * 1e9 +
           (double)(end->tv_nsec - start->tv_nsec);
}

static long bench_calc_moshier(void) {
    double xx[6];
    char serr[256];
    int bodies[] = {SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS,
                    SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO};
    int nbodies = 10;
    int iflag = SEFLG_SPEED | SEFLG_MOSEPH;

    struct timespec start, end;
    long iters = 0;
    double total = 0.0;

    /* warm up */
    for (int i = 0; i < nbodies; i++)
        swe_calc_ut(J2000, bodies[i], iflag, xx, serr);

    clock_gettime(CLOCK_MONOTONIC, &start);
    while (total < 3e9) { /* run for >= 3 seconds */
        for (int i = 0; i < nbodies; i++)
            swe_calc_ut(J2000, bodies[i], iflag, xx, serr);
        iters++;
        clock_gettime(CLOCK_MONOTONIC, &end);
        total = elapsed_ns(&start, &end);
    }

    long ns_per_op = (long)(total / (double)iters);
    printf("calc_moshier_full_chart: %ld ns/iter (%ld iters)\n", ns_per_op, iters);
    return ns_per_op;
}

static long bench_calc_swiss(void) {
    double xx[6];
    char serr[256];
    int bodies[] = {SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS,
                    SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO};
    int nbodies = 10;
    int iflag = SEFLG_SPEED | SEFLG_SWIEPH;

    swe_set_ephe_path("../swisseph/ephe");

    /* test if files are available */
    int rc = swe_calc_ut(J2000, SE_SUN, iflag, xx, serr);
    if (rc < 0) {
        printf("calc_swiss_full_chart: SKIP (files not found)\n");
        return 0;
    }

    struct timespec start, end;
    long iters = 0;
    double total = 0.0;

    clock_gettime(CLOCK_MONOTONIC, &start);
    while (total < 3e9) {
        for (int i = 0; i < nbodies; i++)
            swe_calc_ut(J2000, bodies[i], iflag, xx, serr);
        iters++;
        clock_gettime(CLOCK_MONOTONIC, &end);
        total = elapsed_ns(&start, &end);
    }

    long ns_per_op = (long)(total / (double)iters);
    printf("calc_swiss_full_chart:   %ld ns/iter (%ld iters)\n", ns_per_op, iters);
    return ns_per_op;
}

static long bench_houses_placidus(void) {
    double cusps[13], ascmc[10];

    /* warm up */
    swe_houses(J2000, 47.37, 8.55, 'P', cusps, ascmc);

    struct timespec start, end;
    long iters = 0;
    double total = 0.0;

    clock_gettime(CLOCK_MONOTONIC, &start);
    while (total < 3e9) {
        swe_houses(J2000, 47.37, 8.55, 'P', cusps, ascmc);
        iters++;
        clock_gettime(CLOCK_MONOTONIC, &end);
        total = elapsed_ns(&start, &end);
    }

    long ns_per_op = (long)(total / (double)iters);
    printf("houses_placidus:         %ld ns/iter (%ld iters)\n", ns_per_op, iters);
    return ns_per_op;
}

static long bench_eclipse_search(void) {
    double tret[10];
    char serr[256];
    int iflag = SEFLG_MOSEPH;

    /* warm up */
    swe_sol_eclipse_when_glob(J2000, iflag, 0, tret, 0, serr);

    struct timespec start, end;
    long iters = 0;
    double total = 0.0;

    clock_gettime(CLOCK_MONOTONIC, &start);
    while (total < 3e9) {
        swe_sol_eclipse_when_glob(J2000, iflag, 0, tret, 0, serr);
        iters++;
        clock_gettime(CLOCK_MONOTONIC, &end);
        total = elapsed_ns(&start, &end);
    }

    long ns_per_op = (long)(total / (double)iters);
    printf("sol_eclipse_when_glob:   %ld ns/iter (%ld iters)\n", ns_per_op, iters);
    return ns_per_op;
}

int main(void) {
    printf("C Swiss Ephemeris benchmark (single-thread, libswe.a)\n");
    printf("------------------------------------------------------\n");

    bench_calc_moshier();
    bench_calc_swiss();
    bench_houses_placidus();
    bench_eclipse_search();

    swe_close();
    return 0;
}
