#include "swisseph.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int test_version(void) {
    const char *v = swisseph_version();
    if (v == NULL) {
        fprintf(stderr, "FAIL: swisseph_version returned NULL\n");
        return 1;
    }
    if (strlen(v) == 0) {
        fprintf(stderr, "FAIL: swisseph_version returned empty string\n");
        return 1;
    }
    printf("  version: %s\n", v);
    return 0;
}

static int test_config_default(void) {
    SweConfig config;
    swisseph_config_default(&config);
    if (config.ephemeris_source != 0) {
        fprintf(stderr, "FAIL: default ephemeris_source should be 0 (Moshier), got %d\n",
                config.ephemeris_source);
        return 1;
    }
    if (!isnan(config.tidal_acceleration)) {
        fprintf(stderr, "FAIL: default tidal_acceleration should be NAN\n");
        return 1;
    }
    if (!isnan(config.delta_t_userdef)) {
        fprintf(stderr, "FAIL: default delta_t_userdef should be NAN\n");
        return 1;
    }
    return 0;
}

static int test_new_calc_free(void) {
    SweConfig config;
    swisseph_config_default(&config);

    SweEphemeris *handle = NULL;
    char err[256] = {0};
    int32_t ret = swisseph_new(&config, &handle, err, sizeof(err));
    if (ret != 0) {
        fprintf(stderr, "FAIL: swisseph_new returned %d: %s\n", ret, err);
        return 1;
    }
    if (handle == NULL) {
        fprintf(stderr, "FAIL: swisseph_new returned OK but handle is NULL\n");
        return 1;
    }

    double xx[6] = {0};
    int32_t flags_used = 0;
    memset(err, 0, sizeof(err));
    ret = swisseph_calc_ut(
        handle,
        2451545.0,   /* J2000.0 */
        0,           /* Sun */
        256,         /* SEFLG_SPEED */
        NULL,        /* no geopos */
        NULL,        /* no sid_mode */
        xx,
        &flags_used,
        err,
        sizeof(err)
    );
    if (ret != 0) {
        fprintf(stderr, "FAIL: swisseph_calc_ut returned %d: %s\n", ret, err);
        swisseph_free(handle);
        return 1;
    }

    /* Sun longitude at J2000 should be ~280 degrees */
    if (xx[0] < 270.0 || xx[0] > 290.0) {
        fprintf(stderr, "FAIL: Sun longitude %.6f out of plausible range [270, 290]\n", xx[0]);
        swisseph_free(handle);
        return 1;
    }
    printf("  Sun@J2000 lon=%.6f lat=%.6f dist=%.6f\n", xx[0], xx[1], xx[2]);

    swisseph_free(handle);
    return 0;
}

static int test_invalid_body(void) {
    SweConfig config;
    swisseph_config_default(&config);

    SweEphemeris *handle = NULL;
    char err[256] = {0};
    swisseph_new(&config, &handle, err, sizeof(err));

    double xx[6] = {0};
    memset(err, 0, sizeof(err));
    int32_t ret = swisseph_calc_ut(
        handle,
        2451545.0,
        -999,        /* invalid body */
        256,
        NULL,
        NULL,
        xx,
        NULL,
        err,
        sizeof(err)
    );
    if (ret >= 0) {
        fprintf(stderr, "FAIL: invalid body should return negative, got %d\n", ret);
        swisseph_free(handle);
        return 1;
    }
    if (strlen(err) == 0) {
        fprintf(stderr, "FAIL: err_buf should be non-empty for invalid body\n");
        swisseph_free(handle);
        return 1;
    }
    printf("  invalid body error: %s (code %d)\n", err, ret);

    swisseph_free(handle);
    return 0;
}

static int test_free_null(void) {
    /* Must not crash */
    swisseph_free(NULL);
    return 0;
}

static int test_share(void) {
    SweConfig config;
    swisseph_config_default(&config);

    SweEphemeris *handle = NULL;
    char err[256] = {0};
    int32_t ret = swisseph_new(&config, &handle, err, sizeof(err));
    if (ret != 0) {
        fprintf(stderr, "FAIL: swisseph_new returned %d: %s\n", ret, err);
        return 1;
    }

    /* Share the handle */
    SweEphemeris *shared = NULL;
    memset(err, 0, sizeof(err));
    ret = swisseph_share(handle, &shared, err, sizeof(err));
    if (ret != 0) {
        fprintf(stderr, "FAIL: swisseph_share returned %d: %s\n", ret, err);
        swisseph_free(handle);
        return 1;
    }
    if (shared == NULL) {
        fprintf(stderr, "FAIL: swisseph_share returned OK but shared is NULL\n");
        swisseph_free(handle);
        return 1;
    }

    /* Free original first */
    swisseph_free(handle);

    /* Shared handle must still work */
    double xx[6] = {0};
    int32_t flags_used = 0;
    memset(err, 0, sizeof(err));
    ret = swisseph_calc_ut(shared, 2451545.0, 0, 256, NULL, NULL, xx, &flags_used, err, sizeof(err));
    if (ret != 0) {
        fprintf(stderr, "FAIL: calc_ut on shared handle returned %d: %s\n", ret, err);
        swisseph_free(shared);
        return 1;
    }
    if (xx[0] < 270.0 || xx[0] > 290.0) {
        fprintf(stderr, "FAIL: shared Sun longitude %.6f out of range\n", xx[0]);
        swisseph_free(shared);
        return 1;
    }
    printf("  shared Sun@J2000 lon=%.6f\n", xx[0]);

    swisseph_free(shared);
    return 0;
}

static int test_julday_revjul(void) {
    double jd = swisseph_julday(2000, 1, 1, 12.0, 1);
    if (fabs(jd - 2451545.0) > 1e-10) {
        fprintf(stderr, "FAIL: julday(2000,1,1,12) = %.10f, expected 2451545.0\n", jd);
        return 1;
    }

    int32_t y, m, d;
    double h;
    swisseph_revjul(jd, 1, &y, &m, &d, &h);
    if (y != 2000 || m != 1 || d != 1 || fabs(h - 12.0) > 1e-10) {
        fprintf(stderr, "FAIL: revjul roundtrip: %d-%d-%d %.2f\n", y, m, d, h);
        return 1;
    }
    return 0;
}

static int test_degnorm(void) {
    double r = swisseph_degnorm(-10.0);
    if (fabs(r - 350.0) > 1e-12) {
        fprintf(stderr, "FAIL: degnorm(-10) = %.10f, expected 350.0\n", r);
        return 1;
    }
    return 0;
}

int main(void) {
    int failures = 0;
    printf("C smoke tests:\n");

    printf("  test_version...\n");
    failures += test_version();

    printf("  test_config_default...\n");
    failures += test_config_default();

    printf("  test_new_calc_free...\n");
    failures += test_new_calc_free();

    printf("  test_invalid_body...\n");
    failures += test_invalid_body();

    printf("  test_free_null...\n");
    failures += test_free_null();

    printf("  test_julday_revjul...\n");
    failures += test_julday_revjul();

    printf("  test_degnorm...\n");
    failures += test_degnorm();

    printf("  test_share...\n");
    failures += test_share();

    if (failures > 0) {
        fprintf(stderr, "\n%d test(s) FAILED\n", failures);
        return 1;
    }
    printf("\nAll C smoke tests passed.\n");
    return 0;
}
