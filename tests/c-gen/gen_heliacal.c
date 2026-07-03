/*
 * Generates golden reference data for the heliacal public API: swe_vis_limit_mag.
 * This harness links libswe.a normally (does NOT #include swehel.c).
 *
 * Objects: venus, sirius, moon, mercury. Observer: Cairo (31.25°E, 30.1°N, 30m).
 * Per-object UT instants chosen so the object is above the horizon in the desired
 * sky-brightness regime. Two dates: 2005-Jan-01 (JD 2453371) and 2005-Jan-25
 * (JD 2453395, near full moon for Moon-above-horizon cases).
 *
 * Compile (from repo root):
 *   cd tests/c-gen && cc -O2 -I../../../swisseph -o gen_heliacal \
 *      gen_heliacal.c ../../../swisseph/libswe.a -lm
 * Run (from repo root):
 *   ./tests/c-gen/gen_heliacal > tests/golden-data/heliacal.json
 */

#include <stdio.h>
#include <string.h>
#include "swephexp.h"

/* Cairo observer */
static double dgeo[3] = {31.25, 30.1, 30.0};
static double datm_default[4] = {1013.25, 15.0, 40.0, 40.0};
static double dobs_default[6] = {36.0, 1.0, 0.0, 0.0, 0.0, 0.0};

static char serr_buf[256];

static void emit_case(int *first, double tjd, const char *objname,
                       int helflag, const char *flag_desc, int epheflag,
                       double *dret, int retval) {
    if (!*first) printf(",\n");
    *first = 0;
    printf("  {\"tjd_ut\": %.20e, \"object\": \"%s\", "
           "\"helflag\": %d, \"flag_desc\": \"%s\", \"epheflag\": %d, "
           "\"retval\": %d, "
           "\"dret\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
           tjd, objname, helflag, flag_desc, epheflag, retval,
           dret[0], dret[1], dret[2], dret[3], dret[4], dret[5], dret[6], dret[7]);
}

static void run_case(int *first, double tjd, const char *obj,
                     int helflag, const char *flag_desc, int epheflag) {
    double datm[4], dobs[6], dret[8];
    char objname[256];
    memcpy(datm, datm_default, sizeof(datm));
    memcpy(dobs, dobs_default, sizeof(dobs));
    memset(dret, 0, sizeof(dret));
    strncpy(objname, obj, sizeof(objname)-1);
    objname[sizeof(objname)-1] = '\0';
    int retval = swe_vis_limit_mag(tjd, dgeo, datm, dobs,
                                    objname, helflag, dret, serr_buf);
    if (retval == ERR) {
        fprintf(stderr, "ERR: %s obj=%s jd=%.4f helflag=%d\n",
                serr_buf, obj, tjd, helflag);
        return;
    }
    emit_case(first, tjd, obj, helflag, flag_desc, epheflag, dret, retval);
}

static void run_case_optic(int *first, double tjd, const char *obj,
                           int helflag, const char *flag_desc, int epheflag,
                           double *dobs_custom) {
    double datm[4], dobs[6], dret[8];
    char objname[256];
    memcpy(datm, datm_default, sizeof(datm));
    memcpy(dobs, dobs_custom, 6 * sizeof(double));
    memset(dret, 0, sizeof(dret));
    strncpy(objname, obj, sizeof(objname)-1);
    objname[sizeof(objname)-1] = '\0';
    int retval = swe_vis_limit_mag(tjd, dgeo, datm, dobs,
                                    objname, helflag, dret, serr_buf);
    if (retval == ERR) {
        fprintf(stderr, "ERR: %s obj=%s jd=%.4f helflag=%d\n",
                serr_buf, obj, tjd, helflag);
        return;
    }
    emit_case(first, tjd, obj, helflag, flag_desc, epheflag, dret, retval);
}

int main(void) {
    swe_set_ephe_path("../swisseph/ephe");
    int first = 1;
    printf("{\"vis_limit\": [\n");

    /* ── Main battery: per-object UT instants × SWIEPH ────────────
     * Date 1: 2005-Jan-01 (JD 2453371)
     *   Venus/Mercury: above horizon during daytime (UT 04:00-12:00)
     *   Sirius: above horizon at night (UT 18:00-04:00)
     * Date 2: 2005-Jan-25 (JD 2453395, near full moon 2005-Jan-25)
     *   Moon: above horizon at night
     */

    /* Venus: daytime cases (Sun up, Venus up) */
    double venus_uts[] = {
        2453371.0,      /* 12:00 UT — high Sun, Bday regime */
        2453370.6667,   /* 04:00 UT — morning, near sunrise */
        2453370.75,     /* 06:00 UT — morning, Sun rising */
        2453370.8333,   /* 08:00 UT — morning, Sun well up */
        2453371.5,      /* 00:00 UT Jan 2 — nighttime (Venus below horizon → -2) */
    };
    for (int i = 0; i < 5; i++)
        run_case(&first, venus_uts[i], "venus", SEFLG_SWIEPH, "SWIEPH", SEFLG_SWIEPH);

    /* Sirius: nighttime cases (star visible at night) */
    double sirius_uts[] = {
        2453371.25,     /* 18:00 UT — deep twilight, Sirius up */
        2453371.4167,   /* 22:00 UT — full night */
        2453371.5,      /* 00:00 UT Jan 2 — late night */
        2453371.625,    /* 03:00 UT Jan 2 — pre-dawn */
        2453371.0,      /* 12:00 UT — daytime (Sirius below horizon → -2) */
    };
    for (int i = 0; i < 5; i++)
        run_case(&first, sirius_uts[i], "sirius", SEFLG_SWIEPH, "SWIEPH", SEFLG_SWIEPH);

    /* Moon: near full moon, visible at night */
    double moon_uts[] = {
        2453395.4167,   /* Jan 25 22:00 UT — full moon, night */
        2453395.5,      /* Jan 26 00:00 UT — Moon high at night */
        2453395.25,     /* Jan 25 18:00 UT — evening, Moon rising */
        2453371.4167,   /* Jan 01 22:00 UT — Moon may be up */
        2453371.0,      /* Jan 01 12:00 UT — daytime (Moon may be below → -2) */
    };
    for (int i = 0; i < 5; i++)
        run_case(&first, moon_uts[i], "moon", SEFLG_SWIEPH, "SWIEPH", SEFLG_SWIEPH);

    /* Mercury: daytime cases */
    double mercury_uts[] = {
        2453371.0,      /* 12:00 UT — daytime */
        2453370.75,     /* 06:00 UT — morning */
        2453370.8333,   /* 08:00 UT — morning */
        2453371.4167,   /* 22:00 UT — nighttime (Mercury below horizon → -2) */
    };
    for (int i = 0; i < 4; i++)
        run_case(&first, mercury_uts[i], "mercury", SEFLG_SWIEPH, "SWIEPH", SEFLG_SWIEPH);

    /* ── VISLIM_DARK: force Sun/Moon to -90° ──────────────────── */
    run_case(&first, 2453371.4167, "sirius",
             SEFLG_SWIEPH | SE_HELFLAG_VISLIM_DARK, "SWIEPH_VISLIM_DARK", SEFLG_SWIEPH);
    run_case(&first, 2453371.0, "venus",
             SEFLG_SWIEPH | SE_HELFLAG_VISLIM_DARK, "SWIEPH_VISLIM_DARK", SEFLG_SWIEPH);

    /* ── VISLIM_NOMOON: suppress Moon contribution ───────────── */
    run_case(&first, 2453371.4167, "sirius",
             SEFLG_SWIEPH | SE_HELFLAG_VISLIM_NOMOON, "SWIEPH_VISLIM_NOMOON", SEFLG_SWIEPH);
    run_case(&first, 2453395.4167, "mercury",
             SEFLG_SWIEPH | SE_HELFLAG_VISLIM_NOMOON, "SWIEPH_VISLIM_NOMOON", SEFLG_SWIEPH);

    /* ── OPTICAL_PARAMS: custom dobs (binocular telescope) ───── */
    {
        double dobs_optic[6] = {36.0, 1.0, 1.0, 10.0, 50.0, 0.8};
        run_case_optic(&first, 2453371.4167, "sirius",
                       SEFLG_SWIEPH | SE_HELFLAG_OPTICAL_PARAMS, "SWIEPH_OPTICAL_PARAMS",
                       SEFLG_SWIEPH, dobs_optic);
        run_case_optic(&first, 2453371.0, "venus",
                       SEFLG_SWIEPH | SE_HELFLAG_OPTICAL_PARAMS, "SWIEPH_OPTICAL_PARAMS",
                       SEFLG_SWIEPH, dobs_optic);
    }

    /* ── Forced VISLIM_PHOTOPIC / VISLIM_SCOTOPIC ─────────────── */
    run_case(&first, 2453371.4167, "sirius",
             SEFLG_SWIEPH | SE_HELFLAG_VISLIM_PHOTOPIC, "SWIEPH_VISLIM_PHOTOPIC", SEFLG_SWIEPH);
    run_case(&first, 2453371.4167, "sirius",
             SEFLG_SWIEPH | SE_HELFLAG_VISLIM_SCOTOPIC, "SWIEPH_VISLIM_SCOTOPIC", SEFLG_SWIEPH);

    /* ── MOSEPH duplicates (planets only, no stars) ───────────── */
    run_case(&first, 2453371.0, "venus", SEFLG_MOSEPH, "MOSEPH", SEFLG_MOSEPH);
    run_case(&first, 2453370.75, "venus", SEFLG_MOSEPH, "MOSEPH", SEFLG_MOSEPH);
    run_case(&first, 2453371.4167, "moon", SEFLG_MOSEPH, "MOSEPH", SEFLG_MOSEPH);
    run_case(&first, 2453395.4167, "moon", SEFLG_MOSEPH, "MOSEPH", SEFLG_MOSEPH);
    run_case(&first, 2453371.0, "mercury", SEFLG_MOSEPH, "MOSEPH", SEFLG_MOSEPH);
    run_case(&first, 2453370.8333, "mercury", SEFLG_MOSEPH, "MOSEPH", SEFLG_MOSEPH);

    printf("\n]}\n");
    return 0;
}
