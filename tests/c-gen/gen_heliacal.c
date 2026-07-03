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

    printf("\n],\n");

    /* ── arcvis: swe_topo_arcus_visionis ────────────────────────────
     * Grid: mag × (azi_obj,alt_obj) × (azi_sun,azi_moon,alt_moon) × tjd_ut
     * Rotated (not full cross-product) to keep ~32 cases.
     */
    printf("\"arcvis\": [\n");
    first = 1;
    {
        double mags[] = {-4.6, -1.46, 0.5, 2.0};
        int nmags = 4;
        double obj_geom[][2] = {{120.0, 10.0}, {240.0, 3.0}};
        int nobj = 2;
        double sun_moon[][3] = {{90.0, 270.0, -10.0}, {300.0, 60.0, 20.0}};
        int nsm = 2;
        double tjds[] = {2453371.0, 2451545.0};
        int ntjd = 2;

        for (int im = 0; im < nmags; im++) {
            /* Rotate obj_geom/sun_moon/tjd indices */
            int io = im % nobj;
            int is = im % nsm;
            int it = im % ntjd;
            double datm[4], dobs[6], dret_val;
            char serr[256] = {0};
            memcpy(datm, datm_default, sizeof(datm));
            memcpy(dobs, dobs_default, sizeof(dobs));
            int retval = swe_topo_arcus_visionis(
                tjds[it], dgeo, datm, dobs,
                SEFLG_SWIEPH,
                mags[im],
                obj_geom[io][0], obj_geom[io][1],
                sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                &dret_val, serr);
            if (retval == ERR) {
                fprintf(stderr, "arcvis ERR: %s mag=%.2f\n", serr, mags[im]);
                continue;
            }
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"tjd_ut\": %.20e, \"mag\": %.20e, "
                   "\"azi_obj\": %.20e, \"alt_obj\": %.20e, "
                   "\"azi_sun\": %.20e, \"azi_moon\": %.20e, \"alt_moon\": %.20e, "
                   "\"helflag\": %d, \"retval\": %d, \"dret\": %.20e}",
                   tjds[it], mags[im],
                   obj_geom[io][0], obj_geom[io][1],
                   sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                   SEFLG_SWIEPH, retval, dret_val);
        }
        /* Additional cases: vary tjd and sun/moon geometry for each mag */
        for (int im = 0; im < nmags; im++) {
            int io = (im + 1) % nobj;
            int is = (im + 1) % nsm;
            int it = (im + 1) % ntjd;
            double datm[4], dobs[6], dret_val;
            char serr[256] = {0};
            memcpy(datm, datm_default, sizeof(datm));
            memcpy(dobs, dobs_default, sizeof(dobs));
            int retval = swe_topo_arcus_visionis(
                tjds[it], dgeo, datm, dobs,
                SEFLG_SWIEPH,
                mags[im],
                obj_geom[io][0], obj_geom[io][1],
                sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                &dret_val, serr);
            if (retval == ERR) {
                fprintf(stderr, "arcvis ERR: %s mag=%.2f\n", serr, mags[im]);
                continue;
            }
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"tjd_ut\": %.20e, \"mag\": %.20e, "
                   "\"azi_obj\": %.20e, \"alt_obj\": %.20e, "
                   "\"azi_sun\": %.20e, \"azi_moon\": %.20e, \"alt_moon\": %.20e, "
                   "\"helflag\": %d, \"retval\": %d, \"dret\": %.20e}",
                   tjds[it], mags[im],
                   obj_geom[io][0], obj_geom[io][1],
                   sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                   SEFLG_SWIEPH, retval, dret_val);
        }
        /* Full cross: mag[0] × all combos for broader coverage */
        for (int io = 0; io < nobj; io++)
            for (int is = 0; is < nsm; is++)
                for (int it = 0; it < ntjd; it++) {
                    double datm[4], dobs[6], dret_val;
                    char serr[256] = {0};
                    memcpy(datm, datm_default, sizeof(datm));
                    memcpy(dobs, dobs_default, sizeof(dobs));
                    int retval = swe_topo_arcus_visionis(
                        tjds[it], dgeo, datm, dobs,
                        SEFLG_SWIEPH,
                        mags[0],
                        obj_geom[io][0], obj_geom[io][1],
                        sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                        &dret_val, serr);
                    if (retval == ERR) {
                        fprintf(stderr, "arcvis ERR: %s\n", serr);
                        continue;
                    }
                    if (!first) printf(",\n");
                    first = 0;
                    printf("  {\"tjd_ut\": %.20e, \"mag\": %.20e, "
                           "\"azi_obj\": %.20e, \"alt_obj\": %.20e, "
                           "\"azi_sun\": %.20e, \"azi_moon\": %.20e, \"alt_moon\": %.20e, "
                           "\"helflag\": %d, \"retval\": %d, \"dret\": %.20e}",
                           tjds[it], mags[0],
                           obj_geom[io][0], obj_geom[io][1],
                           sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                           SEFLG_SWIEPH, retval, dret_val);
                }
    }
    printf("\n],\n");

    /* ── helangle: swe_heliacal_angle ──────────────────────────────
     * Same grid minus alt_obj. Output: dret[0..2].
     */
    printf("\"helangle\": [\n");
    first = 1;
    {
        double mags[] = {-4.6, -1.46, 0.5, 2.0};
        int nmags = 4;
        double azi_objs[] = {120.0, 240.0};
        int nazi = 2;
        double sun_moon[][3] = {{90.0, 270.0, -10.0}, {300.0, 60.0, 20.0}};
        int nsm = 2;
        double tjds[] = {2453371.0, 2451545.0};
        int ntjd = 2;

        for (int im = 0; im < nmags; im++) {
            int ia = im % nazi;
            int is = im % nsm;
            int it = im % ntjd;
            double datm[4], dobs[6], dret[3];
            char serr[256] = {0};
            memcpy(datm, datm_default, sizeof(datm));
            memcpy(dobs, dobs_default, sizeof(dobs));
            memset(dret, 0, sizeof(dret));
            int retval = swe_heliacal_angle(
                tjds[it], dgeo, datm, dobs,
                SEFLG_SWIEPH,
                mags[im],
                azi_objs[ia],
                sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                dret, serr);
            if (retval == ERR) {
                fprintf(stderr, "helangle ERR: %s mag=%.2f\n", serr, mags[im]);
                continue;
            }
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"tjd_ut\": %.20e, \"mag\": %.20e, "
                   "\"azi_obj\": %.20e, "
                   "\"azi_sun\": %.20e, \"azi_moon\": %.20e, \"alt_moon\": %.20e, "
                   "\"helflag\": %d, \"retval\": %d, "
                   "\"dret\": [%.20e, %.20e, %.20e]}",
                   tjds[it], mags[im],
                   azi_objs[ia],
                   sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                   SEFLG_SWIEPH, retval,
                   dret[0], dret[1], dret[2]);
        }
        /* Second rotation */
        for (int im = 0; im < nmags; im++) {
            int ia = (im + 1) % nazi;
            int is = (im + 1) % nsm;
            int it = (im + 1) % ntjd;
            double datm[4], dobs[6], dret[3];
            char serr[256] = {0};
            memcpy(datm, datm_default, sizeof(datm));
            memcpy(dobs, dobs_default, sizeof(dobs));
            memset(dret, 0, sizeof(dret));
            int retval = swe_heliacal_angle(
                tjds[it], dgeo, datm, dobs,
                SEFLG_SWIEPH,
                mags[im],
                azi_objs[ia],
                sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                dret, serr);
            if (retval == ERR) {
                fprintf(stderr, "helangle ERR: %s mag=%.2f\n", serr, mags[im]);
                continue;
            }
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"tjd_ut\": %.20e, \"mag\": %.20e, "
                   "\"azi_obj\": %.20e, "
                   "\"azi_sun\": %.20e, \"azi_moon\": %.20e, \"alt_moon\": %.20e, "
                   "\"helflag\": %d, \"retval\": %d, "
                   "\"dret\": [%.20e, %.20e, %.20e]}",
                   tjds[it], mags[im],
                   azi_objs[ia],
                   sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                   SEFLG_SWIEPH, retval,
                   dret[0], dret[1], dret[2]);
        }
        /* Full cross for mag[0] */
        for (int ia = 0; ia < nazi; ia++)
            for (int is = 0; is < nsm; is++)
                for (int it = 0; it < ntjd; it++) {
                    double datm[4], dobs[6], dret[3];
                    char serr[256] = {0};
                    memcpy(datm, datm_default, sizeof(datm));
                    memcpy(dobs, dobs_default, sizeof(dobs));
                    memset(dret, 0, sizeof(dret));
                    int retval = swe_heliacal_angle(
                        tjds[it], dgeo, datm, dobs,
                        SEFLG_SWIEPH,
                        mags[0],
                        azi_objs[ia],
                        sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                        dret, serr);
                    if (retval == ERR) {
                        fprintf(stderr, "helangle ERR: %s\n", serr);
                        continue;
                    }
                    if (!first) printf(",\n");
                    first = 0;
                    printf("  {\"tjd_ut\": %.20e, \"mag\": %.20e, "
                           "\"azi_obj\": %.20e, "
                           "\"azi_sun\": %.20e, \"azi_moon\": %.20e, \"alt_moon\": %.20e, "
                           "\"helflag\": %d, \"retval\": %d, "
                           "\"dret\": [%.20e, %.20e, %.20e]}",
                           tjds[it], mags[0],
                           azi_objs[ia],
                           sun_moon[is][0], sun_moon[is][1], sun_moon[is][2],
                           SEFLG_SWIEPH, retval,
                           dret[0], dret[1], dret[2]);
                }
    }
    printf("\n],\n");

    /* ── pheno: swe_heliacal_pheno_ut ────────────────────────────────
     * Moon × TypeEvent=3 (evening first) at young-crescent evenings,
     * Venus × TypeEvent=1/2, Sirius × TypeEvent=1.
     * Observer: Mecca for Moon, Cairo for planets/stars.
     */
    printf("\"pheno\": [\n");
    first = 1;
    {
        double dgeo_mecca[3] = {39.83, 21.42, 300.0};
        double dgeo_cairo[3] = {31.25, 30.1, 30.0};

        struct pheno_case {
            const char *obj;
            int type_event;
            double tjd;
            double *geo;
            int helflag;
            const char *desc;
        };

        struct pheno_case cases[] = {
            /* Moon evening first at young-crescent evenings, Mecca */
            {"moon", 3, 2453469.229,  dgeo_mecca, SEFLG_SWIEPH, "moon_ef_mecca1"},
            {"moon", 3, 2453498.240,  dgeo_mecca, SEFLG_SWIEPH, "moon_ef_mecca2"},
            /* Moon with HIGH_PRECISION */
            {"moon", 3, 2453469.229,  dgeo_mecca, SEFLG_SWIEPH | SE_HELFLAG_AVKIND_VR, "moon_ef_mecca1_hp"},
            /* Venus morning first, Cairo */
            {"venus", 1, 2453720.0,   dgeo_cairo, SEFLG_SWIEPH, "venus_mf_cairo1"},
            {"venus", 2, 2453000.0,   dgeo_cairo, SEFLG_SWIEPH, "venus_el_cairo1"},
            {"venus", 1, 2453720.0,   dgeo_cairo, SEFLG_SWIEPH | SE_HELFLAG_AVKIND_VR, "venus_mf_cairo1_hp"},
            /* Sirius morning first, Cairo */
            {"sirius", 1, 2453586.0,  dgeo_cairo, SEFLG_SWIEPH, "sirius_mf_cairo"},
            {"sirius", 1, 2453586.0,  dgeo_cairo, SEFLG_SWIEPH | SE_HELFLAG_AVKIND_VR, "sirius_mf_cairo_hp"},
            /* Mercury morning first, Cairo */
            {"mercury", 1, 2453720.0, dgeo_cairo, SEFLG_SWIEPH, "mercury_mf_cairo"},
            /* Venus evening last, Mecca */
            {"venus", 2, 2453500.0,   dgeo_mecca, SEFLG_SWIEPH, "venus_el_mecca"},
            /* Moon morning last, Cairo */
            {"moon", 4, 2453469.0,    dgeo_cairo, SEFLG_SWIEPH, "moon_ml_cairo"},
            /* Mars evening first (early-exit guard path) */
            {"mars", 3, 2453500.0,    dgeo_cairo, SEFLG_SWIEPH, "mars_ef_cairo"},
            /* Jupiter evening first (early-exit guard path) */
            {"jupiter", 3, 2453500.0, dgeo_cairo, SEFLG_SWIEPH, "jupiter_ef_cairo"},
            /* MOSEPH duplicates */
            {"venus", 1, 2453720.0,   dgeo_cairo, SEFLG_MOSEPH, "venus_mf_cairo_moseph"},
            {"moon", 3, 2453469.229,  dgeo_mecca, SEFLG_MOSEPH, "moon_ef_mecca_moseph"},
        };
        int ncases = sizeof(cases) / sizeof(cases[0]);

        for (int i = 0; i < ncases; i++) {
            double datm[4], dobs[6], darr[50];
            char objname[256], serr[256];
            memcpy(datm, datm_default, sizeof(datm));
            memcpy(dobs, dobs_default, sizeof(dobs));
            memset(darr, 0, sizeof(darr));
            strncpy(objname, cases[i].obj, sizeof(objname)-1);
            objname[sizeof(objname)-1] = '\0';
            serr[0] = '\0';

            int retval = swe_heliacal_pheno_ut(
                cases[i].tjd,
                cases[i].geo,
                datm, dobs,
                objname,
                cases[i].type_event,
                cases[i].helflag,
                darr, serr);

            if (retval == ERR) {
                fprintf(stderr, "pheno ERR: %s obj=%s jd=%.4f desc=%s\n",
                        serr, cases[i].obj, cases[i].tjd, cases[i].desc);
                continue;
            }
            if (!first) printf(",\n");
            first = 0;
            printf("  {\"tjd_ut\": %.20e, \"object\": \"%s\", "
                   "\"type_event\": %d, "
                   "\"geo\": [%.20e, %.20e, %.20e], "
                   "\"helflag\": %d, \"desc\": \"%s\", "
                   "\"retval\": %d, \"darr\": [",
                   cases[i].tjd, cases[i].obj,
                   cases[i].type_event,
                   cases[i].geo[0], cases[i].geo[1], cases[i].geo[2],
                   cases[i].helflag, cases[i].desc,
                   retval);
            for (int j = 0; j < 28; j++) {
                if (j > 0) printf(", ");
                printf("%.20e", darr[j]);
            }
            printf("]}");
        }
    }
    printf("\n]}\n");

    return 0;
}
