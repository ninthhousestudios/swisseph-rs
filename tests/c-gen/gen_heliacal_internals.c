/*
 * Generates golden reference data for the heliacal module's pure-math
 * atmospheric extinction and optics layer (swehel.c static functions).
 *
 * Because the functions under test are `static` in swehel.c, we include
 * the source directly — this gives us access to all internals without
 * needing to modify the upstream C.
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -I../../swisseph -o gen_heliacal_internals \
 *      gen_heliacal_internals.c ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_heliacal_internals > ../golden-data/heliacal_internals.json
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

/* Include swehel.c to get access to static functions.
 * Use the full relative path since -I only covers headers. */
#include "../../../swisseph/swehel.c"

/* ── Helpers ──────────────────────────────────────────────────────── */

static int first_item = 1;

static void comma(void) {
    if (!first_item) printf(",");
    first_item = 0;
}

/* ── Extinction battery ──────────────────────────────────────────── */

static void print_extinction(void) {
    /* Rotated grid: ~96 cases */
    double alts[] = {-1, 0, 5, 20, 45, 90};
    double alt_s_vals[] = {-30, -10, 0, 10};
    double sunra_vals[] = {50, 200};
    double lat_vals[] = {30.1, 47.4};
    double height_vals[] = {0, 1500};
    /* Two datm configs */
    double datm_configs[][4] = {
        {1013.25, 15, 40, 40},
        {900, 25, 70, 20}
    };

    int n_alto = sizeof(alts)/sizeof(alts[0]);
    int n_alts = sizeof(alt_s_vals)/sizeof(alt_s_vals[0]);
    int n_sunra = sizeof(sunra_vals)/sizeof(sunra_vals[0]);
    int n_lat = sizeof(lat_vals)/sizeof(lat_vals[0]);
    int n_height = sizeof(height_vals)/sizeof(height_vals[0]);
    int n_datm = sizeof(datm_configs)/sizeof(datm_configs[0]);

    char serr[AS_MAXCH];

    printf("\"extinction\":[");
    first_item = 1;

    for (int ia = 0; ia < n_alto; ia++) {
        double AltO = alts[ia];
        for (int is = 0; is < n_alts; is++) {
            double AltS = alt_s_vals[is];
            for (int ir = 0; ir < n_sunra; ir++) {
                double sunra = sunra_vals[ir];
                /* Rotate through lat/height/datm to reduce the cross-product */
                int idx = ia * n_alts * n_sunra + is * n_sunra + ir;
                int il = idx % n_lat;
                int ih = idx % n_height;
                int id = idx % n_datm;

                double Lat = lat_vals[il];
                double HeightEye = height_vals[ih];
                double datm[4];
                memcpy(datm, datm_configs[id], sizeof(datm));

                /* Compute extinction components */
                int helflag = 0;

                double dm = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
                double kt_val = kt(AltS, sunra, Lat, HeightEye, datm[1], datm[2], datm[3], 4, serr);
                double kR_val = kR(AltS, HeightEye);
                double kOZ_val = kOZ(AltS, sunra, Lat);
                double kW_val = kW(HeightEye, datm[1], datm[2]);
                double ka_val = ka(AltS, sunra, Lat, HeightEye, datm[1], datm[2], datm[3], serr);

                comma();
                printf("{\"AltO\":%.17g,\"AltS\":%.17g,\"sunra\":%.17g,"
                       "\"Lat\":%.17g,\"HeightEye\":%.17g,"
                       "\"datm\":[%.17g,%.17g,%.17g,%.17g],"
                       "\"Deltam\":%.17g,\"kt\":%.17g,\"kR\":%.17g,"
                       "\"kOZ\":%.17g,\"kW\":%.17g,\"ka\":%.17g}\n",
                       AltO, AltS, sunra, Lat, HeightEye,
                       datm[0], datm[1], datm[2], datm[3],
                       dm, kt_val, kR_val, kOZ_val, kW_val, ka_val);
            }
        }
    }
    printf("]");
}

/* ── Airmass battery ─────────────────────────────────────────────── */

static void print_airmass(void) {
    double app_alt_vals[] = {0, 1, 5, 10, 30, 60, 90};
    double press_vals[] = {1013.25, 900};
    int n_alt = sizeof(app_alt_vals)/sizeof(app_alt_vals[0]);
    int n_press = sizeof(press_vals)/sizeof(press_vals[0]);

    printf("\"airmass\":[");
    first_item = 1;

    for (int ia = 0; ia < n_alt; ia++) {
        double AppAltO = app_alt_vals[ia];
        for (int ip = 0; ip < n_press; ip++) {
            double Press = press_vals[ip];
            double airm = Airmass(AppAltO, Press);

            /* Also compute Xext/Xlay at the same zenith distance */
            double zend = (90 - AppAltO) * DEGTORAD;
            if (zend > PI/2) zend = PI/2;
            double xr = Xext(scaleHrayleigh, zend, Press);
            double xw = Xext(scaleHwater, zend, Press);
            double xa = Xext(scaleHaerosol, zend, Press);
            double xoz = Xlay(scaleHozone, zend, Press);

            comma();
            printf("{\"AppAltO\":%.17g,\"Press\":%.17g,"
                   "\"Airmass\":%.17g,\"Xext_rayleigh\":%.17g,"
                   "\"Xext_water\":%.17g,\"Xext_aerosol\":%.17g,"
                   "\"Xlay_ozone\":%.17g}\n",
                   AppAltO, Press, airm, xr, xw, xa, xoz);
        }
    }
    printf("]");
}

/* ── AppAlt battery ──────────────────────────────────────────────── */

static void print_app_alt(void) {
    double alt_vals[] = {-2, 0, 1, 5, 20, 60};
    double temp_vals[] = {0, 15, 30};
    double pres_vals[] = {1013.25, 900};
    int n_alt = sizeof(alt_vals)/sizeof(alt_vals[0]);
    int n_temp = sizeof(temp_vals)/sizeof(temp_vals[0]);
    int n_pres = sizeof(pres_vals)/sizeof(pres_vals[0]);

    printf("\"app_alt\":[");
    first_item = 1;

    for (int ia = 0; ia < n_alt; ia++) {
        double alt = alt_vals[ia];
        for (int it = 0; it < n_temp; it++) {
            double TempE = temp_vals[it];
            for (int ip = 0; ip < n_pres; ip++) {
                double PresE = pres_vals[ip];
                double app = AppAltfromTopoAlt(alt, TempE, PresE, 0);
                double topo = TopoAltfromAppAlt(alt, TempE, PresE);

                comma();
                printf("{\"alt\":%.17g,\"TempE\":%.17g,\"PresE\":%.17g,"
                       "\"AppAltfromTopoAlt\":%.17g,"
                       "\"TopoAltfromAppAlt\":%.17g}\n",
                       alt, TempE, PresE, app, topo);
            }
        }
    }
    printf("]");
}

/* ── Optics battery ──────────────────────────────────────────────── */

static void print_optic(void) {
    double B_vals[] = {1e-5, 1.0, 100.0, 1645.0, 1e5};
    int n_b = sizeof(B_vals)/sizeof(B_vals[0]);

    printf("\"optic\":[");
    first_item = 1;

    struct {
        double dobs[6];
        int helflag;
        const char *label;
    } configs[] = {
        {{36, 1, 0, 0, 0, 0}, 0, "default"},
        {{60, 1, 0, 0, 0, 0}, 0, "age60"},
        {{36, 1, 1, 0, 0, 0}, 0, "binocular"},
        {{36, 1, 1, 10, 50, 0.8}, SE_HELFLAG_OPTICAL_PARAMS, "optical_params"},
    };
    int n_configs = sizeof(configs)/sizeof(configs[0]);

    for (int ib = 0; ib < n_b; ib++) {
        double B = B_vals[ib];
        for (int ic = 0; ic < n_configs; ic++) {
            /* Make a copy and apply defaults */
            double dobs[6];
            memcpy(dobs, configs[ic].dobs, sizeof(dobs));
            double datm[4] = {1013.25, 15, 40, 0};
            double dgeo[3] = {8.55, 47.37, 500};
            int helflag = configs[ic].helflag;
            default_heliacal_parameters(datm, dgeo, dobs, helflag);

            double kX = 0.3; /* representative extinction */
            double cvs = CVA(B, dobs[1], helflag);
            double pd = PupilDia(dobs[0], B);
            double of_intensity = OpticFactor(B, kX, dobs, 2451545.0, "", 0, helflag);
            double of_background = OpticFactor(B, kX, dobs, 2451545.0, "", 1, helflag);

            comma();
            printf("{\"B\":%.17g,\"config\":\"%s\","
                   "\"dobs\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],"
                   "\"helflag\":%d,\"kX\":%.17g,"
                   "\"CVA\":%.17g,\"PupilDia\":%.17g,"
                   "\"OpticFactor_intensity\":%.17g,"
                   "\"OpticFactor_background\":%.17g}\n",
                   B, configs[ic].label,
                   dobs[0], dobs[1], dobs[2], dobs[3], dobs[4], dobs[5],
                   helflag, kX,
                   cvs, pd, of_intensity, of_background);
        }
    }
    printf("]");
}

/* ── Brightness battery ──────────────────────────────────────────── */

static void print_brightness(void) {
    /* Rotated grid: ~50 cases covering the parameter space */
    double alt_o_vals[] = {0, 5, 20, 60};
    double azi_o = 180;
    double alt_s_vals[] = {-30, -12, -5, -1, 2, 10, 30};
    double azi_s_vals[] = {90, 270};
    double alt_m_vals[] = {-10, 0, 30};
    double azi_m = 120;
    double sunra_vals[] = {50, 200};
    double lat = 30.1;
    double height_vals[] = {0, 1500};
    double datm_configs[][4] = {
        {1013.25, 15, 40, 40},
        {900, 25, 70, 20}
    };
    double jdn_vals[] = {2451545.0, 2455197.5};

    int n_alto = sizeof(alt_o_vals)/sizeof(alt_o_vals[0]);
    int n_alts = sizeof(alt_s_vals)/sizeof(alt_s_vals[0]);
    int n_azis = sizeof(azi_s_vals)/sizeof(azi_s_vals[0]);
    int n_altm = sizeof(alt_m_vals)/sizeof(alt_m_vals[0]);
    int n_sunra = sizeof(sunra_vals)/sizeof(sunra_vals[0]);
    int n_height = sizeof(height_vals)/sizeof(height_vals[0]);
    int n_datm = sizeof(datm_configs)/sizeof(datm_configs[0]);
    int n_jdn = sizeof(jdn_vals)/sizeof(jdn_vals[0]);

    char serr[AS_MAXCH];
    int helflag = 0;

    printf("\"brightness\":[");
    first_item = 1;

    int idx = 0;
    for (int ia = 0; ia < n_alto; ia++) {
        double AltO = alt_o_vals[ia];
        for (int is = 0; is < n_alts; is++) {
            double AltS = alt_s_vals[is];
            /* Rotate through other dimensions */
            int izs = idx % n_azis;
            int im = idx % n_altm;
            int ir = idx % n_sunra;
            int ih = idx % n_height;
            int id = idx % n_datm;
            int ij = idx % n_jdn;

            double AziS = azi_s_vals[izs];
            double AltM = alt_m_vals[im];
            double sunra = sunra_vals[ir];
            double HeightEye = height_vals[ih];
            double datm[4];
            memcpy(datm, datm_configs[id], sizeof(datm));
            double JDNDaysUT = jdn_vals[ij];

            double bn_val = Bn(AltO, JDNDaysUT, AltS, sunra, lat, HeightEye, datm, helflag, serr);
            double bm_val = Bm(AltO, azi_o, AltM, azi_m, AltS, AziS, sunra, lat, HeightEye, datm, helflag, serr);
            double btwi_val = Btwi(AltO, azi_o, AltS, AziS, sunra, lat, HeightEye, datm, helflag, serr);
            double bday_val = Bday(AltO, azi_o, AltS, AziS, sunra, lat, HeightEye, datm, helflag, serr);
            double bcity_val = Bcity(0, datm[0]);
            double bsky_val = Bsky(AltO, azi_o, AltM, azi_m, JDNDaysUT, AltS, AziS, sunra, lat, HeightEye, datm, helflag, serr);

            comma();
            printf("{\"AltO\":%.17g,\"AziO\":%.17g,"
                   "\"AltM\":%.17g,\"AziM\":%.17g,"
                   "\"AltS\":%.17g,\"AziS\":%.17g,"
                   "\"sunra\":%.17g,\"Lat\":%.17g,"
                   "\"HeightEye\":%.17g,"
                   "\"datm\":[%.17g,%.17g,%.17g,%.17g],"
                   "\"JDNDaysUT\":%.17g,"
                   "\"Bn\":%.17g,\"Bm\":%.17g,\"Btwi\":%.17g,"
                   "\"Bday\":%.17g,\"Bcity\":%.17g,\"Bsky\":%.17g}\n",
                   AltO, azi_o,
                   AltM, azi_m,
                   AltS, AziS,
                   sunra, lat,
                   HeightEye,
                   datm[0], datm[1], datm[2], datm[3],
                   JDNDaysUT,
                   bn_val, bm_val, btwi_val,
                   bday_val, bcity_val, bsky_val);
            idx++;
        }
    }
    printf("]");
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(void) {
    printf("{");
    print_extinction();
    printf(",");
    print_airmass();
    printf(",");
    print_app_alt();
    printf(",");
    print_optic();
    printf(",");
    print_brightness();
    printf("}\n");
    return 0;
}
