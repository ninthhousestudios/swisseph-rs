/*
 * Generates golden reference data for the corrections module:
 *   1. meff() — copied from sweph.c since it's static
 *   2. swi_aberr_light() — called directly with iflag=0 (position only)
 *   3. Pipeline tests — swe_calc with various correction flag combos
 *
 * Compile:
 *   cc -Wall -I../../swisseph -o gen_corrections gen_corrections.c \
 *      ../../swisseph/libswe.a -lm
 * Run:
 *   ./gen_corrections > ../golden-data/corrections.json
 */

#include <stdio.h>
#include <math.h>
#include "swephexp.h"

/* ---------- meff: copied verbatim from sweph.c (static there) ---------- */

struct meff_ele { double r, m; };
static const struct meff_ele eff_arr[] = {
    {1.000, 1.000000}, {0.990, 0.999979}, {0.980, 0.999940}, {0.970, 0.999881},
    {0.960, 0.999811}, {0.950, 0.999724}, {0.940, 0.999622}, {0.930, 0.999497},
    {0.920, 0.999354}, {0.910, 0.999192}, {0.900, 0.999000}, {0.890, 0.998786},
    {0.880, 0.998535}, {0.870, 0.998242}, {0.860, 0.997919}, {0.850, 0.997571},
    {0.840, 0.997198}, {0.830, 0.996792}, {0.820, 0.996316}, {0.810, 0.995791},
    {0.800, 0.995226}, {0.790, 0.994625}, {0.780, 0.993991}, {0.770, 0.993326},
    {0.760, 0.992598}, {0.750, 0.991770}, {0.740, 0.990873}, {0.730, 0.989919},
    {0.720, 0.988912}, {0.710, 0.987856}, {0.700, 0.986755}, {0.690, 0.985610},
    {0.680, 0.984398}, {0.670, 0.982986}, {0.660, 0.981437}, {0.650, 0.979779},
    {0.640, 0.978024}, {0.630, 0.976182}, {0.620, 0.974256}, {0.610, 0.972253},
    {0.600, 0.970174}, {0.590, 0.968024}, {0.580, 0.965594}, {0.570, 0.962797},
    {0.560, 0.959758}, {0.550, 0.956515}, {0.540, 0.953088}, {0.530, 0.949495},
    {0.520, 0.945741}, {0.510, 0.941838}, {0.500, 0.937790}, {0.490, 0.933563},
    {0.480, 0.928668}, {0.470, 0.923288}, {0.460, 0.917527}, {0.450, 0.911432},
    {0.440, 0.905035}, {0.430, 0.898353}, {0.420, 0.891022}, {0.410, 0.882940},
    {0.400, 0.874312}, {0.390, 0.865206}, {0.380, 0.855423}, {0.370, 0.844619},
    {0.360, 0.833074}, {0.350, 0.820876}, {0.340, 0.808031}, {0.330, 0.793962},
    {0.320, 0.778931}, {0.310, 0.763021}, {0.300, 0.745815}, {0.290, 0.727557},
    {0.280, 0.708234}, {0.270, 0.687583}, {0.260, 0.665741}, {0.250, 0.642597},
    {0.240, 0.618252}, {0.230, 0.592586}, {0.220, 0.565747}, {0.210, 0.537697},
    {0.200, 0.508554}, {0.190, 0.478420}, {0.180, 0.447322}, {0.170, 0.415454},
    {0.160, 0.382892}, {0.150, 0.349955}, {0.140, 0.316691}, {0.130, 0.283565},
    {0.120, 0.250431}, {0.110, 0.218327}, {0.100, 0.186794}, {0.090, 0.156287},
    {0.080, 0.128421}, {0.070, 0.102237}, {0.060, 0.077393}, {0.050, 0.054833},
    {0.040, 0.036361}, {0.030, 0.020953}, {0.020, 0.009645}, {0.010, 0.002767},
    {0.000, 0.000000}
};

static double meff(double r) {
    double f, m;
    int i;
    if (r <= 0) return 0.0;
    if (r >= 1) return 1.0;
    for (i = 0; eff_arr[i].r > r; i++)
        ;
    f = (r - eff_arr[i-1].r) / (eff_arr[i].r - eff_arr[i-1].r);
    m = eff_arr[i-1].m + f * (eff_arr[i].m - eff_arr[i-1].m);
    return m;
}

/* ---------- extern declaration for swi_aberr_light ---------- */

extern void swi_aberr_light(double *xx, double *xe, int32 iflag);

/* ---------- main ---------- */

int main(void) {
    char serr[256];
    swe_set_ephe_path(NULL);
    printf("{\n");

    /* ======== Section 1: meff ======== */
    printf("  \"meff\": [\n");
    double meff_inputs[] = {
        -0.1, 0.0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05,
        0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
        0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
        0.90, 0.95, 0.99, 0.995, 1.0, 1.5
    };
    int nmeff = sizeof(meff_inputs) / sizeof(meff_inputs[0]);
    for (int i = 0; i < nmeff; i++) {
        double r = meff_inputs[i];
        double m = meff(r);
        printf("    {\"r\": %.20e, \"result\": %.20e}", r, m);
        if (i < nmeff - 1) printf(",");
        printf("\n");
    }
    printf("  ],\n");

    /* ======== Section 2: aberr_light (position only, iflag=0) ======== */
    printf("  \"aberr_light\": [\n");

    /* Earth velocity in AU/day — realistic orbital velocity ~29.78 km/s */
    double earth_states[][6] = {
        /* pos (unused by our fn) + vel (AU/day) */
        {0.0, 0.0, 0.0, -0.0172, 0.0, 0.0},           /* vx only */
        {0.0, 0.0, 0.0, 0.0, 0.0172, 0.0},             /* vy only */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0075},             /* vz only */
        {0.0, 0.0, 0.0, -0.01, 0.0149, 0.0003},        /* realistic mix */
        {0.0, 0.0, 0.0, -0.0001, 0.0001, 0.0},         /* near-zero vel */
    };
    int nearth = sizeof(earth_states) / sizeof(earth_states[0]);

    double planet_pos[][6] = {
        /* positions (AU) + speed (unused for iflag=0) */
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},               /* unit x */
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},               /* unit y */
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},               /* unit z */
        {1.0, 2.0, 3.0, 0.0, 0.0, 0.0},               /* diagonal */
        {0.5, 0.0, 0.0, 0.0, 0.0, 0.0},               /* 0.5 AU */
        {5.0, 1.0, -2.0, 0.0, 0.0, 0.0},              /* far + off-axis */
        {30.0, 5.0, 3.0, 0.0, 0.0, 0.0},              /* Neptune distance */
        {0.3, -0.4, 0.1, 0.0, 0.0, 0.0},              /* inner planet */
    };
    int nplanet = sizeof(planet_pos) / sizeof(planet_pos[0]);

    int aberr_count = 0;
    int aberr_total = nearth * nplanet;
    for (int ie = 0; ie < nearth; ie++) {
        for (int ip = 0; ip < nplanet; ip++) {
            double xx[6];
            for (int k = 0; k < 6; k++) xx[k] = planet_pos[ip][k];
            swi_aberr_light(xx, earth_states[ie], 0);
            printf("    {\"input\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
                   "\"earth\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e], "
                   "\"output\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                   planet_pos[ip][0], planet_pos[ip][1], planet_pos[ip][2],
                   planet_pos[ip][3], planet_pos[ip][4], planet_pos[ip][5],
                   earth_states[ie][0], earth_states[ie][1], earth_states[ie][2],
                   earth_states[ie][3], earth_states[ie][4], earth_states[ie][5],
                   xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
            aberr_count++;
            if (aberr_count < aberr_total) printf(",");
            printf("\n");
        }
    }
    printf("  ],\n");

    /* ======== Section 3: Pipeline tests via swe_calc ======== */
    printf("  \"pipeline\": [\n");

    int pipe_bodies[] = { SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN };
    const char *pipe_names[] = { "mercury", "venus", "mars", "jupiter", "saturn" };
    int npbodies = 5;
    double pipe_epochs[] = { 2451545.0, 2460545.0, 2433295.0 };
    int npepochs = 3;

    int base = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_ICRS | SEFLG_J2000
             | SEFLG_EQUATORIAL | SEFLG_XYZ;

    int pipe_count = 0;
    int pipe_total = npbodies * npepochs;
    for (int ib = 0; ib < npbodies; ib++) {
        for (int ie = 0; ie < npepochs; ie++) {
            double tjd = pipe_epochs[ie];
            double true_pos[6], aberr_pos[6], defl_pos[6], both_pos[6];

            swe_calc(tjd, pipe_bodies[ib],
                     base | SEFLG_NOABERR | SEFLG_NOGDEFL, true_pos, serr);
            swe_calc(tjd, pipe_bodies[ib],
                     base | SEFLG_NOGDEFL, aberr_pos, serr);
            swe_calc(tjd, pipe_bodies[ib],
                     base | SEFLG_NOABERR, defl_pos, serr);
            swe_calc(tjd, pipe_bodies[ib],
                     base, both_pos, serr);

            printf("    {\"tjd\": %.20e, \"body\": %d, \"body_name\": \"%s\",\n",
                   tjd, pipe_bodies[ib], pipe_names[ib]);
            printf("     \"true_pos\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e],\n",
                   true_pos[0], true_pos[1], true_pos[2],
                   true_pos[3], true_pos[4], true_pos[5]);
            printf("     \"aberr_pos\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e],\n",
                   aberr_pos[0], aberr_pos[1], aberr_pos[2],
                   aberr_pos[3], aberr_pos[4], aberr_pos[5]);
            printf("     \"defl_pos\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e],\n",
                   defl_pos[0], defl_pos[1], defl_pos[2],
                   defl_pos[3], defl_pos[4], defl_pos[5]);
            printf("     \"both_pos\": [%.20e, %.20e, %.20e, %.20e, %.20e, %.20e]}",
                   both_pos[0], both_pos[1], both_pos[2],
                   both_pos[3], both_pos[4], both_pos[5]);

            pipe_count++;
            if (pipe_count < pipe_total) printf(",");
            printf("\n");
        }
    }
    printf("  ]\n");

    printf("}\n");
    swe_close();
    return 0;
}
