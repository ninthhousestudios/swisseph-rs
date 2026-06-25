pub mod data;

use crate::constants::*;
use crate::flags::CalcFlags;
use crate::math::{normalize_degrees, normalize_radians};
use crate::types::*;
use data::*;

// Max multiplier magnitudes per Delaunay arg: MM, MS, FF, DD, OM
const MAX_MULT: [usize; 5] = [3, 2, 4, 4, 2];

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn nutation(jd: f64, flags: CalcFlags, models: &AstroModels) -> Nutation {
    let is_jplhor = flags.contains(CalcFlags::DPSIDEPS_1980)
        || (flags.contains(CalcFlags::JPLHOR_APPROX)
            && models.jplhora_mode == JplHoraMode::V3
            && jd <= DPSI_DEPS_IAU1980_TJD0_HORIZONS);

    if is_jplhor {
        let mut nut = calc_nutation_iau1980(jd, NutationModel::IAU1980);
        if !flags.contains(CalcFlags::DPSIDEPS_1980) {
            nut.dpsi += DPSI_IAU1980_TJD0 / 3600.0 * DEGTORAD;
            nut.deps += DEPS_IAU1980_TJD0 / 3600.0 * DEGTORAD;
        }
        return nut;
    }

    let model = models.nutation;
    let mut nut = match model {
        NutationModel::IAU1980 | NutationModel::IAUCorr1987 => calc_nutation_iau1980(jd, model),
        NutationModel::IAU2000A | NutationModel::IAU2000B => calc_nutation_iau2000ab(jd, model),
        NutationModel::Woolard => calc_nutation_woolard(jd),
    };

    if matches!(model, NutationModel::IAU2000A | NutationModel::IAU2000B)
        && flags.contains(CalcFlags::JPLHOR_APPROX)
        && models.jplhora_mode == JplHoraMode::V2
    {
        nut.dpsi += -41.7750 / 3600.0 / 1000.0 * DEGTORAD;
        nut.deps += -6.8192 / 3600.0 / 1000.0 * DEGTORAD;
    }

    nut
}

// ---------------------------------------------------------------------------
// IAU 1980 — Delaunay arguments (Seidelmann 1982 / FK5)
// ---------------------------------------------------------------------------

fn delaunay_iau1980(t: f64) -> [f64; 5] {
    let t2 = t * t;
    // C evaluation form: c1*T + c0 + (c3*T + c2)*T2
    let mm = 1717915922.633 * t + 485866.733 + (0.064 * t + 31.310) * t2;
    let ms = 129596581.224 * t + 1287099.804 + (-0.012 * t - 0.577) * t2;
    let ff = 1739527263.137 * t + 335778.877 + (0.011 * t - 13.257) * t2;
    let dd = 1602961601.328 * t + 1072261.307 + (0.019 * t - 6.891) * t2;
    let om = -6962890.539 * t + 450160.280 + (0.008 * t + 7.455) * t2;
    [mm, ms, ff, dd, om].map(|a| normalize_degrees(a / 3600.0) * DEGTORAD)
}

// ---------------------------------------------------------------------------
// Sin/cos multiple precomputation (recurrence)
// ---------------------------------------------------------------------------

fn precompute_sincos(args: &[f64; 5]) -> ([[f64; 4]; 5], [[f64; 4]; 5]) {
    let mut ss = [[0.0_f64; 4]; 5];
    let mut cc = [[0.0_f64; 4]; 5];
    for k in 0..5 {
        let su = args[k].sin();
        let cu = args[k].cos();
        ss[k][0] = su;
        cc[k][0] = cu;
        let mut sv = 2.0 * su * cu;
        let mut cv = cu * cu - su * su;
        ss[k][1] = sv;
        cc[k][1] = cv;
        for i in 2..MAX_MULT[k] {
            let s = su * cv + cu * sv;
            cv = cu * cv - su * sv;
            sv = s;
            ss[k][i] = sv;
            cc[k][i] = cv;
        }
    }
    (ss, cc)
}

// ---------------------------------------------------------------------------
// Sin/cos of a linear combination of arguments
// ---------------------------------------------------------------------------

fn sincos_of_combination(ss: &[[f64; 4]; 5], cc: &[[f64; 4]; 5], mults: &[i16; 5]) -> (f64, f64) {
    let mut sin_acc = 0.0;
    let mut cos_acc = 1.0;
    let mut started = false;
    for m_idx in 0..5 {
        let mut j = mults[m_idx] as i32;
        if j > 100 {
            j = 0;
        }
        if j == 0 {
            continue;
        }
        let abs_j = j.unsigned_abs() as usize;
        let mut su = ss[m_idx][abs_j - 1];
        if j < 0 {
            su = -su;
        }
        let cu = cc[m_idx][abs_j - 1];
        if !started {
            sin_acc = su;
            cos_acc = cu;
            started = true;
        } else {
            let sw = su * cos_acc + cu * sin_acc;
            cos_acc = cu * cos_acc - su * sin_acc;
            sin_acc = sw;
        }
    }
    (sin_acc, cos_acc)
}

// ---------------------------------------------------------------------------
// IAU 1980 nutation (105 standard + 7 Herring correction terms)
// ---------------------------------------------------------------------------

fn calc_nutation_iau1980(jd: f64, model: NutationModel) -> Nutation {
    let t = (jd - J2000) / 36525.0;
    let args = delaunay_iau1980(t);
    let (ss, cc) = precompute_sincos(&args);

    // Dominant OM term (not in NT table)
    let mut dpsi = (-0.01742 * t - 17.1996) * ss[4][0];
    let mut deps = (0.00089 * t + 9.2025) * cc[4][0];

    for row in NT.iter() {
        if row[0] >= 100 && model != NutationModel::IAUCorr1987 {
            continue;
        }
        let mults: [i16; 5] = [row[0], row[1], row[2], row[3], row[4]];
        let (sv, cv) = sincos_of_combination(&ss, &cc, &mults);

        let mut f = row[5] as f64 * 0.0001;
        if row[6] != 0 {
            f += 0.00001 * t * row[6] as f64;
        }
        let mut g = row[7] as f64 * 0.0001;
        if row[8] != 0 {
            g += 0.00001 * t * row[8] as f64;
        }
        if row[0] >= 100 {
            f *= 0.1;
            g *= 0.1;
        }
        if row[0] != 102 {
            dpsi += f * sv;
            deps += g * cv;
        } else {
            dpsi += f * cv;
            deps += g * sv;
        }
    }

    Nutation {
        dpsi: dpsi * DEGTORAD / 3600.0,
        deps: deps * DEGTORAD / 3600.0,
    }
}

// ---------------------------------------------------------------------------
// IAU 2000 — Delaunay arguments (Simon et al. 1994)
// ---------------------------------------------------------------------------

fn delaunay_iau2000(t: f64) -> [f64; 5] {
    let m =
        485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470))));
    let sm =
        1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149))));
    let f =
        335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * 0.00000417)));
    let d =
        1072260.70369 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169))));
    let om =
        450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939))));
    [m, sm, f, d, om].map(|a| normalize_degrees(a / 3600.0) * DEGTORAD)
}

// ---------------------------------------------------------------------------
// IAU 2000 planetary arguments (linear in radians)
// ---------------------------------------------------------------------------

fn planetary_args_iau2000(t: f64) -> [f64; 14] {
    [
        normalize_radians(2.35555598 + 8328.6914269554 * t),
        normalize_radians(6.24006013 + 628.301955 * t),
        normalize_radians(1.627905234 + 8433.466158131 * t),
        normalize_radians(5.198466741 + 7771.3771468121 * t),
        normalize_radians(2.18243920 + -33.757045 * t),
        normalize_radians(4.402608842 + 2608.7903141574 * t),
        normalize_radians(3.176146697 + 1021.3285546211 * t),
        normalize_radians(1.753470314 + 628.3075849991 * t),
        normalize_radians(6.203480913 + 334.0612426700 * t),
        normalize_radians(0.599546497 + 52.9690962641 * t),
        normalize_radians(0.874016757 + 21.3299104960 * t),
        normalize_radians(5.481293871 + 7.4781598567 * t),
        normalize_radians(5.321159000 + 3.8127774000 * t),
        (0.02438175 + 0.00000538691 * t) * t,
    ]
}

// ---------------------------------------------------------------------------
// IAU 2000A/B nutation
// ---------------------------------------------------------------------------

fn calc_nutation_iau2000ab(jd: f64, model: NutationModel) -> Nutation {
    let t = (jd - J2000) / 36525.0;
    let del = delaunay_iau2000(t);
    let [m, sm, f, d, om] = del;

    // Luni-solar nutation
    let inls = if model == NutationModel::IAU2000B {
        NLS_2000B
    } else {
        NLS_COUNT
    };
    let mut dpsi = 0.0_f64;
    let mut deps = 0.0_f64;
    for i in (0..inls).rev() {
        let n = &NLS[i];
        let c = &CLS[i];
        let darg = normalize_radians(
            n[0] as f64 * m
                + n[1] as f64 * sm
                + n[2] as f64 * f
                + n[3] as f64 * d
                + n[4] as f64 * om,
        );
        let (sinarg, cosarg) = (darg.sin(), darg.cos());
        dpsi += (c[0] as f64 + c[1] as f64 * t) * sinarg + c[2] as f64 * cosarg;
        deps += (c[3] as f64 + c[4] as f64 * t) * cosarg + c[5] as f64 * sinarg;
    }
    let mut nutlo_0 = dpsi * O1MAS2DEG;
    let mut nutlo_1 = deps * O1MAS2DEG;

    // Planetary nutation + P03 correction (2000A only)
    if model == NutationModel::IAU2000A {
        let pa = planetary_args_iau2000(t);
        dpsi = 0.0;
        deps = 0.0;
        for i in (0..NPL_COUNT).rev() {
            let n = &NPL[i];
            let ic = &ICPL[i];
            let mut darg = 0.0;
            for j in 0..14 {
                darg += n[j] as f64 * pa[j];
            }
            let darg = normalize_radians(darg);
            let (sinarg, cosarg) = (darg.sin(), darg.cos());
            dpsi += ic[0] as f64 * sinarg + ic[1] as f64 * cosarg;
            deps += ic[2] as f64 * sinarg + ic[3] as f64 * cosarg;
        }
        nutlo_0 += dpsi * O1MAS2DEG;
        nutlo_1 += deps * O1MAS2DEG;

        // P03 precession correction (Capitaine et al. 2005), in microarcseconds
        let dpsi_p03 = -8.1 * om.sin() - 0.6 * (2.0 * f - 2.0 * d + 2.0 * om).sin()
            + t * (47.8 * om.sin()
                + 3.7 * (2.0 * f - 2.0 * d + 2.0 * om).sin()
                + 0.6 * (2.0 * f + 2.0 * om).sin()
                - 0.6 * (2.0 * om).sin());
        let deps_p03 = t * (-25.6 * om.cos() - 1.6 * (2.0 * f - 2.0 * d + 2.0 * om).cos());
        nutlo_0 += dpsi_p03 / (3600.0 * 1_000_000.0);
        nutlo_1 += deps_p03 / (3600.0 * 1_000_000.0);
    }

    Nutation {
        dpsi: nutlo_0 * DEGTORAD,
        deps: nutlo_1 * DEGTORAD,
    }
}

// ---------------------------------------------------------------------------
// Woolard nutation
// ---------------------------------------------------------------------------

fn calc_nutation_woolard(jd: f64) -> Nutation {
    let t = (jd - J1900) / 36525.0;
    let t2 = t * t;
    let frac = |x: f64| -> f64 { x - (x as i64 as f64) };

    let ls = 279.697 + 0.000303 * t2 + 360.0 * frac(100.0021358 * t);
    let ld = 270.434 - 0.001133 * t2 + 360.0 * frac(1336.855231 * t);
    let ms = 358.476 - 0.00015 * t2 + 360.0 * frac(99.99736056000026 * t);
    let md = 296.105 + 0.009192 * t2 + 360.0 * frac(13255523.59 * t);
    let nm = 259.183 + 0.002078 * t2 - 360.0 * frac(5.372616667 * t);

    let tls = 2.0 * ls * DEGTORAD;
    let nm = nm * DEGTORAD;
    let tnm = 2.0 * nm;
    let ms = ms * DEGTORAD;
    let tld = 2.0 * ld * DEGTORAD;
    let md = md * DEGTORAD;

    let dpsi = (-17.2327 - 0.01737 * t) * nm.sin()
        + (-1.2729 - 0.00013 * t) * tls.sin()
        + 0.2088 * tnm.sin()
        - 0.2037 * tld.sin()
        + (0.1261 - 0.00031 * t) * ms.sin()
        + 0.0675 * md.sin()
        - (0.0497 - 0.00012 * t) * (tls + ms).sin()
        - 0.0342 * (tld - nm).sin()
        - 0.0261 * (tld + md).sin()
        + 0.0214 * (tls - ms).sin()
        - 0.0149 * (tls - tld + md).sin()
        + 0.0124 * (tls - nm).sin()
        + 0.0114 * (tld - md).sin();

    let deps = (9.21 + 0.00091 * t) * nm.cos() + (0.5522 - 0.00029 * t) * tls.cos()
        - 0.0904 * tnm.cos()
        + 0.0884 * tld.cos()
        + 0.0216 * (tls + ms).cos()
        + 0.0183 * (tld - nm).cos()
        + 0.0113 * (tld + md).cos()
        - 0.0093 * (tls - ms).cos()
        - 0.0066 * (tls - nm).cos();

    Nutation {
        dpsi: dpsi / 3600.0 * DEGTORAD,
        deps: deps / 3600.0 * DEGTORAD,
    }
}
