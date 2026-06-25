use std::f64::consts::TAU;

use crate::constants::*;
use crate::flags::CalcFlags;
use crate::math::poly_eval;
use crate::types::*;

const OFFSET_EPS_JPLHORIZONS: f64 = 35.95;
const DCOR_EPS_JPL_TJD0: f64 = 2437846.5;

#[rustfmt::skip]
const DCOR_EPS_JPL: [f64; 51] = [
    36.726, 36.627, 36.595, 36.578, 36.640, 36.659, 36.731, 36.765,
    36.662, 36.555, 36.335, 36.321, 36.354, 36.227, 36.289, 36.348, 36.257, 36.163,
    35.979, 35.896, 35.842, 35.825, 35.912, 35.950, 36.093, 36.191, 36.009, 35.943,
    35.875, 35.771, 35.788, 35.753, 35.822, 35.866, 35.771, 35.732, 35.543, 35.498,
    35.449, 35.409, 35.497, 35.556, 35.672, 35.760, 35.596, 35.565, 35.510, 35.394,
    35.385, 35.375, 35.415,
];

pub fn obliquity(jd: f64, flags: CalcFlags, models: &AstroModels) -> Epsilon {
    let t = (jd - J2000) / 36525.0;

    let is_jplhor = flags.contains(CalcFlags::DPSIDEPS_1980)
        || (flags.contains(CalcFlags::JPLHOR_APPROX)
            && models.jplhora_mode == JplHoraMode::V3
            && jd <= DPSI_DEPS_IAU1980_TJD0_HORIZONS);

    let prec_short = models.prec_shortterm;
    let prec_long = models.prec_longterm;

    let eps = if is_jplhor {
        if jd > 2378131.5 && jd < 2525323.5 {
            obliquity_iau1976(t)
        } else {
            obliquity_owen1990(jd) * DEGTORAD
        }
    } else if flags.contains(CalcFlags::JPLHOR_APPROX) && models.jplhora_mode == JplHoraMode::V2 {
        obliquity_iau1976(t)
    } else if prec_short == PrecessionModel::IAU1976 && t.abs() <= PREC_IAU_1976_CTIES {
        obliquity_iau1976(t)
    } else if prec_long == PrecessionModel::IAU1976 {
        obliquity_iau1976(t)
    } else if prec_short == PrecessionModel::IAU2000 && t.abs() <= PREC_IAU_2000_CTIES {
        obliquity_iau2000(t)
    } else if prec_long == PrecessionModel::IAU2000 {
        obliquity_iau2000(t)
    } else if prec_short == PrecessionModel::IAU2006 && t.abs() <= PREC_IAU_2006_CTIES {
        obliquity_iau2006(t)
    } else if prec_long == PrecessionModel::Newcomb {
        obliquity_newcomb(jd)
    } else if prec_long == PrecessionModel::IAU2006 {
        obliquity_iau2006(t)
    } else if prec_long == PrecessionModel::Bretagnon2003 {
        obliquity_bretagnon2003(t)
    } else if prec_long == PrecessionModel::Simon1994 {
        obliquity_simon1994(t)
    } else if prec_long == PrecessionModel::Williams1994 {
        obliquity_williams1994(t)
    } else if prec_long == PrecessionModel::Laskar1986 || prec_long == PrecessionModel::WillEpsLask
    {
        obliquity_laskar1986(t)
    } else if prec_long == PrecessionModel::Owen1990 {
        obliquity_owen1990(jd) * DEGTORAD
    } else {
        // Vondrák 2011 (default)
        let mut eps = obliquity_vondrak2011(jd);
        if flags.contains(CalcFlags::JPLHOR_APPROX) && models.jplhora_mode != JplHoraMode::V2 {
            eps += jpl_eps_correction(jd);
        }
        eps
    };

    Epsilon::new(eps)
}

// ---------------------------------------------------------------------------
// Individual obliquity models (all return radians)
// ---------------------------------------------------------------------------

fn obliquity_iau1976(t: f64) -> f64 {
    poly_eval(&[84381.448, -46.8150, -5.9e-4, 1.813e-3], t) * DEGTORAD / 3600.0
}

fn obliquity_iau2000(t: f64) -> f64 {
    poly_eval(&[84381.406, -46.84024, -5.9e-4, 1.813e-3], t) * DEGTORAD / 3600.0
}

fn obliquity_iau2006(t: f64) -> f64 {
    poly_eval(
        &[
            84381.406, -46.836769, -1.831e-4, 2.0034e-3, -5.76e-7, -4.34e-8,
        ],
        t,
    ) * DEGTORAD
        / 3600.0
}

fn obliquity_bretagnon2003(t: f64) -> f64 {
    poly_eval(
        &[
            84381.40880,
            -46.836051,
            -1.667e-4,
            1.99911e-3,
            -5.23e-7,
            -2.48e-8,
            -3e-11,
        ],
        t,
    ) * DEGTORAD
        / 3600.0
}

fn obliquity_simon1994(t: f64) -> f64 {
    poly_eval(
        &[84381.412, -46.80927, -1.52e-4, 1.9989e-3, -5.1e-7, 2.5e-8],
        t,
    ) * DEGTORAD
        / 3600.0
}

fn obliquity_williams1994(t: f64) -> f64 {
    poly_eval(&[84381.409, -46.833960, -1.74e-4, 2.0e-3, -1.0e-6], t) * DEGTORAD / 3600.0
}

fn obliquity_laskar1986(t: f64) -> f64 {
    let u = t / 10.0;
    poly_eval(
        &[
            84381.448, -468.093, -0.0155, 1.99925, -5.138e-3, -2.4967e-3, -3.905e-5, 7.12e-7,
            2.787e-7, 5.79e-9, 2.45e-10,
        ],
        u,
    ) * (DEGTORAD / 3600.0)
}

fn obliquity_newcomb(jd: f64) -> f64 {
    let tn = (jd - 2396758.0) / 36525.0;
    poly_eval(&[84451.68, -46.837, -0.0085, 0.0017], tn) * DEGTORAD / 3600.0
}

// ---------------------------------------------------------------------------
// Owen 1990 — Chebyshev polynomial evaluation
// ---------------------------------------------------------------------------

const OWEN_T0S: [f64; 5] = [-3392455.5, -470455.5, 2451544.5, 5373544.5, 8295544.5];

#[rustfmt::skip]
const OWEN_EPS0_COEF: [[f64; 10]; 5] = [
    [23.699391439256386, 5.2330816033981775e-1, -5.6259493384864815e-2, -8.2033318431602032e-3, 6.6774163554156385e-4, 2.4931584012812606e-5, -3.1313623302407878e-6, 2.0343814827951515e-7, 2.9182026615852936e-8, -4.1118760893281951e-9],
    [24.124759551704588, -1.2094875596566286e-1, -8.3914869653015218e-2, 3.5357075322387405e-3, 6.4557467824807032e-4, -2.5092064378707704e-5, -1.7631607274450848e-6, 1.3363622791424094e-7, 1.5577817511054047e-8, -2.4613907093017122e-9],
    [23.439103144206208, -4.9386077073143590e-1, -2.3965445283267805e-4, 8.6637485629656489e-3, -5.2828151901367600e-5, -4.3951004595359217e-5, -1.1058785949914705e-6, 6.2431490022621172e-8, 3.4725376218710764e-8, 1.3658853127005757e-9],
    [22.724671295125046, -1.6041813558650337e-1, 7.0646783888132504e-2, 1.4967806745062837e-3, -6.6857270989190734e-4, 5.7578378071604775e-6, 3.3738508454638728e-6, -2.2917813537654764e-7, -2.1019907929218137e-8, 4.3139832091694682e-9],
    [22.914636050333696, 3.2123508304962416e-1, 3.6633220173792710e-2, -5.9228324767696043e-3, -1.882379107379328e-4, 3.2274552870236244e-5, 4.9052463646336507e-7, -5.9064298731578425e-8, -2.0485712675098837e-8, -6.2163304813908160e-10],
];

fn owen_t0_icof(jd: f64) -> (f64, usize) {
    let mut t0 = OWEN_T0S[0];
    let mut icof = 0;
    for i in 1..5 {
        if jd >= (OWEN_T0S[i - 1] + OWEN_T0S[i]) / 2.0 {
            t0 = OWEN_T0S[i];
            icof = i;
        }
    }
    (t0, icof)
}

fn obliquity_owen1990(jd: f64) -> f64 {
    let (t0, icof) = owen_t0_icof(jd);
    let x = (jd - t0) / 36525.0 / 40.0;

    let mut tau = [0.0; 10];
    tau[1] = x;
    for i in 2..=9 {
        tau[i] = x * tau[i - 1];
    }

    let k = [
        1.0,
        tau[1],
        2.0 * tau[2] - 1.0,
        4.0 * tau[3] - 3.0 * tau[1],
        8.0 * tau[4] - 8.0 * tau[2] + 1.0,
        16.0 * tau[5] - 20.0 * tau[3] + 5.0 * tau[1],
        32.0 * tau[6] - 48.0 * tau[4] + 18.0 * tau[2] - 1.0,
        64.0 * tau[7] - 112.0 * tau[5] + 56.0 * tau[3] - 7.0 * tau[1],
        128.0 * tau[8] - 256.0 * tau[6] + 160.0 * tau[4] - 32.0 * tau[2] + 1.0,
        256.0 * tau[9] - 576.0 * tau[7] + 432.0 * tau[5] - 120.0 * tau[3] + 9.0 * tau[1],
    ];

    let coef = &OWEN_EPS0_COEF[icof];
    let mut eps = 0.0;
    for i in 0..10 {
        eps += k[i] * coef[i];
    }
    eps
}

// ---------------------------------------------------------------------------
// Vondrák 2011
// ---------------------------------------------------------------------------

const PEPOL_EPS: [f64; 4] = [84028.206305, 0.3624445, -0.00004039, -0.000000110];

#[rustfmt::skip]
const PEPER_PERIOD: [f64; 10] = [
    409.90, 396.15, 537.22, 402.90, 417.15, 288.92, 4043.00, 306.00, 277.00, 203.00,
];

#[rustfmt::skip]
const PEPER_COS_EPS: [f64; 10] = [
    753.872780, -247.805823, 379.471484, -53.880558, -90.109153,
    -353.600190, -63.115353, -28.248187, 17.703387, 38.911307,
];

#[rustfmt::skip]
const PEPER_SIN_EPS: [f64; 10] = [
    -1704.720302, -862.308358, 447.832178, -889.571909, 190.402846,
    -56.564991, -296.222622, -75.859952, 67.473503, 3.014055,
];

fn obliquity_vondrak2011(jd: f64) -> f64 {
    let t = (jd - J2000) / 36525.0;
    let w = TAU * t;

    let mut q = 0.0;
    for i in 0..10 {
        let a = w / PEPER_PERIOD[i];
        q += a.cos() * PEPER_COS_EPS[i] + a.sin() * PEPER_SIN_EPS[i];
    }

    let mut pw = 1.0;
    for &c in &PEPOL_EPS {
        q += c * pw;
        pw *= t;
    }

    q * (DEGTORAD / 3600.0)
}

// ---------------------------------------------------------------------------
// JPL Horizons obliquity correction
// ---------------------------------------------------------------------------

fn jpl_eps_correction(jd: f64) -> f64 {
    let tofs = (jd - DCOR_EPS_JPL_TJD0) / 365.25;
    let dofs = if tofs < 0.0 {
        DCOR_EPS_JPL[0]
    } else if tofs >= (DCOR_EPS_JPL.len() - 1) as f64 {
        DCOR_EPS_JPL[DCOR_EPS_JPL.len() - 1]
    } else {
        let t0 = tofs as usize;
        (tofs - t0 as f64) * (DCOR_EPS_JPL[t0] - DCOR_EPS_JPL[t0 + 1]) + DCOR_EPS_JPL[t0]
    };
    dofs / (1000.0 * 3600.0) * DEGTORAD
}
