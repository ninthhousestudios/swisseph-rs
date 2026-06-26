use crate::constants::*;
use crate::flags::CalcFlags;
use crate::math::owen_chebyshev_basis;
use crate::obliquity;
use crate::types::*;

const AS2R: f64 = DEGTORAD / 3600.0;
const EPS0_VONDRAK: f64 = 84381.406 * AS2R;

fn horner(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, &c| acc * x + c)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Vondrák 2011: general precession in longitude (pA)
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const PEPS_POL: [[f64; 2]; 4] = [
    [8134.017132, 84028.206305],
    [5043.0520035, 0.3624445],
    [-0.00710733, -0.00004039],
    [0.000000271, -0.000000110],
];

#[rustfmt::skip]
const PEPS_PER: [[f64; 10]; 5] = [
    [409.90, 396.15, 537.22, 402.90, 417.15, 288.92, 4043.00, 306.00, 277.00, 203.00],
    [-6908.287473, -3198.706291, 1453.674527, -857.748557, 1173.231614, -156.981465, 371.836550, -216.619040, 193.691479, 11.891524],
    [753.872780, -247.805823, 379.471484, -53.880558, -90.109153, -353.600190, -63.115353, -28.248187, 17.703387, 38.911307],
    [-2845.175469, 449.844989, -1255.915323, 886.736783, 418.887514, 997.912441, -240.979710, 76.541307, -36.788069, -170.964086],
    [-1704.720302, -862.308358, 447.832178, -889.571909, 190.402846, -56.564991, -296.222622, -75.859952, 67.473503, 3.014055],
];

/// General precession in longitude (pA) via Vondrák 2011 model. Returns radians.
pub fn ldp_peps(jd: f64) -> f64 {
    let t = (jd - J2000) / 36525.0;
    let mut p = 0.0;

    let w = TWOPI * t;
    for i in 0..10 {
        let a = w / PEPS_PER[0][i];
        let (s, c) = a.sin_cos();
        p += c * PEPS_PER[1][i] + s * PEPS_PER[3][i];
    }

    let mut pw = 1.0;
    for i in 0..4 {
        p += PEPS_POL[i][0] * pw;
        pw *= t;
    }

    p * AS2R
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn precess(
    pos: &mut [f64; 3],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    direction: PrecessionDirection,
) {
    if jd == J2000 {
        return;
    }
    let t = (jd - J2000) / 36525.0;
    let prec_long = models.prec_longterm;
    let prec_short = models.prec_shortterm;

    let is_jplhor = flags.contains(CalcFlags::DPSIDEPS_1980)
        || (flags.contains(CalcFlags::JPLHOR_APPROX)
            && models.jplhora_mode == JplHoraMode::V3
            && jd <= DPSI_DEPS_IAU1980_TJD0_HORIZONS);

    if is_jplhor {
        if jd > 2378131.5 && jd < 2525323.5 {
            precess_1(pos, jd, PrecessionModel::IAU1976, direction);
        } else {
            precess_3(pos, jd, flags, PrecessionModel::Owen1990, direction);
        }
    } else if prec_short == PrecessionModel::IAU1976 && t.abs() <= PREC_IAU_1976_CTIES {
        precess_1(pos, jd, PrecessionModel::IAU1976, direction);
    } else if prec_long == PrecessionModel::IAU1976 {
        precess_1(pos, jd, PrecessionModel::IAU1976, direction);
    } else if prec_short == PrecessionModel::IAU2000 && t.abs() <= PREC_IAU_2000_CTIES {
        precess_1(pos, jd, PrecessionModel::IAU2000, direction);
    } else if prec_long == PrecessionModel::IAU2000 {
        precess_1(pos, jd, PrecessionModel::IAU2000, direction);
    } else if prec_short == PrecessionModel::IAU2006 && t.abs() <= PREC_IAU_2006_CTIES {
        precess_1(pos, jd, PrecessionModel::IAU2006, direction);
    } else if prec_long == PrecessionModel::IAU2006 {
        precess_1(pos, jd, PrecessionModel::IAU2006, direction);
    } else if prec_long == PrecessionModel::Bretagnon2003 {
        precess_1(pos, jd, PrecessionModel::Bretagnon2003, direction);
    } else if prec_long == PrecessionModel::Newcomb {
        precess_1(pos, jd, PrecessionModel::Newcomb, direction);
    } else if prec_long == PrecessionModel::Laskar1986 {
        precess_2(
            pos,
            jd,
            flags,
            models,
            PrecessionModel::Laskar1986,
            direction,
        );
    } else if prec_long == PrecessionModel::Simon1994 {
        precess_2(
            pos,
            jd,
            flags,
            models,
            PrecessionModel::Simon1994,
            direction,
        );
    } else if prec_long == PrecessionModel::Williams1994
        || prec_long == PrecessionModel::WillEpsLask
    {
        precess_2(
            pos,
            jd,
            flags,
            models,
            PrecessionModel::Williams1994,
            direction,
        );
    } else if prec_long == PrecessionModel::Owen1990 {
        precess_3(pos, jd, flags, PrecessionModel::Owen1990, direction);
    } else {
        precess_3(pos, jd, flags, PrecessionModel::Vondrak2011, direction);
    }
}

// ---------------------------------------------------------------------------
// precess_1: Euler angle method
// ---------------------------------------------------------------------------

fn precess_1(pos: &mut [f64; 3], jd: f64, model: PrecessionModel, direction: PrecessionDirection) {
    let t = (jd - J2000) / 36525.0;

    let (zz, z, th) = match model {
        PrecessionModel::IAU1976 => {
            let zz = ((0.017998 * t + 0.30188) * t + 2306.2181) * t * DEGTORAD / 3600.0;
            let z = ((0.018203 * t + 1.09468) * t + 2306.2181) * t * DEGTORAD / 3600.0;
            let th = ((-0.041833 * t - 0.42665) * t + 2004.3109) * t * DEGTORAD / 3600.0;
            (zz, z, th)
        }
        PrecessionModel::IAU2000 => {
            let zz = (((((-0.0000002 * t - 0.0000327) * t + 0.0179663) * t + 0.3019015) * t
                + 2306.0809506)
                * t
                + 2.5976176)
                * DEGTORAD
                / 3600.0;
            let z = (((((-0.0000003 * t - 0.000047) * t + 0.0182237) * t + 1.0947790) * t
                + 2306.0803226)
                * t
                - 2.5976176)
                * DEGTORAD
                / 3600.0;
            let th = ((((-0.0000001 * t - 0.0000601) * t - 0.0418251) * t - 0.4269353) * t
                + 2004.1917476)
                * t
                * DEGTORAD
                / 3600.0;
            (zz, z, th)
        }
        PrecessionModel::IAU2006 => {
            let zz = (((((-0.0000003173 * t - 0.000005971) * t + 0.01801828) * t + 0.2988499) * t
                + 2306.083227)
                * t
                + 2.650545)
                * DEGTORAD
                / 3600.0;
            let z = (((((-0.0000002904 * t - 0.000028596) * t + 0.01826837) * t + 1.0927348) * t
                + 2306.077181)
                * t
                - 2.650545)
                * DEGTORAD
                / 3600.0;
            let th = ((((-0.00000011274 * t - 0.000007089) * t - 0.04182264) * t - 0.4294934) * t
                + 2004.191903)
                * t
                * DEGTORAD
                / 3600.0;
            (zz, z, th)
        }
        PrecessionModel::Bretagnon2003 => {
            let zz = ((((((-0.00000000013 * t - 0.0000003040) * t - 0.000005708) * t
                + 0.01801752)
                * t
                + 0.3023262)
                * t
                + 2306.080472)
                * t
                + 2.72767)
                * DEGTORAD
                / 3600.0;
            let z = ((((((-0.00000000005 * t - 0.0000002486) * t - 0.000028276) * t
                + 0.01826676)
                * t
                + 1.0956768)
                * t
                + 2306.076070)
                * t
                - 2.72767)
                * DEGTORAD
                / 3600.0;
            let th = ((((((0.000000000009 * t + 0.00000000036) * t - 0.0000001127) * t
                - 0.000007291)
                * t
                - 0.04182364)
                * t
                - 0.4266980)
                * t
                + 2004.190936)
                * t
                * DEGTORAD
                / 3600.0;
            (zz, z, th)
        }
        PrecessionModel::Newcomb => {
            let mills = 365242.198782;
            let t1 = (J2000 - B1850) / mills;
            let t2 = (jd - B1850) / mills;
            let tn = t2 - t1;
            let tn2 = tn * tn;
            let tn3 = tn2 * tn;
            let z1 = 23035.5548 + 139.720 * t1 + 0.069 * t1 * t1;
            let zz = (z1 * tn + (30.242 - 0.269 * t1) * tn2 + 17.996 * tn3) * (DEGTORAD / 3600.0);
            let z = (z1 * tn + (109.478 - 0.387 * t1) * tn2 + 18.324 * tn3) * (DEGTORAD / 3600.0);
            let th = ((20051.125 - 85.294 * t1 - 0.365 * t1 * t1) * tn
                + (-42.647 - 0.365 * t1) * tn2
                - 41.802 * tn3)
                * (DEGTORAD / 3600.0);
            (zz, z, th)
        }
        _ => return,
    };

    let (sin_th, cos_th) = th.sin_cos();
    let (sin_zz, cos_zz) = zz.sin_cos();
    let (sin_z, cos_z) = z.sin_cos();
    let a = cos_zz * cos_th;
    let b = sin_zz * cos_th;

    let x = match direction {
        PrecessionDirection::J2000ToDate => [
            (a * cos_z - sin_zz * sin_z) * pos[0]
                - (b * cos_z + cos_zz * sin_z) * pos[1]
                - sin_th * cos_z * pos[2],
            (a * sin_z + sin_zz * cos_z) * pos[0]
                - (b * sin_z - cos_zz * cos_z) * pos[1]
                - sin_th * sin_z * pos[2],
            cos_zz * sin_th * pos[0] - sin_zz * sin_th * pos[1] + cos_th * pos[2],
        ],
        PrecessionDirection::DateToJ2000 => [
            (a * cos_z - sin_zz * sin_z) * pos[0]
                + (a * sin_z + sin_zz * cos_z) * pos[1]
                + cos_zz * sin_th * pos[2],
            -(b * cos_z + cos_zz * sin_z) * pos[0]
                - (b * sin_z - cos_zz * cos_z) * pos[1]
                - sin_zz * sin_th * pos[2],
            -sin_th * cos_z * pos[0] - sin_th * sin_z * pos[1] + cos_th * pos[2],
        ],
    };

    *pos = x;
}

// ---------------------------------------------------------------------------
// precess_2: Ecliptic long-term method
// Coefficients in descending order (highest-degree first) matching C source.
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const PA_COF_WILLIAMS: [f64; 10] = [
    -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
    -0.235316, 0.076, 110.5407, 50287.70000,
];
#[rustfmt::skip]
const NODE_COF_WILLIAMS: [f64; 11] = [
    6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10,
    -3.54e-9, -1.8103e-7, 1.26e-7, 7.436169e-5, -0.04207794833, 3.052115282424,
];
#[rustfmt::skip]
const INCL_COF_WILLIAMS: [f64; 11] = [
    1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11,
    -5.4000441e-11, 1.32115526e-9, -6.012e-7, -1.62442e-5, 0.00227850649, 0.0,
];

#[rustfmt::skip]
const PA_COF_SIMON: [f64; 10] = [
    -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
    -0.235316, 0.07732, 111.2022, 50288.200,
];
#[rustfmt::skip]
const NODE_COF_SIMON: [f64; 11] = [
    6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10,
    -3.54e-9, -1.8103e-7, 2.579e-8, 7.4379679e-5, -0.0420782900, 3.0521126906,
];
#[rustfmt::skip]
const INCL_COF_SIMON: [f64; 11] = [
    1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11,
    -5.4000441e-11, 1.32115526e-9, -5.99908e-7, -1.624383e-5, 0.002278492868, 0.0,
];

#[rustfmt::skip]
const PA_COF_LASKAR: [f64; 10] = [
    -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
    -0.235316, 0.07732, 111.1971, 50290.966,
];
#[rustfmt::skip]
const NODE_COF_LASKAR: [f64; 11] = [
    6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 6.3190131e-10,
    -3.48388152e-9, -1.813065896e-7, 2.75036225e-8, 7.4394531426e-5,
    -0.042078604317, 3.052112654975,
];
#[rustfmt::skip]
const INCL_COF_LASKAR: [f64; 11] = [
    1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11,
    -5.4000441e-11, 1.32115526e-9, -5.998737027e-7, -1.6242797091e-5,
    0.002278495537, 0.0,
];

fn precess_2(
    pos: &mut [f64; 3],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    model: PrecessionModel,
    direction: PrecessionDirection,
) {
    let (pa_cof, node_cof, incl_cof): (&[f64], &[f64], &[f64]) = match model {
        PrecessionModel::Laskar1986 => (&PA_COF_LASKAR, &NODE_COF_LASKAR, &INCL_COF_LASKAR),
        PrecessionModel::Simon1994 => (&PA_COF_SIMON, &NODE_COF_SIMON, &INCL_COF_SIMON),
        _ => (&PA_COF_WILLIAMS, &NODE_COF_WILLIAMS, &INCL_COF_WILLIAMS),
    };

    let t = (jd - J2000) / 36525.0;

    // Step 1: rotate from equator to ecliptic at source epoch
    let eps = match direction {
        PrecessionDirection::DateToJ2000 => obliquity::obliquity(jd, flags, models).eps,
        PrecessionDirection::J2000ToDate => obliquity::obliquity(J2000, flags, models).eps,
    };
    let (sin_eps, cos_eps) = eps.sin_cos();
    let mut x = [
        pos[0],
        cos_eps * pos[1] + sin_eps * pos[2],
        -sin_eps * pos[1] + cos_eps * pos[2],
    ];

    // Evaluate polynomials in T = t/10 (thousands of years)
    let t10 = t / 10.0;
    let pa = horner(pa_cof, t10) * ((DEGTORAD / 3600.0) * t10);
    let w = horner(node_cof, t10);
    let incl = horner(incl_cof, t10);

    // Step 2: Z rotation to the node
    let angle = match direction {
        PrecessionDirection::DateToJ2000 => w + pa,
        PrecessionDirection::J2000ToDate => w,
    };
    let (sin_a, cos_a) = angle.sin_cos();
    let tmp = cos_a * x[0] + sin_a * x[1];
    x[1] = -sin_a * x[0] + cos_a * x[1];
    x[0] = tmp;

    // Step 3: X rotation by inclination
    let incl_angle = match direction {
        PrecessionDirection::DateToJ2000 => -incl,
        PrecessionDirection::J2000ToDate => incl,
    };
    let (sin_a, cos_a) = incl_angle.sin_cos();
    let tmp = cos_a * x[1] + sin_a * x[2];
    x[2] = -sin_a * x[1] + cos_a * x[2];
    x[1] = tmp;

    // Step 4: Z rotation back from the node
    let angle = match direction {
        PrecessionDirection::DateToJ2000 => -w,
        PrecessionDirection::J2000ToDate => -w - pa,
    };
    let (sin_a, cos_a) = angle.sin_cos();
    let tmp = cos_a * x[0] + sin_a * x[1];
    x[1] = -sin_a * x[0] + cos_a * x[1];
    x[0] = tmp;

    // Step 5: rotate from ecliptic to equator at target epoch
    let eps = match direction {
        PrecessionDirection::DateToJ2000 => obliquity::obliquity(J2000, flags, models).eps,
        PrecessionDirection::J2000ToDate => obliquity::obliquity(jd, flags, models).eps,
    };
    let (sin_eps, cos_eps) = eps.sin_cos();
    let tmp = cos_eps * x[1] - sin_eps * x[2];
    x[2] = sin_eps * x[1] + cos_eps * x[2];
    x[1] = tmp;

    *pos = x;
}

// ---------------------------------------------------------------------------
// precess_3: Matrix dispatch (Vondrák 2011, Owen 1990)
// ---------------------------------------------------------------------------

fn precess_3(
    pos: &mut [f64; 3],
    jd: f64,
    flags: CalcFlags,
    model: PrecessionModel,
    direction: PrecessionDirection,
) {
    let rp = if model == PrecessionModel::Owen1990 {
        owen_pre_matrix(jd, flags)
    } else {
        pre_pmat(jd)
    };

    let x = match direction {
        PrecessionDirection::J2000ToDate => [
            pos[0] * rp[0][0] + pos[1] * rp[0][1] + pos[2] * rp[0][2],
            pos[0] * rp[1][0] + pos[1] * rp[1][1] + pos[2] * rp[1][2],
            pos[0] * rp[2][0] + pos[1] * rp[2][1] + pos[2] * rp[2][2],
        ],
        PrecessionDirection::DateToJ2000 => [
            pos[0] * rp[0][0] + pos[1] * rp[1][0] + pos[2] * rp[2][0],
            pos[0] * rp[0][1] + pos[1] * rp[1][1] + pos[2] * rp[2][1],
            pos[0] * rp[0][2] + pos[1] * rp[1][2] + pos[2] * rp[2][2],
        ],
    };

    *pos = x;
}

// ---------------------------------------------------------------------------
// Vondrák 2011: ecliptic pole
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const PECL_POL: [[f64; 2]; 4] = [
    [5851.607687, -1600.886300],
    [-0.1189000, 1.1689818],
    [-0.00028913, -0.00000020],
    [0.000000101, -0.000000437],
];

#[rustfmt::skip]
const PECL_PER: [[f64; 8]; 5] = [
    [708.15, 2309.0, 1620.0, 492.2, 1183.0, 622.0, 882.0, 547.0],
    [-5486.751211, -17.127623, -617.517403, 413.44294, 78.614193, -180.732815, -87.676083, 46.140315],
    [-684.66156, 2446.28388, 399.671049, -356.652376, -186.387003, -316.80007, 198.296701, 101.135679],
    [667.66673, -2354.886252, -428.152441, 376.202861, 184.778874, 335.321713, -185.138669, -120.97283],
    [-5523.863691, -549.74745, -310.998056, 421.535876, -36.776172, -145.278396, -34.74445, 22.885731],
];

fn pre_pecl(jd: f64) -> [f64; 3] {
    let t = (jd - J2000) / 36525.0;
    let mut p = 0.0;
    let mut q = 0.0;

    let w = TWOPI * t;
    for i in 0..8 {
        let a = w / PECL_PER[0][i];
        let (s, c) = a.sin_cos();
        p += c * PECL_PER[1][i] + s * PECL_PER[3][i];
        q += c * PECL_PER[2][i] + s * PECL_PER[4][i];
    }

    let mut pw = 1.0;
    for i in 0..4 {
        p += PECL_POL[i][0] * pw;
        q += PECL_POL[i][1] * pw;
        pw *= t;
    }

    p *= AS2R;
    q *= AS2R;

    let z2 = 1.0 - p * p - q * q;
    let z = if z2 < 0.0 { 0.0 } else { z2.sqrt() };
    let (s, c) = EPS0_VONDRAK.sin_cos();
    [p, -q * c - z * s, -q * s + z * c]
}

// ---------------------------------------------------------------------------
// Vondrák 2011: equator pole
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const PEQU_POL: [[f64; 2]; 4] = [
    [5453.282155, -73750.930350],
    [0.4252841, -0.7675452],
    [-0.00037173, -0.00018725],
    [-0.000000152, 0.000000231],
];

#[rustfmt::skip]
const PEQU_PER: [[f64; 14]; 5] = [
    [256.75, 708.15, 274.2, 241.45, 2309.0, 492.2, 396.1, 288.9, 231.1, 1610.0, 620.0, 157.87, 220.3, 1200.0],
    [-819.940624, -8444.676815, 2600.009459, 2755.17563, -167.659835, 871.855056, 44.769698, -512.313065, -819.415595, -538.071099, -189.793622, -402.922932, 179.516345, -9.814756],
    [75004.344875, 624.033993, 1251.136893, -1102.212834, -2660.66498, 699.291817, 153.16722, -950.865637, 499.754645, -145.18821, 558.116553, -23.923029, -165.405086, 9.344131],
    [81491.287984, 787.163481, 1251.296102, -1257.950837, -2966.79973, 639.744522, 131.600209, -445.040117, 584.522874, -89.756563, 524.42963, -13.549067, -210.157124, -44.919798],
    [1558.515853, 7774.939698, -2219.534038, -2523.969396, 247.850422, -846.485643, -1393.124055, 368.526116, 749.045012, 444.704518, 235.934465, 374.049623, -171.33018, -22.899655],
];

fn pre_pequ(jd: f64) -> [f64; 3] {
    let t = (jd - J2000) / 36525.0;
    let mut x = 0.0;
    let mut y = 0.0;

    let w = TWOPI * t;
    for i in 0..14 {
        let a = w / PEQU_PER[0][i];
        let (s, c) = a.sin_cos();
        x += c * PEQU_PER[1][i] + s * PEQU_PER[3][i];
        y += c * PEQU_PER[2][i] + s * PEQU_PER[4][i];
    }

    let mut pw = 1.0;
    for i in 0..4 {
        x += PEQU_POL[i][0] * pw;
        y += PEQU_POL[i][1] * pw;
        pw *= t;
    }

    x *= AS2R;
    y *= AS2R;

    let w2 = x * x + y * y;
    let z = if w2 < 1.0 { (1.0 - w2).sqrt() } else { 0.0 };
    [x, y, z]
}

// ---------------------------------------------------------------------------
// Vondrák 2011: precession matrix from pole vectors
// ---------------------------------------------------------------------------

fn pre_pmat(jd: f64) -> [[f64; 3]; 3] {
    let peqr = pre_pequ(jd);
    let pecl = pre_pecl(jd);

    let v = cross(peqr, pecl);
    let w = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let eqx = [v[0] / w, v[1] / w, v[2] / w];
    let v2 = cross(peqr, eqx);

    [eqx, v2, peqr]
}

// ---------------------------------------------------------------------------
// Owen 1990: precession matrix from Chebyshev angles
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const OWEN_PSIA_COEF: [[f64; 10]; 5] = [
    [-218.57864954903122, 51.752257487741612, 1.3304715765661958e-1, 9.2048123521890745e-2, -6.0877528127241278e-3, -7.0013893644531700e-5, -4.9217728385458495e-5, -1.8578234189053723e-6, 7.4396426162029877e-7, -5.9157528981843864e-9],
    [-111.94350527506128, 55.175558131675861, 4.7366115762797613e-1, -4.7701750975398538e-2, -9.2445765329325809e-3, 7.0962838707454917e-4, 1.5140455277814658e-4, -7.7813159018954928e-7, -2.4729402281953378e-6, -1.0898887008726418e-7],
    [-2.041452011529441e-1, 55.969995858494106, -1.9295093699770936e-1, -5.6819574830421158e-3, 1.1073687302518981e-2, -9.0868489896815619e-5, -1.1999773777895820e-4, 9.9748697306154409e-6, 5.7911493603430550e-7, -2.3647526839778175e-7],
    [111.61366860604471, 56.404525305162447, 4.4403302410703782e-1, 7.1490030578883907e-2, -4.9184559079790816e-3, -1.3912698949042046e-3, -6.8490613661884005e-5, 1.2394328562905297e-6, 1.7719847841480384e-6, 2.4889095220628068e-7],
    [228.40683531269390, 60.056143904919826, 2.9583200718478960e-2, -1.5710838319490748e-1, -7.0017356811600801e-3, 3.3009615142224537e-3, 2.0318123852537664e-4, -6.5840216067828310e-5, -5.9077673352976155e-6, 1.3983942185303064e-6],
];

#[rustfmt::skip]
const OWEN_OMA_COEF: [[f64; 10]; 5] = [
    [25.541291140949806, 2.377889511272162e-1, -3.7337334723142133e-1, 2.4579295485161534e-2, 4.3840999514263623e-3, -3.1126873333599556e-4, -9.8443045771748915e-6, -7.9403103080496923e-7, 1.0840116743893556e-9, 9.2865105216887919e-9],
    [24.429357654237926, -9.5205745947740161e-1, 8.6738296270534816e-2, 3.0061543426062955e-2, -4.1532480523019988e-3, -3.7920928393860939e-4, 3.5117012399609737e-5, 4.6811877283079217e-6, -8.1836046585546861e-8, -6.1803706664211173e-8],
    [23.450465062489337, -9.7259278279739817e-2, 1.1082286925130981e-2, -3.1469883339372219e-2, -1.0041906996819648e-4, 5.6455168475133958e-4, -8.4403910211030209e-6, -3.8269157371098435e-6, 3.1422585261198437e-7, 9.3481729116773404e-9],
    [22.581778052947806, -8.7069701538602037e-1, -9.8140710050197307e-2, 2.6025931340678079e-2, 4.8165322168786755e-3, -1.906558772193363e-4, -4.6838759635421777e-5, -1.6608525315998471e-6, -3.2347811293516124e-8, 2.8104728109642000e-9],
    [21.518861835737142, 2.0494789509441385e-1, 3.5193604846503161e-1, 1.5305977982348925e-2, -7.5015367726336455e-3, -4.0322553186065610e-4, 1.0655320434844041e-4, 7.1792339586935752e-6, -1.603874697543020e-6, -1.613563462813512e-7],
];

#[rustfmt::skip]
const OWEN_CHIA_COEF: [[f64; 10]; 5] = [
    [8.2378850337329404e-1, -3.7443109739678667, 4.0143936898854026e-1, 8.1822830214590811e-2, -8.5978790792656293e-3, -2.8350488448426132e-5, -4.2474671728156727e-5, -1.6214840884656678e-6, 7.8560442001953050e-7, -1.032016641696707e-8],
    [-2.1726062070318606, 7.8470515033132925e-1, 4.4044931004195718e-1, -8.0671247169971653e-2, -8.9672662444325007e-3, 9.2248978383109719e-4, 1.5143472266372874e-4, -1.6387009056475679e-6, -2.4405558979328144e-6, -1.0148113464009015e-7],
    [-4.8518673570735556e-1, 1.0016737299946743e-1, -4.7074888613099918e-1, -5.8604054305076092e-3, 1.4300208240553435e-2, -6.7127991650300028e-5, -1.3703764889645475e-4, 9.0505213684444634e-6, 6.0368690647808607e-7, -2.2135404747652171e-7],
    [-2.0950740076326087, -9.4447359463206877e-1, 4.0940512860493755e-1, 1.0261699700263508e-1, -5.3133241571955160e-3, -1.6634631550720911e-3, -5.9477519536647907e-5, 2.9651387319208926e-6, 1.6434499452070584e-6, 2.3720647656961084e-7],
    [6.3315163285678715e-1, 3.5241082918420464, 2.1223076605364606e-1, -1.5648122502767368e-1, -9.1964075390801980e-3, 3.3896161239812411e-3, 2.1485178626085787e-4, -6.6261759864793735e-5, -5.9257969712852667e-6, 1.3918759086160525e-6],
];

fn owen_pre_matrix(jd: f64, flags: CalcFlags) -> [[f64; 3]; 3] {
    let (icof, k) = owen_chebyshev_basis(jd);

    let mut psia = 0.0;
    let mut oma = 0.0;
    let mut chia = 0.0;
    for i in 0..10 {
        psia += k[i] * OWEN_PSIA_COEF[icof][i];
        oma += k[i] * OWEN_OMA_COEF[icof][i];
        chia += k[i] * OWEN_CHIA_COEF[icof][i];
    }

    if flags.intersects(CalcFlags::DPSIDEPS_1980 | CalcFlags::JPLHOR_APPROX) {
        psia += -0.000018560;
    }

    let eps0 = 84381.448 / 3600.0 * DEGTORAD;
    psia *= DEGTORAD;
    chia *= DEGTORAD;
    oma *= DEGTORAD;

    let (sin_eps0, cos_eps0) = eps0.sin_cos();
    let (sin_chia, cos_chia) = chia.sin_cos();
    let (sin_psia, cos_psia) = psia.sin_cos();
    let (sin_oma, cos_oma) = oma.sin_cos();

    [
        [
            cos_chia * cos_psia + sin_chia * cos_oma * sin_psia,
            (-cos_chia * sin_psia + sin_chia * cos_oma * cos_psia) * cos_eps0
                + sin_chia * sin_oma * sin_eps0,
            (-cos_chia * sin_psia + sin_chia * cos_oma * cos_psia) * sin_eps0
                - sin_chia * sin_oma * cos_eps0,
        ],
        [
            -sin_chia * cos_psia + cos_chia * cos_oma * sin_psia,
            (sin_chia * sin_psia + cos_chia * cos_oma * cos_psia) * cos_eps0
                + cos_chia * sin_oma * sin_eps0,
            (sin_chia * sin_psia + cos_chia * cos_oma * cos_psia) * sin_eps0
                - cos_chia * sin_oma * cos_eps0,
        ],
        [
            sin_oma * sin_psia,
            sin_oma * cos_psia * cos_eps0 - cos_oma * sin_eps0,
            sin_oma * cos_psia * sin_eps0 + cos_oma * cos_eps0,
        ],
    ]
}
