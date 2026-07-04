use crate::constants::*;

const EFF_ARR: [(f64, f64); 101] = [
    (1.000, 1.000000),
    (0.990, 0.999979),
    (0.980, 0.999940),
    (0.970, 0.999881),
    (0.960, 0.999811),
    (0.950, 0.999724),
    (0.940, 0.999622),
    (0.930, 0.999497),
    (0.920, 0.999354),
    (0.910, 0.999192),
    (0.900, 0.999000),
    (0.890, 0.998786),
    (0.880, 0.998535),
    (0.870, 0.998242),
    (0.860, 0.997919),
    (0.850, 0.997571),
    (0.840, 0.997198),
    (0.830, 0.996792),
    (0.820, 0.996316),
    (0.810, 0.995791),
    (0.800, 0.995226),
    (0.790, 0.994625),
    (0.780, 0.993991),
    (0.770, 0.993326),
    (0.760, 0.992598),
    (0.750, 0.991770),
    (0.740, 0.990873),
    (0.730, 0.989919),
    (0.720, 0.988912),
    (0.710, 0.987856),
    (0.700, 0.986755),
    (0.690, 0.985610),
    (0.680, 0.984398),
    (0.670, 0.982986),
    (0.660, 0.981437),
    (0.650, 0.979779),
    (0.640, 0.978024),
    (0.630, 0.976182),
    (0.620, 0.974256),
    (0.610, 0.972253),
    (0.600, 0.970174),
    (0.590, 0.968024),
    (0.580, 0.965594),
    (0.570, 0.962797),
    (0.560, 0.959758),
    (0.550, 0.956515),
    (0.540, 0.953088),
    (0.530, 0.949495),
    (0.520, 0.945741),
    (0.510, 0.941838),
    (0.500, 0.937790),
    (0.490, 0.933563),
    (0.480, 0.928668),
    (0.470, 0.923288),
    (0.460, 0.917527),
    (0.450, 0.911432),
    (0.440, 0.905035),
    (0.430, 0.898353),
    (0.420, 0.891022),
    (0.410, 0.882940),
    (0.400, 0.874312),
    (0.390, 0.865206),
    (0.380, 0.855423),
    (0.370, 0.844619),
    (0.360, 0.833074),
    (0.350, 0.820876),
    (0.340, 0.808031),
    (0.330, 0.793962),
    (0.320, 0.778931),
    (0.310, 0.763021),
    (0.300, 0.745815),
    (0.290, 0.727557),
    (0.280, 0.708234),
    (0.270, 0.687583),
    (0.260, 0.665741),
    (0.250, 0.642597),
    (0.240, 0.618252),
    (0.230, 0.592586),
    (0.220, 0.565747),
    (0.210, 0.537697),
    (0.200, 0.508554),
    (0.190, 0.478420),
    (0.180, 0.447322),
    (0.170, 0.415454),
    (0.160, 0.382892),
    (0.150, 0.349955),
    (0.140, 0.316691),
    (0.130, 0.283565),
    (0.120, 0.250431),
    (0.110, 0.218327),
    (0.100, 0.186794),
    (0.090, 0.156287),
    (0.080, 0.128421),
    (0.070, 0.102237),
    (0.060, 0.077393),
    (0.050, 0.054833),
    (0.040, 0.036361),
    (0.030, 0.020953),
    (0.020, 0.009645),
    (0.010, 0.002767),
    (0.000, 0.000000),
];

/// Effective mass fraction of the Sun for gravitational light deflection, as a function of
/// the impact-parameter ratio `r` (`sin(a)/sin(a_sun)`), interpolated from the `EFF_ARR` table.
/// Port of C `meff`.
pub fn meff(r: f64) -> f64 {
    if r <= 0.0 {
        return 0.0;
    }
    if r >= 1.0 {
        return 1.0;
    }
    let i = EFF_ARR.iter().position(|(ri, _)| *ri <= r).unwrap();
    let (r0, m0) = EFF_ARR[i - 1];
    let (r1, m1) = EFF_ARR[i];
    let f = (r - r0) / (r1 - r0);
    m0 + f * (m1 - m0)
}

/// Apply annual aberration of light to a geocentric position (and, if `has_speed`, speed) in
/// place. `xx` is `[x, y, z, vx, vy, vz]` (AU, AU/day); `earth_vel` is Earth's barycentric
/// velocity (AU/day). Port of C `swi_aberr_light`.
pub fn aberr_light(xx: &mut [f64; 6], earth_vel: &[f64; 3], has_speed: bool) {
    let mut v = [0.0; 3];
    let mut v2 = 0.0;
    for i in 0..3 {
        v[i] = earth_vel[i] / 24.0 / 3600.0 / CLIGHT * AUNIT;
        v2 += v[i] * v[i];
    }

    let orig = [xx[0], xx[1], xx[2]];
    let ru = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
    let b_1 = (1.0 - v2).sqrt();

    let mut f1 = 0.0;
    for i in 0..3 {
        f1 += xx[i] * v[i];
    }
    f1 /= ru;
    let f2 = 1.0 + f1 / (1.0 + b_1);

    for i in 0..3 {
        xx[i] = (b_1 * xx[i] + f2 * ru * v[i]) / (1.0 + f1);
    }

    if has_speed {
        let intv = PLAN_SPEED_INTV;
        let mut u_prev = [0.0; 3];
        for i in 0..3 {
            u_prev[i] = orig[i] - intv * xx[i + 3];
        }
        let ru_prev =
            (u_prev[0] * u_prev[0] + u_prev[1] * u_prev[1] + u_prev[2] * u_prev[2]).sqrt();
        let mut f1_prev = 0.0;
        for i in 0..3 {
            f1_prev += u_prev[i] * v[i];
        }
        f1_prev /= ru_prev;
        let f2_prev = 1.0 + f1_prev / (1.0 + b_1);
        for i in 0..3 {
            let xx2_i = (b_1 * u_prev[i] + f2_prev * ru_prev * v[i]) / (1.0 + f1_prev);
            let dx1 = xx[i] - orig[i];
            let dx2 = xx2_i - u_prev[i];
            xx[i + 3] += (dx1 - dx2) / intv;
        }
    }
}

struct DeflectOutput {
    result: [f64; 3],
    ru: f64,
    u_normalized: [f64; 3],
}

fn deflect_position(u_in: &[f64; 3], earth_pos: &[f64; 3], planet_pos: &[f64; 3]) -> DeflectOutput {
    let ru = (u_in[0] * u_in[0] + u_in[1] * u_in[1] + u_in[2] * u_in[2]).sqrt();
    let re =
        (earth_pos[0] * earth_pos[0] + earth_pos[1] * earth_pos[1] + earth_pos[2] * earth_pos[2])
            .sqrt();
    let rq = (planet_pos[0] * planet_pos[0]
        + planet_pos[1] * planet_pos[1]
        + planet_pos[2] * planet_pos[2])
        .sqrt();

    let u = [u_in[0] / ru, u_in[1] / ru, u_in[2] / ru];
    let e = [earth_pos[0] / re, earth_pos[1] / re, earth_pos[2] / re];
    let q = [planet_pos[0] / rq, planet_pos[1] / rq, planet_pos[2] / rq];

    let uq = u[0] * q[0] + u[1] * q[1] + u[2] * q[2];
    let ue = u[0] * e[0] + u[1] * e[1] + u[2] * e[2];
    let qe = q[0] * e[0] + q[1] * e[1] + q[2] * e[2];

    let sina = (1.0 - ue * ue).sqrt();
    let sin_sunr = SUN_RADIUS / re;
    let meff_fact = if sina < sin_sunr {
        meff(sina / sin_sunr)
    } else {
        1.0
    };

    let g1 = 2.0 * HELGRAVCONST * meff_fact / CLIGHT / CLIGHT / AUNIT / re;
    let g2 = 1.0 + qe;

    let corr = g1 / g2;
    DeflectOutput {
        result: [
            ru * (u[0] + corr * (uq * e[0] - ue * q[0])),
            ru * (u[1] + corr * (uq * e[1] - ue * q[1])),
            ru * (u[2] + corr * (uq * e[2] - ue * q[2])),
        ],
        ru,
        u_normalized: u,
    }
}

/// Apply relativistic gravitational deflection of light by the Sun to a geocentric position
/// (and, if `has_speed`, speed) in place. `xx` is `[x, y, z, vx, vy, vz]` (AU, AU/day);
/// `earth_helio`/`planet_helio` are Earth's/the target body's heliocentric state. Port of C
/// `swi_deflect_light`.
pub fn deflect_light(
    xx: &mut [f64; 6],
    earth_helio: &[f64; 6],
    planet_helio: &[f64; 6],
    has_speed: bool,
) {
    let u_in = [xx[0], xx[1], xx[2]];
    let e_pos = [earth_helio[0], earth_helio[1], earth_helio[2]];
    let q_pos = [planet_helio[0], planet_helio[1], planet_helio[2]];

    let out1 = deflect_position(&u_in, &e_pos, &q_pos);

    if has_speed {
        let dtsp = -DEFL_SPEED_INTV;
        let mut u_pert = [0.0; 3];
        let mut e_pert = [0.0; 3];
        let mut q_pert = [0.0; 3];
        for i in 0..3 {
            u_pert[i] = xx[i] - dtsp * xx[i + 3];
            e_pert[i] = earth_helio[i] - dtsp * earth_helio[i + 3];
            q_pert[i] = u_pert[i] + earth_helio[i] - dtsp * earth_helio[i + 3];
        }
        let out2 = deflect_position(&u_pert, &e_pert, &q_pert);
        for i in 0..3 {
            let dx1 = out1.result[i] - xx[i];
            let dx2 = out2.result[i] - out2.u_normalized[i] * out2.ru;
            xx[i + 3] += (dx1 - dx2) / dtsp;
        }
    }

    xx[0] = out1.result[0];
    xx[1] = out1.result[1];
    xx[2] = out1.result[2];
}
