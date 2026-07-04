use crate::constants::*;
use crate::flags::CalcFlags;
use crate::math::{
    cartesian_to_polar, cartesian_to_polar_with_speed, polar_to_cartesian,
    polar_to_cartesian_with_speed,
};
use crate::types::*;

const DCOR_RA_JPL_TJD0: f64 = 2437846.5;

#[rustfmt::skip]
const DCOR_RA_JPL: [f64; 51] = [
    -51.257, -51.103, -51.065, -51.503, -51.224, -50.796, -51.161, -51.181,
    -50.932, -51.064, -51.182, -51.386, -51.416, -51.428, -51.586, -51.766, -52.038, -52.370,
    -52.553, -52.397, -52.340, -52.676, -52.348, -51.964, -52.444, -52.364, -51.988, -52.212,
    -52.370, -52.523, -52.541, -52.496, -52.590, -52.629, -52.788, -53.014, -53.053, -52.902,
    -52.850, -53.087, -52.635, -52.185, -52.588, -52.292, -51.796, -51.961, -52.055, -52.134,
    -52.165, -52.141, -52.255,
];

// Verbatim frame-bias matrix from the C source; full digits preserved.
#[allow(clippy::excessive_precision)]
#[rustfmt::skip]
const BIAS_IAU2006: [[f64; 3]; 3] = [
    [ 0.99999999999999412,  0.00000007078368695, -0.00000008056214212],
    [-0.00000007078368961,  0.99999999999999700, -0.00000003306427981],
    [ 0.00000008056213978,  0.00000003306428553,  0.99999999999999634],
];

#[rustfmt::skip]
const BIAS_IAU2000: [[f64; 3]; 3] = [
    [ 0.9999999999999942,  0.0000000707827948, -0.0000000805621738],
    [-0.0000000707827974,  0.9999999999999969, -0.0000000330604088],
    [ 0.0000000805621715,  0.0000000330604145,  0.9999999999999962],
];

/// Apply the frame-bias rotation between the J2000 mean-equatorial frame and GCRS/ICRS to
/// `pos` (`[x, y, z, vx, vy, vz]`) in place, per `models.bias` and `direction`. Also applies the
/// JPL Horizons approximation correction when `JPLHOR_APPROX` is set. Port of C `swi_bias`.
pub fn frame_bias(
    pos: &mut [f64; 6],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    direction: FrameTransform,
) {
    if models.bias == BiasModel::None {
        return;
    }

    if flags.contains(CalcFlags::JPLHOR_APPROX) {
        if models.jplhora_mode == JplHoraMode::V2 {
            return;
        }
        if models.jplhora_mode == JplHoraMode::V3 && jd < DPSI_DEPS_IAU1980_TJD0_HORIZONS {
            return;
        }
    }

    let rb = match models.bias {
        BiasModel::IAU2006 => &BIAS_IAU2006,
        _ => &BIAS_IAU2000,
    };

    let has_speed = flags.contains(CalcFlags::SPEED);

    let mut xx = [0.0; 6];
    match direction {
        FrameTransform::J2000ToGcrs => {
            approx_jplhor(pos, jd, flags, models, FrameTransform::J2000ToGcrs);
            for i in 0..3 {
                xx[i] = pos[0] * rb[i][0] + pos[1] * rb[i][1] + pos[2] * rb[i][2];
                if has_speed {
                    xx[i + 3] = pos[3] * rb[i][0] + pos[4] * rb[i][1] + pos[5] * rb[i][2];
                }
            }
        }
        FrameTransform::GcrsToJ2000 => {
            for i in 0..3 {
                xx[i] = pos[0] * rb[0][i] + pos[1] * rb[1][i] + pos[2] * rb[2][i];
                if has_speed {
                    xx[i + 3] = pos[3] * rb[0][i] + pos[4] * rb[1][i] + pos[5] * rb[2][i];
                }
            }
            approx_jplhor(&mut xx, jd, flags, models, FrameTransform::GcrsToJ2000);
        }
    }

    pos[..3].copy_from_slice(&xx[..3]);
    if has_speed {
        pos[3..6].copy_from_slice(&xx[3..6]);
    }
}

fn approx_jplhor(
    x: &mut [f64; 6],
    jd: f64,
    flags: CalcFlags,
    models: &AstroModels,
    direction: FrameTransform,
) {
    if !flags.contains(CalcFlags::JPLHOR_APPROX) {
        return;
    }
    if models.jplhora_mode == JplHoraMode::V2 {
        return;
    }

    let t = (jd - DCOR_RA_JPL_TJD0) / 365.25;
    let dofs = if t < 0.0 {
        DCOR_RA_JPL[0]
    } else if t >= (DCOR_RA_JPL.len() - 1) as f64 {
        DCOR_RA_JPL[DCOR_RA_JPL.len() - 1]
    } else {
        let t0 = t as usize;
        (t - t0 as f64) * (DCOR_RA_JPL[t0] - DCOR_RA_JPL[t0 + 1]) + DCOR_RA_JPL[t0]
    };
    let dofs_rad = dofs / (1000.0 * 3600.0) * DEGTORAD;

    let mut polar = cartesian_to_polar([x[0], x[1], x[2]]);
    match direction {
        FrameTransform::J2000ToGcrs => polar[0] -= dofs_rad,
        FrameTransform::GcrsToJ2000 => polar[0] += dofs_rad,
    }
    let cart = polar_to_cartesian(polar);
    x[..3].copy_from_slice(&cart);
}

// GCRS/ICRS ↔ FK5 rotation (swephlib.c:2292–2333).
// backward=true: FK5 → ICRS (rb · x).
// backward=false: ICRS → FK5 (rb^T · x).
// Writes back all 6 components unconditionally (matching C); since the fixstar
// pipeline forces SPEED internally, has_speed is always true in practice.
#[rustfmt::skip]
const RB: [[f64; 3]; 3] = [
    [0.9999999999999928, 0.0000001110223287, 0.0000000441180557],
    [-0.0000001110223330, 0.9999999999999891, 0.0000000964779176],
    [-0.0000000441180450, -0.0000000964779225, 0.9999999999999943],
];

pub(crate) fn icrs2fk5(x: &mut [f64; 6], has_speed: bool, backward: bool) {
    let mut xx = [0.0f64; 6];
    if backward {
        // FK5 → ICRS: xx = rb · x (row-major)
        for i in 0..3 {
            xx[i] = x[0] * RB[i][0] + x[1] * RB[i][1] + x[2] * RB[i][2];
            if has_speed {
                xx[i + 3] = x[3] * RB[i][0] + x[4] * RB[i][1] + x[5] * RB[i][2];
            }
        }
    } else {
        // ICRS → FK5: xx = rb^T · x (column-major = transpose)
        for i in 0..3 {
            xx[i] = x[0] * RB[0][i] + x[1] * RB[1][i] + x[2] * RB[2][i];
            if has_speed {
                xx[i + 3] = x[3] * RB[0][i] + x[4] * RB[1][i] + x[5] * RB[2][i];
            }
        }
    }
    *x = xx;
}

// FK4 → FK5 RA correction in polar space (swephlib.c:4098–4112).
// Reference: Explanatory Supplement to the Astronomical Almanac, p. 167f.
pub(crate) fn fk4_fk5(xp: &mut [f64; 6], tjd: f64) {
    if xp[0] == 0.0 && xp[1] == 0.0 && xp[2] == 0.0 {
        return;
    }
    let correct_speed = xp[3] != 0.0;
    *xp = cartesian_to_polar_with_speed(*xp);
    // 0.035 is a standalone constant, NOT divided by 36524.2198782
    xp[0] += (0.035 + 0.085 * (tjd - B1950) / 36524.2198782) / 3600.0 * 15.0 * DEGTORAD;
    if correct_speed {
        xp[3] += (0.085 / 36524.2198782) / 3600.0 * 15.0 * DEGTORAD;
    }
    *xp = polar_to_cartesian_with_speed(*xp);
}
