use crate::error::Error;

use super::JplFile;

/// Maximum Chebyshev coefficients per component this implementation supports;
/// sizes the fixed `pc[]`/`vc[]` buffers in [`interp`]. DE441's max ncf is 14.
/// [`super::header::parse_header`] rejects files whose ncf exceeds this.
pub(super) const MAX_NCF: usize = 18;

/// Read record `nr` from the mmap'd file into a Vec<f64>.
/// buf[0] = segment start JD, buf[1] = segment end JD, buf[2..] = coefficients.
fn read_record(file: &JplFile, nr: usize) -> Vec<f64> {
    let h = file.header();
    let order = file.byte_order();
    let bytes = file.bytes();
    let offset = nr * h.irecsz;
    let ncoeffs = h.ncoeffs;
    let mut buf = Vec::with_capacity(ncoeffs);
    for k in 0..ncoeffs {
        let bo = offset + k * 8;
        let b: [u8; 8] = bytes[bo..bo + 8].try_into().unwrap();
        buf.push(order.read_f64(b));
    }
    buf
}

/// JPL Chebyshev interpolation. (swejpl.c:472–591)
///
/// `coeffs` — slice of buf starting at body's first coefficient (0-based, after ipt offset).
/// `t`      — normalized segment time ∈ [0, 1).
/// `intv`   — segment length in days (ss[2]).
/// `ncf`    — Chebyshev coefficients per component per sub-interval.
/// `ncm`    — components (3 for planets/Moon/Sun, 2 for nutations).
/// `na`     — sub-intervals per segment.
///
/// Returns [f64; 6]: positions in [0..ncm], velocities in [ncm..2*ncm].
/// Velocity components are only meaningful when `need_speed` is true.
fn interp(
    coeffs: &[f64],
    t: f64,
    intv: f64,
    ncf: usize,
    ncm: usize,
    na: usize,
    need_speed: bool,
) -> [f64; 6] {
    // Sub-interval selection (swejpl.c:497-504)
    let dt1 = t.floor();
    let temp = na as f64 * t;
    let ni = (temp - dt1) as usize;
    let tc = (temp.rem_euclid(1.0) + dt1) * 2.0 - 1.0;

    // Position Chebyshev polynomials T_0..T_{ncf-1} (swejpl.c:511-533)
    let mut pc = [0.0f64; MAX_NCF];
    pc[0] = 1.0;
    pc[1] = tc;
    let twot = tc * 2.0;
    for i in 2..ncf {
        pc[i] = twot * pc[i - 1] - pc[i - 2];
    }

    let mut pv = [0.0f64; 6];

    // Position: descending dot product per component
    for (c, pv_c) in pv.iter_mut().enumerate().take(ncm) {
        let base = (c + ni * ncm) * ncf;
        let mut sum = 0.0f64;
        let mut j = ncf;
        while j > 0 {
            j -= 1;
            sum += pc[j] * coeffs[base + j];
        }
        *pv_c = sum;
    }

    if need_speed {
        // Velocity derivative polynomials (swejpl.c:540-553)
        let mut vc = [0.0f64; MAX_NCF];
        vc[0] = 0.0;
        vc[1] = 1.0;
        if ncf > 2 {
            vc[2] = twot + twot; // 4*tc
        }
        for i in 3..ncf {
            vc[i] = twot * vc[i - 1] + 2.0 * pc[i - 1] - vc[i - 2];
        }
        let bma = 2.0 * na as f64 / intv;
        for c in 0..ncm {
            let base = (c + ni * ncm) * ncf;
            let mut sum = 0.0f64;
            // Descending from ncf-1 to 1 (j=0 term is always zero since vc[0]=0)
            let mut j = ncf;
            while j > 1 {
                j -= 1;
                sum += vc[j] * coeffs[base + j];
            }
            pv[ncm + c] = sum * bma;
        }
    }

    pv
}

/// Record selection, read, and body interpolation. (swejpl.c:783–851)
///
/// `list[i] > 0` means interpolate body i (bodies 0..9: Mercury..Moon).
/// Returns `(pv[0..13], pvsun)` where `pv[J_EARTH]` = raw EMB (before Earth/Moon
/// decomposition), `pv[J_MOON]` = geocentric Moon.
pub(super) fn state(
    file: &JplFile,
    et: f64,
    list: &[u8; 12],
    do_bary: bool,
    need_speed: bool,
) -> Result<([[f64; 6]; 13], [f64; 6]), Error> {
    let h = file.header();
    let ss = h.ss;
    let ipt = h.ipt;

    if et < ss[0] || et > ss[1] {
        return Err(Error::BeyondEphemerisLimits {
            jd_tt: et,
            start: ss[0],
            end: ss[1],
        });
    }

    // Epoch decomposition (swejpl.c:783–797). Match C's accumulation order:
    // et_fr = s - floor(s) (NOT et - et_mn), so the fractional day is formed from
    // the same operands C uses.
    let s = et - 0.5;
    let et_mn_floor = s.floor();
    let et_fr = s - et_mn_floor;
    let et_mn = et_mn_floor + 0.5;

    // Record number (+2: records 0 and 1 are header and constants)
    let mut nr = ((et_mn - ss[0]) / ss[2]) as i32 + 2;
    if et_mn == ss[1] {
        nr -= 1;
    }

    // Normalized time within segment ∈ [0, 1)
    let t = (et_mn - ((nr - 2) as f64 * ss[2] + ss[0]) + et_fr) / ss[2];

    let intv = ss[2]; // do_km always false: AU + AU/day
    let buf = read_record(file, nr as usize);
    let aufac = 1.0 / h.au;

    // Always interpolate SSBary Sun (ipt[30..32]), ifl=need_speed (swejpl.c:486)
    let sun_off = (ipt[30] - 1) as usize;
    let sun_ncf = ipt[31] as usize;
    let sun_na = ipt[32] as usize;
    let sun_raw = interp(&buf[sun_off..], t, intv, sun_ncf, 3, sun_na, need_speed);
    let mut pvsun = [0.0f64; 6];
    for k in 0..6 {
        pvsun[k] = sun_raw[k] * aufac;
    }

    // Bodies 0..9: Mercury..Moon (swejpl.c:489–500)
    let mut pv = [[0.0f64; 6]; 13];
    for i in 0..10usize {
        if list[i] == 0 {
            continue;
        }
        let b_off = (ipt[i * 3] - 1) as usize;
        let b_ncf = ipt[i * 3 + 1] as usize;
        let b_na = ipt[i * 3 + 2] as usize;
        let raw = interp(&buf[b_off..], t, intv, b_ncf, 3, b_na, need_speed);
        for k in 0..6 {
            pv[i][k] = if i < 9 && !do_bary {
                raw[k] * aufac - pvsun[k]
            } else {
                raw[k] * aufac
            };
        }
    }

    Ok((pv, pvsun))
}
