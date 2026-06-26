use crate::constants::{J2000, STR};
use crate::math::mods3600;

use super::PlantTbl;

const TIMESCALE: f64 = 3652500.0;

const FREQS: [f64; 9] = [
    53810162868.8982,
    21066413643.3548,
    12959774228.3429,
    6890507749.3988,
    1092566037.7991,
    439960985.5372,
    154248119.3933,
    78655032.0744,
    52272245.1795,
];

const PHASES: [f64; 9] = [
    252.25090552 * 3600.0,
    181.97980085 * 3600.0,
    100.46645683 * 3600.0,
    355.43299958 * 3600.0,
    34.35151874 * 3600.0,
    50.07744430 * 3600.0,
    314.05500511 * 3600.0,
    304.34866548 * 3600.0,
    860492.1546,
];

fn sscc(k: usize, arg: f64, n: i8, ss: &mut [[f64; 24]; 9], cc: &mut [[f64; 24]; 9]) {
    let n = n as usize;
    let su = arg.sin();
    let cu = arg.cos();
    ss[k][0] = su;
    cc[k][0] = cu;
    let mut sv = 2.0 * su * cu;
    let mut cv = cu * cu - su * su;
    ss[k][1] = sv;
    cc[k][1] = cv;
    for i in 2..n {
        let s = su * cv + cu * sv;
        cv = cu * cv - su * sv;
        sv = s;
        ss[k][i] = sv;
        cc[k][i] = cv;
    }
}

pub fn moshplan2(jd: f64, table: &PlantTbl) -> [f64; 3] {
    let t = (jd - J2000) / TIMESCALE;

    let mut ss = [[0.0f64; 24]; 9];
    let mut cc = [[0.0f64; 24]; 9];
    for i in 0..9 {
        if table.max_harmonic[i] > 0 {
            let sr = (mods3600(FREQS[i] * t) + PHASES[i]) * STR;
            sscc(i, sr, table.max_harmonic[i], &mut ss, &mut cc);
        }
    }

    let arg = table.arg_tbl;
    let mut pl = 0usize;
    let mut pb = 0usize;
    let mut pr = 0usize;
    let mut p = 0usize;
    let mut sl = 0.0f64;
    let mut sb = 0.0f64;
    let mut sr = 0.0f64;

    loop {
        let np = arg[p];
        p += 1;
        if np < 0 {
            break;
        }

        if np == 0 {
            let nt = arg[p] as usize;
            p += 1;

            let mut cu = table.lon_tbl[pl];
            pl += 1;
            for _ in 0..nt {
                cu = cu * t + table.lon_tbl[pl];
                pl += 1;
            }
            sl += mods3600(cu);

            let mut cu = table.lat_tbl[pb];
            pb += 1;
            for _ in 0..nt {
                cu = cu * t + table.lat_tbl[pb];
                pb += 1;
            }
            sb += cu;

            let mut cu = table.rad_tbl[pr];
            pr += 1;
            for _ in 0..nt {
                cu = cu * t + table.rad_tbl[pr];
                pr += 1;
            }
            sr += cu;
            continue;
        }

        let np = np as usize;
        let mut sv = 0.0f64;
        let mut cv = 0.0f64;
        let mut k1 = false;
        for _ in 0..np {
            let j = arg[p] as i32;
            p += 1;
            let m = arg[p] as usize;
            p += 1;
            if j == 0 {
                continue;
            }
            let k = (j.unsigned_abs() - 1) as usize;
            let su = if j < 0 { -ss[m - 1][k] } else { ss[m - 1][k] };
            let cu = cc[m - 1][k];
            if !k1 {
                sv = su;
                cv = cu;
                k1 = true;
            } else {
                let tmp = su * cv + cu * sv;
                cv = cu * cv - su * sv;
                sv = tmp;
            }
        }

        let nt = arg[p] as usize;
        p += 1;

        let mut cu = table.lon_tbl[pl];
        pl += 1;
        let mut su = table.lon_tbl[pl];
        pl += 1;
        for _ in 0..nt {
            cu = cu * t + table.lon_tbl[pl];
            pl += 1;
            su = su * t + table.lon_tbl[pl];
            pl += 1;
        }
        sl += cu * cv + su * sv;

        let mut cu = table.lat_tbl[pb];
        pb += 1;
        let mut su = table.lat_tbl[pb];
        pb += 1;
        for _ in 0..nt {
            cu = cu * t + table.lat_tbl[pb];
            pb += 1;
            su = su * t + table.lat_tbl[pb];
            pb += 1;
        }
        sb += cu * cv + su * sv;

        let mut cu = table.rad_tbl[pr];
        pr += 1;
        let mut su = table.rad_tbl[pr];
        pr += 1;
        for _ in 0..nt {
            cu = cu * t + table.rad_tbl[pr];
            pr += 1;
            su = su * t + table.rad_tbl[pr];
            pr += 1;
        }
        sr += cu * cv + su * sv;
    }

    [
        STR * sl,
        STR * sb,
        STR * table.distance * sr + table.distance,
    ]
}
