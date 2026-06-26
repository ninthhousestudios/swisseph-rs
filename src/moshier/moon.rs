use crate::constants::{AUNIT, J2000, STR};
use crate::math::mods3600;

use super::moon_tables::*;

pub struct MeanElements {
    pub m: f64,
    pub nf: f64,
    pub mp: f64,
    pub d: f64,
    pub swelp: f64,
}

struct PlanetaryElements {
    ve: f64,
    ea: f64,
    ma: f64,
    ju: f64,
    sa: f64,
}

struct MoonState {
    t: f64,
    ve: f64,
    ea: f64,
    ma: f64,
    ju: f64,
    sa: f64,
    l: f64,
    l1: f64,
    l2: f64,
    l3: f64,
    l4: f64,
    b: f64,
    f_arg: f64,
    moonpol: [f64; 3],
    ss: [[f64; 8]; 5],
    cc: [[f64; 8]; 5],
}

pub fn mean_elements(t: f64) -> MeanElements {
    let t2 = t * t;
    let frac_t = t % 1.0;

    let m = mods3600(129600000.0 * frac_t - 3418.961646 * t + 1287104.76154)
        + ((((((((1.62e-20 * t - 1.0390e-17) * t - 3.83508e-15) * t + 4.237343e-13) * t
            + 8.8555011e-11)
            * t
            - 4.77258489e-8)
            * t
            - 1.1297037031e-5)
            * t
            + 1.4732069041e-4)
            * t
            - 0.552891801772)
            * t2;

    let nf =
        mods3600(1739232000.0 * frac_t + 295263.0983 * t - 2.079419901760e-01 * t + 335779.55755)
            + ((Z[2] * t + Z[1]) * t + Z[0]) * t2;

    let mp =
        mods3600(1717200000.0 * frac_t + 715923.4728 * t - 2.035946368532e-01 * t + 485868.28096)
            + ((Z[5] * t + Z[4]) * t + Z[3]) * t2;

    let d =
        mods3600(1601856000.0 * frac_t + 1105601.4603 * t + 3.962893294503e-01 * t + 1072260.73512)
            + ((Z[8] * t + Z[7]) * t + Z[6]) * t2;

    let swelp =
        mods3600(1731456000.0 * frac_t + 1108372.83264 * t - 6.784914260953e-01 * t + 785939.95571)
            + ((Z[11] * t + Z[10]) * t + Z[9]) * t2;

    MeanElements {
        m,
        nf,
        mp,
        d,
        swelp,
    }
}

fn mean_elements_pl(t: f64, t2: f64) -> PlanetaryElements {
    let ve = mods3600(210664136.4335482 * t + 655127.283046)
        + ((((((((-9.36e-023 * t - 1.95e-20) * t + 6.097e-18) * t + 4.43201e-15) * t
            + 2.509418e-13)
            * t
            - 3.0622898e-10)
            * t
            - 2.26602516e-9)
            * t
            - 1.4244812531e-5)
            * t
            + 0.005871373088)
            * t2;

    let ea = mods3600(129597742.26669231 * t + 361679.214649)
        + ((((((((-1.16e-22 * t + 2.976e-19) * t + 2.8460e-17) * t - 1.08402e-14) * t
            - 1.226182e-12)
            * t
            + 1.7228268e-10)
            * t
            + 1.515912254e-7)
            * t
            + 8.863982531e-6)
            * t
            - 2.0199859001e-2)
            * t2;

    let ma = mods3600(68905077.59284 * t + 1279559.78866) + (-1.043e-5 * t + 9.38012e-3) * t2;

    let ju =
        mods3600(10925660.428608 * t + 123665.342120) + (1.543273e-5 * t - 3.06037836351e-1) * t2;

    let sa = mods3600(4399609.65932 * t + 180278.89694)
        + ((4.475946e-8 * t - 6.874806e-5) * t + 7.56161437443e-1) * t2;

    PlanetaryElements { ve, ea, ma, ju, sa }
}

fn sscc(k: usize, arg: f64, n: usize, ss: &mut [[f64; 8]; 5], cc: &mut [[f64; 8]; 5]) {
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

fn chewm(
    table: &[i16],
    nlines: usize,
    nangles: usize,
    typflg: u8,
    ans: &mut [f64; 3],
    ss: &[[f64; 8]; 5],
    cc: &[[f64; 8]; 5],
) {
    let mut p = 0;
    for _ in 0..nlines {
        let mut sv = 0.0f64;
        let mut cv = 0.0f64;
        let mut k1 = false;
        for m in 0..nangles {
            let j = table[p] as i32;
            p += 1;
            if j != 0 {
                let k = (j.unsigned_abs() - 1) as usize;
                let su = if j < 0 { -ss[m][k] } else { ss[m][k] };
                let cu = cc[m][k];
                if !k1 {
                    sv = su;
                    cv = cu;
                    k1 = true;
                } else {
                    let ff = su * cv + cu * sv;
                    cv = cu * cv - su * sv;
                    sv = ff;
                }
            }
        }
        match typflg {
            1 => {
                let j = table[p] as f64;
                p += 1;
                let k = table[p] as f64;
                p += 1;
                ans[0] += (10000.0 * j + k) * sv;
                let j = table[p] as f64;
                p += 1;
                let k = table[p] as f64;
                p += 1;
                if k != 0.0 {
                    ans[2] += (10000.0 * j + k) * cv;
                }
            }
            2 => {
                let j = table[p] as f64;
                p += 1;
                let k = table[p] as f64;
                p += 1;
                ans[0] += j * sv;
                ans[2] += k * cv;
            }
            3 => {
                let j = table[p] as f64;
                p += 1;
                let k = table[p] as f64;
                p += 1;
                ans[1] += (10000.0 * j + k) * sv;
            }
            4 => {
                let j = table[p] as f64;
                p += 1;
                ans[1] += j * sv;
            }
            _ => unreachable!(),
        }
    }
}

fn moon1(s: &mut MoonState, me: &MeanElements) {
    for i in 0..5 {
        for j in 0..8 {
            s.ss[i][j] = 0.0;
            s.cc[i][j] = 0.0;
        }
    }

    sscc(0, STR * me.d, 6, &mut s.ss, &mut s.cc);
    sscc(1, STR * me.m, 4, &mut s.ss, &mut s.cc);
    sscc(2, STR * me.mp, 4, &mut s.ss, &mut s.cc);
    sscc(3, STR * me.nf, 4, &mut s.ss, &mut s.cc);

    s.moonpol[0] = 0.0;
    s.moonpol[1] = 0.0;
    s.moonpol[2] = 0.0;

    // Phase A: T² series
    chewm(&LRT2, NLRT2, 4, 2, &mut s.moonpol, &s.ss, &s.cc);
    chewm(&BT2, NBT2, 4, 4, &mut s.moonpol, &s.ss, &s.cc);

    // Planetary perturbation terms — Phase A
    let f = 18.0 * s.ve - 16.0 * s.ea;
    s.f_arg = f;

    // Term 1: 18V - 16E - l
    let g = STR * (f - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l = 6.367278 * cg + 12.747036 * sg;
    s.l1 = 23123.70 * cg - 10570.02 * sg;
    s.l2 = Z[12] * cg + Z[13] * sg;
    s.moonpol[2] += 5.01 * cg + 2.72 * sg;

    // Term 2: 10V - 3E - l
    let g = STR * (10.0 * s.ve - 3.0 * s.ea - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.253102 * cg + 0.503359 * sg;
    s.l1 += 1258.46 * cg + 707.29 * sg;
    s.l2 += Z[14] * cg + Z[15] * sg;

    // Term 3: 8V - 13E
    let g = STR * (8.0 * s.ve - 13.0 * s.ea);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.187231 * cg - 0.127481 * sg;
    s.l1 += -319.87 * cg - 18.34 * sg;
    s.l2 += Z[16] * cg + Z[17] * sg;

    let a = 4.0 * s.ea - 8.0 * s.ma + 3.0 * s.ju;

    // Term 4: 4E - 8M + 3J
    let g = STR * a;
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.866287 * cg + 0.248192 * sg;
    s.l1 += 41.87 * cg + 1053.97 * sg;
    s.l2 += Z[18] * cg + Z[19] * sg;

    // Term 5: 4E - 8M + 3J - l
    let g = STR * (a - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.165009 * cg + 0.044176 * sg;
    s.l1 += 4.67 * cg + 201.55 * sg;

    // Term 6: 18V - 16E
    let g = STR * f;
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.330401 * cg + 0.661362 * sg;
    s.l1 += 1202.67 * cg - 555.59 * sg;
    s.l2 += Z[20] * cg + Z[21] * sg;

    // Term 7: 18V - 16E - 2l
    let g = STR * (f - 2.0 * me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.352185 * cg + 0.705041 * sg;
    s.l1 += 1283.59 * cg - 586.43 * sg;

    // Term 8: 2J - 5S
    let g = STR * (2.0 * s.ju - 5.0 * s.sa);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.034700 * cg + 0.160041 * sg;
    s.l2 += Z[22] * cg + Z[23] * sg;

    // Term 9: L - F
    let g = STR * (me.swelp - me.nf);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.000116 * cg + 7.063040 * sg;
    s.l1 += 298.8 * sg;

    // T³ terms
    let sg = (STR * me.m).sin();
    s.l3 = Z[24] * sg;
    s.l4 = 0.0;

    // Radius corrections (T³ relative to T² base)
    let g = STR * (2.0 * me.d - me.m);
    s.moonpol[2] += -0.2655 * g.cos() * s.t;

    let g = STR * (me.m - me.mp);
    s.moonpol[2] += -0.1568 * g.cos() * s.t;

    let g = STR * (me.m + me.mp);
    s.moonpol[2] += 0.1309 * g.cos() * s.t;

    let g = STR * (2.0 * (me.d + me.m) - me.mp);
    s.moonpol[2] += 0.5568 * g.cos() * s.t;

    s.l2 += s.moonpol[0];

    let g = STR * (2.0 * me.d - me.m - me.mp);
    s.moonpol[2] += -0.1910 * g.cos() * s.t;

    s.moonpol[1] *= s.t;
    s.moonpol[2] *= s.t;

    // Phase B: T¹ series
    s.moonpol[0] = 0.0;

    chewm(&BT, NBT, 4, 4, &mut s.moonpol, &s.ss, &s.cc);
    chewm(&LRT, NLRT, 4, 1, &mut s.moonpol, &s.ss, &s.cc);

    // Latitude perturbations — Phase B
    let g = STR * (f - me.mp - me.nf - 2355767.6);
    s.moonpol[1] += -1127.0 * g.sin();

    let g = STR * (f - me.mp + me.nf - 235353.6);
    s.moonpol[1] += -1123.0 * g.sin();

    let g = STR * (s.ea + me.d + 51987.6);
    s.moonpol[1] += 1303.0 * g.sin();

    let g = STR * me.swelp;
    s.moonpol[1] += 342.0 * g.sin();

    // Longitude+speed perturbations — Phase B
    // Term 10: 2V - 3E
    let g = STR * (2.0 * s.ve - 3.0 * s.ea);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.343550 * cg - 0.000276 * sg;
    s.l1 += 105.90 * cg + 336.53 * sg;

    // Term 11: 18V - 16E - 2D
    let g = STR * (f - 2.0 * me.d);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.074668 * cg + 0.149501 * sg;
    s.l1 += 271.77 * cg - 124.20 * sg;

    // Term 12: 18V - 16E - 2D - l
    let g = STR * (f - 2.0 * me.d - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.073444 * cg + 0.147094 * sg;
    s.l1 += 265.24 * cg - 121.16 * sg;

    // Term 13: 18V - 16E + 2D - l
    let g = STR * (f + 2.0 * me.d - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.072844 * cg + 0.145829 * sg;
    s.l1 += 265.18 * cg - 121.29 * sg;

    // Term 14: 18V - 16E + 2(D - l)
    let g = STR * (f + 2.0 * (me.d - me.mp));
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.070201 * cg + 0.140542 * sg;
    s.l1 += 255.36 * cg - 116.79 * sg;

    // Term 15: E + D - F
    let g = STR * (s.ea + me.d - me.nf);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.288209 * cg - 0.025901 * sg;
    s.l1 += -63.51 * cg - 240.14 * sg;

    // Term 16: 2E - 3J + 2D - l
    let g = STR * (2.0 * s.ea - 3.0 * s.ju + 2.0 * me.d - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += 0.077865 * cg + 0.438460 * sg;
    s.l1 += 210.57 * cg + 124.84 * sg;

    // Term 17: E - 2Ma
    let g = STR * (s.ea - 2.0 * s.ma);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.216579 * cg + 0.241702 * sg;
    s.l1 += 197.67 * cg + 125.23 * sg;

    // Term 18: 4E - 8Ma + 3J + l
    let g = STR * (a + me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.165009 * cg + 0.044176 * sg;
    s.l1 += 4.67 * cg + 201.55 * sg;

    // Term 19: 4E - 8Ma + 3J + 2D - l
    let g = STR * (a + 2.0 * me.d - me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.133533 * cg + 0.041116 * sg;
    s.l1 += 6.95 * cg + 187.07 * sg;

    // Term 20: 4E - 8Ma + 3J - 2D + l
    let g = STR * (a - 2.0 * me.d + me.mp);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.133430 * cg + 0.041079 * sg;
    s.l1 += 6.28 * cg + 169.08 * sg;

    // Term 21: 3V - 4E
    let g = STR * (3.0 * s.ve - 4.0 * s.ea);
    let cg = g.cos();
    let sg = g.sin();
    s.l += -0.175074 * cg + 0.003035 * sg;
    s.l1 += 49.17 * cg + 150.57 * sg;

    // Term 22: 2(E + D - l) - 3J + 213534"
    let g = STR * (2.0 * (s.ea + me.d - me.mp) - 3.0 * s.ju + 213534.0);
    s.l1 += 158.4 * g.sin();

    s.l1 += s.moonpol[0];

    let aa = 0.1 * s.t;
    s.moonpol[1] *= aa;
    s.moonpol[2] *= aa;
}

fn moon2(s: &mut MoonState, me: &MeanElements) {
    let f = s.f_arg;

    let g = STR * (2.0 * (s.ea - s.ju + me.d) - me.mp + 648431.172);
    s.l += 1.14307 * g.sin();

    let g = STR * (s.ve - s.ea + 648035.568);
    s.l += 0.82155 * g.sin();

    let g = STR * (3.0 * (s.ve - s.ea) + 2.0 * me.d - me.mp + 647933.184);
    s.l += 0.64371 * g.sin();

    let g = STR * (s.ea - s.ju + 4424.04);
    s.l += 0.63880 * g.sin();

    let g = STR * (me.swelp + me.mp - me.nf + 4.68);
    s.l += 0.49331 * g.sin();

    let g = STR * (me.swelp - me.mp - me.nf + 4.68);
    s.l += 0.4914 * g.sin();

    let g = STR * (me.swelp + me.nf + 2.52);
    s.l += 0.36061 * g.sin();

    let g = STR * (2.0 * s.ve - 2.0 * s.ea + 736.2);
    s.l += 0.30154 * g.sin();

    let g = STR * (2.0 * s.ea - 3.0 * s.ju + 2.0 * me.d - 2.0 * me.mp + 36138.2);
    s.l += 0.28282 * g.sin();

    let g = STR * (2.0 * s.ea - 2.0 * s.ju + 2.0 * me.d - 2.0 * me.mp + 311.0);
    s.l += 0.24516 * g.sin();

    let g = STR * (s.ea - s.ju - 2.0 * me.d + me.mp + 6275.88);
    s.l += 0.21117 * g.sin();

    let g = STR * (2.0 * (s.ea - s.ma) - 846.36);
    s.l += 0.19444 * g.sin();

    let g = STR * (2.0 * (s.ea - s.ju) + 1569.96);
    s.l -= 0.18457 * g.sin();

    let g = STR * (2.0 * (s.ea - s.ju) - me.mp - 55.8);
    s.l += 0.18256 * g.sin();

    let g = STR * (s.ea - s.ju - 2.0 * me.d + 6490.08);
    s.l += 0.16499 * g.sin();

    let g = STR * (s.ea - 2.0 * s.ju - 212378.4);
    s.l += 0.16427 * g.sin();

    let g = STR * (2.0 * (s.ve - s.ea - me.d) + me.mp + 1122.48);
    s.l += 0.16088 * g.sin();

    let g = STR * (s.ve - s.ea - me.mp + 32.04);
    s.l -= 0.15350 * g.sin();

    let g = STR * (s.ea - s.ju - me.mp + 4488.88);
    s.l += 0.14346 * g.sin();

    let g = STR * (2.0 * (s.ve - s.ea + me.d) - me.mp - 8.64);
    s.l += 0.13594 * g.sin();

    let g = STR * (2.0 * (s.ve - s.ea - me.d) + 1319.76);
    s.l += 0.13432 * g.sin();

    let g = STR * (s.ve - s.ea - 2.0 * me.d + me.mp - 56.16);
    s.l -= 0.13122 * g.sin();

    let g = STR * (s.ve - s.ea + me.mp + 54.36);
    s.l -= 0.12722 * g.sin();

    let g = STR * (3.0 * (s.ve - s.ea) - me.mp + 433.8);
    s.l += 0.12539 * g.sin();

    let g = STR * (s.ea - s.ju + me.mp + 4002.12);
    s.l += 0.10994 * g.sin();

    let g = STR * (20.0 * s.ve - 21.0 * s.ea - 2.0 * me.d + me.mp - 317511.72);
    s.l += 0.10652 * g.sin();

    let g = STR * (26.0 * s.ve - 29.0 * s.ea - me.mp + 270002.52);
    s.l += 0.10490 * g.sin();

    let g = STR * (3.0 * s.ve - 4.0 * s.ea + me.d - me.mp - 322765.56);
    s.l += 0.10386 * g.sin();

    // Latitude terms
    let g = STR * (me.swelp + 648002.556);
    s.b = 8.04508 * g.sin();

    let g = STR * (s.ea + me.d + 996048.252);
    s.b += 1.51021 * g.sin();

    let g = STR * (f - me.mp + me.nf + 95554.332);
    s.b += 0.63037 * g.sin();

    let g = STR * (f - me.mp - me.nf + 95553.792);
    s.b += 0.63014 * g.sin();

    let g = STR * (me.swelp - me.mp + 2.9);
    s.b += 0.45587 * g.sin();

    let g = STR * (me.swelp + me.mp + 2.5);
    s.b += -0.41573 * g.sin();

    let g = STR * (me.swelp - 2.0 * me.nf + 3.2);
    s.b += 0.32623 * g.sin();

    let g = STR * (me.swelp - 2.0 * me.d + 2.5);
    s.b += 0.29855 * g.sin();
}

fn moon3(s: &mut MoonState, me: &MeanElements) {
    s.moonpol[0] = 0.0;
    chewm(&LR, NLR, 4, 1, &mut s.moonpol, &s.ss, &s.cc);
    chewm(&MB, NMB, 4, 3, &mut s.moonpol, &s.ss, &s.cc);

    s.l += (((s.l4 * s.t + s.l3) * s.t + s.l2) * s.t + s.l1) * s.t * 1.0e-5;

    s.moonpol[0] = me.swelp + s.l + 1.0e-4 * s.moonpol[0];
    s.moonpol[1] = 1.0e-4 * s.moonpol[1] + s.b;
    s.moonpol[2] = 1.0e-4 * s.moonpol[2] + 385000.52899;
}

fn moon4(s: &mut MoonState) {
    s.moonpol[2] /= AUNIT / 1000.0;
    s.moonpol[0] = STR * mods3600(s.moonpol[0]);
    s.moonpol[1] = STR * s.moonpol[1];
}

pub fn moshmoon2(jd: f64) -> [f64; 3] {
    let t = (jd - J2000) / 36525.0;
    let t2 = t * t;
    let me = mean_elements(t);
    let pl = mean_elements_pl(t, t2);
    let mut s = MoonState {
        t,
        ve: pl.ve,
        ea: pl.ea,
        ma: pl.ma,
        ju: pl.ju,
        sa: pl.sa,
        l: 0.0,
        l1: 0.0,
        l2: 0.0,
        l3: 0.0,
        l4: 0.0,
        b: 0.0,
        f_arg: 0.0,
        moonpol: [0.0; 3],
        ss: [[0.0; 8]; 5],
        cc: [[0.0; 8]; 5],
    };
    moon1(&mut s, &me);
    moon2(&mut s, &me);
    moon3(&mut s, &me);
    moon4(&mut s);
    s.moonpol
}
