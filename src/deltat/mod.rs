pub mod data;

use crate::constants::*;
use crate::context::EphemerisConfig;
use crate::types::*;
use data::*;

// ---------------------------------------------------------------------------
// Tidal acceleration correction (swephlib.c:3143–3151)
// ---------------------------------------------------------------------------

fn adjust_for_tidacc(
    ans: f64,
    y: f64,
    tid_acc: f64,
    tid_acc0: f64,
    adjust_after_1955: bool,
) -> f64 {
    if y < 1955.0 || adjust_after_1955 {
        let b = y - 1955.0;
        ans + (-0.000091) * (tid_acc - tid_acc0) * b * b
    } else {
        ans
    }
}

// ---------------------------------------------------------------------------
// Tidal acceleration resolver (swephlib.c:3198–3240)
// ---------------------------------------------------------------------------

fn resolve_tidal_acceleration(config: &EphemerisConfig) -> f64 {
    if let Some(ta) = config.tidal_acceleration {
        if ta != TIDAL_AUTOMATIC {
            return ta;
        }
    }
    match config.ephemeris_source {
        EphemerisSource::Moshier => TIDAL_DE404,
        // JPL/Swiss: when file backends land, read denum from file header.
        // Until then, fall back to DE431.
        EphemerisSource::Jpl | EphemerisSource::Swiss => TIDAL_DEFAULT,
    }
}

// ---------------------------------------------------------------------------
// Long-term parabola: Morrison & Stephenson (swephlib.c:2841–2846)
// ---------------------------------------------------------------------------

fn deltat_longterm_morrison_stephenson(tjd: f64) -> f64 {
    let y = 2000.0 + (tjd - J2000) / 365.2425;
    let u = (y - 1820.0) / 100.0;
    -20.0 + 32.0 * u * u
}

// ---------------------------------------------------------------------------
// Stephenson & Morrison 1984 (inline in calc_deltat, swephlib.c:2656–2668)
// ---------------------------------------------------------------------------

fn deltat_stephenson_morrison_1984(tjd: f64) -> f64 {
    let y = 2000.0 + (tjd - J2000) / 365.25;
    let ans = if y >= 948.0 {
        let b = 0.01 * (y - 2000.0);
        (23.58 * b + 100.3) * b + 101.6
    } else {
        let b = 0.01 * (y - 2000.0) + 3.75;
        35.0 * b * b + 40.0
    };
    ans / 86400.0
}

// ---------------------------------------------------------------------------
// Stephenson 1997 (swephlib.c:2848–2887)
// ---------------------------------------------------------------------------

fn deltat_stephenson_morrison_1997(tjd: f64, tid_acc: f64) -> f64 {
    let y = 2000.0 + (tjd - J2000) / 365.25;

    if y < -500.0 {
        let b = (y - 1735.0) * 0.01;
        let mut ans = -20.0 + 35.0 * b * b;
        ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);

        if y >= -600.0 {
            let ans2 = adjust_for_tidacc(DT97[0], -500.0, tid_acc, TIDAL_26, false);
            let b3 = (-500.0 - 1735.0) * 0.01;
            let ans3 = adjust_for_tidacc(-20.0 + 35.0 * b3 * b3, y, tid_acc, TIDAL_26, false);
            let dd = ans3 - ans2;
            let b_blend = (y - (-600.0)) * 0.01;
            ans -= dd * b_blend;
        }
        ans / 86400.0
    } else {
        let iy = ((y.floor() - TAB97_START as f64) / TAB97_STEP as f64) as usize;
        let dd = (y - (TAB97_START as f64 + TAB97_STEP as f64 * iy as f64)) / TAB97_STEP as f64;
        let mut ans = DT97[iy] + (DT97[iy + 1] - DT97[iy]) * dd;
        ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
        ans / 86400.0
    }
}

// ---------------------------------------------------------------------------
// Stephenson & Morrison 2004 (swephlib.c:2890–2933)
// ---------------------------------------------------------------------------

fn deltat_stephenson_morrison_2004(tjd: f64, tid_acc: f64) -> f64 {
    let y = 2000.0 + (tjd - J2000) / 365.2425;

    if y < -1000.0 {
        let mut ans = deltat_longterm_morrison_stephenson(tjd);
        ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);

        if y >= -1100.0 {
            let ans2 = adjust_for_tidacc(DT2[0], -1000.0, tid_acc, TIDAL_26, false);
            let tjd0 = (-1000.0 - 2000.0) * 365.2425 + J2000;
            let ans3 = deltat_longterm_morrison_stephenson(tjd0);
            let ans3 = adjust_for_tidacc(ans3, y, tid_acc, TIDAL_26, false);
            let dd = ans3 - ans2;
            let b = (y - (-1100.0)) * 0.01;
            ans -= dd * b;
        }
        ans / 86400.0
    } else {
        // Mixed-year bug: outer routing uses Gregorian, table indexing uses Julian
        let y_jul = 2000.0 + (tjd - 2451557.5) / 365.25;
        let iy = ((y_jul.floor() - TAB2_START as f64) / TAB2_STEP as f64) as usize;
        let dd = (y_jul - (TAB2_START as f64 + TAB2_STEP as f64 * iy as f64)) / TAB2_STEP as f64;
        let mut ans = DT2[iy] + (DT2[iy + 1] - DT2[iy]) * dd;
        ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
        ans / 86400.0
    }
}

// ---------------------------------------------------------------------------
// Stephenson et al. 2016 (swephlib.c:3001–3036)
// ---------------------------------------------------------------------------

fn deltat_stephenson_etc_2016(tjd: f64, tid_acc: f64) -> f64 {
    let y = 2000.0 + (tjd - J2000) / 365.2425;

    let mut dt = None;
    for row in &DTCF16 {
        if tjd < row[0] {
            break;
        }
        if tjd < row[1] {
            let t = (tjd - row[0]) / (row[1] - row[0]);
            dt = Some(row[2] + row[3] * t + row[4] * t * t + row[5] * t * t * t);
            break;
        }
    }

    let dt = match dt {
        Some(v) => v,
        None if y < -720.0 => {
            let t = (y - 1825.0) / 100.0;
            -320.0 + 32.5 * t * t - 179.7337208
        }
        None => {
            let t = (y - 1825.0) / 100.0;
            -320.0 + 32.5 * t * t + 269.4790417
        }
    };

    let dt = adjust_for_tidacc(dt, y, tid_acc, TIDAL_STEPHENSON_2016, true);
    dt / 86400.0
}

// ---------------------------------------------------------------------------
// Espenak & Meeus 2006 (swephlib.c:3038–3084)
// ---------------------------------------------------------------------------

fn deltat_espenak_meeus_1620(tjd: f64, tid_acc: f64) -> f64 {
    let y = 2000.0 + (tjd - J2000) / 365.2425;

    let ans = if y < -500.0 {
        deltat_longterm_morrison_stephenson(tjd)
    } else if y < 500.0 {
        let u = y / 100.0;
        (((((0.0090316521 * u + 0.022174192) * u - 0.1798452) * u - 5.952053) * u + 33.78311) * u
            - 1014.41)
            * u
            + 10583.6
    } else if y < 1600.0 {
        let u = (y - 1000.0) / 100.0;
        (((((0.0083572073 * u - 0.005050998) * u - 0.8503463) * u + 0.319781) * u + 71.23472) * u
            - 556.01)
            * u
            + 1574.2
    } else if y < 1700.0 {
        let u = y - 1600.0;
        120.0 - 0.9808 * u - 0.01532 * u * u + u * u * u / 7129.0
    } else if y < 1800.0 {
        let u = y - 1700.0;
        (((-u / 1174000.0 + 0.00013336) * u - 0.0059285) * u + 0.1603) * u + 8.83
    } else if y < 1860.0 {
        let u = y - 1800.0;
        ((((((0.000000000875 * u - 0.0000001699) * u + 0.0000121272) * u - 0.00037436) * u
            + 0.0041116)
            * u
            + 0.0068612)
            * u
            - 0.332447)
            * u
            + 13.72
    } else if y < 1900.0 {
        let u = y - 1860.0;
        ((((u / 233174.0 - 0.0004473624) * u + 0.01680668) * u - 0.251754) * u + 0.5737) * u + 7.62
    } else if y < 1920.0 {
        let u = y - 1900.0;
        (((-0.000197 * u + 0.0061966) * u - 0.0598939) * u + 1.494119) * u - 2.79
    } else if y < 1941.0 {
        let u = y - 1920.0;
        21.20 + 0.84493 * u - 0.076100 * u * u + 0.0020936 * u * u * u
    } else if y < 1961.0 {
        let u = y - 1950.0;
        29.07 + 0.407 * u - u * u / 233.0 + u * u * u / 2547.0
    } else if y < 1986.0 {
        let u = y - 1975.0;
        45.45 + 1.067 * u - u * u / 260.0 - u * u * u / 718.0
    } else if y < 2005.0 {
        let u = y - 2000.0;
        ((((0.00002373599 * u + 0.000651814) * u + 0.0017275) * u - 0.060374) * u + 0.3345) * u
            + 63.86
    } else {
        0.0
    };

    let ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
    ans / 86400.0
}

// ---------------------------------------------------------------------------
// Tabulated + Bessel interpolation + future extrapolation (swephlib.c:2733–2839)
// ---------------------------------------------------------------------------

fn deltat_aa(tjd: f64, tid_acc: f64, model: DeltaTModel) -> f64 {
    let tabsiz = TABSIZ;
    let tabend = TABSTART as f64 + tabsiz as f64 - 1.0;
    let y = 2000.0 + (tjd - 2451544.5) / 365.25;

    if y <= tabend {
        bessel_interpolation(y, tid_acc)
    } else {
        future_extrapolation(y, tabend, tabsiz, model)
    }
}

fn bessel_interpolation(y: f64, tid_acc: f64) -> f64 {
    let tabsiz = TABSIZ;
    let p_int = y.floor();
    let iy = (p_int - TABSTART as f64) as usize;
    let p = y - p_int;

    let mut ans = DT[iy];

    // First order
    let k = iy + 1;
    if k >= tabsiz {
        let ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
        return ans / 86400.0;
    }
    ans += p * (DT[k] - DT[iy]);

    // Guard for second differences
    if iy < 1 || iy + 2 >= tabsiz {
        let ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
        return ans / 86400.0;
    }

    // First differences with boundary zero-padding
    let mut d = [0.0_f64; 5];
    let k_start = iy as i32 - 2;
    for i in 0..5 {
        let k = k_start + i as i32;
        if k >= 0 && (k + 1) < tabsiz as i32 {
            d[i] = DT[(k + 1) as usize] - DT[k as usize];
        }
    }

    // Second differences
    for i in 0..4 {
        d[i] = d[i + 1] - d[i];
    }
    let mut b = 0.25 * p * (p - 1.0);
    ans += b * (d[1] + d[2]);

    // Guard for third differences
    if iy + 2 >= tabsiz {
        let ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
        return ans / 86400.0;
    }

    // Third differences
    for i in 0..3 {
        d[i] = d[i + 1] - d[i];
    }
    b *= 2.0 / 3.0;
    ans += (p - 0.5) * b * d[1];

    // Guard for fourth differences
    if iy < 2 || iy + 3 > tabsiz {
        let ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
        return ans / 86400.0;
    }

    // Fourth differences
    for i in 0..2 {
        d[i] = d[i + 1] - d[i];
    }
    let b = 0.125 * b * (p + 1.0) * (p - 2.0);
    ans += b * (d[0] + d[1]);

    let ans = adjust_for_tidacc(ans, y, tid_acc, TIDAL_26, false);
    ans / 86400.0
}

fn future_extrapolation(y: f64, tabend: f64, tabsiz: usize, model: DeltaTModel) -> f64 {
    let (mut ans, ans2) = if model == DeltaTModel::StephensonEtc2016 && y < 2500.0 {
        let b = y - 2000.0;
        let a = b * b * b * 121.0 / 30000000.0 + b * b / 1250.0 + b * 521.0 / 3000.0 + 64.0;
        let b2 = tabend - 2000.0;
        let a2 = b2 * b2 * b2 * 121.0 / 30000000.0 + b2 * b2 / 1250.0 + b2 * 521.0 / 3000.0 + 64.0;
        (a, a2)
    } else if model == DeltaTModel::StephensonEtc2016 {
        let b = 0.01 * (y - 2000.0);
        (b * b * 32.5 + 42.5, 0.0)
    } else {
        let b = 0.01 * (y - 1820.0);
        let a = -20.0 + 31.0 * b * b;
        let b2 = 0.01 * (tabend - 1820.0);
        let a2 = -20.0 + 31.0 * b2 * b2;
        (a, a2)
    };

    // 100-year transition blend
    if y <= tabend + 100.0 && !(model == DeltaTModel::StephensonEtc2016 && y >= 2500.0) {
        let ans3 = DT[tabsiz - 1];
        let dd = ans2 - ans3;
        ans += dd * (y - (tabend + 100.0)) * 0.01;
    }

    ans / 86400.0
}

// ---------------------------------------------------------------------------
// Main dispatcher (swephlib.c:2545–2699)
// ---------------------------------------------------------------------------

pub fn calc_deltat(tjd: f64, config: &EphemerisConfig) -> f64 {
    let model = config.astro_models.delta_t;
    let tid_acc = resolve_tidal_acceleration(config);
    let y = 2000.0 + (tjd - J2000) / 365.25;
    let y_greg = 2000.0 + (tjd - J2000) / 365.2425;

    match model {
        DeltaTModel::StephensonEtc2016 if tjd < 2435108.5 => {
            let mut dt = deltat_stephenson_etc_2016(tjd, tid_acc);
            if tjd >= 2434108.5 {
                dt += (1.0 - (2435108.5 - tjd) / 1000.0) * 0.6610218 / 86400.0;
            }
            dt
        }
        DeltaTModel::EspenakMeeus2006 if tjd < 2317746.13090277789 => {
            deltat_espenak_meeus_1620(tjd, tid_acc)
        }
        DeltaTModel::StephensonMorrison2004 if y < 1620.0 => {
            if y < 1600.0 {
                deltat_stephenson_morrison_2004(tjd, tid_acc)
            } else {
                let dd = (y - 1600.0) / 20.0;
                let ans = DT2[26] + dd * (DT[0] - DT2[26]);
                let ans = adjust_for_tidacc(ans, y_greg, tid_acc, TIDAL_26, false);
                ans / 86400.0
            }
        }
        DeltaTModel::Stephenson1997 if y < 1620.0 => {
            if y < 1600.0 {
                deltat_stephenson_morrison_1997(tjd, tid_acc)
            } else {
                let dd = (y - 1600.0) / 20.0;
                let ans = DT97[42] + dd * (DT[0] - DT97[42]);
                let ans = adjust_for_tidacc(ans, y_greg, tid_acc, TIDAL_26, false);
                ans / 86400.0
            }
        }
        DeltaTModel::StephensonMorrison1984 if y < 1620.0 => deltat_stephenson_morrison_1984(tjd),
        _ if y >= 1620.0 => deltat_aa(tjd, tid_acc, model),
        _ => 0.0 / 86400.0,
    }
}
