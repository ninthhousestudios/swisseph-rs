use swisseph::Ephemeris;
use swisseph::flags::CalcFlags;
use swisseph::math::normalize_degrees;
use swisseph::types::{Body, CalendarType};

use crate::args::SweTestArgs;

const ZODIAC_NAMES: [&str; 12] = [
    "ar", "ta", "ge", "cn", "le", "vi", "li", "sc", "sa", "cp", "aq", "pi",
];

const DEGREE_SIGN: &str = "\u{b0}";

const AUNIT_TO_LIGHTYEAR: f64 = 1.0 / 63241.07708427;
const AUNIT_TO_KM: f64 = 1.495978707e8;

bitflags::bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct DmsFlags: u32 {
        const ROUND_SEC  = 1;
        const ROUND_MIN  = 2;
        const ZODIAC     = 4;
        const LZEROES    = 8;
        const EQUATORIAL = 0x0800; // SEFLG_EQUATORIAL bit value
        const ALLOW_361  = 64;
    }
}

pub struct FormatContext<'a> {
    pub name: String,
    pub ipl: i32,
    pub body: Option<Body>,
    pub jd: f64,
    pub tjd_ut: f64,
    pub tjd_tt: f64,
    pub year: i32,
    pub month: i32,
    pub day: i32,
    pub hour: f64,
    pub cal: CalendarType,
    pub is_ut: bool,
    pub data: [f64; 6],
    pub xequ: Option<[f64; 6]>,
    pub xaz: Option<[f64; 3]>,
    pub xcart: Option<[f64; 6]>,
    pub xecart: Option<[f64; 6]>,
    pub hpos: Option<f64>,
    pub hposj: Option<f64>,
    pub armc: Option<f64>,
    pub attr: Option<[f64; 6]>,
    pub args: &'a SweTestArgs,
    pub is_label: bool,
    pub is_house: bool,
    pub is_ayanamsa: bool,
    pub is_first: bool,
}

pub struct FormatNeeds {
    pub equatorial: bool,
    pub azalt: bool,
    pub zenith: bool,
    pub ecl_cartesian: bool,
    pub equ_cartesian: bool,
    pub house_pos: bool,
    pub phenomena: bool,
    pub speed: bool,
}

pub fn scan_format_needs(fmt: &str) -> FormatNeeds {
    let mut needs = FormatNeeds {
        equatorial: false,
        azalt: false,
        zenith: false,
        ecl_cartesian: false,
        equ_cartesian: false,
        house_pos: false,
        phenomena: false,
        speed: false,
    };
    let has_double_s = fmt.contains("SS") || fmt.contains("ss");
    for ch in fmt.chars() {
        match ch {
            'A' | 'a' | 'D' | 'd' | 'Q' | 'm' => needs.equatorial = true,
            'z' => {
                needs.equatorial = true;
                needs.zenith = true;
            }
            'I' | 'i' | 'H' | 'h' | 'K' | 'k' => needs.azalt = true,
            'U' | 'X' => needs.ecl_cartesian = true,
            'u' | 'x' => needs.equ_cartesian = true,
            'G' | 'g' | 'j' => needs.house_pos = true,
            '+' | '-' | '*' | '/' | '=' => needs.phenomena = true,
            'S' | 's' => needs.speed = true,
            _ => {}
        }
    }
    if has_double_s
        && (fmt.contains('A') || fmt.contains('a') || fmt.contains('D') || fmt.contains('d'))
    {
        needs.equatorial = true;
    }
    if has_double_s && (fmt.contains('U') || fmt.contains('X')) {
        needs.ecl_cartesian = true;
    }
    if has_double_s && (fmt.contains('u') || fmt.contains('x')) {
        needs.equ_cartesian = true;
    }
    needs
}

fn dms(xv: f64, flags: DmsFlags, extra_prec: bool) -> String {
    if xv.is_nan() {
        return "nan".into();
    }
    let mut xv = xv;
    if xv >= 360.0 && !flags.contains(DmsFlags::ALLOW_361) {
        xv = 0.0;
    }
    let c = if flags.contains(DmsFlags::EQUATORIAL) {
        "h"
    } else {
        DEGREE_SIGN
    };
    let sgn = if xv < 0.0 {
        xv = -xv;
        -1i32
    } else {
        1
    };
    if flags.contains(DmsFlags::ROUND_MIN) {
        if !flags.contains(DmsFlags::ALLOW_361) {
            xv = normalize_degrees(xv + 0.5 / 60.0);
        }
    } else if flags.contains(DmsFlags::ROUND_SEC) {
        if !flags.contains(DmsFlags::ALLOW_361) {
            xv = normalize_degrees(xv + 0.5 / 3600.0);
        }
    } else if extra_prec {
        xv += 0.000000005 / 3600.0;
    } else {
        xv += 0.00005 / 3600.0;
    }
    let mut s = String::with_capacity(32);
    if flags.contains(DmsFlags::ZODIAC) {
        let mut izod = (xv / 30.0) as usize;
        if izod >= 12 {
            izod = 0;
        }
        xv %= 30.0;
        let kdeg = xv as i32;
        s.push_str(&format!("{:2} {} ", kdeg, ZODIAC_NAMES[izod]));
        xv -= kdeg as f64;
    } else {
        let kdeg = xv as i32;
        s.push_str(&format!(" {:3}{}", kdeg, c));
        xv -= kdeg as f64;
    }
    xv *= 60.0;
    let kmin = xv as i32;
    if flags.contains(DmsFlags::ZODIAC) && flags.contains(DmsFlags::ROUND_MIN) {
        s.push_str(&format!("{:2}", kmin));
    } else {
        s.push_str(&format!("{:2}'", kmin));
    }
    if flags.contains(DmsFlags::ROUND_MIN) {
        if sgn < 0 {
            insert_negative_sign(&mut s);
        }
        return s;
    }
    xv -= kmin as f64;
    xv *= 60.0;
    let ksec = xv as i32;
    if flags.contains(DmsFlags::ROUND_SEC) {
        s.push_str(&format!("{:2}\"", ksec));
    } else {
        s.push_str(&format!("{:2}", ksec));
    }
    if flags.contains(DmsFlags::ROUND_SEC) {
        if sgn < 0 {
            insert_negative_sign(&mut s);
        }
        return s;
    }
    xv -= ksec as f64;
    if extra_prec {
        let k = (xv * 100_000_000.0) as i64;
        s.push_str(&format!(".{:08}", k));
    } else {
        let k = (xv * 10000.0) as i32;
        s.push_str(&format!(".{:04}", k));
    }
    if sgn < 0 {
        insert_negative_sign(&mut s);
    }
    s
}

fn insert_negative_sign(s: &mut String) {
    if let Some(pos) = s.find(|c: char| c.is_ascii_digit())
        && pos > 0
    {
        let byte_pos = s.char_indices().nth(pos - 1).map(|(i, _)| i).unwrap_or(0);
        s.replace_range(byte_pos..byte_pos + 1, "-");
    }
}

fn round_flag(args: &SweTestArgs) -> DmsFlags {
    let mut f = DmsFlags::empty();
    if args.round_sec {
        f |= DmsFlags::ROUND_SEC;
    }
    if args.round_min {
        f |= DmsFlags::ROUND_MIN;
    }
    f
}

fn fmt_decimal(val: f64, extra_prec: bool) -> String {
    if extra_prec {
        format!("{:>11.11}", val)
    } else {
        format!("{:>11.7}", val)
    }
}

fn fmt_decimal_dist(val: f64, extra_prec: bool) -> String {
    if extra_prec {
        format!("{:>18.16}", val)
    } else {
        format!("{:>14.9}", val)
    }
}

fn fmt_cartesian(val: f64, extra_prec: bool) -> String {
    if extra_prec {
        format!("{:.17}", val)
    } else {
        format!("{:>14.9}", val)
    }
}

pub fn format_line(ctx: &FormatContext, eph: &Ephemeris) -> String {
    let fmt = &ctx.args.format;
    let gap = &ctx.args.gap;
    let rflag = round_flag(ctx.args);
    let ep = ctx.args.extra_precision;
    let chars: Vec<char> = fmt.chars().collect();
    let mut parts: Vec<String> = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let ch = chars[i];
        let col = format_char(ch, &chars, i, ctx, eph, rflag, ep, gap);
        if let Some(col) = col {
            parts.push(col);
        }
        // SS/ss: consume the second S/s
        if (ch == 'S' || ch == 's')
            && i + 1 < chars.len()
            && (chars[i + 1] == 'S' || chars[i + 1] == 's')
        {
            i += 1;
        }
        i += 1;
    }
    parts.join(gap)
}

fn format_char(
    ch: char,
    chars: &[char],
    idx: usize,
    ctx: &FormatContext,
    eph: &Ephemeris,
    rflag: DmsFlags,
    ep: bool,
    gap: &str,
) -> Option<String> {
    if ctx.is_label {
        return Some(label_for(ch));
    }
    if (ctx.is_house || ctx.is_ayanamsa)
        && matches!(
            ch,
            'b' | 'B'
                | 'r'
                | 'R'
                | 'x'
                | 'X'
                | 'u'
                | 'U'
                | 'Q'
                | 'n'
                | 'N'
                | 'f'
                | 'F'
                | 'j'
                | '+'
                | '-'
                | '*'
                | '/'
                | '='
        )
    {
        return None;
    }
    if ctx.is_ayanamsa && matches!(ch, 's' | 'S') {
        return None;
    }
    if ctx.args.horizontal && !ctx.is_first && matches!(ch, 'y' | 'Y' | 'J' | 'T' | 't') {
        return None;
    }
    match ch {
        'y' => Some(format!("{}", ctx.year)),
        'Y' => {
            let t2 = swisseph::date::julday(ctx.year, 1, 1, 0.0, ctx.cal);
            let y_frac = (ctx.jd - t2) / 365.0;
            Some(format!("{:.2}", ctx.year as f64 + y_frac))
        }
        'p' => Some(format!("{}", ctx.ipl)),
        'P' => Some(format!("{:<15}", ctx.name)),
        'J' => {
            let y_frac = (ctx.jd - ctx.jd.floor()) * 100.0;
            if y_frac.floor() != y_frac {
                Some(format!("{:.5}", ctx.jd))
            } else {
                Some(format!("{:.2}", ctx.jd))
            }
        }
        'T' => {
            let mut s = format!("{:02}.{:02}.{:04}", ctx.day, ctx.month, ctx.year);
            if ctx.cal == CalendarType::Julian {
                s.push('j');
            }
            if ctx.hour != 0.0 {
                let parts = swisseph::math::split_degrees(
                    ctx.hour,
                    if ctx.args.round_sec {
                        swisseph::flags::SplitDegFlags::ROUND_SEC
                    } else {
                        swisseph::flags::SplitDegFlags::empty()
                    },
                );
                if ctx.args.step_unit == crate::args::StepUnit::Seconds
                    && ctx.args.step_size.abs() < 1.0
                {
                    s.push_str(&format!(
                        " {}:{:02}:{:02}.{:.0}",
                        parts.degrees,
                        parts.minutes,
                        parts.seconds,
                        parts.second_fraction * 100.0,
                    ));
                } else {
                    s.push_str(&format!(
                        " {}:{:02}:{:02}",
                        parts.degrees, parts.minutes, parts.seconds
                    ));
                }
                if ctx.is_ut {
                    s.push_str(" UT");
                } else {
                    s.push_str(" TT");
                }
            }
            Some(s)
        }
        't' => Some(format!(
            "{:02}{:02}{:02}",
            ctx.year % 100,
            ctx.month,
            ctx.day
        )),
        'L' => Some(dms(ctx.data[0], rflag, ep)),
        'l' => {
            if rflag.contains(DmsFlags::ROUND_MIN) {
                Some(format!("{:>6.2}", ctx.data[0]))
            } else {
                Some(fmt_decimal(ctx.data[0], ep))
            }
        }
        'Z' => Some(dms(ctx.data[0], rflag | DmsFlags::ZODIAC, ep)),
        'S' | 's' => format_speed(ch, chars, idx, ctx, rflag, ep, gap),
        'B' => Some(dms(ctx.data[1], rflag, ep)),
        'b' => Some(fmt_decimal(ctx.data[1], ep)),
        'R' => Some(fmt_decimal_dist(ctx.data[2], ep)),
        'r' => {
            if ctx.body == Some(Body::Moon) {
                if let Ok(pheno) = eph.pheno(ctx.tjd_tt, Body::Moon, ctx.args.build_iflag()) {
                    Some(format!("{:>13.5}\"", pheno.0.horizontal_parallax * 3600.0))
                } else {
                    Some(fmt_decimal_dist(ctx.data[2], ep))
                }
            } else {
                Some(fmt_decimal_dist(ctx.data[2], ep))
            }
        }
        'W' => Some(format!("{:>14.9}", ctx.data[2] * AUNIT_TO_LIGHTYEAR)),
        'w' => Some(format!("{:>14.9}", ctx.data[2] * AUNIT_TO_KM)),
        'q' => {
            if ctx.is_label {
                return Some("reldist".into());
            }
            let body = ctx.body?;
            let iflagi = ctx.args.build_iflag()
                & (CalcFlags::SWIEPH
                    | CalcFlags::JPLEPH
                    | CalcFlags::MOSEPH
                    | CalcFlags::HELCTR
                    | CalcFlags::BARYCTR);
            let dar = match eph.orbit_max_min_true_distance(ctx.tjd_tt, body, iflagi) {
                Ok((dmax, dmin, dtrue)) => {
                    if (dmax - dmin).abs() < f64::EPSILON {
                        0
                    } else {
                        ((1.0 - (dtrue - dmin) / (dmax - dmin)) * 1000.0 + 0.5) as i32
                    }
                }
                Err(_) => 0,
            };
            Some(format!("{dar:>5}"))
        }
        'A' => {
            if let Some(ref xequ) = ctx.xequ {
                Some(dms(xequ[0] / 15.0, rflag | DmsFlags::EQUATORIAL, ep))
            } else {
                Some("  --RA--".into())
            }
        }
        'a' => {
            if let Some(ref xequ) = ctx.xequ {
                Some(fmt_decimal(xequ[0], ep))
            } else {
                Some("  --RA--".into())
            }
        }
        'D' => {
            if let Some(ref xequ) = ctx.xequ {
                Some(dms(xequ[1], rflag, ep))
            } else {
                Some("  --dec--".into())
            }
        }
        'd' => {
            if let Some(ref xequ) = ctx.xequ {
                Some(fmt_decimal(xequ[1], ep))
            } else {
                Some("  --dec--".into())
            }
        }
        'I' => {
            if let Some(ref xaz) = ctx.xaz {
                Some(dms(xaz[0], rflag, ep))
            } else {
                Some("  --az--".into())
            }
        }
        'i' => {
            if let Some(ref xaz) = ctx.xaz {
                Some(fmt_decimal(xaz[0], ep))
            } else {
                Some("  --az--".into())
            }
        }
        'H' => {
            if let Some(ref xaz) = ctx.xaz {
                Some(dms(xaz[1], rflag, ep))
            } else {
                Some("  --alt--".into())
            }
        }
        'h' => {
            if let Some(ref xaz) = ctx.xaz {
                Some(fmt_decimal(xaz[1], ep))
            } else {
                Some("  --alt--".into())
            }
        }
        'K' => {
            if let Some(ref xaz) = ctx.xaz {
                Some(dms(xaz[2], rflag, ep))
            } else {
                Some("  --alt--".into())
            }
        }
        'k' => {
            if let Some(ref xaz) = ctx.xaz {
                Some(fmt_decimal(xaz[2], ep))
            } else {
                Some("  --alt--".into())
            }
        }
        'U' | 'X' => {
            if let Some(ref xcart) = ctx.xcart {
                let ar = if ch == 'U' {
                    (xcart[0] * xcart[0] + xcart[1] * xcart[1] + xcart[2] * xcart[2]).sqrt()
                } else {
                    1.0
                };
                Some(format!(
                    "{}{}{}{}{}",
                    fmt_cartesian(xcart[0] / ar, ep),
                    gap,
                    fmt_cartesian(xcart[1] / ar, ep),
                    gap,
                    fmt_cartesian(xcart[2] / ar, ep),
                ))
            } else {
                Some("  --cart--".into())
            }
        }
        'u' | 'x' => {
            if let Some(ref xecart) = ctx.xecart {
                let ar = if ch == 'u' {
                    (xecart[0] * xecart[0] + xecart[1] * xecart[1] + xecart[2] * xecart[2]).sqrt()
                } else {
                    1.0
                };
                Some(format!(
                    "{}{}{}{}{}",
                    fmt_cartesian(xecart[0] / ar, ep),
                    gap,
                    fmt_cartesian(xecart[1] / ar, ep),
                    gap,
                    fmt_cartesian(xecart[2] / ar, ep),
                ))
            } else {
                Some("  --ecart--".into())
            }
        }
        'G' => {
            if let Some(hpos) = ctx.hpos {
                Some(dms(hpos, rflag, ep))
            } else {
                Some("  --hpos--".into())
            }
        }
        'g' => {
            if let Some(hpos) = ctx.hpos {
                Some(fmt_decimal(hpos, ep))
            } else {
                Some("  --hpos--".into())
            }
        }
        'j' => {
            if let Some(hposj) = ctx.hposj {
                Some(fmt_decimal(hposj, ep))
            } else {
                Some("  --hpos--".into())
            }
        }
        'Q' => {
            let mut s = format!("{:<15}", ctx.name);
            s.push_str(&dms(ctx.data[0], rflag, ep));
            s.push_str(&dms(ctx.data[1], rflag, ep));
            s.push_str(&format!("  {:>14.9}", ctx.data[2]));
            s.push_str(&dms(ctx.data[3], rflag, ep));
            s.push_str(&dms(ctx.data[4], rflag, ep));
            s.push_str(&format!("  {:>14.9}\n", ctx.data[5]));
            if let Some(ref xequ) = ctx.xequ {
                s.push_str(&format!("               {}", dms(xequ[0], rflag, ep)));
                s.push_str(&dms(xequ[1], rflag, ep));
                s.push_str(&format!("                {}", dms(xequ[3], rflag, ep)));
                s.push_str(&dms(xequ[4], rflag, ep));
            }
            Some(s)
        }
        '+' => {
            if let Some(ref attr) = ctx.attr {
                if ctx.args.format.contains('l') {
                    Some(fmt_decimal(attr[0], ep))
                } else {
                    Some(dms(attr[0], rflag, ep))
                }
            } else {
                Some("  --phase--".into())
            }
        }
        '-' => {
            if let Some(ref attr) = ctx.attr {
                Some(format!("  {:>14.9}", attr[1]))
            } else {
                Some("  --phase--".into())
            }
        }
        '*' => {
            if let Some(ref attr) = ctx.attr {
                if ctx.args.format.contains('l') {
                    Some(fmt_decimal(attr[2], ep))
                } else {
                    Some(dms(attr[2], rflag, ep))
                }
            } else {
                Some("  --elong--".into())
            }
        }
        '/' => {
            if let Some(ref attr) = ctx.attr {
                Some(dms(attr[3], rflag, ep))
            } else {
                Some("  --diam--".into())
            }
        }
        '=' => {
            if let Some(ref attr) = ctx.attr {
                Some(format!("  {:>6.3}m", attr[4]))
            } else {
                Some("  --magn--".into())
            }
        }
        'V' | 'v' => {
            static HEXA: [i32; 64] = [
                1, 43, 14, 34, 9, 5, 26, 11, 10, 58, 38, 54, 61, 60, 41, 19, 13, 49, 30, 55, 37,
                63, 22, 36, 25, 17, 21, 51, 42, 3, 27, 24, 2, 23, 8, 20, 16, 35, 45, 12, 15, 52,
                39, 53, 62, 56, 31, 33, 7, 4, 29, 59, 40, 64, 47, 6, 46, 18, 48, 57, 32, 50, 28,
                44,
            ];
            let xhds = normalize_degrees(ctx.data[0] - 223.25);
            let ihex = (xhds / 5.625).floor() as usize;
            let iline = ((xhds / 0.9375).floor() as i32) % 6 + 1;
            let igate = if ihex < 64 { HEXA[ihex] } else { 0 };
            if ch == 'V' {
                let pct = ((100.0 * (xhds / 0.9375).fract()).round()) as i32;
                Some(format!("{:2}.{} {:2}%", igate, iline, pct))
            } else {
                Some(format!("{:2}.{}", igate, iline))
            }
        }
        'm' => {
            if let Some(ref xequ) = ctx.xequ {
                if let Some(armc) = ctx.armc {
                    let md = swisseph::math::diff_degrees(xequ[0], armc).abs();
                    Some(fmt_decimal(md, ep))
                } else {
                    Some("  --MD--".into())
                }
            } else {
                Some("  --MD--".into())
            }
        }
        'z' => {
            if let Some(ref xaz) = ctx.xaz {
                let zd = 90.0 - xaz[1];
                Some(fmt_decimal(zd, ep))
            } else {
                Some("  --ZD--".into())
            }
        }
        'N' | 'n' => {
            let body = ctx.body?;
            let imeth = if ch == 'n' {
                swisseph::nodaps::NodApsMethod::MEAN
            } else {
                swisseph::nodaps::NodApsMethod::OSCU
            };
            let result = eph
                .nod_aps(ctx.tjd_tt, body, ctx.args.build_iflag(), imeth)
                .ok()?;
            let use_dms_here = ctx.args.dms;
            let mut s = if use_dms_here {
                dms(result.ascending[0], rflag | DmsFlags::ZODIAC, ep)
            } else {
                fmt_decimal(result.ascending[0], ep)
            };
            s.push_str(gap);
            if use_dms_here {
                s.push_str(&dms(result.descending[0], rflag | DmsFlags::ZODIAC, ep));
            } else {
                s.push_str(&fmt_decimal(result.descending[0], ep));
            }
            Some(s)
        }
        'F' | 'f' => {
            let body = ctx.body?;
            let imeth = if ch == 'f' {
                swisseph::nodaps::NodApsMethod::MEAN
            } else {
                swisseph::nodaps::NodApsMethod::OSCU
            };
            let iflag = ctx.args.build_iflag();
            let result = eph.nod_aps(ctx.tjd_tt, body, iflag, imeth).ok()?;
            let mut s = fmt_decimal(result.perihelion[0], ep);
            s.push_str(gap);
            s.push_str(&fmt_decimal(result.aphelion[0], ep));
            // focal point
            let imeth_foc = imeth | swisseph::nodaps::NodApsMethod::FOPOINT;
            if let Ok(foc) = eph.nod_aps(ctx.tjd_tt, body, iflag, imeth_foc) {
                s.push_str(gap);
                s.push_str(&fmt_decimal(foc.aphelion[0], ep));
            }
            Some(s)
        }
        'c' | 'e' | 'o' => None,
        _ => None,
    }
}

fn format_speed(
    ch: char,
    chars: &[char],
    idx: usize,
    ctx: &FormatContext,
    rflag: DmsFlags,
    ep: bool,
    gap: &str,
) -> Option<String> {
    let has_double_s = idx + 1 < chars.len() && (chars[idx + 1] == 'S' || chars[idx + 1] == 's');
    let has_cart = chars.iter().any(|&c| matches!(c, 'X' | 'U' | 'x' | 'u'));
    if has_double_s || has_cart {
        let mut parts: Vec<String> = Vec::new();
        for &fc in chars {
            let col = match fc {
                'L' | 'Z' => {
                    if ctx.is_label {
                        "lon/day".into()
                    } else {
                        dms(ctx.data[3], rflag | DmsFlags::ALLOW_361, ep)
                    }
                }
                'l' => {
                    if ctx.is_label {
                        "lon/day".into()
                    } else {
                        fmt_decimal(ctx.data[3], ep)
                    }
                }
                'B' => {
                    if ctx.is_label {
                        "lat/day".into()
                    } else {
                        dms(ctx.data[4], rflag | DmsFlags::ALLOW_361, ep)
                    }
                }
                'b' => {
                    if ctx.is_label {
                        "lat/day".into()
                    } else {
                        fmt_decimal(ctx.data[4], ep)
                    }
                }
                'A' => {
                    if ctx.is_label {
                        "RA/day".into()
                    } else if let Some(ref xequ) = ctx.xequ {
                        dms(
                            xequ[3] / 15.0,
                            rflag | DmsFlags::EQUATORIAL | DmsFlags::ALLOW_361,
                            ep,
                        )
                    } else {
                        "  --RA--".into()
                    }
                }
                'a' => {
                    if ctx.is_label {
                        "RA/day".into()
                    } else if let Some(ref xequ) = ctx.xequ {
                        fmt_decimal(xequ[3], ep)
                    } else {
                        "  --RA--".into()
                    }
                }
                'D' => {
                    if ctx.is_label {
                        "dcl/day".into()
                    } else if let Some(ref xequ) = ctx.xequ {
                        dms(xequ[4], rflag | DmsFlags::ALLOW_361, ep)
                    } else {
                        "  --dec--".into()
                    }
                }
                'd' => {
                    if ctx.is_label {
                        "dcl/day".into()
                    } else if let Some(ref xequ) = ctx.xequ {
                        fmt_decimal(xequ[4], ep)
                    } else {
                        "  --dec--".into()
                    }
                }
                'R' | 'r' => {
                    if ctx.is_label {
                        "AU/day".into()
                    } else {
                        fmt_decimal_dist(ctx.data[5], ep)
                    }
                }
                'U' | 'X' => {
                    if let Some(ref xcart) = ctx.xcart {
                        let ar = if fc == 'U' {
                            (xcart[0] * xcart[0] + xcart[1] * xcart[1] + xcart[2] * xcart[2]).sqrt()
                        } else {
                            1.0
                        };
                        format!(
                            "{}{}{}{}{}",
                            fmt_cartesian(xcart[3] / ar, ep),
                            gap,
                            fmt_cartesian(xcart[4] / ar, ep),
                            gap,
                            fmt_cartesian(xcart[5] / ar, ep),
                        )
                    } else {
                        "  --cart--".into()
                    }
                }
                'u' | 'x' => {
                    if let Some(ref xecart) = ctx.xecart {
                        let ar = if fc == 'u' {
                            (xecart[0] * xecart[0] + xecart[1] * xecart[1] + xecart[2] * xecart[2])
                                .sqrt()
                        } else {
                            1.0
                        };
                        format!(
                            "{}{}{}{}{}",
                            fmt_cartesian(xecart[3] / ar, ep),
                            gap,
                            fmt_cartesian(xecart[4] / ar, ep),
                            gap,
                            fmt_cartesian(xecart[5] / ar, ep),
                        )
                    } else {
                        "  --ecart--".into()
                    }
                }
                _ => continue,
            };
            parts.push(col);
        }
        Some(parts.join(gap))
    } else if ch == 'S' {
        Some(dms(ctx.data[3], rflag | DmsFlags::ALLOW_361, ep))
    } else {
        // 's'
        Some(fmt_decimal(ctx.data[3], ep))
    }
}

fn label_for(ch: char) -> String {
    match ch {
        'y' | 'Y' => "year".into(),
        'p' => "obj.nr".into(),
        'P' => format!("{:<15}", "name"),
        'J' => "julday".into(),
        'T' => "date    ".into(),
        't' => "date".into(),
        'L' | 'Z' | 'l' => "long.".into(),
        'S' | 's' => "deg/day".into(),
        'B' => "lat.    ".into(),
        'b' => "lat.    ".into(),
        'R' | 'r' => "distAU   ".into(),
        'W' => "distLY   ".into(),
        'w' => "distkm   ".into(),
        'q' => "reldist".into(),
        'A' | 'a' => "RA      ".into(),
        'D' | 'd' => "decl      ".into(),
        'I' | 'i' => "azimuth".into(),
        'H' | 'h' => "height".into(),
        'K' | 'k' => "hgtApp".into(),
        'U' | 'X' | 'u' | 'x' => "x0".into(),
        'G' | 'g' => "housPos".into(),
        'j' => "houseNr".into(),
        'Q' => "Q".into(),
        '+' => "phase".into(),
        '-' => "phase".into(),
        '*' => "elong".into(),
        '/' => "diamet".into(),
        '=' => "magn".into(),
        'V' | 'v' => "hds".into(),
        'm' => "MD      ".into(),
        'z' => "ZD      ".into(),
        'N' | 'n' => "nodAsc".into(),
        'F' | 'f' => "peri".into(),
        _ => String::new(),
    }
}
