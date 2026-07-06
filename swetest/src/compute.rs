use swisseph::flags::CalcFlags;
use swisseph::types::{Body, CalendarType};
use swisseph::{Ephemeris, EphemerisConfig};

use crate::args::{BodySpec, DiffMode, StepUnit, SweTestArgs, TimeMode};
use crate::format::{self, FormatContext, FormatNeeds};

const VERSION: &str = "0.1.0";
pub(crate) const GREG_BOUNDARY_JD: f64 = 2299161.0; // 15 Oct 1582

pub(crate) fn calendar_for_jd(jd: f64) -> CalendarType {
    if jd < GREG_BOUNDARY_JD {
        CalendarType::Julian
    } else {
        CalendarType::Gregorian
    }
}

fn calendar_for_date(year: i32, month: i32, day: i32) -> CalendarType {
    if year < 1582 || (year == 1582 && (month < 10 || (month == 10 && day < 15))) {
        CalendarType::Julian
    } else {
        CalendarType::Gregorian
    }
}

fn parse_date_string(s: &str) -> (i32, i32, i32) {
    let parts: Vec<&str> = s.split('.').collect();
    match parts.len() {
        3 => {
            let d = parts[0].parse::<i32>().unwrap_or(1);
            let m = parts[1].parse::<i32>().unwrap_or(1);
            let mut y = parts[2].parse::<i32>().unwrap_or(2000);
            if y >= 0 && y < 100 {
                y += 2000;
            }
            (y, m, d)
        }
        2 => {
            let d = parts[0].parse::<i32>().unwrap_or(1);
            let m = parts[1].parse::<i32>().unwrap_or(1);
            (2000, m, d)
        }
        1 => {
            let mut y = parts[0].parse::<i32>().unwrap_or(2000);
            if y >= 0 && y < 100 {
                y += 2000;
            }
            (y, 1, 1)
        }
        _ => (2000, 1, 1),
    }
}

fn parse_time_string(s: &str) -> f64 {
    if s.is_empty() {
        return 0.0;
    }
    let parts: Vec<&str> = if s.contains(':') {
        s.split(':').collect()
    } else {
        s.split('.').collect()
    };
    let h = parts
        .first()
        .and_then(|p| p.parse::<f64>().ok())
        .unwrap_or(0.0);
    if parts.len() == 1 {
        if s.contains(':') {
            return h;
        }
        return s.parse::<f64>().unwrap_or(0.0);
    }
    let m = parts
        .get(1)
        .and_then(|p| p.parse::<f64>().ok())
        .unwrap_or(0.0);
    let sec = parts
        .get(2)
        .and_then(|p| p.parse::<f64>().ok())
        .unwrap_or(0.0);
    h + m / 60.0 + sec / 3600.0
}

pub(crate) struct EpochInfo {
    pub(crate) tjd_ut: f64,
    pub(crate) tjd_tt: f64,
    pub(crate) year: i32,
    pub(crate) month: i32,
    pub(crate) day: i32,
    pub(crate) hour: f64,
    pub(crate) cal: CalendarType,
    pub(crate) is_ut: bool,
}

pub(crate) fn resolve_start_jd(args: &SweTestArgs, config: &EphemerisConfig) -> EpochInfo {
    let thour = parse_time_string(&args.time_input);
    let is_ut = matches!(
        args.time_mode,
        TimeMode::UT | TimeMode::UTC | TimeMode::LMT | TimeMode::LAT
    );

    let (year, month, day, jd_input) = match &args.begin_date {
        None => {
            let jd = swisseph::date::julday(2000, 1, 1, 0.0, CalendarType::Gregorian);
            (2000, 1, 1, jd)
        }
        Some(s) => {
            let trimmed = s.trim();
            if trimmed.starts_with('j') || trimmed.starts_with('J') {
                let jd: f64 = trimmed[1..].parse().unwrap_or(2451545.0);
                let cal = calendar_for_jd(jd);
                let (y, m, d, _h) = swisseph::date::revjul(jd, cal);
                (y, m, d, jd)
            } else {
                let (y, m, d) = parse_date_string(trimmed);
                let cal = calendar_for_date(y, m, d);
                let jd = swisseph::date::julday(y, m, d, 0.0, cal);
                (y, m, d, jd)
            }
        }
    };

    let cal = calendar_for_date(year, month, day);
    let tjd = jd_input + thour / 24.0;

    let (tjd_ut, tjd_tt) = match args.time_mode {
        TimeMode::ET => {
            let dt = swisseph::deltat::calc_deltat(tjd, config);
            (tjd - dt, tjd)
        }
        TimeMode::UT => {
            let dt = swisseph::deltat::calc_deltat(tjd, config);
            (tjd, tjd + dt)
        }
        TimeMode::LMT => {
            let tjd_lmt_ut = tjd - args.geo_longitude / 360.0;
            let dt = swisseph::deltat::calc_deltat(tjd_lmt_ut, config);
            (tjd_lmt_ut, tjd_lmt_ut + dt)
        }
        TimeMode::UTC | TimeMode::LAT => {
            let dt = swisseph::deltat::calc_deltat(tjd, config);
            (tjd, tjd + dt)
        }
    };

    EpochInfo {
        tjd_ut,
        tjd_tt,
        year,
        month,
        day,
        hour: thour,
        cal,
        is_ut,
    }
}

fn step_jd(
    args: &SweTestArgs,
    istep: i32,
    base_year: i32,
    base_month: i32,
    base_day: i32,
    base_hour: f64,
    base_jd: f64,
) -> f64 {
    if istep <= 1 {
        return base_jd;
    }
    let offset = (istep - 1) as f64 * args.step_size;
    match args.step_unit {
        StepUnit::Days => base_jd + offset,
        StepUnit::Years => {
            let y = base_year + (istep - 1) * args.step_size as i32;
            let cal = calendar_for_date(y, base_month, base_day);
            swisseph::date::julday(y, base_month, base_day, base_hour, cal)
        }
        StepUnit::Months => {
            let total_months = base_month + (istep - 1) * args.step_size as i32;
            let y = base_year + (total_months - 1) / 12;
            let m = (total_months - 1) % 12 + 1;
            let cal = calendar_for_date(y, m, base_day);
            swisseph::date::julday(y, m, base_day, base_hour, cal)
        }
        StepUnit::Minutes => base_jd + offset / (24.0 * 60.0),
        StepUnit::Seconds => base_jd + offset / (24.0 * 3600.0),
    }
}

pub(crate) fn format_time(hour: f64) -> String {
    let mut total_sec = (hour * 3600.0).round() as i64;
    if total_sec < 0 {
        total_sec += 86400;
    }
    let h = (total_sec / 3600) % 24;
    let m = (total_sec % 3600) / 60;
    let s = total_sec % 60;
    format!("{h:>2}:{m:02}:{s:02}")
}

fn sidereal_mode_name(sid_mode: i32) -> String {
    swisseph::types::SiderealMode::try_from(sid_mode)
        .ok()
        .and_then(|m| m.name().map(|s| s.to_owned()))
        .unwrap_or_else(|| format!("mode {sid_mode}"))
}

fn print_header(args: &SweTestArgs, eph: &Ephemeris, info: &EpochInfo, iflag: CalcFlags) {
    if !args.with_header {
        return;
    }

    let cal_str = if info.cal == CalendarType::Gregorian {
        "greg."
    } else {
        "jul."
    };
    let time_str = format_time(info.hour);
    let time_label = if info.is_ut { "UT" } else { "ET" };
    println!(
        "date (dmy) {}.{}.{} {cal_str}   {time_str} {time_label}        version {VERSION}",
        info.day, info.month, info.year,
    );

    if info.is_ut {
        let dt_sec = (info.tjd_tt - info.tjd_ut) * 86400.0;
        println!("UT: {:.7}     delta t: {:.6} sec", info.tjd_ut, dt_sec);
        let tt_time = info.hour + dt_sec / 3600.0;
        println!("ET: {:.7}     {}", info.tjd_tt, format_time(tt_time));
    } else {
        println!("ET: {:.7}", info.tjd_tt);
        let dt_sec = (info.tjd_tt - info.tjd_ut) * 86400.0;
        let ut_time = info.hour - dt_sec / 3600.0;
        println!(
            "UT: {:.7}     delta t: {:.6} sec    {}",
            info.tjd_ut,
            dt_sec,
            format_time(ut_time),
        );
    }

    if let Ok(ecl_nut) = eph.calc(info.tjd_tt, Body::EclipticNutation, CalcFlags::empty()) {
        let d = &ecl_nut.data;
        println!("Epsilon (true, mean)   {:.7}   {:.7}", d[0], d[1]);
        println!("Nutation               {:.7}   {:.7}", d[2], d[3]);
    }

    if args.sidereal {
        if let Ok(aya) = eph.get_ayanamsa_ex(info.tjd_tt, iflag) {
            let mode_name = sidereal_mode_name(args.sid_mode);
            println!("Ayanamsa ({mode_name})   {aya:.7}");
        }
    }

    if args.topocentric || args.have_geopos {
        println!(
            "geo. long {:.4}, lat {:.4}, alt {:.1} m",
            args.geo_longitude, args.geo_latitude, args.geo_elevation,
        );
    }

    if args.do_houses {
        let hsys = swisseph::types::HouseSystem::try_from(args.house_system as u8)
            .map(|h| h.name().to_owned())
            .unwrap_or_else(|_| format!("{}", args.house_system));
        println!("Houses system {hsys}");
    }

    println!();
}

pub(crate) fn parse_int_arg(s: &Option<String>) -> Option<i32> {
    s.as_ref().and_then(|v| v.parse::<i32>().ok())
}

pub(crate) fn make_asteroid_body(num: i32) -> Option<Body> {
    swisseph::types::AsteroidId::new(num)
        .ok()
        .map(Body::Asteroid)
}

fn make_plmoon_body(id: i32) -> Option<Body> {
    swisseph::types::PlanetMoonId::new(id)
        .ok()
        .map(Body::PlanetMoon)
}

fn make_fictitious_body(id: i32) -> Option<Body> {
    Body::fictitious(id).ok()
}

pub(crate) fn body_to_ipl(body: Body) -> i32 {
    match body {
        Body::Sun => 0,
        Body::Moon => 1,
        Body::Mercury => 2,
        Body::Venus => 3,
        Body::Mars => 4,
        Body::Jupiter => 5,
        Body::Saturn => 6,
        Body::Uranus => 7,
        Body::Neptune => 8,
        Body::Pluto => 9,
        Body::MeanNode => 10,
        Body::TrueNode => 11,
        Body::MeanApogee => 12,
        Body::OscuApogee => 13,
        Body::Earth => 14,
        Body::Chiron => 15,
        Body::Pholus => 16,
        Body::Ceres => 17,
        Body::Pallas => 18,
        Body::Juno => 19,
        Body::Vesta => 20,
        Body::IntpApogee => 21,
        Body::IntpPerigee => 22,
        Body::EclipticNutation => -1,
        Body::Fictitious(id) => id.raw_id(),
        Body::Asteroid(id) => 10000 + id.mpc_number(),
        Body::PlanetMoon(id) => id.encoded(),
    }
}

pub(crate) fn body_name(eph: &Ephemeris, spec: &BodySpec, args: &SweTestArgs) -> String {
    match spec {
        BodySpec::Planet(body) => eph.get_planet_name(*body),
        BodySpec::Asteroid => {
            if let Some(body) = parse_int_arg(&args.asteroid_number).and_then(make_asteroid_body) {
                eph.get_planet_name(body)
            } else {
                "asteroid ?".into()
            }
        }
        BodySpec::FixedStar => args.star_name.clone().unwrap_or_else(|| "star ?".into()),
        BodySpec::PlanetMoon => {
            if let Some(id) = parse_int_arg(&args.planet_moon) {
                format!("planet moon {id}")
            } else {
                "planet moon ?".into()
            }
        }
        BodySpec::Fictitious => {
            if let Some(body) = parse_int_arg(&args.fictitious).and_then(make_fictitious_body) {
                eph.get_planet_name(body)
            } else {
                "fictitious ?".into()
            }
        }
        BodySpec::EclipticNutation => "Ecl. Nut.".into(),
        BodySpec::Labels => String::new(),
        BodySpec::DeltaT => "Delta T".into(),
        BodySpec::TimeEquation => "Time Equ.".into(),
        BodySpec::SiderealTime => "Sid. Time".into(),
        BodySpec::Ayanamsha => "Ayanamsha".into(),
    }
}

pub(crate) fn resolve_body(spec: &BodySpec, args: &SweTestArgs) -> Option<Body> {
    match spec {
        BodySpec::Planet(body) => Some(*body),
        BodySpec::Asteroid => parse_int_arg(&args.asteroid_number).and_then(make_asteroid_body),
        BodySpec::PlanetMoon => parse_int_arg(&args.planet_moon).and_then(make_plmoon_body),
        BodySpec::Fictitious => parse_int_arg(&args.fictitious).and_then(make_fictitious_body),
        BodySpec::EclipticNutation => Some(Body::EclipticNutation),
        _ => None,
    }
}

fn compute_supplementary(
    eph: &Ephemeris,
    body: Body,
    star: Option<&str>,
    tjd_tt: f64,
    tjd_ut: f64,
    iflag: CalcFlags,
    needs: &FormatNeeds,
    args: &SweTestArgs,
) -> (
    Option<[f64; 6]>,
    Option<[f64; 3]>,
    Option<[f64; 6]>,
    Option<[f64; 6]>,
    Option<f64>,
    Option<f64>,
    Option<f64>,
    Option<[f64; 6]>,
) {
    let mut xequ = None;
    let mut xaz = None;
    let mut xcart = None;
    let mut xecart = None;
    let mut hpos = None;
    let mut hposj = None;
    let mut armc_val = None;
    let mut attr = None;

    if needs.equatorial {
        let iflag2 = iflag | CalcFlags::EQUATORIAL;
        let result = if let Some(s) = star {
            eph.fixstar2(s, tjd_tt, iflag2).ok().map(|(_, r)| r)
        } else {
            eph.calc(tjd_tt, body, iflag2).ok()
        };
        if let Some(r) = result {
            xequ = Some(r.data);
        }
    }

    if needs.azalt {
        let whicheph = iflag & (CalcFlags::SWIEPH | CalcFlags::JPLEPH | CalcFlags::MOSEPH);
        let iflgt = whicheph | CalcFlags::EQUATORIAL | CalcFlags::TOPOCTR;
        let topo_result = if let Some(s) = star {
            eph.fixstar2(s, tjd_tt, iflgt).ok().map(|(_, r)| r)
        } else {
            eph.calc(tjd_tt, body, iflgt).ok()
        };
        if let Some(r) = topo_result {
            let geopos = [args.geo_longitude, args.geo_latitude, args.geo_elevation];
            let xin = [r.data[0], r.data[1]];
            let az_result = eph.azalt(
                tjd_ut,
                swisseph::azalt::AzAltDir::EquToHor,
                geopos,
                args.atmosphere[0],
                args.atmosphere[1],
                0.0,
                xin,
            );
            xaz = Some(az_result);
        }
    }

    if needs.ecl_cartesian {
        let iflag2 = iflag | CalcFlags::XYZ;
        let result = if let Some(s) = star {
            eph.fixstar2(s, tjd_tt, iflag2).ok().map(|(_, r)| r)
        } else {
            eph.calc(tjd_tt, body, iflag2).ok()
        };
        if let Some(r) = result {
            xcart = Some(r.data);
        }
    }

    if needs.equ_cartesian {
        let iflag2 = iflag | CalcFlags::XYZ | CalcFlags::EQUATORIAL;
        let result = if let Some(s) = star {
            eph.fixstar2(s, tjd_tt, iflag2).ok().map(|(_, r)| r)
        } else {
            eph.calc(tjd_tt, body, iflag2).ok()
        };
        if let Some(r) = result {
            xecart = Some(r.data);
        }
    }

    if needs.house_pos {
        let sidt = swisseph::sidereal_time::sidereal_time(tjd_ut, eph.config());
        let armc = swisseph::math::normalize_degrees(sidt * 15.0 + args.geo_longitude);
        armc_val = Some(armc);
        if let Ok(ecl_nut) = eph.calc(tjd_tt, Body::EclipticNutation, CalcFlags::empty()) {
            let xobl = ecl_nut.data[0];
            let hsys = swisseph::types::HouseSystem::try_from(args.house_system as u8)
                .unwrap_or(swisseph::types::HouseSystem::Placidus);
            // Get the ecliptic position for house_pos
            let xsv = if let Some(s) = star {
                eph.fixstar2(s, tjd_tt, iflag).ok().map(|(_, r)| r.data)
            } else {
                eph.calc(tjd_tt, body, iflag).ok().map(|r| r.data)
            };
            if let Some(xsv) = xsv {
                let xpin = [xsv[0], xsv[1]];
                if let Ok(hp) =
                    swisseph::houses::house_pos(armc, args.geo_latitude, xobl, hsys, xpin, None)
                {
                    hposj = Some(hp);
                    if hsys == swisseph::types::HouseSystem::Gauquelin {
                        hpos = Some((hp - 1.0) * 10.0);
                    } else {
                        hpos = Some((hp - 1.0) * 30.0);
                    }
                }
            }
        }
    } else if needs.equatorial {
        // armc needed for meridian distance even without full house_pos
        let sidt = swisseph::sidereal_time::sidereal_time(tjd_ut, eph.config());
        armc_val = Some(swisseph::math::normalize_degrees(
            sidt * 15.0 + args.geo_longitude,
        ));
    }

    if needs.phenomena && star.is_none() {
        if let Ok(pheno) = eph.pheno(tjd_tt, body, iflag) {
            attr = Some([
                pheno.0.phase_angle,
                pheno.0.phase,
                pheno.0.elongation,
                pheno.0.apparent_diameter,
                pheno.0.apparent_magnitude,
                pheno.0.horizontal_parallax,
            ]);
        }
    }

    (xequ, xaz, xcart, xecart, hpos, hposj, armc_val, attr)
}

fn compute_body(
    eph: &Ephemeris,
    spec: &BodySpec,
    args: &SweTestArgs,
    tjd_tt: f64,
    tjd_ut: f64,
    iflag: CalcFlags,
    needs: &FormatNeeds,
    info: &EpochInfo,
) {
    let name = body_name(eph, spec, args);

    match spec {
        BodySpec::Labels => return,
        BodySpec::DeltaT => {
            let dt = swisseph::deltat::calc_deltat(tjd_ut, eph.config());
            let dt_sec = dt * 86400.0;
            println!("{name:<15} {dt_sec:.6} sec");
            return;
        }
        BodySpec::TimeEquation => {
            match eph.time_equ(tjd_ut) {
                Ok(e) => {
                    let e_sec = e * 86400.0;
                    let sign = if e_sec < 0.0 { "-" } else { "" };
                    let abs_sec = e_sec.abs();
                    let m = (abs_sec / 60.0) as i32;
                    let s = abs_sec - m as f64 * 60.0;
                    println!("{name:<15} {sign}{m}m {s:.2}s");
                }
                Err(e) => println!("{name:<15} error: {e}"),
            }
            return;
        }
        BodySpec::SiderealTime => {
            let sidt = swisseph::sidereal_time::sidereal_time(tjd_ut, eph.config());
            println!("{name:<15} {}", format_time(sidt));
            return;
        }
        BodySpec::Ayanamsha => {
            match eph.get_ayanamsa_ex(tjd_tt, iflag) {
                Ok(aya) => println!("{name:<15} {aya:.7}"),
                Err(e) => println!("{name:<15} error: {e}"),
            }
            return;
        }
        _ => {}
    }

    // For bodies that go through the format engine
    let is_fixstar = matches!(spec, BodySpec::FixedStar);
    let star_name = if is_fixstar {
        args.star_name.as_deref()
    } else {
        None
    };

    // Primary computation
    let (calc_name, data) = if is_fixstar {
        if let Some(ref star) = args.star_name {
            match eph.fixstar2(star, tjd_tt, iflag) {
                Ok((canonical, result)) => (canonical, result.data),
                Err(e) => {
                    println!("{name:<15} error: {e}");
                    return;
                }
            }
        } else {
            println!("{name:<15} error: no star name (-xf)");
            return;
        }
    } else if let Some(body) = resolve_body(spec, args) {
        match eph.calc(tjd_tt, body, iflag) {
            Ok(result) => (name.clone(), result.data),
            Err(e) => {
                println!("{name:<15} error: {e}");
                return;
            }
        }
    } else {
        println!("{name:<15} error: invalid body specification");
        return;
    };

    let body = resolve_body(spec, args);

    // Supplementary computations
    let (xequ, xaz, xcart, xecart, hpos, hposj, armc, attr) = if let Some(b) = body {
        compute_supplementary(eph, b, star_name, tjd_tt, tjd_ut, iflag, needs, args)
    } else {
        (None, None, None, None, None, None, None, None)
    };

    let ipl = match spec {
        BodySpec::Planet(b) => body_to_ipl(*b),
        BodySpec::Asteroid => parse_int_arg(&args.asteroid_number).unwrap_or(0),
        BodySpec::PlanetMoon => parse_int_arg(&args.planet_moon).unwrap_or(0),
        BodySpec::Fictitious => parse_int_arg(&args.fictitious).unwrap_or(0),
        _ => 0,
    };

    let ctx = FormatContext {
        name: calc_name,
        ipl,
        body,
        jd: if info.is_ut { tjd_ut } else { tjd_tt },
        tjd_ut,
        tjd_tt,
        year: info.year,
        month: info.month,
        day: info.day,
        hour: info.hour,
        cal: info.cal,
        is_ut: info.is_ut,
        data,
        xequ,
        xaz,
        xcart,
        xecart,
        hpos,
        hposj,
        armc,
        attr,
        args,
        is_label: false,
        is_house: false,
    };

    println!("{}", format::format_line(&ctx, eph));
}

pub fn run(args: &SweTestArgs, eph: &Ephemeris) {
    let iflag = args.build_iflag();
    let config = eph.config();
    let start = resolve_start_jd(args, config);
    let needs = format::scan_format_needs(&args.format);

    let step_count = if args.step_count == 0 {
        20
    } else {
        args.step_count
    };

    for istep in 1..=step_count {
        let tjd_step = step_jd(
            args,
            istep,
            start.year,
            start.month,
            start.day,
            start.hour,
            if start.is_ut {
                start.tjd_ut
            } else {
                start.tjd_tt
            },
        );

        let (tjd_ut, tjd_tt) = match args.time_mode {
            TimeMode::ET => {
                let dt = swisseph::deltat::calc_deltat(tjd_step, config);
                (tjd_step - dt, tjd_step)
            }
            TimeMode::UT | TimeMode::UTC => {
                let dt = swisseph::deltat::calc_deltat(tjd_step, config);
                (tjd_step, tjd_step + dt)
            }
            TimeMode::LMT => {
                let dt = swisseph::deltat::calc_deltat(tjd_step, config);
                (tjd_step, tjd_step + dt)
            }
            TimeMode::LAT => {
                let dt = swisseph::deltat::calc_deltat(tjd_step, config);
                (tjd_step, tjd_step + dt)
            }
        };

        let display_jd = if start.is_ut { tjd_ut } else { tjd_tt };
        let cal = calendar_for_jd(display_jd);
        let (y, m, d, h) = swisseph::date::revjul(display_jd, cal);

        let info = EpochInfo {
            tjd_ut,
            tjd_tt,
            year: y,
            month: m,
            day: d,
            hour: h,
            cal,
            is_ut: start.is_ut,
        };

        if istep == 1 || args.with_header_always {
            print_header(args, eph, &info, iflag);
        }

        if args.do_ayanamsa && !args.sidereal {
            match eph.get_ayanamsa_ex(tjd_tt, iflag) {
                Ok(aya) => {
                    let name = format!("Ayanamsha {}", sidereal_mode_name(args.sid_mode));
                    let mut data = [0.0_f64; 6];
                    data[0] = aya;
                    let ctx = FormatContext {
                        name,
                        ipl: -1,
                        body: None,
                        jd: if info.is_ut { tjd_ut } else { tjd_tt },
                        tjd_ut,
                        tjd_tt,
                        year: info.year,
                        month: info.month,
                        day: info.day,
                        hour: info.hour,
                        cal: info.cal,
                        is_ut: info.is_ut,
                        data,
                        xequ: None,
                        xaz: None,
                        xcart: None,
                        xecart: None,
                        hpos: None,
                        hposj: None,
                        armc: None,
                        attr: None,
                        args,
                        is_label: false,
                        is_house: true,
                    };
                    println!("{}", format::format_line(&ctx, eph));
                }
                Err(e) => println!("Ayanamsha       error: {e}"),
            }
            continue;
        }

        let bodies = args.body_specs();
        for spec in &bodies {
            compute_body(eph, spec, args, tjd_tt, tjd_ut, iflag, &needs, &info);
        }

        if istep < step_count {
            println!();
        }
    }
}
