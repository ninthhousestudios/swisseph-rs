use swisseph::Ephemeris;
use swisseph::flags::{CalcFlags, EclipseFlags, HeliacalFlags, RiseSetFlags};
use swisseph::types::{Body, CalendarType};

use crate::args::{BodySpec, EclipseFilters, SpecialEvent, SweTestArgs};
use crate::compute;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn hms(hours: f64, leading_zeros: bool) -> String {
    let x = hours.abs() + 0.5 / 36000.0;
    let h = x as i32 % 24;
    let mf = (x - x.floor()) * 60.0;
    let m = mf as i32;
    let sf = (mf - m as f64) * 60.0;
    let si = sf as i32;
    let frac = ((sf - si as f64) * 10.0) as i32;
    let sign = if hours < 0.0 { "-" } else { " " };
    if leading_zeros {
        format!("{sign} {h:02}:{m:02}:{si:02}.{frac}")
    } else {
        format!("{sign}{h:>3}:{m:02}:{si:02}.{frac}")
    }
}

fn hms_from_tjd(tjd: f64) -> String {
    let mut x = (tjd + 0.5) % 1.0;
    if x < 0.0 {
        x += 1.0;
    }
    format!("{} ", hms(x * 24.0, true))
}

fn get_gregjul(cal: CalendarType, year: i32) -> &'static str {
    match cal {
        CalendarType::Julian => " jul",
        CalendarType::Gregorian if year < 1700 => " greg",
        _ => "",
    }
}

fn jd_to_date(jd: f64) -> (i32, i32, i32, f64, CalendarType) {
    let cal = compute::calendar_for_jd(jd);
    let (y, m, d, h) = swisseph::date::revjul(jd, cal);
    (y, m, d, h, cal)
}

fn eclipse_ifltype(f: &EclipseFilters) -> EclipseFlags {
    let mut flags = EclipseFlags::empty();
    if f.total {
        flags |= EclipseFlags::TOTAL;
    }
    if f.annular {
        flags |= EclipseFlags::ANNULAR;
    }
    if f.annular_total {
        flags |= EclipseFlags::HYBRID;
    }
    if f.partial {
        flags |= EclipseFlags::PARTIAL;
    }
    if f.penumbral {
        flags |= EclipseFlags::PENUMBRAL;
    }
    if f.central {
        flags |= EclipseFlags::CENTRAL;
    }
    if f.noncentral {
        flags |= EclipseFlags::NONCENTRAL;
    }
    flags
}

fn geopos(args: &SweTestArgs) -> [f64; 3] {
    [args.geo_longitude, args.geo_latitude, args.geo_elevation]
}

fn build_epheflag(args: &SweTestArgs) -> CalcFlags {
    args.build_iflag() & (CalcFlags::SWIEPH | CalcFlags::JPLEPH | CalcFlags::MOSEPH)
}

fn saros_str(series: f64, member: f64) -> String {
    format!("{}/{}", series as i32, member as i32)
}

fn eclipse_type_str(flags: EclipseFlags) -> &'static str {
    if flags.contains(EclipseFlags::TOTAL) {
        "total"
    } else if flags.contains(EclipseFlags::HYBRID) {
        "ann-tot"
    } else if flags.contains(EclipseFlags::ANNULAR) {
        "annular"
    } else if flags.contains(EclipseFlags::PARTIAL) {
        "partial"
    } else if flags.contains(EclipseFlags::PENUMBRAL) {
        "penumbral"
    } else {
        ""
    }
}

fn dms_round_sec(deg: f64) -> String {
    let negative = deg < 0.0;
    let d = deg.abs() + 0.5 / 3600.0;
    let dd = d as i32;
    let mf = (d - dd as f64) * 60.0;
    let mm = mf as i32;
    let sf = (mf - mm as f64) * 60.0;
    let ss = sf as i32;
    let sdeg = if negative {
        format!("{:>4}", -(dd as i32))
    } else {
        format!(" {:>3}", dd)
    };
    format!("{sdeg}\u{00b0}{mm:02}'{ss:>2}\"")
}

fn resolve_first_body(args: &SweTestArgs) -> (Option<Body>, Option<String>) {
    let specs = args.body_specs();
    if let Some(spec) = specs.first() {
        match spec {
            BodySpec::FixedStar => (None, args.star_name.clone()),
            _ => (compute::resolve_body(spec, args), None),
        }
    } else {
        (Some(Body::Sun), None)
    }
}

fn apply_gap(s: &str, args: &SweTestArgs) -> String {
    if args.have_gap_parameter {
        s.replace('\t', &args.gap)
    } else {
        s.to_owned()
    }
}

fn do_print(s: &str, args: &SweTestArgs) {
    print!("{}", apply_gap(s, args));
}

// ---------------------------------------------------------------------------
// Rise/set + meridian transit
// ---------------------------------------------------------------------------

fn rise_set_flags(args: &SweTestArgs) -> RiseSetFlags {
    let mut rsmi = RiseSetFlags::empty();
    if args.no_refrac {
        rsmi |= RiseSetFlags::NO_REFRACTION;
    }
    if args.disc_center {
        rsmi |= RiseSetFlags::DISC_CENTER;
    }
    if args.disc_bottom {
        rsmi |= RiseSetFlags::DISC_BOTTOM;
    }
    if args.hindu {
        rsmi |= RiseSetFlags::HINDU_RISING;
    }
    rsmi
}

fn call_rise_set(args: &SweTestArgs, eph: &Ephemeris, start_ut: f64) {
    let epheflag = build_epheflag(args);
    let rsmi_base = rise_set_flags(args);
    let geo = geopos(args);
    let (body, starname) = resolve_first_body(args);
    let body = body.unwrap_or(Body::Sun);

    let nstep = if args.step_count == 0 {
        1
    } else {
        args.step_count
    };

    if args.with_header {
        do_print(
            &format!(
                "\ngeo. long {:.6}, lat {:.6}, alt {:.6}\n\n",
                args.geo_longitude, args.geo_latitude, args.geo_elevation
            ),
            args,
        );
    }

    let mut t_ut = start_ut;
    let end_ut = start_ut + nstep as f64;
    let mut last_was_empty = false;

    while t_ut <= end_ut {
        let rise_rsmi = RiseSetFlags::RISE | rsmi_base;
        let set_rsmi = RiseSetFlags::SET | rsmi_base;

        let rise_result = eph.rise_trans(
            t_ut,
            body,
            starname.as_deref(),
            epheflag,
            rise_rsmi,
            geo,
            args.atmosphere[0],
            args.atmosphere[1],
        );
        let set_result = eph.rise_trans(
            t_ut,
            body,
            starname.as_deref(),
            epheflag,
            set_rsmi,
            geo,
            args.atmosphere[0],
            args.atmosphere[1],
        );

        let rise_jd = rise_result.ok().map(|r| r.time);
        let set_jd = set_result.ok().map(|r| r.time);

        if rise_jd.is_none() && set_jd.is_none() {
            if !last_was_empty {
                do_print("-  -\n", args);
            }
            last_was_empty = true;
            t_ut += 1.0;
            continue;
        }
        last_was_empty = false;

        let mut line = String::new();
        line.push_str("rise     ");
        if let Some(rjd) = rise_jd {
            let (y, m, d, jut, _) = jd_to_date(rjd);
            line.push_str(&format!("{d:>2}.{m:02}.{y:04}\t{}    ", hms(jut, true)));
        } else {
            line.push_str("         -\t           -    ");
        }
        line.push_str("set      ");
        if let Some(sjd) = set_jd {
            let (y, m, d, jut, _) = jd_to_date(sjd);
            line.push_str(&format!("{d:>2}.{m:02}.{y:04}\t{}    ", hms(jut, true)));
        } else {
            line.push_str("         -\t           -    ");
        }
        if let (Some(rjd), Some(sjd)) = (rise_jd, set_jd) {
            let dt_hours = (sjd - rjd) * 24.0;
            line.push_str(&format!("dt ={}", hms(dt_hours, true)));
        }
        line.push('\n');
        do_print(&line, args);

        t_ut += 1.0;
    }
}

fn call_meridian_transit(args: &SweTestArgs, eph: &Ephemeris, start_ut: f64) {
    let epheflag = build_epheflag(args);
    let geo = geopos(args);
    let (body, starname) = resolve_first_body(args);
    let body = body.unwrap_or(Body::Sun);

    let nstep = if args.step_count == 0 {
        1
    } else {
        args.step_count
    };

    if args.with_header {
        do_print(
            &format!(
                "\ngeo. long {:.6}, lat {:.6}, alt {:.6}\n\n",
                args.geo_longitude, args.geo_latitude, args.geo_elevation
            ),
            args,
        );
    }

    let mut t_ut = start_ut;

    for _ in 0..nstep {
        let mt_result = eph.rise_trans(
            t_ut,
            body,
            starname.as_deref(),
            epheflag,
            RiseSetFlags::MTRANSIT,
            geo,
            args.atmosphere[0],
            args.atmosphere[1],
        );
        let it_search_t = mt_result.as_ref().map(|r| r.time).unwrap_or(t_ut);
        let it_result = eph.rise_trans(
            it_search_t,
            body,
            starname.as_deref(),
            epheflag,
            RiseSetFlags::ITRANSIT,
            geo,
            args.atmosphere[0],
            args.atmosphere[1],
        );

        let mut line = String::new();
        line.push_str("mtransit ");
        if let Ok(r) = &mt_result {
            let (y, m, d, jut, _) = jd_to_date(r.time);
            line.push_str(&format!("{d:>2}.{m:02}.{y:04}\t{}    ", hms(jut, true)));
        } else {
            line.push_str("         -\t           -    ");
        }
        line.push_str("itransit ");
        if let Ok(r) = &it_result {
            let (y, m, d, jut, _) = jd_to_date(r.time);
            line.push_str(&format!("{d:>2}.{m:02}.{y:04}\t{}    ", hms(jut, true)));
        } else {
            line.push_str("         -\t           -    ");
        }
        line.push('\n');
        do_print(&line, args);

        if let Ok(r) = mt_result {
            t_ut = r.time + 0.001;
        } else {
            t_ut += 1.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Solar eclipse
// ---------------------------------------------------------------------------

fn call_solar_eclipse(args: &SweTestArgs, eph: &Ephemeris, start_ut: f64) {
    let ifl = build_epheflag(args);
    let ifltype = eclipse_ifltype(&args.eclipse_filters);
    let geo = geopos(args);
    let nstep = if args.step_count == 0 {
        1
    } else {
        args.step_count
    };

    if args.eclipse_filters.how {
        call_solar_eclipse_how(args, eph, start_ut, ifl, &geo);
        return;
    }

    if args.eclipse_filters.local {
        if args.with_header {
            do_print(
                &format!(
                    "\ngeo. long {:.6}, lat {:.6}, alt {:.6}\n\n",
                    args.geo_longitude, args.geo_latitude, args.geo_elevation
                ),
                args,
            );
        }
        call_solar_eclipse_local(args, eph, start_ut, ifl, ifltype, &geo, nstep);
    } else {
        do_print("\n", args);
        call_solar_eclipse_global(args, eph, start_ut, ifl, ifltype, nstep);
    }
}

fn call_solar_eclipse_how(
    args: &SweTestArgs,
    eph: &Ephemeris,
    tjd_ut: f64,
    ifl: CalcFlags,
    geo: &[f64; 3],
) {
    match eph.sol_eclipse_how(tjd_ut, ifl, *geo) {
        Ok(how) => {
            let typ = eclipse_type_str(how.flags);
            if typ.is_empty() {
                do_print("no solar eclipse\n", args);
            } else {
                do_print(
                    &format!("{typ} solar eclipse: {:.6}\n", how.magnitude),
                    args,
                );
            }
        }
        Err(e) => eprintln!("error: {e}"),
    }
}

fn call_solar_eclipse_local(
    args: &SweTestArgs,
    eph: &Ephemeris,
    start_ut: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    geo: &[f64; 3],
    nstep: i32,
) {
    let mut t_ut = start_ut;
    let mut ii = 0;

    while ii < nstep {
        let result = match eph.sol_eclipse_when_loc(t_ut, ifl, *geo, args.backward) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: {e}");
                return;
            }
        };

        let ecl = result.attr.flags;
        if !ifltype.is_empty() && (ecl & ifltype).is_empty() {
            t_ut = result.time_maximum + (if args.backward { -10.0 } else { 10.0 });
            continue;
        }

        let typ = if ecl.contains(EclipseFlags::TOTAL) {
            "total   "
        } else if ecl.contains(EclipseFlags::ANNULAR) {
            "annular "
        } else {
            "partial "
        };

        let tmax = result.time_maximum;
        let (y, m, d, jut, cal) = jd_to_date(tmax);
        let sgj = get_gregjul(cal, y);
        let saros = saros_str(result.attr.saros_series, result.attr.saros_member);

        let mut s = format!(
            "{typ}solar eclipse\t{d:>2}.{m:02}.{y:04}{sgj}\t{}\t{:.4}/{:.4}/{:.4}\tsaros {saros}\t{tmax:.6}\n",
            hms(jut, true),
            result.attr.nasa_magnitude,
            result.attr.magnitude,
            result.attr.obscuration,
        );

        // Contacts
        let mut cline = String::from("\t");
        for &(t, _vis) in &[
            (result.time_first_contact, EclipseFlags::VISIBLE),
            (result.time_second_contact, EclipseFlags::VISIBLE),
            (result.time_third_contact, EclipseFlags::VISIBLE),
            (result.time_fourth_contact, EclipseFlags::VISIBLE),
        ] {
            if t != 0.0 {
                cline.push_str(&format!("{} ", hms_from_tjd(t)));
            } else {
                cline.push_str("   -         ");
            }
        }
        let dt = swisseph::deltat::calc_deltat(tmax, eph.config()) * 86400.0;
        let cline_trimmed = cline.trim_end();
        s.push_str(cline_trimmed);
        s.push_str(&format!(" dt={dt:.1}\n"));
        do_print(&s, args);

        t_ut = tmax + (if args.backward { -10.0 } else { 10.0 });
        ii += 1;
    }
}

fn call_solar_eclipse_global(
    args: &SweTestArgs,
    eph: &Ephemeris,
    start_ut: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    nstep: i32,
) {
    let mut t_ut = start_ut;

    for _ in 0..nstep {
        let result = match eph.sol_eclipse_when_glob(t_ut, ifl, ifltype, args.backward) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: {e}");
                return;
            }
        };

        let tmax = result.time_maximum;
        let ecl = result.flags;

        let where_result = eph.sol_eclipse_where(tmax, ifl);
        let how_result = if let Ok(ref w) = where_result {
            eph.sol_eclipse_how(tmax, ifl, [w.central_longitude, w.central_latitude, 0.0])
                .ok()
        } else {
            None
        };

        let mut typ = eclipse_type_str(ecl).to_owned();
        if ecl.contains(EclipseFlags::NONCENTRAL) && !ecl.contains(EclipseFlags::PARTIAL) {
            typ.push_str(" non-central");
        }

        let (y, m, d, jut, cal) = jd_to_date(tmax);
        let sgj = get_gregjul(cal, y);
        let attr = how_result.as_ref();
        let path_km = attr.map(|a| a.core_diameter_km).unwrap_or(0.0);
        let nasa_mag = attr.map(|a| a.nasa_magnitude).unwrap_or(0.0);
        let mag = attr.map(|a| a.magnitude).unwrap_or(0.0);
        let obsc = attr.map(|a| a.obscuration).unwrap_or(0.0);
        let saros = attr
            .map(|a| saros_str(a.saros_series, a.saros_member))
            .unwrap_or_default();

        let mut s = format!(
            "{typ} solar\t{d:>2}.{m:02}.{y:04}{sgj}\t{}\t{path_km:.6} km\t{nasa_mag:.4}/{mag:.4}/{obsc:.4}\tsaros {saros}\t{tmax:.6}\n",
            hms(jut, false),
        );

        // Contact times
        let mut cline = String::from("\t");
        for &t in &[
            result.time_begin,
            result.time_totality_begin,
            result.time_totality_end,
            result.time_end,
        ] {
            if t != 0.0 {
                cline.push_str(&format!("{} ", hms_from_tjd(t)));
            } else {
                cline.push_str("   -         ");
            }
        }
        let dt = swisseph::deltat::calc_deltat(tmax, eph.config()) * 86400.0;
        let cline_trimmed = cline.trim_end();
        s.push_str(cline_trimmed);
        s.push_str(&format!(" dt={dt:.1}\n"));

        // Geographic coordinates of maximum
        if let Ok(ref w) = where_result {
            s.push_str(&format!(
                "\t{}\t{}\n",
                dms_round_sec(w.central_longitude),
                dms_round_sec(w.central_latitude),
            ));
        }

        do_print(&s, args);

        t_ut = tmax + (if args.backward { -10.0 } else { 10.0 });
    }
}

// ---------------------------------------------------------------------------
// Lunar eclipse
// ---------------------------------------------------------------------------

fn call_lunar_eclipse(args: &SweTestArgs, eph: &Ephemeris, start_ut: f64) {
    let ifl = build_epheflag(args);
    let ifltype = eclipse_ifltype(&args.eclipse_filters);
    let geo = geopos(args);
    let nstep = if args.step_count == 0 {
        1
    } else {
        args.step_count
    };

    if args.eclipse_filters.how {
        match eph.lun_eclipse_how(start_ut, ifl, geo) {
            Ok(how) => {
                let typ = eclipse_type_str(how.flags);
                if typ.is_empty() {
                    do_print("no lunar eclipse\n", args);
                } else {
                    do_print(
                        &format!("{typ} lunar eclipse: {:.6} o/o\n", how.umbral_magnitude),
                        args,
                    );
                }
            }
            Err(e) => eprintln!("error: {e}"),
        }
        return;
    }

    if args.eclipse_filters.local {
        if args.with_header {
            do_print(
                &format!(
                    "\ngeo. long {:.6}, lat {:.6}, alt {:.6}\n\n",
                    args.geo_longitude, args.geo_latitude, args.geo_elevation
                ),
                args,
            );
        }
        call_lunar_eclipse_local(args, eph, start_ut, ifl, ifltype, &geo, nstep);
    } else {
        do_print("\n", args);
        call_lunar_eclipse_global(args, eph, start_ut, ifl, ifltype, &geo, nstep);
    }
}

fn call_lunar_eclipse_local(
    args: &SweTestArgs,
    eph: &Ephemeris,
    start_ut: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    geo: &[f64; 3],
    nstep: i32,
) {
    let mut t_ut = start_ut;

    for _ in 0..nstep {
        let result = match eph.lun_eclipse_when_loc(t_ut, ifl, *geo, args.backward) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: {e}");
                return;
            }
        };

        let ecl = result.attr.flags;
        if !ifltype.is_empty() && (ecl & ifltype).is_empty() {
            t_ut = result.time_maximum + (if args.backward { -25.0 } else { 25.0 });
            continue;
        }

        let typ = if ecl.contains(EclipseFlags::TOTAL) {
            "total   "
        } else if ecl.contains(EclipseFlags::PENUMBRAL) {
            "penumb. "
        } else {
            "partial "
        };

        let tmax = result.time_maximum;
        let (y, m, d, jut, cal) = jd_to_date(tmax);
        let sgj = get_gregjul(cal, y);
        let saros = saros_str(result.attr.saros_series, result.attr.saros_member);

        let mut s = format!(
            "{typ}lunar eclipse\t{d:>2}.{m:02}.{y:04}{sgj}\t{}\t{:.4}/{:.4}\tsaros {saros}\t{tmax:.6}\n",
            hms(jut, true),
            result.attr.umbral_magnitude,
            result.attr.penumbral_magnitude,
        );

        // Visibility-gated phase times
        let mut cline = String::from("\t");
        let contacts = [
            (result.time_penumbral_begin, EclipseFlags::PENUMBBEG_VISIBLE),
            (result.time_partial_begin, EclipseFlags::PARTBEG_VISIBLE),
            (result.time_totality_begin, EclipseFlags::TOTBEG_VISIBLE),
            (result.time_totality_end, EclipseFlags::TOTEND_VISIBLE),
            (result.time_partial_end, EclipseFlags::PARTEND_VISIBLE),
        ];
        for (t, vis_flag) in contacts {
            if t != 0.0 && ecl.contains(vis_flag) {
                cline.push_str(&format!("  {} ", hms_from_tjd(t)));
            } else {
                cline.push_str("      -         ");
            }
        }
        let dt = swisseph::deltat::calc_deltat(tmax, eph.config()) * 86400.0;
        let cline_trimmed = cline.trim_end();
        s.push_str(cline_trimmed);
        s.push_str(&format!(" dt={dt:.1}\n"));
        do_print(&s, args);

        t_ut = tmax + (if args.backward { -25.0 } else { 25.0 });
    }
}

fn call_lunar_eclipse_global(
    args: &SweTestArgs,
    eph: &Ephemeris,
    start_ut: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    geo: &[f64; 3],
    nstep: i32,
) {
    let mut t_ut = start_ut;

    for _ in 0..nstep {
        let result = match eph.lun_eclipse_when(t_ut, ifl, ifltype, args.backward) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: {e}");
                return;
            }
        };

        let tmax = result.time_maximum;
        let ecl = result.flags;

        let how = eph.lun_eclipse_how(tmax, ifl, *geo).ok();

        let typ = if ecl.contains(EclipseFlags::TOTAL) {
            "total "
        } else if ecl.contains(EclipseFlags::PENUMBRAL) {
            "penumb. "
        } else {
            "partial "
        };

        let (y, m, d, jut, cal) = jd_to_date(tmax);
        let sgj = get_gregjul(cal, y);
        let attr_mag = how.as_ref().map(|h| h.umbral_magnitude).unwrap_or(0.0);
        let attr_dr = how.as_ref().map(|h| h.penumbral_magnitude).unwrap_or(0.0);
        let saros = how
            .as_ref()
            .map(|h| saros_str(h.saros_series, h.saros_member))
            .unwrap_or_default();

        let mut s = format!(
            "{typ}lunar eclipse\t{d:>2}.{m:02}.{y:04}{sgj}\t{}\t{attr_mag:.4}/{attr_dr:.4}\tsaros {saros}\t{tmax:.6}\n",
            hms(jut, true),
        );

        // Contact times (ungated for global)
        let mut cline = String::from("\t");
        let contacts = [
            result.time_penumbral_begin,
            result.time_partial_begin,
            result.time_totality_begin,
            result.time_totality_end,
            result.time_partial_end,
            result.time_penumbral_end,
        ];
        for t in contacts {
            if t != 0.0 {
                cline.push_str(&format!("{} ", hms_from_tjd(t)));
            } else {
                cline.push_str("   -         ");
            }
        }
        let dt = swisseph::deltat::calc_deltat(tmax, eph.config()) * 86400.0;
        let cline_trimmed = cline.trim_end();
        s.push_str(cline_trimmed);
        s.push_str(&format!(" dt={dt:.1}\n"));
        do_print(&s, args);

        t_ut = tmax + (if args.backward { -1.0 } else { 1.0 });
    }
}

// ---------------------------------------------------------------------------
// Occultation
// ---------------------------------------------------------------------------

fn call_occultation(args: &SweTestArgs, eph: &Ephemeris, start_ut: f64) {
    let ifl = build_epheflag(args);
    let ifltype = eclipse_ifltype(&args.eclipse_filters);
    let geo = geopos(args);
    let nstep = if args.step_count == 0 {
        1
    } else {
        args.step_count
    };

    let (mut body, starname) = resolve_first_body(args);
    if body == Some(Body::Moon) {
        body = Some(Body::Mercury);
    }
    let body = body.unwrap_or(Body::Mercury);

    do_print("\n", args);

    if args.eclipse_filters.local {
        if args.with_header {
            do_print(
                &format!(
                    "geo. long {:.6}, lat {:.6}, alt {:.6}\n\n",
                    args.geo_longitude, args.geo_latitude, args.geo_elevation
                ),
                args,
            );
        }
        call_occultation_local(
            args,
            eph,
            start_ut,
            ifl,
            ifltype,
            body,
            starname.as_deref(),
            &geo,
            nstep,
        );
    } else {
        call_occultation_global(
            args,
            eph,
            start_ut,
            ifl,
            ifltype,
            body,
            starname.as_deref(),
            nstep,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn call_occultation_local(
    args: &SweTestArgs,
    eph: &Ephemeris,
    start_ut: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    body: Body,
    starname: Option<&str>,
    geo: &[f64; 3],
    nstep: i32,
) {
    let mut t_ut = start_ut;
    let mut ii = 0;

    while ii < nstep {
        let result = match eph.lun_occult_when_loc(t_ut, body, starname, ifl, *geo, args.backward) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: {e}");
                return;
            }
        };

        let ecl = result.flags;
        if ecl.is_empty() {
            t_ut += if args.backward { -10.0 } else { 10.0 };
            continue;
        }

        if !ifltype.is_empty() && (ecl & ifltype).is_empty() {
            t_ut = result.time_maximum + (if args.backward { -10.0 } else { 10.0 });
            continue;
        }

        let typ = eclipse_type_str(ecl);
        let tmax = result.time_maximum;
        let (y, m, d, jut, _cal) = jd_to_date(tmax);

        let mut s = format!(
            "{typ:<17}{d:>2}.{m:02}.{y:04}\t{}\t{:.6}\t{tmax:.6}\n",
            hms(jut, true),
            result.attr.magnitude,
        );

        // Contacts
        let mut cline = String::from("\t");
        for t in [
            result.time_first_contact,
            result.time_second_contact,
            result.time_third_contact,
            result.time_fourth_contact,
        ] {
            if t != 0.0 {
                cline.push_str(&format!("{} ", hms_from_tjd(t)));
            } else {
                cline.push_str("   -         ");
            }
        }
        let dt = swisseph::deltat::calc_deltat(tmax, eph.config()) * 86400.0;
        let cline_trimmed = cline.trim_end();
        s.push_str(cline_trimmed);
        s.push_str(&format!(" dt={dt:.1}\n"));
        do_print(&s, args);

        t_ut = tmax + (if args.backward { -1.0 } else { 1.0 });
        ii += 1;
    }
}

#[allow(clippy::too_many_arguments)]
fn call_occultation_global(
    args: &SweTestArgs,
    eph: &Ephemeris,
    start_ut: f64,
    ifl: CalcFlags,
    ifltype: EclipseFlags,
    body: Body,
    starname: Option<&str>,
    nstep: i32,
) {
    let mut t_ut = start_ut;
    let mut ii = 0;

    while ii < nstep {
        let result =
            match eph.lun_occult_when_glob(t_ut, body, starname, ifl, ifltype, args.backward) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("error: {e}");
                    return;
                }
            };

        let ecl = result.flags;
        if ecl.is_empty() {
            t_ut += if args.backward { -1.0 } else { 1.0 };
            continue;
        }

        let tmax = result.time_maximum;

        let mut typ = eclipse_type_str(ecl).to_owned();
        if ecl.contains(EclipseFlags::NONCENTRAL) && !ecl.contains(EclipseFlags::PARTIAL) {
            typ.push_str(" non-central");
        }

        let where_result = eph.lun_occult_where(tmax, body, starname, ifl);

        let (y, m, d, jut, _cal) = jd_to_date(tmax);
        let path_km = where_result
            .as_ref()
            .ok()
            .map(|w| w.core_diameter_km)
            .unwrap_or(0.0);

        let mut s = format!(
            "{typ:<17}{d:>2}.{m:02}.{y:04}\t{}\t{path_km:.6} km\t{tmax:.6}\n",
            hms(jut, true),
        );

        // Contact times
        let mut cline = String::from("\t");
        for t in [
            result.time_begin,
            result.time_totality_begin,
            result.time_totality_end,
            result.time_end,
        ] {
            if t != 0.0 {
                cline.push_str(&format!("{} ", hms_from_tjd(t)));
            } else {
                cline.push_str("   -         ");
            }
        }
        let dt = swisseph::deltat::calc_deltat(tmax, eph.config()) * 86400.0;
        let cline_trimmed = cline.trim_end();
        s.push_str(cline_trimmed);
        s.push_str(&format!(" dt={dt:.1}\n"));

        // Geographic coordinates
        if let Ok(ref w) = where_result {
            s.push_str(&format!(
                "\t{}\t{}\n",
                dms_round_sec(w.central_longitude),
                dms_round_sec(w.central_latitude),
            ));
        }

        do_print(&s, args);

        t_ut = tmax + (if args.backward { -1.0 } else { 1.0 });
        ii += 1;
    }
}

// ---------------------------------------------------------------------------
// Heliacal events
// ---------------------------------------------------------------------------

fn call_heliacal(args: &SweTestArgs, eph: &Ephemeris, start_ut: f64) {
    use swisseph::heliacal::HeliacalEventType;

    let epheflag = build_epheflag(args);
    let geo = geopos(args);
    let nstep = if args.step_count == 0 {
        1
    } else {
        args.step_count
    };

    let (body, starname) = resolve_first_body(args);
    let object_name = if let Some(ref s) = starname {
        s.clone()
    } else if let Some(b) = body {
        eph.get_planet_name(b)
    } else {
        "Sun".to_owned()
    };

    let mut datm = args.atmosphere;
    let mut dobs = [0.0_f64; 6];
    dobs[..args.observer_params.len().min(6)]
        .copy_from_slice(&args.observer_params[..args.observer_params.len().min(6)]);

    let mut helflag = HeliacalFlags::from_bits_truncate(args.hel_flag as u32);
    if args.hel_using_av {
        helflag |= HeliacalFlags::AV;
    }

    if args.with_header {
        do_print(
            &format!(
                "\ngeo. long {:.6}, lat {:.6}, alt {:.6}\n\n",
                args.geo_longitude, args.geo_latitude, args.geo_elevation
            ),
            args,
        );
    }

    let event_type = match args.search_flag {
        1 => HeliacalEventType::MorningFirst,
        2 => HeliacalEventType::EveningLast,
        3 => HeliacalEventType::EveningFirst,
        4 => HeliacalEventType::MorningLast,
        _ => HeliacalEventType::MorningFirst,
    };

    let event_name = |et: HeliacalEventType| -> &'static str {
        match et {
            HeliacalEventType::MorningFirst => "heliacal rising ",
            HeliacalEventType::EveningLast => "heliacal setting",
            HeliacalEventType::EveningFirst => "evening first   ",
            HeliacalEventType::MorningLast => "morning last    ",
            HeliacalEventType::AcronymchalRising => "evening rising  ",
            HeliacalEventType::AcronymchalSetting => "morning setting ",
        }
    };

    let mut t_ut = start_ut;

    for _ in 0..nstep {
        let result = match eph.heliacal_ut(
            t_ut,
            &geo,
            &mut datm,
            &mut dobs,
            &object_name,
            event_type,
            epheflag,
            helflag,
        ) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: {e}");
                return;
            }
        };

        let tjd = result.start_visible;
        let (y, m, d, _, _cal) = jd_to_date(tjd);

        if helflag.contains(HeliacalFlags::AV) {
            do_print(
                &format!(
                    "{} {}: {}/{:02}/{:02} {} UT ({:.5})\n",
                    object_name,
                    event_name(event_type),
                    y,
                    m,
                    d,
                    hms_from_tjd(tjd).trim(),
                    tjd,
                ),
                args,
            );
        } else {
            let opt_time = if result.optimum_visibility != 0.0 {
                hms_from_tjd(result.optimum_visibility).trim().to_owned()
            } else {
                "-".to_owned()
            };
            let end_time = if result.end_visible != 0.0 {
                hms_from_tjd(result.end_visible).trim().to_owned()
            } else {
                "-".to_owned()
            };
            let dur = if result.end_visible != 0.0 && result.start_visible != 0.0 {
                (result.end_visible - result.start_visible) * 1440.0
            } else {
                0.0
            };
            do_print(
                &format!(
                    "{} {}: {}/{:02}/{:02} {} UT ({:.5}), opt {opt_time}, end {end_time}, dur {dur:.1} min\n",
                    object_name,
                    event_name(event_type),
                    y,
                    m,
                    d,
                    hms_from_tjd(tjd).trim(),
                    tjd,
                ),
                args,
            );
        }

        t_ut = tjd + 1.0;
    }
}

// ---------------------------------------------------------------------------
// Orbital elements
// ---------------------------------------------------------------------------

fn call_orbital_elements(args: &SweTestArgs, eph: &Ephemeris, tjd_tt: f64) {
    let iflag = args.build_iflag();
    let specs = args.body_specs();

    for spec in &specs {
        let body = match compute::resolve_body(spec, args) {
            Some(b) => b,
            None => continue,
        };

        match eph.get_orbital_elements(tjd_tt, body, iflag) {
            Ok(el) => {
                let peri_jd = el.perihelion_passage;
                let (py, pm, pd, pjut, _pcal) = jd_to_date(peri_jd);
                let sdateperi = format!("{pd:>2}.{pm:02}.{py:04},{}", hms(pjut, true),);

                println!("semiaxis         \t{:.6}", el.semi_major_axis);
                println!("eccentricity     \t{:.6}", el.eccentricity);
                println!("inclination      \t{:.6}", el.inclination);
                println!("asc. node       \t{:.6}", el.ascending_node);
                println!("arg. pericenter  \t{:.6}", el.arg_perihelion);
                println!("pericenter       \t{:.6}", el.perihelion_lon);
                println!("mean longitude   \t{:.6}", el.mean_longitude);
                println!("mean anomaly     \t{:.6}", el.mean_anomaly);
                println!("ecc. anomaly     \t{:.6}", el.eccentric_anomaly);
                println!("true anomaly     \t{:.6}", el.true_anomaly);
                println!(
                    "time pericenter  \t{:.6} {sdateperi}",
                    el.perihelion_passage
                );
                println!("dist. pericenter \t{:.6}", el.perihelion_distance);
                println!("dist. apocenter  \t{:.6}", el.aphelion_distance);
                println!("mean daily motion\t{:.6}", el.mean_daily_motion);
                println!("sid. period (y)  \t{:.6}", el.sidereal_period);
                println!("trop. period (y) \t{:.6}", el.tropical_period);
                println!("synodic cycle (d)\t{:.6}", el.synodic_period);
            }
            Err(e) => {
                let name = compute::body_name(eph, spec, args);
                println!("{name}: error: {e}");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Main dispatch
// ---------------------------------------------------------------------------

pub fn run(args: &SweTestArgs, eph: &Ephemeris) {
    let start = compute::resolve_start_jd(args, eph.config());

    if args.orbital_elements {
        call_orbital_elements(args, eph, start.tjd_tt);
        return;
    }

    if let Some(event) = args.special_event {
        match event {
            SpecialEvent::RiseSet => call_rise_set(args, eph, start.tjd_ut),
            SpecialEvent::MeridianTransit => call_meridian_transit(args, eph, start.tjd_ut),
            SpecialEvent::SolarEclipse => call_solar_eclipse(args, eph, start.tjd_ut),
            SpecialEvent::LunarEclipse => call_lunar_eclipse(args, eph, start.tjd_ut),
            SpecialEvent::Occultation => call_occultation(args, eph, start.tjd_ut),
            SpecialEvent::Heliacal => call_heliacal(args, eph, start.tjd_ut),
        }
    }
}
