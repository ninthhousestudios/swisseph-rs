use std::path::PathBuf;

use swisseph::flags::CalcFlags;
use swisseph::{Body, EphemerisConfig, EphemerisSource, TopoPosition};

const PLSEL_D: &str = "0123456789mtA";
const PLSEL_P: &str = "0123456789mtABCcgDEFGHI";
const PLSEL_H: &str = "JKLMNOPQRSTUVWXYZw";
const PLSEL_A: &str = "0123456789mtABCcgDEFGHIJKLMNOPQRSTUVWXYZw";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeMode {
    ET,
    UT,
    UTC,
    LMT,
    LAT,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepUnit {
    Days,
    Minutes,
    Seconds,
    Years,
    Months,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EphemerisChoice {
    Swiss,
    Moshier,
    Jpl,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecialEvent {
    LunarEclipse,
    SolarEclipse,
    Occultation,
    RiseSet,
    MeridianTransit,
    Heliacal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiffMode {
    None,
    Diff,
    DiffAbs,
    DiffHelio,
}

#[derive(Debug, Clone)]
pub struct SiderealArgs {
    pub mode: i32,
    pub user_t0: f64,
    pub user_ayan: f64,
    pub user_ut: bool,
}

#[derive(Debug, Clone)]
pub struct HouseArgs {
    pub longitude: f64,
    pub latitude: f64,
    pub system: char,
    pub hpos_method: i32,
}

#[derive(Debug, Clone)]
pub struct EclipseFilters {
    pub total: bool,
    pub annular: bool,
    pub annular_total: bool,
    pub partial: bool,
    pub penumbral: bool,
    pub central: bool,
    pub noncentral: bool,
    pub local: bool,
    pub how: bool,
    pub hocal: bool,
}

impl Default for EclipseFilters {
    fn default() -> Self {
        Self {
            total: false,
            annular: false,
            annular_total: false,
            partial: false,
            penumbral: false,
            central: false,
            noncentral: false,
            local: false,
            how: false,
            hocal: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SweTestArgs {
    // Time input
    pub begin_date: Option<String>,
    pub time_input: String,
    pub time_mode: TimeMode,

    // Body selection
    pub planet_selection: String,
    pub asteroid_number: Option<String>,
    pub star_name: Option<String>,
    pub planet_moon: Option<String>,
    pub fictitious: Option<String>,

    // Stepping
    pub step_count: i32,
    pub has_n: bool,
    pub step_size: f64,
    pub step_unit: StepUnit,

    // Observer
    pub geo_longitude: f64,
    pub geo_latitude: f64,
    pub geo_elevation: f64,
    pub have_geopos: bool,
    pub heliocentric: bool,
    pub barycentric: bool,
    pub planetocentric: Option<i32>,

    // Frame/flags
    pub speed: bool,
    pub speed3: bool,
    pub no_speed: bool,
    pub no_aberration: bool,
    pub no_deflection: bool,
    pub no_nutation: bool,
    pub truepos: bool,
    pub j2000: bool,
    pub icrs: bool,
    pub center_of_body: bool,
    pub topocentric: bool,
    pub equatorial: bool,
    pub xyz: bool,
    pub force_iflag: Option<i32>,

    // Ephemeris source
    pub ephemeris: EphemerisChoice,
    pub ephe_dir: Option<String>,
    pub jpl_file: String,

    // Sidereal
    pub sidereal: bool,
    pub do_ayanamsa: bool,
    pub sid_mode: i32,
    pub sid_user_t0: f64,
    pub sid_user_ayan: f64,
    pub sid_user_ut: bool,

    // Houses
    pub do_houses: bool,
    pub house_system: char,
    pub hpos_method: i32,

    // Output format
    pub format: String,
    pub gap: String,
    pub have_gap_parameter: bool,
    pub with_header: bool,
    pub with_header_always: bool,
    pub horizontal: bool,
    pub dms: bool,
    pub round_sec: bool,
    pub round_min: bool,
    pub extra_precision: bool,
    pub short_output: bool,

    // Special events
    pub special_event: Option<SpecialEvent>,
    pub search_flag: i32,
    pub hel_using_av: bool,

    // Eclipse filters
    pub eclipse_filters: EclipseFilters,

    // Rise/set options
    pub no_refrac: bool,
    pub disc_center: bool,
    pub disc_bottom: bool,
    pub hindu: bool,

    // Heliacal
    pub atmosphere: [f64; 4],
    pub observer_params: [f64; 6],

    // Differential
    pub diff_mode: DiffMode,
    pub diff_planet: char,

    // Misc
    pub backward: bool,
    pub direction: i32,
    pub tidal_acc: Option<f64>,
    pub astro_models: Option<String>,
    pub show_file_limit: bool,
    pub chart_link: bool,
    pub use_fixstar2: bool,
    pub orbital_elements: bool,
    pub nutation_output: bool,
    pub hel_flag: i32,

    // testaa
    pub testaa: Option<String>,

    // with_glp (Astrodienst internal)
    pub with_glp: bool,
}

impl Default for SweTestArgs {
    fn default() -> Self {
        Self {
            begin_date: None,
            time_input: String::new(),
            time_mode: TimeMode::ET,

            planet_selection: PLSEL_D.to_owned(),
            asteroid_number: None,
            star_name: None,
            planet_moon: None,
            fictitious: None,

            step_count: 1,
            has_n: false,
            step_size: 1.0,
            step_unit: StepUnit::Days,

            geo_longitude: 0.0,
            geo_latitude: 51.5,
            geo_elevation: 0.0,
            have_geopos: false,
            heliocentric: false,
            barycentric: false,
            planetocentric: None,

            speed: false,
            speed3: false,
            no_speed: false,
            no_aberration: false,
            no_deflection: false,
            no_nutation: false,
            truepos: false,
            j2000: false,
            icrs: false,
            center_of_body: false,
            topocentric: false,
            equatorial: false,
            xyz: false,
            force_iflag: None,

            ephemeris: EphemerisChoice::Swiss,
            ephe_dir: None,
            jpl_file: "de441.eph".to_owned(),

            sidereal: false,
            do_ayanamsa: false,
            sid_mode: 0, // SE_SIDM_FAGAN_BRADLEY
            sid_user_t0: 0.0,
            sid_user_ayan: 0.0,
            sid_user_ut: false,

            do_houses: false,
            house_system: 'P',
            hpos_method: 0,

            format: "PLBRS".to_owned(),
            gap: " ".to_owned(),
            have_gap_parameter: false,
            with_header: true,
            with_header_always: false,
            horizontal: false,
            dms: false,
            round_sec: false,
            round_min: false,
            extra_precision: false,
            short_output: false,

            special_event: None,
            search_flag: 0,
            hel_using_av: false,

            eclipse_filters: EclipseFilters::default(),

            no_refrac: false,
            disc_center: false,
            disc_bottom: false,
            hindu: false,

            atmosphere: [1013.25, 15.0, 40.0, 0.0],
            observer_params: [0.0; 6],

            diff_mode: DiffMode::None,
            diff_planet: '0',

            backward: false,
            direction: 1,
            tidal_acc: None,
            astro_models: None,
            show_file_limit: false,
            chart_link: false,
            use_fixstar2: false,
            orbital_elements: false,
            nutation_output: false,
            hel_flag: 0,

            testaa: None,

            with_glp: false,
        }
    }
}

fn parse_comma_f64s(s: &str, out: &mut [f64]) {
    let s = s.strip_prefix('[').unwrap_or(s);
    for (i, part) in s.split(',').enumerate() {
        if i >= out.len() {
            break;
        }
        if let Ok(v) = part.trim().parse::<f64>() {
            out[i] = v;
        }
    }
}

pub fn parse_args(args: &[String]) -> Result<SweTestArgs, String> {
    let mut a = SweTestArgs::default();
    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];

        if let Some(rest) = arg.strip_prefix("-utc") {
            a.time_mode = TimeMode::UTC;
            if !rest.is_empty() {
                a.time_input = rest.to_owned();
            }
        } else if let Some(rest) = arg.strip_prefix("-ut") {
            a.time_mode = TimeMode::UT;
            if !rest.is_empty() {
                a.time_input = rest.to_owned();
            }
        } else if arg.starts_with("-glp") {
            a.with_glp = true;
        } else if arg.starts_with("-hor") {
            a.horizontal = true;
        } else if arg.starts_with("-head") {
            a.with_header = false;
        } else if arg.starts_with("+head") {
            a.with_header_always = true;
        } else if arg == "-j2000" {
            a.j2000 = true;
        } else if arg == "-icrs" {
            a.icrs = true;
        } else if arg == "-cob" {
            a.center_of_body = true;
        } else if let Some(rest) = arg.strip_prefix("-ay") {
            a.do_ayanamsa = true;
            a.sid_mode = rest.parse::<i32>().unwrap_or(0);
        } else if let Some(rest) = arg.strip_prefix("-sidt0") {
            a.sidereal = true;
            let mode = rest.parse::<i32>().unwrap_or(0);
            a.sid_mode = if mode == 0 { 0 } else { mode };
            a.sid_mode |= 256; // SE_SIDBIT_ECL_T0
        } else if let Some(rest) = arg.strip_prefix("-sidsp") {
            a.sidereal = true;
            let mode = rest.parse::<i32>().unwrap_or(0);
            a.sid_mode = if mode == 0 { 0 } else { mode };
            a.sid_mode |= 512; // SE_SIDBIT_SSY_PLANE
        } else if let Some(rest) = arg.strip_prefix("-sidudef") {
            a.sidereal = true;
            a.sid_mode = 255; // SE_SIDM_USER
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            let parts: Vec<&str> = rest.splitn(3, ',').collect();
            if !parts.is_empty() {
                a.sid_user_t0 = parts[0].parse::<f64>().unwrap_or(0.0);
            }
            if parts.len() > 1 {
                a.sid_user_ayan = parts[1].parse::<f64>().unwrap_or(0.0);
            }
            if rest.contains("jdisut") {
                a.sid_user_ut = true;
            }
        } else if let Some(rest) = arg.strip_prefix("-sidbit") {
            a.sid_mode |= rest.parse::<i32>().unwrap_or(0);
        } else if let Some(rest) = arg.strip_prefix("-sid") {
            a.sidereal = true;
            a.sid_mode = rest.parse::<i32>().unwrap_or(0);
        } else if arg == "-jplhora" {
            // SEFLG_JPLHOR_APPROX — internal, store as force_iflag overlay
        } else if arg == "-jplhor" {
            // SEFLG_JPLHOR — internal
        } else if let Some(rest) = arg.strip_prefix("-j") {
            // -jNNNN → Julian day as begin_date (note: begin_date includes the 'j' prefix)
            a.begin_date = Some(format!("j{rest}"));
        } else if let Some(rest) = arg.strip_prefix("-ejpl") {
            a.ephemeris = EphemerisChoice::Jpl;
            if !rest.is_empty() {
                a.jpl_file = rest.to_owned();
            }
        } else if let Some(rest) = arg.strip_prefix("-edir") {
            if !rest.is_empty() {
                a.ephe_dir = Some(rest.to_owned());
            }
        } else if arg == "-eswe" {
            a.ephemeris = EphemerisChoice::Swiss;
        } else if arg == "-emos" {
            a.ephemeris = EphemerisChoice::Moshier;
        } else if let Some(rest) = arg.strip_prefix("-helflag") {
            a.hel_flag = rest.parse::<i32>().unwrap_or(0);
            if a.hel_flag >= 64 {
                // SE_HELFLAG_AV
                a.hel_using_av = true;
            }
        } else if arg == "-hel" {
            a.heliocentric = true;
        } else if arg == "-bary" {
            a.barycentric = true;
        } else if let Some(rest) = arg.strip_prefix("-house") {
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            let parts: Vec<&str> = rest.split(',').collect();
            if !parts.is_empty() {
                a.geo_longitude = parts[0].parse::<f64>().unwrap_or(0.0);
            }
            if parts.len() > 1 {
                a.geo_latitude = parts[1].parse::<f64>().unwrap_or(0.0);
            }
            if parts.len() > 2 {
                if let Some(c) = parts[2].chars().next() {
                    a.house_system = c;
                }
            }
            a.do_houses = true;
            a.have_geopos = true;
        } else if let Some(rest) = arg.strip_prefix("-hsy") {
            if let Some(c) = rest.chars().next() {
                a.house_system = c;
            }
            if rest.len() > 1 {
                a.hpos_method = rest[1..].parse::<i32>().unwrap_or(0);
            }
            a.have_geopos = true;
        } else if let Some(rest) = arg.strip_prefix("-topo") {
            a.topocentric = true;
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            let mut vals = [0.0_f64; 3];
            parse_comma_f64s(rest, &mut vals);
            a.geo_longitude = vals[0];
            a.geo_latitude = vals[1];
            a.geo_elevation = vals[2];
            a.have_geopos = true;
        } else if let Some(rest) = arg.strip_prefix("-geopos") {
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            let mut vals = [0.0_f64; 3];
            parse_comma_f64s(rest, &mut vals);
            a.geo_longitude = vals[0];
            a.geo_latitude = vals[1];
            a.geo_elevation = vals[2];
            a.have_geopos = true;
        } else if arg == "-true" {
            a.truepos = true;
        } else if arg == "-noaberr" {
            a.no_aberration = true;
        } else if arg == "-nodefl" {
            a.no_deflection = true;
        } else if arg == "-nonut" {
            a.no_nutation = true;
        } else if arg == "-speed3" {
            a.speed3 = true;
        } else if arg == "-speed" {
            a.speed = true;
        } else if arg == "-nospeed" {
            a.no_speed = true;
        } else if let Some(rest) = arg.strip_prefix("-testaa") {
            a.ephemeris = EphemerisChoice::Jpl;
            a.jpl_file = "de200.eph".to_owned();
            match rest {
                "95" => a.begin_date = Some("j2449975.5".to_owned()),
                "96" => a.begin_date = Some("j2450442.5".to_owned()),
                "97" => a.begin_date = Some("j2450482.5".to_owned()),
                _ => {}
            }
            a.format = "PADRu".to_owned();
            a.time_mode = TimeMode::ET;
            a.planet_selection = "3".to_owned();
            a.testaa = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-lmt") {
            a.time_mode = TimeMode::LMT;
            if !rest.is_empty() {
                a.time_input = rest.to_owned();
            }
        } else if arg == "-lat" {
            a.time_mode = TimeMode::LAT;
        } else if arg == "-lim" {
            a.show_file_limit = true;
        } else if arg == "-clink" {
            a.chart_link = true;
        } else if arg == "-lunecl" {
            a.special_event = Some(SpecialEvent::LunarEclipse);
        } else if arg == "-solecl" {
            a.special_event = Some(SpecialEvent::SolarEclipse);
            a.have_geopos = true;
        } else if arg == "-short" {
            a.short_output = true;
        } else if arg == "-occult" {
            a.special_event = Some(SpecialEvent::Occultation);
            a.have_geopos = true;
        } else if arg == "-ep" {
            a.extra_precision = true;
        } else if arg == "-hocal" {
            a.eclipse_filters.hocal = true;
        } else if arg == "-how" {
            a.eclipse_filters.how = true;
        } else if arg == "-total" {
            a.eclipse_filters.total = true;
        } else if arg == "-annular" {
            a.eclipse_filters.annular = true;
        } else if arg == "-anntot" {
            a.eclipse_filters.annular_total = true;
        } else if arg == "-partial" {
            a.eclipse_filters.partial = true;
        } else if arg == "-penumbral" {
            a.eclipse_filters.penumbral = true;
        } else if arg == "-noncentral" {
            a.eclipse_filters.central = false;
            a.eclipse_filters.noncentral = true;
        } else if arg == "-central" {
            a.eclipse_filters.noncentral = false;
            a.eclipse_filters.central = true;
        } else if arg == "-local" {
            a.eclipse_filters.local = true;
        } else if arg == "-rise" {
            a.special_event = Some(SpecialEvent::RiseSet);
            a.have_geopos = true;
        } else if arg == "-norefrac" {
            a.no_refrac = true;
        } else if arg == "-disccenter" {
            a.disc_center = true;
        } else if arg == "-hindu" {
            a.hindu = true;
            a.no_refrac = true;
            a.disc_center = true;
        } else if arg == "-discbottom" {
            a.disc_bottom = true;
        } else if arg == "-metr" {
            a.special_event = Some(SpecialEvent::MeridianTransit);
            a.have_geopos = true;
        } else if let Some(rest) = arg.strip_prefix("-amod") {
            a.astro_models = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-tidacc") {
            a.tidal_acc = rest.parse::<f64>().ok();
        } else if let Some(rest) = arg.strip_prefix("-hev") {
            a.special_event = Some(SpecialEvent::Heliacal);
            a.search_flag = 0;
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            if !rest.is_empty() {
                // Extract numeric part (before any non-numeric suffix like "AV")
                let numeric: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
                if !numeric.is_empty() {
                    a.search_flag = numeric.parse::<i32>().unwrap_or(0);
                }
            }
            if arg.contains("AV") {
                a.hel_using_av = true;
            }
            a.have_geopos = true;
        } else if let Some(rest) = arg.strip_prefix("-at") {
            parse_comma_f64s(rest, &mut a.atmosphere);
        } else if let Some(rest) = arg.strip_prefix("-obs") {
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            let mut vals = [0.0_f64; 2];
            parse_comma_f64s(rest, &mut vals);
            a.observer_params[0] = vals[0];
            a.observer_params[1] = vals[1];
        } else if let Some(rest) = arg.strip_prefix("-opt") {
            let rest = rest.strip_prefix('[').unwrap_or(rest);
            parse_comma_f64s(rest, &mut a.observer_params);
        } else if arg == "-orbel" {
            a.orbital_elements = true;
        } else if arg == "-bwd" {
            a.backward = true;
            a.direction = -1;
        } else if let Some(rest) = arg.strip_prefix("-pc") {
            a.planetocentric = Some(rest.parse::<i32>().unwrap_or(0));
        } else if let Some(rest) = arg.strip_prefix("-p") {
            match rest.chars().next() {
                Some('d') => a.planet_selection = PLSEL_D.to_owned(),
                Some('p') => a.planet_selection = PLSEL_P.to_owned(),
                Some('h') => a.planet_selection = PLSEL_H.to_owned(),
                Some('a') => a.planet_selection = PLSEL_A.to_owned(),
                _ => a.planet_selection = rest.to_owned(),
            }
        } else if let Some(rest) = arg.strip_prefix("-xs") {
            a.asteroid_number = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-xv") {
            a.planet_moon = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-xf") {
            a.star_name = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-xz") {
            a.fictitious = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-x") {
            a.star_name = Some(rest.to_owned());
        } else if arg == "-nut" {
            a.nutation_output = true;
        } else if let Some(rest) = arg.strip_prefix("-n") {
            let n = rest.parse::<i32>().unwrap_or(0);
            a.step_count = if n == 0 { 20 } else { n };
            a.has_n = true;
        } else if let Some(rest) = arg.strip_prefix("-i") {
            let iflag_f = rest.parse::<i32>().unwrap_or(0);
            a.force_iflag = Some(iflag_f);
            if iflag_f & 0x200 != 0 {
                // SEFLG_XYZ
                a.format = "PX".to_owned();
            }
        } else if arg == "-swefixstar2" {
            a.use_fixstar2 = true;
        } else if let Some(rest) = arg.strip_prefix("-s") {
            let last_char = rest.chars().last();
            let numeric_part = match last_char {
                Some('m' | 's' | 'y' | 'o') => &rest[..rest.len() - 1],
                _ => rest,
            };
            a.step_size = numeric_part.parse::<f64>().unwrap_or(1.0);
            match last_char {
                Some('m') => a.step_unit = StepUnit::Minutes,
                Some('s') => a.step_unit = StepUnit::Seconds,
                Some('y') => a.step_unit = StepUnit::Years,
                Some('o') => a.step_unit = StepUnit::Months,
                _ => a.step_unit = StepUnit::Days,
            }
        } else if let Some(rest) = arg.strip_prefix("-b") {
            a.begin_date = Some(rest.to_owned());
        } else if let Some(rest) = arg.strip_prefix("-f") {
            a.format = rest.to_owned();
        } else if let Some(rest) = arg.strip_prefix("-g") {
            a.have_gap_parameter = true;
            if rest.is_empty() {
                a.gap = "\t".to_owned();
            } else {
                a.gap = rest.to_owned();
            }
        } else if arg == "-dms" {
            a.dms = true;
        } else if arg.starts_with("-d") || arg.starts_with("-D") {
            let mode_char = arg.as_bytes()[1];
            let mut rest = &arg[2..];
            if rest.starts_with('h') {
                a.diff_mode = DiffMode::DiffHelio;
                rest = &rest[1..];
            } else if mode_char == b'd' {
                a.diff_mode = DiffMode::Diff;
            } else {
                a.diff_mode = DiffMode::DiffAbs;
            }
            a.diff_planet = rest.chars().next().unwrap_or('0');
        } else if arg == "-roundsec" {
            a.round_sec = true;
        } else if arg == "-roundmin" {
            a.round_min = true;
        } else if let Some(rest) = arg.strip_prefix("-t") {
            if !rest.is_empty() {
                a.time_input.push_str(rest);
            }
        } else if arg.starts_with("-h") || arg.starts_with("-?") {
            return Err("help".to_owned());
        } else {
            return Err(format!("illegal option {arg}"));
        }

        i += 1;
    }
    Ok(a)
}

#[derive(Debug, Clone)]
pub enum BodySpec {
    Planet(Body),
    Asteroid,
    FixedStar,
    PlanetMoon,
    Fictitious,
    EclipticNutation,
    Labels,
    DeltaT,
    TimeEquation,
    SiderealTime,
    Ayanamsha,
}

pub fn letter_to_body(letter: char) -> BodySpec {
    match letter {
        '0' => BodySpec::Planet(Body::Sun),
        '1' => BodySpec::Planet(Body::Moon),
        '2' => BodySpec::Planet(Body::Mercury),
        '3' => BodySpec::Planet(Body::Venus),
        '4' => BodySpec::Planet(Body::Mars),
        '5' => BodySpec::Planet(Body::Jupiter),
        '6' => BodySpec::Planet(Body::Saturn),
        '7' => BodySpec::Planet(Body::Uranus),
        '8' => BodySpec::Planet(Body::Neptune),
        '9' => BodySpec::Planet(Body::Pluto),
        'm' => BodySpec::Planet(Body::MeanNode),
        't' => BodySpec::Planet(Body::TrueNode),
        'A' => BodySpec::Planet(Body::MeanApogee),
        'B' => BodySpec::Planet(Body::OscuApogee),
        'C' => BodySpec::Planet(Body::Earth),
        'c' => BodySpec::Planet(Body::IntpApogee),
        'g' => BodySpec::Planet(Body::IntpPerigee),
        'D' => BodySpec::Planet(Body::Chiron),
        'E' => BodySpec::Planet(Body::Pholus),
        'F' => BodySpec::Planet(Body::Ceres),
        'G' => BodySpec::Planet(Body::Pallas),
        'H' => BodySpec::Planet(Body::Juno),
        'I' => BodySpec::Planet(Body::Vesta),
        c @ 'J'..='Z' => {
            let id = (c as i32) - ('J' as i32) + 40; // SE_CUPIDO=40
            BodySpec::Planet(Body::fictitious(id).unwrap())
        }
        'w' => BodySpec::Planet(Body::fictitious(57).unwrap()), // Waldemath
        's' => BodySpec::Asteroid,
        'v' => BodySpec::PlanetMoon,
        'z' => BodySpec::Fictitious,
        'f' => BodySpec::FixedStar,
        'n' | 'o' => BodySpec::EclipticNutation,
        'e' => BodySpec::Labels,
        'q' => BodySpec::DeltaT,
        'y' => BodySpec::TimeEquation,
        'x' => BodySpec::SiderealTime,
        'b' => BodySpec::Ayanamsha,
        'd' | 'p' | 'h' | 'a' => BodySpec::Labels, // preset selectors, should not appear as individual letters
        _ => BodySpec::Labels,
    }
}

impl SweTestArgs {
    pub fn build_iflag(&self) -> CalcFlags {
        let mut flags = CalcFlags::empty();

        match self.ephemeris {
            EphemerisChoice::Swiss => flags |= CalcFlags::SWIEPH,
            EphemerisChoice::Moshier => flags |= CalcFlags::MOSEPH,
            EphemerisChoice::Jpl => flags |= CalcFlags::JPLEPH,
        }

        if self.speed {
            flags |= CalcFlags::SPEED;
        }
        if self.speed3 {
            flags |= CalcFlags::SPEED3;
        }
        if self.truepos {
            flags |= CalcFlags::TRUEPOS;
        }
        if self.no_aberration {
            flags |= CalcFlags::NOABERR;
        }
        if self.no_deflection {
            flags |= CalcFlags::NOGDEFL;
        }
        if self.no_nutation {
            flags |= CalcFlags::NONUT;
        }
        if self.j2000 {
            flags |= CalcFlags::J2000;
        }
        if self.icrs {
            flags |= CalcFlags::ICRS;
        }
        if self.equatorial {
            flags |= CalcFlags::EQUATORIAL;
        }
        if self.xyz {
            flags |= CalcFlags::XYZ;
        }
        if self.topocentric {
            flags |= CalcFlags::TOPOCTR;
        }
        if self.heliocentric {
            flags |= CalcFlags::HELCTR;
        }
        if self.barycentric {
            flags |= CalcFlags::BARYCTR;
        }
        if self.sidereal {
            flags |= CalcFlags::SIDEREAL;
        }
        if self.center_of_body {
            flags |= CalcFlags::CENTER_BODY;
        }

        if let Some(iflag_f) = self.force_iflag {
            flags = CalcFlags::from_bits_truncate(iflag_f as u32);
        }

        // C's post-parse logic: if format contains S/s/Q and no SPEED3 and not no_speed, add SPEED
        if self.format.contains('S') || self.format.contains('s') || self.format.contains('Q') {
            if !flags.contains(CalcFlags::SPEED3) && !self.no_speed {
                flags |= CalcFlags::SPEED;
            }
        }

        flags
    }

    pub fn to_ephemeris_config(&self) -> EphemerisConfig {
        let mut config = EphemerisConfig {
            ephemeris_source: match self.ephemeris {
                EphemerisChoice::Swiss => EphemerisSource::Swiss,
                EphemerisChoice::Moshier => EphemerisSource::Moshier,
                EphemerisChoice::Jpl => EphemerisSource::Jpl,
            },
            ..EphemerisConfig::default()
        };

        if let Some(ref dir) = self.ephe_dir {
            config.ephe_path = Some(PathBuf::from(dir));
        } else if self.ephemeris != EphemerisChoice::Moshier {
            config.ephe_path = Some(PathBuf::from("."));
        }

        if self.ephemeris == EphemerisChoice::Jpl {
            config.jpl_filename = Some(self.jpl_file.clone());
        }

        if self.topocentric || self.have_geopos {
            if self.topocentric {
                config.topographic = Some(TopoPosition {
                    longitude: self.geo_longitude,
                    latitude: self.geo_latitude,
                    altitude: self.geo_elevation,
                });
            }
        }

        if self.sidereal || self.do_ayanamsa {
            config.set_sidereal_mode(
                if self.sid_mode == 255 {
                    self.sid_mode
                        | if self.sid_user_ut {
                            1024 // SE_SIDBIT_USER_UT
                        } else {
                            0
                        }
                } else {
                    self.sid_mode
                },
                self.sid_user_t0,
                self.sid_user_ayan,
            );
        }

        if let Some(tid_acc) = self.tidal_acc {
            config.tidal_acceleration = Some(tid_acc);
        }

        config
    }

    pub fn body_specs(&self) -> Vec<BodySpec> {
        self.planet_selection.chars().map(letter_to_body).collect()
    }
}
