use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, RiseSetFlags};

fn jd_to_time_str(jd: f64) -> String {
    let (y, m, d, h) = swisseph::date::revjul(jd, swisseph::types::CalendarType::Gregorian);
    let hours = h as u32;
    let minutes = ((h - hours as f64) * 60.0) as u32;
    format!("{y:04}-{m:02}-{d:02} {hours:02}:{minutes:02} UT")
}

fn main() {
    // 2024-Mar-20, observer at Zurich (8.55 E, 47.37 N, 500m)
    let jd_ut = swisseph::date::julday(2024, 3, 20, 0.0, swisseph::types::CalendarType::Gregorian);
    let geopos = [8.55, 47.37, 500.0]; // [lon_east, lat_north, alt_m]
    let atpress = 1013.25; // standard pressure (mbar)
    let attemp = 15.0; // standard temperature (C)

    let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();

    println!("Rise/set times for 2024-Mar-20, Zurich (8.55E, 47.37N, 500m)");
    println!();

    for body in [Body::Sun, Body::Moon] {
        let rise = eph
            .rise_trans(
                jd_ut,
                body,
                None,
                CalcFlags::empty(),
                RiseSetFlags::RISE,
                geopos,
                atpress,
                attemp,
            )
            .unwrap();

        let set = eph
            .rise_trans(
                jd_ut,
                body,
                None,
                CalcFlags::empty(),
                RiseSetFlags::SET,
                geopos,
                atpress,
                attemp,
            )
            .unwrap();

        println!(
            "{:6} rise: {}",
            format!("{body:?}"),
            jd_to_time_str(rise.time)
        );
        println!("{:6}  set: {}", "", jd_to_time_str(set.time));
        println!();
    }
}
