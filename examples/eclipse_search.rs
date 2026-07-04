use swisseph::{CalcFlags, EclipseFlags, Ephemeris, EphemerisConfig};

fn jd_to_date_str(jd: f64) -> String {
    let (y, m, d, h) = swisseph::date::revjul(jd, swisseph::types::CalendarType::Gregorian);
    let hours = h as u32;
    let minutes = ((h - hours as f64) * 60.0) as u32;
    format!("{y:04}-{m:02}-{d:02} {hours:02}:{minutes:02} UT")
}

fn eclipse_type_str(flags: EclipseFlags) -> &'static str {
    if flags.contains(EclipseFlags::TOTAL) {
        "Total"
    } else if flags.contains(EclipseFlags::HYBRID) {
        "Hybrid"
    } else if flags.contains(EclipseFlags::ANNULAR) {
        "Annular"
    } else if flags.contains(EclipseFlags::PARTIAL) {
        "Partial"
    } else if flags.contains(EclipseFlags::PENUMBRAL) {
        "Penumbral"
    } else {
        "Unknown"
    }
}

fn main() {
    let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();

    // Start searching from 2024-Jan-01
    let mut tjd = swisseph::date::julday(2024, 1, 1, 0.0, swisseph::types::CalendarType::Gregorian);

    println!("Next 3 solar eclipses from 2024-Jan-01:");
    println!();
    for _ in 0..3 {
        let ecl = eph
            .sol_eclipse_when_glob(tjd, CalcFlags::empty(), EclipseFlags::empty(), false)
            .unwrap();
        println!(
            "  {} - {} solar eclipse",
            jd_to_date_str(ecl.time_maximum),
            eclipse_type_str(ecl.flags),
        );
        tjd = ecl.time_maximum + 1.0;
    }

    println!();
    println!("Next 3 lunar eclipses from 2024-Jan-01:");
    println!();
    tjd = swisseph::date::julday(2024, 1, 1, 0.0, swisseph::types::CalendarType::Gregorian);
    for _ in 0..3 {
        let ecl = eph
            .lun_eclipse_when(tjd, CalcFlags::empty(), EclipseFlags::empty(), false)
            .unwrap();
        println!(
            "  {} - {} lunar eclipse",
            jd_to_date_str(ecl.time_maximum),
            eclipse_type_str(ecl.flags),
        );
        tjd = ecl.time_maximum + 1.0;
    }
}
