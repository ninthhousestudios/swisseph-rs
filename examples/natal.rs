use swisseph::{
    Body, CalcFlags, DegreeParts, Ephemeris, EphemerisConfig, HouseSystem, SplitDegFlags,
};

const BODIES: [Body; 10] = [
    Body::Sun,
    Body::Moon,
    Body::Mercury,
    Body::Venus,
    Body::Mars,
    Body::Jupiter,
    Body::Saturn,
    Body::Uranus,
    Body::Neptune,
    Body::Pluto,
];

const SIGNS: [&str; 12] = [
    "Ari", "Tau", "Gem", "Can", "Leo", "Vir", "Lib", "Sco", "Sag", "Cap", "Aqu", "Pis",
];

fn format_zodiacal(lon: f64) -> String {
    let parts: DegreeParts =
        swisseph::math::split_degrees(lon, SplitDegFlags::ZODIACAL | SplitDegFlags::ROUND_SEC);
    format!(
        "{:2}{}{}' {:02}\"",
        parts.degrees, SIGNS[parts.sign as usize], parts.minutes, parts.seconds
    )
}

fn main() {
    // 1985-Jul-15 14:30 UT, Zurich (8.55 E, 47.37 N)
    let jd_ut = swisseph::date::julday(1985, 7, 15, 14.5, swisseph::types::CalendarType::Gregorian);

    let eph = Ephemeris::new(EphemerisConfig::default()).unwrap();

    println!("Natal chart: 1985-Jul-15 14:30 UT, Zurich (8.55E, 47.37N)");
    println!("JD(UT) = {jd_ut:.6}");
    println!();
    println!("--- Planetary positions (tropical) ---");
    println!("{:<10} {:>16}", "Body", "Longitude");

    for body in BODIES {
        let result = eph.calc_ut(jd_ut, body, CalcFlags::SPEED).unwrap();
        println!(
            "{:<10} {:>16}",
            format!("{body:?}"),
            format_zodiacal(result.data[0])
        );
    }

    println!();
    println!("--- Placidus house cusps ---");
    let houses = eph
        .houses_ex2(
            jd_ut,
            CalcFlags::empty(),
            47.37,
            8.55,
            HouseSystem::Placidus,
        )
        .unwrap();

    for i in 1..=12 {
        println!("  Cusp {:2}: {}", i, format_zodiacal(houses.cusps[i]));
    }

    println!();
    println!("  Asc: {}", format_zodiacal(houses.ascmc.ascendant));
    println!("  MC:  {}", format_zodiacal(houses.ascmc.mc));
}
