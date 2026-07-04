use swisseph::{Body, CalcFlags, Ephemeris, EphemerisConfig, SiderealMode};

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

fn main() {
    // J2000.0
    let jd_ut = 2451545.0;

    // Tropical ephemeris
    let eph_trop = Ephemeris::new(EphemerisConfig::default()).unwrap();

    // Sidereal ephemeris (Lahiri ayanamsa)
    let mut sid_config = EphemerisConfig::default();
    sid_config.set_sidereal_mode(SiderealMode::Lahiri as i32, 0.0, 0.0);
    let eph_sid = Ephemeris::new(sid_config).unwrap();

    // Get the ayanamsa value at this date
    let ayanamsa = eph_sid.get_ayanamsa_ut(jd_ut, CalcFlags::empty()).unwrap();

    println!("Tropical vs Lahiri Sidereal positions at J2000.0 (JD 2451545.0)");
    println!("Lahiri ayanamsa: {ayanamsa:.6} degrees");
    println!();
    println!(
        "{:<10} {:>12} {:>12} {:>10}",
        "Body", "Tropical", "Sidereal", "Diff"
    );
    println!("{}", "-".repeat(48));

    for body in BODIES {
        let trop = eph_trop.calc_ut(jd_ut, body, CalcFlags::SPEED).unwrap();
        let sid = eph_sid
            .calc_ut(jd_ut, body, CalcFlags::SPEED | CalcFlags::SIDEREAL)
            .unwrap();

        let diff = trop.data[0] - sid.data[0];
        println!(
            "{:<10} {:>12.6} {:>12.6} {:>10.6}",
            format!("{body:?}"),
            trop.data[0],
            sid.data[0],
            diff,
        );
    }
}
