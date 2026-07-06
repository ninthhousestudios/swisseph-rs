mod args;
mod compute;
mod format;

use std::process;

use swisseph::Ephemeris;

use args::parse_args;

fn print_help() {
    eprintln!(
        "\
swetest — Rust port of the Swiss Ephemeris test program
Version: 0.1.0

Command-line flags:
  Input time:
    -b[date]      begin date, e.g. -b1.1.2000 (day.month.year)
    -bj[jd]       begin date as Julian day, e.g. -bj2451545
    -j[jd]        Julian day (same as -bj)
    -t[time]      time of day, e.g. -t12:30 or -t12.5
    -ut[time]     input is UT (Universal Time)
    -utc[time]    input is UTC
    -lmt[time]    input is Local Mean Time
    -lat          input is Local Apparent Time

  Body selection:
    -pd           default planets: 0123456789mtA
    -pp           main planets + main asteroids: 0123456789mtABCcgDEFGHI
    -ph           hypothetical/fictitious: JKLMNOPQRSTUVWXYZw
    -pa           all of the above
    -p[letters]   custom body selection (see planet letters below)
    -xs[n]        asteroid number n
    -xf[name]     fixed star name
    -xv[n]        planetary moon number n
    -xz[n]        fictitious body number n
    -x[name]      fixed star name (same as -xf)

  Planet letters:
    0=Sun  1=Moon  2=Mercury  3=Venus  4=Mars  5=Jupiter
    6=Saturn  7=Uranus  8=Neptune  9=Pluto  m=MeanNode  t=TrueNode
    A=MeanApogee  B=OscuApogee  C=Earth  c=InterpApogee  g=InterpPerigee
    D=Chiron  E=Pholus  F=Ceres  G=Pallas  H=Juno  I=Vesta
    J..Z=Fictitious(Cupido..WhiteMoon)  w=Waldemath
    s=asteroid(use -xs)  f=fixstar(use -xf)  v=planetmoon(use -xv)
    z=fictitious(use -xz)
    e=labels  q=deltaT  y=time_equation  x=sidereal_time  b=ayanamsha

  Stepping:
    -n[count]     number of steps (default 1, 0 → 20)
    -s[step]      step size, suffixes: y=years, o=months, m=minutes, s=seconds

  Observer/frame:
    -topo[lon,lat,elev]   topocentric (degrees, meters)
    -geopos[lon,lat,elev] geographic position (no TOPOCTR flag)
    -hel          heliocentric
    -bary         barycentric
    -pc[n]        planetocentric (center body n)

  Computation flags:
    -speed        compute speed (set automatically by format)
    -speed3       3-point speed
    -nospeed      suppress speed
    -noaberr      no aberration correction
    -nodefl       no gravitational deflection
    -nonut        no nutation
    -true         true/geometric position
    -j2000        J2000 equator/ecliptic
    -icrs         ICRS frame
    -cob          center-of-body flag
    -i[n]         force iflag to n

  Ephemeris source:
    -eswe         Swiss Ephemeris (default)
    -emos         Moshier ephemeris (built-in)
    -ejpl[file]   JPL ephemeris, optional filename
    -edir[path]   ephemeris file directory

  Sidereal:
    -sid[n]       sidereal mode n (0 = Fagan-Bradley)
    -sidt0[n]     sidereal, ecliptic of t0
    -sidsp[n]     sidereal, solar system plane
    -sidudef[t0,ayan0]  user-defined sidereal mode
    -sidbit[n]    sidereal bit flags
    -ay[n]        ayanamsa output for mode n

  Houses:
    -house[lon,lat,sys]   house cusps at position (P=Placidus etc.)
    -hsy[sys]     house system character

  Output format:
    -f[fmt]       format string (default PLBRS)
    -g[sep]       column separator (default space, empty → tab)
    -head         suppress header
    +head         always show header
    -hor          horizontal output
    -dms          degrees/minutes/seconds
    -roundsec     round to seconds
    -roundmin     round to minutes
    -ep           extra precision
    -short        short output

  Special events:
    -solecl       solar eclipse
    -lunecl       lunar eclipse
    -occult       occultation
    -rise         rise/set
    -metr         meridian transit
    -hev[type]    heliacal event
    -orbel        orbital elements

  Eclipse filters:
    -total  -partial  -annular  -anntot  -penumbral
    -central  -noncentral  -local  -how  -hocal

  Rise/set options:
    -norefrac     no atmospheric refraction
    -disccenter   disc center
    -discbottom   disc bottom edge
    -hindu        Hindu rise/set (disc center, no refraction)

  Heliacal:
    -at[p,t,rh,vis]   atmospheric conditions
    -obs[age,SN]      observer parameters
    -opt[age,SN,...]  optic parameters

  Differential:
    -d[planet]    differential (subtract planet from all)
    -D[planet]    differential (absolute value)
    -dh[planet]   differential heliocentric vs geocentric

  Miscellaneous:
    -bwd          backward (direction = -1)
    -tidacc[n]    override tidal acceleration
    -testaa95/96/97  Astronomical Almanac test
    -lim          show file limits
    -clink        chart link output
    -swefixstar2  use fixstar2 interface"
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let parsed = match parse_args(&args) {
        Ok(a) => a,
        Err(e) if e == "help" => {
            print_help();
            process::exit(0);
        }
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    };

    let config = parsed.to_ephemeris_config();

    match Ephemeris::new(config) {
        Ok(eph) => {
            compute::run(&parsed, &eph);
        }
        Err(e) => {
            eprintln!("swetest: failed to construct Ephemeris: {e}");
            process::exit(1);
        }
    }
}
