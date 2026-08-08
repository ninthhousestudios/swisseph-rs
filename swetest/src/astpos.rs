//! `-astpos` option: list named asteroids near an ecliptic longitude.
//!
//! Port of swetest.c's `swe_get_named_ast_list` + `print_asteroids` (the C
//! comment notes these "should move to swephlib.c in next release" — for now
//! they live in the test program, so they live here too, not in the library).

use std::fs;
use std::path::PathBuf;

use swisseph::Ephemeris;
use swisseph::flags::CalcFlags;

use crate::args::SweTestArgs;
use crate::compute::{make_asteroid_body, resolve_start_jd};
use crate::format::{DmsFlags, dms};

const DEGREE_SIGN: &str = "\u{b0}";

/// Read `seasnam.txt` and collect the catalog numbers of *named* asteroids.
///
/// A line names an asteroid when the character at column 8 (0-based) is not a
/// digit: provisional designations (e.g. `1998 SF36`) begin with a year digit
/// there and are skipped. A line too short to have a column 8 is treated as a
/// non-digit and thus included — matching swetest.c's benign out-of-bounds read
/// of `si[8]` on the handful of name-less entries at the top of the file.
pub(crate) fn get_named_ast_list(ephe_dir: &str) -> Result<Vec<i32>, String> {
    let path = PathBuf::from(ephe_dir).join(swisseph::constants::ASTNAMFILE);
    let contents =
        fs::read_to_string(&path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let mut list = Vec::new();
    for line in contents.lines() {
        let col8_is_digit = line.as_bytes().get(8).is_some_and(u8::is_ascii_digit);
        if col8_is_digit {
            continue;
        }
        // C uses atoi(si): skip leading whitespace, then read the leading digits.
        let num: String = line
            .trim_start()
            .chars()
            .take_while(char::is_ascii_digit)
            .collect();
        if let Ok(n) = num.parse::<i32>() {
            list.push(n);
        }
    }
    Ok(list)
}

/// Entry point for `-astpos`. Computes at the begin-date epoch and prints every
/// named asteroid whose ecliptic longitude is within `astpos_orb` of the target.
pub(crate) fn run(args: &SweTestArgs, eph: &Ephemeris) {
    let Some(dref) = args.astpos else {
        return;
    };
    let orb = args.astpos_orb;

    let config = args.to_ephemeris_config();
    let info = resolve_start_jd(args, &config);

    // Same ephe-path resolution as to_ephemeris_config for the Swiss backend.
    let ephe_dir = args.ephe_dir.clone().unwrap_or_else(|| ".".to_owned());
    let list = match get_named_ast_list(&ephe_dir) {
        Ok(l) => l,
        Err(e) => {
            eprintln!(" error in swe_get_named_ast_list(): {e}");
            return;
        }
    };
    let nast = list.len();

    let si = dms(dref, DmsFlags::ZODIAC | DmsFlags::ROUND_SEC, false);
    println!();
    println!("Asteroids near {dref:.6}{DEGREE_SIGN} ({si}) within orb {orb:.3}{DEGREE_SIGN}");
    println!("\t(out of {nast} named asteroids)");
    println!();

    for num in list {
        let Some(body) = make_asteroid_body(num) else {
            continue;
        };
        // C calls swe_calc with iflag=0, whose default ephemeris is Swiss
        // (SEFLG_DEFAULTEPH == SEFLG_SWIEPH) — geocentric tropical ecliptic,
        // ignoring the user's -sid / -topo / -emos flags. Our stateless calc has
        // no implicit default, so SWIEPH is set explicitly. Missing asteroid
        // files return Err and are skipped, matching C's `rc >= 0` guard.
        let Ok(res) = eph.calc(info.tjd_tt, body, CalcFlags::SWIEPH) else {
            continue;
        };
        let mut d = swisseph::math::diff_degrees(dref, res.data[0]);
        if d.abs() > orb {
            continue;
        }
        let m = if d < 0.0 {
            d = -d;
            '-'
        } else {
            ' '
        };
        let pname = eph.get_planet_name(body);
        // C: printf("%.3f%c\t%d\t%-20s\n", d, m, arr[i], pname) — the orb sign
        // char sits flush against the number for easy `sort`-ing.
        println!("{d:.3}{m}\t{num}\t{pname:<20}");
    }
}

#[cfg(test)]
mod tests {
    use super::get_named_ast_list;

    #[test]
    fn named_list_uses_column_8_and_atoi() {
        // Fixed-width seasnam.txt: number in cols 0-5, name from col 8.
        let dir = std::env::temp_dir().join(format!("astpos_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let contents = concat!(
            "000022 \n",          // no name: col 8 absent -> included as 22 (C's benign OOB read)
            "000001  Ceres\n",    // col 8 = 'C' (letter) -> named, atoi = 1
            "025143  Itokawa\n",  // col 8 = 'I' (letter) -> named, atoi = 25143
            "099942  1998 XY1\n", // col 8 = '1' (digit) -> provisional designation, skipped
            "# comment line\n",   // col 8 = 'e' (letter) but atoi finds no digits -> dropped
        );
        let path = dir.join(swisseph::constants::ASTNAMFILE);
        std::fs::write(&path, contents).unwrap();

        let list = get_named_ast_list(dir.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&dir).ok();

        // 22 (no-name, included), 1 (Ceres), 25143 (has a name at col 8);
        // 99942 skipped (digit at col 8); comment line yields atoi 0 -> dropped.
        assert_eq!(list, vec![22, 1, 25143]);
    }
}
