// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

//! Fixed-star catalog — loading, searching, and the 8 built-in ayanamsa
//! reference stars.

use std::collections::HashMap;
use std::path::Path;

use crate::constants::{DEGTORAD, KM_S_TO_AU_CTY};
use crate::error::Error;

/// A single fixed-star catalog record (one line of `sefstars.txt`), with angular quantities
/// already converted to radians and proper motion/parallax/radial velocity to per-century units.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Star {
    /// Traditional star name (e.g. "Sirius"), or empty for Bayer-only records.
    pub name: String,
    /// Bayer/Flamsteed designation (e.g. "alfCMa").
    pub bayer: String,
    /// Normalized search key used for catalog lookup (lowercased name, or `,`-prefixed Bayer).
    pub skey: String,
    /// Reference epoch: `0.0` for ICRS, else a Julian-year epoch (e.g. `1950.0`, `2000.0`).
    pub epoch: f64,
    /// Right ascension at epoch, radians.
    pub ra: f64,
    /// Declination at epoch, radians.
    pub de: f64,
    /// Proper motion in right ascension, radians/century (already divided by `cos(de)`).
    pub ramot: f64,
    /// Proper motion in declination, radians/century.
    pub demot: f64,
    /// Radial velocity, AU/century.
    pub radvel: f64,
    /// Parallax, radians.
    pub parall: f64,
    /// Visual magnitude.
    pub mag: f64,
}

/// In-memory fixed-star catalog, indexed for lookup by sequential number, Bayer designation,
/// and traditional name.
pub struct StarCatalog {
    /// Records with a unique Bayer designation, sorted by `skey` (used for sequential-number
    /// lookup and Bayer-designation search).
    pub bayer_records: Vec<Star>,
    /// Records with a non-empty traditional name (used for name search).
    pub named_records: Vec<Star>,
    by_bayer: HashMap<String, usize>,
    by_name: HashMap<String, usize>,
}

impl StarCatalog {
    /// An empty catalog (no records), returned when `sefstars.txt` is unavailable.
    pub fn empty() -> Self {
        StarCatalog {
            bayer_records: Vec::new(),
            named_records: Vec::new(),
            by_bayer: HashMap::new(),
            by_name: HashMap::new(),
        }
    }

    /// Number of unique Bayer-designated records in the catalog.
    pub fn n_real(&self) -> usize {
        self.bayer_records.len()
    }

    /// Number of named records in the catalog.
    pub fn n_named(&self) -> usize {
        self.named_records.len()
    }

    /// Look up a star by sequential number, Bayer designation (`,`-prefixed, optionally
    /// wildcarded with a trailing `%`), or traditional name.
    pub fn search(&self, input: &str) -> Result<Star, Error> {
        let normalized = format_search_name(input)?;

        // Mode 1: sequential number (first char is ASCII digit)
        if normalized.starts_with(|c: char| c.is_ascii_digit()) {
            let n: usize = normalized
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>()
                .parse()
                .unwrap_or(0);
            if n < 1 || n > self.n_real() {
                return Err(Error::FileFormat(format!(
                    "sequential fixed star number {n} is not available"
                )));
            }
            return Ok(self.bayer_records[n - 1].clone());
        }

        // Wildcard mode: ends with '%'
        if normalized.ends_with('%') {
            let prefix = &normalized[..normalized.len() - 1];
            for star in &self.named_records {
                if star.skey.starts_with(prefix) {
                    return Ok(star.clone());
                }
            }
            return Err(Error::FileFormat(format!("star not found: {input}")));
        }

        // Mode 2: Bayer designation (contains ',')
        if let Some(comma_pos) = normalized.find(',') {
            let bayer_key = &normalized[comma_pos..];
            if let Some(&idx) = self.by_bayer.get(bayer_key) {
                return Ok(self.bayer_records[idx].clone());
            }
            return Err(Error::FileFormat(format!("star not found: {input}")));
        }

        // Mode 3: traditional name
        if let Some(&idx) = self.by_name.get(&normalized) {
            return Ok(self.named_records[idx].clone());
        }

        Err(Error::FileFormat(format!("star not found: {input}")))
    }
}

// Parse the epoch field: "ICRS" → 0.0, "2000" → 2000.0, "1950" → 1950.0.
// sefstars.txt only ever holds these three clean tokens, so a strict full-string
// parse suffices; this is NOT a faithful atof (which would scan a leading float
// from e.g. "1950BESSEL"). Non-numeric → 0.0, matching atof on "ICRS".
fn parse_epoch(s: &str) -> f64 {
    s.trim().parse::<f64>().unwrap_or(0.0)
}

fn cut_string(line: &str) -> Result<Star, Error> {
    let fields: Vec<&str> = line.split(',').map(|f| f.trim()).collect();
    if fields.len() < 14 {
        return Err(Error::FileFormat(format!(
            "data of star '{},{}' incomplete",
            fields.first().unwrap_or(&""),
            fields.get(1).unwrap_or(&""),
        )));
    }

    let name = fields[0].to_string();
    let bayer = fields[1].to_string();
    let epoch = parse_epoch(fields[2]);

    let ra_h: f64 = fields[3].parse().unwrap_or(0.0);
    let ra_m: f64 = fields[4].parse().unwrap_or(0.0);
    let ra_s: f64 = fields[5].parse().unwrap_or(0.0);

    let sde_d = fields[6]; // raw string for sign detection (handles "-00")
    let de_d: f64 = sde_d.parse().unwrap_or(0.0);
    let de_m: f64 = fields[7].parse().unwrap_or(0.0);
    let de_s: f64 = fields[8].parse().unwrap_or(0.0);

    let mut ra_pm: f64 = fields[9].parse().unwrap_or(0.0);
    let mut de_pm: f64 = fields[10].parse().unwrap_or(0.0);
    let mut radvel: f64 = fields[11].parse().unwrap_or(0.0);
    let mut parall: f64 = fields[12].parse().unwrap_or(0.0);
    let mag: f64 = fields[13].parse().unwrap_or(999.99);

    // RA → degrees (hours to degrees)
    let mut ra = (ra_s / 3600.0 + ra_m / 60.0 + ra_h) * 15.0;

    // Dec → degrees; sign from raw field-6 string (strchr(sde_d, '-') in C)
    // Handles "-00" correctly where atof would return 0.0 with wrong sign
    let mut de = if sde_d.contains('-') {
        -de_s / 3600.0 - de_m / 60.0 + de_d
    } else {
        de_s / 3600.0 + de_m / 60.0 + de_d
    };

    // Proper motion: 0.001 arcsec/yr (new sefstars.txt format) → deg/century
    ra_pm = ra_pm / 10.0 / 3600.0;
    de_pm = de_pm / 10.0 / 3600.0;

    // Parallax sign fix (handles historical bug in old catalog entries)
    if parall < 0.0 {
        parall = -parall;
    }
    parall /= 1000.0; // mas → arcsec

    // Radial velocity: km/s → AU/century
    radvel *= KM_S_TO_AU_CTY;

    // Parallax: arcsec → degrees
    if parall > 1.0 {
        parall = 1.0 / parall / 3600.0;
    } else {
        parall /= 3600.0;
    }

    // All angular quantities → radians
    ra *= DEGTORAD;
    de *= DEGTORAD;
    ra_pm *= DEGTORAD;
    de_pm *= DEGTORAD;
    // Catalog stores μα×cos(δ); divide once to recover the pure spherical RA rate
    ra_pm /= de.cos();
    parall *= DEGTORAD;

    Ok(Star {
        name,
        bayer,
        skey: String::new(),
        epoch,
        ra,
        de,
        ramot: ra_pm,
        demot: de_pm,
        radvel,
        parall,
        mag,
    })
}

// Normalize caller input to a canonical search key (port of fixstar_format_search_name).
// Truncates to 40 chars, removes spaces, lowercases the part before the first comma;
// the Bayer part after the comma is left case-sensitive.
fn format_search_name(input: &str) -> Result<String, Error> {
    let end = input.char_indices().nth(40).map_or(input.len(), |(i, _)| i);
    let truncated = &input[..end];
    let no_spaces: String = truncated.chars().filter(|c| *c != ' ').collect();
    let result = if let Some(pos) = no_spaces.find(',') {
        let before = no_spaces[..pos].to_lowercase();
        let after = &no_spaces[pos..];
        format!("{before}{after}")
    } else {
        no_spaces.to_lowercase()
    };
    if result.is_empty() {
        return Err(Error::FileFormat("star name empty".to_string()));
    }
    Ok(result)
}

fn parse_catalog(contents: &str) -> Result<StarCatalog, Error> {
    let mut bayer_records: Vec<Star> = Vec::new();
    let mut named_records: Vec<Star> = Vec::new();
    let mut by_name: HashMap<String, usize> = HashMap::new();
    let mut last_bayer = String::new();

    for line in contents.lines() {
        // Skip blank and whitespace-only lines silently (fixes C's 'data corrupted' bug on blank
        // lines in the legacy swe_fixstar file-scan path; the load_all_fixed_stars path also
        // skips blank lines in C but errors only in the older path).
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }

        let star = cut_string(line)?;

        // Named record (field 0 non-empty)
        if !star.name.is_empty() {
            let skey: String = star
                .name
                .to_lowercase()
                .chars()
                .filter(|c| *c != ' ')
                .collect();
            let mut named = star.clone();
            named.skey = skey.clone();
            let idx = named_records.len();
            named_records.push(named);
            by_name.insert(skey, idx);
        }

        // Bayer record (consecutive dedup)
        if star.bayer == last_bayer {
            continue;
        }
        let bayer_no_spaces: String = star.bayer.chars().filter(|c| *c != ' ').collect();
        let skey = format!(",{bayer_no_spaces}");
        let mut bayer = star.clone();
        bayer.skey = skey;
        bayer_records.push(bayer);
        last_bayer = star.bayer;
    }

    // Sort bayer_records by skey (lexicographic == C strcmp for ASCII);
    // sequential-number lookup uses position in this sorted array.
    bayer_records.sort_by(|a, b| a.skey.cmp(&b.skey));

    let by_bayer: HashMap<String, usize> = bayer_records
        .iter()
        .enumerate()
        .map(|(i, s)| (s.skey.clone(), i))
        .collect();

    Ok(StarCatalog {
        bayer_records,
        named_records,
        by_bayer,
        by_name,
    })
}

/// Hardcoded fallback records for the 8 ayanamsa reference stars (port of get_builtin_star).
/// Called before catalog search so these resolve even when sefstars.txt is unavailable.
/// Note: SgrA* uses epoch "2000" (FK5) not "ICRS" — replicating C exactly (sweph.c:6783).
pub fn builtin_star(input: &str) -> Option<Star> {
    let record = if input.starts_with("spica") || input.starts_with("Spica") {
        "Spica,alVir,ICRS,13,25,11.57937,-11,09,40.7501,-42.35,-30.67,1,13.06,0.97,-10,3672"
    } else if input.contains(",zePsc") || input.starts_with("revati") || input.starts_with("Revati")
    {
        "Revati,zePsc,ICRS,01,13,43.88735,+07,34,31.2745,145,-55.69,15,18.76,5.187,06,174"
    } else if input.contains(",deCnc") || input.starts_with("pushya") || input.starts_with("Pushya")
    {
        // Also matches ",deCnc" for the Sheoran ayanamsa (returns same record)
        "Pushya,deCnc,ICRS,08,44,41.09921,+18,09,15.5034,-17.67,-229.26,17.14,24.98,3.94,18,2027"
    } else if input.contains(",laSco") || input.starts_with("mula") || input.starts_with("Mula") {
        "Mula,laSco,ICRS,17,33,36.52012,-37,06,13.7648,-8.53,-30.8,-3,5.71,1.62,-37,11673"
    } else if input.contains(",SgrA*") {
        "Gal. Center,SgrA*,2000,17,45,40.03599,-29,00,28.1699,-2.755718425,-5.547,0.0,0.125,999.99,0,0"
    } else if input.contains(",GP1958") {
        "Gal. Pole IAU1958,GP1958,1950,12,49,0.0,27,24,0.0,0.0,0.0,0.0,0.0,0.0,0,0"
    } else if input.contains(",GPol") {
        // Matches both GALEQU_TRUE and GALEQU_MULA (same record)
        "Gal. Pole,GPol,ICRS,12,51,36.7151981,27,06,11.193172,0.0,0.0,0.0,0.0,0.0,0,0"
    } else {
        return None;
    };
    cut_string(record).ok()
}

/// Load star catalog from `ephe_path/sefstars.txt`. Returns an empty catalog on any failure;
/// builtin_star() still resolves the 8 ayanamsa reference stars without the file.
pub fn load_catalog(ephe_path: Option<&Path>) -> StarCatalog {
    let Some(dir) = ephe_path else {
        return StarCatalog::empty();
    };
    let path = dir.join("sefstars.txt");
    let contents = match std::fs::read_to_string(&path) {
        Ok(s) => s,
        Err(_) => return StarCatalog::empty(),
    };
    parse_catalog(&contents).unwrap_or_else(|_| StarCatalog::empty())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-12;

    fn check(actual: f64, expected: f64, label: &str) {
        assert!(
            (actual - expected).abs() < EPS,
            "{label}: actual={actual:.15e} expected={expected:.15e} diff={:.3e}",
            (actual - expected).abs()
        );
    }

    #[test]
    fn test_cut_string_sirius() {
        let line =
            "alfCMa,alfCMa,ICRS,06,45,08.91728,-16,42,58.0171,-546.01,-1223.07,-5.50,379.21,-1.46";
        let s = cut_string(line).unwrap();

        assert_eq!(s.name, "alfCMa");
        assert_eq!(s.bayer, "alfCMa");
        assert_eq!(s.epoch, 0.0); // ICRS

        let ra_deg = (8.91728_f64 / 3600.0 + 45.0 / 60.0 + 6.0) * 15.0;
        // field[6] = "-16" contains '-' → negative dec
        let de_deg = -58.0171_f64 / 3600.0 - 42.0 / 60.0 + (-16.0_f64);
        let de_rad = de_deg * PI / 180.0;

        check(s.ra, ra_deg * PI / 180.0, "ra");
        check(s.de, de_rad, "de");
        check(
            s.ramot,
            (-546.01_f64 / 10.0 / 3600.0 * PI / 180.0) / de_rad.cos(),
            "ramot",
        );
        check(s.demot, -1223.07_f64 / 10.0 / 3600.0 * PI / 180.0, "demot");
        check(s.radvel, -5.50_f64 * KM_S_TO_AU_CTY, "radvel");
        // parall = 379.21 mas; < 1000 arcsec after /1000 → /3600 path
        check(
            s.parall,
            379.21_f64 / 1000.0 / 3600.0 * PI / 180.0,
            "parall",
        );
        check(s.mag, -1.46, "mag");
    }

    #[test]
    fn test_cut_string_aldebaran() {
        let line =
            "alfTau,alfTau,ICRS,04,35,55.23907,+16,30,33.4885,63.45,-188.94,54.398,48.94,0.86";
        let s = cut_string(line).unwrap();

        assert_eq!(s.epoch, 0.0);

        let ra_deg = (55.23907_f64 / 3600.0 + 35.0 / 60.0 + 4.0) * 15.0;
        // field[6] = "+16" — no '-', positive dec
        let de_deg = 33.4885_f64 / 3600.0 + 30.0 / 60.0 + 16.0_f64;
        let de_rad = de_deg * PI / 180.0;

        check(s.ra, ra_deg * PI / 180.0, "ra");
        check(s.de, de_rad, "de");
        check(
            s.ramot,
            (63.45_f64 / 10.0 / 3600.0 * PI / 180.0) / de_rad.cos(),
            "ramot",
        );
        check(s.demot, -188.94_f64 / 10.0 / 3600.0 * PI / 180.0, "demot");
        check(s.radvel, 54.398_f64 * KM_S_TO_AU_CTY, "radvel");
        check(s.parall, 48.94_f64 / 1000.0 / 3600.0 * PI / 180.0, "parall");
        check(s.mag, 0.86, "mag");
    }

    #[test]
    fn test_cut_string_sgra_star() {
        // SgrA* from sefstars.txt (ICRS epoch, not the builtin which uses epoch 2000)
        let line = "SgrA*,SgrA*,ICRS,17,45,40.03599,-29,00,28.1699,-2.755718425,-5.547,0.0,0.125,999.99,0,0";
        let s = cut_string(line).unwrap();

        assert_eq!(s.name, "SgrA*");
        assert_eq!(s.bayer, "SgrA*");
        assert_eq!(s.epoch, 0.0); // ICRS

        let ra_deg = (40.03599_f64 / 3600.0 + 45.0 / 60.0 + 17.0) * 15.0;
        let de_deg = -28.1699_f64 / 3600.0 - 0.0 / 60.0 + (-29.0_f64);
        let de_rad = de_deg * PI / 180.0;

        check(s.ra, ra_deg * PI / 180.0, "ra");
        check(s.de, de_rad, "de");
        check(s.radvel, 0.0, "radvel");
        check(s.parall, 0.125_f64 / 1000.0 / 3600.0 * PI / 180.0, "parall");
        check(s.mag, 999.99, "mag");
    }

    #[test]
    fn test_parse_catalog_skips_blank_and_comment() {
        let input = concat!(
            "# comment line\n",
            "alfCMa,alfCMa,ICRS,06,45,08.91728,-16,42,58.0171,-546.01,-1223.07,-5.50,379.21,-1.46\n",
            "\n",
            "   \n",
            "Sirius,alfCMa,ICRS,06,45,08.91728,-16,42,58.0171,-546.01,-1223.07,-5.50,379.21,-1.46\n",
        );
        let cat = parse_catalog(input).unwrap();
        assert_eq!(cat.n_real(), 1); // one unique bayer (alfCMa)
        assert_eq!(cat.n_named(), 2); // "alfCMa" + "Sirius"
    }

    #[test]
    fn test_load_real_catalog() {
        let cat = load_catalog(Some(Path::new("../swisseph/ephe")));
        assert!(
            cat.n_real() > 1000,
            "expected >1000 bayer records, got {}",
            cat.n_real()
        );
        assert!(
            cat.n_named() > 1000,
            "expected >1000 named records, got {}",
            cat.n_named()
        );
    }

    #[test]
    fn test_search_by_name() {
        let cat = load_catalog(Some(Path::new("../swisseph/ephe")));
        let s = cat.search("Sirius").unwrap();
        assert_eq!(s.bayer, "alfCMa");
    }

    #[test]
    fn test_search_by_bayer() {
        let cat = load_catalog(Some(Path::new("../swisseph/ephe")));
        let s = cat.search(",alfCMa").unwrap();
        assert_eq!(s.bayer, "alfCMa");
    }

    #[test]
    fn test_search_sequential() {
        let cat = load_catalog(Some(Path::new("../swisseph/ephe")));
        let s = cat.search("1").unwrap();
        // Sequential lookup returns a bayer record (comma-prefixed skey, sorted)
        assert!(
            s.skey.starts_with(','),
            "skey should be comma-prefixed: {}",
            s.skey
        );
    }

    #[test]
    fn test_builtin_spica() {
        let s = builtin_star("Spica").unwrap();
        assert_eq!(s.epoch, 0.0); // ICRS
    }

    #[test]
    fn test_builtin_sgra_star() {
        let s = builtin_star(",SgrA*").unwrap();
        assert_eq!(s.epoch, 2000.0); // FK5 — intentional C behaviour (sweph.c:6783)
    }

    #[test]
    fn test_builtin_gp1958() {
        let s = builtin_star(",GP1958").unwrap();
        assert_eq!(s.epoch, 1950.0); // FK4 B1950
    }

    #[test]
    fn test_builtin_gpol() {
        let s = builtin_star(",GPol").unwrap();
        assert_eq!(s.epoch, 0.0); // ICRS
    }

    #[test]
    fn test_builtin_none() {
        assert!(builtin_star("unknown_star").is_none());
    }
}
