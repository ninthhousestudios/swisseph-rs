// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

// This file is part of swisseph-rs.
//
// Copyright (c) 2025 Josh Harper <josh@ninthhouse.studio>
//
// swisseph-rs is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

use std::collections::HashSet;
use std::fmt::Write as _;
use std::fs;
use std::io::{self, Write};
use std::process;

use clap::Parser;

const SIMBAD_BASE: &str = "https://simbad.cds.unistra.fr/simbad/sim-id";
const EXPECTED_FIELD_COUNT: usize = 14;

#[derive(Parser)]
#[command(
    name = "make-swe-stars",
    about = "Generate sefstars.txt entries for fixed stars from SIMBAD data."
)]
struct Cli {
    /// Append entries to this file (default: print to stdout).
    #[arg(short, long)]
    output_file: Option<String>,

    /// Read star names from a file (one per line, # comments allowed).
    #[arg(short, long)]
    input_file: Option<String>,

    /// Show what would be generated without writing anything.
    #[arg(short = 'n', long)]
    dry_run: bool,

    /// Verify existing sefstars.txt entries against current SIMBAD data.
    #[arg(long, num_args = 1..)]
    verify: Vec<String>,

    /// Star identifiers (Bayer, Flamsteed, HIP, name).
    stars: Vec<String>,
}

fn main() {
    let cli = Cli::parse();

    if !cli.verify.is_empty() {
        let path = cli.output_file.as_deref().unwrap_or("sefstars.txt");
        verify_against_simbad(&cli.verify, path);
        return;
    }

    let mut star_names: Vec<String> = cli.stars;
    if let Some(ref path) = cli.input_file {
        match read_star_list(path) {
            Ok(names) => star_names.extend(names),
            Err(e) => {
                eprintln!("Error reading input file '{path}': {e}");
                process::exit(1);
            }
        }
    }

    if star_names.is_empty() {
        Cli::parse_from(["make-swe-stars", "--help"]);
        return;
    }

    let mut existing_names = HashSet::new();
    if let Some(ref path) = cli.output_file {
        existing_names = load_existing_names(path);
    }

    let mut all_output = String::new();

    for star in &star_names {
        let parsed = match query_simbad(star) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("Error querying SIMBAD for '{star}': {e}");
                continue;
            }
        };

        let entry_lines = build_entry_lines(&parsed);

        let mut skipped = Vec::new();
        let mut added = Vec::new();
        for line in &entry_lines {
            if line.starts_with('#') {
                added.push(line.as_str());
                continue;
            }
            let first_field = line.split(',').next().unwrap_or("");
            if existing_names.contains(first_field) {
                skipped.push(first_field);
            } else {
                added.push(line.as_str());
                existing_names.insert(first_field.to_owned());
            }
        }

        if !skipped.is_empty() {
            eprintln!(
                "Skipping duplicate names for '{star}': {}",
                skipped.join(", ")
            );
        }

        if !validate_entry(&added, star) {
            continue;
        }

        if cli.dry_run {
            println!("--- {star} ---");
            for line in &added {
                print!("{line}");
            }
        } else {
            for line in &added {
                all_output.push_str(line);
            }
        }
    }

    if cli.dry_run {
        return;
    }

    if let Some(ref path) = cli.output_file {
        if let Err(e) = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .and_then(|mut f| f.write_all(all_output.as_bytes()))
        {
            eprintln!("Error writing to '{path}': {e}");
            process::exit(1);
        }
    } else {
        print!("{all_output}");
    }
}

#[derive(Debug)]
struct SimbadResult {
    trad_name: String,
    nomen_name: String,
    hip_id: String,
    ra_hour: String,
    ra_minute: String,
    ra_sec: String,
    dec_degree: String,
    dec_minute: String,
    dec_sec: String,
    pmra: String,
    pmde: String,
    rad_vel: String,
    parallax: String,
    mag_v: String,
}

fn simbad_url(name: &str) -> String {
    let encoded = name.replace(' ', "+");
    format!(
        "{SIMBAD_BASE}?Ident={encoded}&NbIdent=1&Radius=2&Radius.unit=arcmin\
         &submit=submit%20id&output.format=ASCII"
    )
}

fn query_simbad(name: &str) -> Result<SimbadResult, String> {
    let url = simbad_url(name);
    let body = ureq::get(&url)
        .call()
        .map_err(|e| format!("Network error: {e}"))?
        .body_mut()
        .read_to_string()
        .map_err(|e| format!("Read error: {e}"))?;

    parse_simbad_response(&body)
        .ok_or_else(|| format!("Could not parse SIMBAD response for '{name}'"))
}

fn parse_simbad_response(text: &str) -> Option<SimbadResult> {
    let lines: Vec<&str> = text.lines().collect();

    let mut mag_v = None;
    let mut hip_id = String::from("no hip id");
    let mut trad_name = String::new();
    let mut nomen_name = String::from("noMen");
    let mut ra_hour = None;
    let mut ra_minute = None;
    let mut ra_sec = None;
    let mut dec_degree = None;
    let mut dec_minute = None;
    let mut dec_sec = None;
    let mut pmra = None;
    let mut pmde = None;
    let mut parallax = None;
    let mut rad_vel = None;

    for (n, line) in lines.iter().enumerate() {
        if n > 28 {
            for tail_line in lines.iter().take(lines.len().min(50)).skip(28) {
                if tail_line.contains("HIP") {
                    let parts: Vec<&str> = tail_line.split_whitespace().collect();
                    for (i, element) in parts.iter().enumerate() {
                        if *element == "HIP"
                            && let Some(num) = parts.get(i + 1)
                        {
                            hip_id = format!("HIP {num}");
                        }
                    }
                }
                if tail_line.contains("NAME") {
                    let parts: Vec<&str> = tail_line.split_whitespace().collect();
                    for (i, element) in parts.iter().enumerate() {
                        if *element == "NAME"
                            && let Some(val) = parts.get(i + 1)
                        {
                            trad_name = val.to_string();
                        }
                    }
                }
            }
            break;
        }
        match n {
            2 => trad_name = line.to_string(),
            5 => {
                let tokens: Vec<&str> = line.split(' ').collect();
                if tokens.len() > 3 {
                    nomen_name = format!("{}{}", tokens[2], tokens[3]);
                }
            }
            7 => {
                let icrs: Vec<&str> = line.split(' ').collect();
                if icrs.len() > 7 {
                    ra_hour = Some(icrs[1].to_string());
                    ra_minute = Some(icrs[2].to_string());
                    ra_sec = Some(icrs[3].to_string());
                    dec_degree = Some(icrs[5].to_string());
                    dec_minute = Some(icrs[6].to_string());
                    dec_sec = Some(icrs[7].to_string());
                }
            }
            11 => {
                let pm: Vec<&str> = line.split(' ').collect();
                if pm.len() > 3 {
                    pmra = Some(pm[2].to_string());
                    pmde = Some(pm[3].to_string());
                }
            }
            12 => {
                let para: Vec<&str> = line.split(' ').collect();
                if para.len() > 1 {
                    parallax = Some(para[1].to_string());
                }
            }
            13 => {
                let rv: Vec<&str> = line.split(' ').collect();
                if rv.len() > 2 {
                    rad_vel = Some(rv[2].to_string());
                }
            }
            _ => {}
        }
        if line.contains("Flux V") {
            let flux: Vec<&str> = line.split(' ').collect();
            if flux.len() > 3 {
                mag_v = Some(flux[3].to_string());
            }
        }
    }

    Some(SimbadResult {
        trad_name,
        nomen_name,
        hip_id,
        ra_hour: ra_hour?,
        ra_minute: ra_minute?,
        ra_sec: ra_sec?,
        dec_degree: dec_degree?,
        dec_minute: dec_minute?,
        dec_sec: dec_sec?,
        pmra: pmra?,
        pmde: pmde?,
        rad_vel: rad_vel?,
        parallax: parallax?,
        mag_v: mag_v.unwrap_or_else(|| "0".to_string()),
    })
}

fn build_entry_lines(parsed: &SimbadResult) -> Vec<String> {
    let nomen = &parsed.nomen_name;
    let long_form = nomen_to_long_form(nomen);
    let data_fields = format!(
        "{nomen},ICRS,{},{},{},{},{},{},{},{},{},{},{}",
        parsed.ra_hour,
        parsed.ra_minute,
        parsed.ra_sec,
        parsed.dec_degree,
        parsed.dec_minute,
        parsed.dec_sec,
        parsed.pmra,
        parsed.pmde,
        parsed.rad_vel,
        parsed.parallax,
        parsed.mag_v
    );

    let mut comment = format!("#0# {nomen}, {long_form}");
    let extra_ids: Vec<&str> = [parsed.trad_name.as_str(), parsed.hip_id.as_str()]
        .into_iter()
        .filter(|id| !id.is_empty() && *id != "no hip id")
        .collect();
    for id in &extra_ids {
        write!(comment, ", {id}").unwrap();
    }
    comment.push('\n');

    let mut lines = vec![comment];
    lines.push(format!("{nomen},{data_fields}\n"));
    lines.push(format!("{long_form},{data_fields}\n"));

    for id in &extra_ids {
        lines.push(format!("{id},{data_fields}\n"));
    }

    lines
}

fn validate_entry(lines: &[&str], star_name: &str) -> bool {
    for line in lines {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.trim().split(',').collect();
        if fields.len() < EXPECTED_FIELD_COUNT {
            eprintln!(
                "Warning: entry for '{star_name}' has {} fields (expected {EXPECTED_FIELD_COUNT}), skipping",
                fields.len()
            );
            return false;
        }
    }
    true
}

fn load_existing_names(path: &str) -> HashSet<String> {
    let mut names = HashSet::new();
    let content = match fs::read_to_string(path) {
        Ok(c) => c,
        Err(_) => return names,
    };
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if let Some(first) = line.split(',').next() {
            names.insert(first.to_owned());
        }
    }
    names
}

fn read_star_list(path: &str) -> io::Result<Vec<String>> {
    let content = fs::read_to_string(path)?;
    Ok(content
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(String::from)
        .collect())
}

fn verify_against_simbad(nomen_names: &[String], sefstars_path: &str) {
    let content = match fs::read_to_string(sefstars_path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("File not found: {sefstars_path}: {e}");
            return;
        }
    };

    let mut existing: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.trim().split(',').map(String::from).collect();
        if fields.len() >= EXPECTED_FIELD_COUNT {
            existing.insert(fields[0].clone(), fields);
        }
    }

    for nomen in nomen_names {
        let file_fields = match existing.get(nomen.as_str()) {
            Some(f) => f,
            None => {
                println!("  {nomen}: not found in {sefstars_path}");
                continue;
            }
        };

        let parsed = match query_simbad(nomen) {
            Ok(p) => p,
            Err(e) => {
                println!("  {nomen}: SIMBAD error: {e}");
                continue;
            }
        };

        let simbad_ra = format!("{},{},{}", parsed.ra_hour, parsed.ra_minute, parsed.ra_sec);
        let simbad_dec = format!(
            "{},{},{}",
            parsed.dec_degree, parsed.dec_minute, parsed.dec_sec
        );
        let file_ra = format!("{},{},{}", file_fields[3], file_fields[4], file_fields[5]);
        let file_dec = format!("{},{},{}", file_fields[6], file_fields[7], file_fields[8]);

        if simbad_ra == file_ra && simbad_dec == file_dec {
            println!("  {nomen}: OK (coordinates match)");
        } else {
            println!("  {nomen}: MISMATCH");
            println!("    file:   RA={file_ra}  Dec={file_dec}");
            println!("    simbad: RA={simbad_ra}  Dec={simbad_dec}");
        }

        let simbad_mag = &parsed.mag_v;
        let file_mag = file_fields.get(13).map_or("?", |s| s.as_str());
        if simbad_mag != file_mag {
            println!("    mag V: file={file_mag} simbad={simbad_mag}");
        }
    }
}

fn greek_to_long(abbr: &str) -> Option<&'static str> {
    match abbr {
        "alf" => Some("Alpha"),
        "bet" => Some("Beta"),
        "gam" | "g" => Some("Gamma"),
        "del" | "d" => Some("Delta"),
        "eps" => Some("Epsilon"),
        "zet" => Some("Zeta"),
        "eta" => Some("Eta"),
        "tet" => Some("Theta"),
        "iot" => Some("Iota"),
        "kap" => Some("Kappa"),
        "lam" => Some("Lambda"),
        "mu." => Some("Mu"),
        "nu." => Some("Nu"),
        "ksi" => Some("Xi"),
        "omi" => Some("Omicron"),
        "pi." => Some("Pi"),
        "rho" => Some("Rho"),
        "sig" => Some("Sigma"),
        "tau" => Some("Tau"),
        "ups" => Some("Upsilon"),
        "phi" => Some("Phi"),
        "chi" => Some("Chi"),
        "psi" => Some("Psi"),
        "ome" => Some("Omega"),
        _ => None,
    }
}

fn constellation_to_long(abbr: &str) -> Option<&'static str> {
    match abbr {
        "Ari" => Some("Arietis"),
        "Tau" => Some("Tauri"),
        "Gem" => Some("Geminorum"),
        "Cnc" => Some("Cancri"),
        "Leo" => Some("Leonis"),
        "Vir" => Some("Virginis"),
        "Lib" => Some("Librae"),
        "Sco" => Some("Scorpii"),
        "Oph" => Some("Ophiuci"),
        "Sgr" => Some("Sagittarii"),
        "Cap" => Some("Capricorni"),
        "Aqr" => Some("Aquarii"),
        "And" => Some("Andromedae"),
        "Ant" => Some("Antliae"),
        "Aps" => Some("Apodis"),
        "Ara" => Some("Arae"),
        "Psc" => Some("Piscium"),
        "Eri" => Some("Eridani"),
        "Cae" => Some("Caeli"),
        "Cam" => Some("Camelopardalis"),
        "Cas" => Some("Cassiopeiae"),
        "Cen" => Some("Centauri"),
        "Cep" => Some("Cephei"),
        "UMa" => Some("Ursae Majoris"),
        "UMi" => Some("Ursae Minoris"),
        "Aql" => Some("Aquilae"),
        "Hyd" => Some("Hydrae"),
        "Sct" => Some("Scuti"),
        "Sex" => Some("Sextantis"),
        "Sge" => Some("Sagittae"),
        "Boo" => Some("Bootis"),
        "Dra" => Some("Draconis"),
        "Del" => Some("Delphini"),
        "Dor" => Some("Doradus"),
        "Equ" => Some("Equulei"),
        "For" => Some("Fornacis"),
        "Cyg" => Some("Cygni"),
        "Gru" => Some("Gruis"),
        "Ori" => Some("Orionis"),
        "Cet" => Some("Ceti"),
        "Cha" => Some("Chamaeleontis"),
        "Cir" => Some("Circini"),
        "Col" => Some("Columbae"),
        "Com" => Some("Comae Berenices"),
        "CrB" => Some("Coronae Borealis"),
        "CrA" => Some("Coronae Australis"),
        "TCrB" => Some("TCoronae Borealis"),
        "Crt" => Some("Crateris"),
        "Cru" => Some("Crucis"),
        "Crv" => Some("Corvi"),
        "CVn" => Some("Canum Venaticorum"),
        "CMa" => Some("Canis Majoris"),
        "CMi" => Some("Canis Minoris"),
        "Aur" => Some("Aurigae"),
        "Car" => Some("Carinae"),
        "Lyr" => Some("Lyrae"),
        "Lep" => Some("Leporis"),
        "Men" => Some("Mensae"),
        "Mic" => Some("Microscopii"),
        "Mon" => Some("Monocerotis"),
        "Mus" => Some("Muscae"),
        "Nor" => Some("Normae"),
        "Oct" => Some("Octantis"),
        "Ind" => Some("Indi"),
        "Pav" => Some("Pavonis"),
        "Peg" => Some("Pegasi"),
        "Phe" => Some("Phoenicis"),
        "LMi" => Some("Leonis Minoris"),
        "Lup" => Some("Lupi"),
        "Lyn" => Some("Lyncis"),
        "Ser" => Some("Serpentis"),
        "Tel" => Some("Telescopii"),
        "TrA" => Some("Trianguli Australis"),
        "Tri" => Some("Trianguli"),
        "Tuc" => Some("Tucanae"),
        "Her" => Some("Herculis"),
        "Hor" => Some("Horologii"),
        "Hya" => Some("Hydrae"),
        "Hyi" => Some("Hydri"),
        "Lac" => Some("Lacertae"),
        "Per" => Some("Persei"),
        "Pic" => Some("Pictoris"),
        "PsA" => Some("Piscis Austrini"),
        "Pup" => Some("Puppis"),
        "Pyx" => Some("Pyxidis"),
        "Ret" => Some("Reticuli"),
        "Scl" => Some("Sculptoris"),
        "Vel" => Some("Velorum"),
        "Vol" => Some("Volantis"),
        "Vul" => Some("Vulpeculae"),
        "VC" => Some("Virgo Cluster"),
        "M" => Some("Messier Object"),
        "NGC" => Some("New General Catalogue"),
        "HIP" => Some("Hipparcos Catalogue"),
        "HR" => Some("Bright Star Catalogue"),
        "HD" => Some("Henry Draper Catalogue"),
        _ => None,
    }
}

fn nomen_to_long_form(nomen: &str) -> String {
    let nomen = nomen.strip_prefix(',').unwrap_or(nomen);

    // Flamsteed: starts with digits, e.g. "48Lib"
    if nomen.starts_with(|c: char| c.is_ascii_digit()) {
        let split = nomen
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(nomen.len());
        let number = &nomen[..split];
        let constellation = &nomen[split..];
        if let Some(long) = constellation_to_long(constellation) {
            return format!("{long}{number}");
        }
        return nomen.to_string();
    }

    // Special catalogue prefixes: VC, HD, HR, HIP, NGC, M
    let special = ["NGC", "HIP", "HR", "HD", "VC"];
    for prefix in special {
        if let Some(rest) = nomen.strip_prefix(prefix)
            && let Some(long) = constellation_to_long(prefix)
        {
            if rest.is_empty() {
                return long.to_string();
            }
            return format!("{long} {}", rest.trim());
        }
    }
    // Messier: starts with M but not "mu."
    if nomen.starts_with('M') && !nomen.starts_with("mu.") {
        let rest = &nomen[1..];
        if let Some(long) = constellation_to_long("M") {
            return format!("{long} {}", rest.trim());
        }
    }

    // Bayer: 3-char greek + constellation, optionally with a number suffix on the greek
    // e.g. "alfTau", "eps01Ori"
    if nomen.len() >= 4 {
        let greek_abbr = &nomen[..3];
        let remainder = &nomen[3..];

        // Check for a numeric suffix between greek and constellation (e.g. "01" in "eps01Ori")
        let num_end = remainder
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(remainder.len());
        let number = &remainder[..num_end];
        let constellation = &remainder[num_end..];

        if greek_abbr.chars().all(|c| c.is_lowercase() || c == '.')
            && let Some(greek_long) = greek_to_long(greek_abbr)
            && let Some(const_long) = constellation_to_long(constellation)
        {
            if number.is_empty() {
                return format!("{greek_long} {const_long}");
            }
            return format!("{greek_long} {const_long} {number}");
        }
    }

    nomen.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bayer_designation() {
        assert_eq!(nomen_to_long_form("alfTau"), "Alpha Tauri");
        assert_eq!(nomen_to_long_form("betOri"), "Beta Orionis");
        assert_eq!(nomen_to_long_form("gamSgr"), "Gamma Sagittarii");
        assert_eq!(nomen_to_long_form("zetUMa"), "Zeta Ursae Majoris");
    }

    #[test]
    fn bayer_with_number() {
        assert_eq!(nomen_to_long_form("eps01Ori"), "Epsilon Orionis 01");
    }

    #[test]
    fn flamsteed_designation() {
        assert_eq!(nomen_to_long_form("48Lib"), "Librae48");
        assert_eq!(nomen_to_long_form("61Cyg"), "Cygni61");
    }

    #[test]
    fn hip_designation() {
        assert_eq!(nomen_to_long_form("HIP12345"), "Hipparcos Catalogue 12345");
    }

    #[test]
    fn messier_object() {
        assert_eq!(nomen_to_long_form("M31"), "Messier Object 31");
    }

    #[test]
    fn ngc_catalogue() {
        assert_eq!(nomen_to_long_form("NGC1234"), "New General Catalogue 1234");
    }

    #[test]
    fn leading_comma_stripped() {
        assert_eq!(nomen_to_long_form(",alfTau"), "Alpha Tauri");
    }

    #[test]
    fn unknown_passes_through() {
        assert_eq!(nomen_to_long_form("Foobar"), "Foobar");
    }

    #[test]
    fn mu_dot_not_messier() {
        assert_eq!(nomen_to_long_form("mu.Sgr"), "Mu Sagittarii");
    }

    #[test]
    fn validate_good_entry() {
        let line = "alfTau,alfTau,ICRS,4,35,55.2,16,30,33,62,-189,54,48,0.85\n";
        assert!(validate_entry(&[line], "alfTau"));
    }

    #[test]
    fn validate_short_entry() {
        let line = "alfTau,alfTau,ICRS,4,35\n";
        assert!(!validate_entry(&[line], "alfTau"));
    }

    #[test]
    fn comments_skip_validation() {
        assert!(validate_entry(&["# this is a comment\n"], "test"));
    }
}
