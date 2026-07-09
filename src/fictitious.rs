// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Ninth House Studios LLC

//! Fictitious / hypothetical planet elements and orbital mechanics — Uranian
//! (Hamburg School) bodies, Waldemath Black Moon, and custom `seorbel.txt` entries.

use std::path::Path;

use crate::constants::*;
use crate::error::Error;
use crate::flags::CalcFlags;
use crate::math::{normalize_degrees, normalize_radians, rotate_x};
use crate::obliquity::obliquity;
use crate::precession::precess;
use crate::types::{AstroModels, PrecessionDirection};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const KGAUSS_GEO: f64 = 0.0000298122353216;

// ---------------------------------------------------------------------------
// Built-in table — Neely-revised Uranian + 7 shared rows (§1)
// ---------------------------------------------------------------------------

struct BuiltinRow {
    epoch: f64,
    equinox: f64,
    mano: f64,
    sema: f64,
    ecce: f64,
    parg: f64,
    node: f64,
    incl: f64,
    name: &'static str,
}

// Only the Neely (SE_NEELY) set — the non-Neely #else branch is dead code in
// every standard C build and has no golden-test signal to validate against.
const BUILTIN_TABLE: [BuiltinRow; NFICT_ELEM as usize] = [
    // 0..7: Witte/Sieggruen planets, refined by James Neely
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 163.7409,
        sema: 40.99837,
        ecce: 0.00460,
        parg: 171.4333,
        node: 129.8325,
        incl: 1.0833,
        name: "Cupido",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 27.6496,
        sema: 50.66744,
        ecce: 0.00245,
        parg: 148.1796,
        node: 161.3339,
        incl: 1.0500,
        name: "Hades",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 165.1232,
        sema: 59.21436,
        ecce: 0.00120,
        parg: 299.0440,
        node: 0.0000,
        incl: 0.0000,
        name: "Zeus",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 169.0193,
        sema: 64.81960,
        ecce: 0.00305,
        parg: 208.8801,
        node: 0.0000,
        incl: 0.0000,
        name: "Kronos",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 138.0533,
        sema: 70.29949,
        ecce: 0.00000,
        parg: 0.0000,
        node: 0.0000,
        incl: 0.0000,
        name: "Apollon",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 351.3350,
        sema: 73.62765,
        ecce: 0.00000,
        parg: 0.0000,
        node: 0.0000,
        incl: 0.0000,
        name: "Admetos",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 55.8983,
        sema: 77.25568,
        ecce: 0.00000,
        parg: 0.0000,
        node: 0.0000,
        incl: 0.0000,
        name: "Vulkanus",
    },
    BuiltinRow {
        epoch: J1900,
        equinox: J1900,
        mano: 165.5163,
        sema: 83.66907,
        ecce: 0.00000,
        parg: 0.0000,
        node: 0.0000,
        incl: 0.0000,
        name: "Poseidon",
    },
    // 8..14: non-Neely-variant rows (same in both builds)
    BuiltinRow {
        epoch: 2368547.66,
        equinox: 2431456.5,
        mano: 0.0,
        sema: 77.775,
        ecce: 0.3,
        parg: 0.7,
        node: 0.0,
        incl: 0.0,
        name: "Isis-Transpluto",
    },
    BuiltinRow {
        epoch: 1856113.380954,
        equinox: 1856113.380954,
        mano: 0.0,
        sema: 234.8921,
        ecce: 0.981092,
        parg: 103.966,
        node: -44.567,
        incl: 158.708,
        name: "Nibiru",
    },
    BuiltinRow {
        epoch: 2374696.5,
        equinox: J2000,
        mano: 0.0,
        sema: 101.2,
        ecce: 0.411,
        parg: 208.5,
        node: 275.4,
        incl: 32.4,
        name: "Harrington",
    },
    BuiltinRow {
        epoch: 2395662.5,
        equinox: 2395662.5,
        mano: 34.05,
        sema: 36.15,
        ecce: 0.10761,
        parg: 284.75,
        node: 0.0,
        incl: 0.0,
        name: "Leverrier",
    },
    BuiltinRow {
        epoch: 2395662.5,
        equinox: 2395662.5,
        mano: 24.28,
        sema: 37.25,
        ecce: 0.12062,
        parg: 299.11,
        node: 0.0,
        incl: 0.0,
        name: "Adams",
    },
    BuiltinRow {
        epoch: 2425977.5,
        equinox: 2425977.5,
        mano: 281.0,
        sema: 43.0,
        ecce: 0.202,
        parg: 204.9,
        node: 0.0,
        incl: 0.0,
        name: "Lowell",
    },
    BuiltinRow {
        epoch: 2425977.5,
        equinox: 2425977.5,
        mano: 48.95,
        sema: 55.1,
        ecce: 0.31,
        parg: 280.1,
        node: 100.0,
        incl: 15.0,
        name: "Pickering",
    },
];

// ---------------------------------------------------------------------------
// T-term polynomial expression (§3.1)
// ---------------------------------------------------------------------------

/// Pre-parsed representation of a check_t_terms expression.
/// Evaluates to: `constant + coeffs[0]*T + coeffs[1]*T^2 + coeffs[2]*T^3 + coeffs[3]*T^4`
/// where T = t_days / 36525.
#[derive(Clone, Debug)]
struct ElementExpr {
    constant: f64,
    coeffs: [f64; 4],
    has_t_terms: bool,
}

impl ElementExpr {
    fn from_constant(val: f64) -> Self {
        Self {
            constant: val,
            coeffs: [0.0; 4],
            has_t_terms: false,
        }
    }

    fn eval(&self, t_days: f64) -> f64 {
        let t = t_days / 36525.0;
        self.constant
            + self.coeffs[0] * t
            + self.coeffs[1] * t * t
            + self.coeffs[2] * t * t * t
            + self.coeffs[3] * t * t * t * t
    }
}

/// Parse a check_t_terms expression string into an ElementExpr.
///
/// The C algorithm (swemplan.c:916-967) accumulates terms as running products
/// of numeric literals and T-power lookups, committed on `+`, `-`, or end.
///
/// C's tt[] array: tt[0]=tt[1]=T^1, tt[2]=T^2, tt[3]=T^3, tt[4]=T^4.
/// The "T0" quirk (tt[0]=T^1, not T^0=1) is preserved by mapping indices 0,1
/// both to T^1 in the polynomial output.
fn parse_t_terms(input: &str) -> Result<ElementExpr, Error> {
    let has_t_terms = input.contains('+') || input.contains('-');

    let mut constant = 0.0_f64;
    let mut coeffs = [0.0_f64; 4];
    let mut fac = 1.0_f64;
    let mut t_power = 0_usize; // accumulated T power for current term (0 = constant)
    let mut z = 0_usize;

    let bytes = input.as_bytes();
    let mut i = 0;

    loop {
        // Skip spaces/tabs
        while i < bytes.len() && (bytes[i] == b' ' || bytes[i] == b'\t') {
            i += 1;
        }

        if i >= bytes.len() || bytes[i] == b'+' || bytes[i] == b'-' {
            // Commit current term
            if z > 0 {
                if t_power == 0 {
                    constant += fac;
                } else if (1..=4).contains(&t_power) {
                    coeffs[t_power - 1] += fac;
                }
            }

            if i >= bytes.len() {
                break;
            }

            let sign = if bytes[i] == b'-' { -1.0 } else { 1.0 };
            fac = sign;
            t_power = 0;
            i += 1;
        } else {
            // Skip multiplication signs and whitespace
            while i < bytes.len() && (bytes[i] == b'*' || bytes[i] == b' ' || bytes[i] == b'\t') {
                i += 1;
            }

            if i < bytes.len() && (bytes[i] == b't' || bytes[i] == b'T') {
                i += 1;
                if i >= bytes.len()
                    || bytes[i] == b'+'
                    || bytes[i] == b'-'
                    || bytes[i] == b' '
                    || bytes[i] == b'\t'
                    || bytes[i] == b'*'
                {
                    // Bare T → tt[0] = T^1
                    t_power += 1;
                } else {
                    // Tn → tt[n], mapped to effective power
                    let start = i;
                    while i < bytes.len() && bytes[i].is_ascii_digit() {
                        i += 1;
                    }
                    let idx: usize = std::str::from_utf8(&bytes[start..i])
                        .unwrap_or("0")
                        .parse()
                        .unwrap_or(0);
                    // C quirk: tt[0]=tt[1]=T^1, tt[2]=T^2, tt[3]=T^3, tt[4]=T^4
                    let power = if idx <= 1 { 1 } else { idx };
                    if power > 4 {
                        return Err(Error::FileFormat(format!(
                            "T-term power {idx} out of range in expression"
                        )));
                    }
                    t_power += power;
                }
            } else if i < bytes.len() {
                // Numeric literal
                let start = i;
                while i < bytes.len()
                    && (bytes[i].is_ascii_digit()
                        || bytes[i] == b'.'
                        || bytes[i] == b'e'
                        || bytes[i] == b'E'
                        || ((bytes[i] == b'+' || bytes[i] == b'-')
                            && i > start
                            && (bytes[i - 1] == b'e' || bytes[i - 1] == b'E')))
                {
                    i += 1;
                }
                let num_str = std::str::from_utf8(&bytes[start..i]).unwrap_or("0");
                let val: f64 = num_str.parse().map_err(|_| {
                    Error::FileFormat(format!("invalid number '{num_str}' in T-term expression"))
                })?;
                fac *= val;
            }
        }
        z += 1;
    }

    Ok(ElementExpr {
        constant,
        coeffs,
        has_t_terms,
    })
}

// ---------------------------------------------------------------------------
// Equinox representation
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Equinox {
    Fixed(f64),
    Date,
}

// ---------------------------------------------------------------------------
// Parsed orbital element row
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct ParsedRow {
    epoch: f64,
    equinox: Equinox,
    mano: ElementExpr,
    sema: ElementExpr,
    ecce: ElementExpr,
    parg: ElementExpr,
    node: ElementExpr,
    incl: ElementExpr,
    name: String,
    is_geo: bool,
    /// Built-in rows skip degnorm on angle columns — C's built-in path does
    /// plain `value * DEGTORAD` without `swe_degnorm` (swemplan.c:720-730).
    from_builtin: bool,
}

// ---------------------------------------------------------------------------
// Resolved elements (after T-term evaluation for a specific tjd)
// ---------------------------------------------------------------------------

/// Fictitious body's Keplerian orbital elements after T-term evaluation at a
/// specific Julian day.
pub struct ResolvedElements {
    /// Reference epoch (Julian day) the elements are given for.
    pub tjd0: f64,
    /// Equinox (Julian day) the elements are referred to.
    pub tequ: f64,
    /// Mean anomaly at `tjd0` (degrees).
    pub mano: f64,
    /// Semi-major axis (AU).
    pub sema: f64,
    /// Orbital eccentricity.
    pub ecce: f64,
    /// Argument (or longitude) of perihelion (degrees).
    pub parg: f64,
    /// Longitude of the ascending node (degrees).
    pub node: f64,
    /// Orbital inclination (degrees).
    pub incl: f64,
    /// Body name as given in the catalog row.
    pub name: String,
    /// Whether the elements are geocentric (vs. heliocentric).
    pub is_geo: bool,
}

// ---------------------------------------------------------------------------
// Fictitious catalog
// ---------------------------------------------------------------------------

/// Catalog of fictitious/hypothetical planet orbital elements, either the
/// built-in Neely-revised Uranian table or one parsed from a `seorbel.txt`-style file.
pub struct FictitiousCatalog {
    rows: Vec<ParsedRow>,
    from_file: bool,
}

impl FictitiousCatalog {
    /// Constructs the catalog from the built-in Neely-revised Uranian planet
    /// table (no file I/O).
    pub fn builtin() -> Self {
        let rows = BUILTIN_TABLE
            .iter()
            .map(|r| ParsedRow {
                epoch: r.epoch,
                equinox: Equinox::Fixed(r.equinox),
                mano: ElementExpr::from_constant(r.mano),
                sema: ElementExpr::from_constant(r.sema),
                ecce: ElementExpr::from_constant(r.ecce),
                parg: ElementExpr::from_constant(r.parg),
                node: ElementExpr::from_constant(r.node),
                incl: ElementExpr::from_constant(r.incl),
                name: r.name.to_string(),
                is_geo: false,
                from_builtin: true,
            })
            .collect();
        Self {
            rows,
            from_file: false,
        }
    }

    /// Number of catalog entries.
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Whether this catalog was loaded from an external file (vs. the built-in table).
    pub fn from_file(&self) -> bool {
        self.from_file
    }
}

/// Parse epoch or equinox field: J2000/J1900/B1950 sentinels, JDATE (equinox only),
/// or a literal JD number.
fn parse_epoch_equinox(field: &str, allow_jdate: bool) -> Result<Equinox, Error> {
    let trimmed = field.trim();
    if trimmed.len() >= 5 {
        let lower: String = trimmed[..5].to_ascii_lowercase();
        if lower == "j2000" {
            return Ok(Equinox::Fixed(J2000));
        }
        if lower == "b1950" {
            return Ok(Equinox::Fixed(B1950));
        }
        if lower == "j1900" {
            return Ok(Equinox::Fixed(J1900));
        }
        if allow_jdate && trimmed.len() >= 5 {
            let jdate_lower: String = trimmed[..5].to_ascii_lowercase();
            if jdate_lower == "jdate" {
                return Ok(Equinox::Date);
            }
        }
    }
    // Check for invalid j/b prefixed values
    let first = trimmed.as_bytes().first().copied().unwrap_or(0);
    if first == b'j' || first == b'J' || first == b'b' || first == b'B' {
        return Err(Error::FileFormat(format!(
            "invalid epoch/equinox '{trimmed}'"
        )));
    }
    let val: f64 = trimmed
        .parse()
        .map_err(|_| Error::FileFormat(format!("invalid epoch/equinox number '{trimmed}'")))?;
    Ok(Equinox::Fixed(val))
}

/// Parse seorbel.txt contents into a catalog.
fn parse_orbel_file(contents: &str) -> Result<FictitiousCatalog, Error> {
    let mut rows = Vec::new();

    for line in contents.lines() {
        // Strip leading whitespace
        let trimmed = line.trim_start();

        // Skip blank lines and comments
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with('\r')
            || trimmed.starts_with('\n')
        {
            continue;
        }

        // Truncate at inline comment
        let content = if let Some(pos) = trimmed.find('#') {
            &trimmed[..pos]
        } else {
            trimmed
        };

        // Split on commas
        let fields: Vec<&str> = content.split(',').collect();
        if fields.len() < 9 {
            return Err(Error::FileFormat(format!(
                "nine elements required, got {} in line: {content}",
                fields.len()
            )));
        }

        // Column 1: epoch
        let epoch_eq = parse_epoch_equinox(fields[0], false)?;
        let epoch = match epoch_eq {
            Equinox::Fixed(v) => v,
            Equinox::Date => {
                return Err(Error::FileFormat("JDATE not allowed for epoch".to_string()));
            }
        };

        // Column 2: equinox
        let equinox = parse_epoch_equinox(fields[1], true)?;

        // Columns 3-8: elements via check_t_terms
        let mano = parse_t_terms(fields[2].trim()).map_err(|_| {
            Error::FileFormat(format!(
                "mean anomaly value invalid: '{}'",
                fields[2].trim()
            ))
        })?;
        let sema = parse_t_terms(fields[3].trim()).map_err(|_| {
            Error::FileFormat(format!("semi-axis value invalid: '{}'", fields[3].trim()))
        })?;
        let ecce = parse_t_terms(fields[4].trim()).map_err(|_| {
            Error::FileFormat(format!("eccentricity invalid: '{}'", fields[4].trim()))
        })?;
        let parg = parse_t_terms(fields[5].trim()).map_err(|_| {
            Error::FileFormat(format!(
                "perihelion argument value invalid: '{}'",
                fields[5].trim()
            ))
        })?;
        let node = parse_t_terms(fields[6].trim()).map_err(|_| {
            Error::FileFormat(format!("node value invalid: '{}'", fields[6].trim()))
        })?;
        let incl = parse_t_terms(fields[7].trim()).map_err(|_| {
            Error::FileFormat(format!("inclination value invalid: '{}'", fields[7].trim()))
        })?;

        // Column 9: name
        let name = fields[8].trim().to_string();

        // Column 10 (optional): geo flag
        let is_geo = if fields.len() > 9 {
            fields[9].to_ascii_lowercase().contains("geo")
        } else {
            false
        };

        rows.push(ParsedRow {
            epoch,
            equinox,
            mano,
            sema,
            ecce,
            parg,
            node,
            incl,
            name,
            is_geo,
            from_builtin: false,
        });
    }

    Ok(FictitiousCatalog {
        rows,
        from_file: true,
    })
}

/// Load fictitious-body catalog from `ephe_path/seorbel.txt`.
/// Falls back to the built-in 15-row table only if the file is absent.
/// A present but malformed file returns an error — matching C, which only
/// falls back when the file cannot be opened (swemplan.c:707-713).
pub fn load_fictitious_catalog(ephe_path: Option<&Path>) -> Result<FictitiousCatalog, Error> {
    let Some(dir) = ephe_path else {
        return Ok(FictitiousCatalog::builtin());
    };
    let path = dir.join(FICTFILE);
    let contents = match std::fs::read_to_string(&path) {
        Ok(s) => s,
        Err(_) => return Ok(FictitiousCatalog::builtin()),
    };
    parse_orbel_file(&contents)
}

// ---------------------------------------------------------------------------
// Element resolution
// ---------------------------------------------------------------------------

/// Resolve orbital elements for a given body row index and Julian day.
///
/// `ipl` is the 0-based row index (already `raw_id - SE_FICT_OFFSET`).
/// Returns fully evaluated elements with angles in radians.
pub fn resolve_elements(
    catalog: &FictitiousCatalog,
    ipl: usize,
    tjd: f64,
) -> Result<ResolvedElements, Error> {
    if ipl >= catalog.rows.len() {
        if !catalog.from_file && ipl >= NFICT_ELEM as usize {
            return Err(Error::FileFormat(format!(
                "error no elements for fictitious body no {:.0}",
                ipl as f64
            )));
        }
        return Err(Error::FileFormat(format!(
            "elements for planet {:.0} not found",
            ipl as f64
        )));
    }

    let row = &catalog.rows[ipl];

    // Epoch
    let mut tjd0 = row.epoch;
    let tt = tjd - tjd0;

    // Equinox
    let tequ = match &row.equinox {
        Equinox::Fixed(v) => *v,
        Equinox::Date => tjd,
    };

    // Mean anomaly — if T-terms present, override epoch to tjd (§3 special case)
    let mano_raw = row.mano.eval(tt);
    if row.mano.has_t_terms {
        tjd0 = tjd;
    }
    // C's built-in path (swemplan.c:720-730) does plain `value * DEGTORAD`
    // without swe_degnorm; the file path applies degnorm then DEGTORAD.
    let to_rad = |val: f64| -> f64 {
        if row.from_builtin {
            val * DEGTORAD
        } else {
            normalize_degrees(val) * DEGTORAD
        }
    };
    let mano = to_rad(mano_raw);

    // Semi-axis
    let sema = row.sema.eval(tt);
    if sema <= 0.0 {
        return Err(Error::FileFormat(format!(
            "semi-axis value invalid ({sema})"
        )));
    }

    // Eccentricity
    let ecce = row.ecce.eval(tt);
    if !(0.0..1.0).contains(&ecce) {
        return Err(Error::FileFormat(format!(
            "eccentricity invalid ({ecce}), no parabolic or hyperbolic orbits allowed"
        )));
    }

    // Angular elements
    let parg = to_rad(row.parg.eval(tt));
    let node = to_rad(row.node.eval(tt));
    let incl = to_rad(row.incl.eval(tt));

    Ok(ResolvedElements {
        tjd0,
        tequ,
        mano,
        sema,
        ecce,
        parg,
        node,
        incl,
        name: row.name.clone(),
        is_geo: row.is_geo,
    })
}

/// Look up the human-readable name for a fictitious body.
pub fn fictitious_name(catalog: &FictitiousCatalog, ipl: usize) -> String {
    if ipl < catalog.rows.len() {
        catalog.rows[ipl].name.clone()
    } else {
        "name not found".to_string()
    }
}

// ---------------------------------------------------------------------------
// Kepler equation solver (swephlib.c:4065-4096)
// ---------------------------------------------------------------------------

/// Solve the Kepler equation E - e*sin(E) = M for the eccentric anomaly E.
///
/// Port of swi_kepler: fixed-point iteration for e < 0.4, Newton's method
/// for e >= 0.4. Convergence tolerance 1e-12 radians, no iteration cap.
pub fn kepler(mut e: f64, m: f64, ecce: f64) -> f64 {
    let mut de = 1.0_f64;
    if ecce < 0.4 {
        while de > 1e-12 {
            let e0 = e;
            e = m + ecce * e0.sin();
            de = (e - e0).abs();
        }
    } else {
        while de > 1e-12 {
            let e0 = e;
            let x = (m + ecce * e0.sin() - e0) / (1.0 - ecce * e0.cos());
            de = x.abs();
            if de < 1e-2 {
                e = e0 + x;
            } else {
                e = normalize_radians(e0 + x);
                de = (e - e0).abs();
            }
        }
    }
    e
}

// ---------------------------------------------------------------------------
// swi_osc_el_plan — elements → J2000 equatorial barycentric state vector (§4)
// ---------------------------------------------------------------------------

/// Compute heliocentric (or geocentric) J2000-equatorial-barycentric position
/// and velocity for a fictitious body from its Keplerian orbital elements.
///
/// `ipl` is the 0-based row index (already `raw_id - SE_FICT_OFFSET`).
/// `xearth` / `xsun` are 6-vectors (position + velocity) for barycentric Earth
/// and Sun respectively, already computed by the caller.
pub fn osc_el_plan(
    tjd: f64,
    catalog: &FictitiousCatalog,
    ipl: usize,
    xearth: &[f64; 6],
    xsun: &[f64; 6],
    models: &AstroModels,
) -> Result<[f64; 6], Error> {
    let elem = resolve_elements(catalog, ipl, tjd)?;

    // §4.2 Daily motion (deg/day → rad/day)
    // Evaluation order preserved: ((constant * DEGTORAD) / sema) / sqrt(sema)
    let sema_sqrt = elem.sema.sqrt();
    let mut dmot = (0.9856076686 * DEGTORAD) / elem.sema / sema_sqrt;
    if elem.is_geo {
        dmot /= SUN_EARTH_MRAT.sqrt();
    }

    // Gaussian constant
    let k = if elem.is_geo {
        KGAUSS_GEO / sema_sqrt
    } else {
        KGAUSS / sema_sqrt
    };

    // §4.3 Gaussian P/Q/R rotation vectors
    let cosnode = elem.node.cos();
    let sinnode = elem.node.sin();
    let cosincl = elem.incl.cos();
    let sinincl = elem.incl.sin();
    let cosparg = elem.parg.cos();
    let sinparg = elem.parg.sin();

    let pqr = [
        cosparg * cosnode - sinparg * cosincl * sinnode, // P.x
        -sinparg * cosnode - cosparg * cosincl * sinnode, // Q.x
        sinincl * sinnode,                               // R.x (unused)
        cosparg * sinnode + sinparg * cosincl * cosnode, // P.y
        -sinparg * sinnode + cosparg * cosincl * cosnode, // Q.y
        -sinincl * cosnode,                              // R.y (unused)
        sinparg * sinincl,                               // P.z
        cosparg * sinincl,                               // Q.z
        cosincl,                                         // R.z (unused)
    ];

    // §4.4 Kepler equation
    let m = normalize_radians(elem.mano + (tjd - elem.tjd0) * dmot);
    let mut e_anom = m;

    // High-eccentricity initial-guess refinement (only Nibiru, ecce=0.981092)
    if elem.ecce > 0.975 {
        let mut m2 = e_anom * RADTODEG;
        let m_180_or_0 = if m2 > 150.0 && m2 < 210.0 {
            m2 -= 180.0;
            180.0
        } else {
            0.0
        };
        if m2 > 330.0 {
            m2 -= 360.0;
        }
        let msgn;
        if m2 < 0.0 {
            m2 = -m2;
            msgn = -1.0;
        } else {
            msgn = 1.0;
        }
        if m2 < 30.0 {
            m2 *= DEGTORAD;
            let alpha = (1.0 - elem.ecce) / (4.0 * elem.ecce + 0.5);
            let beta = m2 / (8.0 * elem.ecce + 1.0);
            // C bug: pow(x, 1/3) is pow(x, 0) = 1.0 due to integer division.
            // Reproduce exactly: zeta is always 1.0. (swemplan.c:638)
            let zeta = (beta + (beta * beta + alpha * alpha).sqrt()).powf(0.0);
            let mut sigma = zeta - alpha / 2.0;
            sigma -= 0.078 * sigma.powi(5) / (1.0 + elem.ecce);
            e_anom = msgn * (m2 + elem.ecce * (3.0 * sigma - 4.0 * sigma.powi(3))) + m_180_or_0;
        }
    }

    e_anom = kepler(e_anom, m, elem.ecce);

    // §4.5 Position/velocity in orbital plane, then rotate to ecliptic
    let cose = e_anom.cos();
    let sine = e_anom.sin();
    let fac = ((1.0 - elem.ecce) * (1.0 + elem.ecce)).sqrt();
    let rho = 1.0 - elem.ecce * cose;

    // Orbital-plane coords (z components are zero for 2-body Kepler)
    let ox = elem.sema * (cose - elem.ecce);
    let oy = elem.sema * fac * sine;
    let ovx = -k * sine / rho;
    let ovy = k * fac * cose / rho;

    // Rotate to ecliptic (equinox of elements)
    let mut xp = [
        pqr[0] * ox + pqr[1] * oy,
        pqr[3] * ox + pqr[4] * oy,
        pqr[6] * ox + pqr[7] * oy,
        pqr[0] * ovx + pqr[1] * ovy,
        pqr[3] * ovx + pqr[4] * ovy,
        pqr[6] * ovx + pqr[7] * ovy,
    ];

    // §4.6 Ecliptic → equatorial via obliquity at tequ
    let oe = obliquity(elem.tequ, CalcFlags::empty(), models);
    let pos = rotate_x([xp[0], xp[1], xp[2]], -oe.eps);
    let vel = rotate_x([xp[3], xp[4], xp[5]], -oe.eps);
    xp[0..3].copy_from_slice(&pos);
    xp[3..6].copy_from_slice(&vel);

    // Precess to J2000 if equinox != J2000
    if elem.tequ != J2000 {
        let mut pos3 = [xp[0], xp[1], xp[2]];
        let mut vel3 = [xp[3], xp[4], xp[5]];
        precess(
            &mut pos3,
            elem.tequ,
            CalcFlags::empty(),
            models,
            PrecessionDirection::DateToJ2000,
        );
        precess(
            &mut vel3,
            elem.tequ,
            CalcFlags::empty(),
            models,
            PrecessionDirection::DateToJ2000,
        );
        xp[0..3].copy_from_slice(&pos3);
        xp[3..6].copy_from_slice(&vel3);
    }

    // Barycentric shift: geocentric elements → add Earth, heliocentric → add Sun
    let anchor = if elem.is_geo { xearth } else { xsun };
    for i in 0..6 {
        xp[i] += anchor[i];
    }

    Ok(xp)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_t_terms_constant() {
        let expr = parse_t_terms("163.7409").unwrap();
        assert!(!expr.has_t_terms);
        assert_eq!(expr.constant, 163.7409);
        assert_eq!(expr.coeffs, [0.0; 4]);
    }

    #[test]
    fn test_parse_t_terms_linear() {
        let expr = parse_t_terms("252.8987988 + 707550.7341 * T").unwrap();
        assert!(expr.has_t_terms);
        assert_eq!(expr.constant, 252.8987988);
        assert!((expr.coeffs[0] - 707550.7341).abs() < 1e-10);
        assert_eq!(expr.coeffs[1], 0.0);
    }

    #[test]
    fn test_parse_t_terms_subtraction() {
        let expr = parse_t_terms("47.787931-1670.056*T").unwrap();
        assert!(expr.has_t_terms);
        assert!((expr.constant - 47.787931).abs() < 1e-10);
        assert!((expr.coeffs[0] - (-1670.056)).abs() < 1e-10);
    }

    #[test]
    fn test_parse_t_terms_waldemath_complex() {
        let expr = parse_t_terms("70.3407215 + 109023.2634989 * T").unwrap();
        assert!(expr.has_t_terms);
        assert!((expr.constant - 70.3407215).abs() < 1e-10);
        assert!((expr.coeffs[0] - 109023.2634989).abs() < 1e-10);
    }

    #[test]
    fn test_parse_t_terms_eval() {
        let expr = parse_t_terms("252.8987988 + 707550.7341 * T").unwrap();
        // At t_days = 36525 (T=1), result should be 252.8987988 + 707550.7341
        let result = expr.eval(36525.0);
        assert!((result - (252.8987988 + 707550.7341)).abs() < 1e-8);
    }

    #[test]
    fn test_builtin_catalog() {
        let catalog = FictitiousCatalog::builtin();
        assert_eq!(catalog.rows.len(), 15);
        assert!(!catalog.from_file);

        // Verify Cupido
        let elem = resolve_elements(&catalog, 0, J2000).unwrap();
        assert_eq!(elem.name, "Cupido");
        assert!(!elem.is_geo);
        assert_eq!(elem.tjd0, J1900);
    }

    #[test]
    fn test_bodies_55_58_need_file() {
        let catalog = FictitiousCatalog::builtin();
        // Bodies 55-58 (indices 15-18) should fail without file
        assert!(resolve_elements(&catalog, 15, J2000).is_err());
        assert!(resolve_elements(&catalog, 16, J2000).is_err());
        assert!(resolve_elements(&catalog, 17, J2000).is_err());
        assert!(resolve_elements(&catalog, 18, J2000).is_err());
    }

    #[test]
    fn test_kepler_zero_eccentricity() {
        // For e=0, E=M immediately
        let m = 1.5;
        let e = kepler(m, m, 0.0);
        assert!((e - m).abs() < 1e-12);
    }

    #[test]
    fn test_kepler_low_eccentricity() {
        let m = 1.0;
        let ecce = 0.1;
        let e = kepler(m, m, ecce);
        // Verify Kepler equation: E - e*sin(E) = M
        let residual = (e - ecce * e.sin() - m).abs();
        assert!(residual < 1e-12);
    }

    #[test]
    fn test_kepler_high_eccentricity() {
        let m = 0.5;
        let ecce = 0.9;
        let e = kepler(m, m, ecce);
        let residual = (e - ecce * e.sin() - m).abs();
        assert!(residual < 1e-12);
    }

    #[test]
    fn test_parse_seorbel_file() {
        let ephe_path = std::path::Path::new("../swisseph/ephe");
        if !ephe_path.join(FICTFILE).exists() {
            return;
        }
        let contents = std::fs::read_to_string(ephe_path.join(FICTFILE)).unwrap();
        let catalog = parse_orbel_file(&contents).unwrap();
        // Should have at least 19 rows (the standard set)
        assert!(catalog.rows.len() >= 19);
        assert!(catalog.from_file);

        // Verify Cupido is row 0
        assert_eq!(catalog.rows[0].name, "Cupido");

        // Verify Kronos sema from FILE is 64.81690 (not 64.81960)
        let kronos_sema = catalog.rows[3].sema.eval(0.0);
        assert!((kronos_sema - 64.81690).abs() < 1e-10);

        // Verify Vulcan (row 15) has T-terms in mano
        assert!(catalog.rows[15].mano.has_t_terms);

        // Verify White Moon (row 16) is geocentric
        assert!(catalog.rows[16].is_geo);

        // Verify Waldemath (row 18) is geocentric
        assert!(catalog.rows[18].is_geo);
    }
}
