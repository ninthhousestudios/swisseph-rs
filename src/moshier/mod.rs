pub mod moon;
pub mod moon_tables;
pub mod planets;
pub mod tables;

pub use tables::PLANETS;

pub struct PlantTbl {
    pub max_harmonic: [i8; 9],
    pub max_power_of_t: i8,
    pub arg_tbl: &'static [i8],
    pub lon_tbl: &'static [f64],
    pub lat_tbl: &'static [f64],
    pub rad_tbl: &'static [f64],
    pub distance: f64,
}

#[cfg(test)]
mod tests {
    use super::PlantTbl;
    use super::tables::*;

    #[test]
    fn element_counts_match_c() {
        // Total terms from c-ref-moshier.md line 114–124
        // arg_tbl lengths verified against swemptab.h array sizes
        let expected: &[(&str, &PlantTbl, usize)] = &[
            ("mer404", &MER404, 130),
            ("ven404", &VEN404, 108),
            ("ear404", &EAR404, 135),
            ("mar404", &MAR404, 201),
            ("jup404", &JUP404, 142),
            ("sat404", &SAT404, 215),
            ("ura404", &URA404, 177),
            ("nep404", &NEP404, 59),
            ("plu404", &PLU404, 173),
        ];

        for (name, tbl, total_terms) in expected {
            // Walk arg_tbl to count terms and verify coefficient consumption
            let mut p = 0;
            let mut term_count = 0usize;
            let mut lon_consumed = 0usize;
            let mut lat_consumed = 0usize;
            let mut rad_consumed = 0usize;

            while p < tbl.arg_tbl.len() {
                let np = tbl.arg_tbl[p];
                p += 1;
                if np < 0 {
                    break;
                }
                term_count += 1;
                if np == 0 {
                    // Polynomial term
                    let nt = tbl.arg_tbl[p] as usize;
                    p += 1;
                    lon_consumed += nt + 1;
                    lat_consumed += nt + 1;
                    rad_consumed += nt + 1;
                } else {
                    // Periodic term: np argument pairs, then nt
                    p += (np as usize) * 2;
                    let nt = tbl.arg_tbl[p] as usize;
                    p += 1;
                    lon_consumed += 2 * (nt + 1);
                    lat_consumed += 2 * (nt + 1);
                    rad_consumed += 2 * (nt + 1);
                }
            }

            assert_eq!(
                term_count, *total_terms,
                "{name}: expected {total_terms} terms, got {term_count}"
            );
            assert_eq!(
                lon_consumed,
                tbl.lon_tbl.len(),
                "{name}: lon_tbl length mismatch (consumed {lon_consumed}, actual {})",
                tbl.lon_tbl.len()
            );
            assert_eq!(
                lat_consumed,
                tbl.lat_tbl.len(),
                "{name}: lat_tbl length mismatch (consumed {lat_consumed}, actual {})",
                tbl.lat_tbl.len()
            );
            assert_eq!(
                rad_consumed,
                tbl.rad_tbl.len(),
                "{name}: rad_tbl length mismatch (consumed {rad_consumed}, actual {})",
                tbl.rad_tbl.len()
            );
        }
    }

    #[test]
    fn planets_array_order() {
        // C planets[] order: mer, ven, ear, mar, jup, sat, ura, nep, plu
        let expected_distances = [
            MER404.distance,
            VEN404.distance,
            EAR404.distance,
            MAR404.distance,
            JUP404.distance,
            SAT404.distance,
            URA404.distance,
            NEP404.distance,
            PLU404.distance,
        ];
        for (i, dist) in expected_distances.iter().enumerate() {
            assert_eq!(PLANETS[i].distance, *dist, "PLANETS[{i}] distance mismatch");
        }
    }
}
