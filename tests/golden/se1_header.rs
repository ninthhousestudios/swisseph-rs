use serde::Deserialize;
use std::path::PathBuf;
use swisseph::sweph_file::SwissEphFile;

#[derive(Deserialize)]
struct PlanetCase {
    body_id: i32,
    iflg: u32,
    ncoe: usize,
    neval: usize,
    rmax: f64,
    dseg: f64,
    tfstart: f64,
    tfend: f64,
    lndx0: usize,
    nndx: usize,
    telem: f64,
    prot: f64,
    qrot: f64,
    dprot: f64,
    dqrot: f64,
    peri: f64,
    dperi: f64,
    #[serde(default)]
    refep: Option<Vec<f64>>,
}

#[derive(Deserialize)]
struct FileCase {
    filename: String,
    version: i32,
    file_type: String,
    tfstart: f64,
    tfend: f64,
    denum: i32,
    #[allow(dead_code)]
    byte_order: String,
    planets: Vec<PlanetCase>,
}

#[derive(Deserialize)]
struct GoldenData {
    files: Vec<FileCase>,
}

fn ephe_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe")
        .join(name)
}

#[test]
fn header_fields_match_c() {
    let data: GoldenData = serde_json::from_str(
        &std::fs::read_to_string(super::golden_data_path("se1_header.json")).unwrap(),
    )
    .unwrap();

    for fc in &data.files {
        let path = ephe_path(&fc.filename);
        let eph = SwissEphFile::open(&path)
            .unwrap_or_else(|e| panic!("failed to open {}: {e}", fc.filename));
        let h = eph.header();

        assert_eq!(h.version, fc.version, "{}: version", fc.filename);
        super::assert_f64_exact(
            &format!("{}: file tfstart", fc.filename),
            fc.tfstart,
            h.time_range.0,
        );
        super::assert_f64_exact(
            &format!("{}: file tfend", fc.filename),
            fc.tfend,
            h.time_range.1,
        );
        assert_eq!(h.denum, fc.denum, "{}: denum", fc.filename);

        let expected_type = match fc.file_type.as_str() {
            "planet" => swisseph::sweph_file::FileType::Planet,
            "moon" => swisseph::sweph_file::FileType::Moon,
            "main_asteroid" => swisseph::sweph_file::FileType::MainAsteroid,
            "asteroid" => swisseph::sweph_file::FileType::Asteroid,
            "planetary_moon" => swisseph::sweph_file::FileType::PlanetaryMoon,
            other => panic!("unknown file type: {other}"),
        };
        assert_eq!(h.file_type, expected_type, "{}: file_type", fc.filename);

        assert_eq!(
            eph.planets().len(),
            fc.planets.len(),
            "{}: planet count",
            fc.filename
        );

        for pc in &fc.planets {
            let pd = eph
                .planet_data(pc.body_id)
                .unwrap_or_else(|| panic!("missing body_id {} in {}", pc.body_id, fc.filename));
            let label = format!("{}:body{}", fc.filename, pc.body_id);
            assert_eq!(pd.iflg, pc.iflg, "{label}: iflg");
            assert_eq!(pd.ncoe, pc.ncoe, "{label}: ncoe");
            assert_eq!(pd.neval, pc.neval, "{label}: neval");
            super::assert_f64_exact(&format!("{label}: rmax"), pc.rmax, pd.rmax);
            super::assert_f64_exact(&format!("{label}: dseg"), pc.dseg, pd.dseg);
            super::assert_f64_exact(&format!("{label}: tfstart"), pc.tfstart, pd.tfstart);
            super::assert_f64_exact(&format!("{label}: tfend"), pc.tfend, pd.tfend);
            assert_eq!(pd.lndx0, pc.lndx0, "{label}: lndx0");
            assert_eq!(pd.nndx, pc.nndx, "{label}: nndx");
            super::assert_f64_exact(&format!("{label}: telem"), pc.telem, pd.telem);
            super::assert_f64_exact(&format!("{label}: prot"), pc.prot, pd.prot);
            super::assert_f64_exact(&format!("{label}: qrot"), pc.qrot, pd.qrot);
            super::assert_f64_exact(&format!("{label}: dprot"), pc.dprot, pd.dprot);
            super::assert_f64_exact(&format!("{label}: dqrot"), pc.dqrot, pd.dqrot);
            super::assert_f64_exact(&format!("{label}: peri"), pc.peri, pd.peri);
            super::assert_f64_exact(&format!("{label}: dperi"), pc.dperi, pd.dperi);
            match (&pc.refep, &pd.refep) {
                (Some(expected), Some(actual)) => {
                    assert_eq!(expected.len(), actual.len(), "{label}: refep length");
                    for (j, (e, a)) in expected.iter().zip(actual).enumerate() {
                        super::assert_f64_exact(&format!("{label}: refep[{j}]"), *e, *a);
                    }
                }
                (None, None) => {}
                (Some(_), None) => panic!("{label}: expected refep but got None"),
                (None, Some(_)) => panic!("{label}: expected no refep but got Some"),
            }
        }
    }
}

#[test]
fn byte_order_detection_all_files() {
    let ephe_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("swisseph")
        .join("ephe");
    let mut count = 0;
    let mut stack = vec![ephe_dir.clone()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "se1") {
                SwissEphFile::open(&path)
                    .unwrap_or_else(|e| panic!("failed to open {}: {e}", path.display()));
                count += 1;
            }
        }
    }
    assert!(count > 0, "no .se1 files found in {}", ephe_dir.display());
}

#[test]
fn body_file_id_mapping() {
    use swisseph::Body;
    use swisseph::sweph_file::body_file_id;

    assert_eq!(body_file_id(Body::Sun), Some(0));
    assert_eq!(body_file_id(Body::Moon), Some(1));
    assert_eq!(body_file_id(Body::Mercury), Some(2));
    assert_eq!(body_file_id(Body::Venus), Some(3));
    assert_eq!(body_file_id(Body::Mars), Some(4));
    assert_eq!(body_file_id(Body::Jupiter), Some(5));
    assert_eq!(body_file_id(Body::Saturn), Some(6));
    assert_eq!(body_file_id(Body::Uranus), Some(7));
    assert_eq!(body_file_id(Body::Neptune), Some(8));
    assert_eq!(body_file_id(Body::Pluto), Some(9));
    assert_eq!(body_file_id(Body::Earth), Some(0));
    assert_eq!(body_file_id(Body::MeanNode), None);
    assert_eq!(body_file_id(Body::TrueNode), None);
}
