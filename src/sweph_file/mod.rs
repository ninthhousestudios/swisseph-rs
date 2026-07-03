mod evaluate;
mod parse;
mod segment;
pub mod types;

pub use evaluate::evaluate_body;

use std::path::Path;

use memmap2::Mmap;

use crate::error::Error;
use crate::types::Body;

pub use types::{
    AsteroidMeta, ByteOrder, FileHeader, FileType, PlanetFileData, SEI_FLG_HELIO, SEI_SUNBARY,
};

pub struct SwissEphFile {
    mmap: Mmap,
    header: FileHeader,
    planets: Vec<PlanetFileData>,
}

impl SwissEphFile {
    pub fn open(path: &Path) -> Result<Self, Error> {
        let file_type = detect_file_type(path)?;
        let file =
            std::fs::File::open(path).map_err(|_| Error::FileNotFound(path.to_path_buf()))?;
        // SAFETY: the caller must ensure the file is not truncated, replaced, or
        // modified by another process while this mapping is live. Ephemeris .se1
        // files are static data installed once and never mutated at runtime.
        let mmap = unsafe { Mmap::map(&file) }
            .map_err(|e| Error::FileFormat(format!("mmap failed: {e}")))?;
        let (header, planets) = parse::parse_file(&mmap, file_type)?;
        Ok(Self {
            mmap,
            header,
            planets,
        })
    }

    pub fn header(&self) -> &FileHeader {
        &self.header
    }

    pub fn planet_data(&self, body_id: i32) -> Option<&PlanetFileData> {
        self.planets.iter().find(|p| p.body_id == body_id)
    }

    pub fn planets(&self) -> &[PlanetFileData] {
        &self.planets
    }

    pub fn bytes(&self) -> &[u8] {
        &self.mmap
    }
}

fn detect_file_type(path: &Path) -> Result<FileType, Error> {
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
    if stem.starts_with("sepl") {
        Ok(FileType::Planet)
    } else if stem.starts_with("semo") {
        Ok(FileType::Moon)
    } else if stem.starts_with("seas") {
        Ok(FileType::MainAsteroid)
    } else if stem.starts_with("sepm") {
        Ok(FileType::PlanetaryMoon)
    } else if (stem.starts_with("se") && stem.len() > 2 && stem.as_bytes()[2].is_ascii_digit())
        || (stem.starts_with('s') && stem.len() > 1 && stem.as_bytes()[1].is_ascii_digit())
    {
        Ok(FileType::Asteroid)
    } else {
        Err(Error::FileFormat(format!(
            "unrecognized SE1 file type: {}",
            path.display()
        )))
    }
}

/// Map public Body enum to the body ID used in SE1 file ipl[] arrays.
/// Returns None for bodies not stored in SE1 files (mean nodes, fictitious, etc.).
/// Note: Body::Sun and Body::Earth both map to 0 (the EMB entry in planet files).
pub fn body_file_id(body: Body) -> Option<i32> {
    match body {
        Body::Sun | Body::Earth => Some(0),
        Body::Moon => Some(1),
        Body::Mercury => Some(2),
        Body::Venus => Some(3),
        Body::Mars => Some(4),
        Body::Jupiter => Some(5),
        Body::Saturn => Some(6),
        Body::Uranus => Some(7),
        Body::Neptune => Some(8),
        Body::Pluto => Some(9),
        _ => None,
    }
}

fn asteroid_file_candidates(dir: &Path, mpc: i32) -> [std::path::PathBuf; 4] {
    let base = if mpc > 99999 {
        format!("s{mpc:06}")
    } else {
        format!("se{mpc:05}")
    };
    let subdir = format!("ast{}", mpc / 1000);
    [
        dir.join(&subdir).join(format!("{base}.se1")),
        dir.join(&subdir).join(format!("{base}s.se1")),
        dir.join(format!("{base}.se1")),
        dir.join(format!("{base}s.se1")),
    ]
}

pub fn open_asteroid_file(dir: &Path, mpc: i32) -> Result<SwissEphFile, Error> {
    let candidates = asteroid_file_candidates(dir, mpc);
    for path in &candidates {
        match SwissEphFile::open(path) {
            Ok(f) => return Ok(f),
            Err(Error::FileNotFound(_)) => continue,
            Err(e) => return Err(e),
        }
    }
    Err(Error::FileNotFound(candidates[0].clone()))
}

pub fn open_ephemeris_files(dir: &Path, prefix: &str) -> Result<Vec<SwissEphFile>, Error> {
    let mut files = Vec::new();
    let entries = std::fs::read_dir(dir).map_err(|_| Error::FileNotFound(dir.to_path_buf()))?;
    for entry in entries {
        let entry = entry.map_err(|_| Error::FileNotFound(dir.to_path_buf()))?;
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if name_str.starts_with(prefix) && name_str.ends_with(".se1") {
            files.push(SwissEphFile::open(&entry.path())?);
        }
    }
    files.sort_by(|a, b| {
        a.header()
            .time_range
            .0
            .partial_cmp(&b.header().time_range.0)
            .unwrap()
    });
    Ok(files)
}

/// Select the ephemeris file for `jd`. `files` must be sorted ascending by
/// file-level `time_range.0` (as `open_ephemeris_files` guarantees).
///
/// Picks the latest-starting file whose file-level tfstart is ≤ `jd` and whose
/// per-planet range covers `jd`. This mirrors C's `swi_gen_filename` logic: the
/// file named for a given epoch is the one whose century boundary is the largest
/// that does not exceed the epoch. Using `<=` (not strict `<`) matches C's
/// behavior at exact file boundaries like jd=2378496.5 (1800-Jan-1 = sepl_18's
/// tfstart): C opens sepl_18 for the main position at jd, then switches to sepl_12
/// for the retarded epoch (jd-dt, which falls before sepl_18's tfstart). Callers
/// that need the retarded-time file should call this function separately with the
/// retarded jd.
pub fn find_file_for_jd(files: &[SwissEphFile], body_id: i32, jd: f64) -> Option<&SwissEphFile> {
    files.iter().rev().find(|f| {
        let (file_start, _) = f.header().time_range;
        file_start <= jd
            && f.planet_data(body_id)
                .is_some_and(|pd| jd >= pd.tfstart && jd <= pd.tfend)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn _assert_send_sync() {
        fn assert_send<T: Send>() {}
        fn assert_sync<T: Sync>() {}
        assert_send::<SwissEphFile>();
        assert_sync::<SwissEphFile>();
    }

    fn ephe_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ephe")
    }

    #[test]
    fn asteroid_eros_se1() {
        let path = ephe_dir().join("ast0/se00433s.se1");
        if !path.exists() {
            return;
        }
        let f = SwissEphFile::open(&path).unwrap();
        assert_eq!(f.header().file_type, FileType::Asteroid);
        assert_eq!(f.planets()[0].body_id, 10433);
        assert_eq!(f.header().time_range, (2268922.5, 2488522.5));
        let meta = f.header().asteroid.as_ref().unwrap();
        assert_eq!(meta.h, 10.38);
        assert_eq!(meta.g, 0.15);
        assert_eq!(meta.name, "Eros");
    }

    #[test]
    fn asteroid_eris_6digit_se1() {
        let path = ephe_dir().join("ast136/s136199s.se1");
        if !path.exists() {
            return;
        }
        let f = SwissEphFile::open(&path).unwrap();
        assert_eq!(f.header().file_type, FileType::Asteroid);
        assert_eq!(f.planets()[0].body_id, 146199);
        let meta = f.header().asteroid.as_ref().unwrap();
        assert!(meta.name.contains("Eris"), "name was: {}", meta.name);
    }

    #[test]
    fn open_asteroid_file_eros() {
        let dir = ephe_dir();
        if !dir.join("ast0/se00433s.se1").exists() {
            return;
        }
        let f = open_asteroid_file(&dir, 433).unwrap();
        assert_eq!(f.planets()[0].body_id, 10433);
    }

    #[test]
    fn open_asteroid_file_missing() {
        let dir = ephe_dir();
        let result = open_asteroid_file(&dir, 99999);
        assert!(result.is_err());
    }

    #[test]
    fn detect_file_type_6digit() {
        let path = Path::new("/tmp/s136108s.se1");
        assert_eq!(detect_file_type(path).unwrap(), FileType::Asteroid);
    }

    #[test]
    fn ephemeris_new_with_asteroids() {
        let dir = ephe_dir();
        if !dir.join("ast0/se00433s.se1").exists() {
            return;
        }
        let config = crate::config::EphemerisConfig {
            ephe_path: Some(dir),
            asteroid_numbers: vec![433],
            ..Default::default()
        };
        let eph = crate::context::Ephemeris::new(config);
        assert!(eph.is_ok());
    }

    #[test]
    fn ephemeris_new_missing_asteroid_errors() {
        let dir = ephe_dir();
        let config = crate::config::EphemerisConfig {
            ephe_path: Some(dir),
            asteroid_numbers: vec![99999],
            ..Default::default()
        };
        let eph = crate::context::Ephemeris::new(config);
        assert!(eph.is_err());
    }

    #[test]
    fn ephemeris_new_asteroid_numbers_without_path_errors() {
        let config = crate::config::EphemerisConfig {
            asteroid_numbers: vec![433],
            ..Default::default()
        };
        let eph = crate::context::Ephemeris::new(config);
        assert!(eph.is_err());
    }
}
