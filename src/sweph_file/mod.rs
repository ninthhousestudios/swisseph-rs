mod parse;
pub mod types;

use std::path::Path;

use memmap2::Mmap;

use crate::error::Error;
use crate::types::Body;

pub use types::{ByteOrder, FileHeader, FileType, PlanetFileData};

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
    } else if stem.starts_with("se") && stem.len() > 2 && stem.as_bytes()[2].is_ascii_digit() {
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

#[cfg(test)]
mod tests {
    use super::SwissEphFile;

    fn _assert_send_sync() {
        fn assert_send<T: Send>() {}
        fn assert_sync<T: Sync>() {}
        assert_send::<SwissEphFile>();
        assert_sync::<SwissEphFile>();
    }
}
