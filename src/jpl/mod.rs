//! JPL DE ephemeris (`.bin`) file reader and Chebyshev evaluation.
//!
//! Low-level internals; exposed for golden tests and advanced use.

mod header;
mod interp;

pub use header::{ByteOrder, JplHeader};

use std::path::Path;

use memmap2::Mmap;

use crate::error::Error;

// JPL body indices (swejpl.h:68–83). Used as `ntarg`/`ncent` in `jpl_pleph`
// and as slot indices in the internal `pv[13]` array.
/// JPL body index for Mercury.
pub const J_MERCURY: i32 = 0;
/// JPL body index for Venus.
pub const J_VENUS: i32 = 1;
/// JPL body index for Earth.
pub const J_EARTH: i32 = 2;
/// JPL body index for Mars.
pub const J_MARS: i32 = 3;
/// JPL body index for Jupiter.
pub const J_JUPITER: i32 = 4;
/// JPL body index for Saturn.
pub const J_SATURN: i32 = 5;
/// JPL body index for Uranus.
pub const J_URANUS: i32 = 6;
/// JPL body index for Neptune.
pub const J_NEPTUNE: i32 = 7;
/// JPL body index for Pluto.
pub const J_PLUTO: i32 = 8;
/// JPL body index for the Moon (geocentric).
pub const J_MOON: i32 = 9;
/// JPL body index for the Sun.
pub const J_SUN: i32 = 10;
/// JPL body index for the Solar System Barycenter.
pub const J_SBARY: i32 = 11;
/// JPL body index for the Earth-Moon Barycenter.
pub const J_EMB: i32 = 12;
/// JPL body index for nutations.
pub const J_NUT: i32 = 13;
/// JPL body index for lunar librations.
pub const J_LIB: i32 = 14;

/// A memory-mapped, parsed JPL DE ephemeris file.
pub struct JplFile {
    mmap: Mmap,
    header: JplHeader,
}

impl JplFile {
    /// Open and parse the JPL DE file at `path`, memory-mapping its contents.
    pub fn open(path: &Path) -> Result<Self, Error> {
        let file =
            std::fs::File::open(path).map_err(|_| Error::FileNotFound(path.to_path_buf()))?;
        // SAFETY: the caller must ensure the file is not truncated, replaced, or
        // modified by another process while this mapping is live. JPL DE files
        // are static data installed once and never mutated at runtime.
        let mmap = unsafe { Mmap::map(&file) }
            .map_err(|e| Error::FileFormat(format!("mmap failed: {e}")))?;
        let header = header::parse_header(&mmap)?;
        header::validate_file_length(&mmap, &header)?;
        Ok(Self { mmap, header })
    }

    /// Return the parsed file header.
    pub fn header(&self) -> &JplHeader {
        &self.header
    }

    /// Return the raw memory-mapped file bytes.
    pub fn bytes(&self) -> &[u8] {
        &self.mmap
    }

    /// Return the file's byte order.
    pub fn byte_order(&self) -> ByteOrder {
        self.header.byte_order
    }
}

/// Return the state of body `ntarg` relative to `ncent` in barycentric equatorial
/// J2000/ICRF. Units: AU (position), AU/day (velocity). (swejpl.c:362–449)
///
/// `ntarg` / `ncent`: J_* constants defined in this module.
/// `need_speed`: when false, velocity components of the result are zero.
pub fn jpl_pleph(
    file: &JplFile,
    et: f64,
    ntarg: i32,
    ncent: i32,
    need_speed: bool,
) -> Result<[f64; 6], Error> {
    let val: u8 = if need_speed { 2 } else { 1 };
    let mut list = [0u8; 12];

    // Populate list[] for target and center. Sun and SBARY need no list entry
    // (Sun comes from pvsun always; SBARY is the zero origin). (swejpl.c:374–416)
    for &body in &[ntarg, ncent] {
        match body {
            b if b == J_MOON => {
                list[J_MOON as usize] = val;
                list[J_EARTH as usize] = val;
            }
            b if b == J_EARTH => {
                list[J_EARTH as usize] = val;
                list[J_MOON as usize] = val;
            }
            b if b == J_EMB => {
                list[J_EARTH as usize] = val;
            }
            b if (0..10).contains(&b) => {
                list[b as usize] = val;
            }
            _ => {} // J_SUN, J_SBARY: no list entry needed
        }
    }

    let (mut pv, pvsun) = interp::state(file, et, &list, true, need_speed)?;

    // Post-state assembly (swejpl.c:418–447).
    // Order matters: copy EMB slot before Earth/Moon decomposition alters pv[J_EARTH].
    if ntarg == J_SUN || ncent == J_SUN {
        pv[J_SUN as usize] = pvsun;
    }
    if ntarg == J_SBARY || ncent == J_SBARY {
        pv[J_SBARY as usize] = [0.0; 6];
    }
    if ntarg == J_EMB || ncent == J_EMB {
        pv[J_EMB as usize] = pv[J_EARTH as usize];
    }

    // Earth/Moon decomposition: pv[J_EARTH] from state() is the EMB (barycentric
    // Earth-Moon Barycenter); pv[J_MOON] is geocentric Moon.
    let is_earth_moon_pair =
        (ntarg == J_EARTH && ncent == J_MOON) || (ntarg == J_MOON && ncent == J_EARTH);

    if is_earth_moon_pair {
        // Moon is already geocentric; zero Earth so the result is the raw geocentric Moon.
        pv[J_EARTH as usize] = [0.0; 6];
    } else {
        let emrat = file.header().emrat;
        if list[J_EARTH as usize] > 0 {
            // EMB → barycentric Earth: Earth = EMB - Moon_geo / (emrat + 1)
            let moon = pv[J_MOON as usize];
            for k in 0..6 {
                pv[J_EARTH as usize][k] -= moon[k] / (emrat + 1.0);
            }
        }
        if list[J_MOON as usize] > 0 {
            // Geocentric Moon → barycentric Moon: Moon_bary = Moon_geo + Earth_bary
            let earth = pv[J_EARTH as usize];
            for k in 0..6 {
                pv[J_MOON as usize][k] += earth[k];
            }
        }
    }

    let mut rrd = [0.0f64; 6];
    for k in 0..6 {
        rrd[k] = pv[ntarg as usize][k] - pv[ncent as usize][k];
    }

    Ok(rrd)
}

#[cfg(test)]
mod tests {
    use super::JplFile;

    fn _assert_send_sync() {
        fn assert_send<T: Send>() {}
        fn assert_sync<T: Sync>() {}
        assert_send::<JplFile>();
        assert_sync::<JplFile>();
    }
}
