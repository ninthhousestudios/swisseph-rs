mod header;

pub use header::{ByteOrder, JplHeader};

use std::path::Path;

use memmap2::Mmap;

use crate::error::Error;

pub struct JplFile {
    mmap: Mmap,
    header: JplHeader,
}

impl JplFile {
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

    pub fn header(&self) -> &JplHeader {
        &self.header
    }

    pub fn bytes(&self) -> &[u8] {
        &self.mmap
    }

    pub fn byte_order(&self) -> ByteOrder {
        self.header.byte_order
    }
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
