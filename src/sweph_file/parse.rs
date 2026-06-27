use crate::error::Error;

use super::types::{
    ByteOrder, ENDIAN_TEST, FileHeader, FileType, PlanetFileData, SE_AST_OFFSET, SE_PLMOON_OFFSET,
    SEI_FLG_ELLIPSE,
};

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    order: ByteOrder,
}

impl<'a> Reader<'a> {
    fn new(data: &'a [u8], pos: usize, order: ByteOrder) -> Self {
        Self { data, pos, order }
    }

    fn ensure(&self, n: usize) -> Result<(), Error> {
        if self.pos + n > self.data.len() {
            Err(Error::FileFormat("unexpected end of file".into()))
        } else {
            Ok(())
        }
    }

    fn read_u8(&mut self) -> Result<u8, Error> {
        self.ensure(1)?;
        let v = self.data[self.pos];
        self.pos += 1;
        Ok(v)
    }

    fn read_i16(&mut self) -> Result<i16, Error> {
        self.ensure(2)?;
        let bytes: [u8; 2] = self.data[self.pos..self.pos + 2].try_into().unwrap();
        self.pos += 2;
        Ok(self.order.read_i16(bytes))
    }

    fn read_i32(&mut self) -> Result<i32, Error> {
        self.ensure(4)?;
        let bytes: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(self.order.read_i32(bytes))
    }

    fn read_f64(&mut self) -> Result<f64, Error> {
        self.ensure(8)?;
        let bytes: [u8; 8] = self.data[self.pos..self.pos + 8].try_into().unwrap();
        self.pos += 8;
        Ok(self.order.read_f64(bytes))
    }

    fn skip(&mut self, n: usize) -> Result<(), Error> {
        self.ensure(n)?;
        self.pos += n;
        Ok(())
    }
}

fn find_crlf(data: &[u8], start: usize) -> Result<usize, Error> {
    for i in start..data.len().saturating_sub(1) {
        if data[i] == b'\r' && data[i + 1] == b'\n' {
            return Ok(i);
        }
    }
    Err(Error::FileFormat("missing CRLF in text header".into()))
}

fn parse_version(line: &[u8]) -> Result<i32, Error> {
    let s = std::str::from_utf8(line)
        .map_err(|_| Error::FileFormat("invalid UTF-8 in header".into()))?;
    let digit_start = s
        .find(|c: char| c.is_ascii_digit())
        .ok_or_else(|| Error::FileFormat("no version number in header".into()))?;
    let digit_end = s[digit_start..]
        .find(|c: char| !c.is_ascii_digit())
        .map(|i| digit_start + i)
        .unwrap_or(s.len());
    s[digit_start..digit_end]
        .parse::<i32>()
        .map_err(|_| Error::FileFormat("invalid version number".into()))
}

fn detect_byte_order(data: &[u8], offset: usize) -> Result<ByteOrder, Error> {
    if offset + 4 > data.len() {
        return Err(Error::FileFormat("file too short for endian test".into()));
    }
    let bytes: [u8; 4] = data[offset..offset + 4].try_into().unwrap();
    if u32::from_le_bytes(bytes) == ENDIAN_TEST {
        Ok(ByteOrder::Little)
    } else if u32::from_be_bytes(bytes) == ENDIAN_TEST {
        Ok(ByteOrder::Big)
    } else {
        Err(Error::FileFormat("invalid endian test value".into()))
    }
}

pub(super) fn parse_file(
    data: &[u8],
    file_type: FileType,
) -> Result<(FileHeader, Vec<PlanetFileData>), Error> {
    // --- Text header ---
    // Individual asteroid files (SEI_FILE_ANY_AST) have a 4th MPC elements line.
    // All other types (planet, moon, main asteroid) have 3 text lines.
    let num_text_lines: usize = match file_type {
        FileType::Asteroid => 4,
        _ => 3,
    };

    let first_crlf = find_crlf(data, 0)?;
    let version = parse_version(&data[..first_crlf])?;

    let mut pos = 0;
    for _ in 0..num_text_lines {
        let crlf = find_crlf(data, pos)?;
        pos = crlf + 2;
    }
    let binary_start = pos;

    // --- Binary header ---
    let order = detect_byte_order(data, binary_start)?;
    let mut r = Reader::new(data, binary_start + 4, order);

    let file_len = r.read_i32()? as usize;
    if file_len != data.len() {
        return Err(Error::FileFormat(format!(
            "file length mismatch: header says {file_len}, actual {}",
            data.len()
        )));
    }

    let denum = r.read_i32()?;
    let tfstart = r.read_f64()?;
    let tfend = r.read_f64()?;

    let nplan_raw = r.read_i16()? as i32;
    let (nbytes_ipl, nplan) = if nplan_raw > 256 {
        (4usize, (nplan_raw % 256) as usize)
    } else {
        (2usize, nplan_raw as usize)
    };
    if nplan == 0 || nplan > 50 {
        return Err(Error::FileFormat(format!("invalid planet count: {nplan}")));
    }

    let mut ipl = Vec::with_capacity(nplan);
    for _ in 0..nplan {
        let id = if nbytes_ipl == 4 {
            r.read_i32()?
        } else {
            r.read_i16()? as i32
        };
        ipl.push(id);
    }

    if file_type == FileType::Asteroid {
        r.skip(30)?;
    }

    // CRC32 field (skip — not validating)
    r.skip(4)?;

    // Physical constants: 5 doubles (read and discard)
    r.skip(40)?;

    // --- Per-planet metadata ---
    let mut planets = Vec::with_capacity(nplan);
    for i in 0..nplan {
        let ipli = ipl[i];

        let lndx0 = r.read_i32()? as usize;
        let iflg = r.read_u8()? as u32;
        let ncoe = r.read_u8()? as usize;

        let rmax_raw = r.read_i32()?;
        let rmax = if ipli >= SE_PLMOON_OFFSET
            && ipli < SE_AST_OFFSET
            && ((ipli % 100) == 99 || (ipli - 9000) / 100 == 4)
        {
            rmax_raw as f64 / 1_000_000.0
        } else {
            rmax_raw as f64 / 1000.0
        };

        let orbit_tfstart = r.read_f64()?;
        let orbit_tfend = r.read_f64()?;
        let dseg = r.read_f64()?;
        let telem = r.read_f64()?;
        let prot = r.read_f64()?;
        let dprot = r.read_f64()?;
        let qrot = r.read_f64()?;
        let dqrot = r.read_f64()?;
        let peri = r.read_f64()?;
        let dperi = r.read_f64()?;

        let nndx = ((orbit_tfend - orbit_tfstart + 0.1) / dseg) as usize;

        let refep = if iflg & SEI_FLG_ELLIPSE != 0 {
            let count = ncoe * 2;
            let mut coeffs = Vec::with_capacity(count);
            for _ in 0..count {
                coeffs.push(r.read_f64()?);
            }
            Some(coeffs)
        } else {
            None
        };

        planets.push(PlanetFileData {
            body_id: ipli,
            iflg,
            ncoe,
            neval: ncoe,
            rmax,
            dseg,
            tfstart: orbit_tfstart,
            tfend: orbit_tfend,
            lndx0,
            nndx,
            telem,
            prot,
            qrot,
            dprot,
            dqrot,
            peri,
            dperi,
            refep,
        });
    }

    let header = FileHeader {
        version,
        file_type,
        time_range: (tfstart, tfend),
        denum,
        byte_order: order,
    };

    Ok((header, planets))
}
