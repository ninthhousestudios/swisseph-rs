use crate::error::Error;

use super::types::{
    AsteroidMeta, ByteOrder, ENDIAN_TEST, FileHeader, FileType, PlanetFileData, SE_AST_OFFSET,
    SE_PLMOON_OFFSET, SEI_FLG_ELLIPSE,
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

    fn read_bytes(&mut self, n: usize) -> Result<&'a [u8], Error> {
        self.ensure(n)?;
        let slice = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(slice)
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

fn atof_prefix(bytes: &[u8]) -> f64 {
    let mut i = 0;
    while i < bytes.len() && bytes[i] == b' ' {
        i += 1;
    }
    if i >= bytes.len() {
        return 0.0;
    }
    let start = i;
    if bytes[i] == b'+' || bytes[i] == b'-' {
        i += 1;
    }
    let has_digits_before = i < bytes.len() && bytes[i].is_ascii_digit();
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i < bytes.len() && bytes[i] == b'.' {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
    }
    if !has_digits_before && (i == start || i == start + 1) {
        return 0.0;
    }
    if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
        let j = i + 1;
        let mut k = j;
        if k < bytes.len() && (bytes[k] == b'+' || bytes[k] == b'-') {
            k += 1;
        }
        if k < bytes.len() && bytes[k].is_ascii_digit() {
            i = k;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
        }
    }
    let s = std::str::from_utf8(&bytes[start..i]).unwrap_or("");
    s.parse::<f64>().unwrap_or(0.0)
}

pub(super) fn parse_mpc_elements(line: &[u8]) -> (f64, f64, f64) {
    let mut sp = 0;
    while sp < line.len() && line[sp] == b' ' {
        sp += 1;
    }
    while sp < line.len() && line[sp].is_ascii_digit() {
        sp += 1;
    }
    if sp < line.len() {
        sp += 1;
    }
    let i = sp;

    let h = if 35 + i < line.len() {
        atof_prefix(&line[35 + i..])
    } else {
        0.0
    };
    let mut g = if 42 + i < line.len() {
        atof_prefix(&line[42 + i..])
    } else {
        0.0
    };
    if g == 0.0 {
        g = 0.15;
    }

    let diam_raw = if 51 + i + 7 <= line.len() {
        atof_prefix(&line[51 + i..51 + i + 7])
    } else if 51 + i < line.len() {
        atof_prefix(&line[51 + i..])
    } else {
        0.0
    };
    let diameter_km = if diam_raw == 0.0 {
        1329.0 / 0.15_f64.sqrt() * 10f64.powf(-0.2 * h)
    } else {
        diam_raw
    };

    (h, g, diameter_km)
}

fn extract_asteroid_name(line: &[u8], ipl0: i32, name_field_30: &[u8]) -> String {
    let mut sp = 0;
    while sp < line.len() && line[sp] == b' ' {
        sp += 1;
    }
    while sp < line.len() && line[sp].is_ascii_digit() {
        sp += 1;
    }
    if sp < line.len() {
        sp += 1;
    }
    let i = sp;
    let lastnam = 19;

    let sastnam_len = (lastnam + i).min(line.len());
    let sastnam = &line[..sastnam_len];

    let mut j = 4usize;
    while j < 10 && j < sastnam.len() && sastnam[j] != b' ' {
        j += 1;
    }

    let sastno = &sastnam[..j.min(sastnam.len())];
    let num_str = std::str::from_utf8(sastno).unwrap_or("").trim();
    let parsed_num = num_str.parse::<i32>().unwrap_or(0);

    let name_bytes = if parsed_num == ipl0 - SE_AST_OFFSET || parsed_num == ipl0 {
        let start = (j + 1).min(sastnam.len());
        let end = (start + lastnam).min(sastnam.len());
        &sastnam[start..end]
    } else {
        name_field_30
    };

    let mut name = String::from_utf8_lossy(name_bytes).into_owned();
    name = name.trim_end().to_string();
    if let Some(pos) = name.find("  ") {
        name.truncate(pos);
    }
    name
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
        FileType::Asteroid | FileType::PlanetaryMoon => 4,
        _ => 3,
    };

    let first_crlf = find_crlf(data, 0)?;
    let version = parse_version(&data[..first_crlf])?;

    let mut pos = 0;
    let mut mpc_line: Option<&[u8]> = None;
    for line_idx in 0..num_text_lines {
        let crlf = find_crlf(data, pos)?;
        if line_idx == 3 && file_type == FileType::Asteroid {
            mpc_line = Some(&data[pos..crlf]);
        }
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

    // Asteroid name + 30-byte field handling (c-ref-asteroid.md §3.3)
    let asteroid = if let Some(line) = mpc_line {
        let (h, g, diameter_km) = parse_mpc_elements(line);
        let ipl0 = ipl.first().copied().unwrap_or(0);
        let name_field_30 = r.read_bytes(30)?;
        let name = extract_asteroid_name(line, ipl0, name_field_30);
        Some(AsteroidMeta {
            h,
            g,
            diameter_km,
            name,
        })
    } else if file_type == FileType::PlanetaryMoon {
        r.skip(30)?;
        None
    } else {
        None
    };

    // CRC32 field (skip — not validating)
    r.skip(4)?;

    // Physical constants: 5 doubles (read and discard)
    r.skip(40)?;

    // --- Per-planet metadata ---
    let mut planets = Vec::with_capacity(nplan);
    for &ipli in &ipl {
        let lndx0 = r.read_i32()? as usize;
        let iflg = r.read_u8()? as u32;
        let ncoe = r.read_u8()? as usize;

        let rmax_raw = r.read_i32()?;
        let rmax = if (SE_PLMOON_OFFSET..SE_AST_OFFSET).contains(&ipli)
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
        asteroid,
    };

    Ok((header, planets))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atof_prefix() {
        assert_eq!(atof_prefix(b"10.38"), 10.38);
        assert_eq!(atof_prefix(b"  0.15"), 0.15);
        assert_eq!(atof_prefix(b""), 0.0);
        assert_eq!(atof_prefix(b"   "), 0.0);
        assert_eq!(atof_prefix(b"-3.5e2"), -350.0);
        assert_eq!(atof_prefix(b"42xyz"), 42.0);
    }

    #[test]
    fn test_parse_mpc_elements_eros() {
        let line =
            b"000433 Eros               L.H. Wasserman  10.38  0.15                    4   0 ";
        let (h, g, diam) = parse_mpc_elements(line);
        assert_eq!(h, 10.38);
        assert_eq!(g, 0.15);
        let expected_diam = 1329.0 / 0.15_f64.sqrt() * 10f64.powf(-0.2 * 10.38);
        assert_eq!(diam, expected_diam);
    }
}
