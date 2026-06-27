use crate::error::Error;

use super::types::{ByteOrder, PlanetFileData};

pub(super) struct SegmentData {
    pub coeffs: Vec<f64>,
    pub tseg0: f64,
    #[allow(dead_code)]
    pub tseg1: f64,
}

pub(super) fn unpack_segment(
    data: &[u8],
    planet: &PlanetFileData,
    jd: f64,
    order: ByteOrder,
) -> Result<SegmentData, Error> {
    let iseg = ((jd - planet.tfstart) / planet.dseg) as usize;
    let tseg0 = planet.tfstart + iseg as f64 * planet.dseg;
    let tseg1 = tseg0 + planet.dseg;

    let fpos = read_segment_offset(data, planet.lndx0 + iseg * 3, order)?;
    let ncoe = planet.ncoe;
    let rmax = planet.rmax;
    let mut coeffs = vec![0.0; ncoe * 3];
    let mut pos = fpos;

    for icoord in 0..3usize {
        let idbl_base = icoord * ncoe;
        let mut idbl = idbl_base;

        let (nsizes, nsize, header_len) = parse_coord_header(data, pos)?;
        pos += header_len;

        let nco: usize = nsize[..nsizes].iter().sum();
        if nco > ncoe {
            return Err(Error::FileFormat(format!(
                "nco {nco} exceeds ncoe {ncoe} for coord {icoord}"
            )));
        }

        for i in 0..nsizes {
            let count = nsize[i];
            if count == 0 {
                continue;
            }
            match i {
                0..=3 => {
                    let width = 4 - i;
                    pos = unpack_integer_coeffs(
                        data,
                        pos,
                        order,
                        width,
                        count,
                        rmax,
                        &mut coeffs,
                        &mut idbl,
                    )?;
                }
                4 => {
                    pos = unpack_nibble_coeffs(
                        data,
                        pos,
                        order,
                        count,
                        rmax,
                        &mut coeffs,
                        &mut idbl,
                    )?;
                }
                5 => {
                    pos = unpack_quarter_coeffs(
                        data,
                        pos,
                        order,
                        count,
                        rmax,
                        &mut coeffs,
                        &mut idbl,
                    )?;
                }
                _ => unreachable!(),
            }
        }
    }

    Ok(SegmentData {
        coeffs,
        tseg0,
        tseg1,
    })
}

fn read_segment_offset(data: &[u8], pos: usize, order: ByteOrder) -> Result<usize, Error> {
    if pos + 3 > data.len() {
        return Err(Error::FileFormat("segment index out of bounds".into()));
    }
    let mut buf = [0u8; 4];
    match order {
        ByteOrder::Big => {
            buf[1] = data[pos];
            buf[2] = data[pos + 1];
            buf[3] = data[pos + 2];
            Ok(i32::from_be_bytes(buf) as usize)
        }
        ByteOrder::Little => {
            buf[0] = data[pos];
            buf[1] = data[pos + 1];
            buf[2] = data[pos + 2];
            Ok(i32::from_le_bytes(buf) as usize)
        }
    }
}

fn parse_coord_header(data: &[u8], pos: usize) -> Result<(usize, [usize; 6], usize), Error> {
    if pos + 2 > data.len() {
        return Err(Error::FileFormat("coord header out of bounds".into()));
    }
    let c0 = data[pos];
    let c1 = data[pos + 1];

    if c0 & 128 != 0 {
        if pos + 4 > data.len() {
            return Err(Error::FileFormat(
                "extended coord header out of bounds".into(),
            ));
        }
        let c2 = data[pos + 2];
        let c3 = data[pos + 3];
        let nsize = [
            (c1 >> 4) as usize,
            (c1 & 0x0f) as usize,
            (c2 >> 4) as usize,
            (c2 & 0x0f) as usize,
            (c3 >> 4) as usize,
            (c3 & 0x0f) as usize,
        ];
        Ok((6, nsize, 4))
    } else {
        let nsize = [
            (c0 >> 4) as usize,
            (c0 & 0x0f) as usize,
            (c1 >> 4) as usize,
            (c1 & 0x0f) as usize,
            0,
            0,
        ];
        Ok((4, nsize, 2))
    }
}

fn ensure(data: &[u8], pos: usize, n: usize) -> Result<(), Error> {
    if pos + n > data.len() {
        Err(Error::FileFormat("unexpected end of segment data".into()))
    } else {
        Ok(())
    }
}

fn read_integer(data: &[u8], pos: usize, width: usize, order: ByteOrder) -> u32 {
    match (width, order) {
        (4, ByteOrder::Little) => {
            u32::from_le_bytes([data[pos], data[pos + 1], data[pos + 2], data[pos + 3]])
        }
        (4, ByteOrder::Big) => {
            u32::from_be_bytes([data[pos], data[pos + 1], data[pos + 2], data[pos + 3]])
        }
        (3, ByteOrder::Little) => u32::from_le_bytes([data[pos], data[pos + 1], data[pos + 2], 0]),
        (3, ByteOrder::Big) => u32::from_be_bytes([0, data[pos], data[pos + 1], data[pos + 2]]),
        (2, ByteOrder::Little) => u16::from_le_bytes([data[pos], data[pos + 1]]) as u32,
        (2, ByteOrder::Big) => u16::from_be_bytes([data[pos], data[pos + 1]]) as u32,
        (1, _) => data[pos] as u32,
        _ => unreachable!(),
    }
}

fn unpack_integer_coeffs(
    data: &[u8],
    mut pos: usize,
    order: ByteOrder,
    width: usize,
    count: usize,
    rmax: f64,
    coeffs: &mut [f64],
    idbl: &mut usize,
) -> Result<usize, Error> {
    ensure(data, pos, width * count)?;
    for _ in 0..count {
        let raw = read_integer(data, pos, width, order);
        // C order: (int_val) / 1e9 * rmax / 2
        coeffs[*idbl] = if raw & 1 != 0 {
            -(((raw.wrapping_add(1)) / 2) as f64 / 1e9 * rmax / 2.0)
        } else {
            (raw / 2) as f64 / 1e9 * rmax / 2.0
        };
        pos += width;
        *idbl += 1;
    }
    Ok(pos)
}

fn unpack_nibble_coeffs(
    data: &[u8],
    mut pos: usize,
    order: ByteOrder,
    count: usize,
    rmax: f64,
    coeffs: &mut [f64],
    idbl: &mut usize,
) -> Result<usize, Error> {
    let nbytes = (count + 1) / 2;
    ensure(data, pos, nbytes)?;
    let mut j = 0;
    for _ in 0..nbytes {
        let byte = read_integer(data, pos, 1, order);
        pos += 1;
        let mut lval = byte;
        let mut o: u32 = 16;
        for _ in 0..2 {
            if j >= count {
                break;
            }
            // C order: int_val * rmax / 2 / 1e9
            if lval & o != 0 {
                coeffs[*idbl] = -(((lval.wrapping_add(o)) / o / 2) as f64 * rmax / 2.0 / 1e9);
            } else {
                coeffs[*idbl] = (lval / o / 2) as f64 * rmax / 2.0 / 1e9;
            }
            lval %= o;
            o /= 16;
            j += 1;
            *idbl += 1;
        }
    }
    Ok(pos)
}

fn unpack_quarter_coeffs(
    data: &[u8],
    mut pos: usize,
    order: ByteOrder,
    count: usize,
    rmax: f64,
    coeffs: &mut [f64],
    idbl: &mut usize,
) -> Result<usize, Error> {
    let nbytes = (count + 3) / 4;
    ensure(data, pos, nbytes)?;
    let mut j = 0;
    for _ in 0..nbytes {
        let byte = read_integer(data, pos, 1, order);
        pos += 1;
        let mut lval = byte;
        let mut o: u32 = 64;
        for _ in 0..4 {
            if j >= count {
                break;
            }
            // C order: int_val * rmax / 2 / 1e9
            if lval & o != 0 {
                coeffs[*idbl] = -(((lval.wrapping_add(o)) / o / 2) as f64 * rmax / 2.0 / 1e9);
            } else {
                coeffs[*idbl] = (lval / o / 2) as f64 * rmax / 2.0 / 1e9;
            }
            lval %= o;
            o /= 4;
            j += 1;
            *idbl += 1;
        }
    }
    Ok(pos)
}
