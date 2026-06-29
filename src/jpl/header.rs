use crate::error::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteOrder {
    Little,
    Big,
}

impl ByteOrder {
    pub fn read_i16(self, bytes: [u8; 2]) -> i16 {
        match self {
            Self::Little => i16::from_le_bytes(bytes),
            Self::Big => i16::from_be_bytes(bytes),
        }
    }

    pub fn read_i32(self, bytes: [u8; 4]) -> i32 {
        match self {
            Self::Little => i32::from_le_bytes(bytes),
            Self::Big => i32::from_be_bytes(bytes),
        }
    }

    pub fn read_f64(self, bytes: [u8; 8]) -> f64 {
        match self {
            Self::Little => f64::from_le_bytes(bytes),
            Self::Big => f64::from_be_bytes(bytes),
        }
    }
}

pub(super) struct Reader<'a> {
    pub(super) data: &'a [u8],
    pub(super) pos: usize,
    pub(super) order: ByteOrder,
}

impl<'a> Reader<'a> {
    pub(super) fn new(data: &'a [u8], pos: usize, order: ByteOrder) -> Self {
        Self { data, pos, order }
    }

    pub(super) fn ensure(&self, n: usize) -> Result<(), Error> {
        if self.pos + n > self.data.len() {
            Err(Error::FileFormat("unexpected end of file".into()))
        } else {
            Ok(())
        }
    }

    pub(super) fn read_i32(&mut self) -> Result<i32, Error> {
        self.ensure(4)?;
        let bytes: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(self.order.read_i32(bytes))
    }

    pub(super) fn read_f64(&mut self) -> Result<f64, Error> {
        self.ensure(8)?;
        let bytes: [u8; 8] = self.data[self.pos..self.pos + 8].try_into().unwrap();
        self.pos += 8;
        Ok(self.order.read_f64(bytes))
    }

    #[allow(dead_code)]
    pub(super) fn skip(&mut self, n: usize) -> Result<(), Error> {
        self.ensure(n)?;
        self.pos += n;
        Ok(())
    }

    #[allow(dead_code)]
    pub(super) fn seek(&mut self, pos: usize) {
        self.pos = pos;
    }
}

pub struct JplHeader {
    pub byte_order: ByteOrder,
    pub ss: [f64; 3],
    pub au: f64,
    pub emrat: f64,
    pub denum: i32,
    pub ncon: i32,
    pub ipt: [i32; 39],
    pub ksize: usize,
    pub irecsz: usize,
    pub ncoeffs: usize,
}

/// Detect byte order by reading ss[2] (segment length in days) at file offset 2668.
/// JPL files have no magic number; endianness is detected by value plausibility.
/// ss[2] must lie in [1.0, 200.0] days. (swejpl.c:217–226)
fn detect_byte_order(data: &[u8]) -> Result<ByteOrder, Error> {
    if data.len() < 2676 {
        return Err(Error::FileFormat("file too short for JPL header".into()));
    }
    // ss[2] starts at offset 2652 + 2*8 = 2668
    let b: [u8; 8] = data[2668..2676].try_into().unwrap();
    let ss2_le = f64::from_le_bytes(b);
    if (1.0..=200.0).contains(&ss2_le) {
        return Ok(ByteOrder::Little);
    }
    let ss2_be = f64::from_be_bytes(b);
    if (1.0..=200.0).contains(&ss2_be) {
        return Ok(ByteOrder::Big);
    }
    Err(Error::FileFormat(
        "cannot detect byte order: JPL ss[2] implausible in both endiannesses".into(),
    ))
}

/// Compute ksize from ipt[]. (swejpl.c:275–291)
fn compute_ksize(ipt: &[i32; 39]) -> Result<usize, Error> {
    let mut kmx = 0i32;
    let mut khi = 0usize;
    for i in 0..13usize {
        if ipt[i * 3] > kmx {
            kmx = ipt[i * 3];
            khi = i + 1; // 1-based
        }
    }
    if khi == 0 {
        return Err(Error::FileFormat(
            "JPL ipt table has all-zero offsets".into(),
        ));
    }
    // khi==12 → nutations (2 components); all others → 3 components
    let nd = if khi == 12 { 2i32 } else { 3i32 };
    let k1 = khi * 3; // 1-based, so ipt[k1-3] = ipt[(khi-1)*3] = body's starting offset
    let ksize_raw = (ipt[k1 - 3] + nd * ipt[k1 - 2] * ipt[k1 - 1] - 1) * 2;
    let mut ksize = ksize_raw as usize;
    if ksize == 1546 {
        ksize = 1652; // DE102: padded to match DE200 record size
    }
    if !(1000..=5000).contains(&ksize) {
        return Err(Error::FileFormat(format!(
            "computed JPL ksize {ksize} out of range [1000, 5000]"
        )));
    }
    Ok(ksize)
}

/// Parse the JPL DE file header from the mmap'd bytes. (swejpl.c:189–328, 668–730)
///
/// Record 0 layout (byte offsets):
///   0–251:   title (252 bytes)
///   252–2651: constant names (2400 bytes)
///   2652–2675: ss[3] (3 × f64)
///   2676–2679: ncon (i32)
///   2680–2687: au (f64)
///   2688–2695: emrat (f64)
///   2696–2839: ipt[0..35] (36 × i32)
///   2840–2843: numde (i32)
///   2844–2855: lpt[0..2] → ipt[36..38] (3 × i32)
pub(super) fn parse_header(data: &[u8]) -> Result<JplHeader, Error> {
    let byte_order = detect_byte_order(data)?;

    // Start reading at offset 2652 (after title + cnam blocks)
    let mut r = Reader::new(data, 2652, byte_order);

    let ss = [r.read_f64()?, r.read_f64()?, r.read_f64()?];
    let ncon = r.read_i32()?;
    let au = r.read_f64()?;
    let emrat = r.read_f64()?;

    let mut ipt = [0i32; 39];
    for item in ipt[..36].iter_mut() {
        *item = r.read_i32()?;
    }
    let denum = r.read_i32()?;
    for item in ipt[36..].iter_mut() {
        *item = r.read_i32()?; // lpt → ipt[36..38]
    }

    // Validate ss[0] and ss[1] plausibility (swejpl.c:228–236)
    const JD_MIN: f64 = -5_583_942.0;
    const JD_MAX: f64 = 9_025_909.0;
    if !(JD_MIN..=JD_MAX).contains(&ss[0]) {
        return Err(Error::FileFormat(format!(
            "JPL ss[0]={} outside plausibility range",
            ss[0]
        )));
    }
    if !(JD_MIN..=JD_MAX).contains(&ss[1]) {
        return Err(Error::FileFormat(format!(
            "JPL ss[1]={} outside plausibility range",
            ss[1]
        )));
    }

    let ksize = compute_ksize(&ipt)?;
    let irecsz = 4 * ksize;
    let ncoeffs = ksize / 2;

    Ok(JplHeader {
        byte_order,
        ss,
        au,
        emrat,
        denum,
        ncon,
        ipt,
        ksize,
        irecsz,
        ncoeffs,
    })
}

/// Validate the file length against the expected size derived from the header.
/// Accepts the exact expected size or expected + one extra record. (swejpl.c:732–762)
pub(super) fn validate_file_length(data: &[u8], header: &JplHeader) -> Result<(), Error> {
    let nseg = ((header.ss[1] - header.ss[0]) / header.ss[2]).round() as usize;

    // Sum coefficient doubles for all 13 bodies (nutations have 2 components, rest 3)
    let mut expected_doubles: usize = 0;
    for i in 0..13usize {
        let k = if i == 11 { 2usize } else { 3usize };
        let ncf = header.ipt[i * 3 + 1] as usize;
        let na = header.ipt[i * 3 + 2] as usize;
        expected_doubles += ncf * na * k * nseg;
    }
    expected_doubles += 2 * nseg; // buf[0] and buf[1] (segment start/end JD) per record

    let expected_bytes = expected_doubles * 8 + 2 * header.ksize * 4;
    let actual = data.len();

    if actual != expected_bytes && actual != expected_bytes + header.ksize * 4 {
        return Err(Error::FileFormat(format!(
            "JPL file length mismatch: expected {} or {}, got {}",
            expected_bytes,
            expected_bytes + header.ksize * 4,
            actual,
        )));
    }

    Ok(())
}
