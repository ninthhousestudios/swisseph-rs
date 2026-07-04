//! Send + Sync guarantee for `Ephemeris`.
//!
//! The future C-ABI/Dart FFI layer (docs/ffi-design.md) hands a `*mut Ephemeris` across
//! isolate boundaries and calls into one shared instance from many threads with zero
//! synchronization. That is only sound if `Ephemeris: Send + Sync` with no interior
//! mutability — which the stateless architecture guarantees today. These tests pin that
//! guarantee: the compile-time assertion breaks the build if a non-Sync field is ever
//! added, and the concurrency test asserts that parallel shared-reference use produces
//! bitwise-identical results to serial computation.

use swisseph::{Body, CalcFlags, CalcResult, Ephemeris, EphemerisConfig};

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn ephemeris_is_send_sync() {
    assert_send_sync::<Ephemeris>();
    assert_send_sync::<EphemerisConfig>();
}

const PLANETS: [Body; 10] = [
    Body::Sun,
    Body::Moon,
    Body::Mercury,
    Body::Venus,
    Body::Mars,
    Body::Jupiter,
    Body::Saturn,
    Body::Uranus,
    Body::Neptune,
    Body::Pluto,
];

/// 100 epochs spanning 1900–2037, well inside Moshier validity.
fn epochs() -> impl Iterator<Item = f64> {
    (0..100).map(|i| 2415020.5 + i as f64 * 500.0)
}

/// f64 payloads compared as raw bits: the guarantee is bitwise determinism,
/// not epsilon closeness.
fn to_bits(r: &CalcResult) -> ([u64; 6], u32) {
    (r.data.map(f64::to_bits), r.flags_used.bits())
}

fn compute_all(eph: &Ephemeris) -> Vec<([u64; 6], u32)> {
    let flags = CalcFlags::MOSEPH | CalcFlags::SPEED;
    let mut out = Vec::with_capacity(100 * PLANETS.len());
    for jd in epochs() {
        for body in PLANETS {
            let r = eph
                .calc_ut(jd, body, flags)
                .unwrap_or_else(|e| panic!("calc_ut({jd}, {body:?}) failed: {e}"));
            out.push(to_bits(&r));
        }
    }
    out
}

/// 8 threads sharing one `&Ephemeris` must each produce results bitwise-identical
/// to a serially computed baseline — the "zero synchronization, identical answers"
/// contract in executable form.
#[test]
fn shared_ephemeris_is_bitwise_deterministic_across_threads() {
    let eph = Ephemeris::new(EphemerisConfig::default()).expect("Ephemeris::new");
    let baseline = compute_all(&eph);

    std::thread::scope(|s| {
        let handles: Vec<_> = (0..8).map(|_| s.spawn(|| compute_all(&eph))).collect();
        for (i, h) in handles.into_iter().enumerate() {
            let results = h.join().expect("thread panicked");
            assert_eq!(
                results, baseline,
                "thread {i} diverged bitwise from the serial baseline"
            );
        }
    });
}
