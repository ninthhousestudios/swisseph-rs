use criterion::{Criterion, black_box, criterion_group, criterion_main};
use std::hint::black_box as std_bb;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
use swisseph::CalcFlags;
use swisseph::nutation;
use swisseph::types::{AstroModels, NutationModel};
use swisseph::{Body, EclipseFlags, Ephemeris, EphemerisConfig, HouseSystem};

const J2000: f64 = 2451545.0;

const BODIES: [Body; 10] = [
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

fn moshier_eph() -> Ephemeris {
    Ephemeris::new(EphemerisConfig::default()).unwrap()
}

fn swiss_eph() -> Option<Ephemeris> {
    use swisseph::types::EphemerisSource;
    let ephe_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../swisseph/ephe");
    if !ephe_path.exists() {
        return None;
    }
    let config = EphemerisConfig {
        ephemeris_source: EphemerisSource::Swiss,
        ephe_path: Some(ephe_path),
        ..Default::default()
    };
    Ephemeris::new(config).ok()
}

fn bench_calc_moshier(c: &mut Criterion) {
    let eph = moshier_eph();
    c.bench_function("calc_moshier_full_chart", |b| {
        b.iter(|| {
            for &body in &BODIES {
                let _ = black_box(eph.calc_ut(black_box(J2000), body, CalcFlags::SPEED));
            }
        });
    });
}

fn bench_calc_swiss(c: &mut Criterion) {
    let Some(eph) = swiss_eph() else {
        eprintln!("SKIP: calc_swiss — Swiss Ephemeris files not found at ../swisseph/ephe");
        return;
    };
    c.bench_function("calc_swiss_full_chart", |b| {
        b.iter(|| {
            for &body in &BODIES {
                let _ = black_box(eph.calc_ut(black_box(J2000), body, CalcFlags::SPEED));
            }
        });
    });
}

fn bench_houses_placidus(c: &mut Criterion) {
    let eph = moshier_eph();
    c.bench_function("houses_placidus", |b| {
        b.iter(|| {
            let _ = black_box(eph.houses_ex2(
                black_box(J2000),
                CalcFlags::empty(),
                black_box(47.37),
                black_box(8.55),
                HouseSystem::Placidus,
            ));
        });
    });
}

fn bench_eclipse_search(c: &mut Criterion) {
    let eph = moshier_eph();
    c.bench_function("sol_eclipse_when_glob", |b| {
        b.iter(|| {
            let _ = black_box(eph.sol_eclipse_when_glob(
                black_box(J2000),
                CalcFlags::empty(),
                EclipseFlags::empty(),
                false,
            ));
        });
    });
}

fn bench_nutation_iau2000b(c: &mut Criterion) {
    let models = AstroModels {
        nutation: NutationModel::IAU2000B,
        ..Default::default()
    };
    c.bench_function("nutation_iau2000b", |b| {
        b.iter(|| {
            black_box(nutation::nutation(
                black_box(J2000),
                CalcFlags::empty(),
                &models,
            ));
        });
    });
}

fn bench_moshier_moon(c: &mut Criterion) {
    c.bench_function("moshier_moon_series", |b| {
        b.iter(|| {
            black_box(swisseph::moshier::moon::moshmoon2(black_box(J2000)));
        });
    });
}

fn bench_thread_scaling(c: &mut Criterion) {
    let eph = Arc::new(moshier_eph());
    let total_epochs = 10_000;

    for n_threads in [1, 2, 4, 8] {
        c.bench_function(&format!("thread_scaling_{n_threads}t"), |b| {
            b.iter_custom(|iters| {
                let mut total = std::time::Duration::ZERO;
                for _ in 0..iters {
                    let eph = Arc::clone(&eph);
                    let start = Instant::now();
                    std::thread::scope(|s| {
                        let chunk = total_epochs / n_threads;
                        for t in 0..n_threads {
                            let eph = &eph;
                            s.spawn(move || {
                                let base_jd = J2000 + (t * chunk) as f64;
                                for i in 0..chunk {
                                    let jd = base_jd + i as f64;
                                    for &body in &BODIES {
                                        let _ =
                                            std_bb(eph.calc_ut(std_bb(jd), body, CalcFlags::SPEED));
                                    }
                                }
                            });
                        }
                    });
                    total += start.elapsed();
                }
                total
            });
        });
    }
}

criterion_group!(
    benches,
    bench_calc_moshier,
    bench_calc_swiss,
    bench_houses_placidus,
    bench_eclipse_search,
    bench_nutation_iau2000b,
    bench_moshier_moon,
    bench_thread_scaling,
);
criterion_main!(benches);
