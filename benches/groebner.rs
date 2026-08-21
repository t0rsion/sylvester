//! Wall-time benchmarks for both backends on three small systems at
//! p = 32003, to catch regressions between runs. No cell here is certified
//! or parallel.

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use sylvester::{Backend, ComputeOptions, Ideal, PolynomialRing};

fn ideal(modulus: u64, variables: &[&str], generators: &[&str]) -> Ideal {
    let ring = PolynomialRing::prime_field(modulus, variables.iter().copied())
        .expect("the modulus is prime");
    let polynomials: Vec<_> = generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the text parses"))
        .collect();
    ring.ideal(polynomials)
        .expect("the polynomials share a ring")
}

fn cyclic_3(modulus: u64) -> Ideal {
    ideal(
        modulus,
        &["x", "y", "z"],
        &["x + y + z", "x*y + y*z + z*x", "x*y*z - 1"],
    )
}

fn cyclic_4(modulus: u64) -> Ideal {
    ideal(
        modulus,
        &["x", "y", "z", "w"],
        &[
            "x + y + z + w",
            "x*y + y*z + z*w + w*x",
            "x*y*z + y*z*w + z*w*x + w*x*y",
            "x*y*z*w - 1",
        ],
    )
}

fn katsura_4(modulus: u64) -> Ideal {
    ideal(
        modulus,
        &["a", "b", "c", "d", "e"],
        &[
            "a + 2*b + 2*c + 2*d + 2*e - 1",
            "2*a*b - b",
            "2*a*c + b^2 - c",
            "2*a*d + 2*b*c - d",
            "2*a*e + 2*b*d + c^2 - e",
        ],
    )
}

fn bench_backends(c: &mut Criterion) {
    let modulus = 32003;
    let systems = [
        ("cyclic-3", cyclic_3(modulus)),
        ("cyclic-4", cyclic_4(modulus)),
        ("katsura-4", katsura_4(modulus)),
    ];

    for (name, ideal) in &systems {
        for backend in [Backend::Matrix, Backend::Classic] {
            let label = format!("{name}/{}", backend_name(backend));
            c.bench_function(&label, |b| {
                b.iter(|| {
                    let basis = ideal
                        .groebner_basis(ComputeOptions::new().backend(backend))
                        .expect("a run without a budget cannot stop early");
                    black_box(basis.len())
                })
            });
        }
    }
}

fn backend_name(backend: Backend) -> &'static str {
    match backend {
        Backend::Matrix => "matrix",
        Backend::Classic => "classic",
    }
}

criterion_group!(benches, bench_backends);
criterion_main!(benches);
