use criterion::{criterion_group, criterion_main, Criterion};

extern crate num;
extern crate num_bigint;
extern crate num_traits;
use num_bigint::{BigInt, ToBigInt};

use babyjubjub_rs::{utils, Point};

fn criterion_benchmark(c: &mut Criterion) {
    let x: BigInt = BigInt::parse_bytes(
        b"17777552123799933955779906779655732241715742912184938656739573121738514868268",
        10,
    )
    .unwrap();
    c.bench_function("modulus", |b| {
        b.iter(|| utils::modulus(&x, &babyjubjub_rs::Q))
    });

    let p: Point = Point {
        x: BigInt::parse_bytes(
            b"17777552123799933955779906779655732241715742912184938656739573121738514868268",
            10,
        )
        .unwrap(),
        y: BigInt::parse_bytes(
            b"2626589144620713026669568689430873010625803728049924121243784502389097019475",
            10,
        )
        .unwrap(),
    };
    let q = p.clone();

    c.bench_function("add", |b| b.iter(|| p.add(&q)));

    c.bench_function("mul_scalar_small", |b| {
        b.iter(|| p.mul_scalar(&3.to_bigint().unwrap()))
    });
    let r: BigInt = BigInt::parse_bytes(
        b"2626589144620713026669568689430873010625803728049924121243784502389097019475",
        10,
    )
    .unwrap();
    c.bench_function("mul_scalar", |b| b.iter(|| p.mul_scalar(&r)));

    let sk = babyjubjub_rs::new_key();
    let pk = sk.public().unwrap();
    let msg = 5.to_bigint().unwrap();
    c.bench_function("sign_poseidon", |b| {
        b.iter(|| sk.sign_poseidon(msg.clone()))
    });
    let sig = sk.sign_poseidon(msg.clone()).unwrap();
    c.bench_function("verify_poseidon", |b| {
        b.iter(|| babyjubjub_rs::verify_poseidon(pk.clone(), sig.clone(), msg.clone()))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
