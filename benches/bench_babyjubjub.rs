use criterion::{criterion_group, criterion_main, Criterion};

extern crate rand;
#[macro_use]
extern crate ff;
use ff::*;

extern crate num;
extern crate num_bigint;
use num_bigint::BigInt;

use babyjubjub_rs::{utils, Point};

fn criterion_benchmark(c: &mut Criterion) {
    // let x: BigInt = BigInt::parse_bytes(
    //     b"17777552123799933955779906779655732241715742912184938656739573121738514868268",
    //     10,
    // )
    // .unwrap();
    // c.bench_function("modulus", |b| {
    //     b.iter(|| utils::modulus(&x, &babyjubjub_rs::Q))
    // });

    let p: Point = Point {
        x: babyjubjub_rs::Fr::from_str(
            "17777552123799933955779906779655732241715742912184938656739573121738514868268",
        )
        .unwrap(),
        y: babyjubjub_rs::Fr::from_str(
            "2626589144620713026669568689430873010625803728049924121243784502389097019475",
        )
        .unwrap(),
        z: babyjubjub_rs::Fr::one(),
    };
    let q = p.clone();

    c.bench_function("add", |b| b.iter(|| p.add(&q)));
    let r: BigInt = BigInt::parse_bytes(b"3", 10).unwrap();
    c.bench_function("mul_scalar_small", |b| b.iter(|| p.mul_scalar(&r)));
    let r: BigInt = BigInt::parse_bytes(
        b"2626589144620713026669568689430873010625803728049924121243784502389097019475",
        10,
    )
    .unwrap();
    c.bench_function("mul_scalar", |b| b.iter(|| p.mul_scalar(&r)));

    // c.bench_function("compress", |b| b.iter(|| p.compress()));
    // let p_comp = p.compress();
    // c.bench_function("decompress", |b| {
    //     b.iter(|| babyjubjub_rs::decompress_point(p_comp))
    // });

    // let sk = babyjubjub_rs::new_key();
    // let pk = sk.public().unwrap();
    // let msg = 5.to_bigint().unwrap();
    // c.bench_function("sign_poseidon", |b| {
    //     b.iter(|| sk.sign_poseidon(msg.clone()))
    // });
    // let sig = sk.sign_poseidon(msg.clone()).unwrap();
    // c.bench_function("verify_poseidon", |b| {
    //     b.iter(|| babyjubjub_rs::verify_poseidon(pk.clone(), sig.clone(), msg.clone()))
    // });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
