use criterion::{criterion_group, criterion_main, Criterion};

extern crate rand;
#[macro_use]
extern crate ff;
use ff::*;

extern crate num;
extern crate num_bigint;
use num_bigint::{BigInt, ToBigInt};

use babyjubjub_rs::{utils, Point};

fn criterion_benchmark(c: &mut Criterion) {
    let p: Point = Point {
        x: babyjubjub_rs::Fr::from_str(
            "17777552123799933955779906779655732241715742912184938656739573121738514868268",
        )
        .unwrap(),
        y: babyjubjub_rs::Fr::from_str(
            "2626589144620713026669568689430873010625803728049924121243784502389097019475",
        )
        .unwrap(),
    };
    let q = p.clone();

    let p_projective = p.projective();
    let q_projective = q.projective();

    c.bench_function("add", |b| b.iter(|| p_projective.add(&q_projective)));
    let r: BigInt = BigInt::parse_bytes(b"3", 10).unwrap();
    c.bench_function("mul_scalar_small", |b| b.iter(|| p.mul_scalar(&r)));
    let r: BigInt = BigInt::parse_bytes(
        b"2626589144620713026669568689430873010625803728049924121243784502389097019475",
        10,
    )
    .unwrap();
    c.bench_function("mul_scalar", |b| b.iter(|| p.mul_scalar(&r)));

    c.bench_function("point compress", |b| b.iter(|| p.compress()));
    let p_comp = p.compress();
    c.bench_function("point decompress", |b| {
        b.iter(|| babyjubjub_rs::decompress_point(p_comp))
    });

    let sk = babyjubjub_rs::new_key();
    let pk = sk.public().unwrap();
    let msg = 5.to_bigint().unwrap();
    c.bench_function("sign", |b| b.iter(|| sk.sign(msg.clone())));
    let sig = sk.sign(msg.clone()).unwrap();
    c.bench_function("verify", |b| {
        b.iter(|| babyjubjub_rs::verify(pk.clone(), sig.clone(), msg.clone()))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
