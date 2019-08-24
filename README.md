# babyjubjub-rs [![Crates.io](https://img.shields.io/crates/v/babyjubjub-rs.svg)](https://crates.io/crates/babyjubjub-rs) [![Build Status](https://travis-ci.org/arnaucube/babyjubjub-rs.svg?branch=master)](https://travis-ci.org/arnaucube/babyjubjub-rs)
BabyJubJub elliptic curve implementation in Rust

Uses MiMC7 hash function: https://github.com/arnaucube/mimc-rs

## Warning
Doing this in my free time to get familiar with Rust, do not use in production.

- [x] point addition
- [x] point scalar multiplication
- [x] eddsa keys generation
- [x] eddsa signature
- [x] eddsa signature verification
- [ ] {point, pk, signature} compress&decompress parsers




### References
- JubJub curve explanation: https://z.cash/technology/jubjub/
	- Rust: https://github.com/zkcrypto/jubjub
	- Python: https://github.com/daira/jubjub
- BabyJubJub curve:
	- C++ https://github.com/barryWhiteHat/baby_jubjub_ecc
	- Javascript & Circom: https://github.com/iden3/circomlib
	- Go https://github.com/iden3/go-iden3-crypto
