extern crate num;
extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::{One, Zero};

pub fn modulus(a: &BigInt, m: &BigInt) -> BigInt {
    ((a % m) + m) % m
}

pub fn mod_inverse0(a: &BigInt, q: &BigInt) -> BigInt {
    let mut mn = (q.clone(), a.clone());
    let mut xy: (BigInt, BigInt) = (Zero::zero(), One::one());

    let big_zero: BigInt = Zero::zero();
    while mn.1 != big_zero {
        xy = (xy.1.clone(), xy.0 - (mn.0.clone() / mn.1.clone()) * xy.1);
        mn = (mn.1.clone(), modulus(&mn.0, &mn.1));
    }

    while xy.0 < Zero::zero() {
        xy.0 = modulus(&xy.0, q);
    }
    xy.0
}

/*
pub fn mod_inverse1(a0: BigInt, m0: BigInt) -> BigInt {
    if m0 == One::one() {
        return One::one();
    }

    let (mut a, mut m, mut x0, mut inv): (BigInt, BigInt, BigInt, BigInt) =
        (a0, m0.clone(), Zero::zero(), One::one());

    while a > One::one() {
        inv = inv - (&a / m.clone()) * x0.clone();
        a = a % m.clone();
        std::mem::swap(&mut a, &mut m);
        std::mem::swap(&mut x0, &mut inv);
    }

    if inv < Zero::zero() {
        inv += m0.clone()
    }
    inv
}

pub fn mod_inverse2(a: BigInt, q: BigInt) -> BigInt {
    let mut aa: BigInt = a;
    let mut qq: BigInt = q;
    if qq < Zero::zero() {
        qq = -qq;
    }
    if aa < Zero::zero() {
        aa = -aa;
    }
    let d = num::Integer::gcd(&aa, &qq);
    if d != One::one() {
        println!("ERR no mod_inv");
    }
    let res: BigInt;
    if d < Zero::zero() {
        res = d + qq;
    } else {
        res = d;
    }
    res
}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_inverse() {
        let a = BigInt::parse_bytes(b"123456789123456789123456789123456789123456789", 10).unwrap();
        let b = BigInt::parse_bytes(b"12345678", 10).unwrap();
        assert_eq!(
            mod_inverse0(&a, &b),
            BigInt::parse_bytes(b"641883", 10).unwrap()
        );
    }
}
