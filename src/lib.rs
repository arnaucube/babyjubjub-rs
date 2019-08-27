#[macro_use]
extern crate arrayref;
extern crate generic_array;
extern crate mimc_rs;
extern crate num;
extern crate num_bigint;
extern crate num_traits;
extern crate rand;

use blake2::{Blake2b, Digest};
use mimc_rs::Mimc7;

use num_bigint::{BigInt, RandBigInt, Sign, ToBigInt};
use num_traits::{One, Zero};

use generic_array::GenericArray;

mod utils;

#[derive(Clone, Debug)]
pub struct Point {
    pub x: BigInt,
    pub y: BigInt,
}

pub struct Signature {
    r_b8: Point,
    s: BigInt,
}

pub struct PrivateKey {
    bbjj: Babyjubjub,
    key: BigInt,
}

impl PrivateKey {
    pub fn public(&self) -> Point {
        // https://tools.ietf.org/html/rfc8032#section-5.1.5
        let pk = &self.bbjj.mul_scalar(self.bbjj.b8.clone(), self.key.clone());
        pk.clone()
    }

    pub fn sign(&self, msg: BigInt) -> Signature {
        // https://tools.ietf.org/html/rfc8032#section-5.1.6
        let mut hasher = Blake2b::new();
        let (_, sk_bytes) = self.key.to_bytes_be();
        hasher.input(sk_bytes);
        let mut h = hasher.result(); // h: hash(sk)
                                     // s: h[32:64]
        let s = GenericArray::<u8, generic_array::typenum::U32>::from_mut_slice(&mut h[32..64]);
        let (_, msg_bytes) = msg.to_bytes_be();
        let r_bytes = utils::concatenate_arrays(s, &msg_bytes);
        let mut r = BigInt::from_bytes_be(Sign::Plus, &r_bytes[..]);
        r = utils::modulus(&r, &self.bbjj.sub_order);
        let r8: Point = self.bbjj.mul_scalar(self.bbjj.b8.clone(), r.clone());
        let a = &self.public();

        let hm_input = vec![r8.x.clone(), r8.y.clone(), a.x.clone(), a.y.clone(), msg];
        let mimc7 = Mimc7::new();
        let hm = mimc7.hash(hm_input);

        let mut s = &self.key << 3;
        s = hm * s;
        s = r + s;
        s = s % &self.bbjj.sub_order;

        Signature {
            r_b8: r8.clone(),
            s: s,
        }
    }
}

pub struct Babyjubjub {
    d: BigInt,
    a: BigInt,
    q: BigInt,
    b8: Point,
    // order: BigInt,
    sub_order: BigInt,
}

impl Babyjubjub {
    pub fn new() -> Babyjubjub {
        let d: BigInt = BigInt::parse_bytes(b"168696", 10).unwrap();
        let a: BigInt = BigInt::parse_bytes(b"168700", 10).unwrap();
        let q: BigInt = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();
        let b8: Point = Point {
            x: BigInt::parse_bytes(
                b"5299619240641551281634865583518297030282874472190772894086521144482721001553",
                10,
            )
            .unwrap(),
            y: BigInt::parse_bytes(
                b"16950150798460657717958625567821834550301663161624707787222815936182638968203",
                10,
            )
            .unwrap(),
        };
        let order: BigInt = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088614511777268538073601725287587578984328",
            10,
        )
        .unwrap();
        let sub_order: BigInt = &order >> 3;

        Babyjubjub {
            d: d,
            a: a,
            q: q,
            b8: b8,
            // order: order,
            sub_order: sub_order,
        }
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        // x = (x1*y2+y1*x2)/(c*(1+d*x1*x2*y1*y2))
        // y = (y1*y2-x1*x2)/(c*(1-d*x1*x2*y1*y2))

        // x = (x1 * y2 + y1 * x2) / (1 + d * x1 * y1 * y2)
        let one: BigInt = One::one();
        let x_num: BigInt = &p.x * &q.y + &p.y * &q.x;
        let x_den: BigInt = &one + &self.d * &p.x * &q.x * &p.y * &q.y;
        let x_den_inv = utils::modinv(&x_den, &self.q);
        let x: BigInt = utils::modulus(&(&x_num * &x_den_inv), &self.q);

        // y = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)
        let y_num = &p.y * &q.y - &self.a * &p.x * &q.x;
        let y_den = utils::modulus(&(&one - &self.d * &p.x * &q.x * &p.y * &q.y), &self.q);
        let y_den_inv = utils::modinv(&y_den, &self.q);
        let y: BigInt = utils::modulus(&(&y_num * &y_den_inv), &self.q);

        Point { x: x, y: y }
    }

    pub fn mul_scalar(&self, p: Point, n: BigInt) -> Point {
        // TODO use & in p and n to avoid clones on function call
        let mut r: Point = Point {
            x: Zero::zero(),
            y: One::one(),
        };
        let mut rem: BigInt = n;
        let mut exp: Point = p;

        let zero: BigInt = Zero::zero();
        let one: BigInt = One::one();
        while rem != zero {
            let is_odd = &rem & &one == one;
            if is_odd == true {
                r = self.add(&r, &exp);
            }
            exp = self.add(&exp, &exp);
            rem = rem >> 1;
        }
        r.x = utils::modulus(&r.x, &self.q);
        r.y = utils::modulus(&r.y, &self.q);
        r
    }

    pub fn compress_point(&self, p: &Point) -> [u8; 32] {
        let mut r: [u8; 32];
        let (_, y_bytes) = p.y.to_bytes_le();
        r = *array_ref!(y_bytes, 0, 32);
        if &p.x > &(&self.q >> 1) {
            r[31] = r[31] | 0x80;
        }
        r
    }

    pub fn decompress_point(&self, bb: [u8; 32]) -> Point {
        // https://tools.ietf.org/html/rfc8032#section-5.2.3
        let mut sign: bool = false;
        let mut b = bb.clone();
        if b[31] & 0x80 != 0x00 {
            sign = true;
            b[31] = b[31] & 0x7F;
        }
        let y: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[..]);
        if y >= self.q {
            // println!("ERROR0");
        }
        let one: BigInt = One::one();

        // x^2 = (1 - y^2) / (a - d * y^2) (mod p)
        let mut x: BigInt = utils::modulus(
            &((one - utils::modulus(&(&y * &y), &self.q))
                * utils::modinv(
                    &utils::modulus(
                        &(&self.a - utils::modulus(&(&self.d * (&y * &y)), &self.q)),
                        &self.q,
                    ),
                    &self.q,
                )),
            &self.q,
        );
        x = utils::modsqrt(&x, &self.q);

        if (sign && x >= Zero::zero()) || (!sign && x < Zero::zero()) {
            x = x * -1.to_bigint().unwrap();
        }
        x = utils::modulus(&x, &self.q);
        Point { x: x, y: y }
    }

    pub fn compress_sig(&self, sig: &Signature) -> [u8; 64] {
        let mut b: Vec<u8> = Vec::new();
        b.append(&mut self.compress_point(&sig.r_b8).to_vec());
        // let (_, mut s_bytes) = sig.s.to_bytes_le();
        let (_, mut s_bytes) = sig.s.to_bytes_le();
        println!("sbytes LENGTH {:?}", s_bytes.len());
        // let mut s_32bytes: [u8; 32] = [0; 32];
        // s_32bytes[..].copy_from_slice(&s_bytes[..]);
        b.append(&mut s_bytes);
        let mut r: [u8; 64] = [0; 64];
        // r = *array_ref!(b, 0, 64);
        // r.copy_from_slice(&b[..]);
        println!("b LENGTH {:?}", b.len());
        // if b.len() < 64 {
        //     // let diff = 64 - b.len();
        //     println!("less than 64, add padding");
        //     let e: [u8; 1] = [0];
        //     b.append(&mut e.to_vec());
        // }
        r.copy_from_slice(&b[..]);
        println!("r {:?}", r.len());
        r
    }

    pub fn decompress_sig(&self, b: &[u8; 64]) -> Signature {
        let r_b8_bytes: [u8; 32] = *array_ref!(b[..32], 0, 32);
        let s: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[32..]);
        let r_b8 = &self.decompress_point(r_b8_bytes);
        Signature {
            r_b8: r_b8.clone(),
            s: s,
        }
    }

    pub fn new_key(&self) -> PrivateKey {
        // https://tools.ietf.org/html/rfc8032#section-5.1.5
        let mut rng = rand::thread_rng();
        let sk_raw = rng.gen_biguint(1024).to_bigint().unwrap();

        let mut hasher = Blake2b::new();
        let (_, sk_raw_bytes) = sk_raw.to_bytes_be();
        hasher.input(sk_raw_bytes);
        let mut h = hasher.result();

        h[0] = h[0] & 0xF8;
        h[31] = h[31] & 0x7F;
        h[31] = h[31] | 0x40;

        let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);

        let bbjj_new = Babyjubjub {
            d: self.d.clone(),
            a: self.a.clone(),
            q: self.q.clone(),
            b8: self.b8.clone(),
            sub_order: self.sub_order.clone(),
        };
        PrivateKey {
            bbjj: bbjj_new,
            key: sk,
        }
    }

    pub fn verify(&self, pk: Point, sig: Signature, msg: BigInt) -> bool {
        let hm_input = vec![
            sig.r_b8.x.clone(),
            sig.r_b8.y.clone(),
            pk.x.clone(),
            pk.y.clone(),
            msg,
        ];
        let mimc7 = Mimc7::new();
        let hm = mimc7.hash(hm_input);
        let l = self.mul_scalar(self.b8.clone(), sig.s);
        let r = self.add(&sig.r_b8, &self.mul_scalar(pk, 8.to_bigint().unwrap() * hm));
        if l.x == r.x && l.y == r.y {
            return true;
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate rustc_hex;
    use rustc_hex::ToHex;

    #[test]
    fn test_add_same_point() {
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
        let q: Point = Point {
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
        let bbjj = Babyjubjub::new();
        let res = bbjj.add(&p, &q);
        assert_eq!(
            res.x.to_string(),
            "6890855772600357754907169075114257697580319025794532037257385534741338397365"
        );
        assert_eq!(
            res.y.to_string(),
            "4338620300185947561074059802482547481416142213883829469920100239455078257889"
        );
    }
    #[test]
    fn test_add_different_points() {
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
        let q: Point = Point {
            x: BigInt::parse_bytes(
                b"16540640123574156134436876038791482806971768689494387082833631921987005038935",
                10,
            )
            .unwrap(),
            y: BigInt::parse_bytes(
                b"20819045374670962167435360035096875258406992893633759881276124905556507972311",
                10,
            )
            .unwrap(),
        };
        let bbjj = Babyjubjub::new();
        let res = bbjj.add(&p, &q);
        assert_eq!(
            res.x.to_string(),
            "7916061937171219682591368294088513039687205273691143098332585753343424131937"
        );
        assert_eq!(
            res.y.to_string(),
            "14035240266687799601661095864649209771790948434046947201833777492504781204499"
        );
    }

    #[test]
    fn test_mul_scalar() {
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
        let bbjj = Babyjubjub::new();
        let res_m = bbjj.mul_scalar(p.clone(), 3.to_bigint().unwrap());
        let res_a = bbjj.add(&p, &p);
        let res_a = bbjj.add(&res_a, &p);
        assert_eq!(res_m.x, res_a.x);
        assert_eq!(
            res_m.x.to_string(),
            "19372461775513343691590086534037741906533799473648040012278229434133483800898"
        );
        assert_eq!(
            res_m.y.to_string(),
            "9458658722007214007257525444427903161243386465067105737478306991484593958249"
        );

        let n = BigInt::parse_bytes(
            b"14035240266687799601661095864649209771790948434046947201833777492504781204499",
            10,
        )
        .unwrap();
        let res2 = bbjj.mul_scalar(p.clone(), n);
        assert_eq!(
            res2.x.to_string(),
            "17070357974431721403481313912716834497662307308519659060910483826664480189605"
        );
        assert_eq!(
            res2.y.to_string(),
            "4014745322800118607127020275658861516666525056516280575712425373174125159339"
        );
    }

    #[test]
    fn test_point_compress_decompress() {
        let bbjj = Babyjubjub::new();
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
        let p_comp = bbjj.compress_point(&p);
        assert_eq!(
            p_comp[..].to_hex(),
            "53b81ed5bffe9545b54016234682e7b2f699bd42a5e9eae27ff4051bc698ce85"
        );
        let p2 = bbjj.decompress_point(p_comp);
        assert_eq!(p.x, p2.x);
        assert_eq!(p.y, p2.y);
    }

    #[test]
    fn test_new_key_sign_verify0() {
        let bbjj = Babyjubjub::new();
        let sk = bbjj.new_key();
        let pk = sk.public();
        let msg = 5.to_bigint().unwrap();
        let sig = sk.sign(msg.clone());
        let v = bbjj.verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_new_key_sign_verify1() {
        let bbjj = Babyjubjub::new();
        let sk = bbjj.new_key();
        let pk = sk.public();
        let msg = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();
        let sig = sk.sign(msg.clone());
        let v = bbjj.verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_signature_compress_decompress() {
        let bbjj = Babyjubjub::new();
        let sk = bbjj.new_key();
        let pk = sk.public();
        let msg = 5.to_bigint().unwrap();
        let sig = sk.sign(msg.clone());

        let compressed_sig = bbjj.compress_sig(&sig);
        println!("compressedsig {:?}", compressed_sig.to_hex());
        let decompressed_sig = bbjj.decompress_sig(&compressed_sig);
        assert_eq!(&sig.r_b8.x, &decompressed_sig.r_b8.x);
        assert_eq!(&sig.r_b8.y, &decompressed_sig.r_b8.y);
        assert_eq!(&sig.s, &decompressed_sig.s);

        let v = bbjj.verify(pk, decompressed_sig, msg);
        assert_eq!(v, true);
    }
}
