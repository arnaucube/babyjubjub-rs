// WARNING still updating the code, it works, but is still in process the refactor.

extern crate rand;
#[macro_use]
extern crate ff;
use ff::*;

use poseidon_rs::Poseidon;
pub type Fr = poseidon_rs::Fr; // alias

#[macro_use]
extern crate arrayref;
extern crate generic_array;
extern crate num;
extern crate num_bigint;
extern crate num_traits;

extern crate rand6;
use rand6::Rng;

use blake2::{Blake2b, Digest};
use std::cmp::min;

use num_bigint::{BigInt, RandBigInt, RandomBits, Sign, ToBigInt};
use num_traits::{One, Zero};

use generic_array::GenericArray;

pub mod utils;

#[macro_use]
extern crate lazy_static;

lazy_static! {
    static ref D: Fr = Fr::from_str("168696").unwrap();
    static ref D_big: BigInt = BigInt::parse_bytes(b"168696", 10).unwrap();
    static ref A: Fr = Fr::from_str("168700").unwrap();
    static ref A_big: BigInt = BigInt::parse_bytes(b"168700", 10).unwrap();
    pub static ref Q: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",10
    )
    .unwrap();
    // pub static ref Q: Fr = Fr::from_str(
    //     "21888242871839275222246405745257275088548364400416034343698204186575808495617"
    // )
    // .unwrap();
    static ref B8: Point = Point {
        x: Fr::from_str(
            "5299619240641551281634865583518297030282874472190772894086521144482721001553",
        )
        .unwrap(),
        y: Fr::from_str(
            "16950150798460657717958625567821834550301663161624707787222815936182638968203",
        )
        .unwrap(),
            // z: Fr::one(),
    };
    static ref ORDER: Fr = Fr::from_str(
        "21888242871839275222246405745257275088614511777268538073601725287587578984328",
    )
    .unwrap();

    // SUBORDER = ORDER >> 3
    // static ref SUBORDER: Fr = Fr::from_str(
    //     "i2736030358979909402780800718157159386076813972158567259200215660948447373041",
    // )
    // .unwrap();
    static ref SUBORDER: BigInt = &BigInt::parse_bytes(
                b"21888242871839275222246405745257275088614511777268538073601725287587578984328",
                        10,
                            )
            .unwrap()
                    >> 3;
    static ref poseidon: poseidon_rs::Poseidon = Poseidon::new();
}

#[derive(Clone, Debug)]
pub struct PointProjective {
    pub x: Fr,
    pub y: Fr,
    pub z: Fr,
}

impl PointProjective {
    pub fn affine(&self) -> Point {
        if self.z.is_zero() {
            return Point {
                x: Fr::zero(),
                y: Fr::zero(),
            };
        }

        let mut zinv = self.z.inverse().unwrap();
        let mut x = self.x;
        x.mul_assign(&zinv);
        let mut y = self.y;
        y.mul_assign(&zinv);

        Point {
            x: x.clone(),
            y: y.clone(),
        }
    }
    pub fn add(&self, q: &PointProjective) -> Result<PointProjective, String> {
        // add-2008-bbjlp https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp
        let mut a = self.z;
        a.mul_assign(&q.z);
        let mut b = a;
        b.square();
        let mut c = self.x;
        c.mul_assign(&q.x);
        let mut d = self.y;
        d.mul_assign(&q.y);
        let mut e = D.clone();
        e.mul_assign(&c);
        e.mul_assign(&d);
        let mut f = b;
        f.sub_assign(&e);
        let mut g = b;
        g.add_assign(&e);
        let mut x1y1 = self.x;
        x1y1.add_assign(&self.y);
        let mut x2y2 = q.x;
        x2y2.add_assign(&q.y);
        let mut aux = x1y1;
        aux.mul_assign(&x2y2);
        aux.sub_assign(&c);
        aux.sub_assign(&d);
        let mut x3 = a;
        x3.mul_assign(&f);
        x3.mul_assign(&aux);
        let mut ac = A.clone();
        ac.mul_assign(&c);
        let mut dac = d;
        dac.sub_assign(&ac);
        let mut y3 = a;
        y3.mul_assign(&g);
        y3.mul_assign(&dac);
        let mut z3 = f;
        z3.mul_assign(&g);

        Ok(PointProjective {
            x: x3.clone(),
            y: y3.clone(),
            z: z3.clone(),
        })
    }
}

#[derive(Clone, Debug)]
pub struct Point {
    pub x: Fr,
    pub y: Fr,
}

impl Point {
    pub fn projective(&self) -> PointProjective {
        PointProjective {
            x: self.x.clone(),
            y: self.y.clone(),
            z: Fr::one(),
        }
    }

    pub fn mul_scalar(&self, n: &BigInt) -> Result<Point, String> {
        let mut r: PointProjective = PointProjective {
            x: Fr::zero(),
            y: Fr::one(),
            z: Fr::one(),
        };
        let mut exp: PointProjective = self.projective();
        let (_, b) = n.to_bytes_le();
        for i in 0..n.bits() {
            if test_bit(&b, i) {
                r = r.add(&exp)?;
            }
            exp = exp.add(&exp)?;
        }
        Ok(r.affine())
    }

    pub fn compress(&self) -> [u8; 32] {
        let p = &self;
        let mut r: [u8; 32] = [0; 32];
        let x_big = BigInt::parse_bytes(to_hex(&p.x).as_bytes(), 16).unwrap();
        let y_big = BigInt::parse_bytes(to_hex(&p.y).as_bytes(), 16).unwrap();
        let (_, y_bytes) = y_big.to_bytes_le();
        let len = min(y_bytes.len(), r.len());
        r[..len].copy_from_slice(&y_bytes[..len]);
        if &x_big > &(&Q.clone() >> 1) {
            r[31] = r[31] | 0x80;
        }
        r
    }

    pub fn equals(&self, p: Point) -> bool {
        if self.x == p.x && self.y == p.y {
            return true;
        }
        false
    }
}

pub fn test_bit(b: &Vec<u8>, i: usize) -> bool {
    return b[i / 8] & (1 << (i % 8)) != 0;
}

pub fn decompress_point(bb: [u8; 32]) -> Result<Point, String> {
    // https://tools.ietf.org/html/rfc8032#section-5.2.3
    let mut sign: bool = false;
    let mut b = bb.clone();
    if b[31] & 0x80 != 0x00 {
        sign = true;
        b[31] = b[31] & 0x7F;
    }
    let y: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[..]);
    if y >= Q.clone() {
        return Err("y outside the Finite Field over R".to_string());
    }
    let one: BigInt = One::one();

    // x^2 = (1 - y^2) / (a - d * y^2) (mod p)
    let den = utils::modinv(
        &utils::modulus(
            &(&A_big.clone() - utils::modulus(&(&D_big.clone() * (&y * &y)), &Q)),
            &Q,
        ),
        &Q,
    )?;
    let mut x: BigInt = utils::modulus(&((one - utils::modulus(&(&y * &y), &Q)) * den), &Q);
    x = utils::modsqrt(&x, &Q)?;

    if sign && !(&x > &(&Q.clone() >> 1)) || (!sign && (&x > &(&Q.clone() >> 1))) {
        x = x * -1.to_bigint().unwrap();
    }
    x = utils::modulus(&x, &Q);
    let x_fr: Fr = Fr::from_str(&x.to_string()).unwrap();
    let y_fr: Fr = Fr::from_str(&y.to_string()).unwrap();
    Ok(Point { x: x_fr, y: y_fr })
}

#[derive(Debug, Clone)]
pub struct Signature {
    r_b8: Point,
    s: BigInt,
}

impl Signature {
    pub fn compress(&self) -> [u8; 64] {
        let mut b: Vec<u8> = Vec::new();
        b.append(&mut self.r_b8.compress().to_vec());
        let (_, s_bytes) = self.s.to_bytes_le();
        let mut s_32bytes: [u8; 32] = [0; 32];
        let len = min(s_bytes.len(), s_32bytes.len());
        s_32bytes[..len].copy_from_slice(&s_bytes[..len]);
        b.append(&mut s_32bytes.to_vec());
        let mut r: [u8; 64] = [0; 64];
        r[..].copy_from_slice(&b[..]);
        r
    }
}

pub fn decompress_signature(b: &[u8; 64]) -> Result<Signature, String> {
    let r_b8_bytes: [u8; 32] = *array_ref!(b[..32], 0, 32);
    let s: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[32..]);
    let r_b8 = decompress_point(r_b8_bytes);
    match r_b8 {
        Result::Err(err) => return Err(err.to_string()),
        Result::Ok(res) => Ok(Signature {
            r_b8: res.clone(),
            s: s,
        }),
    }
}

pub struct PrivateKey {
    key: BigInt,
}

impl PrivateKey {
    pub fn public(&self) -> Result<Point, String> {
        // https://tools.ietf.org/html/rfc8032#section-5.1.5
        let pk = B8.mul_scalar(&self.key)?;
        Ok(pk.clone())
    }

    pub fn sign(&self, msg: BigInt) -> Result<Signature, String> {
        if msg > Q.clone() {
            return Err("msg outside the Finite Field".to_string());
        }
        // https://tools.ietf.org/html/rfc8032#section-5.1.6
        let mut hasher = Blake2b::new();
        let (_, sk_bytes) = self.key.to_bytes_be();
        hasher.input(sk_bytes);
        let mut h = hasher.result(); // h: hash(sk)
                                     // s: h[32:64]
        let (_, msg_bytes) = msg.to_bytes_be();
        let msgFr: Fr = Fr::from_str(&msg.to_string()).unwrap();

        let s = GenericArray::<u8, generic_array::typenum::U32>::from_mut_slice(&mut h[32..64]);
        let r_bytes = utils::concatenate_arrays(s, &msg_bytes);
        let mut r = BigInt::from_bytes_be(Sign::Plus, &r_bytes[..]);
        r = utils::modulus(&r, &SUBORDER);
        let r8: Point = B8.mul_scalar(&r)?;
        let a = &self.public()?;

        let hm_input = vec![r8.x.clone(), r8.y.clone(), a.x.clone(), a.y.clone(), msgFr];
        let hm = poseidon.hash(hm_input)?;

        let mut s = &self.key << 3;
        let hmB = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
        s = hmB * s;
        s = r + s;
        s = s % &SUBORDER.clone();

        Ok(Signature {
            r_b8: r8.clone(),
            s: s,
        })
    }

    pub fn sign_schnorr(&self, m: BigInt) -> Result<(Point, BigInt), String> {
        // random r
        let mut rng = rand6::thread_rng();
        let k = rng.gen_biguint(1024).to_bigint().unwrap();

        // r = k·G
        let r = B8.mul_scalar(&k)?;

        // h = H(x, r, m)
        let pk = &self.public()?;
        let h = schnorr_hash(&pk, m, &r)?;

        // s= k+x·h
        let s = k + &self.key * &h;
        Ok((r, s))
    }
}

pub fn schnorr_hash(pk: &Point, msg: BigInt, c: &Point) -> Result<BigInt, String> {
    if msg > Q.clone() {
        return Err("msg outside the Finite Field".to_string());
    }
    let msgFr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    let hm_input = vec![pk.x.clone(), pk.y.clone(), c.x.clone(), c.y.clone(), msgFr];
    let h = poseidon.hash(hm_input)?;
    println!("h {:?}", h.to_string());
    let hB = BigInt::parse_bytes(to_hex(&h).as_bytes(), 16).unwrap();
    Ok(hB)
}

pub fn verify_schnorr(pk: Point, m: BigInt, r: Point, s: BigInt) -> Result<bool, String> {
    // sG = s·G
    let sg = B8.mul_scalar(&s)?;

    // r + h · x
    let h = schnorr_hash(&pk, m, &r)?;
    let pk_h = pk.mul_scalar(&h)?;
    let right = r.projective().add(&pk_h.projective())?;

    Ok(sg.equals(right.affine()))
}

pub fn new_key() -> PrivateKey {
    // https://tools.ietf.org/html/rfc8032#section-5.1.5
    let mut rng = rand6::thread_rng();
    let sk_raw = rng.gen_biguint(1024).to_bigint().unwrap();

    let mut hasher = Blake2b::new();
    let (_, sk_raw_bytes) = sk_raw.to_bytes_be();
    hasher.input(sk_raw_bytes);
    let mut h = hasher.result();

    h[0] = h[0] & 0xF8;
    h[31] = h[31] & 0x7F;
    h[31] = h[31] | 0x40;

    let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);

    PrivateKey { key: sk }
}

pub fn verify(pk: Point, sig: Signature, msg: BigInt) -> bool {
    if msg > Q.clone() {
        return false;
    }
    let msgFr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    let hm_input = vec![
        sig.r_b8.x.clone(),
        sig.r_b8.y.clone(),
        pk.x.clone(),
        pk.y.clone(),
        msgFr,
    ];
    let hm = match poseidon.hash(hm_input) {
        Result::Err(_) => return false,
        Result::Ok(hm) => hm,
    };
    let l = match B8.mul_scalar(&sig.s) {
        Result::Err(_) => return false,
        Result::Ok(l) => l,
    };
    let hmB = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
    let r = match sig.r_b8.projective().add(
        &pk.mul_scalar(&(8.to_bigint().unwrap() * hmB))
            .unwrap()
            .projective(),
    ) {
        Result::Err(_) => return false,
        Result::Ok(r) => r,
    };
    l.equals(r.affine())
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate rustc_hex;
    use rustc_hex::{FromHex, ToHex};

    #[test]
    fn test_add_same_point() {
        let p: PointProjective = PointProjective {
            x: Fr::from_str(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
            z: Fr::one(),
        };
        let q: PointProjective = PointProjective {
            x: Fr::from_str(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
            z: Fr::one(),
        };
        let res = p.add(&q).unwrap().affine();
        assert_eq!(
            res.x,
            Fr::from_str(
                "6890855772600357754907169075114257697580319025794532037257385534741338397365"
            )
            .unwrap()
        );
        assert_eq!(
            res.y,
            Fr::from_str(
                "4338620300185947561074059802482547481416142213883829469920100239455078257889"
            )
            .unwrap()
        );
    }
    #[test]
    fn test_add_different_points() {
        let p: PointProjective = PointProjective {
            x: Fr::from_str(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
            z: Fr::one(),
        };
        let q: PointProjective = PointProjective {
            x: Fr::from_str(
                "16540640123574156134436876038791482806971768689494387082833631921987005038935",
            )
            .unwrap(),
            y: Fr::from_str(
                "20819045374670962167435360035096875258406992893633759881276124905556507972311",
            )
            .unwrap(),
            z: Fr::one(),
        };
        let res = p.add(&q).unwrap().affine();
        assert_eq!(
            res.x,
            Fr::from_str(
                "7916061937171219682591368294088513039687205273691143098332585753343424131937"
            )
            .unwrap()
        );
        assert_eq!(
            res.y,
            Fr::from_str(
                "14035240266687799601661095864649209771790948434046947201833777492504781204499"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_mul_scalar() {
        let p: Point = Point {
            x: Fr::from_str(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
        };
        let res_m = p.mul_scalar(&3.to_bigint().unwrap()).unwrap();
        let res_a = p.projective().add(&p.projective()).unwrap();
        let res_a = res_a.add(&p.projective()).unwrap().affine();
        assert_eq!(res_m.x, res_a.x);
        assert_eq!(
            res_m.x,
            Fr::from_str(
                "19372461775513343691590086534037741906533799473648040012278229434133483800898"
            )
            .unwrap()
        );
        assert_eq!(
            res_m.y,
            Fr::from_str(
                "9458658722007214007257525444427903161243386465067105737478306991484593958249"
            )
            .unwrap()
        );

        let n = BigInt::parse_bytes(
            b"14035240266687799601661095864649209771790948434046947201833777492504781204499",
            10,
        )
        .unwrap();
        let res2 = p.mul_scalar(&n).unwrap();
        assert_eq!(
            res2.x,
            Fr::from_str(
                "17070357974431721403481313912716834497662307308519659060910483826664480189605"
            )
            .unwrap()
        );
        assert_eq!(
            res2.y,
            Fr::from_str(
                "4014745322800118607127020275658861516666525056516280575712425373174125159339"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_new_key_sign_verify_0() {
        let sk = new_key();
        let pk = sk.public().unwrap();
        let msg = 5.to_bigint().unwrap();
        let sig = sk.sign(msg.clone()).unwrap();
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_new_key_sign_verify_1() {
        let sk = new_key();
        let pk = sk.public().unwrap();
        let msg = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();
        let sig = sk.sign(msg.clone()).unwrap();
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_point_compress_decompress() {
        let p: Point = Point {
            x: Fr::from_str(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
        };
        let p_comp = p.compress();
        assert_eq!(
            p_comp[..].to_hex(),
            "53b81ed5bffe9545b54016234682e7b2f699bd42a5e9eae27ff4051bc698ce85"
        );
        let p2 = decompress_point(p_comp).unwrap();
        assert_eq!(p.x, p2.x);
        assert_eq!(p.y, p2.y);
    }

    #[test]
    fn test_point_decompress0() {
        let y_bytes_raw = "b5328f8791d48f20bec6e481d91c7ada235f1facf22547901c18656b6c3e042f"
            .from_hex()
            .unwrap();
        let mut y_bytes: [u8; 32] = [0; 32];
        y_bytes.copy_from_slice(&y_bytes_raw);
        let p = decompress_point(y_bytes).unwrap();

        let expected_px_raw = "b86cc8d9c97daef0afe1a4753c54fb2d8a530dc74c7eee4e72b3fdf2496d2113"
            .from_hex()
            .unwrap();
        let mut e_px_bytes: [u8; 32] = [0; 32];
        e_px_bytes.copy_from_slice(&expected_px_raw);
        let expected_px: Fr =
            Fr::from_str(&BigInt::from_bytes_le(Sign::Plus, &e_px_bytes).to_string()).unwrap();
        assert_eq!(&p.x, &expected_px);
    }

    #[test]
    fn test_point_decompress1() {
        let y_bytes_raw = "70552d3ff548e09266ded29b33ce75139672b062b02aa66bb0d9247ffecf1d0b"
            .from_hex()
            .unwrap();
        let mut y_bytes: [u8; 32] = [0; 32];
        y_bytes.copy_from_slice(&y_bytes_raw);
        let p = decompress_point(y_bytes).unwrap();

        let expected_px_raw = "30f1635ba7d56f9cb32c3ffbe6dca508a68c7f43936af11a23c785ce98cb3404"
            .from_hex()
            .unwrap();
        let mut e_px_bytes: [u8; 32] = [0; 32];
        e_px_bytes.copy_from_slice(&expected_px_raw);
        let expected_px: Fr =
            Fr::from_str(&BigInt::from_bytes_le(Sign::Plus, &e_px_bytes).to_string()).unwrap();
        assert_eq!(&p.x, &expected_px);
    }

    #[test]
    fn test_point_decompress_loop() {
        for _ in 0..5 {
            let random_bytes = rand6::thread_rng().gen::<[u8; 32]>();
            let sk_raw: BigInt = BigInt::from_bytes_le(Sign::Plus, &random_bytes[..]);
            let mut hasher = Blake2b::new();
            let (_, sk_raw_bytes) = sk_raw.to_bytes_be();
            hasher.input(sk_raw_bytes);
            let mut h = hasher.result();

            h[0] = h[0] & 0xF8;
            h[31] = h[31] & 0x7F;
            h[31] = h[31] | 0x40;

            let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);
            let point = B8.mul_scalar(&sk).unwrap();
            let cmp_point = point.compress();
            let dcmp_point = decompress_point(cmp_point).unwrap();

            assert_eq!(&point.x, &dcmp_point.x);
            assert_eq!(&point.y, &dcmp_point.y);
        }
    }

    #[test]
    fn test_signature_compress_decompress() {
        let sk = new_key();
        let pk = sk.public().unwrap();

        for i in 0..5 {
            let msg_raw = "123456".to_owned() + &i.to_string();
            let msg = BigInt::parse_bytes(msg_raw.as_bytes(), 10).unwrap();
            let sig = sk.sign(msg.clone()).unwrap();

            let compressed_sig = sig.compress();
            let decompressed_sig = decompress_signature(&compressed_sig).unwrap();
            assert_eq!(&sig.r_b8.x, &decompressed_sig.r_b8.x);
            assert_eq!(&sig.r_b8.y, &decompressed_sig.r_b8.y);
            assert_eq!(&sig.s, &decompressed_sig.s);

            let v = verify(pk.clone(), decompressed_sig, msg);
            assert_eq!(v, true);
        }
    }

    #[test]
    fn test_schnorr_signature() {
        let sk = new_key();
        let pk = sk.public().unwrap();

        let msg = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();
        let (s, e) = sk.sign_schnorr(msg.clone()).unwrap();
        println!("s {:?}", s.x.to_string());
        println!("s {:?}", s.y.to_string());
        println!("e {:?}", e.to_string());
        let verification = verify_schnorr(pk, msg, s, e).unwrap();
        assert_eq!(true, verification);
    }
}
