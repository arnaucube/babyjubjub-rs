// BabyJubJub elliptic curve implementation in Rust.
// For LICENSE check https://github.com/arnaucube/babyjubjub-rs

use ff::*;
use std::fmt::Error;
use num::Num;
use std::fmt;
use serde::{Serialize, ser::SerializeStruct, de::Visitor, de::MapAccess, Deserialize, Deserializer};
use poseidon_rs::Poseidon;
pub type Fr = poseidon_rs::Fr; // alias

extern crate rand_new;
extern crate rand;
extern crate ff;

// Create a new primefield for the subgroup defined by the base point, order Fl:
#[derive(PrimeField)]
#[PrimeFieldModulus = "2736030358979909402780800718157159386076813972158567259200215660948447373041"]
#[PrimeFieldGenerator = "7"] // TODO: double check this is a valid generator!!!
pub struct Fl(FpRepr);

use arrayref::array_ref;

#[cfg(not(feature = "aarch64"))]
use blake_hash::Digest; // compatible version with Blake used at circomlib

#[cfg(feature = "aarch64")]
extern crate blake; // compatible version with Blake used at circomlib

use std::{cmp::min, str::FromStr};

use num_bigint::{BigInt, RandBigInt, Sign, ToBigInt};
use num_traits::One;

use generic_array::GenericArray;

pub mod utils;

use lazy_static::lazy_static;

lazy_static! {
    static ref D: Fr = Fr::from_str("168696").unwrap();
    static ref D_BIG: BigInt = BigInt::parse_bytes(b"168696", 10).unwrap();
    static ref A: Fr = Fr::from_str("168700").unwrap();
    static ref A_BIG: BigInt = BigInt::parse_bytes(b"168700", 10).unwrap();
    pub static ref Q: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",10
    )
        .unwrap();

    pub static ref G: Point = Point {
        x: Fr::from_str(
                "995203441582195749578291179787384436505546430278305826713579947235728471134",
            )
            .unwrap(),
            y: Fr::from_str(
                "5472060717959818805561601436314318772137091100104008585924551046643952123905",
            )
            .unwrap(),
    };

    pub static ref B8: Point = Point {
        x: Fr::from_str(
               "5299619240641551281634865583518297030282874472190772894086521144482721001553",
           )
            .unwrap(),
            y: Fr::from_str(
                "16950150798460657717958625567821834550301663161624707787222815936182638968203",
            )
            .unwrap(),
    };

    pub static ref O: Point = Point {
        x: Fr::zero(),
        y: Fr::one()
    };

    // Order expressed as a bigint (it is larger than the r modulus of the Fr struct)
    pub static ref ORDER: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088614511777268538073601725287587578984328",
        10,
    ).unwrap();

    // SUBORDER = ORDER >> 3
    pub static ref SUBORDER: BigInt = &BigInt::parse_bytes(
        b"21888242871839275222246405745257275088614511777268538073601725287587578984328",
        10,
    ).unwrap()
    >> 3;
    pub static ref POSEIDON: poseidon_rs::Poseidon = Poseidon::new();
    
    // MAX_MSG is maximum message length that can be encoded into a point
    pub static ref MAX_MSG: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",10
    )
        .unwrap()
        >> 10;

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

        let zinv = self.z.inverse().unwrap();
        let mut x = self.x;
        x.mul_assign(&zinv);
        let mut y = self.y;
        y.mul_assign(&zinv);

        Point { x, y }
    }

    #[allow(clippy::many_single_char_names)]
    pub fn add(&self, q: &PointProjective) -> PointProjective {
        // add-2008-bbjlp https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#addition-add-2008-bbjlp
        let mut a = self.z;
        a.mul_assign(&q.z);
        let mut b = a;
        b.square();
        let mut c = self.x;
        c.mul_assign(&q.x);
        let mut d = self.y;
        d.mul_assign(&q.y);
        let mut e = *D;
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
        let mut ac = *A;
        ac.mul_assign(&c);
        let mut dac = d;
        dac.sub_assign(&ac);
        let mut y3 = a;
        y3.mul_assign(&g);
        y3.mul_assign(&dac);
        let mut z3 = f;
        z3.mul_assign(&g);

        PointProjective {
            x: x3,
            y: y3,
            z: z3,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Point {
    pub x: Fr,
    pub y: Fr,
}

pub trait ToDecimalString {
    fn to_dec_string(&self) -> String;
}
impl ToDecimalString for Fr {
    fn to_dec_string(&self) -> String {
        let s = self.to_string();
        let hex_str = s[5..s.len()-1].to_string();
        BigInt::from_str_radix(&hex_str, 16).unwrap().to_string()
    }
}

impl ToDecimalString for Fl {
    fn to_dec_string(&self) -> String {
        let s = self.to_string();
        let hex_str = s[5..s.len()-1].to_string();
        BigInt::from_str_radix(&hex_str, 16).unwrap().to_string()
    }
}
pub trait FrBigIntConversion<T> {
    fn from_bigint(bi: &BigInt) -> T;
    fn to_bigint(&self) -> BigInt;
}

impl FrBigIntConversion<Fr> for Fr {
    // Note: this could probably be more efficient by converting bigint to raw repr to Fr
    fn from_bigint(bi: &BigInt) -> Fr {
        Fr::from_str(&bi.to_string()).unwrap()
    }
    fn to_bigint(&self) -> BigInt {
        BigInt::from_str(&self.to_dec_string()).unwrap()
    }
}

impl FrBigIntConversion<Fl> for Fl {
    // Note: this could probably be more efficient by converting bigint to raw repr to Fr
    fn from_bigint(bi: &BigInt) -> Fl {
        Fl::from_str(&bi.to_string()).unwrap()
    }
    fn to_bigint(&self) -> BigInt {
        BigInt::from_str(&self.to_dec_string()).unwrap()
    }
}

impl Serialize for Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: serde::Serializer {
                let mut state = serializer.serialize_struct("Point", 2)?;
                state.serialize_field("x", &self.x.to_dec_string())?;
                state.serialize_field("y", &self.y.to_dec_string())?;
                state.end()
    }
}


impl<'de> Deserialize<'de> for Point {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>
    {
        deserializer.deserialize_map(PointVisitor)
    }
}
// For deserialization:
struct PointVisitor;

impl<'de> Visitor<'de> for PointVisitor {
    type Value = Point;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(formatter, "a map with keys 'x' and 'y'")
    }

    fn visit_map<M>(self, mut map: M) -> Result<Self::Value, M::Error>
    where
        M: MapAccess<'de>
    {
        let mut x = None;
        let mut y = None;

        while let Some(k) = map.next_key::<String>()? {
            if k == "x" {
                let x_: String = map.next_value()?;
                x = Fr::from_str(&x_);
            }
            else if k == "y" {
                let y_: String = map.next_value()?;
                y = Fr::from_str(&y_);
            }
            else {
                return Err(serde::de::Error::custom(&format!("Invalid key: {}", k)));
            }
        }

        if x.is_none() || y.is_none() {
            return Err(serde::de::Error::custom("Missing first or second"));
        }

        Ok(Point { x: x.unwrap(), y: y.unwrap() })
    }
}

// For addition shorthand
impl Point {
    pub fn projective(&self) -> PointProjective {
        PointProjective {
            x: self.x,
            y: self.y,
            z: Fr::one(),
        }
    }

    pub fn add(&self, another_point: &Point) -> Point {
        self.projective().add(&another_point.projective()).affine()
    }
    pub fn neg(&self) -> Point {
        let mut x_inverse = Fr::zero();
        x_inverse.sub_assign(&self.x);
        Point {
            x: x_inverse,
            y: self.y
        }
    }

    pub fn mul_scalar(&self, n: &BigInt) -> Point {
        let mut r: PointProjective = PointProjective {
            x: Fr::zero(),
            y: Fr::one(),
            z: Fr::one(),
        };
        let mut exp: PointProjective = self.projective();
        let (_, b) = n.to_bytes_le();
        for i in 0..n.bits() {
            if test_bit(&b, i.try_into().unwrap()) {
                r = r.add(&exp);
            }
            exp = exp.add(&exp);
        }
        r.affine()
    }

    pub fn compress(&self) -> [u8; 32] {
        let p = &self;
        let mut r: [u8; 32] = [0; 32];
        let x_big = BigInt::parse_bytes(to_hex(&p.x).as_bytes(), 16).unwrap();
        let y_big = BigInt::parse_bytes(to_hex(&p.y).as_bytes(), 16).unwrap();
        let (_, y_bytes) = y_big.to_bytes_le();
        let len = min(y_bytes.len(), r.len());
        r[..len].copy_from_slice(&y_bytes[..len]);
        if x_big > (&Q.clone() >> 1) {
            r[31] |= 0x80;
        }
        r
    }

    pub fn equals(&self, p: Point) -> bool {
        if self.x == p.x && self.y == p.y {
            return true;
        }
        false
    }

    /// Koblitz decoding method, adapted for this curve: 
        /// message m must be < r/(2^10) //using 2^10=1024 instead of Koblitz' parameter of 1000, so arithmetic can happen easily mod 2. This shouldn't have any (negative) impact: https://crypto.stackexchange.com/questions/103132/encoding-a-message-as-points-on-elliptic-curve
        /// Try finding a point with y value m*1024+0, m*1024+1, .... m*1024+5617 (5617 are last four digits of prime r)
        /// There is an approximately 1/(2^1024) chance no point will be encodable,
        /// since each y value has probability of about 1/2 of being on the curve
    pub fn from_msg_vartime(msg: &BigInt) -> Option<Point> {
        #[allow(non_snake_case)]
        let ACC_UNDER = 1024; // Last four digits of prime r. MAX_MSG * 1024 + ACC_UNDER = r
        assert!(msg <= &MAX_MSG);
        let mut acc: u16 = 0;
        // Start with message * 10000 as x coordinate
        let mut x: Fr = Fr::from_str(&msg.to_str_radix(10)).unwrap();
        let mut y: Fr; 
        x.mul_assign(&Fr::from_str(&ACC_UNDER.to_string()).unwrap());
        let one = Fr::one();

        while acc < ACC_UNDER {
            // If x is on curve, calculate what y^2 should be, by (ax^2 - 1) / (dx^2 - 1)
            let mut x2 = x.clone();
            x2.mul_assign(&x);

            // Numerator will be (ax^2 - 1) and denominator will be (dx^2 - 1)
            let mut numerator = x2;
            let mut denominator = x2.clone();

            numerator.mul_assign(&A);
            denominator.mul_assign(&D);

            numerator.sub_assign(&one);
            denominator.sub_assign(&one);

            // If the point is on the curve, numerator/denominator will be y^2. Check whether numerator/denominator is a quadratic residue:
            numerator.mul_assign(&denominator.inverse().unwrap()); // Note: this is no longer a numerator since it was divided in this step
            
            if let LegendreSymbol::QuadraticResidue = numerator.legendre() {
                y = numerator.sqrt().unwrap();
                let pt = Point {x, y};
                if pt.in_subgroup() {
                    return Some(Point {x, y})
                }
            } 
            acc += 1;
            x.add_assign(&one);
        }
        return None
    }

    /// Converts a point to a message by dividing its x by 1024 (a.k.a. right-shifting by 10)
    pub fn to_msg(&self) -> Fr {
        let mut msg = self.x.clone().into_repr();
        msg.shr(10);
        Fr::from_repr(msg).unwrap()

    }
    /// Checks that a point is on the BabyJubJub curve. Does not check the point is in the correct subgroup.
    pub fn on_curve(&self) -> bool {
        let mut x2 = self.x.clone();
        let mut y2 = self.y.clone();
        x2.mul_assign(&self.x);
        y2.mul_assign(&self.y);
        // compute left hand side ax^2+y^2
        let mut lhs = x2.clone();
        lhs.mul_assign(&A);
        lhs.add_assign(&y2);
        // compute right hand side: x^2*y^2*d+1
        let mut rhs = x2.clone();
        rhs.mul_assign(&y2);
        rhs.mul_assign(&D);
        rhs.add_assign(&Fr::one());

        lhs.eq(&rhs)
    }

    /// Checks that a point's order is equal to the subgroup order. Does not check the point is on the curve.
    pub fn in_subgroup(&self) -> bool {
        let should_be_zero = self.mul_scalar(&SUBORDER);
        should_be_zero.equals({
            Point { x: Fr::zero(), y: Fr::one() }
        })
    }

    pub fn from_xy_strings(x: String, y: String) -> Point {
        Point {
            x: Fr::from_str(&x).unwrap(),
            y: Fr::from_str(&y).unwrap()
        }
    }
}

pub fn test_bit(b: &[u8], i: usize) -> bool {
    b[i / 8] & (1 << (i % 8)) != 0
}

pub fn decompress_point(bb: [u8; 32]) -> Result<Point, String> {
    // https://tools.ietf.org/html/rfc8032#section-5.2.3
    let mut sign: bool = false;
    let mut b = bb;
    if b[31] & 0x80 != 0x00 {
        sign = true;
        b[31] &= 0x7F;
    }
    let y: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[..]);
    if y >= Q.clone() {
        return Err("y outside the Finite Field over R".to_string());
    }
    let one: BigInt = One::one();

    // x^2 = (1 - y^2) / (a - d * y^2) (mod p)
    let den = utils::modinv(
        &utils::modulus(
            &(&A_BIG.clone() - utils::modulus(&(&D_BIG.clone() * (&y * &y)), &Q)),
            &Q,
        ),
        &Q,
    )?;
    let mut x: BigInt = utils::modulus(&((one - utils::modulus(&(&y * &y), &Q)) * den), &Q);
    x = utils::modsqrt(&x, &Q)?;

    if sign && (x <= (&Q.clone() >> 1)) || (!sign && (x > (&Q.clone() >> 1))) {
        x *= -(1.to_bigint().unwrap());
    }
    x = utils::modulus(&x, &Q);
    let x_fr: Fr = Fr::from_str(&x.to_string()).unwrap();
    let y_fr: Fr = Fr::from_str(&y.to_string()).unwrap();
    Ok(Point { x: x_fr, y: y_fr })
}

#[cfg(not(feature = "aarch64"))]
fn blh(b: &[u8]) -> Vec<u8> {
    let hash = blake_hash::Blake512::digest(b);
    hash.to_vec()
}

#[cfg(feature = "aarch64")]
fn blh(b: &[u8]) -> Vec<u8> {
    let mut hash = [0; 64];
    blake::hash(512, b, &mut hash).unwrap();
    hash.to_vec()
}

#[derive(Debug, Clone, Serialize)]
pub struct Signature {
    pub r_b8: Point,
    pub s: BigInt,
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
        Result::Err(err) => Err(err),
        Result::Ok(res) => Ok(Signature { r_b8: res, s }),
    }
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ElGamalEncryption {
    pub c1: Point,
    pub c2: Point
}

pub struct PrivateKey {
    pub key: [u8; 32],
}

impl PrivateKey {
    pub fn import(b: Vec<u8>) -> Result<PrivateKey, String> {
        if b.len() != 32 {
            return Err(String::from("imported key can not be bigger than 32 bytes"));
        }
        let mut sk: [u8; 32] = [0; 32];
        sk.copy_from_slice(&b[..32]);
        Ok(PrivateKey { key: sk })
    }

    pub fn scalar_key(&self) -> BigInt {
        // not-compatible with circomlib implementation, but using Blake2b
        // let mut hasher = Blake2b::new();
        // hasher.update(sk_raw_bytes);
        // let mut h = hasher.finalize();

        // compatible with circomlib implementation
        let hash: Vec<u8> = blh(&self.key);
        let mut h: Vec<u8> = hash[..32].to_vec();

        // prune buffer following RFC 8032
        // https://tools.ietf.org/html/rfc8032#page-13
        h[0] &= 0xF8;
        h[31] &= 0x7F;
        h[31] |= 0x40;

        let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);
        sk >> 3
    }

    pub fn public(&self) -> Point {
        B8.mul_scalar(&self.scalar_key())
    }

    pub fn sign(&self, msg: BigInt) -> Result<Signature, String> {
        if msg > Q.clone() {
            return Err("msg outside the Finite Field".to_string());
        }
        // let (_, sk_bytes) = self.key.to_bytes_le();
        // let mut hasher = Blake2b::new();
        // hasher.update(sk_bytes);
        // let mut h = hasher.finalize(); // h: hash(sk), s: h[32:64]
        let mut h: Vec<u8> = blh(&self.key);

        let (_, msg_bytes) = msg.to_bytes_le();
        let mut msg32: [u8; 32] = [0; 32];
        msg32[..msg_bytes.len()].copy_from_slice(&msg_bytes[..]);
        let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();

        // https://tools.ietf.org/html/rfc8032#section-5.1.6
        let s = GenericArray::<u8, generic_array::typenum::U32>::from_mut_slice(&mut h[32..64]);
        let r_bytes = utils::concatenate_arrays(s, &msg32);
        let r_hashed: Vec<u8> = blh(&r_bytes);
        let mut r = BigInt::from_bytes_le(Sign::Plus, &r_hashed[..]);
        r = utils::modulus(&r, &SUBORDER);
        let r_b8: Point = B8.mul_scalar(&r);
        let a = &self.public();

        let hm_input = vec![r_b8.x, r_b8.y, a.x, a.y, msg_fr];
        let hm = POSEIDON.hash(hm_input)?;

        let mut s = &self.scalar_key() << 3;
        let hm_b = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
        s = hm_b * s;
        s = r + s;
        s %= &SUBORDER.clone();

        Ok(Signature { r_b8, s: s })
    }

    #[allow(clippy::many_single_char_names)]
    pub fn sign_schnorr(&self, m: BigInt) -> Result<(Point, BigInt), String> {
        // random r
        let mut rng = rand_new::thread_rng();
        let k = rng.gen_biguint(1024).to_bigint().unwrap();

        // r = k·G
        let r = B8.mul_scalar(&k);

        // h = H(x, r, m)
        let pk = self.public();
        let h = schnorr_hash(&pk, m, &r)?;

        // s= k+x·h
        let sk_scalar = self.scalar_key();
        let s = k + &sk_scalar * &h;
        Ok((r, s))
    }

    pub fn decrypt_elgamal(&self, encrypted_point: ElGamalEncryption) -> Point {
        // Make sure inputs aren't bad (i imagine this check could be skipped for performance reasons, but it seems a sanity check here would be helpful)
        assert!(encrypted_point.c1.on_curve(), "Error: C1 is not on the curve!");
        assert!(encrypted_point.c1.in_subgroup(), "Error: C1 is not in the subgroup!");
        assert!(encrypted_point.c2.on_curve(), "Error: C2 is not on the curve!");
        assert!(encrypted_point.c2.in_subgroup(), "Error: C2 is not in the subgroup!");
        
        let shared_secret = encrypted_point.c1.mul_scalar(&self.scalar_key());
        // Subtract the shared secret
        encrypted_point.c2.add(&shared_secret.neg())
    }

}

// encrypts msg to public key to_pubkey, using some random scalar nonce
pub fn encrypt_elgamal(to_pubkey: &Point, nonce: &BigInt, msg: &Point) -> ElGamalEncryption {
    let shared_secret = to_pubkey.mul_scalar(&nonce);
    let public_nonce = B8.mul_scalar(&nonce);
    // let msg_point = point_for_msg(msg);
    let msg_plus_secret = msg.add(&shared_secret);

    ElGamalEncryption {
        c1: public_nonce,
        c2: msg_plus_secret
    }
}

pub fn schnorr_hash(pk: &Point, msg: BigInt, c: &Point) -> Result<BigInt, String> {
    if msg > Q.clone() {
        return Err("msg outside the Finite Field".to_string());
    }
    let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    let hm_input = vec![pk.x, pk.y, c.x, c.y, msg_fr];
    let h = POSEIDON.hash(hm_input)?;
    let h_b = BigInt::parse_bytes(to_hex(&h).as_bytes(), 16).unwrap();
    Ok(h_b)
}

pub fn verify_schnorr(pk: Point, m: BigInt, r: Point, s: BigInt) -> Result<bool, String> {
    // sG = s·G
    let sg = B8.mul_scalar(&s);

    // r + h · x
    let h = schnorr_hash(&pk, m, &r)?;
    let pk_h = pk.mul_scalar(&h);
    let right = r.add(&pk_h);

    Ok(sg.equals(right))
}

pub fn new_key() -> PrivateKey {
    // https://tools.ietf.org/html/rfc8032#section-5.1.5
    let mut rng = rand_new::thread_rng();
    let sk_raw = rng.gen_biguint(1024).to_bigint().unwrap();
    let (_, sk_raw_bytes) = sk_raw.to_bytes_be();
    PrivateKey::import(sk_raw_bytes[..32].to_vec()).unwrap()
}

pub fn verify(pk: Point, sig: Signature, msg: BigInt) -> bool {
    if msg > Q.clone() {
        return false;
    }
    let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    let hm_input = vec![sig.r_b8.x, sig.r_b8.y, pk.x, pk.y, msg_fr];
    let hm = match POSEIDON.hash(hm_input) {
        Result::Err(_) => return false,
        Result::Ok(hm) => hm,
    };
    let l = B8.mul_scalar(&sig.s);
    let hm_b = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
    let r = sig
        .r_b8
        .add(&pk.mul_scalar(&(8.to_bigint().unwrap() * hm_b)));
    l.equals(r)
}


pub struct DLEQProof {
    pub A: Point,
    pub B: Point,
    pub xA: Point,
    pub xB: Point,
    pub challenge: Fl,
    pub response: Fl,
}
impl DLEQProof {
     // Prove x*A = variable xA and x*B = variable xB
    pub fn new(x: Fl, point_A: Point, point_B: Point) -> Result<DLEQProof, Error> {
        let x_bigint = x.to_bigint();

        let modulus = SUBORDER.clone();
        // TODO: better error handling (not assert), make it more efficient too:
        assert!(x_bigint < modulus);

        let k_bigint = rand_new::thread_rng().gen_biguint(512).to_bigint().unwrap() % &modulus;
        let k = Fl::from_bigint(&k_bigint);

        let xA = point_A.mul_scalar(&x_bigint);
        let xB = point_B.mul_scalar(&x_bigint);

        let kA = point_A.mul_scalar(&k_bigint);
        let kB = point_B.mul_scalar(&k_bigint);
        
        let challenge = DLEQProof::get_challenge(&point_A, &point_B, &xA, &xB, &kA, &kB);
        let mut challenge_x = challenge.clone(); challenge_x.mul_assign(&x);
        // response = k - challenge * x;
        let mut response = k.clone(); response.sub_assign(&challenge_x);

        Ok(DLEQProof {
            A: point_A,
            B: point_B,
            xA: xA,
            xB: xB,
            challenge: challenge,
            response: response,
        })
    }

    pub fn verify(&self) -> bool {
        // This should equal kA
        let kA_ = 
            self.A.mul_scalar(&self.response.to_bigint()).add(
                &self.xA.mul_scalar(&self.challenge.to_bigint())
            );
        
        // This should equal kB
        let kB_ = 
            self.B.mul_scalar(&self.response.to_bigint()).add(
                &self.xB.mul_scalar(&self.challenge.to_bigint())
        );

        let challenge = DLEQProof::get_challenge(&self.A, &self.B, &self.xA, &self.xB, &kA_, &kB_);

        return challenge == self.challenge;

    }
    // Generates randomness for DLEQ Fiat-Shamir transform
    fn get_challenge(point_A: &Point, point_B: &Point, xA: &Point, xB: &Point, kA: &Point, kB: &Point) -> Fl {
        // This could probably be faster if we neglect either the x or y coordinate, but rn doing both to be safe until i study this more
        let inputs: Vec<Fr> = vec![
                point_A.x,point_A.y, 
                xA.x,xA.y,
                point_B.x,point_B.y,
                xB.x,xB.y,
                kA.x,kA.y,
                kB.x,kB.y
        ];
        let input: Vec<u8> = inputs.iter().map(|f| {
            let as_hex = to_hex(f);
            hex::decode(as_hex).unwrap()
        })
        .flatten()
        .collect();
        
        let challenge_hash = blh(&input);
        Fl::from_bigint(&BigInt::from_bytes_be(Sign::Plus, &challenge_hash.as_slice()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::hex;
    use rand_new::Rng;
    use num_traits::FromPrimitive;

    #[test]
    fn test_dleq() {
        let a = BigInt::parse_bytes(b"123456789", 10).unwrap();
        let point_A = B8.mul_scalar(&a);
        let b = BigInt::parse_bytes(b"55555555555555512345678944444444444", 10).unwrap();
        let point_B = B8.mul_scalar(&b);

        let mut rng = rand_new::thread_rng();
        let x = Fl::from_bigint(&rng.gen_biguint(512).to_bigint().unwrap());
        let proof = DLEQProof::new(x, point_A, point_B).unwrap();
        assert!(proof.verify());
    }
    #[test]
    fn test_on_curve() {
        assert_eq!(B8.on_curve(), true);
        assert_eq!(B8.mul_scalar(&12345.to_bigint().unwrap()).on_curve(), true);

        let some_point = Point { x: Fr::from_str("1234").unwrap(), y: Fr::from_str("5678").unwrap() };
        assert_eq!(some_point.on_curve(), false);
    }

    #[test]
    fn test_in_subgroup() {
        assert_eq!(B8.in_subgroup(), true);
        assert_eq!(B8.mul_scalar(&12345.to_bigint().unwrap()).in_subgroup(), true);

        assert_eq!(G.in_subgroup(), false);
        assert_eq!(G.mul_scalar(&5.to_bigint().unwrap()).in_subgroup(), false);
        assert_eq!(G.mul_scalar(&7.to_bigint().unwrap()).in_subgroup(), false);
        
        assert_eq!(G.mul_scalar(&8.to_bigint().unwrap()).in_subgroup(), true);
        assert_eq!(G.mul_scalar(&16.to_bigint().unwrap()).in_subgroup(), true);
        assert_eq!(G.mul_scalar(&8000.to_bigint().unwrap()).in_subgroup(), true);

        


        
    }
    #[test]
    fn test_from_msg_vartime() {
        let msg = 123456789.to_bigint().unwrap();
        assert!(Point::from_msg_vartime(&msg).unwrap().on_curve());

        // Try with some more random numbers -- it's extremely unlikely to get lucky will with valid points 20 times in a row if it's not always producing valid points
        for n in 0..20 {
            let m = rand_new::thread_rng().gen_bigint_range(&0.to_bigint().unwrap() , &MAX_MSG);
            assert!(Point::from_msg_vartime(&m).unwrap().on_curve());
            assert!(Point::from_msg_vartime(&m).unwrap().in_subgroup());
        }

    }

    #[test]
    fn test_from_to_msg() {        
        // Convert from msg to point back to msg and make sure it works:
        let msg = 123456789.to_bigint().unwrap();
        assert!(
            Point::from_msg_vartime(&msg).unwrap().to_msg()
            .eq(
            &Fr::from_str(&msg.to_string()).unwrap()
            )
        );

        // Try with some more random numbers -- it's extremely unlikely to get lucky will with valid points 20 times in a row if it's not always producing valid points
        for n in 0..20 {
            let m = rand_new::thread_rng().gen_bigint_range(&0.to_bigint().unwrap() , &MAX_MSG);
                                                                                                                                                                                       assert!(
                Point::from_msg_vartime(&m).unwrap().to_msg()
                .eq(
                &Fr::from_str(&m.to_string()).unwrap()
                )
            );
        }
    }
    #[test]
    fn test_neg() {
        let some_point = B8.mul_scalar(&BigInt::from_u8(0x69).unwrap());
        let another_point = B8.mul_scalar(&BigInt::from_u8(0x99).unwrap());
        let mut some_point_x_inverse = Fr::zero();
        some_point_x_inverse.sub_assign(&some_point.x);
        // assert_eq!(some_point_x_inverse, some_point.x.inverse().unwrap());
        assert!(some_point.equals(
            some_point.add(&another_point).add(
            &another_point.neg())
        ));

    }
    #[test]
    fn test_elgamal() {
        let some_point = B8.mul_scalar(&BigInt::from_u8(0x69).unwrap());
        let some_privkey = new_key();
        let some_point_encrypted = encrypt_elgamal(
            &some_privkey.public(), 
            &BigInt::parse_bytes(b"ABCDEF123456789", 16).unwrap(), 
            &some_point
        );
        let some_point_encrypted_decrypted = some_privkey.decrypt_elgamal(some_point_encrypted);

        assert_eq!(some_point.x, some_point_encrypted_decrypted.x);
        assert_eq!(some_point.y, some_point_encrypted_decrypted.y);
    }
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
        let res = p.add(&q).affine();
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
        let res = p.add(&q).affine();
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
        let res_m = p.mul_scalar(&3.to_bigint().unwrap());
        let res_a = p.add(&p);
        let res_a = res_a.add(&p);
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
        let res2 = p.mul_scalar(&n);
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
        let pk = sk.public();
        let msg = 5.to_bigint().unwrap();
        let sig = sk.sign(msg.clone()).unwrap();
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_new_key_sign_verify_1() {
        let sk = new_key();
        let pk = sk.public();
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
            hex::encode(p_comp),
            "53b81ed5bffe9545b54016234682e7b2f699bd42a5e9eae27ff4051bc698ce85"
        );
        let p2 = decompress_point(p_comp).unwrap();
        assert_eq!(p.x, p2.x);
        assert_eq!(p.y, p2.y);
    }

    #[test]
    fn test_point_decompress0() {
        let y_bytes_raw =
            hex::decode("b5328f8791d48f20bec6e481d91c7ada235f1facf22547901c18656b6c3e042f")
                .unwrap();
        let mut y_bytes: [u8; 32] = [0; 32];
        y_bytes.copy_from_slice(&y_bytes_raw);
        let p = decompress_point(y_bytes).unwrap();

        let expected_px_raw =
            hex::decode("b86cc8d9c97daef0afe1a4753c54fb2d8a530dc74c7eee4e72b3fdf2496d2113")
                .unwrap();
        let mut e_px_bytes: [u8; 32] = [0; 32];
        e_px_bytes.copy_from_slice(&expected_px_raw);
        let expected_px: Fr =
            Fr::from_str(&BigInt::from_bytes_le(Sign::Plus, &e_px_bytes).to_string()).unwrap();
        assert_eq!(&p.x, &expected_px);
    }

    #[test]
    fn test_point_decompress1() {
        let y_bytes_raw =
            hex::decode("70552d3ff548e09266ded29b33ce75139672b062b02aa66bb0d9247ffecf1d0b")
                .unwrap();
        let mut y_bytes: [u8; 32] = [0; 32];
        y_bytes.copy_from_slice(&y_bytes_raw);
        let p = decompress_point(y_bytes).unwrap();

        let expected_px_raw =
            hex::decode("30f1635ba7d56f9cb32c3ffbe6dca508a68c7f43936af11a23c785ce98cb3404")
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
            let random_bytes = rand_new::thread_rng().gen::<[u8; 32]>();
            let sk_raw: BigInt = BigInt::from_bytes_le(Sign::Plus, &random_bytes[..]);
            let (_, sk_raw_bytes) = sk_raw.to_bytes_be();
            let mut h: Vec<u8> = blh(&sk_raw_bytes);

            h[0] = h[0] & 0xF8;
            h[31] = h[31] & 0x7F;
            h[31] = h[31] | 0x40;

            let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);
            let point = B8.mul_scalar(&sk);
            let cmp_point = point.compress();
            let dcmp_point = decompress_point(cmp_point).unwrap();

            assert_eq!(&point.x, &dcmp_point.x);
            assert_eq!(&point.y, &dcmp_point.y);
        }
    }

    #[test]
    fn test_signature_compress_decompress() {
        let sk = new_key();
        let pk = sk.public();

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
        let pk = sk.public();

        let msg = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();
        let (s, e) = sk.sign_schnorr(msg.clone()).unwrap();
        let verification = verify_schnorr(pk, msg, s, e).unwrap();
        assert_eq!(true, verification);
    }

    #[test]
    fn test_circomlib_testvector() {
        let sk_raw_bytes =
            hex::decode("0001020304050607080900010203040506070809000102030405060708090001")
                .unwrap();

        // test blake compatible with circomlib implementation
        let h: Vec<u8> = blh(&sk_raw_bytes);
        assert_eq!(hex::encode(h), "c992db23d6290c70ffcc02f7abeb00b9d00fa8b43e55d7949c28ba6be7545d3253882a61bd004a236ef1cdba01b27ba0aedfb08eefdbfb7c19657c880b43ddf1");

        // test private key
        let sk = PrivateKey::import(
            hex::decode("0001020304050607080900010203040506070809000102030405060708090001")
                .unwrap(),
        )
        .unwrap();
        assert_eq!(
            sk.scalar_key().to_string(),
            "6466070937662820620902051049739362987537906109895538826186780010858059362905"
        );

        // test public key
        let pk = sk.public();
        assert_eq!(
            pk.x.to_string(),
            "Fr(0x1d5ac1f31407018b7d413a4f52c8f74463b30e6ac2238220ad8b254de4eaa3a2)"
        );
        assert_eq!(
            pk.y.to_string(),
            "Fr(0x1e1de8a908826c3f9ac2e0ceee929ecd0caf3b99b3ef24523aaab796a6f733c4)"
        );

        // test signature & verification
        let msg = BigInt::from_bytes_le(Sign::Plus, &hex::decode("00010203040506070809").unwrap());
        let sig = sk.sign(msg.clone()).unwrap();
        assert_eq!(
            sig.r_b8.x.to_string(),
            "Fr(0x192b4e51adf302c8139d356d0e08e2404b5ace440ef41fc78f5c4f2428df0765)"
        );
        assert_eq!(
            sig.r_b8.y.to_string(),
            "Fr(0x2202bebcf57b820863e0acc88970b6ca7d987a0d513c2ddeb42e3f5d31b4eddf)"
        );
        assert_eq!(
            sig.s.to_string(),
            "1672775540645840396591609181675628451599263765380031905495115170613215233181"
        );
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }
}
