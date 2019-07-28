extern crate num;
extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::{One, Zero};

mod utils;

#[derive(Clone, Debug)]
pub struct Point {
    pub x: BigInt,
    pub y: BigInt,
}

pub struct Babyjubjub {
    d: BigInt,
    a: BigInt,
    q: BigInt,
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
        Babyjubjub { d: d, a: a, q: q }
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        // x = (x1*y2+y1*x2)/(c*(1+d*x1*x2*y1*y2))
        // y = (y1*y2-x1*x2)/(c*(1-d*x1*x2*y1*y2))

        // x = (x1 * y2 + y1 * x2) / (1 + d * x1 * y1 * y2)
        let one: BigInt = One::one();
        let x_num: BigInt = &p.x * &q.y + &p.y * &q.x;
        let x_den: BigInt = &one + &self.d * &p.x * &q.x * &p.y * &q.y;
        let x_den_inv = utils::mod_inverse0(&x_den, &self.q);
        // let x_den_inv = utils::mod_inverse1(x_den, self.q.clone());
        // let x_den_inv = utils::mod_inverse2(x_den, self.q.clone());
        let x: BigInt = utils::modulus(&(&x_num * &x_den_inv), &self.q);

        // y = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)
        let y_num = &p.y * &q.y - &self.a * &p.x * &q.x;
        let y_den = utils::modulus(&(&one - &self.d * &p.x * &q.x * &p.y * &q.y), &self.q);
        let y_den_inv = utils::mod_inverse0(&y_den, &self.q);
        // let y_den_inv = utils::mod_inverse1(y_den, self.q.clone());
        // let y_den_inv = utils::mod_inverse2(y_den, self.q.clone());
        let y: BigInt = utils::modulus(&(&y_num * &y_den_inv), &self.q);

        Point { x: x, y: y }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

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
}
