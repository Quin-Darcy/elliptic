mod utils;
use rand::Rng;

// Given a prime p and coefficients (A, B) which define an elliptic curve E: x^3 + Ax^2 + B,
// this function returns a point on this curve.
pub fn get_curve_point(curve_coeffs: (u128, u128), p: u128) -> (u128, u128) {
    let mut x: u128;
    let mut y: u128;
    let mut z: u128;

    loop {
        x = rand::thread_rng().gen_range(0..p);
        z = (utils::tools::power(x, 3_u128, p) + curve_coeffs.0*x + curve_coeffs.1) % p;
        loop {
            if utils::tools::quadratic_residue(z, p) {
                break;
            } else {
                x = rand::thread_rng().gen_range(0..p);
                z = (utils::tools::power(x, 3_u128, p) + curve_coeffs.0*x + curve_coeffs.1) % p;
            }
        }
        y = utils::tools::mod_sqrt(z, p);
        if y != 1 && y != 0{
            break;
        }
    }
    (x, y)
}

pub fn order(point: (u128, u128), curve_coeffs: (u128, u128), p: u128) -> u128 {
    let mut i: u128 = 0;
    let mut sum: (u128, u128) = (0, 0);
    let mut x: u128 = point.0.clone();
    let mut y: u128 = point.1.clone();
    
    while sum != (p, p) {
        sum = utils::tools::add((x, y), point, curve_coeffs.0, p);
        x = sum.0;
        y = sum.1;
        i += 1;
    }
    i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inv_test() {
        assert_eq!(utils::tools::inv(3_u128, 13_u128), 9_u128);
    }

    #[test]
    fn mod_sub_test() {
        assert_eq!(utils::tools::mod_sub(5_u128, 11_u128, 13_u128), 7_u128);
    }

    #[test]
    fn add_test() {
        assert_eq!(utils::tools::add((9_u128, 7_u128), (2_u128, 10_u128), 3_u128, 13_u128), (12_u128, 11_u128));
    }

    #[test]
    fn mod_sqrt_test() {
        let root: u128 = utils::tools::mod_sqrt(12_u128, 13_u128);
        assert_eq!(utils::tools::power(root, 2_u128, 13_u128), 12_u128);
    }

    #[test]
    fn quadriatic_residue_test() {
        assert_eq!(utils::tools::quadratic_residue(12_u128, 13_u128), true);
    }

    #[test]
    fn tonelli_shanks_test() {
        assert_eq!(utils::tools::tonelli_shanks(12_u128, 13_u128), 5_u128);
    }

    #[test]
    fn coeffs_test() {
        let coeffs: (u128, u128) = utils::tools::get_coeffs(13_u128);
        let term: u128 = 4*utils::tools::power(coeffs.0, 3_u128, 13_u128)+27*utils::tools::power(coeffs.1, 2_u128, 13_u128);
        assert_eq!(term % 13_u128 == 0, false);
    }

    #[test]
    fn curve_point_test() {
        let curve_coeffs: (u128, u128) = utils::tools::get_coeffs(13_u128);
        let point: (u128, u128) = get_curve_point(curve_coeffs, 13_u128);
        let lhs: u128 = utils::tools::power(point.1, 2_u128, 13_u128);
        let rhs: u128 = utils::tools::power(point.0, 3_u128, 13_u128)+curve_coeffs.0*point.0+curve_coeffs.1;
        assert_eq!(lhs, rhs);
    }
}
