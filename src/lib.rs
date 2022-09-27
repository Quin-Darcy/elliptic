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
    fn modular_inverse_test() {
        assert_eq!(3*utils::tools::inv(3, 13) % 13, 1);
        assert_eq!(5*utils::tools::inv(5, 41) % 41, 1);
    }

    #[test]
    fn modular_subtraction_test() {
        assert_eq!(utils::tools::mod_sub(5, 11, 13), 7);
        assert_eq!(utils::tools::mod_sub(9, 2, 41), 7);
    }

    #[test]
    fn quadriatic_residue_test() {
        assert_eq!(utils::tools::quadratic_residue(12, 13), true);
    }

    #[test]
    fn tonelli_shanks_test() {
        let root: u128 = utils::tools::tonelli_shanks(5, 41);
        assert_eq!(root == 13 || root == 28, true);
    }

    #[test]
    fn modular_square_root_test() {
        let root: u128 = utils::tools::mod_sqrt(12, 13);
        assert_eq!(utils::tools::power(root, 2, 13), 12);
    }

    #[test]
    fn elliptic_curve_addition_test() {
        assert_eq!(utils::tools::add((9, 7), (2, 10), 3, 13), (12, 11));
        assert_eq!(utils::tools::add((17, 11), (39, 24), 17, 41), (10, 36));
    }

    #[test]
    fn elliptic_curve_coefficients_test() {
        let coeffs: (u128, u128) = utils::tools::get_coeffs(13);
        let term: u128 = 4*utils::tools::power(coeffs.0, 3, 13)
            +27*utils::tools::power(coeffs.1, 2, 13);
        
        assert_eq!(term % 13 == 0, false);
    }

    #[test]
    fn elliptic_curve_point_test() {
        let curve_coeffs: (u128, u128) = utils::tools::get_coeffs(13);
        let point: (u128, u128) = get_curve_point(curve_coeffs, 13);
        let lhs: u128 = utils::tools::power(point.1, 2, 13);
        let rhs: u128 = utils::tools::power(point.0, 3, 13)
            +curve_coeffs.0*point.0+curve_coeffs.1;
        
        assert_eq!(lhs, rhs % 13);
    }
}
