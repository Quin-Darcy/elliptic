#![allow(non_snake_case)]
mod utils;
use rand::Rng;


pub fn ec_add(P: (u128, u128), Q: (u128, u128), a: u128, p: u128) -> (u128, u128) {
    let x1: u128 = P.0; let y1: u128 = P.1;
    let x2: u128 = Q.0; let y2: u128 = Q.1;
    let x3: u128; let y3: u128;
    let lambda: u128;

    if x1 == p || y1 == p {
        return (x2, y2);
    } 
    if x2 == p || y2 == p {
        return (x1, y1);
    }
    if x1 == x2 && y1 == p-y2 {
        return (p, p);
    }

    if x1 == x2 && y1 == y2 {
        lambda = ((3*utils::power(x1, 2, p)+a)*utils::mod_inv(2*y1, p)) % p;
    } else {
        lambda = (utils::mod_sub(y2, y1, p)*utils::mod_inv(utils::mod_sub(x2, x1, p), p)) % p;
    }

    x3 = utils::mod_sub(utils::power(lambda, 2, p), x1+x2, p); 
    y3 = utils::mod_sub(lambda*utils::mod_sub(x1, x3, p), y1, p);

    (x3, y3)
}

// Returns valid coefficients for an elliptic curve E: y^2 = x^3 + Ax + B
// such that 4A^3 + 27B^2 != 0 (mod p).
pub fn get_ec_coeffs(p: u128) -> (u128, u128) {
    let mut a: u128 = rand::thread_rng().gen_range(0..p);
    let mut b: u128 = rand::thread_rng().gen_range(0..p);

    loop {
        if (4*a.pow(3_u32) + 27*b.pow(2_u32)) % p != 0 {
            break;
        } else {
            a = rand::thread_rng().gen_range(0..p);
            b = rand::thread_rng().gen_range(0..p);
        }
    }
    (a, b)
}

// Given a prime p and coefficients (A, B) which define an elliptic curve E: x^3 + Ax^2 + B,
// this function returns a point on this curve.
pub fn get_curve_point(curve_coeffs: (u128, u128), p: u128) -> (u128, u128) {
    let mut x: u128;
    let mut y: u128;
    let mut z: u128;

    loop {
        x = rand::thread_rng().gen_range(0..p);
        z = (utils::power(x, 3_u128, p) + curve_coeffs.0*x + curve_coeffs.1) % p;
        loop {
            if utils::quadratic_residue(z, p) {
                break;
            } else {
                x = rand::thread_rng().gen_range(0..p);
                z = (utils::power(x, 3_u128, p) + curve_coeffs.0*x + curve_coeffs.1) % p;
            }
        }
        y = utils::mod_sqrt(z, p);
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
        sum = ec_add((x, y), point, curve_coeffs.0, p);
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
        assert_eq!(3*utils::mod_inv(3, 13) % 13, 1);
        assert_eq!(5*utils::mod_inv(5, 41) % 41, 1);
    }

    #[test]
    fn modular_subtraction_test() {
        assert_eq!(utils::mod_sub(5, 11, 13), 7);
        assert_eq!(utils::mod_sub(9, 2, 41), 7);
    }

    #[test]
    fn quadriatic_residue_test() {
        assert_eq!(utils::quadratic_residue(12, 13), true);
    }

    #[test]
    fn tonelli_shanks_test() {
        let root: u128 = utils::tonelli_shanks(5, 41);
        assert_eq!(root == 13 || root == 28, true);
    }

    #[test]
    fn modular_square_root_test() {
        let root: u128 = utils::mod_sqrt(12, 13);
        assert_eq!(utils::power(root, 2, 13), 12);
    }

    #[test]
    fn elliptic_curve_addition_test() {
        assert_eq!(ec_add((9, 7), (2, 10), 3, 13), (12, 11));
        assert_eq!(ec_add((17, 11), (39, 24), 17, 41), (10, 36));
    }

    #[test]
    fn elliptic_curve_coefficients_test() {
        let coeffs: (u128, u128) = get_ec_coeffs(13);
        let term: u128 = 4*utils::power(coeffs.0, 3, 13)
            +27*utils::power(coeffs.1, 2, 13);
        
        assert_eq!(term % 13 == 0, false);
    }

    #[test]
    fn elliptic_curve_point_test() {
        let curve_coeffs: (u128, u128) = get_ec_coeffs(13);
        let point: (u128, u128) = get_curve_point(curve_coeffs, 13);
        let lhs: u128 = utils::power(point.1, 2, 13);
        let rhs: u128 = utils::power(point.0, 3, 13)
            +curve_coeffs.0*point.0+curve_coeffs.1;
        
        assert_eq!(lhs, rhs % 13);
    }
}
