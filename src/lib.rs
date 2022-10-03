#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
use num_traits::{Zero, One};
use num_bigint::{BigUint, RandBigInt};
mod utils;


fn ec_scalar_mult(k: &BigUint, P: &(BigUint, BigUint), a: &BigUint, p: &BigUint) -> (BigUint, BigUint) {
    let bit_str: String = format!("{:b}", k).to_string();
    let bit_chars: Vec<char> = bit_str.chars().collect::<Vec<_>>();
    let mut q: (BigUint, BigUint) = P.clone();
    
    let l: usize = bit_chars.len();
    for i in 1..l {
        q = ec_add(&q, &q, a, p);
        if bit_chars[i] == '1' {
            q = ec_add(&q, P, a, p);
        }
    }
    q
}

fn ec_add(P: &(BigUint, BigUint), Q: &(BigUint, BigUint), a: &BigUint, p: &BigUint) -> (BigUint, BigUint) {
    let x1: &BigUint = &P.0; let y1: &BigUint = &P.1;
    let x2: &BigUint = &Q.0; let y2: &BigUint = &Q.1;
    let x3: BigUint; let y3: BigUint;
    let lambda: BigUint;

    if x1 == p || y1 == p {
        return ((*x2).clone(), (*y2).clone());
    } 
    if x2 == p || y2 == p {
        return ((*x1).clone(), (*y1).clone());
    }
    if x1 == x2 && y1 == &(p-y2) {
        return ((*p).clone(), (*p).clone());
    }

    if x1 == x2 && y1 == y2 {
        lambda = ((3_u32*x1.modpow(&BigUint::from(2_u32), &p)+a)*utils::mod_inv(&(BigUint::from(2_u32)*y1), p)) % p; 
    } else {
        lambda = (utils::mod_sub(&y2, &y1, p)*utils::mod_inv(&utils::mod_sub(&x2, &x1, p), p)) % p; 
    }

    x3 = utils::mod_sub(&lambda.modpow(&BigUint::from(2_u32), p), &(x1+x2), p); 
    y3 = utils::mod_sub(&(lambda*utils::mod_sub(&x1, &x3, p)), &y1, p);

    (x3, y3)
}

pub fn get_ec_coeffs(p: &BigUint) -> (BigUint, BigUint) {
    let mut a: BigUint = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);
    let mut b: BigUint = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);

    loop {
        let coeff_cond: BigUint = 4_u32*a.modpow(&BigUint::from(3_u32), &p)
                                    +27_u32*b.modpow(&BigUint::from(2_u32), &p);
        
        if coeff_cond % p != Zero::zero() {
            break;
        } else {
            a = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);
            b = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);
        }
    }
    (a, b)
}

pub fn get_ec_point(curve_coeffs: &(BigUint, BigUint), p: &BigUint) -> (BigUint, BigUint) {
    let mut x: BigUint;
    let mut y: BigUint;
    let mut z: BigUint;

    loop {
        x = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);
        z = (x.modpow(&BigUint::from(3_u32), &p)+&curve_coeffs.0*&x+&curve_coeffs.1) % p; 
        loop {
            if utils::quadratic_residue(&z, &p) {
                break;
            } else {
                x = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);
                z = (x.modpow(&BigUint::from(3_u32), &p)+&curve_coeffs.0*&x+&curve_coeffs.1) % p;
            }
        }
        y = utils::mod_sqrt(&z, &p);
        if y != One::one() && y != Zero::zero() {
            break;
        }
    }
    (x, y)
}

pub fn ec_point_order(P: &(BigUint, BigUint), a: &BigUint, p: &BigUint) -> BigUint {
    let mut scalar: BigUint = BigUint::from(1_u32);
    let mut current_multiple: (BigUint, BigUint);
    let zero: (BigUint, BigUint) = (p.clone(), p.clone());

    loop {
        current_multiple = ec_scalar_mult(&scalar, P, a, p);
        if current_multiple == zero || scalar >= 2_u32*p+1_u32 {
            break;
        } else {
            scalar += 1_u32;
        }
    }
    scalar
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modular_inverse_test() {
        assert_eq!(3_u32*utils::mod_inv(&BigUint::from(3_u32), &BigUint::from(13_u32)) % 13_u32, One::one());
        assert_eq!(5_u32*utils::mod_inv(&BigUint::from(5_u32), &BigUint::from(41_u32)) % 41_u32, One::one());
    }

    #[test]
    fn modular_subtraction_test() {
        assert_eq!(utils::mod_sub(&BigUint::from(5_u32), &BigUint::from(11_u32), &BigUint::from(13_u32)), BigUint::from(7_u32));
        assert_eq!(utils::mod_sub(&BigUint::from(9_u32), &BigUint::from(2_u32), &BigUint::from(41_u32)), BigUint::from(7_u32));
    }

    #[test]
    fn quadriatic_residue_test() {
        assert_eq!(utils::quadratic_residue(&BigUint::from(12_u32), &BigUint::from(13_u32)), true);
    }

    #[test]
    fn tonelli_shanks_test() {
        let root: BigUint = utils::tonelli_shanks(&BigUint::from(5_u32), &BigUint::from(41_u32));
        assert_eq!(root == BigUint::from(13_u32) || root == BigUint::from(28_u32), true);
    }

    #[test]
    fn modular_square_root_test() {
        let root: BigUint = utils::mod_sqrt(&BigUint::from(12_u32), &BigUint::from(13_u32));
        assert_eq!(root.modpow(&BigUint::from(2_u32), &BigUint::from(13_u32)), BigUint::from(12_u32));
    }

    #[test]
    fn elliptic_curve_addition_test() {
        let tup1: (BigUint, BigUint) = (BigUint::from(9_u32), BigUint::from(7_u32));
        let tup2: (BigUint, BigUint) = (BigUint::from(2_u32), BigUint::from(10_u32));
        let tup3: (BigUint, BigUint) = (BigUint::from(17_u32), BigUint::from(11_u32));
        let tup4: (BigUint, BigUint) = (BigUint::from(39_u32), BigUint::from(24_u32));

        assert_eq!(ec_add(&tup1, &tup2, &BigUint::from(3_u32), &BigUint::from(13_u32)), (BigUint::from(12_u32), BigUint::from(11_u32)));
        assert_eq!(ec_add(&tup3, &tup4, &BigUint::from(17_u32), &BigUint::from(41_u32)), (BigUint::from(10_u32), BigUint::from(36_u32)));
    }

    #[test]
    fn elliptic_curve_scalar_multiplication_test() {
        let point: (BigUint, BigUint) = (BigUint::from(2_u32), BigUint::from(10_u32));
        let scalar: BigUint = BigUint::from(5_u32);
        let coeff: BigUint = BigUint::from(3_u32);
        let p: BigUint = BigUint::from(13_u32);

        assert_eq!(ec_scalar_mult(&scalar, &point, &coeff, &p), (BigUint::from(1_u32), BigUint::from(5_u32)));
    }

    #[test]
    fn elliptic_curve_coefficients_test() {
        let coeffs: (BigUint, BigUint) = get_ec_coeffs(&BigUint::from(13_u32));
        let term: BigUint = 4_u32*(coeffs.0).modpow(&BigUint::from(3_u32), &BigUint::from(13_u32))
                            +27_u32*(coeffs.1).modpow(&BigUint::from(2_u32), &BigUint::from(13_u32));
        //let term: u128 = 4*utils::power(coeffs.0, 3, 13)
        //    +27*utils::power(coeffs.1, 2, 13);
        
        assert_eq!((term % 13_u32).is_zero(), false);
    }

    #[test]
    fn elliptic_curve_point_test() {
        let curve_coeffs: (BigUint, BigUint) = get_ec_coeffs(&BigUint::from(13_u32));
        let point: (BigUint, BigUint) = get_ec_point(&curve_coeffs, &BigUint::from(13_u32));
        let lhs: BigUint = (point.1).modpow(&BigUint::from(2_u32), &BigUint::from(13_u32));
        let rhs: BigUint = (point.0).modpow(&BigUint::from(3_u32), &BigUint::from(13_u32))
                            +curve_coeffs.0*point.0+curve_coeffs.1;
        
        assert_eq!(lhs, rhs % 13_u32);
    }

    #[test]
    fn elliptic_curve_point_order_test() {
        let p: BigUint = BigUint::from(41_u32);
        let curve_coeffs: (BigUint, BigUint) = get_ec_coeffs(&p);
        let point: (BigUint, BigUint) = get_ec_point(&curve_coeffs, &p);
        let order: BigUint = ec_point_order(&point, &curve_coeffs.0, &p);
        let zero: (BigUint, BigUint) = (p.clone(), p.clone());

        assert_eq!(ec_scalar_mult(&order, &point, &curve_coeffs.0, &p), zero);
    }
}
