#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
use rand::Rng;
use num_traits::{Zero, One};
use num_traits::cast::ToPrimitive;
use num_bigint::{BigUint, RandBigInt};
mod utils;


pub fn ec_add(P: &(BigUint, BigUint), Q: &(BigUint, BigUint), a: &BigUint, p: &BigUint) -> (BigUint, BigUint) {
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
