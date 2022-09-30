#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
use rand::Rng;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{Zero, One};
use num_traits::cast::ToPrimitive;


pub fn mod_inv(a: &BigUint, p: &BigUint) -> BigUint {
    a.modpow(&(p-2_u32), &p)
}

pub fn mod_sub(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    (a+(p-(b % p))) % p
}

pub fn quadratic_residue(z: &BigUint, p: &BigUint) -> bool {
    let exp: BigUint = (p-1_u32)/2_u32;
    let legendre: BigUint = z.modpow(&exp, p);

    if legendre == Zero::zero() || legendre == One::one() {
        return true;
    } else {
        return false;
    }
}

pub fn tonelli_shanks(n: &BigUint, p: &BigUint) -> BigUint {
    let q: BigUint;
    let mut z: BigUint;
    let mut s: u32 = 0;
    let mut d: BigUint = p-1_u32;

    // Express p-1 = q2^s
    while &d % 2_u32 == Zero::zero() {
        d = d / 2_u32;
        s += 1;
    }
    q = (p-1_u32) / 2_u32.pow(s);

    // Obtain a quadratic non-residue, z.
    loop {
        z = rand::thread_rng().gen_biguint_range(&Zero::zero(), &p);
        if !quadratic_residue(&z, &p) {
            break;
        }
    }

    let mut m: BigUint = BigUint::from(s.clone());
    let mut c: BigUint = z.modpow(&q, &p);
    let mut t: BigUint = n.modpow(&q, &p);
    let mut r: BigUint = n.modpow(&((&q+1_u32)/2_u32), &p);
    let mut b: BigUint;

    loop {
        if t == Zero::zero() {
            return Zero::zero();
        } else if t == One::one() {
            return r;
        } else {
            let mut i: u32 = 0;
            while t.modpow(&BigUint::from(2_u32.pow(i)), &p) != One::one() {
                i += 1;
            }
            if BigUint::from(i) == m {
                return One::one();
            }

            b = c.modpow(&BigUint::from(2_u32).modpow(&(&m-BigUint::from(i)-1_u32), &p), &p);
            m = BigUint::from(i);
            c = b.modpow(&BigUint::from(2_u32), &p);
            t = t*b.modpow(&BigUint::from(2_u32), &p) % p;
            r = r*&b % p; 
        }
    }
}

pub fn mod_sqrt(z: &BigUint, p: &BigUint) -> BigUint {
    // Multiples of p only have 0 as its roots
    if z.is_zero() {
        return Zero::zero();
    }

    // If p = 3(mod 4), then the square root is calculated directly
    if p % 4_u32 == BigUint::from(3_u32) {
        return z.modpow(&((p+1_u32)/4_u32), &p); 
    }

    // If p = 1(mod 4), we must use Tonelli-Shanks Algorithm
    tonelli_shanks(&z, &p)
}
