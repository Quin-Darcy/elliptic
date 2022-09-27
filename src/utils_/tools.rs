#![allow(non_snake_case)]
use rand::Rng;
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;


pub fn power(base: u128, exp: u128, p: u128) -> u128 {
    BigUint::from(base).modpow(&BigUint::from(exp), 
                               &BigUint::from(p)).to_u128() .unwrap()
}

pub fn inv(a: u128, p: u128) -> u128 {
    power(a, p-2, p)
}

pub fn mod_sub(a: u128, b: u128, p: u128) -> u128 {
    (a+(p-(b % p))) % p
}

pub fn add(P: (u128, u128), Q: (u128, u128), a: u128, p: u128) -> (u128, u128) {
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
        lambda = ((3*power(x1, 2, p)+a)*inv(2*y1, p)) % p;
    } else {
        lambda = (mod_sub(y2, y1, p)*inv(mod_sub(x2, x1, p), p)) % p;
    }

    x3 = mod_sub(power(lambda, 2, p), x1+x2, p); 
    y3 = mod_sub(lambda*mod_sub(x1, x3, p), y1, p);

    (x3, y3)
}

// This function returns one of the square roots of z over F_p.
pub fn mod_sqrt(z: u128, p: u128) -> u128 {
    // Multiples of p only have 0 as its roots
    if z == 0 {
        return 0;
    }

    // If p = 3(mod 4), then the square root is calculated directly
    if p % 4 == 3 {
        return power(z, (p+1)/4, p);
    }

    // If p = 1(mod 4), we must use Tonelli-Shanks Algorithm
    tonelli_shanks(z, p)
}

// Uses Euler's criterion to determine if given number z is a quadratic
// residue in the field F_p, i.e., there exists n in F_p s.t. n^2 = z (mod p).
pub fn quadratic_residue(z: u128, p: u128) -> bool {
    let legendre: u128 = power(z, (p-1)/2, p);

    if legendre == 0 || legendre == 1 {
        return true;
    } else {
        return false;
    }
}

// See https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
pub fn tonelli_shanks(n: u128, p: u128) -> u128 {
    let q: u128;
    let mut z: u128;
    let mut s: u32 = 0;
    let mut d: u128 = p-1;

    // Express p-1 = q2^s
    while d % 2 == 0 {
        d >>= 1;
        s += 1;
    }
    q = (p-1) / 2_u128.pow(s);

    // Obtain a quadratic non-residue, qnr.
    loop {
        z = rand::thread_rng().gen_range(0..p);
        if !quadratic_residue(z, p) {
            break;
        }
    }

    let mut m: u128 = s.clone() as u128;
    let mut c: u128 = power(z, q, p);
    let mut t: u128 = power(n, q, p);
    let mut r: u128 = power(n, (q+1)/2, p);
    let mut b: u128;

    loop {
        if t == 0 {
            return 0;
        } else if t == 1 {
            return r;
        } else {
            let mut i: u32 = 0;
            while power(t, 2_u128.pow(i), p) != 1 {
                i += 1;
            }
            if (i as u128) == m {
                return 1;
            }

            b = power(c, power(2_u128, m-(i as u128)-1, p), p);
            m = i as u128;
            c = power(b, 2_u128, p);
            t = t*power(b, 2_u128, p) % p;
            r = r*b % p;
        }
    }
}

// Returns valid coefficients for an elliptic curve E: y^2 = x^3 + Ax + B
// such that 4A^3 + 27B^2 != 0 (mod p).
pub fn get_coeffs(p: u128) -> (u128, u128) {
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
