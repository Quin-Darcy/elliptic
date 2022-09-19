use rand::Rng;
use num_bigint::{BigUint, RandBigInt};
use num_traits::cast::ToPrimitive;
use num_traits::identities::{Zero};


pub fn gen_curve(n: BigUint) -> Vec<[BigUint; 2]> {
    let p: BigUint = n.clone();
    let coeffs: [BigUint; 2] = get_coeffs(&p);
    let squares: Vec<Vec<BigUint>> = get_squares(&p);
    
    get_curve(&p, coeffs, squares)
}

fn get_coeffs(p: &BigUint) -> [BigUint; 2] {
    let mut a: BigUint;
    let mut b: BigUint;
    let mut c: BigUint;
    let coeffs: [BigUint; 2];
    loop {
        a = rand::thread_rng().gen_biguint_range(&BigUint::zero(), &p.clone());
        b = rand::thread_rng().gen_biguint_range(&BigUint::zero(), &p.clone());
        c = (4_u32)*BigUint::pow(&a, 3_u32) + (27_u32)*BigUint::pow(&b, 2_u32);

        if !(c % p.clone()).is_zero() {
            coeffs = [a, b];
            return coeffs;
        }
    }
}

fn get_squares(p: &BigUint) -> Vec<Vec<BigUint>> {
    let mut index: usize;
    let mut squares: Vec<Vec<BigUint>> = Vec::new();
    let mut i: BigUint = BigUint::zero();
    
    while i < *p {
        squares.push(Vec::new());
        i += BigUint::from(1_u32);
    }
    i = BigUint::zero();
    
    while i < *p {
        index = BigUint::modpow(&i, &BigUint::from(2_u32), p).to_usize().unwrap();
        squares[index].push(i.clone());
        i += BigUint::from(1_u32);
    }
    squares
}

fn get_curve(p: &BigUint, coeffs: [BigUint; 2], squares: Vec<Vec<BigUint>>) -> Vec<[BigUint; 2]> {
    let mut index: usize;
    let mut curve: Vec<[BigUint; 2]> = Vec::new();
    let mut x: BigUint = BigUint::zero();

    while x < *p {
        index = ((x.pow(3) + coeffs[0].clone()*&x + coeffs[1].clone()) % p).to_usize().unwrap();
        if squares[index].len() > 0 {
            for y in &squares[index] {
                curve.push([x.clone(), (*y).clone()]);
            }
        }
        x += BigUint::from(1_u32);
    }
    curve
}
