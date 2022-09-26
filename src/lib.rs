mod utils;

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
    let mut lambda: u128;
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



