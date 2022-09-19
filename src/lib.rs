use rand::Rng;

pub fn gen_curve(n: &u128) -> Vec<[u128; 2]> {
    let p: u128 = (*n).clone();
    let coeffs: [u128; 2] = get_coeffs(&p);
    let squares: Vec<Vec<u128>> = get_squares(&p);

    get_curve(&p, coeffs, squares)
}

fn get_coeffs(p: &u128) -> [u128; 2] {
    let mut a: u128;
    let mut b: u128;
    let mut c: u128;
    let coeffs: [u128; 2];
    loop {
        a = rand::thread_rng().gen_range(0..*p);
        b = rand::thread_rng().gen_range(0..*p);
        c = 4*a.pow(3) + 27*b.pow(2);

        if c % p != 0 {
            coeffs = [a, b];
            return coeffs;
        }
    }
}

fn get_squares(p: &u128) -> Vec<Vec<u128>> {
    let mut index: usize;
    let mut squares: Vec<Vec<u128>> = Vec::new();
    for _i in 0..*p {
        squares.push(Vec::new());
    }
    for x in 0..*p {
        index = ((x as u128).pow(2) % p) as usize;
        squares[index].push(x as u128);
    }
    squares
}

fn get_curve(p: &u128, coeffs: [u128; 2], squares: Vec<Vec<u128>>) -> Vec<[u128; 2]> {
    let mut index: usize;
    let mut curve: Vec<[u128; 2]> = Vec::new();
    for x in 0..*p {
        index = (((x as u128).pow(3) + coeffs[0]*(x as u128) + coeffs[1]) % p) as usize;
        if squares[index].len() > 0 {
            for y in &squares[index] {
                curve.push([x as u128, *y]);
            }
        }
    }
    curve
}



