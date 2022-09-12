use rand::Rng;

pub fn gen_curve(n: &u32) -> Vec<[u32; 2]> {
    let p: u32 = (*n).clone();
    let coeffs: [u32; 2] = get_coeffs(&p);
    let squares: Vec<Vec<u32>> = get_squares(&p);

    get_curve(&p, coeffs, squares)
}

fn get_coeffs(p: &u32) -> [u32; 2] {
    let mut a: u32;
    let mut b: u32;
    let mut c: u32;
    let coeffs: [u32; 2];
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

fn get_squares(p: &u32) -> Vec<Vec<u32>> {
    let mut index: usize;
    let mut squares: Vec<Vec<u32>> = Vec::new();
    for _i in 0..*p {
        squares.push(Vec::new());
    }
    for x in 0..*p {
        index = ((x as u32).pow(2) % p) as usize;
        squares[index].push(x as u32);
    }
    squares
}

fn get_curve(p: &u32, coeffs: [u32; 2], squares: Vec<Vec<u32>>) -> Vec<[u32; 2]> {
    let mut index: usize;
    let mut curve: Vec<[u32; 2]> = Vec::new();
    for x in 0..*p {
        index = (((x as u32).pow(3) + coeffs[0]*(x as u32) + coeffs[1]) % p) as usize;
        if squares[index].len() > 0 {
            for y in &squares[index] {
                curve.push([x as u32, *y]);
            }
        }
    }
    curve
}



