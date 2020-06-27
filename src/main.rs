fn main() {
    const IMAGE_WIDTH: u16 = 256;
    const IMAGE_HEIGHT: u16 = 256;

    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");
    
    for j in (0..IMAGE_HEIGHT).rev() {
        for i in 0..IMAGE_WIDTH {
            let r = f64::from(i) / f64::from(IMAGE_WIDTH-1);
            let g = f64::from(j) / f64::from(IMAGE_HEIGHT-1);
            let b = 0.25f64;

            let ir = (r * 256.0) as u8;
            let ig = (g * 256.0) as u8;
            let ib = (b * 256.0) as u8;

            println!("{} {} {}", ir, ig, ib);
        }
    }
}
