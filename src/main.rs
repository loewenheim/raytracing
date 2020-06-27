use raytracing::Color;
use raytracing::geometry::{Vec3, Point3, Ray};

fn main() {
    const IMAGE_WIDTH: u16 = 256;
    const IMAGE_HEIGHT: u16 = 256;

    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");
    
    for j in (0..IMAGE_HEIGHT).rev() {
        eprint!("\rScanlines remaining: {} ", j);
        for i in 0..IMAGE_WIDTH {
            let r = f64::from(i) / f64::from(IMAGE_WIDTH-1);
            let g = f64::from(j) / f64::from(IMAGE_HEIGHT-1);
            let b = 0.25;

            let color = Color(Vec3::new(r, g, b));
            println!("{}", color);
        }
    }
    eprintln!("\nDone.");
}
