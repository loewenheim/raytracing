use raytracing::color::Color;
use raytracing::geometry::{Point3, Ray, Vec3};

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u16 = 384;
    const IMAGE_HEIGHT: u16 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u16;

    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");

    let viewport_height = 2.0;
    let viewport_width = ASPECT_RATIO * viewport_height;
    let focal_length = 1.0;

    let origin = Point3::default();
    let horizontal = Vec3::new(viewport_width, 0.0, 0.0);
    let vertical = Vec3::new(0.0, viewport_height, 0.0);
    let lower_left_corner =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, focal_length);

    for j in (0..IMAGE_HEIGHT).rev() {
        eprint!("\rScanlines remaining: {} ", j);
        for i in 0..IMAGE_WIDTH {
            let u = f64::from(i) / f64::from(IMAGE_WIDTH - 1);
            let v = f64::from(j) / f64::from(IMAGE_HEIGHT - 1);

            let ray = Ray {
                origin,
                direction: (lower_left_corner + horizontal * u + vertical * v) - origin,
            };

            println!("{}", ray.color());
        }
    }
    eprintln!("\nDone.");
}
