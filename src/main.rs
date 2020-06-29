use image::{ImageBuffer, RgbImage};
use rand::distributions::{Distribution, Uniform};
use raytracing::camera::*;
use raytracing::geometry::{Intersection, Point3, Sphere};
use raytracing::light::*;

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 384;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_DEPTH: usize = 50;

    let camera = Camera::default();

    let between = Uniform::from(0.0..1.0);
    let mut rng = rand::thread_rng();

    let mut world = Vec::new();
    world.push(Box::new(Sphere {
        center: Point3([0.0, 0.0, -1.0]),
        radius: 0.5,
    }) as Box<dyn Intersection>);

    world.push(Box::new(Sphere {
        center: Point3([0.0, -100.5, -1.0]),
        radius: 100.0,
    }) as Box<dyn Intersection>);

    let image: RgbImage = ImageBuffer::from_fn(IMAGE_WIDTH, IMAGE_HEIGHT, |i, j| {
        (0..SAMPLES_PER_PIXEL)
            .map(|_| {
                let u = (f64::from(i) + between.sample(&mut rng)) / f64::from(IMAGE_WIDTH - 1);
                let v = (f64::from(IMAGE_HEIGHT - j) + between.sample(&mut rng))
                    / f64::from(IMAGE_HEIGHT - 1);
                ray_color(&camera.ray(u, v), &world, &mut rng, MAX_DEPTH)
            })
            .sum::<Color>()
            .map(|c| c.sqrt())
            .into()
    });

    image.save("test.png").unwrap();
}
