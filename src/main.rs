use image::{ImageBuffer, RgbImage};
use rand::distributions::{Distribution, Uniform};
use raytracing::camera::*;
use raytracing::geometry::{Point3, Sphere, Vec3};
use raytracing::light::*;
use raytracing::materials::{Dielectric, Lambertian, Metal};
use raytracing::{MaterialObject, Reflective};

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 640;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_DEPTH: usize = 50;

    let origin = Point3([3.0, 3.0, 2.0]);
    let looking_at = Point3([0.0, 0.0, -1.0]);
    let looking_direction = looking_at - origin;
    let vup = Vec3([0.0, 1.0, 0.0]);
    let vfov = 20.0;
    let focus_dist = looking_direction.norm();
    let aperture = 2.0;
    let camera = Camera::new(
        origin,
        looking_direction,
        vup,
        vfov,
        ASPECT_RATIO,
        aperture,
        focus_dist,
    );

    let between = Uniform::from(0.0..1.0);
    let mut rng = rand::thread_rng();
    let mut world = Vec::new();

    let sphere1 = MaterialObject::new(
        Sphere {
            center: Point3([0.0, -100.5, -1.0]),
            radius: 100.0,
        },
        Lambertian {
            albedo: Color::new(0.8, 0.0, 0.8),
        },
    );

    let sphere2 = MaterialObject::new(
        Sphere {
            center: Point3([0.0, 0.0, -1.0]),
            radius: 0.5,
        },
        Lambertian {
            albedo: Color::new(0.1, 0.2, 0.5),
        },
    );

    let sphere3 = MaterialObject::new(
        Sphere {
            center: Point3([1.0, 0.0, -1.0]),
            radius: 0.5,
        },
        Metal {
            albedo: Color::new(0.8, 0.6, 0.2),
            fuzz: 0.0,
        },
    );

    let sphere4 = MaterialObject::new(
        Sphere {
            center: Point3([-1.0, 0.0, -1.0]),
            radius: 0.5,
        },
        Dielectric { ref_index: 1.5 },
    );

    let sphere5 = MaterialObject::new(
        Sphere {
            center: Point3([0.0, 0.5, -3.0]),
            radius: 1.0,
        },
        Metal {
            albedo: Color::new(0.5, 0.5, 0.5),
            fuzz: 0.02,
        },
    );

    world.push(Box::new(sphere1) as Box<dyn Reflective<_>>);
    world.push(Box::new(sphere2) as Box<dyn Reflective<_>>);
    world.push(Box::new(sphere3) as Box<dyn Reflective<_>>);
    world.push(Box::new(sphere4) as Box<dyn Reflective<_>>);
    world.push(Box::new(sphere5) as Box<dyn Reflective<_>>);

    let image: RgbImage = ImageBuffer::from_fn(IMAGE_WIDTH, IMAGE_HEIGHT, |i, j| {
        (0..SAMPLES_PER_PIXEL)
            .map(|_| {
                let u = (f64::from(i) + between.sample(&mut rng)) / f64::from(IMAGE_WIDTH - 1);
                let v = (f64::from(IMAGE_HEIGHT - j) + between.sample(&mut rng))
                    / f64::from(IMAGE_HEIGHT - 1);
                ray_color(&camera.ray(u, v, &mut rng), &world, &mut rng, MAX_DEPTH)
            })
            .sum::<Color>()
            .map(|c| c.sqrt())
            .into()
    });

    image.save("test.png").unwrap();
}
