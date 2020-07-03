use image::{ImageBuffer, RgbImage};
use rayon::prelude::*;
use raytracing::camera::{Camera, CameraOptions};
use raytracing::color::Color;
use raytracing::geometry::{Point3, Shape, Vec3};
use raytracing::materials::Material;
use raytracing::{pixel, random_world, BvhNode, Object};

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 1280;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_DEPTH: usize = 50;

    let looking_from = Point3([13.0, 2.0, 3.0]);
    let looking_at = Point3([0.0, 0.0, 0.0]);
    let direction = looking_at - looking_from;

    let cam_options = CameraOptions {
        origin: looking_from,
        direction,
        vup: Vec3([0.0, 1.0, 0.0]),
        vfov: 20.0,
        focus_distance: 10.0,
        aperture: 0.1,
        aspect_ratio: ASPECT_RATIO,
    };

    let camera = Camera::new(cam_options);

    let mut rng = rand::thread_rng();

    let mut world = random_world(&mut rng);

    let sphere1 = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 1.0, 0.0]),
            radius: 1.0,
        },

        material: Material::glass(),
    };

    let sphere2 = Object {
        shape: Shape::Sphere {
            center: Point3([-4.0, 1.0, 0.0]),
            radius: 1.0,
        },

        material: Material::Lambertian {
            albedo: Color::new(0.4, 0.2, 0.1),
        },
    };

    let sphere3 = Object {
        shape: Shape::Sphere {
            center: Point3([4.0, 1.0, 0.0]),
            radius: 1.0,
        },

        material: Material::Metal {
            albedo: Color::new(0.7, 0.6, 0.5),
            fuzz: 0.0,
        },
    };

    world.push(sphere1);
    world.push(sphere2);
    world.push(sphere3);

    let world = BvhNode::create(&mut world, &mut rng);

    let pixels: Vec<[u8; 3]> = (0..IMAGE_HEIGHT)
        .into_par_iter()
        .rev()
        .flat_map(move |j| {
            let world = world.clone();
            (0..IMAGE_WIDTH).into_par_iter().map(move |i| {
                pixel(
                    &camera,
                    &world,
                    (i, j),
                    (IMAGE_WIDTH, IMAGE_HEIGHT),
                    (0.0, 1.0),
                    SAMPLES_PER_PIXEL,
                    MAX_DEPTH,
                )
            })
        })
        .collect();

    let pixels = pixels.concat();

    let image: RgbImage = ImageBuffer::from_raw(IMAGE_WIDTH, IMAGE_HEIGHT, pixels).unwrap();

    image.save("test.png").unwrap();
}
