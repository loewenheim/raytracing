use image::{ImageBuffer, RgbImage};
use rand::distributions::{Distribution, Uniform};
use raytracing::camera::{Camera, CameraOptions};
use raytracing::geometry::{Point3, Shape, Vec3};
use raytracing::light::*;
use raytracing::materials::Material;
use raytracing::{Object, World};

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 640;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_DEPTH: usize = 50;

    let looking_from = Point3([3.0, 3.0, 2.0]);
    let looking_at = Point3([0.0, 0.0, -1.0]);
    let direction = looking_at - looking_from;

    let cam_options = CameraOptions {
        origin: looking_from,
        direction,
        vup: Vec3([0.0, 1.0, 0.0]),
        vfov: 20.0,
        focus_distance: direction.norm(),
        aperture: 0.0,
        aspect_ratio: ASPECT_RATIO,
    };

    let camera = Camera::new(cam_options);

    let between = Uniform::from(0.0..1.0);
    let mut rng = rand::thread_rng();
    let mut world = World::new();

    let ground = Object {
        shape: Shape::Plane {
            normal: Vec3([0.0, 1.0, 0.0]).normed(),
            offset: -0.5,
        },
        material: Material::Lambertian {
            albedo: Color::new(0.8, 0.0, 0.8),
        },
    };

    let ceiling = Object {
        shape: Shape::Plane {
            normal: Vec3([0.0, -1.0, 0.0]).normed(),
            offset: 2.0,
        },
        material: Material::Lambertian {
            albedo: Color::new(0.0, 0.0, 0.0),
        },
    };

    let sphere2 = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 0.0, -1.0]),
            radius: 0.5,
        },
        material: Material::Lambertian {
            albedo: Color::new(0.1, 0.2, 0.5),
        },
    };

    let sphere3 = Object {
        shape: Shape::Sphere {
            center: Point3([1.2, 0.0, -1.0]),
            radius: 0.5,
        },
        material: Material::Metal {
            albedo: Color::new(0.8, 0.6, 0.2),
            fuzz: 0.0,
        },
    };

    let sphere4 = Object {
        shape: Shape::Sphere {
            center: Point3([-1.2, 0.0, -1.0]),
            radius: 0.5,
        },
        material: Material::Dielectric {
            refraction_index: 1.5,
        },
    };

    let sphere5 = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 0.5, -3.0]),
            radius: 1.0,
        },
        material: Material::Metal {
            albedo: Color::new(0.5, 0.5, 0.5),
            fuzz: 0.02,
        },
    };

    world.push(ground);
    world.push(ceiling);
    world.push(sphere2);
    world.push(sphere3);
    world.push(sphere4);
    world.push(sphere5);

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
