use image::{ImageBuffer, RgbImage};
use perlin_noise::PerlinNoise;
use raytracing::camera::{Camera, CameraOptions};
use raytracing::geometry::{Axes, Point3, Shape, Vec3};
use raytracing::materials::Material;
use raytracing::textures::Texture;
use raytracing::{pixels, ImageOptions, Object, World};
use std::sync::Arc;

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 1280;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_DEPTH: usize = 50;

    let mut rng = rand::thread_rng();

    let image_options = ImageOptions {
        width: IMAGE_WIDTH,
        height: IMAGE_HEIGHT,
        samples_per_pixel: SAMPLES_PER_PIXEL,
        max_depth: MAX_DEPTH,
    };

    let camera = Camera::new(CameraOptions {
        aperture: 0.0,
        aspect_ratio: ASPECT_RATIO,
        focus_distance: 10.0,
        looking_at: Point3([0.0, 0.0, 0.0]),
        origin: Point3([13.0, 2.0, 3.0]),
        shutter_open: 0.0,
        shutter_close: 0.01,
        vfov: 20.0,
        vup: Vec3([0.0, 1.0, 0.0]),
    });

    let mut objects = Vec::new();

    let perlin = Arc::new(PerlinNoise::new());
    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([0.0, -1000.0, 0.0]),
            radius: 1000.0,
        },
        material: Material::Lambertian {
            texture: Texture::Noise {
                noise: perlin.clone(),
                scale: 1.0,
            },
        },
    });

    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 2.0, 0.0]),
            radius: 2.0,
        },
        material: Material::Lambertian {
            texture: Texture::Noise {
                noise: perlin.clone(),
                scale: 1.0,
            },
        },
    });

    let diff_light = Material::DiffuseLight {
        emit: Texture::SolidColor(Vec3([4.0, 4.0, 4.0])),
    };

    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 7.0, 0.0]),
            radius: 2.0,
        },

        material: diff_light.clone(),
    });

    objects.push(Object {
        shape: Shape::Rectangle {
            lower_left: (3.0, 1.0),
            upper_right: (5.0, 3.0),
            height: -2.0,
            axes: Axes::XY,
        },

        material: diff_light,
    });

    let world = World::new(objects, Vec3([0.0, 0.0, 0.0]), &mut rng);
    let image: RgbImage = ImageBuffer::from_raw(
        IMAGE_WIDTH,
        IMAGE_HEIGHT,
        pixels(&camera, &world, image_options),
    )
    .unwrap();

    image.save("test.png").unwrap();
}
