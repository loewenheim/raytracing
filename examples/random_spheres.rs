use image::{ImageBuffer, RgbImage};
use indicatif::{ProgressBar, ProgressStyle};
use rand::Rng;
use raytracing::camera::{Camera, CameraOptions};
use raytracing::geometry::{Point3, Shape, Vec3};
use raytracing::materials::Material;
use raytracing::textures::Texture;
use raytracing::{pixels, BvhNode, ImageOptions, Object, Optics, World};
use std::sync::{Arc, Mutex};

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 640;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_REFLECTIONS: usize = 50;
    let progress_bar = ProgressBar::new((IMAGE_WIDTH * IMAGE_HEIGHT).into());
    progress_bar.set_style(ProgressStyle::default_bar().template("{wide_bar} {percent:2}%"));
    let progress_bar = Arc::new(Mutex::new(progress_bar));

    let image_options = ImageOptions {
        width: IMAGE_WIDTH,
        height: IMAGE_HEIGHT,
        samples_per_pixel: SAMPLES_PER_PIXEL,
        max_reflections: MAX_REFLECTIONS,
    };

    let camera = Camera::new(CameraOptions {
        origin: Point3([13.0, 2.0, 3.0]),
        looking_at: Point3([0.0, 0.0, 0.0]),
        vup: Vec3([0.0, 1.0, 0.0]),
        vfov: 20.0,
        focus_distance: 10.0,
        aperture: 0.1,
        aspect_ratio: ASPECT_RATIO,
        shutter_open: 0.0,
        shutter_close: 0.5,
    });

    let mut rng = rand::thread_rng();

    let mut objects: Vec<Box<dyn Optics>> = Vec::new();

    objects.push(Box::new(Object {
        shape: Box::new(Shape::Sphere {
            center: Point3([0.0, -1000.0, 0.0]),
            radius: 1000.0,
        }),

        material: Material::Lambertian {
            texture: Texture::Checkered {
                even: Box::new(Texture::SolidColor(Vec3([0.2, 0.3, 0.1]))),
                odd: Box::new(Texture::SolidColor(Vec3([0.9, 0.9, 0.9]))),
            },
        },
    }));

    for a in -11..11 {
        for b in -11..11 {
            let center = Point3([
                a as f64 + 0.9 * rng.gen::<f64>(),
                0.2,
                b as f64 + 0.9 * rng.gen::<f64>(),
            ]);

            if center.dist(&Point3([4.0, 0.2, 0.0])) > 0.9 {
                let roll: f64 = rng.gen();
                let (center2, material) = if roll < 0.8 {
                    let center2 = center + Vec3([0.0, rng.gen_range(0.0, 0.5), 0.0]);
                    let texture = Texture::SolidColor(
                        Vec3::random(0.0..1.0, &mut rng) + Vec3::random(0.0..1.0, &mut rng),
                    );
                    (center2, Material::Lambertian { texture })
                } else if roll < 0.95 {
                    let albedo = Vec3::random(0.5..1.0, &mut rng);
                    let fuzz = rng.gen_range(0.0, 0.5);

                    (center, Material::Metal { albedo, fuzz })
                } else {
                    (center, Material::glass())
                };

                let shape = Box::new(Shape::sphere(0.2, (0.0, center), (1.0, center2)));
                objects.push(Box::new(Object { shape, material }));
            }
        }
    }

    objects.push(Box::new(Object {
        shape: Box::new(Shape::Sphere {
            center: Point3([0.0, 1.0, 0.0]),
            radius: 1.0,
        }),

        material: Material::glass(),
    }));

    objects.push(Box::new(Object {
        shape: Box::new(Shape::Sphere {
            center: Point3([-4.0, 1.0, 0.0]),
            radius: 1.0,
        }),

        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.4, 0.2, 0.1])),
        },
    }));

    objects.push(Box::new(Object {
        shape: Box::new(Shape::Sphere {
            center: Point3([4.0, 1.0, 0.0]),
            radius: 1.0,
        }),

        material: Material::Metal {
            albedo: Vec3([0.7, 0.6, 0.5]),
            fuzz: 0.0,
        },
    }));

    let world = World {
        objects: Box::new(BvhNode::create(objects, &mut rng)),
        background_color: Vec3([1.0, 1.0, 1.0]),
    };

    let image: RgbImage = ImageBuffer::from_raw(
        IMAGE_WIDTH,
        IMAGE_HEIGHT,
        pixels(&camera, &world, image_options, Arc::clone(&progress_bar)),
    )
    .unwrap();

    image.save("test.png").unwrap();
}
