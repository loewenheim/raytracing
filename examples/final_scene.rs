use image::{ImageBuffer, RgbImage};
use rand::Rng;
use raytracing::camera::{Camera, CameraOptions};
use raytracing::geometry::{Point3, Shape, Vec3};
use raytracing::materials::Material;
use raytracing::textures::Texture;
use raytracing::{pixels, ImageOptions, Object, Scene};

fn main() {
    const ASPECT_RATIO: f64 = 1.0;
    const IMAGE_WIDTH: u32 = 500;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_REFLECTIONS: usize = 50;
    const BOXES_PER_SIDE: usize = 20;

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

    let mut objects = Vec::new();

    for i in 0..BOXES_PER_SIDE {
        for j in 0..BOXES_PER_SIDE {
            let w = 100.0;
            let x0 = -1000.0 + i as f64 * w;
            let y0 = 0.0;
            let z0 = -1000.0 + j as f64 * w;
            let x1 = x0 + w;
            let y1 = rng.gen_range(1.0, 101.0);
            let z1 = z0 + w;

            objects.push(Object {
                shape: Shape::new_box(Point3([x0, y0, z0]), Point3([x1, y1, z1])),
                material: Material::Lambertian {
                    texture: Texture::SolidColor(Vec3([0.48, 0.83, 0.53])),
                },
            });
        }
    }

    // light
    objects.push(Object {
        shape: Shape::rectangle(Point3([123.0, 554.0, 147.0]), Point3([423.0, 554.0, 412.0])),
        material: Material::DiffuseLight {
            emit: Texture::SolidColor(Vec3([7.0, 7.0, 7.0])),
        },
    });

    // moving sphere
    objects.push(Object {
        shape: Shape::sphere(
            50.0,
            (0.0, Point3([400.0, 400.0, 200.0])),
            (1.0, Point3([430.0, 400.0, 200.0])),
        ),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.7, 0.3, 0.1])),
        },
    });

    // glass sphere
    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([260.0, 150.0, 45.0]),
            radius: 50.0,
        },

        material: Material::glass(),
    });

    // metal sphere
    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 150.0, 145.0]),
            radius: 50.0,
        },

        material: Material::Metal {
            albedo: Vec3([0.8, 0.8, 0.9]),
            fuzz: 10.0,
        },
    });

    // glass sphere with blue smoke inside
    let sphere = Shape::Sphere {
        center: Point3([360.0, 150.0, 145.0]),
        radius: 70.0,
    };

    objects.push(Object {
        shape: sphere,
        material: Material::glass(),
    });

    objects.push(Object {
        shape: Shape::ConstantMedium {
            boundary: Box::new(sphere),
            density: 0.2,
        },

        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.2, 0.4, 0.9])),
        },
    });

    // huge ball full of white mist
    objects.push(Object {
        shape: Shape::ConstantMedium {
            boundary: Box::new(Shape::Sphere {
                center: Point3([0.0, 0.0, 0.0]),
                radius: 5000.0,
            }),
            density: 0.0001,
        },

        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([1.0, 1.0, 1.0])),
        },
    });

    // earth
    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([400.0, 200.0, 400.0]),
            radius: 100.0,
        },

        material: Material::Lambertian {
            texture: Texture::image("earthmap.jpg"),
        },
    });

    // noise sphere
    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([220.0, 280.0, 300.0]),
            radius: 80.0,
        },

        material: Material::Lambertian {
            texture: Texture::Noise {
                noise: std::sync::Arc::new(noise::Perlin::new()),
                scale: 0.1,
            },
        },
    });

    for _ in 0..1000 {
        objects.push(Object {
            shape: Shape::Sphere {
                center: Point3::random(0.0..165.0, &mut rng),
                radius: 10.0,
            },

            material: Material::Lambertian {
                texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
            },
        });
    }

    objects.push(Object {
        shape: Shape::Sphere {
            center: Point3([0.0, -1000.0, 0.0]),
            radius: 1000.0,
        },

        material: Material::Lambertian {
            texture: Texture::Checkered {
                even: Box::new(Texture::SolidColor(Vec3([0.2, 0.3, 0.1]))),
                odd: Box::new(Texture::SolidColor(Vec3([0.9, 0.9, 0.9]))),
            },
        },
    });

    let scene = Scene::new(objects, Vec3([1.0, 1.0, 1.0]), &mut rng);

    let image: RgbImage = ImageBuffer::from_raw(
        IMAGE_WIDTH,
        IMAGE_HEIGHT,
        pixels(&camera, &scene, image_options),
    )
    .unwrap();

    image.save("test.png").unwrap();
}
