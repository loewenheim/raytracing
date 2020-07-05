use image::{ImageBuffer, RgbImage};
use raytracing::camera::{Camera, CameraOptions};
use raytracing::geometry::{Axes, Point3, Shape, Vec3};
use raytracing::materials::Material;
use raytracing::textures::Texture;
use raytracing::{pixels, ImageOptions, Object, World};

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
        focus_distance: 1.0,
        looking_at: Point3([278.0, 278.0, 0.0]),
        origin: Point3([278.0, 278.0, -800.0]),
        shutter_open: 0.0,
        shutter_close: 0.01,
        vfov: 40.0,
        vup: Vec3([0.0, 1.0, 0.0]),
    });

    let mut objects = Vec::new();

    objects.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::YZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 555.0,
        }
        .flipped(),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.12, 0.45, 0.15])),
        },
    });

    objects.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::YZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 0.0,
        },
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.65, 0.05, 0.05])),
        },
    });

    objects.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XZ,
            lower_left: (213.0, 227.0),
            upper_right: (343.0, 332.0),
            height: 554.0,
        },
        material: Material::DiffuseLight {
            emit: Texture::SolidColor(Vec3([15.0, 15.0, 15.0])),
        },
    });

    objects.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 0.0,
        }
        .flipped(),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    });

    objects.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 555.0,
        },
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    });

    objects.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XY,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 555.0,
        }
        .flipped(),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
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
