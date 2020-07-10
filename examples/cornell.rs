use image::{ImageBuffer, RgbImage};
use indicatif::{ProgressBar, ProgressStyle};
use raytracing::camera::{Camera, CameraOptions};
use raytracing::geometry::{Axis, Geometry, MyBox, Point3, Rotate, Shape, Vec3};
use raytracing::materials::Material;
use raytracing::textures::Texture;
use raytracing::{pixels, BvhNode, ImageOptions, Object, Optics, World};
use std::sync::{Arc, Mutex};

fn main() {
    const ASPECT_RATIO: f64 = 1.0;
    const IMAGE_WIDTH: u32 = 750;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: usize = 100;
    const MAX_REFLECTIONS: usize = 50;

    let progress_bar = ProgressBar::new((IMAGE_WIDTH * IMAGE_HEIGHT).into());
    progress_bar.set_style(ProgressStyle::default_bar().template("{wide_bar} {percent:2}%"));
    let progress_bar = Arc::new(Mutex::new(progress_bar));

    let mut rng = rand::thread_rng();

    let image_options = ImageOptions {
        width: IMAGE_WIDTH,
        height: IMAGE_HEIGHT,
        samples_per_pixel: SAMPLES_PER_PIXEL,
        max_reflections: MAX_REFLECTIONS,
    };

    let camera = Camera::new(CameraOptions {
        aperture: 0.0,
        aspect_ratio: ASPECT_RATIO,
        focus_distance: 10.0,
        looking_at: Point3([278.0, 278.0, 0.0]),
        origin: Point3([278.0, 278.0, -800.0]),
        shutter_open: 0.0,
        shutter_close: 1.0,
        vfov: 40.0,
        vup: Vec3([0.0, 1.0, 0.0]),
    });

    let mut objects: Vec<Box<dyn Optics>> = Vec::new();

    // left wall
    objects.push(Box::new(Object {
        shape: Box::new(Shape::rectangle(
            Point3([0.0, 0.0, 0.0]),
            Point3([0.0, 555.0, 555.0]),
        )),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.65, 0.05, 0.05])),
        },
    }));

    // right wall
    objects.push(Box::new(Object {
        shape: Box::new(
            Shape::rectangle(Point3([555.0, 0.0, 0.0]), Point3([555.0, 555.0, 555.0])).flip(),
        ),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.12, 0.45, 0.15])),
        },
    }));

    // floor
    objects.push(Box::new(Object {
        shape: Box::new(Shape::rectangle(
            Point3([0.0, 0.0, 0.0]),
            Point3([555.0, 0.0, 555.0]),
        )),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    }));

    // ceiling
    objects.push(Box::new(Object {
        shape: Box::new(
            Shape::rectangle(Point3([0.0, 555.0, 0.0]), Point3([555.0, 555.0, 555.0])).flip(),
        ),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    }));

    // back wall
    objects.push(Box::new(Object {
        shape: Box::new(
            Shape::rectangle(Point3([0.0, 0.0, 555.0]), Point3([555.0, 555.0, 555.0])).flip(),
        ),
        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    }));

    // first box
    objects.push(Box::new(Object {
        shape: Box::new(
            MyBox::new(Point3([0.0, 0.0, 0.0]), Point3([165.0, 330.0, 165.0]))
                .rotate(Axis::Y, 15.0)
                .translate(Vec3([265.0, 0.0, 295.0])),
        ),

        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    }));

    // second box
    objects.push(Box::new(Object {
        shape: Box::new(
            MyBox::new(Point3([0.0, 0.0, 0.0]), Point3([165.0, 165.0, 165.0]))
                .rotate(Axis::Y, -18.0)
                .translate(Vec3([130.0, 0.0, 65.0])),
        ),

        material: Material::Lambertian {
            texture: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    }));

    // light
    objects.push(Box::new(Object {
        shape: Box::new(Shape::rectangle(
            Point3([213.0, 554.0, 227.0]),
            Point3([343.0, 554.0, 332.0]),
        )),
        material: Material::DiffuseLight {
            emit: Texture::SolidColor(Vec3([15.0, 15.0, 15.0])),
        },
    }));

    let world = World {
        objects: Box::new(BvhNode::create(objects, &mut rng)),
        background_color: Vec3([0.0, 0.0, 0.0]),
    };
    let image: RgbImage = ImageBuffer::from_raw(
        IMAGE_WIDTH,
        IMAGE_HEIGHT,
        pixels(&camera, &world, image_options, Arc::clone(&progress_bar)),
    )
    .unwrap();

    image.save("test.png").unwrap();
}
