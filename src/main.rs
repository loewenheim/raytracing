use image::{ImageBuffer, RgbImage};
use raytracing::geometry::{Intersection, Point3, Ray, Sphere, Vec3};

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 384;
    const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;

    let viewport_height = 2.0;
    let viewport_width = ASPECT_RATIO * viewport_height;
    let focal_length = 1.0;

    let origin = Point3::default();
    let horizontal = Vec3([viewport_width, 0.0, 0.0]);
    let vertical = Vec3([0.0, viewport_height, 0.0]);
    let lower_left_corner =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3([0.0, 0.0, focal_length]);

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
        let u = f64::from(i) / f64::from(IMAGE_WIDTH - 1);
        let v = f64::from(IMAGE_HEIGHT - j) / f64::from(IMAGE_HEIGHT - 1);

        let ray = Ray {
            origin,
            direction: (lower_left_corner + horizontal * u + vertical * v) - origin,
        };

        ray.color(&world)
    });

    image.save("test.png").unwrap();
}
