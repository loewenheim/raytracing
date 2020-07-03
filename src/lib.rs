pub mod geometry;
pub mod materials;
pub mod textures;

use camera::{Camera, CameraOptions};
use geometry::{Axes, Boundable, BoundingBox, Intersection, Point3, Ray, Shape, Vec3};
use materials::Material;
use perlin_noise::PerlinNoise;
use rand::Rng;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::sync::Arc;
use textures::Texture;

pub fn color(v: Vec3) -> [u8; 3] {
    [
        (v[0].max(0.0).min(1.0) * 255.999) as u8,
        (v[1].max(0.0).min(1.0) * 255.999) as u8,
        (v[2].max(0.0).min(1.0) * 255.999) as u8,
    ]
}
#[derive(Debug, Clone)]
pub enum BvhNode<'a, T: Boundable> {
    Leaf {
        bounding_box: Option<BoundingBox>,
        shape: &'a T,
    },

    Branch {
        bounding_box: Option<BoundingBox>,
        left: Box<BvhNode<'a, T>>,
        right: Box<BvhNode<'a, T>>,
    },
}

impl<'a, T: Boundable> BvhNode<'a, T> {
    pub fn create<R: Rng + ?Sized>(objects: &'a mut [T], rng: &mut R) -> Self {
        assert!(!objects.is_empty());
        let n = objects.len();

        if n == 1 {
            let shape = &objects[0];
            let bounding_box = shape.bound();
            Self::Leaf {
                shape,
                bounding_box,
            }
        } else {
            let axis: usize = rng.gen_range(0, 3);

            objects.sort_by(|o1, o2| match (o1.bound(), o2.bound()) {
                (Some(bb1), Some(bb2)) => bb1.min[axis]
                    .partial_cmp(&bb2.min[axis])
                    .unwrap_or(Ordering::Equal),
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
                (None, None) => Ordering::Equal,
            });
            let (objects_left, objects_right) = objects.split_at_mut(n / 2);
            let left = Self::create(objects_left, rng);
            let right = Self::create(objects_right, rng);
            let bounding_box =
                if let (Some(lb), Some(rb)) = (left.bounding_box(), right.bounding_box()) {
                    Some(lb + rb)
                } else {
                    None
                };

            Self::Branch {
                left: Box::new(left),
                right: Box::new(right),
                bounding_box,
            }
        }
    }

    pub fn bounding_box(&self) -> Option<BoundingBox> {
        match *self {
            Self::Leaf { bounding_box, .. } => bounding_box,
            Self::Branch { bounding_box, .. } => bounding_box,
        }
    }

    pub fn bounding_boxes(&self) -> Vec<Option<BoundingBox>> {
        match self {
            Self::Leaf { bounding_box, .. } => vec![*bounding_box],
            Self::Branch {
                bounding_box,
                left,
                right,
            } => {
                let mut l = left.bounding_boxes();
                let mut r = right.bounding_boxes();
                let mut result = vec![*bounding_box];
                result.append(&mut l);
                result.append(&mut r);
                result
            }
        }
    }
}

impl<'a> BvhNode<'a, Object> {
    pub fn scatter<R: Rng + ?Sized>(
        &self,
        ray: &Ray,
        tmin: f64,
        tmax: f64,
        time: f64,
        rng: &mut R,
    ) -> Option<RayHit> {
        self.scatter_(ray, tmin, tmax, time, rng).map(|x| x.1)
    }

    fn scatter_<R: Rng + ?Sized>(
        &self,
        ray: &Ray,
        tmin: f64,
        mut tmax: f64,
        time: f64,
        rng: &mut R,
    ) -> Option<(f64, RayHit)> {
        if let Some(bb) = self.bounding_box() {
            if !bb.hit(ray, tmin, tmax) {
                return None;
            }
        }
        match self {
            Self::Leaf { shape, .. } => (*shape).scatter(ray, tmin, tmax, time, rng),
            Self::Branch { left, right, .. } => {
                let mut scattered = None;
                if let Some((t, new_scattered)) = (*left).scatter_(ray, tmin, tmax, time, rng) {
                    scattered = Some((t, new_scattered));
                    tmax = t;
                }
                if let Some((t, new_scattered)) = (*right).scatter_(ray, tmin, tmax, time, rng) {
                    scattered = Some((t, new_scattered));
                }

                scattered
            }
        }
    }
}

pub fn default_world<R: Rng + ?Sized>(rng: &mut R) -> Vec<Object> {
    let mut world = random_world(rng);
    let sphere1 = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 1.0, 0.0]),
            radius: 1.0,
        },

        material: Material::glass(),
    };

    //let sphere2 = Object {
    //    shape: Shape::Sphere {
    //        center: Point3([-4.0, 1.0, 0.0]),
    //        radius: 1.0,
    //    },

    //    material: Material::Lambertian {
    //        albedo: Texture::SolidColor(Color::new(0.4, 0.2, 0.1)),
    //    },
    //};
    let sphere2 = earth(Point3([4.0, 1.0, 0.0]), 1.0);

    let sphere3 = Object {
        shape: Shape::Sphere {
            center: Point3([-4.0, 1.0, 0.0]),
            radius: 1.0,
        },

        material: Material::Metal {
            albedo: Vec3([0.7, 0.6, 0.5]),
            fuzz: 0.0,
        },
    };

    world.push(sphere1);
    world.push(sphere2);
    world.push(sphere3);

    world
}

pub fn two_perlin_spheres<R: Rng + ?Sized>(rng: &mut R) -> Vec<Object> {
    let mut world = Vec::new();

    let perlin = Arc::new(PerlinNoise::new());
    let sphere1 = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, -1000.0, 0.0]),
            radius: 1000.0,
        },
        material: Material::Lambertian {
            albedo: Texture::Noise(perlin.clone()),
        },
    };

    let sphere2 = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 2.0, 0.0]),
            radius: 2.0,
        },
        material: Material::Lambertian {
            albedo: Texture::Noise(perlin.clone()),
        },
    };

    let diff_light = Material::DiffuseLight {
        emit: Texture::SolidColor(Vec3([4.0, 4.0, 4.0])),
    };

    let light_sphere = Object {
        shape: Shape::Sphere {
            center: Point3([0.0, 7.0, 0.0]),
            radius: 2.0,
        },

        material: diff_light.clone(),
    };

    let light_rectangle = Object {
        shape: Shape::Rectangle {
            lower_left: (3.0, 1.0),
            upper_right: (5.0, 3.0),
            height: -2.0,
            axes: Axes::XY,
        },

        material: diff_light,
    };

    world.push(sphere1);
    world.push(sphere2);
    world.push(light_sphere);
    world.push(light_rectangle);

    world
}

pub fn cornell_box(aspect_ratio: f64) -> (Vec<Object>, Camera) {
    let mut world = Vec::new();

    world.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::YZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 555.0,
        }
        .flipped(),
        material: Material::Lambertian {
            albedo: Texture::SolidColor(Vec3([0.12, 0.45, 0.15])),
        },
    });

    world.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::YZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 0.0,
        },
        material: Material::Lambertian {
            albedo: Texture::SolidColor(Vec3([0.65, 0.05, 0.05])),
        },
    });

    world.push(Object {
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

    world.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 0.0,
        }
        .flipped(),
        material: Material::Lambertian {
            albedo: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    });

    world.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XZ,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 555.0,
        },
        material: Material::Lambertian {
            albedo: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    });

    world.push(Object {
        shape: Shape::Rectangle {
            axes: Axes::XY,
            lower_left: (0.0, 0.0),
            upper_right: (555.0, 555.0),
            height: 555.0,
        }
        .flipped(),
        material: Material::Lambertian {
            albedo: Texture::SolidColor(Vec3([0.73, 0.73, 0.73])),
        },
    });

    let look_from = Point3([278.0, 278.0, -800.0]);
    let look_at = Point3([278.0, 278.0, 0.0]);
    let vup = Vec3([0.0, 1.0, 0.0]);
    let focus_distance = 1.0;
    let aperture = 0.0;
    let vfov = 40.0;

    (
        world,
        Camera::new(CameraOptions {
            origin: look_from,
            direction: look_at - look_from,
            vup,
            focus_distance,
            aperture,
            vfov,
            aspect_ratio,
        }),
    )
}

pub fn earth(center: Point3, radius: f64) -> Object {
    let texture = Texture::image("earthmap.jpg");
    Object {
        shape: Shape::Sphere { center, radius },

        material: Material::Lambertian { albedo: texture },
    }
}

pub fn random_world<R: Rng + ?Sized>(rng: &mut R) -> Vec<Object> {
    let mut world = Vec::new();

    let ground = Object {
        //shape: Shape::Plane {
        //    normal: Vec3([0.0, 1.0, 0.0]).normed(),
        //    offset: 0.0,
        //},
        shape: Shape::Sphere {
            center: Point3([0.0, -1000.0, 0.0]),
            radius: 1000.0,
        },

        material: Material::Lambertian {
            albedo: Texture::Checkered {
                even: Box::new(Texture::SolidColor(Vec3([0.2, 0.3, 0.1]))),
                odd: Box::new(Texture::SolidColor(Vec3([0.9, 0.9, 0.9]))),
            },
        },
    };

    world.push(ground);

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
                    let albedo = Texture::SolidColor(
                        Vec3::random(0.0..1.0, rng) + Vec3::random(0.0..1.0, rng),
                    );
                    (center2, Material::Lambertian { albedo })
                } else if roll < 0.95 {
                    let albedo = Vec3::random(0.5..1.0, rng);
                    let fuzz = rng.gen_range(0.0, 0.5);

                    (center, Material::Metal { albedo, fuzz })
                } else {
                    (center, Material::glass())
                };

                let shape = Shape::sphere(0.2, (0.0, center), (1.0, center2));
                world.push(Object { shape, material })
            }
        }
    }

    world
}

#[derive(Clone)]
pub struct Object {
    pub shape: Shape,
    pub material: Material,
}

impl Object {
    pub fn scatter<R: Rng + ?Sized>(
        &self,
        r: &Ray,
        tmin: f64,
        tmax: f64,
        time: f64,
        rng: &mut R,
    ) -> Option<(f64, RayHit)> {
        self.shape
            .intersect(r, tmin, tmax, time)
            .as_ref()
            .map(|p| (p.t, self.material.scatter(p, rng)))
    }
}

impl Boundable for Object {
    fn bound(&self) -> Option<BoundingBox> {
        self.shape.bound()
    }
}

pub enum RayHit {
    Scattered { ray: Ray, attenuation: Vec3 },

    Emitted(Vec3),
}

pub fn pixel<'a>(
    camera: &Camera,
    world: &BvhNode<Object>,
    background_color: Vec3,
    (row, column): (u32, u32),
    (width, height): (u32, u32),
    (open_time, close_time): (f64, f64),
    samples: usize,
    max_depth: usize,
) -> [u8; 3] {
    let color_vec: Vec3 = (0..samples)
        .into_par_iter()
        .map(|_| {
            let mut rng = rand::thread_rng();
            let u = (f64::from(row) + rng.gen::<f64>()) / f64::from(width - 1);
            let v = (f64::from(column) + rng.gen::<f64>()) / f64::from(height - 1);
            let time = rng.gen_range(open_time, close_time);
            ray_color(
                &camera.ray(u, v, &mut rng),
                &world,
                background_color,
                time,
                &mut rng,
                max_depth,
            )
        })
        .sum::<Vec3>()
        / samples as f64;
    color(color_vec)
}

pub fn ray_color<R>(
    ray: &Ray,
    world: &BvhNode<Object>,
    background_color: Vec3,
    time: f64,
    rng: &mut R,
    depth: usize,
) -> Vec3
where
    R: Rng + ?Sized,
{
    if depth == 0 {
        return Vec3([0.0, 0.0, 0.0]);
    }
    match world.scatter(ray, 0.001, f64::INFINITY, time, rng) {
        None => background_color,

        Some(RayHit::Emitted(color)) => color,
        Some(RayHit::Scattered { attenuation, ray }) => {
            attenuation * ray_color(&ray, world, background_color, time, rng, depth - 1)
        }
    }
}

pub mod camera {
    use super::geometry::*;
    use rand::distributions::Distribution;
    use rand::Rng;

    #[derive(Debug, Clone, Copy)]
    pub struct Camera {
        origin: Point3,
        horizontal: Vec3,
        vertical: Vec3,
        lower_left_corner: Point3,
        vfov: f64,
        aspect_ratio: f64,
        onb: Onb,
        lens_radius: f64,
    }

    impl Camera {
        pub fn new(
            CameraOptions {
                origin,
                aperture,
                direction,
                vfov,
                focus_distance,
                aspect_ratio,
                vup,
            }: CameraOptions,
        ) -> Self {
            let h = (vfov / 2.0).to_radians().tan();
            let viewport_height = 2.0 * h;
            let viewport_width = viewport_height * aspect_ratio;

            let onb = Onb::from23(vup, -direction);

            let horizontal = *onb[0] * viewport_width * focus_distance;
            let vertical = *onb[1] * viewport_height * focus_distance;
            let lower_left_corner =
                origin - horizontal / 2.0 - vertical / 2.0 - *onb[2] * focus_distance;
            let lens_radius = aperture / 2.0;

            Self {
                origin,
                horizontal,
                vertical,
                lower_left_corner,
                vfov,
                aspect_ratio,
                onb,
                lens_radius,
            }
        }

        pub fn ray<R: Rng + ?Sized>(&self, s: f64, t: f64, rng: &mut R) -> Ray {
            let rd = UnitDisc.sample(rng) * self.lens_radius;
            let offset = *self.onb[0] * rd[0] + *self.onb[1] * rd[1];
            let origin = self.origin + offset;
            let direction =
                (self.lower_left_corner + self.horizontal * s + self.vertical * t) - origin;
            Ray { origin, direction }
        }
    }

    pub struct CameraOptions {
        pub origin: Point3,
        pub aperture: f64,
        pub focus_distance: f64,
        pub aspect_ratio: f64,
        pub direction: Vec3,
        pub vfov: f64,
        pub vup: Vec3,
    }
}
