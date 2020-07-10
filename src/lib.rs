pub mod geometry;
pub mod materials;
pub mod textures;

use camera::Camera;
use geometry::{BoundingBox, Geometry, IntersectionPoint, Ray, Vec3};
use indicatif::ProgressBar;
use materials::Material;
use rand::rngs::ThreadRng;
use rand::Rng;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::sync::{Arc, Mutex};

pub trait Optics: Geometry {
    fn scatter(
        &self,
        ray: &Ray,
        tmin: f64,
        tmax: f64,
        time: f64,
        rng: &mut ThreadRng,
    ) -> Option<RayHit> {
        self.scatter_(ray, tmin, tmax, time, rng).map(|(_, rh)| rh)
    }

    fn scatter_(
        &self,
        ray: &Ray,
        tmin: f64,
        tmax: f64,
        time: f64,
        rng: &mut ThreadRng,
    ) -> Option<(f64, RayHit)>;
}

pub fn pixels(
    camera: &Camera,
    world: &World,
    image_options: ImageOptions,
    progress_bar: Arc<Mutex<ProgressBar>>,
) -> Vec<u8> {
    (0..image_options.height)
        .into_par_iter()
        .rev()
        .flat_map(move |j| {
            let progress_bar = Arc::clone(&progress_bar);
            (0..image_options.width).into_par_iter().map(move |i| {
                let progress_bar = Arc::clone(&progress_bar);
                pixel(&camera, &world, (i, j), image_options, progress_bar)
            })
        })
        .collect::<Vec<[u8; 3]>>()
        .concat()
}

fn pixel<'a>(
    camera: &Camera,
    world: &World,
    (row, column): (u32, u32),
    ImageOptions {
        width,
        height,
        samples_per_pixel,
        max_reflections,
    }: ImageOptions,
    progress_bar: Arc<Mutex<ProgressBar>>,
) -> [u8; 3] {
    let v: Vec3 = (0..samples_per_pixel)
        .into_par_iter()
        .map(|_| {
            let mut rng = rand::thread_rng();
            let u = (f64::from(row) + rng.gen::<f64>()) / f64::from(width - 1);
            let v = (f64::from(column) + rng.gen::<f64>()) / f64::from(height - 1);
            ray_color(
                &camera.ray(u, v, &mut rng),
                world,
                camera.random_time(&mut rng),
                &mut rng,
                max_reflections,
            )
        })
        .sum::<Vec3>()
        / samples_per_pixel as f64;

    let progress_bar = progress_bar.lock().unwrap();

    progress_bar.inc(1);
    [
        (v[0].max(0.0).min(1.0) * 255.999) as u8,
        (v[1].max(0.0).min(1.0) * 255.999) as u8,
        (v[2].max(0.0).min(1.0) * 255.999) as u8,
    ]
}

fn ray_color(
    ray: &Ray,
    world: &World,
    time: f64,
    rng: &mut ThreadRng,
    max_reflections: usize,
) -> Vec3 {
    let mut ray = *ray;
    let mut result = Vec3([1.0, 1.0, 1.0]);

    for _ in 0..max_reflections {
        match world.objects.scatter(&ray, 0.001, f64::INFINITY, time, rng) {
            None => return result * world.background_color,
            Some(RayHit::Emitted(color)) => return result * color,
            Some(RayHit::Scattered {
                attenuation,
                ray: new_ray,
            }) => {
                result *= attenuation;
                ray = new_ray;
            }
        }
    }

    Vec3([0.0, 0.0, 0.0])
}

#[derive(Debug, Clone, Copy)]
pub struct ImageOptions {
    pub height: u32,
    pub width: u32,
    pub samples_per_pixel: usize,
    pub max_reflections: usize,
}

pub struct World {
    pub objects: Box<dyn Optics>,
    pub background_color: Vec3,
}

#[derive(Debug)]
pub enum BvhNode {
    Leaf {
        bounding_box: Option<BoundingBox>,
        object: Box<dyn Optics>,
    },

    Branch {
        bounding_box: Option<BoundingBox>,
        left: Box<BvhNode>,
        right: Box<BvhNode>,
    },
}

impl BvhNode {
    pub fn create<R: Rng + ?Sized>(objects: Vec<Box<dyn Optics>>, rng: &mut R) -> Self {
        let mut objects: Vec<Option<Box<dyn Optics>>> = objects.into_iter().map(Some).collect();
        Self::create_(&mut objects, rng)
    }

    fn create_<R: Rng + ?Sized>(objects: &mut [Option<Box<dyn Optics>>], rng: &mut R) -> Self {
        assert!(!objects.is_empty());
        let n = objects.len();

        if n == 1 {
            let object = objects[0].take().unwrap();
            let bounding_box = object.bound();
            Self::Leaf {
                object,
                bounding_box,
            }
        } else {
            let axis: usize = rng.gen_range(0, 3);

            objects.sort_by(|o1, o2| {
                match (o1.as_ref().unwrap().bound(), o2.as_ref().unwrap().bound()) {
                    (Some(bb1), Some(bb2)) => bb1.min[axis]
                        .partial_cmp(&bb2.min[axis])
                        .unwrap_or(Ordering::Equal),
                    (None, Some(_)) => Ordering::Less,
                    (Some(_), None) => Ordering::Greater,
                    (None, None) => Ordering::Equal,
                }
            });
            let (objects_left, objects_right) = objects.split_at_mut(n / 2);
            let left = Self::create_(objects_left, rng);
            let right = Self::create_(objects_right, rng);
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

    fn bounding_box(&self) -> Option<BoundingBox> {
        match *self {
            Self::Leaf { bounding_box, .. } => bounding_box,
            Self::Branch { bounding_box, .. } => bounding_box,
        }
    }
}

impl Optics for BvhNode {
    fn scatter_(
        &self,
        ray: &Ray,
        tmin: f64,
        mut tmax: f64,
        time: f64,
        rng: &mut ThreadRng,
    ) -> Option<(f64, RayHit)> {
        if let Some(bb) = self.bounding_box() {
            if !bb.hit(ray, tmin, tmax) {
                return None;
            }
        }
        match self {
            Self::Leaf { object, .. } => (*object).scatter_(ray, tmin, tmax, time, rng),
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

impl Geometry for BvhNode {
    fn intersect(
        &self,
        ray: &Ray,
        tmin: f64,
        mut tmax: f64,
        time: f64,
    ) -> Option<IntersectionPoint> {
        if let Some(bb) = self.bounding_box() {
            if !bb.hit(ray, tmin, tmax) {
                return None;
            }
        }

        match self {
            Self::Leaf { object, .. } => object.intersect(ray, tmin, tmax, time),
            Self::Branch { left, right, .. } => {
                let mut ip = None;

                if let Some(res) = (*left).intersect(ray, tmin, tmax, time) {
                    ip = Some(res);
                    tmax = res.t;
                }

                if let Some(res) = (*right).intersect(ray, tmin, tmax, time) {
                    ip = Some(res);
                }

                ip
            }
        }
    }

    fn bound(&self) -> Option<BoundingBox> {
        self.bounding_box()
    }
}

#[derive(Debug)]
pub struct Object {
    pub shape: Box<dyn Geometry + Sync + Send>,
    pub material: Material,
}

impl Optics for Object {
    fn scatter_(
        &self,
        r: &Ray,
        tmin: f64,
        tmax: f64,
        time: f64,
        rng: &mut ThreadRng,
    ) -> Option<(f64, RayHit)> {
        self.shape
            .intersect(r, tmin, tmax, time)
            .as_ref()
            .map(|p| (p.t, self.material.scatter(p, rng)))
    }
}

impl Geometry for Object {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64, time: f64) -> Option<IntersectionPoint> {
        self.shape.intersect(ray, tmin, tmax, time)
    }
    fn bound(&self) -> Option<BoundingBox> {
        self.shape.bound()
    }
}

pub enum RayHit {
    Scattered { ray: Ray, attenuation: Vec3 },
    Emitted(Vec3),
}

pub mod camera {
    use super::geometry::*;
    use rand::distributions::Distribution;
    use rand::Rng;

    #[derive(Debug, Clone, Copy)]
    pub struct Camera {
        aspect_ratio: f64,
        horizontal: Vec3,
        lens_radius: f64,
        lower_left_corner: Point3,
        onb: Onb,
        origin: Point3,
        shutter_close: f64,
        shutter_open: f64,
        vertical: Vec3,
        vfov: f64,
    }

    impl Camera {
        pub fn new(
            CameraOptions {
                aperture,
                aspect_ratio,
                focus_distance,
                looking_at,
                origin,
                shutter_close,
                shutter_open,
                vfov,
                vup,
            }: CameraOptions,
        ) -> Self {
            let h = (vfov / 2.0).to_radians().tan();
            let viewport_height = 2.0 * h;
            let viewport_width = viewport_height * aspect_ratio;

            let onb = Onb::from23(vup, origin - looking_at);

            let horizontal = *onb[0] * viewport_width * focus_distance;
            let vertical = *onb[1] * viewport_height * focus_distance;
            let lower_left_corner =
                origin - horizontal / 2.0 - vertical / 2.0 - *onb[2] * focus_distance;
            let lens_radius = aperture / 2.0;

            Self {
                aspect_ratio,
                horizontal,
                lens_radius,
                lower_left_corner,
                onb,
                origin,
                shutter_close,
                shutter_open,
                vertical,
                vfov,
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

        pub fn random_time<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
            rng.gen_range(self.shutter_open, self.shutter_close)
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CameraOptions {
        pub aperture: f64,
        pub aspect_ratio: f64,
        pub focus_distance: f64,
        pub looking_at: Point3,
        pub origin: Point3,
        pub shutter_close: f64,
        pub shutter_open: f64,
        pub vfov: f64,
        pub vup: Vec3,
    }
}
