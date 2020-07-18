pub mod geometry;
pub mod materials;
pub mod textures;

use camera::Camera;
use geometry::{Boundable, BoundingBox, Intersection, Ray, Shape, Vec3};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use materials::Material;
use rand::Rng;
use rayon::prelude::*;
use std::cmp::Ordering;

/// Computes a vector of byte values representing an image (3 bytes per pixel).
/// # Arguments
/// * `camera` - The camera the scene is viewed from
/// * `scene` - The scene being viewed
/// * `image_options` - Settings for the generated image
pub fn pixels(camera: &Camera, scene: &Scene, image_options: ImageOptions) -> Vec<u8> {
    let bar = ProgressBar::new(image_options.height as _);
    bar.set_style(ProgressStyle::default_bar().template("Lines: {wide_bar} {percent:2} %"));
    (0..image_options.height)
        .into_par_iter()
        .rev()
        .progress_with(bar)
        .flat_map(move |j| {
            (0..image_options.width)
                .into_par_iter()
                .map(move |i| pixel(&camera, &scene, (i, j), image_options))
        })
        .collect::<Vec<[u8; 3]>>()
        .concat()
}

/// Computes RGB byte values for one image pixel.
/// # Arguments
/// * `camera` - The camera the scene is viewed from
/// * `scene` - The scene being viewed
/// * `(row, column)` - The coordinates of the pixel
/// * `image_options` - Settings for the generated image
fn pixel<'a>(
    camera: &Camera,
    scene: &Scene,
    (row, column): (u32, u32),
    ImageOptions {
        width,
        height,
        samples_per_pixel,
        max_reflections,
    }: ImageOptions,
) -> [u8; 3] {
    let v: Vec3 = (0..samples_per_pixel)
        .into_par_iter()
        .map(|_| {
            let mut rng = rand::thread_rng();
            let u = (f64::from(row) + rng.gen::<f64>()) / f64::from(width - 1);
            let v = (f64::from(column) + rng.gen::<f64>()) / f64::from(height - 1);
            ray_color(
                &camera.ray(u, v, &mut rng),
                scene,
                camera.random_time(&mut rng),
                &mut rng,
                max_reflections,
            )
        })
        .sum::<Vec3>()
        / samples_per_pixel as f64;

    [
        (v[0].max(0.0).min(1.0).sqrt() * 255.999) as u8,
        (v[1].max(0.0).min(1.0).sqrt() * 255.999) as u8,
        (v[2].max(0.0).min(1.0).sqrt() * 255.999) as u8,
    ]
}

/// Computes the color of a ray.
/// # Arguments
/// * `ray` - The ray to color
/// * `scene` - The scene being viewed
/// * `time` - The time at which the ray is emitted
/// * `rng` - A random number generator
/// * `max_reflections` - The maximum number of time a ray can be reflected or scattered
fn ray_color<R>(ray: &Ray, scene: &Scene, time: f64, rng: &mut R, max_reflections: usize) -> Vec3
where
    R: Rng + ?Sized,
{
    let mut ray = *ray;
    let mut result = Vec3([1.0, 1.0, 1.0]);

    for _ in 0..max_reflections {
        match scene.objects.scatter(&ray, 0.001, f64::INFINITY, time, rng) {
            None => return result * scene.background_color,
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

/// Settings for generating an image.
#[derive(Debug, Clone, Copy)]
pub struct ImageOptions {
    /// The image's height
    pub height: u32,
    /// The images width
    pub width: u32,
    /// How many rays to sample and average for every pixel
    pub samples_per_pixel: usize,
    /// The maximum number of times a ray can be reflected or scattered
    pub max_reflections: usize,
}

#[derive(Clone, Debug)]
pub struct Scene {
    objects: BvhNode<Object>,
    background_color: Vec3,
}

impl Scene {
    /// Creates a new scene. The vector of objects is converted to a BVH tree
    pub fn new<R: Rng + ?Sized>(objects: Vec<Object>, background_color: Vec3, rng: &mut R) -> Self {
        Self {
            objects: BvhNode::create(objects, rng),
            background_color,
        }
    }
}

/// A Node in a BVH (bounded volume hierarchy) tree
#[derive(Debug, Clone)]
enum BvhNode<T: Boundable> {
    /// A leaf, containing an object and possibly a bounding box
    Leaf {
        bounding_box: Option<BoundingBox>,
        object: T,
    },

    /// A branch, containing a left and right [`BVHNode`] and possibly
    /// a bounding box
    Branch {
        bounding_box: Option<BoundingBox>,
        left: Box<BvhNode<T>>,
        right: Box<BvhNode<T>>,
    },
}

impl<T: Boundable> BvhNode<T> {
    fn create<R: Rng + ?Sized>(objects: Vec<T>, rng: &mut R) -> Self {
        Self::create_(&mut objects.into_iter().map(Some).collect::<Vec<_>>(), rng)
    }

    /// Creates a BVH tree from a slice of objects by randomly choosing an axis
    /// splitting the slice in half along that axis, and recursing
    fn create_<R: Rng + ?Sized>(objects: &mut [Option<T>], rng: &mut R) -> Self {
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

impl BvhNode<Object> {
    /// Computes the scattering of a ray hitting the BVH tree.
    /// We first check whether the tree's bounding box is hit.
    /// If so, we recurse further into the tree.
    /// Only when we arrive at a leaf do we actually check whether
    /// an object is hit.
    fn scatter<R: Rng + ?Sized>(
        &self,
        ray: &Ray,
        tmin: f64,
        tmax: f64,
        time: f64,
        rng: &mut R,
    ) -> Option<RayHit> {
        self.scatter_(ray, tmin, tmax, time, rng).map(|x| x.1)
    }

    #[doc(hidden)]
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
            Self::Leaf { object, .. } => (*object).scatter(ray, tmin, tmax, time, rng),
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

/// A combination of geometry and material
#[derive(Clone, Debug)]
pub struct Object {
    /// The physical shape of the object
    pub shape: Shape,
    /// The material of the object
    pub material: Material,
}

impl Object {
    /// Computes the result of a ray hitting the object and
    /// interacting with the material
    ///
    /// # Arguments
    /// * `r` - The ray
    /// * `time` - The time at which the ray is emitted
    /// * `rng` - A random number generator
    fn scatter<R: Rng + ?Sized>(
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

/// The result of a ray hitting an [`Object`]
///
/// [`Object`]: ./struct.Object.html
pub enum RayHit {
    /// The ray was scattered, resulting in a new ray and
    /// an attenuation value
    Scattered { ray: Ray, attenuation: Vec3 },
    /// The ray hit a light source of a given color
    Emitted(Vec3),
}

pub mod camera {
    use super::geometry::*;
    use rand::distributions::Distribution;
    use rand::Rng;

    /// A representation of a camera
    #[derive(Debug, Clone, Copy)]
    pub struct Camera {
        /// The aspect ratio of the camera's viewport
        aspect_ratio: f64,
        /// The vector from the lower left corner of the
        /// viewport to the lower right corner
        horizontal: Vec3,
        /// The radius of the lens
        lens_radius: f64,
        /// Coordinates of the lower left corner of the viewport
        lower_left_corner: Point3,
        /// An orthonormal base consisting of
        /// * the camera's "right" direction,
        /// * the camera's "up" direction,
        /// * the camera's "back" (away from the viewport) direction
        onb: Onb,
        /// The point the camera is situated at
        origin: Point3,
        /// The time at which the shutter closes
        shutter_close: f64,
        /// The time at which the shutter opens
        shutter_open: f64,
        /// The vector from the lower left corner of the
        /// viewport to the upper left corner
        vertical: Vec3,
        /// The camera's vertical field of view in degrees
        vfov: f64,
    }

    impl Camera {
        /// Creates a new camera from a [`CameraOptions`] struct
        ///
        /// [`CameraOptions`]: ./struct.CameraOptions.html
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

        /// Computes the ray from the camera's origin to the point (s, t)
        /// in viewport coordinates. The ray is randomized slightly to allow
        /// for random sampling
        pub(crate) fn ray<R: Rng + ?Sized>(&self, s: f64, t: f64, rng: &mut R) -> Ray {
            let rd = InsideUnitDisc.sample(rng) * self.lens_radius;
            let offset = *self.onb[0] * rd[0] + *self.onb[1] * rd[1];
            let origin = self.origin + offset;
            let direction =
                (self.lower_left_corner + self.horizontal * s + self.vertical * t) - origin;
            Ray { origin, direction }
        }

        /// A random time between the shutter opening and closing
        pub(crate) fn random_time<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
            rng.gen_range(self.shutter_open, self.shutter_close)
        }
    }

    /// A struct capturing the options needed to create a camera
    #[derive(Debug, Clone, Copy)]
    pub struct CameraOptions {
        /// The radius of the camera's aperture
        pub aperture: f64,
        /// The aspect ratio of the camera's viewport
        pub aspect_ratio: f64,
        /// The distance at which things are in focus
        pub focus_distance: f64,
        /// A point the camera is looking straight at
        pub looking_at: Point3,
        /// The point the camera is situated at
        pub origin: Point3,
        /// The time at which the shutter closes
        pub shutter_close: f64,
        /// The time at which the shutter opens
        pub shutter_open: f64,
        /// The camera's vertical field of view in degrees
        pub vfov: f64,
        /// A vector describing the camera's "up" direction.
        pub vup: Vec3,
    }
}
