use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::f64::consts::PI;
use std::iter::Sum;
use std::ops::{
    Add, AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

pub trait Intersection {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64, time: f64) -> Option<IntersectionPoint>;
}

pub trait Boundable {
    fn bound(&self) -> Option<BoundingBox>;
}

#[derive(Clone, Copy, Debug)]
pub enum Shape {
    Sphere {
        center: Point3,
        radius: f64,
    },

    Plane {
        normal: UnitVec3,
        offset: f64,
    },

    MovingSphere {
        ray: Ray,
        radius: f64,
        start_time: f64,
        end_time: f64,
    },
}

impl Shape {
    pub fn sphere(
        radius: f64,
        (start_time, start_center): (f64, Point3),
        (end_time, end_center): (f64, Point3),
    ) -> Self {
        assert!(end_time > start_time);

        if start_center == end_center {
            Self::Sphere {
                center: start_center,
                radius,
            }
        } else {
            Self::MovingSphere {
                start_time,
                end_time,
                radius,
                ray: Ray {
                    origin: start_center,
                    direction: end_center - start_center,
                },
            }
        }
    }
}

impl Intersection for Shape {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64, time: f64) -> Option<IntersectionPoint> {
        match *self {
            Shape::Sphere { center, radius } => {
                let Ray {
                    origin: o,
                    direction: dir,
                } = ray;
                let oc = *o - center;
                let a = dir.norm_squared();
                let half_b = oc.dot(*dir);
                let c = oc.norm_squared() - radius * radius;
                let discriminant = half_b * half_b - a * c;

                if discriminant < 0.0 {
                    None
                } else {
                    let t = (-half_b - discriminant.sqrt()) / a;
                    if t > tmin && t < tmax {
                        let point = ray.at(t);
                        let normal = (point - center).normed();
                        let (u, v) = {
                            let p = (point - center) / radius;
                            let phi = p[2].atan2(p[0]);
                            let theta = p[1].asin();
                            let u = 1.0 - (phi + PI) / (2.0 * PI);
                            let v = (theta + PI / 2.0) / PI;
                            (u, v)
                        };
                        Some(IntersectionPoint {
                            t,
                            point,
                            normal,
                            in_vec: dir.normed(),
                            surface_coordinates: (u, v),
                        })
                    } else {
                        None
                    }
                }
            }

            Shape::Plane { normal, offset } => {
                if ray.direction.dot(*normal).abs() <= 1e-10 {
                    None
                } else {
                    let u = Vec3(ray.origin.0);
                    let t = (offset - u.dot(*normal)) / ray.direction.dot(*normal);

                    if t >= tmin && t < tmax {
                        let intersection_point = IntersectionPoint {
                            t,
                            in_vec: ray.direction.normed(),
                            normal,
                            point: ray.at(t),
                            surface_coordinates: (0.0, 0.0),
                        };
                        Some(intersection_point)
                    } else {
                        None
                    }
                }
            }

            Shape::MovingSphere {
                ray: r,
                radius,
                start_time,
                ..
            } => {
                let sphere = Shape::Sphere {
                    center: r.at(time - start_time),
                    radius,
                };

                sphere.intersect(ray, tmin, tmax, time)
            }
        }
    }
}

impl Boundable for Shape {
    fn bound(&self) -> Option<BoundingBox> {
        match *self {
            Shape::Plane { .. } => None,
            Shape::Sphere { center, radius } => {
                let min = center - Vec3([1.0, 1.0, 1.0]) * radius;
                let max = center + Vec3([1.0, 1.0, 1.0]) * radius;

                Some(BoundingBox { min, max })
            }
            Shape::MovingSphere {
                ray,
                radius,
                start_time,
                end_time,
            } => Some(
                (&Shape::Sphere {
                    center: ray.at(0.0),
                    radius,
                })
                    .bound()
                    .unwrap()
                    + (&Shape::Sphere {
                        center: ray.at(end_time - start_time),
                        radius,
                    })
                        .bound()
                        .unwrap(),
            ),
        }
    }
}

impl<T: Boundable> Boundable for Vec<T> {
    fn bound(&self) -> Option<BoundingBox> {
        if self.is_empty() {
            return None;
        }

        let mut result = None;
        for shape in self.iter() {
            match (result, shape.bound()) {
                (Some(old_box), Some(new_box)) => result = Some(old_box + new_box),
                (None, Some(new_box)) => result = Some(new_box),
                (_, None) => return None,
            }
        }

        result
    }
}
/// Contains information about the intersection of a ray
/// with an object: the actual point, the normal vector
/// at that point, and the parameter of the ray.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    pub point: Point3,
    pub normal: UnitVec3,
    pub t: f64,
    pub in_vec: UnitVec3,
    pub surface_coordinates: (f64, f64),
}

impl IntersectionPoint {
    pub fn face(&self) -> Face {
        if self.in_vec.dot(*self.normal) < 0.0 {
            Face::Front
        } else {
            Face::Back
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Face {
    Front,
    Back,
}

/// A point in three-dimensional space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3(pub [f64; 3]);

impl Point3 {
    pub fn dist(&self, other: &Point3) -> f64 {
        (*self - *other).norm()
    }
}

impl Deref for Point3 {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Point3 {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Default for Point3 {
    fn default() -> Self {
        Self([0.0, 0.0, 0.0])
    }
}

/// A vector in three-dimensional space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3(pub [f64; 3]);

impl Vec3 {
    /// Returns the (Euclidean) norm of the vector.
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// Returns the square of the (Euclidean) norm of the vector.
    pub fn norm_squared(&self) -> f64 {
        self[0].powi(2) + self[1].powi(2) + self[2].powi(2)
    }

    /// Returns the dot (scalar) product of this vector with another.
    pub fn dot(self, other: Self) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }

    /// Returns the cross (vector) product of this vector with another.
    pub fn cross(self, other: Self) -> Self {
        Self([
            self[1] * other[2] - self[2] * other[1],
            self[2] * other[0] - self[0] * other[2],
            self[0] * other[1] - self[1] * other[0],
        ])
    }

    /// Returns a unit vector with the same direction as this.
    pub fn normed(&self) -> UnitVec3 {
        UnitVec3(*self / self.norm())
    }
}

impl Deref for Vec3 {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Vec3 {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Default for Vec3 {
    fn default() -> Self {
        Self([0.0, 0.0, 0.0])
    }
}

impl Add<Vec3> for Vec3 {
    type Output = Self;

    fn add(self, other: Vec3) -> Self::Output {
        Self([self[0] + other[0], self[1] + other[1], self[2] + other[2]])
    }
}

impl AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        self[0] += other[0];
        self[1] += other[1];
        self[2] += other[2];
    }
}

impl Add<Vec3> for Point3 {
    type Output = Self;

    fn add(self, vec: Vec3) -> Self::Output {
        Self([self[0] + vec[0], self[1] + vec[1], self[2] + vec[2]])
    }
}

impl AddAssign<Vec3> for Point3 {
    fn add_assign(&mut self, vec: Vec3) {
        self[0] += vec[0];
        self[1] += vec[1];
        self[2] += vec[2];
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self([-self[0], -self[1], -self[2]])
    }
}

impl Sub<Vec3> for Vec3 {
    type Output = Self;

    fn sub(self, other: Vec3) -> Self::Output {
        Self([self[0] - other[0], self[1] - other[1], self[2] - other[2]])
    }
}

impl SubAssign<Vec3> for Vec3 {
    fn sub_assign(&mut self, other: Vec3) {
        self[0] -= other[0];
        self[1] -= other[1];
        self[2] -= other[2];
    }
}

impl Sub<Vec3> for Point3 {
    type Output = Self;

    fn sub(self, vec: Vec3) -> Self::Output {
        Self([self[0] - vec[0], self[1] - vec[1], self[2] - vec[2]])
    }
}

impl SubAssign<Vec3> for Point3 {
    fn sub_assign(&mut self, vec: Vec3) {
        self[0] -= vec[0];
        self[1] -= vec[1];
        self[2] -= vec[2];
    }
}

impl Sub<Point3> for Point3 {
    type Output = Vec3;

    fn sub(self, other: Self) -> Self::Output {
        Vec3([self[0] - other[0], self[1] - other[1], self[2] - other[2]])
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self::Output {
        Self([self[0] * scalar, self[1] * scalar, self[2] * scalar])
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, scalar: f64) {
        self[0] *= scalar;
        self[1] *= scalar;
        self[2] *= scalar;
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, scalar: f64) -> Self::Output {
        Self([self[0] / scalar, self[1] / scalar, self[2] / scalar])
    }
}

impl DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, scalar: f64) {
        self[0] /= scalar;
        self[1] /= scalar;
        self[2] /= scalar;
    }
}

impl Sum<Vec3> for Vec3 {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Vec3>,
    {
        iter.fold(Vec3([0.0, 0.0, 0.0]), Add::add)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct UnitVec3(Vec3);

impl Deref for UnitVec3 {
    type Target = Vec3;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Neg for UnitVec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Onb([UnitVec3; 3]);

impl Onb {
    pub fn from23(v: Vec3, w: Vec3) -> Self {
        let w = w.normed();
        let u = v.cross(*w).normed();
        let v = w.cross(*u).normed();

        Self([u, v, w])
    }
}

impl Deref for Onb {
    type Target = [UnitVec3; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
/// A ray in three-dimensional space, i.e., a set of the form
/// {A + tb | t ∈ ℝ+}, where A is a [Point3] and b a [Vec3].
#[derive(Clone, Copy, Debug)]
pub struct Ray {
    pub origin: Point3,
    pub direction: Vec3,
}

impl Ray {
    /// Returns the point origin + t * direction.
    pub fn at(&self, t: f64) -> Point3 {
        if t <= 0.0 {
            self.origin
        } else {
            self.origin + self.direction * t
        }
    }
}

pub struct InsideUnitSphere;

impl Distribution<Vec3> for InsideUnitSphere {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        let range = Uniform::from(-1.0..1.0);
        let mut vec = Vec3([range.sample(rng), range.sample(rng), range.sample(rng)]);

        while vec.norm_squared() >= 1.0 {
            vec = Vec3([range.sample(rng), range.sample(rng), range.sample(rng)]);
        }

        vec
    }
}

impl Distribution<Point3> for InsideUnitSphere {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Point3 {
        let p: Vec3 = self.sample(rng);
        Point3(p.0)
    }
}
pub struct OnUnitSphere;

impl Distribution<UnitVec3> for OnUnitSphere {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> UnitVec3 {
        let phi = rng.gen_range(0.0, 2.0 * std::f64::consts::PI);
        let z = rng.gen_range(-1.0, 1.0);
        let r = 1.0 - z * z;

        let vec = Vec3([r * phi.cos(), r * phi.sin(), z]);

        vec.normed()
    }
}

impl Distribution<Point3> for OnUnitSphere {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Point3 {
        let p: UnitVec3 = self.sample(rng);
        Point3((p.0).0)
    }
}

pub fn random_unit_vector<R: Rng + ?Sized>(rng: &mut R) -> UnitVec3 {
    OnUnitSphere.sample(rng)
}

pub struct UnitDisc;

impl Distribution<Vec3> for UnitDisc {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        let range = Uniform::from(-1.0..1.0);
        let mut p = Vec3([range.sample(rng), range.sample(rng), 0.0]);

        while p.norm() >= 1.0 {
            p = Vec3([range.sample(rng), range.sample(rng), 0.0]);
        }

        p
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Interval {
    Empty,
    Nonempty(f64, f64),
}

impl Interval {
    pub fn new(start: f64, end: f64) -> Self {
        if start > end {
            Interval::Empty
        } else {
            Interval::Nonempty(start, end)
        }
    }
}

fn intersect_interval(ray: &Ray, interval: &Interval, coord: usize) -> Interval {
    assert!(coord < 3);
    match interval {
        Interval::Empty => Interval::Empty,
        Interval::Nonempty(start, end) => {
            let origin = ray.origin[coord];
            let direction = ray.direction[coord];

            if direction == 0.0 {
                return Interval::Empty;
            }

            let (t0, t1) = ((start - origin) / direction, (end - origin) / direction);

            Interval::new(t0.min(t1), t0.max(t1))
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    pub min: Point3,
    pub max: Point3,
}

impl BoundingBox {
    pub fn hit(&self, ray: &Ray, mut tmin: f64, mut tmax: f64) -> bool {
        for i in 0..3 {
            let d = ray.direction[i];
            let o = ray.origin[i];
            let mut t0 = (self.min[i] - o) / d;
            let mut t1 = (self.max[i] - o) / d;
            if d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            tmin = t0.max(tmin);
            tmax = t1.min(tmax);

            if tmax <= tmin {
                return false;
            }
        }

        true
    }
}

impl Add<BoundingBox> for BoundingBox {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let new_min = Point3([
            self.min[0].min(other.min[0]),
            self.min[1].min(other.min[1]),
            self.min[2].min(other.min[2]),
        ]);

        let new_max = Point3([
            self.max[0].max(other.max[0]),
            self.max[1].max(other.max[1]),
            self.max[2].max(other.max[2]),
        ]);

        Self {
            min: new_min,
            max: new_max,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn inside_unit_sphere() {
        let mut rng = rand::thread_rng();

        for _ in 0..20 {
            let v: Vec3 = InsideUnitSphere.sample(&mut rng);

            assert!(v.norm() < 1.0);
        }
    }

    #[test]
    fn on_unit_sphere() {
        let mut rng = rand::thread_rng();

        for _ in 0..20 {
            let v: UnitVec3 = OnUnitSphere.sample(&mut rng);

            assert!((v.norm() - 1.0).abs() <= 1e-10);
        }
    }
}
