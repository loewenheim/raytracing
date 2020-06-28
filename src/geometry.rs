
use super::Rgb;
use std::ops::{
    Add, AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// A trait for things that can be intersected by rays (such as spheres).
pub trait Intersection {
    /// Computes the intersection point of the given
    /// ray with this object.
    /// Returns None if the ray misses.
    fn intersection(&self, ray: &Ray, tmin: f64, tmax: f64) -> Option<IntersectionPoint>;
}

/// Contains information about the intersection of a ray
/// with an object: the actual point, the normal vector
/// at that point, and the parameter of the ray.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    pub point: Point3,
    pub normal: Vec3,
    pub t: f64,
    pub face: Face,
}

#[derive(Debug, Clone, Copy)]
pub enum Face {
    Front,
    Back,
}

/// A point in three-dimensional space.
#[derive(Clone, Copy, Debug)]
pub struct Point3(pub [f64; 3]);

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
#[derive(Clone, Copy, Debug)]
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
    pub fn unit(&self) -> Self {
        *self / self.norm()
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

impl Into<Rgb> for Vec3 {
    fn into(self) -> Rgb {
        let [red, green, blue] = self.0;
        let convert = |x| (x * 255.999) as u8;
        image::Rgb([convert(red), convert(green), convert(blue)])
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

/// A ray in three-dimensional space, i.e., a set of the form
/// {A + tb | t ∈ ℝ+}, where A is a [Point3] and b a [Vec3].
#[derive(Clone, Copy)]
pub struct Ray {
    pub origin: Point3,
    pub direction: Vec3,
}

impl Ray {
    /// Returns the point origin + t * direction.
    pub fn at(&self, t: f64) -> Point3 {
        assert!(t >= 0.0);
        self.origin + self.direction * t
    }

    /// Rgbs the point the ray hits.
    pub fn color<I: Intersection>(&self, world: &I) -> super::Rgb {
        match world.intersection(self, 0.0, f64::INFINITY) {
            Some(IntersectionPoint { normal, .. }) => {
                ((normal + Vec3([1.0, 1.0, 1.0])) * 0.5).into()
            }
            None => {
                let unit = self.direction.unit();
                let t = 0.5 * (unit[1] + 1.0);
                let blue = Vec3([0.5, 0.7, 1.0]);
                let white = Vec3([1.0, 1.0, 1.0]);

                (blue * t + white * (1.0 - t)).into()
            }
        }
    }
}

/// A sphere in three-dimensional space.
#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    pub center: Point3,
    pub radius: f64,
}

impl Intersection for Sphere {
    fn intersection(&self, ray: &Ray, tmin: f64, tmax: f64) -> Option<IntersectionPoint> {
        let Sphere { center, radius: r } = self;
        let Ray {
            origin: o,
            direction: dir,
        } = ray;
        let oc = *o - *center;
        let a = dir.norm_squared();
        let half_b = oc.dot(*dir);
        let c = oc.norm_squared() - r.powi(2);
        let discriminant = half_b.powi(2) - a * c;

        if discriminant < 0.0 {
            None
        } else {
            let t = (-half_b - discriminant.sqrt()) / a;
            if t > tmin && t < tmax {
                let point = ray.at(t);
                let normal = (point - *center) / *r;
                let face = if dir.dot(normal) < 0.0 {
                    Face::Front
                } else {
                    Face::Back
                };
                Some(IntersectionPoint {
                    t,
                    point,
                    normal,
                    face,
                })
            } else {
                None
            }
        }
    }
}

impl Intersection for Vec<Box<dyn Intersection>> {
    fn intersection(&self, ray: &Ray, tmin: f64, mut tmax: f64) -> Option<IntersectionPoint> {
        let mut intersection_point = None;

        for object in self.iter() {
            if let Some(new_ip) = object.intersection(ray, tmin, tmax) {
                intersection_point = Some(new_ip);
                tmax = new_ip.t;
            }
        }

        intersection_point
    }
}
