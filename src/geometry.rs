use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::iter::Sum;
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
    pub in_vec: Vec3,
}

impl IntersectionPoint {
    pub fn random_scatter<R: Rng + ?Sized>(&self, rng: &mut R) -> Ray {
        let Self { point, normal, .. } = *self;
        let sphere = Sphere {
            center: point + normal,
            radius: 1.0,
        };
        let target = On(sphere).sample(rng);
        Ray {
            origin: point,
            direction: target - point,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Face {
    Front,
    Back,
}

/// A point in three-dimensional space.
#[derive(Clone, Copy, Debug)]
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
    pub fn normed(&self) -> Self {
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
        assert!(t > 0.0);
        self.origin + self.direction * t
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
                let normal = (point - *center).normed();
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
                    in_vec: *dir,
                })
            } else {
                None
            }
        }
    }
}

pub struct Inside(Sphere);

impl Distribution<Point3> for Inside {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Point3 {
        let range = Uniform::from(-self.0.radius..self.0.radius);
        let mut vec = Vec3([range.sample(rng), range.sample(rng), range.sample(rng)]);

        while vec.norm() >= self.0.radius {
            vec = Vec3([range.sample(rng), range.sample(rng), range.sample(rng)]);
        }

        self.0.center + vec
    }
}

pub struct On(Sphere);

impl Distribution<Point3> for On {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Point3 {
        let phi = rng.gen_range(0.0, 2.0 * std::f64::consts::PI);
        let z = rng.gen_range(-1.0, 1.0);
        let r = 1.0 - z * z;

        let vec = Vec3([r * phi.cos(), r * phi.sin(), z]) * self.0.radius;

        self.0.center + vec
    }
}

pub fn random_unit_vector<R: Rng + ?Sized>(rng: &mut R) -> Vec3 {
    let point = On(Sphere {
        center: Point3::default(),
        radius: 1.0,
    })
    .sample(rng);

    Vec3(point.0)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn inside_sphere_random_point() {
        let mut rng = rand::thread_rng();

        let sphere = Sphere {
            center: Point3([1.0, 0.0, 0.0]),
            radius: 2.0,
        };

        for _ in 0..20 {
            let p = Inside(sphere).sample(&mut rng);

            assert!(p.dist(&sphere.center) < sphere.radius);
        }
    }

    #[test]
    fn on_sphere_random_point() {
        let mut rng = rand::thread_rng();

        let sphere = Sphere {
            center: Point3([1.0, 0.0, 0.0]),
            radius: 2.0,
        };

        for _ in 0..20 {
            let p = On(sphere).sample(&mut rng);

            assert!((p.dist(&sphere.center) - sphere.radius).abs() / sphere.radius <= 0.15);
        }
    }
}
