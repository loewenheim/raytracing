mod shape;

use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::iter::Sum;
use std::ops::{
    Add, AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, Neg, Range, Sub, SubAssign,
};

pub use shape::{Boundable, BoundingBox, Face, Intersection, IntersectionPoint, Shape};

pub trait Rotate {
    fn rotate(self, axis: Axis, angle: f64) -> Self;
}

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Debug, Clone, Copy)]
pub enum Axes {
    XY,
    XZ,
    YZ,
}

impl Axes {
    fn coords(&self) -> (usize, usize, usize) {
        match self {
            Axes::XY => (0, 1, 2),
            Axes::XZ => (0, 2, 1),
            Axes::YZ => (1, 2, 0),
        }
    }
}

/// A point in three-dimensional space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3(pub [f64; 3]);

impl Point3 {
    pub fn dist(&self, other: &Point3) -> f64 {
        (*self - *other).norm()
    }
}

impl Rotate for Point3 {
    fn rotate(self, axis: Axis, angle: f64) -> Self {
        let mut out = self;
        let angle = angle.to_radians();
        let cos = angle.cos();
        let sin = angle.sin();
        let (c0, c1) = match axis {
            Axis::X => (1, 2),
            Axis::Y => (2, 0),
            Axis::Z => (0, 1),
        };

        out[c0] = self[c0] * cos - self[c1] * sin;
        out[c1] = self[c0] * sin + self[c1] * cos;

        out
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

    pub fn random<R: Rng + ?Sized>(range: Range<f64>, rng: &mut R) -> Self {
        let range = Uniform::from(range);
        Vec3([range.sample(rng), range.sample(rng), range.sample(rng)])
    }
}

impl Rotate for Vec3 {
    fn rotate(self, axis: Axis, angle: f64) -> Self {
        let mut out = self;
        let angle = angle.to_radians();
        let cos = angle.cos();
        let sin = angle.sin();
        let (c0, c1) = match axis {
            Axis::X => (1, 2),
            Axis::Y => (2, 0),
            Axis::Z => (0, 1),
        };

        out[c0] = self[c0] * cos - self[c1] * sin;
        out[c1] = self[c0] * sin + self[c1] * cos;

        out
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

impl Mul<Vec3> for Vec3 {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        Self([self[0] * other[0], self[1] * other[1], self[2] * other[2]])
    }
}

impl MulAssign<Vec3> for Vec3 {
    fn mul_assign(&mut self, other: Self) {
        self[0] *= other[0];
        self[1] *= other[1];
        self[2] *= other[2];
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

impl UnitVec3 {
    pub fn rotate(self, axis: Axis, angle: f64) -> Self {
        Self(self.0.rotate(axis, angle))
    }
}

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

pub struct InsideUnitDisc;

impl Distribution<Vec3> for InsideUnitDisc {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        let range = Uniform::from(-1.0..1.0);
        let mut p = Vec3([range.sample(rng), range.sample(rng), 0.0]);

        while p.norm() >= 1.0 {
            p = Vec3([range.sample(rng), range.sample(rng), 0.0]);
        }

        p
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

    #[test]
    fn rotated_sphere() {
        let sphere = Shape::Sphere {
            center: Point3([0.0, 0.0, 0.0]),
            radius: 2.0,
        }
        .rotate(Axis::Y, 45.0);

        let ray = Ray {
            origin: Point3([10.0, 0.0, 0.0]),
            direction: Vec3([-2.0, 0.0, 0.0]),
        };

        assert!(sphere.intersect(&ray, 0.0, 10.0, 0.0).is_some())
    }

    #[test]
    fn rotated_box() {
        let my_box = Shape::new_box(Point3([0.0, 0.0, 0.0]), Point3([10.0, 10.0, 10.0]))
            .rotate(Axis::Y, 45.0);
        let ray = Ray {
            origin: Point3([-10.0, 5.0, 12.0]),
            direction: Vec3([1.0, 0.0, -1.0]),
        };
        let ip = my_box.intersect(&ray, 0.0, 20.0, 0.0);

        println!("{:?}", ip);

        assert!(ip.is_some());
    }
}
