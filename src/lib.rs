pub mod color {
    use std::fmt;

    /// A color with fractional RGB values (intended to be in [0, 1]).
    #[derive(Clone, Copy, Debug)]
    pub struct Color {
        pub red: f64,
        pub green: f64,
        pub blue: f64,
    }

    impl Color {
        /// Standard white color.
        pub fn white() -> Self {
            Self {
                red: 1.0,
                green: 1.0,
                blue: 1.0,
            }
        }

        /// Standard red color.
        pub fn red() -> Self {
            Self {
                red: 1.0,
                green: 0.0,
                blue: 0.0,
            }
        }

        /// Blends this color with another. A `ratio` of 1.0 results in `self`,
        /// 0.0 results in `other`.
        pub fn blend(&self, other: &Self, ratio: f64) -> Self {
            assert!(ratio >= 0.0);
            assert!(ratio <= 1.0);
            Self {
                red: ratio * self.red + (1.0 - ratio) * other.red,
                green: ratio * self.green + (1.0 - ratio) * other.green,
                blue: ratio * self.blue + (1.0 - ratio) * other.blue,
            }
        }
    }

    impl fmt::Display for Color {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let factor = 255.999;
            write!(
                f,
                "{} {} {}",
                (self.red * factor) as u8,
                (self.green * factor) as u8,
                (self.blue * factor) as u8
            )
        }
    }

    impl From<super::geometry::Vec3> for Color {
        fn from(vec: super::geometry::Vec3) -> Self {
            Self {
                red: vec[0],
                green: vec[1],
                blue: vec[2],
            }
        }
    }
}

pub mod geometry {
    use super::color::Color;
    use std::ops::{
        Add, AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    };

    /// A point in three-dimensional space.
    #[derive(Clone, Copy, Debug)]
    pub struct Point3([f64; 3]);

    impl Point3 {
        /// Returns a new point with the given three coordinates.
        pub fn new(x: f64, y: f64, z: f64) -> Self {
            Self([x, y, z])
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
    pub struct Vec3([f64; 3]);

    impl Vec3 {
        /// Returns a new vector with the given three components.
        pub fn new(x: f64, y: f64, z: f64) -> Self {
            Self([x, y, z])
        }

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
            Self::new(
                self[1] * other[2] - self[2] * other[1],
                self[2] * other[0] - self[0] * other[2],
                self[0] * other[1] - self[1] * other[0],
            )
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

    impl Add<Vec3> for Vec3 {
        type Output = Self;

        fn add(self, other: Vec3) -> Self::Output {
            Self::new(self[0] + other[0], self[1] + other[1], self[2] + other[2])
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
            Self::new(self[0] + vec[0], self[1] + vec[1], self[2] + vec[2])
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
            Self::new(-self[0], -self[1], -self[2])
        }
    }

    impl Sub<Vec3> for Vec3 {
        type Output = Self;

        fn sub(self, other: Vec3) -> Self::Output {
            Self::new(self[0] - other[0], self[1] - other[1], self[2] - other[2])
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
            Self::new(self[0] - vec[0], self[1] - vec[1], self[2] - vec[2])
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
            Vec3::new(self[0] - other[0], self[1] - other[1], self[2] - other[2])
        }
    }

    impl Mul<f64> for Vec3 {
        type Output = Self;
        fn mul(self, scalar: f64) -> Self::Output {
            Self::new(self[0] * scalar, self[1] * scalar, self[2] * scalar)
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
            Self::new(self[0] / scalar, self[1] / scalar, self[2] / scalar)
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

        /// Colors the point the ray hits.
        pub fn color(&self) -> super::color::Color {
            let sphere = Sphere {
                center: Point3::new(0.0, 0.0, -1.0),
                radius: 0.5,
            };

            match intersection(self, &sphere) {
                Some(t) => {
                    let n = (self.at(t) - Point3::new(0.0, 0.0, -1.0)).unit();
                    Color::from((n + Vec3::new(1.0, 1.0, 1.0)) * 0.5)
                }
                None => {
                    let unit = self.direction.unit();
                    let t = 0.5 * (unit[1] + 1.0);
                    let blue = Color {
                        red: 0.5,
                        green: 0.7,
                        blue: 1.0,
                    };

                    blue.blend(&Color::white(), t)
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

    /// Computes the intersection point of the given
    /// ray with the given sphere (as a ray parameter).
    /// Returns None if the ray misses the sphere.
    pub fn intersection(
        Ray {
            origin: o,
            direction: dir,
        }: &Ray,
        Sphere {
            center: c,
            radius: r,
        }: &Sphere,
    ) -> Option<f64> {
        let oc = *o - *c;
        let a = dir.dot(*dir);
        let b = 2.0 * oc.dot(*dir);
        let c = oc.dot(oc) - r.powi(2);
        let discriminant = b.powi(2) - 4.0 * a * c;

        if discriminant < 0.0 {
            None
        } else {
            Some((-b - discriminant.sqrt()) / (2.0 * a))
        }
    }
}
