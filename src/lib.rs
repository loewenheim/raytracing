pub mod geometry;

use geometry::{Intersection, IntersectionPoint, Ray};
use light::Color;
use materials::Material;
use rand::Rng;

pub trait Reflective<R: Rng + ?Sized> {
    fn reflect(&self, r: &Ray, tmin: f64, tmax: f64, rng: &mut R) -> Option<Scattered>;
}

pub struct MaterialObject<I, M, R>
where
    I: Intersection,
    M: Material<R>,
    R: Rng + ?Sized,
{
    shape: I,
    material: M,
    phantom: std::marker::PhantomData<R>,
}

impl<I, M, R> MaterialObject<I, M, R>
where
    I: Intersection,
    M: Material<R>,
    R: Rng + ?Sized,
{
    pub fn new(shape: I, material: M) -> Self {
        Self {
            shape,
            material,
            phantom: std::marker::PhantomData,
        }
    }
}

pub struct Scattered {
    ray: Ray,
    attenuation: Color,
    intersection_point: IntersectionPoint,
}

impl<I, M, R> Reflective<R> for MaterialObject<I, M, R>
where
    I: Intersection,
    M: Material<R>,
    R: Rng + ?Sized,
{
    fn reflect(&self, r: &Ray, tmin: f64, tmax: f64, rng: &mut R) -> Option<Scattered> {
        self.shape.intersection(r, tmin, tmax).as_ref().and_then(|p| {
            self.material.scatter(p, rng).map(|(attenuation, ray)|
            Scattered {
                ray,
                attenuation,
                intersection_point: *p,
            })
        })
    }
}

impl<R> Reflective<R> for Vec<Box<dyn Reflective<R>>>
where
    R: Rng + ?Sized,
{
    fn reflect(&self, ray: &Ray, tmin: f64, mut tmax: f64, rng: &mut R) -> Option<Scattered> {
        let mut scattered = None;

        for object in self.iter() {
            if let Some(new_scattered) = object.reflect(ray, tmin, tmax, rng) {
                tmax = new_scattered.intersection_point.t;
                scattered = Some(new_scattered);
            }
        }

        scattered
    }
}

pub mod materials {
    use super::geometry::{random_unit_vector, IntersectionPoint, Ray, Vec3, Face};
    use super::light::Color;
    use rand::Rng;

    pub trait Material<R: Rng + ?Sized> {
        fn scatter(&self, p: &IntersectionPoint, rng: &mut R) -> Option<(Color, Ray)>;
    }

    pub struct Lambertian {
        pub albedo: Color,
    }

    impl<R: Rng + ?Sized> Material<R> for Lambertian {
        fn scatter(
            &self,
            IntersectionPoint { normal, point, .. }: &IntersectionPoint,
            rng: &mut R,
        ) -> Option<(Color, Ray)> {
            let direction = random_unit_vector(rng) + *normal;
            let scattered = Ray {
                origin: *point,
                direction,
            };

            Some((self.albedo, scattered))
        }
    }

    pub struct Metal {
        pub albedo: Color,
        pub fuzz: f64,
    }

    impl<R: Rng + ?Sized> Material<R> for Metal {
        fn scatter(
            &self,
            IntersectionPoint {
                normal,
                point,
                in_vec,
                ..
            }: &IntersectionPoint,
            rng: &mut R,
        ) -> Option<(Color, Ray)> {
            let reflected = reflect(in_vec, normal);
            let scattered = reflected + random_unit_vector(rng) * self.fuzz;
            if scattered.dot(*normal) >= 0.0 {
            Some((
                self.albedo,
                Ray {
                    origin: *point,
                    direction: scattered,
                },
            ))
            } else {
                None
            }
        }
    }

    pub struct Dielectric {
        pub ref_index: f64,
    }

    impl<R: Rng + ?Sized> Material<R> for Dielectric {
        fn scatter(
            &self,
            IntersectionPoint {
                normal,
                point,
                in_vec,
                face,
                ..
            }: &IntersectionPoint,
            rng: &mut R,
        ) -> Option<(Color, Ray)> {
            let eta_ratio = match face {
                Face::Front => 1.0/self.ref_index,
                Face::Back => self.ref_index,
            };

            let in_vec_unit = in_vec.normed();
            let cos_theta = -in_vec_unit.dot(*normal).min(1.0);
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

            let direction = if eta_ratio * sin_theta > 1.0 || rng.gen::<f64>() < schlick(cos_theta, eta_ratio) {
                reflect(&in_vec_unit, normal)
            } else {
                refract(&in_vec_unit, normal, eta_ratio)
            };

            Some((Color::new(1.0, 1.0, 1.0), Ray {
                origin: *point,
                direction,
            }))
        }
    }

    fn reflect(&v: &Vec3, &n: &Vec3) -> Vec3 {
        v - n * 2.0 * v.dot(n)
    }

    fn refract(&v: &Vec3, &n: &Vec3, eta_ratio: f64) -> Vec3 {
        let cos_theta = -v.dot(n);
        let out_par = (v + n * cos_theta) * eta_ratio;
        let out_perp = n * -((1.0 - out_par.norm_squared()).sqrt());
        out_par + out_perp
    }

    fn schlick(cosine: f64, ref_index: f64) -> f64 {
        let mut r0 = (1.0 - ref_index) / (1.0 + ref_index);
        r0 = r0 * r0;
        r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
    }
}

pub mod light {
    use super::geometry::*;
    use super::{Reflective, Scattered};
    use rand::Rng;
    use std::iter::Sum;
    use std::ops::{Deref, Mul};

    pub fn ray_color<I, R>(ray: &Ray, world: &I, rng: &mut R, depth: usize) -> Color
    where
        I: Reflective<R>,
        R: Rng + ?Sized,
    {
        if depth == 0 {
            return Color::new(0.0, 0.0, 0.0);
        }
        match world.reflect(ray, 0.001, f64::INFINITY, rng) {
            Some(Scattered {
                ray, attenuation, ..
            }) => attenuation * ray_color(&ray, world, rng, depth - 1),
            None => {
                let unit = ray.direction.normed();
                let t = 0.5 * (unit[1] + 1.0);
                let blue = Color::new(0.5, 0.7, 1.0);
                let white = Color::new(1.0, 1.0, 1.0);

                blue.mix(&white, t)
            }
        }
    }

    pub type Rgb = image::Rgb<u8>;

    #[derive(Debug, Clone, Copy)]
    pub struct Color([f64; 3]);

    impl Sum<Color> for Color {
        fn sum<I: Iterator<Item = Color>>(iter: I) -> Self {
            let mut count: usize = 0;
            let mixed = iter.fold([0.0f64, 0.0, 0.0], |acc, c| {
                count += 1;
                [acc[0] + c.0[0], acc[1] + c.0[1], acc[2] + c.0[2]]
            });

            Color([
                mixed[0] / count as f64,
                mixed[1] / count as f64,
                mixed[2] / count as f64,
            ])
        }
    }

    impl Color {
        pub fn new(red: f64, green: f64, blue: f64) -> Self {
            assert!(
                red >= 0.0
                    && red <= 1.0
                    && green >= 0.0
                    && green <= 1.0
                    && blue >= 0.0
                    && blue <= 1.0,
                format!("{:?} is not a legal color vetcor", (red, green, blue))
            );
            Self([red, green, blue])
        }

        pub fn mix(&self, other: &Self, t: f64) -> Self {
            assert!(t >= 0.0 && t <= 1.0, format!("{} is not in [0, 1]", t));

            let [r1, g1, b1] = self.0;
            let [r2, g2, b2] = other.0;

            Self::new(
                t * r1 + (1.0 - t) * r2,
                t * g1 + (1.0 - t) * g2,
                t * b1 + (1.0 - t) * b2,
            )
        }

        pub fn map<F: FnMut(f64) -> f64>(self, mut f: F) -> Self {
            let [r, g, b] = self.0;
            Color::new(f(r), f(g), f(b))
        }
    }

    impl Into<Rgb> for Color {
        fn into(self) -> Rgb {
            let [red, green, blue] = self.0;
            let convert = |x| (x * 255.999) as u8;
            image::Rgb([convert(red), convert(green), convert(blue)])
        }
    }

    impl Mul<f64> for Color {
        type Output = Self;
        fn mul(self, scalar: f64) -> Self::Output {
            Self([self[0] * scalar, self[1] * scalar, self[2] * scalar])
        }
    }

    impl Mul<Color> for Color {
        type Output = Self;
        fn mul(self, other: Self) -> Self::Output {
            Self([self[0] * other[0], self[1] * other[1], self[2] * other[2]])
        }
    }

    impl Deref for Color {
        type Target = [f64; 3];

        fn deref(&self) -> &Self::Target {
            &self.0
        }
    }
}

pub mod camera {
    use super::geometry::*;

    #[derive(Debug, Clone, Copy)]
    pub struct Camera {
        origin: Point3,
        horizontal: Vec3,
        vertical: Vec3,
        lower_left_corner: Point3,
    }

    impl Camera {
        pub fn ray(&self, u: f64, v: f64) -> Ray {
            let direction =
                (self.lower_left_corner + self.horizontal * u + self.vertical * v) - self.origin;
            Ray {
                origin: self.origin,
                direction,
            }
        }
    }

    impl Default for Camera {
        fn default() -> Self {
            camera(1.0)
                .aspect_ratio(16.0 / 9.0)
                .viewport_height(2.0)
                .build()
        }
    }

    pub struct CameraBuilder {
        aspect_ratio: Option<f64>,
        viewport_height: Option<f64>,
        viewport_width: Option<f64>,
        focal_length: f64,
    }

    impl CameraBuilder {
        pub fn aspect_ratio(mut self, aspect_ratio: f64) -> Self {
            match (self.viewport_height, self.viewport_width) {
                (Some(_), Some(_)) => panic!("Viewport height and width are already set"),
                (Some(h), None) => self.viewport_width = Some(h * aspect_ratio),
                (None, Some(w)) => self.viewport_height = Some(w / aspect_ratio),
                (None, None) => {}
            }
            self.aspect_ratio = Some(aspect_ratio);
            self
        }

        pub fn viewport_height(mut self, viewport_height: f64) -> Self {
            match (self.aspect_ratio, self.viewport_width) {
                (Some(_), Some(_)) => panic!("Aspect ratio and viewport width are already set"),
                (Some(ar), None) => self.viewport_width = Some(viewport_height * ar),
                (None, Some(w)) => self.aspect_ratio = Some(w / viewport_height),
                (None, None) => {}
            }
            self.viewport_height = Some(viewport_height);
            self
        }

        pub fn viewport_width(mut self, viewport_width: f64) -> Self {
            match (self.aspect_ratio, self.viewport_height) {
                (Some(_), Some(_)) => panic!("Aspect ratio and viewport height are already set"),
                (Some(ar), None) => self.viewport_height = Some(viewport_width / ar),
                (None, Some(h)) => self.aspect_ratio = Some(viewport_width / h),
                (None, None) => {}
            }
            self.viewport_width = Some(viewport_width);
            self
        }

        pub fn build(self) -> Camera {
            let origin = Point3::default();
            let horizontal = Vec3([
                self.viewport_width.expect(
                    "You need to set two of viewport height, viewport width, and aspect ratio",
                ),
                0.0,
                0.0,
            ]);
            let vertical = Vec3([
                0.0,
                self.viewport_height.expect(
                    "You need to set two of viewport height, viewport width, and aspect ratio",
                ),
                0.0,
            ]);
            let lower_left_corner =
                origin - horizontal / 2.0 - vertical / 2.0 - Vec3([0.0, 0.0, self.focal_length]);

            Camera {
                origin,
                horizontal,
                vertical,
                lower_left_corner,
            }
        }
    }
    pub fn camera(focal_length: f64) -> CameraBuilder {
        CameraBuilder {
            aspect_ratio: None,
            viewport_height: None,
            viewport_width: None,
            focal_length,
        }
    }
}
