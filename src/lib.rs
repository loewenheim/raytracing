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
        self.shape
            .intersection(r, tmin, tmax)
            .as_ref()
            .and_then(|p| {
                self.material
                    .scatter(p, rng)
                    .map(|(attenuation, ray)| Scattered {
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
        u: Vec3,
        v: Vec3,
        w: Vec3,
        lens_radius: f64,
    }

    impl Camera {
        pub fn new(
            look_from: Point3,
            looking: Vec3,
            vup: Vec3,
            vfov: f64,
            aspect_ratio: f64,
            aperture: f64,
            focus_dist: f64,
        ) -> Self {
            let h = (vfov / 2.0).to_radians().tan();
            let viewport_height = 2.0 * h;
            let viewport_width = viewport_height * aspect_ratio;

            let w = (-looking).normed();
            let u = vup.cross(w).normed();
            let v = w.cross(u);

            let horizontal = u * viewport_width * focus_dist;
            let vertical = v * viewport_height * focus_dist;
            let lower_left_corner = look_from - horizontal / 2.0 - vertical / 2.0 - w * focus_dist;
            let lens_radius = aperture / 2.0;

            Self {
                origin: look_from,
                horizontal,
                vertical,
                lower_left_corner,
                vfov,
                aspect_ratio,
                u,
                v,
                w,
                lens_radius,
            }
        }

        pub fn ray<R: Rng + ?Sized>(&self, s: f64, t: f64, rng: &mut R) -> Ray {
            let rd = UnitDisc.sample(rng) * self.lens_radius;
            let offset = self.u * rd[0] + self.v * rd[1];
            let origin = self.origin + offset;
            let direction =
                (self.lower_left_corner + self.horizontal * s + self.vertical * t) - origin;
            Ray { origin, direction }
        }
    }
}
