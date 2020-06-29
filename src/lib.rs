pub mod geometry;

pub mod light {
    use super::geometry::*;
    use rand::Rng;
    use std::iter::Sum;
    use std::ops::{Deref, Mul};

    pub fn ray_color<I, R>(ray: &Ray, world: &I, rng: &mut R, depth: usize) -> Color
    where
        I: Intersection,
        R: Rng + ?Sized,
    {
        if depth == 0 {
            return Color::new(0.0, 0.0, 0.0);
        }
        match world.intersection(ray, 0.0, f64::INFINITY) {
            Some(intersection_point) => {
                let ray = intersection_point.random_scatter(rng);
                ray_color(&ray, world, rng, depth - 1) * 0.5
            }
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
