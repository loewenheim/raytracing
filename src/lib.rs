pub type Rgb = image::Rgb<u8>;

pub mod geometry;

pub mod light {
    use super::geometry::*;
    use rand::Rng;

    pub fn color_vec<I, R>(ray: &Ray, world: &I, rng: &mut R, depth: usize) -> Vec3
    where
        I: Intersection,
        R: Rng + ?Sized,
    {
        if depth == 0 {
            return Vec3::default();
        }
        match world.intersection(ray, 0.0, f64::INFINITY) {
            Some(intersection_point) => {
                let ray = intersection_point.random_scatter(rng);
                color_vec(&ray, world, rng, depth - 1)
            }
            None => {
                let unit = ray.direction.normed();
                let t = 0.5 * (unit[1] + 1.0);
                let blue = Vec3([0.5, 0.7, 1.0]);
                let white = Vec3([1.0, 1.0, 1.0]);

                blue * t + white * (1.0 - t)
            }
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
