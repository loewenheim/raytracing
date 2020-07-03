use super::geometry::{Point3, Vec3};
use image::RgbImage;
use perlin_noise::PerlinNoise;
use std::sync::Arc;

#[derive(Clone)]
pub enum Texture {
    SolidColor(Vec3),
    Checkered {
        odd: Box<Texture>,
        even: Box<Texture>,
    },
    Noise(Arc<PerlinNoise>),
    Image(RgbImage),
}

impl Texture {
    pub fn image<P: AsRef<std::path::Path>>(path: P) -> Self {
        let image = image::open(path).unwrap().as_rgb8().unwrap().to_owned();
        Self::Image(image)
    }

    pub fn color_at(&self, u: f64, v: f64, p: &Point3) -> Vec3 {
        match self {
            Self::SolidColor(color) => *color,
            Self::Checkered { odd, even } => {
                if (10.0 * p[0]).sin() * (10.0 * p[1]).sin() * (10.0 * p[2]).sin() < 0.0 {
                    odd.color_at(u, v, p)
                } else {
                    even.color_at(u, v, p)
                }
            }
            Self::Noise(perlin) => Vec3([1.0, 1.0, 1.0]) * perlin.get3d(p.0),
            Self::Image(image) => {
                let u = u.max(0.0).min(1.0);
                let v = 1.0 - v.max(0.0).min(1.0);

                let i = (u * image.width() as f64) as u32;
                let j = (v * image.height() as f64) as u32;

                let i = i.min(image.width() - 1);
                let j = j.min(image.height() - 1);

                let pixel = image.get_pixel(i, j);

                Vec3([
                    pixel[0] as f64 / 255.0,
                    pixel[1] as f64 / 255.0,
                    pixel[2] as f64 / 255.0,
                ])
            }
        }
    }
}
