
use super::geometry::{random_unit_vector, Face, IntersectionPoint, Ray, Vec3};
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
        let direction = *random_unit_vector(rng) + **normal;
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
        let scattered = reflected + *random_unit_vector(rng) * self.fuzz;
        if scattered.dot(**normal) >= 0.0 {
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
    fn scatter(&self, p: &IntersectionPoint, rng: &mut R) -> Option<(Color, Ray)> {
        let eta_ratio = match p.face() {
            Face::Front => 1.0 / self.ref_index,
            Face::Back => self.ref_index,
        };

        let in_vec_unit = p.in_vec.normed();
        let cos_theta = -in_vec_unit.dot(*p.normal).min(1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        let direction =
            if eta_ratio * sin_theta > 1.0 || rng.gen::<f64>() < schlick(cos_theta, eta_ratio) {
                reflect(&in_vec_unit, &p.normal)
            } else {
                refract(&in_vec_unit, &p.normal, eta_ratio)
            };

        Some((
            Color::new(1.0, 1.0, 1.0),
            Ray {
                origin: p.point,
                direction,
            },
        ))
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
