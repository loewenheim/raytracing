use super::geometry::{random_unit_vector, Face, IntersectionPoint, Ray, UnitVec3, Vec3};
use super::light::Color;
use super::Scattered;
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub enum Material {
    Lambertian { albedo: Color },

    Metal { albedo: Color, fuzz: f64 },

    Dielectric { refraction_index: f64 },
}

impl Material {
    pub fn glass() -> Self {
        Material::Dielectric {
            refraction_index: 1.5,
        }
    }

    pub fn scatter<R: Rng + ?Sized>(
        &self,
        p: &IntersectionPoint,
        rng: &mut R,
    ) -> Option<Scattered> {
        match *self {
            Material::Lambertian { albedo } => {
                let direction = *random_unit_vector(rng) + *p.normal;
                let scattered = Ray {
                    origin: p.point,
                    direction,
                };

                Some(Scattered {
                    attenuation: albedo,
                    ray: scattered,
                })
            }

            Material::Metal { albedo, fuzz } => {
                let reflected = reflect(&p.in_vec, &p.normal);
                let scattered = *reflected + *random_unit_vector(rng) * fuzz;
                if scattered.dot(*p.normal) >= 0.0 {
                    Some(Scattered {
                        attenuation: albedo,
                        ray: Ray {
                            origin: p.point,
                            direction: scattered,
                        },
                    })
                } else {
                    None
                }
            }

            Material::Dielectric { refraction_index } => {
                let eta_ratio = match p.face() {
                    Face::Front => 1.0 / refraction_index,
                    Face::Back => refraction_index,
                };

                let in_vec_unit = p.in_vec.normed();
                let cos_theta = -in_vec_unit.dot(*p.normal).min(1.0);
                let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

                let direction = if eta_ratio * sin_theta > 1.0
                    || rng.gen::<f64>() < schlick(cos_theta, eta_ratio)
                {
                    *reflect(&in_vec_unit, &p.normal)
                } else {
                    refract(&in_vec_unit, &p.normal, eta_ratio)
                };

                Some(Scattered {
                    attenuation: Color::new(1.0, 1.0, 1.0),
                    ray: Ray {
                        origin: p.point,
                        direction,
                    },
                })
            }
        }
    }
}

fn reflect(&v: &UnitVec3, &n: &UnitVec3) -> UnitVec3 {
    (*v - *n * 2.0 * (*v).dot(*n)).normed()
}

fn refract(&v: &UnitVec3, &n: &UnitVec3, eta_ratio: f64) -> Vec3 {
    let cos_theta = -v.dot(*n);
    let out_par = (*v + *n * cos_theta) * eta_ratio;
    let out_perp = *n * -((1.0 - out_par.norm_squared()).sqrt());
    out_par + out_perp
}

fn schlick(cosine: f64, refraction_index: f64) -> f64 {
    let mut r0 = (1.0 - refraction_index) / (1.0 + refraction_index);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}
