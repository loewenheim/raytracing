use super::geometry::{random_unit_vector, Face, IntersectionPoint, Ray, UnitVec3, Vec3};
use crate::textures::Texture;
use crate::RayHit;
use rand::Rng;

#[derive(Clone, Debug)]
pub enum Material {
    Lambertian { texture: Texture },

    Metal { albedo: Vec3, fuzz: f64 },

    Dielectric { albedo: Vec3, refraction_index: f64 },

    DiffuseLight { emit: Texture },
}

impl Material {
    pub fn glass() -> Self {
        Material::Dielectric {
            refraction_index: 1.5,
            albedo: Vec3([1.0, 1.0, 1.0]),
        }
    }

    pub fn scatter<R: Rng + ?Sized>(&self, p: &IntersectionPoint, rng: &mut R) -> RayHit {
        match self {
            Material::Lambertian { texture } => {
                let direction = *random_unit_vector(rng) + *p.normal;
                let scattered = Ray {
                    origin: p.point,
                    direction,
                };

                RayHit::Scattered {
                    attenuation: texture.color_at(
                        p.surface_coordinates.0,
                        p.surface_coordinates.1,
                        &p.point,
                    ),
                    ray: scattered,
                }
            }

            Material::Metal { albedo, fuzz } => {
                let reflected = reflect(&p.in_vec, &p.normal);
                let scattered = *reflected + *random_unit_vector(rng) * *fuzz;
                if scattered.dot(*p.normal) >= 0.0 {
                    RayHit::Scattered {
                        attenuation: *albedo,
                        ray: Ray {
                            origin: p.point,
                            direction: scattered,
                        },
                    }
                } else {
                    RayHit::Emitted(Vec3([0.0, 0.0, 0.0]))
                }
            }

            Material::Dielectric {
                albedo,
                refraction_index,
            } => {
                let eta_ratio = match p.face() {
                    Face::Front => 1.0 / refraction_index,
                    Face::Back => *refraction_index,
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

                RayHit::Scattered {
                    attenuation: *albedo,
                    ray: Ray {
                        origin: p.point,
                        direction,
                    },
                }
            }

            Material::DiffuseLight { emit } => RayHit::Emitted(emit.color_at(
                p.surface_coordinates.0,
                p.surface_coordinates.1,
                &p.point,
            )),
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
