use super::{Axes, Axis, Point3, Ray, Rotate, UnitVec3, Vec3};
use rand::Rng;
use std::f64::consts::PI;
use std::ops::Add;

pub trait Intersection {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64, time: f64) -> Option<IntersectionPoint>;
}

impl<I: Intersection> Intersection for Vec<I> {
    fn intersect(
        &self,
        ray: &Ray,
        tmin: f64,
        mut tmax: f64,
        time: f64,
    ) -> Option<IntersectionPoint> {
        let mut intersection_point = None;
        for i in self.iter() {
            if let Some(ip) = i.intersect(ray, tmin, tmax, time) {
                tmax = ip.t;
                intersection_point = Some(ip);
            }
        }

        intersection_point
    }
}

pub trait Boundable {
    fn bound(&self) -> Option<BoundingBox>;
}

impl<T: Boundable> Boundable for Vec<T> {
    fn bound(&self) -> Option<BoundingBox> {
        if self.is_empty() {
            return None;
        }

        let mut result = None;
        for shape in self.iter() {
            match (result, shape.bound()) {
                (Some(old_box), Some(new_box)) => result = Some(old_box + new_box),
                (None, Some(new_box)) => result = Some(new_box),
                (_, None) => return None,
            }
        }

        result
    }
}

#[derive(Clone, Debug)]
pub enum Shape {
    Sphere {
        center: Point3,
        radius: f64,
    },

    Plane {
        normal: UnitVec3,
        offset: f64,
    },

    MovingSphere {
        ray: Ray,
        radius: f64,
        start_time: f64,
        end_time: f64,
    },

    Rectangle {
        axes: Axes,
        lower_left: (f64, f64),
        upper_right: (f64, f64),
        height: f64,
    },

    ConstantMedium {
        boundary: Box<Shape>,
        density: f64,
    },

    Box {
        min: Point3,
        max: Point3,
        sides: Vec<Shape>,
    },

    Flipped(Box<Shape>),

    Translated {
        inner: Box<Shape>,
        offset: Vec3,
    },

    Rotated {
        axis: Axis,
        inner: Box<Shape>,
        angle: f64,
    },
}

impl Shape {
    pub fn flipped(self) -> Self {
        Self::Flipped(Box::new(self))
    }

    pub fn translated(self, offset: Vec3) -> Self {
        Self::Translated {
            inner: Box::new(self),
            offset,
        }
    }

    pub fn rotated(self, axis: Axis, angle: f64) -> Self {
        self.rotate(axis, angle)
    }

    pub fn rectangle(lower_left: Point3, upper_right: Point3) -> Self {
        let axes = if lower_left[0] == upper_right[0] {
            Axes::YZ
        } else if lower_left[1] == upper_right[1] {
            Axes::XZ
        } else {
            Axes::XY
        };

        let (p1, p2, o) = axes.coords();

        let height = lower_left[o];
        let lower_left = (lower_left[p1], lower_left[p2]);
        let upper_right = (upper_right[p1], upper_right[p2]);

        Self::Rectangle {
            lower_left,
            upper_right,
            axes,
            height,
        }
    }

    pub fn new_box(min: Point3, max: Point3) -> Self {
        let mut sides = Vec::new();

        let Point3([x0, y0, z0]) = min;
        let Point3([x1, y1, z1]) = max;

        sides.push(
            Shape::Rectangle {
                height: x0,
                lower_left: (y0, z0),
                upper_right: (y1, z1),
                axes: Axes::YZ,
            }
            .flipped(),
        );

        sides.push(Shape::Rectangle {
            height: x1,
            lower_left: (y0, z0),
            upper_right: (y1, z1),
            axes: Axes::YZ,
        });

        sides.push(
            Shape::Rectangle {
                height: y0,
                lower_left: (x0, z0),
                upper_right: (x1, z1),
                axes: Axes::XZ,
            }
            .flipped(),
        );

        sides.push(Shape::Rectangle {
            height: y1,
            lower_left: (x0, z0),
            upper_right: (x1, z1),
            axes: Axes::XZ,
        });

        sides.push(
            Shape::Rectangle {
                height: z0,
                lower_left: (x0, y0),
                upper_right: (x1, y1),
                axes: Axes::XY,
            }
            .flipped(),
        );

        sides.push(Shape::Rectangle {
            height: z1,
            lower_left: (x0, y0),
            upper_right: (x1, y1),
            axes: Axes::XY,
        });

        Self::Box { min, max, sides }
    }

    pub fn sphere(
        radius: f64,
        (start_time, start_center): (f64, Point3),
        (end_time, end_center): (f64, Point3),
    ) -> Self {
        assert!(end_time > start_time);

        if start_center == end_center {
            Self::Sphere {
                center: start_center,
                radius,
            }
        } else {
            Self::MovingSphere {
                start_time,
                end_time,
                radius,
                ray: Ray {
                    origin: start_center,
                    direction: end_center - start_center,
                },
            }
        }
    }
}

impl Intersection for Shape {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64, time: f64) -> Option<IntersectionPoint> {
        match self {
            Shape::Sphere { center, radius } => {
                let Ray {
                    origin: o,
                    direction: dir,
                } = ray;
                let oc = *o - *center;
                let a = dir.norm_squared();
                let half_b = oc.dot(*dir);
                let c = oc.norm_squared() - radius * radius;
                let discriminant = half_b * half_b - a * c;

                if discriminant < 0.0 {
                    None
                } else {
                    let t = (-half_b - discriminant.sqrt()) / a;
                    if t > tmin && t < tmax {
                        let point = ray.at(t);
                        let normal = (point - *center).normed();
                        let (u, v) = {
                            let p = (point - *center) / *radius;
                            let phi = p[2].atan2(p[0]);
                            let theta = p[1].asin();
                            let u = 1.0 - (phi + PI) / (2.0 * PI);
                            let v = (theta + PI / 2.0) / PI;
                            (u, v)
                        };
                        Some(IntersectionPoint {
                            t,
                            point,
                            normal,
                            in_vec: dir.normed(),
                            surface_coordinates: (u, v),
                        })
                    } else {
                        None
                    }
                }
            }

            Shape::Plane { normal, offset } => {
                if ray.direction.dot(**normal).abs() <= 1e-10 {
                    None
                } else {
                    let u = Vec3(ray.origin.0);
                    let t = (offset - u.dot(**normal)) / ray.direction.dot(**normal);

                    if t >= tmin && t < tmax {
                        let intersection_point = IntersectionPoint {
                            t,
                            in_vec: ray.direction.normed(),
                            normal: *normal,
                            point: ray.at(t),
                            surface_coordinates: (0.0, 0.0),
                        };
                        Some(intersection_point)
                    } else {
                        None
                    }
                }
            }

            Shape::MovingSphere {
                ray: r,
                radius,
                start_time,
                ..
            } => {
                let sphere = Shape::Sphere {
                    center: r.at(time - start_time),
                    radius: *radius,
                };

                sphere.intersect(ray, tmin, tmax, time)
            }

            Shape::Rectangle {
                lower_left: (x0, y0),
                upper_right: (x1, y1),
                height,
                axes,
            } => {
                let (p1, p2, o) = axes.coords();

                let t = (height - ray.origin[o]) / ray.direction[o];
                if t < tmin || t > tmax {
                    return None;
                }

                let point = ray.at(t);
                if point[p1] < *x0 || point[p1] > *x1 || point[p2] < *y0 || point[p2] > *y1 {
                    return None;
                }

                let u = (point[p1] - x0) / (x1 - x0);
                let v = (point[p2] - y0) / (y1 - y0);

                let mut normal = Vec3([0.0, 0.0, 0.0]);
                normal[o] = 1.0;

                Some(IntersectionPoint {
                    point,
                    surface_coordinates: (u, v),
                    normal: normal.normed(),
                    t,
                    in_vec: ray.direction.normed(),
                })
            }

            Self::Flipped(inner) => {
                inner
                    .intersect(ray, tmin, tmax, time)
                    .map(|ip| IntersectionPoint {
                        normal: -ip.normal,
                        ..ip
                    })
            }

            Self::Box { sides, .. } => sides.intersect(ray, tmin, tmax, time),
            Self::Translated { inner, offset } => {
                let moved_ray = Ray {
                    origin: ray.origin - *offset,
                    ..*ray
                };

                inner
                    .intersect(&moved_ray, tmin, tmax, time)
                    .map(|ip| IntersectionPoint {
                        point: ip.point + *offset,
                        ..ip
                    })
            }

            Self::Rotated { axis, angle, inner } => {
                let rotated_ray = Ray {
                    origin: ray.origin.rotate(*axis, -angle),
                    direction: ray.direction.rotate(*axis, -angle),
                };

                inner
                    .intersect(&rotated_ray, tmin, tmax, time)
                    .map(|ip| IntersectionPoint {
                        point: ip.point.rotate(*axis, *angle),
                        normal: ip.normal.rotate(*axis, *angle),
                        in_vec: ip.in_vec.rotate(*axis, *angle),
                        ..ip
                    })
            }

            Self::ConstantMedium { boundary, density } => boundary
                .intersect(ray, f64::NEG_INFINITY, f64::INFINITY, time)
                .and_then(|mut ip1| {
                    boundary
                        .intersect(ray, ip1.t + 0.001, f64::INFINITY, time)
                        .and_then(|mut ip2| {
                            let mut rng = rand::thread_rng();
                            ip1.t = ip1.t.max(tmin);
                            ip2.t = ip2.t.min(tmax);

                            if ip1.t >= ip2.t {
                                return None;
                            }

                            ip1.t = ip1.t.max(0.0);

                            let length = ray.direction.norm();
                            let distance_inside_boundary = (ip2.t - ip1.t) * length;
                            let hit_distance = -1.0 / density * rng.gen::<f64>().ln();

                            if hit_distance > distance_inside_boundary {
                                return None;
                            }

                            let t = ip1.t + hit_distance / length;
                            let point = ray.at(t);

                            let normal = Vec3([1.0, 0.0, 0.0]).normed();

                            Some(IntersectionPoint {
                                point,
                                normal,
                                t,
                                in_vec: ray.direction.normed(),
                                surface_coordinates: (0.0, 0.0),
                            })
                        })
                }),
        }
    }
}

impl Boundable for Shape {
    fn bound(&self) -> Option<BoundingBox> {
        match self {
            Shape::Plane { .. } => None,
            Shape::Sphere { center, radius } => {
                let min = *center - Vec3([1.0, 1.0, 1.0]) * *radius;
                let max = *center + Vec3([1.0, 1.0, 1.0]) * *radius;

                Some(BoundingBox { min, max })
            }
            Shape::MovingSphere {
                ray,
                radius,
                start_time,
                end_time,
            } => Some(
                (&Shape::Sphere {
                    center: ray.at(0.0),
                    radius: *radius,
                })
                    .bound()
                    .unwrap()
                    + (&Shape::Sphere {
                        center: ray.at(end_time - start_time),
                        radius: *radius,
                    })
                        .bound()
                        .unwrap(),
            ),

            Shape::Rectangle {
                lower_left: (x0, y0),
                upper_right: (x1, y1),
                height,
                axes,
            } => {
                let (p1, p2, o) = axes.coords();

                let mut min = Point3::default();
                let mut max = Point3::default();

                min[p1] = *x0;
                min[p2] = *y0;
                min[o] = height - 0.0001;

                max[p1] = *x1;
                max[p2] = *y1;
                max[o] = height + 0.0001;

                Some(BoundingBox { min, max })
            }

            Self::Flipped(inner) => inner.bound(),

            Self::Box { min, max, .. } => Some(BoundingBox {
                min: *min,
                max: *max,
            }),

            Self::Translated { inner, offset } => {
                inner.bound().map(|BoundingBox { min, max }| BoundingBox {
                    min: min + *offset,
                    max: max + *offset,
                })
            }

            Self::Rotated { inner, axis, angle } => inner.bound().map(|bb| {
                let mut min = Point3([f64::INFINITY, f64::INFINITY, f64::INFINITY]);
                let mut max = Point3([f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY]);

                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..2 {
                            let x = i as f64 * bb.max[0] + (1.0 - i as f64) * bb.min[0];
                            let y = j as f64 * bb.max[1] + (1.0 - j as f64) * bb.min[1];
                            let z = k as f64 * bb.max[2] + (1.0 - k as f64) * bb.min[2];

                            let tester = Vec3([x, y, z]).rotate(*axis, *angle);

                            for c in 0..3 {
                                min[c] = min[c].min(tester[c]);
                                max[c] = max[c].max(tester[c]);
                            }
                        }
                    }
                }

                BoundingBox { min, max }
            }),

            Self::ConstantMedium { boundary, .. } => boundary.bound(),
        }
    }
}

impl Rotate for Shape {
    fn rotate(self, axis: Axis, angle: f64) -> Self {
        Self::Rotated {
            inner: Box::new(self),
            angle,
            axis,
        }
    }
}

/// Contains information about the intersection of a ray
/// with an object: the actual point, the normal vector
/// at that point, and the parameter of the ray.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    pub point: Point3,
    pub normal: UnitVec3,
    pub t: f64,
    pub in_vec: UnitVec3,
    pub surface_coordinates: (f64, f64),
}

impl IntersectionPoint {
    pub fn face(&self) -> Face {
        if self.in_vec.dot(*self.normal) < 0.0 {
            Face::Front
        } else {
            Face::Back
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Face {
    Front,
    Back,
}

#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    pub min: Point3,
    pub max: Point3,
}

impl BoundingBox {
    pub fn hit(&self, ray: &Ray, mut tmin: f64, mut tmax: f64) -> bool {
        for i in 0..3 {
            let d = ray.direction[i];
            let o = ray.origin[i];
            let mut t0 = (self.min[i] - o) / d;
            let mut t1 = (self.max[i] - o) / d;
            if d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            tmin = t0.max(tmin);
            tmax = t1.min(tmax);

            if tmax <= tmin {
                return false;
            }
        }

        true
    }
}

impl Add<BoundingBox> for BoundingBox {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let new_min = Point3([
            self.min[0].min(other.min[0]),
            self.min[1].min(other.min[1]),
            self.min[2].min(other.min[2]),
        ]);

        let new_max = Point3([
            self.max[0].max(other.max[0]),
            self.max[1].max(other.max[1]),
            self.max[2].max(other.max[2]),
        ]);

        Self {
            min: new_min,
            max: new_max,
        }
    }
}
