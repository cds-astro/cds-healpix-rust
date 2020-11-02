// use rand::{thread_rng, Rng};
use crate::SMALLER_EDGE2OPEDGE_DIST;
use super::coo3d::*;

static R_MAX: f64 = SMALLER_EDGE2OPEDGE_DIST[0];

#[derive(Debug)]
pub struct Cone {
  center: UnitVect3,
  radius: f64, // in radians
}

impl Cone {

  pub fn new(center: UnitVect3, radius: f64) -> Cone {
    Cone { center, radius }
  }

  /// Returns the smallest cone enclosing the given input points.
  /// Retuns `None` if the radius is > to ~48.7 deg (since for a higher value we have to test the 
  /// 12 base cells). We use this since the algorithm is not made to work with nonreflex cones (i.e.
  /// if input points are distributed on more than an hemisphere).
  /// For a more general algorithm working with reflex cones, one should have a look at the 
  /// very similar "minimum enclosing sphere" problem (which is probably more efficient since
  /// the code here use angular distances and thus computes trigonometric functions).
  /// We use the algorithm of 
  /// Barequet & Elber (2005) "Optimal bounding cones of vectors in three dimensions",
  /// [see here](https://www.sciencedirect.com/science/article/pii/S0020019004002911?via%3Dihub).
  pub fn mec<T: Vec3 + UnitVec3>(points: &[T]) -> Option<Cone> {
    match points.len() {
      0 | 1 => None,
      2 => {
        let cone = mec_2(&points[0], &points[1]);
        if cone.radius > R_MAX { None } else { Some(cone) }
      },
      3 => {
        let cone = mec_3(&points[0], &points[1], &points[2]);
        if cone.radius > R_MAX { None } else { Some(cone) }
      },
      _ => mec_n(points),
    }
  }

  /// Returns a bounding cone (not the smalles one), i.e. a cone containing all the given points.
  pub fn bounding_cone<T: UnitVec3>(points: &[T]) -> Cone {
    // Compute the gravity center direction
    let (x, y, z) = points.iter().fold((0.0, 0.0, 0.0), |(x, y, z), p| {
      (x + p.x(), y + p.y(), z + p.z())
    });
    let norm = (pow2(x) + pow2(y) + pow2(z)).sqrt();
    let center = UnitVect3::new(x / norm, y / norm, z / norm);
    // Compute dmax
    let d2_max = points.iter()
      .map(|p| center.squared_euclidean_dist(p))
      .fold(0.0, f64::max);
    Cone::new(center, 2.0 * (0.5 * d2_max.sqrt()).asin())
  }

  pub fn radius(&self) -> f64 {
    self.radius
  }

  pub fn center(&self) -> &UnitVect3 {
    &self.center
  }

  /// Returns `true` if the given point is inside the cone.
  pub fn contains<T: Vec3 + UnitVec3>(&self, point: &T) -> bool {
    self.center.ang_dist(point) <= self.radius
  }

}


/// Returns the minimum enclosing cone, i.e. the cone containig the two given points and having
/// the smallest possible radius. In this trivial case, the diameter of the cone is the arc (ab).
#[inline]
fn mec_2<T: Vec3 + UnitVec3>(a: &T, b: &T) -> Cone {
  let radius = 0.5 * a.ang_dist(b);
  let center = a.arc_center(b);
  Cone::new(center, radius)
}


/// Returns the Minimum Enclosing Cone, i.e. the cone containing the three given points and having
/// the smallest possible radius.
pub(crate) fn mec_3<T: Vec3 + UnitVec3>(a: &T, b: &T, c: &T) -> Cone {
  let da = b.ang_dist(c);
  let db = a.ang_dist(c);
  let dc = a.ang_dist(b);
  let (mut center, mut radius) = if da > db && da > dc {
    (b.arc_center(c), 0.5 * da)
  } else if db > dc {
    (a.arc_center(c), 0.5 * db)
  } else {
    (a.arc_center(b), 0.5 * dc)
  };
  if c.ang_dist(&center) > radius {
    radius = circum_radius(da, db, dc);
    center = circum_center(a, b, c, radius);
  }
  Cone::new(center, radius)
}


// We remove the dependance to rand
fn mec_n<T: Vec3 + UnitVec3>(points: &[T]) -> Option<Cone> {
  // let mut rng = thread_rng();
  let points: Vec<&T> = points.iter().map(|v| v).collect();
  // rng.shuffle(points[0..points.len()]);
  let mut cone = mec_2(points[0], points[1]);
  if cone.radius > R_MAX { return None; }
  for i in 2..points.len() {
    if !cone.contains(points[i]) {
      // rng.shuffle(points[0..i]); // try without this and compare performances
      cone = mec_2(points[0], points[i]);
      if cone.radius > R_MAX { return None; }
      for j in 1..i {
        if !cone.contains(points[j]) {
          // rng.shuffle(points[0..j]); // try without this and compare performances
          cone = mec_2(points[j], points[i]);
          if cone.radius > R_MAX { return None; }
          for k in 0..j {
            if !cone.contains(points[k]) {
              cone = mec_3(points[k], points[j], points[i]);
              if cone.radius > R_MAX { return None; }
            }
          }
        }
      }
    }
  }
  return Some(cone);
}

/// Returns the angular radius (in radians) of the circumcircle of a
/// spherical triangle of given side lengths a, b and c.
/// # Inputs:
/// - `a` first size length (in radians)
/// - `b` second size length (in radians)
/// - `c` third size length (in radians)
/// # Output: 
/// - the circumcircle angular radius (in radians)
#[inline]
fn circum_radius(a: f64, b: f64, c: f64) -> f64 {
  let sin_half_a = (0.5 * a).sin();
  let sin_half_b = (0.5 * b).sin();
  let sin_half_c = (0.5 * c).sin();
  let n = pow2(sin_half_a * sin_half_b * sin_half_c);
  let d = 0.25
    * (sin_half_a + sin_half_b + sin_half_c)
    * (sin_half_a + sin_half_b - sin_half_c)
    * (sin_half_a - sin_half_b + sin_half_c)
    * (sin_half_b + sin_half_c - sin_half_a);
  (n / d).sqrt().asin()
}

/// Computes the center on the unit sphere of the circumcircle of radius r of a spherical triangle
/// of given vertices a, b and c.
/// # Inputs:
/// - `a` first vertex
/// - `b` second vertex
/// - `c` third vertex
/// - `radius` spherical radius of the circumcircle
/// # Output: 
/// - the center on the unit sphere of the circumcircle of radius r of a spherical triangle of given
///   vertices a, b and c.
fn circum_center<T: Vec3 + UnitVec3>(a: &T, b: &T, c: &T, radius: f64) -> UnitVect3 {
  let e = 1.0 - 0.5 * pow2(radius);
  // Simple cramer resolution of AX = E (here A --> X)
  let d = a.x() * (b.y() * c.z() - c.y() * b.z())
    - b.x() * (a.y() * c.z() - c.y() * a.z())
    + c.x() * (a.y() * b.z() - b.y() * a.z());
  let x = b.y() * c.z() - c.y() * b.z()
    - a.y() * c.z() + c.y() * a.z()
    + a.y() * b.z() - b.y() * a.z();
  let y = a.x() * (c.z() - b.z()) - b.x() * (c.z() - a.z()) + c.x() * (b.z() - a.z());
  let z = a.x() * (b.y() - c.y()) - b.x() * (a.y() - c.y()) + c.x() * (a.y() - b.y());
  UnitVect3::new((e * x) / d, (e * y) / d, (e * z) / d)
}


#[inline]
fn pow2(x: f64) -> f64 {
  x * x
}
