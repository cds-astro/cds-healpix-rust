//! This module contains elements of linear algebra and spherical geometry needed to solve
//! astrometrical problems.
//!
use std::ops::{Add, Sub};

use crate::Customf64;

pub trait HasLocalRotMatrix {
  fn local_rot_matrix(&self) -> M3x3;
}

#[derive(PartialEq, Debug)]
pub struct Coo {
  /// longitude in `[0, 2pi[`, in `rad`
  pub lon: f64, // usual units: deg
  /// latitude in `[-pi/2, pi/2]`, in `rad`
  pub lat: f64, // usual units: deg
}

impl Coo {
  pub fn new(lon: f64, lat: f64) -> Self {
    Self { lon, lat }
  }

  pub fn from_deg(lon: f64, lat: f64) -> Self {
    Coo {
      lon: lon.to_radians(),
      lat: lat.to_radians(),
    }
  }

  pub fn from_xyz(xyz: &XYZ) -> Self {
    Self::from(xyz.0, xyz.1, xyz.2)
  }

  /// # Output
  /// * `(lon, lat)`, in degrees
  pub fn to_deg(&self) -> (f64, f64) {
    (self.lon.to_degrees(), self.lat.to_degrees())
  }

  /// Create spherical coordinates from cartesian coordinates.
  /// # Remark
  /// The cartesian vector do not have to be normalized.
  pub fn from(x: f64, y: f64, z: f64) -> Self {
    const TWICE_PI: f64 = 2.0 * std::f64::consts::PI;
    let mut lon = y.atan2(x);
    if lon < 0.0_f64 {
      lon += TWICE_PI;
    } else if lon == TWICE_PI {
      lon = 0.0;
    }
    let lat = z.atan2((x.pow2() + y.pow2()).sqrt());
    Coo { lon, lat }
  }

  /// Provides the vector $\vec{r0} = (x, y, z)$, see Eq. 1.2.11 in
  /// [The Hipparcos and Tycho Catalogues](https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf)
  /// in radians.
  /// Remark: the vector is normalized (unit vector)
  pub fn xyz(&self) -> XYZ {
    XYZ::from_coo(self.lon, self.lat)
  }
}

impl HasLocalRotMatrix for Coo {
  /// Components of the 3x3 rotation matrix transforming a vector into the
  /// reference frame to a vector into the local frame (i.e. the frame in
  /// which the position of the projection center is (1, 0, 0).
  ///
  /// The 3 transposed vector composing this matrix correspond to
  /// $\vec{r_0}$, $\vec{p_0}$, $\vec{q_0}$ of Eq. 1.2.16 in
  /// [The Hipparcos and Tycho Catalogues](https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf)
  ///
  /// # Remark 1:
  /// * $r_{22} =  \cos(lon)$
  /// * $r_{21} = -\sin(lon)$
  /// * $r_{33} =  \cos(lat)$
  /// * $r_{13} =  \sin(lat)$
  ///
  /// # Remark 2
  /// If you need to compute the `xyz` coordinates of the vector then you
  /// may choose to compute the matrix from the Euclidean coordinates not
  /// to have to call twice the trigonometric functions (perform first a
  /// bench to see if the compiler is not optimizing by caching the trigonometric
  /// function results the code with Euclidean coordinate is more complex since
  /// division by 0 have to be tested).
  fn local_rot_matrix(&self) -> M3x3 {
    let (sin_lon, cos_lon) = self.lon.sin_cos();
    let (sin_lat, cos_lat) = self.lat.sin_cos();
    M3x3::new(
      XYZt::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat), // r0
      XYZt::new(-sin_lon, cos_lon, 0.0),                        // p0
      XYZt::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat), // q0
    )
  }
}

/// Can be called on array ref
/// ```rust
/// use cdshealpix::nested::map::astrometry::math::scalar;
///
/// let a1: [f64; 3] = [1.0, 2.0, 3.0];
/// let a2: [f64; 3] = [2.0, 3.0, 4.0];
/// assert_eq!(scalar(&a1, &a2), 20.0);
/// ```
pub fn scalar<'a, A, B>(lhs: A, rhs: B) -> f64
where
  A: IntoIterator<Item = &'a f64>,
  B: IntoIterator<Item = &'a f64>,
{
  lhs
    .into_iter()
    .zip(rhs.into_iter())
    .map(|(a, b)| a * b)
    .sum()
}

pub trait HasXYZ {
  fn x(&self) -> f64;
  fn y(&self) -> f64;
  fn z(&self) -> f64;

  fn to_coo(&self) -> Coo {
    Coo::from(self.x(), self.y(), self.z())
  }

  fn norm(&self) -> f64 {
    self.squared_norm().sqrt()
  }

  fn squared_norm(&self) -> f64 {
    self.x().pow2() + self.y().pow2() + self.z().pow2()
  }

  fn scalar<T: HasXYZ>(&self, rhs: &T) -> f64 {
    self.x() * rhs.x() + self.y() * rhs.y() + self.z() * rhs.z()
  }

  /// Scalar product of the y and z coordinates ignoring the x coordinate).
  /// Used when we know in advance that either `self.x()` or `rhs.x()` equals 0
  fn scalar_yz<T: HasXYZ>(&self, rhs: &T) -> f64 {
    self.y() * rhs.y() + self.z() * rhs.z()
  }
}

/// 3 dimentional vector
#[derive(Debug, Clone)]
pub struct XYZ(pub f64, pub f64, pub f64);

impl XYZ {
  pub fn new(x: f64, y: f64, z: f64) -> XYZ {
    XYZ(x, y, z)
  }

  pub fn from_coo(lon: f64, lat: f64) -> XYZ {
    let (sin_lon, cos_lon) = lon.sin_cos();
    let (sin_lat, cos_lat) = lat.sin_cos();
    XYZ::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
  }

  pub fn time(mut self, cte: f64) -> Self {
    self.0 *= cte;
    self.1 *= cte;
    self.2 *= cte;
    self
  }
}

impl<T: HasXYZ> Add<T> for XYZ {
  type Output = Self;
  fn add(mut self, rhs: T) -> Self::Output {
    self.0 += rhs.x();
    self.1 += rhs.y();
    self.2 += rhs.z();
    self
  }
}

impl<T: HasXYZ> Sub<T> for XYZ {
  type Output = Self;
  fn sub(mut self, rhs: T) -> Self::Output {
    self.0 -= rhs.x();
    self.1 -= rhs.y();
    self.2 -= rhs.z();
    self
  }
}

impl HasXYZ for XYZ {
  fn x(&self) -> f64 {
    self.0
  }
  fn y(&self) -> f64 {
    self.1
  }
  fn z(&self) -> f64 {
    self.2
  }
}
impl HasXYZ for &XYZ {
  fn x(&self) -> f64 {
    self.0
  }
  fn y(&self) -> f64 {
    self.1
  }
  fn z(&self) -> f64 {
    self.2
  }
}

impl HasLocalRotMatrix for XYZ {
  fn local_rot_matrix(&self) -> M3x3 {
    let n_xy = self.x().hypot(self.y());
    let norm = n_xy.hypot(self.z());
    let (sin_lon, cos_lon, sin_lat, cos_lat) = if n_xy.eq0() {
      (0.0, 1.0, 1.0, 0.0)
    } else {
      (
        self.y() / n_xy,
        self.x() / n_xy,
        self.z() / norm,
        n_xy / norm,
      )
    };
    if norm.eq0() {
      Default::default()
    } else if (norm - 1.0).eq0() {
      M3x3(
        XYZt(self.x(), self.y(), self.z()),                    // r0
        XYZt(-sin_lon, cos_lon, 0.0),                          // p0
        XYZt(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat), // q0
      )
    } else {
      M3x3(
        XYZt(self.x() / norm, self.y() / norm, sin_lat),
        XYZt(-sin_lon, cos_lon, 0.0),
        XYZt(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat),
      )
    }
  }
}

/// Transpose of a 3 dimensional vector.
#[derive(Debug)]
pub struct XYZt(pub f64, pub f64, pub f64);

impl XYZt {
  pub fn new(x: f64, y: f64, z: f64) -> XYZt {
    XYZt(x, y, z)
  }

  pub fn from_coo(lon: f64, lat: f64) -> XYZt {
    let (sin_lon, cos_lon) = lon.sin_cos();
    let (sin_lat, cos_lat) = lat.sin_cos();
    XYZt::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
  }

  pub fn time(mut self, cte: f64) -> Self {
    self.0 *= cte;
    self.1 *= cte;
    self.2 *= cte;
    self
  }
}

impl HasXYZ for XYZt {
  fn x(&self) -> f64 {
    self.0
  }
  fn y(&self) -> f64 {
    self.1
  }
  fn z(&self) -> f64 {
    self.2
  }
}

#[derive(Debug)]
pub struct M3x3(pub XYZt, pub XYZt, pub XYZt);

/// The defalt returns the identity matrix
impl Default for M3x3 {
  fn default() -> Self {
    Self(
      XYZt(1.0, 0.0, 0.0),
      XYZt(0.0, 1.0, 0.0),
      XYZt(0.0, 0.0, 1.0),
    )
  }
}

impl M3x3 {
  pub fn new(row_x: XYZt, row_y: XYZt, row_z: XYZt) -> M3x3 {
    M3x3(row_x, row_y, row_z)
  }

  /// Rotation matrix around the `x-axis` so that the new frame is obtain from a rotation
  /// of angle `rotation_angle_rad` around the `x-axis`.
  pub fn from_rotx(rotation_angle_rad: f64) -> Self {
    let (s, c) = rotation_angle_rad.sin_cos();
    Self::new(
      XYZt(1_f64, 0_f64, 0_f64),
      XYZt(0_f64, c, -s),
      XYZt(0_f64, s, c),
    )
  }

  /// Rotation matrix around the `y-axis` so that the new frame is obtain from a rotation
  /// of angle `rotation_angle_rad` around the `y-axis`.
  pub fn from_roty(rotation_angle_rad: f64) -> Self {
    let (s, c) = rotation_angle_rad.sin_cos();
    Self::new(
      XYZt(c, 0_f64, s),
      XYZt(0_f64, 1_f64, 0_f64),
      XYZt(-s, 0_f64, c),
    )
  }

  /// Rotation matrix around the `z-axis` so that the new frame is obtain from a rotation
  /// of angle `rotation_angle_rad` around the `z-axis`.
  pub fn from_rotz(rotation_angle_rad: f64) -> Self {
    let (s, c) = rotation_angle_rad.sin_cos();
    Self::new(
      XYZt(c, -s, 0_f64),
      XYZt(s, c, 0_f64),
      XYZt(0_f64, 0_f64, 1_f64),
    )
  }

  /// Returns the rotation matrix corresponding to transform the old frame in the new frame
  /// obtained by:
  /// * a first rotation of `angle_z_rad` radians around the `z-axis`, leading to `x'y'z'`
  /// * a second rotation of `angle_yp_rad` radians around the `y'-axis`, leading to `x''y''z''`
  /// * a thrid rotation of `angle_zpp_rad` radians around the `z''-axis`, leading to `XYZ`
  #[allow(dead_code)]
  fn from_zyz_euler_intrinsic(angle_z_rad: f64, angle_yp_rad: f64, angle_zpp_rad: f64) -> Self {
    // z1-y′-z2″ (intrinsic rotations) or z2-y-z1 (extrinsic rotations)
    let (s1, c1) = angle_z_rad.sin_cos();
    let (s2, c2) = angle_yp_rad.sin_cos();
    let (s3, c3) = angle_zpp_rad.sin_cos();
    Self::new(
      XYZt(c1 * c2 * c3 - s1 * s3, -c3 * s1 - c1 * c2 * s3, c1 * s2),
      XYZt(c1 * s3 + c2 * c3 * s1, c1 * c3 - c2 * s1 * s3, s1 * s2),
      XYZt(-c3 * s2, s2 * s3, c2),
    )
  }

  pub fn xx(&self) -> f64 {
    self.row_x().x()
  }
  pub fn xy(&self) -> f64 {
    self.row_x().y()
  }
  pub fn xz(&self) -> f64 {
    self.row_x().z()
  }
  pub fn yx(&self) -> f64 {
    self.row_y().x()
  }
  pub fn yy(&self) -> f64 {
    self.row_y().y()
  }
  pub fn yz(&self) -> f64 {
    self.row_y().z()
  }
  pub fn zx(&self) -> f64 {
    self.row_z().x()
  }
  pub fn zy(&self) -> f64 {
    self.row_z().y()
  }
  pub fn zz(&self) -> f64 {
    self.row_z().z()
  }

  pub fn col_x(&self) -> XYZ {
    XYZ(self.row_x().x(), self.row_y().x(), self.row_z().x())
  }

  pub fn col_y(&self) -> XYZ {
    XYZ(self.row_x().y(), self.row_y().y(), self.row_z().y())
  }

  pub fn col_z(&self) -> XYZ {
    XYZ(self.row_x().z(), self.row_y().z(), self.row_z().z())
  }

  pub fn row_x(&self) -> &XYZt {
    &self.0
  }
  pub fn row_y(&self) -> &XYZt {
    &self.1
  }
  pub fn row_z(&self) -> &XYZt {
    &self.2
  }

  pub fn scale(mut self, cte: f64) -> M3x3 {
    self.0 = self.0.time(cte);
    self.1 = self.1.time(cte);
    self.2 = self.2.time(cte);
    self
  }

  pub fn time(&self, rhs: &M3x3) -> M3x3 {
    let rx = self.row_x();
    let ry = self.row_y();
    let rz = self.row_z();
    let cx = rhs.col_x();
    let cy = rhs.col_y();
    let cz = rhs.col_z();
    M3x3(
      XYZt(rx.scalar(&cx), rx.scalar(&cy), rx.scalar(&cz)),
      XYZt(ry.scalar(&cx), ry.scalar(&cy), ry.scalar(&cz)),
      XYZt(rz.scalar(&cx), rz.scalar(&cy), rz.scalar(&cz)),
    )
  }

  pub fn time_transpose_of(&self, rhs: &M3x3) -> M3x3 {
    let rx = self.row_x();
    let ry = self.row_y();
    let rz = self.row_z();
    let cx = rhs.row_x();
    let cy = rhs.row_y();
    let cz = rhs.row_z();
    M3x3(
      XYZt(rx.scalar(cx), rx.scalar(cy), rx.scalar(cz)),
      XYZt(ry.scalar(cx), ry.scalar(cy), ry.scalar(cz)),
      XYZt(rz.scalar(cx), rz.scalar(cy), rz.scalar(cz)),
    )
  }

  pub fn transpose_time(&self, rhs: &M3x3) -> M3x3 {
    let rx = self.col_x();
    let ry = self.col_y();
    let rz = self.col_z();
    let cx = rhs.col_x();
    let cy = rhs.col_y();
    let cz = rhs.col_z();
    M3x3(
      XYZt(rx.scalar(&cx), rx.scalar(&cy), rx.scalar(&cz)),
      XYZt(ry.scalar(&cx), ry.scalar(&cy), ry.scalar(&cz)),
      XYZt(rz.scalar(&cx), rz.scalar(&cy), rz.scalar(&cz)),
    )
  }

  /// Transform the input matrix in the global frame to a matrix in the local frame.
  pub fn rotate_matrix(&self, m_global: &M3x3) -> M3x3 {
    self.time(&m_global.time_transpose_of(&self))
  }

  /// Transform the input matrix in the local frame to a matrix in th global frame.
  pub fn unrotate_matrix(&self, m_local: &M3x3) -> M3x3 {
    self.transpose_time(&m_local.time(&self))
  }

  /// From global to local
  /// # Input
  /// - `v`: vector in the global frame
  /// # Output
  /// - vector in the local frame
  pub fn rotate(&self, v: &XYZ) -> XYZ {
    XYZ(
      self.row_x().scalar(v),
      self.row_y().scalar(v),
      self.row_z().scalar(v),
    )
  }

  /// From global to local
  /// # Input
  /// - `v`: vector in the global frame
  /// # Output
  /// - `(y, z)` the `y` and `z` coordiantes in local frame
  pub fn rotate_proj_yz(&self, v: &XYZ) -> (f64, f64) {
    (self.row_y().scalar(v), self.row_z().scalar(v))
  }

  /// From local to global
  /// <=> rotate with the transpose of the matrix
  /// # Input
  /// - `v`: vector in the local frame
  /// # Output
  /// - vector in the global frame
  pub fn unrotate(&self, v: &XYZ) -> XYZ {
    XYZ(
      self.col_x().scalar(v),
      self.col_y().scalar(v),
      self.col_z().scalar(v),
    )
  }

  pub fn to_global(&self, local_rot_matrix: &M3x3) -> M3x3 {
    local_rot_matrix.unrotate_matrix(&self)
  }
}

#[cfg(test)]
mod tests {

  use super::M3x3;

  #[test]
  fn test_from_zyz_euler_intrinsic() {
    // Angular distance between (0, 90) and (15 * (12 + 49/60), 27.4)
    /*let z1 = (15.0 * (12.0_f64 + 49.0_f64 / 60.0_f64)).to_radians();
    let y2 = -27.4_f64.to_radians();
    let z3 = 123.0_f64.to_radians();*/

    // FK4 to GAL
    let z1 = 192.25_f64.to_radians();
    let y2 = 62.6_f64.to_radians(); // (90.0_f64 - 27.4_f64) because we rotate the North pole, not he vernal point, else - 27.4 deg
    let z3 = 57_f64.to_radians(); // (180.0_f64 - 123.0_f64).to_radians();

    let m33 = M3x3::from_zyz_euler_intrinsic(z1, y2, z3);
    println!("{:?}", m33);
    println!("");
    let mz1 = M3x3::from_rotz(z1);
    let my2 = M3x3::from_roty(y2);
    let mz3 = M3x3::from_rotz(z3);
    let m33_v2 = mz1.time(&my2.time(&mz3));
    println!("{:?}", m33_v2);

    println!("----------------------------");

    // echo "cos(d2r(192.25)) * cos(d2r(62.6)) * cos(d2r(57)) - sin(d2r(192.25)) * sin(d2r(57))" | bc -l
    // -0.06698873941515072831
    /*
      -0.0669887394151507166,   0.4927284660753235682,  -0.8676008111514348650,
      -0.8727557658519926782,  -0.4503469580199613469,  -0.1883746017229204474,
      -0.4835389146321842459,   0.7445846332830310861,   0.4601997847838516236
    */

    // echo "cos(d2r(192.85948)) * cos(d2r(62.87175)) * cos(d2r(57.06808)) - sin(d2r(192.85948)) * sin(d2r(57.06808))" | bc -l

    /*
    import math
    math.cos(math.radians(192.85948)) * math.cos(math.radians(62.87175)) * math.cos(math.radians(57.06808)) - math.sin(math.radians(192.85948)) * math.sin(math.radians(57.06808))

    math.sin(math.radians(192.85948))
    -0.22256070189765312
    */

    // FK5 to GAL
    let z1 = 192.85948_f64.to_radians();
    let y2 = 62.87175_f64.to_radians(); // (90_f64 - (27.0 + 7.0 / 60.0 + 41.704/3600.0)).to_radians(); = 62.87175
    let z3 = 57.06808_f64.to_radians(); // 90.0 - 32.93192

    println!("");
    println!("{:?}", z1.sin_cos());
    println!("{:?}", y2.sin_cos());
    println!("{:?}", z3.sin_cos());
    println!("");
    println!("{}", z1.cos() * y2.cos() * z3.cos() - z1.sin() * z3.sin());

    let m33 = M3x3::from_zyz_euler_intrinsic(z1, y2, z3);
    println!("{:?}", m33);
    println!("");
    let mz1 = M3x3::from_rotz(z1);
    let my2 = M3x3::from_roty(y2);
    let mz3 = M3x3::from_rotz(z3);
    let m33_v2 = mz1.time(&my2.time(&mz3));
    println!("{:?}", m33_v2);
  }
}
