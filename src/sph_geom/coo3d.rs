
use std::f64::{EPSILON};
use std::f64::consts::{PI};

use super::super::Customf64;

const TWO_PI: f64 = 2.0 * PI;
const HALF_PI: f64 = 0.5 * PI;

// see https://www.nalgebra.org/


// Euclidean coordinates

pub trait Vec3 {
  
  fn new(x: f64, y: f64, z: f64) -> Self where Self: Sized;
  
  fn x(&self) -> f64;
  fn y(&self) -> f64;
  fn z(&self) -> f64;
  
  fn norm(&self) -> f64 { 
    self.squared_norm().sqrt()
  }
  
  fn squared_norm(&self) -> f64 {
    // squared_norm_of(Vec3::x(self), Vec3::y(self), Vec3::z(self))
    squared_norm_of(self.x(), self.y(), self.z())
  }
  
  fn dot_product<V: Vec3>(&self, other: &V) -> f64 {
    // Vec3::x(self) * vec3.x() + Vec3::y(self) * vec3.y() + Vec3::z(self) * vec3.z()
    self.x() * other.x() + self.y() * other.y() + self.z() * other.z()
  }
  
  fn lonlat(&self) -> (f64, f64) {
    // lonlat_of(Vec3::x(self), Vec3::y(self), Vec3::z(self))
    lonlat_of(self.x(), self.y(), self.z())
  }
  
  fn opposite(&self) -> Self where Self: Sized{
    Self::new(-self.x(), -self.y(), -self.z())
  }

  fn squared_euclidean_dist<V: Vec3>(&self, other: &V) -> f64 {
    pow2(self.x() - other.x()) + pow2(self.y() - other.y()) + pow2(self.z() - other.z())
  }

  fn euclidean_dist<V: Vec3>(&self, other: &V) -> f64 {
    self.squared_euclidean_dist(other).sqrt()
  }

  #[inline]
  fn normalized(&self) -> UnitVect3 {
    let norm = self.norm();
    UnitVect3 {
      x: self.x() / norm,
      y: self.y() / norm,
      z: self.z() / norm,
     }
  }

  /*#[inline]
  fn normalized_opposite(&self) -> UnitVect3 {
    let norm = self.norm();
    UnitVect3 {
      x: -self.x() / norm,
      y: -self.y() / norm,
      z: -self.z() / norm,
    }
  }*/
}

// Apply to all references to a type that implements the Vec3 trait
impl<'a, T> Vec3 for &'a T where T: Vec3 {
  
  fn new(_x: f64, _y: f64, _z: f64) -> Self { // Voir si ca marche en pratique
    panic!("Method must be defined for each implementor!");
  }
  
  #[inline]
  fn x(&self) -> f64 { Vec3::x(*self) }

  #[inline]
  fn y(&self) -> f64 { Vec3::y(*self) }

  #[inline]
  fn z(&self) -> f64 { Vec3::z(*self) }
}

// impl<'a, T> Vec3 for &'a mut T where T: Vec3 {}




pub trait UnitVec3: Vec3 {

  /*#[inline]
  fn check_is_unit(&self) {
    assert!(UnitVect3D::is_unit_from_squared_norm(self.squared_norm()));
  }*/
  
  fn cross_prod_norm<T: Vec3 + UnitVec3>(&self, other: &T) -> f64 {
    let nx = self.y() * other.z() - self.z() * other.y();
    let ny = self.z() * other.x() - self.x() * other.z();
    let nz = self.x() * other.y() - self.y() * other.x();
    (nx.pow2() + ny.pow2() + nz.pow2()).sqrt()
  }
  
  /// Compute the angular distance between this vector and the other given vector
  fn ang_dist<T: Vec3 + UnitVec3>(&self, other: &T) -> f64 {
    let cos = self.dot_product(other);
    let sin = self.cross_prod_norm(other);
    debug_assert!(sin >= 0.0);
    sin.atan2(cos)
    /* As noticed by M. Reinecke, the folowing formula is numerically unstable for angles near PI.
    let half_eucl = 0.5 * self.euclidean_dist(other);
    2.0 * half_eucl.asin()*/
    // One can use also use the Vincenty formula from (lon_a, lat_a), (lon_b, lat_b)
  }
  
  /// Returns the angular
  fn arc_center<T: Vec3 + UnitVec3>(&self, other: &T) -> UnitVect3 {
    let one_over_twice_norm = 1.0 / (1.0 + self.dot_product(other));
    if one_over_twice_norm.is_infinite() {
      UnitVect3 {x: 1.0, y: 0.0, z: 0.0} // any Unit vector is ok
    } else {
      // Check for numerical inaccuracy in cases where one_over_twice_norm 
      // is very large but not infinite?
      UnitVect3 {
        x: one_over_twice_norm * (self.x() + other.x()),
        y: one_over_twice_norm * (self.y() + other.y()),
        z: one_over_twice_norm * (self.z() + other.z()),
      }
    }
  }
  
  fn to_struct(&self) -> UnitVect3 {
    UnitVect3{x: self.x(), y: self.y(), z: self.z() }
  }

  /*fn to_opposite(&self) -> UnitVect3 {
    UnitVect3{x: -self.x(), y: -self.y(), z: -self.z() }
  }*/
}

#[inline]
fn pow2(x: f64) -> f64 {
  x * x
}

/* Commented beacause not used so far
#[inline]
pub fn norm_of(x: f64, y: f64, z: f64) -> f64 {
  squared_norm_of(x, y, z).sqrt()
}
*/

#[inline]
pub fn squared_norm_of(x: f64, y: f64, z: f64) -> f64 {
  x.pow2() + y.pow2() + z.pow2()
}

/* Commented beacause not used so far
#[inline]
pub fn is_unit(x: f64, y: f64, z: f64) -> bool {
  is_unit_from_squared_norm(x.pow2() + y.pow2() + z.pow2())
}*/

#[inline]
pub fn is_unit_from_norm(norm: f64) -> bool {
  (norm - 1.0_f64).abs() <= EPSILON
}

#[inline]
pub fn is_unit_from_squared_norm(squared_norm: f64) -> bool {
  is_unit_from_norm(squared_norm)
}


/*#[inline]
pub fn dot_product(v1: &Vec3, v2: &Vec3) -> f64 {
  v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z()
}*/

#[inline]
pub fn dot_product<T1, T2>(v1: &T1, v2: &T2) -> f64 
  where T1: Vec3, T2: Vec3 {
  v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z()
}

#[inline]
pub fn cross_product<T1, T2>(v1: T1, v2: T2) -> Vect3
  where T1: Vec3, T2: Vec3 {
  Vect3::new(
    v1.y() * v2.z() - v1.z() * v2.y(),
    v1.z() * v2.x() - v1.x() * v2.z(),
    v1.x() * v2.y() - v1.y() * v2.x()
  )
}

/* Commented beacause not used so far
#[inline]
pub fn cross_product_opt<T1, T2>(o1: Option<T1>, o2: Option<T2>) -> Option<Vect3>
  where T1: Vec3, T2: Vec3 {
  if o1.is_none() || o2.is_none() {
    None
  } else {
    Some(cross_product(o1.unwrap(), o2.unwrap()))
  }
}*/

/*#[inline]
WE CAN COMPUTE THE coo OF NON UNIT VECTS!!
pub fn lonlat_of(mut x: f64, mut y: f64, mut z: f64) -> (f64, f64) {
  let squared_norm = squared_norm_of(x, y, z);
  if !is_unit_from_squared_norm(squared_norm) {
    let norm = squared_norm.sqrt();
    x /= norm;
    y /= norm;
    z /= norm;
  }
  lonlat_of_unsafe(x, y, z)
}*/

#[inline]
pub fn lonlat_of(x: f64, y: f64, z: f64) -> (f64, f64) {
  let mut lon = y.atan2(x);
  if lon < 0.0_f64 {
    lon += TWO_PI;
  } else if lon == TWO_PI {
    lon = 0.0;
  }
  let lat = z.atan2((x.pow2() + y.pow2()).sqrt());
  debug_assert!(0.0 <= lon && lon <= TWO_PI);
  debug_assert!(-HALF_PI <= lat && lat < HALF_PI);
  (lon, lat)
}

#[derive(Debug)]
pub struct Vect3 {
  x: f64,
  y: f64,
  z: f64,
}

/*impl Vect3 {
  pub fn new(x: f64, y: f64, z: f64) -> Vect3 {
    Vect3{ x, y, z }
  }
}*/

impl Vec3 for Vect3 {
  #[inline]
  fn new(x: f64, y: f64, z: f64) -> Vect3 { // Self where Self: Sized;
    Vect3{ x, y, z }
  }
  #[inline]
  fn x(&self) -> f64 { self.x }

  #[inline]
  fn y(&self) -> f64 { self.y }

  #[inline]
  fn z(&self) -> f64 { self.z }
}

// pub struct UnitVect3 (f64, f64, f64);
#[derive(Debug)]
pub struct UnitVect3 {
  x: f64,
  y: f64,
  z: f64,
}

impl UnitVect3 {
  /*pub fn new(x: f64, y: f64, z: f64) -> UnitVect3 {
    let norm2 = squared_norm_of(x, y, z);
    if is_unit_from_squared_norm(norm2) {
      UnitVect3::new_unsafe(x, y, z)
    } else {
      let norm = norm2.sqrt();
      UnitVect3::new_unsafe(x / norm, y / norm, z / norm)
    }
    
  }*/
  
  #[inline]
  pub fn new_unsafe(x: f64, y: f64, z: f64) -> UnitVect3 {
    UnitVect3{ x, y, z }
  }
  
  #[inline]
  pub fn lonlat(&self) -> LonLat {
    let mut lon = f64::atan2(self.y(), self.x());
    if lon < 0.0_f64 {
      lon += TWO_PI;
    }
    let lat = f64::atan2(self.z(), (self.x.pow2() + self.y.pow2()).sqrt());
    LonLat {lon, lat}
  }
}

impl Vec3 for UnitVect3 {
  
  fn new(x: f64, y: f64, z: f64) -> UnitVect3 {
    let norm2 = squared_norm_of(x, y, z);
    if is_unit_from_squared_norm(norm2) {
      UnitVect3::new_unsafe(x, y, z)
    } else {
      let norm = norm2.sqrt();
      UnitVect3::new_unsafe(x / norm, y / norm, z / norm)
    }
  }
  
  #[inline]
  fn x(&self) -> f64 { self.x }

  #[inline]
  fn y(&self) -> f64 { self.y }

  #[inline]
  fn z(&self) -> f64 { self.z }
}


impl UnitVec3 for UnitVect3 {

}

// Geographic coordinates

pub trait LonLatT {
  fn lon(&self) -> f64;
  fn lat(&self) -> f64;
  fn vec3(&self) -> UnitVect3 {
    vec3_of(self.lon(), self.lat())
  }
}

// Apply to all references to a type that implements the Vec3 trait
impl<'a, T> LonLatT for &'a T where T: LonLatT {
  #[inline]
  fn lon(&self) -> f64 { LonLatT::lon(*self) }

  #[inline]
  fn lat(&self) -> f64  { LonLatT::lat(*self) }

  #[inline]
  fn vec3(&self) -> UnitVect3 {
    LonLatT::vec3(*self)
  }
}

// pub const fn vec3_of(lon: f64, lat: f64) -> UnitVect3 {
pub fn vec3_of(lon: f64, lat: f64) -> UnitVect3 {
  let (sin_lon, cos_lon) = lon.sin_cos();
  let (sin_lat, cos_lat) = lat.sin_cos();
  UnitVect3 {
    x: cos_lat * cos_lon,
    y: cos_lat * sin_lon,
    z: sin_lat,
  }
} 

#[derive(Debug)]
pub struct LonLat {
  pub lon: f64, 
  pub lat: f64,
}

impl LonLatT for LonLat {
  #[inline]
  fn lon(&self) -> f64 {
    self.lon
  }
  #[inline]
  fn lat(&self) -> f64 {
    self.lat
  }
}

// Specific Coo3D

#[derive(Debug)]
pub struct Coo3D {
  x: f64,
  y: f64,
  z: f64,
  lon: f64,
  lat: f64,
}

impl Coo3D {

  pub fn from<T: UnitVec3>(v: T) -> Coo3D {
    Coo3D::from_vec3(v.x(), v.y(), v.z())
  }
  
  pub fn from_vec3(x: f64, y: f64, z: f64) -> Coo3D {
    let (lon, lat) = lonlat_of(x, y, z);
    Coo3D {x, y, z, lon, lat}
  }
  
  /// lon and lat in radians
  pub fn from_sph_coo(lon: f64, lat: f64) -> Coo3D {
    let v = vec3_of(lon, lat);
    if lon < 0.0 || TWO_PI <= lon || lat < -HALF_PI || HALF_PI < lat {
      let (new_lon, new_lat) = lonlat_of(v.x(), v.y(), v.z());
      Coo3D {x: v.x(), y: v.y(), z: v.z(), lon: new_lon, lat: new_lat}
    } else {
      Coo3D {x: v.x(), y: v.y(), z: v.z(), lon, lat}
    }
  }
  
}

impl LonLatT for Coo3D {

  #[inline]
  fn lon(&self) -> f64 { self.lon }

  #[inline]
  fn lat(&self) -> f64 { self.lat }
}

impl Vec3 for Coo3D {
  
  #[inline]
  fn new(x: f64, y: f64, z: f64) -> Coo3D {
    Coo3D::from_vec3(x, y, z)
  }
  
  #[inline]
  fn x(&self) -> f64 { self.x }
  
  #[inline]
  fn y(&self) -> f64 { self.y }
  
  #[inline]
  fn z(&self) -> f64 { self.z }
}

impl UnitVec3 for Coo3D { }


