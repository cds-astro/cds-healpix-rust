
use super::super::Customf64;

#[derive(Debug)]
pub struct Ellipse {
  // Derived quantities to speed up computations
  sigx2: f64,
  sigy2: f64,
  rho_sigx_sigy: f64,
  // Derived quantities
  one_over_det: f64,
}

impl Ellipse {
  
  /// Create a new ellipse for the given oriented ellipse parameters.
  /// # Inputs
  /// - `a` ellipse semi-major axis
  /// - `b` ellipse semi-minor axis
  /// - `(sin(theta), cos(theta)` in which `theta` is the counterclockwise angle between the x-axis 
  ///                             and the semi-major axis = `theta_radians.sin_cos()`
  pub fn from_oriented(a: f64, b: f64, theta_sin_cos: (f64, f64)) -> Ellipse {
    let a2 = a.pow2();
    let b2 = b.pow2();
    let (sin_t, cos_t) = theta_sin_cos;
    let sin2_t = sin_t.pow2();
    let cos2_t = cos_t.pow2();
    let sigx2 = a2 * cos2_t + b2 * sin2_t;
    let sigy2 = a2 * sin2_t + b2 * cos2_t;
    let rho_sigx_sigy = cos_t * sin_t * (a2 - b2);
    let det = sigx2 * sigy2 - rho_sigx_sigy.pow2();
    Ellipse {
      sigx2,
      sigy2,
      rho_sigx_sigy,
      one_over_det: 1.0 / det,
    }
  }

  pub fn from_cor_matrix(sig_x: f64, sig_y: f64, rho: f64) -> Ellipse {
    let sigx2 = sig_x * sig_x;
    let sigy2 = sig_y * sig_y;
    let rho_sigx_sigy = rho * sig_x * sig_y;
    let det = sigx2 * sigy2 - rho_sigx_sigy.pow2();
    Ellipse {
      sigx2,
      sigy2,
      rho_sigx_sigy,
      one_over_det: 1.0 / det,
    }
  }
  
  pub fn from_cov_matrix(sigx2: f64, sigy2: f64, rho_sigx_sigy: f64) -> Ellipse {
    let det = sigx2 * sigy2 - rho_sigx_sigy.pow2();
    Ellipse {
      sigx2,
      sigy2,
      rho_sigx_sigy,
      one_over_det: 1.0 / det,
    }
  }
  
  fn extended_by_circle(&self, sig: f64) -> Ellipse {
    Ellipse::from_cov_matrix(
      self.sigx2 + sig.pow2(),
      self.sigy2 + sig.pow2(),
      self.rho_sigx_sigy,
    )
  }
  
  fn extended(&self, other: &Ellipse) -> Ellipse {
    Ellipse::from_cov_matrix(
      self.sigx2 + other.sigx2,
      self.sigy2 + other.sigy2,
      self.rho_sigx_sigy + other.rho_sigx_sigy,
    )
  }
  
  pub fn mahalanobis_distance(&self, x: f64, y: f64) -> f64 {
    let x2 = x.pow2();
    let y2 = y.pow2();
    self.one_over_det * (x2 * self.sigy2 - (self.rho_sigx_sigy * x * y).twice() + y2 * self.sigx2)
  }
  
  pub fn contains(&self, x: f64, y: f64) -> bool {
    self.mahalanobis_distance(x, y) <= 1.0
  }
  
  pub fn overlap(&self, x: f64, y: f64, other: &Ellipse) -> bool {
    self.extended(other).contains(x, y)
  }
  
  pub fn path_along_edge(a: f64, b: f64, theta: f64, half_num_points: usize) -> Box<[(f64, f64)]> {
    let step = (2.0 * a) / (half_num_points as f64);
    let (sin_t, cos_t) = theta.sin_cos();
    let mut v = Vec::<(f64, f64)>::with_capacity(half_num_points << 1);
    for i in 0..half_num_points {
      let x = a - (i as f64) * step;
      let y = b * (1.0 - (x / a).pow2()).sqrt();
      v.push((x, y));
    }
    for i in 0..half_num_points {
      let (x, y) = v[i];
      v.push((-x, -y));
    }
    v.into_boxed_slice()
  }
  
}