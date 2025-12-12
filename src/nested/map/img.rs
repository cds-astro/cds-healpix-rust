use std::{
  cmp::Ordering,
  collections::HashSet,
  error::Error,
  f64::consts::PI,
  fs::File,
  io::{self, BufWriter, Write},
  ops::RangeInclusive,
  path::Path,
};

pub use colorous::Gradient;
use mapproj::{
  img2celestial::Img2Celestial, img2proj::ReversedEastPngImgXY2ProjXY, pseudocyl::mol::Mol,
  CanonicalProjection, CenteredProjection, ImgXY, LonLat,
};

use crate::nested::{
  self,
  map::{
    astrometry::{gal::Galactic, math::Coo},
    mom::{Mom, ZUniqHashT},
    skymap::SkyMap,
    HHash,
  },
};

/// Wanted image coordinate frame.
pub enum ImgCooFrame {
  Equatorial,
  Galactic,
}

pub enum PosConversion {
  SameMapAndImg,
  EqMap2GalImg,
  GalMap2EqImg,
}
impl PosConversion {
  pub fn convert_img_pos_to_map_pos(&self) -> fn(f64, f64) -> (f64, f64) {
    match self {
      PosConversion::SameMapAndImg => |l, b| (l, b),
      PosConversion::EqMap2GalImg => Self::gal2eq(),
      PosConversion::GalMap2EqImg => Self::eq2gal(),
    }
  }
  pub fn convert_map_pos_to_img_pos(&self) -> fn(f64, f64) -> (f64, f64) {
    match self {
      PosConversion::SameMapAndImg => |l, b| (l, b),
      PosConversion::EqMap2GalImg => Self::eq2gal(),
      PosConversion::GalMap2EqImg => Self::gal2eq(),
    }
  }
  pub fn gal2eq() -> fn(f64, f64) -> (f64, f64) {
    |lon, lat| {
      const GAL: Galactic = Galactic::new_for_icrs_including_fk5_icrs_offsets_from_mignard2002();
      let coo_eq = Coo::new(lon, lat);
      let coo_gal = GAL.coo_gal2eq(&coo_eq);
      (coo_gal.lon, coo_gal.lat)
    }
  }
  pub fn eq2gal() -> fn(f64, f64) -> (f64, f64) {
    |lon, lat| {
      const GAL: Galactic = Galactic::new_for_icrs_including_fk5_icrs_offsets_from_mignard2002();
      let coo_gal = Coo::new(lon, lat);
      let coo_eq = GAL.coo_eq2gal(&coo_gal);
      (coo_eq.lon, coo_eq.lat)
    }
  }
}

pub trait Val: PartialOrd + Copy {
  fn to_f64(&self) -> f64;
  fn to_bits_repr(&self) -> u64;
}

impl Val for u8 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for u16 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for u32 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for u64 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self
  }
}

impl Val for i8 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for i16 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for i32 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for i64 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    *self as u64
  }
}

impl Val for f32 {
  fn to_f64(&self) -> f64 {
    *self as f64
  }
  fn to_bits_repr(&self) -> u64 {
    self.to_bits() as u64
  }
}

impl Val for f64 {
  fn to_f64(&self) -> f64 {
    *self
  }
  fn to_bits_repr(&self) -> u64 {
    self.to_bits()
  }
}

pub enum ColorMapFunctionType {
  Linear,
  LinearLog,
  LinearPow2,
  LinearSqrt,
  LinearAsinh,
  Log,
  Pow2,
  Sqrt,
  Asinh,
}

///
/// * `min`: minimum value corresponding to 0 in the ColorMap
/// * `max`: maximum value corresponding to 1 in the ColorMap
#[derive(Debug, Copy, Clone)]
pub enum ColorMapFunction {
  Linear {
    min: f64,
    max: f64,
  },
  // Use formula: 0.5*log(x+0.01)+1 returning value in [0, 1], x max must be 0.99
  LinearLog {
    min: f64,
    max: f64,
  },
  LinearPow2 {
    min: f64,
    max: f64,
  },
  LinearSqrt {
    min: f64,
    max: f64,
  },
  LinearAsinh {
    min: f64,
    max: f64,
  },
  Log {
    min: f64,
    max: f64,
    fmin: f64,
    fmax: f64,
  },
  Pow2 {
    min: f64,
    max: f64,
    fmin: f64,
    fmax: f64,
  },
  Sqrt {
    min: f64,
    max: f64,
    fmin: f64,
    fmax: f64,
  },
  Asinh {
    min: f64,
    max: f64,
    fmin: f64,
    fmax: f64,
  },
}
impl ColorMapFunction {
  pub fn new_from(cm_type: ColorMapFunctionType, min: f64, max: f64) -> Self {
    match cm_type {
      ColorMapFunctionType::Linear => Self::new_linear(min, max),
      ColorMapFunctionType::LinearLog => Self::new_linear_log(min, max),
      ColorMapFunctionType::LinearPow2 => Self::new_linear_pow2(min, max),
      ColorMapFunctionType::LinearSqrt => Self::new_linear_sqrt(min, max),
      ColorMapFunctionType::LinearAsinh => Self::new_linear_asinh(min, max),
      ColorMapFunctionType::Log => Self::new_log(min, max),
      ColorMapFunctionType::Pow2 => Self::new_pow2(min, max),
      ColorMapFunctionType::Sqrt => Self::new_sqrt(min, max),
      ColorMapFunctionType::Asinh => Self::new_asinh(min, max),
    }
  }
  pub fn new_linear(min: f64, max: f64) -> Self {
    ColorMapFunction::Linear { min, max }
  }
  pub fn new_linear_log(min: f64, max: f64) -> Self {
    ColorMapFunction::LinearLog { min, max }
  }
  pub fn new_linear_pow2(min: f64, max: f64) -> Self {
    ColorMapFunction::LinearPow2 { min, max }
  }
  pub fn new_linear_sqrt(min: f64, max: f64) -> Self {
    ColorMapFunction::LinearSqrt { min, max }
  }
  pub fn new_linear_asinh(min: f64, max: f64) -> Self {
    ColorMapFunction::LinearAsinh { min, max }
  }

  pub fn new_log(min: f64, max: f64) -> Self {
    ColorMapFunction::Log {
      min,
      max,
      fmin: (min + 0.01).ln(),
      fmax: (max + 0.01).ln(),
    }
  }
  pub fn new_pow2(min: f64, max: f64) -> Self {
    ColorMapFunction::Pow2 {
      min,
      max,
      fmin: min * min,
      fmax: max * max,
    }
  }
  pub fn new_sqrt(min: f64, max: f64) -> Self {
    ColorMapFunction::Sqrt {
      min,
      max,
      fmin: min.sqrt(),
      fmax: max.sqrt(),
    }
  }
  pub fn new_asinh(min: f64, max: f64) -> Self {
    ColorMapFunction::Asinh {
      min,
      max,
      fmin: min.asinh(),
      fmax: max.asinh(),
    }
  }

  /// Returns a value in `[0, 1]`.
  pub fn value(&self, value: f64) -> f64 {
    match self {
      ColorMapFunction::Linear { min, max } => Self::linear_val(*min, *max, value, |v| v),
      ColorMapFunction::LinearLog { min, max } =>
      // 1 + log10(0.01) / 2 = 0 since log10(10^-2) = -2 * log10(10) = -2 * ln(10)/ln(10)
      // 1 + log10(1.00) / 2 = 1 since log10(1) = ln(1)/ln(10) and ln(1) = 0
      {
        Self::linear_val(*min, *max, value, |v| 1.0 + 0.5 * (0.99 * v + 0.01).log10())
      }
      ColorMapFunction::LinearPow2 { min, max } => Self::linear_val(*min, *max, value, |v| v * v),
      ColorMapFunction::LinearSqrt { min, max } => {
        Self::linear_val(*min, *max, value, |v| v.sqrt())
      }
      ColorMapFunction::LinearAsinh { min, max } => {
        Self::linear_val(*min, *max, value, |v| v.asinh() / 1.0_f64.asinh())
      }
      ColorMapFunction::Log {
        min,
        max,
        fmin,
        fmax,
      } => Self::val(*min, *max, *fmin, *fmax, value, |v| (v + 0.01).ln()),
      ColorMapFunction::Pow2 {
        min,
        max,
        fmin,
        fmax,
      } => Self::val(*min, *max, *fmin, *fmax, value, |v| v * v),
      ColorMapFunction::Sqrt {
        min,
        max,
        fmin,
        fmax,
      } => Self::val(*min, *max, *fmin, *fmax, value, |v| v.sqrt()),
      ColorMapFunction::Asinh {
        min,
        max,
        fmin,
        fmax,
      } => Self::val(*min, *max, *fmin, *fmax, value, |v| v.asinh()),
    }
  }
  /// Reverse of the value, i.e. returns a value in `[1, 0]`.
  pub fn value_rev(&self, value: f64) -> f64 {
    1.0 - self.value(value)
  }

  // Function F takes a values in [0, 1] and must return a value in [0, 1].
  fn linear_val<F>(min: f64, max: f64, value: f64, f: F) -> f64
  where
    F: Fn(f64) -> f64,
  {
    if value > max {
      1.0
    } else if value < min {
      0.0
    } else {
      let a = (value - min) / (max - min);
      f(a)
    }
  }
  fn val<F>(min: f64, max: f64, fmin: f64, fmax: f64, value: f64, f: F) -> f64
  where
    F: Fn(f64) -> f64,
  {
    if value > max {
      1.0
    } else if value < min {
      0.0
    } else {
      (f(value) - fmin) / (fmax - fmin)
    }
  }
}

/// Returns an RGBA array (each pixel is made of 4 successive u8: RGBA) using the Mollweide projection.
///
/// # Params
/// * `skymap`: the skymap to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
pub fn to_skymap_img_default<'a, S>(
  skymap: &'a S,
  img_size: (u16, u16),
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Option<Gradient>,
  color_map_func_type: Option<ColorMapFunctionType>,
) -> Result<Vec<u8>, Box<dyn Error>>
where
  S: SkyMap<'a>,
  S::ValueType: Val,
{
  let proj = Mol::new();
  let color_map = color_map.unwrap_or(colorous::TURBO);
  let color_map_func_type = color_map_func_type.unwrap_or(ColorMapFunctionType::Linear);
  to_skymap_img(
    skymap,
    img_size,
    proj,
    proj_center,
    proj_bounds,
    pos_convert,
    color_map,
    color_map_func_type,
  )
}

/// Returns an RGBA array (each pixel is made of 4 successive u8: RGBA).
///
/// # Params
/// * `skymap`: the skymap to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj`: a projection, if different from Mollweide;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
/// * `pos_convert`: to handle a different coordinate system between the skymap and the image.
/// * `color_map`:
/// * `color_map_func`: the color map fonction, build it using:
// TODO: take in input the value coding NULL
pub fn to_skymap_img<'a, P, S>(
  skymap: &'a S,
  img_size: (u16, u16),
  proj: P,
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Gradient,
  color_map_func_type: ColorMapFunctionType,
) -> Result<Vec<u8>, Box<dyn Error>>
where
  P: CanonicalProjection,
  S: SkyMap<'a>,
  S::ValueType: Val,
{
  const N_VAL_MAX: usize = 25; // 255 distinct color values
  let depth = skymap.depth();

  // We remove from value 0s and NaNs
  let mut iter = skymap.values().filter(|v| v == v && v.to_f64() != 0.0);
  let (min, max) = if let Some(first_value) = iter.next() {
    let mut value_set = HashSet::with_capacity(N_VAL_MAX);
    value_set.insert(first_value.to_bits_repr());
    let (min, max, value_set) = iter.fold(
      (first_value, first_value, value_set),
      |(min, max, mut value_set), val| {
        if value_set.len() < N_VAL_MAX {
          value_set.insert(val.to_bits_repr());
        }
        // unwrap() are ok since we removed NaN above
        (
          std::cmp::min_by(val, min, |a, b| a.partial_cmp(b).unwrap()),
          std::cmp::max_by(val, max, |a, b| a.partial_cmp(b).unwrap()),
          value_set,
        )
      },
    );
    let (mut min, max) = (min.to_f64(), max.to_f64());
    // We do this so that the min value does not lead to the 0.0 value in the color map,
    // if we hace less than N_VAL_MAX distinct values.
    if min == max {
      debug_assert_eq!(value_set.len(), 1);
      min -= 1.0;
    } else if (1..=N_VAL_MAX).contains(&value_set.len()) {
      if min > 0.0 {
        min = 0.0_f64.max(min - (max - min) / value_set.len() as f64)
      } else {
        min -= (max - min) / value_set.len() as f64;
      }
    }
    (min, max)
  } else {
    // To avoid / by 0
    (0.0, 1.0)
  };

  let color_map_func = ColorMapFunction::new_from(color_map_func_type, min, max);

  let (size_x, size_y) = img_size;
  let mut v: Vec<u8> = Vec::with_capacity((size_x as usize * size_y as usize) << 2);

  let (proj_range_x, proj_range_y) = proj_bounds.unwrap_or((
    proj
      .bounds()
      .x_bounds()
      .as_ref()
      .cloned()
      .unwrap_or_else(|| -PI..=PI),
    proj
      .bounds()
      .y_bounds()
      .as_ref()
      .cloned()
      .unwrap_or_else(|| -PI..=PI),
  ));

  let img2proj =
    ReversedEastPngImgXY2ProjXY::from((size_x, size_y), (&proj_range_x, &proj_range_y));
  let mut img2cel = Img2Celestial::new(img2proj, CenteredProjection::new(proj));
  if let Some((lon, lat)) = proj_center {
    img2cel.set_proj_center_from_lonlat(&LonLat::new(lon, lat));
  }

  let hpx = nested::get(depth);

  let pos_convert = pos_convert.unwrap_or(PosConversion::SameMapAndImg);
  let mappos2imgpos = pos_convert.convert_map_pos_to_img_pos();
  let imgpos2mappos = pos_convert.convert_img_pos_to_map_pos();

  for y in 0..size_y {
    for x in 0..size_x {
      if let Some(lonlat) = img2cel.img2lonlat(&ImgXY::new(x as f64, y as f64)) {
        let (lon, lat) = imgpos2mappos(lonlat.lon(), lonlat.lat());
        let idx = hpx.hash(lon, lat);
        let val = skymap.get(S::HashType::from_u64(idx));
        // num_traits::cast::ToPrimitive contains to_f64
        let color = color_map.eval_continuous(color_map_func.value(val.to_f64()));
        v.push(color.r);
        v.push(color.g);
        v.push(color.b);
        v.push(255);
      } else {
        // Not in the proj area
        v.push(255);
        v.push(255);
        v.push(255);
        v.push(0);
      }
    }
  }
  // But, in case of sparse map with cells smaller than img pixels
  let color0 = color_map.eval_continuous(0.0);
  for (idx, val) in skymap.entries().filter_map(|(idx, val)| {
    let val = val.to_f64();
    if val == val && val != 0.0 {
      Some((idx.to_u64(), val))
    } else {
      None
    }
  }) {
    let (lon_rad, lat_rad) = nested::center(depth, idx);
    let (lon_rad, lat_rad) = mappos2imgpos(lon_rad, lat_rad);
    if let Some(xy) = img2cel.lonlat2img(&LonLat::new(lon_rad, lat_rad)) {
      let ix = xy.x() as u16;
      let iy = xy.y() as u16;
      if ix < img_size.0 && iy < img_size.1 {
        // <<2 <=> *4
        let from = (xy.y() as usize * size_x as usize + ix as usize) << 2;
        // If pixel was associated to value 0.0
        if v[from] == color0.r && v[from + 1] == color0.g && v[from + 2] == color0.b {
          let color = color_map.eval_continuous(color_map_func.value(val));
          v[from] = color.r;
          v[from + 1] = color.g;
          v[from + 2] = color.b;
          v[from + 3] = 255;
        }
      }
    }
  }
  Ok(v)
}

/// Returns an RGBA array (each pixel is made of 4 successive u8: RGBA) using the Mollweide projection.
///
/// # Params
/// * `mom`: the MOM to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
pub fn to_mom_img_default<'a, M>(
  mom: &'a M,
  img_size: (u16, u16),
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Option<Gradient>,
  color_map_func_type: Option<ColorMapFunctionType>,
) -> Result<Vec<u8>, Box<dyn Error>>
where
  M: Mom<'a>,
  M::ValueType: Val,
{
  let proj = Mol::new();
  let color_map = color_map.unwrap_or(colorous::TURBO);
  let color_map_func_type = color_map_func_type.unwrap_or(ColorMapFunctionType::Linear);
  to_mom_img(
    mom,
    img_size,
    proj,
    proj_center,
    proj_bounds,
    pos_convert,
    color_map,
    color_map_func_type,
  )
}

/// Returns an RGBA array (each pixel is made of 4 successive u8: RGBA).
///
/// # Params
/// * `mom`: the MOM to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj`: a projection, if different from Mollweide;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
/// * `pos_convert`: to handle a different coordinate system between the skymap and the image.
/// * `color_map`:
/// * `color_map_func`: the color map fonction, build it using:
// TODO: take in input the value coding NULL
pub fn to_mom_img<'a, P, M>(
  mom: &'a M,
  img_size: (u16, u16),
  proj: P,
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Gradient,
  color_map_func_type: ColorMapFunctionType,
) -> Result<Vec<u8>, Box<dyn Error>>
where
  P: CanonicalProjection,
  M: Mom<'a>,
  M::ValueType: Val,
{
  let depth = mom.depth_max();

  // Unwrap is ok here since we already tested 'n_cell' which must be > 0.
  let mut iter = mom.values_copy();
  let first_value = iter.next().unwrap();
  let (min, max) = iter.fold((first_value, first_value), |(min, max), val| {
    (
      std::cmp::min_by(val, min, |a, b| {
        a.partial_cmp(b).unwrap_or_else(move || {
          if a != a {
            Ordering::Greater
          } else {
            Ordering::Less
          }
        })
      }),
      std::cmp::max_by(val, max, |a, b| {
        a.partial_cmp(b).unwrap_or_else(move || {
          if a != a {
            Ordering::Greater
          } else {
            Ordering::Less
          }
        })
      }),
    )
  });
  let (min, max) = (min.to_f64(), max.to_f64());
  let color_map_func = ColorMapFunction::new_from(color_map_func_type, min, max);

  let (size_x, size_y) = img_size;
  let mut v: Vec<u8> = Vec::with_capacity((size_x as usize * size_y as usize) << 2);

  let (proj_range_x, proj_range_y) = proj_bounds.unwrap_or((
    proj
      .bounds()
      .x_bounds()
      .as_ref()
      .cloned()
      .unwrap_or_else(|| -PI..=PI),
    proj
      .bounds()
      .y_bounds()
      .as_ref()
      .cloned()
      .unwrap_or_else(|| -PI..=PI),
  ));

  let img2proj =
    ReversedEastPngImgXY2ProjXY::from((size_x, size_y), (&proj_range_x, &proj_range_y));
  let mut img2cel = Img2Celestial::new(img2proj, CenteredProjection::new(proj));
  if let Some((lon, lat)) = proj_center {
    img2cel.set_proj_center_from_lonlat(&LonLat::new(lon, lat));
  }

  let hpx = nested::get(depth);

  let pos_convert = pos_convert.unwrap_or(PosConversion::SameMapAndImg);
  let mappos2imgpos = pos_convert.convert_map_pos_to_img_pos();
  let imgpos2mappos = pos_convert.convert_img_pos_to_map_pos();

  for y in 0..size_y {
    for x in 0..size_x {
      if let Some(lonlat) = img2cel.img2lonlat(&ImgXY::new(x as f64, y as f64)) {
        let (lon, lat) = imgpos2mappos(lonlat.lon(), lonlat.lat());
        let idx = hpx.hash(lon, lat);
        let color = if let Some((_, val)) = mom.get_copy_of_cell_containing_unsafe(
          M::ZUniqHType::to_zuniq(depth, M::ZUniqHType::from_u64(idx)),
        ) {
          color_map.eval_continuous(color_map_func.value(val.to_f64()))
        } else {
          color_map.eval_continuous(color_map_func.value(0.0))
        };
        v.push(color.r);
        v.push(color.g);
        v.push(color.b);
        v.push(255);
      } else {
        // Not in the proj area
        v.push(255);
        v.push(255);
        v.push(255);
        v.push(0);
      }
    }
  }
  // But, in case of sparse map with cells smaller than img pixels
  for (depth, idx, val) in mom.entries_copy().filter_map(|(zuniq, val)| {
    let val = val.to_f64();
    if val > 0.0 {
      let (depth, idx) = M::ZUniqHType::from_zuniq(zuniq);
      Some((depth, idx.to_u64(), val))
    } else {
      None
    }
  }) {
    let (lon_rad, lat_rad) = nested::center(depth, idx);
    let (lon_rad, lat_rad) = mappos2imgpos(lon_rad, lat_rad);
    if let Some(xy) = img2cel.lonlat2img(&LonLat::new(lon_rad, lat_rad)) {
      let ix = xy.x() as u16;
      let iy = xy.y() as u16;
      if ix < img_size.0 && iy < img_size.1 {
        let from = (xy.y() as usize * size_x as usize + ix as usize) << 2; // <<2 <=> *4
        if v[from] == 0 {
          let color = color_map.eval_continuous(color_map_func.value(val));
          v[from] = color.r;
          v[from + 1] = color.g;
          v[from + 2] = color.b;
          v[from + 3] = 128;
        }
      }
    }
  }
  Ok(v)
}

/// # Params
/// * `skymap`: the skymap to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj`: a projection, if different from Mollweide;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
/// * `writer`: the writer in which the image is going to be written
pub fn to_skymap_png<'a, P, S, W>(
  skymap: &'a S,
  img_size: (u16, u16),
  proj: Option<P>,
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Option<Gradient>,
  color_map_func_type: Option<ColorMapFunctionType>,
  writer: W,
) -> Result<(), Box<dyn Error>>
where
  P: CanonicalProjection,
  S: SkyMap<'a>,
  S::ValueType: Val,
  W: Write,
{
  let (xsize, ysize) = img_size;
  let data = if let Some(proj) = proj {
    let color_map = color_map.unwrap_or(colorous::TURBO);
    let color_map_func_type = color_map_func_type.unwrap_or(ColorMapFunctionType::Linear);
    to_skymap_img(
      skymap,
      img_size,
      proj,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
    )
  } else {
    to_skymap_img_default(
      skymap,
      img_size,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
    )
  }?;
  let mut encoder = png::Encoder::new(writer, xsize as u32, ysize as u32); // Width is 2 pixels and height is 1.
  encoder.set_color(png::ColorType::Rgba);
  encoder.set_depth(png::BitDepth::Eight);
  let mut writer = encoder.write_header()?;
  writer.write_image_data(&data).expect("Wrong encoding");
  Ok(())
}

/// # Params
/// * `skymap`: the skymap to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj`: a projection, if different from Mollweide;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
/// * `path`: the path of th PNG file to be written.
/// * `view`: set to true to visualize the saved image.
#[cfg(not(target_arch = "wasm32"))]
pub fn to_skymap_png_file<'a, S, P, A>(
  skymap_implicit: &'a S,
  img_size: (u16, u16),
  proj: Option<P>,
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Option<Gradient>,
  color_map_func_type: Option<ColorMapFunctionType>,
  path: A,
  view: bool,
) -> Result<(), Box<dyn Error>>
where
  P: CanonicalProjection,
  S: SkyMap<'a>,
  S::ValueType: Val,
  A: AsRef<Path>,
{
  // Brackets are important to be sure the file is closed before trying to open it.
  {
    let file = File::create(path.as_ref())?;
    let mut writer = BufWriter::new(file);
    to_skymap_png(
      skymap_implicit,
      img_size,
      proj,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
      &mut writer,
    )?;
  }
  if view {
    show_with_default_app(path.as_ref().to_string_lossy().as_ref())?;
  }
  Ok(())
}

/// # Params
/// * `mom`: the MOM to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj`: a projection, if different from Mollweide;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
/// * `writer`: the writer in which the image is going to be written
pub fn to_mom_png<'a, P, M, W>(
  mom: &'a M,
  img_size: (u16, u16),
  proj: Option<P>,
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Option<Gradient>,
  color_map_func_type: Option<ColorMapFunctionType>,
  writer: W,
) -> Result<(), Box<dyn Error>>
where
  P: CanonicalProjection,
  M: Mom<'a>,
  M::ValueType: Val,
  W: Write,
{
  let (xsize, ysize) = img_size;
  let data = if let Some(proj) = proj {
    let color_map = color_map.unwrap_or(colorous::TURBO);
    let color_map_func_type = color_map_func_type.unwrap_or(ColorMapFunctionType::Linear);
    to_mom_img(
      mom,
      img_size,
      proj,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
    )
  } else {
    to_mom_img_default(
      mom,
      img_size,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
    )
  }?;
  let mut encoder = png::Encoder::new(writer, xsize as u32, ysize as u32); // Width is 2 pixels and height is 1.
  encoder.set_color(png::ColorType::Rgba);
  encoder.set_depth(png::BitDepth::Eight);
  let mut writer = encoder.write_header()?;
  writer.write_image_data(&data).expect("Wrong encoding");
  Ok(())
}

/// # Params
/// * `mom`: the MOM to be print;
/// * `size`: the `(X, Y)` number of pixels in the image;
/// * `proj`: a projection, if different from Mollweide;
/// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
///   if different from `(0, 0)`;
/// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
///   which depends on the projection. For unbounded projections, de default value
///   is `(-PI..PI, -PI..PI)`.
/// * `path`: the path of th PNG file to be written.
/// * `view`: set to true to visualize the saved image.
#[cfg(not(target_arch = "wasm32"))]
pub fn to_mom_png_file<'a, M, P, A>(
  mom: &'a M,
  img_size: (u16, u16),
  proj: Option<P>,
  proj_center: Option<(f64, f64)>,
  proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
  pos_convert: Option<PosConversion>,
  color_map: Option<Gradient>,
  color_map_func_type: Option<ColorMapFunctionType>,
  path: A,
  view: bool,
) -> Result<(), Box<dyn Error>>
where
  P: CanonicalProjection,
  M: Mom<'a>,
  M::ValueType: Val,
  A: AsRef<Path>,
{
  // Brackets are important to be sure the file is closed before trying to open it.
  {
    let file = File::create(path.as_ref())?;
    let mut writer = BufWriter::new(file);
    to_mom_png(
      mom,
      img_size,
      proj,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
      &mut writer,
    )?;
  }
  if view {
    show_with_default_app(path.as_ref().to_string_lossy().as_ref())?;
  }
  Ok(())
}

// Adapted from https://github.com/igiagkiozis/plotly/blob/master/plotly/src/plot.rs
#[cfg(target_os = "linux")]
pub fn show_with_default_app(path: &str) -> Result<(), io::Error> {
  use std::process::Command;
  Command::new("xdg-open").args([path]).output()?;
  // .map_err(|e| e.into())?;
  Ok(())
}

#[cfg(target_os = "macos")]
pub fn show_with_default_app(path: &str) -> Result<(), io::Error> {
  use std::process::Command;
  Command::new("open").args(&[path]).output()?;
  Ok(())
}

#[cfg(target_os = "windows")]
pub fn show_with_default_app(path: &str) -> Result<(), io::Error> {
  use std::process::Command;
  Command::new("cmd")
    .arg("/C")
    .arg(format!(r#"start {}"#, path))
    .output()?;
  Ok(())
}
