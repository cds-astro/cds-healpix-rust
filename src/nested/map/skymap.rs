use std::{
  error::Error,
  fs::File,
  io::{BufReader, BufWriter, Read, Seek, Write},
  iter::{Enumerate, Map},
  marker::PhantomData,
  ops::{Deref, RangeInclusive},
  path::Path,
  slice::Iter,
};

use colorous::Gradient;
use mapproj::CanonicalProjection;
use num_traits::ToBytes;

use crate::{n_hash, nested::map::HHash};

#[cfg(not(target_arch = "wasm32"))]
use super::img::show_with_default_app;
use super::{
  fits::{error::FitsError, read::from_fits_skymap, write::write_implicit_skymap_fits},
  img::{to_skymap_png, ColorMapFunctionType, PosConversion},
};

/// Trait marking the type of the values writable in a FITS skymap.
pub trait SkyMapValue: ToBytes + Clone {
  /// FITS size, in bytes, of a value.
  fn fits_naxis1() -> u8;
  /// FITS type of the value
  fn fits_tform() -> &'static str;
}

impl SkyMapValue for u8 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "B"
  }
}
impl SkyMapValue for i16 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "I"
  }
}
impl SkyMapValue for i32 {
  fn fits_naxis1() -> u8 {
    4
  }
  fn fits_tform() -> &'static str {
    "J"
  }
}
impl SkyMapValue for i64 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "K"
  }
}
impl SkyMapValue for u32 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "J"
  }
}
impl SkyMapValue for u64 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "K"
  }
}
impl SkyMapValue for f32 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "E"
  }
}
impl SkyMapValue for f64 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }
  fn fits_tform() -> &'static str {
    "D"
  }
}

pub trait SkyMap<'a> {
  /// Type of the HEALPix hash value (mainly `u32` or `u64`).
  type HashType: HHash;
  /// Type of the value associated to each HEALPix cell.
  type ValueType: 'a + SkyMapValue;
  /// Type of the iterator iterating on the skymap values.
  type ValuesIt: Iterator<Item = &'a Self::ValueType>;
  /// Type of the iterator iterating on the skymap borrowed entries.
  /// WARNING: we are so far stucked with iterator on ranges,
  /// e.g `(0..n_cell).iter().zip(...)`, since it relies on the `Step` trait
  /// which requires `nightly builds`.
  /// In the case of `implicit` skymaps, a solution is to use `enumerate`.
  type EntriesIt: Iterator<Item = (Self::HashType, &'a Self::ValueType)>;
  // Type of the iterator iterating over the skymap owned entries.
  // type DrainEntriesIt: Iterator<Item = (Self::HashType, Self::ValueType)>;

  /// Depth (<=> HEALPix order) of the skymap.
  fn depth(&self) -> u8;

  /// Tells whether the map is implicit or not.
  /// If implicit, method `values` and `entries` will return as many items as the number
  /// of HEALPix cell at the map HEALPix depth.
  fn is_implicit(&self) -> bool;

  /// Returns the number of elements in the skymap.
  /// For implicit skymaps, the number of elements equals the number of HEALPix cells at the
  /// skymap depth/order.
  fn len(&self) -> usize;

  /// Returns the value associated with the HEALPix cell of given hash number.
  fn get(&self, hash: Self::HashType) -> &Self::ValueType;

  /// Returns all values associated with HEALPix cells, ordered by increasing cell hash number.
  fn values(&'a self) -> Self::ValuesIt;

  /// Returns all entries, i.e. HEALPix cell hash / value tuples, ordered by increasing cell hash number.
  fn entries(&'a self) -> Self::EntriesIt;

  // In case we want to build mom from complex type that are costly to clone.
  // fn drain_entries(self) -> Self::DrainEntriesIt;
}

pub struct ImplicitSkyMapArray<H: HHash, V: SkyMapValue> {
  depth: u8,
  values: Box<[V]>,
  _htype: PhantomData<H>,
}
impl<'a, H: HHash, V: SkyMapValue + 'a> ImplicitSkyMapArray<H, V> {
  /// WARNING: we assume that the coherency between the depth and the number of elements in the
  ///array has already been tested.
  pub fn new(depth: u8, values: Box<[V]>) -> Self {
    assert_eq!(
      n_hash(depth) as usize,
      values.deref().len(),
      "Wrong implicit skymap size. Epecgted: {}. Actual: {}.",
      n_hash(depth),
      values.len()
    );
    Self {
      depth,
      values,
      _htype: PhantomData,
    }
  }
}
impl<'a, H: HHash, V: SkyMapValue + 'a> SkyMap<'a> for ImplicitSkyMapArray<H, V> {
  type HashType = H;
  type ValueType = V;
  type ValuesIt = Iter<'a, Self::ValueType>;
  type EntriesIt = Map<Enumerate<Self::ValuesIt>, fn((usize, &V)) -> (H, &V)>;

  fn depth(&self) -> u8 {
    self.depth
  }

  fn len(&self) -> usize {
    self.values.len()
  }

  fn is_implicit(&self) -> bool {
    true
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    &self.values.deref()[hash.as_()]
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.values.deref().iter()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self
      .values
      .deref()
      .iter()
      .enumerate()
      .map(move |(h, v)| (H::from_usize(h), v))
  }
}

pub enum SkyMapEnum {
  ImplicitU64U8(ImplicitSkyMapArray<u64, u8>),
  ImplicitU64I16(ImplicitSkyMapArray<u64, i16>),
  ImplicitU64I32(ImplicitSkyMapArray<u64, i32>),
  ImplicitU64I64(ImplicitSkyMapArray<u64, i64>),
  ImplicitU64F32(ImplicitSkyMapArray<u64, f32>),
  ImplicitU64F64(ImplicitSkyMapArray<u64, f64>),
}

impl SkyMapEnum {
  #[cfg(not(target_arch = "wasm32"))]
  pub fn from_fits_file<P: AsRef<Path>>(path: P) -> Result<Self, FitsError> {
    File::open(path)
      .map_err(FitsError::Io)
      .map(BufReader::new)
      .and_then(SkyMapEnum::from_fits)
  }

  pub fn from_fits<R: Read + Seek>(reader: BufReader<R>) -> Result<Self, FitsError> {
    from_fits_skymap(reader)
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    match &self {
      Self::ImplicitU64U8(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64I16(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64I32(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64I64(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64F32(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64F64(s) => write_implicit_skymap_fits(writer, s.values.deref()),
    }
  }

  #[cfg(not(target_arch = "wasm32"))]
  pub fn to_skymap_png_file<P: CanonicalProjection, W: AsRef<Path>>(
    &self,
    img_size: (u16, u16),
    proj: Option<P>,
    proj_center: Option<(f64, f64)>,
    proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
    pos_convert: Option<PosConversion>,
    color_map: Option<Gradient>,
    color_map_func_type: Option<ColorMapFunctionType>,
    path: W,
    view: bool,
  ) -> Result<(), Box<dyn Error>> {
    File::create(path.as_ref())
      .map_err(|e| e.into())
      .map(BufWriter::new)
      .and_then(|mut writer| {
        self.to_skymap_png(
          img_size,
          proj,
          proj_center,
          proj_bounds,
          pos_convert,
          color_map,
          color_map_func_type,
          &mut writer,
        )
      })
      .and_then(|()| {
        if view {
          show_with_default_app(path.as_ref().to_string_lossy().as_ref()).map_err(|e| e.into())
        } else {
          Ok(())
        }
      })
  }

  pub fn to_skymap_png<P: CanonicalProjection, W: Write>(
    &self,
    img_size: (u16, u16),
    proj: Option<P>,
    proj_center: Option<(f64, f64)>,
    proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
    pos_convert: Option<PosConversion>,
    color_map: Option<Gradient>,
    color_map_func_type: Option<ColorMapFunctionType>,
    writer: W,
  ) -> Result<(), Box<dyn Error>> {
    match &self {
      Self::ImplicitU64U8(s) => to_skymap_png(
        s,
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      Self::ImplicitU64I16(s) => to_skymap_png(
        s,
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      Self::ImplicitU64I32(s) => to_skymap_png(
        s,
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      Self::ImplicitU64I64(s) => to_skymap_png(
        s,
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      Self::ImplicitU64F32(s) => to_skymap_png(
        s,
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      Self::ImplicitU64F64(s) => to_skymap_png(
        s,
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
    }
  }
}
