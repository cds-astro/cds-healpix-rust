//! For the Skymap FITS file convention, see [this document](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf).

use std::{
  error::Error,
  fs::File,
  io::{BufReader, BufWriter, Read, Seek, Write},
  ops::{Add, AddAssign, Deref, RangeInclusive},
  path::Path,
};

use colorous::Gradient;
use log::warn;
use num_traits::{ToBytes, Zero};

use mapproj::CanonicalProjection;

#[cfg(not(target_arch = "wasm32"))]
use super::img::show_with_default_app;
use super::{
  fits::{error::FitsError, read::from_fits_skymap, write::write_implicit_skymap_fits},
  img::{to_skymap_png, ColorMapFunctionType, PosConversion},
};
use crate::nested::map::fits::write::write_explicit_skymap_fits;
use crate::nested::map::skymap::explicit::ExplicitSkyMapBTree;
use crate::nested::map::{
  skymap::implicit::{ImplicitCountMap, ImplicitCountMapU32, ImplicitDensityMap, ImplicitSkyMapArray},
  HHash,
};
use crate::nested::n_hash;

pub mod explicit;
pub mod implicit;

/// Trait marking the type of the values writable in a FITS skymap.
pub trait SkyMapValue: ToBytes + Add + AddAssign + Clone + Zero {
  /// FITS size, in bytes, of a value.
  const FITS_NAXIS1: u8 = size_of::<Self>() as u8;
  /// FITS TFORM type of the value
  const FITS_TFORM: &'static str;
  /// Non standard value used in FITS to support various datatypes (including unsigned integer, ...).
  const FITS_DATATYPE: &'static str;
}

impl SkyMapValue for u8 {
  const FITS_TFORM: &'static str = "B";
  const FITS_DATATYPE: &'static str = "u8";
}
impl SkyMapValue for i16 {
  const FITS_TFORM: &'static str = "B";
  const FITS_DATATYPE: &'static str = "i16";
}
impl SkyMapValue for i32 {
  const FITS_TFORM: &'static str = "J";
  const FITS_DATATYPE: &'static str = "i32";
}
impl SkyMapValue for i64 {
  const FITS_TFORM: &'static str = "K";
  const FITS_DATATYPE: &'static str = "i64";
}
impl SkyMapValue for u32 {
  const FITS_TFORM: &'static str = "J";
  const FITS_DATATYPE: &'static str = "u32";
}
impl SkyMapValue for u64 {
  const FITS_TFORM: &'static str = "K";
  const FITS_DATATYPE: &'static str = "u64";
}
impl SkyMapValue for f32 {
  const FITS_TFORM: &'static str = "E";
  const FITS_DATATYPE: &'static str = "f32";
}
impl SkyMapValue for f64 {
  const FITS_TFORM: &'static str = "D";
  const FITS_DATATYPE: &'static str = "f64";
}

// Make a struct:
// SkyMapMap<S, F>
// where
//   F: FnMut(S::ValueType) -> B,
// {
//   skymap: S,
//   f: F
// }
//
// And implement SkyMap on it.

#[allow(clippy::len_without_is_empty)]
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
  /// Type of iterator iterating on owned entries.
  type OwnedEntriesIt: Iterator<Item = (Self::HashType, Self::ValueType)>;

  /// Depth (<=> HEALPix order) of the skymap.
  fn depth(&self) -> u8;

  /// Tells whether the map is implicit or not.
  /// If implicit, method `values` and `entries` will return as many items as the number
  /// of HEALPix cell at the map HEALPix depth.
  fn is_implicit(&self) -> bool;

  fn implicit_byte_size(&self) -> usize {
    size_of::<Self::ValueType>()
      * if self.is_implicit() {
        self.len()
      } else {
        n_hash(self.depth()) as usize
      }
  }

  /// We assume that cell having a value = 0 are not written.
  fn explicit_byte_size(&'a self) -> usize {
    (size_of::<Self::HashType>() + size_of::<Self::ValueType>())
      * if self.is_implicit() {
        let mut n = 0;
        for v in self.values() {
          if !v.is_zero() {
            n += 1;
          }
        }
        n
      } else {
        self.len()
      }
  }

  /// Returns the "best" representation ot be used to serialize (in FITS) this index.
  /// The "best" is determined from the size, in byte, of both representation.
  /// # Params
  /// * `impl_over_expl_limit_ratio` limit on the ratio of the implicit byte size
  /// over the explicit byte size.
  ///     + above the limit, the explicit representation is chosen.
  ///     + below the limit, the implicit representation is chosen.
  /// (This allows to put a size limit on the explicit representation: looking to a value in 'explicit'
  /// requires a binary search, which is longer than the direct access allowed by `implicit`).
  fn is_implicit_the_best_representation(&'a self, impl_over_expl_limit_ratio: f64) -> bool {
    (self.implicit_byte_size() as f64)
      < self.explicit_byte_size() as f64 * impl_over_expl_limit_ratio
  }

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

  /// In case we want to build a map from complex type that are costly to clone.
  fn owned_entries(self) -> Self::OwnedEntriesIt;
}

#[derive(Debug)]
pub enum SkyMapEnum {
  // Implicit
  ImplicitU64U8(ImplicitSkyMapArray<u64, u8>),
  ImplicitU64I16(ImplicitSkyMapArray<u64, i16>),
  ImplicitU64I32(ImplicitSkyMapArray<u64, i32>),
  ImplicitU64I64(ImplicitSkyMapArray<u64, i64>),
  ImplicitU64F32(ImplicitSkyMapArray<u64, f32>),
  ImplicitU64F64(ImplicitSkyMapArray<u64, f64>),
  // Explicit U32
  ExplicitU32U8(ExplicitSkyMapBTree<u32, u8>),
  ExplicitU32I16(ExplicitSkyMapBTree<u32, i16>),
  ExplicitU32I32(ExplicitSkyMapBTree<u32, i32>),
  ExplicitU32I64(ExplicitSkyMapBTree<u32, i64>),
  ExplicitU32F32(ExplicitSkyMapBTree<u32, f32>),
  ExplicitU32F64(ExplicitSkyMapBTree<u32, f64>),
  // Explicit U64
  ExplicitU64U8(ExplicitSkyMapBTree<u64, u8>),
  ExplicitU64I16(ExplicitSkyMapBTree<u64, i16>),
  ExplicitU64I32(ExplicitSkyMapBTree<u64, i32>),
  ExplicitU64I64(ExplicitSkyMapBTree<u64, i64>),
  ExplicitU64F32(ExplicitSkyMapBTree<u64, f32>),
  ExplicitU64F64(ExplicitSkyMapBTree<u64, f64>),
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

  /// Returns the "best" representation ot be used to serialize (in FITS) this index.
  /// The "best" is determined from the size, in byte, of both representation.
  /// # Params
  /// * `impl_over_expl_limit_ratio` limit on the ratio of the implicit byte size
  /// over the explicit byte size.
  ///     + above the limit, the explicit representation is chosen.
  ///     + below the limit, the implicit representation is chosen.
  /// (This allows to put a size limit on the explicit representation: looking to a value in 'explicit'
  /// requires a binary search, which is longer than the direct access allowed by `implicit`).
  pub fn is_implicit_the_best_representation(&self, impl_over_expl_limit_ratio: f64) -> bool {
    match &self {
      Self::ImplicitU64U8(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ImplicitU64I16(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ImplicitU64I32(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ImplicitU64I64(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ImplicitU64F32(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ImplicitU64F64(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU32U8(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU32I16(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU32I32(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU32I64(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU32F32(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU32F64(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU64U8(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU64I16(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU64I32(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU64I64(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU64F32(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
      Self::ExplicitU64F64(e) => e.is_implicit_the_best_representation(impl_over_expl_limit_ratio),
    }
  }

  pub fn to_implicit(self) -> Self {
    match self {
      Self::ImplicitU64U8(_) => self,
      Self::ImplicitU64I16(_) => self,
      Self::ImplicitU64I32(_) => self,
      Self::ImplicitU64I64(_) => self,
      Self::ImplicitU64F32(_) => self,
      Self::ImplicitU64F64(_) => self,
      Self::ExplicitU32U8(e) => Self::ImplicitU64U8(e.to_implicit_map().into()),
      Self::ExplicitU32I16(e) => Self::ImplicitU64I16(e.to_implicit_map().into()),
      Self::ExplicitU32I32(e) => Self::ImplicitU64I32(e.to_implicit_map().into()),
      Self::ExplicitU32I64(e) => Self::ImplicitU64I64(e.to_implicit_map().into()),
      Self::ExplicitU32F32(e) => Self::ImplicitU64F32(e.to_implicit_map().into()),
      Self::ExplicitU32F64(e) => Self::ImplicitU64F64(e.to_implicit_map().into()),
      Self::ExplicitU64U8(e) => Self::ImplicitU64U8(e.to_implicit_map()),
      Self::ExplicitU64I16(e) => Self::ImplicitU64I16(e.to_implicit_map()),
      Self::ExplicitU64I32(e) => Self::ImplicitU64I32(e.to_implicit_map()),
      Self::ExplicitU64I64(e) => Self::ImplicitU64I64(e.to_implicit_map()),
      Self::ExplicitU64F32(e) => Self::ImplicitU64F32(e.to_implicit_map()),
      Self::ExplicitU64F64(e) => Self::ImplicitU64F64(e.to_implicit_map()),
    }
  }

  pub fn to_explicit(self) -> Self {
    match self {
      Self::ImplicitU64U8(e) => Self::ExplicitU64U8(e.to_explicit_map()),
      Self::ImplicitU64I16(e) => Self::ExplicitU64I16(e.to_explicit_map()),
      Self::ImplicitU64I32(e) => Self::ExplicitU64I32(e.to_explicit_map()),
      Self::ImplicitU64I64(e) => Self::ExplicitU64I64(e.to_explicit_map()),
      Self::ImplicitU64F32(e) => Self::ExplicitU64F32(e.to_explicit_map()),
      Self::ImplicitU64F64(e) => Self::ExplicitU64F64(e.to_explicit_map()),
      Self::ExplicitU32U8(_) => self,
      Self::ExplicitU32I16(_) => self,
      Self::ExplicitU32I32(_) => self,
      Self::ExplicitU32I64(_) => self,
      Self::ExplicitU32F32(_) => self,
      Self::ExplicitU32F64(_) => self,
      Self::ExplicitU64U8(_) => self,
      Self::ExplicitU64I16(_) => self,
      Self::ExplicitU64I32(_) => self,
      Self::ExplicitU64I64(_) => self,
      Self::ExplicitU64F32(_) => self,
      Self::ExplicitU64F64(_) => self,
    }
  }

  pub fn to_count_map(self) -> Result<ImplicitCountMap, String> {
    match self {
      Self::ImplicitU64I32(skymap) => Ok(
        ImplicitSkyMapArray::new(skymap.depth, unsafe {
          std::mem::transmute::<Box<[i32]>, Box<[u32]>>(skymap.values)
        })
        .into(),
      ),
      Self::ImplicitU64I64(skymap) => {
        warn!("Transforms i64 skymap into u32 skymap without checking possible overflow!");
        let values = skymap.values.iter().map(|&v| v as u32).collect();
        Ok(ImplicitSkyMapArray::new(skymap.depth, values).into())
      }
      Self::ImplicitU64F32(skymap) => {
        warn!("Transforms f32 skymap into u32 skymap assuming integer values!");
        let values = skymap.values.iter().map(|&v| v as u32).collect();
        Ok(ImplicitSkyMapArray::new(skymap.depth, values).into())
      }
      Self::ImplicitU64F64(skymap) => {
        warn!("Transforms f64 skymap into u32 skymap assuming integer values!");
        let values = skymap.values.iter().map(|&v| v as u32).collect();
        Ok(ImplicitSkyMapArray::new(skymap.depth, values).into())
      }
      _ => Err(String::from("Unable to convert to count map.")),
    }
  }

  pub fn to_count_map_u32(self) -> Result<ImplicitCountMapU32, String> {
    match self {
      Self::ImplicitU64I32(skymap) => Ok(
        ImplicitSkyMapArray::new(skymap.depth, unsafe {
          std::mem::transmute::<Box<[i32]>, Box<[u32]>>(skymap.values)
        })
        .into(),
      ),
      _ => Err(String::from("Unable to convert to count map.")),
    }
  }

  pub fn to_dens_map(self) -> Result<ImplicitDensityMap, String> {
    match self {
      Self::ImplicitU64F64(skymap) => {
        Ok(ImplicitSkyMapArray::new(skymap.depth, skymap.values).into())
      }
      _ => Err(String::from("Unable to convert to count map.")),
    }
  }

  pub fn par_add(self, rhs: Self) -> Result<Self, FitsError> {
    match (self, rhs) {
      (Self::ImplicitU64I32(l), Self::ImplicitU64I32(r)) => Ok(Self::ImplicitU64I32(l.par_add(r))),
      (Self::ImplicitU64F32(l), Self::ImplicitU64F32(r)) => Ok(Self::ImplicitU64F32(l.par_add(r))),
      (Self::ImplicitU64F64(l), Self::ImplicitU64F64(r)) => Ok(Self::ImplicitU64F64(l.par_add(r))),
      _ => Err(FitsError::Custom {
        msg: format!("Skymap type not supported or not the same in 'add' operation"),
      }),
    }
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    match &self {
      // Implicit
      Self::ImplicitU64U8(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64I16(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64I32(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64I64(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64F32(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      Self::ImplicitU64F64(s) => write_implicit_skymap_fits(writer, s.values.deref()),
      // Explicit U32
      Self::ExplicitU32U8(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU32I16(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU32I32(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU32I64(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU32F32(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU32F64(s) => write_explicit_skymap_fits(writer, s),
      // Explicit U64
      Self::ExplicitU64U8(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU64I16(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU64I32(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU64I64(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU64F32(s) => write_explicit_skymap_fits(writer, s),
      Self::ExplicitU64F64(s) => write_explicit_skymap_fits(writer, s),
    }
  }

  pub fn to_fits_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits(BufWriter::new(file)))
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
      // Explicit U32
      Self::ExplicitU32U8(s) => to_skymap_png(
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
      Self::ExplicitU32I16(s) => to_skymap_png(
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
      Self::ExplicitU32I32(s) => to_skymap_png(
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
      Self::ExplicitU32I64(s) => to_skymap_png(
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
      Self::ExplicitU32F32(s) => to_skymap_png(
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
      Self::ExplicitU32F64(s) => to_skymap_png(
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
      // Explicit U64
      Self::ExplicitU64U8(s) => to_skymap_png(
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
      Self::ExplicitU64I16(s) => to_skymap_png(
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
      Self::ExplicitU64I32(s) => to_skymap_png(
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
      Self::ExplicitU64I64(s) => to_skymap_png(
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
      Self::ExplicitU64F32(s) => to_skymap_png(
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
      Self::ExplicitU64F64(s) => to_skymap_png(
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

#[cfg(test)]
mod tests {

  fn init_logger() {
    let log_level = log::LevelFilter::max();
    // let log_level = log::LevelFilter::Error;

    let _ = env_logger::builder()
      // Include all events in tests
      .filter_level(log_level)
      // Ensure events are captured by `cargo test`
      .is_test(true)
      // Ignore errors initializing the logger if tests race to configure it
      .try_init();
  }

  /*  Test only on personal computer
  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_xmm_slew_dens() {
    use std::{fs::read_to_string, time::SystemTime};
    use log::debug;
    use mapproj::pseudocyl::mol::Mol;
    use crate::nested::map::{
      img::{to_mom_png_file, to_skymap_png_file, ColorMapFunctionType, PosConversion},
      skymap::{CountMap, CountMapU32, DensityMap},
    };

    let path = "local_resources/xmmsl3_241122_posonly.csv";
    let content = read_to_string(path).unwrap();
    let img_size = (1366, 768);
    let depth = 8;

    let it = content.lines().skip(1).map(|row| {
      let (l, b) = row.split_once(',').unwrap();
      (
        l.parse::<f64>().unwrap().to_radians(),
        b.parse::<f64>().unwrap().to_radians(),
      )
    });
    let dens_map = DensityMap::from_positions(depth, it);
    to_skymap_png_file::<'_, _, Mol, _>(
      &dens_map.0,
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), // LinearLog
      "local_resources/xmmsl3_241122.dens_map.png",
      false,
    )
    .unwrap();
    // println!("{:?}", dens_map);
    let dens_mom = dens_map.to_chi2_mom();
    to_mom_png_file::<'_, _, Mol, _>(
      &dens_mom,
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
      "local_resources/xmmsl3_241122.dens_mom.png",
      false,
    )
    .unwrap();
  }
  */

  /*  Test only on personal computer
  #[test]
  #[cfg(all(target_os = "linux", not(target_arch = "wasm32")))]
  fn test_xmm_slew_count() {
    use std::{fs::read_to_string, time::SystemTime};
    use log::debug;
    use mapproj::pseudocyl::mol::Mol;
    use crate::nested::map::{
      img::{to_mom_png_file, to_skymap_png_file, ColorMapFunctionType, PosConversion},
      skymap::{CountMap, CountMapU32, DensityMap},
    };

    let path = "local_resources/xmmsl3_241122_posonly.csv";
    let content = read_to_string(path).unwrap();
    let img_size = (1366, 768);
    let depth = 7;

    let it = content.lines().skip(1).map(|row| {
      let (l, b) = row.split_once(',').unwrap();
      (
        l.parse::<f64>().unwrap().to_radians(),
        b.parse::<f64>().unwrap().to_radians(),
      )
    });
    let dens_map = CountMap::from_positions(depth, it);
    to_skymap_png_file::<'_, _, Mol, _>(
      &dens_map.0,
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), // LinearLog
      "local_resources/xmmsl3_241122.dens_map_from_counts.png",
      false,
    )
    .unwrap();

    // println!("{:?}", dens_map);
    let dens_mom = dens_map.to_chi2_mom();
    to_mom_png_file::<'_, _, Mol, _>(
      &dens_mom,
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
      "local_resources/xmmsl3_241122.dens_mom_from_counts.png",
      false,
    )
    .unwrap();
  }
  */

  /* Test only on personal computer
  #[test]
  #[cfg(all(target_os = "linux", not(target_arch = "wasm32")))]
  fn test_countmap_par() {
    use std::{fs::read_to_string, time::SystemTime};
    use log::debug;
    use mapproj::pseudocyl::mol::Mol;
    use crate::nested::map::{
      img::{to_mom_png_file, to_skymap_png_file, ColorMapFunctionType, PosConversion},
      skymap::{CountMap, CountMapU32, DensityMap},
    };

    let n_threads = Some(4);
    let path = "./local_resources/input11.csv"; // 3.5 GB file, depth=4 chunk_size=2M => total time = 8.57 s, i.e 400 MB/s
    let depth = 6; // Test also with 10
    let has_header = false;
    let chunk_size = 2_000_000;
    // Init logger
    init_logger();
    // Build thread pool
    let mut pool_builder = rayon::ThreadPoolBuilder::new();
    if let Some(n_threads) = n_threads {
      pool_builder = pool_builder.num_threads(n_threads);
    }
    let thread_pool = pool_builder.build().unwrap();

    let tstart = SystemTime::now();
    let count_map = CountMapU32::from_csv_file_par(
      path,
      1,
      2,
      Some(','),
      has_header,
      depth,
      chunk_size,
      &thread_pool,
    )
    .unwrap();
    debug!(
      "Count map computed in {} ms",
      SystemTime::now()
        .duration_since(tstart)
        .unwrap_or_default()
        .as_millis()
    );
  }*/
}
