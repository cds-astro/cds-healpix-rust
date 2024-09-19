use std::{
  error::Error,
  fs::File,
  io::{BufReader, BufWriter, Read, Seek, Write},
  ops::{Deref, RangeInclusive},
  path::Path,
};

use colorous::Gradient;
use mapproj::CanonicalProjection;

// Make a trait SkyMap: Iterable<Item=(u64, V)> { depth(), get(hash) }
// => Make a map from a FITS file without loading data in memory (bot implicit and explicit)

use super::{
  fits::{error::FitsError, read::from_fits_skymap, write::write_implicit_skymap_fits},
  img::{show_with_default_app, to_png, ColorMapFunctionType, PosConversion},
};

// Implicit map
pub struct SkyMap {
  pub depth: u8,
  pub values: SkyMapArray,
}

pub enum SkyMapArray {
  U8(Box<[u8]>),
  I16(Box<[i16]>),
  I32(Box<[i32]>),
  I64(Box<[i64]>),
  F32(Box<[f32]>),
  F64(Box<[f64]>),
}

impl SkyMap {
  pub fn from_fits_file<P: AsRef<Path>>(path: P) -> Result<Self, FitsError> {
    File::open(path)
      .map_err(FitsError::Io)
      .map(BufReader::new)
      .and_then(SkyMap::from_fits)
  }
  pub fn from_fits<R: Read + Seek>(reader: BufReader<R>) -> Result<Self, FitsError> {
    from_fits_skymap(reader)
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    match &self.values {
      SkyMapArray::U8(b) => write_implicit_skymap_fits(writer, b.deref()),
      SkyMapArray::I16(b) => write_implicit_skymap_fits(writer, b.deref()),
      SkyMapArray::I32(b) => write_implicit_skymap_fits(writer, b.deref()),
      SkyMapArray::I64(b) => write_implicit_skymap_fits(writer, b.deref()),
      SkyMapArray::F32(b) => write_implicit_skymap_fits(writer, b.deref()),
      SkyMapArray::F64(b) => write_implicit_skymap_fits(writer, b.deref()),
    }
  }

  pub fn to_png_file<P: CanonicalProjection, W: AsRef<Path>>(
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
        self.to_png(
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

  pub fn to_png<P: CanonicalProjection, W: Write>(
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
    match &self.values {
      SkyMapArray::U8(b) => to_png(
        b.deref(),
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      SkyMapArray::I16(b) => to_png(
        b.deref(),
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      SkyMapArray::I32(b) => to_png(
        b.deref(),
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      SkyMapArray::I64(b) => to_png(
        b.deref(),
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      SkyMapArray::F32(b) => to_png(
        b.deref(),
        img_size,
        proj,
        proj_center,
        proj_bounds,
        pos_convert,
        color_map,
        color_map_func_type,
        writer,
      ),
      SkyMapArray::F64(b) => to_png(
        b.deref(),
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
