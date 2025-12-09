use std::{error::Error, path::PathBuf};

use clap::{Args, Subcommand};
use log::error;

use mapproj::{
  cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
  hybrid::hpx::Hpx,
  pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
  zenithal::{air::Air, arc::Arc, feye::Feye, sin::Sin, stg::Stg, tan::Tan, zea::Zea},
};

use hpxlib::nested::{
  iter::{
    sortedcoo2sortedhash::SortedHashIt, sortedhash2hashcount::SortedHash2HashCountIncludingZeroIt,
  },
  map::{
    img::{ColorMapFunctionType, PosConversion, Val},
    mom::{
      impls::{bslice::FITSMom, zvec::MomVecImpl},
      new_chi2_count_merger_no_depth_threshold, new_chi2_count_merger_with_depth_threshold,
      Mom as MomTrait, ViewableMom, WritableMom,
    },
    skymap::SkyMapValue,
  },
};

use crate::{
  input::{
    hash::{HashItOperation, HashListInput},
    pos::{PosItOperation, PosListInput},
  },
  map::ColorMapFunction,
  view::Mode,
};

/// Create and manipulate HEALPix count and density MOMs.
#[derive(Debug, Subcommand)]
pub enum Mom {
  // Chi2Count From sorted HPX / sorted coo
  //From(MomFrom),    // HATS,  ProbDens (threshold merge), CountMap, DensMap (chi2merge, threshold merge), Sorted iterator of pos (to interpret has counts)
  Count(Count),
  Op(Operation), // add or mult: we can decide wether the split method return 4 time the initial value or divide the value in 4!
  Convert(Convert),
  View(View),
}
impl Mom {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match self {
      // Self::From(e) => e.exec(),
      Self::Count(e) => e.exec(),
      Self::Op(e) => e.exec(),
      Self::Convert(e) => e.exec(),
      Self::View(e) => e.exec(),
    }
  }
}

/// Build a chi2 count MOM from a list of pre-ordered HEALPix indices or positions.
#[derive(Debug, Args)]
pub struct Count {
  /// Healpix maximum depth (aka order) of the count MOM.
  #[clap(value_parser = clap::value_parser!(u8).range(0..=29))]
  depth: u8,
  /// Completeness of the chi2 distribution of 3 degrees of freedom.
  #[clap(short = 't', long, default_value_t = 16.266)]
  threshold: f64,
  /// Depth threshold
  #[clap(short = 'd', long, value_name = "DEPTH_MIN")]
  depth_threshold: Option<u8>,
  /// Count MOM destination FITS file.
  output: PathBuf,
  #[clap(subcommand)]
  input: PosOrHash,
}
impl Count {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    self.input.exec(
      self.depth,
      self.threshold,
      self.depth_threshold,
      self.output,
    )
  }
}

#[derive(Debug, Subcommand)]
enum PosOrHash {
  /// Create a chi2 MOM from a list of HEALPix ordered positions
  Pos {
    #[clap(subcommand)]
    input: PosListInput,
  },
  /// Create a chi2 MOM from a list of ordered HEALPix hash (aka ipix) at the maximum MOM depth
  Hash {
    #[clap(subcommand)]
    input: HashListInput,
  },
}
impl PosOrHash {
  pub fn exec(
    self,
    depth: u8,
    threshold: f64,
    depth_threshold: Option<u8>,
    output: PathBuf,
  ) -> Result<(), Box<dyn Error>> {
    match self {
      Self::Pos { input } => build_count_mom_from_pos(depth, threshold, depth_threshold, input),
      Self::Hash { input } => build_count_mom_from_hash(depth, threshold, depth_threshold, input),
    }
    .and_then(|mom| mom.to_fits_file(output, "DENSITY").map_err(|e| e.into()))
  }
}
fn build_count_mom_from_pos(
  depth: u8,
  threshold: f64,
  depth_threshold: Option<u8>,
  input: PosListInput,
) -> Result<MomVecImpl<u64, u32>, Box<dyn Error>> {
  struct Op {
    depth: u8,
    threshold: f64,
    depth_threshold: Option<u8>,
  }
  impl PosItOperation for Op {
    type R = MomVecImpl<u64, u32>;
    fn exec<I>(self, pos_it: I) -> Result<Self::R, Box<dyn Error>>
    where
      I: Iterator<Item = Result<(f64, f64), Box<dyn Error>>>,
    {
      let it = pos_it.filter_map(|r| match r {
        Ok(lonlat) => Some(lonlat),
        Err(e) => {
          error!("Error reading/parsing line: {:?}", e);
          None
        }
      });
      let hash_count_it = SortedHash2HashCountIncludingZeroIt::new_from_fallible(
        self.depth,
        SortedHashIt::new(self.depth, it),
      );
      match self.depth_threshold {
        None => MomVecImpl::<u64, u32>::from_hpx_sorted_entries_fallible(
          self.depth,
          hash_count_it,
          new_chi2_count_merger_no_depth_threshold(self.threshold),
        ),
        Some(depth_threshold) => MomVecImpl::<u64, u32>::from_hpx_sorted_entries_fallible(
          self.depth,
          hash_count_it,
          new_chi2_count_merger_with_depth_threshold(self.threshold, depth_threshold),
        ),
      }
      .map_err(|e| e.into())
    }
  }
  input.exec(Op {
    depth,
    threshold,
    depth_threshold,
  })
}
fn build_count_mom_from_hash(
  depth: u8,
  threshold: f64,
  depth_threshold: Option<u8>,
  input: HashListInput,
) -> Result<MomVecImpl<u64, u32>, Box<dyn Error>> {
  struct Op {
    depth: u8,
    threshold: f64,
    depth_threshold: Option<u8>,
  }
  impl HashItOperation for Op {
    type R = MomVecImpl<u64, u32>;
    fn exec<I>(self, hash_it: I) -> Result<Self::R, Box<dyn Error>>
    where
      I: Iterator<Item = Result<u64, Box<dyn Error>>>,
    {
      let it = hash_it.filter_map(|r| match r {
        Ok(hash) => Some(hash),
        Err(e) => {
          error!("Error reading/parsing line: {:?}", e);
          None
        }
      });
      let hash_count_it = SortedHash2HashCountIncludingZeroIt::new_from_infallible(self.depth, it);
      match self.depth_threshold {
        None => MomVecImpl::<u64, u32>::from_hpx_sorted_entries_fallible(
          self.depth,
          hash_count_it,
          new_chi2_count_merger_no_depth_threshold(self.threshold),
        ),
        Some(depth_threshold) => MomVecImpl::<u64, u32>::from_hpx_sorted_entries_fallible(
          self.depth,
          hash_count_it,
          new_chi2_count_merger_with_depth_threshold(self.threshold, depth_threshold),
        ),
      }
      .map_err(|e| e.into())
    }
  }
  input.exec(Op {
    depth,
    threshold,
    depth_threshold,
  })
}

#[derive(Debug, Subcommand)]
pub enum Conversion {
  /// Transforms a count MOM into a density MOM.
  Count2dens,
  /// Transforms a FITS MOM into a BINTABLE FITS MOM file.
  Bintable,
  // gw2fits
  /*/// Transforms a density map into a MOM based on a chi2 merge algorithm.
  Dens2chi2mom {
    /// Completeness of the chi2 distribution of 3 degrees of freedom.
    #[clap(default_value_t = 16.266)]
    threshold: f64,
  },*/
}

/// Map conversion (counts -> density, ...).
#[derive(Debug, Args)]
pub struct Convert {
  /// Convert from a MOM type to another MOM type
  #[clap(subcommand)]
  conversion: Conversion,
  /// Path of the input map FITS file.
  #[clap(value_name = "IN_FILE")]
  input: PathBuf,
  /// Path of the output FITS file.
  #[clap(value_name = "OUT_FILE")]
  output: PathBuf,
}
impl Convert {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let mom = FITSMom::from_fits_file(self.input)?;
    match self.conversion {
      Conversion::Count2dens => match mom {
        FITSMom::U32U32(e) => MomVecImpl::<u32, f64>::from_counts_to_densities(&e.get_mom())
          .to_fits_file(self.output, format!("{}_DENSITY", e.get_value_name()))
          .map_err(|e| e.into()),
        FITSMom::U64U32(e) => MomVecImpl::<u64, f64>::from_counts_to_densities(&e.get_mom())
          .to_fits_file(self.output, format!("{}_DENSITY", e.get_value_name()))
          .map_err(|e| e.into()),
        _ => Err(String::from("Input MOM probably not a count MOM (value type not u32).").into()),
      },
      Conversion::Bintable => mom.to_fits_bintable_file(self.output).map_err(|e| e.into()),
    }
  }
}

#[derive(Debug, Subcommand)]
pub enum OpType {
  /// Multiply left values by right values by a constant and post-chi2 merge.
  /// Values must be f32 or f64 and are assumed to be densities (i.e. remain the same splitting a cell).
  MultAndChi2Merge {
    /// Constant to be used in the multiplication.
    #[clap(short, long, value_name = "CTE", default_value_t = 1.0)]
    cte: f64,
    /// Name assigned to the computed value in the output file.
    #[clap(short = 'n', long, value_name = "NAME", default_value = "MULT")]
    value_name: String,
    /// Completeness of the chi2 distribution of 3 degrees of freedom below which cells are post-merged.
    #[clap(short, long, default_value_t = 16.266)]
    threshold: f64,
    /// Depth threshold
    #[clap(short = 't', long, value_name = "DEPTH_MIN")]
    depth_threshold: Option<u8>,
  },
}

/// Map operations (add, ...).
#[derive(Debug, Args)]
pub struct Operation {
  /// Operation to be performed
  #[clap(subcommand)]
  op: OpType,
  /// Path of the left input MOM FITS file.
  #[clap(value_name = "LHS_FILE")]
  left: PathBuf,
  /// Path of the right input MOM FITS file.
  #[clap(value_name = "RHS_FILE")]
  right: PathBuf,
  /// Path of the output FITS file.
  #[clap(value_name = "OUT_FILE")]
  output: PathBuf,
}
impl Operation {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let l = FITSMom::from_fits_file(self.left)?;
    let r = FITSMom::from_fits_file(self.right)?;
    match self.op {
      OpType::MultAndChi2Merge {
        cte,
        value_name,
        threshold,
        depth_threshold,
      } => match (l, r) {
        (FITSMom::U32U32(_l), FITSMom::U32U32(_r)) => {
          Err(String::from("Method not available for integers").into())
        }
        (FITSMom::U32F32(l), FITSMom::U32F32(r)) => {
          MomVecImpl::lhs_time_rhs_time_cte_and_chi2merge_assuming_densities(
            l.get_mom(),
            r.get_mom(),
            cte as f32,
            threshold,
            depth_threshold,
          )
          .to_fits_file(self.output, value_name)
          .map_err(|e| e.into())
        }
        (FITSMom::U32F64(l), FITSMom::U32F64(r)) => {
          MomVecImpl::lhs_time_rhs_time_cte_and_chi2merge_assuming_densities(
            l.get_mom(),
            r.get_mom(),
            cte,
            threshold,
            depth_threshold,
          )
          .to_fits_file(self.output, value_name)
          .map_err(|e| e.into())
        }
        (FITSMom::U64U32(_l), FITSMom::U64U32(_r)) => {
          Err(String::from("Method not available for integers").into())
        }
        (FITSMom::U64F32(l), FITSMom::U64F32(r)) => {
          MomVecImpl::lhs_time_rhs_time_cte_and_chi2merge_assuming_densities(
            l.get_mom(),
            r.get_mom(),
            cte as f32,
            threshold,
            depth_threshold,
          )
          .to_fits_file(self.output, value_name)
          .map_err(|e| e.into())
        }
        (FITSMom::U64F64(l), FITSMom::U64F64(r)) => {
          MomVecImpl::lhs_time_rhs_time_cte_and_chi2merge_assuming_densities(
            l.get_mom(),
            r.get_mom(),
            cte,
            threshold,
            depth_threshold,
          )
          .to_fits_file(self.output, value_name)
          .map_err(|e| e.into())
        }
        _ => Err(String::from("MOM must be of same type.").into()),
      },
    }
  }
}

/// Save a PNG file representing the MOM and visualize it.
#[derive(Debug, Args)]
pub struct View {
  /// Path of the input MOM FITS file.
  #[clap(value_name = "FITS_FILE")]
  input: PathBuf,

  /// Path of the output PNG file
  #[clap(value_name = "PNG_FILE")]
  output: PathBuf,
  /// Generate the output PNG without showing it.
  #[clap(short = 's', long = "silent")]
  hide: bool,

  /// Color map function
  #[clap(short, long, value_enum, default_value_t = ColorMapFunction::LinearLog)]
  color_map_fn: ColorMapFunction,

  /// Image generation parameters
  #[clap(subcommand)]
  mode: Mode,
}
impl View {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    fn to_png_file<'a, V: SkyMapValue + Val + 'a, M: MomTrait<'a, ValueType = V>>(
      mom: &'a M,
      mode: Mode,
      output: PathBuf,
      hide: bool,
      color_map_fn: ColorMapFunction,
    ) -> Result<(), Box<dyn Error>> {
      let color_map_fn_type = color_map_fn.get();
      match mode {
        Mode::AllSky { y_size } => mom.to_mom_png_file::<_, _>(
          (y_size << 1, y_size),
          Some(Mol::new()),
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(color_map_fn_type),
          output,
          !hide,
        ),
        Mode::Custom {
          proj,
          center_lon,
          center_lat,
          img_size_x,
          img_size_y,
          proj_bounds_x,
          proj_bounds_y,
        } => {
          let img_size = (img_size_x, img_size_y);
          let center = Some((center_lon.to_radians(), center_lat.to_radians()));
          let proj_bounds =
            if let (Some(proj_bounds_x), Some(proj_bounds_y)) = (proj_bounds_x, proj_bounds_y) {
              Some((proj_bounds_x.0, proj_bounds_y.0))
            } else {
              None
            };
          match proj.as_str() {
            // Cylindrical
            "car" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Car::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "cea" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Cea::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(ColorMapFunctionType::LinearLog),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "cyp" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Cyp::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "mer" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Mer::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            // Hybrid
            "hpx" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Hpx::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            // Psudocylindrical
            "ait" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Ait::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "mol" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Mol::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "par" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Par::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "sfl" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Sfl::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            // Zenithal
            "air" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Air::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "arc" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Arc::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "feye" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Feye::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "sin" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Sin::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "stg" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Stg::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "tan" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Tan::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            "zea" => mom
              .to_mom_png_file::<_, _>(
                img_size,
                Some(Zea::new()),
                center,
                proj_bounds,
                Some(PosConversion::EqMap2GalImg),
                None,
                Some(color_map_fn_type),
                output,
                !hide,
              )
              .map_err(|e| e.into()),
            _ => Err(format!("Unknown projection '{}'", &proj).into()),
          }
        }
      }
    }

    match FITSMom::from_fits_file(self.input)? {
      FITSMom::U32U32(fits_mom) => to_png_file(
        &fits_mom.get_mom(),
        self.mode,
        self.output,
        self.hide,
        self.color_map_fn,
      ),
      FITSMom::U32F32(fits_mom) => to_png_file(
        &fits_mom.get_mom(),
        self.mode,
        self.output,
        self.hide,
        self.color_map_fn,
      ),
      FITSMom::U32F64(fits_mom) => to_png_file(
        &fits_mom.get_mom(),
        self.mode,
        self.output,
        self.hide,
        self.color_map_fn,
      ),
      FITSMom::U64U32(fits_mom) => to_png_file(
        &fits_mom.get_mom(),
        self.mode,
        self.output,
        self.hide,
        self.color_map_fn,
      ),
      FITSMom::U64F32(fits_mom) => to_png_file(
        &fits_mom.get_mom(),
        self.mode,
        self.output,
        self.hide,
        self.color_map_fn,
      ),
      FITSMom::U64F64(fits_mom) => to_png_file(
        &fits_mom.get_mom(),
        self.mode,
        self.output,
        self.hide,
        self.color_map_fn,
      ),
    }
  }
}
