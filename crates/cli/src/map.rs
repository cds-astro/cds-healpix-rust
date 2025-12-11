use std::{error::Error, marker::PhantomData, path::PathBuf};

use clap::{Args, Subcommand, ValueEnum};
use colorous;
use log::{error, warn};

use mapproj::{
  cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
  hybrid::hpx::Hpx,
  pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
  zenithal::{air::Air, arc::Arc, feye::Feye, sin::Sin, stg::Stg, tan::Tan, zea::Zea},
};

use super::{
  get_thread_pool,
  input::pos::{PosCsvConsumed, PosItOperation, PosListInput},
  view::Mode,
};
use hpxlib::nested::map::img::Gradient;
use hpxlib::nested::map::{
  img::{ColorMapFunctionType, PosConversion},
  mom::WritableMom,
  skymap::{explicit::ExplicitCountMap, implicit::ImplicitCountMapU32, SkyMap, SkyMapEnum},
  HHash,
};

/// Create and manipulate HEALPix count and density maps.
#[derive(Debug, Subcommand)]
pub enum Map {
  Count(Count),
  Dens(Dens),
  Op(Operation), // Make a MOM from MAP1 and MAP2 (either left or right reach max N sources)
  Convert(Convert),
  View(View),
}
impl Map {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match self {
      Self::Count(e) => e.exec(),
      Self::Dens(e) => e.exec(),
      Self::Op(e) => e.exec(),
      Self::Convert(e) => e.exec(),
      Self::View(e) => e.exec(),
    }
  }
}

/// Build a count map from a list of positions.
#[derive(Debug, Args)]
pub struct Count {
  /// Healpix depth (aka order) of the count map.
  #[clap(value_parser = clap::value_parser!(u8).range(0..=29))]
  depth: u8,
  /// Count map destination FITS file.
  output: PathBuf,
  #[arg(short = 'e', long)]
  /// Use in-memory `explicit` representation instead of `implicit` (adapted for maps covering less that half the sky).
  explicit: bool,
  #[arg(short = 'r', long = "ratio")]
  /// Limit on the ratio of the implicit over the explicit byte sizes for the FITS serialisation.
  /// Above the limit, the explicit representation is chosen.
  /// Below the limit, the implicit representation is chosen.
  /// If unset, the in-memory representation is chosen.
  implicit_over_explicit_ratio: Option<f64>,
  #[command(subcommand)]
  input: PosListInput,
}
impl Count {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    if self.explicit {
      if self.depth < 13 {
        build_count_map_explicit::<u32>(self.depth, self.input).and_then(|map| {
          match self.implicit_over_explicit_ratio {
            Some(threshold) if map.is_implicit_the_best_representation(0, threshold) => {
              map.into_implicit_skymap(0).to_fits_file(self.output)
            }
            _ => map.to_fits_file(self.output),
          }
          .map_err(|e| e.into())
        })
      } else {
        build_count_map_explicit::<u64>(self.depth, self.input).and_then(|map| {
          match self.implicit_over_explicit_ratio {
            Some(threshold) if map.is_implicit_the_best_representation(0, threshold) => {
              map.into_implicit_skymap(0).to_fits_file(self.output)
            }
            _ => map.to_fits_file(self.output),
          }
          .map_err(|e| e.into())
        })
      }
    } else {
      build_count_map_implicit(self.depth, self.input).and_then(|map| {
        match self.implicit_over_explicit_ratio {
          Some(threshold) if !map.is_implicit_the_best_representation(0, threshold) => {
            map.into_explicit_skymap(0).to_fits_file(self.output)
          }
          _ => map.to_fits_file(self.output),
        }
        .map_err(|e| e.into())
      })
    }
  }
}

/// Build a density map from a list of positions.
#[derive(Debug, Args)]
pub struct Dens {
  /// Healpix depth (aka order) of the count map. Max 12 for in-memory implicit maps.
  depth: u8,
  /// Density map destination FITS file.
  output: PathBuf,
  #[arg(short = 'e', long)]
  /// Use in-memory `explicit` representation instead of `implicit` (adapted for maps covering less that half the sky).
  explicit: bool,
  #[arg(short = 'r', long = "ratio")]
  /// Limit on the ratio of the implicit over the explicit byte sizes for the FITS serialisation.
  /// Above the limit, the explicit representation is chosen.
  /// Below the limit, the implicit representation is chosen.
  /// If unset, the in-memory representation is chosen.
  implicit_over_explicit_ratio: Option<f64>,
  #[command(subcommand)]
  input: PosListInput,
}
impl Dens {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    if self.explicit {
      if self.depth < 13 {
        build_count_map_explicit::<u32>(self.depth, self.input).and_then(|map| {
          let map = map.to_dens_map_par();
          match self.implicit_over_explicit_ratio {
            Some(threshold) if map.is_implicit_the_best_representation(0.0, threshold) => {
              map.into_implicit_skymap(0.0).to_fits_file(self.output)
            }
            _ => map.to_fits_file(self.output),
          }
          .map_err(|e| e.into())
        })
      } else {
        build_count_map_explicit::<u64>(self.depth, self.input).and_then(|map| {
          let map = map.to_dens_map_par();
          match self.implicit_over_explicit_ratio {
            Some(threshold) if map.is_implicit_the_best_representation(0.0, threshold) => {
              warn!("No implicit representation available for DensityMap of depth >= 13!");
              // map.into_implicit_skymap().to_fits_file(self.output)
              map.to_fits_file(self.output)
            }
            _ => map.to_fits_file(self.output),
          }
          .map_err(|e| e.into())
        })
      }
    } else {
      build_count_map_implicit(self.depth, self.input).and_then(|map| {
        map
          .to_dens_map_par() // WARNING: this is performed in the global thread pool!!
          .to_fits_file(self.output)
          .map_err(|e| e.into())
      })
    }
  }
}

fn build_count_map_implicit(
  depth: u8,
  input: PosListInput,
) -> Result<ImplicitCountMapU32, Box<dyn Error>> {
  match input {
    PosListInput::List(e) => {
      struct Op {
        depth: u8,
      }
      impl PosItOperation for Op {
        type R = ImplicitCountMapU32;
        fn exec<I>(self, pos_it: I) -> Result<Self::R, Box<dyn Error>>
        where
          I: Iterator<Item = Result<(f64, f64), Box<dyn Error>>>,
        {
          let it = pos_it.filter_map(|r| match r {
            Ok(t) => Some(t),
            Err(e) => {
              error!("Error reading/parsing line: {:?}", e);
              None
            }
          });
          Ok(ImplicitCountMapU32::from_positions(self.depth, it))
        }
      }
      e.exec_op(Op { depth })
    }
    PosListInput::Csv(PosCsvConsumed {
      input,
      lon,
      lat,
      delimiter,
      header,
      parallel,
      chunk_size,
    }) => {
      let thread_pool = get_thread_pool(parallel);
      if input == PathBuf::from(r"-") {
        ImplicitCountMapU32::from_csv_stdin_par(
          lon - 1,
          lat - 1,
          Some(delimiter),
          header,
          depth,
          chunk_size,
          &thread_pool,
        )
      } else {
        ImplicitCountMapU32::from_csv_file_par(
          input,
          lon - 1,
          lat - 1,
          Some(delimiter),
          header,
          depth,
          chunk_size,
          &thread_pool,
        )
      }
      .map_err(|e| e.into())
    }
  }
}

fn build_count_map_explicit<H: HHash>(
  depth: u8,
  input: PosListInput,
) -> Result<ExplicitCountMap<H>, Box<dyn Error>> {
  match input {
    PosListInput::List(e) => {
      struct Op<H: HHash> {
        depth: u8,
        _phantom: PhantomData<H>,
      }
      impl<H: HHash> PosItOperation for Op<H> {
        type R = ExplicitCountMap<H>;
        fn exec<I>(self, pos_it: I) -> Result<Self::R, Box<dyn Error>>
        where
          I: Iterator<Item = Result<(f64, f64), Box<dyn Error>>>,
        {
          let it = pos_it.filter_map(|r| match r {
            Ok(t) => Some(t),
            Err(e) => {
              error!("Error reading/parsing line: {:?}", e);
              None
            }
          });
          Ok(ExplicitCountMap::from_positions(self.depth, it))
        }
      }
      e.exec_op(Op {
        depth,
        _phantom: PhantomData,
      })
    }
    PosListInput::Csv(PosCsvConsumed {
      input,
      lon,
      lat,
      delimiter,
      header,
      parallel,
      chunk_size,
    }) => {
      let thread_pool = get_thread_pool(parallel);
      if input == PathBuf::from(r"-") {
        ExplicitCountMap::from_csv_stdin_par(
          lon - 1,
          lat - 1,
          Some(delimiter),
          header,
          depth,
          chunk_size,
          &thread_pool,
        )
      } else {
        ExplicitCountMap::from_csv_file_par(
          input,
          lon - 1,
          lat - 1,
          Some(delimiter),
          header,
          depth,
          chunk_size,
          &thread_pool,
        )
      }
      .map_err(|e| e.into())
    }
  }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum OpType {
  Add,
}

/// Map operations (add, ...).
#[derive(Debug, Args)]
pub struct Operation {
  /// Operation to be performed
  #[clap(value_enum)]
  op: OpType,
  /// Path of the left input map FITS file.
  #[clap(value_name = "LHS_FILE")]
  left: PathBuf,
  /// Path of the right input map FITS file.
  #[clap(value_name = "RHS_FILE")]
  right: PathBuf,
  /// Path of the output FITS file.
  #[clap(value_name = "OUT_FILE")]
  output: PathBuf,
  /// Set the number of threads [default: use all available threads]
  #[clap(long, value_name = "N")]
  parallel: Option<usize>,
}
impl Operation {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let l = SkyMapEnum::from_fits_file(self.left)?;
    let r = SkyMapEnum::from_fits_file(self.right)?;
    match self.op {
      OpType::Add => {
        let thread_pool = get_thread_pool(self.parallel);
        thread_pool
          .install(|| l.par_add(r, Default::default()))
          .and_then(|res| res.to_fits_file(self.output))
          .map_err(|e| e.into())
      }
    }
  }
}

// #[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
#[derive(Debug, Subcommand)]
pub enum Conversion {
  /// Transforms a count map into a density map.
  Count2dens,
  /// Transforms a count map into a MOM based on a chi2 merge algorithm.
  Count2chi2mom {
    /// Completeness of the chi2 distribution of 3 degrees of freedom.
    #[clap(default_value_t = 16.266)]
    threshold: f64,
    /// Depth threshold
    #[clap(short, long, value_name = "DEPTH_MIN")]
    depth_threshold: Option<u8>,
  },
  Count2mom {
    /// Upper value on a cell count
    #[clap(value_name = "COUNT_MAX")]
    threshold: u32,
  },
  /// Transforms a density map into a MOM based on a chi2 merge algorithm.
  Dens2chi2mom {
    /// Completeness of the chi2 distribution of 3 degrees of freedom.
    #[clap(default_value_t = 16.266)]
    threshold: f64,
    /// Depth threshold (i.e. smallest depth)
    #[clap(short, long, value_name = "DEPTH_MIN")]
    depth_threshold: Option<u8>,
  },
}

/// Map conversion (counts -> density, ...).
#[derive(Debug, Args)]
pub struct Convert {
  /// Convert from a map type to another map type
  #[clap(subcommand)]
  conversion: Conversion,
  /// Path of the input map FITS file.
  #[clap(value_name = "IN_FILE")]
  input: PathBuf,
  /// Path of the output FITS file.
  #[clap(value_name = "OUT_FILE")]
  output: PathBuf,
  /// Set the number of threads used for 'count2dens' [default: use all available threads]
  #[clap(long, value_name = "N")]
  parallel: Option<usize>,
}
impl Convert {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let skymap = SkyMapEnum::from_fits_file(self.input)?;
    match self.conversion {
      Conversion::Count2dens => {
        let count_map = skymap.to_count_map()?;
        let thread_pool = get_thread_pool(self.parallel);
        let dens_map = thread_pool.install(|| count_map.to_dens_map_par());
        dens_map.to_fits_file(self.output)
      }
      Conversion::Count2chi2mom {
        threshold,
        depth_threshold,
      } => {
        let count_map = skymap.to_count_map()?;
        let mom = count_map.to_chi2_mom(threshold, depth_threshold);
        mom.to_fits_file(self.output, "count")
      }
      Conversion::Count2mom { threshold } => {
        let count_map = skymap.to_count_map()?;
        count_map
          .to_upper_count_threshold_mom(threshold)
          .to_fits_file(self.output, "count")
      }
      Conversion::Dens2chi2mom {
        threshold,
        depth_threshold,
      } => {
        let dens_map = skymap.to_dens_map()?;
        let mom = dens_map.to_chi2_mom(threshold, depth_threshold);
        mom.to_fits_file(self.output, "density")
      }
    }
    .map_err(|e| e.into())
  }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum ColorMap {
  Turbo,
  Viridis,
  Inferno,
  Magma,
  Plasma,
  Cividis,
  Warm,
  Cool,
  CubeHelix,
  YellowGreenBlue,
  YellowOrangeRed,
}
impl ColorMap {
  pub fn get(&self) -> Gradient {
    match self {
      Self::Turbo => colorous::TURBO,
      Self::Viridis => colorous::VIRIDIS,
      Self::Inferno => colorous::INFERNO,
      Self::Magma => colorous::MAGMA,
      Self::Plasma => colorous::PLASMA,
      Self::Cividis => colorous::CIVIDIS,
      Self::Warm => colorous::WARM,
      Self::Cool => colorous::COOL,
      Self::CubeHelix => colorous::CUBEHELIX,
      Self::YellowGreenBlue => colorous::YELLOW_GREEN_BLUE,
      Self::YellowOrangeRed => colorous::YELLOW_ORANGE_RED,
    }
  }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum ColorMapFunction {
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
impl ColorMapFunction {
  pub fn get(&self) -> ColorMapFunctionType {
    match self {
      Self::Linear => ColorMapFunctionType::Linear,
      Self::LinearLog => ColorMapFunctionType::LinearLog,
      Self::LinearPow2 => ColorMapFunctionType::LinearPow2,
      Self::LinearSqrt => ColorMapFunctionType::LinearSqrt,
      Self::LinearAsinh => ColorMapFunctionType::LinearAsinh,
      Self::Log => ColorMapFunctionType::Log,
      Self::Pow2 => ColorMapFunctionType::Pow2,
      Self::Sqrt => ColorMapFunctionType::Sqrt,
      Self::Asinh => ColorMapFunctionType::Asinh,
    }
  }
}

/// Save a PNG file representing the map and visualize it.
#[derive(Debug, Args)]
pub struct View {
  /// Path of the input map FITS file.
  #[clap(value_name = "FITS_FILE")]
  input: PathBuf,

  /// Path of the output PNG file
  #[clap(value_name = "PNG_FILE")]
  output: PathBuf,
  /// Generate the output PNG without showing it.
  #[clap(short = 's', long = "silent")]
  hide: bool,

  /// Show the MAP in equatorial instead of Galactic
  #[clap(short = 'e', long = "equatorial")]
  equatorial: bool,

  /// Color map
  #[clap(long, value_enum, default_value_t = ColorMap::Turbo)]
  color_map: ColorMap,
  /// Color map function
  #[clap(short, long, value_enum, default_value_t = ColorMapFunction::LinearLog)]
  color_map_fn: ColorMapFunction,

  /// Image generation parameters
  #[clap(subcommand)]
  mode: Mode,
}
impl View {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let skymap = SkyMapEnum::from_fits_file(self.input)?;
    let pos_conv = if self.equatorial {
      None
    } else {
      Some(PosConversion::EqMap2GalImg)
    };
    let color_map = Some(self.color_map.get());
    let color_map_fn_type = Some(self.color_map_fn.get());
    match self.mode {
      Mode::AllSky { y_size } => skymap.to_skymap_png_file::<_, _>(
        (y_size << 1, y_size),
        Some(Mol::new()),
        None,
        None,
        pos_conv,
        color_map,
        color_map_fn_type,
        self.output,
        !self.hide,
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
          "car" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Car::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "cea" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Cea::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "cyp" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Cyp::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "mer" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Mer::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          // Hybrid
          "hpx" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Hpx::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          // Psudocylindrical
          "ait" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Ait::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "mol" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Mol::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "par" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Par::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "sfl" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Sfl::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          // Zenithal
          "air" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Air::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "arc" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Arc::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "feye" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Feye::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "sin" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Sin::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "stg" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Stg::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "tan" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Tan::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          "zea" => skymap
            .to_skymap_png_file::<_, _>(
              img_size,
              Some(Zea::new()),
              center,
              proj_bounds,
              pos_conv,
              color_map,
              color_map_fn_type,
              self.output,
              !self.hide,
            )
            .map_err(|e| e.into()),
          _ => Err(format!("Unknown projection '{}'", &proj).into()),
        }
      }
    }
  }
}
