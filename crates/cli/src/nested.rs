use std::{
  error::Error,
  io::{stdout, StdoutLock, Write},
};

use clap::{Args, Subcommand};

use hpxlib::{
  best_starting_depth,
  compass_point::MainWind,
  has_best_starting_depth, n_hash,
  nested::{
    bilinear_interpolation, center, get, hash, hash_with_dxdy, neighbours, to_uniq, to_uniq_ivoa,
    to_zuniq, vertices,
  },
  DEPTH_MAX,
};

use crate::input::{depthash::DepthHashInput, hash::HashInput, pos::PosInput};

#[derive(Debug, Subcommand)]
pub enum Nested {
  /// Prints the table of Healpix layers properties (depth, nside, ncells, ...)
  Table,
  #[clap(name = "bestdepth")]
  BestStartingDepth(BestStartingDepth),
  Hash(Hash),
  Bilinear(Bilinear),
  Center(Center),
  Vertices(Vertices),
  // PathAlongCellEdge() ??
  // Grid ??
  Neighbours(Neighbours),
  #[clap(name = "touniq")]
  ToUniq(ToUniq),
  #[clap(name = "toring", subcommand)]
  ToRing(ToRing),
  //#[clap(name = "cov")]
  //Coverage(Coverage),
}

impl Nested {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match self {
      Self::Table => Self::print_table(),
      Self::BestStartingDepth(e) => e.exec(),
      Self::Hash(e) => e.exec(),
      Self::Bilinear(e) => e.exec(),
      Self::Center(e) => e.exec(),
      Self::Vertices(e) => e.exec(),
      Self::Neighbours(e) => e.exec(),
      Self::ToUniq(e) => e.exec(),
      Self::ToRing(e) => e.exec(),
      // Self::Coverage(e) => e.exec(),
    }
  }

  fn print_table() -> Result<(), Box<dyn Error>> {
    const ARCMIN: f64 = 60.0;
    const ARCSEC: f64 = 3_600.0;
    const MAS: f64 = 3_600_000.0;
    const UMAS: f64 = 3_600_000_000.0;
    let mut handle = std::io::stdout().lock();
    writeln!(&mut handle, "HEALPix NESTED layers info:")?;
    writeln!(
      &mut handle,
      "{:>5} {:>10} {:>20} {:>12}",
      "depth", "nside", "ncells", "cellSize"
    )?;
    for depth in 0..=DEPTH_MAX {
      let ncells = n_hash(depth);
      let area_rad2 = std::f64::consts::PI / ((ncells >> 2) as f64);
      let side_deg = area_rad2.sqrt().to_degrees();
      if side_deg < 1.0 / MAS {
        writeln!(
          &mut handle,
          "   {:2} {:10} {:20} {:8.4} μas",
          depth,
          1_u64 << depth,
          ncells,
          side_deg * UMAS
        )?;
      } else if side_deg < 1.0 / ARCSEC {
        writeln!(
          &mut handle,
          "   {:2} {:10} {:20} {:8.4} mas",
          depth,
          1_u64 << depth,
          ncells,
          side_deg * MAS
        )?;
      } else if side_deg < 1.0 / ARCMIN {
        writeln!(
          &mut handle,
          "   {:2} {:10} {:20} {:8.4} ″  ",
          depth,
          1_u64 << depth,
          ncells,
          side_deg * ARCSEC
        )?;
      } else if side_deg < 1.0 {
        writeln!(
          &mut handle,
          "   {:2} {:10} {:20} {:8.4} ′  ",
          depth,
          1_u64 << depth,
          ncells,
          side_deg * ARCMIN
        )?;
      } else {
        writeln!(
          &mut handle,
          "   {:2} {:10} {:20} {:8.4} °  ",
          depth,
          1_u64 << depth,
          ncells,
          side_deg
        )?;
      }
    }
    Ok(())
  }
}

/// For a cone of given radius, returns the largest depth at which the cone covers max 9 cells.
#[derive(Debug, Args)]
pub struct BestStartingDepth {
  /// Cone radius, in decimal degrees.
  #[clap(value_name = "RADIUS")]
  radius: f64,
  /// Input radius in radians.
  #[clap(long, conflicts_with = "arcmin")]
  radian: bool,
  /// Input radius in arcminutes.
  #[clap(long, conflicts_with = "arcsec")]
  arcmin: bool,
  /// Input radius in arcseconds.
  #[clap(long, conflicts_with = "radian")]
  arcsec: bool,
}
impl BestStartingDepth {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let radius_rad = if self.radian {
      self.radius
    } else if self.arcmin {
      (self.radius / 60.0).to_radians()
    } else if self.arcsec {
      (self.radius / 3600.0).to_radians()
    } else {
      self.radius.to_radians()
    };
    if has_best_starting_depth(radius_rad) {
      let bsd = best_starting_depth(radius_rad);
      write!(stdout(), "{}\n", bsd)
    } else {
      write!(stdout(), "0\n")
    }
    .map_err(|e| e.into())
  }
}

/// Get the Healpix hash value (aka ipix or cell number) of given equatorial position(s)
#[derive(Debug, Args)]
pub struct Hash {
  /// Depth (aka order) of the Healpix layer
  depth: u8,
  /// Additionally provides the offsets (dx, dy) ∊ ([0..1[, [0..1[) of the position inside the cell.
  #[clap(long)]
  with_dxdy: bool,
  /// Compute the hash using the HEALPix UNIQ notation
  #[clap(long, conflicts_with = "suniq")]
  uniq: bool,
  /// Compute the hash using the MSB sentinel bit UNIQ notation
  #[clap(long, conflicts_with = "zuniq")]
  suniq: bool,
  /// Compute the hash using the LSB sentinel bit UNIQ notation
  #[clap(long, conflicts_with = "uniq")]
  zuniq: bool,
  #[command(subcommand)]
  pos: PosInput,
}
impl Hash {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let width = if self.zuniq {
      20
    } else {
      1 + n_hash(self.depth).ilog10() as usize
    };
    let htype = (self.uniq as u8) | ((self.suniq as u8) << 1) | ((self.zuniq as u8) << 2);
    let tr = |h: u64| -> u64 {
      match htype {
        0 => h,
        1 => to_uniq_ivoa(self.depth, h),
        2 => to_uniq(self.depth, h),
        4 => to_zuniq(self.depth, h),
        _ => unreachable!(),
      }
    };
    if !self.with_dxdy {
      self.pos.exec(
        move |lon_rad: f64,
              lat_rad: f64,
              write: &mut StdoutLock<'static>,
              _sep: char|
              -> Result<(), Box<dyn Error>> {
          let h = hash(self.depth, lon_rad, lat_rad);
          write!(write, "{:width$}", tr(h)).map_err(|e| e.into())
        },
      )
    } else {
      self.pos.exec(
        move |lon_rad: f64,
              lat_rad: f64,
              write: &mut StdoutLock<'static>,
              sep: char|
              -> Result<(), Box<dyn Error>> {
          let (h, dx, dy) = hash_with_dxdy(self.depth, lon_rad, lat_rad);
          write!(
            write,
            "{:width$}{}{:17.15}{}{:17.15}",
            tr(h),
            sep,
            dx,
            sep,
            dy,
          )
          .map_err(|e| e.into())
        },
      )
    }
  }
}

/// Get the bilinear interpolation parameters, 4x(hash, weight), of given equatorial position(s)
#[derive(Debug, Args)]
pub struct Bilinear {
  /// Depth (aka order) of the Healpix layer
  depth: u8,
  #[command(subcommand)]
  pos: PosInput,
}
impl Bilinear {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    let width = 1 + n_hash(self.depth).ilog10() as usize;
    self.pos.exec(
      |lon_rad: f64,
       lat_rad: f64,
       write: &mut StdoutLock<'static>,
       sep: char|
       -> Result<(), Box<dyn Error>> {
        let [(h1, w1), (h2, w2), (h3, w3), (h4, w4)] =
          bilinear_interpolation(self.depth, lon_rad, lat_rad);
        write!(
          write,
          "{:width$}{}{:17.15}{}{:width$}{}{:17.15}{}{:width$}{}{:17.15}{}{:width$}{}{:17.15}",
          h1, sep, w1, sep, h2, sep, w2, sep, h3, sep, w3, sep, h4, sep, w4
        )
        .map_err(|e| e.into())
      },
    )
  }
}

/// Get the longitude and latitude (in degrees) of the center of the cell of given hash value.
#[derive(Debug, Args)]
pub struct Center {
  /// Depth (aka order) of the Healpix layer
  depth: u8,
  #[command(subcommand)]
  hash: HashInput,
}
impl Center {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    self.hash.exec(
      |hash: u64, write: &mut StdoutLock<'static>, sep: char| -> Result<(), Box<dyn Error>> {
        let (lon, lat) = center(self.depth, hash);
        write!(
          write,
          "{:017.13}{}{:+017.13}",
          lon.to_degrees(),
          sep,
          lat.to_degrees()
        )
        .map_err(|e| e.into())
      },
    )
  }
}

/// Get the longitude and latitude (in degrees) of the 4 vertices (S, W, N, E) of the cell of given hash value.
#[derive(Debug, Args)]
pub struct Vertices {
  /// Depth (aka order) of the Healpix layer
  depth: u8,
  #[command(subcommand)]
  hash: HashInput,
}
impl Vertices {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    self.hash.exec(
      |hash: u64, write: &mut StdoutLock<'static>, sep: char| -> Result<(), Box<dyn Error>> {
        let [(l1, b1), (l2, b2), (l3, b3), (l4, b4)] = vertices(self.depth, hash);
        write!(
          write,
          "{:016.12}{}{:+016.12}{}{:016.12}{}{:+016.12}{}{:016.12}{}{:+016.12}{}{:016.12}{}{:+016.12}",
          l1.to_degrees(),
          sep,
          b1.to_degrees(),
          sep,
          l2.to_degrees(),
          sep,
          b2.to_degrees(),
          sep,
          l3.to_degrees(),
          sep,
          b3.to_degrees(),
          sep,
          l4.to_degrees(),
          sep,
          b4.to_degrees()
        )
          .map_err(|e| e.into())
      },
    )
  }
}

/// Get the hash value (aka ipix or cell number) of the neighbours of the cell of given hash value (S, SE, E, SW, NE, W, NW, N).
#[derive(Debug, Args)]
pub struct Neighbours {
  /// Depth (aka order) of the Healpix layer
  depth: u8,
  #[command(subcommand)]
  hash: HashInput,
}
impl Neighbours {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    self.hash.exec(
      |hash: u64, write: &mut StdoutLock<'static>, sep: char| -> Result<(), Box<dyn Error>> {
        let neigh = neighbours(self.depth, hash, false);
        if let Some(h) = neigh.get(MainWind::S) {
          write!(write, "{}", h)?;
        }
        for direction in [
          MainWind::SE,
          MainWind::E,
          MainWind::SW,
          MainWind::NE,
          MainWind::W,
          MainWind::NW,
          MainWind::N,
        ] {
          match neigh.get(direction) {
            Some(h) => write!(write, "{}{}", sep, h),
            None => write!(write, "{}", sep),
          }?;
        }
        Ok(())
      },
    )
  }
}

/// Transform regular hash values (aka ipix or cell numbers) into the UNIQ representation.
#[derive(Debug, Args)]
pub struct ToUniq {
  /// Transform into the MSB sentinel bit UNIQ notation instead of the regular Healpix UNIQ (the SUNIQ notation is valid for other MOC quantities).
  #[clap(long, conflicts_with = "zuniq")]
  suniq: bool,
  /// Transform into the LSB sentinel bit UNIQ notation instead of the regular Healpix UNIQ (ZUNIQ natural ordering follow the z-order curve: the order does not depend on the depth of each cell).
  #[clap(long)]
  zuniq: bool,
  #[command(subcommand)]
  fixed_or_variable: FixedOrVariableDepth,
}
impl ToUniq {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match (self.fixed_or_variable, self.suniq, self.zuniq) {
      (FixedOrVariableDepth::FixedDepth(FixedDepth { depth, input }), false, false) => input.exec(
        |hash: u64, write: &mut StdoutLock<'static>, _sep: char| -> Result<(), Box<dyn Error>> {
          write!(write, "{}", to_uniq_ivoa(depth, hash)).map_err(|e| e.into())
        },
      ),
      (FixedOrVariableDepth::FixedDepth(FixedDepth { depth, input }), true, false) => input.exec(
        |hash: u64, write: &mut StdoutLock<'static>, _sep: char| -> Result<(), Box<dyn Error>> {
          write!(write, "{}", to_uniq(depth, hash)).map_err(|e| e.into())
        },
      ),
      (FixedOrVariableDepth::FixedDepth(FixedDepth { depth, input }), false, true) => input.exec(
        |hash: u64, write: &mut StdoutLock<'static>, _sep: char| -> Result<(), Box<dyn Error>> {
          write!(write, "{}", to_zuniq(depth, hash)).map_err(|e| e.into())
        },
      ),
      (FixedOrVariableDepth::VarDepth(e), false, false) => e.exec(
        |depth: u8,
         hash: u64,
         write: &mut StdoutLock<'static>,
         _sep: char|
         -> Result<(), Box<dyn Error>> {
          write!(write, "{}", to_uniq_ivoa(depth, hash)).map_err(|e| e.into())
        },
      ),
      (FixedOrVariableDepth::VarDepth(e), true, false) => e.exec(
        |depth: u8,
         hash: u64,
         write: &mut StdoutLock<'static>,
         _sep: char|
         -> Result<(), Box<dyn Error>> {
          write!(write, "{}", to_uniq(depth, hash)).map_err(|e| e.into())
        },
      ),
      (FixedOrVariableDepth::VarDepth(e), false, true) => e.exec(
        |depth: u8,
         hash: u64,
         write: &mut StdoutLock<'static>,
         _sep: char|
         -> Result<(), Box<dyn Error>> {
          write!(write, "{}", to_zuniq(depth, hash)).map_err(|e| e.into())
        },
      ),
      _ => unreachable!(), // since with use 'conflicts_with'
    }
  }
}

#[derive(Debug, Subcommand)]
pub enum FixedOrVariableDepth {
  /// From fixed depth hash values.
  #[clap(name = "singledepth")]
  FixedDepth(FixedDepth),
  /// From hash values at various depth.
  #[clap(name = "multidepth", subcommand)]
  VarDepth(DepthHashInput),
}

/// Transform regular hash values (aka ipix or cell numbers) from the NESTED to the RING scheme.
#[derive(Debug, Subcommand)]
pub enum ToRing {
  /// Transform fixed depth hash values from the NESTED to the RING scheme.
  #[clap(name = "singledepth")]
  FixedDepth(FixedDepth),
  /// Transform hash values at various depth from the NESTED to the RING scheme.
  #[clap(name = "multidepth", subcommand)]
  VarDepth(DepthHashInput),
}
impl ToRing {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match self {
      Self::FixedDepth(FixedDepth { depth, input }) => {
        let layer = get(depth);
        input.exec(
          |hash: u64, write: &mut StdoutLock<'static>, _sep: char| -> Result<(), Box<dyn Error>> {
            write!(write, "{}", layer.to_ring(hash)).map_err(|e| e.into())
          },
        )
      }
      Self::VarDepth(e) => e.exec(
        |depth: u8,
         hash: u64,
         write: &mut StdoutLock<'static>,
         _sep: char|
         -> Result<(), Box<dyn Error>> {
          write!(write, "{}", get(depth).to_ring(hash)).map_err(|e| e.into())
        },
      ),
    }
  }
}

#[derive(Debug, Args)]
pub struct FixedDepth {
  /// Depth (aka order) of the Healpix layer
  depth: u8,
  #[command(subcommand)]
  input: HashInput,
}
