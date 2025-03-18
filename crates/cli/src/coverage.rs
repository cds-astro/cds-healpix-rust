use std::{
  error::Error,
  f64::consts::{FRAC_PI_2, PI},
  fmt::{self, Display, Formatter},
  fs::{self, File},
  io::{BufWriter, Write},
  num::ParseFloatError,
  path::PathBuf,
  str::FromStr,
};

use clap::{Args, Subcommand, ValueEnum};
use nom::{
  error::{convert_error, VerboseError},
  Err,
};
use thiserror::Error;

use stc_s::{
  space::{
    common::{
      region::{BoxParams, CircleParams, ConvexParams, EllipseParams, PolygonParams},
      FillFrameRefposFlavor, Flavor, Frame, FromPosToVelocity, SpaceUnit,
    },
    position::Position,
    positioninterval::PositionInterval,
  },
  visitor::{impls::donothing::VoidVisitor, CompoundVisitor, SpaceVisitor, StcVisitResult},
  Stc,
};

use hpxlib::{
  nested::{
    bmoc::BMOC, box_coverage, cone_coverage_approx_custom, elliptical_cone_coverage_custom,
    polygon_coverage, ring_coverage_approx_custom, zone_coverage,
  },
  sph_geom::frame::RefToLocalRotMatrix,
  TWICE_PI,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum BMocOutputType {
  /// CSV serialization.
  CSV,
  // /// ASCII serialization.
  // ASCII,
  /// FITS serialization.
  FITS,
  /// FITS BINTABLE serialization.
  BINTABLE,
}
impl Display for BMocOutputType {
  fn fmt(&self, f: &mut Formatter) -> fmt::Result {
    match self {
      Self::CSV => f.write_str("csv"),
      // Self::ASCII => f.write_str("ascii"),
      Self::FITS => f.write_str("fits"),
      Self::BINTABLE => f.write_str("bintable"),
    }
  }
}

/// Compute the list of NESTED Healpix cells in a given sky region
/// (see also moc-cli).
#[derive(Debug, Args)]
pub struct Coverage {
  /// Maximum depth (aka order) of the Healpix cells in the coverage.
  depth: u8,
  /// Output type
  #[clap(short = 't', long = "out-type", default_value_t = BMocOutputType::CSV)]
  output_type: BMocOutputType,
  /// Write in the given output file instead of stdout (mandatory for FITS).
  #[clap(short, long, value_name = "FILE")]
  output_file: Option<PathBuf>,
  #[clap(subcommand)]
  cov: CoverageType,
}

impl Coverage {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match (self.output_file, self.output_type) {
      (None, BMocOutputType::CSV) => {
        let writer = std::io::stdout().lock();
        self.cov.exec(self.depth, self.output_type, writer)
      }
      (None, out_type) => Err(format!("Sysout output not posible for: {}", out_type).into()),
      (Some(path), out_type) => {
        let writer = File::create(path).map(BufWriter::new)?;
        self.cov.exec(self.depth, out_type, writer)
      }
    }
  }
}

// MAKE A BMOC SERIALIZATION?!!
// CSV: hierach / flat
// FITS?!

/*
BMOC output with options: hierarch/flat/with or without flag/...

*/

#[derive(Debug, Subcommand)]
pub enum CoverageType {
  #[clap(name = "cone")]
  /// Cone coverage
  Cone {
    /// Longitude of the cone center (in degrees)
    lon_deg: f64,
    /// Latitude of the cone center (in degrees)
    #[clap(allow_negative_numbers = true)]
    lat_deg: f64,
    /// Radius of the cone (in degrees)
    r_deg: f64,
    /// Extra depth at which the computations are made to be more precise
    #[clap(short = 'd', long, default_value_t = 2)]
    delta_depth: u8,
  },
  #[clap(name = "ellipse")]
  ///Elliptical cone coverage
  EllipticalCone {
    /// Longitude of the elliptical cone center (in degrees)
    lon_deg: f64,
    /// Latitude of the elliptical cone center (in degrees)
    #[clap(allow_negative_numbers = true)]
    lat_deg: f64,
    /// Elliptical cone semi-major axis (in degrees)
    a_deg: f64,
    /// Elliptical cone semi-minor axis (in degrees)
    b_deg: f64,
    /// Elliptical cone position angle (in degrees)
    pa_deg: f64,
    /// Extra depth at which the computations are made to be more precise
    #[clap(short = 'd', long, default_value_t = 2)]
    delta_depth: u8,
  },
  #[clap(name = "ring")]
  /// Ring coverage
  Ring {
    /// Longitude of the ring center (in degrees)
    lon_deg: f64,
    /// Latitude of the ring center (in degrees)
    #[clap(allow_negative_numbers = true)]
    lat_deg: f64,
    /// External radius of the ring (in degrees)
    r_max_deg: f64,
    /// Internal radius of the ring (in degrees)
    r_min_deg: f64,
    /// Extra depth at which the computations are made to be more precise
    #[clap(short = 'd', long, default_value_t = 2)]
    delta_depth: u8,
  },
  #[clap(name = "zone")]
  /// Zone coverage
  Zone {
    /// Longitude min, in degrees
    lon_deg_min: f64,
    /// Latitude min, in degrees
    #[clap(allow_negative_numbers = true)]
    lat_deg_min: f64,
    /// Longitude max, in degrees
    lon_deg_max: f64,
    /// Latitude max, in degrees
    #[clap(allow_negative_numbers = true)]
    lat_deg_max: f64,
  },
  #[clap(name = "box")]
  /// Box coverage
  Box {
    /// Longitude of the box center, in degrees
    lon_deg: f64,
    /// Latitude of the box center, in degrees
    #[clap(allow_negative_numbers = true)]
    lat_deg: f64,
    /// Semi-major axis of the box, in degrees
    a_deg: f64,
    /// Semi-minor axis of the box, in degrees
    b_deg: f64,
    /// Position angle of the box, in degrees
    pa_deg: f64,
  },
  #[clap(name = "polygon")]
  /// Polygon coverage
  Polygon {
    /// List of vertices: "(lon,lat),(lon,lat),...,(lon,lat)" in degrees
    vertices_deg: Vertices,
  },
  #[clap(name = "stcs")]
  /// STC-S region coverage.
  Stcs {
    /// STC-S region
    stcs: String,
  },
  #[clap(name = "stcsfile")]
  /// STC-S region provided in a file coverage.
  StcsFile {
    #[clap(value_name = "FILE")]
    /// Path of the file containing the STC-S region
    file_path: PathBuf,
  },
}

impl CoverageType {
  pub fn exec<W: Write>(
    self,
    depth: u8,
    out_type: BMocOutputType,
    writer: W,
  ) -> Result<(), Box<dyn Error>> {
    match self {
      Self::Cone {
        lon_deg,
        lat_deg,
        r_deg,
        delta_depth,
      } => {
        let lon = lon_deg2rad(lon_deg)?;
        let lat = lat_deg2rad(lat_deg)?;
        let r = r_deg.to_radians();
        if r <= 0.0 || PI <= r {
          Err::<_, Box<dyn Error>>(String::from("Radius must be in ]0, pi[").into())
        } else {
          Ok(cone_coverage_approx_custom(depth, delta_depth, lon, lat, r))
        }
      }
      Self::EllipticalCone {
        lon_deg,
        lat_deg,
        a_deg,
        b_deg,
        pa_deg,
        delta_depth,
      } => {
        let lon = lon_deg2rad(lon_deg)?;
        let lat = lat_deg2rad(lat_deg)?;
        let a = a_deg.to_radians();
        let b = b_deg.to_radians();
        let pa = pa_deg.to_radians();
        if a <= 0.0 || FRAC_PI_2 <= a {
          Err(String::from("Semi-major axis must be in ]0, pi/2]").into())
        } else if b <= 0.0 || a <= b {
          Err(String::from("Semi-minor axis must be in ]0, a[").into())
        } else if pa <= 0.0 || PI <= pa {
          Err(String::from("Position angle must be in [0, pi[").into())
        } else {
          Ok(elliptical_cone_coverage_custom(
            depth,
            delta_depth,
            lon,
            lat,
            a,
            b,
            pa,
          ))
        }
      }
      Self::Ring {
        lon_deg,
        lat_deg,
        r_max_deg,
        r_min_deg,
        delta_depth,
      } => {
        let lon = lon_deg2rad(lon_deg)?;
        let lat = lat_deg2rad(lat_deg)?;
        let r_min = r_min_deg.to_radians();
        let r_max = r_max_deg.to_radians();
        if r_min <= 0.0 || r_max <= 0.0 || PI <= r_min || PI <= r_max {
          Err(String::from("Radius must be in ]0, pi[").into())
        } else {
          Ok(ring_coverage_approx_custom(
            depth,
            delta_depth,
            lon,
            lat,
            r_min,
            r_max,
          ))
        }
      }
      Self::Zone {
        lon_deg_min,
        lat_deg_min,
        lon_deg_max,
        lat_deg_max,
      } => {
        let lon_min = lon_deg2rad(lon_deg_min)?;
        let lat_min = lat_deg2rad(lat_deg_min)?;
        let lon_max = lon_deg2rad(lon_deg_max)?;
        let lat_max = lat_deg2rad(lat_deg_max)?;
        Ok(zone_coverage(depth, lon_min, lat_min, lon_max, lat_max))
      }
      Self::Box {
        lon_deg,
        lat_deg,
        a_deg,
        b_deg,
        pa_deg,
      } => {
        let lon = lon_deg2rad(lon_deg)?;
        let lat = lat_deg2rad(lat_deg)?;
        let a = a_deg.to_radians();
        let b = b_deg.to_radians();
        let pa = pa_deg.to_radians();
        if a <= 0.0 || FRAC_PI_2 <= a {
          Err(String::from("Semi-major axis must be in ]0, pi/2]").into())
        } else if b <= 0.0 || a <= b {
          Err(String::from("Semi-minor axis must be in ]0, a[").into())
        } else if pa <= 0.0 || FRAC_PI_2 <= pa {
          Err(String::from("Position angle must be in [0, pi[").into())
        } else {
          Ok(box_coverage(depth, lon, lat, a, b, pa))
        }
      }
      Self::Polygon { vertices_deg } => Ok(polygon_coverage(
        depth,
        vertices_deg
          .list
          .iter()
          .map(|(lon_deg, lat_deg)| {
            let lon = lon_deg2rad(*lon_deg)?;
            let lat = lat_deg2rad(*lat_deg)?;
            Ok((lon, lat))
          })
          .collect::<Result<Vec<(f64, f64)>, String>>()?
          .as_slice(),
        true,
      )),
      Self::Stcs { stcs } => {
        let stcs_query = stcs2stcsquery(stcs.as_str()).map_err(|e| e.to_string())?;
        Ok(stcs_query.to_bmoc(depth))
      }
      Self::StcsFile { file_path } => {
        let stcs_query = fs::read_to_string(&file_path)
          .map_err(|e| format!("Error reading file '{:?}': {:?}", file_path, e))
          .and_then(|stcs| stcs2stcsquery(stcs.as_str()).map_err(|e| e.to_string()))?;
        Ok(stcs_query.to_bmoc(depth))
      }
    }
    .and_then(|bmoc| {
      match out_type {
        BMocOutputType::CSV => bmoc.to_csv(writer),
        BMocOutputType::FITS => bmoc.to_fits(writer),
        BMocOutputType::BINTABLE => bmoc.to_bintable_fits(writer),
      }
      .map_err(|e| e.into())
    })
  }
}

#[derive(Debug, Clone)]
pub struct Vertices {
  // (ra0,dec0),(ra1,dec1),...,(ran,decn)
  list: Vec<(f64, f64)>,
}

impl FromStr for Vertices {
  type Err = ParseFloatError;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    let list: Vec<f64> = s
      .replace(['(', ')'], "")
      .split(',')
      .map(|t| str::parse::<f64>(t.trim()))
      .collect::<Result<Vec<f64>, _>>()?;
    Ok(Vertices {
      list: list
        .iter()
        .step_by(2)
        .zip(list.iter().skip(1).step_by(2))
        .map(|(lon, lat)| (*lon, *lat))
        .collect(),
    })
  }
}

fn lon_deg2rad(lon_deg: f64) -> Result<f64, String> {
  let lon = lon_deg.to_radians();
  if !(0.0..TWICE_PI).contains(&lon) {
    Err(String::from("Longitude must be in [0, 2pi["))
  } else {
    Ok(lon)
  }
}

fn lat_deg2rad(lat_deg: f64) -> Result<f64, String> {
  let lat = lat_deg.to_radians();
  if !(-FRAC_PI_2..=FRAC_PI_2).contains(&lat) {
    Err(String::from("Latitude must be in [-pi/2, pi/2]"))
  } else {
    Ok(lat)
  }
}

#[derive(Debug)]
pub struct Cone {
  /// Center longitude, in `[0, 2\pi[` radians
  lon: f64,
  /// Center latitude in `[-\pi/2, \pi/2]` radians
  lat: f64,
  /// Radius, in radians
  r: f64,
}

impl Cone {
  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2[`
  /// * if `r` not in `]0, \pi[`
  pub fn new(lon: f64, lat: f64, r: f64) -> Self {
    assert!(0.0 <= lon && lon < TWICE_PI);
    assert!(-FRAC_PI_2 <= lat && lat <= FRAC_PI_2);
    assert!(0.0 < r && r < PI);
    Self { lon, lat, r }
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    cone_coverage_approx_custom(depth, 3, self.lon, self.lat, self.r)
  }
}

#[derive(Debug)]
pub struct EllipticalCone {
  /// Center longitude, in `[0, 2\pi[` radians
  lon: f64,
  /// Center latitude in `[-\pi/2, \pi/2]` radians
  lat: f64,
  /// Semi-major axis, in radians
  a: f64,
  /// Semi-minor axis, in radians
  b: f64,
  /// Positional angle (east-of-north), in `[0, \pi[` radians
  pa: f64,
}

impl EllipticalCone {
  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2]`
  /// * if `a` not in `]0, \pi/2]` or \`]0, \pi/2]`
  /// * if `b` not in `]0, a[`
  /// * if `pa` not in `[0, \pi[`
  pub fn new(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> EllipticalCone {
    assert!(0.0 <= lon && lon < TWICE_PI);
    assert!(-FRAC_PI_2 <= lat && lat <= FRAC_PI_2);
    assert!(0.0 < a && a <= FRAC_PI_2);
    assert!(0.0 < b && b < a);
    assert!(0.0 < pa && pa < PI);
    EllipticalCone { lon, lat, a, b, pa }
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    elliptical_cone_coverage_custom(depth, 3, self.lon, self.lat, self.a, self.b, self.pa)
  }
}

#[derive(Debug)]
pub struct Zone {
  /// Minimal longitude (inclusive), in `[0, 2\pi[` radians
  lon_min: f64,
  /// Minimal latitude (inclusive if positive, else exclusive) in `[-\pi/2, \pi/2[` radians
  lat_min: f64,
  /// Maximal longitude (exclusive), in `[0, 2\pi]` radians
  lon_max: f64,
  /// Maximal latitude (exlusive if positive, else inclusive) in `]-\pi/2, \pi/2]` radians
  lat_max: f64,
}

impl Zone {
  /// # Remarks
  /// * If `lon_min > lon_max` then we consider that the zone crosses the primary meridian.
  /// * The north pole is included only if `lon_min == 0 && lat_max == pi/2`
  /// # Panics
  /// * if `lon_min` or `lon_max` not in `[0, 2\pi[`
  /// * if `lat_min` or `lat_max` not in `[-\pi/2, \pi/2[`
  /// * `lat_min >= lat_max`
  pub fn new(lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64) -> Zone {
    assert!((0.0..TWICE_PI).contains(&lon_min) && 0.0 < lon_max && lon_max <= TWICE_PI);
    assert!(
      (-FRAC_PI_2..FRAC_PI_2).contains(&lat_min) && -FRAC_PI_2 < lat_max && lat_max <= FRAC_PI_2
    );
    assert!(lat_min <= lat_max);
    Zone {
      lon_min,
      lat_min,
      lon_max,
      lat_max,
    }
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    zone_coverage(
      depth,
      self.lon_min,
      self.lat_min,
      self.lon_max,
      self.lat_max,
    )
  }
}

#[derive(Debug)]
pub struct Polygon {
  vertices: Vec<(f64, f64)>,
}

impl Polygon {
  pub fn new(vertices: Vec<(f64, f64)>) -> Self {
    assert!(vertices.len() > 2);
    Self { vertices }
  }

  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2]`
  /// * if `a` not in `]0, \pi/2]` or \`]0, \pi/2]`
  /// * if `b` not in `]0, a[`
  /// * if `pa` not in `[0, \pi[`
  pub fn from_box(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> Polygon {
    assert!((0.0..TWICE_PI).contains(&lon));
    assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&lat));
    assert!(0.0 < a && a <= FRAC_PI_2);
    assert!(0.0 < b && b < a);
    assert!((0.0..PI).contains(&pa));
    // Compute spherical coordinates
    let frame_rotation = RefToLocalRotMatrix::from_center(lon, lat);
    // By application of the Thales theorem, the new point has the property:
    //   sin(new_lat) / sin(dlat) = (cos(dlon) * cos(new_lat)) / cos(dlat)
    // With (imagine looking a the triangle in the (x, z) plane:
    // * old_x = cos(lat)
    // * new_x = cos(dlon) * cos(new_lat)
    // * old_z = sin(dlat)
    // * new_z = sin(new_lat)
    // Leading to:
    //   tan(new_lat) = cos(dlon) * tan(dlat)
    let lon = a;
    let (sin_lon, cos_lon) = lon.sin_cos();
    let lat = (cos_lon * b.tan()).atan();
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_pa, cos_pa) = pa.sin_cos();
    // Rotation by the position angle
    // - upper right (before rotation by PA)
    let (x1, y1, z1) = (cos_lon * cos_lat, sin_lon * cos_lat, sin_lat);
    // - apply rotation (sin and cos are revere since theta = pi/2 - pa)
    let (y2, z2) = (y1 * sin_pa - z1 * cos_pa, y1 * cos_pa + z1 * sin_pa);
    let mut vertices = Vec::with_capacity(4);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - lower right (before rotation by PA) = (y1, -z1)
    let (y2, z2) = (y1 * sin_pa + z1 * cos_pa, y1 * cos_pa - z1 * sin_pa);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - lower left (before rotation by PA) = (-y1, -z1)
    let (y2, z2) = (-y1 * sin_pa + z1 * cos_pa, -y1 * cos_pa - z1 * sin_pa);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - upper left (before rotation by PA) = (-y1, z1)
    let (y2, z2) = (-y1 * sin_pa - z1 * cos_pa, -y1 * cos_pa + z1 * sin_pa);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    Polygon::new(vertices)
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    polygon_coverage(depth, self.vertices.as_slice(), true)
  }
}

#[derive(Debug)]
pub enum StcsElem {
  AllSky,
  Cone(Cone),
  Ellipse(EllipticalCone),
  Zone(Zone),
  Polygon(Polygon), // Box => Polygon
  Not {
    elem: Box<StcsElem>,
  },
  Union {
    elems: Vec<StcsElem>,
  },
  Intersection {
    elems: Vec<StcsElem>,
  },
  Difference {
    left: Box<StcsElem>,
    right: Box<StcsElem>,
  },
}
impl StcsElem {
  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    match self {
      Self::AllSky => BMOC::new_allsky(depth),
      Self::Cone(elem) => elem.to_bmoc(depth),
      Self::Ellipse(elem) => elem.to_bmoc(depth),
      Self::Zone(elem) => elem.to_bmoc(depth),
      Self::Polygon(elem) => elem.to_bmoc(depth),
      Self::Not { elem } => elem.to_bmoc(depth).not(),
      Self::Union { elems } => elems
        .iter()
        .map(|elem| elem.to_bmoc(depth))
        .reduce(|acc, curr| acc.or(&curr))
        .unwrap_or(BMOC::new_empty(depth)),
      Self::Intersection { elems } => elems
        .iter()
        .map(|elem| elem.to_bmoc(depth))
        .reduce(|acc, curr| acc.and(&curr))
        .unwrap_or(BMOC::new_empty(depth)),
      Self::Difference { left, right } => left.to_bmoc(depth).xor(&right.to_bmoc(depth)),
    }
  }
}

#[derive(Debug)]
pub struct StcQuery {
  stc: StcsElem,
}
impl StcQuery {
  pub fn to_bmoc(self, depth: u8) -> BMOC {
    self.stc.to_bmoc(depth)
  }
}

#[derive(Error, Debug)]
pub enum Stc2StcQueryError {
  #[error("Frame other than ICRS not supported (yet). Found: {found:?}")]
  FrameIsNotICRS { found: Frame },
  #[error("Flavor other than Spher2 not supported (yet). Found: {found:?}")]
  FlavorIsNotSpher2 { found: Flavor },
  #[error("Units ther than 'deg' not (yet?!) supported. Found: {found:?}")]
  UnitsNotSupported { found: Vec<SpaceUnit> },
  #[error("Convex shape not (yet?!) supported.")]
  ConvexNotSupported,
  #[error("Simple position not supported.")]
  SimplePositionNotSupported,
  #[error("Empty position range not supported.")]
  EmptyPositionRangeNotSupported,
  #[error("Position interval not supported.")]
  PositionIntervalNotSupported,
  #[error("invalid header (expected {expected:?}, found {found:?})")]
  WrongNumberOfParams { expected: u8, found: u8 },
  #[error("Longitude value out of bounds. Expected: [0, 360[. Actual: {value:?}")]
  WrongLongitude { value: f64 },
  #[error("Latitude value out of bounds. Expected: [-90, 90[. Actual: {value:?}")]
  WrongLatitude { value: f64 },
  #[error("STC-S string parsing not complete. Remaining: {rem:?}")]
  ParseHasRemaining { rem: String },
  #[error("STC-S string parsing incomplete: {msg:?}")]
  ParseIncomplete { msg: String },
  #[error("STC-S string parsing error: {msg:?}")]
  ParseFailure { msg: String },
  #[error("STC-S string parsing failure: {msg:?}")]
  ParseError { msg: String },
  #[error("No space sub-phrase found in STC-S string")]
  NoSpaceFound,
  #[error("Custom error: {msg:?}")]
  Custom { msg: String },
}

#[derive(Debug, Clone)]
struct Stc2StcQuery;

impl CompoundVisitor for Stc2StcQuery {
  type Value = StcsElem;
  type Error = Stc2StcQueryError;

  fn visit_allsky(&mut self) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::AllSky)
  }

  fn visit_circle(&mut self, circle: &CircleParams) -> Result<Self::Value, Self::Error> {
    // Get params
    let lon_deg = circle
      .center()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty circle longitude"),
      })?;
    let lat_deg = circle
      .center()
      .get(1)
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty circle latitude"),
      })?;
    let radius_deg = circle.radius();
    // Convert params
    let lon = lon_deg2rad_stcs(*lon_deg)?;
    let lat = lat_deg2rad_stcs(*lat_deg)?;
    let r = radius_deg.to_radians();
    if r <= 0.0 || PI <= r {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Radius out of bounds in  CIRCLE. Expected: ]0, 180[. Actual: {}.",
          r
        ),
      })
    } else {
      Ok(StcsElem::Cone(Cone::new(lon, lat, r)))
    }
  }

  fn visit_ellipse(&mut self, ellipse: &EllipseParams) -> Result<Self::Value, Self::Error> {
    // Get params
    let lon_deg = ellipse
      .center()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse longitude"),
      })?;
    let lat_deg = ellipse
      .center()
      .get(1)
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse latitude"),
      })?;
    let a_deg = ellipse.radius_a();
    let b_deg = ellipse.radius_b();
    let pa_deg = ellipse.pos_angle();
    // Convert params
    let lon = lon_deg2rad_stcs(*lon_deg)?;
    let lat = lat_deg2rad_stcs(*lat_deg)?;
    let a = a_deg.to_radians();
    let b = b_deg.to_radians();
    let pa = pa_deg.to_radians();
    if a <= 0.0 || FRAC_PI_2 <= a {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Semi-major axis out of bounds. Expected: ]0, 90[. Actual: {}.",
          a_deg
        ),
      })
    } else if b <= 0.0 || a <= b {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Semi-minor axis out of bounds. Expected: ]0, {}[. Actual: {}.",
          a_deg, b_deg
        ),
      })
    } else if pa <= 0.0 || PI <= pa {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Position angle out of bounds. Expected: [0, 180[. Actual: {}.",
          pa_deg
        ),
      })
    } else {
      Ok(StcsElem::Ellipse(EllipticalCone::new(lon, lat, a, b, pa)))
    }
  }

  fn visit_box(&mut self, skybox: &BoxParams) -> Result<Self::Value, Self::Error> {
    // Get params
    let lon_deg = skybox
      .center()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse longitude"),
      })?;
    let lat_deg = skybox
      .center()
      .get(1)
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse latitude"),
      })?;
    let mut a_deg = skybox
      .bsize()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty bsize on latitude"),
      })?;
    let mut b_deg = skybox
      .bsize()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty bsize on longitude"),
      })?;
    let mut pa_deg = skybox.bsize().first().copied().unwrap_or(90.0);
    if a_deg < b_deg {
      std::mem::swap(&mut b_deg, &mut a_deg);
      pa_deg = 90.0 - pa_deg;
    }
    // Convert params
    let lon = lon_deg2rad_stcs(*lon_deg)?;
    let lat = lat_deg2rad_stcs(*lat_deg)?;
    let a = a_deg.to_radians();
    let b = b_deg.to_radians();
    let pa = pa_deg.to_radians();
    if a <= 0.0 || FRAC_PI_2 <= a {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Box semi-major axis out of bounds. Expected: ]0, 90[. Actual: {}.",
          a_deg
        ),
      })
    } else if b <= 0.0 || a <= b {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Box semi-minor axis out of bounds. Expected: ]0, {}[. Actual: {}.",
          a_deg, b_deg
        ),
      })
    } else if !(0.0..PI).contains(&pa) {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Position angle out of bounds. Expected: [0, 180[. Actual: {}.",
          pa_deg
        ),
      })
    } else {
      Ok(StcsElem::Polygon(Polygon::from_box(lon, lat, a, b, pa)))
    }
  }

  fn visit_polygon(&mut self, polygon: &PolygonParams) -> Result<Self::Value, Self::Error> {
    let vertices_deg = polygon.vertices();
    let vertices = vertices_deg
      .iter()
      .step_by(2)
      .zip(vertices_deg.iter().skip(1).step_by(2))
      .map(|(lon_deg, lat_deg)| {
        let lon = lon_deg2rad_stcs(*lon_deg)?;
        let lat = lat_deg2rad_stcs(*lat_deg)?;
        Ok((lon, lat))
      })
      .collect::<Result<Vec<(f64, f64)>, Stc2StcQueryError>>()?;
    // Do something to use the right convention (depending on vertex order)!!
    Ok(StcsElem::Polygon(Polygon::new(vertices)))
  }

  fn visit_convex(&mut self, _convex: &ConvexParams) -> Result<Self::Value, Self::Error> {
    Err(Stc2StcQueryError::ConvexNotSupported)
  }

  fn visit_not(&mut self, elem: Self::Value) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::Not {
      elem: Box::new(elem),
    })
  }

  fn visit_union(&mut self, elems: Vec<Self::Value>) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::Union { elems })
  }

  fn visit_intersection(&mut self, elems: Vec<Self::Value>) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::Intersection { elems })
  }

  fn visit_difference(
    &mut self,
    left: Self::Value,
    right: Self::Value,
  ) -> Result<Self::Value, Self::Error> {
    // Warning: we interpret 'difference' as being a 'symmetrical difference', i.e. xor (not minus)
    Ok(StcsElem::Difference {
      left: Box::new(left),
      right: Box::new(right),
    })
  }
}

impl SpaceVisitor for Stc2StcQuery {
  type Value = StcsElem;
  type Error = Stc2StcQueryError;
  type C = Self;

  fn new_compound_visitor(
    &self,
    fill_frame_refpos_flavor: &FillFrameRefposFlavor,
    from_pos_to_velocity: &FromPosToVelocity,
  ) -> Result<Self, Self::Error> {
    // Check ICRS frame
    let frame = fill_frame_refpos_flavor.frame();
    if frame != Frame::ICRS {
      return Err(Stc2StcQueryError::FrameIsNotICRS { found: frame });
    }
    // Check SPHER2 flavor
    let flavor = fill_frame_refpos_flavor.flavor();
    if let Some(flavor) = flavor {
      if flavor != Flavor::Spher2 {
        return Err(Stc2StcQueryError::FlavorIsNotSpher2 { found: flavor });
      }
    }
    // Check units
    let opt_units = from_pos_to_velocity.unit().cloned();
    if let Some(units) = opt_units {
      for unit in units.iter().cloned() {
        if unit != SpaceUnit::Deg {
          return Err(Stc2StcQueryError::UnitsNotSupported { found: units });
        }
      }
    }
    Ok(self.clone())
  }

  fn visit_position_simple(self, _: &Position) -> Result<Self::Value, Self::Error> {
    Err(Stc2StcQueryError::SimplePositionNotSupported)
  }

  fn visit_position_interval(
    self,
    interval: &PositionInterval,
  ) -> Result<Self::Value, Self::Error> {
    // We use compound visitor only to check interval parameters
    self.new_compound_visitor(&interval.pre, &interval.post)?;
    let corners = interval
      .lo_hi_limits
      .iter()
      .step_by(2)
      .zip(interval.lo_hi_limits.iter().skip(1).step_by(2))
      .map(|(lon_deg, lat_deg)| {
        let lon = lon_deg2rad_stcs(*lon_deg)?;
        let lat = lat_deg2rad_stcs(*lat_deg)?;
        Ok((lon, lat))
      })
      .collect::<Result<Vec<(f64, f64)>, Stc2StcQueryError>>()?;
    let corners_it = corners
      .iter()
      .cloned()
      .step_by(2)
      .zip(corners.iter().cloned().skip(1).step_by(2));
    let mut zones = corners_it
      .map(|((ra_min, dec_min), (ra_max, dec_max))| {
        StcsElem::Zone(Zone::new(ra_min, dec_min, ra_max, dec_max))
      })
      .collect::<Vec<StcsElem>>();
    match zones.len() {
      0 => Err(Stc2StcQueryError::EmptyPositionRangeNotSupported),
      1 => Ok(zones.drain(..).last().unwrap()), // Unwrap ok since we tested the number of elements
      _ => Ok(StcsElem::Union { elems: zones }),
    }
  }

  fn visit_allsky(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_circle(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_ellipse(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_box(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_polygon(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_convex(self, _: StcsElem) -> Result<Self::Value, Self::Error> {
    unreachable!() // because an error is raised before calling this
  }

  fn visit_not(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_union(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_intersection(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_difference(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }
}

fn lon_deg2rad_stcs(lon_deg: f64) -> Result<f64, Stc2StcQueryError> {
  let mut lon = lon_deg.to_radians();
  if lon == TWICE_PI {
    lon = 0.0;
  }
  if !(0.0..TWICE_PI).contains(&lon) {
    Err(Stc2StcQueryError::WrongLongitude { value: lon_deg })
  } else {
    Ok(lon)
  }
}

fn lat_deg2rad_stcs(lat_deg: f64) -> Result<f64, Stc2StcQueryError> {
  let lat = lat_deg.to_radians();
  if !(-FRAC_PI_2..=FRAC_PI_2).contains(&lat) {
    Err(Stc2StcQueryError::WrongLatitude { value: lat_deg })
  } else {
    Ok(lat)
  }
}

/// Create new StcsQuery from the given STC-S string.
///
/// # WARNING
/// * `DIFFERENCE` is interpreted as a symmetrical difference (it is a `MINUS` in the STC standard)
/// * `Polygon` do not follow the STC-S standard: here self-intersecting polygons are supported
/// * No implicit conversion: the STC-S will be rejected if
///     + the frame is different from `ICRS`
///     + the flavor is different from `Spher2`
///     + the units are different from `degrees`
/// * Time, Spectral and Redshift sub-phrases are ignored
///
/// # Params
/// * `ascii_stcs`: lthe STC-S string
///
/// # Output
/// - The new StcsQuery (or an error)
pub fn stcs2stcsquery(stcs_ascii: &str) -> Result<StcQuery, Stc2StcQueryError> {
  match Stc::parse::<VerboseError<&str>>(stcs_ascii.trim()) {
    Ok((rem, stcs)) => {
      if !rem.is_empty() {
        return Err(Stc2StcQueryError::ParseHasRemaining {
          rem: rem.to_string(),
        });
      }
      let StcVisitResult { space, .. } =
        stcs.accept(VoidVisitor, Stc2StcQuery, VoidVisitor, VoidVisitor);
      match space {
        None => Err(Stc2StcQueryError::NoSpaceFound),
        Some(space_res) => space_res.map(|stc| StcQuery { stc }),
      }
    }
    Err(err) => Err(match err {
      Err::Incomplete(_) => Stc2StcQueryError::ParseIncomplete {
        msg: String::from("Incomplete parsing."),
      },
      Err::Error(e) => Stc2StcQueryError::ParseIncomplete {
        msg: convert_error(stcs_ascii, e),
      },
      Err::Failure(e) => Stc2StcQueryError::ParseIncomplete {
        msg: convert_error(stcs_ascii, e),
      },
    }),
  }
}
