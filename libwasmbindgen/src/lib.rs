//! Interface for accesing the lib in WebAssembly.
//! The Javascript name (js_name) have been chosen with T. Boch

use std::collections::HashMap;

extern crate cdshealpix;
extern crate wasm_bindgen;

// see https://rustwasm.github.io/wasm-bindgen/
use wasm_bindgen::prelude::*;

use cdshealpix::nested::{hash, center, neighbours, vertices, cone_coverage_approx_flat, cone_coverage_approx, polygon_coverage};
use cdshealpix::nested::bmoc::*;
use cdshealpix::compass_point::{MainWind, MainWindMap};

use cdshealpix::sph_geom::Polygon;
use cdshealpix::sph_geom::coo3d::{Coo3D, LonLat};

const DEPHT_MAX: u8 = 24;

/// Returns the max order one can use in Javascipts.
/// Equals 24 since we store ou HEALPix indices on a Javascript number which is a double.
/// A double has a mantissa of 52 bits = 24 * 2 (2 bits per level) + 4 (4 bits to code base cells).
#[wasm_bindgen(js_name = getOrderMax)]
pub fn get_depth_max() -> u8 {
  DEPHT_MAX
}

#[wasm_bindgen(js_name = PolygonMap)]
pub struct PolygonMap {
  polygons: HashMap<u32, Polygon>,
}

#[wasm_bindgen(js_class = PolygonMap)]
impl PolygonMap {

  pub fn new() -> PolygonMap {
    PolygonMap {
      polygons: Default::default(),
    }
  }

  /// Add a polygon of given identifier to the map, the provided liste is:
  /// [ra_v1_deg, dec_v1_deg, ra_v2_deg, de_v2_deg, ..., ra_vn_deg, de_vn_deg]
  #[wasm_bindgen(js_name = addPolygon)]
  pub fn add_polygon(&mut self, id: u32, vertices_coos: Box<[f64]>) {
    let mut vertices: Vec<LonLat> = Vec::with_capacity(vertices_coos.len() >> 1);
    for i in (0..vertices_coos.len()).step_by(2) {
      vertices.push(LonLat { 
          lon: vertices_coos[i].to_radians(), 
          lat: vertices_coos[i + 1].to_radians(),
      });
    }
    self.polygons.insert(id, Polygon::new(vertices.into_boxed_slice()));
  }

  /// Remove the polygon of given identifier to the map.
  #[wasm_bindgen(js_name = rmPolygon)]
  pub fn remove_polygon(&mut self, id: u32) {
    self.polygons.remove(&id);
  }
  
  /// Returns the list of the identifiers of the polygon containing the given point.
  #[wasm_bindgen(js_name = polygonContaining)]
  pub fn polygon_containing(&self, ra_deg: f64, de_deg: f64) -> Box<[u32]> {
    let coo: Coo3D = Coo3D::from_sph_coo(ra_deg.to_radians(), de_deg.to_radians());
    let mut res: Vec<u32> = Default::default();
    for (id, polygon) in self.polygons.iter() {
      if polygon.contains(&coo) {
        res.push(*id);
      }
    }
    res.into_boxed_slice()
  }

  /*pub fn from(n_vertices: Box<[usize]>, vertices_coos: Box<[f64]>) -> PolygonSet {
    let mut polygons: Vec<Polygon> = Vec::with_capacity(n_vertices.len());
    let mut cumul_size = 0;
    for size in n_vertices.iter() {
    }
    PolygonSet {
      polygons
    }
  }*/


}

/// Returns the cell number (hash value) in the NESTED scheme associated with the given position
/// on the unit sphere
/// # Inputs
/// - `order`: the order of the HEALPix hash we want in output
/// - `lon`: longitude in degrees, support reasonably large positive and negative values
///          producing accurate results with a naive range reduction like modulo 360
///          (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
/// - `lat`: latitude in degrees, must be in `[-90, 90]`
/// # Output
/// - the cell number (hash value) associated with the given position on the unit sphere,
///   in `[0, 12*nside^2[`
/// # Panics
///   If `lat` **not in** `[-90, 90]`, this method panics.
#[wasm_bindgen(js_name = lonlatToNested)]
pub fn wasm_hash(depth: u8, lon: f64, lat: f64) -> f64 {
  hash(depth, lon.to_radians(), lat.to_radians()) as f64
}

/// Returns the cell numbers (hash value) in the NESTED scheme associated with the given positions 
/// on the unit sphere
/// # Inputs
/// - `order`: the order of the HEALPix hashes we want in output
/// - `coords`: an array storing consecutively the coordinates of the positions we look for the
///             cell numbers: `[lon_1, lat_1, lon_2, lat_2, ..., lon_n, lat_n]` 
///   - `lon`: longitude in degrees, support reasonably large positive and negative values
///            producing accurate results with a naive range reduction like modulo 2*pi
///            (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
///   - `lat`: latitude in degrees, must be in `[-90, 90]`
/// # Output
/// - the cell number (hash value) associated with the given position on the unit sphere,
///   in `[0, 12*nside^2[`
/// # Panics
///   If `lat` **not in** `[-90, 90]`, this method panics.
#[wasm_bindgen(js_name = multiLonlatToNested)]
pub fn wasm_hash_multi(depth: u8, coords: &[f64]) -> Box<[f64]> {
  let mut vec: Vec<f64> = Vec::with_capacity(coords.len() >> 1);
  for i in (0..coords.len()).step_by(2) {
    vec.push(hash(depth, coords[i].to_radians(), coords[i + 1].to_radians()) as f64);
  }
  vec.into_boxed_slice()
}
// [[64; 2]] a tester

/*
/// Defines a set of spherical coordinates, in radians.
#[wasm_bindgen]
pub struct LonLat {
  #[wasm_bindgen(readonly)]
  pub lon: f64,
  #[wasm_bindgen(readonly)]
  pub lat: f64,
}

impl LonLat {
  pub fn from_vals(lon: f64, lat: f64) -> LonLat {
    LonLat { lon, lat}
  }

  pub fn from_tuple(lonlat: (f64, f64)) -> LonLat {
    LonLat::from_vals(lonlat.0, lonlat.1)
  }
}*/

/// Compute the position on the unit sphere of the center (in the Euclidean projection plane)
/// of the cell associated to the given cell number.
/// # Input
/// - `order`: the order of the cell
/// - `icell`: the cell number value of the cell we look for the unprojected center, in the NESTED scheme
/// # Output
/// - `[lon, lat]` in degrees, the unprojected position (on the unit sphere) of the center of 
///   the cell in the Euclidean plane
///   - `lon`, longitude in `[0, 360]` degrees;
///   - `lat`, latitude in `[-90, 90]` degrees. 
/// 
/// # Panics
/// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
#[wasm_bindgen(js_name = nestedCenter)]
pub fn wasm_center(depth: u8, hash: f64) -> Box<[f64]> {
  let (lon, lat) = center(depth, hash as u64);
  Box::new([lon.to_degrees(), lat.to_degrees()])
}

/*
/// This object stores the 4 vertices of an HEALPix cell.
/// There is one vertex per cardinal point, i.e. at South, East, North and West.
#[wasm_bindgen]
pub struct CellVertices {
  lons: [f64; 4],
  lats: [f64; 4],
}

#[wasm_bindgen]
impl CellVertices {
  /// Returns the longitude of the south vertex, in radians
  #[wasm_bindgen(js_name = getSouthLon)]
  pub fn get_south_lon(&self) -> f64 {
    self.lons[0]
  }
  /// Returns the latitude of the south vertex, in radians
  #[wasm_bindgen(js_name = getSouthLat)]
  pub fn get_south_lat(&self) -> f64 {
    self.lats[0]
  }
  /// Returns the longitude of the east vertex, in radians
  #[wasm_bindgen(js_name = getEastLon)]
  pub fn get_east_lon(&self) -> f64 {
    self.lons[1]
  }
  /// Returns the latitude of the east vertex, in radians
  #[wasm_bindgen(js_name = getEastLat)]
  pub fn get_east_lat(&self) -> f64 {
    self.lats[1]
  }
  /// Returns the longitude of the north vertex, in radians
  #[wasm_bindgen(js_name = getNorthLon)]
  pub fn get_north_lon(&self) -> f64 {
    self.lons[2]
  }
  /// Returns the latitude of the north vertex, in radians
  #[wasm_bindgen(js_name = getNorthLat)]
  pub fn get_north_lat(&self) -> f64 {
    self.lats[2]
  }
  /// Returns the longitude of the west vertex, in radians
  #[wasm_bindgen(js_name = getWestLon)]
  pub fn get_west_lon(&self) -> f64 {
    self.lons[3]
  }
  /// Returns the latitude of the west vertex, in radians
  #[wasm_bindgen(js_name = getWestLat)]
  pub fn get_west_lat(&self) -> f64 {
    self.lats[3]
  }
}


/// Computes the location on the unit sphere of the 4 vertices of the given HEALPix cell
/// (define by its depth and number).
/// # Inputs
/// - `depth` the depth of the cell we look for the vertices
/// - `hash` the nested number of the cell we look for the vertices
#[wasm_bindgen(js_name = nestedVertices)]
pub fn wasm_vertices(depth: u8, hash: f64) -> CellVertices  {
  let [(s_lon, s_lat), (e_lon, e_lat), (n_lon, n_lat), (w_lon, w_lat)] = vertices(depth, hash as u64);
  CellVertices {
    lons: [s_lon, e_lon, n_lon, w_lon],
    lats: [s_lat, e_lat, n_lat, w_lat],
  }
}

/// Computes the location on the unit sphere of the 4 vertices of the given HEALPix cell
/// (define by its depth and number).
/// The coordinates are provided in radians and are stored in the returned array as following:
/// `[south_lon, south_lat, east_lon, east_lat, north_lon, north_lat, west_lon, west_lat]`
/// # Inputs
/// - `depth` the depth of the cell we look for the vertices
/// - `hash` the nested number of the cell we look for the vertices
#[wasm_bindgen(js_name = hpxVerticesFlat)]
pub fn wasm_vertices_flat(depth: u8, hash: f64) -> Box<[f64]>  {
  let [(s_lon, s_lat), (e_lon, e_lat), (n_lon, n_lat), (w_lon, w_lat)] = vertices(depth, hash as u64);
  let flat: [f64; 8] = [s_lon, s_lat, e_lon, e_lat, n_lon, n_lat, w_lon, w_lat];
  // let flat: [[f64; 2]; 4] = [[s_lon, s_lat], [e_lon, e_lat], [n_lon, n_lat], [w_lon, w_lat]];
  Box::new(flat)
}

*/

/// Computes the location on the unit sphere of the 4 vertices of the given HEALPix cell
/// (define by its depth and number).
/// # Inputs
/// - `order` the order of the cell we look for the vertices
/// - `icell`: the cell number value of the cell we look for the unprojected center, in the NESTED scheme
/// # Output
/// - array containing the longitudes and latitudes (in degrees) of the vertices in the following order:
///   `[SouthLon, SouthLat, EastLon, EastLat, NoethLon, NorthLat, WestLon, WestLat]`
#[wasm_bindgen(js_name = nestedVertices)]
pub fn wasm_vertices(depth: u8, hash: f64) -> Box<[f64]>  {
  let [(s_lon, s_lat), (e_lon, e_lat), (n_lon, n_lat), (w_lon, w_lat)] = vertices(depth, hash as u64);
  Box::new([
    s_lon.to_degrees(), s_lat.to_degrees(),
    e_lon.to_degrees(), e_lat.to_degrees(),
    n_lon.to_degrees(), n_lat.to_degrees(),
    w_lon.to_degrees(), w_lat.to_degrees(),
  ])
}

/*#[wasm_bindgen] CA NE MARCHE PAS!!
pub fn wasm_vertices_struct(depth: u8, hash: f64) -> Box<[LonLat]>  {
  let [s, e, n, w] = vertices(depth, hash as u64);
  Box::new([
    LonLat::from_tuple(s), 
    LonLat::from_tuple(e), 
    LonLat::from_tuple(n), 
    LonLat::from_tuple(w)
  ])
}*/


// Neighbours

/// Contains the list of the neighbours of the `center` HEALPix cell.
/// In general, a cell has 8 neighbours.
/// But the top and bottom corners of equatorial base cells do not have 
/// north and south neighbours respectively.
/// The left most and right most corners of the polar caps base cells do not have
/// west and east corners respectively.
#[wasm_bindgen]
pub struct Neighbours {
  #[wasm_bindgen(readonly)]
  pub south: Option<f64>,
  #[wasm_bindgen(readonly)]
  pub southeast: f64,
  #[wasm_bindgen(readonly)]
  pub east: Option<f64>,
  #[wasm_bindgen(readonly)]
  pub southwest: f64,
  #[wasm_bindgen(readonly)]
  pub center: f64,
  #[wasm_bindgen(readonly)]
  pub norhteast: f64,
  #[wasm_bindgen(readonly)]
  pub west: Option<f64>,
  #[wasm_bindgen(readonly)]
  pub norhtwest: f64,
  #[wasm_bindgen(readonly)]
  pub north: Option<f64>,
}

impl Neighbours {
  fn from_map(neig_map: MainWindMap<u64>) -> Neighbours {
    Neighbours {
      south: neig_map.get(MainWind::S).map(|&h| h as f64),
      southeast: *neig_map.get(MainWind::SE).unwrap() as f64,
      east: neig_map.get(MainWind::E).map(|&h| h as f64),
      southwest: *neig_map.get(MainWind::SW).unwrap() as f64,
      center: *neig_map.get(MainWind::C).unwrap() as f64,
      norhteast: *neig_map.get(MainWind::NE).unwrap() as f64,
      west: neig_map.get(MainWind::W).map(|&h| h as f64),
      norhtwest: *neig_map.get(MainWind::NW).unwrap() as f64,
      north: neig_map.get(MainWind::N).map(|&h| h as f64),
    }
  }
}

/// Returns the hash values of all the neighbour cells of the cell of given cell number.
/// The given cell itself is included.
/// # Inputs
/// - `order` the order of the cell we look for the neighbours
/// - `icell` the nested number of the cell we look for the neighbours
#[wasm_bindgen(js_name = nestedNeighbours)]
pub fn wasm_neighbours(depth: u8, hash: f64) -> Neighbours {
  Neighbours::from_map(neighbours(depth, hash as u64, true))
}


// Bmoc cell

/// Represents a BMOC object resulting from a cone search or a polygon query. 
/// For practical reasons of data exchange between Javascript and WebAssembly, each elements
/// of a cell (depth, hash, flag) are stored in separate arrays.
/// (I firs wanted to return an array of Cell, but so far it is not possible with wasm-bindgen)
#[wasm_bindgen(js_name = BMOC)]
pub struct ColBMOC {
  #[wasm_bindgen(readonly)]
  pub n_cells: u32,
  #[wasm_bindgen(readonly)]
  pub depth_max: u8,
  depths:  Box<[u8]>,
  hashs:  Box<[f64]>,
  is_full_flags:  Box<[bool]>,
}

/// Represents a BMOC cell: its depth, number and flag telling if the cell is fully/partially covered.
#[wasm_bindgen]
pub struct Cell {
  /// The order of the HEALPix cell
  #[wasm_bindgen(readonly)]
  pub order: u8,
  /// The nested HEALPix cell number 
  #[wasm_bindgen(readonly)]
  pub icell: f64,
  /// The flag telling if the cell if fully or partially covered
  #[wasm_bindgen(readonly)]
  pub isfull: bool,
}

#[wasm_bindgen(js_class = BMOC)]
impl ColBMOC {
  /// Returns the number of cells (of various depth) in the BMOC
  #[wasm_bindgen(js_name = getSize)]
  pub fn len(&self) -> u32 {
    self.n_cells
  }
  /// Returns the maximal depth of the BMOC
  #[wasm_bindgen(js_name = getDepth)]
  pub fn get_moc_depth(&self) -> u8 {
    self.depth_max
  }
  /// Utility method replacing the calls to `getCellDepth`, `getCellHash` and `getCellFlag`
  /// # Inputs
  /// - `i` the index of the cell in the BMOC
  #[wasm_bindgen(js_name = getCell)]
  pub fn get_cell(&self, icell: u32) -> Cell {
    let icell = icell as usize;
    Cell {
      order: self.depths[icell],
      icell: self.hashs[icell],
      isfull: self.is_full_flags[icell],
    }
  }
  
  /// Returns the depth of cell number `i` in the BMOC
  /// # Inputs
  /// - `i` the index of the cell in the BMOC
  #[wasm_bindgen(js_name = getCellDepth)]
  pub fn get_cell_depth(&self, icell: u32) -> u8 {
    self.depths[icell as usize]
  }
  /// Returns the hash value (or pixel number) of the cell number `i` in the BMOC
  /// # Inputs
  /// - `i` the index of the cell in the BMOC
  #[wasm_bindgen(js_name = getCellHash)]
  pub fn get_cell_hash(&self, icell: u32) -> f64 {
    self.hashs[icell as usize]
  }
  /// Returns the status flag of the cell number `i` in the BMOC 
  /// (`true`: cell fully covered; `false`: cell partially covered)
  /// # Inputs
  /// - `i` the index of the cell in the BMOC
  #[wasm_bindgen(js_name = getCellFlag)]
  pub fn get_cell_flag(&self, icell: u32) -> bool {
    self.is_full_flags[icell as usize]
  }
}

/*impl BMocCell {
  fn new(raw_value: u64, depth_max: u8) -> BMocCell {
    // Extract the flag
    let is_full = (raw_value | 1u64) as u8;
    // Remove the flag bit, then divide by 2 (2 bits per level)
    let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
    // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
    let hash = raw_value >> (2 + (delta_depth << 1));
    let depth = depth_max - delta_depth;
    BMocCell { raw_value: raw_value as f64, depth, hash: hash as f64, is_full }
  }
}*/

// Cone

/// Returns the BMOC covered by the given cone.
/// # Inputs
/// - `order`: the maximum order of the cell in the returned BMOC
/// - `cone_lon`: cone center longitude, in degrees
/// - `cone_lat`: cone center latitude in degrees, must be in `[-90, 90]`
/// - `cone_radius`: radius of the cone, in degrees
/// # Outputs
/// - the BMOC covered by the given cone
#[wasm_bindgen(js_name = nestedQueryConeBMOC)]
pub fn wasm_cone_coverage(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> ColBMOC {
  bmoc2colbmoc(cone_coverage_approx(depth, cone_lon.to_radians(), cone_lat.to_radians(), cone_radius.to_radians()))
}

/// Returns the flat list of HEALPix cell of given order covered by the given cone.
/// # Inputs
/// - `order`: the order of the returned cells
/// - `cone_lon`: cone center longitude, in degrees
/// - `cone_lat`: cone center latitude in degrees, must be in `[-90, 90]`
/// - `cone_radius`: radius of the cone, in degrees
/// # Outputs
/// - the flat list of HEALPix cell of given order covered by the given cone
#[wasm_bindgen(js_name = nestedQueryCone)]
pub fn wasm_cone_coverage_flat(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> Box<[f64]>  {
  let a: Vec<f64> = cone_coverage_approx_flat(depth, cone_lon.to_radians(), cone_lat.to_radians(), cone_radius.to_radians())
    .iter().map(|&h| h as f64).collect();
  a.into_boxed_slice()
}


/// Returns the BMOC covered by the given polygon.
/// # Inputs
/// - `order`: the maximum order of the cell in the returned BMOC
/// - `vertices_coos`: the flat list of polygon vertices coordinates, in radians
///    `[v1.lon, v1.lat, v2.lon, v2.lat, ..., vN.lon, vN.lat]`
/// # Outputs
/// - the BMOC covered by the given polygon
#[wasm_bindgen(js_name = nestedQueryPolygonBMOC)]
pub fn wasm_polygon_coverage(depth: u8, vertices_coos: Box<[f64]>) -> ColBMOC {
  let mut vertices:Vec<(f64, f64)> = Vec::with_capacity(vertices_coos.len() >> 1);
  for i in (0..vertices_coos.len()).step_by(2) {
    vertices.push((vertices_coos[i], vertices_coos[i+1]));
  }
  bmoc2colbmoc(polygon_coverage(depth, &vertices.into_boxed_slice(), true))
}

/// Returns the BMOC covered by the given polygon.
/// # Inputs
/// - `order`: the order of the returned cells
/// - `vertices_coos`: the flat list of polygon vertices coordinates, in radians
///    `[v1.lon, v1.lat, v2.lon, v2.lat, ..., vN.lon, vN.lat]`
/// # Outputs
/// - the BMOC covered by the given polygon
#[wasm_bindgen(js_name = nestedQueryPolygon)]
pub fn wasm_polygon_coverage_flat(depth: u8, vertices_coos: Box<[f64]>) -> Box<[f64]> {
  let mut vertices:Vec<(f64, f64)> = Vec::with_capacity(vertices_coos.len() >> 1);
  for i in (0..vertices_coos.len()).step_by(2) {
    vertices.push((vertices_coos[i], vertices_coos[i+1]));
  }
  let bmoc = polygon_coverage(depth, &vertices.into_boxed_slice(), true);
  let flat_iter: BMOCFlatIter = bmoc.flat_iter();
  let mut flat: Vec<f64> = Vec::with_capacity(flat_iter.deep_size());
  for c in flat_iter {
    flat.push(c as f64);
  }
  flat.into_boxed_slice()
}

fn bmoc2colbmoc(bmoc: BMOC) -> ColBMOC {
  let n_elems = bmoc.entries.len();
  let mut depths: Vec<u8> = Vec::with_capacity(n_elems);
  let mut hashs: Vec<f64> = Vec::with_capacity(n_elems);
  let mut flags: Vec<bool> = Vec::with_capacity(n_elems);
  for raw_value in bmoc.entries.iter() {
    let cell = bmoc.from_raw_value(*raw_value);
    depths.push(cell.depth);
    hashs.push(cell.hash as f64);
    flags.push(cell.is_full);
  }
  ColBMOC {
    n_cells: n_elems as u32,
    depth_max: bmoc.get_depth_max(),
    depths: depths.into_boxed_slice(),
    hashs: hashs.into_boxed_slice(),
    is_full_flags: flags.into_boxed_slice(),
  }
}

/*
Found here: https://users.rust-lang.org/t/flattening-a-vector-of-tuples/11409/2
fn flatten_coords(coords: &[(i32, i32)]) -> Vec<i32> {
  let mut result = data.to_vec();
  unsafe {
    result.set_len(data.len() * 2);
    std::mem::transmute(result)
  }
}*/
