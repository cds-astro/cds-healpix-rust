use std::{
  cmp::Ordering,
  collections::{
    btree_map::{IntoIter as BTreeMapIntoIter, Iter as BTreeMapIter, Values},
    BTreeMap,
  },
  f64::consts::PI,
  fs::File,
  io::{stdin, BufRead, BufReader, BufWriter, Error as IoError, Write},
  iter::{Map, Peekable},
  ops::{Add, AddAssign},
  path::Path,
};

use rayon::{
  iter::{IntoParallelRefIterator, ParallelIterator},
  slice::ParallelSlice,
  ThreadPool,
};

use crate::nested::map::skymap::DegradableBySumming;
use crate::{
  n_hash,
  nested::{
    get,
    map::{
      fits::{error::FitsError, write::write_explicit_skymap_fits},
      mom::{
        impls::zvec::MomVecImpl, new_chi2_density_ref_merger_no_depth_threshold,
        new_chi2_density_ref_merger_with_depth_threshold, Mom,
      },
      skymap::{
        implicit::{
          ImplicitCountMap, ImplicitCountMapU32, ImplicitDensityMap, ImplicitSkyMapArray,
        },
        SkyMap, SkyMapValue,
      },
      HHash,
    },
    sort::get_hpx_opt,
  },
};

/// `BTreeMap` based structure: slow loading (tree creation), but fast queries.
/// To read from FITS file, an alternative is a sorted array with binary search.
#[derive(Debug)]
pub struct ExplicitSkyMapBTree<H: HHash, V: SkyMapValue> {
  depth: u8,
  entries: BTreeMap<H, V>,
  _zero: V,
}
impl<H: HHash, V: SkyMapValue> ExplicitSkyMapBTree<H, V> {
  pub fn new(depth: u8, entries: BTreeMap<H, V>) -> Self {
    Self {
      depth,
      entries,
      _zero: V::zero(),
    }
  }
  pub fn into_implicit_map(self, null_value: V) -> ImplicitSkyMapArray<H, V> {
    // We may implement an "Implicit iterator" to write a very large map on the disk
    // (with an implicit map possibly larger that the RAM).
    let depth = self.depth;
    let mut values = vec![null_value; n_hash(self.depth) as usize];
    for (k, v) in self.entries {
      values[k.to_usize().unwrap()] = v;
    }
    ImplicitSkyMapArray::new(depth, values.into_boxed_slice())
  }
}
impl<V: SkyMapValue> ExplicitSkyMapBTree<u32, V> {
  pub fn into_u64_hash(self) -> ExplicitSkyMapBTree<u64, V> {
    ExplicitSkyMapBTree::new(
      self.depth,
      self
        .entries
        .into_iter()
        .map(|(k, v)| (k as u64, v))
        .collect(),
    )
  }
}

impl<H: HHash, V: SkyMapValue> Add for ExplicitSkyMapBTree<H, V> {
  type Output = Self;

  fn add(mut self, rhs: Self) -> Self::Output {
    let l = self.entries.into_iter().peekable();
    let r = rhs.entries.into_iter().peekable();
    struct SortedIt<H: HHash, V: SkyMapValue> {
      l: Peekable<BTreeMapIntoIter<H, V>>,
      r: Peekable<BTreeMapIntoIter<H, V>>,
    }
    impl<H: HHash, V: SkyMapValue> Iterator for SortedIt<H, V> {
      type Item = (H, V);
      fn next(&mut self) -> Option<Self::Item> {
        match (self.l.peek(), self.r.peek()) {
          (Some((hl, _)), Some((hr, _))) => match Ord::cmp(hl, hr) {
            Ordering::Equal => {
              let (hl, vl) = self.l.next().unwrap();
              let (_, vr) = self.r.next().unwrap();
              Some((hl, vl + vr))
            }
            Ordering::Less => self.l.next(),
            Ordering::Greater => self.r.next(),
          },
          (Some(_), None) => self.l.next(),
          (None, Some(_)) => self.r.next(),
          (None, None) => None,
        }
      }
      fn size_hint(&self) -> (usize, Option<usize>) {
        let (min_l, max_l) = self.l.size_hint();
        let (min_r, max_r) = self.r.size_hint();
        (
          usize::min(min_l, min_r),
          if let (Some(max_l), Some(max_r)) = (max_l, max_r) {
            Some(max_l + max_r) // check larger the usize to returnn None...
          } else {
            None
          },
        )
      }
    }
    self.entries = SortedIt { l, r }.collect();
    self
  }
}
impl<H: HHash, V: SkyMapValue + Send + Sync + AddAssign> ExplicitSkyMapBTree<H, V> {
  pub fn par_add(mut self, rhs: Self) -> Self {
    // No par algo so far... TODO: find one!
    self = self + rhs;
    self
  }
}
impl<'a, H: HHash, V: SkyMapValue + 'a> SkyMap<'a> for ExplicitSkyMapBTree<H, V> {
  type HashType = H;
  type ValueType = V;
  type ValuesIt = Values<'a, H, V>;
  type EntriesIt = Map<BTreeMapIter<'a, H, V>, fn((&'a H, &'a V)) -> (H, &'a V)>;
  type OwnedEntriesIt = BTreeMapIntoIter<H, V>;

  fn depth(&self) -> u8 {
    self.depth
  }

  fn is_implicit(&self) -> bool {
    false
  }

  fn len(&self) -> usize {
    self.entries.len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    match self.entries.get(&hash) {
      Some(v) => v,
      None => &self._zero,
    }
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.entries.values()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self.entries.iter().map(|(h, v)| (*h, v))
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    BTreeMap::into_iter(self.entries)
  }
}

impl<H: HHash, V: SkyMapValue + Send + Sync + AddAssign> DegradableBySumming
  for ExplicitSkyMapBTree<H, V>
{
  type Degraded = Self;

  fn degrade_sum(self, new_depth: u8) -> Self::Degraded {
    if new_depth >= self.depth {
      self
    } else {
      let twice_dd = ((self.depth - new_depth) << 1) as usize;
      let mut entries = BTreeMap::new();
      let mut it = self.owned_entries();
      if let Some((k, v)) = it.next() {
        let mut prev_k = k >> twice_dd;
        let mut prev_v = v;
        for (k, v) in it {
          let curr_k = k >> twice_dd;
          if curr_k == prev_k {
            prev_v += v;
          } else {
            entries.insert(prev_k, prev_v);
            prev_k = curr_k;
            prev_v = v;
          }
        }
        entries.insert(prev_k, prev_v);
      }
      Self::new(new_depth, entries)
    }
  }
}

/// SkyMap implementation use to store counts, rely on u32 HEALPix index instead of u64.
#[derive(Debug)]
pub struct ExplicitCountMap<H: HHash>(ExplicitSkyMapBTree<H, u32>);
impl<H: HHash> From<ExplicitSkyMapBTree<H, u32>> for ExplicitCountMap<H> {
  fn from(value: ExplicitSkyMapBTree<H, u32>) -> Self {
    Self(value)
  }
}
impl<H: HHash> ExplicitCountMap<H> {
  pub fn as_explicit_skymap_btree(&self) -> &ExplicitSkyMapBTree<H, u32> {
    &self.0
  }
  pub fn into_explicit_skymap_btree(self) -> ExplicitSkyMapBTree<H, u32> {
    self.0
  }
  /// Build a count skymap from an iterator over HEALPix cells at the given depth.
  pub fn from_hash_values<I>(depth: u8, hash_values_it: I) -> Self
  where
    I: Iterator<Item = H>,
  {
    let mut btree = BTreeMap::new();
    for h in hash_values_it {
      // Most maps will have a lot of values in a same cell, so we don't mind a slow first insert
      // and favor fast updates.
      match btree.get_mut(&h) {
        Some(c) => *c += 1,
        None => assert!(btree.insert(h, 1).is_none()),
      }
    }
    Self(ExplicitSkyMapBTree::new(depth, btree))
  }

  // From Sorted Hash Values?!!

  /// Build a count skymap from an iterator over position ( (ra, dec), in radian).
  pub fn from_positions<I>(depth: u8, pos_it_rad: I) -> Self
  where
    I: Iterator<Item = (f64, f64)>,
  {
    let layer = get(depth);
    Self::from_hash_values(
      depth,
      pos_it_rad.map(move |(l, b)| H::from_u64(layer.hash(l, b))),
    )
  }

  pub fn from_csv_stdin_par(
    ilon: usize,
    ilat: usize,
    separator: Option<char>,
    has_header: bool,
    depth: u8,
    chunk_size: usize,
    thread_pool: &ThreadPool,
  ) -> Result<Self, IoError> {
    Self::from_csv_gen_par(
      BufReader::new(stdin()),
      ilon,
      ilat,
      separator,
      has_header,
      depth,
      chunk_size,
      thread_pool,
    )
  }

  pub fn from_csv_file_par<P: AsRef<Path>>(
    path: P,
    ilon: usize,
    ilat: usize,
    separator: Option<char>,
    has_header: bool,
    depth: u8,
    chunk_size: usize,
    thread_pool: &ThreadPool,
  ) -> Result<Self, IoError> {
    Self::from_csv_gen_par(
      BufReader::new(File::open(&path)?),
      ilon,
      ilat,
      separator,
      has_header,
      depth,
      chunk_size,
      thread_pool,
    )
  }

  pub fn from_csv_gen_par<R: BufRead + Send>(
    bufr: R,
    ilon: usize,
    ilat: usize,
    separator: Option<char>,
    has_header: bool,
    depth: u8,
    chunk_size: usize,
    thread_pool: &ThreadPool,
  ) -> Result<Self, IoError> {
    let mut it = bufr.lines().peekable();
    // Handle starting comments
    while let Some(Ok(_)) = it.next_if(|res| {
      res
        .as_ref()
        .map(|line| line.starts_with('#'))
        .unwrap_or(false)
    }) {}
    // Handle header line
    if has_header {
      it.next().transpose()?;
    }
    // Ok, go!
    Self::from_csv_it_par(it, ilon, ilat, separator, depth, chunk_size, thread_pool)
  }

  /// Compute the count map from an iterator iterating on raw CSV rows.
  /// # Params
  /// * `it` the iterator an rows
  /// * `ilon` index of the column containing the longitude, starting at 0
  /// * `ilat` index of the column containing the latitude, starting at 0
  /// * `separator` file separator (',' if None)
  /// * `depth` HEALPix depth (or order) of the map
  /// * `chunk_size` number of rows to be processed in parallel (the memory will hold twice this number
  ///   since a chunk is read while another chunk is processed)
  /// * `thread_pool` the thread pool in which the process will be executed
  /// # Note
  /// The parallel processing chosen resort on a single sequential read.
  /// Although it is not the fastest option with SSDs, it should ensure reasonably good performances
  /// with both SSDs and HDDs.
  /// # Warning
  /// Ensure that you have already removed the comment and the possible header line,
  /// E.g. using:
  /// ```rust,ignore
  ///     let mut it = it.peekable();
  ///     // Handle starting comments
  ///     while let Some(Ok(line)) = it.next_if(|res| {
  ///       res
  ///         .as_ref()
  ///         .map(|line| line.starts_with('#'))
  ///         .unwrap_or(false)
  ///     }) { }
  ///     // Handle header line
  ///     if has_header {
  ///       it.next().transpose()?;
  ///     }
  /// ```
  pub fn from_csv_it_par<I>(
    mut it: I,
    ilon: usize,
    ilat: usize,
    separator: Option<char>,
    depth: u8,
    chunk_size: usize,
    thread_pool: &ThreadPool,
  ) -> Result<Self, IoError>
  where
    I: Iterator<Item = Result<String, IoError>> + Send,
  {
    let separator = separator.unwrap_or(',');
    let layer = get(depth);
    let n_thread = thread_pool.current_num_threads();
    let hpx = get_hpx_opt(ilon, ilat, separator, layer);
    let count_map_fn = |chunk: Vec<String>| {
      chunk
        .par_chunks((chunk.len() / (n_thread << 2)).max(10_000))
        .map(|elems| {
          ExplicitCountMap::<H>::from_hash_values(
            depth,
            elems.iter().filter_map(&hpx).map(H::from_u64),
          )
        })
        .reduce_with(|mapl, mapr| mapl.par_add(mapr))
    };
    fn load_n<I: Iterator<Item = Result<String, IoError>>>(
      chunk_size: usize,
      it: &mut I,
    ) -> Result<Vec<String>, IoError> {
      it.take(chunk_size).collect()
    }

    let mut chunk = load_n(chunk_size, it.by_ref())?;
    let mut count_map = Self(ExplicitSkyMapBTree::new(depth, BTreeMap::new()));
    while !chunk.is_empty() {
      let (next_chunk, new_count_map) = thread_pool.join(
        || load_n(chunk_size, it.by_ref()),
        || {
          if let Some(local_count_map) = count_map_fn(chunk) {
            count_map.par_add(local_count_map)
          } else {
            count_map
          }
        },
      );
      chunk = next_chunk?;
      count_map = new_count_map;
    }
    Ok(count_map)
  }

  pub fn to_dens_map(&self) -> ExplicitDensityMap<H> {
    let depth = self.depth();
    let one_over_area = (n_hash(depth) >> 2) as f64 / PI;
    ExplicitDensityMap(ExplicitSkyMapBTree::<H, f64>::new(
      depth,
      self
        .0
        .entries
        .iter()
        .map(|(h, count)| (h.clone(), *count as f64 * one_over_area))
        .collect(),
    ))
  }

  /// For custom thread pool, use `thread_pool.install(|| to_dens_map_par())`
  pub fn to_dens_map_par(&self) -> ExplicitDensityMap<H> {
    let depth = self.depth();
    let one_over_area = (n_hash(depth) >> 2) as f64 / PI;
    ExplicitDensityMap(ExplicitSkyMapBTree::<H, f64>::new(
      depth,
      self
        .0
        .entries
        .par_iter()
        .map(|(h, count)| (h.clone(), *count as f64 * one_over_area))
        .collect(),
    ))
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    write_explicit_skymap_fits(writer, self.as_explicit_skymap_btree())
  }

  pub fn to_fits_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits(BufWriter::new(file)))
  }
}
impl ExplicitCountMap<u32> {
  pub fn into_implicit_skymap(self, null_value: u32) -> ImplicitCountMapU32 {
    self.0.into_implicit_map(null_value).into()
  }
}
impl ExplicitCountMap<u64> {
  pub fn into_implicit_skymap(self, null_value: u32) -> ImplicitCountMap {
    self.0.into_implicit_map(null_value).into()
  }
}

impl<H: HHash> Add for ExplicitCountMap<H> {
  type Output = Self;

  fn add(mut self, rhs: Self) -> Self::Output {
    /*let l = self.0.entries.into_iter().peekable();
    let r = rhs.0.entries.into_iter().peekable();
    struct SortedIt<H: HHash> {
      l: Peekable<BTreeMapIntoIter<H, u32>>,
      r: Peekable<BTreeMapIntoIter<H, u32>>,
    }
    impl<H: HHash> Iterator for SortedIt<H> {
      type Item = (H, u32);
      fn next(&mut self) -> Option<Self::Item> {
        match (self.l.peek(), self.r.peek()) {
          (Some((hl, _)), Some((hr, _))) => match Ord::cmp(hl, hr) {
            Ordering::Equal => {
              let (hl, vl) = self.l.next().unwrap();
              let (_, vr) = self.r.next().unwrap();
              Some((hl, vl + vr))
            }
            Ordering::Less => self.l.next(),
            Ordering::Greater => self.r.next(),
          },
          (Some(_), None) => self.l.next(),
          (None, Some(_)) => self.r.next(),
          (None, None) => None,
        }
      }
      fn size_hint(&self) -> (usize, Option<usize>) {
        let (min_l, max_l) = self.l.size_hint();
        let (min_r, max_r) = self.r.size_hint();
        (
          usize::min(min_l, min_r),
          if let (Some(max_l), Some(max_r)) = (max_l, max_r) {
            Some(max_l + max_r) // check larger the usize to returnn None...
          } else {
            None
          },
        )
      }
    }
    self.0.entries = SortedIt { l, r }.collect();*/
    self.0 = self.0 + rhs.0;
    self
  }
}
impl<H: HHash> ExplicitCountMap<H> {
  pub fn par_add(mut self, rhs: Self) -> Self {
    self.0 = self.0.par_add(rhs.0);
    self
  }
}
impl<'a, H: HHash> SkyMap<'a> for ExplicitCountMap<H> {
  type HashType = H;
  type ValueType = u32;
  type ValuesIt = Values<'a, H, u32>;
  type EntriesIt = Map<BTreeMapIter<'a, H, u32>, fn((&'a H, &'a u32)) -> (H, &'a u32)>;
  type OwnedEntriesIt = BTreeMapIntoIter<H, u32>;

  fn depth(&self) -> u8 {
    self.as_explicit_skymap_btree().depth()
  }

  fn is_implicit(&self) -> bool {
    false
  }

  fn len(&self) -> usize {
    self.as_explicit_skymap_btree().len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    self.as_explicit_skymap_btree().get(hash)
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.as_explicit_skymap_btree().values()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self.as_explicit_skymap_btree().entries()
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    self.into_explicit_skymap_btree().owned_entries()
  }
}

/// SkyMap implementation use to store densities, rely on HEALPix index instead of u64.
#[derive(Debug)]
pub struct ExplicitDensityMap<H: HHash>(ExplicitSkyMapBTree<H, f64>);
impl<H: HHash> From<ExplicitSkyMapBTree<H, f64>> for ExplicitDensityMap<H> {
  fn from(value: ExplicitSkyMapBTree<H, f64>) -> Self {
    Self(value)
  }
}
impl<H: HHash> ExplicitDensityMap<H> {
  pub fn as_explicit_skymap_btree(&self) -> &ExplicitSkyMapBTree<H, f64> {
    &self.0
  }
  pub fn into_explicit_skymap_btree(self) -> ExplicitSkyMapBTree<H, f64> {
    self.0
  }
  /// Build a density skymap from an iterator over position ( (ra, dec), in radian).
  /// # Panics
  /// * if `depth > 12`.
  pub fn from_positions<I>(depth: u8, pos_it_rad: I) -> Self
  where
    I: Iterator<Item = (f64, f64)>,
  {
    let layer = get(depth);
    let one_over_cell_area = layer.n_hash as f64 / (4.0 * PI);
    let mut btree = BTreeMap::new();
    for h in pos_it_rad.map(move |(l, b)| H::from_u64(layer.hash(l, b))) {
      // Most maps will have a lot of values in a same cell, so we don't mind a slow first insert
      // and favor fast updates.
      match btree.get_mut(&h) {
        Some(c) => *c += one_over_cell_area,
        None => assert!(btree.insert(h, one_over_cell_area).is_none()),
      }
    }
    Self(ExplicitSkyMapBTree::new(depth, btree))
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    write_explicit_skymap_fits(writer, self.as_explicit_skymap_btree())
  }

  pub fn to_fits_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits(BufWriter::new(file)))
  }

  // to_png
}

impl ExplicitDensityMap<u32> {
  /// # Params
  /// * `chi2_of_3dof_threshold`: threshold on the value of the chi square distribution with 3
  /// degrees of freedom below which we consider the 4 values of 4 sibling cells as coming
  /// from the same normal distribution which mean and variance comes from a poisson distribution.
  /// Here a few typical values corresponding the the given completeness:
  ///     + Completeness = 90.0% =>  6.251
  ///     + Completeness = 95.0% =>  7.815
  ///     + Completeness = 97.5% =>  9.348
  ///     + Completeness = 99.0% => 11.345
  ///     + Completeness = 99.9% => 16.266
  /// * `depth_threshold`: threshold on `depth` to avoid making to low resolution cells, i.e MOM minimum depth
  pub fn to_chi2_mom(
    &self,
    chi2_of_3dof_threshold: f64,
    depth_threshold: Option<u8>,
  ) -> MomVecImpl<u32, f64> {
    // WARNING: result will bve different from to_chi2_mom on implicit map because
    // no value is different from value = 0.
    match depth_threshold {
      None => MomVecImpl::from_skymap_ref(
        &self.0,
        new_chi2_density_ref_merger_no_depth_threshold(chi2_of_3dof_threshold),
      ),
      Some(depth_threshold) => MomVecImpl::from_skymap_ref(
        &self.0,
        new_chi2_density_ref_merger_with_depth_threshold(chi2_of_3dof_threshold, depth_threshold),
      ),
    }
  }
}
impl ExplicitDensityMap<u64> {
  /// # Params
  /// * `chi2_of_3dof_threshold`: threshold on the value of the chi square distribution with 3
  /// degrees of freedom below which we consider the 4 values of 4 sibling cells as coming
  /// from the same normal distribution which mean and variance comes from a poisson distribution.
  /// Here a few typical values corresponding the the given completeness:
  ///     + Completeness = 90.0% =>  6.251
  ///     + Completeness = 95.0% =>  7.815
  ///     + Completeness = 97.5% =>  9.348
  ///     + Completeness = 99.0% => 11.345
  ///     + Completeness = 99.9% => 16.266
  /// * `depth_threshold`: threshold on `depth` to avoid making to low resolution cells
  pub fn to_chi2_mom(
    &self,
    chi2_of_3dof_threshold: f64,
    depth_threshold: Option<u8>,
  ) -> MomVecImpl<u64, f64> {
    // WARNING: result will bve different from to_chi2_mom on implicit map because
    // no value is different from value = 0.
    match depth_threshold {
      None => MomVecImpl::from_skymap_ref(
        &self.0,
        new_chi2_density_ref_merger_no_depth_threshold(chi2_of_3dof_threshold),
      ),
      Some(depth_threshold) => MomVecImpl::from_skymap_ref(
        &self.0,
        new_chi2_density_ref_merger_with_depth_threshold(chi2_of_3dof_threshold, depth_threshold),
      ),
    }
  }
}
impl ExplicitDensityMap<u32> {
  pub fn into_implicit_skymap(self, null_value: f64) -> ImplicitDensityMap {
    self.0.into_implicit_map(null_value).into()
  }
}
/*impl ExplicitDensityMap<u64> {
  pub fn into_implicit_skymap(self) -> ImplicitDensityMap {
    self.0.into_implicit_map().into()
  }
}*/

impl<'a, H: HHash> SkyMap<'a> for ExplicitDensityMap<H> {
  type HashType = H;
  type ValueType = f64;
  type ValuesIt = Values<'a, H, f64>;
  type EntriesIt = Map<BTreeMapIter<'a, H, f64>, fn((&'a H, &'a f64)) -> (H, &'a f64)>;
  type OwnedEntriesIt = BTreeMapIntoIter<H, f64>;

  fn depth(&self) -> u8 {
    self.as_explicit_skymap_btree().depth()
  }

  fn is_implicit(&self) -> bool {
    false
  }

  fn len(&self) -> usize {
    self.as_explicit_skymap_btree().len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    self.as_explicit_skymap_btree().get(hash)
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.as_explicit_skymap_btree().values()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self.as_explicit_skymap_btree().entries()
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    self.into_explicit_skymap_btree().owned_entries()
  }
}
