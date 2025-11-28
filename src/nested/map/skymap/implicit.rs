use std::{
  collections::BTreeMap,
  f64::consts::PI,
  fs::File,
  io::{stdin, BufRead, BufReader, BufWriter, Error as IoError, Write},
  iter::{Enumerate, Map},
  marker::PhantomData,
  ops::{Add, AddAssign, Deref},
  path::Path,
  slice::Iter,
  vec::IntoIter,
};

use itertools::Itertools;
use rayon::{
  prelude::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelIterator, ParallelSlice, ParallelSliceMut,
  },
  ThreadPool,
};

use super::{SkyMap, SkyMapValue};
use crate::{
  n_hash,
  nested::{
    get,
    map::{
      fits::{error::FitsError, write::write_implicit_skymap_fits},
      mom::{
        impls::zvec::MomVecImpl, new_chi2_count_ref_merger, new_chi2_density_ref_merger, Mom,
        ZUniqHashT,
      },
      skymap::explicit::{ExplicitCountMap, ExplicitSkyMapBTree},
      HHash,
    },
    sort::get_hpx_opt,
  },
};

#[derive(Debug)]
pub struct ImplicitSkyMapArray<H: HHash, V: SkyMapValue> {
  pub depth: u8,
  pub values: Box<[V]>,
  _htype: PhantomData<H>,
}
impl<V: SkyMapValue> From<ImplicitSkyMapArray<u32, V>> for ImplicitSkyMapArray<u64, V> {
  fn from(value: ImplicitSkyMapArray<u32, V>) -> Self {
    Self::new(value.depth, value.values)
  }
}

impl<H: HHash, V: SkyMapValue> ImplicitSkyMapArray<H, V> {
  /// WARNING: we assume that the coherency between the depth and the number of elements in the
  ///array has already been tested.
  pub fn new(depth: u8, values: Box<[V]>) -> Self {
    assert_eq!(
      n_hash(depth) as usize,
      values.deref().len(),
      "Wrong implicit skymap size. Expected: {}. Actual: {}.",
      n_hash(depth),
      values.len()
    );
    Self {
      depth,
      values,
      _htype: PhantomData,
    }
  }

  pub fn into_explicit_map(self, null_value: V) -> ExplicitSkyMapBTree<H, V> {
    let depth = self.depth;
    let btreemap: BTreeMap<H, V> = self
      .owned_entries()
      .filter(move |(_k, v)| null_value.ne(v))
      .collect();
    ExplicitSkyMapBTree::new(depth, btreemap)
  }
}
impl<H: HHash, V: SkyMapValue + Send + Sync + AddAssign> ImplicitSkyMapArray<H, V> {
  pub fn par_add(mut self, rhs: Self) -> Self {
    self
      .values
      .as_parallel_slice_mut()
      .par_iter_mut()
      .zip_eq(rhs.values.into_vec().into_par_iter())
      .for_each(|(l, r)| l.add_assign(r));
    self
  }
}
impl<'a, H: HHash, V: SkyMapValue + 'a> SkyMap<'a> for ImplicitSkyMapArray<H, V> {
  type HashType = H;
  type ValueType = V;
  type ValuesIt = Iter<'a, Self::ValueType>;
  type EntriesIt = Map<Enumerate<Self::ValuesIt>, fn((usize, &V)) -> (H, &V)>;
  type OwnedEntriesIt = Map<Enumerate<IntoIter<Self::ValueType>>, fn((usize, V)) -> (H, V)>;

  fn depth(&self) -> u8 {
    self.depth
  }

  fn is_implicit(&self) -> bool {
    true
  }

  fn len(&self) -> usize {
    self.values.len()
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

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    Box::into_iter(self.values)
      .enumerate()
      .map(move |(h, v)| (H::from_usize(h), v))
  }
}

#[derive(Debug)]
pub struct ImplicitSkyMapArrayRef<'a, H: HHash, V: SkyMapValue> {
  depth: u8,
  values: &'a [V],
  _htype: PhantomData<H>,
}
impl<'a, H: HHash, V: SkyMapValue + 'a> ImplicitSkyMapArrayRef<'a, H, V> {
  /// WARNING: we assume that the coherency between the depth and the number of elements in the
  ///array has already been tested.
  pub fn new(depth: u8, values: &'a [V]) -> Self {
    assert_eq!(
      n_hash(depth) as usize,
      values.len(),
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

  /*pub fn par_add(mut self, rhs: Self, pool: &ThreadPool) -> Self {
    pool.install(|| {
      self
        .values
        .par_iter_mut()
        .zip_eq(rhs.values.into_par_iter())
        .for_each(|(l, r)| *l += r)
    });
    self
  }*/
}
impl<'a, H: HHash, V: SkyMapValue + Clone + 'a> SkyMap<'a> for ImplicitSkyMapArrayRef<'a, H, V> {
  type HashType = H;
  type ValueType = V;
  type ValuesIt = Iter<'a, Self::ValueType>;
  type EntriesIt = Map<Enumerate<Self::ValuesIt>, fn((usize, &V)) -> (H, &V)>;
  type OwnedEntriesIt = Map<Enumerate<Self::ValuesIt>, fn((usize, &V)) -> (H, V)>;

  fn depth(&self) -> u8 {
    self.depth
  }

  fn is_implicit(&self) -> bool {
    true
  }

  fn len(&self) -> usize {
    self.values.len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    &self.values[hash.as_()]
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.values.iter()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self
      .values
      .iter()
      .enumerate()
      .map(move |(h, v)| (H::from_usize(h), v))
  }

  // Make a owned_values method!!

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    self
      .values
      .iter()
      .enumerate()
      .map(move |(h, v)| (H::from_usize(h), v.clone()))
  }
}

/// SkyMap implementation use to store counts.
#[derive(Debug)]
pub struct ImplicitCountMap(ImplicitSkyMapArray<u64, u32>);
impl From<ImplicitSkyMapArray<u64, u32>> for ImplicitCountMap {
  fn from(value: ImplicitSkyMapArray<u64, u32>) -> Self {
    Self(value)
  }
}
impl ImplicitCountMap {
  pub fn as_implicit_skymap_array(&self) -> &ImplicitSkyMapArray<u64, u32> {
    &self.0
  }
  pub fn into_implicit_skymap_array(self) -> ImplicitSkyMapArray<u64, u32> {
    self.0
  }
  pub fn into_explicit_skymap(self, null_value: u32) -> ExplicitCountMap<u64> {
    self.0.into_explicit_map(null_value).into()
  }
  /// Build a count skymap from an iterator over HEALPix cells at the given depth.
  /// # Panics
  /// * if `depth > 12`.
  pub fn from_hash_values<I>(depth: u8, hash_values_it: I) -> Self
  where
    I: Iterator<Item = u32>,
  {
    assert!(
      depth < 13,
      "Wrong count map input depth. Expected: < 13. Actual: {}",
      depth
    );
    let mut counts = vec![0_u32; n_hash(depth) as usize].into_boxed_slice();
    for h in hash_values_it {
      counts[h as usize] += 1;
    }
    Self(ImplicitSkyMapArray::new(depth, counts))
  }

  /// Build a count skymap from an iterator over position ( (ra, dec), in radian).
  /// # Panics
  /// * if `depth > 12`.
  pub fn from_positions<I>(depth: u8, pos_it_rad: I) -> Self
  where
    I: Iterator<Item = (f64, f64)>,
  {
    let layer = get(depth);
    Self::from_hash_values(depth, pos_it_rad.map(|(l, b)| layer.hash(l, b) as u32))
  }

  pub fn par_add(mut self, rhs: Self) -> Self {
    self.0 = self.0.par_add(rhs.0);
    self
  }

  pub fn to_dens_map(&self) -> ImplicitDensityMap {
    let depth = self.depth();
    let on_over_area = (n_hash(depth) >> 2) as f64 / PI;
    ImplicitDensityMap(ImplicitSkyMapArray::new(
      depth,
      self
        .0
        .values
        .iter()
        .map(|count| *count as f64 * on_over_area)
        .collect::<Vec<f64>>()
        .into_boxed_slice(),
    ))
  }

  /// For custom thread pool, use `thread_pool.install(|| to_dens_map_par())`
  pub fn to_dens_map_par(&self) -> ImplicitDensityMap {
    let depth = self.depth();
    let on_over_area = (n_hash(depth) >> 2) as f64 / PI;
    ImplicitDensityMap(ImplicitSkyMapArray::new(
      depth,
      self
        .0
        .values
        .par_iter()
        .map(|count| *count as f64 * on_over_area)
        .collect::<Vec<f64>>()
        .into_boxed_slice(),
    ))
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    write_implicit_skymap_fits(writer, self.as_implicit_skymap_array().values.deref())
  }

  pub fn to_fits_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits(BufWriter::new(file)))
  }

  /// Merge sibling cells till the sum of their counts remains lower than
  /// the given threshold.
  pub fn to_upper_count_threshold_mom(&self, upper_count_threshold: u32) -> MomVecImpl<u64, u32> {
    let merger = |_depth: u8, _hash: u64, [n0, n1, n2, n3]: [&u32; 4]| -> Option<u32> {
      let sum = *n0 + *n1 + *n2 + *n3;
      if sum < upper_count_threshold {
        Some(sum)
      } else {
        None
      }
    };
    MomVecImpl::from_skymap_ref(&self.0, merger)
  }

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
    /*let chi2_merger = |_depth: u8, _hash: u64, [n0, n1, n2, n3]: [&u32; 4]| -> Option<u32> {
      // With Poisson distribution:
      // * mu_i = source density in cell i
      // * sigma_i = sqrt(mu_i)
      // weight_i = 1 / sigma_i^2
      // mu_e = weighted_mean = ( sum_{1=1}^4 weight_i * mu_i ) / ( sum_{1=1}^4 weight_i )
      //                      = 4 / ( sum_{1=1}^4 1/mu_i )
      // V_e^{-1} = sum_{1=1}^4 1/mu_i
      // Applying Pineau et al. 2017:
      // => sum_{1=1}^4 (mu_i - mu_e)^2 / mu_i = ... = (sum_{1=1}^4 mu_i) - 4 * mu_e
      // Normal law product of 4 1D normal laws and apply Pineau 2017 to find the above equation:
      // 1/sqrt(2 * pi) * exp[ -1/2 * sum_{i=1}^4 ( (x - mu_i)/sqrt(sigma_i) )^2] / sqrt(prod_{i=1}^4 sigma_i)

      let mu0 = *n0 as f64;
      let mu1 = *n1 as f64;
      let mu2 = *n2 as f64;
      let mu3 = *n3 as f64;

      let sum = mu0 + mu1 + mu2 + mu3;
      let weighted_var_inv =
        1.0 / mu0.max(1.0) + 1.0 / mu1.max(1.0) + 1.0 / mu2.max(1.0) + 1.0 / mu3.max(1.0);
      let weighted_mean = 4.0 / weighted_var_inv;
      let chi2_of_3dof = sum - 4.0 * weighted_mean;

      // chi2 3 dof:
      // 90.0% =>  6.251
      // 95.0% =>  7.815
      // 97.5% =>  9.348
      // 99.0% => 11.345
      // 99.9% => 16.266
      if chi2_of_3dof < chi2_of_3dof_threshold {
        Some(*n0 + *n1 + *n2 + *n3)
      } else {
        None
      }
    };*/
    let chi2_merger = new_chi2_count_ref_merger(chi2_of_3dof_threshold);
    let mom = match depth_threshold {
      None => MomVecImpl::from_skymap_ref(&self.0, chi2_merger),
      Some(depth_threshold) => MomVecImpl::from_skymap_ref(
        &self.0,
        |depth: u8, hash: u64, na: [&u32; 4]| -> Option<u32> {
          if depth >= depth_threshold {
            chi2_merger(depth, hash, na)
          } else {
            None
          }
        },
      ),
    };
    // Create a new MOM transforming number of sources into densities.
    MomVecImpl::from_map(mom, |z, v| {
      v as f64 / (4.0 * PI / (n_hash(u64::depth_from_zuniq(z))) as f64)
    })
  }

  // to_png
  // to_fits
}
impl Add for ImplicitCountMap {
  type Output = Self;

  fn add(mut self, rhs: Self) -> Self::Output {
    self
      .0
      .values
      .iter_mut()
      .zip_eq(rhs.0.values.iter())
      .for_each(|(l, r)| *l += r);
    self
  }
}
impl<'a> SkyMap<'a> for ImplicitCountMap {
  type HashType = u64;
  type ValueType = u32;
  type ValuesIt = Iter<'a, Self::ValueType>;
  type EntriesIt = Map<Enumerate<Self::ValuesIt>, fn((usize, &u32)) -> (u64, &u32)>;
  type OwnedEntriesIt = Map<Enumerate<IntoIter<Self::ValueType>>, fn((usize, u32)) -> (u64, u32)>;

  fn depth(&self) -> u8 {
    self.as_implicit_skymap_array().depth
  }

  fn is_implicit(&self) -> bool {
    true
  }

  fn len(&self) -> usize {
    self.as_implicit_skymap_array().values.len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    &self.as_implicit_skymap_array().values.deref()[hash as usize]
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.as_implicit_skymap_array().values.deref().iter()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self
      .as_implicit_skymap_array()
      .values
      .deref()
      .iter()
      .enumerate()
      .map(move |(h, v)| (u64::from_usize(h), v))
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    Box::into_iter(self.into_implicit_skymap_array().values)
      .enumerate()
      .map(move |(h, v)| (u64::from_usize(h), v))
  }
}

/// SkyMap implementation use to store counts, rely on u32 HEALPix index instead of u64.
#[derive(Debug)]
pub struct ImplicitCountMapU32(ImplicitSkyMapArray<u32, u32>);
impl From<ImplicitSkyMapArray<u32, u32>> for ImplicitCountMapU32 {
  fn from(value: ImplicitSkyMapArray<u32, u32>) -> Self {
    Self(value)
  }
}
impl ImplicitCountMapU32 {
  pub fn as_implicit_skymap_array(&self) -> &ImplicitSkyMapArray<u32, u32> {
    &self.0
  }
  pub fn into_implicit_skymap_array(self) -> ImplicitSkyMapArray<u32, u32> {
    self.0
  }
  pub fn into_explicit_skymap(self, null_value: u32) -> ExplicitCountMap<u32> {
    self.0.into_explicit_map(null_value).into()
  }
  /// Build a count skymap from an iterator over HEALPix cells at the given depth.
  /// # Panics
  /// * if `depth > 12`.
  pub fn from_hash_values<I>(depth: u8, hash_values_it: I) -> Self
  where
    I: Iterator<Item = u32>,
  {
    assert!(
      depth < 13,
      "Wrong count map input depth. Expected: < 13. Actual: {}",
      depth
    );
    let mut counts = vec![0_u32; n_hash(depth) as usize].into_boxed_slice();
    for h in hash_values_it {
      counts[h as usize] += 1;
    }
    Self(ImplicitSkyMapArray::new(depth, counts))
  }

  /// Build a count skymap from an iterator over position ( (ra, dec), in radian).
  /// # Panics
  /// * if `depth > 12`.
  pub fn from_positions<I>(depth: u8, pos_it_rad: I) -> Self
  where
    I: Iterator<Item = (f64, f64)>,
  {
    assert!(
      depth < 13,
      "Wrong count map input depth. Expected: < 13. Actual: {}",
      depth
    );
    let layer = get(depth);
    let mut counts = vec![0_u32; layer.n_hash as usize].into_boxed_slice();
    for (l, b) in pos_it_rad {
      counts[layer.hash(l, b) as usize] += 1;
    }
    Self(ImplicitSkyMapArray::new(depth, counts))
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
    /*let mut it = BufReader::new(File::open(&path)?).lines().peekable();
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
    Self::from_csv_it_par(it, ilon, ilat, separator, depth, chunk_size, thread_pool)*/
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
    let n_hash = layer.n_hash as usize;
    let n_thread = thread_pool.current_num_threads();
    let hpx = get_hpx_opt(ilon, ilat, separator, layer);
    let count_map_fn = |chunk: Vec<String>| {
      chunk
        .par_chunks((chunk.len() / (n_thread << 2)).max(10_000))
        .map(|elems| {
          ImplicitCountMapU32::from_hash_values(
            depth,
            elems.iter().filter_map(&hpx).map(|h| h as u32),
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
    let mut count_map = Self(ImplicitSkyMapArray::new(
      depth,
      vec![0_u32; n_hash].into_boxed_slice(),
    ));
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

  pub fn par_add(mut self, rhs: Self) -> Self {
    /* self
    .0
    .values
    .par_iter_mut()
    .zip_eq(rhs.0.values.into_par_iter())
    .for_each(|(l, r)| *l += r)*/
    self.0 = self.0.par_add(rhs.0);
    self
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    write_implicit_skymap_fits(writer, self.as_implicit_skymap_array().values.deref())
  }

  pub fn to_fits_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(&path)
      .map_err(|err| FitsError::IoWithPath {
        path: path.as_ref().to_string_lossy().into(),
        err,
      })
      .and_then(|file| self.to_fits(BufWriter::new(file)))
  }

  pub fn to_dens_map(&self) -> ImplicitDensityMap {
    let depth = self.depth();
    let one_over_area = (n_hash(depth) >> 2) as f64 / PI;
    ImplicitDensityMap(ImplicitSkyMapArray::new(
      depth,
      self
        .0
        .values
        .iter()
        .map(|count| *count as f64 * one_over_area)
        .collect::<Vec<f64>>()
        .into_boxed_slice(),
    ))
  }

  pub fn to_dens_map_par(&self) -> ImplicitDensityMap {
    let depth = self.depth();
    let one_over_area = (n_hash(depth) >> 2) as f64 / PI;
    ImplicitDensityMap(ImplicitSkyMapArray::new(
      depth,
      self
        .0
        .values
        .par_iter()
        .map(|count| *count as f64 * one_over_area)
        .collect::<Vec<f64>>()
        .into_boxed_slice(),
    ))
  }

  /// # Params
  /// * `chi2_of_3dof_threshold`: threshold on the value of the chi square distribution with 3
  /// degrees of freedom below which we consider the 4 values of 4 sibling cells as coming
  /// from the same normal distribution which mean and variance comes from a poisson distribution.
  /// Here a few typical values corresponding the the given completeness:
  /// * Completeness = 90.0% =>  6.251
  /// * Completeness = 95.0% =>  7.815
  /// * Completeness = 97.5% =>  9.348
  /// * Completeness = 99.0% => 11.345
  /// * Completeness = 99.9% => 16.266
  pub fn to_chi2_mom(&self, chi2_of_3dof_threshold: f64) -> MomVecImpl<u32, u32> {
    let chi2_merger = |_depth: u8, _hash: u32, [n0, n1, n2, n3]: [&u32; 4]| -> Option<u32> {
      // With Poisson distribution:
      // * mu_i = source density in cell i
      // * sigma_i = sqrt(mu_i)
      // weight_i = 1 / sigma_i^2
      // mu_e = weighted_mean = ( sum_{1=1}^4 weight_i * mu_i ) / ( sum_{1=1}^4 weight_i )
      //                      = 4 / ( sum_{1=1}^4 1/mu_i )
      // V_e^{-1} = sum_{1=1}^4 1/mu_i
      // Applying Pineau et al. 2017:
      // => sum_{1=1}^4 (mu_i - mu_e)^2 / mu_i = ... = (sum_{1=1}^4 mu_i) - 4 * mu_e
      // Normal law product of 4 1D normal laws and apply Pineau 2017 to find the above equation:
      // 1/sqrt(2 * pi) * exp[ -1/2 * sum_{i=1}^4 ( (x - mu_i)/sqrt(sigma_i) )^2] / sqrt(prod_{i=1}^4 sigma_i)

      let mu0 = *n0 as f64;
      let mu1 = *n1 as f64;
      let mu2 = *n2 as f64;
      let mu3 = *n3 as f64;

      let sum = mu0 + mu1 + mu2 + mu3;
      let weighted_var_inv =
        1.0 / mu0.max(1.0) + 1.0 / mu1.max(1.0) + 1.0 / mu2.max(1.0) + 1.0 / mu3.max(1.0);
      let weighted_mean = 4.0 / weighted_var_inv;
      let chi2_of_3dof = sum - 4.0 * weighted_mean;

      // chi2 3 dof:
      // 90.0% =>  6.251
      // 95.0% =>  7.815
      // 97.5% =>  9.348
      // 99.0% => 11.345
      // 99.9% => 16.266
      if chi2_of_3dof < chi2_of_3dof_threshold {
        Some(*n0 + *n1 + *n2 + *n3)
      } else {
        None
      }
    };
    /*let mom = MomVecImpl::from_skymap_ref(&self.0, chi2_merger);
    // Create a new MOM transforming number of sources into densities.
    MomVecImpl::from(mom, |z, v| {
      v as f64 / (4.0 * PI / (n_hash(u32::depth_from_zuniq(z))) as f64)
    })*/
    MomVecImpl::from_skymap_ref(&self.0, chi2_merger)
  }

  // to_png
  // to_fits
}
impl Add for ImplicitCountMapU32 {
  type Output = Self;

  fn add(mut self, rhs: Self) -> Self::Output {
    self
      .0
      .values
      .iter_mut()
      .zip_eq(rhs.0.values.iter())
      .for_each(|(l, r)| *l += r);
    self
  }
}
impl<'a> SkyMap<'a> for ImplicitCountMapU32 {
  type HashType = u32;
  type ValueType = u32;
  type ValuesIt = Iter<'a, Self::ValueType>;
  type EntriesIt = Map<Enumerate<Self::ValuesIt>, fn((usize, &u32)) -> (u32, &u32)>;
  type OwnedEntriesIt = Map<Enumerate<IntoIter<Self::ValueType>>, fn((usize, u32)) -> (u32, u32)>;

  fn depth(&self) -> u8 {
    self.as_implicit_skymap_array().depth
  }

  fn is_implicit(&self) -> bool {
    true
  }

  fn len(&self) -> usize {
    self.as_implicit_skymap_array().values.len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    &self.as_implicit_skymap_array().values.deref()[hash as usize]
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.as_implicit_skymap_array().values.deref().iter()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self
      .as_implicit_skymap_array()
      .values
      .deref()
      .iter()
      .enumerate()
      .map(move |(h, v)| (u32::from_usize(h), v))
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    Box::into_iter(self.into_implicit_skymap_array().values) // a normal into_iter() calls the slice impl
      .enumerate()
      .map(move |(h, v)| (u32::from_usize(h), v))
  }
}

/// SkyMap implementation use to store densities.
#[derive(Debug)]
pub struct ImplicitDensityMap(ImplicitSkyMapArray<u32, f64>);
impl From<ImplicitSkyMapArray<u32, f64>> for ImplicitDensityMap {
  fn from(value: ImplicitSkyMapArray<u32, f64>) -> Self {
    Self(value)
  }
}
impl ImplicitDensityMap {
  pub fn as_implicit_skymap_array(&self) -> &ImplicitSkyMapArray<u32, f64> {
    &self.0
  }
  pub fn into_implicit_skymap_array(self) -> ImplicitSkyMapArray<u32, f64> {
    self.0
  }
  /// Build a count skymap from an iterator over position ( (ra, dec), in radian).
  /// # Panics
  /// * if `depth > 12`.
  pub fn from_positions<I>(depth: u8, pos_it_rad: I) -> Self
  where
    I: Iterator<Item = (f64, f64)>,
  {
    assert!(
      depth < 13,
      "Wrong count map input depth. Expected: < 13. Actual: {}",
      depth
    );
    let layer = get(depth);
    let mut densities = vec![0_f64; layer.n_hash as usize].into_boxed_slice();
    let one_over_cell_area = layer.n_hash as f64 / (4.0 * PI);
    for (l, b) in pos_it_rad {
      densities[layer.hash(l, b) as usize] += one_over_cell_area;
    }
    Self(ImplicitSkyMapArray::new(depth, densities))
  }

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
  ) -> MomVecImpl<u32, f64> {
    let chi2_merger = new_chi2_density_ref_merger(chi2_of_3dof_threshold);
    match depth_threshold {
      None => MomVecImpl::from_skymap_ref(&self.0, chi2_merger),
      Some(depth_threshold) => MomVecImpl::from_skymap_ref(
        &self.0,
        |depth: u8, hash: u32, [n0, n1, n2, n3]: [&f64; 4]| -> Option<f64> {
          if depth >= depth_threshold {
            chi2_merger(depth, hash, [n0, n1, n2, n3])
          } else {
            None
          }
        },
      ),
    }
  }

  pub fn to_fits<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    write_implicit_skymap_fits(writer, self.as_implicit_skymap_array().values.deref())
  }

  pub fn to_fits_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits(BufWriter::new(file)))
  }

  // to_png
}
