/*
Strategy:
* create 48 directories if 68 files each (=> 3072 files)
* create a 'cumul count' at a given order at a given depth (e.g. 8)
* spilt too large files in xx ranges computed from 'cumul counts'
* read sorting each file when reading it
*/

use std::{
  cmp::Ordering,
  error::Error,
  fs::{self, read_dir, remove_dir, remove_file, File, OpenOptions},
  io::{
    self, BufRead, BufReader, BufWriter, Error as IoError, ErrorKind as IoErrorKind, Lines, Seek,
    Write,
  },
  iter::{Peekable, Zip},
  marker::{Send, Sync},
  num::ParseIntError,
  ops::Range,
  path::{Path, PathBuf},
  time::{Duration, SystemTime},
  vec::IntoIter,
};

use bincode::{
  self,
  config::{FixintEncoding, WithOtherIntEncoding},
  DefaultOptions, Options,
};
use log::{debug, error, info, warn};
use rayon::{prelude::ParallelSliceMut, ThreadPool};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::nested::{
  get,
  map::skymap::{CountMapU32, SkyMap},
  n_hash,
};

/// Type defining a row that can be sorted internally.
pub trait IntSortable: Send {}
/// All types that implement `Send` are `IntSortable`.
impl<A: Send> IntSortable for A {}

/// Type defining a row than can be sorted externally, i.e. need to be serialized and deserialized.
/// Keep it mind that to ser/de  `&[u8]` or `Vec<u8>`, using
/// [serde_bytes](https://docs.rs/serde_bytes/latest/serde_bytes/) can boost performances.
pub trait ExtSortable: IntSortable + Serialize + DeserializeOwned {}
/// All types that implement `IntSortable`, `Serialize` and `DeserializeOwned` are `ExtSortable`.
impl<A: IntSortable + Serialize + DeserializeOwned> ExtSortable for A {}

// Declare a SortError ?

/// Sort internally (i.e. in memory) the given array of elements according to the order 29 HEALPix
/// cell computed using the given `hpx29` function.
/// If the provided input parameter `n_thread` is `None`, the maximum number of thread available
/// on the machine is used. Use `n_thread = Some(1)` for a monothreaded execution (without having
/// to instantiate a thread pool).
/// # Note
/// Having a costly `hpx29` method is ok: internally it will be called only once per row
/// (not multiple times durng the sort process).
pub fn hpx_internal_sort<T, F>(mut elems: &mut [T], hpx29: F, n_threads: Option<usize>)
where
  T: IntSortable,
  F: Fn(&T) -> u64 + Sync,
{
  match n_threads {
    Some(1) => elems.sort_by_cached_key(hpx29),
    _ => get_thread_pool(n_threads).install(|| elems.par_sort_by_cached_key(&hpx29)),
  }
}

/// Returns a pool of thread containing:
/// * the maximum number of threads available on the machine, if `n_threads` is `None`
/// * the given number of threads in the `n_threads` option if it is `Some`.
fn get_thread_pool(n_threads: Option<usize>) -> ThreadPool {
  let mut pool_builder = rayon::ThreadPoolBuilder::new();
  if let Some(n_threads) = n_threads {
    pool_builder = pool_builder.num_threads(n_threads);
  }
  pool_builder.build().unwrap()
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SimpleExtSortInfo {
  /// Maximum number of elements per chunk (i.e. per file).
  n_elems_per_chunk: u32,
  /// Number of used threads.
  n_threads: Option<usize>,
  /// Total number of rows.
  n_tot: u64,
  /// Healpix depth of the ranges in filenames.
  depth: u8,
  /// Healpix ranges at `depth`, one range per chunk/file, together with the number of rows in the chunk/file.
  ordered_ranges_counts: Vec<(Range<u32>, u32)>,
}
impl SimpleExtSortInfo {
  fn new(
    n_elems_per_chunk: u32,
    n_threads: Option<usize>,
    n_tot: u64,
    depth: u8,
    ordered_ranges_counts: Vec<(Range<u32>, u32)>,
  ) -> Self {
    Self {
      n_elems_per_chunk,
      n_threads,
      n_tot,
      depth,
      ordered_ranges_counts,
    }
  }
}

pub struct SimpleExtSortParams {
  /// Directory containing the temporary directories/files.
  tmp_dir: PathBuf,
  /// Number of items to be sorted at once in memory, correspond to the maximum number of items
  /// in a temporary file at second pass.
  n_elems_per_chunk: u32,
  /// Provide a given number of threads, else use all threads available.
  n_threads: Option<usize>,
}
impl Default for SimpleExtSortParams {
  fn default() -> Self {
    Self {
      tmp_dir: PathBuf::from(".sort_tmp/"),
      n_elems_per_chunk: 10_000_000,
      n_threads: None,
    }
  }
}
impl SimpleExtSortParams {
  const ALL_FILENAME: &'static str = "hpxsort.all.unsorted.bin";
  const PREFIX: &'static str = "hpxsort.";
  const INFO_FILENAME: &'static str = "hpxsort.info.toml";
  const COUNT_MAP_FILENAME: &'static str = "hpxsort.countmap.fits";
  const SUCCESS_FILENAME: &'static str = "hpxsort.success";
  const SUFFIX: &'static str = ".unsorted.bin";

  fn clean<P: AsRef<Path>>(tmp_dir: P) -> Result<(), IoError> {
    let mut path = PathBuf::from(tmp_dir.as_ref());
    debug!("PATH: {:?}", &path);
    assert!(path.file_name().is_some());
    // Remove SUCCESS first (since we are going to remove files, we have to remove the SUCCESS flag first)
    path.push(Self::SUCCESS_FILENAME);
    debug!("Remove file {:?} if exists.", &path);
    fs::exists(&path).and_then(|exists| if exists { remove_file(&path) } else { Ok(()) })?;
    // Remove count map file
    path.set_file_name(Self::COUNT_MAP_FILENAME);
    debug!("Remove file {:?} if exists.", &path);
    fs::exists(&path).and_then(|exists| if exists { remove_file(&path) } else { Ok(()) })?;
    // Remove info file
    path.set_file_name(Self::INFO_FILENAME);
    debug!("Remove file {:?} if exists.", &path);
    fs::exists(&path).and_then(|exists| if exists { remove_file(&path) } else { Ok(()) })?;
    // Remove all file if exists
    path.set_file_name(Self::ALL_FILENAME);
    debug!("Remove file {:?} if exists.", &path);
    fs::exists(&path).and_then(|exists| if exists { remove_file(&path) } else { Ok(()) })?;
    // Remove temporary files
    path.pop();
    for file in read_dir(&path)? {
      let file_path = file?.path();
      if file_path.is_file()
        && file_path
          .file_name()
          .and_then(|os_str| os_str.to_str())
          .and_then(|name| SimpleExtSortParams::parse_tmp_file_name(name))
          .is_some()
      {
        debug!("Remove file {:?} if exists.", &file_path);
        remove_file(file_path)?;
      }
    }
    // Finally remove the directory if empty
    let tmp_dir_ref = tmp_dir.as_ref();
    // debug_assert_eq!(tmp_dir_ref.as_os_str(), path.as_os_str());
    debug!("Remove dir {:?} if exists and is empty.", tmp_dir_ref);
    if fs::exists(tmp_dir_ref).and_then(move |exists| {
      if exists {
        read_dir(tmp_dir_ref).map(|mut it| it.next().is_none())
      } else {
        Ok(false)
      }
    })? {
      remove_dir(tmp_dir_ref)?;
    } else {
      warn!("Unable to remove directory {:?}", tmp_dir_ref);
    }
    Ok(())
  }

  fn set_tmp_dir(mut self, tmp_dir: PathBuf) -> Self {
    self.tmp_dir = tmp_dir;
    self
  }

  fn set_n_elems_per_chunk(mut self, n_elems_per_chunk: u32) -> Self {
    self.n_elems_per_chunk = n_elems_per_chunk;
    self
  }

  fn set_n_threads(mut self, n_threads: usize) -> Self {
    self.n_threads = Some(n_threads);
    self
  }

  fn from_dir_and_info(tmp_dir: PathBuf, info: &SimpleExtSortInfo) -> Self {
    Self {
      tmp_dir,
      n_elems_per_chunk: info.n_elems_per_chunk,
      n_threads: info.n_threads.clone(),
    }
  }

  fn create_tmp_dir(&self) -> Result<(), IoError> {
    debug!("Create tmp dir: {}", &self.tmp_dir.to_string_lossy());
    fs::create_dir_all(&self.tmp_dir)
  }

  /// Open or create, in write mode (replacing previous content if any).
  fn create_file_all(&self) -> Result<File, IoError> {
    let mut path = self.tmp_dir.clone();
    path.push(Self::ALL_FILENAME);
    debug!(
      "Create or open to overwrite file: {}",
      path.to_string_lossy()
    );
    OpenOptions::new()
      .append(false)
      .write(true)
      .create(true)
      .open(path)
  }
  fn open_file_all(&self) -> Result<File, IoError> {
    let mut path = self.tmp_dir.clone();
    path.push(Self::ALL_FILENAME);
    File::open(path)
  }

  fn create_file_or_open_to_append(&self, depth: u8, range: &Range<u32>) -> Result<File, IoError> {
    let n_digits = (1 + n_hash(depth).ilog10()) as usize;
    let filename = format!(
      "{}d{}f{:0n_digits$}t{:0n_digits$}{}",
      Self::PREFIX,
      depth,
      range.start,
      range.end,
      Self::SUFFIX
    );
    let mut path = self.tmp_dir.clone();
    path.push(filename);
    debug!("Create or open to append file: {}", path.to_string_lossy());
    OpenOptions::new().append(true).create(true).open(path)
  }

  // Rea/Write info file

  fn create_info_file_path(&self) -> PathBuf {
    let mut path = self.tmp_dir.clone();
    path.push(Self::INFO_FILENAME);
    path
  }

  fn create_info_file_path_gen(path: &PathBuf) -> PathBuf {
    let mut path = path.clone();
    path.push(Self::INFO_FILENAME);
    path
  }

  fn write_info(
    &self,
    n_tot: u64,
    depth: u8,
    ordered_ranges_counts: Vec<(Range<u32>, u32)>,
  ) -> Result<(), Box<dyn Error>> {
    let info = SimpleExtSortInfo::new(
      self.n_elems_per_chunk,
      self.n_threads.clone(),
      n_tot,
      depth,
      ordered_ranges_counts,
    );
    let content = toml::to_string_pretty(&info)?;
    let path = self.create_info_file_path();
    fs::write(&path, content)
      .map_err(|e| format!("Error writing file {}: {:?}.", path.to_string_lossy(), e).into())
  }

  fn read_info(&self) -> Result<SimpleExtSortInfo, Box<dyn Error>> {
    let path = self.create_info_file_path();
    fs::read_to_string(&path)
      .map_err(|e| format!("Error reading file {}: {:?}.", path.to_string_lossy(), e).into())
      .and_then(|content| toml::from_str(&content).map_err(|e| e.into()))
  }

  fn read_info_gen(path: &PathBuf) -> Result<SimpleExtSortInfo, Box<dyn Error>> {
    let path = Self::create_info_file_path_gen(path);
    fs::read_to_string(&path)
      .map_err(|e| format!("Error reading file {}: {:?}.", path.to_string_lossy(), e).into())
      .and_then(|content| toml::from_str(&content).map_err(|e| e.into()))
  }

  fn create_countmap_file_path(&self) -> PathBuf {
    let mut path = self.tmp_dir.clone();
    path.push(Self::COUNT_MAP_FILENAME);
    path
  }

  // Read/Write SUCESS file

  fn create_ok_file_path(&self) -> PathBuf {
    let mut path = self.tmp_dir.clone();
    path.push(Self::SUCCESS_FILENAME);
    path
  }
  fn write_ok(&self) -> Result<(), Box<dyn Error>> {
    let path = self.create_ok_file_path();
    fs::write(&path, "")
      .map_err(|e| format!("Error writing file {}: {:?}.", path.to_string_lossy(), e).into())
  }
  fn ok_file_exists(&self) -> Result<bool, Box<dyn Error>> {
    let path = self.create_ok_file_path();
    fs::exists(&path).map_err(|e| {
      format!(
        "Error cheking for file {}: {:?}.",
        path.to_string_lossy(),
        e
      )
      .into()
    })
  }

  //

  /// Returns the depth and the range encoded in the filename.
  fn parse_tmp_file_name(name: &str) -> Option<(u8, Range<u32>)> {
    name
      .strip_prefix(format!("{}d", Self::PREFIX).as_str())
      .and_then(|s| s.strip_suffix(Self::SUFFIX))
      .and_then(|s| s.split_once('f'))
      .and_then(|(depth, s)| s.split_once('t').map(|(from, to)| (depth, from, to)))
      .and_then(|(depth, from, to)| {
        match (depth.parse::<u8>(), from.parse::<u32>(), to.parse::<u32>()) {
          (Ok(d), Ok(f), Ok(t)) => Some((d, f..t)),
          _ => None,
        }
      })
  }

  /// Returns a list of file `path`, `depth` dans HEALPix `range` sorted by increasing range starting
  /// hash value.
  fn get_ordered_files_in_tmp_dir(&self) -> Result<Vec<(PathBuf, u8, Range<u32>)>, IoError> {
    read_dir(&self.tmp_dir).map(|it| {
      let mut path_depth_range_vec: Vec<(PathBuf, u8, Range<u32>)> = it
        .filter_map(|res_dir_entry| {
          res_dir_entry.ok().and_then(|dir_entry| {
            let path = dir_entry.path();
            if path.is_file() {
              path
                .file_name()
                .and_then(|os_str| os_str.to_str())
                .and_then(|name| SimpleExtSortParams::parse_tmp_file_name(name))
                .map(|(depth, range)| (path, depth, range))
            } else {
              None
            }
          })
        })
        .collect();
      // TODO: check that ranges are non-overlapping and depth is always the same?
      path_depth_range_vec.sort_by(|(_, _, rl), (_, _, rr)| rl.start.cmp(&rr.start));
      path_depth_range_vec.into_iter().collect()
    })
  }
}

///
/// # Params
/// * `it`: the iterator over all rows to be sorted
/// * `count_map`: the skymap containing the number of rows inside each cell. It **must**
///   be consistent with the rows we get from the input iterator `it`.
///   If not known, a first pass is required, see [hpx_external_sort]{#hpx_external_sort}.
/// * `hpx29`: the function computing HEALPix hash at depth 29 from a row reference.
/// * `sort_param`: parameters specific to the sort algorithm
///
/// # Info
/// * Why taking in input an iterator of `Result`?
///     + because we expect large iteration that do no hold in memory, meaning that reading
///       from a file (or the network) is involved, and thus an error may occur during the iteration.
///
/// # Algo
/// This external sort is optimal: the input data is written and then read exactly once.
/// This is possible only providing its HEALPix distribution.
/// It probably means a first pass is required to compute the HEALPix distribution.
/// When reading from a file, we have to read twice the input file:
/// * first to compute the HEALPix cells distribution;
/// * second to provide this method input iterator.
/// When reading from an iterator (distribution unknown), we additionally have to
/// write the input (possibly pre-sorting it) before re-reading it.
///
/// # Warning
/// Fails if the HEALPix cell distribution in the input iterator does not fit the one provided
/// in the input `map`.
///
pub fn hpx_external_sort_with_knowledge<'a, T, E, I, S, F>(
  it: I,
  count_map: &'a S,
  hpx29: F,
  sort_params: Option<SimpleExtSortParams>,
) -> Result<impl Iterator<Item = Result<T, Box<dyn Error>>>, Box<dyn Error>>
where
  T: ExtSortable,
  E: Error + 'static,
  I: Iterator<Item = Result<T, E>>,
  S: SkyMap<'a, HashType = u32, ValueType = u32>,
  F: Fn(&T) -> u64 + Sync,
{
  let params = sort_params.unwrap_or_default();
  let tmp_dir = params.tmp_dir.clone();
  hpx_external_sort_with_knowledge_write_tmp(it, count_map, &hpx29, Some(params))
    .and_then(|()| hpx_external_sort_with_knowledge_read_tmp(hpx29, tmp_dir))
}

/// See [hpx_external_sort_with_knowledge](#hpx_external_sort_with_knowledge).
/// This method perform the first part of the operation, i.e. creating and writing the (temporary)
/// file structure on the disk.
/// The second part consits in reading (and sorting in memory) the temporary files iteratively.
///
/// We separate the two parts to:
/// * be able to only perform the read operation if an error occurs while provessing the output iterator.
/// * allow for file-by-file sorting and serialization to create a structure on the disk (possibly indexing each file)
///   for direct use.
pub fn hpx_external_sort_with_knowledge_write_tmp<'a, T, E, I, S, F>(
  mut it: I,
  count_map: &'a S,
  hpx29: F,
  sort_params: Option<SimpleExtSortParams>,
) -> Result<(), Box<dyn Error>>
where
  T: ExtSortable,
  E: Error + 'static,
  I: Iterator<Item = Result<T, E>>,
  S: SkyMap<'a, HashType = u32, ValueType = u32>,
  F: Fn(&T) -> u64 + Sync,
{
  // Get general params
  let params = sort_params.unwrap_or_default();
  let depth = count_map.depth();
  let dd = 29 - depth;
  let twice_dd = dd << 1;
  // Config bincode
  let bincode = get_bincode();
  // .allow_trailing_bytes();
  // Create the temporary directory
  params.create_tmp_dir()?;
  // Init thread pool
  let pool = get_thread_pool(params.n_threads);
  // Get 'chunks', i.e. files corresponding to (ordered) HEALPix ranges containing a maximum of `n_max` rows.
  let n_max = params.n_elems_per_chunk;
  let (n_tot, mut ranges_counts) = skymap2ranges(count_map, n_max);
  params.write_info(n_tot, depth, ranges_counts.clone())?;

  let n_max = n_max as usize;
  let n_chunks = ranges_counts.len();
  // Now iterates on rows and sort them before writing into 'chunks' files.
  let mut entries_buf: Vec<T> = Vec::with_capacity(n_max);
  let mut n_tot_it = 0_u64;
  loop {
    // Populate the row buffer from the iterator.
    // I don't want ot map and collect to be able to reuse the same temporary Vec.
    // Without `Result`, we could have used directly `entries_buf.extend((&mut it).take(n_max));`
    let tstart = SystemTime::now();
    // Here, instead of adding into the existing vector, get  n new vector from a sender/receiver!!
    // to continue loading data while the parllel sort occurs!!
    for res in (&mut it).take(n_max) {
      entries_buf.push(res?);
    }
    // MAKE AN ALGO WRITING/READING NEXT WHILE SORTING
    let n = entries_buf.len();
    debug!(
      "Row chunk of size {} populated to be sorted in {} ms",
      n,
      SystemTime::now()
        .duration_since(tstart)
        .unwrap_or_default()
        .as_millis()
    );
    // In case n_tot % n_max == 0, quit the loop
    if n == 0 {
      break;
    }
    // Parallel sort
    let tstart = SystemTime::now();
    pool.install(|| entries_buf.par_sort_by_cached_key(&hpx29));
    debug!(
      "Parallel sort performed in {} ms",
      SystemTime::now()
        .duration_since(tstart)
        .unwrap_or_default()
        .as_millis()
    );

    // Write in files
    let mut entries_view = entries_buf.as_mut_slice();
    // Find first and last range, then iter_mut and this sub_slice
    // * unwrap() is ok for first and last since we 'break' if the buffer is empty.
    // * unwrap() is ok on binary_search since ranges cover the full hash range.
    let first_h = (hpx29(&entries_view.first().unwrap()) >> twice_dd) as u32;
    let last_h = (hpx29(&entries_view.last().unwrap()) >> twice_dd) as u32;
    let rstart = ranges_counts
      .binary_search_by(get_range_binsearch(first_h))
      .unwrap();
    let rend = ranges_counts
      .binary_search_by(get_range_binsearch(last_h))
      .unwrap();

    for (range, count) in &mut ranges_counts[rstart..=rend] {
      /*debug!(
        "Range [{}, {}]; count before: {}",
        range.start, range.end, count
      );*/
      let to = entries_view.partition_point(|row| {
        let h = (hpx29(row) >> twice_dd) as u32;
        range.contains(&h)
      });
      if to > 0 {
        let (to_be_writen, remaining) = entries_view.split_at_mut(to);
        entries_view = remaining;
        let tstart = SystemTime::now();
        let file = params.create_file_or_open_to_append(depth, range)?;
        let mut bufw = BufWriter::new(file);
        for row in to_be_writen {
          bincode.serialize_into(&mut bufw, row)?;
        }
        debug!(
          "Range file [{}, {}), {} row written in {} ms.",
          range.start,
          range.end,
          to,
          SystemTime::now()
            .duration_since(tstart)
            .unwrap_or_default()
            .as_millis()
        );

        let to = to as u32;
        if to > *count {
          // TODO: remove tmp files?
          return Err(
            format!(
              "Wrong number of counts for range [{}, {}). N remaining: {}. N to be removed: {}.",
              range.start, range.end, count, &to
            )
            .into(),
          );
        }
        *count -= to;
        /*debug!(
          "Range [{}, {}]; count after: {}",
          range.start, range.end, count
        );*/
      }
    }
    n_tot_it += n as u64;
    entries_buf.clear();
    // assert_eq!(entries_buf.len(), 0);
    if n < n_max {
      break;
    }
  }
  // Post write operation checks
  // * check ntot
  if n_tot_it != n_tot {
    // TODO: remove tmp files?
    return Err(
      format!(
        "Wrong number of written rows. Expected: {}. Actual: {}",
        n_tot, n_tot_it
      )
      .into(),
    );
  }
  // * check healpix ranges counts (all 0)
  for (range, count) in &ranges_counts {
    if *count != 0 {
      // TODO: remove tmp files?
      return Err(
        format!(
          "Wrong number of written rows for range [{}, {}). N remaining: {}",
          range.start, range.end, count
        )
        .into(),
      );
    }
  }
  // write the count map??!
  params.write_ok().map_err(|e| e.into())
}

/// See [hpx_external_sort_with_knowledge](#hpx_external_sort_with_knowledge).
/// This method perform the second part of the operation, i.e. read the temporary writen files
/// and sort them on-the-fly to provide a sorted iterator over all entries.
pub fn hpx_external_sort_with_knowledge_read_tmp<T, F>(
  hpx29: F,
  tmp_dir: PathBuf,
) -> Result<impl Iterator<Item = Result<T, Box<dyn Error>>>, Box<dyn Error>>
where
  T: ExtSortable,
  F: Fn(&T) -> u64 + Sync,
{
  // Get params
  let info = SimpleExtSortParams::read_info_gen(&tmp_dir)?;
  // debug!("Loaded info file: {:?}", &info);
  let params = SimpleExtSortParams::from_dir_and_info(tmp_dir, &info);
  // Get list of written files, ordered according to the range.
  let ordered_files = params.get_ordered_files_in_tmp_dir()?;
  // debug!("Ordered files: {:?}", &ordered_files);
  if ordered_files.is_empty() {
    return Err(format!("No tmp file found in {:?}", &params.tmp_dir).into());
  }
  // Init thread pool
  let thread_pool = get_thread_pool(params.n_threads);

  fn load_file<TT: ExtSortable>(nrows: usize, path: PathBuf) -> Result<Vec<TT>, Box<dyn Error>> {
    let path_str = path.to_string_lossy().to_string();
    debug!("Start loading file {}", path_str);
    let tstart = SystemTime::now();
    let file = File::open(path)?;
    let file_len = file.metadata()?.len();
    let mut bufr = BufReader::new(file);
    let bincode = get_bincode();
    let mut rows = Vec::with_capacity(nrows);
    for _ in 0..nrows {
      rows.push(bincode.deserialize_from(&mut bufr)?);
    }
    let pos = bufr.stream_position()?;
    if pos != file_len {
      Err(
        format!(
          "File position '{}' does not match file len '{}'.",
          pos, file_len
        )
        .into(),
      )
    } else {
      debug!(
        "Read file {} in {} ms",
        path_str,
        SystemTime::now()
          .duration_since(tstart)
          .unwrap_or_default()
          .as_millis()
      );
      Ok(rows)
    }
  }
  // load and sort and put in a MPSC ??

  // Unwrap ok since we tested that `ordered_files` is not empty.
  let mut ordered_files_counts_it = ordered_files.into_iter().zip(info.ordered_ranges_counts);
  /*let (file_path, depth, range) = ordered_files_it.next().unwrap();
  let i_nrows = info.ordered_ranges_counts.binary_search_by_key(&range, |&(range, count)| range)
    .map_err(|_| format!("Range {:?} not found in file list.", &range).into())?;*/
  let ((file_path, _depth, lrange), (rrange, nrows)) = ordered_files_counts_it.next().unwrap();
  assert_eq!(
    rrange, lrange,
    "File range difference from counts range: {:?} != {:?}.",
    &rrange, &lrange
  );

  let mut file_content = load_file(nrows as usize, file_path)?;
  thread_pool.install(|| file_content.par_sort_by_cached_key(&hpx29));

  struct GlobalIt<TT: ExtSortable, FF: Fn(&TT) -> u64 + Sync> {
    /// Pool of thread used for in-memmory sort.
    thread_pool: ThreadPool,
    /// Method to extract/compute the HEALPix order29 index from a row.
    hpx29: FF,
    /// Iterates on the files and ragen info: ((PathBuf, Depth, Range), (Range, NRows)).
    ordered_files_counts_it: Zip<IntoIter<(PathBuf, u8, Range<u32>)>, IntoIter<(Range<u32>, u32)>>,
    /// Iterates on the sorted rows of a file.
    rows_it: IntoIter<TT>,
  }
  impl<TT: ExtSortable, FF: Fn(&TT) -> u64 + Sync> Iterator for GlobalIt<TT, FF> {
    type Item = Result<TT, Box<dyn Error>>;

    fn next(&mut self) -> Option<Self::Item> {
      match self.rows_it.next() {
        Some(next_row) => Some(Ok(next_row)),
        None => match self.ordered_files_counts_it.next() {
          Some(((path, _depth, lrange), (rrange, nrows))) => {
            assert_eq!(
              rrange, lrange,
              "File range difference from counts range: {:?} != {:?}.",
              &rrange, &lrange
            );
            // Here, instead of loading the file, get it from a sender/receiver!!
            match load_file(nrows as usize, path.clone()) {
              Ok(mut file_content) => {
                let tstart = SystemTime::now();
                self
                  .thread_pool
                  .install(|| file_content.par_sort_by_cached_key(&self.hpx29));
                debug!(
                  "Loaded file {} sorted in {} ms",
                  path.to_string_lossy(),
                  SystemTime::now()
                    .duration_since(tstart)
                    .unwrap_or_default()
                    .as_millis()
                );
                self.rows_it = file_content.into_iter();
                self.next()
              }
              Err(e) => Some(Err(e)),
            }
          }
          None => None,
        },
      }
    }
  }

  Ok(GlobalIt {
    thread_pool,
    hpx29,
    ordered_files_counts_it,
    rows_it: file_content.into_iter(),
  })
}

// DefaultOptions
fn get_bincode() -> WithOtherIntEncoding<DefaultOptions, FixintEncoding> {
  DefaultOptions::new().with_fixint_encoding()
  // .allow_trailing_bytes()
}

/// Method telling if a range is smaller, contains or is greater than the given input hash value.
fn get_range_binsearch(hash: u32) -> impl Fn(&(Range<u32>, u32)) -> Ordering {
  move |(target, _)| {
    if target.end <= hash {
      Ordering::Less
    } else if hash < target.start {
      Ordering::Greater
    } else {
      Ordering::Equal
    }
  }
}

/// Returns ordered ranges at the skymap depth containing a sum of `counts` of maximum `threshold`,
/// together with this sum of counts each range contains.
/// The first value of the returned tuple is the total number of counts in the skymap.
fn skymap2ranges<'a, S>(counts: &'a S, threshold: u32) -> (u64, Vec<(Range<u32>, u32)>)
where
  S: SkyMap<'a, HashType = u32, ValueType = u32>,
{
  // 1000 chosen assuming a total count of 10x10^9 and a threshold of 10x10^6.
  let mut range_count_vec: Vec<(Range<u32>, u32)> = Vec::with_capacity(1000);
  let mut n_tot = 0_u64;

  let depth = counts.depth();

  let mut start = 0_u32;
  let mut cumul_count = 0_u32;

  for (hash, count) in counts.entries() {
    let sum = cumul_count + count;
    if sum <= threshold {
      cumul_count = sum;
    } else {
      range_count_vec.push((start..hash, cumul_count));
      n_tot += cumul_count as u64;
      start = hash;
      cumul_count = *count;
    }
  }
  range_count_vec.push((start..n_hash(depth) as u32, cumul_count));
  n_tot += cumul_count as u64;

  (n_tot, range_count_vec)
}

/// # Params
/// * `iterator`: the iterator on depth 29 healpix cells computed from the rows, and used to compute the count map
/// * `iterable`: object returning an iterator on the rows that are actually sorted and written in output
/// # Warning
/// The input `iterator` and `iterable` **must** provide the same count map!
/// # Note 1
/// Both the input and the output are row oriented, another mechanism should be considered
/// to sort a column oriented format.
/// Also, only basic row oriented formats are supported (since we serialize row by row and not
/// globally providing a serializer). It could be improved in a future, more elaborate, version.
/// # Note 2
/// For simple case, hy not use a `Clonable IntoIter`?
/// I just found [this post](https://orxfun.github.io/orxfun-notes/#/missing-iterable-traits-2024-12-13)
/// that seems to tackle a similar problem here (iterating several time over a same Collection).
pub fn hpx_external_sort<'a, T, E, I, J, F>(
  iterator: I, // replace this by another iterator on a different type of rows + hpx29: then compute in parallel (parsing + ContMap + reduce/fold countmap)
  // FARIE METHOD PARALLEL COUNT_MAP SEPAREE Et L'APPELER ICI!!
  iterable: J,
  hpx29: F,
  depth: u8,
  sort_params: Option<SimpleExtSortParams>,
) -> Result<impl Iterator<Item = Result<T, Box<dyn Error>>>, Box<dyn Error>>
where
  T: ExtSortable,
  E: Error + 'static,
  I: Iterator<Item = Result<u64, E>>,
  J: IntoIterator<Item = Result<T, E>>,
  F: Fn(&T) -> u64 + Sync,
{
  debug!("Starts computing the count map of depth {}...", depth);
  let tstart = SystemTime::now();
  let params = sort_params.unwrap_or_default();
  let twice_dd = (29 - depth) << 1;
  let count_map = CountMapU32::from_hash_values(
    depth,
    iterator
      .enumerate()
      .filter_map(|(irow, row_res)| match row_res {
        Ok(row_hash) => Some((row_hash >> twice_dd) as u32),
        Err(e) => {
          error!("Error at line {}, line ignored: {:?}", irow, e);
          None
        }
      }),
  );
  debug!(
    "... count map of depth {} computed in {} ms.",
    depth,
    SystemTime::now()
      .duration_since(tstart)
      .unwrap_or_default()
      .as_millis()
  );
  let tstart = SystemTime::now();
  params.create_tmp_dir()?;
  count_map.to_fits_file(params.create_countmap_file_path())?;
  debug!(
    "Count map of writen in {} ms.",
    SystemTime::now()
      .duration_since(tstart)
      .unwrap_or_default()
      .as_millis()
  );
  hpx_external_sort_with_knowledge(iterable.into_iter(), &count_map, hpx29, Some(params))
}

pub fn hpx_external_sort_stream<'a, T, E, I, F>(
  stream: I,
  hpx29: F,
  depth: u8,
  sort_params: Option<SimpleExtSortParams>,
) -> Result<impl Iterator<Item = Result<T, Box<dyn Error>>>, Box<dyn Error>>
where
  T: ExtSortable,
  E: Error + 'static,
  I: Iterator<Item = Result<T, E>>,
  F: Fn(&T) -> u64 + Sync,
{
  let params = sort_params.unwrap_or_default();
  let twice_dd = (29 - depth) << 1;
  let bincode = get_bincode();
  /// Write tmp file containing all rows, and compute the count map
  let count_map = {
    let tmp_file_all = params
      .create_tmp_dir()
      .and_then(|()| params.create_file_all())?;
    let mut bufw = BufWriter::new(tmp_file_all);
    CountMapU32::from_hash_values(
      depth,
      stream
        .enumerate()
        .filter_map(|(irow, row_res)| match row_res {
          Ok(row) => {
            let hash = (hpx29(&row) >> twice_dd) as u32;
            match bincode.serialize_into(&mut bufw, &row) {
              Ok(()) => Some(hash),
              Err(e) => {
                error!("Error writing line {}, line ignored: {:?}", irow, e);
                None
              }
            }
          }
          Err(e) => {
            error!("Error reading line {}, line ignored: {:?}", irow, e);
            None
          }
        }),
    )
  };
  /// Get an iterator over all written rows and apply the second part of the algorithm.
  let n_rows = count_map.values().map(|count| *count as u64).sum();
  let mut bufr = params.open_file_all().map(BufReader::new)?;
  hpx_external_sort_with_knowledge(
    (0..n_rows).map(move |_| bincode.deserialize_from(&mut bufr)),
    &count_map,
    hpx29,
    Some(params),
  )
}

/// Apply the HEALPix external sort on a CSV file.
/// # Params
/// * `input_path` path of the file to be sorted.
/// * `output_path` path of a file in which the result has to be put.
/// * `output_overwrite` do not fail if the output path already exists, but overwrite the file content
///   (WARNING: setting this flag to `true` is discouraged).
/// * `ilon`: index of the longitude column (starting at 0)
/// * `ilat`: index of the latitude column (starting at 0)
/// * `has_header`: tells that the first non-commented line is a header line (to be ignored in the
///   sort process, but to be copied in output)
/// * `separator`: ASCII value of the separator (`b','` if None)
/// * `depth`: depth used to compute the count map and for the file ranges
/// * `sort_params`: common sorting parameters
/// * `clean`: remove temporary files if the operation succeed
/// # Note
/// For small files, we advice to read the full content of the in memory and use the
/// [hpx_internal_sort](#hpx_internal_sort) method.
pub fn hpx_external_sort_csv_file<IN: AsRef<Path>, OUT: AsRef<Path>>(
  input_path: IN,
  output_path: OUT,
  output_overwrite: bool,
  ilon: usize,
  ilat: usize,
  has_header: bool,
  separator: Option<char>,
  depth: u8,
  sort_params: Option<SimpleExtSortParams>,
  clean: bool,
) -> Result<(), Box<dyn Error>> {
  let separator = separator.unwrap_or(',');
  let layer29 = get(29);
  let hpx29 = move |s: &String| {
    let cols = s.split(separator).collect::<Vec<&str>>();
    match (cols[ilon].parse::<f64>(), cols[ilat].parse::<f64>()) {
      (Ok(lon), Ok(lat)) => layer29.hash(lon.to_radians(), lat.to_radians()),
      _ => {
        error!("Error parsing coordinates at line: {}. Hash set to 0.", s);
        0
      }
    }
  };
  let mut line_res_it = BufReader::new(File::open(&input_path)?).lines().peekable();
  let mut bufw = BufWriter::new(if output_overwrite {
    OpenOptions::new()
      .write(true)
      .create(true)
      .open(output_path)
  } else {
    OpenOptions::new()
      .write(true)
      .create_new(true)
      .open(output_path)
  }?);
  // Handle starting comments
  while let Some(Ok(line)) = line_res_it.next_if(|res| {
    res
      .as_ref()
      .map(|line| line.starts_with('#'))
      .unwrap_or(false)
  }) {
    bufw.write_all(line.as_bytes())?;
  }
  // Handle header line
  if has_header {
    if let Some(header) = line_res_it.next().transpose()? {
      bufw.write_all(header.as_bytes())?;
    }
  }
  // Create an iterable from the original input files (yes, we need to read it twice)
  struct IntoIterStruct<P: AsRef<Path>> {
    input_path: P,
    has_header: bool,
  }
  impl<P: AsRef<Path>> IntoIterator for IntoIterStruct<P> {
    type Item = Result<String, IoError>;
    type IntoIter = Peekable<Lines<BufReader<File>>>;

    fn into_iter(self) -> Self::IntoIter {
      let file = match File::open(self.input_path.as_ref()) {
        Ok(file) => file,
        Err(e) => panic!(
          "Unable to re-open file {}: {:?}",
          self.input_path.as_ref().to_string_lossy(),
          e
        ),
      };
      let mut line_res_it = BufReader::new(file).lines().peekable();
      // Consume starting comments
      while line_res_it
        .next_if(|res| {
          res
            .as_ref()
            .map(|line| line.starts_with('#'))
            .unwrap_or(false)
        })
        .is_some()
      {}
      // Comsume the header line if any
      if self.has_header {
        line_res_it.next();
      }
      line_res_it
    }
  }
  // Get the sorted iterator
  let tmp_path = sort_params
    .as_ref()
    .map(|p| p.tmp_dir.clone())
    .unwrap_or(SimpleExtSortParams::default().tmp_dir);
  let mut sorted_it = hpx_external_sort(
    line_res_it.map(|row_res| row_res.map(|row| hpx29(&row))),
    IntoIterStruct {
      input_path,
      has_header,
    },
    hpx29,
    depth,
    sort_params,
  )?;
  debug!("Starts writing sorted rows in output file...");
  let tstart = SystemTime::now();
  // Write the results
  for row_res in sorted_it {
    row_res.and_then(|row| writeln!(bufw, "{}", row).map_err(|e| e.into()))?;
  }
  debug!(
    "... writing sorted rows in {} ms",
    SystemTime::now()
      .duration_since(tstart)
      .unwrap_or_default()
      .as_millis()
  );
  if clean {
    SimpleExtSortParams::clean(tmp_path).map_err(|e| e.into())
  } else {
    Ok(())
  }
}

/// Create an HEALPix ordered CSV file containing one row per HEALPix hash values at the given
/// input depth.
/// Each row contains the healpix hash value, the longitude and latitude of the center of the
/// HEALPix cell and additional dummy columns adding to 34 bytes (including the coma).
pub fn create_test_file(depth: u8, path: &str) -> Result<(), IoError> {
  debug!("Starts writing test file...");
  let tstart = SystemTime::now();
  let n_hash = n_hash(depth);
  let width = (n_hash.ilog10() + 1) as usize; // => 50_331_648 rows
  let dummmy_cols = "abc,def,ghi,jkl,mno,pqr,stu,vwx,yz"; // 34 bytes
  let layer = get(depth);

  let mut bufw = BufWriter::new(File::create(path)?);

  for h in 0..n_hash {
    let (lon, lat) = layer.center(h as u64);
    writeln!(
      bufw,
      "{:0width$},{:014.10},{:+014.10},{}",
      h,
      lon.to_degrees(),
      lat.to_degrees(),
      dummmy_cols
    )?;
  }
  debug!(
    "... done in {} ms",
    SystemTime::now()
      .duration_since(tstart)
      .unwrap_or_default()
      .as_millis()
  );
  Ok(())
}

// create_test_file("hpx_sort_result.csv")
// Use the 'shuf' linux command to shuffle the file
// Run the process
// Compare the output file with the generated file
// Compare with linux sort!!

#[cfg(test)]
mod tests {
  use super::*;
  use std::process::{Command, ExitStatus};

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

  #[cfg(target_os = "linux")]
  #[test]
  fn testok_sort() {
    init_logger();
    // Run with:
    // > cargo test testok_sort --release -- --nocapture
    // Compare with:
    //    time sort --buffer-size=70M --parallel 4 input.csv -T sort_tmp > toto.res

    let depth_file = 5; //10;
    let depth_sort = 4;
    fs::create_dir_all("./local_resources/").unwrap();
    create_test_file(depth_file, "./local_resources/test.csv").unwrap();
    Command::new("bash")
      .arg("-c")
      .arg("shuf ./local_resources/test.csv -o ./local_resources/input.csv")
      .output()
      .expect("failed to execute process");
    hpx_external_sort_csv_file(
      "./local_resources/input.csv",
      "./local_resources/output.csv",
      true,
      1,
      2,
      false,
      Some(','),
      depth_sort,
      Some(SimpleExtSortParams {
        tmp_dir: PathBuf::from("./local_resources/sort_tmp/"),
        n_elems_per_chunk: 1000, // 1_000_000
        n_threads: Some(4),
      }),
      true,
    )
    .unwrap();
    let out = Command::new("bash")
      .arg("-c")
      .arg("diff ./local_resources/test.csv ./local_resources/output.csv")
      .output()
      .expect("failed to execute process");
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
    Command::new("bash")
      .arg("-c")
      .arg("rm -r ./local_resources/test.csv ./local_resources/output.csv")
      .output()
      .expect("failed to execute process");
  }

  /*
  #[cfg(target_os = "linux")]
  #[test]
  fn testok_bigsort() {
    init_logger();
    // Run with:
    // > cargo test testok_sort --release -- --nocapture
    // Compare with:
    //   time sort --buffer-size=70M --parallel 4 input.csv -T sort_tmp > toto.res

    let depth_file = 11;
    let depth_sort = 4;
    fs::create_dir_all("/data/pineau/sandbox/").unwrap();
    create_test_file(depth_file, "/data/pineau/sandbox/test.csv").unwrap();
    Command::new("bash")
      .arg("-c")
      .arg("shuf /data/pineau/sandbox/test.csv -o /data/pineau/sandbox/input.csv")
      .output()
      .expect("failed to execute process");
    hpx_external_sort_csv_file(
      "/data/pineau/sandbox/input.csv",
      "/data/pineau/sandbox/output.csv",
      true,
      1,
      2,
      false,
      Some(','),
      depth_sort,
      Some(SimpleExtSortParams {
        tmp_dir: PathBuf::from("/data/pineau/sandbox/sort_tmp/"),
        n_elems_per_chunk: 5_000_000,
        n_threads: Some(4),
      }),
      false,
    )
    .unwrap();
    let out = Command::new("bash")
      .arg("-c")
      .arg("diff /data/pineau/sandbox/test.csv /data/pineau/sandbox/output.csv")
      .output()
      .expect("failed to execute process");
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
    Command::new("bash")
    .arg("-c")
    .arg("rm -r /data/pineau/sandbox/test.csv /data/pineau/sandbox/output.csv")
    .output()
    .expect("failed to execute process");
  }*/
}
