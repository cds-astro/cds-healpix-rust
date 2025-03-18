use rayon::ThreadPool;

pub mod coverage;
pub mod hcidx;
pub mod input;
pub mod map;
pub mod mom;
pub mod nested;
pub mod qhcidx;
pub mod sort;
pub mod view;

/// Returns a pool of thread containing:
/// * the maximum number of threads available on the machine, if `n_threads` is `None`
/// * the given number of threads in the `n_threads` option if it is `Some`.
pub fn get_thread_pool(n_threads: Option<usize>) -> ThreadPool {
  let mut pool_builder = rayon::ThreadPoolBuilder::new();
  if let Some(n_threads) = n_threads {
    pool_builder = pool_builder.num_threads(n_threads);
  }
  pool_builder.build().unwrap()
}
