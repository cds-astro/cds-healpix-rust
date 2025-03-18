extern crate libc;

//We could have used:
// - https://github.com/bluejekyll/pg-extend-rs
//   - https://bluejekyll.github.io/blog/rust/2018/12/27/announcing-pg-extend.html
// - https://github.com/jeff-davis/postgres-extension.rs
// See PostgreSQL code here: https://github.com/postgres/postgres/

// use libc::c_float;
use libc::c_double;

extern crate cdshealpix;
//use cdshealpix;

#[no_mangle]
pub extern fn nside(depth: u8) -> u32 {
  cdshealpix::nside(depth)
}

#[no_mangle]
pub extern fn nest_hash(depth: u8, lon_deg: c_double, lat_deg: c_double) -> u64 {
  cdshealpix::nested::hash(depth, lon_deg.to_radians(), lat_deg.to_radians())
}

#[no_mangle]
pub extern fn nest_center(depth: u8, icell: u64, res_ptr: *mut c_double) {
  assert!(!res_ptr.is_null(), "Null pointer in center()");
  let res: &mut[c_double] = unsafe { std::slice::from_raw_parts_mut(res_ptr, 2) }; 
  let (lon, lat) = cdshealpix::nested::center(depth, icell);
  res[0] = lon.to_degrees();
  res[1] = lat.to_degrees();
}





