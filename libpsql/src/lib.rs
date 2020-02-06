extern crate libc;

//We could have used:
// - https://github.com/bluejekyll/pg-extend-rs
//   - https://bluejekyll.github.io/blog/rust/2018/12/27/announcing-pg-extend.html
// - https://github.com/jeff-davis/postgres-extension.rs
// See PostgreSQL code here: https://github.com/postgres/postgres/

// use libc::int8_t;
// use libc::int16_t;
// use libc::int32_t;
// use libc::int64_t;

use libc::uint8_t;
// use libc::uint16_t;
use libc::uint32_t;
use libc::uint64_t;

// use libc::c_float;
use libc::c_double;

extern crate cdshealpix;
//use cdshealpix;

#[no_mangle]
pub extern fn nside(depth: uint8_t) -> uint32_t {
    cdshealpix::nside(depth)
}

#[no_mangle]
pub extern fn hash(depth: uint8_t, lon_deg: c_double, lat_deg: c_double) -> uint64_t {
    cdshealpix::nested::hash(depth, lon_deg.to_radians(), lat_deg.to_radians())
}
