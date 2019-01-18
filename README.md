
<meta charset="utf-8"/>

# `cdshealpix-rust`

**[CDS](http://cdsweb.u-strasbg.fr) implementation of the HEALPix tesselation in Rust and modules to generate libraries in WebAssemnly, Python, ...**

[![](https://meritbadge.herokuapp.com/cdshealpix)](https://crates.io/crates/cdshealpix)
[![](https://img.shields.io/crates/d/cdshealpix.svg)](https://crates.io/crates/cdshealpix)
[![API Documentation on docs.rs](https://docs.rs/cdshealpix/badge.svg)](https://docs.rs/cdshealpix/0.1.1/cdshealpix/)

About
-----

For informations on HEALPix in general, see:
 * The [official web site](https://healpix.jpl.nasa.gov/)
 * The [Wikipedia page](https://en.wikipedia.org/wiki/HEALPix)
 * The two main related papers: [Gorsky (2005)](http://adsabs.harvard.edu/abs/2005ApJ...622..759G) and [Calabretta (2007)](http://adsabs.harvard.edu/abs/2007MNRAS.381..865C)

See also the [official page](https://healpix.sourceforge.io/) containing GPL v2 codes in Fortran, C++, Java, IDL, Python, ...

(Help me to add links to other HEALPix resources and codes).

This library is mainly a port of a part of the CDS Java library available [here](https://github.com/cds-astro/cds-healpix-java). 

Features
--------

 * Supports the **HEALix Nested scheme**
 * Supports approximated `cone` and `polygon` queries
 * Supports `BMOC` (MOC with a flag telling if a cell is fully or partially covered by a surface) as a result of `cone` and `polygon` queries

Missing Features
----------------

 * Not supported
   * RING scheme
   * Spherical Harmonics computations
   * (Help me fill this)
 * Not yet implemented
   * Exact cone and polygon solution
   * Logical operations on MOC/BMOC
   * MOC creation from arbritrary input

Examples
--------

Compute the cell number of a given position on the unit-sphere at a given HEALPix depth.

```rust
use cdshealpix::{nside};
use cdshealpix::nested::{get_or_create, Layer};


let depth = 12_u8;
let lon = 12.5_f64.to_radians();
let lat = 89.99999_f64.to_radians();

let nested_d12 = get_or_create(depth);
let nside = nside(depth) as u_64;
let expected_cell_number = nside * nside - 1

assert_eq!(expected_cell_number, nested_d12.hash(lon, lat));
```

Get the spherical coorinates of the 4 vertices of a given cell at a given depth:

```rust
use cdshealpix::nested::{get_or_create, Layer};

let depth = 12_u8;
let celll_number= 10_u64;

let nested_d12 = get_or_create(depth);

let [
  (lon_south, lat_south), 
  (lon_east,  lat_east), 
  (lon_north, lat_north), 
  (lon_west,  lat_west)
] = nested_d12.vertices(cell_number);

```

Get a hierarchical view (a [MOC](http://www.ivoa.net/documents/MOC/)) on the cells overlapped by a given cone:

```rust
use cdshealpix::nested::{get_or_create, Layer};

let depth = 6_u8;
let nested_d6 = get_or_create(depth);

let lon = 13.158329_f64.to_radians();
let lat = -72.80028_f64.to_radians();
let radius = 5.64323_f64.to_radians();
 
let moc = nested3.cone_overlap_approx(lon, lat, radius);
```

Standalone
----------

The code source of the very beginning of a standalone exec can be found in `src/bin.rs`.

WebAssembly
-----------

To build and use the WebAssembly (and Javascript) files, the `libwasmbingen` directory.
We rely on [wasm-bingen](https://github.com/rustwasm/wasm-bindgen).


Python
------

See the `libpython` directory containing so far a very basic python script showing how to use the library from python.
We are currently working on making a clean Python wrapper and generating Python Wheels for a simple install through `pip`. 


License
-------

Like most projects in Rust, this project is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or
   http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or
   http://opensource.org/licenses/MIT)

at your option.


Contribution
------------

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this project by you, as defined in the Apache-2.0 license,
shall be dual licensed as above, without any additional terms or conditions.


Disclaimer
----------

This is my very first project in Rust so please be indulgent and find better examples if you want to learn Rust.

