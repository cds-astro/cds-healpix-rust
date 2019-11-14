
<meta charset="utf-8"/>

# `cdshealpix-rust`

**[CDS](http://cdsweb.u-strasbg.fr) implementation of the HEALPix tesselation in Rust and modules to generate libraries in WebAssembly, Python, ...**

[![](https://meritbadge.herokuapp.com/cdshealpix)](https://crates.io/crates/cdshealpix)
[![](https://img.shields.io/crates/d/cdshealpix.svg)](https://crates.io/crates/cdshealpix)
[![API Documentation on docs.rs](https://docs.rs/cdshealpix/badge.svg)](https://docs.rs/cdshealpix/)

About
-----

This library is an implementation in Rust of the HEALPix tesselation.
This implementation has been made by the Strasbourg astronomical Data Centre (*Centre de Données astronomique de Strasbourg*, [CDS](http://cdsweb.u-strasbg.fr)).

Initially, it is a port of a part of the CDS Java library available [here](https://github.com/cds-astro/cds-healpix-java),
but improvement have been added while porting the code.

For informations on HEALPix in general, see:
 * The [official web site](https://healpix.jpl.nasa.gov/)
 * The [Wikipedia page](https://en.wikipedia.org/wiki/HEALPix)
 * The two main reference papers: [Gorski (2005)](http://adsabs.harvard.edu/abs/2005ApJ...622..759G) and [Calabretta (2007)](http://adsabs.harvard.edu/abs/2007MNRAS.381..865C)

Official implementations, are available [here](https://healpix.sourceforge.io/). It contains GPL v2 codes in Fortran, C++, Java, IDL, Python, ...

Other independant HEALPix implementations:
 * [Astropy-healpix](https://github.com/astropy/astropy-healpix) python wrapper using a C code (C code by Dustin Lang, python wrapper by Thomas Robitaille and others)
 * [Javascript/Typescript](https://github.com/michitaro/healpix) implementation by Koike Michitaro
 * [Julia](https://github.com/ziotom78/Healpix.jl) implementation by Maurizio Tomasi
 * [C](https://sourceforge.net/projects/healpix/files/healpix_bare_1.0/) "official" core functionalities implementation in BSD by Martin Reinecke
 * ... (Help me to add links to other HEALPix resources and codes).

Warning
-------

For best performances on your specific hardware, you can compile using:
```bash
RUSTFLAGS='-C target-cpu=native' cargo build --release
```
This uses BMI2 instructions PDEP and PEXT, if supported by your processor, for bit interleaving.

However, the implementaion of those instructions on **AMD Ryzen processors** are **extremely slow** (20x slower than a lookup table, 
doubling the `hash` computation time)! 
You can test it usingi:
```bash
RUSTFLAGS='-C target-cpu=native' cargo bench
```
If the result of `ZOrderCurve/BMI` is slower thatn `ZOrderCurve/LUPT`, compile without the `native` support:
```bash
cargo build --release
```

Features
--------

 * Supports the **HEALix Nested scheme**
     + Supports approximated `cone` and `elliptical cone` coverage plus **exact** `polygon` coverage queries
     + Supports `BMOC` (MOC with a flag telling if a cell is fully or partially covered by a surface) as a result of `cone`, `polygon` ot `elliptical cone` coverage queries
     + Supports logical operations on `BMOCs` and `BMOC` creation from a list of cell number at a given depth
 * Supports the **HEALPix Ring scheme** with **any** NSIDE (i.e. not necessarilly powers of 2)

Missing Features
----------------

 * Not supported
   * Polygon and ellipse in the RING scheme
   * Spherical Harmonics computations
   * (Help me fill this)
 * Not yet implemented
   * Exact cone and ellipse solution (but using the `custom` approx methods, one can handle the rate of false positives)  
   * Cone query in the RING scheme

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
let cell_number= 10_u64;

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
 
let moc = nested_d6.cone_overlap_approx(lon, lat, radius);
for cell in moc.into_iter() {
    println!("cell: {:?}", cell);
}
```

Standalone
----------

(Not on crates.io, but on github) 
The code source of the very beginning of a standalone exec can be found in `cli/src/bin.rs`.

WebAssembly
-----------

(Not on crates.io, but on github) 
To build and use the WebAssembly (and Javascript) files, the `libwasmbingen` directory.
We rely on [wasm-bingen](https://github.com/rustwasm/wasm-bindgen).


Python
------

(Not on crates.io, but on github) 
See the `libpython` directory containing a very first integration in python  using [CFFI](https://cffi.readthedocs.io/en/latest/).

For a clean Python wrapper and associated Wheels, see Matthieu Baumann's project [cds-healpix-python](https://github.com/cds-astro/cds-healpix-python/).
To use the library in python, install it through `pip` (examples are provided on github [cds-healpix-python](https://github.com/cds-astro/cds-healpix-python/)):
```bash
pip install cdshealpix
```

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

It a first code in Rust, feel free to give some advice/feedback.

