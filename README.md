<meta charset="utf-8"/>

# `cdshealpix-rust`

**[CDS](http://cdsweb.u-strasbg.fr) implementation of the HEALPix tesselation in Rust and modules to generate libraries
in WebAssembly, Python, ...**

[![](https://img.shields.io/crates/v/cdshealpix.svg)](https://crates.io/crates/cdshealpix)
[![](https://img.shields.io/crates/d/cdshealpix.svg)](https://crates.io/crates/cdshealpix)
[![API Documentation on docs.rs](https://docs.rs/cdshealpix/badge.svg)](https://docs.rs/cdshealpix/)
![Rust](https://github.com/cds-astro/cds-healpix-rust/workflows/Rust/badge.svg)


About
-----

This library is an implementation in Rust of the HEALPix tesselation.
This implementation has been made by the Strasbourg astronomical Data Centre (*Centre de Données astronomique de
Strasbourg*, [CDS](http://cdsweb.u-strasbg.fr)).

It is used in:

* [Aladin Lite V3](https://github.com/cds-astro/aladin-lite)
* [The CDS HEALPix Python](https://github.com/cds-astro/cds-healpix-python) package
* [HPXCli](https://github.com/cds-astro/cds-healpix-rust/tree/master/crates/cli) for command line HEALPix related
  manipulations
* [The CDS MOC library in Rust](https://github.com/cds-astro/cds-moc-rust) used in:
    + [MOCPy](https://github.com/cds-astro/mocpy), a Python wrapper to manipulate MOCs;
    + [MOCli](https://github.com/cds-astro/cds-moc-rust/tree/main/crates/cli) a standalone command line tool to
      manipulate MOCs on linux, MacOS and Windows;
    + [MOCWasm](https://github.com/cds-astro/cds-moc-rust/tree/main/crates/wasm), a WASM library to manipulate MOCs from
      web browsers;
    + [MOCSet](https://github.com/cds-astro/cds-moc-rust/tree/main/crates/set), a standalone command line tool to build,
      update and query a persistent set of MOCs.

* CDS internal developments
* *Please help me fill in this list*

Initially, it is a port of a part of the CDS Java library
available [here](https://github.com/cds-astro/cds-healpix-java),
but improvement have been added while porting the code and new features are added.

For information on HEALPix in general, see:

* The [official web site](https://healpix.jpl.nasa.gov/)
* The [Wikipedia page](https://en.wikipedia.org/wiki/HEALPix)
* The two main reference papers: [Gorski (2005)](http://adsabs.harvard.edu/abs/2005ApJ...622..759G)
  and [Calabretta (2007)](http://adsabs.harvard.edu/abs/2007MNRAS.381..865C)

Official implementations, are available [here](https://healpix.sourceforge.io/). It contains GPL v2 codes in Fortran,
C++, Java, IDL, Python, ...

Other independent HEALPix implementations:

* [Astropy-healpix](https://github.com/astropy/astropy-healpix) python wrapper using a C code (C code by Dustin Lang,
  python wrapper by Thomas Robitaille and others)
* [Javascript/Typescript](https://github.com/michitaro/healpix) implementation by Koike Michitaro
* [Julia](https://github.com/ziotom78/Healpix.jl) implementation by Maurizio Tomasi
* [C](https://sourceforge.net/projects/healpix/files/healpix_bare_1.0/) "official" core functionalities implementation
  in BSD by Martin Reinecke
* *Please Help me adding links to other HEALPix resources and codes*

Warning
-------

For best performances on your specific hardware, you can compile using:

```bash
RUSTFLAGS='-C target-cpu=native' cargo build --release
```

This uses BMI2 instructions PDEP and PEXT, if supported by your processor, for bit interleaving.

However, the implementaion of those instructions on **AMD Ryzen processors** are **extremely slow** (20x slower than a
lookup table,
doubling the `hash` computation time)!
You can test it using:

```bash
RUSTFLAGS='-C target-cpu=native' cargo bench
```

If the result of `ZOrderCurve/BMI` is slower thatn `ZOrderCurve/LUPT`, compile without the `native` support:

```bash
cargo build --release
```

Target 32 bit on a 64 bit linux
-------------------------------

```rust
rustup target install i686-unknown-linux-gnu
sudo apt-get install gcc-multilib
RUSTFLAGS=' - C target-cpu=native' cargo build - - target=i686-unknown-linux-gnu - - release
```

Features
--------

* Supports the **HEALix Nested scheme**
    + Supports approximated `cone` and `elliptical cone` coverage plus **exact** `polygon` coverage queries
    + Supports `BMOC` (MOC with a flag telling if a cell is fully or partially covered by a surface) as a result of
      `cone`, `polygon` ot `elliptical cone` coverage queries
    + Supports logical operations on `BMOCs` and `BMOC` creation from a list of cell number at a given depth
    + Supports implicit HEALPix density maps (single column so far), with PNG image creation for density maps
    + Supports (non standard) HEALPix multi-order maps (MOM, single columns so far), with PNG image creation for density
      MOMs
    + Supports HEALPix external sort
    + Supports sorted indexation
* Supports the **HEALPix Ring scheme** with **any** NSIDE (i.e. not necessarilly powers of 2)

Missing Features
----------------

* Not supported
    * Polygon and ellipse in the RING scheme
    * Spherical Harmonics computations
    * *Help me fill this*
* Not yet implemented
    * Exact cone and ellipse solution (but using the `custom` approx methods, one can handle the rate of false
      positives)
    * Cone query in the RING scheme

Examples
--------

Compute the cell number of a given position on the unit-sphere at a given HEALPix depth.

```rust
use cdshealpix::nside;
use cdshealpix::nested::{get, Layer};


let depth = 12_u8;
let lon = 12.5_f64.to_radians();
let lat = 89.99999_f64.to_radians();

let nested_d12 = get(depth);
let nside = nside(depth) as u64;
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
use cdshealpix::nested::{get, Layer};

let depth = 6_u8;
let nested_d6 = get(depth);

let lon = 13.158329_f64.to_radians();
let lat = - 72.80028_f64.to_radians();
let radius = 5.64323_f64.to_radians();

let moc = nested_d6.cone_overlap_approx(lon, lat, radius);
for cell in moc.into_iter() {
println ! ("cell: {:?}", cell);
}
```

Standalone
----------

See [HPXCli](https://github.com/cds-astro/cds-healpix-rust/tree/master/crates/cli).


WebAssembly
-----------

(Not really maintained so far: if you need it, please let us know!)  
To build and use the WebAssembly (and Javascript) files, the `libwasmbingen` directory.
We rely on [wasm-bingen](https://github.com/rustwasm/wasm-bindgen).


Python
------

See [cdshealpix python](https://github.com/cds-astro/cds-healpix-python/) available
on [pypi](https://pypi.org/project/cdshealpix/):

```bash
pip install -U cdshealpix
```

ToDo list
---------

* [ ] Modify elliptical cone: compute distance to both foci
* [ ] Implement the exact cone solution

Acknowledgements
----------------

If you use this code and work in a scientific public domain
(especially astronomy), please acknowledge its usage and the
[CDS](https://en.wikipedia.org/wiki/Centre_de_donn%C3%A9es_astronomiques_de_Strasbourg)
who developed it.
It may help us in promoting our work to our financiers.


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

This code started has my first Rust code.

