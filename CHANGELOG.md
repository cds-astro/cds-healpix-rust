# `cdshealpix` Change Log

## 0.3.0 

Released 2019-06-05.

## Added

* Added in the NESTED scheme
    + `sph_coo`: spherical coordinates from an HEALPix cell + offsets `(dx, dy)`
    + `hash_with_dxdy`: hash value together with the `(dx, dy)` offsets
    + `path_along_cell_side`: provides points along a cell side
    + `path_along_cell_edge`: provides points along a cell edge
    + `bilinear_interpolation`: bilinear interpolation on the nested scheme
    + set `center_of_projected_cell` public
* Support of ring in NESTED:
    + `to_ring`: convert a NESTED index into a RING index
    + `from_ring`: convert a RING index into a NESTED index
* Starts supporting the RING scheme
   + `hash`: compute the RING number from (lon, lat) coordiantes
   + `center`: get the coordiates of the center of a RING cell
   + `center_of_projected_cell`: like center, but coordinates are gicen in the Euclidean projection plane
* All
   + `base_cell_from_proj_coo`: experiment to be tested 

### Bug correction

* Fix polygon potential bug (see method `is_in_lon_range`)

--------------------------------------------------------------------------------


## 0.2.0

Released 2019-05-10.

### Info

* We have started to run astropy-healpix tests on cdshealpix in the python wrapper project

### Added

* Methods to get internal/external edges of a cell a deeper depth

### Bug correction

* Fix elliptical cone
* Fix numerical precision of hash (cell number from coordinates) in NE border of North Polar Cap base cells 
* Fix numerical precision of angular distances computation near from PI

--------------------------------------------------------------------------------


## 0.1.6

Released 2019-03-14.

### Bug correction

* Elliptical cone: fix SIN projection conditions
* Elliptical cone: better handle large ellipses (when ellipse semi-major axis + bounding cone radius > PI/2)
* WARNING: still bugged in crates.io, but fixed on github!

--------------------------------------------------------------------------------


## 0.1.5

Released 2019-03-06.


### Bug correction

* Now ensures that lon in [0, 2pi[ and lat in [-pi/2, +pi/2] in Coo3D

### Added

* add support to elliptical cones

--------------------------------------------------------------------------------


## 0.1.4

Released 2019-04-01.

### Bug correction

* fix documentation error with katex on doc.rs
* first debug of the exact solution of `query_polygon`

## WARNING

* BMOC logical operators still to be tested with `is_full` flags possibly set to `false`
* More test are needed on the exact polygon algo

--------------------------------------------------------------------------------

## 0.1.3

Released 2019-01-31.

### Bug correction

* fix BMOC logical operations (`not`, `and`, `or`, `xor`)
* fix method `MOCBuilderUnsafe::to_lower_depth` in module `cdshealpix::nested::bmoc`: first cell was ignored when not part of a larger cell.

### Added

* first version of the exacy polygon algorithm
* `flat iterator returning cells (containing flags) in  BMOC
* `BMOCBuilderFixedDepth` to build a BMOC from a list of cells at the MOC depth
* posibility to add LaTeX formula in the documentation (using the crate [katex-doc](https://crates.io/crates/katex-doc)).
* add sub-direcotry containing tests on MOCs (in a sub-module not to add the
  `serde` and `serde-json` dependencies in the man crate) 
* add sub-directory containing an example of usage of cdshealpix-python

### API changes

* `(polygon|cone)_overlap.*` changed to `(polygon|cone)_coverage.*`
* `polygon_coverage_approx` changed to `polygon_coverage` with a boolean to select approx / exact solution

### WARNING

* BMOC logical operators still to be tested with `is_full` flags possibly set to `false`
* Exact polygon algo still to be tested, probably bugged so far!! (I know, creating branches would be cleaner)

--------------------------------------------------------------------------------

## 0.1.2

Released 2019-01-24.

### Changes

* allow compilation with `stable` (but `nighlty` needed to actually run tests and benches)
* put the CLI in a sub-project to remove structops dependencies

### Bug correction

* `is_flag` always set to false in BMOC now fixed

### Added

* `to_range` method in module `cdsxhealpix::nested`
* `to_ranges` and `from_ranges` method in module `cdshealpix::nested::bmoc`
* logical operators on BMOC (`not`, `and`, `or`, `xor`)  **BUT THEY STILL HAVE TO BE TESTED**
* Fix memory leak in Python when calling `bmoc_free` (by Matthieu Baumann)

--------------------------------------------------------------------------------

## 0.1.1

Released 2019-01-18.


