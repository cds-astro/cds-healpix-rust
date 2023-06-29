# `cdshealpix` Change Log

## 0.6.6

Released 2023-06-29

### Add

* Method `is_empty` and traits `Copy, Clone, Eq, ...` on `OrdinalSet` and `CardinalSet`

--------------------------------------------------------------------------------



## 0.6.5

Released 2023-06-28

### Add

* Add `OrdinalSet` in module `compass_point`
* Update dependencies versions
    + update base64 code

### Change

* Use cargo fmt with option `--config tab_spaces=2` to reformat the code

### Bug correction

* fix numerical issue in 'intersect_small_circle' when great circle arc is small
  (<1 arcsec) by resorting to locally flat approximation


--------------------------------------------------------------------------------


## 0.6.4

Released 2022-10-17

### Add

* Add latitude check in `bilinear_interpolation` 
* Add computation of lat=cte small circle intersection with polygon edges (great circle arcs)
* Add Sph geom polygon/parallel intersection
* Update dependencies
* Remove a few warnings in test and bench
* Make clippy happier

--------------------------------------------------------------------------------



## 0.6.3

Released 2022-04-11

### Add

* Add `into_flat_iter` in BMOC (to be able to flatten an iterator of BMOC `flat_iter`) 
* Remove a few warnings

--------------------------------------------------------------------------------


## 0.6.2

Released 2022-04-01

### Bug correction

* In the computation of constants in `ConstantsC2V`, 
  impacting the `largest_center_to_vertex_distance` family methods: 
  wrong contant 4/pi instead of pi/4 leading to overestimating the 
  `largest_center_to_vertex_distance`

--------------------------------------------------------------------------------


## 0.6.1

Released 2022-03-22

### Bug correction

* Method `box_coverage`: relax a too restrictive constraint preventing a == b

--------------------------------------------------------------------------------


## 0.6.0

Released 2022-01-25

### Add

* Method `cone_coverage_fullin`: returns HEALPix nested cells fully covered by the given cone
* Method `cone_coverage_centers`: returns HEALPix nested cells having their center in the given cone
* Method `ring_coverage_approx`: returns HEALPix nested cells overlapped by a ring
* Method `ring_coverage_approx_custom`: returns HEALPix nested cells overlapped by a ring (with a custom delta depth for better precision)

--------------------------------------------------------------------------------


## 0.5.5

Released 2020-08-05

### Add

* Add `box_coverage`

### Bug correction

* In polygon: contrary to the Java code (which was correct) a cell was added to the 
  moc if its 4 vertices where in the polygon (without testing possible intersection with the polygon).
  This may lead to wrong coverages in the case of self-intersecting polygons
  (see `nested::testok_polygone_exact_fxp`).
* Fix benches  

--------------------------------------------------------------------------------


## 0.5.4

Released 2020-11-02

### Modified

* Default behaviour for polygons

### Add

* Add `zone_coverage`.
* Add `custom_polygon_coverage` (to be able to select the default area or its complementary)

--------------------------------------------------------------------------------


## 0.5.3

Released 2020-06-25.

### Changed

* Remove dead code (not all dead code, but large parts).
* Fix compilation warnings
* Fix clippy warnings

--------------------------------------------------------------------------------


## 0.5.2

Released 2020-06-22.

### Bug correction

* When Newton-Raphson method fails when looking for "special points" in polygons 
  (due to divergence or too slow convergence) the method failed.
  Now returns None (i.e. no special point found).

--------------------------------------------------------------------------------


## 0.5.1

Released 2020-05-25.

### Bug correction

* When polygons are very large, previous code sometimes returned the complementary
  polygon. We now decided that for large polygons the gravity center should be inside
  the polygon.


--------------------------------------------------------------------------------


## 0.5.0

Released 2020-03-04.

### Added

* Add a MOC module for testing purpose
    + Add MOC compression/decompression iterators
    + Add logical MOC operations (not, and, or) taking iterators and returning iterators
    + Add Cells to Range and Range to Cells conversion (iterator based too)
* Add the method `test_coo` in `BMOC`

--------------------------------------------------------------------------------


## 0.4.1

Released 2020-02-11.

### Changed

* Transform methods in `const`, add `dyn` for trait objects.
* Add an enum implementing the `ZOrderCurve` trait to allow...
* ... the usage of a `const fn` to compute `Layers` structs at compile time.
* Replace a few constants by the ones defined in the Rust standard library
 
### Added

* `to_uniq`, `from_uniq`, `to_uniq_ivo` and `from_uniq_ivoa` to handle 
  uniq hash notation (i.e. uniq value for all possible (depth, hash) tuples).
* Add BMOC lossy compression/decompression


### Bug correction

* Fix `polygon_coverage` bug due to a bug in `great_circle_arcs_are_overlapping_in_lon`
  when a great circle crosses the RA=0 meridian.


--------------------------------------------------------------------------------


## 0.4.0 

Released 2019-11-14.

### Changed

* Do not requires `nighlty` any more (replace built-in bench by Criterion)
* Change `hash` internals (improve performance and fix specific cases)

### Bug correction

* Polygon specific case (by fixing the hash method)
* `hash`: fix very specific cases (lon = n * PI/2) in the South polar cap

--------------------------------------------------------------------------------


## 0.3.2 

Released 2019-09-05.

### Added

* `hash_with_dxdy` algo to be ported in glsl (WebGL) for AladinLite: is takes
  in input a vector (x, y, z) made of 3 single precision floats.
  Most graphics cards being limited to single precision, the method is able to
  provided indices at a maximum depth of 14. At depth 13, the precision
  on dx and dy is better than 1/512 (<=> images at order 13 + 9).

### Bug correction

* Fix NESTED `bilinear_interpolation`: ~1/4 of cases in which 2 out of 4 HEALPix
 cell numbers where swapped

--------------------------------------------------------------------------------


## 0.3.1 

Released 2019-08-19.

### Bug correction

* Fix NESTED `hash` rare panic due to numerical inaccuracies

--------------------------------------------------------------------------------


## 0.3.0 

Released 2019-06-25.

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
   + `hash`: compute the RING number from (lon, lat) coordinates
   + `hash_with_dxdy`: compute the RING number from (lon, lat) coordinates and additionally provide the offsets (dx, dy) in the cell
   + `center`: get the coordiates of the center of a RING cell
   + `center_of_projected_cell`: like center, but coordinates are gicen in the Euclidean projection plane
   + `sph_coo`: get a coordinate on the sphere from a cell and a position (dx, dy) in the cell
   + `vertices`: provide the 4 vertices of a given cell 
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


