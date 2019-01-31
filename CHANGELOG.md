# `cdshealpix` Change Log



--------------------------------------------------------------------------------

## 0.1.4

Released 2019-01-31.

### Bug correction

* method `MOCBuilderUnsafe::to_lower_depth` in module `cdshealpix::nested::bmoc`: first cell was ignored when not part of a larger cell.

### Added

* first version of the exacy polygon algorithm
* `flat iterator returning cells (containing flags) in  BMOC
* posibility to add LaTeX formula in the documentation (using the crate [katex-doc](https://crates.io/crates/katex-doc)).

### API changes

* `(polygon|cone)_overlap.*` changed to `(polygon|cone)_coverage.*`
* `polygon_coverage_approx` changed to `polygon_coverage` with a boolean to select approx / exact solution

### WARNING

* Exact polygon algo still to be tested, probably bugged so far!!

--------------------------------------------------------------------------------

## 0.1.3

Released 2019-01-27.

### Bug correction

* fix BMOC logical operations (`not`, `and`, `or`, `xor`)

### Added

* `BMOCBuilderFixedDepth` to build a BMOC from a list of cells at the MOC depth
* add sub-direcotry containing tests on MOCs (in a sub-module not to add the
  `serde` and `serde-json` dependencies in the man crate) 

### WARNING

* BMOC logical operators still to be tested with `is_full` flags possibly set to `false`

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



