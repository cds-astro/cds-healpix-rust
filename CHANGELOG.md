# `cdshealpix` Change Log



--------------------------------------------------------------------------------

## 0.1.3

Released 2019-xx-xx.

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



