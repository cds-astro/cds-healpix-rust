# `cdshealpix` Change Log

--------------------------------------------------------------------------------

## 0.1.2

Released 2019-01-24.

### Changes

* allow compilation with `stable` (but `nighlty` needed to actually run tests and benches)
* put the CLI in a sub-project to remove structops dependencies

### Bug correction

* `is_flag` always set to false in BMOC 

### Added

* `to_range` method in module `cdsxhealpix::nested`
* `to_ranges` and `from_ranges` method in module `cdshealpix::nested::bmoc`
* logical operators on BMOC (`not`, `and`, `or`, `xor`)  **BUT THEY STILL HAVE TO BE TESTED**
* Fix memory leak in Python when calling `bmoc_free` (by Matthieu Baumann)

--------------------------------------------------------------------------------

## 0.1.1

Released 2019-01-18.



