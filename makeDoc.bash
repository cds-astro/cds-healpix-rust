#!/bin/bash

from="resources/4doc"
dest="target/doc/cdshealpix"

#RUSTFLAGS='-C target-feature=+bmi,+bmi2' cargo doc
cargo doc

cp ${from}/hpx_proj.png ${dest}/hpx_proj.png
cp ${from}/d_center_vertex.png ${dest}/d_center_vertex.png

#RUSTFLAGS='-C target-feature=+bmi,+bmi2' cargo doc --open
cargo doc --open

