#!/bin/bash

from="resources/4doc"
dest="target/doc/cdshealpix"

#RUSTFLAGS='-C target-feature=+bmi,+bmi2' cargo doc
cargo doc

cp ${from}/hpx_proj.png ${dest}/hpx_proj.png
cp ${from}/d_center_vertex.png ${dest}/d_center_vertex.png
cp ${from}/external_edge.png ${dest}/nested/external_edge.png

#RUSTFLAGS='-C target-feature=+bmi,+bmi2' cargo doc --open
RUSTDOCFLAGS="--html-in-header katex.html" cargo doc --no-deps --open

