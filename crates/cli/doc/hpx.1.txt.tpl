moc(1)
=====

Name
----
hpx - HEALPix command line tool


Synopsis
--------

*hpx* _SUBCMD_ _SUBCMDPARAMS_

*hpx* *--version*

*hpx* *--help*

*hpx* _SUBCMD_ *--help*

*hpx* _SUBCMD_ _SUBCMDPARAMS_ *--help*

*command* | *hpx* _SUBCMD_  _SUBCMDPARAM_


SUBCMD
------

_proj_::
  Computes the projected coordinates (x, y) ∊ ([0..8[, [-2..2]) of input equatorial coordinates

_unproj_::
  Computes the equatorial coordinates of Healpix projected coordinates (x, y) ∊ ([0..8[, [-2..2])
  
_nested_::
  Operations in the HEALPix NESTED scheme

_sort_::
  Sorts a CSV file by order 29 HEALPix NESTED indices, uses external sort to support huge files
  
_hcidx_::
  Create an index on an HEALPix NESTED sorted CSV file to then quickly retrieve rows in a given HEALPix cell

_qhcidx_::
  Query an HEALPix sorted and indexed CSV file (see the 'hcidx' command)

_map_::
  Create and manipulate HEALPix count and density maps

_mom_::
  Create and manipulate HEALPix count and density maps

_cov_::
  Compute the list of NESTED Healpix cells in a given sky region (see also moc-cli)


Examples
--------

hpx proj value 0.0 0.0

hpx unproj value 0.0  0.0

echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' | hpx nested center 0 list -d , | hpx proj csv -d , --paste

hpx nested table

hpx nested bestdepth --arcsec 1.0

hpx nested hash 0 value 0.0 +1.0

echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' | hpx nested center 0 list -d , | hpx nested hash --with-dxdy --zuniq 0 list -d ,

hpx nested bilinear 3 value 123.456789 +12.3456789

hpx nested neighbours 0 value 5

hpx cov 9 cone --delta-depth 3 123.456789 -12.3456789 0.2 \
  | egrep -v "^#" \
  | hpx nested touniq --zuniq singledepth 9 csv -d , -f 2 --header --paste --paste-header "zuniq"

hpx sort --header -d , -o hipmain.sorted.csv  hipmain.csv

hpx hcidx --header -d , --depth 5 -o hipmain.sorted.hci.fits hipmain.sorted.csv

hpx qhcidx hipmain.sorted.hci.fits cell 0 1

echo "Creates a first density maps from a generated sequence of points..."
for n in $(seq 0 47); do for i in $(seq 0 ${n}); do echo "$n"; done; done \
  | hpx nested center 1 list \
  | hpx map count 2 count.d2.l.fits list
echo "... and view it."
hpx map view --silent --color-map-fn linear count.d2.l.fits count.d2.l.png allsky 200

echo "Add the two maps..."
hpx map op add count.d2.l.fits count.d2.r.fits count.d2.m.fits
echo "... and view the result."
hpx map view --silent --color-map-fn linear count.d2.m.fits count.d2.m.png allsky 200

echo "Build a MOM from the density map..."
hpx map convert --parallel 2 hip.d5.dens.fits hip.d5.dens.mom.fits dens2chi2mom
echo "... and view it."
hpx mom view --silent --color-map-fn linear hip.d5.dens.mom.fits hip.d5.dens.mom.png allsky 200
echo "Convert the MOM into in a BINTABLE MOM (check you can read it in, e.g.,  TOPCAT)"
hpx mom convert hip.d5.dens.mom.fits hip.d5.dens.mom.bintable.fits bintable

hpx cov 10 cone --delta-depth 3 123.456789 +12.3456789 0.2
hpx cov $(($(hpx nested bestdepth 0.2)+3)) cone --delta-depth 3 123.456789 +12.3456789 0.2

hpx cov 4 zone 0.0 -10.0 20.0 10.0

hpx cov 5 polygon 10.0,-10.0,10.0,+20.0,20.0,+20.0,20.0,-10.0

hpx cov 9 stcs "Circle ICRS 147.6 69.9 0.4"

hpx cov 13 stcs "Difference ICRS (
  Polygon 272.536719 -19.461249 272.542612 -19.476380 272.537389 -19.491509 272.540192 -19.499823
    272.535455 -19.505218 272.528024 -19.505216 272.523437 -19.500298 272.514082 -19.503376
    272.502271 -19.500966 272.488647 -19.490390  272.481932 -19.490913 272.476737 -19.486589
    272.487633 -19.455645 272.500386 -19.444996 272.503003 -19.437557 272.512303 -19.432436
    272.514132 -19.423973 272.522103 -19.421523 272.524511 -19.413250 272.541021 -19.400024
    272.566264 -19.397500 272.564202 -19.389111 272.569055 -19.383210 272.588186 -19.386539
    272.593376 -19.381832 272.596327 -19.370541 272.624911 -19.358915 272.629256 -19.347842
    272.642277 -19.341020 272.651322 -19.330424 272.653174 -19.325079 272.648903 -19.313708
    272.639616 -19.311098 272.638128 -19.303083 272.632705 -19.299839 272.627971 -19.289408
    272.628226 -19.276293 272.633750 -19.270590 272.615109 -19.241810 272.614704 -19.221196
    272.618224 -19.215572 272.630809 -19.209945 272.633540 -19.198681 272.640711 -19.195292
    272.643028 -19.186751 272.651477 -19.182729 272.649821 -19.174859 272.656782 -19.169272
    272.658933 -19.161883 272.678012 -19.159481 272.689173 -19.176982 272.689395 -19.183512
    272.678006 -19.204016 272.671112 -19.206598 272.664854 -19.203523 272.662760 -19.211156
    272.654435 -19.214434 272.652969 -19.222085 272.656724 -19.242136 272.650071 -19.265092
    272.652868 -19.274296 272.660871 -19.249462 272.670041 -19.247807 272.675533 -19.254935
    272.673291 -19.273917 272.668710 -19.279245 272.671460 -19.287043 272.667507 -19.293933
    272.669261 -19.300601 272.663969 -19.307130 272.672626 -19.308954 272.675225 -19.316490
    272.657188 -19.349105 272.657638 -19.367455 272.662447 -19.372035 272.662232 -19.378566
    272.652479 -19.386871 272.645819 -19.387933 272.642279 -19.398277 272.629282 -19.402739
    272.621487 -19.398197 272.611782 -19.405716 272.603367 -19.404667 272.586162 -19.422703
    272.561792 -19.420008 272.555815 -19.413012 272.546500 -19.415611 272.537427 -19.424213
    272.533081 -19.441402
  Union (
    Polygon 272.511081 -19.487278 272.515300 -19.486595 272.517029 -19.471442
            272.511714 -19.458837 272.506430 -19.459001 272.496401 -19.474322 272.504821 -19.484924
    Polygon 272.630446 -19.234210 272.637274 -19.248542 272.638942 -19.231476 272.630868 -19.226364
  )
)"


DESCRIPTION
-----------

HEALPix manipulation on the command line.


VERSION
-------
{VERSION}


HOMEPAGE
--------
https://github.com/cds-astro/cds-healpix-rust

Please report bugs and feature requests in the issue tracker.


AUTHORS
-------
F.-X. Pineau <francois-xavier.pineau@astro.unistra.fr>


