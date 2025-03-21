#!/bin/bash


#1:expected filename; 2: actual filename
assert_eq(){
  local diff=$(diff $1 $2)
  if [[ ${diff} != "" ]]; then
    echo "KO!"
    echo "$1 != $2; diff (<: expected, >: actual):"
    echo "${diff}"
    exit -1;
  fi
}

# Stands for 'generic test'
# 1: name of the method to be called
# 2: description of the method
gtest(){
  local name="$1"
  local expect="${name}.expect"
  local actual="${name}.actual"
  echo -n "* ${name}: $2... "
  $1 > ${actual}
  [[ ! -f ${expect} ]] && { cp ${actual} ${expect}; }
  assert_eq ${expect} ${actual}
  echo "ok!"
  rm ${actual}
}

##########################################################
# Test 'hpx nested center list' (required to test 'porj' #
##########################################################

test_nested_center_list(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list
}
gtest "test_nested_center_list" "Get centers of the 12 base cells"

test_nested_center_list_csv(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d ,
}
gtest "test_nested_center_list_csv" "Get centers of the 12 base cells, with ',' output separator"


###########################################################
# Test 'hpx proj ...' (requires 'hpx nested center list') #
###########################################################

test_proj_value(){
  hpx proj value 0.0 0.0
  hpx proj value 0.0 +45.0
  hpx proj value 0.0 -45.0
}
gtest "test_proj_value" "Projection of a single value"

test_proj_list(){
    echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list \
    | hpx proj list
}
gtest "test_proj_list" "Projections of the center of the 12 base cells in a list"

test_proj_csv(){
    echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , \
    | hpx proj csv -d , 
}
gtest "test_proj_csv" "Projections of the center of the 12 base cells in csv"

test_proj_csv_paste(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , \
    | hpx proj csv -d , --paste
}
gtest "test_proj_csv_paste" "Projections of the center of the 12 base cells, with paste option"

test_proj_csv_head(){
    ( echo "center_lon,center_lat" && echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , ) \
    | hpx proj csv -d , --header
}
gtest "test_proj_csv_head" "Projections of the center of the 12 base cells, with header option"

test_proj_csv_head_paste(){
    ( echo "center_lon,center_lat" && echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , ) \
    | hpx proj csv -d , --header --paste --paste-header "proj_x,proj_y"
}
gtest "test_proj_csv_head_paste" "Projections of the center of the 12 base cells, with header and paste options"


#########################
# Test 'hpx unproj ...' #
#########################

test_unproj_value(){
  hpx unproj value 0.0  0.0
  hpx unproj value 0.0 +1.0
  hpx unproj value 0.0 -1.0
}
gtest "test_unproj_value" "Unproject single position"

test_unproj_list(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list \
    | hpx proj list \
    | hpx unproj list
}
gtest "test_unproj_list" "Unproject the centers of the 12 base cells"

test_unproj_csv_head_paste(){
    echo "d0h 0 1 2 3 4 5 6 7 8 9 10 11" \
    | tr ' ' '\n' \
    | hpx nested center 0 csv -d , --header --paste --paste-header "center_lon,center_lat" \
    | hpx proj   csv -d , -l 2 -b 3 --header --paste --paste-header "proj_x,proj_y" \
    | hpx unproj csv -d , -x 4 -y 5 --header --paste --paste-header "unproj_lon,unproj_lat"
}
gtest "test_unproj_csv_head_paste" "Join unprojections of the center of 12 base cells, with header"


##########################################################
# Test 'hpx nested ...' (hpx center list already tested) #
##########################################################

test_nested_table(){
  hpx nested table 
}
gtest "test_nested_table" "Get NESTED table"

test_nested_bestdepth(){
  echo "Best depth for 1 rad"
  hpx nested bestdepth 57.295779513
  hpx nested bestdepth --radian 1.0
  hpx nested bestdepth --arcmin 3437.7467707
  hpx nested bestdepth --arcsec 206264.80624
  echo "Best depth for 1 deg"
  hpx nested bestdepth 1.0
  hpx nested bestdepth --radian 0.01745329
  hpx nested bestdepth --arcmin 60.0
  hpx nested bestdepth --arcsec 3600.0
  echo "Best depth for 1 arcmin"
  hpx nested bestdepth 0.0166666666
  hpx nested bestdepth --radian 0.00029
  hpx nested bestdepth --arcmin 1.0
  hpx nested bestdepth --arcsec 60.0
  echo "Best depth for 1 arcsec"
  hpx nested bestdepth 0.0002777777
  hpx nested bestdepth --radian 0.000004847733
  hpx nested bestdepth --arcmin 0.01666666666
  hpx nested bestdepth --arcsec 1.0
  echo "Best depth for 1 mas"
  hpx nested bestdepth --arcsec 0.001
}
gtest "test_nested_bestdepth" "Get best starting depth for various radii"

test_nested_hash_value(){
  hpx nested hash         0 value 0.0 0.0
  hpx nested hash         0 value 0.0 +1.0
  hpx nested hash         0 value 0.0 -1.0
  hpx nested hash --uniq  0 value 0.0 0.0
  hpx nested hash --suniq 0 value 0.0 0.0
  hpx nested hash --zuniq 0 value 0.0 0.0
  hpx nested hash --with-dxdy         0 value 0.0 0.0
  hpx nested hash --with-dxdy --uniq  0 value 0.0 0.0
  hpx nested hash --with-dxdy --suniq 0 value 0.0 0.0
  hpx nested hash --with-dxdy --zuniq 0 value 0.0 0.0
}
gtest "test_nested_hash_value" "Compute hash value of (0.0, 0.0) with various options"

test_nested_hash_list(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list \
    | hpx nested hash 0 list
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list \
    | hpx nested hash --with-dxdy --uniq  0 list
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , \
    | hpx nested hash --with-dxdy --zuniq 0 list -d ,
}
gtest "test_nested_hash_list" "Compute hash values from a list"


test_nested_hash_csv(){
  ( echo "center_lon,center_lat" && echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
  | hpx nested center 0 list -d , ) \
  | hpx nested hash 14 csv -d , --header --paste --paste-header "hpx14"
}
gtest "test_nested_hash_csv" "Add hash values to a CSV file with header and paste"

test_nested_bilinear_value(){
  hpx nested bilinear 0 value 0.0 0.0
  hpx nested bilinear 1 value 0.0 0.0
  hpx nested bilinear 3 value 123.456789 +12.3456789
  hpx nested bilinear 3 value 123.456789 -12.3456789
}
gtest "test_nested_bilinear_value" "Compute bilinear values"

test_nested_bilinera_list(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list \
    | hpx nested bilinear 1 list
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , \
    | hpx nested bilinear 1 list -d ,
}
gtest "test_nested_bilinera_list" "Compute bilinera values from lists"


test_nested_bilinera_csv(){
  echo "d0h 0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 csv -d , --header --paste --paste-header "clon_deg,clat_deg"  \
    | hpx nested bilinear 1 csv -d , -l 2 -b 3 --header
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 csv -d , --paste  \
    | hpx nested bilinear 1 csv -d , -l 2 -b 3 --paste
  echo "d0h 0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 csv -d , --header --paste --paste-header "clon_deg,clat_deg"  \
    | hpx nested bilinear 1 csv -d , -l 2 -b 3 --header --paste --paste-header "h1,w1,h2,w2,h3,w3,h4,w4"
}
gtest "test_nested_bilinera_csv" "Compute bilinera values from CSVs"

test_nested_vertices_value(){
  hpx nested vertices 0 value 0
  hpx nested vertices 0 value 1
  hpx nested vertices 0 value 2
  hpx nested vertices 0 value 3
  hpx nested vertices 0 value 4
  hpx nested vertices 0 value 5
  hpx nested vertices 0 value 6
  hpx nested vertices 0 value 7
  hpx nested vertices 0 value 8
  hpx nested vertices 0 value 9
  hpx nested vertices 0 value 10
  hpx nested vertices 0 value 11
}
gtest "test_nested_vertices_value" "Compute base cells vertices"

test_nested_vertices_list(){
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested vertices 0 list -d '|'
}
gtest "test_nested_vertices_list" "Compute base cells vertices from a list"

test_nested_vertices_csv(){
  ( echo "center_lon,center_lat" && echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
    | hpx nested center 0 list -d , ) \
    | hpx nested hash 0 csv -d , --header --paste --paste-header "hpx0" \
    | hpx nested vertices 0 csv -d ',' -f 3 --header --paste --paste-header "slon,slat,wlon,wlat,elon,elat,nlon,nlat"
}
gtest "test_nested_vertices_csv" "Compute base cells vertices from a csv"


test_nested_neigh_value(){
  hpx nested neighbours 0 value 0
  hpx nested neighbours 0 value 1
  hpx nested neighbours 0 value 2
  hpx nested neighbours 0 value 3
  hpx nested neighbours 0 value 4
  hpx nested neighbours 0 value 5
  hpx nested neighbours 0 value 6
  hpx nested neighbours 0 value 7
  hpx nested neighbours 0 value 8
  hpx nested neighbours 0 value 9
  hpx nested neighbours 0 value 10
  hpx nested neighbours 0 value 11
}
gtest "test_nested_neigh_value" "Get base cells neighbours"

test_nested_neigh_list(){
  echo "S,SE,E,SW,NE,W,NW,N" && \
  echo "0 1 2 3 4 5 6 7 8 9 10 11" | tr ' ' '\n' \
  | hpx nested neighbours 0 list -d ,
}
gtest "test_nested_neigh_list" "Get base cells neighbours from a list"

test_nested_tozuniq_singledepth_csv(){
  hpx cov 9 cone --delta-depth 3 123.456789 -12.3456789 0.2 \
  | egrep -v "^#" \
  | hpx nested touniq --zuniq singledepth 9 csv -d , -f 2 --header --paste --paste-header "zuniq"
}
gtest "test_nested_tozuniq_singledepth_csv" "Compute zuniq from CSV, singledepth"

test_nested_tozuniq_multidepth_csv(){
  hpx cov 10 cone --delta-depth 3 123.456789 -12.3456789 0.2 \
  | egrep -v "^#" \
  | hpx nested touniq --zuniq multidepth csv -d , --header --paste --paste-header "zuniq" 
}
gtest "test_nested_tozuniq_multidepth_csv" "Comute zuniq from CSV, multidepth"

# test toring? (very similar to uniq...)


#######################
# Test 'hpx sort ...' #
#######################

test_sort(){
  # Sort a file
  hpx sort --header -d , -o hipmain.sorted.csv  hipmain.csv
  # Ensure is sorted (explicitely computing the HEALPi index):
  [[ "$(hpx nested hash 15 csv --header -d, --paste hipmain.sorted.csv | sort -c -n -t , -k 3)" != "" ]] && { echo "File not sorted!"; exit 1; }
  # Ensures we get the same result by computing the HEALPix index and using the linux sort...
  diff <(hpx nested hash 15 csv --header -d, --paste hipmain.sorted.csv) \
       <(hpx nested hash 15 csv --header -d, --paste hipmain.csv | sort -n -t , -k 3)
  [[ "$?" != "0" ]] && { echo "Diff on sorted files is not empty!"; exit 1; } 
}
test_sort

#########################
# Test 'hpx hcidx ...'  #
# Test 'hpx qhcidx ...' #
#########################

test_hcidx_qhicid(){
  hpx hcidx --header -d , --depth 5 -o hipmain.sorted.hci.fits hipmain.sorted.csv
  diff <(tail -n +2 hipmain.sorted.csv) \
       <(for i in $(seq 0 11); do hpx qhcidx hipmain.sorted.hci.fits cell 0 $i; done | egrep "^[0-9]")
  [[ "$?" != "0" ]] && { echo "Diff on qhcidx is not empt!"; exit 1; } 
}
test_hcidx_qhicid


######################
# Test 'hpx map ...' #
######################

test_map_count(){
  # Creates a first density maps from a generated sequence of points...
  for n in $(seq 0 47); do for i in $(seq 0 ${n}); do echo "$n"; done; done \
    | hpx nested center 1 list \
    | hpx map count 2 count.d2.l.fits list
  # ... and view it.
  hpx map view --silent --color-map-fn linear count.d2.l.fits count.d2.l.png allsky 200
 
  # Creates a second density maps from a generated sequence of points...
  for n in $(seq 0 191); do for i in $(seq 0 $((191-${n}))); do echo "$n"; done; done \
    | hpx nested center 2 list \
    | hpx map count 2 count.d2.r.fits list
  # ... and view it.
  hpx map view --silent --color-map-fn linear count.d2.r.fits count.d2.r.png allsky 200

  # Add the two maps...
  hpx map op add count.d2.l.fits count.d2.r.fits count.d2.m.fits
  # ... and view the result.
  hpx map view --silent --color-map-fn linear count.d2.m.fits count.d2.m.png allsky 200

  # Compare results
  assert_eq count.d2.l.png.expect count.d2.l.png
  assert_eq count.d2.r.png.expect count.d2.r.png
  assert_eq count.d2.m.png.expect count.d2.m.png
}	
test_map_count



######################
# Test 'hpx mom ...' #
######################

test_mom(){
  # Create a density map from a catalogue...
  hpx map dens 5 hip.d5.dens.fits csv -d , --header hipmain.csv
  # ... and view it.
  hpx map view --silent --color-map-fn linear hip.d5.dens.fits hip.d5.dens.png allsky 200
  # Build a MOM from the density map...
  hpx map convert --parallel 2 hip.d5.dens.fits hip.d5.dens.mom.fits dens2chi2mom
  # ... and view it.
  hpx mom view --silent --color-map-fn linear hip.d5.dens.mom.fits hip.d5.dens.mom.png allsky 200
  # Convert the MOM into in a BINTABLE MOM (check you can read it in, e.g.,  TOPCAT)
  hpx mom convert hip.d5.dens.mom.fits hip.d5.dens.mom.bintable.fits bintable

  # Compare results
  assert_eq hip.d5.dens.png.expect hip.d5.dens.png
  assert_eq hip.d5.dens.mom.png.expect hip.d5.dens.mom.png
}

test_mom



######################
# Test 'hpx cov ...' #
######################

test_cov_cone(){
  local r="0.2"
  hpx cov $(($(hpx nested bestdepth $r)+3)) cone --delta-depth 3 123.456789 +12.3456789 $r
  hpx cov $(($(hpx nested bestdepth $r)+3)) cone --delta-depth 3 123.456789 -12.3456789 $r
}
gtest "test_cov_cone" "Cone BMOC"

test_cov_ellipse(){
  local a="0.2"
  local b="0.1"
  hpx cov $(($(hpx nested bestdepth $a)+3)) ellipse --delta-depth 3 123.456789 +12.3456789 $a $b 23.0
  hpx cov $(($(hpx nested bestdepth $a)+3)) ellipse --delta-depth 3 123.456789 -12.3456789 $a $b 23.0
}
gtest "test_cov_ellipse" "Ellipse BMOC"

test_cov_ring(){
  local R="0.2"
  local r="0.1"
  hpx cov $(($(hpx nested bestdepth $r)+3)) ring --delta-depth 3 123.456789 +12.3456789 $R $r
  hpx cov $(($(hpx nested bestdepth $r)+3)) ring --delta-depth 3 123.456789 -12.3456789 $R $r
}
gtest "test_cov_ring" "Ring BMOC"

test_cov_zone(){
  hpx cov 4 zone 0.0 -10.0 20.0 10.0
}
gtest "test_cov_zone" "Zone BMOC"

test_cov_box(){
  local a="0.2"
  local b="0.1"
  hpx cov $(($(hpx nested bestdepth $a)+3)) box 123.456789 +12.3456789 $a $b 23.0
  hpx cov $(($(hpx nested bestdepth $a)+3)) box 123.456789 -12.3456789 $a $b 23.0
}
gtest "test_cov_box" "Box BMOC"

test_cov_polygon(){
  hpx cov 5 polygon 10.0,-10.0,10.0,+20.0,20.0,+20.0,20.0,-10.0
}
gtest "test_cov_polygon" "Polygon BMOC"

test_cov_stcs(){
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
}
gtest "test_cov_stcs" "STC-S BMOCs"

test_cov_fits(){
  hpx cov 9 --out-type fits --output-file cone.bmoc.fits cone 080.8942 -69.7561 0.2
  # Also compute center to check we do are aroudn the center of the cone
  hpx cov 9 --out-type csv convert cone.bmoc.fits \
    | tail -n +2 | hpx nested center 9 csv -d , -f 2 --header --paste --paste-header "lon,lat"
}
gtest "test_cov_fits" "Make a FITS BMOC and convert it into CSV"

