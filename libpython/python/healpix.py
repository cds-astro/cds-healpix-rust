import math
import numpy as np

import datetime
import time
from pprint import pprint

from cffi import FFI

def _as_u64_array(array):
    """ Cast np.long64 array to a pointer to long 64s."""
    return ffi.cast("unsigned long int*", array.ctypes.data)

ffi = FFI()
ffi.cdef("""

   // Returns the cell number (hash value) associated with the given position on the unit sphere
   // in the HEALPix NESTED scheme.
   // Inputs:
   // - depth: HEALPix depth, must be in [0, 24?]
   // - lon: longitude, in radians
   // - lat: latitude, in radians, must be in [-pi/2, pi/2]
   // Output:
   // - the nested cell number if the given position a thte given depth
   unsigned long int hpx_hash(unsigned char depth, double lon, double lat);

   // Returns the cell number (hash value) associated with the given positions on the unit sphere 
   // in the HEALPix NESTED scheme.
   // Inputs
   // - depth: HEALPix depth, must be in [0, 24?]
   // - n_elems: number of positions submited (= size of the `coords` array divided by 2 = size of `result`)
   // - coords: coordinates of the position, in radians, following [lon1, lat1, lon2, lat2, ..., lonN, latN]
   // - result: array storing the result [hash1, hash2, ..., hashN]
   // Output:
   // - no output: the result is stored in the `result` array
   // We use the `result` array so that the memory is managed by Python and do not have to be free
   // by an explicit call to a specific free function.
   void* hpx_hash_multi(unsigned short depth, int n_elems, double* coords, unsigned long int* result);   
  
     
   void* hpx_center(unsigned char depth, unsigned long int hash, double* lon, double* lat);
     
   void* hpx_center_multi(unsigned char depth, int n_elems, unsigned long int* hash_ptr, double* res_ptr);
   
   void* hpx_vertices(unsigned char depth, unsigned long int hash, double* res_ptr);
     
   void* hpx_neighbours(unsigned char depth, unsigned long int hash, long int* res_ptr);
     
   // Structure storing a BMOC cell informations, i.e.
   // the cell depth, the cell number and a flag telling if the cell is fylly (1) or partially (0) covered
   typedef struct {
       unsigned char depth;
       unsigned long int hash;
       unsigned char flag;
   } bmoc_cell;
  
   // BMOC: simple ordrered array of BMOC cells
   typedef struct {
     int ncells;
     bmoc_cell* cells;
   } bmoc;

   // Free the BMOC memory owned by Rust
   void bmoc_free(bmoc* bmoc);
   

   bmoc* hpx_query_cone_approx(unsigned char depth, double lon, double lat, double radius);

   bmoc* hpx_query_cone_approx_custom(unsigned char depth, unsigned chardelta_depth, double lon, double lat, double radius);

   bmoc* hpx_query_polygon_approx(unsigned char depth, int n_vertices, double* vertices_coords);


""")

# C = ffi.dlopen("../target/release/libcdshealpix_ffi.so")
C = ffi.dlopen("/home/pineau/IdeaProjects/rust-healpix/libpython/target/release/libcdshealpix_ffi.so")

# SIMPLE HASH
print("Hash: ")
a = datetime.datetime.now()
h = C.hpx_hash(6, 0.5, 0.4)
print str(datetime.datetime.now() - a)
print(h)

a = datetime.datetime.now()
h = C.hpx_hash(6, 0.72, -0.1)
print str(datetime.datetime.now() - a)
print(h)

a = datetime.datetime.now()
h = C.hpx_hash(8, 1.2, -0.12)
print str(datetime.datetime.now() - a)
print(h)

# MULTI-HASH: NOT NUMPY VERSION
res = ffi.new("unsigned long int[2]")
print("Hashes: ")
a = datetime.datetime.now()
C.hpx_hash_multi(6, 2, [0.5, 0.4, 0.5, -0.4], res)
print str(datetime.datetime.now() - a)
print(res)
print(res[0])
print(res[1])

# MUTLI-HASH: NUMPY VERSION
#res2 = np.array([0,0])
res2 = np.zeros((3,), dtype=np.uint64)
print("Hashes: ")
print(res2)
res2_ptr = _as_u64_array(res2)
C.hpx_hash_multi(10, 3, [0.0, 0.0, 0.0, 0.5, 0.25, 0.25], res2_ptr)
for i in range(3):
  print(res2[i])

print("##############")
print("#### CONE ####")

a = datetime.datetime.now()
cone = C.hpx_query_cone_approx(3, math.radians(36.80105218), math.radians(56.78028536), math.radians(14.93))
print str(datetime.datetime.now() - a)
print("N cells: ")
print(cone.ncells)

a = datetime.datetime.now()
cone = C.hpx_query_cone_approx(8, math.radians(36.80105218), math.radians(56.78028536), math.radians(14.93))
print str(datetime.datetime.now() - a)
print("N cells: ")
print(cone.ncells)
#for i in range(cone.ncells):
#  print "cell1: ", cone.cells[i].depth, cone.cells[i].hash, cone.cells[i].flag

cone = C.hpx_query_cone_approx_custom(3, 2, math.radians(36.80105218), math.radians(56.78028536), math.radians(14.93))
print("N cells: ")
print(cone.ncells)
for i in range(cone.ncells):
  print "cell1: ", cone.cells[i].depth, cone.cells[i].hash, cone.cells[i].flag



print("#################")
print("#### POLYGON ####")

poly = C.hpx_query_polygon_approx(3, 3, [0.0, 0.0, 0.0, 0.5, 0.25, 0.25])
print(poly)
print("N cells: ")
print(poly.ncells)

for i in range(poly.ncells):
  print("cell1: ")
  print(poly.cells[i].depth)
  print(poly.cells[i].hash)
  print(poly.cells[i].flag)

#print("cell2: ")
#print(poly.cells[1].depth)
#print(poly.cells[1].hash)
#print(poly.cells[1].flag)

lons = [83.71315909, 83.71378887, 83.71297292, 83.71233919] 
lats = [-5.44217436,-5.44298864, -5.44361751, -5.4428033]
coords = []
    
for i in range(len(lons)):
    coords.append(lons[i] * np.pi/180.0)
    coords.append(lats[i] * np.pi/180.0)

a = datetime.datetime.now()
poly = C.hpx_query_polygon_approx(20, len(lats), coords)
# time.sleep(0.01)
print str(datetime.datetime.now() - a)

ipix_d = {'17': [], '18': [], '19': [], '20': []}
for i in range(poly.ncells):
    depth = poly.cells[i].depth
    ipix = poly.cells[i].hash

    ipix_d[str(depth)].append(ipix)

pprint(ipix_d)

a = datetime.datetime.now()
C.bmoc_free(poly)
print str(datetime.datetime.now() - a)



