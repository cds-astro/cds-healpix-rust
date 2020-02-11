#ifndef PG_CDS_HEALPIX
#define PG_CDS_HEALPIX

#include <string.h>
#include <stdint.h>

#include "postgres.h"
#include "fmgr.h"
#include "math.h"
#include <utils/array.h>
#include <catalog/pg_type.h>
#include <utils/lsyscache.h>
#include "utils/builtins.h"
#include <limits.h>

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif


////////////////////////////////////////////////////////////////////////////////
//
// RUST FUNCTIONS DECLARATION (see details in lib.rs)
//

/** [ nside(ushort depth) --> uint ] */
extern uint32_t nside(uint8_t);

/** [ nest_hash(ushort, float, float) --> ulong ] */
extern uint64_t nest_hash(uint8_t, float8, float8);

/** [ nest_center(ushort, ulong, float*) ] */
extern void nest_center(uint8_t, uint64_t, float8*);


////////////////////////////////////////////////////////////////////////////////
//
// ARGUMENT TESTING methods
//

/** [ check_hpx_depth(ushort depth) ] */
void check_hpx_depth(uint8_t);

/** [ check_lon(float lon_deg) ] */
void check_lon(float8);

/** [ check_lat(float lat_deg) ] */
void check_lat(float8);

/** [ check_hpx_hash(ushort depth, ulong hash) ] */
void check_hpx_hash(uint8_t, uint64_t);


////////////////////////////////////////////////////////////////////////////////
//
// HEALPIX FUNCTION FOR SQL
//

/** [ hpx_nside(long) --> long ] */
Datum hpx_nside(PG_FUNCTION_ARGS);

/** [ hpx_hash(long, float, float) --> long ] */
Datum hpx_hash(PG_FUNCTION_ARGS);

/** [ hpx_center(long depth, long hash) --> float* ] */
Datum hpx_center(PG_FUNCTION_ARGS);

#endif
