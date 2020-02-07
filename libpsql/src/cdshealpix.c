#include "cdshealpix.h"



////////////////////////////////////////////////////////////////////////////////
//
// ARGUMENT TESTING methods
//

/**
 * Check the given HEALPix depth.
 *
 * If not inside [0,29], report an error in the PostgreSQL backend.
 * In such case, the query execution is immediately stopped.
 *
 * @param depth  The HEALPix depth to check.
 */
void check_hpx_depth(uint8_t depth) {
    if (depth > 29) {
        ereport(ERROR,
	      (
	        errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
	        errmsg("Incorrect HEALPix depth: %d!", depth),
	        //, errdetail("value %d not inside [0, 29]", depth),
            errhint("An HEALPix depth must be inside [0,29].")
	      )
	    );
    }
}

/**
 * Check the given latitude (in degrees).
 *
 * If not inside [-90,90], report an error in the PostgreSQL backend.
 * In such case, the query execution is immediately stopped.
 *
 * @param lat_deg  The latitude to check.
 */
void check_lat(float8 lat_deg) {
    if (lat_deg > 90 || lat_deg < -90) {
        ereport(ERROR,
	      (
	        errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
	        errmsg("Incorrect latitude: %2f!", lat_deg),
	        //errdetail(" %f not inside [-90, 90]", lat_deg),
	        errhint("A latitude must be inside [-90,90] degrees.")
	      )
	    );
    }
}

/**
 * Check the given HEALPix cell hash/index.
 *
 * If the given hash is not inside [0,12*nside^2], report an error in the
 * PostgreSQL backend. In such case, the query execution is immediately stopped.
 *
 * @param depth  HEALPix depth.
 * @param hash   The HEALPix cell hash/index to check.
 */
void check_hpx_hash(uint8_t depth, uint64_t hash) {
    uint64_t nsideValue = (uint64_t) nside(depth);
    uint64_t max_hash = 12*nsideValue*nsideValue - 1;
    if (hash < 0 || hash > max_hash) {
        ereport(ERROR,
	      (
	        errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
	        errmsg("Incorrect HEALPix cell hash: %lu!", hash),
	        //, errdetail("value %d not inside [0, 29]", depth),
            errhint("An HEALPix cell hash at depth %d must be inside [0,%lu].", depth, max_hash)
	      )
	    );
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// NSIDE method
//

PG_FUNCTION_INFO_V1(hpx_nside);

/**
 * [ hpx_nside(long) --> long ]
 *
 * Compute the HEALPix nside for the given HEALPix depth.
 *
 * An error is reported to the PostgreSQL backend if the given depth is
 * incorrect (see check_hpx_index(...)).
 *
 * @param depth  HEALPix depth.
 *
 * @return The corresponding nside.
 *
 * @depends check_hpx_depth(long)
 */
Datum hpx_nside(PG_FUNCTION_ARGS) {
    // Get argument (healpix depth):
    uint8_t depth = (uint8_t) PG_GETARG_INT32(0);
    // Check it:
    check_hpx_depth(depth);
    // Compute the nside:
    PG_RETURN_UINT32( nside(depth) );
}


////////////////////////////////////////////////////////////////////////////////
//
// HASH method
//

PG_FUNCTION_INFO_V1(hpx_hash);

/**
 * [ hpx_hash(long, float, float) --> long ]
 *
 * Compute the HEALPix cell hash (or index) for the given longitude and
 * latitude at the given HEALPix depth.
 *
 * An error is reported to the PostgreSQL backend if the given depth or
 * latitude is incorrect (see check_hpx_depth(...) and check_latitude(...)).
 *
 * @param depth    HEALPix depth.
 * @param lon_deg  Longitude (in degrees).
 * @param lat_deg  Lattitude (in degrees)
 *
 * @return  The corresponding HEALPix hash/index.
 *
 * @depends check_hpx_depth(long)
 * @depends check_hpx_latitude(float)
 */
Datum hpx_hash(PG_FUNCTION_ARGS) {
    // Get arguments (depth, lon, lat):
    uint8_t depth = (uint8_t) PG_GETARG_INT32(0);
    float8 lon_deg = PG_GETARG_FLOAT8(1);
    float8 lat_deg = PG_GETARG_FLOAT8(2);
    // Check them:
    check_hpx_depth(depth);
    check_lat(lat_deg);
    // Compute the HEALPix cell hash:
    PG_RETURN_INT64( (int64_t) nest_hash(depth, lon_deg, lat_deg) );
}


////////////////////////////////////////////////////////////////////////////////
//
// CENTER method
//

PG_FUNCTION_INFO_V1(hpx_center);

/**
 * [ hpx_center(long, long) --> float* ]
 *
 * Compute the center of the specified HEALPix cell.
 *
 * An error is reported to the PostgreSQL backend if the given depth is
 * incorrect (see check_hpx_index(...)).
 *
 * @param depth  HEALPix depth.
 * @param hash   Hash/Index of an HEALPix cell at the given depth.
 *
 * @return  An array of 2 float values: [ lon_deg, lat_deg ].
 *
 * @depends check_hpx_depth(long)
 */
Datum hpx_center(PG_FUNCTION_ARGS) {
    float8* coordsRust;
    Datum*  coordsPSQL;

    // Get arguments (depth, hpx hash):
    uint8_t depth = (uint8_t) PG_GETARG_INT32(0);
    uint64_t hash = PG_GETARG_INT64(1);

    // Check them:
    check_hpx_depth(depth);
    check_hpx_hash(depth, hash);

    // Compute the center of the specified HEALPix cell:
    coordsRust = palloc0(sizeof(float8) * 2);
    nest_center(depth, hash, coordsRust);

    // Cast the result to Datum array (for Postgres):
    coordsPSQL = palloc0(sizeof(Datum) * 2);
    coordsPSQL[0] = Float8GetDatum(coordsRust[0]);
    coordsPSQL[1] = Float8GetDatum(coordsRust[1]);
    pfree(coordsRust);
    
    // Return the final PostgreSQL array object:
    //   note: construct_array returns an ArrayType* ; its parameters are:
    //           - pointer on the array (supposed to be a Datum*)
    //           - number of array items
    //           - OID of the items' datatype
    //           - size of this datatype
    //           - whether items should be passed by value (true) or reference (false)
    //           - items alignment ('d' for double, 'i' for integer, ...)
    PG_RETURN_ARRAYTYPE_P( construct_array(coordsPSQL, 2, FLOAT8OID, sizeof(float8), true, 'd') );
}



