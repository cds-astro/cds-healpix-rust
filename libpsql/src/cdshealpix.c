#include <string.h>
#include <stdint.h>

#include "postgres.h"
#include "fmgr.h"
#include "utils/builtins.h"
#include <limits.h>

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

// See:
// - https://www.postgresql.org/docs/current/xfunc-c.html
// - https://www.postgresql.org/docs/current/sql-createaggregate.html

// see psql dataytpes:
// - https://docs.postgresql.fr/10/xfunc-c.html

// see https://doxygen.postgresql.org/fmgr_8h.html

//////////////////
// NSIDE method //

extern uint32_t nside(uint32_t);

PG_FUNCTION_INFO_V1(hpx_nside);

Datum hpx_nside(PG_FUNCTION_ARGS) {
    uint8_t depth = (uint8_t) PG_GETARG_INT32(0);
    if (depth > 29) {
        ereport(ERROR,
	      (
	        errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
	        errmsg("wrong depth value"),
	        errdetail("value %d not in [0, 29]", depth),
	        errhint("change the depth value")
	      )
	    );
    }
    PG_RETURN_UINT32(nside(depth));
}

/////////////////
// HASH method //

extern uint64_t hash(uint8_t, float8, float8);

PG_FUNCTION_INFO_V1(hpx_hash);

Datum hpx_hash(PG_FUNCTION_ARGS) {
    uint8_t depth = (uint8_t) PG_GETARG_INT32(0);
    float8 lon_deg = PG_GETARG_FLOAT8(1);
    float8 lat_deg = PG_GETARG_FLOAT8(2);
    if (depth > 29) {
        ereport(ERROR,
	     (
	       errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
	       errmsg("wrong depth value"),
	       errdetail("value %d not in [0, 29]", depth),
	       errhint("change the depth value")
	     )
	  );
    }
    PG_RETURN_INT64((int64_t) hash(depth, lon_deg, lat_deg));
}


