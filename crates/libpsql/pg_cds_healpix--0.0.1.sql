
CREATE FUNCTION hpx_nside(integer) RETURNS integer
  AS '$libdir/cdshealpix', 'hpx_nside'
  LANGUAGE C
  STRICT
  IMMUTABLE;

CREATE FUNCTION hpx_hash(integer, double precision, double precision) RETURNS bigint
  AS '$libdir/cdshealpix', 'hpx_hash'
  LANGUAGE C
  STRICT
  IMMUTABLE;

CREATE FUNCTION hpx_center(integer, bigint) RETURNS float8[]
  AS '$libdir/cdshealpix', 'hpx_center'
  LANGUAGE C
  STRICT
  IMMUTABLE;

