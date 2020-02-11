CREATE EXTENSION pg_cds_healpix;

-- TEST HPX_NSIDE

SELECT hpx_nside(0);  -- expected: 1
SELECT hpx_nside(5);  -- expected: 32
SELECT hpx_nside(13); -- expected: 8192
SELECT hpx_nside(NULL); -- expected: NULL
SELECT hpx_nside(30); -- expected: ERROR!
SELECT hpx_nside(-1); -- expected: ERROR!

-- TEST HPX_HASH

SELECT hpx_hash(0, 0, 0);                   -- expected: 4
SELECT hpx_hash(15, 250, -90);              -- expected: 10737418240
SELECT hpx_hash(NULL, 0, 0);                -- expected: NULL
SELECT hpx_hash(0, NULL, 0);                -- expected: NULL
SELECT hpx_hash(0, 0, NULL);                -- expected: NULL
SELECT hpx_hash(30, 0, 0);                  -- expected: ERROR!
SELECT hpx_hash(-1, 0, 0);                  -- expected: ERROR!
SELECT hpx_hash(0, 0, 91);                  -- expected: ERROR!
SELECT hpx_hash(0, 0, -91);                 -- expected: ERROR!
SELECT hpx_hash(0, 'NaN'::float8, 0);       -- expected: ERROR!
SELECT hpx_hash(0, 'Infinity'::float8, 0);  -- expected: ERROR!
SELECT hpx_hash(0, '-Infinity'::float8, 0); -- expected: ERROR!
SELECT hpx_hash(0, 0, 'NaN'::float8);       -- expected: ERROR!
SELECT hpx_hash(0, 0, 'Infinity'::float8);  -- expected: ERROR!
SELECT hpx_hash(0, 0, '-Infinity'::float8); -- expected: ERROR!

-- TEST HXP_CENTER

SELECT hpx_center(0, 4);            -- expected: {0,0}
SELECT hpx_center(15, 10737418240); -- expected: {225,-89.9985723325173}
SELECT hpx_center(NULL, 4);         -- expected: NULL
SELECT hpx_center(0, NULL);         -- expected: NULL
SELECT hpx_center(30, 4);           -- expected: ERROR!
SELECT hpx_center(-1, 4);           -- expected: ERROR!
SELECT hpx_center(0, 13);           -- expected: ERROR!

