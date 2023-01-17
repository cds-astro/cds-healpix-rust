
# PSQL extension for HEALPix

This PSQL extension installs HEALPix functions in a PostgreSQL database.

## Available functions


* **`hpx_nside(INTEGER)->INTEGER`** :
  compute the Nside for the given HEALPix depth.
  
  _Example:_
  
  ```
  psql=# SELECT hpx_nside(4);
  
  hpx_nside
  -----------
          16
  (1 row)
  ```
  
  _Warning: If the given depth is not inside [0,29], the query stops immediately
   with an error message._
  
* **`hpx_hash(INTEGER, DOUBLE PRECISION, DOUBLE PRECISION)->BIGINT`** :
  compute the hash/index of the HEALPix cell including the given position at the
  given depth.
  
  _Example:_
  
  ```
  psql=# SELECT hpx_hash(6, 0.0, 0.0);
   hpx_hash
  ----------
      19456
  (1 row)
  ```
  
  _Warning: If the given depth is not inside [0,29] or if the given latitude is
   not inside [-90,90], the query stops immediately with an error message._

* **`hpx_center(INTEGER, BIGINT)->DOUBLE PRECISION[]`** :
  compute the center of the specified HEALPix cell.
  
  _Example:_
  
  ```
  psql=# SELECT hpx_center(6, 1234);
            hpx_center          
  ------------------------------
   {69.609375,34.2288663278126}
  (1 row)
  ```
  
 _Warning: If the given depth is not inside [0,29] or if the given hash is not
  inside [0,12*Nside^2[, the query stops immediately with an error message._

 
## Installation

```bash
RUSTFLAGS='-C target-cpu=native' cargo build --release
make
make install
```

Then, to use this HEALPix extension in the database, you must first run the
following SQL instruction (only once):

```sql
CREATE EXTENSION pg_cds_healpix;
```


## Regression tests

If the installation is successful, the extension can be tested by running
several regression tests with the following command:

```bash
make installcheck
```

If you get an error like below:

```
============== dropping database "contrib_regression" ==============

psql: FATAL: role "xxxx" does not exist ... ...
```

...where `xxxx` is the name of the current user, try providing the name of a
user having permissions to create a database on your PostgreSQL installation
(e.g. the `postgres` user):

```bash
PGUSER=postgres make installcheck
```


## Help about PostgreSQL extensions

* **Tutorials:**
  - <http://big-elephants.com/2015-10/writing-postgres-extensions-part-i/>
  - <https://www.postgresql.org/docs/9.4/static/extend-pgxs.html>

* **Some PostgreSQL documentation:**
  - [C-Language Functions in PostgreSQL extension API](https://www.postgresql.org/docs/current/xfunc-c.html)
  - [PSQL datatypes](https://docs.postgresql.fr/10/xfunc-c.html)
  - [fmgr.h](https://doxygen.postgresql.org/fmgr_8h.html)
  - [Examples of array manipulation](https://github.com/pjungwir/aggs_for_arrays)
  - [CREATE AGGREGATE](https://www.postgresql.org/docs/current/sql-createaggregate.html)

* **Hints:**
  - To find out where is $libdir, execute the following command:
    ```bash
    pg_config --pkglibdir
    ```
