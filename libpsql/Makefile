
EXTENSION = pg_cds_healpix
DATA = pg_cds_healpix--0.0.1.sql
REGRESS = simpletests

MODULE_big = cdshealpix
OBJS = src/cdshealpix.o target/release/libcdshealpix_libpsql.a

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

#clean:
#	@rm src/cdshealpix.o src/cdshealpix.so

