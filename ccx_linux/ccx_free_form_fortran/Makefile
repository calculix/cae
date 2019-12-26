
DIR=../SPOOLES.2.2

CFLAGS = -Wall -O2 -I $(DIR) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -O2 -ffree-form

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.16.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

LIBS = \
       $(DIR)/spooles.a \
	../ARPACK/libarpack_INTEL.a \
       -lpthread -lm -lc

ccx_2.16: $(OCCXMAIN) ccx_2.16.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.16.c; $(FC) -Wall -O2 -o $@ $(OCCXMAIN) ccx_2.16.a $(LIBS)

ccx_2.16.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
