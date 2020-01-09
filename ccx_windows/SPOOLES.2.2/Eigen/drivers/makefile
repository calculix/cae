include ../../Make.inc

#-----------------------------------------------------------------------
#
#  set suffix rule *.c --> *.o
#
.c.o : 
	$(PURIFY) $(CC) -c $(CFLAGS) $(MPI_INCLUDE_DIR) $<
#
#-----------------------------------------------------------------------
#
#  redefine CC for mpi on sgi
#
# CCmpi = /opt/mpi/bin/mpicc
#
#-----------------------------------------------------------------------

F77 = f77

DRIVERS = testSerial testMT testMPI

#-----------------------------------------------------------------------
# environment on sparcstation 20
#
# locations of eigensolver libraries
#
LANCZOS_ST_LIB  = /home/pattersn/LANCZOS/lib/Serial_lncz.a
LANCZOS_MT_LIB  = /home/pattersn/LANCZOS/lib/MT_lncz.a
LANCZOS_MPI_LIB = /home/pattersn/LANCZOS/lib/MPI_lncz.a
#
# locations of spooles libraries
#
SPOOLES_DIR     = ../..
SPOOLES_LIB     = $(SPOOLES_DIR)/spooles.a
SPOOLES_MT_LIB  = $(SPOOLES_DIR)/MT/src/spoolesMT.a
SPOOLES_MPI_LIB = $(SPOOLES_DIR)/MPI/src/spoolesMPI.a
#
# libraries needed for the eigensolver
#
LIBS_FOR_ST     = $(LANCZOS_ST_LIB) ../srcST/Bridge.a $(SPOOLES_LIB) -lm
LIBS_FOR_MT     = $(LANCZOS_MT_LIB) ../srcMT/BridgeMT.a \
                  $(SPOOLES_MT_LIB) $(SPOOLES_LIB) ${THREAD_LIBS} -lm
LIBS_FOR_MPI    = $(LANCZOS_MPI_LIB) ../srcMPI/BridgeMPI.a \
                  $(SPOOLES_MPI_LIB) $(SPOOLES_LIB) ${MPI_LIBS} -lm
#-----------------------------------------------------------------------
#
# dependencies for Bridge objects
# 
STlibs :
	cd ../srcST ; make Bridge.a
MTlibs :
	cd ../srcMT ; make BridgeMT.a
MPIlibs :
	cd ../srcMPI ; make BridgeMPI.a
#-----------------------------------------------------------------------

all_drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

testSerial : testSerial.o STlibs ../../spooles.a
	${PURIFY} ${F77} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS_FOR_ST}

testMT : testMT.o MTlibs ../../spooles.a
	${PURIFY} ${F77} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS_FOR_MT}

testMPI : testMPI.o MPIlibs ../../spooles.a
	${PURIFY} ${F77} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS_FOR_MPI}
