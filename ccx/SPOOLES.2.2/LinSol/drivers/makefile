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
#  redefine CC for mpi
#
# CCmpi = /opt/mpi/bin/mpicc
#
#-----------------------------------------------------------------------


DRIVERS = testWrapper testWrapperMT testWrapperMPI
 
STlib :
	cd ../srcST ; make Bridge.a

MTlib :
	cd ../srcMT ; make BridgeMT.a
	cd ../../MT/src ; make spoolesMT.a

MPIlib :
	cd ../srcMPI ; make BridgeMPI.a
	cd ../../MPI/src ; make spoolesMPI.a

LIBS_FOR_ST = ../srcST/Bridge.a ../../spooles.a -lm

LIBS_FOR_MT = ../srcMT/BridgeMT.a ../../MT/src/spoolesMT.a \
              ../../spooles.a ${THREAD_LIBS} -lm

LIBS_FOR_MPI = ../srcMPI/BridgeMPI.a ../../MPI/src/spoolesMPI.a \
               ../../spooles.a ${MPI_LIBS} -lm

drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

testWrapper : STlib testWrapper.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS_FOR_ST}

testWrapperMT : MTlib testWrapperMT.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS_FOR_MT}

testWrapperMPI : MPIlib testWrapperMPI.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS_FOR_MPI}
