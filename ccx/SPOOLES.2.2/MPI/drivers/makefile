include ../../Make.inc
 
#-----------------------------------------------------------------------
#
#  set suffix rule *.c --> *.o
#
.c.o : 
	$(PURIFY) $(CC) -c $(CFLAGS) $(MPI_INCLUDE_DIR) $<
#

#-----------------------------------------------------------------------

LIBS = ../src/spoolesMPI.a ../../spooles.a $(MPI_LIBS) -lm

DRIVERS = allInOneMPI patchAndGoMPI testGridMPI testIV_allgather \
          testIVL_allgather testIVL_alltoall testSplitDenseMtx \
          testSplitInpMtx testSymbFac testGather testMMM \
          testIVL_Bcast testGraph_Bcast \
          testScatterDenseMtx testScatterInpMtx
 
libs :
	cd ../src ; make spoolesMPI.a 

drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
allInOneMPI : allInOneMPI.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
patchAndGoMPI : patchAndGoMPI.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testGather : testGather.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testGraph_Bcast : testGraph_Bcast.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testGridMPI : testGridMPI.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testMMM : testMMM.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testScatterDenseMtx : testScatterDenseMtx.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${LIBS} 
 
testScatterInpMtx : testScatterInpMtx.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${LIBS} 
 
testSymbFac : testSymbFac.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${LIBS} 
 
testSplitDenseMtx : testSplitDenseMtx.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testSplitInpMtx : testSplitInpMtx.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
 
testIV_allgather : testIV_allgather.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
testIVL_allgather : testIVL_allgather.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
testIVL_alltoall : testIVL_alltoall.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
testIVL_Bcast : testIVL_Bcast.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

