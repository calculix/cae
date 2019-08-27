include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm

DRIVERS = testIO \
          test_scalevec \
          test_solve \
          test_solveH \
          test_solveT \
          test_solveupd \
          test_solveupdH \
          test_solveupdT \
          test_sort

drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
testIO : testIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_scalevec : test_scalevec.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_solve : test_solve.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_solveH : test_solveH.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_solveT : test_solveT.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_solveupd : test_solveupd.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_solveupdH : test_solveupdH.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_solveupdT : test_solveupdT.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_sort : test_sort.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
