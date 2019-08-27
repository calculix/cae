include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = checkComponents \
          compressGraph \
          expandGraph \
          mkGridGraph \
          testIO \
          testIsSymmetric \
          testWirebasket\
          writeMetisFile
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

checkComponents : checkComponents.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

compressGraph : compressGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

expandGraph : expandGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

mkGridGraph : mkGridGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testChacoIO : testChacoIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testIO : testIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testIsSymmetric : testIsSymmetric.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testWirebasket : testWirebasket.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

writeMetisFile : writeMetisFile.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

writeAIJ : writeAIJ.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
