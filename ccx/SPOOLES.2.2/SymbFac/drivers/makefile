include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = testSymbfacInpMtx \
          testSymbfacGraph
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

testSymbfacGraph : testSymbfacGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testSymbfacInpMtx : testSymbfacInpMtx.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
