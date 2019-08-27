include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a ${THREAD_LIBS} -lm

DRIVERS = testSimple testGrid
 
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
testSimple : testSimple.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}
 
testGrid : testGrid.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}
