include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm

DRIVERS = testGrid testQRgrid testFactor
 
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
 
testQRgrid : testQRgrid.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o \
                        -o $@ ${LIBS} ${THREAD_LIBS}
 
testGrid : testGrid.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o \
                        -o $@ ${LIBS} ${THREAD_LIBS}
 
testFactor : testFactor.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o \
                        -o $@ ${LIBS} ${THREAD_LIBS}
