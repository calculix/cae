include ../../Make.inc
 
#-----------------------------------------------------------------------
 
DRIVERS = test_norms testQR

LIBS = ../../spooles.a -lm
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}
 
testQR : testQR.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_norms : test_norms.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIF_GCC_VERSION} ${LIBS}
