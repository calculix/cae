include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm
 
DRIVERS = testDrand
 
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
testDrand : testDrand.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

