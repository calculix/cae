include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm
 
DRIVERS = testIO testDM extractBPG
 
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}

testIO : testIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
testDM : testDM.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
extractBPG : extractBPG.o  ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

