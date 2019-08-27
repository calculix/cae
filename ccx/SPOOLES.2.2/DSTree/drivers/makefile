include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = testIO writeStagesIV testDomWeightStages
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

testIO : testIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testDomWeightStages : testDomWeightStages.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

writeStagesIV : writeStagesIV.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
