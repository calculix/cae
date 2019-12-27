include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = testDDviaFishnet \
          testTwoSetViaBKL \
          testSmoothBisector \
          testRBviaDDsep \
          mkDSTree \
          testRBviaDDsep2
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

mkDSTree : mkDSTree.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testSmoothBisector : testSmoothBisector.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testTwoSetViaBKL : testTwoSetViaBKL.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testDDviaFishnet : testDDviaFishnet.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testRBviaDDsep2 : testRBviaDDsep2.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testRBviaDDsep : testRBviaDDsep.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
