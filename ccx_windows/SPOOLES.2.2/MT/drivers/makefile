include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../src/spoolesMT.a ../../spooles.a ${THREAD_LIBS} -lm

DRIVERS = testMMM \
          testGridMT \
          testQRgridMT \
          allInOneMT \
          patchAndGoMT \
          testFactorMT

libs :
	cd ../src ; make spoolesMT.a

drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

allInOneMT : allInOneMT.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testFactorMT : testFactorMT.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 

testMMM : testMMM.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 

testQRgridMT : testQRgridMT.o libs ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS} 
 
testGridMT : testGridMT.o libs ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS} 
 
patchAndGoMT : patchAndGoMT.o libs ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

