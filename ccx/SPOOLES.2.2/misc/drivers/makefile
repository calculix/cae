include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm
LIBS = ../../spooles.a -lpthread -lm

DRIVERS = allInOne drawGraph patchAndGo QRallInOne testOrderViaND \
          testOrderViaMS testOrderViaMMD testNDperm testSemi \
          testOptPart testOrderViaBestOfNDandMS
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

allInOne : allInOne.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

patchAndGo : patchAndGo.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

QRallInOne : QRallInOne.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

drawGraph : drawGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testNDperm : testNDperm.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testOptPart : testOptPart.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testOrderViaBestOfNDandMS : testOrderViaBestOfNDandMS.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testOrderViaMMD : testOrderViaMMD.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testOrderViaMS : testOrderViaMS.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testOrderViaND : testOrderViaND.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testSemi : testSemi.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
