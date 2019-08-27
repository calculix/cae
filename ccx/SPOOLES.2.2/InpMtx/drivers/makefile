include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = adjToGraph \
          createGraph \
          createGraphForATA \
          getProfile \
          mkLaplacianMtx \
          mkNaturalFactorMtx \
          readAIJ \
          readAIJ2 \
          testExtract \
          testFullAdj \
          testFullAdj2 \
          testGMMM \
          testGMVM \
          testHBIO \
          testHBIO2 \
          testInput \
          testIO \
          testMMM \
          testR2D \
          weightedAdjToGraph
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

adjToGraph : adjToGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

createGraph : createGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

createGraphForATA : createGraphForATA.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

getProfile : getProfile.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

mkLaplacianMtx : mkLaplacianMtx.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

mkNaturalFactorMtx : mkNaturalFactorMtx.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

readAIJ : readAIJ.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

readAIJ2 : readAIJ2.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testExtract : testExtract.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testFullAdj : testFullAdj.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testFullAdj2 : testFullAdj2.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testGMMM : testGMMM.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testGMVM : testGMVM.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testHBIO : testHBIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testHBIO2 : testHBIO2.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testInput : testInput.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testIO : testIO.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testMMM : testMMM.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

testR2D : testR2D.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

weightedAdjToGraph : weightedAdjToGraph.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
