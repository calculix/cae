include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = \
       ../src/ILUMtx.a \
       ../../InpMtx/src/InpMtx.a \
       ../../DenseMtx/src/DenseMtx.a \
       ../../A2/src/A2.a \
       ../../EGraph/src/EGraph.a \
       ../../ETree/src/ETree.a \
       ../../Graph/src/Graph.a \
       ../../Tree/src/Tree.a \
       ../../Coords/src/Coords.a \
       ../../IVL/src/IVL.a \
       ../../IV/src/IV.a \
       ../../DV/src/DV.a \
       ../../ZV/src/ZV.a \
       ../../Utilities/src/Utilities.a \
       ../../Drand/src/Drand.a

DRIVERS = testFactor testSolve
 
libs :
	cd ../src ; make ILUMtx.a
	cd ../../InpMtx/src ; make InpMtx.a
	cd ../../DenseMtx/src ; make DenseMtx.a
	cd ../../A2/src ; make A2.a
	cd ../../EGraph/src ; make EGraph.a
	cd ../../Graph/src ; make Graph.a
	cd ../../ETree/src ; make ETree.a
	cd ../../Tree/src ; make Tree.a
	cd ../../Coords/src ; make Coords.a
	cd ../../IVL/src ; make IVL.a
	cd ../../IV/src ; make IV.a
	cd ../../DV/src ; make DV.a
	cd ../../ZV/src ; make ZV.a
	cd ../../Utilities/src ; make Utilities.a
	cd ../../Drand/src ; make Drand.a

drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

testFactor : testFactor.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} -lm

testSolve : testSolve.o libs ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} -lm
