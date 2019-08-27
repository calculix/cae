include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS =	../src/Iter.a ${THREAD_LIBS} ../../spooles.a -lm
LIBS =	${THREAD_LIBS} ../../spooles.a -lm


DRIVERS = iter test_colDotProduct test_colGenAxpy test_frobNorm\
          test_twoNormOfColumn test_colCopy test_mmm

drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}

clean_test:
	- rm -f test_colDotProduct test_colGenAxpy test_frobNorm\
	  test_twoNormOfColumn test_colCopy \
	  test_colDotProduct.o test_colGenAxpy.o test_frobNorm.o\
	  test_twoNormOfColumn.o test_colCopy.o 
 
iter    : iter.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

test_colDotProduct:  test_colDotProduct.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

test_colGenAxpy:  test_colGenAxpy.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

test_frobNorm:  test_frobNorm.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

test_twoNormOfColumn:  test_twoNormOfColumn.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

test_colCopy:  test_colCopy.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

test_mmm: test_mmm.o ../../spooles.a
	${PURIFY} ${CC} $(PURIFY_GCC_VERSION) $@.o -o $@ ${LIBS}

