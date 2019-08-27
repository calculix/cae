include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm

DRIVERS = test_addChevron \
          test_assmbChv \
          test_copyEntriesToVector\
          test_copyBigEntriesToVector\
          test_factor \
          test_findPivot \
          test_maxabs \
          test_swap \
          test_r1upd \
          test_r2upd \
          test_update
          
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
test_assmbChv : test_assmbChv.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_factor : test_factor.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_findPivot : test_findPivot.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_maxabs : test_maxabs.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_r2upd : test_r2upd.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_r1upd : test_r1upd.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_addChevron : test_addChevron.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_copyBigEntriesToVector : test_copyBigEntriesToVector.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}

test_copyEntriesToVector : test_copyEntriesToVector.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_swap : test_swap.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}
 
test_update : test_update.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ ${PURIFY_GCC_VERSION} ${LIBS}

