include ../../Make.inc

#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm
 
DRIVERS = test_sort test_sortUpAndCompress

drivers : ${DRIVERS}
 
clean :
	- rm -f $.a *.o ${DRIVERS}
 
test_sort : test_sort.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
 
test_sortUpAndCompress : test_sortUpAndCompress.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
