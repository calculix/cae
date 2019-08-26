include ../../Make.inc

#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm
 
DRIVERS = test_hash 
 
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
test_hash : test_hash.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
