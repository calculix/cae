include ../../Make.inc
 
#-----------------------------------------------------------------------
 
LIBS = ../../spooles.a -lm
 
DRIVERS = testIO mk9PCoords
 
drivers : ${DRIVERS}
 
clean :
	- rm -f *.a *.o ${DRIVERS}
 
mk9PCoords : mk9PCoords.o ../../spooles.a
	${PURIFY} $(CC) $@.o -o $@ $(PURIFY_GCC_VERSION) $(LIBS)
 
testIO : testIO.o ../../spooles.a
	${PURIFY} $(CC) $@.o -o $@ $(PURIFY_GCC_VERSION) $(LIBS)

