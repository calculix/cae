include ../../Make.inc

#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = drawTree 
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

drawTree : drawTree.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS} 
