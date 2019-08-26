include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = testIO \
          mkGraph \
          mkGridEGraph
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

mkGraph : mkGraph.o  ../../spooles.a
	$(PURIFY) $(CC) $@.o -o $@ $(PURIFY_GCC_VERSION) $(LIBS)

mkGridEGraph : mkGridEGraph.o ../../spooles.a
	$(PURIFY) $(CC) $@.o -o $@ $(PURIFY_GCC_VERSION)$(LIBS)

testIO : testIO.o ../../spooles.a
	$(PURIFY) $(CC) $@.o -o $@ $(PURIFY_GCC_VERSION) $(LIBS)
