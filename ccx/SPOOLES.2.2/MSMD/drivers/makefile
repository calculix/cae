include ../../Make.inc
 
#-----------------------------------------------------------------------

LIBS = ../../spooles.a -lm

DRIVERS = orderViaMMD \
          orderViaMS \
          orderViaND \
          orderViaStages 
 
drivers : ${DRIVERS}

clean :
	- rm -f *.a *.o ${DRIVERS}

orderViaStages : orderViaStages.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

orderViaMMD : orderViaMMD.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

orderViaND : orderViaND.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}

orderViaMS : orderViaMS.o ../../spooles.a
	${PURIFY} ${CC} $@.o -o $@ $(PURIFY_GCC_VERSION) ${LIBS}
