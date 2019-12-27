include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = GPart
 
$(OBJ).a : \
      $(OBJ).a(DDviaFishnet.o) \
      $(OBJ).a(DDviaProjection.o) \
      $(OBJ).a(DDsepInfo.o) \
      $(OBJ).a(RBviaDDsep.o) \
      $(OBJ).a(TwoSetViaBKL.o) \
      $(OBJ).a(basics.o) \
      $(OBJ).a(domSegMap.o) \
      $(OBJ).a(identifyWideSep.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(makeYCmap.o) \
      $(OBJ).a(smoothBisector.o) \
      $(OBJ).a(smoothBy2layers.o) \
      $(OBJ).a(smoothYSep.o) \
      $(OBJ).a(split.o) \
      $(OBJ).a(util.o)
	$(RANLIB) $(OBJ).a
	@echo $(OBJ).a is now up-to-date
 
clean :
	- rm -f *.a *.o

makeLib :
	perl ../../makeLib > makeG
	make -f makeG
	rm -f makeG

updateLib :
	perl ../../updLib > makeG
	make -f makeG
	rm -f makeG

