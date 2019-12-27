include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = InpMtx
 
$(OBJ).a : \
      $(OBJ).a(adjForATA.o) \
      $(OBJ).a(basics.o) \
      $(OBJ).a(extract.o) \
      $(OBJ).a(fullAdj.o) \
      $(OBJ).a(gmmm.o) \
      $(OBJ).a(gmvm.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(input.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(IO.o) \
      $(OBJ).a(map.o) \
      $(OBJ).a(mvm.o) \
      $(OBJ).a(mvmVector.o) \
      $(OBJ).a(permute.o) \
      $(OBJ).a(profile.o) \
      $(OBJ).a(support.o) \
      $(OBJ).a(util.o)
	@ $(RANLIB) $(OBJ).a
	@ echo $(OBJ).a is now up-to-date
 
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

