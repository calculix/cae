include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = SubMtx

$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(initRandom.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(IO.o) \
      $(OBJ).a(scalevec.o) \
      $(OBJ).a(solve.o) \
      $(OBJ).a(solveH.o) \
      $(OBJ).a(solveT.o) \
      $(OBJ).a(solveupd.o) \
      $(OBJ).a(solveupdH.o) \
      $(OBJ).a(solveupdT.o) \
      $(OBJ).a(sort.o) \
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

