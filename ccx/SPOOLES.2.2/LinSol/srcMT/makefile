include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = BridgeMT

$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(factor.o) \
      $(OBJ).a(factorSetup.o) \
      $(OBJ).a(info.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(setparams.o) \
      $(OBJ).a(setup.o) \
      $(OBJ).a(solve.o) \
      $(OBJ).a(solveSetup.o)
	$(RANLIB) $(OBJ).a
	@echo $(OBJ).a is now up-to-date

clean :
	- rm -f *.a *.o
