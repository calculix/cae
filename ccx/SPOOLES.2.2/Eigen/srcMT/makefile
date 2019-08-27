include ../../Make.inc

#-----------------------------------------------------------------------
 
OBJ = BridgeMT
 
$(OBJ).a : \
      $(OBJ).a(CleanupMT.o) \
      $(OBJ).a(FactorMT.o) \
      $(OBJ).a(MatMulMT.o) \
      $(OBJ).a(SetupMT.o) \
      $(OBJ).a(SolveMT.o)
	$(RANLIB) $(OBJ).a
	@echo $(OBJ).a is now up-to-date
 
clean :
	- rm -f *.a *.o

