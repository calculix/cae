include ../../Make.inc

#-----------------------------------------------------------------------
 
OBJ = Bridge
 
$(OBJ).a : \
      $(OBJ).a(Cleanup.o) \
      $(OBJ).a(Factor.o) \
      $(OBJ).a(MatMul.o) \
      $(OBJ).a(Setup.o) \
      $(OBJ).a(Solve.o)
	$(RANLIB) $(OBJ).a
	@echo $(OBJ).a is now up-to-date
 
clean :
	- rm -f *.a *.o

