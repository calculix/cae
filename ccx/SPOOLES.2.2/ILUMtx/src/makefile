include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = ILUMtx
 
$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(factor.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(IO.o) \
      $(OBJ).a(misc.o) \
      $(OBJ).a(solve.o)
	@ $(RANLIB) $(OBJ).a
	@ echo $(OBJ).a is now up-to-date
 
clean :
	- rm -f *.a *.o
