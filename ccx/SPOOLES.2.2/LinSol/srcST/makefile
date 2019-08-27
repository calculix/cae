include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Bridge

$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(factor.o) \
      $(OBJ).a(info.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(setparams.o) \
      $(OBJ).a(setup.o) \
      $(OBJ).a(solve.o)
	$(RANLIB) $(OBJ).a
	@echo $(OBJ).a is now up-to-date

clean :
	- rm -f *.a *.o
