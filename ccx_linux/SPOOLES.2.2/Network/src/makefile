include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Network
 
$(OBJ).a : \
      $(OBJ).a(addArc.o) \
      $(OBJ).a(augmentPath.o) \
      $(OBJ).a(basics.o) \
      $(OBJ).a(findAugmentingPath.o) \
      $(OBJ).a(findMaxFlow.o) \
      $(OBJ).a(findMincut.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(IO.o)
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

