include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Chv

$(OBJ).a : \
      $(OBJ).a(assemble.o) \
      $(OBJ).a(basics.o) \
      $(OBJ).a(copy.o) \
      $(OBJ).a(factor.o) \
      $(OBJ).a(findPivot.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(IO.o) \
      $(OBJ).a(search.o) \
      $(OBJ).a(swap.o) \
      $(OBJ).a(update.o) \
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

