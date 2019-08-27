include ../../Make.inc

#-----------------------------------------------------------------------
 
OBJ = IV
 
$(OBJ).a : \
      $(OBJ).a(IO.o) \
      $(OBJ).a(basics.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(instance.o) \
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

