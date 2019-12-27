include ../../Make.inc

#-----------------------------------------------------------------------
 
OBJ = PatchAndGoInfo
 
$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(init.o)
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

