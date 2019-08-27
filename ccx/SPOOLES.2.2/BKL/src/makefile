include ../../Make.inc
 
OBJ = BKL
 
$(OBJ).a : \
           $(OBJ).a(basics.o) \
           $(OBJ).a(evalfcn.o) \
           $(OBJ).a(exhSearch.o) \
           $(OBJ).a(fidmat.o) \
           $(OBJ).a(init.o) \
           $(OBJ).a(util.o)
	$(RANLIB) $(OBJ).a
	@echo lib $(OBJ).a is now up-to-date
 
 
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

