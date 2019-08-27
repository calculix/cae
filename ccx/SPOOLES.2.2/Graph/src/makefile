include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Graph
 
$(OBJ).a : \
   $(OBJ).a(IO.o) \
   $(OBJ).a(basics.o) \
   $(OBJ).a(compress.o) \
   $(OBJ).a(equivMap.o) \
   $(OBJ).a(expand.o) \
   $(OBJ).a(fillFromOffsets.o) \
   $(OBJ).a(init.o) \
   $(OBJ).a(setListsFromOffsets.o) \
   $(OBJ).a(util.o) \
   $(OBJ).a(wirebasket.o)
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

