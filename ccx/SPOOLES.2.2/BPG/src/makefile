include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = BPG
 
$(OBJ).a : \
   $(OBJ).a(basics.o) \
   $(OBJ).a(DM.o) \
   $(OBJ).a(init.o) \
   $(OBJ).a(IO.o) \
   $(OBJ).a(makeGraphs.o) \
   $(OBJ).a(maxFlow.o) \
   $(OBJ).a(pseudo.o)
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

