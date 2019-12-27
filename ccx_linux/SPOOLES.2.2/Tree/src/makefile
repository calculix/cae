include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Tree
 
$(OBJ).a : \
   $(OBJ).a(basics.o) \
   $(OBJ).a(compress.o) \
   $(OBJ).a(draw.o) \
   $(OBJ).a(getCoords.o) \
   $(OBJ).a(init.o) \
   $(OBJ).a(instance.o) \
   $(OBJ).a(IO.o) \
   $(OBJ).a(justify.o) \
   $(OBJ).a(maximizeGain.o) \
   $(OBJ).a(metrics.o) \
   $(OBJ).a(perms.o) \
   $(OBJ).a(permute.o) \
   $(OBJ).a(subtree.o) \
   $(OBJ).a(util.o)
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

