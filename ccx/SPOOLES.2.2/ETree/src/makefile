include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = ETree
 
$(OBJ).a : \
   $(OBJ).a(IO.o) \
   $(OBJ).a(basics.o) \
   $(OBJ).a(compress.o) \
   $(OBJ).a(init.o) \
   $(OBJ).a(initFromSubtree.o) \
   $(OBJ).a(instance.o) \
   $(OBJ).a(justify.o) \
   $(OBJ).a(maps.o) \
   $(OBJ).a(metrics.o) \
   $(OBJ).a(ms.o) \
   $(OBJ).a(permute.o) \
   $(OBJ).a(semi.o) \
   $(OBJ).a(storage.o) \
   $(OBJ).a(transform.o) \
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

