include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = MSMD
 
$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(cleanReachSet.o) \
      $(OBJ).a(eliminate.o) \
      $(OBJ).a(fillPerms.o) \
      $(OBJ).a(findInodes.o) \
      $(OBJ).a(frontETree.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(makeSchurComplement.o) \
      $(OBJ).a(MSMDinfo.o) \
      $(OBJ).a(MSMDvtx.o) \
      $(OBJ).a(order.o) \
      $(OBJ).a(update.o)
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

