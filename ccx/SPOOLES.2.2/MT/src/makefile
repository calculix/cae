include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = spoolesMT
 
$(OBJ).a : \
      $(OBJ).a(factorMT.o) \
      $(OBJ).a(mvm.o) \
      $(OBJ).a(QRfactorMT.o) \
      $(OBJ).a(QRsolveMT.o) \
      $(OBJ).a(solveMT.o)
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

