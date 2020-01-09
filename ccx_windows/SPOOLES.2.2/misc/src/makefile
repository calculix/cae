include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = misc
 
$(OBJ).a : \
      $(OBJ).a(drawGraphEPS.o) \
      $(OBJ).a(fpgrids.o) \
      $(OBJ).a(localND.o) \
      $(OBJ).a(ND.o) \
      $(OBJ).a(ND2.o) \
      $(OBJ).a(mkNDlinsys.o) \
      $(OBJ).a(mkNDlinsysQR.o) \
      $(OBJ).a(orderViaBestOfNDandMS.o) \
      $(OBJ).a(orderViaMMD.o) \
      $(OBJ).a(orderViaMS.o) \
      $(OBJ).a(orderViaND.o)
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

