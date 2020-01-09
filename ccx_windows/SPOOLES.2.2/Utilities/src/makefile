include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Utilities
 
$(OBJ).a : \
   $(OBJ).a(axpy.o) \
   $(OBJ).a(CV.o) \
   $(OBJ).a(DV.o) \
   $(OBJ).a(FV.o) \
   $(OBJ).a(IV.o) \
   $(OBJ).a(PCV.o) \
   $(OBJ).a(PDV.o) \
   $(OBJ).a(PFV.o) \
   $(OBJ).a(PIV.o) \
   $(OBJ).a(IP.o) \
   $(OBJ).a(I2OP.o) \
   $(OBJ).a(iohb.o) \
   $(OBJ).a(newsort.o) \
   $(OBJ).a(sortAndCompress.o) \
   $(OBJ).a(ZV.o)
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

