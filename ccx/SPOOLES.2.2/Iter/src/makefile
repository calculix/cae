include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = Iter

$(OBJ).a : \
	$(OBJ).a(bicgstabr.o) \
	$(OBJ).a(bicgstabl.o) \
	$(OBJ).a(mlbicgstabr.o) \
	$(OBJ).a(mlbicgstabl.o) \
	$(OBJ).a(pcgr.o)\
	$(OBJ).a(pcgl.o)\
	$(OBJ).a(tfqmrr.o)\
	$(OBJ).a(tfqmrl.o)\
	$(OBJ).a(zbicgstabr.o) \
	$(OBJ).a(zbicgstabl.o) \
	$(OBJ).a(zmlbicgstabr.o) \
	$(OBJ).a(zmlbicgstabl.o) \
	$(OBJ).a(zpcgl.o)\
	$(OBJ).a(zpcgr.o)\
	$(OBJ).a(bgmresl.o)\
	$(OBJ).a(bgmresr.o)\
	$(OBJ).a(ztfqmrr.o)\
	$(OBJ).a(ztfqmrl.o)\
	$(OBJ).a(DenseMtx_mmm.o) \
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
