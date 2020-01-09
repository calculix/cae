include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = FrontMtx
 

$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(factor.o) \
      $(OBJ).a(factorUtil.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(IO.o) \
      $(OBJ).a(loadEntries.o) \
      $(OBJ).a(permute.o) \
      $(OBJ).a(postponed.o) \
      $(OBJ).a(postProcess.o) \
      $(OBJ).a(QRfactor.o) \
      $(OBJ).a(QRsolve.o) \
      $(OBJ).a(QRutil.o) \
      $(OBJ).a(solve.o) \
      $(OBJ).a(solveUtil.o) \
      $(OBJ).a(split.o) \
      $(OBJ).a(storeFront.o) \
      $(OBJ).a(update.o) \
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

