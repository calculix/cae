include ../../Make.inc
 
#-----------------------------------------------------------------------
 
OBJ = A2

$(OBJ).a : \
      $(OBJ).a(basics.o) \
      $(OBJ).a(copyEntriesToVector.o) \
      $(OBJ).a(init.o) \
      $(OBJ).a(instance.o) \
      $(OBJ).a(IO.o) \
      $(OBJ).a(makeStaircase.o) \
      $(OBJ).a(norms.o) \
      $(OBJ).a(QRreduce.o) \
      $(OBJ).a(sort.o) \
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

