include ../../Make.inc

#-----------------------------------------------------------------------
#
#  set suffix rule *.c --> *.a
#
.c.a : 
	$(PURIFY) $(CC) -c $(CFLAGS) $(MPI_INCLUDE_DIR) $*.c -o $*.o
	$(AR) $(ARFLAGS) $(OBJ).a $*.o
	rm $*.o

#-----------------------------------------------------------------------
 
OBJ = BridgeMPI
 
$(OBJ).a : \
      $(OBJ).a(JimMatMulMPI.o) \
      $(OBJ).a(JimSolveMPI.o) \
      $(OBJ).a(CleanupMPI.o) \
      $(OBJ).a(FactorMPI.o) \
      $(OBJ).a(MatMulMPI.o) \
      $(OBJ).a(SetupMPI.o) \
      $(OBJ).a(SolveMPI.o)
	$(RANLIB) $(OBJ).a
	@echo $(OBJ).a is now up-to-date
 
clean :
	- rm -f *.a *.o

