include ../../Make.inc

#-----------------------------------------------------------------------
#
#  set suffix rule *.c --> *.a
#
.c.a : 
	$(PURIFY) $(CC) -c $(CFLAGS) $(MPI_INCLUDE_DIR) $*.c -o $*.o
	$(AR) $(ARFLAGS) $(OBJ).a $*.o
	rm $*.o
#
#-----------------------------------------------------------------------

OBJ = spoolesMPI

$(OBJ).a : \
      $(OBJ).a(aggListMPI.o) \
      $(OBJ).a(colmapMPI.o) \
      $(OBJ).a(DenseMtx_gather.o) \
      $(OBJ).a(DenseMtx_scatterAdd.o) \
      $(OBJ).a(ETree_Bcast.o) \
      $(OBJ).a(factorMPI.o) \
      $(OBJ).a(fullAdjMPI.o) \
      $(OBJ).a(Graph_Bcast.o) \
      $(OBJ).a(IVallgather.o) \
      $(OBJ).a(IVLallgather.o) \
      $(OBJ).a(IVL_alltoall.o) \
      $(OBJ).a(IVL_Bcast.o) \
      $(OBJ).a(makeSendRecvIVLs.o) \
      $(OBJ).a(MMM.o) \
      $(OBJ).a(postProcess.o) \
      $(OBJ).a(rowmapMPI.o) \
      $(OBJ).a(solveMPI.o) \
      $(OBJ).a(splitDenseMtx.o) \
      $(OBJ).a(splitFrontMtx.o) \
      $(OBJ).a(splitInpMtx.o) \
      $(OBJ).a(splitPencil.o) \
      $(OBJ).a(symbfacMPI.o) \
      $(OBJ).a(utilities.o)
	@ $(RANLIB) $(OBJ).a
	@ echo $(OBJ).a is now up-to-date

clean :
	- rm -f *.a *.o
