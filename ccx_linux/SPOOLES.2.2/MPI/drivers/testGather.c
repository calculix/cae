/*  testGather.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------------------------
   this program tests the DenseMtx_MPI_gatherRows() method.

   we create a distributed DenseMtx object X, 
   which is partitioned among the processors.
   the entries of X have a particular form so we can check our results.

   each processor then defines a set of rows for its Y matrix,
   and figures out where those rows come from, and puts this
   information into a recvIVL object.

   using the IVL_MPI_alltoall() method, the processors then
   construct their individual sendIVL objects that tell them
   which rows of their owned X matrix they need to send to
   the other processors.

   the lists in sendIVL are then made local w.r.t. X.
   the lists in recvIVL are then made local w.r.t. Y.

   the processors then fill the entries in their Y matrices
   using the DenseMtx_MPI_gatherRows() method.

   the processors then check that their entries in Y are correct.

   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
{
char        *buffer ;
DenseMtx    *X, *Y ;
double      error, gerror, imag, real, t1, t2, value ;
Drand       drand ;
int         count, iirow, iproc, irow, inc1, inc2, jcol, length, myid, 
            msglvl, ncol, ncolX, ncolY, nproc, nrow, nrowX, nrowY, 
            seed, tag, size, type ;
int         *colindX, *colindY, *counts, *head, *link, *list, *map, 
            *rowidsY, *rowindX, *rowindY, *temp ;
int         stats[4], tstats[4] ;
IV          *mapIV ;
IVL         *recvIVL, *sendIVL ;
FILE        *msgFile ;
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
if ( argc != 9 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow ncol inc1 inc2 seed "
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- type of entries"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow    -- number of rows"
"\n    ncol    -- number of columns"
"\n    inc1    -- row increment"
"\n    inc2    -- column increment"
"\n    seed    -- random number seed"
"\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else {
   length = strlen(argv[2]) + 1 + 4 ;
   buffer = CVinit(length, '\0') ;
   sprintf(buffer, "%s.%d", argv[2], myid) ;
   if ( (msgFile = fopen(buffer, "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], argv[2]) ;
      return(0) ;
   }
   CVfree(buffer) ;
}
type = atoi(argv[3]) ;
nrow = atoi(argv[4]) ;
ncol = atoi(argv[5]) ;
inc1 = atoi(argv[6]) ;
inc2 = atoi(argv[7]) ;
seed = atoi(argv[8]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl     -- %d" 
        "\n msgFile    -- %s" 
        "\n type       -- %d" 
        "\n nrow       -- %d" 
        "\n ncol       -- %d" 
        "\n inc1       -- %d" 
        "\n inc2       -- %d" 
        "\n seed       -- %d" 
        "\n",
        argv[0], msglvl, argv[2], type, nrow, ncol, 
        inc1, inc2, seed) ;
fflush(msgFile) ;
if ( inc1 != 1 && inc2 != 1 ) {
   fprintf(stderr, "\n fatal error, inc1 = %d, inc2 = %d", inc1, inc2) ;
   exit(-1) ;
}
/*
   ----------------------------------
   set up the random number generator
   ----------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_setSeed(&drand, seed) ;
/*
   -----------------------------------------------------------
   generate the mapIV object that maps rows of X to processors
   -----------------------------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, nrow, NULL) ;
map = IV_entries(mapIV) ;
Drand_setUniform(&drand, 0, nproc) ;
Drand_fillIvector(&drand, nrow, map) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n mapIV") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
}
/*
   ------------------------------------------------
   make sure that each processor has some rows of X
   ------------------------------------------------
*/
counts = IVinit(nproc, 0) ;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   counts[map[irow]]++ ;
}
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( counts[iproc] == 0 ) {
      fprintf(stderr, "\n proc %d has no rows", iproc) ;
      MPI_Finalize() ;
      exit(-1) ;
   }
}
nrowX = counts[myid] ;
fprintf(msgFile, "\n nrowX = %d", nrowX) ;
IVfree(counts) ;
/*
   -----------------------------------------------------
   generate this processor's portion of the dense matrix
   -----------------------------------------------------
*/
X = DenseMtx_new() ;
if ( inc1 == 1 ) {
   DenseMtx_init(X, type, 1, -1, nrowX, ncol, 1, nrowX) ;
} else if ( inc2 == 1 ) {
   DenseMtx_init(X, type, 1, -1, nrowX, ncol, ncol, 1) ;
}
DenseMtx_zero(X) ;
DenseMtx_rowIndices(X, &nrowX, &rowindX) ;
for ( irow = iirow = 0 ; irow < nrow ; irow++ ) {
   if ( map[irow] == myid ) {
      rowindX[iirow++] = irow ;
   }
}
DenseMtx_columnIndices(X, &ncolX, &colindX) ;
IVramp(ncolX, colindX, 0, 1) ;
if ( DENSEMTX_IS_REAL(X) ) {
   for ( iirow = 0 ; iirow < nrowX ; iirow++ ) {
      irow = rowindX[iirow] ;
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         DenseMtx_setRealEntry(X, iirow, jcol, irow + nrow*jcol) ;
      }
   }
} else {
   for ( iirow = 0 ; iirow < nrowX ; iirow++ ) {
      irow = rowindX[iirow] ;
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         DenseMtx_setComplexEntry(X, iirow, jcol, 
                             irow + nrow*jcol, 2*(irow + nrow*jcol)) ;
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n X") ;
   DenseMtx_writeForHumanEye(X, msgFile) ;
}
/*
   ---------------------------------
   generate a vector of rowids for Y
   ---------------------------------
*/
{ int ii ; double value ;
for ( ii = 0 ; ii <= myid ; ii++ ) {
   value = Drand_value(&drand) ;
}
}
Drand_setUniform(&drand, 1, nrow+1) ;
nrowY = (int) Drand_value(&drand) ;
fprintf(msgFile, "\n nrowY = %d", nrowY) ;
rowidsY = IVinit(nrowY, -1) ;
temp    = IVinit(nrow, -1) ;
IVramp(nrow, temp, 0, 1) ;
IVshuffle(nrow, temp, seed + 2*myid + 1) ;
IVcopy(nrowY, rowidsY, temp) ;
IVqsortUp(nrowY, rowidsY) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rowids for Y") ;
   IVfprintf(msgFile, nrowY, rowidsY) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   set up the recvIVL object, list iproc contains 
   the rows this processor needs from iproc
   ----------------------------------------------
*/
recvIVL = IVL_new() ;
IVL_init1(recvIVL, IVL_CHUNKED, nproc) ;
head = IVinit(nproc, -1) ;
link = IVinit(nrow,  -1) ;
for ( iirow = 0 ; iirow < nrowY ; iirow++ ) {
   irow = rowidsY[iirow] ;
   iproc = map[irow] ;
   link[irow] = head[iproc] ;
   head[iproc] = irow ;
}
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   count = 0 ;
   for ( irow = head[iproc] ; irow != -1 ; irow = link[irow] ) {
      temp[count++] = irow ;
   }
   IVqsortUp(count, temp) ;
   IVL_setList(recvIVL, iproc, count, temp) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n recvIVL") ;
   IVL_writeForHumanEye(recvIVL, msgFile) ;
}
/*
   ----------------------------------------------
   set up the sendIVL object, list iproc contains
   the rows this processor will send to iproc
   ----------------------------------------------
*/
stats[0] = stats[1] = stats[2] = stats[3] = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
tag = 1 ;
MARKTIME(t1) ;
sendIVL = IVL_MPI_alltoall(recvIVL, NULL, 
                          stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : IVL_MPI_alltoall", t2 - t1) ;
fprintf(msgFile, "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n sendIVL") ;
   IVL_writeForHumanEye(sendIVL, msgFile) ;
}
/*
   ----------------------------------------
   make the lists in sendIVL local w.r.t. X
   ----------------------------------------
*/
for ( irow = 0 ; irow < nrowX ; irow++ ) {
   temp[rowindX[irow]] = irow ;
}
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   IVL_listAndSize(sendIVL, iproc, &size, &list) ;
   for ( iirow = 0 ; iirow < size ; iirow++ ) {
      list[iirow] = temp[list[iirow]] ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n sendIVL with local lists") ;
   IVL_writeForHumanEye(sendIVL, msgFile) ;
}
/*
   ----------------------------------------
   make the lists in recvIVL local w.r.t. Y
   ----------------------------------------
*/
for ( irow = 0 ; irow < nrowY ; irow++ ) {
   temp[rowidsY[irow]] = irow ;
}
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   IVL_listAndSize(recvIVL, iproc, &size, &list) ;
   for ( iirow = 0 ; iirow < size ; iirow++ ) {
      list[iirow] = temp[list[iirow]] ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n recvIVL with local lists") ;
   IVL_writeForHumanEye(recvIVL, msgFile) ;
}
/*
   -------------------
   create the Y matrix
   -------------------
*/
Y = DenseMtx_new() ;
if ( inc1 == 1 ) {
   DenseMtx_init(Y, type, 1, -1, nrowY, ncol, 1, nrowY) ;
} else {
   DenseMtx_init(Y, type, 1, -1, nrowY, ncol, ncol, 1) ;
}
DenseMtx_zero(Y) ;
DenseMtx_rowIndices(Y, &nrowY, &rowindY) ;
IVcopy(nrowY, rowindY, rowidsY) ;
DenseMtx_columnIndices(Y, &ncolY, &colindY) ;
IVramp(ncolY, colindY, 0, 1) ;
/*
   ----------------------------------------------
   gather the entries of Y from the distributed X
   ----------------------------------------------
*/
stats[0] = stats[1] = stats[2] = stats[3] = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
tag = 1 ;
MARKTIME(t1) ;
DenseMtx_MPI_gatherRows(Y, X, sendIVL, recvIVL, stats, msglvl,
                        msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : DenseMtx_MPI_gatherRows", t2 - t1) ;
fprintf(msgFile, "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Y") ;
   DenseMtx_writeForHumanEye(Y, msgFile) ;
}
/*
   ---------------------------------
   check that the local Y is correct
   ---------------------------------
*/
if ( DENSEMTX_IS_REAL(Y) ) {
   for ( iirow = 0, error = 0.0 ; iirow < nrowY ; iirow++ ) {
      irow = rowindY[iirow] ;
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         DenseMtx_realEntry(Y, iirow, jcol, &real) ;
         value = real - (irow + nrow*jcol) ;
         error += value*value ;
      }
   }
} else {
   for ( iirow = 0, error = 0.0 ; iirow < nrowY ; iirow++ ) {
      irow = rowindY[iirow] ;
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         DenseMtx_complexEntry(Y, iirow, jcol, &real, &imag) ;
         value = real - (irow + nrow*jcol) ;
         error += value*value ;
         value = imag - 2*(irow + nrow*jcol) ;
         error += value*value ;
      }
   }
}
fprintf(msgFile, "\n\n local error  = %12.4e", error) ;
MPI_Reduce((void *) &error, (void *) &gerror, 1, MPI_DOUBLE,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global error = %12.4e", gerror) ;
   fflush(msgFile) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
IV_free(mapIV) ;
DenseMtx_free(X) ;
DenseMtx_free(Y) ;
IVfree(rowidsY) ;
IVfree(temp) ;
IVfree(head) ;
IVfree(link) ;
IVL_free(recvIVL) ;
IVL_free(sendIVL) ;

MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
