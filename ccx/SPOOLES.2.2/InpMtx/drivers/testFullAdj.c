/*  testFullAdj.c  */

#include "../InpMtx.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   generate a random InpMtx object for a matrix A and
   return an IVL object with the structure of A + A^T

   created -- 97nov05, cca
   ---------------------------------------------------
*/
{
Drand    *drand ;
int      msglvl, ient, ii, irow, isHere, jcol, nent, nedgesMissing,
         nvtx, seed, size ;
int      *colids, *list, *rowids ;
InpMtx   *inpmtxA ;
FILE     *msgFile ;
IVL      *adjIVL ;

if ( argc != 6 ) {
   fprintf(stdout, 
      "\n\n usage : testFullAdj msglvl msgFile nvtx nent seed"
      "\n    msglvl  -- message level"
      "\n    msgFile -- message file"
      "\n    nvtx    -- number of rows and columns"
      "\n    nent    -- bound on number of entries"
      "\n    seed    -- random number seed"
      "\n") ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
nvtx = atoi(argv[3]) ;
nent = atoi(argv[4]) ;
seed = atoi(argv[5]) ;
fprintf(msgFile, 
        "\n testIO "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n nvtx     -- %d" 
        "\n nent     -- %d" 
        "\n seed     -- %d" 
        "\n",
        msglvl, argv[2], nvtx, nent, seed) ;
fflush(msgFile) ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
Drand_setUniform(drand, 0, nvtx) ;
/*
   -----------------------------------
   initialize the InpMtx object for A
   -----------------------------------
*/
inpmtxA = InpMtx_new() ;
InpMtx_init(inpmtxA, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, nent, nvtx) ;
/*
   ----------------------
   load with random edges
   ----------------------
*/
rowids = IVinit(nent, -1) ;
colids = IVinit(nent, -1) ;
Drand_fillIvector(drand, nent, rowids) ;
Drand_fillIvector(drand, nent, colids) ;
for ( ient = 0 ; ient < nent ; ient++ ) {
   irow = rowids[ient] ;
   jcol = colids[ient] ;
   if ( msglvl > 0 ) {
      fprintf(msgFile, "\n loading (%5d,%5d)", irow, jcol) ;
      fflush(msgFile) ;
   }
   InpMtx_inputEntry(inpmtxA, irow, jcol) ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after loading raw data") ;
   InpMtx_writeForHumanEye(inpmtxA, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   sort, compress and change to vector form
   ----------------------------------------
*/
InpMtx_sortAndCompress(inpmtxA) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after sort and compress") ;
   InpMtx_writeForHumanEye(inpmtxA, msgFile) ;
   fflush(msgFile) ;
}
InpMtx_convertToVectors(inpmtxA) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after convert to vectors") ;
   InpMtx_writeForHumanEye(inpmtxA, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   get the full adjacency structure of A + A^T
   -------------------------------------------
*/
adjIVL = InpMtx_fullAdjacency(inpmtxA) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n full adjacency IVL object") ;
   IVL_writeForHumanEye(adjIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------
   check that each (irow,jcol) is in the full adjacency
   ----------------------------------------------------
*/
for ( ient = 0, nedgesMissing = 0 ; ient < nent ; ient++ ) {
   irow = rowids[ient] ;
   jcol = colids[ient] ;
   IVL_listAndSize(adjIVL, irow, &size, &list) ;
   for ( ii = 0, isHere = 0 ; ii < size ; ii++ ) {
      if ( list[ii] == jcol ) {
         isHere = 1 ;
         break ;
      }
   }
   if ( isHere != 1 ) {
      fprintf(stderr, "\n fatal error, (%d,%d) not in adjIVL",
              irow, jcol) ;
      nedgesMissing++ ;
   }
   IVL_listAndSize(adjIVL, jcol, &size, &list) ;
   for ( ii = 0, isHere = 0 ; ii < size ; ii++ ) {
      if ( list[ii] == irow ) {
         isHere = 1 ;
         break ;
      }
   }
   if ( isHere != 1 ) {
      fprintf(stderr, "\n fatal error, (%d,%d) not in adjIVL",
              irow, jcol) ;
      nedgesMissing++ ;
   }
}
fprintf(msgFile, "\n %d edges missing from adjIVL", nedgesMissing) ;
/*
   -------------
   free the data
   -------------
*/
IVfree(colids) ;
IVfree(rowids) ;
InpMtx_free(inpmtxA) ;
IVL_free(adjIVL) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
