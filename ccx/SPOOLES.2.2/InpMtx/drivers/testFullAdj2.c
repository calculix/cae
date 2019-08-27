/*  testFullAdj2.c  */

#include "../InpMtx.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------------
   generate random InpMtx object for matrices A and B
   return an IVL object with the structure of (A+B) + (A+B)^T

   created -- 97nov05, cca
   ----------------------------------------------------------
*/
{
Drand    *drand ;
int      msglvl, ient, ii, irow, isHere, jcol, nedgesMissing,
         nentA, nentB, nvtx, seed, size ;
int      *colidsA, *colidsB, *list, *rowidsA, *rowidsB ;
InpMtx   *inpmtxA, *inpmtxB ;
FILE     *msgFile ;
IVL      *adjIVL ;

if ( argc != 7 ) {
   fprintf(stdout, 
      "\n\n usage : testFullAdj2 msglvl msgFile nvtx nentA nentB seed"
      "\n    msglvl  -- message level"
      "\n    msgFile -- message file"
      "\n    nvtx    -- number of rows and columns"
      "\n    nentA   -- bound on number of entries in A"
      "\n    nentB   -- bound on number of entries in B"
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
nvtx  = atoi(argv[3]) ;
nentA = atoi(argv[4]) ;
nentB = atoi(argv[5]) ;
seed  = atoi(argv[6]) ;
fprintf(msgFile, 
        "\n testIO "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n nvtx     -- %d" 
        "\n nentA    -- %d" 
        "\n nentB    -- %d" 
        "\n seed     -- %d" 
        "\n",
        msglvl, argv[2], nvtx, nentA, nentB, seed) ;
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
InpMtx_init(inpmtxA, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, nentA, nvtx) ;
/*
   ----------------------
   load with random edges
   ----------------------
*/
rowidsA = IVinit(nentA, -1) ;
colidsA = IVinit(nentA, -1) ;
Drand_fillIvector(drand, nentA, rowidsA) ;
Drand_fillIvector(drand, nentA, colidsA) ;
for ( ient = 0 ; ient < nentA ; ient++ ) {
   irow = rowidsA[ient] ;
   jcol = colidsA[ient] ;
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
   -----------------------------------
   initialize the InpMtx object for B
   -----------------------------------
*/
inpmtxB = InpMtx_new() ;
InpMtx_init(inpmtxB, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, nentB, nvtx) ;
/*
   ----------------------
   load with random edges
   ----------------------
*/
rowidsB = IVinit(nentB, -1) ;
colidsB = IVinit(nentB, -1) ;
Drand_fillIvector(drand, nentB, rowidsB) ;
Drand_fillIvector(drand, nentB, colidsB) ;
for ( ient = 0 ; ient < nentB ; ient++ ) {
   irow = rowidsB[ient] ;
   jcol = colidsB[ient] ;
   if ( msglvl > 0 ) {
      fprintf(msgFile, "\n loading (%5d,%5d)", irow, jcol) ;
      fflush(msgFile) ;
   }
   InpMtx_inputEntry(inpmtxB, irow, jcol) ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after loading raw data") ;
   InpMtx_writeForHumanEye(inpmtxB, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   sort, compress and change to vector form
   ----------------------------------------
*/
InpMtx_sortAndCompress(inpmtxB) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after sort and compress") ;
   InpMtx_writeForHumanEye(inpmtxB, msgFile) ;
   fflush(msgFile) ;
}
InpMtx_convertToVectors(inpmtxB) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after convert to vectors") ;
   InpMtx_writeForHumanEye(inpmtxB, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------------
   get the full adjacency structure of (A+B) + (A+B)^T
   ---------------------------------------------------
*/
fprintf(msgFile, "\n inpmtxA = %p, inpmtxB = %p", inpmtxA, inpmtxB) ;
adjIVL = InpMtx_fullAdjacency2(inpmtxA, inpmtxB) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n full adjacency IVL object") ;
   IVL_writeForHumanEye(adjIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   check that each (irow,jcol) from A is in the full adjacency
   -----------------------------------------------------------
*/
for ( ient = 0, nedgesMissing = 0 ; ient < nentA ; ient++ ) {
   irow = rowidsA[ient] ;
   jcol = colidsA[ient] ;
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
fprintf(msgFile, "\n %d edges of A missing from adjIVL", 
        nedgesMissing) ;
/*
   -----------------------------------------------------------
   check that each (irow,jcol) from B is in the full adjacency
   -----------------------------------------------------------
*/
for ( ient = 0, nedgesMissing = 0 ; ient < nentB ; ient++ ) {
   irow = rowidsB[ient] ;
   jcol = colidsB[ient] ;
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
fprintf(msgFile, "\n %d edges of B missing from adjIVL", 
        nedgesMissing) ;
/*
   -------------
   free the data
   -------------
*/
IVfree(colidsA) ;
IVfree(rowidsA) ;
InpMtx_free(inpmtxA) ;
IVfree(colidsB) ;
IVfree(rowidsB) ;
InpMtx_free(inpmtxB) ;
IVL_free(adjIVL) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
