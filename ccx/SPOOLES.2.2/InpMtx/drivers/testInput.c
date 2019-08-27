/*  testInput.c  */

#include "../../timings.h"
#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------
   test the InpMtx input methods for Peter Schartz

   created -- 98sep04, cca
   -----------------------------------------------
*/
{
double   estpar, growth, t1, t2 ;
double   *entries ;
FILE     *msgFile ;
int      count, ii, irow, maxsize, msglvl, nent, neqns, n1, n2, n3,
         size, size1, size2, type ;
int      *indices, *indices1, *indices2, *list ;
InpMtx   *mtxA ;
IVL      *adjIVL, *fullIVL, *lowerIVL ;

if ( argc != 9 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile n1 n2 n3 estpar growth"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    type     -- type of entries"
      "\n       0 -- indices only"
      "\n       1 -- real entries"
      "\n       2 -- complex entries"
      "\n    n1       -- # of grid points in first direction"
      "\n    n2       -- # of grid points in second direction"
      "\n    n3       -- # of grid points in third direction"
      "\n    estpar   -- estimation for nent"
      "\n    growth   -- growth factor"
      "\n", argv[0]) ;
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
type   = atoi(argv[3]) ;
n1     = atoi(argv[4]) ;
n2     = atoi(argv[5]) ;
n3     = atoi(argv[6]) ;
estpar = atof(argv[7]) ;
growth = atof(argv[8]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n type     -- %d" 
        "\n n1       -- %d" 
        "\n n2       -- %d" 
        "\n n3       -- %d" 
        "\n estpar   -- %f" 
        "\n growth   -- %f" 
        "\n",
        argv[0], msglvl, argv[2], type, n1, n2, n3, estpar, growth) ;
fflush(msgFile) ;
if ( n1 <= 0 || n2 <= 0 || n3 <= 0 || estpar < 0.0 || growth <= 1.0 ) {
   fprintf(stderr, "\n fatal error in testInput, bad input\n") ;
   exit(-1) ;
}
/*
   -----------------------------------
   set up the grid adjacency structure
   -----------------------------------
*/
neqns = n1 * n2 * n3 ;
MARKTIME(t1) ;
if ( n1 == 1 ) {
   adjIVL = IVL_make9P(n2, n3, 1) ;
} else if ( n2 == 1 ) {
   adjIVL = IVL_make9P(n1, n3, 1) ;
} else if ( n3 == 1 ) {
   adjIVL = IVL_make9P(n1, n2, 1) ;
} else {
   adjIVL = IVL_make27P(n1, n2, n3, 1) ;
}
nent = IVL_tsize(adjIVL) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : make full adjacency, %d entries", 
        t2 - t1, nent) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n full adjacency structure, %d entries", nent) ;
   IVL_writeForHumanEye(adjIVL, msgFile) ;
}
/*
   ----------------------------------
   make the lower adjacency structure
   ----------------------------------
*/
MARKTIME(t1) ;
lowerIVL = IVL_new() ;
IVL_init1(lowerIVL, IVL_CHUNKED, neqns) ;
list = IVinit(neqns, -1) ;
for ( irow = 0 ; irow < neqns ; irow++ ) {
   IVL_listAndSize(adjIVL, irow, &size, &indices) ;
   for ( ii = count = 0 ; ii < size ; ii++ ) {
      if ( indices[ii] >= irow ) {
         list[count++] = indices[ii] ;
      }
   }
   IVL_setList(lowerIVL, irow, count, list) ;
}
nent = IVL_tsize(lowerIVL) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : make lower adjacency, %d entries", 
        t2 - t1, nent) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n lower adjacency structure, %d entries", nent);
   IVL_writeForHumanEye(adjIVL, msgFile) ;
}
/*
   ---------------------------------------------------
   create a vector to hold entries,
   its size is the maximum size of the lower adjacency
   ---------------------------------------------------
*/
maxsize = IVL_maxListSize(adjIVL) ;
entries = DVinit(2*maxsize, 0.0) ;
/*
   ----------------------------
   initialize the InpMtx object
   ----------------------------
*/
MARKTIME(t1) ;
mtxA = InpMtx_new() ;
InpMtx_init(mtxA, INPMTX_BY_COLUMNS, type, estpar*nent, 0) ;
InpMtx_setResizeMultiple(mtxA, growth) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : initialize InpMtx", t2 - t1) ;

/*
   ----------------
   load the columns
   ----------------
*/
MARKTIME(t1) ;
if ( INPMTX_IS_INDICES_ONLY(mtxA) ) {
   for ( irow = 0 ; irow < neqns ; irow++ ) {
      IVL_listAndSize(lowerIVL, irow, &size, &indices) ;
      InpMtx_inputColumn(mtxA, irow, size, indices) ;
   }
} else if ( INPMTX_IS_REAL_ENTRIES(mtxA) ) {
   for ( irow = 0 ; irow < neqns ; irow++ ) {
      IVL_listAndSize(lowerIVL, irow, &size, &indices) ;
      InpMtx_inputRealColumn(mtxA, irow, size, indices, entries) ;
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(mtxA) ) {
   for ( irow = 0 ; irow < neqns ; irow++ ) {
      IVL_listAndSize(lowerIVL, irow, &size, &indices) ;
      InpMtx_inputComplexColumn(mtxA, irow, size, indices, entries) ;
   }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : load entries by columns", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n mtxA") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
}
/*
   -----------------------------
   sort and compress the entries
   -----------------------------
*/
MARKTIME(t1) ;
InpMtx_sortAndCompress(mtxA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : sort and compress", t2 - t1) ;
/*
   -------------------
   set the vector mode
   -------------------
*/
MARKTIME(t1) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : convert to vectors", t2 - t1) ;
/*
   --------------------------------------
   construct the full adjacency structure
   --------------------------------------
*/
MARKTIME(t1) ;
fullIVL = InpMtx_fullAdjacency(mtxA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : construct the full adjacency", 
        t2 - t1) ;
/*
   -----------------------------------------
   compare the two full adjacency structures
   -----------------------------------------
*/
for ( irow = 0 ; irow < neqns ; irow++ ) {
   IVL_listAndSize(adjIVL,  irow, &size1, &indices1) ;
   IVL_listAndSize(fullIVL, irow, &size2, &indices2) ;
   if ( size1 != size2 ) {
      fprintf(msgFile, "\n\n error, irow %d, size1 %d, size2 %d",
              irow, size1, size2) ;
      exit(-1) ;
   }
   for ( ii = 0 ; ii < size1 ; ii++ ) {
      if ( indices1[ii] != indices2[ii] ) {
         fprintf(msgFile, "\n\n error, irow %d", irow) ;
         fprintf(msgFile, "\n indices1") ;
         IVfprintf(msgFile, size1, indices1) ;
         fprintf(msgFile, "\n indices2") ;
         IVfprintf(msgFile, size1, indices2) ;
         exit(-1) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVL_free(adjIVL) ;
IVL_free(lowerIVL) ;
IVL_free(fullIVL) ;
InpMtx_free(mtxA) ;
DVfree(entries) ;
IVfree(list) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
