/*  mkNaturalFactorMtx.c  */

#include "../InpMtx.h"
#include "../../EGraph.h"
#include "../../Coords.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   create a natural factor matrix.
   for a n1 x n2 grid there are (n1-1)*(n2-1) elements
      and four rows per element, so the matrix is
      4*(n1-1)*(n2-1) x n1*n2
   for a n1 x n2 x n3 grid there are (n1-1)*(n2-1)*(n3-1) elements
      and eight rows per element, so the matrix is
      8*(n1-1)*(n2-1)*(n3-1) x n1*n2*n3

   created -- 97sep19, cca
   ---------------------------------------------------------------
*/
void
main ( int argc, char *argv[] ) 
{
char       *outFileName ;
InpMtx    *inpmtx ;
double     *entries ;
double     t1, t2 ;
Drand      *drand ;
EGraph     *egraph ;
FILE       *msgFile ;
int        ielem, irow, jrow, msglvl, nelem, nvtx, n1, n2, n3, seed, size ;
int        *indices ;
/*
   ---------------
   check the input
   ---------------
*/
if ( argc != 8 ) {
   fprintf(stdout, 
        "\n\n usage : %s msglvl msgFile n1 n2 n3 seed outFile "
        "\n    msglvl  -- message level"
        "\n    msgFile -- message file"
        "\n    n1      -- number of grid points in the first direction"
        "\n    n2      -- number of grid points in the second direction"
        "\n    n3      -- number of grid points in the third direction"
        "\n    seed    -- random number seed"
        "\n    outFile -- file to contain the InpMtx object"
        "\n               must be *.dinpmtxb or *.dinpmtxf"
        "\n", argv[0]) ;
   return ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return ;
}
n1   = atoi(argv[3]) ;
n2   = atoi(argv[4]) ;
n3   = atoi(argv[5]) ;
seed = atoi(argv[6]) ;
outFileName = argv[7] ;
/*
   --------------
   echo the input
   --------------
*/
fprintf(msgFile, "\n input to %s"
        "\n msglvl  = %d"
        "\n msgFile = %s"
        "\n n1      = %d"
        "\n n2      = %d"
        "\n n3      = %d"
        "\n seed    = %d"
        "\n outFile = %s"
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, seed, outFileName) ;
fflush(msgFile) ;
/*
   ------------------------
   create the EGraph object
   ------------------------
*/
MARKTIME(t1) ;
if ( n3 == 1 ) {
   egraph = EGraph_make9P(n1, n2, 1) ;
   entries = DVinit(4, 0.0) ;
} else {
   egraph = EGraph_make27P(n1, n2, n3, 1) ;
   entries = DVinit(8, 0.0) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : create egraph ", t2 - t1) ;
if ( msglvl > 2 ) {
   EGraph_writeForHumanEye(egraph, msgFile) ;
   fflush(msgFile) ;
}
nvtx  = egraph->nvtx  ;
nelem = egraph->nelem ;
/*
   -------------------------------
   create the random number object
   -------------------------------
*/
drand = Drand_new() ;
Drand_init(drand) ;
Drand_setUniform(drand, -1.0, 1.0) ;
Drand_setSeed(drand, seed) ;
/*
   -------------------------
   create the InpMtx object
   -------------------------
*/
MARKTIME(t1) ;
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_ROWS, SPOOLES_REAL, 0, 0) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : initialize the InpMtx object", 
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n InpMtx after initialization") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
for ( ielem = 0, jrow = 0 ; ielem < nelem ; ielem++ ) {
   IVL_listAndSize(egraph->adjIVL, ielem, &size, &indices) ;
   if ( n3 == 1 ) {
      for ( irow = 0 ; irow < 4 ; irow++, jrow++ ) {
         Drand_fillDvector(drand, size, entries) ;
         InpMtx_inputRealRow(inpmtx, jrow, size, indices, entries) ;
      }
   } else {
      for ( irow = 0 ; irow < 8 ; irow++, jrow++ ) {
         Drand_fillDvector(drand, size, entries) ;
         InpMtx_inputRealRow(inpmtx, jrow, size, indices, entries) ;
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n InpMtx object for a %d x %d x %d grid",
           n1, n2, n3) ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   optionally write out the InpMtx object to a file
   -------------------------------------------------
*/
if ( strcmp(outFileName, "none") != 0 ) {
   InpMtx_writeToFile(inpmtx, outFileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
EGraph_free(egraph) ;
InpMtx_free(inpmtx) ;
Drand_free(drand) ;
DVfree(entries) ;


fprintf(msgFile, "\n") ;

return ; }

/*--------------------------------------------------------------------*/
