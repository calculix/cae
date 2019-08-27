/*  mkLaplacianMtx.c  */

#include "../InpMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   create a laplacian matrix for a regular grid
   for a n1 x n2 grid we use the following stencil
      [ -1 -1 -1 ]
      [ -1  8 -1 ]
      [ -1 -1 -1 ]
   for a n1 x n2 x n3 grid we use the following stencil
      [ -1 -1 -1 ] [ -1 -1 -1 ] [ -1 -1 -1 ]
      [ -1 -1 -1 ] [ -1 27 -1 ] [ -1 -1 -1 ]
      [ -1 -1 -1 ] [ -1 -1 -1 ] [ -1 -1 -1 ]

   created -- 98nov06, cca
   ---------------------------------------------------------------
*/
void
main ( int argc, char *argv[] ) 
{
char       *outFileName ;
double     *dvec ;
double     t1, t2 ;
FILE       *msgFile ;
InpMtx     *inpmtx ;
int        ient, ii, msglvl, nent, nvtx, n1, n2, n3, 
           size, stencil, v, vsize, w ;
int        *ivec1, *ivec2, *vadj ;
IVL        *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( argc != 7 ) {
   fprintf(stdout, 
        "\n\n usage : %s msglvl msgFile n1 n2 n3 outFile "
        "\n    msglvl  -- message level"
        "\n    msgFile -- message file"
        "\n    n1      -- number of grid points in the first direction"
        "\n    n2      -- number of grid points in the second direction"
        "\n    n3      -- number of grid points in the third direction"
        "\n    outFile -- file to contain the InpMtx object"
        "\n               must be *.inpmtxb or *.inpmtxf"
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
outFileName = argv[6] ;
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
        "\n outFile = %s"
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, outFileName) ;
fflush(msgFile) ;
nvtx = n1 * n2 * n3 ;
/*
   ----------------------------------------
   create the grid graph's adjacency object
   ----------------------------------------
*/
if ( n1 == 1 ) {
   adjIVL = IVL_make9P(n2, n3, 1) ;
   stencil = 9 ;
} else if ( n2 == 1 ) {
   adjIVL = IVL_make9P(n1, n3, 1) ;
   stencil = 9 ;
} else if ( n3 == 1 ) {
   adjIVL = IVL_make9P(n1, n2, 1) ;
   stencil = 9 ;
} else {
   adjIVL = IVL_make27P(n1, n2, n3, 1) ;
   stencil = 27 ;
}
nent = adjIVL->tsize ;
/*
   -------------------------
   create the InpMtx object
   -------------------------
*/
MARKTIME(t1) ;
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_ROWS, SPOOLES_REAL, nent, 0) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : initialize the InpMtx object", 
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n InpMtx after initialization") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
dvec  = InpMtx_dvec(inpmtx) ;
if ( stencil == 9 ) {
   for ( v = ient = 0 ; v < nvtx ; v++ ) {
      IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) == v ) {
            ivec1[ient] =  v  ;
            ivec2[ient] =  w  ;
            dvec[ient]  = 8.0 ;
            ient++ ;
         } else if ( w > v ) {
            ivec1[ient] =   v  ;
            ivec2[ient] =   w  ;
            dvec[ient]  = -1.0 ;
            ient++ ;
         }
      }
   }
} else {
   for ( v = ient = 0 ; v < nvtx ; v++ ) {
      IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) == v ) {
            ivec1[ient] =   v  ;
            ivec2[ient] =   w  ;
            dvec[ient]  = 27.0 ;
            ient++ ;
         } else if ( w > v ) {
            ivec1[ient] =   v  ;
            ivec2[ient] =   w  ;
            dvec[ient]  = -1.0 ;
            ient++ ;
         }
      }
   }
}
inpmtx->nent = ient ;
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
IVL_free(adjIVL) ;
InpMtx_free(inpmtx) ;

fprintf(msgFile, "\n") ;

return ; }

/*--------------------------------------------------------------------*/
