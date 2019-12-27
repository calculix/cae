/*  testFactor.c  */

#include "../ILUMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------------------
   generate a random matrix and test the drop tolerance factorization.
   the output is a matlab file to test correctness.

   created -- 98oct08, cca
 ---------------------------------------------------------------------
*/
{
double   droptol, nops, t1, t2 ;
Drand    *drand ;
int      msglvl, neqns, nitem, rc, seed, symflag, type ;
InpMtx   *A ;
ILUMtx   *mtxLDU ;
FILE     *matlabFile, *msgFile ;

if ( argc != 10 ) {
   fprintf(stdout, 
      "\n\n %% usage : %s msglvl msgFile type symflag "
      "\n            neqns nitem seed droptol matlabFile"
      "\n %%    msglvl  -- message level"
      "\n %%    msgFile -- message file"
      "\n %%    type    -- type of matrix entries"
      "\n %%       1 -- real"
      "\n %%       2 -- complex"
      "\n %%    symflag -- symmetry flag"
      "\n %%       0 -- symmetric"
      "\n %%       1 -- hermitian"
      "\n %%       2 -- nonsymmetric"
      "\n %%    neqns   -- number of equations"
      "\n %%    nitem   -- number of items"
      "\n %%    seed    -- random number seed"
      "\n %%    droptol -- drop tolerance"
      "\n %%    matlabFile -- matlab file name"
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
type    = atoi(argv[3]) ;
symflag = atoi(argv[4]) ;
neqns   = atoi(argv[5]) ;
nitem   = atoi(argv[6]) ;
seed    = atoi(argv[7]) ;
droptol = atof(argv[8]) ;
if ( strcmp(argv[9], "stdout") == 0 ) {
   matlabFile = stdout ;
} else if ( strcmp(argv[2], argv[9]) == 0 ) {
   matlabFile = msgFile ;
} else if ( strcmp("none", argv[9]) != 0 ) {
   if ( (matlabFile = fopen(argv[9], "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], argv[9]) ;
      return(-1) ;
   }
} else {
   matlabFile = NULL ;
}
fprintf(msgFile, 
        "\n %% %s "
        "\n %% msglvl     -- %d" 
        "\n %% msgFile    -- %s" 
        "\n %% type       -- %d" 
        "\n %% symflag    -- %d" 
        "\n %% neqns      -- %d" 
        "\n %% nitem      -- %d" 
        "\n %% seed       -- %d" 
        "\n %% droptol    -- %e"
        "\n %% matlabFile -- %s"
        "\n",
        argv[0], msglvl, argv[2], type, symflag, 
        neqns, nitem, seed, droptol, argv[9]) ;
fflush(msgFile) ;
if ( type != SPOOLES_REAL && type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n invalid value %d for type\n", type) ;
   exit(-1) ;
}
switch ( symflag ) {
case SPOOLES_SYMMETRIC    :
case SPOOLES_HERMITIAN    :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n invalid value %d for symflag\n", symflag) ;
   exit(-1) ;
   break ;
}
if ( symflag == SPOOLES_HERMITIAN && type == SPOOLES_REAL ) {
   fprintf(stderr, "\n symflag is hermitian and type is real") ;
   exit(-1) ;
}
if ( neqns <= 0 || nitem <= 0 ) {
   fprintf(stderr, "\n invalid value: neqns = %d, nitem = %d",
           neqns, nitem) ;
   exit(-1) ;
}
/*
   ----------------------------
   initialize the matrix object
   ----------------------------
*/
A = InpMtx_new() ;
InpMtx_randomMatrix(A, type, INPMTX_BY_CHEVRONS, INPMTX_BY_VECTORS,
                    neqns, neqns, symflag, 1, nitem, seed) ;
/*
   -------------------------------------------
   write the assembled matrix to a matlab file
   -------------------------------------------
*/
if ( matlabFile != NULL ) {
   InpMtx_writeForMatlab(A, "A", matlabFile) ;
   if ( symflag == SPOOLES_SYMMETRIC ) {
      fprintf(matlabFile,
              "\n   for k = 1:%d"
              "\n      for j = k+1:%d"
              "\n         A(j,k) = A(k,j) ;"
              "\n      end"
              "\n   end", neqns, neqns) ;
   } else if ( symflag == SPOOLES_HERMITIAN ) {
      fprintf(matlabFile,
              "\n   for k = 1:%d"
              "\n      for j = k+1:%d"
              "\n         A(j,k) = ctranspose(A(k,j)) ;"
              "\n      end"
              "\n   end", neqns, neqns) ;
   }
   InpMtx_changeCoordType(A, INPMTX_BY_CHEVRONS) ;
   InpMtx_changeStorageMode(A, INPMTX_BY_VECTORS) ;
/*
   -----------------------------------------------
   compute the incomplete factorization via matlab
   -----------------------------------------------
*/
   fprintf(msgFile, "\n\n droptol = %24.16e ;", droptol) ;
   fprintf(msgFile, "\n\n [ L, D, U ] = ilu(A, droptol) ;") ;
}
/*
   ------------------------------------
   compute the incomplete factorization
   ------------------------------------
*/
mtxLDU = ILUMtx_new() ;
ILUMtx_init(mtxLDU, neqns, type, symflag, 
            SPOOLES_BY_COLUMNS, SPOOLES_BY_ROWS) ;
nops = 0.0 ;
MARKTIME(t1) ;
rc = ILUMtx_factor(mtxLDU, A, droptol, &nops, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n %% CPU %8.3f : compute ILU, %.0f operations, %.2f mflops\n",
        t2 - t1, nops, 1.e-6*nops/(t2 - t1)) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error return %d from ILUMtx_factor()", rc) ;
   return(rc) ;
}
if ( matlabFile != NULL ) {
   rc = ILUMtx_writeForMatlab(mtxLDU, "Lhat", "Dhat", "Uhat", 
                              matlabFile) ;
/*
   -----------------
   compute the error
   -----------------
*/
   fprintf(matlabFile, "\n\n errorL = max(max(abs(Lhat - L))) ;") ;
   fprintf(matlabFile, "\n\n errorD = max(max(abs(Dhat - D))) ;") ;
   fprintf(matlabFile, "\n\n errorU = max(max(abs(Uhat - U))) ;") ;
   fprintf(matlabFile, "\n\n error = [ errorL errorD errorU ]") ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
InpMtx_free(A) ;
ILUMtx_free(mtxLDU) ;

fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
