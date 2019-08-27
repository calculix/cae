/*  test_assmbChv.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------
   test the Chv_assembleChv() method.

   created -- 98apr18, cca
   ------------------------------------
*/
{
Chv     *chvI, *chvJ ;
double   imag, real, t1, t2 ;
double   *entriesI, *entriesJ ;
Drand    *drand ;
FILE     *msgFile ;
int      ierr, ii, irow, jcol,
         lastcol, msglvl, ncolI, ncolJ, nDI, nDJ, nentI, nentJ, 
         nrowI, nrowJ, nUI, nUJ, seed, symflag, type ;
int      *colindI, *colindJ, *rowindI, *rowindJ, *temp ;

if ( argc != 10 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile nDJ nUJ nDI nUI type symflag seed "
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    nDJ     -- # of rows and columns in the (1,1) block"
"\n    nUJ     -- # of columns in the (1,2) block"
"\n    nDI     -- # of rows and columns in the (1,1) block"
"\n    nUI     -- # of columns in the (1,2) block"
"\n    type    -- entries type"
"\n       1 --> real"
"\n       2 --> complex"
"\n    symflag -- symmetry flag"
"\n       0 --> symmetric"
"\n       1 --> hermitian"
"\n       2 --> nonsymmetric"
"\n    seed    -- random number seed"
"\n", argv[0]) ;
   return(0) ;
}
if ( (msglvl = atoi(argv[1])) < 0 ) {
   fprintf(stderr, "\n message level must be positive\n") ;
   exit(-1) ;
}
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n unable to open file %s\n", argv[2]) ;
   return(-1) ;
}
nDJ     = atoi(argv[3]) ;
nUJ     = atoi(argv[4]) ;
nDI     = atoi(argv[5]) ;
nUI     = atoi(argv[6]) ;
type    = atoi(argv[7]) ;
symflag = atoi(argv[8]) ;
seed    = atoi(argv[9]) ;
if (  nDJ <= 0 || nUJ < 0 
   || nDI <= 0 || nUI < 0 
   || nDI >= nDJ || (nDI + nUI) >= (nDJ + nUJ)
   || nUI >= (nDJ + nUJ - nDI)
   || (  symflag != SPOOLES_SYMMETRIC
      && symflag != SPOOLES_HERMITIAN
      && symflag != SPOOLES_NONSYMMETRIC) ) {
   fprintf(stderr, "\n invalid input"
      "\n nDJ = %d, nUJ = %d, nDI = %d, nUI = %d, symflag = %d\n",
           nDJ, nUJ, nDI, nUI, symflag) ;
   exit(-1) ;
}
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
Drand_init(drand) ;
Drand_setSeed(drand, seed) ;
Drand_setUniform(drand, -1.0, 1.0) ;
/*
   ----------------------------
   initialize the ChvJ object
   ----------------------------
*/
MARKTIME(t1) ;
chvJ = Chv_new() ;
Chv_init(chvJ, 0, nDJ, nUJ, nUJ, type, symflag) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize chv object",
        t2 - t1) ;
fflush(msgFile) ;
Chv_columnIndices(chvJ, &ncolJ, &colindJ) ;
temp = IVinit(2*(nDJ+nUJ), -1) ;
IVramp(2*(nDJ+nUJ), temp, 0, 1) ;
IVshuffle(2*(nDJ+nUJ), temp, ++seed) ;
IVcopy(ncolJ, colindJ, temp) ;
IVfree(temp) ;
IVqsortUp(ncolJ, colindJ) ;
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
   Chv_rowIndices(chvJ, &nrowJ, &rowindJ) ;
   IVcopy(nrowJ, rowindJ, colindJ) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n %% column indices") ;
   IVfprintf(msgFile, ncolJ, colindJ) ;
}
lastcol = colindJ[ncolJ-1] ;
nentJ = Chv_nent(chvJ) ;
entriesJ = Chv_entries(chvJ) ;
if ( CHV_IS_REAL(chvJ) ) {
   Drand_fillDvector(drand, nentJ, entriesJ) ;
} else if ( CHV_IS_COMPLEX(chvJ) ) {
   Drand_fillDvector(drand, 2*nentJ, entriesJ) ;
}
if ( CHV_IS_HERMITIAN(chvJ) ) {
/*
   ---------------------------------------------------------
   hermitian example, set imaginary part of diagonal to zero
   ---------------------------------------------------------
*/
   for ( irow = 0 ; irow < nDJ ; irow++ ) {
      Chv_complexEntry(chvJ, irow, irow, &real, &imag) ;
      Chv_setComplexEntry(chvJ, irow, irow, real, 0.0) ;
   }
}
/*
   ---------------------------
   initialize the ChvI object
   ---------------------------
*/
chvI = Chv_new() ;
Chv_init(chvI, 0, nDI, nUI, nUI, type, symflag) ;
Chv_columnIndices(chvI, &ncolI, &colindI) ;
temp = IVinit(ncolJ, -1) ;
IVramp(ncolJ, temp, 0, 1) ;
while ( 1 ) {
   IVshuffle(ncolJ, temp, ++seed) ;
   IVqsortUp(ncolI, temp) ;
   if ( temp[0] < nDJ ) {
      break ;
   }
}
for ( ii = 0 ; ii < ncolI ; ii++ ) {
   colindI[ii] = colindJ[temp[ii]] ;
}
IVfree(temp) ;
if ( CHV_IS_NONSYMMETRIC(chvI) ) {
   Chv_rowIndices(chvI, &nrowI, &rowindI) ;
   IVcopy(nrowI, rowindI, colindI) ;
}
nentI = Chv_nent(chvI) ;
entriesI = Chv_entries(chvI) ;
if ( CHV_IS_REAL(chvI) ) {
   Drand_fillDvector(drand, nentI, entriesI) ;
} else if ( CHV_IS_COMPLEX(chvI) ) {
   Drand_fillDvector(drand, 2*nentI, entriesI) ;
}
if ( CHV_IS_HERMITIAN(chvI) ) {
/*
   ---------------------------------------------------------
   hermitian example, set imaginary part of diagonal to zero
   ---------------------------------------------------------
*/
   for ( irow = 0 ; irow < nDI ; irow++ ) {
      Chv_complexEntry(chvI, irow, irow, &real, &imag) ;
      Chv_setComplexEntry(chvI, irow, irow, real, 0.0) ;
   }
}
/*
   --------------------------------------------------
   write out the two chevron objects to a matlab file
   --------------------------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n a = zeros(%d,%d) ;", lastcol+1, lastcol+1) ;
   Chv_writeForMatlab(chvJ, "a", msgFile) ;
   fprintf(msgFile, "\n b = zeros(%d,%d) ;", lastcol+1, lastcol+1) ;
   Chv_writeForMatlab(chvI, "b", msgFile) ;
}
/*
   ---------------------------------------------
   assemble the chvI object into the chvJ object
   ---------------------------------------------
*/
Chv_assembleChv(chvJ, chvI) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n %% after assembly") ;
   fprintf(msgFile, "\n c = zeros(%d,%d) ;", lastcol+1, lastcol+1) ;
   Chv_writeForMatlab(chvJ, "c", msgFile) ;
}
/*
   -----------------
   compute the error
   -----------------
*/
fprintf(msgFile, "\n max(max(abs(c - (b + a))))") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Chv_free(chvJ) ;
Chv_free(chvI) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
