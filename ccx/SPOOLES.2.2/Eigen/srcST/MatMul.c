/*  MatMul.c  */

#include "../Bridge.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_MatMul = 0 ;
static double time_MatMul = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose --- to compute a matrix-vector multiply y[] = C * x[]
     where C is the identity, A or B (depending on *pprbtype).

   *pprbtype -- problem type
      *pprbtype = 1 --> vibration problem, matrix is A
      *pprbtype = 2 --> buckling problem, matrix is B
      *pprbtype = 3 --> matrix is identity, y[] = x[]
   *pnrows -- # of rows in x[]
   *pncols -- # of columns in x[]

   created -- 98aug11, cca & jcp
   -------------------------------------------------------------
*/
void 
MatMul ( 
   int      *pnrows, 
   int      *pncols, 
   double   x[], 
   double   y[],
   int      *pprbtype,
   void     *data
) {
int   ncols, nent, nrows ;
#if MYDEBUG > 0
double   t1, t2 ;
MARKTIME(t1) ;
count_MatMul++ ;
fprintf(stdout, "\n (%d) MatMul()", count_MatMul) ;
fflush(stdout) ;
#endif

nrows = *pnrows ;
ncols = *pncols ;
nent  = nrows*ncols ;
if ( *pprbtype == 3 ) {
/*
    --------------------------
    ... matrix is the identity
    --------------------------
*/
   DVcopy(nent, y, x) ;
   return;
} else {
   Bridge     *bridge = (Bridge *) data ; 
   DenseMtx   *X, *Y ;
   double     alpha[2] = {1.0, 0.0} ;
/*
   ---------------------------------
   setup x and y as DenseMtx objects
   ---------------------------------
*/
   X = bridge->X ;
   DenseMtx_init(X, SPOOLES_REAL, 0, 0, nrows, ncols, 1, nrows);
   DVcopy (nent, DenseMtx_entries(X), x);
   if ( bridge->msglvl > 2 ) {
      fprintf(bridge->msgFile, "\n inside MatMul, X") ;
      DenseMtx_writeForHumanEye(X, bridge->msgFile) ;
      fflush(bridge->msgFile) ;
   }
   Y = bridge->Y ;
   DenseMtx_init(Y, SPOOLES_REAL, 0, 0, nrows, ncols, 1, nrows);
   DenseMtx_zero(Y);
   if ( *pprbtype == 1 ) {
/*
       ---------------------------------------
       ... vibration case matrix is 'm' or 'b'
       ---------------------------------------
*/
       InpMtx_sym_mmm(bridge->B, Y, alpha, X);
   } else {
/*
       --------------------------------------
       ... buckling case matrix is 'k' or 'a'
       --------------------------------------
*/
       InpMtx_sym_mmm(bridge->A, Y, alpha, X);
   }
   if ( bridge->msglvl > 2 ) {
      fprintf(bridge->msgFile, "\n inside MatMul, Y") ;
      DenseMtx_writeForHumanEye(Y, bridge->msgFile) ;
      fflush(bridge->msgFile) ;
   }
/*
   --------------------------------
   copy solution into output vector
   --------------------------------
*/
   DVcopy (nent, y, DenseMtx_entries(Y) );
}
#if MYDEBUG > 0
MARKTIME(t2) ;
time_MatMul += t2 - t1 ;
fprintf(stdout, ", %8.3f seconds, %8.3f total time", 
        t2 - t1, time_MatMul) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
