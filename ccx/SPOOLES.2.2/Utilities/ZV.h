/*  ZV.h  */

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to return the absolute value of (areal,aimag)
 
   created -- 98jan24, cca
   --------------------------------------------------------
*/
double
Zabs (
   double  real,
   double  imag
) ;
/*
   ---------------------------------------------------
   purpose -- given (areal,aimag), 
              compute (breal,bimag) = 1/(areal,aimag), 
 
   put breal into *pbreal
   put bimag into *pbimag
 
   created -- 98jan23, cca
   ---------------------------------------------------
*/
int
Zrecip ( 
   double   areal,
   double   aimag,
   double   *preal,
   double   *pimag
) ;
/*
   ---------------------------------------------------------------------
   given [ (areal,aimag) (breal,bimag) ]
         [ (creal,cimag) (dreal,dimag) ]
   compute [ (ereal,eimag) (freal,fimag) ]
           [ (greal,gimag) (hreal,himag) ]
   where
   I = [ (areal,aimag) (breal,bimag) ] * [ (ereal,eimag) (freal,fimag) ]
       [ (creal,cimag) (dreal,dimag) ]   [ (greal,gimag) (hreal,himag) ]
 
   note, any of {pereal, peimag, pfreal, pfimag, pgreal, pgimag,
   phreal, phimag} can be NULL. if not NULL, their value is filled.
   this is useful when the input matrix is symmetric or hermitian.
 
   created -- 98jan23, cca
   ---------------------------------------------------------------------
*/
int
Zrecip2 (
   double areal,
   double aimag,
   double breal,
   double bimag,
   double creal,
   double cimag,
   double dreal,
   double dimag,
   double *pereal,
   double *peimag,
   double *pfreal,
   double *pfimag,
   double *pgreal,
   double *pgimag,
   double *phreal,
   double *phimag
) ;
/*
   ------------------------------------------------
   purpose -- to allocate, initialize and
              return a complex vector x[]
 
   n -- length of the complex vector (2*n double's)
   
   x[ii] = (real, imag) for 0 <= ii < n
 
   created -- 98jan23, cca
   ------------------------------------------------
*/
double *
ZVinit (
   int      n,
   double   real,
   double   imag
) ;
/*
   -------------------------------------------
   purpose -- to perform a complex dot product
 
   (*prdot,*pidot) = y^T x
 
   where y and x are complex
 
   created -- 98apr15, cca
   -------------------------------------------
*/
void
ZVdotU ( 
   int      size, 
   double   y[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) ;
/*
   ------------------------------------------------------
   purpose -- to perform a conjugated complex dot product
 
   (*prdot,*pidot) = conjugate(y^T) x
 
   where y and x are complex
 
   created -- 98apr15, cca
   ------------------------------------------------------
*/
void
ZVdotC ( 
   int      size, 
   double   y[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) ;
/*
   ---------------------------------------------------
   purpose -- to perform a indexed complex dot product
 
   (*prdot,*pidot) = sum_k y[index[k]]*x[k]
 
   where y and x are complex
 
   created -- 98apr15, cca
   ---------------------------------------------------
*/
void
ZVdotiU ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) ;
/*
   -------------------------------------------------------------
   purpose -- to perform a indexed conjugate complex dot product
 
   (*prdot,*pidot) = sum_k conjugate(y[index[k]])*x[k]
 
   where y and x are complex
 
   created -- 98apr15, cca
   -------------------------------------------------------------
*/
void
ZVdotiC ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) ;
/*
   ------------------------------------
   purpose -- to perform a complex axpy
 
   y := y + (areal, aimag) * x
 
   where y and x are complex
 
   created -- 98jan22, cca
   ------------------------------------
*/
void
ZVaxpy (
   int      size,
   double   y[],
   double   areal,
   double   aimag,
   double   x[]
) ;
/*
   -----------------------------------------------------
   purpose -- to perform a complex axpy with two vectors
 
   z := z + (areal, aimag)*x + (breal, bimag)*y
 
   where y and x are complex
 
   created -- 98jan23, cca
   ----------------------------------------------------
*/
void
ZVaxpy2 ( 
   int      size, 
   double   z[], 
   double   areal, 
   double   aimag, 
   double   x[],
   double   breal, 
   double   bimag, 
   double   y[]
) ;
/*
   ------------------------------------------------------------
   purpose -- to scale a double complex vector by (xreal,ximag)
 
   y := y * (areal, aimag)
 
   created -- 98jan22, cca
   ------------------------------------------------------------
*/
void
ZVscale (
   int      size,
   double   y[],
   double   areal,
   double   aimag
) ;
/*
   ------------------------------------------------------------
   purpose -- to scale a double complex vector by a 2x2 matrix
 
   [ y0 ] := [ d00 d01 ] [ y0 ]
   [ y1 ]    [ d10 d11 ] [ y1 ]
 
   created -- 98jan23, cca
   ------------------------------------------------------------
*/
void
ZVscale2 ( 
   int      size, 
   double   y0[], 
   double   y1[], 
   double   d00r,
   double   d00i,
   double   d01r,
   double   d01i,
   double   d10r,
   double   d10i,
   double   d11r,
   double   d11i
) ;
/*
   -----------------------------------------------
   purpose --- to print a complex vector to a file
 
   created -- 98jan23, cca
   -----------------------------------------------
*/
void
ZVfprintf (
   FILE     *fp,
   int      size,
   double   y[]
) ;
/*
   ----------------------------------
   return the minimum absolute value
   of the entries in a complex vector
 
   created -- 98jan23, cca
   ----------------------------------
*/
double
ZVminabs (
   int      size,
   double   x[]
) ;
/*
   ----------------------------------
   return the maximum absolute value
   of the entries in a complex vector
 
   created -- 98jan23, cca
   ----------------------------------
*/
double
ZVmaxabs (
   int      size,
   double   x[]
) ;
/*
   ----------------------------------
   copy a complex vector into another
   y[] := x[]
 
   created -- 98jan23, cca
   ----------------------------------
*/
void
ZVcopy (
   int      size,
   double   y[],
   double   x[]
) ;
/*
   --------------------------------------
   subtract a complex vector from another
   y[] := x[]
 
   created -- 98may25, cca
   --------------------------------------
*/
void
ZVsub (
   int      size,
   double   y[],
   double   x[]
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to gather y[*] = x[index[*]]
 
   created -- 98apr15, cca
   ---------------------------------------
*/
void
ZVgather ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to scatter y[index[*]] = x[*]
 
   created -- 98apr15, cca
   ----------------------------------------
*/
void
ZVscatter (
   int      size,
   double   y[],
   int      index[],
   double   x[]
) ;
/*
   -----------------------------
   purpose -- to zero the vector
 
   y := 0
 
   where y is complex
 
   created -- 98apr25, cca
   -----------------------------
*/
void
ZVzero ( 
   int      size, 
   double   y[] 
) ;
/*--------------------------------------------------------------------*/
