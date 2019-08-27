/*  axpy.c  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] + alpha[2] * x2[]
   y1[] = y1[] + alpha[3] * x0[] + alpha[4] * x1[] + alpha[5] * x2[]
   y2[] = y2[] + alpha[6] * x0[] + alpha[7] * x1[] + alpha[8] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
DVaxpy33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
register double   a00, a01, a02, a10, a11, a12, a20, a21, a22 ;
register double   x0i, x1i, x2i ;
int               ii ;

a00 = alpha[0] ; a01 = alpha[1] ; a02 = alpha[2] ;
a10 = alpha[3] ; a11 = alpha[4] ; a12 = alpha[5] ;
a20 = alpha[6] ; a21 = alpha[7] ; a22 = alpha[8] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ; x1i = x1[ii] ; x2i = x2[ii] ;
   y0[ii] += a00*x0i + a01*x1i + a02*x2i ;
   y1[ii] += a10*x0i + a11*x1i + a12*x2i ;
   y2[ii] += a20*x0i + a21*x1i + a22*x2i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] 
   y1[] = y1[] + alpha[2] * x0[] + alpha[3] * x1[] 
   y2[] = y2[] + alpha[4] * x0[] + alpha[5] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------
*/
void
DVaxpy32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[]
) {
register double   a00, a01, a10, a11, a20, a21 ;
register double   x0i, x1i ;
int               ii ;

a00 = alpha[0] ; a01 = alpha[1] ;
a10 = alpha[2] ; a11 = alpha[3] ;
a20 = alpha[4] ; a21 = alpha[5] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ; x1i = x1[ii] ;
   y0[ii] += a00*x0i + a01*x1i ;
   y1[ii] += a10*x0i + a11*x1i ;
   y2[ii] += a20*x0i + a21*x1i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   y0[] = y0[] + alpha[0] * x0[] 
   y1[] = y1[] + alpha[1] * x0[] 
   y2[] = y2[] + alpha[2] * x0[] 

   created -- 98dec10, cca
   -----------------------------
*/
void
DVaxpy31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[]
) {
register double   a00, a10, a20 ;
register double   x0i ;
int               ii ;

a00 = alpha[0] ;
a10 = alpha[1] ;
a20 = alpha[2] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ;
   y0[ii] += a00*x0i ;
   y1[ii] += a10*x0i ;
   y2[ii] += a20*x0i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] + alpha[2] * x2[]
   y1[] = y1[] + alpha[3] * x0[] + alpha[4] * x1[] + alpha[5] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
DVaxpy23 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
register double   a00, a01, a02, a10, a11, a12 ;
register double   x0i, x1i, x2i ;
int               ii ;

a00 = alpha[0] ; a01 = alpha[1] ; a02 = alpha[2] ;
a10 = alpha[3] ; a11 = alpha[4] ; a12 = alpha[5] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ; x1i = x1[ii] ; x2i = x2[ii] ;
   y0[ii] += a00*x0i + a01*x1i + a02*x2i ;
   y1[ii] += a10*x0i + a11*x1i + a12*x2i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] 
   y1[] = y1[] + alpha[2] * x0[] + alpha[3] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------
*/
void
DVaxpy22 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[]
) {
register double   a00, a01, a10, a11 ;
register double   x0i, x1i ;
int               ii ;

a00 = alpha[0] ; a01 = alpha[1] ;
a10 = alpha[2] ; a11 = alpha[3] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ; x1i = x1[ii] ;
   y0[ii] += a00*x0i + a01*x1i ;
   y1[ii] += a10*x0i + a11*x1i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   y0[] = y0[] + alpha[0] * x0[] 
   y1[] = y1[] + alpha[1] * x0[] 

   created -- 98dec10, cca
   -----------------------------
*/
void
DVaxpy21 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[]
) {
register double   a00, a10 ;
register double   x0i ;
int               ii ;

a00 = alpha[0] ;
a10 = alpha[1] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ;
   y0[ii] += a00*x0i ;
   y1[ii] += a10*x0i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] + alpha[2] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
DVaxpy13 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
register double   a00, a01, a02 ;
register double   x0i, x1i, x2i ;
int               ii ;

a00 = alpha[0] ; a01 = alpha[1] ; a02 = alpha[2] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ; x1i = x1[ii] ; x2i = x2[ii] ;
   y0[ii] += a00*x0i + a01*x1i + a02*x2i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------
*/
void
DVaxpy12 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[]
) {
register double   a00, a01 ;
register double   x0i, x1i ;
int               ii ;

a00 = alpha[0] ; a01 = alpha[1] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ; x1i = x1[ii] ;
   y0[ii] += a00*x0i + a01*x1i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   y0[] = y0[] + alpha[0] * x0[] 

   created -- 98dec10, cca
   -----------------------------
*/
void
DVaxpy11 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[]
) {
register double   a00 ;
register double   x0i ;
int               ii ;

a00 = alpha[0] ;
for ( ii = 0 ; ii < n ; ii++ ) {
   x0i = x0[ii] ;
   y0[ii] += a00*x0i ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] 
               + alpha[2:3] * x1[] + alpha[4:5] * x2[]
   y1[] = y1[] + alpha[6:7] * x0[] 
               + alpha[8:9] * x1[] + alpha[10:11] * x2[]
   y2[] = y2[] + alpha[12:13] * x0[] 
               + alpha[14:15] * x1[] + alpha[16:17] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
ZVaxpy33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
register double   ar00, ar01, ar02, ar10, ar11, ar12, ar20, ar21, ar22 ;
register double   ai00, ai01, ai02, ai10, ai11, ai12, ai20, ai21, ai22 ;
register double   xr0i, xr1i, xr2i ;
register double   xi0i, xi1i, xi2i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar01 = alpha[ 2] ; ai01 = alpha[ 3] ; 
ar02 = alpha[ 4] ; ai02 = alpha[ 5] ; 
ar10 = alpha[ 6] ; ai10 = alpha[ 7] ; 
ar11 = alpha[ 8] ; ai11 = alpha[ 9] ; 
ar12 = alpha[10] ; ai12 = alpha[11] ; 
ar20 = alpha[12] ; ai20 = alpha[13] ; 
ar21 = alpha[14] ; ai21 = alpha[15] ; 
ar22 = alpha[16] ; ai22 = alpha[17] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 
   xr1i = x1[2*ii] ; xi1i = x1[2*ii+1] ; 
   xr2i = x2[2*ii] ; xi2i = x2[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) 
              +  (ar01*xr1i - ai01*xi1i)
              +  (ar02*xr2i - ai02*xi2i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) 
              +  (ar01*xi1i + ai01*xr1i)
              +  (ar02*xi2i + ai02*xr2i) ;
   y1[2*ii]   += (ar10*xr0i - ai10*xi0i) 
              +  (ar11*xr1i - ai11*xi1i)
              +  (ar12*xr2i - ai12*xi2i) ;
   y1[2*ii+1] += (ar10*xi0i + ai10*xr0i) 
              +  (ar11*xi1i + ai11*xr1i)
              +  (ar12*xi2i + ai12*xr2i) ;
   y2[2*ii]   += (ar20*xr0i - ai20*xi0i) 
              +  (ar21*xr1i - ai21*xi1i)
              +  (ar22*xr2i - ai22*xi2i) ;
   y2[2*ii+1] += (ar20*xi0i + ai20*xr0i) 
              +  (ar21*xi1i + ai21*xr1i)
              +  (ar22*xi2i + ai22*xr2i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] + alpha[2:3] * x1[] 
   y1[] = y1[] + alpha[4:5] * x0[] + alpha[6:7] * x1[] 
   y2[] = y2[] + alpha[8:9] * x0[] + alpha[10:11] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------------
*/
void
ZVaxpy32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[]
) {
register double   ar00, ar01, ar10, ar11, ar20, ar21 ;
register double   ai00, ai01, ai10, ai11, ai20, ai21 ;
register double   xr0i, xr1i ;
register double   xi0i, xi1i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar01 = alpha[ 2] ; ai01 = alpha[ 3] ; 
ar10 = alpha[ 4] ; ai10 = alpha[ 5] ; 
ar11 = alpha[ 6] ; ai11 = alpha[ 7] ; 
ar20 = alpha[ 8] ; ai20 = alpha[ 9] ; 
ar21 = alpha[10] ; ai21 = alpha[11] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 
   xr1i = x1[2*ii] ; xi1i = x1[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) 
              +  (ar01*xr1i - ai01*xi1i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) 
              +  (ar01*xi1i + ai01*xr1i) ;
   y1[2*ii]   += (ar10*xr0i - ai10*xi0i) 
              +  (ar11*xr1i - ai11*xi1i) ;
   y1[2*ii+1] += (ar10*xi0i + ai10*xr0i) 
              +  (ar11*xi1i + ai11*xr1i) ;
   y2[2*ii]   += (ar20*xr0i - ai20*xi0i) 
              +  (ar21*xr1i - ai21*xi1i) ;
   y2[2*ii+1] += (ar20*xi0i + ai20*xr0i) 
              +  (ar21*xi1i + ai21*xr1i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] 
   y1[] = y1[] + alpha[2:3] * x0[] 
   y2[] = y2[] + alpha[4:5] * x0[] 

   created -- 98dec10, cca
   -------------------------------
*/
void
ZVaxpy31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[]
) {
register double   ar00, ar10, ar20 ;
register double   ai00, ai10, ai20 ;
register double   xr0i ;
register double   xi0i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar10 = alpha[ 2] ; ai10 = alpha[ 3] ; 
ar20 = alpha[ 4] ; ai20 = alpha[ 5] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) ;
   y1[2*ii]   += (ar10*xr0i - ai10*xi0i) ;
   y1[2*ii+1] += (ar10*xi0i + ai10*xr0i) ;
   y2[2*ii]   += (ar20*xr0i - ai20*xi0i) ;
   y2[2*ii+1] += (ar20*xi0i + ai20*xr0i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] 
               + alpha[2:3] * x1[] + alpha[4:5] * x2[]
   y1[] = y1[] + alpha[6:7] * x0[] 
               + alpha[8:9] * x1[] + alpha[10:11] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
ZVaxpy23 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
register double   ar00, ar01, ar02, ar10, ar11, ar12 ;
register double   ai00, ai01, ai02, ai10, ai11, ai12 ;
register double   xr0i, xr1i, xr2i ;
register double   xi0i, xi1i, xi2i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar01 = alpha[ 2] ; ai01 = alpha[ 3] ; 
ar02 = alpha[ 4] ; ai02 = alpha[ 5] ; 
ar10 = alpha[ 6] ; ai10 = alpha[ 7] ; 
ar11 = alpha[ 8] ; ai11 = alpha[ 9] ; 
ar12 = alpha[10] ; ai12 = alpha[11] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 
   xr1i = x1[2*ii] ; xi1i = x1[2*ii+1] ; 
   xr2i = x2[2*ii] ; xi2i = x2[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) 
              +  (ar01*xr1i - ai01*xi1i)
              +  (ar02*xr2i - ai02*xi2i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) 
              +  (ar01*xi1i + ai01*xr1i)
              +  (ar02*xi2i + ai02*xr2i) ;
   y1[2*ii]   += (ar10*xr0i - ai10*xi0i) 
              +  (ar11*xr1i - ai11*xi1i)
              +  (ar12*xr2i - ai12*xi2i) ;
   y1[2*ii+1] += (ar10*xi0i + ai10*xr0i) 
              +  (ar11*xi1i + ai11*xr1i)
              +  (ar12*xi2i + ai12*xr2i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] + alpha[2:3] * x1[] 
   y1[] = y1[] + alpha[4:5] * x0[] + alpha[6:7] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------------
*/
void
ZVaxpy22 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[]
) {
register double   ar00, ar01, ar10, ar11 ;
register double   ai00, ai01, ai10, ai11 ;
register double   xr0i, xr1i ;
register double   xi0i, xi1i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar01 = alpha[ 2] ; ai01 = alpha[ 3] ; 
ar10 = alpha[ 4] ; ai10 = alpha[ 5] ; 
ar11 = alpha[ 6] ; ai11 = alpha[ 7] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 
   xr1i = x1[2*ii] ; xi1i = x1[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) 
              +  (ar01*xr1i - ai01*xi1i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) 
              +  (ar01*xi1i + ai01*xr1i) ;
   y1[2*ii]   += (ar10*xr0i - ai10*xi0i) 
              +  (ar11*xr1i - ai11*xi1i) ;
   y1[2*ii+1] += (ar10*xi0i + ai10*xr0i) 
              +  (ar11*xi1i + ai11*xr1i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] 
   y1[] = y1[] + alpha[2:3] * x0[] 

   created -- 98dec10, cca
   -------------------------------
*/
void
ZVaxpy21 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[]
) {
register double   ar00, ar10 ;
register double   ai00, ai10 ;
register double   xr0i ;
register double   xi0i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar10 = alpha[ 2] ; ai10 = alpha[ 3] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) ;
   y1[2*ii]   += (ar10*xr0i - ai10*xi0i) ;
   y1[2*ii+1] += (ar10*xi0i + ai10*xr0i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] 
               + alpha[2:3] * x1[] + alpha[4:5] * x2[]

   created -- 98dec10, cca
   ---------------------------------------------------
*/
void
ZVaxpy13 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
register double   ar00, ar01, ar02 ;
register double   ai00, ai01, ai02 ;
register double   xr0i, xr1i, xr2i ;
register double   xi0i, xi1i, xi2i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar01 = alpha[ 2] ; ai01 = alpha[ 3] ; 
ar02 = alpha[ 4] ; ai02 = alpha[ 5] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 
   xr1i = x1[2*ii] ; xi1i = x1[2*ii+1] ; 
   xr2i = x2[2*ii] ; xi2i = x2[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) 
              +  (ar01*xr1i - ai01*xi1i)
              +  (ar02*xr2i - ai02*xi2i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) 
              +  (ar01*xi1i + ai01*xr1i)
              +  (ar02*xi2i + ai02*xr2i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] + alpha[2:3] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------------
*/
void
ZVaxpy12 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[]
) {
register double   ar00, ar01 ;
register double   ai00, ai01 ;
register double   xr0i, xr1i ;
register double   xi0i, xi1i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
ar01 = alpha[ 2] ; ai01 = alpha[ 3] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 
   xr1i = x1[2*ii] ; xi1i = x1[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) 
              +  (ar01*xr1i - ai01*xi1i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) 
              +  (ar01*xi1i + ai01*xr1i) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] 

   created -- 98dec10, cca
   -------------------------------
*/
void
ZVaxpy11 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[]
) {
register double   ar00 ;
register double   ai00 ;
register double   xr0i ;
register double   xi0i ;
int               ii ;

ar00 = alpha[ 0] ; ai00 = alpha[ 1] ; 
for ( ii = 0 ; ii < n ; ii++ ) {
   xr0i = x0[2*ii] ; xi0i = x0[2*ii+1] ; 

   y0[2*ii]   += (ar00*xr0i - ai00*xi0i) ;
   y0[2*ii+1] += (ar00*xi0i + ai00*xr0i) ;
}
return ; }

/*--------------------------------------------------------------------*/
