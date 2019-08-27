/*  instance.c  */

#include "../A2.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   return the number of rows in the array

   created -- 98may01, cca
   --------------------------------------
*/
int
A2_nrow ( 
   A2   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_nrow(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->n1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   return the number of columns in the array

   created -- 98may01, cca
   -----------------------------------------
*/
int
A2_ncol ( 
   A2   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_ncol(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->n2) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   return the first increment

   created -- 98may01, cca
   --------------------------
*/
int
A2_inc1 ( 
   A2  *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_inc1(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->inc1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   return the second increment

   created -- 98may01, cca
   ---------------------------
*/
int
A2_inc2 ( 
   A2  *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_inc2(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->inc2) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   return a pointer to the entries

   created -- 98may01, cca
   -------------------------------
*/
double *
A2_entries ( 
   A2      *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_entries(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->entries) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   return a pointer to the first entry in a row

   created -- 98may01, cca
   --------------------------------------------
*/
double *
A2_row ( 
   A2   *mtx,
   int   irow
) {
double   *row ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_row(%p,%d)"
           "\n bad input\n", mtx, irow) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in A2_row(%p,%d)"
           "\n bad structure, entries is NULL\n", mtx, irow) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_row(%p,%d)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, irow, mtx->n1) ;
   exit(-1) ;
}
if ( A2_IS_REAL(mtx) ) {
   row = mtx->entries + irow*mtx->inc1 ;
} else if ( A2_IS_COMPLEX(mtx) ) {
   row = mtx->entries + 2*irow*mtx->inc1 ;
} else {
   fprintf(stderr, "\n fatal error in A2_row(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX",
           mtx, irow, mtx->type) ;
   exit(-1) ;
}
return(row) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   return a pointer to the first entry in a column

   created -- 98may01, cca
   -----------------------------------------------
*/
double *
A2_column ( 
   A2   *mtx,
   int   jcol
) {
double   *col ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_column(%p,%d)"
           "\n bad input\n", mtx, jcol) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in A2_column(%p,%d)"
           "\n bad structure, entries is NULL\n", mtx, jcol) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_column(%p,%d)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, jcol, jcol, mtx->n2) ;
   exit(-1) ;
}
if ( A2_IS_REAL(mtx) ) {
   col = mtx->entries + jcol*mtx->inc2 ;
} else if ( A2_IS_COMPLEX(mtx) ) {
   col = mtx->entries + 2*jcol*mtx->inc2 ;
} else {
   fprintf(stderr, "\n fatal error in A2_col(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX",
           mtx, jcol, mtx->type) ;
   exit(-1) ;
}
return(col) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   fill *pValue with the entry in (irow, jcol)

   created -- 98may01, cca
   -------------------------------------------
*/
void
A2_realEntry ( 
   A2       *mtx,
   int      irow,
   int      jcol,
   double   *pValue
) {
int   loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pValue == NULL ) {
   fprintf(stderr, "\n fatal error in A2_realEntry(%p,%d,%d,%p)"
           "\n bad input\n", mtx, irow, jcol, pValue) ;
   exit(-1) ;
}
if ( ! A2_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_realEntry(%p,%d,%d,%p)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, irow, jcol, pValue, mtx->type) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in A2_realEntry(%p,%d,%d,%p)"
           "\n bad structure, entries is NULL\n", 
           mtx, irow, jcol, pValue) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_realEntry(%p,%d,%d,%p)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, jcol, pValue, irow, mtx->n1) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_realEntry(%p,%d,%d,%p)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, irow, jcol, pValue, jcol, mtx->n2) ;
   exit(-1) ;
}
loc = irow*mtx->inc1 + jcol*mtx->inc2 ;
*pValue = mtx->entries[loc] ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   fill (*pReal,*pImag) with the entry in (irow, jcol)

   created -- 98may01, cca
   ---------------------------------------------------
*/
void
A2_complexEntry ( 
   A2       *mtx,
   int      irow,
   int      jcol,
   double   *pReal,
   double   *pImag
) {
int   loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pReal == NULL || pImag == NULL ) {
   fprintf(stderr, "\n fatal error in A2_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad input\n", mtx, irow, jcol, pReal, pImag) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, jcol, pReal, pImag, mtx->type) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in A2_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad structure, entries is NULL\n", 
           mtx, irow, jcol, pReal, pImag) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, jcol, pReal, pImag, irow, mtx->n1) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, irow, jcol, pReal, pImag, jcol, mtx->n2) ;
   exit(-1) ;
}
loc = 2*(irow*mtx->inc1 + jcol*mtx->inc2) ;
*pReal = mtx->entries[loc] ;
*pImag = mtx->entries[loc+1] ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   set the entry in (irow, jcol) to be value

   created -- 98may01, cca
   -----------------------------------------
*/
void
A2_setRealEntry ( 
   A2       *mtx,
   int      irow,
   int      jcol,
   double   value
) {
int   loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_setRealEntry(%p,%d,%d,%f)"
           "\n bad input\n", mtx, irow, jcol, value) ;
   exit(-1) ;
}
if ( ! A2_IS_REAL(mtx) ) {
   fprintf(stderr, 
           "\n fatal error in A2_setRealEntry(%p,%d,%d,%f)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, irow, jcol, value, mtx->type) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_setRealEntry(%p,%d,%d,%f)"
           "\n bad structure, entries is NULL\n", 
           mtx, irow, jcol, value) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, 
           "\n fatal error in A2_setRealEntry(%p,%d,%d,%f)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, jcol, value, irow, mtx->n1) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, 
           "\n fatal error in A2_setRealEntry(%p,%d,%d,%f)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, irow, jcol, value, jcol, mtx->n2) ;
   exit(-1) ;
}
loc = irow*mtx->inc1 + jcol*mtx->inc2 ;
mtx->entries[loc] = value ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   set the entry in (irow, jcol) to be (real,imag)

   created -- 98may01, cca
   -----------------------------------------------
*/
void
A2_setComplexEntry ( 
   A2       *mtx,
   int      irow,
   int      jcol,
   double   real,
   double   imag
) {
int   loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_setComplexEntry(%p,%d,%d,%f,%f)"
           "\n bad input\n", mtx, irow, jcol, real, imag) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, 
           "\n fatal error in A2_setComplexEntry(%p,%d,%d,%f,%f)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, jcol, real, imag, mtx->type) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_setComplexEntry(%p,%d,%d,%f,%f)"
           "\n bad structure, entries is NULL\n", 
           mtx, irow, jcol, real, imag) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, 
           "\n fatal error in A2_setComplexEntry(%p,%d,%d,%f,%f)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, jcol, real, imag, irow, mtx->n1) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, 
           "\n fatal error in A2_setComplexEntry(%p,%d,%d,%f,%f)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, irow, jcol, real, imag, jcol, mtx->n2) ;
   exit(-1) ;
}
loc = 2*(irow*mtx->inc1 + jcol*mtx->inc2) ;
mtx->entries[loc]   = real ;
mtx->entries[loc+1] = imag ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   fill pointers to the matrix first entry
   in row irow and column jcol

   created -- 98may01, cca
   ---------------------------------------
*/
void
A2_pointerToRealEntry ( 
   A2       *mtx,
   int      irow,
   int      jcol,
   double   **ppValue
) {
int   loc  ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || ppValue == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToRealEntry(%p,%d,%d,%p)"
           "\n bad input\n", mtx, irow, jcol, ppValue) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToRealEntry(%p,%d,%d,%p)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, jcol, ppValue, mtx->type) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToRealEntry(%p,%d,%d,%p)"
           "\n bad structure, entries is NULL\n", 
           mtx, irow, jcol, ppValue) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToRealEntry(%p,%d,%d,%p)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, jcol, ppValue, irow, mtx->n1) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToRealEntry(%p,%d,%d,%p)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, irow, jcol, ppValue, jcol, mtx->n2) ;
   exit(-1) ;
}
loc = irow*mtx->inc1 + jcol*mtx->inc2 ;
*ppValue = mtx->entries + loc ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   fill pointers to the matrix first entry
   in row irow and column jcol

   created -- 98may01, cca
   ---------------------------------------
*/
void
A2_pointerToComplexEntry ( 
   A2       *mtx,
   int      irow,
   int      jcol,
   double   **ppReal,
   double   **ppImag
) {
int   loc  ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || ppReal == NULL || ppImag == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad input\n", mtx, irow, jcol, ppReal, ppImag) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, jcol, ppReal, ppImag, mtx->type) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad structure, entries is NULL\n", 
           mtx, irow, jcol, ppReal, ppImag) ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad input, irow = %d, n1 = %d\n", 
           mtx, irow, jcol, ppReal, ppImag, irow, mtx->n1) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, 
           "\n fatal error in A2_pointerToComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad input, jcol = %d, n2 = %d\n", 
           mtx, irow, jcol, ppReal, ppImag, jcol, mtx->n2) ;
   exit(-1) ;
}
loc = 2*(irow*mtx->inc1 + jcol*mtx->inc2) ;
*ppReal = mtx->entries + loc ;
*ppImag = mtx->entries + loc + 1 ;

return ; }

/*--------------------------------------------------------------------*/
