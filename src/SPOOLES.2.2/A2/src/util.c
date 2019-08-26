/*  util.c  */

#include "../A2.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 98may01, cca
   ----------------------------------------------
*/
int
A2_sizeOf (
   A2   *mtx
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_sizeOf(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_sizeOf(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
if ( A2_IS_REAL(mtx) ) {
   nbytes = sizeof(struct _A2) + mtx->nowned*sizeof(double) ;
} else if ( A2_IS_COMPLEX(mtx) ) {
   nbytes = sizeof(struct _A2) + 2*mtx->nowned*sizeof(double) ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   shift the base of the entries and adjust dimensions

   mtx(0:n1-rowoff-1,0:n2-coloff-1) = mtx(rowoff:n1-1,coloff:n2-1) 

   created -- 98may01, cca
   ---------------------------------------------------------------
*/
void
A2_shiftBase (
   A2   *mtx,
   int   rowoff,
   int   coloff
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_shiftbase(%p,%d,%d)"
           "\n bad input\n", mtx, rowoff, coloff) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_shiftBase(%p,%d,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, rowoff, coloff, mtx->type) ;
   exit(-1) ;
}
mtx->n1 -= rowoff ;
mtx->n2 -= coloff ;
if ( A2_IS_REAL(mtx) ) {
   mtx->entries += rowoff*mtx->inc1 + coloff*mtx->inc2 ;
} else if ( A2_IS_COMPLEX(mtx) ) {
   mtx->entries += 2*(rowoff*mtx->inc1 + coloff*mtx->inc2) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   returns 1 if the storage is row major, otherwise returns zero.

   created -- 98may01, cca
   --------------------------------------------------------------
*/
int
A2_rowMajor ( 
   A2   *mtx 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_rowMajor(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_rowMajor(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
if ( mtx->inc2 == 1 ) {
   return(1) ;
} else {
   return(0) ;
} }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   returns 1 if the storage is column major, otherwise returns zero.

   created -- 98may01, cca
   -----------------------------------------------------------------
*/
int
A2_columnMajor ( 
   A2   *mtx 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_columnMajor(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_columnMajor(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
if ( mtx->inc1 == 1 ) {
   return(1) ;
} else {
   return(0) ;
} }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   transpose the matrix
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_transpose (
   A2   *mtx
) {
int   inc1, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_transpose(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_transpose(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
n1        = mtx->n1   ;
mtx->n1   = mtx->n2   ;
mtx->n2   = n1        ;
inc1      = mtx->inc1 ;
mtx->inc1 = mtx->inc2 ;
mtx->inc2 = inc1      ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   extract row[*] = mtx(irow,*)

   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractRow ( 
   A2      *mtx, 
   double   row[],
   int      irow 
) {
double   *entries ;
int      inc2, j, k, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || row == NULL || mtx->entries == NULL
   || irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_extractRow(%p,%p,%d)"
           "\n bad input\n", mtx, row, irow) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_extractRow(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, row, irow, mtx->type) ;
   exit(-1) ;
}
k       = irow * mtx->inc1 ;
n2      = mtx->n2   ;
inc2    = mtx->inc2 ;
entries = mtx->entries ;
if ( A2_IS_REAL(mtx) ) {
   for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
      row[j] = entries[k] ;
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
      row[2*j]   = entries[2*k] ;
      row[2*j+1] = entries[2*k+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   extract col[*] = mtx(*,jcol)

   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractColumn ( 
   A2      *mtx, 
   double   col[],
   int      jcol 
) {
double   *entries ;
int      i, inc1, k, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || col == NULL || mtx->entries == NULL
   || jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_extractColumn(%p,%p,%d)"
           "\n bad input\n", mtx, col, jcol) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_extractColumn(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, col, jcol, mtx->type) ;
   exit(-1) ;
}
k       = jcol * mtx->inc2 ;
n1      = mtx->n1   ;
inc1    = mtx->inc1 ;
entries = mtx->entries ;
if ( A2_IS_REAL(mtx) ) {
   for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
      col[i] = entries[k] ;
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
      col[2*i]   = entries[2*k]   ;
      col[2*i+1] = entries[2*k+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set mtx(irow,*) = y[*]

   created -- 98may01, cca
   -----------------------
*/
void
A2_setRow ( 
   A2      *mtx, 
   double   row[],
   int      irow 
) {
double   *entries ;
int      inc2, j, k, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || row == NULL || irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_setRow(%p,%p,%d)"
           "\n bad input\n", mtx, row, irow) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_setRow(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, row, irow, mtx->type) ;
   exit(-1) ;
}
k       = irow * mtx->inc1 ;
n2      = mtx->n2   ;
inc2    = mtx->inc2 ;
entries = mtx->entries ;
if ( A2_IS_REAL(mtx) ) {
   for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
      entries[k] = row[j] ;
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
      entries[2*k]   = row[2*j]   ;
      entries[2*k+1] = row[2*j+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set mtx(*,jcol) = y[*]

   created -- 98may01, cca
   -----------------------
*/
void
A2_setColumn ( 
   A2      *mtx, 
   double   col[],
   int      jcol 
) {
double   *entries ;
int      inc1, i, k, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || col == NULL || jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_setColumn(%p,%p,%d)"
           "\n bad input\n", mtx, col, jcol) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_setColumn(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, col, jcol, mtx->type) ;
   exit(-1) ;
}
k       = jcol * mtx->inc2 ;
n1      = mtx->n1   ;
inc1    = mtx->inc1 ;
entries = mtx->entries ;
if ( A2_IS_REAL(mtx) ) {
   for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
      entries[k] = col[i] ;
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
      entries[2*k]   = col[2*i]   ;
      entries[2*k+1] = col[2*i+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   extract row[*] = mtx(irow,*)

   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractRowDV ( 
   A2   *mtx, 
   DV    *rowDV,
   int   irow 
) {
double   *entries, *row ;
int      inc2, j, k, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || rowDV == NULL || mtx->entries == NULL
   || irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_extractRowDV(%p,%p,%d)"
           "\n bad input\n", mtx, rowDV, irow) ;
   exit(-1) ;
}
if ( ! A2_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_extractRowDV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, rowDV, irow, mtx->type) ;
   exit(-1) ;
}
if ( DV_size(rowDV) < (n2 = mtx->n2) ) {
   DV_setSize(rowDV, n2) ;
}
row     = DV_entries(rowDV) ;
k       = irow * mtx->inc1 ;
inc2    = mtx->inc2 ;
entries = mtx->entries ;
for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
   row[j] = entries[k] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   extract row[*] = mtx(irow,*)

   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractRowZV ( 
   A2   *mtx, 
   ZV    *rowZV,
   int   irow 
) {
double   *entries, *row ;
int      inc2, j, k, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || rowZV == NULL || mtx->entries == NULL
   || irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_extractRowZV(%p,%p,%d)"
           "\n bad input\n", mtx, rowZV, irow) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_extractRowZV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, rowZV, irow, mtx->type) ;
   exit(-1) ;
}
if ( ZV_size(rowZV) < (n2 = mtx->n2) ) {
   ZV_setSize(rowZV, n2) ;
}
row     = ZV_entries(rowZV) ;
k       = irow * mtx->inc1 ;
inc2    = mtx->inc2 ;
entries = mtx->entries ;
for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
   row[2*j]   = entries[2*k]   ;
   row[2*j+1] = entries[2*k+1] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   extract col[*] = mtx(*,jcol)

   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractColumnDV ( 
   A2   *mtx, 
   DV    *colDV,
   int   jcol 
) {
double   *entries, *col ;
int      i, inc1, k, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || colDV == NULL || mtx->entries == NULL
   || jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_extractColumnDV(%p,%p,%d)"
           "\n bad input\n", mtx, colDV, jcol) ;
   exit(-1) ;
}
if ( ! A2_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_extractColumnDV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, colDV, jcol, mtx->type) ;
   exit(-1) ;
}
if ( DV_size(colDV) < (n1 = mtx->n1) ) {
   DV_setSize(colDV, n1) ;
}
col     = DV_entries(colDV) ;
k       = jcol * mtx->inc2 ;
inc1    = mtx->inc1 ;
entries = mtx->entries ;
for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
   col[i] = entries[k] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   extract col[*] = mtx(*,jcol)

   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractColumnZV ( 
   A2   *mtx, 
   ZV    *colZV,
   int   jcol 
) {
double   *entries, *col ;
int      i, inc1, k, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || colZV == NULL || mtx->entries == NULL
   || jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_extractColumnZV(%p,%p,%d)"
           "\n bad input\n", mtx, colZV, jcol) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_extractColumnZV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, colZV, jcol, mtx->type) ;
   exit(-1) ;
}
if ( ZV_size(colZV) < (n1 = mtx->n1) ) {
   ZV_setSize(colZV, n1) ;
}
col     = ZV_entries(colZV) ;
k       = jcol * mtx->inc2 ;
inc1    = mtx->inc1 ;
entries = mtx->entries ;
for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
   col[2*i]   = entries[2*k]   ;
   col[2*i+1] = entries[2*k+1] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set mtx(irow,*) = y[*]

   created -- 98may01, cca
   -----------------------
*/
void
A2_setRowDV ( 
   A2      *mtx, 
   DV       *rowDV,
   int      irow 
) {
double   *entries, *row ;
int      inc2, j, k, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || rowDV == NULL || DV_size(rowDV) != (n2 = mtx->n2)
   || irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_setRowDV(%p,%p,%d)"
           "\n bad input\n", mtx, rowDV, irow) ;
   exit(-1) ;
}
if ( ! A2_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_setRowDV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, rowDV, irow, mtx->type) ;
   exit(-1) ;
}
k       = irow * mtx->inc1 ;
inc2    = mtx->inc2 ;
entries = mtx->entries ;
row     = DV_entries(rowDV) ;
for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
   entries[k] = row[j] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set mtx(irow,*) = y[*]

   created -- 98may01, cca
   -----------------------
*/
void
A2_setRowZV ( 
   A2    *mtx, 
   ZV    *rowZV,
   int   irow 
) {
double   *entries, *row ;
int      inc2, j, k, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || rowZV == NULL || ZV_size(rowZV) != (n2 = mtx->n2)
   || irow < 0 || irow >= mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_setRowZV(%p,%p,%d)"
           "\n bad input\n", mtx, rowZV, irow) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_setRowZV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, rowZV, irow, mtx->type) ;
   exit(-1) ;
}
k       = irow * mtx->inc1 ;
inc2    = mtx->inc2 ;
entries = mtx->entries ;
row     = ZV_entries(rowZV) ;
for ( j = 0 ; j < n2 ; j++, k += inc2 ) {
   entries[2*k]   = row[2*j]   ;
   entries[2*k+1] = row[2*j+1] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set mtx(*,jcol) = y[*]

   created -- 98may01, cca
   -----------------------
*/
void
A2_setColumnDV ( 
   A2      *mtx, 
   DV       *colDV,
   int      jcol 
) {
double   *col, *entries ;
int      inc1, i, k, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || colDV == NULL || DV_size(colDV) != (n1 = mtx->n1)
   || jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_setColumnDV(%p,%p,%d)"
           "\n bad input\n", mtx, colDV, jcol) ;
   exit(-1) ;
}
if ( ! A2_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_setColumnDV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, colDV, jcol, mtx->type) ;
   exit(-1) ;
}
k       = jcol * mtx->inc2 ;
inc1    = mtx->inc1 ;
entries = mtx->entries ;
col     = DV_entries(colDV) ;
for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
   entries[k] = col[i] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set mtx(*,jcol) = y[*]

   created -- 98may01, cca
   -----------------------
*/
void
A2_setColumnZV ( 
   A2      *mtx, 
   ZV       *colZV,
   int      jcol 
) {
double   *col, *entries ;
int      inc1, i, k, n1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || colZV == NULL || ZV_size(colZV) != (n1 = mtx->n1)
   || jcol < 0 || jcol >= mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_setColumnZV(%p,%p,%d)"
           "\n bad input\n", mtx, colZV, jcol) ;
   exit(-1) ;
}
if ( ! A2_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in A2_setColumnZV(%p,%p,%d)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, colZV, jcol, mtx->type) ;
   exit(-1) ;
}
k       = jcol * mtx->inc2 ;
inc1    = mtx->inc1 ;
entries = mtx->entries ;
col     = ZV_entries(colZV) ;
for ( i = 0 ; i < n1 ; i++, k += inc1 ) {
   entries[2*k]   = col[2*i]   ;
   entries[2*k+1] = col[2*i+1] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   fill the matrix with uniform random numbers in [lower, upper]

   created -- 98may01, cca
   -------------------------------------------------------------
*/
void
A2_fillRandomUniform (
   A2       *a,
   double   lower,
   double   upper,
   int      seed
) {
double   *entries ;
int      i, inc1, inc2, j, loc, n1, n2 ;
Drand    drand ;
/*
   ---------------
   check the input
   ---------------
*/
if ( a == NULL
   || (n1 = a->n1) <= 0
   || (n2 = a->n2) <= 0
   || (inc1 = a->inc1) <= 0
   || (inc2 = a->inc2) <= 0
   || (entries = a->entries) == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_fillRandomUniform(%p,%f,%f,%d)"
           "\n bad input\n",
           a, lower, upper, seed) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_fillRandomUniform(%p,%f,%f,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, lower, upper, seed, a->type) ;
   exit(-1) ;
}
/*
   ----------------
   fill the entries
   ----------------
*/
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setUniform(&drand, lower, upper) ;
for ( j = 0 ; j < n2 ; j++ ) {
   for ( i = 0 ; i < n1 ; i++ ) {
      loc = i*inc1 + j*inc2 ;
      if ( A2_IS_REAL(a) ) {
         entries[loc] = Drand_value(&drand) ;
      } else if ( A2_IS_COMPLEX(a) ) {
         entries[2*loc]   = Drand_value(&drand) ;
         entries[2*loc+1] = Drand_value(&drand) ;
      }
   }
} 
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   fill the matrix with normal(0,1) random numbers

   created -- 98may01, cca
   -----------------------------------------------
*/
void
A2_fillRandomNormal (
   A2      *a,
   double   mean,
   double   variance,
   int      seed
) {
double   *entries ;
int      i, inc1, inc2, j, loc, n1, n2 ;
Drand    drand ;
/*
   ---------------
   check the input
   ---------------
*/
if ( a == NULL
   || (n1 = a->n1) <= 0
   || (n2 = a->n2) <= 0
   || (inc1 = a->inc1) <= 0
   || (inc2 = a->inc2) <= 0
   || (entries = a->entries) == NULL ) {
   fprintf(stderr, "\n fatal error in A2_fillRandomNormal(%p,%d)"
           "\n bad input\n",
           a, seed) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_fillRandomNormal(%p,%f,%f,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, mean, variance, seed, a->type) ;
   exit(-1) ;
}
/*
   ----------------
   fill the entries
   ----------------
*/
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setUniform(&drand, mean, variance) ;
for ( j = 0 ; j < n2 ; j++ ) {
   for ( i = 0 ; i < n1 ; i++ ) {
      loc = i*inc1 + j*inc2 ;
      if ( A2_IS_REAL(a) ) {
         entries[loc] = Drand_value(&drand) ;
      } else if ( A2_IS_COMPLEX(a) ) {
         entries[2*loc]   = Drand_value(&drand) ;
         entries[2*loc+1] = Drand_value(&drand) ;
      }
   }
} 

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   fill the matrix with the identity matrix

   created -- 98may01, cca
   ----------------------------------------
*/
void
A2_fillWithIdentity (
   A2   *a
) {
double   *entries ;
int      ii, inc, inc1, inc2, j, n ;
/*
   ---------------
   check the input
   ---------------
*/
if ( a == NULL
   || (n = a->n1) <= 0
   || n != a->n2
   || (inc1 = a->inc1) <= 0
   || (inc2 = a->inc2) <= 0
   || (inc1 != 1 && inc2 != 1)
   || (entries = a->entries) == NULL ) {
   fprintf(stderr, "\n fatal error in A2_fillWithIdentity(%p)"
           "\n bad input\n", a) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_fillWithIdentity(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, a->type) ;
   exit(-1) ;
}
inc = (inc1 == 1) ? inc2 : inc1 ;
A2_zero(a) ;
for ( j = 0, ii = 0 ; j < n ; j++, ii += inc + 1 ) {
   if ( A2_IS_REAL(a) ) {
      entries[ii] = 1.0 ;
   } else if ( A2_IS_COMPLEX(a) ) {
      entries[2*ii] = 1.0 ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   fill the matrix with zeros

   created -- 98may01, cca
   --------------------------
*/
void
A2_zero (
   A2   *a
) {
double   *entries ;
int      i, inc1, inc2, j, loc, n1, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( a == NULL
   || (n1 = a->n1) <= 0
   || (n2 = a->n2) <= 0
   || (inc1 = a->inc1) <= 0
   || (inc2 = a->inc2) <= 0
   || (entries = a->entries) == NULL ) {
   fprintf(stderr, "\n fatal error in A2_zero(%p)"
           "\n bad input\n", a) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_zero(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, a->type) ;
   exit(-1) ;
}
for ( j = 0 ; j < n2 ; j++ ) {
   for ( i = 0 ; i < n1 ; i++ ) {
      loc =i*inc1 + j*inc2 ;
      if ( A2_IS_REAL(a) ) {
         entries[loc] = 0.0 ;
      } else if ( A2_IS_COMPLEX(a) ) {
         entries[2*loc]   = 0.0 ;
         entries[2*loc+1] = 0.0 ;
      }
   }
} 

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   copy one matrix into another
      A := B

   created  -- 98may01, cca
   ----------------------------
*/
void
A2_copy (
   A2   *A,
   A2   *B
) {
double   *entA, *entB ;
int      inc1A, inc1B, inc2A, inc2B, irow, jcol, locA, locB,
         ncol, ncolA, ncolB, nrow, nrowA, nrowB ;
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL
   || (nrowA = A->n1) < 0
   || (ncolA = A->n2) < 0
   || (inc1A = A->inc1) <= 0
   || (inc2A = A->inc2) <= 0
   || (entA = A->entries) == NULL
   || B == NULL
   || (nrowB = B->n1) < 0
   || (ncolB = B->n2) < 0
   || (inc1B = B->inc1) <= 0
   || (inc2B = B->inc2) <= 0 
   || (entB = B->entries) == NULL ) {
   fprintf(stderr, "\n fatal error in A2_copy(%p,%p)"
           "\n bad input\n", A, B) ;
   if ( A != NULL ) {
      fprintf(stderr, "\n\n first A2 object") ;
      A2_writeStats(A, stderr) ;
   }
   if ( B != NULL ) {
      fprintf(stderr, "\n\n second A2 object") ;
      A2_writeStats(B, stderr) ;
   }
   exit(-1) ;
}
if ( ! (A2_IS_REAL(A) || A2_IS_COMPLEX(A)) ) {
   fprintf(stderr, "\n fatal error in A2_copy(%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           A, B, A->type) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(B) || A2_IS_COMPLEX(B)) ) {
   fprintf(stderr, "\n fatal error in A2_copy(%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           A, B, B->type) ;
   exit(-1) ;
}
if ( A->type != B->type ) {
   fprintf(stderr, "\n fatal error in A2_copy(%p,%p)"
           "\n A's type %d, B's type = %d, must be the same\n",
           A, B, A->type, B->type) ;
   exit(-1) ;
}
nrow = (nrowA <= nrowB) ? nrowA : nrowB ;
ncol = (ncolA <= ncolB) ? ncolA : ncolB ;
if ( A2_IS_REAL(A) ) {
   if ( inc1A == 1 && inc1B == 1 ) {
      double   *colA = entA, *colB = entB ;
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            colA[irow] = colB[irow] ;
         }
         colA += inc2A ;
         colB += inc2B ;
      }
   } else if ( inc2A == 1 && inc2B == 1 ) {
      double   *rowA = entA, *rowB = entB ;
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            rowA[jcol] = rowB[jcol] ;
         }
         rowA += 2*inc1A ;
      }
   } else {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            locA = irow*inc1A + jcol*inc2A ;
            locB = irow*inc1B + jcol*inc2B ;
            entA[locA] = entB[locB] ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(A) ) {
   if ( inc1A == 1 && inc1B == 1 ) {
      double   *colA = entA, *colB = entB ;
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            colA[2*irow]   = colB[2*irow]   ;
            colA[2*irow+1] = colB[2*irow+1] ;
         }
         colA += 2*inc2A ;
         colB += 2*inc2B ;
      }
   } else if ( inc2A == 1 && inc2B == 1 ) {
      double   *rowA = entA, *rowB = entB ;
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            rowA[2*jcol]   = rowB[2*jcol]   ;
            rowA[2*jcol+1] = rowB[2*jcol+1] ;
         }
         rowA += 2*inc1A ;
         rowB += 2*inc1B ;
      }
   } else {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            locA = irow*inc1A + jcol*inc2A ;
            locB = irow*inc1B + jcol*inc2B ;
            entA[2*locA]   = entB[2*locB]   ;
            entA[2*locA+1] = entB[2*locB+1] ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------- 
   subtract one matrix from another 

   A := A - B
   
   created -- 98may01, cca
   ----------------------------
*/
void
A2_sub (
   A2   *A,
   A2   *B
) {
double   *entA, *entB ;
int      inc1A, inc1B, inc2A, inc2B, irow, jcol, locA, locB,
         ncol, ncolA, ncolB, nrow, nrowA, nrowB ;
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL
   || B == NULL
   || (nrowA = A->n1) <= 0
   || (ncolA = A->n2) <= 0
   || (inc1A = A->inc1) <= 0
   || (inc2A = A->inc2) <= 0
   || (nrowB = B->n1) <= 0
   || (ncolB = B->n2) <= 0
   || (inc1B = B->inc1) <= 0
   || (inc2B = B->inc2) <= 0 
   || (entA = A->entries) == NULL
   || (entB = B->entries) == NULL ) {
   fprintf(stderr, "\n fatal error in A2_sub(%p,%p)"
           "\n bad input\n", A, B) ;
   if ( A != NULL ) {
      fprintf(stderr, "\n\n first A2 object") ;
      A2_writeStats(A, stderr) ;
   }
   if ( B != NULL ) {
      fprintf(stderr, "\n\n second A2 object") ;
      A2_writeStats(B, stderr) ;
   }
   exit(-1) ;
}
if ( ! (A2_IS_REAL(A) || A2_IS_COMPLEX(A)) ) {
   fprintf(stderr, "\n fatal error in A2_sub(%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           A, B, A->type) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(B) || A2_IS_COMPLEX(B)) ) {
   fprintf(stderr, "\n fatal error in A2_sub(%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           A, B, B->type) ;
   exit(-1) ;
}
if ( A->type != B->type ) {
   fprintf(stderr, "\n fatal error in A2_sub(%p,%p)"
           "\n A's type %d, B's type = %d, must be the same\n",
           A, B, A->type, B->type) ;
   exit(-1) ;
}
/*
fprintf(stdout, "\n debug : A") ;
A2_writeForHumanEye(A, stdout) ;
fprintf(stdout, "\n debug : B") ;
A2_writeForHumanEye(B, stdout) ;
*/
nrow = (nrowA <= nrowB) ? nrowA : nrowB ;
ncol = (ncolA <= ncolB) ? ncolA : ncolB ;
if ( A2_IS_REAL(A) ) {
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         locA = irow*inc1A + jcol*inc2A ;
         locB = irow*inc1B + jcol*inc2B ;
         entA[locA] -= entB[locB] ;
      }
   }
} else if ( A2_IS_COMPLEX(A) ) {
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         locA = irow*inc1A + jcol*inc2A ;
         locB = irow*inc1B + jcol*inc2B ;
         entA[2*locA]   -= entB[2*locB]   ;
         entA[2*locA+1] -= entB[2*locB+1] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   swap two rows of the matrix

   created -- 98may01, cca
   ---------------------------
*/
void
A2_swapRows (
   A2   *a,
   int   irow1,
   int   irow2
) {
double   temp ;
double   *row1, *row2 ;
int      inc2, j, k, n2 ;
/*
   -----------
   check input
   -----------
*/
if (  a == NULL 
   || irow1 < 0 || irow1 >= a->n1
   || irow2 < 0 || irow2 >= a->n1 ) {
   fprintf(stderr, 
           "\n fatal error in A2_swapRows(%p,%d,%d)"
           "\n bad input\n", a, irow1, irow2) ;
   exit(-1) ;
}
if (  a->n1   <= 0
   || a->inc1 <= 0
   || (n2 = a->n2) <= 0
   || (inc2 = a->inc2) <= 0
   || a->entries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_swapRows(%p,%d,%d)"
           "\n bad structure\n", a, irow1, irow2) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_swapRows(%p,%d,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, irow1, irow2, a->type) ;
   exit(-1) ;
}
if ( irow1 == irow2 ) {
   return ;
}
if ( A2_IS_REAL(a) ) {
   row1 = a->entries + irow1*a->inc1 ;
   row2 = a->entries + irow2*a->inc1 ;
   if ( inc2 == 1 ) {
      for ( j = 0 ; j < n2 ; j++ ) {
         temp    = row1[j] ;
         row1[j] = row2[j] ;
         row2[j] = temp    ;
      }
   } else {
      for ( j = 0, k = 0 ; j < n2 ; j++, k += inc2 ) {
         temp    = row1[k] ;
         row1[k] = row2[k] ;
         row2[k] = temp    ;
      }
   }
} else if ( A2_IS_COMPLEX(a) ) {
   row1 = a->entries + 2*irow1*a->inc1 ;
   row2 = a->entries + 2*irow2*a->inc1 ;
   if ( inc2 == 1 ) {
      for ( j = 0 ; j < n2 ; j++ ) {
         temp        = row1[2*j]   ;
         row1[2*j]   = row2[2*j]   ;
         row2[2*j]   = temp        ;
         temp        = row1[2*j+1] ;
         row1[2*j+1] = row2[2*j+1] ;
         row2[2*j+1] = temp        ;
      }
   } else {
      for ( j = 0, k = 0 ; j < n2 ; j++, k += inc2 ) {
         temp        = row1[2*k]   ;
         row1[2*k]   = row2[2*k]   ;
         row2[2*k]   = temp        ;
         temp        = row1[2*k+1] ;
         row1[2*k+1] = row2[2*k+1] ;
         row2[2*k+1] = temp        ;
      }
   }
}
return ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------
   swap two columns of the matrix

   created -- 98may01, cca
   ------------------------------
*/
void
A2_swapColumns (
   A2   *a,
   int   jcol1,
   int   jcol2
) {
double   temp ;
double   *col1, *col2 ;
int      i, inc1, k, n1 ;
/*
   -----------
   check input
   -----------
*/
if (  a == NULL
   || jcol1 < 0 || jcol1 >= a->n2
   || jcol2 < 0 || jcol2 >= a->n2 ) {
   fprintf(stderr,
           "\n fatal error in A2_swapColumns(%p,%d,%d)"
           "\n bad input\n", a, jcol1, jcol2) ;
   exit(-1) ;
}
if (  (n1 = a->n1)   <= 0
   || (inc1 = a->inc1) <= 0
   || a->n2 <= 0
   || a->inc2 <= 0
   || a->entries == NULL ) {
   fprintf(stderr,
           "\n fatal error in A2_swapColumns(%p,%d,%d)"
           "\n bad structure\n", a, jcol1, jcol2) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_swapColumns(%p,%d,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, jcol1, jcol2, a->type) ;
   exit(-1) ;
}
if ( jcol1 == jcol2 ) {
   return ;
}
if ( A2_IS_REAL(a) ) {
   col1 = a->entries + jcol1*a->inc2 ;
   col2 = a->entries + jcol2*a->inc2 ;
   if ( inc1 == 1 ) {
      for ( i = 0 ; i < n1 ; i++ ) {
         temp    = col1[i] ;
         col1[i] = col2[i] ;
         col2[i] = temp    ;
      }
   } else {
      for ( i = 0, k = 0 ; i < n1 ; i++, k += inc1 ) {
         temp    = col1[k] ;
         col1[k] = col2[k] ;
         col2[k] = temp    ;
      }
   }
} else if ( A2_IS_COMPLEX(a) ) {
   col1 = a->entries + 2*jcol1*a->inc2 ;
   col2 = a->entries + 2*jcol2*a->inc2 ;
   if ( inc1 == 1 ) {
      for ( i = 0 ; i < n1 ; i++ ) {
         temp        = col1[2*i]   ;
         col1[2*i]   = col2[2*i]   ;
         col2[2*i]   = temp        ;
         temp        = col1[2*i+1] ;
         col1[2*i+1] = col2[2*i+1] ;
         col2[2*i+1] = temp        ;
      }
   } else {
      for ( i = 0, k = 0 ; i < n1 ; i++, k += inc1 ) {
         temp        = col1[2*k]   ;
         col1[2*k]   = col2[2*k]   ;
         col2[2*k]   = temp        ;
         temp        = col1[2*k+1] ;
         col1[2*k+1] = col2[2*k+1] ;
         col2[2*k+1] = temp        ;
      }
   }
}
return ; }
 
/*--------------------------------------------------------------------*/
