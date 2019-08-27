/*  norms.c  */

#include "../A2.h"

#define CAUTIOUS 1

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   return the entry of maximum magnitude
 
   created -- 98apr15, cca
   -------------------------------------
*/
double
A2_maxabs (
   A2   *a
) {
double   maxval, val ;
double   *entries, *row ;
int      inc1, inc2, irow, jcol, kk, n1, n2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  a == NULL 
   || (n1 = a->n1) < 0
   || (n2 = a->n2) < 0
   || (inc1 = a->inc1) < 0
   || (inc2 = a->inc2) < 0 ) {
   fprintf(stderr, "\n fatal error in A2_maxabs(%p)"
           "\n bad input\n", a ) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(a) || A2_IS_COMPLEX(a)) ) {
   fprintf(stderr, "\n fatal error in A2_maxabs(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           a, a->type) ;
   exit(-1) ;
}
entries = a->entries ;
maxval = 0.0 ;
row = entries ;
if ( A2_IS_REAL(a) ) {
   for ( irow = 0 ; irow < n1 ; irow++ ) {
      for ( jcol = 0, kk = 0 ; jcol < n2 ; jcol++, kk += inc2 ) {
         val = fabs(row[kk]) ;
         if ( maxval < val ) {
            maxval = val ;
         }
      }
      row += inc1 ;
   }
} else if ( A2_IS_COMPLEX(a) ) {
   for ( irow = 0 ; irow < n1 ; irow++ ) {
      for ( jcol = 0, kk = 0 ; jcol < n2 ; jcol++, kk += inc2 ) {
         val = Zabs(row[2*kk], row[2*kk+1]) ;
         if ( maxval < val ) {
            maxval = val ;
         }
      }
      row += inc1 ;
   }
}
return(maxval) ; }
   
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   return the frobenius norm of the matrix

   created -- 98apr15, cca
   ---------------------------------------
*/
double
A2_frobNorm (
   A2   *mtx
) {
double   norm ;
int      ncol, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in A2_frobNorm(%p) "
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_frobNorm(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
if ( (nrow = mtx->n1) <= 0 || (ncol = mtx->n2) <= 0 ) {
   return(0.0) ;
}
norm  = 0 ;
if ( A2_IS_REAL(mtx) ) {
   if ( mtx->inc1 == 1 ) {
      double   *col ;
      int      inc2 = mtx->inc2, irow, jcol ;
   
      for ( jcol = 0, col = mtx->entries ; 
            jcol < ncol ; 
            jcol++, col += inc2 ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            norm += col[irow]*col[irow] ;
         }
      }
   } else if ( mtx->inc2 == 1 ) {
      double   *row ;
      int      inc1 = mtx->inc1, irow, jcol ;
   
      for ( irow = 0, row = mtx->entries ; 
            irow < nrow ; 
            irow++, row += inc1 ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            norm += row[jcol]*row[jcol] ;
         }
      }
   } else {
      double   *entries = mtx->entries ;
      int      inc1 = mtx->inc1, inc2 = mtx->inc2, irow, jcol, loc ;
      
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            loc = irow*inc1 + jcol*inc2 ;
            norm += entries[loc]*entries[loc] ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   if ( mtx->inc1 == 1 ) {
      double   *col ;
      int      inc2 = mtx->inc2, irow, jcol ;
   
      for ( jcol = 0, col = mtx->entries ; 
            jcol < ncol ; 
            jcol++, col += 2*inc2 ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            norm += col[2*irow]*col[2*irow]
                  + col[2*irow+1]*col[2*irow+1] ;
         }
      }
   } else if ( mtx->inc2 == 1 ) {
      double   *row ;
      int      inc1 = mtx->inc1, irow, jcol ;
   
      for ( irow = 0, row = mtx->entries ; 
            irow < nrow ; 
            irow++, row += 2*inc1 ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            norm += row[2*jcol]*row[2*jcol] 
                  + row[2*jcol+1]*row[2*jcol+1] ;
         }
      }
   } else {
      double   *entries = mtx->entries ;
      int      inc1 = mtx->inc1, inc2 = mtx->inc2, irow, jcol, loc ;
      
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            loc = irow*inc1 + jcol*inc2 ;
            norm += entries[2*loc]*entries[2*loc]
                  + entries[2*loc+1]*entries[2*loc+1] ;
         }
      }
   }
}
norm = sqrt(norm) ;

return(norm) ; }
   
/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   return the one-norm of the matrix

   created -- 98apr15, cca
   ---------------------------------
*/
double
A2_oneNorm (
   A2   *mtx
) {
double   norm ;
int      ncol, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_oneNorm(%p) "
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_oneNorm(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
if ( (nrow = mtx->n1) <= 0 || (ncol = mtx->n2) <= 0 ) {
   return(0.0) ;
}
norm = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   if ( mtx->inc1 == 1 ) {
      double   sum ;
      double   *col ;
      int      inc2 = mtx->inc2, irow, jcol ;
   
      for ( jcol = 0, col = mtx->entries ; 
            jcol < ncol ; 
            jcol++, col += inc2 ) {
         sum = 0.0 ;
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            sum += fabs(col[irow]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   } else {
      double   sum ;
      double   *col ;
      int      inc1 = mtx->inc1, inc2 = mtx->inc2, irow, jcol, kk ;
   
      for ( jcol = 0, col = mtx->entries ; 
            jcol < ncol ; 
            jcol++, col += inc2 ) {
         sum = 0.0 ;
         for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
            sum += fabs(col[kk]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   if ( mtx->inc1 == 1 ) {
      double   sum ;
      double   *col ;
      int      inc2 = mtx->inc2, irow, jcol ;
   
      for ( jcol = 0, col = mtx->entries ; 
            jcol < ncol ; 
            jcol++, col += 2*inc2 ) {
         sum = 0.0 ;
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            sum += Zabs(col[2*irow], col[2*irow+1]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   } else {
      double   sum ;
      double   *col ;
      int      inc1 = mtx->inc1, inc2 = mtx->inc2, irow, jcol, kk ;
   
      for ( jcol = 0, col = mtx->entries ; 
            jcol < ncol ; 
            jcol++, col += 2*inc2 ) {
         sum = 0.0 ;
         for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
            sum += Zabs(col[2*kk], col[2*kk+1]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   }
}   
return(norm) ; }
      
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   return the infinity-norm of the matrix

   created -- 98apr15, cca
   --------------------------------------
*/
double
A2_infinityNorm (
   A2   *mtx
) {
double   norm ;
int      ncol, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_infinityNorm(%p) "
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_infinityNorm(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
if ( (nrow = mtx->n1) <= 0 || (ncol = mtx->n2) <= 0 ) {
   return(0.0) ;
}
norm = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   if ( mtx->inc2 == 1 ) {
      double   sum ;
      double   *row = mtx->entries ;
      int      inc1 = mtx->inc1, irow, jcol ;
   
      for ( irow = 0 ; irow < nrow ; irow++, row += inc1 ) {
         for ( jcol = 0, sum = 0.0 ; jcol < ncol ; jcol++ ) {
            sum += fabs(row[jcol]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   } else {
      double   sum ;
      double   *row = mtx->entries ;
      int      inc1 = mtx->inc1, inc2 = mtx->inc2, irow, jcol, kk ;
   
      for ( irow = 0 ; irow < nrow ; irow++, row += inc1 ) {
         sum = 0.0 ;
         for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
            sum += fabs(row[kk]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   if ( mtx->inc2 == 1 ) {
      double   sum ;
      double   *row ;
      int      inc1 = mtx->inc1, irow, jcol ;
   
      for ( irow = 0, row = mtx->entries ; 
         irow < nrow ; 
            irow++, row += 2*inc1 ) {
         sum = 0.0 ;
         for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
            sum += Zabs(row[2*jcol], row[2*jcol+1]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   } else {
      double   sum ;
      double   *row ;
      int      inc1 = mtx->inc1, inc2 = mtx->inc2, irow, jcol, kk ;
   
      for ( irow = 0, row = mtx->entries ; 
            irow < nrow ; 
            irow++, row += 2*inc1 ) {
         sum = 0.0 ;
         for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
            sum += Zabs(row[2*kk], row[2*kk+1]) ;
         }
         if ( norm < sum ) {
            norm = sum ;
         }
      }
   }
}
return(norm) ; }
   
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the one-norm of column jcol

   created -- 98apr15, cca
   ----------------------------------
*/
double
A2_oneNormOfColumn (
   A2   *mtx,
   int   jcol 
) {
double   sum ;
double   *col ;
int      inc1, irow, kk, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || mtx->entries == NULL 
   || jcol < 0 || jcol > mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_oneNormOfColumn(%p,%d)"
           "\n bad input\n", mtx, jcol) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_oneNormOfColumn(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, jcol, mtx->type) ;
   exit(-1) ;
}
nrow = mtx->n1 ;
sum  = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   col  = mtx->entries + jcol*mtx->inc2 ;
   if ( (inc1 = mtx->inc1) == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         sum += fabs(col[irow]) ;
      }
   } else {
      for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
         sum += fabs(col[kk]) ;
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   col  = mtx->entries + 2*jcol*mtx->inc2 ;
   if ( (inc1 = mtx->inc1) == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         sum += Zabs(col[2*irow], col[2*irow+1]) ;
      }
   } else {
      for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
         sum += Zabs(col[2*kk], col[2*kk+1]) ;
      }
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the two-norm of column jcol

   created -- 98apr15, cca
   ----------------------------------
*/
double
A2_twoNormOfColumn (
   A2   *mtx,
   int   jcol 
) {
double   sum ;
double   *col ;
int      inc1, irow, kk, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || mtx->entries == NULL 
   || jcol < 0 || jcol > mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_twoNormOfColumn(%p,%d)"
           "\n bad input\n", mtx, jcol) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_twoNormOfColumn(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, jcol, mtx->type) ;
   exit(-1) ;
}
nrow = mtx->n1 ;
sum  = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   col  = mtx->entries + jcol*mtx->inc2 ;
   if ( (inc1 = mtx->inc1) == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         sum += col[irow]*col[irow] ;
      }
   } else {
      for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
         sum += col[kk]*col[kk] ;
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   col  = mtx->entries + 2*jcol*mtx->inc2 ;
   if ( (inc1 = mtx->inc1) == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         sum += col[2*irow]*col[2*irow] + col[2*irow+1]*col[2*irow+1] ;
      }
   } else {
      for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
         sum += col[2*kk]*col[2*kk] + col[2*kk+1]*col[2*kk+1] ;
      }
   }
}
sum = sqrt(sum) ;

return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   return the infinity-norm of column jcol

   created -- 98apr15, cca
   ---------------------------------------
*/
double
A2_infinityNormOfColumn (
   A2   *mtx,
   int   jcol 
) {
double   norm, val ;
double   *col ;
int      inc1, irow, kk, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || mtx->entries == NULL 
   || jcol < 0 || jcol > mtx->n2 ) {
   fprintf(stderr, "\n fatal error in A2_infinityNormOfColumn(%p,%d)"
           "\n bad input\n", mtx, jcol) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_infinityNormOfColumn(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, jcol, mtx->type) ;
   exit(-1) ;
}
nrow = mtx->n1 ;
norm = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   col  = mtx->entries + jcol*mtx->inc2 ;
   if ( (inc1 = mtx->inc1) == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         val = fabs(col[irow]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   } else {
      for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
         val = fabs(col[kk]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   col  = mtx->entries + 2*jcol*mtx->inc2 ;
   if ( (inc1 = mtx->inc1) == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         val = Zabs(col[2*irow], col[2*irow+1]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   } else {
      for ( irow = 0, kk = 0 ; irow < nrow ; irow++, kk += inc1 ) {
         val = Zabs(col[2*kk], col[2*kk+1]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   }
}
return(norm) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   return the one-norm of row irow

   created -- 98apr15, cca
   -------------------------------
*/
double
A2_oneNormOfRow (
   A2   *mtx,
   int   irow 
) {
double   sum ;
double   *row ;
int      inc2, jcol, kk, ncol ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || mtx->entries == NULL 
   || irow < 0 || irow > mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_oneNormOfRow(%p,%d)"
           "\n bad input\n", mtx, irow) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_oneNormOfRow(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, irow, mtx->type) ;
   exit(-1) ;
}
ncol = mtx->n2 ;
sum  = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   row  = mtx->entries + irow*mtx->inc1 ;
   if ( (inc2 = mtx->inc2) == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         sum += fabs(row[jcol]) ;
      }
   } else {
      for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
         sum += fabs(row[kk]) ;
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   row  = mtx->entries + 2*irow*mtx->inc1 ;
   if ( (inc2 = mtx->inc2) == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         sum += Zabs(row[2*jcol], row[2*jcol+1]) ;
      }
   } else {
      for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
         sum += Zabs(row[2*kk], row[2*kk+1]) ;
      }
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   return the two-norm of row irow

   created -- 98apr15, cca
   -------------------------------
*/
double
A2_twoNormOfRow (
   A2   *mtx,
   int   irow 
) {
double   sum ;
double   *row ;
int      inc2, jcol, kk, ncol ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || mtx->entries == NULL 
   || irow < 0 || irow > mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_twoNormOfRow(%p,%d)"
           "\n bad input\n", mtx, irow) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_twoNormOfRow(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, irow, mtx->type) ;
   exit(-1) ;
}
ncol = mtx->n2 ;
sum  = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   row  = mtx->entries + irow*mtx->inc1 ;
   if ( (inc2 = mtx->inc2) == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         sum += row[jcol]*row[jcol] ;
      }
   } else {
      for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
         sum += row[kk]*row[kk] ;
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   row  = mtx->entries + 2*irow*mtx->inc1 ;
   if ( (inc2 = mtx->inc2) == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         sum += row[2*jcol]*row[2*jcol] + row[2*jcol+1]*row[2*jcol+1] ;
      }
   } else {
      for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
         sum += row[2*kk]*row[2*kk] + row[2*kk+1]*row[2*kk+1] ;
      }
   }
}
sum = sqrt(sum) ;

return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   return the infinity-norm of row irow

   created -- 98apr15, cca
   ------------------------------------
*/
double
A2_infinityNormOfRow (
   A2   *mtx,
   int   irow 
) {
double   norm, val ;
double   *row ;
int      inc2, jcol, kk, ncol ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || mtx->entries == NULL 
   || irow < 0 || irow > mtx->n1 ) {
   fprintf(stderr, "\n fatal error in A2_infinityNormOfRow(%p,%d)"
           "\n bad input\n", mtx, irow) ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_infinityNormOfRow(%p,%d)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, irow, mtx->type) ;
   exit(-1) ;
}
ncol = mtx->n2 ;
norm = 0.0 ;
if ( A2_IS_REAL(mtx) ) {
   row  = mtx->entries + irow*mtx->inc1 ;
   if ( (inc2 = mtx->inc2) == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         val = fabs(row[jcol]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   } else {
      for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
         val = fabs(row[kk]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   row  = mtx->entries + 2*irow*mtx->inc1 ;
   if ( (inc2 = mtx->inc2) == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         val = Zabs(row[2*jcol], row[2*jcol+1]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   } else {
      for ( jcol = 0, kk = 0 ; jcol < ncol ; jcol++, kk += inc2 ) {
         val = Zabs(row[2*kk], row[2*kk+1]) ;
         if ( norm < val ) {
            norm = val ;
         }
      }
   }
}
return(norm) ; }

/*--------------------------------------------------------------------*/
