/*  DV.c  */

#include "../Utilities.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to compute y[*] := y[*] + x[*]

   created -- 95sep22, cca
   -----------------------------------------
*/
void
DVadd ( 
   int      size, 
   double   y[], 
   double   x[] 
) {
if ( size <= 0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in DVadd"
           "\n invalid input, size = %d, y = %p, x = %p\n",
           size, y, x) ;
   exit(-1) ; 
} else {
   int i ;
   for ( i = 0 ; i < size ; i++ ) {
      y[i] += x[i] ; 
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to compute y[*] := y[*] + alpha * x[*]

   created -- 95sep22, cca
   -------------------------------------------------
*/
void
DVaxpy ( 
   int      size, 
   double   y[], 
   double   alpha, 
   double   x[] 
) {
if ( size < 0 || alpha == 0.0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in DVaxpy"
           "\n invalid input, size = %d, y = %p, alpha = %f, x = %p\n",
           size, y, alpha, x) ;
   exit(-1) ; 
} else {
   int i ;
   for ( i = 0 ; i < size ; i++ ) {
      y[i] += alpha * x[i] ; 
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to compute y[index[*]] := y[index[*]] + alpha * x[*]
 
   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
void
DVaxpyi ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   alpha, 
   double   x[] 
) {
if ( size <= 0 || alpha == 0.0 ) {
   return ;
} else if ( y == NULL || index == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in DVaxpyi, invalid input"
           "\n size = %d, y = %p, index = %p, alpha = %f, x = %p",
           size, y, index, alpha, x) ;
   exit(-1) ; 
} else {
   int i ;
   for ( i = 0 ; i < size ; i++ ) {
      y[index[i]] += alpha * x[i] ; 
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   purpose -- to copy y[*] := x[*]

   created -- 95sep22, cca
   -------------------------------
*/
void
DVcopy ( 
   int      size, 
   double   y[], 
   double   x[] 
) {
if ( size <= 0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in DVcopy, invalid input"
           "\n size = %d, y = %p, x = %p\n", size, y, x) ;
   exit(-1) ; 
} else {
   int i ;
   for ( i = 0 ; i < size ; i++ ) {
      y[i] = x[i] ; 
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- given the pair of arrays (x1[],y1[]), 
              create a pair of arrays (x2[],y2[]) whose
              entries are pairwise chosen from (x1[],y1[])
              and whose distribution is an approximation.

   return value -- the size of the (x2[],y2[]) arrays

   created -- 95sep22, cca
   -------------------------------------------------------
*/
int
DVcompress ( 
   int      size1, 
   double   x1[], 
   double   y1[],
   int      size2, 
   double   x2[], 
   double   y2[] 
) {
double   delta, dx, dy, path, totalPath ;
double   *ds ;
int      i, j ;
/*
   --------------------
   check the input data
   --------------------
*/
if ( size1 <= 0 || size2 <= 0 ) {
   return(0) ;
} else if ( x1 == NULL || y1 == NULL || x2 == NULL || y2 == NULL ) {
   fprintf(stderr, "\n fatal error in DVcompress, invalid data"
           "\n size1 = %d, x1 = %p, y1 = %p"
           "\n size2 = %d, x2 = %p, y2 = %p",
           size1, x1, y1, size2, x2, y2) ;
   exit(-1) ; 
}
/*
   ----------------------------------------
   compute the path length and its segments
   ----------------------------------------
*/
ds = DVinit(size1, 0.0) ;
for ( j = 1 ; j < size1 ; j++ ) {
   dx = x1[j] - x1[j-1] ;
   dy = y1[j] - y1[j-1] ;
   ds[j-1] = sqrt(dx*dx + dy*dy) ;
}
totalPath = DVsum(size1, ds) ;
delta = totalPath / (size2-2) ;
#if MYDEBUG > 0
fprintf(stdout, "\n totalPath = %12.4e, delta = %12.4e, ds",
        totalPath, delta) ;
DVfprintf(stdout, size1, ds) ;
#endif
/*
   ---------------------
   fill the second array
   ---------------------
*/
i = 0 ;
x2[i] = x1[i] ;
y2[i] = y1[i] ;
i++ ;
path  = 0. ;
for ( j = 1 ; j < size1 - 1 ; j++ ) {
   path += ds[j-1] ;
#if MYDEBUG > 0
   fprintf(stdout, "\n j %d, path %12.4e", j, path) ;
#endif
   if ( path >= delta ) {
#if MYDEBUG > 0
      fprintf(stdout, ", accepted") ;
#endif
      x2[i] = x1[j] ;
      y2[i] = y1[j] ;
      i++ ;
      path  = 0. ;
   }
}
x2[i] = x1[size1-1] ;
y2[i] = y1[size1-1] ;
i++ ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
DVfree(ds) ;

return(i) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to copy sum{y[*] * x[*]}

   created -- 95sep22, cca
   -----------------------------------
*/
double
DVdot ( 
   int      size, 
   double   y[], 
   double   x[] 
) {
double   sum = 0.0 ;
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in DVdot, invalid data"
              "\n size = %d, y = %p, x = %p\n", size, y, x) ;
      exit(-1) ;
   } else {
      int      i ;
      for ( i = 0, sum = 0. ; i < size ; i++ ) {
         sum += y[i] * x[i] ; 
      }
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to perform a indexed dot product
 
   sum = sum_k y[index[k]]*x[k]
 
   where y and x are real
 
   created -- 98may02, cca
   --------------------------------------------
*/
double
DVdoti ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[]
) {
double   sum ;
int      ii  ;
 
if (  size < 0 || y == NULL || index == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in DVdoti(%d,%p,%p,%p)"
           "\n bad input\n", size, y, index, x) ;
   exit(-1) ;
} 
for ( ii = 0, sum = 0.0 ; ii < size ; ii++ ) {
   sum += y[index[ii]] * x[ii] ;
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to fill a double vector with a value

   created -- 95sep22, cca
   -----------------------------------------------
*/
void
DVfill ( 
   int      size, 
   double   y[], 
   double   dval 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVfill, invalid data"
              "\n size = %d, y = %p, dval = %f\n", size, y, dval) ;
      exit(-1) ;
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = dval ; 
      }
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to print out a double vector

   created -- 95sep22, cca
   -----------------------------------------
*/
void
DVfprintf ( 
   FILE     *fp, 
   int      size, 
   double   y[]
) {
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVfprintf, invalid input"
              "\n fp = %p, size = %d, y = %p\n", fp, size, y) ;
      exit(-1) ;
   } else {
      int    i ;
      for ( i = 0 ; i < size ; i++ ) {
         if ( i % 6 == 0 ) fprintf(fp, "\n ") ;
         fprintf(fp, "%12.4e", y[i]) ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to free storage taken by a double vector.
              note : y[] must have been returned by DVinit.

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
DVfree ( 
   double   y[] 
) {
if ( y != NULL ) {
   FREE(y) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to read in a double vector
              return value -- # of entries read

   created -- 95sep22, cca
   --------------------------------------------
*/
int
DVfscanf ( 
   FILE     *fp, 
   int      size, 
   double   y[] 
) {
int    i  = 0, rc ;
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVfscanf, invalid input"
              "\n fp = %p, size = %d, y = %p\n", fp, size, y) ;
      exit(-1) ;
   }
   for ( i = 0 ; i < size ; i++ ) {
      if ( fscanf(fp, " %le", y + i)  != 1 ) {
         break ; 
      } 
   }
}
return(i) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to gather y[*] = x[index[*]]

   created -- 95sep22, cca
   ---------------------------------------
*/
void
DVgather ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVgather, invalid input"
              "\n size = %d, y = %p, x = %p, index = %p\n",
              size, y, x, index) ;
      exit(-1) ;
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = x[index[i]] ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to gather add y[*] += x[index[*]] and zero x[index[*]]

   created -- 95sep22, cca
   -----------------------------------------------------------------
*/
void
DVgatherAddZero ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVgatherAddZero, invalid input"
              "\n size = %d, y = %p, x = %p, index = %p\n",
              size, y, x, index) ;
      exit(-1) ;
   } else {
      int   i, j ;
      for ( i = 0 ; i < size ; i++ ) {
         j    =  index[i] ; 
         y[i] += x[j]     ; 
         x[j] =  0.0      ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to gather y[*] = x[index[*]] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------
*/
void
DVgatherZero ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVgatherZero, invalid input"
              "\n size = %d, y = %p, x = %p, index = %p\n",
              size, y, x, index) ;
      exit(-1) ;
   } else {
      int   i, j ;
      for ( i = 0 ; i < size ; i++ ) {
         j    = index[i] ; 
         y[i] = x[j]     ; 
         x[j] = 0.0      ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a double array with size entries
              and fill with value dval
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
double *
DVinit ( 
   int      size, 
   double   dval 
) {
double   *y = NULL ;
if ( size > 0 ) {
   y = DVinit2(size) ;
   DVfill(size, y, dval) ;
}
return(y) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a double array with size entries
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
double *
DVinit2 ( 
   int   size 
) {
double   *y = NULL ;
if ( size > 0 ) {
   ALLOCATE(y, double, size) ;
}
return(y) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[index[*]] := y[*]

   created -- 95sep22, cca
   ------------------------------
*/
void
DVinvPerm ( 
   int      size, 
   double   y[], 
   int      index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVinvPerm, invalid data"
              "\n size = %d, y = %p, index = %p", size, y, index) ;
      exit(-1) ;
   } else {
      double   *x ;
      int      i ;
      x = DVinit2(size) ;
      DVcopy(size, x, y) ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] = x[i] ; 
      }
      DVfree(x) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of maximum value,
              *ploc contains its index 

   created -- 95sep22, cca
   ------------------------------------------------------
*/
double
DVmax ( 
   int      size, 
   double   y[], 
   int      *ploc 
) {
double   maxval = 0.0 ;
int      loc = -1 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVmax, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
      exit(-1) ;
   } else {
      int      i ;
      maxval = y[0] ;
      loc    = 0    ;
      for ( i = 1 ; i < size ; i++ ) {
         if ( maxval < y[i] ) {
            maxval = y[i] ; 
            loc    = i    ; 
         } 
      }
      *ploc = loc ;
   }
}
*ploc = loc ;

return(maxval) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of maximum absolute value,
              *ploc contains its index 

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
double
DVmaxabs ( 
   int      size, 
   double   y[], 
   int      *ploc 
) {
double   maxval = 0.0 ;
int      loc = -1 ;

if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVmaxabs, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
      exit(-1) ;
   } else {
      int      i   ;
      double   val ;
      maxval = (y[0] >= 0.0) ? y[0] : -y[0] ;
      loc    = 0    ;
      for ( i = 1 ; i < size ; i++ ) {
         val = (y[i] >= 0.0) ? y[i] : -y[i] ;
         if ( maxval < val ) {
            maxval = val ; 
            loc    = i   ; 
         } 
      }
      *ploc = loc ;
   }
}
*ploc = loc ;

return(maxval) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of minimum value,
              *ploc contains its index 

   created -- 95sep22, cca
   ------------------------------------------------------
*/
double
DVmin ( 
   int      size, 
   double   y[], 
   int      *ploc 
) {
double   minval = 0.0 ;
int      loc = -1 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVmin, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
      exit(-1) ;
   } else {
      int      i ;
      minval = y[0] ;
      loc    = 0    ;
      for ( i = 1 ; i < size ; i++ ) {
         if ( minval > y[i] ) {
            minval = y[i] ; 
            loc    = i    ; 
         } 
      }
      *ploc = loc ;
   }
}
*ploc = loc ;

return(minval) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of minimum absolute value,
              *ploc contains its index 

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
double
DVminabs ( 
   int      size, 
   double   y[], 
   int      *ploc 
) {
double   minval = 0.0 ;
int      loc = -1 ;

if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVminabs, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
      exit(-1) ;
   } else {
      int      i   ;
      double   val ;
      minval = (y[0] >= 0.0) ? y[0] : -y[0] ;
      loc    = 0    ;
      for ( i = 1 ; i < size ; i++ ) {
         val = (y[i] >= 0.0) ? y[i] : -y[i] ;
         if ( minval > val ) {
            minval = val ; 
            loc    = i   ; 
         } 
      }
      *ploc = loc ;
   }
}
*ploc = loc ;

return(minval) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[*] := y[index[*]]

   created -- 95sep22, cca
   ------------------------------
*/
void
DVperm ( 
   int      size, 
   double   y[], 
   int      index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVperm, invalid data"
              "\n size = %d, y = %p, index = %p\n", size, y, index) ;
      exit(-1) ;
   } else {
      double   *x ;
      int      i ;
      x = DVinit2(size) ;
      DVcopy(size, x, y) ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = x[index[i]] ; 
      }
      DVfree(x) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to fill a double vector with a ramp function

   created -- 95sep22, cca
   -------------------------------------------------------
*/
void
DVramp ( 
   int      size, 
   double   y[], 
   double   start, 
   double   inc 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVramp, invalid input"
              "\n size = %d, y = %p, start = %f, inc = %f\n",
              size, y, start, inc) ;
      exit(-1) ;
   } else {
      int      i ;
      double   val ;
      for ( i = 0, val = start ; i < size ; i++, val += inc ) {
         y[i] = val ; 
      }
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to compute y[*] := y[*] - x[*]

   created -- 95sep22, cca
   -----------------------------------------
*/
void
DVsub ( 
   int      size, 
   double   y[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in DVsub, invalid input"
              "\n size = %d, y = %p, x = %p", size, y, x) ;
      exit(-1) ;
   } else {
      int i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] -= x[i] ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to scale a double vector by alpha

   created -- 95sep22, cca
   --------------------------------------------
*/
void
DVscale ( 
   int      size, 
   double   y[], 
   double   alpha 
) {
if ( size > 0 && alpha != 1.0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVscale, invalid data"
              "\n size = %d, y = %p, alpha = %f\n",
              size, y, alpha) ;
      exit(-1) ;
   } else {
      int i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] *= alpha ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to scatter y[index[*]] = x[*]

   created -- 95sep22, cca
   ----------------------------------------
*/
void
DVscatter ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVscatter, invalid data"
              "\n size = %d, y = %p, index = %p, x = %p\n",
              size, y, index, x) ;
      exit(-1) ; 
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] = x[i] ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to scatter add y[index[*]] += x[*]

   created -- 96aug17, cca
   ---------------------------------------------
*/
void
DVscatterAdd ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVscatterAdd, invalid data"
              "\n size = %d, y = %p, index = %p, x = %p\n",
              size, y, index, x) ;
      exit(-1) ; 
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] += x[i] ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to scatter add y[index[*]] += x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
DVscatterAddZero ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVscatterAddZero, invalid data"
              "\n size = %d, y = %p, index = %p, x = %p\n",
              size, y, index, x) ;
      exit(-1) ; 
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] += x[i] ; 
         x[i] = 0.0 ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to scatter y[index[*]] = x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------
*/
void
DVscatterZero ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in DVscatterZero, invalid data"
              "\n size = %d, y = %p, index = %p, x = %p\n",
              size, y, index, x) ;
      exit(-1) ; 
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] = x[i] ; 
         x[i] = 0.0 ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to return the sum of a double vector

   created -- 95sep22, cca
   -----------------------------------------------
*/
double
DVsum ( 
   int      size, 
   double   y[] 
) {
double   sum = 0.0 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVsum, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int      i ;
      for ( i = 0, sum = 0. ; i < size ; i++ ) {
         sum += y[i] ; 
      }
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the sum of the absolute values 
              of the entries in a double vector

   created -- 95sep22, cca
   ---------------------------------------------------
*/
double
DVsumabs ( 
   int      size, 
   double   y[] 
) {
double   sum = 0.0 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVsumabs, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int      i ;
      for ( i = 0, sum = 0. ; i < size ; i++ ) {
         sum += ((y[i] >= 0.0) ? y[i] : -y[i]) ; 
      }
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- to swap y[*] and x[*]

   created -- 95sep22, cca
   --------------------------------
*/
void
DVswap ( 
   int      size, 
   double   y[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in DVswap, invalid data"
              "\n size = %d, y = %p, x = %p", size, y, x) ;
      exit(-1) ; 
   } else {
      double   temp ;
      int      i ;
      for ( i = 0 ; i < size ; i++ ) {
         temp = x[i] ;
         x[i] = y[i] ;
         y[i] = temp ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   purpose -- to zero a double vector

   created -- 95sep22, cca
   ----------------------------------
*/
void
DVzero ( 
   int      size, 
   double   y[] 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVzero, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = 0. ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to permute an integer vector,
              procedure uses the Drand object

   input --

      size -- size of the vector
      y    -- vector to be permuted
      seed -- seed for the random number generator
              if seed <= 0, simple return

   created -- 95sep22, cca
   -------------------------------------------------
*/
void
DVshuffle ( 
   int      size, 
   double   y[], 
   int      seed 
) {
if ( size > 0 || seed <= 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in DVshuffle, invalid data"
              "\n size = %d, y = %p, seed = %d\n", size, y, seed) ;
      exit(-1) ; 
   } else {
      Drand    drand ;
      double   temp ;
      int      i, j ;

      Drand_setDefaultFields(&drand) ;
      Drand_setSeed(&drand, seed) ;
      for ( i = 0 ; i < size ; i++ ) {
         j = (int) (size * Drand_value(&drand)) ;
         temp = y[i] ;
         y[i] = y[j] ;
         y[j] = temp ;
      }
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to scale a double vector by a 2x2 matrix
 
   [ x ] := [ a b ] [ x ]
   [ y ]    [ c d ] [ y ]
 
   created -- 98jan23, cca
   ---------------------------------------------------
*/
void
DVscale2 ( 
   int      size, 
   double   x[], 
   double   y[], 
   double   a,
   double   b,
   double   c,
   double   d
) {
double   xi, yi ;
int      ii ;
 
if ( size < 0 || x == NULL || y == NULL ) {
   fprintf(stderr, "\n fatal error in DVscale2(%d,%p,%p,...)"
           "\n bad input\n", size, x, y) ;
   exit(-1) ;
} 
for ( ii = 0 ; ii < size ; ii++ ) {
   xi    = x[ii] ;
   yi    = y[ii] ;
   x[ii] = a*xi + b*yi ;
   y[ii] = c*xi + d*yi ;
}
return ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to perform a axpy with two vectors
 
   z := z + a*x + b*y
 
   where y and x are double vectors
 
   created -- 98jan23, cca
   --------------------------------------------
*/
void
DVaxpy2 (
   int      size,
   double   z[],
   double   a,
   double   x[],
   double   b,
   double   y[]
) {
double   xi, yi ;
int      ii ;
 
if ( size < 0 || y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in DVaxpy2()"
           "\n bad input\n") ;
   exit(-1) ;
}
for ( ii = 0 ; ii < size ; ii++ ) {
   xi = x[ii] ;
   yi = y[ii] ;
   z[ii] += a*xi + b*yi ;
}
return ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]
      sums[5] = row1[*] * col2[*]
      sums[6] = row2[*] * col0[*]
      sums[7] = row2[*] * col1[*]
      sums[8] = row2[*] * col2[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot33 (
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    col1[],
   double    col2[],
   double    sums[]
) {
register double   s00, s01, s02, s10, s11, s12, s20, s21, s22 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL || row1 == NULL || row2 == NULL
   || col0 == NULL || col1 == NULL || col2 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot33(%d,%p,%p,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, row1, row2, col0, col1, col2, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 9 dot products
   --------------------------
*/
s00 = s01 = s02 = s10 = s11 = s12 = s20 = s21 = s22 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, r1i, r2i, c0i, c1i, c2i ;

   r0i = row0[i]  ; r1i = row1[i]  ; r2i = row2[i] ;
   c0i = col0[i]  ; c1i = col1[i]  ; c2i = col2[i] ;
   s00 += r0i*c0i ; s01 += r0i*c1i ; s02 += r0i*c2i ;
   s10 += r1i*c0i ; s11 += r1i*c1i ; s12 += r1i*c2i ;
   s20 += r2i*c0i ; s21 += r2i*c1i ; s22 += r2i*c2i ;
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; sums[1] = s01 ; sums[2] = s02 ;
sums[3] = s10 ; sums[4] = s11 ; sums[5] = s12 ;
sums[6] = s20 ; sums[7] = s21 ; sums[8] = s22 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]
      sums[5] = row1[*] * col2[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot23 (
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    col1[],
   double    col2[],
   double    sums[]
) {
register double   s00, s01, s02, s10, s11, s12 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL || row1 == NULL 
   || col0 == NULL || col1 == NULL || col2 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot23(%d,%p,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, row1, col0, col1, col2, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 6 dot products
   --------------------------
*/
s00 = s01 = s02 = s10 = s11 = s12 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, r1i, c0i, c1i, c2i ;

   r0i = row0[i]  ; r1i = row1[i]  ; 
   c0i = col0[i]  ; c1i = col1[i]  ; c2i = col2[i] ;
   s00 += r0i*c0i ; s01 += r0i*c1i ; s02 += r0i*c2i ;
   s10 += r1i*c0i ; s11 += r1i*c1i ; s12 += r1i*c2i ;
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; sums[1] = s01 ; sums[2] = s02 ;
sums[3] = s10 ; sums[4] = s11 ; sums[5] = s12 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot13 (
   int       n,
   double    row0[],
   double    col0[],
   double    col1[],
   double    col2[],
   double    sums[]
) {
register double   s00, s01, s02 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL 
   || col0 == NULL || col1 == NULL || col2 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot13(%d,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, col0, col1, col2, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 3 dot products
   --------------------------
*/
s00 = s01 = s02 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, c0i, c1i, c2i ;

   r0i = row0[i]  ; 
   c0i = col0[i]  ; c1i = col1[i]  ; c2i = col2[i] ;
   s00 += r0i*c0i ; s01 += r0i*c1i ; s02 += r0i*c2i ;
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; sums[1] = s01 ; sums[2] = s02 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row1[*] * col0[*]
      sums[3] = row1[*] * col1[*]
      sums[4] = row2[*] * col0[*]
      sums[5] = row2[*] * col1[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot32 (
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    col1[],
   double    sums[]
) {
register double   s00, s01, s10, s11, s20, s21 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL || row1 == NULL || row2 == NULL
   || col0 == NULL || col1 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot32(%d,%p,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, row1, row2, col0, col1, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 6 dot products
   --------------------------
*/
s00 = s01 = s10 = s11 = s20 = s21 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, r1i, r2i, c0i, c1i ;

   r0i = row0[i]  ; r1i = row1[i]  ; r2i = row2[i] ;
   c0i = col0[i]  ; c1i = col1[i]  ; 
   s00 += r0i*c0i ; s01 += r0i*c1i ; 
   s10 += r1i*c0i ; s11 += r1i*c1i ; 
   s20 += r2i*c0i ; s21 += r2i*c1i ; 
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; sums[1] = s01 ; 
sums[2] = s10 ; sums[3] = s11 ; 
sums[4] = s20 ; sums[5] = s21 ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row1[*] * col0[*]
      sums[3] = row1[*] * col1[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot22 (
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    col1[],
   double    sums[]
) {
register double   s00, s01, s10, s11 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL || row1 == NULL 
   || col0 == NULL || col1 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot22(%d,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, row1, col0, col1, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 4 dot products
   --------------------------
*/
s00 = s01 = s10 = s11 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, r1i, c0i, c1i ;

   r0i = row0[i]  ; r1i = row1[i]  ; 
   c0i = col0[i]  ; c1i = col1[i]  ; 
   s00 += r0i*c0i ; s01 += r0i*c1i ; 
   s10 += r1i*c0i ; s11 += r1i*c1i ; 
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; sums[1] = s01 ; 
sums[2] = s10 ; sums[3] = s11 ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot12 (
   int       n,
   double    row0[],
   double    col0[],
   double    col1[],
   double    sums[]
) {
register double   s00, s01 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL 
   || col0 == NULL || col1 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot12(%d,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, col0, col1, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 2 dot products
   --------------------------
*/
s00 = s01 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, c0i, c1i ;

   r0i = row0[i]  ; 
   c0i = col0[i]  ; c1i = col1[i]  ; 
   s00 += r0i*c0i ; s01 += r0i*c1i ; 
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; sums[1] = s01 ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row1[*] * col0[*]
      sums[2] = row2[*] * col0[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot31 (
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    sums[]
) {
register double   s00, s10, s20 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL || row1 == NULL || row2 == NULL
   || col0 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot31(%d,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, row1, row2, col0, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 3 dot products
   --------------------------
*/
s00 = s10 = s20 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, r1i, r2i, c0i ;

   r0i = row0[i]  ; r1i = row1[i]  ; r2i = row2[i] ;
   c0i = col0[i]  ; 
   s00 += r0i*c0i ; 
   s10 += r1i*c0i ; 
   s20 += r2i*c0i ; 
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; 
sums[1] = s10 ; 
sums[2] = s20 ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[2] = row1[*] * col0[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot21 (
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    sums[]
) {
register double   s00, s10 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL || row1 == NULL 
   || col0 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot21(%d,%p,%p,%p,%p)"
           "\n bad input\n",
           n, row0, row1, col0, sums) ;
   exit(-1) ;
}
/*
   --------------------------
   compute the 2 dot products
   --------------------------
*/
s00 = s10 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, r1i, c0i ;

   r0i = row0[i]  ; r1i = row1[i]  ; 
   c0i = col0[i]  ; 
   s00 += r0i*c0i ; 
   s10 += r1i*c0i ; 
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; 
sums[1] = s10 ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a single dot product

      sums[0] = row0[*] * col0[*]

   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot11 (
   int       n,
   double    row0[],
   double    col0[],
   double    sums[]
) {
register double   s00 ;
register int      i ;
/*
   ---------------
   check the input
   ---------------
*/
if (  sums == NULL 
   || row0 == NULL 
   || col0 == NULL ) {
   fprintf(stderr, "\n fatal error in DVdot11(%d,%p,%p,%p)"
           "\n bad input\n",
           n, row0, col0, sums) ;
   exit(-1) ;
}
/*
   -------------------------
   compute the 1 dot product
   -------------------------
*/
s00 = 0.0 ;
for ( i = 0 ; i < n ; i++ ) {
   register double r0i, c0i ;

   r0i = row0[i]  ; 
   c0i = col0[i]  ; 
   s00 += r0i*c0i ;
}
/*
   ----------------------
   store the dot products
   ----------------------
*/
sums[0] = s00 ; 

return ; }

/*--------------------------------------------------------------------*/
