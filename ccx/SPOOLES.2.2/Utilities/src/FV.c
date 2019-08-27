/*  FV.c  */

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
FVadd ( 
   int     size, 
   float   y[], 
   float   x[] 
) {
if ( size <= 0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in FVadd"
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
FVaxpy ( 
   int     size, 
   float   y[], 
   float   alpha, 
   float   x[] 
) {
if ( size < 0 || alpha == 0.0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in FVaxpy"
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
FVaxpyi ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   alpha, 
   float   x[] 
) {
if ( size <= 0 || alpha == 0.0 ) {
   return ;
} else if ( y == NULL || index == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in FVaxpyi, invalid input"
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
FVcopy ( 
   int     size, 
   float   y[], 
   float   x[] 
) {
if ( size <= 0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in FVcopy, invalid input"
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
FVcompress ( 
   int     size1, 
   float   x1[], 
   float   y1[],
   int     size2, 
   float   x2[], 
   float   y2[] 
) {
float   delta, dx, dy, path, totalPath ;
float   *ds ;
int     i, j ;
/*
   --------------------
   check the input data
   --------------------
*/
if ( size1 <= 0 || size2 <= 0 ) {
   return(0) ;
} else if ( x1 == NULL || y1 == NULL || x2 == NULL || y2 == NULL ) {
   fprintf(stderr, "\n fatal error in FVcompress, invalid data"
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
ds = FVinit(size1, 0.0) ;
for ( j = 1 ; j < size1 ; j++ ) {
   dx = x1[j] - x1[j-1] ;
   dy = y1[j] - y1[j-1] ;
   ds[j-1] = sqrt(dx*dx + dy*dy) ;
}
totalPath = FVsum(size1, ds) ;
delta = totalPath / (size2-2) ;
#if MYDEBUG > 0
fprintf(stdout, "\n totalPath = %12.4e, delta = %12.4e, ds",
        totalPath, delta) ;
FVfprintf(stdout, size1, ds) ;
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
FVfree(ds) ;

return(i) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to copy sum{y[*] * x[*]}

   created -- 95sep22, cca
   -----------------------------------
*/
float
FVdot ( 
   int     size, 
   float   y[], 
   float   x[] 
) {
float   sum = 0.0 ;
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in FVdot, invalid data"
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
   -----------------------------------------------
   purpose -- to fill a float vector with a value

   created -- 95sep22, cca
   -----------------------------------------------
*/
void
FVfill ( 
   int     size, 
   float   y[], 
   float   dval 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVfill, invalid data"
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
   purpose -- to print out a float vector

   created -- 95sep22, cca
   -----------------------------------------
*/
void
FVfprintf ( 
   FILE    *fp, 
   int     size, 
   float   y[]
) {
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVfprintf, invalid input"
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
   purpose -- to free storage taken by a float vector.
              note : y[] must have been returned by FVinit.

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
FVfree ( 
   float   y[] 
) {
if ( y != NULL ) {
   FREE(y) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to read in a float vector
              return value -- # of entries read

   created -- 95sep22, cca
   --------------------------------------------
*/
int
FVfscanf ( 
   FILE    *fp, 
   int     size, 
   float   y[] 
) {
int    i  = 0, rc ;
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVfscanf, invalid input"
              "\n fp = %p, size = %d, y = %p\n", fp, size, y) ;
      exit(-1) ;
   } else {
      for ( i = 0 ; i < size ; i++ ) {
         if ( (rc = fscanf(fp, " %f", y + i)) != 1 ) {
            fprintf(stderr, 
                    "\n fatal error in FVfscanf(%p,%d,%p), rc = %d", 
                    fp, size, y, rc) ;
            break ; 
         } 
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
FVgather ( 
   int     size, 
   float   y[], 
   float   x[], 
   int     index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVgather, invalid input"
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
FVgatherAddZero ( 
   int     size, 
   float   y[], 
   float   x[], 
   int     index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVgatherAddZero, invalid input"
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
FVgatherZero ( 
   int     size, 
   float   y[], 
   float   x[], 
   int     index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVgatherZero, invalid input"
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
   purpose -- allocate a float array with size entries
              and fill with value dval
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
float *
FVinit ( 
   int     size, 
   float   dval 
) {
float   *y = NULL ;
if ( size > 0 ) {
   y = FVinit2(size) ;
   FVfill(size, y, dval) ;
}
return(y) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a float array with size entries
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
float *
FVinit2 ( 
   int   size 
) {
float   *y = NULL ;
if ( size > 0 ) {
   ALLOCATE(y, float, size) ;
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
FVinvPerm ( 
   int     size, 
   float   y[], 
   int     index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVinvPerm, invalid data"
              "\n size = %d, y = %p, index = %p", size, y, index) ;
      exit(-1) ;
   } else {
      float   *x ;
      int      i ;
      x = FVinit2(size) ;
      FVcopy(size, x, y) ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] = x[i] ; 
      }
      FVfree(x) ;
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
float
FVmax ( 
   int     size, 
   float   y[], 
   int     *ploc 
) {
float   maxval = 0.0 ;
int      loc = -1 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVmax, invalid data"
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
float
FVmaxabs ( 
   int     size, 
   float   y[], 
   int     *ploc 
) {
float   maxval = 0.0 ;
int      loc = -1 ;

if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVmaxabs, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
      exit(-1) ;
   } else {
      int      i   ;
      float   val ;
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
float
FVmin ( 
   int     size, 
   float   y[], 
   int     *ploc 
) {
float   minval = 0.0 ;
int      loc = -1 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVmin, invalid data"
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
float
FVminabs ( 
   int     size, 
   float   y[], 
   int     *ploc 
) {
float   minval = 0.0 ;
int      loc = -1 ;

if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVminabs, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
      exit(-1) ;
   } else {
      int      i   ;
      float   val ;
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
FVperm ( 
   int     size, 
   float   y[], 
   int     index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVperm, invalid data"
              "\n size = %d, y = %p, index = %p\n", size, y, index) ;
      exit(-1) ;
   } else {
      float   *x ;
      int      i ;
      x = FVinit2(size) ;
      FVcopy(size, x, y) ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = x[index[i]] ; 
      }
      FVfree(x) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to fill a float vector with a ramp function

   created -- 95sep22, cca
   -------------------------------------------------------
*/
void
FVramp ( 
   int     size, 
   float   y[], 
   float   start, 
   float   inc 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVramp, invalid input"
              "\n size = %d, y = %p, start = %f, inc = %f\n",
              size, y, start, inc) ;
      exit(-1) ;
   } else {
      int      i ;
      float   val ;
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
FVsub ( 
   int     size, 
   float   y[], 
   float   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in FVsub, invalid input"
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
   purpose -- to scale a float vector by alpha

   created -- 95sep22, cca
   --------------------------------------------
*/
void
FVscale ( 
   int     size, 
   float   y[], 
   float   alpha 
) {
if ( size > 0 && alpha != 1.0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVscale, invalid data"
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
FVscatter ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVscatter, invalid data"
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
   -----------------------------------------------------------
   purpose -- to scatter add y[index[*]] += x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
FVscatterAddZero ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVscatterAddZero, invalid data"
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
FVscatterZero ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in FVscatterZero, invalid data"
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
   purpose -- to return the sum of a float vector

   created -- 95sep22, cca
   -----------------------------------------------
*/
float
FVsum ( 
   int     size, 
   float   y[] 
) {
float   sum = 0.0 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVsum, invalid data"
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
              of the entries in a float vector

   created -- 95sep22, cca
   ---------------------------------------------------
*/
float
FVsumabs ( 
   int     size, 
   float   y[] 
) {
float   sum = 0.0 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVsumabs, invalid data"
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
FVswap ( 
   int     size, 
   float   y[], 
   float   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in FVswap, invalid data"
              "\n size = %d, y = %p, x = %p", size, y, x) ;
      exit(-1) ; 
   } else {
      float   temp ;
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
   purpose -- to zero a float vector

   created -- 95sep22, cca
   ----------------------------------
*/
void
FVzero ( 
   int     size, 
   float   y[] 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVzero, invalid data"
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
FVshuffle ( 
   int     size, 
   float   y[], 
   int     seed 
) {
if ( size > 0 || seed <= 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in FVshuffle, invalid data"
              "\n size = %d, y = %p, seed = %d\n", size, y, seed) ;
      exit(-1) ; 
   } else {
      float   temp ;
      int     i, j ;
      Drand   drand ;

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
