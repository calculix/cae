/*  IV.c  */

#include "../Utilities.h"
#include "../../Drand.h"

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
IVcompress ( 
   int   size1, 
   int   x1[], 
   int   y1[],
   int   size2,  
   int   x2[],  
   int   y2[] 
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
   fprintf(stderr, "\n fatal error in IVcompress, invalid data"
           "\n size1 = %d, x1 = %p, y1 = %p"
           "\n size2 = %d, x2 = %p, y2 = %p\n",
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
   ds[j-1] = sqrt((double) (dx*dx + dy*dy)) ;
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
   -------------------------------
   purpose -- to copy y[*] := x[*]

   created -- 95sep22, cca
   -------------------------------
*/
void
IVcopy ( 
   int   size, 
   int   y[], 
   int   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in IVcopy, invalid data"
              "\n size = %d, y = %p, x = %p", size, y, x) ;
      exit(-1) ;
   } else {
      int i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = x[i] ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to fill a int vector with a value

   created -- 95sep22, cca
   -----------------------------------------------
*/
void
IVfill ( 
   int   size, 
   int   y[], 
   int   ival 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVfill, invalid data"
              "\n size = %d, y = %p, ival = %d\n",
              size, y, ival) ;
      exit(-1) ;
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = ival ; 
      }
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- to print out a int vector

   created -- 95sep22, cca
   ------------------------------------
*/
void
IVfprintf ( 
   FILE   *fp, 
   int    size, 
   int    y[] 
) {
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVfprintf, invalid data"
              "\n fp = %p, size = %d, y = %p\n",
              fp, size, y) ;
      exit(-1) ;
   } else {
      int    i ;
      for ( i = 0 ; i < size ; i++ ) {
         if ( i % 16 == 0 ) fprintf(fp, "\n") ;
         fprintf(fp, " %4d", y[i]) ; 
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- to write out an integer vector with eighty column lines

   input --

      fp     -- file pointer, must be formatted and write access
      size   -- length of the vector
      y      -- integer vector
      column -- present column
      pierr  -- pointer to int to hold return value, 
                should be 1 if any print was successful,
                if fprintf() failed, then ierr = -1
  
   return value -- present column

   created -- 95sep22, cca
   mods    -- 95sep29, cca
      added ierr argument to handle error returns from fprintf()
   -------------------------------------------------------------------
*/
int
IVfp80 ( 
   FILE   *fp, 
   int    size, 
   int    y[], 
   int    column,
   int    *pierr
) {
*pierr = 1 ;
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVfp80, invalid input"
              "\n fp = %p, size = %d, y = %p, column = %d\n",
              fp, size, y, column) ;
      exit(-1) ;
   } else {
      int    i, inum, nchar ;
      for ( i = 0 ; i < size ; i++ ) {
         inum = y[i] ;
         if ( inum < 0 ) {
            inum = -inum ; 
            nchar = 2 ; 
         } else if ( inum == 0 ) {
            nchar = 2 ; 
         } else {
            nchar = 1 ; 
         }
         while ( inum > 0 ) {
            nchar++ ;
            inum /= 10 ;
         }
         if ( (column += nchar) >= 80 ) {
            fprintf(fp, "\n") ;
            column = nchar ; 
         }
         if ( (*pierr = fprintf(fp, " %d", y[i])) < 0 ) {
/*
           --------------------
           error with fprintf()
           --------------------
*/
           break ;
         }
      }
   }
}
return(column) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- to free the storage taken by an integer vector.
              note : y must have been returned by IVinit

   created -- 95sep22, cca
   ----------------------------------------------------------
*/
void
IVfree ( 
   int y[] 
) {
if ( y != NULL ) {
   FREE(y) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to read in an int vector

   created -- 95sep22, cca
   -----------------------------------
*/
int
IVfscanf ( 
   FILE   *fp, 
   int    size, 
   int    y[] 
) {
int    i = 0 ;
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVfscanf, invalid data"
              "\n fp = %p, size = %d, y = %p\n", fp, size, y) ;
      exit(-1) ;
   } else {
      for ( i = 0 ; i < size ; i++ ) {
         if ( fscanf(fp, " %d", y + i) != 1 ) {
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
IVgather ( 
   int   size, 
   int   y[], 
   int   x[], 
   int   index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in IVgather, invalid data"
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
   ---------------------------------------------------
   purpose -- allocate a int array with size entries
              and fill with value ival
   
   return value -- a pointer to the start of the array
 
   created : 95sep22, cca
   ---------------------------------------------------
*/
int *
IVinit ( 
   int   size, 
   int   ival 
) {
int   *y = NULL ;
if ( size > 0 ) {
   y = IVinit2(size) ;
   IVfill(size, y, ival) ;
}
return(y) ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- allocate a int array with size entries
   
   return value -- a pointer to the start of the array
 
   created : 95sep22, cca
   ---------------------------------------------------
*/
int *
IVinit2 ( 
   int   size 
) {
int   *y = NULL ;
if ( size > 0 ) {
   ALLOCATE(y, int, size) ;
}
return(y) ; }
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- allocate an int array and fill with the inverse of
              y[]. note, y[] must be a permutation vector.

   created : 95sep22, cca
   -------------------------------------------------------------
*/
int *
IVinverse ( 
   int   size, 
   int   y[] 
) {
int   *x = NULL ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVinverse, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int   i, j ;
      x = IVinit(size, -1) ;
      for ( i = 0 ; i < size ; i++ ) {
         j = y[i] ;
         if ( j < 0 || j >= size || x[j] != -1 ) {
            fprintf(stderr, 
                    "\n fatal error in IVinverse"
                    "\n y[%d] = %d, value out-of-range or repeated",
                    i, j) ;
            exit(-1) ;
         }
         x[j] = i ; 
      }
   }
}
return(x) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[index[*]] := y[*]

   created -- 95sep22, cca
   ------------------------------
*/
void
IVinvPerm ( 
   int   size, 
   int   y[], 
   int   index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in IVinvPerm, invalid data"
              "\n size = %d, y = %p, index = %p\n", size, y, index) ;
      exit(-1) ;
   } else {
      int   *x ;
      int    i ;
      x = IVinit2(size) ;
      IVcopy(size, x, y) ;
      for ( i = 0 ; i < size ; i++ ) {
         y[index[i]] = x[i] ;
      }
      IVfree(x) ;
   }
}
 
return ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to return the first entry of maximum size,
              *ploc contains its index 

   created -- 95sep22, cca
   -----------------------------------------------------
*/
int
IVmax ( 
   int   size, 
   int   y[], 
   int   *ploc 
) {
int   maxval = 0  ;
int   loc    = -1 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVmax, invalid data"
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
int
IVmaxabs (
   int   size,
   int   y[],
   int   *ploc
) {
int   maxval =  0 ;
int   loc    = -1 ;
 
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVmaxabs, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
     exit(-1) ;
   } else {
      int   i, val ;
      maxval = (y[0] >= 0) ? y[0] : -y[0] ;
      loc    = 0    ;
      for ( i = 1 ; i < size ; i++ ) {
         val = (y[i] >= 0) ? y[i] : -y[i] ;
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
int
IVmin (
   int   size,
   int   y[],
   int   *ploc
) {
int   minval =  0 ;
int   loc    = -1 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVmin, invalid data"
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
int
IVminabs (
   int   size,
   int   y[],
   int   *ploc
) {
int   minval =  0 ;
int   loc    = -1 ;
 
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVminabs, invalid data"
              "\n size = %d, y = %p, ploc = %p\n", size, y, ploc) ;
     exit(-1) ;
   } else {
      int   i, val ;
      minval = (y[0] >= 0) ? y[0] : -y[0] ;
      loc    = 0    ;
      for ( i = 1 ; i < size ; i++ ) {
         val = (y[i] >= 0) ? y[i] : -y[i] ;
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
IVperm ( 
   int   size, 
   int   y[], 
   int   index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in IVperm, invalid data"
              "\n size = %d, y = %p, index = %p\n",
              size, y, index) ;
      exit(-1) ;
   } else {
      int   *x ;
      int   i  ;
      x = IVinit2(size) ;
      IVcopy(size, x, y) ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = x[index[i]] ; 
      }
      IVfree(x) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to fill a int vector with a ramp function

   created -- 95sep22, cca
   ----------------------------------------------------
*/
void
IVramp ( 
   int   size, 
   int   y[], 
   int   start, 
   int   inc 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVramp, invalid data"
              "\n size = %d, y = %p, start = %d, inc = %d\n",
              size, y, start, inc) ;
      exit(-1) ;
   } else {
      int   i, j ;
      for ( i = 0, j = start ; i < size ; i++, j += inc ) {
         y[i] = j ; 
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
IVscatter ( 
   int   size, 
   int   y[], 
   int   index[], 
   int   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in IVscatter, invalid data"
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
   -----------------------------------------------
   purpose -- to return the sum of a int vector

   created -- 95sep22, cca
   -----------------------------------------------
*/
int
IVsum ( 
   int   size, 
   int   y[] 
) {
int   sum = 0 ;

if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVsum, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int   i ;
      for ( i = 0, sum = 0 ; i < size ; i++ ) {
         sum += y[i] ; 
      }
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the sum of the absolute values 
              of the entries in an int vector

   created -- 95sep22, cca
   ---------------------------------------------------
*/
int
IVsumabs ( 
   int   size, 
   int   y[] 
) {
int   sum = 0 ;
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVsumabs, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int   i, sum ;
      for ( i = 0, sum = 0 ; i < size ; i++ ) {
         sum += ((y[i] >= 0) ? y[i] : -y[i]) ; 
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
IVswap ( 
   int   size, 
   int   y[], 
   int   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL ) {
      fprintf(stderr, "\n fatal error in IVswap, invalid data"
              "\n size = %d, y = %p, x = %p\n", size, y, x) ;
      exit(-1) ;
   } else {
      int   i, temp ;
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
   purpose -- to zero a int vector

   created -- 95sep22, cca
   ----------------------------------
*/
void
IVzero ( 
   int   size, 
   int   y[] 
) {
if ( size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVzero, invalid data"
              "\n size = %d, y = %p\n", size, y) ;
      exit(-1) ;
   } else {
      int   i ;
      for ( i = 0 ; i < size ; i++ ) {
         y[i] = 0 ; 
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
IVshuffle ( 
   int   size, 
   int   y[], 
   int   seed 
) {
if ( size > 0 && seed > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in IVshuffle, invalid data"
              "\n size = %d, y = %p, seed = %d\n", size, y, seed) ;
     exit(-1) ;
   } else {
      int      i, j, temp ;
      double   value ;
      Drand    drand ;

      Drand_setDefaultFields(&drand) ;
      Drand_setSeed(&drand, seed) ;
      Drand_setUniform(&drand, 0.0, 1.0) ;
      for ( i = 0 ; i < size ; i++ ) {
         value = Drand_value(&drand) ;
         j = (int) (size * value) ;
         temp = y[i] ;
         y[i] = y[j] ;
         y[j] = temp ;
      }
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   locate an instance of target in the vector ivec[size].
   we assume that ivec[] is sorted in ascending order
   so we can use binary search to locate target.

   return value --
      -1  -- if target not found in ivec[]
      loc -- where target = ivec[loc]

   created -- 96may27, cca
   ------------------------------------------------------
*/
int
IVlocateViaBinarySearch (
   int   size,
   int   ivec[],
   int   target
) {
int   left, mid, right ;

if ( size <= 0 || ivec == NULL ) {
   return(-1) ;
}
left  = 0 ;
right = size - 1 ;
if ( target < ivec[left] || ivec[right] < target ) {
   return(-1) ;
} else if ( target == ivec[left] ) {
   return(left) ;
} else if ( target == ivec[right] ) {
   return(right) ;
} else {
   while ( right > left + 1 ) {
      mid = (left + right)/2 ;
      if ( target < ivec[mid] ) {
         right = mid ;
      } else if ( target > ivec[mid] ) {
         left = mid ;
      } else {
         return(mid) ;
      }
   }
}
return(-1) ; }

/*--------------------------------------------------------------------*/
