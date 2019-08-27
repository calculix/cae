/*  LocalND.c  */

#include "../misc.h"

#define MYDEBUG 0
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   compute an old-to-new ordering for 
   local nested dissection in two dimensions

   n1       -- number of grid points in first direction
   n2       -- number of grid points in second direction
   p1       -- number of domains in first direction
   p2       -- number of domains in second direction
   dsizes1  -- domain sizes in first direction, size p1
               if NULL, then we construct our own
   dsizes2  -- domain sizes in second direction, size p2
               if NULL, then we construct our own
   oldToNew -- old-to-new permutation vector

   note : the following must hold
      n1 > 0, n2 >0, n1 >= 2*p1 - 1, n2 >= 2*p2 - 1, p2 > 1
      sum(dsizes1) = n1 - p1 + 1 and sum(dsizes2) = n2 - p2 + 1

   created -- 95nov16, cca
   ------------------------------------------------------------
*/
void
localND2D ( 
   int   n1, 
   int   n2, 
   int   p1, 
   int   p2, 
   int   dsizes1[], 
   int   dsizes2[], 
   int   oldToNew[] 
) {
int   i, idom, ij, isw, j, jdom, jsw, length1, length2, 
      m, m1, m2, msize, now, nvtx ;
int   *length1s, *length2s, *isws, *jsws, *temp ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n1 <= 0 || n2 <= 0 || 2*p1 - 1 > n1 || 2*p2 - 1 > n2 
   || oldToNew == NULL ) {
   fprintf(stderr, "\n fatal error in localND2D(%d,%d,%d,%d,%p,%p,%p)"
           "\n bad input\n",
           n1, n2, p1, p2, dsizes1, dsizes2, oldToNew) ;
   exit(-1) ; 
}
if ( p2 <= 1 ) {
   fprintf(stderr, "\n fatal error in localND2D(%d,%d,%d,%d,%p,%p,%p)"
           "\n p2 = %d, must be > 1", 
           n1, n2, p1, p2, dsizes1, dsizes2, oldToNew, p2) ;
   exit(-1) ; 
}
if ( dsizes1 != NULL && IVsum(p1, dsizes1) != n1 - p1 + 1 ) {
   fprintf(stderr, "\n fatal error in localND2D(%d,%d,%d,%d,%p,%p,%p)"
           "\n IVsum(p1, dsizes1) = %d != %d = n1 - p1 + 1 ",
           n1, n2, p1, p2, dsizes1, dsizes2, oldToNew, 
           IVsum(p1, dsizes1), n1 - p1 + 1) ;
   return ; 
}
if ( dsizes1 != NULL && IVmin(p1, dsizes1, &i) <= 0 ) {
   fprintf(stderr, "\n fatal error in localND2D(%d,%d,%d,%d,%p,%p,%p)"
           "\n IVmin(p1, dsizes1) = %d must be > 0",
           n1, n2, p1, p2, dsizes1, dsizes2, oldToNew, 
           IVmin(p1, dsizes1, &i)) ;
   return ; 
}
if ( dsizes2 != NULL && IVsum(p2, dsizes2) != n2 - p2 + 1 ) {
   fprintf(stderr, "\n fatal error in localND2D(%d,%d,%d,%d,%p,%p,%p)"
           "\n IVsum(p2, dsizes2) = %d != %d = n2 - p2 + 1 ",
           n1, n2, p1, p2, dsizes1, dsizes2, oldToNew, 
           IVsum(p2, dsizes2), n2 - p2 + 1) ;
   return ; 
}
if ( dsizes2 != NULL && IVmin(p2, dsizes2, &i) <= 0 ) {
   fprintf(stderr, "\n fatal error in localND2D(%d,%d,%d,%d,%p,%p,%p)"
           "\n IVmin(p2, dsizes2) = %d must be > 0",
           n1, n2, p1, p2, dsizes1, dsizes2, oldToNew, 
           IVmin(p2, dsizes2, &i)) ;
   return ; 
}
nvtx = n1*n2 ;
/*
   ----------------------------------
   construct the domain sizes vectors
   ----------------------------------
*/
if ( dsizes1 == NULL ) {
   length1s = IVinit(p1, 0) ;
   length1  = (n1 - p1 + 1) / p1 ;
   m1       = (n1 - p1 + 1) % p1 ;
   for ( i = 0 ; i < m1 ; i++ ) {
      length1s[i] = length1 + 1 ;
   }
   for ( ; i < p1 ; i++ ) {
      length1s[i] = length1 ;
   }
} else {
   length1s = dsizes1 ;
}
if ( dsizes2 == NULL ) {
   length2s = IVinit(p2, 0) ;
   length2  = (n2 - p2 + 1) / p2 ;
   m2       = (n2 - p2 + 1) % p2 ;
   for ( i = 0 ; i < m2 ; i++ ) {
      length2s[i] = length2 + 1 ;
   }
   for ( ; i < p2 ; i++ ) {
      length2s[i] = length2 ;
   }
} else {
   length2s = dsizes2 ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n inside localND2D") ;
fprintf(stdout, "\n n1 = %d, n2 = %d, p1 = %d, p2 = %d", 
        n1, n2, p1, p2) ;
fprintf(stdout, "\n length1s[%d] = ", p1) ;
IVfp80(stdout, p1, length1s, 12) ;
fprintf(stdout, "\n length2s[%d] = ", p2) ;
IVfp80(stdout, p2, length2s, 12) ;
#endif
/*
   ---------------------------------------
   determine the first and last domain ids 
   and the array of southwest points
   ---------------------------------------
*/
isws = IVinit(p1, -1) ;
for ( idom = 0, isw = 0 ; idom < p1 ; idom++ ) {
   isws[idom] = isw ;
   isw += length1s[idom] + 1 ;
}
jsws = IVinit(p2, -1) ;
for ( jdom = 0, jsw = 0 ; jdom < p2 ; jdom++ ) {
   jsws[jdom] = jsw ;
   jsw += length2s[jdom] + 1 ;
}
#if MYDEBUG > 1
fprintf(stdout, "\n isws[%d] = ", p1) ;
IVfp80(stdout, p1, isws, 12) ;
fprintf(stdout, "\n jsws[%d] = ", p2) ;
IVfp80(stdout, p2, jsws, 12) ;
#endif
/*
   ----------------------------------------------------------------
   create a temporary permutation vector for the domains' orderings
   ----------------------------------------------------------------
*/
msize = IVmax(p1, length1s, &i) * IVmax(p2, length2s, &i) ;
temp  = IVinit(msize, -1) ;
/*
   ------------------------
   fill in the domain nodes
   ------------------------
*/
now = 0 ;
for ( jdom = 0; jdom < p2 ; jdom++ ) {
   jsw     = jsws[jdom] ;
   length2 = length2s[jdom] ;
   for ( idom = 0 ; idom < p1 ; idom++ ) {
      length1 = length1s[idom] ;
      isw     = isws[idom] ;
      mkNDperm(length1, length2, 1, temp, 0, length1-1, 
               0, length2-1, 0, 0) ;
      for ( m = 0 ; m < length1*length2 ; m++ ) {
         ij = temp[m] ;
         i  = isw + (ij % length1) ;
         j  = jsw + (ij / length1) ;
         ij = i + j*n1 ;
         oldToNew[ij] = now++ ;
      }
   }
}
#if MYDEBUG > 2
fprintf(stdout, "\n old-to-new after domains are numbered") ;
fp2DGrid(n1, n2, oldToNew, stdout) ;
#endif
/*
   ---------------------------------
   fill in the lower separator nodes
   ---------------------------------
*/
for ( jdom = 0 ; jdom < (p2/2) ; jdom++ ) {
   jsw     = jsws[jdom] ;
   length2 = length2s[jdom] ;
   for ( idom = 0 ; idom < p1 ; idom++ ) {
      isw     = isws[idom] ;
      length1 = length1s[idom] ;
      if ( isw > 0 ) {
         i = isw - 1 ;
         for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
            ij = i + j*n1 ;
            oldToNew[ij] = now++ ;
         }
      }
      if ( isw > 0 && jsw > 0 ) {
         i = isw - 1 ;
         j = jsw - 1 ;
         ij = i + j*n1 ;
         oldToNew[ij] = now++ ;
      }
      if ( jsw > 0 ) {
         j = jsw - 1 ;
         for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
            ij = i + j*n1 ;
            oldToNew[ij] = now++ ;
         }
      }
   }
}
#if MYDEBUG > 2
fprintf(stdout, "\n after the lower separators filled in") ;
fp2DGrid(n1, n2, oldToNew, stdout) ;
#endif
/*
   ---------------------------------
   fill in the upper separator nodes
   ---------------------------------
*/
for ( jdom = p2 - 1 ; jdom >= (p2/2) ; jdom-- ) {
   jsw     = jsws[jdom] ;
   length2 = length2s[jdom] ;
   for ( idom = p1 - 1 ; idom >= 0 ; idom-- ) {
      isw     = isws[idom] ;
      length1 = length1s[idom] ;
      if ( isw + length1 < n1 ) {
         i = isw + length1 ;
         for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
            ij = i + j*n1 ;
            oldToNew[ij] = now++ ;
         }
      }
      if ( isw + length1 < n1 && jsw + length2 < n2 ) {
         i = isw + length1 ;
         j = jsw + length2 ;
         ij = i + j*n1 ;
         oldToNew[ij] = now++ ;
      }
      if ( jsw + length2 < n2 ) {
         j = jsw + length2 ;
         for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
            ij = i + j*n1 ;
            oldToNew[ij] = now++ ;
         }
      }
   }
}
#if MYDEBUG > 2
fprintf(stdout, "\n after the upper separators filled in") ;
fp2DGrid(n1, n2, oldToNew, stdout) ;
#endif
/*
   -------------------------------
   fill in the top level separator
   -------------------------------
*/
m1 = p2 / 2 ;
for ( jdom = 0, j = 0 ; jdom < m1 ; jdom++ ) {
   j += length2s[jdom] + 1 ;
}
j-- ;
for ( i = 0 ; i < n1 ; i++ ) { 
   ij = i + j*n1 ;
   oldToNew[ij] = now++ ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( dsizes1 == NULL ) {
   IVfree(length1s) ;
}
if ( dsizes2 == NULL ) {
   IVfree(length2s) ;
}
IVfree(isws) ;
IVfree(jsws) ;
IVfree(temp) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   compute an old-to-new ordering for 
   local nested dissection in three dimensions

   n1       -- number of grid points in first direction
   n2       -- number of grid points in second direction
   n3       -- number of grid points in third direction
   p1       -- number of domains in first direction
   p2       -- number of domains in second direction
   p3       -- number of domains in third direction
   dsizes1  -- domain sizes in first direction, size p1
               if NULL, then we construct our own
   dsizes2  -- domain sizes in second direction, size p2
               if NULL, then we construct our own
   dsizes3  -- domain sizes in third direction, size p3
               if NULL, then we construct our own
   oldToNew -- old-to-new permutation vector

   note : the following must hold
      n1 > 0, n2 >0, n3 > 0,
      n1 >= 2*p1 - 1, n2 >= 2*p2 - 1, n3 >= 2*p3 - 1, p3 > 1
      sum(dsizes1) = n1 - p1 + 1, sum(dsizes2) = n2 - p2 + 1
      sum(dsizes3) = n3 - p3 + 1

   created -- 95nov16, cca
   ------------------------------------------------------------
*/
void
localND3D ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   p1, 
   int   p2, 
   int   p3,
   int   dsizes1[], 
   int   dsizes2[], 
   int   dsizes3[], 
   int   oldToNew[] 
) {
int   i, idom, ijk, isw, j, jdom, jsw, k, kdom, ksw, 
      length1, length2, length3, m, m1, m2, m3, msize, now, nvtx ;
int   *length1s, *length2s, *length3s, *isws, *jsws, *ksws, *temp ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n1 <= 0 || n2 <= 0 || n3 <= 0
   || 2*p1 - 1 > n1 || 2*p2 - 1 > n2 || 2*p3 - 1 > n3 ) {
   fprintf(stderr, "\n error in input data") ;
   return ; 
}
if ( p3 <= 1 ) {
   fprintf(stderr, "\n p3 must be > 1") ;
   return ; 
}
if ( oldToNew == NULL ) {
   fprintf(stderr, "\n oldToNew = NULL") ;
   return ; 
}
if ( dsizes1 != NULL && IVsum(p1, dsizes1) != n1 - p1 + 1 ) {
   fprintf(stderr, "\n IVsum(p1, dsizes1) != n1 - p1 + 1 ") ;
   return ; 
}
if ( dsizes2 != NULL && IVsum(p2, dsizes2) != n2 - p2 + 1 ) {
   fprintf(stderr, "\n IVsum(p2, dsizes2) != n2 - p2 + 1 ") ;
   return ; 
}
if ( dsizes3 != NULL && IVsum(p3, dsizes3) != n3 - p3 + 1 ) {
   fprintf(stderr, "\n IVsum(p3, dsizes3) != n3 - p3 + 1 ") ;
   return ; 
}
if ( dsizes1 != NULL && IVmin(p1, dsizes1, &i) <= 0 ) {
   fprintf(stderr, "\n IVmin(p1, dsizes1) must be > 0") ;
   return ; 
}
if ( dsizes2 != NULL && IVmin(p2, dsizes2, &i) <= 0 ) {
   fprintf(stderr, "\n IVmin(p2, dsizes2) must be > 0") ;
   return ; 
}
if ( dsizes3 != NULL && IVmin(p3, dsizes3, &i) <= 0 ) {
   fprintf(stderr, "\n IVmin(p3, dsizes3) must be > 0") ;
   return ; 
}
nvtx = n1*n2*n3 ;
/*
   ----------------------------------
   construct the domain sizes vectors
   ----------------------------------
*/
if ( dsizes1 == NULL ) {
   length1s = IVinit(p1, 0) ;
   length1 = (n1 - p1 + 1) / p1 ;
   m1 = (n1 - p1 + 1) % p1 ;
   for ( i = 0 ; i < m1 ; i++ ) {
      length1s[i] = length1 + 1 ;
   }
   for ( ; i < p1 ; i++ ) {
      length1s[i] = length1 ;
   }
} else {
   length1s = dsizes1 ;
}
if ( dsizes2 == NULL ) {
   length2s = IVinit(p2, 0) ;
   length2 = (n2 - p2 + 1) / p2 ;
   m2 = (n2 - p2 + 1) % p2 ;
   for ( i = 0 ; i < m2 ; i++ ) {
      length2s[i] = length2 + 1 ;
   }
   for ( ; i < p2 ; i++ ) {
      length2s[i] = length2 ;
   }
} else {
   length2s = dsizes2 ;
}
if ( dsizes3 == NULL ) {
   length3s = IVinit(p3, 0) ;
   length3 = (n3 - p3 + 1) / p3 ;
   m3 = (n3 - p3 + 1) % p3 ;
   for ( i = 0 ; i < m3 ; i++ ) {
      length3s[i] = length3 + 1 ;
   }
   for ( ; i < p3 ; i++ ) {
      length3s[i] = length3 ;
   }
} else {
   length3s = dsizes3 ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n inside localND3D") ;
fprintf(stdout, 
        "\n n1 = %d, n2 = %d, n3 = %d, p1 = %d, p2 = %dm p3 = %d", 
        n1, n2, n3, p1, p2, p3) ;
fprintf(stdout, "\n length1s[%d] = ", p1) ;
IVfp80(stdout, p1, length1s, 12) ;
fprintf(stdout, "\n length2s[%d] = ", p2) ;
IVfp80(stdout, p2, length2s, 12) ;
fprintf(stdout, "\n length3s[%d] = ", p3) ;
IVfp80(stdout, p3, length3s, 13) ;
#endif
/*
   ---------------------------------------
   determine the first and last domain ids 
   and the array of southwest points
   ---------------------------------------
*/
isws = IVinit(p1, -1) ;
for ( idom = 0, isw = 0 ; idom < p1 ; idom++ ) {
   isws[idom] = isw ;
   isw += length1s[idom] + 1 ;
}
jsws = IVinit(p2, -1) ;
for ( jdom = 0, jsw = 0 ; jdom < p2 ; jdom++ ) {
   jsws[jdom] = jsw ;
   jsw += length2s[jdom] + 1 ;
}
ksws = IVinit(p3, -1) ;
for ( kdom = 0, ksw = 0 ; kdom < p3 ; kdom++ ) {
   ksws[kdom] = ksw ;
   ksw += length3s[kdom] + 1 ;
}
#if MYDEBUG > 1
fprintf(stdout, "\n isws[%d] = ", p1) ;
IVfp80(stdout, p1, isws, 12) ;
fprintf(stdout, "\n jsws[%d] = ", p2) ;
IVfp80(stdout, p2, jsws, 12) ;
fprintf(stdout, "\n ksws[%d] = ", p3) ;
IVfp80(stdout, p3, ksws, 12) ;
#endif
/*
   ----------------------------------------------------------------
   create a temporary permutation vector for the domains' orderings
   ----------------------------------------------------------------
*/
msize = IVmax(p1, length1s, &i) * IVmax(p2, length2s, &i) 
                                * IVmax(p3, length3s, &k) ;
temp  = IVinit(msize, -1) ;
/*
   ------------------------
   fill in the domain nodes
   ------------------------
*/
now = 0 ;
for ( kdom = 0 ; kdom < p3 ; kdom++ ) {
   ksw     = ksws[kdom] ;
   length3 = length3s[kdom] ;
   for ( jdom = 0 ; jdom < p2 ; jdom++ ) {
      jsw     = jsws[jdom] ;
      length2 = length2s[jdom] ;
      for ( idom = 0 ; idom < p1 ; idom++ ) {
         isw     = isws[idom] ;
         length1 = length1s[idom] ;
/*
fprintf(stdout, "\n domain (%d,%d,%d), size %d x %d x %d",
        idom, jdom, kdom, length1, length2, length3) ;
fprintf(stdout, "\n (isw, jsw, ksw) = (%d, %d, %d)",
        isw, jsw, ksw) ;
*/
         mkNDperm(length1, length2, length3, temp, 
                  0, length1-1, 0, length2-1, 0, length3-1) ;
         for ( m = 0 ; m < length1*length2*length3 ; m++ ) {
            ijk = temp[m] ;
/*
fprintf(stdout, "\n    m = %d, ijk = %d", m, ijk) ;
*/
            k   = ksw + ijk / (length1*length2) ;
            ijk = ijk % (length1*length2) ;
            j   = jsw + ijk / length1 ;
            i   = isw + ijk % length1 ;
/*
fprintf(stdout, ", (i, j, k) = (%d, %d, %d)", i, j, k) ;
*/
            ijk = i + j*n1 + k*n1*n2 ;
            oldToNew[ijk] = now++ ;
         }
      }
   }
}
#if MYDEBUG > 2
fprintf(stdout, "\n old-to-new after domains are numbered") ;
fp3DGrid(n1, n2, n3, oldToNew, stdout) ;
#endif
/*
   ---------------------------------
   fill in the lower separator nodes
   ---------------------------------
*/
for ( kdom = 0 ; kdom < (p3/2) ; kdom++ ) {
   ksw  = ksws[kdom] ;
   length3   = length3s[kdom] ;
   for ( jdom = 0 ; jdom < p2 ; jdom++ ) {
      jsw  = jsws[jdom] ;
      length2   = length2s[jdom] ;
      for ( idom = 0 ; idom < p1 ; idom++ ) {
         isw  = isws[idom] ;
         length1   = length1s[idom] ;
/*
   -------
   3 faces
   -------
*/
         if ( isw > 0 ) {
            i = isw - 1 ;
            for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
               for ( k = ksw ; k <= ksw + length3 - 1 ; k++ ) {
                  ijk = i + j*n1 + k*n1*n2 ;
                  oldToNew[ijk] = now++ ;
               }
            }
         }
         if ( jsw > 0 ) {
            j = jsw - 1 ;
            for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
               for ( k = ksw ; k <= ksw + length3 - 1 ; k++ ) {
                  ijk = i + j*n1 + k*n1*n2 ;
                  oldToNew[ijk] = now++ ;
               }
            }
         }
         if ( ksw > 0 ) {
            k = ksw - 1 ;
            for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
               for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
                  ijk = i + j*n1 + k*n1*n2 ;
                  oldToNew[ijk] = now++ ;
               }
            }
         }
/*
         -----------
         three edges
         -----------
*/
         if ( isw > 0 && jsw > 0 ) {
            i = isw - 1 ;
            j = jsw - 1 ;
            for ( k = ksw ; k <= ksw + length3 - 1 ; k++ ) {
               ijk = i + j*n1 + k*n1*n2 ;
               oldToNew[ijk] = now++ ;
            }
         }
         if ( isw > 0 && ksw > 0 ) {
            i = isw - 1 ;
            k = ksw - 1 ;
            for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
               ijk = i + j*n1 + k*n1*n2 ;
               oldToNew[ijk] = now++ ;
            }
         }
         if ( jsw > 0 && ksw > 0 ) {
            j = jsw - 1 ;
            k = ksw - 1 ;
            for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
               ijk = i + j*n1 + k*n1*n2 ;
               oldToNew[ijk] = now++ ;
            }
         }
/*
         ----------------
         one corner point
         ----------------
*/
         if ( isw > 0 && jsw > 0 && ksw > 0 ) {
            i = isw - 1 ;
            j = jsw - 1 ;
            k = ksw - 1 ;
            ijk = i + j*n1 + k*n1*n2 ;
            oldToNew[ijk] = now++ ;
         }
      }
   }
}
#if MYDEBUG > 2
fprintf(stdout, "\n after the lower separators filled in") ;
fp2DGrid(n1, n2, oldToNew, stdout) ;
#endif
/*
   ---------------------------------
   fill in the upper separator nodes
   ---------------------------------
*/
for ( kdom = p3 - 1 ; kdom >= (p3/2) ; kdom-- ) {
   ksw  = ksws[kdom] ;
   length3   = length3s[kdom] ;
   for ( jdom = p2 - 1 ; jdom >= 0 ; jdom-- ) {
      jsw  = jsws[jdom] ;
      length2   = length2s[jdom] ;
      for ( idom = p1 - 1 ; idom >= 0 ; idom-- ) {
         isw  = isws[idom] ;
         length1   = length1s[idom] ;
/*
         -------
         3 faces
         -------
*/
         if ( isw + length1 < n1 ) {
            i = isw + length1 ;
            for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
               for ( k = ksw ; k <= ksw + length3 - 1 ; k++ ) {
                  ijk = i + j*n1 + k*n1*n2 ;
                  oldToNew[ijk] = now++ ;
               }
            }
         }
         if ( jsw + length2 < n2 ) {
            j = jsw + length2 ;
            for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
               for ( k = ksw ; k <= ksw + length3 - 1 ; k++ ) {
                  ijk = i + j*n1 + k*n1*n2 ;
                  oldToNew[ijk] = now++ ;
               }
            }
         }
         if ( ksw + length3 < n3 ) {
            k = ksw + length3 ;
            for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
               for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
                  ijk = i + j*n1 + k*n1*n2 ;
                  oldToNew[ijk] = now++ ;
               }
            }
         }
/*
         -----------
         three edges
         -----------
*/
         if ( isw + length1 < n1 && jsw + length2 < n2 ) {
            i = isw + length1 ;
            j = jsw + length2 ;
            for ( k = ksw ; k <= ksw + length3 - 1 ; k++ ) {
               ijk = i + j*n1 + k*n1*n2 ;
               oldToNew[ijk] = now++ ;
            }
         }
         if ( isw + length1 < n1 && ksw + length3 < n3 ) {
            i = isw + length1 ;
            k = ksw + length3 ;
            for ( j = jsw ; j <= jsw + length2 - 1 ; j++ ) {
               ijk = i + j*n1 + k*n1*n2 ;
               oldToNew[ijk] = now++ ;
            }
         }
         if ( jsw + length2 < n2 && ksw + length3 < n3 ) {
            j = jsw + length2 ;
            k = ksw + length3 ;
            for ( i = isw ; i <= isw + length1 - 1 ; i++ ) {
               ijk = i + j*n1 + k*n1*n2 ;
               oldToNew[ijk] = now++ ;
            }
         }
/*
         ----------------
         one corner point
         ----------------
*/
         if ( isw + length1 < n1 && jsw + length2 < n2 
                                 && ksw + length3 < n3 ) {
            i = isw + length1 ;
            j = jsw + length2 ;
            k = ksw + length3 ;
            ijk = i + j*n1 + k*n1*n2 ;
            oldToNew[ijk] = now++ ;
         }
      }
   }
}
#if MYDEBUG > 2
fprintf(stdout, "\n after the upper separators filled in") ;
fp2DGrid(n1, n2, oldToNew, stdout) ;
#endif
/*
   -------------------------------
   fill in the top level separator
   -------------------------------
*/
m1 = p3 / 2 ;
for ( kdom = 0, k = 0 ; kdom < m1 ; kdom++ ) {
   k += length3s[kdom] + 1 ;
}
k-- ;
for ( j = 0 ; j < n2 ; j++ ) { 
   for ( i = 0 ; i < n1 ; i++ ) { 
      ijk = i + j*n1 + k*n1*n2 ;
      oldToNew[ijk] = now++ ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( dsizes1 == NULL ) {
   IVfree(length1s) ;
}
if ( dsizes2 == NULL ) {
   IVfree(length2s) ;
}
if ( dsizes3 == NULL ) {
   IVfree(length3s) ;
}
IVfree(isws) ;
IVfree(jsws) ;
IVfree(ksws) ;
IVfree(temp) ;

return ; }

/*--------------------------------------------------------------------*/
