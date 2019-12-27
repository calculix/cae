/*  ND2.c  */

#include "../misc.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   this file contains procedures to create a nested dissection 
   permutation in one, two or three dimensions.

   note, all separators are double wide

   created -- 96feb01, cca
   ---------------------------------------------------------------
*/

#define DEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------
   static procedures
   -----------------
*/
static void   SplitX ( int n1, int n2, int n3, int newToOld[], 
                       int west, int east, int south, int north,
                       int bottom, int top ) ;
static void   SplitY ( int n1, int n2, int n3, int newToOld[], 
                       int west, int east, int south, int north,
                       int bottom, int top ) ;
static void   SplitZ ( int n1, int n2, int n3, int newToOld[], 
                       int west, int east, int south, int north,
                       int bottom, int top ) ;
static int    WhichCut ( int west, int east, int south, int north,
                         int bottom, int top ) ;
static int    MidPoint ( int left, int right, int glob_center ) ;

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- main procedure. if the input region is m1 x m2 x m3,
              where 0 < m1, m2, m3 < 2, then
              the node is put into the permutation vector.
              otherwise the region is split into three pieces,
              two subregions and a separator, and recursive 
              calls are made to order the subregions

   input --

      n1         -- number of points in the first direction
      n2         -- number of points in the second direction
      n3         -- number of points in the third direction
      newToOld -- pointer to the permutation vector
      west       -- west coordinate
      east       -- east coordinate
      south      -- south coordinate
      north      -- north coordinate
      bottom     -- bottom coordinate
      top        -- top coordinate

   created -- 95nov15, cca
   ------------------------------------------------------------
*/
void
mkNDperm2 ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   newToOld[], 
   int   west, 
   int   east, 
   int   south, 
   int   north, 
   int   bottom, 
   int   top 
) {
int   count, i, j, k ;
if ( n1 <= 0 || n2 <= 0 || n3 <= 0 || newToOld == NULL
   || west   < 0 || east  >= n1
   || south  < 0 || north >= n2
   || bottom < 0 || top   >= n3 ) {
   fprintf(stderr, 
           "\n fatal error in mkND2perm(%d,%d,%d,%p,%d,%d,%d,%d,%d,%d)"
           "\n bad input data\n",
           n1, n2, n3, newToOld, 
           west, east, south, north, bottom, top) ;
   exit(-1) ;
}

if ( east - west <= 1 && north - south <= 1 && top - bottom <= 1 ) {
   count = 0 ;
   for ( i = west ; i <= east ; i++ ) {
      for ( j = south ; j <= north ; j++ ) {
         for ( k = bottom ; k <= top ; k++ ) {
#     if DEBUG > 0
            fprintf(stdout, "\n ND : ordering %d",
                    i + j * n1 + k * n1 * n2) ;
#     endif
            newToOld[count++] = i + j * n1 + k * n1 * n2 ; 
         }
      }
   }
} else {
   switch ( WhichCut(west, east, south, north, bottom, top) ) {
   case 1 : 
#     if DEBUG > 0
      fprintf(stdout, "\n ND : calling SplitX(%d,%d,%d,%d,%d,%d)",
              west, east, south, north, bottom, top) ; 
#     endif
      SplitX(n1, n2, n3, newToOld, 
             west, east, south, north, bottom, top) ; 
      break ;
   case 2 : 
#     if DEBUG > 0
      fprintf(stdout, "\n ND : calling SplitY(%d,%d,%d,%d,%d,%d)",
              west, east, south, north, bottom, top) ; 
#     endif
      SplitY(n1, n2, n3, newToOld, 
             west, east, south, north, bottom, top) ; 
      break ;
   case 3 : 
#     if DEBUG > 0
      fprintf(stdout, "\n ND : calling SplitZ(%d,%d,%d,%d,%d,%d)",
              west, east, south, north, bottom, top) ; 
#     endif
      SplitZ(n1, n2, n3, newToOld, 
             west, east, south, north, bottom, top) ; 
      break ; } }

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to decide which way to partition a subregion

   created -- 95nov16, cca
   -------------------------------------------------------
*/
static int
WhichCut ( 
   int   west, 
   int   east, 
   int   south, 
   int   north, 
   int   bottom, 
   int   top 
) {
int   d1, d2, d3 ;

d1 = east - west + 1 ;
d2 = north - south + 1 ;
d3 = top - bottom + 1 ;
#if DEBUG > 0
   fprintf(stdout, "\n WhichCut : (d1,d2,d3) = (%d,%d,%d)", d1,d2,d3) ;
#endif
if ( d1 > d2 && d1 > d3 ) {
   return(1) ; }
else if ( d2 > d1 && d2 > d3 ) {
   return(2) ; }
else if ( d3 > d1 && d3 > d2 ) {
   return(3) ; }
else if ( d1 == d2 && d1 > d3 ) {
   return(1) ; }
else if ( d1 == d3 && d1 > d2 ) {
   return(1) ; }
else if ( d2 == d3 && d2 > d1 ) {
   return(2) ; }
else {
   return(3) ; } }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to split a subregion with a Y-Z plane
   ------------------------------------------------
*/
static void
SplitX ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   newToOld[], 
   int   west,  
   int   east,  
   int   south,  
   int   north,  
   int   bottom,  
   int   top 
) {
int   depth, east1, east2, height, j, k, k1, k2, k3, m, mid, 
      west1, west2, width1, width2 ;

#  if DEBUG > 0
   fprintf(stdout, "\n inside SplitX(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom, top) ;
#  endif
mid    = MidPoint(west, east, n1/2) ;
west1  = west    ;
east1  = mid - 1 ;
west2  = mid + 2 ;
east2  = east    ;
width1 = east1 - west1  + 1 ;
width2 = east2 - west2  + 1 ;
height = north - south  + 1 ;
depth  = top   - bottom + 1 ;
k1 = 0 ;
k2 = k1 + width1 * height * depth ;
k3 = k2 + width2 * height * depth ;
if ( width1 > 0 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitX : 1. calling ND(%d,%d,%d,%d,%d,%d",
           west1, east1, south, north, bottom, top) ; 
#  endif
   mkNDperm2(n1, n2, n3, newToOld + k1, 
      west1, east1, south, north, bottom, top) ; 
}
if ( width2 > 0 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitX : 2. calling ND(%d,%d,%d,%d,%d,%d",
           west2, east2, south, north, bottom, top) ; 
#  endif
   mkNDperm2(n1, n2, n3, newToOld + k2, 
      west2, east2, south, north, bottom, top) ; 
}
m = k3 ;
for ( k = bottom ; k <= top ; k++ ) {
   for ( j = south ; j <= north ; j++ ) {
#  if DEBUG > 0
      fprintf(stdout, "\n newToOld[%d] = %d",
              m, mid +     j * n1 + k * n1 * n2) ;
      fprintf(stdout, "\n newToOld[%d] = %d",
              m+1, mid + 1 + j * n1 + k * n1 * n2) ;
#  endif
      newToOld[m++] = mid +     j * n1 + k * n1 * n2 ;
      newToOld[m++] = mid + 1 + j * n1 + k * n1 * n2 ; 
   } 
}
#  if DEBUG > 0
   fprintf(stdout, "\n leaving SplitX(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom, top) ;
#  endif

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to split a subregion with a X-Z plane

   created -- 95nov16, cca
   ------------------------------------------------
*/
static void
SplitY ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   newToOld[], 
   int   west, 
   int   east, 
   int   south, 
   int   north,
   int   bottom, 
   int   top 
) {
int   depth, height1, height2, i, k, k1, k2, k3, m, mid, 
      north1, north2, south1, south2, width ;

#  if DEBUG > 0
   fprintf(stdout, "\n inside SplitY(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom, top) ;
#  endif
mid     = MidPoint(south, north, n2/2) ;
south1  = south ;
north1  = mid - 1 ;
south2  = mid + 2 ;
north2  = north ;
#  if DEBUG > 0
   fprintf(stdout, 
      "\n mid = %d, south1 = %d, north1 = %d, south2 = %d, north2 = %d",
      mid, south1, north1, south2, north2) ;
#  endif
width   = east   - west   + 1 ;
height1 = north1 - south1 + 1 ;
height2 = north2 - south2 + 1 ;
depth   = top    - bottom + 1 ;
k1 = 0 ;
k2 = k1 + width * height1 * depth ;
k3 = k2 + width * height2 * depth ;
if ( height1 > 0 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitY : calling ND(%d:%d,%d:%d,%d:%d)",
           west, east, south1, north1, bottom, top) ;
#  endif
   mkNDperm2(n1, n2, n3, newToOld + k1, 
      west, east, south1, north1, bottom, top) ; 
}
if ( height2 > 0 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitY : calling ND(%d:%d,%d:%d,%d:%d)",
           west, east, south2, north2, bottom, top) ;
#  endif
   mkNDperm2(n1, n2, n3, newToOld + k2, 
      west, east, south2, north2, bottom, top) ; 
}
m = k3 ;
for ( k = bottom ; k <= top ; k++ ) {
   for ( i = west ; i <= east ; i++ ) {
#     if DEBUG > 0
      fprintf(stdout, "\n SplitY : newToOld[%d] = %d",
              m, i + mid * n1 + k * n1 * n2) ;
      fprintf(stdout, "\n SplitY : newToOld[%d] = %d",
              m+1, i + (mid + 1) * n1 + k * n1 * n2) ;
#     endif
      newToOld[m++] = i + mid*n1       + k*n1*n2 ;
      newToOld[m++] = i + (mid + 1)*n1 + k*n1*n2 ; 
   } 
}
#  if DEBUG > 0
   fprintf(stdout, "\n leaving SplitY(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom, top) ;
#  endif

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to split a subregion with a X-Y plane

   created -- 95nov16, cca
   ------------------------------------------------
*/
static void
SplitZ ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   newToOld[], 
   int   west, 
   int   east, 
   int   south, 
   int   north,
   int   bottom, 
   int   top 
) {
int   bottom1, bottom2, depth1, depth2, height, i, j, k1, k2, k3, 
      m, mid, top1, top2, width ;

#  if DEBUG > 0
   fprintf(stdout, "\n inside SplitZ(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom, top) ;
#  endif
mid = MidPoint(bottom, top, n2/2) ;
bottom1 = bottom  ;
top1    = mid - 1 ;
bottom2 = mid + 2 ;
top2    = top     ;
width   = east  - west    + 1 ;
height  = north - south   + 1 ;
depth1  = top1  - bottom1 + 1 ;
depth2  = top2  - bottom2 + 1 ;
k1 = 0 ;
k2 = k1 + width * height * depth1 ;
k3 = k2 + width * height * depth2 ;
if ( depth1 > 0 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitZ : 1. calling ND(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom1, top1) ;
#  endif
   mkNDperm2(n1, n2, n3, newToOld + k1, 
      west, east, south, north, bottom1, top1) ; 
}
if ( depth2 > 0 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitZ : 2. calling ND(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom2, top2) ;
#  endif
   mkNDperm2(n1, n2, n3, newToOld + k2, 
      west, east, south, north, bottom2, top2) ; 
}
m = k3 ;
for ( j = south ; j <= north ; j++ ) {
   for ( i = west ; i <= east ; i++ ) {
#  if DEBUG > 0
      fprintf(stdout, "\n newToOld[%d] = %d",
              m, i + j*n1 + mid*n1*n2) ;
      fprintf(stdout, "\n newToOld[%d] = %d",
              m+1, i + j*n1 + (mid+1)*n1*n2) ;
#  endif
      newToOld[m++] = i + j*n1 + mid*n1*n2 ; 
      newToOld[m++] = i + j*n1 + (mid + 1)*n1*n2 ; 
   }  
}
#  if DEBUG > 0
   fprintf(stdout, "\n leaving SplitZ(%d:%d,%d:%d,%d:%d)",
           west, east, south, north, bottom, top) ;
#  endif

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the midpoint of a region
   
   input --

      left        -- left coordinate
      right       -- right coordinate
      glob_center -- glob_center coordinate
   ----------------------------------------------
*/
static int 
MidPoint ( 
   int   left, 
   int   right, 
   int   glob_center 
) {
int   mid ;

/*
mid = (left + right)/2 ;
if ( (left + right) % 2 == 0 ) {
   return(mid) ; }
else if ( mid >= glob_center ) {
   return(mid) ; }
else {
   return(mid+1) ; }
*/
mid = left + (right - left - 1)/2 ;

return(mid) ; }

/*--------------------------------------------------------------------*/
