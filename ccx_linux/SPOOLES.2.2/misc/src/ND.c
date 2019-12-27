/*  ND.c  */

#include "../misc.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   this file contains procedures to create a nested dissection 
   permutation in one, two or three dimensions.
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
   purpose -- main procedure. if the input region is 1 x 1 x 1,
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
mkNDperm ( 
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
if ( n1 <= 0 || n2 <= 0 || n3 <= 0 || newToOld == NULL
   || west   < 0 || east  >= n1
   || south  < 0 || north >= n2
   || bottom < 0 || top   >= n3 ) {
   fprintf(stderr, 
           "\n fatal error in mkNDperm(%d,%d,%d,%p,%d,%d,%d,%d,%d,%d)"
           "\n bad input data\n",
           n1, n2, n3, newToOld, 
           west, east, south, north, bottom, top) ;
   exit(-1) ;
}
if ( west == east && south == north && bottom == top ) {
#     if DEBUG > 0
      fprintf(stdout, "\n ND : ordering %d",
              west + south * n1 + bottom * n1 * n2) ;
#     endif
   newToOld[0] = west + south * n1 + bottom * n1 * n2 ; }
else {
   switch ( WhichCut(west, east, south, north, bottom, top) ) {
   case 1 : 
#     if DEBUG > 0
      fprintf(stdout, "\n ND : calling Split9X") ;
#     endif
      SplitX(n1, n2, n3, newToOld, 
             west, east, south, north, bottom, top) ; 
      break ;
   case 2 : 
#     if DEBUG > 0
      fprintf(stdout, "\n ND : calling Split9Y") ;
#     endif
      SplitY(n1, n2, n3, newToOld, 
             west, east, south, north, bottom, top) ; 
      break ;
   case 3 : 
#     if DEBUG > 0
      fprintf(stdout, "\n ND : calling Split9Z") ;
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
int   east1, east2, j, k, k1, k2, k3, m, m1, west1, west2 ;

m1 = MidPoint(west, east, n1/2) ;
west1 = west ;
east1 = m1 - 1 ;
west2 = m1 + 1 ;
east2 = east ;
k1 = 0 ;
k2 = k1 + (m1 - west) * (north - south + 1) * (top - bottom + 1) ;
k3 = k2 + (east - m1) * (north - south + 1) * (top - bottom + 1) ;
if ( m1 > west ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitX : 1. calling ND") ;
#  endif
   mkNDperm(n1, n2, n3, newToOld + k1, 
      west1, east1, south, north, bottom, top) ; }
if ( east > m1 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitX : 2. calling ND") ;
#  endif
   mkNDperm(n1, n2, n3, newToOld + k2, 
      west2, east2, south, north, bottom, top) ; }
#  if DEBUG > 0
fprintf(stdout, "\n SplitX : 3. calling ND") ;
fprintf(stdout, "\n m1 = %d, south = %d, north = %d",
        m1, south, north) ;
#  endif
m = k3 ;
for ( k = bottom ; k <= top ; k++ ) {
   for ( j = south ; j <= north ; j++ ) {
      newToOld[m++] = m1 + j * n1 + k * n1 * n2 ; } }

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
int   i, k, k1, k2, k3, m, m1, north1, north2, south1, south2 ;

m1 = MidPoint(south, north, n2/2) ;
south1 = south ;
north1 = m1 - 1 ;
south2 = m1 + 1 ;
north2 = north ;
k1 = 0 ;
k2 = k1 + (m1 - south) * (east - west + 1) * (top - bottom + 1) ;
k3 = k2 + (north - m1) * (east - west + 1) * (top - bottom + 1) ;
if ( m1 > south ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitY : calling ND(%d:%d,%d:%d,%d:%d)",
           west, east, south1, north1, bottom, top) ;
#  endif
   mkNDperm(n1, n2, n3, newToOld + k1, 
      west, east, south1, north1, bottom, top) ; }
if ( north > m1 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitY : calling ND(%d:%d,%d:%d,%d:%d)",
           west, east, south2, north2, bottom, top) ;
#  endif
   mkNDperm(n1, n2, n3, newToOld + k2, 
      west, east, south2, north2, bottom, top) ; }
#  if DEBUG > 0
   fprintf(stdout, "\n SplitY : ordering (%d:%d,%d:%d,%d:%d)",
           west, east, m1, m1, bottom, top) ;
#  endif
m = k3 ;
for ( k = bottom ; k <= top ; k++ ) {
   for ( i = west ; i <= east ; i++ ) {
#     if DEBUG > 0
      fprintf(stdout, "\n SplitY : newToOld[%d] = %d",
              m, i + m1 * n1 + k * n1 * n2) ;
#     endif
      newToOld[m++] = i + m1 * n1 + k * n1 * n2 ; } }

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
int   bottom1, bottom2, i, j, k1, k2, k3, m, m1, top1, top2 ;

m1 = MidPoint(bottom, top, n2/2) ;
bottom1 = bottom ;
top1    = m1 - 1 ;
bottom2 = m1 + 1 ;
top2    = top    ;
k1 = 0 ;
k2 = k1 + (m1 - bottom) * (east - west + 1) * (north - south + 1) ;
k3 = k2 + (top - m1) * (east - west + 1) * (north - south + 1) ;
if ( m1 > bottom ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitZ : 1. calling ND") ;
#  endif
   mkNDperm(n1, n2, n3, newToOld + k1, 
      west, east, south, north, bottom1, top1) ; }
if ( top > m1 ) {
#  if DEBUG > 0
   fprintf(stdout, "\n SplitZ : 2. calling ND") ;
#  endif
   mkNDperm(n1, n2, n3, newToOld + k2, 
      west, east, south, north, bottom2, top2) ; }
#  if DEBUG > 0
fprintf(stdout, "\n SplitZ : 3. calling ND") ;
#  endif
m = k3 ;
for ( j = south ; j <= north ; j++ ) {
   for ( i = west ; i <= east ; i++ ) {
      newToOld[m++] = i + j * n1 + m1 * n1 * n2 ; } }

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

mid = (left + right)/2 ;
if ( (left + right) % 2 == 0 ) {
   return(mid) ; }
else if ( mid >= glob_center ) {
   return(mid) ; }
else {
   return(mid+1) ; } }

/*--------------------------------------------------------------------*/
