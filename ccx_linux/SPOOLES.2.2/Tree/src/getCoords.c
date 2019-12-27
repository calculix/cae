/*  getCoords.c  */

#include "../Tree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to get simple x[] and y[] coordinates 
              for the tree vertices

   return values --
      1 -- normal return
     -1 -- tree is NULL
     -2 -- heightflag is invalid
     -3 -- coordflag is invalid
     -4 -- xDV is NULL
     -5 -- yDV is NULL

   created -- 99jan07, cca
   ------------------------------------------------
*/
int
Tree_getSimpleCoords (
   Tree   *tree,
   char   heightflag,
   char   coordflag,
   DV     *xDV,
   DV     *yDV
) {
double   *x, *y ;
int      count, I, J, n, nleaves ;
int      *fch, *par, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n error in Tree_getSimpleCoords()"
           "\n tree is NULL\n") ;
   return(-1) ;
}
if ( heightflag != 'D' && heightflag != 'H' ) {
   fprintf(stderr, "\n error in Tree_getSimpleCoords()"
           "\n invalid heightflag = %c\n", heightflag) ;
   return(-2) ;
}
if ( coordflag != 'C' && coordflag != 'P' ) {
   fprintf(stderr, "\n error in Tree_getSimpleCoords()"
           "\n invalid coordflag = %c\n", coordflag) ;
   return(-3) ;
}
if ( xDV == NULL ) {
   fprintf(stderr, "\n error in Tree_getSimpleCoords()"
           "\n xDV is NULL\n") ;
   return(-4) ;
}
if ( yDV == NULL ) {
   fprintf(stderr, "\n error in Tree_getSimpleCoords()"
           "\n yDV is NULL\n") ;
   return(-5) ;
}
n   = tree->n   ;
par = tree->par ;
fch = tree->fch ;
sib = tree->sib ;
DV_setSize(xDV, n) ;
DV_setSize(yDV, n) ;
x = DV_entries(xDV) ;
y = DV_entries(yDV) ;
switch ( heightflag ) {
case 'D' : {
   int   J, K, maxdepth ;

   for ( J = Tree_preOTfirst(tree), maxdepth = 0 ;
         J != -1 ;
         J = Tree_preOTnext(tree, J) ) {
      if ( (K = par[J]) == -1 ) {
         y[J] = 0.0 ;
      } else {
         y[J] = y[K] + 1.0 ;
      }
      if ( maxdepth < y[J] ) {
         maxdepth = y[J] ;
      }
   }
   if ( coordflag == 'C' ) {
      for ( J = 0 ; J < n ; J++ ) {
         y[J] = maxdepth - y[J] ;
      }
   }
   } break ;
case 'H' : {
   int   height, I, J, maxheight ;

   for ( J = Tree_postOTfirst(tree), maxheight = 0 ;
         J != -1 ;
         J = Tree_postOTnext(tree, J) ) {
      if ( (I = fch[J]) == -1 ) {
         y[J] = 0.0 ;
      } else {
         height = y[I] ;
         for ( I = sib[I] ; I != -1 ; I = sib[I] ) {
            if ( height < y[I] ) {
               height = y[I] ;
            }
         }
         y[J] = height + 1.0 ;
      }
      if ( maxheight < y[J] ) {
         maxheight = y[J] ;
      }
   }
   if ( coordflag == 'P' ) {
      for ( J = 0 ; J < n ; J++ ) {
         y[J] = maxheight - y[J] ;
      }
   }
   } break ;
default :
   break ;
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n y") ;
   DV_writeForHumanEye(yDV, stdout) ;
#endif
DV_zero(xDV) ;
nleaves = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( fch[J] == -1 ) {
      x[J] = nleaves++ ;
   } else {
      for ( I = fch[J], count = 0 ; I != -1 ; I = sib[I] ) {
         x[J] += x[I] ;
         count++ ;
      }
      x[J] /= count ;
   }
}
if ( coordflag == 'C' ) {
   for ( J = 0 ; J < n ; J++ ) {
      x[J] = x[J] / nleaves ;
   }
} else {
   double   r, theta ;

   for ( J = 0 ; J < n ; J++ ) {
      theta = 6.283185 * x[J] / nleaves ;
      r     = y[J] ;
      x[J]  = r * cos(theta) ;
      y[J]  = r * sin(theta) ;
   }
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n x") ;
   DV_writeForHumanEye(xDV, stdout) ;
#endif
return(1) ; }

/*--------------------------------------------------------------------*/
