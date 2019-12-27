/*  init.C  */

#include "../Coords.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   initialize the object, storage freed if repeat call

   created -- 95dec17, cca
   ---------------------------------------------------
*/
void
Coords_init ( 
   Coords   *coords,
   int      type, 
   int      ndim, 
   int      ncoor 
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  coords == NULL 
   || !(type == COORDS_BY_TUPLE || type == COORDS_BY_COORD)
   || ndim < 1 || ncoor <= 0 ) {
   fprintf(stderr, "\n fatal error in Coords_init(%p,%d,%d,%d)"
           "\n bad input\n", coords, type, ndim, ncoor) ;
   exit(-1) ;
}
/*
   --------------
   clear the data
   --------------
*/
Coords_clearData(coords) ;
/*
   -------------------
   set the data fields
   -------------------
*/
coords->type  = type ;
coords->ndim  = ndim ;
coords->ncoor = ncoor ;
coords->coors = FVinit(ndim*ncoor, 0.0) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- initialize a 9-point operator on a n1 x n2 grid

   input --

      bbox  -- bounding box
         bbox[0] -- x southwest
         bbox[1] -- y southwest
         bbox[2] -- x northeast
         bbox[3] -- y northeast
      Type  -- type of object
      n1    -- # of grid points in first direction
      n2    -- # of grid points in second direction
      ncomp -- # of components per grid point

   created -- 95dec17, cca
   ----------------------------------------------------------
*/
void
Coords_init9P ( 
   Coords   *coords,
   float    bbox[], 
   int      type, 
   int      n1, 
   int      n2, 
   int      ncomp 
) {
int      i, idof, j, k ;
float    deltax, deltay ;
float    *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || bbox == NULL 
     || !(type == COORDS_BY_TUPLE || type == COORDS_BY_COORD)
     || n1 <= 0 || n2 <= 0 || ncomp <= 0 ) {
   fprintf(stderr, "\n fatal error in Coords_init9P(%p,%p,%d,%d,%d,%d)"
           "\n bad input\n", coords, bbox, type, n1, n2, ncomp) ;
   exit(-1) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Coords_init(coords, type, 2, n1*n2) ;
/*
   --------------------
   fill the coordinates
   --------------------
*/
deltax = (bbox[2] - bbox[0]) / (n1-1) ;
deltay = (bbox[3] - bbox[1]) / (n2-1) ;
switch ( type ) {
case COORDS_BY_TUPLE :
   x = coords->coors ;
   y = coords->coors + 1 ;
   for ( j = 0, idof = 0 ; j < n2 ; j++ ) {
      for ( i = 0 ; i < n1 ; i++ ) {
         for ( k = 0 ; k < ncomp ; k++, idof += 2 ) {
            x[idof] = bbox[0] + deltax*i ;
            y[idof] = bbox[1] + deltay*j ;
         }
      }
   }
   break ;
case COORDS_BY_COORD :
   x = coords->coors ;
   y = coords->coors + n1*n2 ;
   for ( j = 0, idof = 0 ; j < n2 ; j++ ) {
      for ( i = 0 ; i < n1 ; i++ ) {
         for ( k = 0 ; k < ncomp ; k++, idof++ ) {
            x[idof] = bbox[0] + deltax*i ;
            y[idof] = bbox[1] + deltay*j ;
         }
      }
   }
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- initialize a 27-point operator on a n1 x n2 x n3 grid

   input --

      bbox  -- bounding box
         bbox[0] -- x southwest bottom
         bbox[1] -- y southwest bottom
         bbox[2] -- z southwest bottom
         bbox[3] -- x northeast top
         bbox[4] -- y northeast top
         bbox[5] -- z northeast top
      Type  -- type of object
      n1    -- # of grid points in first direction
      n2    -- # of grid points in second direction
      n3    -- # of grid points in third direction
      ncomp -- # of components per grid point

   created -- 95dec17, cca
   ----------------------------------------------------------------
*/
void
Coords_init27P ( 
   Coords   *coords,
   float    bbox[], 
   int      type, 
   int      n1, 
   int      n2, 
   int      n3, 
   int      ncomp 
) {
int      i, idof, j, k, l ;
float    deltax, deltay, deltaz ;
float    *x, *y, *z ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || bbox == NULL 
     || !(type == COORDS_BY_TUPLE || type == COORDS_BY_COORD)
     || n1 <= 0 || n2 <= 0 || n3 <= 0 || ncomp <= 0 ) {
   fprintf(stderr, 
           "\n fatal error in Coords_init27P(%p,%p,%d,%d,%d,%d,%d)"
           "\n bad input\n", coords, bbox, type, n1, n2, n3, ncomp) ;
   exit(-1) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Coords_init(coords, type, 3, n1*n2*n3) ;
/*
   --------------------
   fill the coordinates
   --------------------
*/
deltax = (bbox[3] - bbox[0]) / (n1-1) ;
deltay = (bbox[4] - bbox[1]) / (n2-1) ;
deltaz = (bbox[5] - bbox[2]) / (n3-1) ;
switch ( type ) {
case COORDS_BY_TUPLE :
   x = coords->coors ;
   y = coords->coors + 1 ;
   z = coords->coors + 2 ;
   for ( k = 0, idof = 0 ; k < n3 ; k++ ) {
      for ( j = 0 ; j < n2 ; j++ ) {
         for ( i = 0 ; i < n1 ; i++ ) {
            for ( l = 0 ; l < ncomp ; l++, idof += 3 ) {
               x[idof] = bbox[0] + deltax*i ;
               y[idof] = bbox[1] + deltay*j ;
               z[idof] = bbox[2] + deltaz*k ;
            }
         }
      }
   }
   break ;
case COORDS_BY_COORD :
   x = coords->coors ;
   y = coords->coors + coords->ncoor ;
   z = coords->coors + 2*coords->ncoor ;
   for ( k = 0, idof = 0 ; k < n3 ; k++ ) {
      for ( j = 0 ; j < n2 ; j++ ) {
         for ( i = 0 ; i < n1 ; i++ ) {
            for ( l = 0 ; l < ncomp ; l++, idof++ ) {
               x[idof] = bbox[0] + deltax*i ;
               y[idof] = bbox[1] + deltay*j ;
               z[idof] = bbox[2] + deltaz*k ;
            }
         }
      }
   }
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
