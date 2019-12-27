/*  basics.C  */

#include "../Coords.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95dec17, cca
   ----------------------------------------------
*/
int
Coords_sizeOf ( 
   Coords   *coords
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if (  coords == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_sizeof(%p)"
           "\n bad input\n", coords) ;
   exit(-1) ;
}
nbytes = sizeof(struct _Coords) ;
if ( coords->ndim > 0 && coords->ncoor > 0 ) {
   nbytes += coords->ndim * coords->ncoor * sizeof(float) ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   return the minimum coordinate in the dim-direction

   created -- 95dec17, cca
   --------------------------------------------------
*/
float
Coords_min ( 
   Coords   *coords,
   int      dim 
) {
float   minval ;
float   *coors ;
int     i, j, ncoor, ndim ;

if ( coords == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_min(%p,%d)"
           "\n bad input \n", coords, dim) ;
   exit(-1) ;
} else if ( !(coords->type == COORDS_BY_TUPLE
           || coords->type == COORDS_BY_COORD) ) {
   fprintf(stderr, "\n fatal error in Coords_min(%p,%d)"
           "\n coords->type = %d",
           coords, dim, coords->type) ;
   exit(-1) ;
} else if ( (ndim = coords->ndim) < 1 ) {
   fprintf(stderr, "\n fatal error in Coords_min(%p,%d)"
           "\n coords->ndim = %d",
           coords, dim, coords->ndim) ;
   exit(-1) ;
} else if ( (ncoor = coords->ncoor) < 1 ) {
   fprintf(stderr, "\n fatal error in Coords_min(%p,%d)"
           "\n coords->ncoor = %d",
           coords, dim, coords->ncoor) ;
   exit(-1) ;
} else if ( (coors = coords->coors) == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_min(%p,%d)"
           "\n coords->coords = %p",
           coords, dim, coords->coors) ;
   exit(-1) ;
} else if ( dim <= 0 || dim > coords->ndim ) {
   fprintf(stderr, "\n fatal error in Coords_min(%p,%d)"
           "\n bad input value, dim %d, ndim %d", 
           coords, dim, dim, ndim) ;
   exit(-1) ;
}
switch ( coords->type ) {
case COORDS_BY_TUPLE :
   i      = dim - 1 ;
   minval = coors[i] ;
   for ( j = 1, i += ndim ; j < ncoor ; j++, i += ndim ) {
      if ( coors[i] < minval ) {
         minval = coors[i] ;
      }
   }
   break ;
case COORDS_BY_COORD :
   i      = (dim - 1)*ncoor ;
   minval = coors[i] ;
   for ( i++ ; i < ncoor ; i++ ) {
      if ( coors[i] < minval ) {
         minval = coors[i] ;
      }
   }
   break ;
}
return(minval) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   return the maximum coordinate in the dim-direction

   created -- 95dec17, cca
   --------------------------------------------------
*/
float
Coords_max ( 
   Coords   *coords,
   int      dim 
) {
float   maxval ;
float   *coors ;
int     i, j, ncoor, ndim ;

if ( coords == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_max(%p,%d)"
           "\n bad input \n", coords, dim) ;
   exit(-1) ;
} else if ( !(coords->type == COORDS_BY_TUPLE
           || coords->type == COORDS_BY_COORD) ) {
   fprintf(stderr, "\n fatal error in Coords_max(%p,%d)"
           "\n coords->type = %d",
           coords, dim, coords->type) ;
   exit(-1) ;
} else if ( (ndim = coords->ndim) < 1 ) {
   fprintf(stderr, "\n fatal error in Coords_max(%p,%d)"
           "\n coords->ndim = %d",
           coords, dim, coords->ndim) ;
   exit(-1) ;
} else if ( (ncoor = coords->ncoor) < 1 ) {
   fprintf(stderr, "\n fatal error in Coords_max(%p,%d)"
           "\n coords->ncoor = %d",
           coords, dim, coords->ncoor) ;
   exit(-1) ;
} else if ( (coors = coords->coors) == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_max(%p,%d)"
           "\n coords->coords = %p",
           coords, dim, coords->coors) ;
   exit(-1) ;
} else if ( dim <= 0 || dim > coords->ndim ) {
   fprintf(stderr, "\n fatal error in Coords_max(%p,%d)"
           "\n bad input value, dim %d, ndim %d", 
           coords, dim, dim, ndim) ;
   exit(-1) ;
}
switch ( coords->type ) {
case COORDS_BY_TUPLE :
   i      = dim - 1 ;
   maxval = coors[i] ;
   for ( j = 1, i += ndim ; j < ncoor ; j++, i += ndim ) {
      if ( coors[i] > maxval ) {
         maxval = coors[i] ;
      }
   }
   break ;
case COORDS_BY_COORD :
   i      = (dim - 1)*ncoor ;
   maxval = coors[i] ;
   for ( i++ ; i < ncoor ; i++ ) {
      if ( coors[i] > maxval ) {
         maxval = coors[i] ;
      }
   }
   break ;
}
return(maxval) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   1 <= idim  <= ndim
   0 <= icoor < ncoor

   created -- 95dec17, cca
   -----------------------
*/
float
Coords_value ( 
   Coords   *coords,
   int      idim, 
   int      icoor 
) {
float   val ;
/*
   ---------------
   check the input
   ---------------
*/
if (  coords == NULL || idim <= 0 
   || !(coords->type == COORDS_BY_TUPLE
         || coords->type == COORDS_BY_COORD)
   || icoor < 0 || coords->ncoor <= icoor || coords->coors == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_value(%p,%d,%d)"
           "\n bad input or object\n", coords, idim, icoor) ;
   exit(-1) ;
}
switch ( coords->type ) {
case COORDS_BY_TUPLE :
   val = coords->coors[idim - 1 + coords->ndim*icoor] ;
   break ;
case COORDS_BY_COORD :
   val = coords->coors[icoor + coords->ncoor*(idim - 1)] ;
   break ;
}
return(val) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   1 <= idim  <= ndim
   0 <= icoor < ncoor

   created -- 95dec17, cca
   -----------------------
*/
void
Coords_setValue ( 
   Coords   *coords,
   int      idim, 
   int      icoor, 
   float    val 
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  coords == NULL || idim <= 0 
   || !(coords->type == COORDS_BY_TUPLE
         || coords->type == COORDS_BY_COORD)
   || icoor < 0 || coords->ncoor <= icoor || coords->coors == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_setValue(%p,%d,%d,%f)"
           "\n bad input or object\n", coords, idim, icoor, val) ;
   exit(-1) ;
}
switch ( coords->type ) {
case COORDS_BY_TUPLE :
   coords->coors[idim - 1 + coords->ndim*icoor] = val ;
   break ;
case COORDS_BY_COORD :
   coords->coors[icoor + coords->ncoor*(idim - 1)] = val ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
