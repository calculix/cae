/*  Coords.h  */

#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   object to hold coordinates of a graph

   type  -- coordinates type
      1 -- use tuple, x(icoor, idim) = coors[idim + ndim*icoor]
      2 -- use bycoord, x(icoor, idim) = coors[icoor + ncoor*idim]
   ndim  -- number of dimensions
   ncoor -- number of coordinates
   coors -- coordinate array
   ---------------------------------------------------------------
*/
typedef struct _Coords   Coords ;
struct _Coords {
   int     type   ;
   int     ndim   ;
   int     ncoor  ;
   float   *coors ;
} ;
#define COORDS_BY_TUPLE 1
#define COORDS_BY_COORD 2
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor

   created -- 95dec17, cca
   -----------------------
*/
Coords *
Coords_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields

   created -- 95dec17, cca
   -----------------------
*/
void
Coords_setDefaultFields (
   Coords   *coords
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 95dec17, cca
   --------------------------------------------------
*/
void
Coords_clearData ( 
   Coords   *coords 
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 95dec17, cca
   ------------------------------------------
*/
Coords *
Coords_free ( 
   Coords   *coords 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
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
) ;
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
) ;
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
   int      Type, 
   int      n1, 
   int      n2, 
   int      n3, 
   int      ncomp 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95dec17, cca
   ----------------------------------------------
*/
int
Coords_sizeOf ( 
   Coords   *coords
) ;
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
) ;
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
) ;
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------
   purpose -- to read in an Coords object from a file

   input --

      fn -- filename, must be *.coordsb or *.coordsf

   return value -- 1 if success, 0 if failure

   created -- 95dec19, cca
   ---------------------------------------------------
*/
int
Coords_readFromFile ( 
   Coords    *coords, 
   char   *fn 
) ;
/*
   ---------------------------------------------------------
   purpose -- to read an Coords object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95dec19, cca
   ---------------------------------------------------------
*/
int
Coords_readFromFormattedFile ( 
   Coords    *coords, 
   FILE   *fp 
) ;
/*
   ------------------------------------------------------
   purpose -- to read an Coords object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95dec19, cca
   ------------------------------------------------------
*/
int
Coords_readFromBinaryFile ( 
   Coords    *coords, 
   FILE   *fp 
) ;
/*
   ----------------------------------------------
   purpose -- to write an Coords object to a file

   input --

      fn -- filename
        *.coordsb -- binary
        *.coordsf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   ----------------------------------------------
*/
int
Coords_writeToFile ( 
   Coords    *coords, 
   char   *fn 
) ;
/*
   --------------------------------------------------------
   purpose -- to write an Coords object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   --------------------------------------------------------
*/
int
Coords_writeToFormattedFile ( 
   Coords    *coords, 
   FILE   *fp 
) ;
/*
   -----------------------------------------------------
   purpose -- to write an Coords object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   -----------------------------------------------------
*/
int
Coords_writeToBinaryFile ( 
   Coords    *coords, 
   FILE   *fp 
) ;
/*
   ----------------------------------------------------
   purpose -- to write an Coords object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   ----------------------------------------------------
*/
int
Coords_writeForHumanEye ( 
   Coords    *coords, 
   FILE   *fp 
) ;
/*
   ------------------------------------------------------------
   purpose -- to write out the statistics for the Coords object

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   ------------------------------------------------------------
*/
int
Coords_writeStats ( 
   Coords    *coords, 
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
