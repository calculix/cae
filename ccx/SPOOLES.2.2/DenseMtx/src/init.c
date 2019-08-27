/*  init.c  */

#include "../DenseMtx.h"

#define MYDEBUG 1

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   return the number of bytes needed for an object of this size
 
   created -- 98may02, cca
   ------------------------------------------------------------
*/
int
DenseMtx_nbytesNeeded (
   int   type,
   int   nrow,
   int   ncol
) {
int   nbytes, ndouble, nint ;
/*
   ---------------
   check the input
   ---------------
*/
if (  nrow < 0 || ncol < 0 ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_nbytesNeeded(%d,%d,%d)"
           "\n bad input\n", type, nrow, ncol) ;
   exit(-1) ;
}
nint = 7 + nrow + ncol ;
if ( type == SPOOLES_REAL ) {
   ndouble = nrow*ncol ;
} else if ( type == SPOOLES_COMPLEX ) {
   ndouble = 2*nrow*ncol ;
} else {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_nbytesNeeded(%d,%d,%d)"
           "\n bad type %d\n", type, nrow, ncol, type) ;
   exit(-1) ;
}
if ( sizeof(int) == sizeof(double) ) {
   nbytes = nint*sizeof(int) + ndouble*sizeof(double) ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   nbytes = ((nint + 1)/2 + ndouble)*sizeof(double) ;
} else {
   fprintf(stderr, "\n error in DenseMtx_nbytesNeeded(%d,%d)"
           "\n sizeof(int) = %d, sizeof(double) = %d",
           nrow, ncol, sizeof(int), sizeof(double)) ;
   exit(-1) ;
}
return(nbytes) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   return the number of bytes in the workspace owned by this object
 
   created -- 98may02, cca
   ----------------------------------------------------------------
*/
int
DenseMtx_nbytesInWorkspace (
   DenseMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_nbytesInWorkspace(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(sizeof(double)*DV_maxsize(&mtx->wrkDV)) ; }
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   set the number of bytes in the workspace owned by this object
 
   created -- 98may02, cca
   -------------------------------------------------------------
*/
void
DenseMtx_setNbytesInWorkspace (
   DenseMtx   *mtx,
   int        nbytes
) {
if ( mtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_setNbytesInWorkspace(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
DV_setSize(&mtx->wrkDV, nbytes/sizeof(double)) ;

return ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- set the fields of the object

   created -- 98may02, cca
   ---------------------------------------
*/
void
DenseMtx_setFields (
   DenseMtx   *mtx,
   int         type,
   int         rowid,
   int         colid,
   int         nrow,
   int         ncol,
   int         inc1,
   int         inc2
) {
double   *dbuffer ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || nrow < 0 || ncol < 0 
   || !((inc1 == ncol && inc2 == 1) || (inc1 == 1 && inc2 == nrow)) ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_setFields(%p,%d,%d,%d,%d,%d,%d)"
           "\n bad input\n", 
           mtx, rowid, colid, nrow, ncol, inc1, inc2) ;
   exit(-1) ;
}
dbuffer = DV_entries(&mtx->wrkDV) ;
ibuffer = (int *) dbuffer ;
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
mtx->type  = ibuffer[0] = type  ;
mtx->rowid = ibuffer[1] = rowid ;
mtx->colid = ibuffer[2] = colid ;
mtx->nrow  = ibuffer[3] = nrow  ;
mtx->ncol  = ibuffer[4] = ncol  ;
mtx->inc1  = ibuffer[5] = inc1  ;
mtx->inc2  = ibuffer[6] = inc2  ;
/*
   -------------------
   set up the pointers
   -------------------
*/
mtx->rowind = ibuffer + 7 ;
mtx->colind = mtx->rowind + nrow ;
if ( sizeof(int) == sizeof(double) ) {
   mtx->entries = dbuffer + 7 + nrow + ncol ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   mtx->entries = dbuffer + (8 + nrow + ncol)/2 ;
}
/*
fprintf(stdout, 
        "\n rowind - ibuffer = %d" 
        "\n colind - rowind  = %d" 
        "\n entries - dbuffer = %d", 
        mtx->rowind - ibuffer,
        mtx->colind - mtx->rowind,
        mtx->entries - dbuffer) ;
*/

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   purpose -- basic initializer

   created -- 98may02, cca
   ----------------------------
*/
void
DenseMtx_init (
   DenseMtx   *mtx,
   int        type,
   int        rowid,
   int        colid,
   int        nrow,
   int        ncol,
   int        inc1,
   int        inc2
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || nrow < 0 || ncol < 0 
   || !((inc1 == ncol && inc2 == 1) || (inc1 == 1 && inc2 == nrow)) ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_init(%p,%d,%d,%d,%d,%d,%d)"
           "\n bad input\n", 
           mtx, rowid, colid, nrow, ncol, inc1, inc2) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in DenseMtx_init(%p,%d,%d,%d,%d,%d,%d,%d)"
           "\n bad type %d\n", 
           mtx, type, rowid, colid, nrow, ncol, inc1, inc2, type) ;
   exit(-1) ;
   break ;
}
/*
   -------------------------------------------------------
   get and set the number of bytes needed in the workspace
   -------------------------------------------------------
*/
nbytes = DenseMtx_nbytesNeeded(type, nrow, ncol) ;
DenseMtx_setNbytesInWorkspace(mtx, nbytes) ;
/*
   --------------
   set the fields
   --------------
*/
DenseMtx_setFields(mtx, type, rowid, colid, nrow, ncol, inc1, inc2) ;
if ( nrow > 0 ) {
   IVramp(nrow, mtx->rowind, 0, 1) ;
}
if ( ncol > 0 ) {
   IVramp(ncol, mtx->colind, 0, 1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- initialize the object from its working storage
              used when the object is a MPI message

   created -- 98may02, cca
   ---------------------------------------------------------
*/
void
DenseMtx_initFromBuffer (
   DenseMtx   *mtx
) {
int   *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_initFromBuffer(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
ibuffer   = (int *) DV_entries(&mtx->wrkDV) ;
DenseMtx_setFields(mtx, ibuffer[0], ibuffer[1], ibuffer[2], 
                   ibuffer[3], ibuffer[4], ibuffer[5], ibuffer[6]) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- initializer with pointers

   created -- 98may02, cca
   ------------------------------------
*/
void
DenseMtx_initWithPointers (
   DenseMtx   *mtx,
   int         type,
   int         rowid,
   int         colid,
   int         nrow,
   int         ncol,
   int         inc1,
   int         inc2,
   int         *rowind,
   int         *colind,
   double      *entries
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || nrow <= 0 || ncol <= 0 || inc1 < 0 || inc2 < 0
   || (inc1 != 1 && inc2 != 1) 
   || entries == NULL || colind == NULL || rowind == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_initWithPointers()"
           "\n mtx = %p, rowid = %d, colid = %d"
           "\n nrow = %d, ncol = %d, inc1 = %d, inc2 = %d"
           "\n rowind = %p, colind = %p, entries = %p "
           "\n bad input\n",
           mtx, rowid, colid, nrow, ncol, inc1, inc2, 
           rowind, colind, entries) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n fatal error in DenseMtx_initWithPointers()"
           "\n bad type %d\n", type) ;
   break ;
}

/*
   ---------------------
   set the scalar fields
   ---------------------
*/
mtx->type  = type  ;
mtx->rowid = rowid ;
mtx->colid = colid ;
mtx->nrow  = nrow  ;
mtx->ncol  = ncol  ;
mtx->inc1  = inc1  ;
mtx->inc2  = inc2  ;
/*
   --------------------------
   set up the working storage
   --------------------------
*/
mtx->rowind  = rowind  ;
mtx->colind  = colind  ;
mtx->entries = entries ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   this method initializes a A2 object 
   to point into the entries

   created -- 98may02, cca
   -----------------------------------
*/
void
DenseMtx_setA2 (
   DenseMtx   *mtx,
   A2         *a2
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || a2 == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setZA2(%p,%p)"
           "\n bad input\n", mtx, a2) ;
   exit(-1) ;
}
/*
if ( DENSEMTX_IS_REAL(mtx) ) {
   A2_init(a2, SPOOLES_REAL, mtx->nrow, mtx->ncol, mtx->inc1, mtx->inc2,
           mtx->entries) ;
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   A2_init(a2, SPOOLES_COMPLEX, mtx->nrow, mtx->ncol, mtx->inc1, mtx->inc2,
           mtx->entries) ;
}
*/
A2_init(a2, mtx->type, mtx->nrow, mtx->ncol, 
        mtx->inc1, mtx->inc2, mtx->entries) ;
return ; }

/*--------------------------------------------------------------------*/
