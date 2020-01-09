/*  init.c  */

#include "../SubMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   return the number of bytes needed for an object of this size
 
   created -- 98may01, cca
   ------------------------------------------------------------
*/
int
SubMtx_nbytesNeeded (
   int   type,
   int   mode,
   int   nrow,
   int   ncol,
   int   nent
) {
int   nbytes, ndouble, nint ;
/*
   ---------------
   check the input
   ---------------
*/
if (  nrow <= 0 || ncol <= 0 || nent < 0 ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_nbytesNeeded(%d,%d,%d,%d,%d)"
           "\n bad input\n", type, mode, nrow, ncol, nent) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in SubMtx_nbytesNeeded(%d,%d,%d,%d,%d)"
           "\n bad type\n", type, mode, nrow, ncol, nent) ;
   exit(-1) ;
}
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
case SUBMTX_SPARSE_TRIPLES :
case SUBMTX_DENSE_SUBROWS :
case SUBMTX_DENSE_SUBCOLUMNS :
case SUBMTX_DIAGONAL :
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in SubMtx_nbytesNeeded(%d,%d,%d,%d,%d)"
           "\n bad mode\n", type, mode, nrow, ncol, nent) ;
   exit(-1) ;
}
/*
   --------------------------------------------
   determine the number of integers and doubles 
   --------------------------------------------
*/
nint    = 7 + nrow + ncol ;
if ( type == SPOOLES_REAL ) {
   ndouble = nent ;
} else if ( type == SPOOLES_COMPLEX ) {
   ndouble = 2*nent ;
}
switch ( mode ) {
case SUBMTX_SPARSE_ROWS :
   nint += nrow + nent ;
   break ;
case SUBMTX_SPARSE_COLUMNS :
   nint += ncol + nent ;
   break ;
case SUBMTX_SPARSE_TRIPLES :
   nint += 2*nent ;
   break ;
case SUBMTX_DENSE_SUBROWS :
   nint += 2*nrow ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   nint += 2*ncol ;
   break ;
case SUBMTX_DIAGONAL :
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   nint += ncol ;
   break ;
default :
   break ;
}
/*
   -----------------------------
   determine the number of bytes
   -----------------------------
*/
if ( sizeof(int) == sizeof(double) ) {
   nbytes = nint*sizeof(int) + ndouble*sizeof(double) ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   nbytes = ((nint + 1)/2 + ndouble)*sizeof(double) ;
} else {
   fprintf(stderr, "\n error in SubMtx_nbytesNeeded(%d,%d)"
           "\n sizeof(int) = %d, sizeof(double) = %d",
           nrow, ncol, sizeof(int), sizeof(double)) ;
   exit(-1) ;
}
return(nbytes) ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   return the number of bytes in use in
   the workspace owned by this object
 
   created -- 98may01, cca
   ------------------------------------
*/
int
SubMtx_nbytesInUse (
   SubMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_nbytesInUse(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(sizeof(double)*DV_size(&mtx->wrkDV)) ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   return a pointer to the workspace owned by this object
 
   created -- 98may01, cca
   ------------------------------------------------------
*/
void *
SubMtx_workspace (
   SubMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_workspace(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return((void *)DV_entries(&mtx->wrkDV)) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   return the number of bytes in the workspace owned by this object
 
   created -- 98may01, cca
   ----------------------------------------------------------------
*/
int
SubMtx_nbytesInWorkspace (
   SubMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_nbytesInWorkspace(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(sizeof(double)*DV_maxsize(&mtx->wrkDV)) ; }
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   set the number of bytes in the workspace owned by this object
 
   created -- 98may01, cca
   -------------------------------------------------------------
*/
void
SubMtx_setNbytesInWorkspace (
   SubMtx   *mtx,
   int    nbytes
) {
if ( mtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_setNbytesInWorkspace(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
DV_setSize(&mtx->wrkDV, nbytes/sizeof(double)) ;

return ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- set the fields of the object

   created -- 98may01, cca
   ---------------------------------------
*/
void
SubMtx_setFields (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_setFields()"
           "\n mtx is NULL\n") ;
   exit(-1) ;
}
if (  nrow <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_setFields()"
           "\n nrow = %d <= 0\n", nrow) ;
   exit(-1) ;
}
if (  ncol <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_setFields()"
           "\n ncol = %d <= 0\n", ncol) ;
   exit(-1) ;
}
if (  nrow <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_setFields()"
           "\n nent = %d <= 0\n", nent) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_setFields()"
           "\n invalid type %d", type) ;
   exit(-1) ;
}
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_DIAGONAL :
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
case SUBMTX_SPARSE_TRIPLES :
case SUBMTX_DENSE_SUBROWS :
case SUBMTX_DENSE_SUBCOLUMNS :
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_setFields()"
           "\n invalid mode %d", mode) ;
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
mtx->mode  = ibuffer[1] = mode  ;
mtx->rowid = ibuffer[2] = rowid ;
mtx->colid = ibuffer[3] = colid ;
mtx->nrow  = ibuffer[4] = nrow  ;
mtx->ncol  = ibuffer[5] = ncol  ;
mtx->nent  = ibuffer[6] = nent  ;
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_DIAGONAL :
   nint = 7 + mtx->nrow + mtx->ncol ;
   break ;
case SUBMTX_SPARSE_ROWS :
   nint = 7 + mtx->nrow + mtx->ncol + mtx->nrow + mtx->nent ;
   break ;
case SUBMTX_SPARSE_COLUMNS :
   nint = 7 + mtx->nrow + mtx->ncol + mtx->ncol + mtx->nent ;
   break ;
case SUBMTX_SPARSE_TRIPLES :
   nint = 7 + mtx->nrow + mtx->ncol + mtx->nent + mtx->nent ;
   break ;
case SUBMTX_DENSE_SUBROWS :
   nint = 7 + mtx->nrow + mtx->ncol + mtx->nrow + mtx->nrow ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   nint = 7 + mtx->nrow + mtx->ncol + mtx->ncol + mtx->ncol ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   nint = 7 + mtx->nrow + mtx->ncol + mtx->nrow ;
   break ;
}
if ( sizeof(int) == sizeof(double) ) {
   mtx->entries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   mtx->entries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   purpose -- basic initializer

   created -- 98may01, cca
   ----------------------------
*/
void
SubMtx_init (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent
) {
int   nbytes ;
int   *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_init()"
           "\n mtx is NULL\n") ;
   exit(-1) ;
}
if (  nrow <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_init()"
           "\n nrow = %d <= 0\n", nrow) ;
   exit(-1) ;
}
if (  ncol <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_init()"
           "\n ncol = %d <= 0\n", ncol) ;
   exit(-1) ;
}
if (  nrow <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_init()"
           "\n nent = %d <= 0\n", nent) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_init()"
           "\n invalid type %d", type) ;
   exit(-1) ;
}
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_DIAGONAL :
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
case SUBMTX_SPARSE_TRIPLES :
case SUBMTX_DENSE_SUBROWS :
case SUBMTX_DENSE_SUBCOLUMNS :
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_init()"
           "\n invalid mode %d", mode) ;
   exit(-1) ;
}
/*
   -------------------------------------------------------
   get and set the number of bytes needed in the workspace
   -------------------------------------------------------
*/
nbytes = SubMtx_nbytesNeeded(type, mode, nrow, ncol, nent) ;
SubMtx_setNbytesInWorkspace(mtx, nbytes) ;
DVzero(nbytes/sizeof(double), (double *) SubMtx_workspace(mtx)) ;
/*
   --------------
   set the fields
   --------------
*/
SubMtx_setFields(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
IVramp(nrow, rowind, 0, 1) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
IVramp(ncol, colind, 0, 1) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- initialize the object from its working storage
              used when the object is a MPI message

   created -- 98may01, cca
   ---------------------------------------------------------
*/
void
SubMtx_initFromBuffer (
   SubMtx   *mtx
) {
int   *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_initFromBuffer(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
ibuffer   = (int *) DV_entries(&mtx->wrkDV) ;
SubMtx_setFields(mtx, ibuffer[0], ibuffer[1], ibuffer[2], 
                 ibuffer[3], ibuffer[4], ibuffer[5], ibuffer[6]) ;

return ; }

/*--------------------------------------------------------------------*/
