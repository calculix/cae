/*  instanceRead.c  */
 
#include "../InpMtx.h"
 
/*--------------------------------------------------------------------*/
/*
   -----------------------
   returns coordinate type

   created -- 98jan28, cca
   -----------------------
*/
int   
InpMtx_coordType (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_coordType(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->coordType) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   returns storage mode

   created -- 98jan28, cca
    ----------------------
*/
int   
InpMtx_storageMode (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_storageMode(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->storageMode) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   returns input mode

   created -- 98jan28, cca
   -----------------------
*/
int   
InpMtx_inputMode (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_inputMode(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->inputMode) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   returns inpmtx->mxnent

   created -- 98jan28, cca
   -----------------------
*/
int   
InpMtx_maxnent (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_maxnent(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->maxnent) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   returns inpmtx->nent

   created -- 98jan28, cca
   -----------------------
*/
int   
InpMtx_nent (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nent(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->nent) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   returns inpmtx->mxnvector

   created -- 98jan28, cca
   -------------------------
*/
int   
InpMtx_maxnvector (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_maxnvector(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->maxnvector) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   returns inpmtx->nvector

   created -- 98jan28, cca
   -----------------------
*/
int   
InpMtx_nvector (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nvector(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->nvector) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   returns inpmtx->resizeMultiple

   created -- 98jan28, cca
   ------------------------------
*/
double   
InpMtx_resizeMultiple (
   InpMtx  *inpmtx
) {
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_resizeMultiple(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(inpmtx->resizeMultiple) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   returns pointer to ivec1[] vector

   created -- 98jan28, cca
   ---------------------------------
*/
int *   
InpMtx_ivec1 (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_ivec1(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(IV_entries(&inpmtx->ivec1IV)) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   returns pointer to ivec2[] vector

   created -- 98jan28, cca
   ---------------------------------
*/
int *   
InpMtx_ivec2 (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_ivec2(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(IV_entries(&inpmtx->ivec2IV)) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   returns pointer to dvec[] vector

   created -- 98jan28, cca
   --------------------------------
*/
double *   
InpMtx_dvec (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_dvec(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(DV_entries(&inpmtx->dvecDV)) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   returns pointer to sizes[] vector

   created -- 98jan28, cca
   ---------------------------------
*/
int *   
InpMtx_sizes (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_sizes(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(IV_entries(&inpmtx->sizesIV)) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   returns pointer to vecids[] vector

   created -- 98jan28, cca
   ----------------------------------
*/
int *   
InpMtx_vecids (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_vecids(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(IV_entries(&inpmtx->vecidsIV)) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   returns pointer to offsets[] vector

   created -- 98jan28, cca
   -----------------------------------
*/
int *   
InpMtx_offsets (
   InpMtx  *inpmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_offsets(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
return(IV_entries(&inpmtx->offsetsIV)) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   retrieve requested vector
   set *pnent to # of entries
       *pindices to address of first index

   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_vector (
   InpMtx   *inpmtx,
   int      id,
   int      *pnent,
   int      **pindices
) {
int   loc, off ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_vector(%p,%d,%p,%p)"
           "\n bad input\n", inpmtx, id, pnent, pindices) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_BY_VECTORS(inpmtx) ) {
   fprintf(stderr, "\n fatal error in InpMtx_vector(%p,%d,%p,%p)"
           "\n bad input\n", inpmtx, id, pnent, pindices) ;
   exit(-1) ;
}
if (   pnent == NULL || pindices == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_vector(%p,%d,%p,%p)"
           "\n NULL input, pnent = %p, pindices = %p",
	   inpmtx, id, pnent, pindices, pnent, pindices) ;
   exit(-1) ;
}
/*
   -------------------------------
   find the location of the vector
   -------------------------------
*/
loc = IV_findValueAscending(&inpmtx->vecidsIV, id) ;
if ( loc == -1 ) {
/*
   ------------------------------------------------------------
   vector is not present, set size to zero and pointers to NULL
   ------------------------------------------------------------
*/
   *pnent    =   0  ;
   *pindices = NULL ;
} else {
/*
  --------------------------------------------------------------
  fill *pnent with the number of entries in vector id.
  fill *pindices with the pointer to the first index in vector id.
  fill *pentries with the pointer to the first entry in vector id.
  --------------------------------------------------------------
*/
   *pnent    = inpmtx->sizesIV.vec[loc]   ;
   off       = inpmtx->offsetsIV.vec[loc] ;
   *pindices = inpmtx->ivec2IV.vec + off  ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   retrieve requested vector
   set *pnent to # of entries
       *pindices to address of first index
       *pentries to address of first entry

   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_realVector (
   InpMtx   *inpmtx,
   int      id,
   int      *pnent,
   int      **pindices,
   double   **pentries
) {
int   loc, off ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_realVector(%p,%d,%p,%p,%p)"
           "\n bad input\n", inpmtx, id, pnent, pindices, pentries) ;
   exit(-1) ;
}
if ( !INPMTX_IS_BY_VECTORS(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_realVector(%p,%d,%p,%p,%p)"
           "\n storageMode must be INPMTX_BY_VECTORS\n", 
           inpmtx, id, pnent, pindices, pentries) ;
   exit(-1) ;
}
if (   pnent == NULL || pindices == NULL || pentries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_realVector(%p,%d,%p,%p,%p)"
           "\n NULL input, pnent = %p, pindices = %p, pentries = %p",
	   inpmtx, id, pnent, pindices, pentries,
	   pnent, pindices, pentries) ;
   exit(-1) ;
}
/*
   -------------------------------
   find the location of the vector
   -------------------------------
*/
loc = IV_findValueAscending(&inpmtx->vecidsIV, id) ;
if ( loc == -1 ) {
/*
   ------------------------------------------------------------
   vector is not present, set size to zero and pointers to NULL
   ------------------------------------------------------------
*/
   *pnent    =   0  ;
   *pindices = NULL ;
   *pentries = NULL ;
} else {
/*
  --------------------------------------------------------------
  fill *pnent with the number of entries in vector id.
  fill *pindices with the pointer to the first index in vector id.
  fill *pentries with the pointer to the first entry in vector id.
  --------------------------------------------------------------
*/
   *pnent    = inpmtx->sizesIV.vec[loc]   ;
   off       = inpmtx->offsetsIV.vec[loc] ;
   *pindices = inpmtx->ivec2IV.vec + off  ;
   *pentries = inpmtx->dvecDV.vec  + off  ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   retrieve requested vector
   set *pnent to # of entries
       *pindices to address of first index
       *pentries to address of first entry

   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_complexVector (
   InpMtx   *inpmtx,
   int      id,
   int      *pnent,
   int      **pindices,
   double   **pentries
) {
int   loc, off ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_complexVector(%p,%d,%p,%p,%p)"
           "\n bad input\n", inpmtx, id, pnent, pindices, pentries) ;
   exit(-1) ;
}
if ( !INPMTX_IS_BY_VECTORS(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_complexVector(%p,%d,%p,%p,%p)"
           "\n storage mode muse be INPMTX_BY_VECTORS\n", 
           inpmtx, id, pnent, pindices, pentries) ;
   exit(-1) ;
}
if (   pnent == NULL || pindices == NULL || pentries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_complexVector(%p,%d,%p,%p,%p)"
           "\n NULL input, pnent = %p, pindices = %p, pentries = %p",
	   inpmtx, id, pnent, pindices, pentries,
	   pnent, pindices, pentries) ;
   exit(-1) ;
}
/*
   -------------------------------
   find the location of the vector
   -------------------------------
*/
loc = IV_findValueAscending(&inpmtx->vecidsIV, id) ;
if ( loc == -1 ) {
/*
   ------------------------------------------------------------
   vector is not present, set size to zero and pointers to NULL
   ------------------------------------------------------------
*/
   *pnent    =   0  ;
   *pindices = NULL ;
   *pentries = NULL ;
} else {
/*
  --------------------------------------------------------------
  fill *pnent with the number of entries in vector id.
  fill *pindices with the pointer to the first index in vector id.
  fill *pentries with the pointer to the first entry in vector id.
  --------------------------------------------------------------
*/
   *pnent    = inpmtx->sizesIV.vec[loc]   ;
   off       = inpmtx->offsetsIV.vec[loc] ;
   *pindices = inpmtx->ivec2IV.vec + off  ;
   *pentries = inpmtx->dvecDV.vec + 2*off ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   sets the maximum numnber of entries.  this methods resizes the
   ivec1[], ivece2[] and dvec[] vectors if newmaxnent != maxnent

   created -- 98jan28, cca
   --------------------------------------------------------------
*/
void   
InpMtx_setMaxnent (
   InpMtx  *inpmtx, 
   int      newmaxnent
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || newmaxnent < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_setMaxnent(%p, %d)"
           "\n bad input\n", inpmtx, newmaxnent) ;
   exit(-1) ;
}
if ( newmaxnent != inpmtx->maxnent ) {
  IV_setMaxsize(&(inpmtx->ivec1IV), newmaxnent) ;
  IV_setMaxsize(&(inpmtx->ivec2IV), newmaxnent) ;
  if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
    DV_setMaxsize(&(inpmtx->dvecDV), newmaxnent) ;
  } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
    DV_setMaxsize(&(inpmtx->dvecDV), 2*newmaxnent) ;
  }
}
inpmtx->maxnent = newmaxnent ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   set the present number of entries

   created -- 98jan28, cca
   --------------------------------
*/
void
InpMtx_setNent (
   InpMtx   *inpmtx,
   int       newnent
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || newnent < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_setNent(%p,%d)"
           "\n bad input\n", inpmtx, newnent) ;
   exit(-1) ;
}
if ( inpmtx->maxnent < newnent ) {
/*
   -------------------------------------------------------
   newnent requested is more than maxnent, set new maxnent 
   -------------------------------------------------------
*/
   InpMtx_setMaxnent(inpmtx, newnent) ;
}
inpmtx->nent = newnent ;
IV_setSize(&inpmtx->ivec1IV, newnent) ;
IV_setSize(&inpmtx->ivec2IV, newnent) ;
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
    DV_setSize(&(inpmtx->dvecDV), newnent) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
    DV_setSize(&inpmtx->dvecDV, 2*newnent) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   sets the maximum number of vectors. 
   if newmaxnent != maxnent then this methods resizes 
   the vecids[], sizes[] and offsets[] vectors 

   created  -- 98jan28, cca
   --------------------------------------------------
*/
void   
InpMtx_setMaxnvector (
   InpMtx  *inpmtx, 
   int      newmaxnvector
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || newmaxnvector < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_Maxnvector(%p, %d)"
           "\n bad input\n", inpmtx, newmaxnvector) ;
   exit(-1) ;
}
if ( newmaxnvector != inpmtx->maxnvector ) {
  IV_setMaxsize(&(inpmtx->vecidsIV),  newmaxnvector) ;
  IV_setMaxsize(&(inpmtx->sizesIV),   newmaxnvector) ;
  IV_setMaxsize(&(inpmtx->offsetsIV), newmaxnvector) ;
}
inpmtx->maxnvector = newmaxnvector ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   set the present number of vectors

   created -- 98jan28, cca
   ---------------------------------
*/
void   
InpMtx_setNvector (
   InpMtx   *inpmtx, 
   int       newnvector
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || newnvector < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_setNvector(%p, %d)"
           "\n bad input\n", inpmtx, newnvector) ;
   exit(-1) ;
}
if ( newnvector > inpmtx->maxnvector ) {
  InpMtx_setMaxnvector(inpmtx, newnvector) ;
}
inpmtx->nvector = newnvector ;
IV_setSize(&inpmtx->vecidsIV,  newnvector) ;
IV_setSize(&inpmtx->sizesIV,   newnvector) ;
IV_setSize(&inpmtx->offsetsIV, newnvector) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   sets inpmtx->resizeMultiple

   created -- 98jan28, cca
    ---------------------------
*/
void   
InpMtx_setResizeMultiple (
   InpMtx   *inpmtx, 
   double    resizeMultiple
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || resizeMultiple < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_setResizeMultiple(%p,%f)"
           "\n bad input\n", inpmtx, resizeMultiple) ;
   exit(-1) ;
}
inpmtx->resizeMultiple = resizeMultiple ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   sets coordType of InpMtx structure to
   allow user to define custom coordinate type.
   Note, new type must be > 3.

   created -- 98jan28, cca
    --------------------------------------------
*/
void   
InpMtx_setCoordType (
   InpMtx  *inpmtx, 
   int      type
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || type <= 3 ) {
   fprintf(stderr, "\n fatal error in InpMtx_setCoordType(%p,%d)"
           "\n bad input\n", inpmtx, type) ;
   if ( type <= 3 ) 
     fprintf(stderr, "\n fatal error in InpMtx_setCoordType"
           "\n reserved coordinate type %d \n", type) ;
      
   exit(-1) ;
}
inpmtx->coordType = type ;

return ; }    

/*--------------------------------------------------------------------*/
