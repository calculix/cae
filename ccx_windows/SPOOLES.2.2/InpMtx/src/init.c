/*  init.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   initialize the object

   coordType -- coordinate type, input supported for types 1, 2 and 3
      1 -- row triples (i, j, a_{i,j})
      2 -- column triples (i, j, a_{j,i})
      3 -- chevron triples 
           (i, k, a_{i,i+k}) if k >= 0
           (i, k, a_{i-k,i}) if k < 0
           i is the chevron, k is the offset
      4 -- custom coordinate, e.g., one could store (I, k, a_{i,j})
           where I is the front where a_{i,j} will be assembled
           and k is the offset into the vector that holds the 
           entries in the front
   inputMode -- mode for input
      0 --> indices only
      1 --> indices and real entries
      2 --> indices and complex entries
   maxnent  -- upper bound on the number of entries,
      also equal to the workspace, so if the assembly includes
      overlapping data, give enough elbow room for efficiency.
   maxnvector -- upper bound on the number of vectors to be supported. 
      this may not be known ahead of time (e.g., the number of vectors 
      may be the number of fronts which is not known before the
      ordering is done and front tree constructed).

   created -- 98jan28, cca
   ------------------------------------------------------------------
*/
void
InpMtx_init (
  InpMtx   *inpmtx,
  int      coordType,
  int      inputMode,
  int      maxnent,
  int      maxnvector
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_init(%p,%d,%d,%d,%d)"
           "\n inpmtx is NULL \n", 
           inpmtx, coordType, inputMode, maxnent, maxnvector) ;
   exit(-1) ;
}
if ( ! ( INPMTX_IS_BY_ROWS(inpmtx)
      || INPMTX_IS_BY_COLUMNS(inpmtx)
      || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_init(%p,%d,%d,%d,%d)"
           "\n bad coordType \n", 
           inpmtx, coordType, inputMode, maxnent, maxnvector) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_INDICES_ONLY(inpmtx)
   ||   INPMTX_IS_REAL_ENTRIES(inpmtx)
   ||   INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_init(%p,%d,%d,%d,%d)"
           "\n bad inputMode \n", 
           inpmtx, coordType, inputMode, maxnent, maxnvector) ;
   exit(-1) ;
}
if ( maxnent < 0 || maxnvector < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_init(%p,%d,%d,%d,%d)"
           "\n maxnent = %d, maxnvector = %d \n", 
           inpmtx, coordType, inputMode, maxnent, maxnvector,
           maxnent, maxnvector) ;
   exit(-1) ;
}
/*
   ------------------------
   clear the data structure
   ------------------------
*/
InpMtx_clearData(inpmtx) ;
/*
   -------------------------------------------------------------------
   if maxnent > 0 then allocate the ivec1[], ivec2[] and dvec[] vectors
   -------------------------------------------------------------------
*/
inpmtx->coordType = coordType ;
inpmtx->inputMode = inputMode ;
if ( maxnent > 0 ) {
   InpMtx_setMaxnent(inpmtx, maxnent) ;
}
/*
   ---------------------------------------
   if nvector > 0 initialize the vecids[], 
   sizes[] and offsets[] vectors.
   ---------------------------------------
*/
if ( maxnvector > 0 ) {
   InpMtx_setMaxnvector(inpmtx, maxnvector) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   change the coordinate type

   created -- 98jan28, cca
   --------------------------
*/
void
InpMtx_changeCoordType (
   InpMtx   *inpmtx,
   int       newType
) {
int   chev, col, i, nent, off, oldType, row, temp ;
int   *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_changeCoordType(%p,%d)"
           "\n bad input\n", inpmtx, newType) ;
   exit(-1) ;
}
if (  newType != INPMTX_BY_ROWS && newType != INPMTX_BY_COLUMNS
   && newType != INPMTX_BY_CHEVRONS && newType != INPMTX_CUSTOM ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_changeCoordType(%p,%d)"
           "\n bad new coordType\n", inpmtx, newType) ;
   exit(-1) ;
}
if ( ! ( INPMTX_IS_BY_ROWS(inpmtx)
      || INPMTX_IS_BY_COLUMNS(inpmtx)
      || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_changeCoordType(%p,%d)"
           "\n bad existing coordType\n", inpmtx, newType) ;
   exit(-1) ;
}
oldType = inpmtx->coordType ;
if ( oldType == newType ) {
   return ;
}
if ( newType == INPMTX_CUSTOM ) {
/*
   ----------------------------------------------------
   custom type, just set the coordType field and return
   ----------------------------------------------------
*/
   inpmtx->coordType = INPMTX_CUSTOM ;
   return ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( oldType == INPMTX_BY_ROWS ) {
   if ( newType == INPMTX_BY_COLUMNS ) {
/*
      ----------------------------------
      row coordType --> column coordType
      ----------------------------------
*/
      for ( i = 0 ; i < nent ; i++ ) {
         temp = ivec1[i] ;
         ivec1[i] = ivec2[i] ;
         ivec2[i] = temp ;
      }
      inpmtx->coordType   = INPMTX_BY_COLUMNS ;
      inpmtx->storageMode = INPMTX_RAW_DATA   ;
   } else if ( newType == INPMTX_BY_CHEVRONS ) {
/*
      -----------------------------------
      row coordType --> chevron coordType
      -----------------------------------
*/
      for ( i = 0 ; i < nent ; i++ ) {
         row = ivec1[i] ;
         col = ivec2[i] ;
         if ( row <= col ) {
            ivec1[i] = row ;
            ivec2[i] = col - row ;
         } else {
            ivec1[i] = col ;
            ivec2[i] = col - row ;
         }
      }
      inpmtx->coordType   = INPMTX_BY_CHEVRONS ;
      inpmtx->storageMode = INPMTX_RAW_DATA    ;
   }
} else if ( oldType == INPMTX_BY_COLUMNS ) {
   if ( newType == INPMTX_BY_ROWS ) {
/*
      ----------------------------------
      column coordType --> row coordType
      ----------------------------------
*/
      for ( i = 0 ; i < nent ; i++ ) {
         temp = ivec1[i] ;
         ivec1[i] = ivec2[i] ;
         ivec2[i] = temp ;
      }
      inpmtx->coordType   = INPMTX_BY_ROWS  ;
      inpmtx->storageMode = INPMTX_RAW_DATA ;
   } else if ( newType == INPMTX_BY_CHEVRONS ) {
/*
      --------------------------------------
      column coordType --> chevron coordType
      --------------------------------------
*/
      for ( i = 0 ; i < nent ; i++ ) {
         col = ivec1[i] ;
         row = ivec2[i] ;
         if ( row <= col ) {
            ivec1[i] = row ;
            ivec2[i] = col - row ;
         } else {
            ivec1[i] = col ;
            ivec2[i] = col - row ;
         }
      }
      inpmtx->coordType   = INPMTX_BY_CHEVRONS ;
      inpmtx->storageMode = INPMTX_RAW_DATA    ;
   }
} else if ( oldType == INPMTX_BY_CHEVRONS ) {
   if ( newType == INPMTX_BY_ROWS ) {
/*
      -----------------------------------
      chevron coordType --> row coordType
      -----------------------------------
*/
      for ( i = 0 ; i < nent ; i++ ) {
         chev = ivec1[i] ;
         off  = ivec2[i] ;
         if ( off >= 0 ) {
            row = chev ;
            col = chev + off ;
         } else {
            col = chev ;
            row = chev - off ;
         }
         ivec1[i] = row ;
         ivec2[i] = col ;
      }
      inpmtx->coordType   = INPMTX_BY_ROWS  ;
      inpmtx->storageMode = INPMTX_RAW_DATA ;
   } else if ( newType == INPMTX_BY_COLUMNS ) {
/*
      --------------------------------------
      chevron coordType --> column coordType
      --------------------------------------
*/
      for ( i = 0 ; i < nent ; i++ ) {
         chev = ivec1[i] ;
         off  = ivec2[i] ;
         if ( off >= 0 ) {
            row = chev ;
            col = chev + off ;
         } else {
            col = chev ;
            row = chev - off ;
         }
         ivec1[i] = col ;
         ivec2[i] = row ;
      }
      inpmtx->coordType   = INPMTX_BY_COLUMNS ;
      inpmtx->storageMode = INPMTX_RAW_DATA   ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   change the storage mode

   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_changeStorageMode (
   InpMtx   *inpmtx,
   int       newMode
) {
int   oldMode ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_changeStorageMode(%p,%d)"
           "\n inpmtx is NULL\n", inpmtx, newMode) ;
   exit(-1) ;
}
if (  newMode != INPMTX_RAW_DATA && newMode != INPMTX_SORTED
   && newMode != INPMTX_BY_VECTORS ) {
   fprintf(stderr, "\n fatal error in InpMtx_changeStorageMode(%p,%d)"
           "\n bad new storage mode\n", inpmtx, newMode) ;
   exit(-1) ;
}
oldMode = inpmtx->storageMode ;
if ( oldMode == newMode ) {
   return ;
}
if ( oldMode == INPMTX_RAW_DATA ) {
   if ( newMode == INPMTX_SORTED ) {
/*
      -------------------------------------
      raw triples --> sorted and compressed 
      -------------------------------------
*/
      InpMtx_sortAndCompress(inpmtx) ;
      inpmtx->storageMode = INPMTX_SORTED ;
   } else if ( newMode == INPMTX_BY_VECTORS ) {
/*
      -----------------------
      raw triples --> vectors
      -----------------------
*/
      InpMtx_sortAndCompress(inpmtx) ;
      InpMtx_convertToVectors(inpmtx) ;
      inpmtx->storageMode = INPMTX_BY_VECTORS ;
   }
} else if ( oldMode == INPMTX_SORTED ) {
   if ( newMode == INPMTX_RAW_DATA ) {
/*
      --------------------------------------
      sorted and compressed --> raw triples, 
      simply set the storageMode field
      --------------------------------------
*/
      inpmtx->storageMode = INPMTX_RAW_DATA ;
   } else if ( newMode == INPMTX_BY_VECTORS ) {
/*
      ---------------------------------
      sorted and compressed --> vectors
      ---------------------------------
*/
      InpMtx_convertToVectors(inpmtx) ;
      inpmtx->storageMode = INPMTX_BY_VECTORS ;
   }
} else if ( oldMode == INPMTX_BY_VECTORS ) {
   if ( newMode == INPMTX_RAW_DATA || newMode == INPMTX_SORTED ) {
/*
      --------------------------------
      vectors --> triples, 
      simply set the storageMode field
      --------------------------------
*/
      inpmtx->storageMode = newMode ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
