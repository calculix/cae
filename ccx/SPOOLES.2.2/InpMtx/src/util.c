/*  util.c  */

#include "../InpMtx.h"

#define MYDEBUG  0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   given the data is in raw triples,
   sort and compress the data

   created  -- 98jan28, cca
   modified -- 98sep04, cca
      test to see if the sort is necessary 
   ---------------------------------------
*/
void
InpMtx_sortAndCompress (
   InpMtx   *inpmtx
) {
int      ient, nent, sortMustBeDone ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_sortAndCompress(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if (  INPMTX_IS_SORTED(inpmtx) 
   || INPMTX_IS_BY_VECTORS(inpmtx) 
   || (nent = inpmtx->nent) == 0 ) {
   inpmtx->storageMode = INPMTX_SORTED ;
   return ;
}
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
sortMustBeDone = 0 ;
for ( ient = 1 ; ient < nent ; ient++ ) {
   if ( ivec1[ient-1] > ivec1[ient] 
      || (   ivec1[ient-1] == ivec1[ient] 
          && ivec2[ient-1] > ivec2[ient] ) ) {
      sortMustBeDone = 1 ;
      break ;
   }
}
if ( sortMustBeDone == 1 ) {
   if ( INPMTX_IS_INDICES_ONLY(inpmtx) ) {
      inpmtx->nent = IV2sortUpAndCompress(nent, ivec1, ivec2) ;
   } else if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
      double   *dvec = InpMtx_dvec(inpmtx) ;
      inpmtx->nent = IV2DVsortUpAndCompress(nent, ivec1, ivec2, dvec) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
      double   *dvec = InpMtx_dvec(inpmtx) ;
      inpmtx->nent = IV2ZVsortUpAndCompress(nent, ivec1, ivec2, dvec) ;
   }
}
inpmtx->storageMode = INPMTX_SORTED ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   convert from sorted and compressed triples to vector

   created -- 98jan28, cca
   ----------------------------------------------------
*/
void
InpMtx_convertToVectors (
   InpMtx   *inpmtx 
) {
int      *ivec1, *ivec2, *offsets, *sizes, *vecids ;
int      first, id, ient, jj, nent, nvector, value ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_convertToVectors(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( INPMTX_IS_BY_VECTORS(inpmtx) || (nent = inpmtx->nent) == 0 ) {
   inpmtx->storageMode = INPMTX_BY_VECTORS ;
   return ;
}
if ( INPMTX_IS_RAW_DATA(inpmtx) ) {
   InpMtx_sortAndCompress(inpmtx) ;
}
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
/*
   -----------------------------------
   find the number of distinct vectors
   -----------------------------------
*/
value = -1 ;
nvector = 0 ;
for ( ient = 0 ; ient < nent ; ient++ ) {
   if ( ivec1[ient] >= 0 && value != ivec1[ient] ) {
      value = ivec1[ient] ;
      nvector++ ;
#if MYDEBUG > 0
      fprintf(stdout, "\n ient %d, value %d, nvector %d", 
              ient, value, nvector) ;
      fflush(stdout) ;
#endif
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n %d vectors", nvector) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------------------------------------
   adjust the sizes of the sizes[] and offsets[] arrays if necessary
   -----------------------------------------------------------------
*/
InpMtx_setNvector(inpmtx, nvector) ;
if ( nvector <= 0 ) {
/*
   -----------------------------
   matrix has no entries, return
   -----------------------------
*/
   inpmtx->storageMode = INPMTX_RAW_DATA ;
   InpMtx_setNent(inpmtx, 0) ;
   return ;
}
vecids  = InpMtx_vecids(inpmtx)  ;
sizes   = InpMtx_sizes(inpmtx)   ;
offsets = InpMtx_offsets(inpmtx) ;
/*
   ------------------------------------------------------------
   set the vector sizes and offsets
   note: we skip all entries whose first coordinate is negative
   ------------------------------------------------------------
*/
for ( first = 0 ; first < nent ; first++ ) {
   if ( ivec1[first] >= 0 ) {
      break ;
   }
}
id = 0 ;
if ( first < nent ) {
   value = ivec1[first] ;
   for ( jj = first + 1 ; jj < nent ; jj++ ) {
      if ( ivec1[jj] != value ) {
         vecids[id]  = value      ;
         sizes[id]   = jj - first ;
         offsets[id] = first      ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n vecids[%d] = %d, sizes[%d] = %d, offsets[%d] = %d",
                 id, vecids[id], id, sizes[id], id, offsets[id]) ;
         fflush(stdout) ;
#endif
         first       = jj         ;
         value       = ivec1[jj]  ;
         id++ ;
      }
   }
   vecids[id]  = value      ;
   sizes[id]   = jj - first ;
   offsets[id] = first      ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n vecids[%d] = %d, sizes[%d] = %d, offsets[%d] = %d",
           id, vecids[id], id, sizes[id], id, offsets[id]) ;
   fflush(stdout) ;
#endif
}
inpmtx->storageMode = INPMTX_BY_VECTORS ;
#if MYDEBUG > 0
fprintf(stdout, "\n vecidsIV") ;
IV_writeForHumanEye(&inpmtx->vecidsIV, stdout) ;
fprintf(stdout, "\n sizesIV") ;
IV_writeForHumanEye(&inpmtx->sizesIV, stdout) ;
fprintf(stdout, "\n offsetsIV") ;
IV_writeForHumanEye(&inpmtx->offsetsIV, stdout) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   drop off-diagonal entries

   created -- 98jan28, cca
   -------------------------
*/
void
InpMtx_dropOffdiagonalEntries (
   InpMtx   *inpmtx
) {
double   *dvec ;
int      count, ii, nent ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_dropOffdiagonalEntries(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_dropOffdiagonalEntries(%p)"
           "\n bad coordType = %d\n", inpmtx, inpmtx->coordType) ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
count = 0 ;
if (  INPMTX_IS_REAL_ENTRIES(inpmtx)
   || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   dvec = InpMtx_dvec(inpmtx) ;
}
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec1[ii] == ivec2[ii] ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count]   = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec1[ii] == ivec2[ii] ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count]   = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec2[ii] == 0 ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count]   = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
}
inpmtx->nent = count ;
IV_setSize(&inpmtx->ivec1IV, count) ;
IV_setSize(&inpmtx->ivec2IV, count) ;
if (  INPMTX_IS_REAL_ENTRIES(inpmtx)
   || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   DV_setSize(&inpmtx->dvecDV, count) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   drop entries in the lower triangle

   created -- 98jan28, cca
   ----------------------------------
*/
void
InpMtx_dropLowerTriangle (
   InpMtx   *inpmtx
) {
double   *dvec ;
int      count, ii, nent ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_dropLowerTriangle(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_dropLowerTriangle(%p)"
           "\n bad coordType \n", inpmtx) ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
count = 0 ;
if (  INPMTX_IS_REAL_ENTRIES(inpmtx)
   || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   dvec = InpMtx_dvec(inpmtx) ;
}
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
/*
fprintf(stdout, "\n ii = %d, (%d,%d)", ii, ivec1[ii], ivec2[ii]) ;
*/
      if ( ivec1[ii] <= ivec2[ii] ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
/*
fprintf(stdout, "\n    ivec1[%d] = %d, ivec2[%d] = %d",
        count, ivec1[count], count, ivec2[count]) ;
*/
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count] = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec1[ii] >= ivec2[ii] ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count] = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec2[ii] >= 0 ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count] = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
}
inpmtx->nent = count ;
IV_setSize(&inpmtx->ivec1IV, count) ;
IV_setSize(&inpmtx->ivec2IV, count) ;
if (  INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   DV_setSize(&inpmtx->dvecDV, count) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   DV_setSize(&inpmtx->dvecDV, 2*count) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   drop entries in the upper triangle

   created -- 98jan28, cca
   ----------------------------------
*/
void
InpMtx_dropUpperTriangle (
   InpMtx   *inpmtx
) {
double   *dvec ;
int      count, ii, nent ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_dropUpperTriangle(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_dropUpperTriangle(%p)"
           "\n bad coordType \n", inpmtx) ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
count = 0 ;
if (  INPMTX_IS_REAL_ENTRIES(inpmtx)
   || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   dvec = InpMtx_dvec(inpmtx) ;
}
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec1[ii] >= ivec2[ii] ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count] = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec1[ii] <= ivec2[ii] ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count] = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec2[ii] <= 0 ) {
         ivec1[count] = ivec1[ii] ;
         ivec2[count] = ivec2[ii] ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
            dvec[count] = dvec[ii]   ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
            dvec[2*count]   = dvec[2*ii]   ;
            dvec[2*count+1] = dvec[2*ii+1] ;
         }
         count++ ;
      }
   }
}
inpmtx->nent = count ;
IV_setSize(&inpmtx->ivec1IV, count) ;
IV_setSize(&inpmtx->ivec2IV, count) ;
if (  INPMTX_IS_REAL_ENTRIES(inpmtx)
   || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   DV_setSize(&inpmtx->dvecDV, count) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   map entries into the lower triangle

   created -- 98jan28, cca
   -----------------------------------
*/
void
InpMtx_mapToLowerTriangle (
   InpMtx   *inpmtx
) {
int      col, ii, nent, row ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToLowerTriangle(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToLowerTriangle(%p)"
           "\n bad coordType\n", inpmtx) ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( (row = ivec1[ii]) < (col = ivec2[ii]) ) {
         ivec1[ii] = col ;
         ivec2[ii] = row ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( (row = ivec2[ii]) < (col = ivec1[ii]) ) {
         ivec1[ii] = row ;
         ivec2[ii] = col ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec2[ii] > 0 ) {
         ivec2[ii] = -ivec2[ii] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   map entries into the upper triangle

   created -- 98jan28, cca
   -----------------------------------
*/
void
InpMtx_mapToUpperTriangle (
   InpMtx   *inpmtx
) {
int      col, ii, nent, row ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToUpperTriangle(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToUpperTriangle(%p)"
           "\n bad coordType = %d, must be 1, 2, or 3\n", 
           inpmtx, inpmtx->coordType) ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( (row = ivec1[ii]) > (col = ivec2[ii]) ) {
         ivec1[ii] = col ;
         ivec2[ii] = row ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( (row = ivec2[ii]) > (col = ivec1[ii]) ) {
         ivec1[ii] = row ;
         ivec2[ii] = col ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec2[ii] < 0 ) {
         ivec2[ii] = -ivec2[ii] ;
      }
   }
}
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   map entries into the upper triangle
   for a hermitian matrix

   created -- 98jan28, cca
   -----------------------------------
*/
void
InpMtx_mapToUpperTriangleH (
   InpMtx   *inpmtx
) {
double   *dvec ;
int      col, ii, nent, row ;
int      *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToUpperTriangle(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToUpperTriangle(%p)"
           "\n bad coordType = %d, must be 1, 2, or 3\n", 
           inpmtx, inpmtx->coordType) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapToUpperTriangleH(%p)"
           "\n input mode is not complex\n", inpmtx) ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
dvec  = InpMtx_dvec(inpmtx) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( (row = ivec1[ii]) > (col = ivec2[ii]) ) {
         ivec1[ii]    = col ;
         ivec2[ii]    = row ;
         dvec[2*ii+1] = -dvec[2*ii+1] ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( (row = ivec2[ii]) > (col = ivec1[ii]) ) {
         ivec1[ii] = row ;
         ivec2[ii] = col ;
         dvec[2*ii+1] = -dvec[2*ii+1] ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( ivec2[ii] < 0 ) {
         ivec2[ii] = -ivec2[ii] ;
         dvec[2*ii+1] = -dvec[2*ii+1] ;
      }
   }
}
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- compute the checksums of the indices and entries

      sums[0] = sum_{ii=0}^{nent} abs(ivec1[ii])
      sums[1] = sum_{ii=0}^{nent} abs(ivec2[ii])
      if real or complex entries then
         sums[2] = sum_{ii=0}^{nent} magnitudes of entries
      endif

   created -- 98may16, cca
   -----------------------------------------------------------
*/
void
InpMtx_checksums (
   InpMtx   *inpmtx,
   double   sums[]
) {
int   ient, nent ;
int   *ivec1, *ivec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_checksums(%p,%p)"
           "\n bad input\n", inpmtx, sums) ;
   exit(-1) ;
}
switch ( inpmtx->inputMode ) {
case INPMTX_INDICES_ONLY :
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_checksums(%p,%p)"
           "\n bad inputMode\n", inpmtx, sums) ;
   exit(-1) ;
}
sums[0] = sums[1] = sums[2] = 0.0 ;
if ( (nent = InpMtx_nent(inpmtx)) <= 0 ) {
   return ;
}
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
for ( ient = 0 ; ient < nent ; ient++ ) {
   sums[0] += abs(ivec1[ient]) ;
   sums[1] += abs(ivec2[ient]) ;
}
switch ( inpmtx->inputMode ) {
case INPMTX_INDICES_ONLY :
   break ;
case SPOOLES_REAL : {
   double *dvec  = InpMtx_dvec(inpmtx)  ;
   for ( ient = 0 ; ient < nent ; ient++ ) {
      sums[2] += fabs(dvec[ient]) ;
   }
   } break ;
case SPOOLES_COMPLEX : {
   double *dvec  = InpMtx_dvec(inpmtx)  ;
   for ( ient = 0 ; ient < nent ; ient++ ) {
      sums[2] += Zabs(dvec[2*ient], dvec[2*ient+1]) ;
   }
   } break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to create an InpMtx object filled with random entries

   input --

      mtx         -- matrix object, if NULL, it is created
      inputMode   -- input mode for the object,
                     indices only, real or complex entries
      coordType   -- coordinate type for the object,
                     by rows, by columns or by chevrons
      storageMode -- storage mode for the object,
                     raw data, sorted or by vectors
      nrow        -- # of rows
      ncol        -- # of columns
      symflag     -- symmetry flag for the matrix,
                     symmetric, hermitian or nonsymmetric
      nonzerodiag -- if 1, entries are placed on the diagonal
      nitem       -- # of items to be placed into the matrix
      seed        --  random number seed

   return value ---
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- bad input mode
     -3 -- bad coordinate type
     -4 -- bad storage mode
     -5 -- nrow or ncol <= 0
     -6 -- bad symmetry flag
     -7 -- hermitian matrix but not complex
     -8 -- symmetric or hermitian matrix but nrow != ncol
     -9 -- nitem < 0
   ----------------------------------------------------------------
*/
int
InpMtx_randomMatrix (
   InpMtx   *mtx,
   int      inputMode,
   int      coordType,
   int      storageMode,
   int      nrow,
   int      ncol,
   int      symflag,
   int      nonzerodiag,
   int      nitem,
   int      seed
) {
double   *dvec ;
Drand    *drand ;
int      col, ii, neqns, row ;
int      *colids, *rowids ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
switch ( inputMode ) {
case INPMTX_INDICES_ONLY :
case SPOOLES_REAL        :
case SPOOLES_COMPLEX     :
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n bad input mode %d\n", inputMode) ;
   return(-2) ;
   break ;
}
switch ( coordType ) {
case INPMTX_BY_ROWS     :
case INPMTX_BY_COLUMNS  :
case INPMTX_BY_CHEVRONS :
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n bad coordinate type %d\n", coordType) ;
   return(-3) ;
   break ;
}
switch ( storageMode ) {
case INPMTX_RAW_DATA   :
case INPMTX_SORTED     :
case INPMTX_BY_VECTORS :
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n bad storage mode%d\n", storageMode) ;
   return(-4) ;
   break ;
}
if ( nrow <= 0 || ncol <= 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n nrow = %d, ncol = %d\n", nrow, ncol) ;
   return(-5) ;
}
switch ( symflag ) {
case SPOOLES_SYMMETRIC    :
case SPOOLES_HERMITIAN    :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n bad symmetry flag%d\n", symflag) ;
   return(-6) ;
   break ;
}
if ( symflag == SPOOLES_HERMITIAN && inputMode != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n symmetryflag is Hermitian, requires complex type\n") ;
   return(-7) ;
}
if ( (symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN)
  && nrow != ncol ) {
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n symmetric or hermitian matrix, nrow %d, ncol%d\n",
           nrow, ncol) ;
   return(-8) ;
}
if ( nitem < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_randomMatrix"
           "\n nitem = %d\n", nitem) ;
   return(-9) ;
}
/*--------------------------------------------------------------------*/
neqns = (nrow <= ncol) ? nrow : ncol ;
if ( nonzerodiag == 1 ) {
   nitem += neqns ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
InpMtx_init(mtx, INPMTX_BY_ROWS, inputMode, nitem, 0) ;
/*
   ----------------
   fill the triples
   ----------------
*/
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
rowids = IVinit(nitem, -1) ;
colids = IVinit(nitem, -1) ;
if ( nonzerodiag == 1 ) {
   IVramp(neqns, rowids, 0, 1) ;
   Drand_setUniform(drand, 0, nrow) ;
   Drand_fillIvector(drand, nitem - neqns, rowids + neqns) ;
   Drand_setUniform(drand, 0, ncol) ;
   IVramp(neqns, colids, 0, 1) ;
   Drand_fillIvector(drand, nitem - neqns, colids + neqns) ;
} else {
   Drand_setUniform(drand, 0, nrow) ;
   Drand_fillIvector(drand, nitem, rowids) ;
   Drand_setUniform(drand, 0, ncol) ;
   Drand_fillIvector(drand, nitem, colids) ;
}
if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
   for ( ii = 0 ; ii < nitem ; ii++ ) {
      if ( (row = rowids[ii]) > (col = colids[ii]) ) {
         rowids[ii] = col ;
         colids[ii] = row ;
      }
   }
}
if ( inputMode == SPOOLES_REAL ) {
   dvec = DVinit(nitem, 0.0) ;
   Drand_setUniform(drand, 0.0, 1.0) ;
   Drand_fillDvector(drand, nitem, dvec) ;
} else if ( inputMode == SPOOLES_COMPLEX ) {
   dvec = DVinit(2*nitem, 0.0) ;
   Drand_setUniform(drand, 0.0, 1.0) ;
   Drand_fillDvector(drand, 2*nitem, dvec) ;
   if ( symflag == SPOOLES_HERMITIAN ) {
      for ( ii = 0 ; ii < nitem ; ii++ ) {
         if ( rowids[ii] == colids[ii] ) {
            dvec[2*ii+1] = 0.0 ;
         }
      }
   }
} else {
   dvec = NULL ;
}
/*
   ----------------
   load the triples
   ----------------
*/
switch ( inputMode ) {
case INPMTX_INDICES_ONLY :
   InpMtx_inputTriples(mtx, nitem, rowids, colids) ;
   break ;
case SPOOLES_REAL :
   InpMtx_inputRealTriples(mtx, nitem, rowids, colids, dvec) ;
   break ;
case SPOOLES_COMPLEX :
   InpMtx_inputComplexTriples(mtx, nitem, rowids, colids, dvec) ;
   break ;
}
/*
   ----------------------------------------
   set the coordinate type and storage mode
   ----------------------------------------
*/
InpMtx_changeCoordType(mtx, coordType) ;
InpMtx_changeStorageMode(mtx, storageMode) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Drand_free(drand) ;
IVfree(rowids) ;
IVfree(colids) ;
if ( dvec != NULL ) {
   DVfree(dvec) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   determine the range of the matrix, 
   i.e., the minimum and maximum rows and columns

   if pmincol != NULL then *pmincol = minimum column id
   if pmaxcol != NULL then *pmaxcol = maximum column id
   if pminrow != NULL then *pminrow = minimum row id
   if pmaxrow != NULL then *pmaxrow = maximum row id

   return value ---
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- no entries in the matrix
     -3 -- invalid coordinate type
   
   created -- 98oct15, cca
   ----------------------------------------------------
*/
int
InpMtx_range (
   InpMtx   *mtx,
   int      *pmincol,
   int      *pmaxcol,
   int      *pminrow,
   int      *pmaxrow
) {
int   maxcol, maxrow, mincol, minrow, nent ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_range()"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
if ( (nent = mtx->nent) <= 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_range()"
           "\n %d entries\n", nent) ;
   return(-2) ;
}
if ( INPMTX_IS_BY_ROWS(mtx) ) {
   int   *rowind = InpMtx_ivec1(mtx) ;
   int   *colind = InpMtx_ivec2(mtx) ;
   int   col, ii, row ;

   minrow = maxrow = rowind[0] ;
   mincol = maxcol = colind[0] ;
   for ( ii = 1 ; ii < nent ; ii++ ) {
      row = rowind[ii] ;
      col = colind[ii] ;
      if ( minrow > row ) {
         minrow = row ;
      } else if ( maxrow < row ) {
         maxrow = row ;
      }
      if ( mincol > col ) {
         mincol = col ;
      } else if ( maxcol < col ) {
         maxcol = col ;
      }
   }
 } else if ( INPMTX_IS_BY_COLUMNS(mtx) ) {
   int   *rowind = InpMtx_ivec2(mtx) ;
   int   *colind = InpMtx_ivec1(mtx) ;
   int   col, ii, row ;

   minrow = maxrow = rowind[0] ;
   mincol = maxcol = colind[0] ;
   for ( ii = 1 ; ii < nent ; ii++ ) {
      row = rowind[ii] ;
      col = colind[ii] ;
      if ( minrow > row ) {
         minrow = row ;
      } else if ( maxrow < row ) {
         maxrow = row ;
      }
      if ( mincol > col ) {
         mincol = col ;
      } else if ( maxcol < col ) {
         maxcol = col ;
      }
   }
 } else if ( INPMTX_IS_BY_CHEVRONS(mtx) ) {
   int   *chvind  = InpMtx_ivec1(mtx) ;
   int   *offsets = InpMtx_ivec2(mtx) ;
   int   col, ii, off, row ;

   if ( (off = offsets[0]) >= 0 ) {
      row = chvind[0], col = row + off ;
   } else {
      col = chvind[0], row = col - off ;
   }
   minrow = maxrow = row ;
   mincol = maxcol = col ;
   for ( ii = 1 ; ii < nent ; ii++ ) {
      if ( (off = offsets[ii]) >= 0 ) {
         row = chvind[ii], col = row + off ;
      } else {
         col = chvind[ii], row = col - off ;
      }
      if ( minrow > row ) {
         minrow = row ;
      } else if ( maxrow < row ) {
         maxrow = row ;
      }
      if ( mincol > col ) {
         mincol = col ;
      } else if ( maxcol < col ) {
         maxcol = col ;
      }
   }
} else {
   fprintf(stderr, "\n fatal error in InpMtx_range()"
           "\n invalid coordType %d\n", mtx->coordType) ;
   return(-3) ;
}
if ( pmincol != NULL ) { *pmincol = mincol ; }
if ( pmaxcol != NULL ) { *pmaxcol = maxcol ; }
if ( pminrow != NULL ) { *pminrow = minrow ; }
if ( pmaxrow != NULL ) { *pmaxrow = maxrow ; }

return(1) ; }

/*--------------------------------------------------------------------*/
