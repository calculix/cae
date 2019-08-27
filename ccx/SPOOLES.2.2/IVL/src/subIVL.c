/*  subIVL.c  */

#include "../IVL.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- initialize subIVL from ivl
     if keeplistIV is not NULL then
        keeplistIV contains the lists of ivl to be placed in subIVL
     else
        all lists are placed in subIVL
     endif
     if keepentriesIV is not NULL then
        keepentriesIV contains the entries in the lists to be kept
     else
        all entries in the kept lists are placed in subIVL
     endif

   return value ---
      1 -- normal return
     -1 -- subIVL is NULL
     -2 -- ivl is NULL
     -3 -- keeplistIV is invalid

   created -- 98oct16, cca
   ----------------------------------------------------------------
*/
int
IVL_initFromSubIVL (
   IVL   *subIVL,
   IVL   *ivl,
   IV    *keeplistIV,
   IV    *keepentriesIV
) {
int   count, ii, ikeep, ilist, maxlistsize, maxval, nkeep, nkeepent,
      nlist, size, val ;
int   *keepent, *keeplist, *list, *map, *temp ;
/*
   ---------------
   check the input
   ---------------
*/
if ( subIVL == NULL ) {
   fprintf(stderr, "\n error in IVL_initFromSubIVL()"
           "\n subIVL is NULL\n") ;
   return(-1) ;
}
if ( ivl == NULL ) {
   fprintf(stderr, "\n error in IVL_initFromSubIVL()"
           "\n ivl is NULL\n") ;
   return(-2) ;
}
nlist = ivl->nlist ;
if ( keeplistIV != NULL ) {
   IV_sizeAndEntries(keeplistIV, &nkeep, &keeplist) ;
   if ( nkeep < 0 || keeplist == NULL ) {
      fprintf(stderr, "\n error in IVL_initFromSubIVL()"
              "\n invalid keeplistIV, nkeep %d, keeplist %p\n",
              nkeep, keeplist) ;
      return(-3) ;
   }
   for ( ii = 0 ; ii < nkeep ; ii++ ) {
      if ( (val = keeplist[ii]) < 0 || val >= nlist ) {
         fprintf(stderr, "\n error in IVL_initFromSubIVL()"
                 "\n invalid keeplistIV, keeplist[%d] = %d, nlist %d\n",
                 ii, val, nlist) ;
         return(-3) ;
      }
   }
} else {
   nkeep = nlist ;
   keeplist = IVinit(nkeep, -1) ;
   IVramp(nkeep, keeplist, 0, 1) ;
}
if ( keepentriesIV != NULL ) {
   IV_sizeAndEntries(keepentriesIV, &nkeepent, &keepent) ;
   maxval = IVL_max(ivl) ;
   if ( maxval >= 0 ) {
      map = IVinit(1 + maxval, -1) ;
      for ( ii = 0 ; ii < nkeepent ; ii++ ) {
         if ( (val = keepent[ii]) >= 0 ) {
            map[val] = ii ;
         }
      }
      maxlistsize = IVL_maxListSize(ivl) ;
      temp = IVinit(maxlistsize, -1) ;
   } else {
      map = NULL ;
   }
} else {
   map = NULL ;
}
/*
   ----------------------------
   initialize the subIVL object
   ----------------------------
*/
IVL_init1(subIVL, IVL_CHUNKED, nkeep) ;
/*
   -----------------------------------
   fill the lists of the subIVL object
   -----------------------------------
*/
for ( ikeep = 0 ; ikeep < nkeep ; ikeep++ ) {
   ilist = keeplist[ikeep] ;
   IVL_listAndSize(ivl, ilist, &size, &list) ;
   if ( map == NULL ) {
      IVL_setList(subIVL, ikeep, size, list) ;
   } else {
      for ( ii = count = 0 ; ii < size ; ii++ ) {
         if ( (val = map[list[ii]]) != -1 ) {
            temp[count++] = val ;
         }
      }
      IVL_setList(subIVL, ikeep, count, temp) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( keeplistIV == NULL ) {
   IVfree(keeplist) ;
}
if ( map != NULL ) {
   IVfree(map) ;
   IVfree(temp) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
