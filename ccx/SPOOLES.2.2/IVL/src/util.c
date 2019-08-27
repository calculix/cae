/*  util.c  */

#include "../IVL.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95sep22, cca
   ----------------------------------------------
*/
int
IVL_sizeOf ( 
   IVL   *ivl 
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_sizeOf(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
nbytes = sizeof(struct _IVL) ;
if ( ivl->nlist > 0 ) {
   nbytes += ivl->nlist * (sizeof(int) + sizeof(int *)) ;
   if ( ivl->type == IVL_SOLO ) {
         nbytes += IVsum(ivl->nlist, ivl->sizes) * sizeof(int) ;
   } else {
      Ichunk   *chunk ;
      for ( chunk = ivl->chunk ; 
            chunk != NULL ; 
            chunk = chunk->next ) {
         nbytes += sizeof(Ichunk) + chunk->size * sizeof(int) ;
      }
   }
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the minimum entry in the lists

   created -- 95sep22, cca
   ---------------------------------------------------
*/
int
IVL_min ( 
   IVL   *ivl 
) {
int   first, i, ilist, minval, nlist, size, val ;
int   *ent ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || (nlist = ivl->nlist) <= 0 ) {
   fprintf(stderr, "\n fatal error in IVL_min(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
first  =  1 ;
minval = -1 ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &ent) ;
   if ( size > 0 ) {
      val = IVmin(size, ent, &i) ;
      if ( first == 1 ) {
         minval = val ;
         first  =  0  ;
      } else if ( minval > val ) {
         minval = val ;
      }
   }
}
return(minval) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the maximum entry in the lists

   created -- 95sep22, cca
   ---------------------------------------------------
*/
int
IVL_max ( 
   IVL   *ivl 
) {
int   first, i, ilist, maxval, nlist, size, val ;
int   *ent ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || (nlist = ivl->nlist) <= 0 ) {
   fprintf(stderr, "\n fatal error in IVL_max(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
first  =  1 ;
maxval = -1 ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &ent) ;
   if ( size > 0 ) {
      val = IVmax(size, ent, &i) ;
      if ( first == 1 ) {
         maxval = val ;
         first  =  0  ;
      } else if ( maxval < val ) {
         maxval = val ;
      }
   }
}
return(maxval) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   return the maximum list size 

   created -- 95sep22, cca
   ----------------------------
*/
int
IVL_maxListSize (
   IVL   *ivl
) {
int   ilist, maxsize, nlist, size ;
int   *ent ;
/*
   -------------------
   check for bad input
   -------------------
*/
if ( ivl == NULL || (nlist = ivl->nlist) <= 0 ) {
   fprintf(stderr, "\n fatal error in IVL_maxListSize(%p)"
           "\n bad input", ivl) ;
   exit(-1) ;
}
for ( ilist = 0, maxsize = 0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &ent) ;
   if ( maxsize < size ) {
      maxsize = size ;
   }
}
return(maxsize) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   return the sum of all the lists

   created -- 95sep29, cca
   -------------------------------
*/
int
IVL_sum (
   IVL   *ivl
) {
int   j, jsize, sum ;
int   *jind ;

if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_sum(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
sum = 0 ;
for ( j = 0 ; j < ivl->nlist ; j++ ) {
   IVL_listAndSize(ivl, j, &jsize, &jind) ;
   if ( jsize > 0 ) {
      sum = sum + IVsum(jsize, jind) ;
   }
}
return(sum) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   sort the adjacency lists in ascending order

   created -- 95sep22, cca
   -------------------------------------------
*/
void
IVL_sortUp ( 
   IVL   *ivl 
) {
int   ilist, nlist, size ;
int   *ent ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || (nlist = ivl->nlist) < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_sortUp(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}

for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &ent) ;
   if ( size > 0 ) {
      IVqsortUp(size, ent) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   create an equivalence map 
   if ( map[j] == map[k] ) then
      the lists for j and k are identical
   endif
   NOTE : each empty list is mapped to a different value

   return value -- pointer to map[] vector

   created -- 95mar15, cca
   -----------------------------------------------------
*/
int *
IVL_equivMap1 (
   IVL   *ivl
) {
int   first, ierr, ii, itest, jtest, nlist, nlist2, ntest, nv2, sum, 
      v, v2, vsize, w, wsize ;
int   *chksum, *issorted, *map, *mark, *vadj, *vids, *wadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || (nlist = ivl->nlist) < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_equivMap(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
if ( nlist == 0 ) {
   return(NULL) ;
}
/*
   --------------
   initialize map
   --------------
*/
map = IVinit(nlist, -1) ;
nlist2 = 0 ;
/*
   ---------------------------------
   sort the lists by their checksums
   ---------------------------------
*/
vids   = IVinit(nlist, -1) ;
chksum = IVinit(nlist,  -1) ;
for ( v = 0, ntest = 0 ; v < nlist ; v++ ) {
   IVL_listAndSize(ivl, v, &vsize, &vadj) ;
   if ( vsize > 0 ) {
/*
      ---------------------------------------------
      list is not empty, store list id and checksum
      ---------------------------------------------
*/
      for ( ii = 0, sum = 0 ; ii < vsize ; ii++ ) {
         sum += vadj[ii] ;
      }
      vids[ntest]   =  v  ;
      chksum[ntest] = sum ;
      ntest++ ;
   } else {
/*
      ------------------------------
      list is empty, map to new list
      ------------------------------
*/
      map[v] = nlist2++ ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n before sort, vids") ;
IVfp80(stdout, ntest, vids, 80, &ierr) ;
fprintf(stdout, "\n before sort, chksum") ;
IVfp80(stdout, ntest, chksum, 80, &ierr) ;
fflush(stdout) ;
#endif
IV2qsortUp(ntest, chksum, vids) ;
#if MYDEBUG > 0
fprintf(stdout, "\n after sort, vids") ;
IVfp80(stdout, ntest, vids, 80, &ierr) ;
fprintf(stdout, "\n after sort, chksum") ;
IVfp80(stdout, ntest, chksum, 80, &ierr) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------------------------------------
   loop over the nonempty lists in the order of increasing checksums
   -----------------------------------------------------------------
*/
issorted = IVinit(nlist, -1) ;
for ( itest = 0 ; itest < ntest ; itest++ ) {
   v = vids[itest] ;
   if ( map[v] == -1 ) {
/*
      -------------------------------------------------
      list v has not been found to be indistinguishable
      to any other list, map it to a new list
      -------------------------------------------------
*/
      map[v] = nlist2++ ;
#if MYDEBUG > 0
      fprintf(stdout, "\n setting map[%d] = %d, chksum[%d] = %d", 
              v, map[v], itest, chksum[itest]) ;
      fflush(stdout) ;
#endif
/*
      --------------------------------------
      loop over lists with the same checksum
      --------------------------------------
*/
      IVL_listAndSize(ivl, v, &vsize, &vadj) ;
      first = 1 ;
      for ( jtest = itest + 1 ; jtest < ntest ; jtest++ ) {
         w = vids[jtest] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n    comparing with %d, chksum[%d] = %d",
              w, jtest, chksum[jtest]) ;
         fflush(stdout) ;
#endif
         if ( chksum[itest] != chksum[jtest] ) {
/*
            --------------------------------------------------
            checksums are not equal, list v cannot be the same
            as any following list, break out of test loop
            --------------------------------------------------
*/
            break ;
         } else {
/*
            -----------------------------------------------------
            lists v and w have the same checksum
            if the list sizes are the same then compare the lists
            -----------------------------------------------------
*/
            IVL_listAndSize(ivl, w, &wsize, &wadj) ;
#if MYDEBUG > 0
         fprintf(stdout, "\n    vsize = %d, wsize = %d", vsize, wsize) ;
         fflush(stdout) ;
#endif
            if ( vsize == wsize ) {
               if ( issorted[v] != 1 ) {
#if MYDEBUG > 0
                  fprintf(stdout, "\n    sorting list for %d", v) ;
                  fflush(stdout) ;
#endif
                  issorted[v] = 1 ;
                  IVqsortUp(vsize, vadj) ;
               }
               if ( issorted[w] != 1 ) {
#if MYDEBUG > 0
                  fprintf(stdout, "\n    sorting list for %d", w) ;
                  fflush(stdout) ;
#endif
                  issorted[w] = 1 ;
                  IVqsortUp(wsize, wadj) ;
               }
               for ( ii = 0 ; ii < vsize ; ii++ ) {
                  if ( vadj[ii] != wadj[ii] ) {
                     break ;
                  }
               }
               if ( ii == vsize ) {
/*
                  ----------------------------------
                  lists are identical, set map for w
                  ----------------------------------
*/
#if MYDEBUG > 0
                  fprintf(stdout, "\n    lists are identical") ;
                  fflush(stdout) ;
#endif
                  map[w] = map[v] ;
               }
            }
         }
      }
   }
}
IVfree(issorted) ;
IVfree(chksum)   ;
IVfree(vids)     ;
#if MYDEBUG > 0
fprintf(stdout, "\n initial map") ;
IVfp80(stdout, nlist, map, 80, &ierr) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------------
   now modify the map to ensure 
   if v2 < w2 then
      min { v | map[v] = v2 } < min { w | map[w] = w2 }
   endif
   ----------------------------------------------------
*/
mark = IVinit(nlist2, -1) ;
for ( v = 0, nv2 = 0 ; v < nlist ; v++ ) {
   v2 = map[v] ;
   if ( mark[v2] == -1 ) {
      mark[v2] = nv2++ ;
   }
   map[v] = mark[v2] ;
}
IVfree(mark) ;
#if MYDEBUG > 0
fprintf(stdout, "\n final map") ;
IVfp80(stdout, nlist, map, 80, &ierr) ;
fflush(stdout) ;
#endif

return(map) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   create an equivalence map 
   if ( map[j] == map[k] ) then
      the lists for j and k are identical
   endif
   NOTE : each empty list is mapped to a different value

   return value -- pointer to map IV object

   created -- 96mar15, cca
   -----------------------------------------------------
*/
IV *
IVL_equivMap2 (
   IVL   *ivl
) {
int   *map ;
IV    *mapIV ;

if ( (map = IVL_equivMap1(ivl)) == NULL ) {
   mapIV = NULL ;
} else { 
   mapIV = IV_new() ;
   IV_init2(mapIV, ivl->nlist, ivl->nlist, 1, map) ;
}
return(mapIV) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   overwrite each list with new entries

   created -- 96oct03, cca
   ------------------------------------
*/
void
IVL_overwrite (
   IVL   *ivl,
   IV    *oldToNewIV
) {
int   ii, ilist, nlist, range, size ;
int   *list, *oldToNew ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || oldToNewIV == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_overwrite(%p,%p)"
           "\n bad input\n", ivl, oldToNewIV) ;
   exit(-1) ;
}
oldToNew = IV_entries(oldToNewIV) ;
range    = IV_size(oldToNewIV) ;
nlist    = ivl->nlist ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &list) ;
   for ( ii = 0 ; ii < size ; ii++ ) {
      if ( 0 <= list[ii] && list[ii] < range ) {
         list[ii] = oldToNew[list[ii]] ;
      }
   } 
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   given an IVL object and a map from old list entries 
   to new list entries, create a new IVL object whose
   new list entries contain no duplicates.

   return value -- pointer to new IVL object

   created -- 96nov07, cca
   ---------------------------------------------------
*/
IVL *
IVL_mapEntries (
   IVL   *ivl,
   IV    *mapIV
) {
int   count, ierr, ii, ilist, maxsize, nlist, range, size, value ;
int   *list, *map, *newlist ;
IVL   *newIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || mapIV == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_mapEntries(%p,%p)"
           "\n bad input\n", ivl, mapIV) ;
   exit(-1) ;
}
nlist = ivl->nlist ;
range = IV_size(mapIV) ;
map   = IV_entries(mapIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n nlist = %d, range = %d, map = %p",
        nlist, range, map) ;
#endif
if (  nlist <= 0 || range < 0 || map == NULL ) {
   return(NULL) ;
}
/*
   -------------------------
   create the new IVL object
   -------------------------
*/
newIVL = IVL_new();
IVL_init1(newIVL, IVL_CHUNKED, nlist) ;
/*
   -------------------
   loop over the lists
   -------------------
*/
maxsize = IVL_maxListSize(ivl) ;
newlist = IVinit(maxsize, -1) ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &list) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n list %d :", ilist) ;
   IVfp80(stdout, size, list, 10, &ierr) ;
#endif
   for ( ii = 0, count = 0 ; ii < size ; ii++ ) {
      if ( 0 <= list[ii] && list[ii] < range ) {
/*
         -----------------------------------------
         old entry is in range, store mapped value
         -----------------------------------------
*/
#if MYDEBUG > 0
         fprintf(stdout, "\n    newlist[%d] = map[%d] = %d",
                 count, list[ii], map[list[ii]]) ;
#endif
         newlist[count++] = map[list[ii]] ;
      }
   }
   if ( count > 0 ) {
/*
      ------------------------------------
      sort the new list in ascending order
      ------------------------------------
*/
#if MYDEBUG > 0
      fprintf(stdout, "\n    unsorted list %d :", ilist) ;
      IVfp80(stdout, count, newlist, 10, &ierr) ;
#endif
      IVqsortUp(count, newlist) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    sorted list %d :", ilist) ;
      IVfp80(stdout, count, newlist, 10, &ierr) ;
#endif
/*
      -----------------------
      purge duplicate entries
      -----------------------
*/
      size  = count ;
      value = -1 ;
      for ( ii = count = 0 ; ii < size ; ii++ ) {
         if ( newlist[ii] != value ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n    keeping entry %d", newlist[ii]) ;
#endif
            newlist[count++] = newlist[ii] ;
            value = newlist[ii] ;
         }
      }
   }
/*
   ----------------------------------
   set the list in the new IVL object
   ----------------------------------
*/
   IVL_setList(newIVL, ilist, count, newlist) ;
}
IVfree(newlist) ;

return(newIVL) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   IVL object ivl1 absorbs the lists and data from IVL object ivl2.
   the lists in ivl2 are mapped into lists in ivl1 using the mapIV
   IV object.

   created -- 96dec06, cca
   ----------------------------------------------------------------
*/
void
IVL_absorbIVL (
   IVL   *ivl1,
   IVL   *ivl2,
   IV    *mapIV
) {
Ichunk   *chunk ;
int      ilist, jlist, nlist2, size ;
int      *ivec, *map ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl1 == NULL || ivl2 == NULL || mapIV == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_absorbIVL(%p,%p,%p)"
           "\n bad input\n", ivl1, ivl2, mapIV) ;
   exit(-1) ;
}
if ( (map = IV_entries(mapIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_absorbIVL(%p,%p,%p)"
           "\n IV_entries(mapIV) is NULL\n", ivl1, ivl2, mapIV) ;
   exit(-1) ;
}
/*
   --------------------------------------------
   check that the sizes of ivl2 and mapIV agree
   --------------------------------------------
*/
if ( IV_size(mapIV) != (nlist2 = ivl2->nlist) ) {
   fprintf(stderr, "\n fatal error in IVL_absorbIVL(%p,%p,%p)"
           "\n ivl2->nlist = %d, IV_size(mapIV) = %d\n", 
           ivl1, ivl2, mapIV, nlist2, IV_size(mapIV)) ;
   exit(-1) ;
}
/*
   -------------------------------
   for each list in ivl2
      get size and pointer
      get mapped list in ivl1
      set size and pointer in ivl1
   -------------------------------
*/
for ( ilist = 0 ; ilist < nlist2 ; ilist++ ) {
   IVL_listAndSize(ivl2, ilist, &size, &ivec) ;
   if ( (jlist = map[ilist]) >= 0 ) {
      IVL_setPointerToList(ivl1, jlist, size, ivec) ;
   }
}
if ( (chunk = ivl2->chunk) != NULL ) {
/*
   -------------------------------------------
   move the chunks of memory from ivl2 to ivl1
   -------------------------------------------
*/
   while ( chunk->next != NULL ) {
      chunk = chunk->next ;
   }
   chunk->next = ivl1->chunk ;
   ivl1->chunk = ivl2->chunk ;
   ivl2->chunk = NULL ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- expand the lists in an IVL object.

   this method was created in support of a symbolic factorization.
   an IVL object is constructed using a compressed graph.
   it must be expanded to reflect the compressed graph.
   the number of lists does not change (there is one list per front)
   but the size of each list may change. so we create a new IVL
   object that contains entries for the uncompressed graph.

   created -- 97feb13, cca
   -----------------------------------------------------------------
*/
IVL *
IVL_expand (
   IVL   *ivl,
   IV    *eqmapIV
) {
int   count, ii, ilist, jj, kk, nlist1, nlist2, nvtx, size ;
int   *head, *link, *list, *map, *temp ;
IVL   *ivl2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || eqmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_expand(%p,%p)"
           "\n bad input\n", ivl, eqmapIV) ;
   exit(-1) ;
}
nlist1 = ivl->nlist ;
/*
   -------------------------------
   get the head[]/link[] structure
   -------------------------------
*/
IV_sizeAndEntries(eqmapIV, &nlist2, &map) ;
nvtx = 1 + IV_max(eqmapIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n nlist1 = %d, nlist2 = %d", nlist1, nlist2) ;
fflush(stdout) ;
#endif
head = IVinit(nvtx,   -1) ;
link = IVinit(nlist2, -1) ;
for ( ii = nlist2 - 1 ; ii >= 0 ; ii-- ) {
   if ( (jj = map[ii]) < 0 || jj >= nvtx ) {
      fprintf(stderr, "\n fatal error in IVL_expand(%p,%p)"
              "\n nlist1 = %d, nvtx = %d, map[%d] = %d\n",
              ivl, eqmapIV, nlist1, nvtx, ii, jj) ;
      exit(-1) ;
   }
   link[ii] = head[jj] ;
   head[jj] = ii ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n head/link structure created") ;
fflush(stdout) ;
#endif
/*
   ---------------------------
   allocate the new IVL object
   ---------------------------
*/
ivl2 = IVL_new() ;
IVL_init1(ivl2, IVL_CHUNKED, nlist1) ;
temp = IVinit(nlist2, -1) ;
for ( ilist = 0 ; ilist < nlist1 ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &list) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n working on list %d", ilist) ;
   IVfprintf(stdout, size, list) ;
   fflush(stdout) ;
#endif
   for ( ii = 0, count = 0 ; ii < size ; ii++ ) {
      jj = list[ii] ;
      for ( kk = head[jj] ; kk != -1 ; kk = link[kk] ) {
         temp[count++] = kk ;
      }
   }
   IVL_setList(ivl2, ilist, count, temp) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head) ;
IVfree(link) ;
IVfree(temp) ;

return(ivl2) ; }

/*--------------------------------------------------------------------*/
