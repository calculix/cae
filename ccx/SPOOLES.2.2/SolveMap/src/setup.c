/*  setup.c  */

#include "../SolveMap.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to set up the linked lists for the forward solve
              local to process myid

   created -- 98mar19, cca
   -----------------------------------------------------------
*/
IP **
SolveMap_forwardSetup (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int   count, J, K, loc, nblock, nfront ;
int   *colids, *map, *rowids ;
IP    *ip, *nextip ;
IP    **heads ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_forwardSetup(%p,%d)"
           "\n solvemap is NULL\n") ;
   exit(-1) ;
}
if ( myid < 0 || myid >= solvemap->nproc ) {
   fprintf(stderr, "\n fatal error in SolveMap_forwardSetup(%p,%d)"
           "\n myid %d, solvemap->nproc %d\n", myid, solvemap->nproc) ;
   exit(-1) ;
}
if ( solvemap == NULL || myid < 0 || myid >= solvemap->nproc ) {
   fprintf(stderr, "\n fatal error in SolveMap_forwardSetup(%p,%d)"
           "\n bad input\n", solvemap, myid) ;
   exit(-1) ;
}
nfront = solvemap->nfront ;
if ( solvemap->symmetryflag == 2 ) {
   nblock = solvemap->nblockLower ;
   map    = solvemap->mapLower    ;
   rowids = solvemap->rowidsLower ;
   colids = solvemap->colidsLower ;
} else {
   nblock = solvemap->nblockUpper ;
   map    = solvemap->mapUpper    ;
   rowids = solvemap->colidsUpper ;
   colids = solvemap->rowidsUpper ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside SolveMap_forwardSetup()") ;
   fprintf(msgFile, ", %d blocks", nblock) ;
   fprintf(msgFile, "\n map") ;
   IVfprintf(msgFile, nblock, map) ;
   fprintf(msgFile, "\n rowids") ;
   IVfprintf(msgFile, nblock, rowids) ;
   fprintf(msgFile, "\n colids") ;
   IVfprintf(msgFile, nblock, colids) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   count the number of necessary IP objects
   ----------------------------------------
*/
for ( loc = count = 0 ; loc < nblock ; loc++ ) {
   if ( map[loc] == myid ) {
      count++ ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n count = %d", count) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   allocate the IP objects and vector of pointers
   ----------------------------------------------
*/
ALLOCATE(heads, struct _IP *, nfront + 2) ;
for ( J = 0 ; J < nfront ; J++ ) {
   heads[J] = NULL ;
}
heads[nfront] = NULL ;
if ( count > 0 ) {
   heads[nfront+1] = IP_init(count, IP_FORWARD) ;
} else {
   heads[nfront+1] = NULL ;
}
if ( count > 0 ) {
/*
   ----------------
   set up the lists
   ----------------
*/
   for ( loc = 0, nextip = heads[nfront + 1] ; loc < nblock ; loc++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n map[%d] = %d", loc, map[loc]) ;
         fflush(msgFile) ;
      }
      if ( map[loc] == myid ) {
         ip = nextip ;
         nextip = ip->next ;
         J = colids[loc] ;
         K = rowids[loc] ;
         ip->val = J ;
         ip->next = heads[K] ;
         heads[K] = ip ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", linking (K,J) = (%d,%d)" 
                    ", heads[%d] = %p", K, J, K, heads[K]) ;
            fflush(msgFile) ;
         }
      }
   }
}
return(heads) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to set up the linked lists for the backward solve
              local to process myid

   created -- 98mar19, cca
   ------------------------------------------------------------
*/
IP **
SolveMap_backwardSetup (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int   count, J, K, loc, nblock, nfront ;
int   *colids, *map, *rowids ;
IP    *ip, *nextip ;
IP    **heads ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || myid < 0 || myid >= solvemap->nproc ) {
   fprintf(stderr, "\n fatal error in SolveMap_backwardSetup(%p,%d)"
           "\n bad input\n", solvemap, myid) ;
   exit(-1) ;
}
nfront = solvemap->nfront ;
nblock = solvemap->nblockUpper ;
map    = solvemap->mapUpper    ;
rowids = solvemap->rowidsUpper ;
colids = solvemap->colidsUpper ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n nfront %d, nblock %d, map %p, rowids %p, colids %p",
           nfront, nblock, map, rowids, colids) ;
   fprintf(msgFile, "\n\n inside SolveMap_backwardSetup()") ;
   fprintf(msgFile, ", %d blocks", nblock) ;
   fflush(msgFile) ;
   fprintf(msgFile, "\n map = %p", map) ;
   fflush(msgFile) ;
   IVfprintf(msgFile, nblock, map) ;
   fflush(msgFile) ;
   fprintf(msgFile, "\n rowids = %p", rowids) ;
   fflush(msgFile) ;
   IVfprintf(msgFile, nblock, rowids) ;
   fflush(msgFile) ;
   fprintf(msgFile, "\n colids = %p", colids) ;
   fflush(msgFile) ;
   IVfprintf(msgFile, nblock, colids) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   count the number of necessary IP objects
   ----------------------------------------
*/
for ( loc = count = 0 ; loc < nblock ; loc++ ) {
   if ( map[loc] == myid ) {
      count++ ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n count = %d", count) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   allocate the IP objects and vector of pointers
   ----------------------------------------------
*/
ALLOCATE(heads, struct _IP *, nfront + 2) ;
for ( J = 0 ; J < nfront ; J++ ) {
   heads[J] = NULL ;
}
heads[nfront] = NULL ;
if ( count > 0 ) {
   heads[nfront+1] = IP_init(count, IP_FORWARD) ;
} else {
   heads[nfront+1] = NULL ;
}
if ( count > 0 ) {
/*
   ----------------
   set up the lists
   ----------------
*/
   for ( loc = 0, nextip = heads[nfront + 1] ; loc < nblock ; loc++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n map[%d] = %d", loc, map[loc]) ;
         fflush(msgFile) ;
      }
      if ( map[loc] == myid ) {
         ip = nextip ;
         nextip = ip->next ;
         J = rowids[loc] ;
         K = colids[loc] ;
         ip->val = K ;
         ip->next = heads[J] ;
         heads[J] = ip ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", linking (J,K) = (%d,%d)" 
                    ", heads[%d] = %p", J, K, K, heads[J]) ;
            fflush(msgFile) ;
         }
      }
   }
}
return(heads) ; }

/*--------------------------------------------------------------------*/
