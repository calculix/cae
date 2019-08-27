/*  transform.c  */

#include "../ETree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   transform an ETree object by 
   (1) merging small fronts into larger fronts
       using the ETree_mergeFrontsOne() method
   (2) merging small fronts into larger fronts
       using the ETree_mergeFrontsAll() method
   (3) merging small fronts into larger fronts
       using the ETree_mergeFrontsAny() method
   (4) split a large front into a chain of smaller fronts
       using the ETree_splitFronts() method

   created  -- 96jun27, cca
   ------------------------------------------------------
*/
ETree *
ETree_transform (
   ETree   *etree,
   int     vwghts[],
   int     maxzeros,
   int     maxfrontsize,
   int     seed
) {
ETree   *etree2 ;
int     nfront, nvtx ;
IV      *nzerosIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || maxfrontsize <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_transform(%p,%p,%d,%d,%d)"
           "\n bad input\n", etree, vwghts, maxzeros, maxfrontsize, 
           seed) ;
   exit(-1) ;
}
nzerosIV = IV_new();
IV_init(nzerosIV, nfront, NULL) ;
IV_fill(nzerosIV, 0) ;
/*
   --------------------------
   first, merge only children
   --------------------------
*/
etree2 = ETree_mergeFrontsOne(etree, maxzeros, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   --------------------------
   second, merge all children
   --------------------------
*/
etree2 = ETree_mergeFrontsAll(etree, maxzeros, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   -------------------------
   third, merge any children
   -------------------------
*/
etree2 = ETree_mergeFrontsAny(etree, maxzeros, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   -----------------------------------
   fourth, split large interior fronts
   -----------------------------------
*/
etree2 = ETree_splitFronts(etree, vwghts, maxfrontsize, seed) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(nzerosIV) ;

return(etree) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   transform an ETree object by 
   (1) merging small fronts into larger fronts
       using the ETree_mergeFrontsOne() method
   (2) merging small fronts into larger fronts
       using the ETree_mergeFrontsAll() method
   (3) split a large front into a chain of smaller fronts
       using the ETree_splitFronts() method

   created  -- 96jun27, cca
   ------------------------------------------------------
*/
ETree *
ETree_transform2 (
   ETree   *etree,
   int     vwghts[],
   int     maxzeros,
   int     maxfrontsize,
   int     seed
) {
ETree   *etree2 ;
int     nfront, nvtx ;
IV      *nzerosIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || maxfrontsize <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_transform2(%p,%p,%d,%d,%d)"
           "\n bad input\n", etree, vwghts, maxzeros, maxfrontsize, 
           seed) ;
   exit(-1) ;
}
nzerosIV = IV_new();
IV_init(nzerosIV, nfront, NULL) ;
IV_fill(nzerosIV, 0) ;
/*
   --------------------------
   first, merge only children
   --------------------------
*/
etree2 = ETree_mergeFrontsOne(etree, maxzeros, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   --------------------------
   second, merge all children
   --------------------------
*/
etree2 = ETree_mergeFrontsAll(etree, maxzeros, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   -----------------------------------
   fourth, split large interior fronts
   -----------------------------------
*/
etree2 = ETree_splitFronts(etree, vwghts, maxfrontsize, seed) ;
ETree_free(etree) ;
etree = etree2 ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(nzerosIV) ;

return(etree) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- merge the front tree allowing only chains of nodes to
      merge that create at most maxzeros zero entries inside a front

   return -- 
      IV object that has the old front to new front map

   created -- 98jan29, cca
   --------------------------------------------------------------------
*/
ETree *
ETree_mergeFrontsOne (
   ETree   *etree,
   int     maxzeros,
   IV      *nzerosIV
) {
ETree   *etree2 ;
int     costJ, J, K, nfront, nvtx, nnew ;
int     *bndwghts, *fch, *map, *nodwghts, *nzeros, *rep, *sib, *temp ;
IV      *mapIV ;
Tree    *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL || nzerosIV == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_mergeFrontsOne(%p,%d,%p)"
           "\n bad input\n", etree, maxzeros, nzerosIV) ;
   exit(-1) ;
}
if ( IV_size(nzerosIV) != nfront ) {
   fprintf(stderr, "\n fatal error in ETree_mergeFrontsOne(%p,%d,%p)"
           "\n size(nzerosIV) = %d, nfront = %d\n", 
           etree, maxzeros, nzerosIV, IV_size(nzerosIV), nfront) ;
   exit(-1) ;
}
nzeros = IV_entries(nzerosIV) ;
tree     = etree->tree ;
fch      = ETree_fch(etree) ;
sib      = ETree_sib(etree) ;
/*
   ----------------------
   set up working storage
   ----------------------
*/
nodwghts = IVinit(nfront, 0) ;
IVcopy(nfront, nodwghts, ETree_nodwghts(etree)) ;
bndwghts = ETree_bndwghts(etree) ;
rep = IVinit(nfront, -1) ;
IVramp(nfront, rep, 0, 1) ;
/*
   ------------------------------------------
   perform a post-order traversal of the tree
   ------------------------------------------
*/
for ( K = Tree_postOTfirst(tree) ;
      K != -1 ;
      K = Tree_postOTnext(tree, K) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n ##### visiting front %d", K) ;
   fflush(stdout) ;
#endif
   if ( (J = fch[K]) != -1 && sib[J] == -1 ) {
      costJ = nodwghts[J]*(nodwghts[K] + bndwghts[K] - bndwghts[J]) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n nzeros[%d] = %d, costJ = %d",
              J, nzeros[J], costJ) ;
      fflush(stdout) ;
#endif
      if ( nzeros[J] + costJ <= maxzeros ) {
         rep[J] = K ;
         nodwghts[K] += nodwghts[J] ;
         nzeros[K] = nzeros[J] + costJ ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n merging %d into %d, |%d| = %d, nzeros = %d",
                 J, K, K, nodwghts[K], nzeros[K]) ;
         fflush(stdout) ;
#endif
      }
   }
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n whoa, finished") ;
   fflush(stdout) ;
#endif
/*
   -------------------------------------------------
   take the map from fronts to representative fronts
   and make the map from old fronts to new fronts
   -------------------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, nfront, NULL) ;
map   = IV_entries(mapIV) ;
for ( J = 0, nnew = 0 ; J < nfront ; J++ ) {
   if ( rep[J] == J ) {
      map[J] = nnew++ ;
   } else {
      K = J ;
      while ( rep[K] != K ) {
         K = rep[K] ;
      }
      rep[J] = K ;
   }
}
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (K = rep[J]) != J ) {
      map[J] = map[K] ;
   }
}
/*
   -------------------------------
   get the compressed ETree object
   -------------------------------
*/
etree2 = ETree_compress(etree, mapIV) ;
/*
   -------------------------
   remap the nzeros[] vector
   -------------------------
*/
temp = IVinit(nfront, NULL) ;
IVcopy(nfront, temp, nzeros) ;
IV_setSize(nzerosIV, nnew) ;
nzeros = IV_entries(nzerosIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( rep[J] == J ) {
      nzeros[map[J]] = temp[J] ;
   }
}
IVfree(temp) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(nodwghts) ;
IVfree(rep)      ;
IV_free(mapIV)   ;

return(etree2) ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- merge the front tree allowing a parent 
              to absorb all children when that creates 
              at most maxzeros zero entries inside a front

   return -- 
      IV object that has the old front to new front map

   created -- 98jan29, cca
   -------------------------------------------------------
*/
ETree *
ETree_mergeFrontsAll (
   ETree   *etree,
   int     maxzeros,
   IV      *nzerosIV
) {
ETree   *etree2 ;
int     cost, J, Jall, K, KandBnd, nfront, nvtx, nnew ;
int     *bndwghts, *fch, *map, *nodwghts, *nzeros, *rep, *sib, *temp ;
IV      *mapIV ;
Tree    *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL || nzerosIV == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_mergeFrontsAll(%p,%d,%p)"
           "\n bad input\n", etree, maxzeros, nzerosIV) ;
   if ( etree != NULL ) {
      fprintf(stderr, "\n nfront = %d, nvtx = %d",
              etree->nfront, etree->nvtx) ;
   }
   exit(-1) ;
}
if ( IV_size(nzerosIV) != nfront ) {
   fprintf(stderr, "\n fatal error in ETree_mergeFrontsAll(%p,%d,%p)"
           "\n size(nzerosIV) = %d, nfront = %d\n", 
           etree, maxzeros, nzerosIV, IV_size(nzerosIV), nfront) ;
   exit(-1) ;
}
nzeros = IV_entries(nzerosIV) ;
/*
   ----------------------
   set up working storage
   ----------------------
*/
tree     = etree->tree ;
fch      = ETree_fch(etree) ;
sib      = ETree_sib(etree) ;
nodwghts = IVinit(nfront, 0) ;
IVcopy(nfront, nodwghts, ETree_nodwghts(etree)) ;
bndwghts = ETree_bndwghts(etree) ;
rep = IVinit(nfront, -1) ;
IVramp(nfront, rep, 0, 1) ;
/*
   ------------------------------------------
   perform a post-order traversal of the tree
   ------------------------------------------
*/
for ( K = Tree_postOTfirst(tree) ;
      K != -1 ;
      K = Tree_postOTnext(tree, K) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n ##### visiting front %d", K) ;
   fflush(stdout) ;
#endif
   if ( (J = fch[K]) != -1 ) {
      KandBnd = nodwghts[K] + bndwghts[K] ;
      Jall = 0 ;
      cost = 2*nzeros[K] ;
      for ( J = fch[K] ; J != -1 ; J = sib[J] ) {
         Jall += nodwghts[J] ;
         cost -= nodwghts[J]*nodwghts[J] ;
         cost += 2*nodwghts[J]*(KandBnd - bndwghts[J]) ;
         cost += 2*nzeros[J] ;
      }
      cost += Jall*Jall ;
      cost = cost/2 ;
#if MYDEBUG > 0
      fprintf(stdout, "\n cost = %d", cost) ;
      fflush(stdout) ;
#endif
      if ( cost <= maxzeros ) {
         for ( J = fch[K] ; J != -1 ; J = sib[J] ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n merging %d into %d", J, K) ;
            fflush(stdout) ;
#endif
            rep[J] = K ;
            nodwghts[K] += nodwghts[J] ;
         }
         nzeros[K] = cost ;
      }
   }
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n whoa, finished") ;
   fflush(stdout) ;
#endif
/*
   -------------------------------------------------
   take the map from fronts to representative fronts
   and make the map from old fronts to new fronts
   -------------------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, nfront, NULL) ;
map   = IV_entries(mapIV) ;
for ( J = 0, nnew = 0 ; J < nfront ; J++ ) {
   if ( rep[J] == J ) {
      map[J] = nnew++ ;
   } else {
      K = J ;
      while ( rep[K] != K ) {
         K = rep[K] ;
      }
      rep[J] = K ;
   }
}
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (K = rep[J]) != J ) {
      map[J] = map[K] ;
   }
}
/*
   -------------------------------
   get the compressed ETree object
   -------------------------------
*/
etree2 = ETree_compress(etree, mapIV) ;
/*
   -------------------------
   remap the nzeros[] vector
   -------------------------
*/
temp = IVinit(nfront, NULL) ;
IVcopy(nfront, temp, nzeros) ;
IV_setSize(nzerosIV, nnew) ;
nzeros = IV_entries(nzerosIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( rep[J] == J ) {
      nzeros[map[J]] = temp[J] ;
   }
}
IVfree(temp) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(nodwghts) ;
IVfree(rep)      ;
IV_free(mapIV)   ;

return(etree2) ; }
   
/*--------------------------------------------------------------------*/
/*
   ---------------------------
   static prototype definition
   ---------------------------
*/
static void visitAny ( int K, int par[], int fch[], int sib[], 
                       int nodwghts[], int bndwghts[], int map[], 
                       int cost[], int nzeros[], int maxzeros) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- merge the front tree allowing at most
              maxzeros zero entries inside a front

   return -- 
      IV object that has the old front to new front map

   created -- 96jun23, cca
   modified -- 97dec18, cca
      bug fixed that incorrectly counted the number of zeros in a front
   --------------------------------------------------------------------
*/
ETree *
ETree_mergeFrontsAny (
   ETree   *etree,
   int     maxzeros,
   IV      *nzerosIV
) {
ETree   *etree2 ;
int     J, K, nfront, nvtx, nnew ;
int     *bndwghts, *cost, *fch, *map, *nodwghts, 
        *nzeros, *par, *place, *rep, *sib, *temp ;
IV      *mapIV ;
Tree    *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_mergeFrontsAny(%p,%d)"
           "\n bad input\n", etree, maxzeros) ;
   exit(-1) ;
}
if ( IV_size(nzerosIV) != nfront ) {
   fprintf(stderr, "\n fatal error in ETree_mergeFrontsAny(%p,%d,%p)"
           "\n size(nzerosIV) = %d, nfront = %d\n", 
           etree, maxzeros, nzerosIV, IV_size(nzerosIV), nfront) ;
   exit(-1) ;
}
nzeros = IV_entries(nzerosIV) ;
tree     = etree->tree ;
nodwghts = IVinit(nfront, 0) ;
bndwghts = IVinit(nfront, 0) ;
par = IVinit(nfront, -1) ;
fch = IVinit(nfront, -1) ;
sib = IVinit(nfront, -1) ;
IVcopy(nfront, par, tree->par) ;
IVcopy(nfront, fch, tree->fch) ;
IVcopy(nfront, sib, tree->sib) ;
IVcopy(nfront, nodwghts, IV_entries(etree->nodwghtsIV)) ;
IVcopy(nfront, bndwghts, IV_entries(etree->bndwghtsIV)) ;
/*
   ----------------------
   set up working storage
   ----------------------
*/
rep = IVinit(nfront, -1) ;
IVramp(nfront, rep, 0, 1) ;
cost   = IVinit(nfront, 0) ;
/*
   ------------------------------------------
   perform a post-order traversal of the tree
   ------------------------------------------
*/
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n ##### visiting front %d", J) ;
   fflush(stdout) ;
#endif
   visitAny(J, par, fch, sib, nodwghts, bndwghts, 
            rep, cost, nzeros, maxzeros) ;
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n whoa, finished") ;
   fflush(stdout) ;
#endif
/*
   -------------------------------------------------
   take the map from fronts to representative fronts
   and make the map from old fronts to new fronts
   -------------------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, nfront, NULL) ;
map   = IV_entries(mapIV) ;
place = IVinit(nfront, -1) ;
for ( J = 0, nnew = 0 ; J < nfront ; J++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n rep[%d] = %d", J, rep[J]) ;
   fflush(stdout) ;
#endif
   if ( rep[J] != J ) {
      K = J ;
      while ( rep[K] != K ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n    rep[%d] = %d", K, rep[K]) ;
      fflush(stdout) ;
#endif
         K = rep[K] ;
      }
      rep[J] = K ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting rep[%d] = %d", J, rep[J]) ;
      fflush(stdout) ;
#endif
   } else {
      place[J] = nnew++ ;
   }
}
for ( J = 0 ; J < nfront ; J++ ) {
   K = rep[J] ;
   map[J] = place[K] ;
}
/*
   -------------------------------
   get the compressed ETree object
   -------------------------------
*/
etree2 = ETree_compress(etree, mapIV) ;
/*
   -------------------------
   remap the nzeros[] vector
   -------------------------
*/
temp = IVinit(nfront, NULL) ;
IVcopy(nfront, temp, nzeros) ;
IV_setSize(nzerosIV, nnew) ;
nzeros = IV_entries(nzerosIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( rep[J] == J ) {
      nzeros[map[J]] = temp[J] ;
   }
}
IVfree(temp) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(par)      ;
IVfree(fch)      ;
IVfree(sib)      ;
IVfree(nodwghts) ;
IVfree(bndwghts) ;
IVfree(rep)      ;
IVfree(cost)     ;
IVfree(place)    ;
IV_free(mapIV)   ;

return(etree2) ; }
   
/*--------------------------------------------------------------------*/
static void
visitAny ( 
   int    K,
   int    par[],
   int    fch[],
   int    sib[],
   int    nodwghts[],
   int    bndwghts[],
   int    rep[],
   int    cost[],
   int    nzeros[],
   int    maxzeros
) {
int   bestJ, firstI, J, lastI, nextJ, prevJ ;

if ( fch[K] == -1 ) {
   return ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n inside visitAny(%d), nzeros(%d) = %d", 
        K, K, nzeros[K]) ;
#endif
/*
   ----------------------------------
   find the child with the least cost
   ----------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n    nodwght %d, bndwght %d",
           nodwghts[K], bndwghts[K]) ;
#endif
bestJ = -1 ;
for ( J = fch[K] ; J != -1 ; J = sib[J] ) {
   cost[J] = nzeros[J] 
           + nodwghts[J] * (nodwghts[K] + bndwghts[K] - bndwghts[J]) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n    child %d, nodwght %d, bndwght %d, cost %d",
           J, nodwghts[J], bndwghts[J], cost[J]) ;
#endif
   if (  bestJ == -1 
      || cost[J] < cost[bestJ]
      || (cost[J] == cost[bestJ] && nodwghts[J] < nodwghts[bestJ]) ) {
      bestJ = J ;
   }
}
if (  (cost[bestJ] + nzeros[K] > maxzeros)
   || (sib[fch[K]] != -1 && par[K] == -1)   ) {
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n    no merge: cost[%d] + nzeros[%d] = %d + %d > %d",
           bestJ, K, cost[bestJ], nzeros[K], maxzeros) ;
#endif
/*
   --------------------------------
   no child can be absorbed, return
   --------------------------------
*/
   return ;
} 
#if MYDEBUG > 0
fprintf(stdout, "\n    merging child %d into %d", bestJ, K) ;
#endif
/*
   -------------------------
   absorb child bestJ into K
   -------------------------
*/
for ( J = fch[K], prevJ = -1 ; J != bestJ ; J = sib[J] ) {
   prevJ = J ;
}
nextJ = sib[bestJ] ;
#if MYDEBUG > 0
fprintf(stdout, "\n    previous sibling = %d, next sibling = %d",
        prevJ, nextJ) ;
#endif
if ( (firstI = fch[bestJ]) == -1 ) {
   if ( prevJ == -1 ) {
      fch[K] = nextJ ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting fch[%d] = %d", K, fch[K]) ;
#endif
   } else {
      sib[prevJ] = nextJ ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting sib[%d] = %d", prevJ, sib[prevJ]) ;
#endif
   }
} else {
   firstI = fch[bestJ] ;
   par[firstI] = K ;
#if MYDEBUG > 0
   fprintf(stdout, "\n    setting par[%d] = %d", firstI, par[firstI]) ;
#endif
   if ( (lastI = sib[firstI]) != -1 ) {
      while ( sib[lastI] != -1 ) {
         par[lastI] = K ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n    setting par[%d] = %d", lastI, par[lastI]) ;
#endif
         lastI = sib[lastI] ;
      }
      par[lastI] = K ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting par[%d] = %d", lastI, par[lastI]) ;
#endif
   }
   if ( prevJ == -1 ) {
      fch[K] = firstI ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting fch[%d] = %d", K, fch[K]) ;
#endif
   } else {
      sib[prevJ] = firstI ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting sib[%d] = %d", prevJ, sib[prevJ]) ;
#endif
   }
   if ( lastI != -1 ) {
      sib[lastI] = nextJ ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    setting sib[%d] = %d", lastI, sib[lastI]) ;
#endif
   }
}
rep[bestJ]  =  K ;
nodwghts[K] += nodwghts[bestJ] ;
nzeros[K]   += cost[bestJ] ;
#if MYDEBUG > 0
fprintf(stdout, "\n    setting rep[%d] = %d", bestJ, rep[bestJ]) ;
fprintf(stdout, "\n    setting nodwghts[%d] = %d", K, nodwghts[K]) ;
fprintf(stdout, "\n    setting nzeros[%d] = %d", K, nzeros[K]) ;
#endif
/*
   -------------
   visit K again
   -------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n\n ### visiting front %d", K) ;
   fflush(stdout) ;
#endif
visitAny(K, par, fch, sib, nodwghts, bndwghts, 
         rep, cost, nzeros, maxzeros) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   expand an ETree object by splitting a large front 
   into a chain of smaller fronts.

   created -- 96jun27, cca
   -------------------------------------------------
*/
ETree *
ETree_splitFronts (
   ETree   *etree,
   int     vwghts[],
   int     maxfrontsize,
   int     seed
) {
ETree   *etree2 ;
int     count, front, ii, I, Inew, J, Jnew, nbnd, newsize, nint, nfront,
        nfront2, nsplit, nvtx, prev, size, sizeJ, v, vwght ;
int     *bndwghts, *fch, *head, *indices, *link, *newbndwghts, *newmap, 
        *newnodwghts, *newpar, *nodwghts, *roots, *sib, *vtxToFront ;
Tree    *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || maxfrontsize <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_splitFronts(%p,%p,%d,%d)"
           "\n bad input\n", etree, vwghts, maxfrontsize, seed) ;
   exit(-1) ;
}
tree       = etree->tree ;
fch        = tree->fch ;
sib        = tree->sib ;
nodwghts   = IV_entries(etree->nodwghtsIV) ;
bndwghts   = IV_entries(etree->bndwghtsIV) ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   --------------------------
   set up the working storage
   --------------------------
*/
newpar      = IVinit(nvtx,   -1) ;
roots       = IVinit(nfront, -1) ;
newmap      = IVinit(nvtx,   -1) ;
newnodwghts = IVinit(nvtx,   -1) ;
newbndwghts = IVinit(nvtx,   -1) ;
head        = IVinit(nfront, -1) ;
link        = IVinit(nvtx,   -1) ;
indices     = IVinit(nvtx,   -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   link[v] = head[front] ;
   head[front] = v ;
}
/*
   ------------------------------------------------
   execute a post-order traversal of the front tree
   ------------------------------------------------
*/
nfront2 = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   sizeJ = 0 ;
   for ( v = head[J], count = 0 ; v != -1 ; v = link[v] ) {
      indices[count++] = v ;
      vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
      sizeJ += vwght ;
   }
   if ( sizeJ != nodwghts[J] ) {
      fprintf(stderr, "\n fatal error in ETree_splitFronts(%p,%p,%d,%d)"
             "\n J = %d, sizeJ = %d, nodwght = %d\n", 
             etree, vwghts, maxfrontsize, seed, J, sizeJ, nodwghts[J]) ;
      exit(-1) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n\n checking out front %d, size %d", J, sizeJ) ;
#endif
   if ( sizeJ <= maxfrontsize || fch[J] == -1 ) {
/*
      -------------------------------------------
      this front is small enough (or is a domain)
      -------------------------------------------
*/
      Jnew = nfront2++ ;
      for ( ii = 0 ; ii < count ; ii++ ) {
         v = indices[ii] ;
         newmap[v] = Jnew ;
#if MYDEBUG > 1
            fprintf(stdout, "\n   mapping vertex %d into new front %d",
                    v, Jnew) ;
#endif
      }
      for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
         Inew = roots[I] ;
         newpar[Inew] = Jnew ;
      }
      newnodwghts[Jnew] = nodwghts[J] ;
      newbndwghts[Jnew] = bndwghts[J] ;
      roots[J] = Jnew ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    front is small enough, Jnew = %d", Jnew) ;
#endif
   } else {
/*
      ------------------------------------------
      this front is too large, split into pieces 
      whose size differs by one vertex
      ------------------------------------------
*/
      nsplit  = (sizeJ + maxfrontsize - 1)/maxfrontsize ;
      newsize = sizeJ / nsplit ;
      if ( sizeJ % nsplit != 0 ) {
         newsize++ ;
      }
#if MYDEBUG > 0
      fprintf(stdout, 
         "\n    front is too large, %d target fronts, target size = %d",
         nsplit, newsize) ;
#endif
      prev    = -1 ;
      nint    = nodwghts[J] ;
      nbnd    = nint + bndwghts[J] ;
      if ( seed > 0 ) {
         IVshuffle(count, indices, seed) ;
      }
      ii = 0 ;
      while ( ii < count ) {
         Jnew = nfront2++ ;
         size = 0 ;
         while ( ii < count ) {
            v = indices[ii] ;
            vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
#if MYDEBUG > 0
            fprintf(stdout, 
                "\n   ii = %d, v = %d, vwght = %d, size = %d",
                ii, v, vwght, size) ;
#endif
/*
   -----------------------------------------------
   97aug28, cca
   bug fix. front is created even if it is too big
   -----------------------------------------------
*/
            if ( newsize >= size + vwght || size == 0 ) {
               newmap[v] = Jnew ;
               size += vwght ;
#if MYDEBUG > 0
               fprintf(stdout, 
                "\n   mapping vertex %d into new front %d, size = %d",
                v, Jnew, size) ;
#endif
               ii++ ;
            } else {
               break ;
            }
         }
         if ( prev == -1 ) {
            for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
               Inew = roots[I] ;
               newpar[Inew] = Jnew ;
            }
         } else {
            newpar[prev] = Jnew ;
         }
         prev = Jnew ;
         newnodwghts[Jnew] = size ;
         nbnd = nbnd - size ;
         newbndwghts[Jnew] = nbnd ;
#if MYDEBUG > 0
         fprintf(stdout, "\n    new front %d, size %d, bnd %d",
                 Jnew, newnodwghts[Jnew], newbndwghts[Jnew]) ;
#endif
      }
      roots[J] = Jnew ;
   }
}
/*
   ---------------------------
   create the new ETree object
   ---------------------------
*/
etree2 = ETree_new() ;
ETree_init1(etree2, nfront2, nvtx) ;
IVcopy(nfront2, etree2->tree->par, newpar) ;
Tree_setFchSibRoot(etree2->tree) ;
IVcopy(nvtx, IV_entries(etree2->vtxToFrontIV), newmap) ;
IVcopy(nfront2, IV_entries(etree2->nodwghtsIV), newnodwghts) ;
IVcopy(nfront2, IV_entries(etree2->bndwghtsIV), newbndwghts) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(newpar) ;
IVfree(roots)  ;
IVfree(newmap) ;
IVfree(newnodwghts) ;
IVfree(newbndwghts) ;
IVfree(head) ;
IVfree(link) ;
IVfree(indices) ;

return(etree2) ; }

/*--------------------------------------------------------------------*/
