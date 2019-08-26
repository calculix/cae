/*  stages.c  */

#include "../DSTree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   create a stages vector for a nested dissection ordering

   created -- 96mar10, cca
   -------------------------------------------------------
*/
IV *
DSTree_NDstages (
   DSTree   *dstree
) {
int    ndomsep, nvtx, v ;
int    *hmetric, *map, *stages ;
IV     *hmetricIV, *mapIV, *stagesIV, *vmetricIV ;
Tree   *tree ;
/*
   --------------
   check the data
   --------------
*/
if (  dstree == NULL 
   || (tree = dstree->tree) == NULL
   || (ndomsep = tree->n) < 1
   || (mapIV = dstree->mapIV) == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_NDstages(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
if ( map == NULL || nvtx < 1 ) {
   fprintf(stderr, "\n fatal error in DSTree_NDstages(%p)"
           "\n bad mapIV object\n", dstree) ;
   exit(-1) ;
}
/*
   ----------------------------------
   get the height metric for the tree
   ----------------------------------
*/
vmetricIV = IV_new() ;
IV_init(vmetricIV, ndomsep, NULL) ;
IV_fill(vmetricIV, 1) ;
hmetricIV = Tree_setHeightImetric(tree, vmetricIV) ;
hmetric = IV_entries(hmetricIV) ;
/*
   ------------------------------------
   set the stages for nested dissection
   ------------------------------------
*/
stagesIV = IV_new() ;
IV_init(stagesIV, nvtx, NULL) ;
stages = IV_entries(stagesIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   stages[v] = hmetric[map[v]] - 1 ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(vmetricIV) ;
IV_free(hmetricIV) ;

return(stagesIV) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   create a stages vector for a ``half'' nested dissection ordering

   created -- 96mar10, cca
   ----------------------------------------------------------------
*/
IV *
DSTree_ND2stages (
   DSTree   *dstree
) {
int    ndomsep, nvtx, v ;
int    *hmetric, *map, *stages ;
IV     *hmetricIV, *mapIV, *stagesIV, *vmetricIV ;
Tree   *tree ;
/*
   --------------
   check the data
   --------------
*/
if (  dstree == NULL 
   || (tree = dstree->tree) == NULL
   || (ndomsep = tree->n) < 1
   || (mapIV = dstree->mapIV) == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_ND2stages(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
if ( map == NULL || nvtx < 1 ) {
   fprintf(stderr, "\n fatal error in DSTree_ND2stages(%p)"
           "\n bad mapIV object\n", dstree) ;
   exit(-1) ;
}
/*
   ----------------------------------
   get the height metric for the tree
   ----------------------------------
*/
vmetricIV = IV_new() ;
IV_init(vmetricIV, ndomsep, NULL) ;
IV_fill(vmetricIV, 1) ;
hmetricIV = Tree_setHeightImetric(tree, vmetricIV) ;
hmetric = IV_entries(hmetricIV) ;
/*
   ------------------------------------
   set the stages for nested dissection
   ------------------------------------
*/
stagesIV = IV_new() ;
IV_init(stagesIV, nvtx, NULL) ;
stages = IV_entries(stagesIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   stages[v] = hmetric[map[v]] - 1 ;
   if ( stages[v] > 0 ) {
      stages[v] = (1 + stages[v])/2 ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(vmetricIV) ;
IV_free(hmetricIV) ;

return(stagesIV) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   create a stages vector for a two-level multisection ordering

   created -- 96mar10, cca
   ------------------------------------------------------------
*/
IV *
DSTree_MS2stages (
   DSTree   *dstree
) {
int    ndomsep, nvtx, v ;
int    *hmetric, *map, *stages ;
IV     *hmetricIV, *mapIV, *stagesIV, *vmetricIV ;
Tree   *tree ;
/*
   --------------
   check the data
   --------------
*/
if (  dstree == NULL 
   || (tree = dstree->tree) == NULL
   || (ndomsep = tree->n) < 1
   || (mapIV = dstree->mapIV) == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_MS2stages(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
if ( map == NULL || nvtx < 1 ) {
   fprintf(stderr, "\n fatal error in DSTree_MS2stages(%p)"
           "\n bad mapIV object\n", dstree) ;
   exit(-1) ;
}
/*
   ----------------------------------
   get the height metric for the tree
   ----------------------------------
*/
vmetricIV = IV_new() ;
IV_init(vmetricIV, ndomsep, NULL) ;
IV_fill(vmetricIV, 1) ;
hmetricIV = Tree_setHeightImetric(tree, vmetricIV) ;
hmetric = IV_entries(hmetricIV) ;
/*
   -----------------------------------------
   set the stages for two-level multisection
   -----------------------------------------
*/
stagesIV = IV_new() ;
IV_init(stagesIV, nvtx, NULL) ;
stages = IV_entries(stagesIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   stages[v] = hmetric[map[v]] - 1 ;
   if ( stages[v] > 0 ) {
      stages[v] = 1 ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(vmetricIV) ;
IV_free(hmetricIV) ;

return(stagesIV) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   create a stages vector for a three-level multisection ordering

   created -- 96mar10, cca
   --------------------------------------------------------------
*/
IV *
DSTree_MS3stages (
   DSTree   *dstree
) {
int    ndomsep, nstage, nvtx, v ;
int    *hmetric, *map, *stages ;
IV     *hmetricIV, *mapIV, *stagesIV, *vmetricIV ;
Tree   *tree ;
/*
   --------------
   check the data
   --------------
*/
if (  dstree == NULL 
   || (tree = dstree->tree) == NULL
   || (ndomsep = tree->n) < 1
   || (mapIV = dstree->mapIV) == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_MS3stages(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
if ( map == NULL || nvtx < 1 ) {
   fprintf(stderr, "\n fatal error in DSTree_MS3stages(%p)"
           "\n bad mapIV object\n", dstree) ;
   exit(-1) ;
}
/*
   ----------------------------------
   get the height metric for the tree
   ----------------------------------
*/
vmetricIV = IV_new() ;
IV_init(vmetricIV, ndomsep, NULL) ;
IV_fill(vmetricIV, 1) ;
hmetricIV = Tree_setHeightImetric(tree, vmetricIV) ;
hmetric = IV_entries(hmetricIV) ;
nstage = IV_max(hmetricIV) ;
/*
   -------------------------------------------
   set the stages for three-level multisection
   -------------------------------------------
*/
stagesIV = IV_new() ;
IV_init(stagesIV, nvtx, NULL) ;
stages = IV_entries(stagesIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   stages[v] = hmetric[map[v]] - 1 ;
   if ( stages[v] > 0 ) {
      if ( 2*stages[v] > nstage ) {
         stages[v] = 2 ;
      } else {
         stages[v] = 1 ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(vmetricIV) ;
IV_free(hmetricIV) ;

return(stagesIV) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   create a stages vector via cutoff on domain weights

   created -- 97jun12, cca
   ---------------------------------------------------
*/
IV *
DSTree_stagesViaDomainWeight (
   DSTree   *dstree,
   int      *vwghts,
   DV       *cutoffDV
) {
double   totvwght ;
double   *cutoffs, *nodewghts, *subtreewghts ;
DV       *nodewghtDV, *subtreeDV ;
int      ireg, istage, jstage, ndomsep, nstage, nvtx, v ;
int      *map, *mark, *regmap, *stages ;
IV       *mapIV, *stagesIV ;
Tree     *tree ;
/*
   --------------
   check the data
   --------------
*/
if (  dstree == NULL 
   || (tree = dstree->tree) == NULL
   || (ndomsep = tree->n) < 1
   || (mapIV = dstree->mapIV) == NULL 
   || cutoffDV == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DSTree_stagesViaDomainWeight(%p,%p,%p)"
           "\n bad input\n", dstree, vwghts, cutoffDV) ;
   exit(-1) ;
}
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
if ( map == NULL || nvtx < 1 ) {
   fprintf(stderr, 
           "\n fatal error in DSTree_stagesViaDomainWeight(%p,%p,%p)"
           "\n bad mapIV object\n", dstree, vwghts, cutoffDV) ;
   exit(-1) ;
}
DV_sizeAndEntries(cutoffDV, &nstage, &cutoffs) ;
if ( cutoffs == NULL || nstage < 1 ) {
   fprintf(stderr, 
           "\n fatal error in DSTree_stagesViaDomainWeight(%p,%p,%p)"
           "\n bad cutoffDV object\n", dstree, vwghts, cutoffDV) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n %d stages", nstage) ;
DVfprintf(stdout, nstage, cutoffs) ;
#endif
/*
   --------------------------------
   get the node metric for the tree
   --------------------------------
*/
nodewghtDV = DV_new() ;
DV_init(nodewghtDV, ndomsep, NULL) ;
DV_fill(nodewghtDV, 0.0) ;
nodewghts = DV_entries(nodewghtDV) ;
totvwght  = 0.0 ;
if ( vwghts == NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      nodewghts[map[v]]++ ;
   }
   totvwght = nvtx ;
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      nodewghts[map[v]] += vwghts[v] ;
      totvwght += vwghts[v] ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n node metric") ;
DV_writeForHumanEye(nodewghtDV, stdout) ;
#endif
/*
   ----------------------------------
   get the subtree metric for the tree
   ----------------------------------
*/
subtreeDV = Tree_setSubtreeDmetric(tree, nodewghtDV) ;
subtreewghts = DV_entries(subtreeDV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n subtree metric before scale") ;
DV_writeForHumanEye(subtreeDV, stdout) ;
#endif
for ( ireg = 0 ; ireg < ndomsep ; ireg++ ) {
   subtreewghts[ireg] /= totvwght ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n subtree metric after scale") ;
DV_writeForHumanEye(subtreeDV, stdout) ;
#endif
/*
   ----------------------------
   mark all stages with support
   ----------------------------
*/
mark = IVinit(nstage, -1) ;
for ( ireg = 0 ; ireg < ndomsep ; ireg++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n region %d, subtree weight %.4f",
           ireg, subtreewghts[ireg]) ;
   fflush(stdout) ;
#endif
   for ( istage = 0 ; istage < nstage - 1 ; istage++ ) {
      if (  cutoffs[istage] <= subtreewghts[ireg] 
         && subtreewghts[ireg] < cutoffs[istage+1] ) {
         mark[istage] = 1 ;
#if MYDEBUG > 0
         fprintf(stdout, ", found in stage %d", istage) ;
         fflush(stdout) ;
#endif
         break ;
      }
   }
   if ( istage == nstage - 1 ) {
      mark[nstage - 1] = 1 ;
   }
}
/*
   ------------------
   slide cutoffs down
   ------------------
*/
for ( istage = jstage = 0 ; istage < nstage ; istage++ ) {
   if ( mark[istage] == 1 ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n stage %d marked", istage) ;
      fflush(stdout) ;
#endif
      cutoffs[jstage++] = cutoffs[istage] ;
   }
}
nstage = jstage ;
/*
   ----------------------------------
   get the map from regions to stages
   ----------------------------------
*/
regmap = IVinit(ndomsep, -1) ;
for ( ireg = 0 ; ireg < ndomsep ; ireg++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n region %d, subtree weight %.4f",
           ireg, subtreewghts[ireg]) ;
   fflush(stdout) ;
#endif
   for ( istage = 0 ; istage < nstage - 1 ; istage++ ) {
      if (  cutoffs[istage] <= subtreewghts[ireg] 
         && subtreewghts[ireg] < cutoffs[istage+1] ) {
#if MYDEBUG > 0
         fprintf(stdout, ", found in stage %d", istage) ;
         fflush(stdout) ;
#endif
         regmap[ireg] = istage ;
         break ;
      }
   }
   if ( istage == nstage - 1 ) {
      regmap[ireg] = nstage - 1 ;
#if MYDEBUG > 0
      fprintf(stdout, ", found in stage %d", nstage - 1) ;
      fflush(stdout) ;
#endif
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n region to stage map") ;
IVfp80(stdout, ndomsep, regmap, 80, &ierr) ;
#endif
/*
   --------------
   set the stages
   --------------
*/
stagesIV = IV_new() ;
IV_init(stagesIV, nvtx, NULL) ;
stages = IV_entries(stagesIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   stages[v] = regmap[map[v]] ;
}
#if MYDEBUG > 0
stageWeights = DVinit(nstage, 0.0) ;
for ( ireg = 0 ; ireg < ndomsep ; ireg++ ) {
   stageWeights[regmap[ireg]] += nodewghts[ireg] ;
}
fprintf(stdout, "\n\n stageWeights, sum = %.4f", 
        DVsum(nstage, stageWeights)) ;
DVfprintf(stdout, nstage, stageWeights) ;
for ( istage = nstage - 2 ; istage >= 0 ; istage-- ) {
   stageWeights[istage] += stageWeights[istage+1] ;
}
fprintf(stdout, "\n\n stageWeights, sum = %.4f", 
        DVsum(nstage, stageWeights)) ;
for ( istage = 0 ; istage < nstage ; istage++ ) {
   stageWeights[istage] /= totvwght ;
}
fprintf(stdout, "\n\n stageWeights, sum = %.4f", 
        DVsum(nstage, stageWeights)) ;
DVfprintf(stdout, nstage, stageWeights) ;
DVfree(stageWeights) ;
#endif
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(nodewghtDV) ;
DV_free(subtreeDV) ;
IVfree(regmap) ;
IVfree(mark) ;

return(stagesIV) ; }

/*--------------------------------------------------------------------*/
