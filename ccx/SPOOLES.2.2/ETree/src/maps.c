/*  maps.c  */

#include "../ETree.h"

#define MYDEBUG  0

#define SUBTREE_SIZE   1
#define SUBTREE_OPS    2
#define SUBTREE_DEFINITION SUBTREE_OPS

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a wrap map of the front tree.

   created -- 96dec12, cca
   --------------------------------------------------------------
*/
IV *
ETree_wrapMap (
   ETree   *frontETree,
   int     type,
   int     symflag,
   DV      *cumopsDV
) {
double   *cumops, *forwardOps ;
DV       *forwardOpsDV ;
int      jthread, J, nfront, nthread ;
int      *owners ;
IV       *ownersIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || cumopsDV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_wrapMap(%p,%p)"
           "\n bad input\n", frontETree, cumopsDV) ;
   exit(-1) ;
}
tree = frontETree->tree ;
DV_sizeAndEntries(cumopsDV, &nthread, &cumops) ;
DV_zero(cumopsDV) ;
/*
   ---------------------------------
   get a vector of forward op counts
   ---------------------------------
*/
forwardOpsDV = ETree_forwardOps(frontETree, type, symflag) ;
DV_sizeAndEntries(forwardOpsDV, &nfront, &forwardOps) ;
/*
   -------------------
   load the map vector
   -------------------
*/
ownersIV = IV_new() ;
IV_init(ownersIV, nfront, NULL) ;
owners = IV_entries(ownersIV) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   jthread = J % nthread ;
   owners[J] = jthread ;
   cumops[jthread] += forwardOps[J] ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(forwardOpsDV) ;

return(ownersIV) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a balanced map of the front tree.
   the fronts are visited in the post-order traversal.

   created -- 96dec12, cca
   ----------------------------------------------------------------
*/
IV *
ETree_balancedMap (
   ETree   *frontETree,
   int     type,
   int     symflag,
   DV      *cumopsDV
) {
double   minops ;
double   *cumops, *forwardOps ;
DV       *forwardOpsDV ;
int      ithread, jthread, J, nfront, nthread ;
int      *owners ;
IV       *ownersIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || cumopsDV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_balancedMap(%p,%p)"
           "\n bad input\n", frontETree, cumopsDV) ;
   exit(-1) ;
}
tree = frontETree->tree ;
DV_sizeAndEntries(cumopsDV, &nthread, &cumops) ;
DV_zero(cumopsDV) ;
/*
   ---------------------------------
   get a vector of forward op counts
   ---------------------------------
*/
forwardOpsDV = ETree_forwardOps(frontETree, type, symflag) ;
DV_sizeAndEntries(forwardOpsDV, &nfront, &forwardOps) ;
/*
   -------------------
   load the map vector
   -------------------
*/
ownersIV = IV_new() ;
IV_init(ownersIV, nfront, NULL) ;
owners = IV_entries(ownersIV) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   jthread = 0 ;
   minops  = cumops[0] ;
   for ( ithread = 1 ; ithread < nthread ; ithread++ ) {
      if ( minops > cumops[ithread] ) {
         jthread = ithread ;
         minops  = cumops[ithread] ;
      }
   }
   owners[J] = jthread ;
   cumops[jthread] += forwardOps[J] ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(forwardOpsDV) ;

return(ownersIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   this method constructs and returns an IV object 
   that holds the map from fronts to threads for a 
   "subtree-subset" map of the front tree.

   created -- 97jan15, cca
   -----------------------------------------------
*/
IV *
ETree_subtreeSubsetMap (
   ETree   *frontETree,
   int     type,
   int     symflag,
   DV      *cumopsDV
) {
double   offset, total ;
double   *cumops, *forwardOps, *tmetric ;
DV       *forwardOpsDV, *tmetricDV ;
int      I, J, mthread, nfront, nthread, q, qmin ;
int      *fch, *firsts, *lasts, *owners, *par, *sib ;
IV       *ownersIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || cumopsDV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_subtreeSubsetMap(%p,%p)"
           "\n bad input\n", frontETree, cumopsDV) ;
   exit(-1) ;
}
tree = frontETree->tree ;
fch  = tree->fch ;
par  = tree->par ;
sib  = tree->sib ;
DV_sizeAndEntries(cumopsDV, &nthread, &cumops) ;
DV_zero(cumopsDV) ;
/*
   ---------------------------------
   get a vector of forward op counts
   ---------------------------------
*/
forwardOpsDV = ETree_forwardOps(frontETree, type, symflag) ;
DV_sizeAndEntries(forwardOpsDV, &nfront, &forwardOps) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n forwardOpsDV") ;
DV_writeForHumanEye(forwardOpsDV, stdout) ;
fflush(stdout) ;
#endif
/*
   --------------------------------
   get the subtree metric DV object
   --------------------------------
*/
tmetricDV = Tree_setSubtreeDmetric(tree, forwardOpsDV) ;
tmetric = DV_entries(tmetricDV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n tmetricDV") ;
DV_writeForHumanEye(tmetricDV, stdout) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------
   left justify the tree to make the 
   assignment algorithm work correctly
   -----------------------------------
*/
ETree_leftJustifyD(frontETree, tmetricDV) ;
/*
   -----------------------------------------------------
   fill two vectors that hold the first and last threads
   that are eligible to own a front 
   -----------------------------------------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n\n pre-order traversal to determine eligible sets") ;
fflush(stdout) ;
#endif
firsts = IVinit(nfront, -1) ;
lasts  = IVinit(nfront, -1) ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n visiting front %d", J) ;
   fflush(stdout) ;
#endif
   if ( par[J] == -1 ) {
      firsts[J] = 0 ;
      lasts[J]  = nthread - 1 ;
#if MYDEBUG > 0
      fprintf(stdout, ", root front") ;
      fflush(stdout) ;
#endif
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n first = %d, last = %d", firsts[J], lasts[J]) ;
   fflush(stdout) ;
#endif
   if ( fch[J] != -1 ) {
      mthread = lasts[J] - firsts[J] + 1 ;
      total   = tmetric[J] - forwardOps[J] ;
#if MYDEBUG > 0
      fprintf(stdout, "\n mthread = %d, total = %.0f", mthread, total) ;
      fflush(stdout) ;
#endif
      for ( I = fch[J], offset = 0.0 ; I != -1 ; I = sib[I] ) {
         firsts[I] = firsts[J] + (int) (mthread*offset/total) ;
#if MYDEBUG > 0
         fprintf(stdout, "\n child %d, offset = %.0f, firsts[%d] = %d",
                 I, offset, I, firsts[I]) ;
         fflush(stdout) ;
#endif
         offset += tmetric[I] ;
         lasts[I] = firsts[J] + (int) (mthread*offset/total) - 1 ;
         if ( lasts[I] < firsts[I] ) {
            lasts[I] = firsts[I] ;
         }
#if MYDEBUG > 0
         fprintf(stdout, "\n child %d, offset = %.0f, lasts[%d] = %d",
                 I, offset, I, lasts[I]) ;
         fflush(stdout) ;
#endif
      }
   }
}
/*
   ---------------------------------------------------------------
   now fill the map IV object and cumops[*] vector with a
   balanced map using the candidate sets via a postorder traversal
   ---------------------------------------------------------------
*/
ownersIV = IV_new() ;
IV_init(ownersIV, nfront, NULL) ;
owners = IV_entries(ownersIV) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n front %d, firsts[%d] = %d, lasts[%d] = %d",
           J, J, firsts[J], J, lasts[J]) ;
   fflush(stdout) ;
#endif
   qmin = firsts[J] ;
   for ( q = firsts[J] + 1 ; q <= lasts[J] ; q++ ) {
      if ( cumops[qmin] > cumops[q] ) {
         qmin = q ;
      }
   }
   owners[J] = qmin ;
   cumops[qmin] += forwardOps[J] ;
#if MYDEBUG > 0
   fprintf(stdout, ", owners[%d] = %d, cumops[%d] = %.0f",
           J, owners[J], qmin, cumops[qmin]) ;
   fflush(stdout) ;
#endif
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(forwardOpsDV) ;
DV_free(tmetricDV) ;
IVfree(firsts) ;
IVfree(lasts) ;

return(ownersIV) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a domain decomposition 
   balanced map of the front tree.
   the domains are mapped to threads using a balanced map,
   and the schur complement fronts are mapped to threads 
   using a balanced map, but the two balanced maps are independent.

   created -- 97jan17, cca
   ----------------------------------------------------------------
*/
IV *
ETree_ddMap (
   ETree    *frontETree,
   int     type,
   int     symflag,
   DV       *cumopsDV,
   double   cutoff
) {
double   minops ;
double   *cumops, *domainops, *forwardOps, *schurops, *tmetric ;
DV       *forwardOpsDV, *tmetricDV ;
int      ithread, jthread, J, K, ndom, nfront, nthread, root ;
int      *ms, *owners, *par, *rootAnc ;
IV       *msIV, *ownersIV, *rootAncIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || cumopsDV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_ddMap(%p,%p,%f)"
           "\n bad input\n", frontETree, cumopsDV, cutoff) ;
   exit(-1) ;
}
nfront = frontETree->nfront ;
tree   = frontETree->tree ;
par    = tree->par ;
DV_sizeAndEntries(cumopsDV, &nthread, &cumops) ;
DV_zero(cumopsDV) ;
/*
   ---------------------------------
   get a vector of forward op counts
   ---------------------------------
*/
forwardOpsDV = ETree_forwardOps(frontETree, type, symflag) ;
DV_sizeAndEntries(forwardOpsDV, &nfront, &forwardOps) ;
#if MYDEBUG > 0
fprintf(stdout, "\n forwardOpsDV") ;
DV_writeForHumanEye(forwardOpsDV, stdout) ;
fflush(stdout) ;
#endif
#if SUBTREE_DEFINITION == SUBTREE_SIZE
/*
   -----------------------------
   get a vector of subtree sizes
   -----------------------------
*/
{ 
DV   *tempDV ;
IV   *tempIV ;
tempIV = ETree_nvtxMetric(frontETree) ;
fprintf(stdout, "\n\n nvtx metric") ;
IV_writeForHumanEye(tempIV, stdout) ;
tempDV = DV_new() ;
for ( J = 0 ; J < nfront ; J++ ) {
   DV_setEntry(tempDV, J, (double) IV_entry(tempIV, J)) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n double nvtx metric") ;
DV_writeForHumanEye(tempDV, stdout) ;
#endif
tmetricDV = Tree_setSubtreeDmetric(tree, tempDV) ;
IV_free(tempIV) ;
DV_free(tempDV) ;
}
#endif
#if SUBTREE_DEFINITION == SUBTREE_OPS
tmetricDV = Tree_setSubtreeDmetric(tree, forwardOpsDV) ;
#endif
/*
   ------------------------
   get a multisector vector
   ------------------------
*/
msIV = IV_new() ;
IV_init(msIV, nfront, NULL) ;
IV_fill(msIV, 0) ;
ms = IV_entries(msIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n double nvtx subtree metric") ;
DV_writeForHumanEye(tmetricDV, stdout) ;
#endif
tmetric   = DV_entries(tmetricDV) ;
cutoff = cutoff * DV_max(tmetricDV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( tmetric[J] < cutoff ) {
      ms[J] = 1 ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n msIV") ;
IV_writeForHumanEye(msIV, stdout) ;
fflush(stdout) ;
#endif
/*
   --------------------------------------------
   create a rootAnc vector, 
   if J is in a domain then
      rootAnc[J] is the root node of the domain
   else
      rootAnc[J] is the root node of the tree
   endif
   --------------------------------------------
*/
rootAncIV = IV_new() ;
IV_init(rootAncIV, nfront, NULL) ;
rootAnc   = IV_entries(rootAncIV) ;
for ( J = nfront - 1 ; J >= 0 ; J-- ) {
   if ( (K = par[J]) == -1 || ms[J] != ms[K] ) {
      rootAnc[J] = J ;
   } else {
      rootAnc[J] = rootAnc[K] ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n rootAncIV") ;
IV_writeForHumanEye(rootAncIV, stdout) ;
fflush(stdout) ;
#endif
/*
   ------------------------------
   initialize the ownersIV object
   ------------------------------
*/
ownersIV = IV_new() ;
IV_init(ownersIV, nfront, NULL) ;
owners = IV_entries(ownersIV) ;
/*
   --------------------------------------------------
   assign the domains to threads using a balanced map
   --------------------------------------------------
*/
domainops = DVinit(nthread, 0.0) ;
root = -1 ;
ndom =  0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( ms[J] == 1 ) {
      if ( root != rootAnc[J] ) {
         ndom++ ;
         root    = rootAnc[J] ;
         jthread = 0 ;
         minops  = domainops[0] ;
         for ( ithread = 1 ; ithread < nthread ; ithread++ ) {
            if ( minops > domainops[ithread] ) {
               jthread = ithread ;
               minops  = domainops[ithread] ;
            }
         }
      }
      owners[J] = jthread ;
      domainops[jthread] += forwardOps[J] ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n %d domains", ndom) ;
fprintf(stdout, "\n domainops") ;
DVfprintf(stdout, nthread, domainops) ;
fflush(stdout) ;
#endif
/*
   ------------------------------------------------------------------
   assign the schur complement fronts to threads using a balanced map
   ------------------------------------------------------------------
*/
schurops = DVinit(nthread, 0.0) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( ms[J] == 0 ) {
      jthread = 0 ;
      minops  = schurops[0] ;
      for ( ithread = 1 ; ithread < nthread ; ithread++ ) {
         if ( minops > schurops[ithread] ) {
            jthread = ithread ;
            minops  = schurops[ithread] ;
         }
      }
      owners[J] = jthread ;
      schurops[jthread] += forwardOps[J] ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n schurops") ;
DVfprintf(stdout, nthread, schurops) ;
fflush(stdout) ;
#endif
/*
   -------------------------------------
   fill the cumulative operations vector
   -------------------------------------
*/
for ( jthread = 0 ; jthread < nthread ; jthread++ ) {
   cumops[jthread] = domainops[jthread] + schurops[jthread] ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n cumops") ;
DVfprintf(stdout, nthread, cumops) ;
fflush(stdout) ;
#endif
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(forwardOpsDV) ;
DV_free(tmetricDV) ;
IV_free(msIV) ;
IV_free(rootAncIV) ;
DVfree(domainops) ;
DVfree(schurops) ;

return(ownersIV) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a domain decomposition 
   balanced map of the front tree.
   the domains are mapped to threads using a balanced map,
   and the schur complement fronts are mapped to threads 
   using a balanced map, but the two balanced maps are independent.

   created -- 97jan17, cca
   ----------------------------------------------------------------
*/
IV *
ETree_ddMapNew (
   ETree   *frontETree,
   int     type,
   int     symflag,
   IV      *msIV,
   DV      *cumopsDV
) {
double   minops ;
double   *cumops, *domprios, *domops, *forwardOps, 
         *schurprios, *schurops ;
DV       *forwardOpsDV ;
int      domid, idom, ireg, ithread, jthread, J, K, nbndJ, 
         ndom, nfront, nJ, nschur, nthread, 
         nvtx, v ;
int      *bndwghts, *domids, *dommap, *frontToRegion, *ms, 
         *nodwghts, *owners, *par, *schurids, *vtxToFront ;
IV       *ownersIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || cumopsDV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_ddMapNew(%p,%p,%p)"
           "\n bad input\n", frontETree, msIV, cumopsDV) ;
   exit(-1) ;
}
nfront     = ETree_nfront(frontETree) ;
nvtx       = ETree_nvtx(frontETree) ;
tree       = ETree_tree(frontETree) ;
vtxToFront = ETree_vtxToFront(frontETree) ;
nodwghts   = ETree_nodwghts(frontETree) ;
bndwghts   = ETree_bndwghts(frontETree) ;
par        = ETree_par(frontETree) ;
DV_sizeAndEntries(cumopsDV, &nthread, &cumops) ;
DV_zero(cumopsDV) ;
ms = IV_entries(msIV) ;
ownersIV = IV_new() ;
IV_init(ownersIV, nfront, NULL) ;
owners = IV_entries(ownersIV) ;
/*
   --------------------------------------
   compute the frontToRegion[] map vector
   --------------------------------------
*/
frontToRegion = IVinit(nfront, -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   J = vtxToFront[v] ;
   frontToRegion[J] = ms[v] ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n frontToRegion[] from msIV") ;
IVfprintf(stdout, nfront, frontToRegion) ;
#endif
ndom = 0 ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   if ( frontToRegion[J] != 0 ) {
      if ( (K = par[J]) != -1 ) {
         if ( frontToRegion[K] == 0 ) {
            frontToRegion[J] = ++ndom ;
         } else {
            frontToRegion[J] = frontToRegion[K] ;
         }
      } else {
         frontToRegion[J] = ++ndom ;
      }
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n frontToRegion[], separate domains") ;
IVfprintf(stdout, nfront, frontToRegion) ;
#endif
/*
   ---------------------------------
   get a vector of forward op counts
   ---------------------------------
*/
forwardOpsDV = ETree_forwardOps(frontETree, type, symflag) ;
forwardOps = DV_entries(forwardOpsDV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n forwardOpsDV") ;
DV_writeForHumanEye(forwardOpsDV, stdout) ;
#endif
/*
   ----------------------------------------
   for each domain, compute the forward ops
   ----------------------------------------
*/
domprios = DVinit(ndom+1, 0.0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (ireg = frontToRegion[J]) > 0 ) {
     domprios[ireg] += forwardOps[J] ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n domprios[], sum = %.0f",
        DVsum(ndom, domprios+1)) ;
DVfprintf(stdout, ndom+1, domprios) ;
#endif
/*
   --------------------------------------------------------
   sort the domains by descending order of their operations
   --------------------------------------------------------
*/
domids = IVinit(ndom, -1) ;
IVramp(ndom, domids, 1, 1) ;
#if MYDEBUG > 0
fprintf(stdout, "\n before sort: domids") ;
IVfprintf(stdout, ndom, domids) ;
fprintf(stdout, "\n before sort: domprios") ;
DVfprintf(stdout, ndom, domprios+1) ;
#endif
DVIVqsortDown(ndom, domprios + 1, domids) ;
#if MYDEBUG > 0
fprintf(stdout, "\n after sort: domids") ;
IVfprintf(stdout, ndom, domids) ;
fprintf(stdout, "\n after sort: domprios") ;
DVfprintf(stdout, ndom, domprios+1) ;
#endif
/*
   ----------------------------
   assign domains to processors
   ----------------------------
*/
dommap = IVinit(ndom+1, -1) ;
domops = DVinit(nthread, 0.0) ;
for ( idom = 0 ; idom < ndom ; idom++ ) {
   domid = domids[idom] ;
   jthread = 0 ;
   minops = domops[0] ;
   for ( ithread = 1 ; ithread < nthread ; ithread++ ) {
      if ( domops[ithread] < minops ) {
         jthread = ithread ;
         minops  = domops[ithread] ;
      }
   }
   dommap[domid] = jthread ;
   domops[jthread] += domprios[idom+1] ;
#if MYDEBUG > 0
   fprintf(stdout, 
  "\n assigning domain %d with ops %.0f to thread %d, new ops = %.0f", 
           domid, domprios[idom+1], jthread, domops[jthread]) ;
#endif
}
#if MYDEBUG > 0
fprintf(stdout, "\n domops[]") ;
DVfprintf(stdout, nthread, domops) ;
#endif
/*
   -------------------------------------
   assign fronts in domains to processes
   -------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (idom = frontToRegion[J]) > 0 ) {
      owners[J] = dommap[idom] ;
   }
}
/*
   -----------------------------------------------------
   now compute priorities of the schur complement fronts
   -----------------------------------------------------
*/
schurprios = DVinit(nfront, 0.0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontToRegion[J] == 0 ) {
      nJ    = nodwghts[J] ;
      nbndJ = bndwghts[J] ;
      schurprios[J] = nJ*nJ*(nJ + nbndJ) ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n schur front ops") ;
DVfprintf(stdout, nfront, schurprios) ;
#endif
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   if ( frontToRegion[J] == 0 ) {
      if ( (K = par[J]) != -1 ) {
         schurprios[J] += schurprios[K] ;
      }
   }
}
/*
   --------------------------------------------------
   sort the priorities of the schur complement fronts
   in descending order of their priorities
   --------------------------------------------------
*/
schurids = IVinit(nfront, -1) ;
for ( J = nschur = 0 ; J < nfront ; J++ ) {
   if ( frontToRegion[J] == 0 ) {
      schurids[nschur]   = J ;
      schurprios[nschur] = schurprios[J] ;
      nschur++ ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n before sort: schurids") ;
IVfprintf(stdout, nschur, schurids) ;
fprintf(stdout, "\n before sort: schurprios") ;
DVfprintf(stdout, nschur, schurprios) ;
#endif
DVIVqsortDown(nschur, schurprios, schurids) ;
#if MYDEBUG > 0
fprintf(stdout, "\n after sort: schurids") ;
IVfprintf(stdout, nschur, schurids) ;
fprintf(stdout, "\n after sort: schurprios") ;
DVfprintf(stdout, nschur, schurprios) ;
#endif
/*
   ---------------------------------
   assign schur fronts to processors
   ---------------------------------
*/
schurops = DVinit(nthread, 0.0) ;
for ( ireg = 0 ; ireg < nschur ; ireg++ ) {
   J = schurids[ireg] ;
   jthread = 0 ;
   minops = schurops[0] ;
   for ( ithread = 1 ; ithread < nthread ; ithread++ ) {
      if ( schurops[ithread] < minops ) {
         jthread = ithread ;
         minops  = schurops[ithread] ;
      }
   }
   owners[J] = jthread ;
   schurops[jthread] += forwardOps[J] ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n assigning schur front %d to thread %d, new ops = %.0f", 
           J, jthread, schurops[jthread]) ;
#endif
}
/*
   ------------------------------
   fill the cumulative ops vector
   ------------------------------
*/
for ( jthread = 0 ; jthread < nthread ; jthread++ ) {
   cumops[jthread] = domops[jthread] + schurops[jthread] ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n domops[]") ;
DVfprintf(stdout, nthread, domops) ;
fprintf(stdout, "\n schurops[]") ;
DVfprintf(stdout, nthread, schurops) ;
fprintf(stdout, "\n sum(domops) = %.0f", DVsum(nthread, domops)) ;
fprintf(stdout, "\n sum(schurops) = %.0f", DVsum(nthread, schurops)) ;
fprintf(stdout, "\n sum(cumops) = %.0f", DV_sum(cumopsDV)) ;
#endif
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(frontToRegion) ;
IVfree(domids) ;
IVfree(dommap) ;
IVfree(schurids) ;
DV_free(forwardOpsDV) ;
DVfree(domprios) ;
DVfree(domops) ;
DVfree(schurprios) ;
DVfree(schurops) ;

return(ownersIV) ; }

/*--------------------------------------------------------------------*/
