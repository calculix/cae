/*  init.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   initialize the object given the 
   number of fronts and vertices

   created -- 96jun23, cca
   -------------------------------
*/
void
ETree_init1 (
   ETree   *etree,
   int     nfront,
   int     nvtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || nfront < 0 || nvtx < nfront ) {
   fprintf(stderr, "\n fatal error in ETree_init1(%p,%d,%d)"
           "\n bad input\n", etree, nfront, nvtx) ;
   exit(-1) ;
}
ETree_clearData(etree) ;
etree->nfront = nfront ;
etree->nvtx   = nvtx   ;
etree->tree   = Tree_new() ;
Tree_init1(etree->tree, nfront) ;
etree->nodwghtsIV = IV_new() ;
IV_init(etree->nodwghtsIV, nfront, NULL) ;
IV_fill(etree->nodwghtsIV, 0) ;
etree->bndwghtsIV = IV_new() ;
IV_init(etree->bndwghtsIV, nfront, NULL) ;
IV_fill(etree->bndwghtsIV, 0) ;
etree->vtxToFrontIV = IV_new() ;
IV_init(etree->vtxToFrontIV, nvtx, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   initialize the ETree object from a graph

   created -- 95nov15, cca
   ----------------------------------------
*/
void
ETree_initFromGraph (
   ETree   *etree,
   Graph   *g
) {
int   ii, nfront, nvtx, u, v, vsize ;
int   *bndwghts, *mark, *nodwghts, *par, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || g == NULL || (nvtx = g->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_initFromGraph(%p,%p)"
           "\n bad input\n", etree, g) ;
   exit(-1) ;
}
/*
   -----------------
   set up the object
   -----------------
*/
nfront = nvtx ;
ETree_init1(etree, nfront, nvtx) ;
nodwghts = IV_entries(etree->nodwghtsIV) ;
bndwghts = IV_entries(etree->bndwghtsIV) ;
par      = etree->tree->par ;
IV_ramp(etree->vtxToFrontIV, 0, 1) ;
/*
   --------------------------------------------------------
   fill the parent, node weight and boundary weight vectors
   --------------------------------------------------------
*/
if ( g->vwghts == NULL ) {
   IVfill(nvtx, nodwghts, 1) ;
} else {
   IVcopy(nvtx, nodwghts, g->vwghts) ;
}
mark = IVinit(nvtx, -1) ;
IVramp(nvtx, mark, 0, 1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   Graph_adjAndSize(g, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      u = vadj[ii] ;
      while ( u < v && mark[u] != v ) {
         bndwghts[u] += nodwghts[v] ;
         if ( mark[u] == u ) {
            par[u] = v ;
         }
         mark[u] = v ;
         u = par[u] ;
      }
   }
}
IVfree(mark) ;
/*
   ------------------------------------
   set the fch[], sub[] and root fields
   ------------------------------------
*/
Tree_setFchSibRoot(etree->tree) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   initialize the ETree object from a graph and permutation

   created -- 95nov15, cca
   --------------------------------------------------------
*/
void
ETree_initFromGraphWithPerms (
   ETree   *etree,
   Graph   *g,
   int     newToOld[],
   int     oldToNew[]
) {
int   ii, nfront, nvtx, unew, uold, vold, vnew, vsize ;
int   *bndwghts, *mark, *nodwghts, *par, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || g == NULL || (nvtx = g->nvtx) <= 0 
   || newToOld == NULL || oldToNew == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_initFromGraph(%p,%p,%p,%p)"
           "\n bad input\n", etree, g, newToOld, oldToNew) ;
   exit(-1) ;
}
nfront = nvtx ;
/*
   --------------------------------------------------
   check that the permutation is really a permutation
   --------------------------------------------------
*/
for ( unew = 0 ; unew < nvtx ; unew++ ) {
   if (    (uold = newToOld[unew]) < 0 
        || uold >= nvtx 
        || oldToNew[uold] != unew ) {
      fprintf(stderr, 
           "\n fatal error in ETree_initFromGraphWithPerms(%p,%p,%p,%p)"
           "\n uold = %d, unew = %d",
         etree, g, newToOld, oldToNew, uold, unew) ;
      if ( 0 <= uold && uold < nvtx ) {
         fprintf(stderr, "\n oldToNew[%d] = %d", uold, oldToNew[uold]) ;
      }
      if ( 0 <= unew && unew < nvtx ) {
         fprintf(stderr, "\n newToOld[%d] = %d", unew, newToOld[unew]) ;
      }
      exit(-1) ;
   }
}
/*
   -----------------
   set up the object
   -----------------
*/
nfront = nvtx ;
ETree_init1(etree, nfront, nvtx) ;
nodwghts = IV_entries(etree->nodwghtsIV) ;
bndwghts = IV_entries(etree->bndwghtsIV) ;
par      = etree->tree->par ;
IVcopy(nvtx, IV_entries(etree->vtxToFrontIV), oldToNew) ;
/*
   --------------------------------------------------------
   fill the parent, node weight and boundary weight vectors
   --------------------------------------------------------
*/
if ( g->vwghts == NULL ) {
   IVfill(nvtx, nodwghts, 1) ;
} else {
   for ( vold = 0 ; vold < nvtx ; vold++ ) {
      nodwghts[oldToNew[vold]] = g->vwghts[vold] ;
   }
}
mark = IVinit(nvtx, -1) ;
IVramp(nvtx, mark, 0, 1) ;
for ( vnew = 0 ; vnew < nvtx ; vnew++ ) {
   vold = newToOld[vnew] ;
   Graph_adjAndSize(g, vold, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      uold = vadj[ii] ;
      if ( uold < nvtx ) {
         unew = oldToNew[uold] ;
         while ( unew < vnew && mark[unew] != vnew ) {
            bndwghts[unew] += nodwghts[vnew] ;
            if ( mark[unew] == unew ) {
               par[unew] = vnew ;
            }
            mark[unew] = vnew ;
            unew = par[unew] ;
         }
      }
   }
}
IVfree(mark) ;
/*
   ------------------------------------
   set the fch[], sub[] and root fields
   ------------------------------------
*/
Tree_setFchSibRoot(etree->tree) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- initialize the front tree for a dense matrix
   
   n -- size of the matrix
   option -- mapping option
      1 --> have all fronts (save the last) contain the same
            number of vertices
      2 --> have all fronts have roughly equal numbers of entries

   created -- 96aug19, cca
   --------------------------------------------------------------
*/
void
ETree_initFromDenseMatrix (
   ETree   *etree,
   int     n,
   int     option,
   int     param
) {
int   b, bnd, first, front, ii, last, nent, nfront, target ;
int   *bndwghts, *nodwghts, *par, *vtxToFront ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL || n <= 0 
   || option < 1 || option > 2 || param <= 0 ) {
   fprintf(stderr, 
           "\n fatal error in ETree_initFromDenseMatrix(%p,%d,%d,%d)"
           "\n bad input\n", etree, n, option, param) ;
   exit(-1) ;
}
ETree_clearData(etree) ;
if ( option == 1 ) {
/*
   ------------------------------------------------------------------  
   first option, let all front be size b except for possibly the root
   ------------------------------------------------------------------  
*/
   b = param ;
   nfront = (n + b - 1)/b ;
   ETree_init1(etree, nfront, n) ;
   nodwghts   = IV_entries(etree->nodwghtsIV) ;
   bndwghts   = IV_entries(etree->bndwghtsIV) ;
   vtxToFront = IV_entries(etree->vtxToFrontIV) ;
   for ( ii = 0 ; ii < n ; ii++ ) {
      vtxToFront[ii] = ii / b ;
   }
   for ( ii = 0, bnd = n ; ii < nfront ; ii++ ) {
      nodwghts[ii] = (b <= bnd) ? b : bnd ;
      bnd -= nodwghts[ii] ;
      bndwghts[ii] = bnd ;
   }
} else if ( option == 2 ) {
/*
   ------------------------------------------------------------
   second option, let each front have the same number of entries
   ------------------------------------------------------------
*/
   target = param ;
   first  = 0 ;
   bnd    = n - 1 ;
   nfront = 0 ;
   while ( first < n ) {
      nent = 2*(n - first) - 1 ;
      last = first + 1 ;
      while ( last < n && (nent + 2*(n - last) - 1) <= target ) {
         nent += 2*(n - last) - 1 ;
         last++ ;
      }
      first = last ;
      nfront++ ;
   }
   ETree_init1(etree, nfront, n) ;
   nodwghts   = IV_entries(etree->nodwghtsIV) ;
   bndwghts   = IV_entries(etree->bndwghtsIV) ;
   vtxToFront = IV_entries(etree->vtxToFrontIV) ;
   first  = 0 ;
   bnd    = n - 1 ;
   front  = 0 ;
   while ( first < n ) {
      nent = 2*(n - first) - 1 ;
      vtxToFront[first] = front ;
      last = first + 1 ;
      while ( last < n && (nent + 2*(n - last) - 1) <= target ) {
         vtxToFront[last] = front ;
         nent += 2*(n - last) - 1 ;
         last++ ;
      }
      last-- ;
      fprintf(stdout, "\n front = %d, first = %d, last = %d, nent = %d",
              front, first, last, nent) ;
      nodwghts[front] = last - first + 1 ;
      bndwghts[front] = n - last - 1 ;
      first = last + 1 ;
      front++ ;
   }
}
par = etree->tree->par ;
IVramp(nfront, par, 1, 1) ;
par[nfront-1] = -1 ;
Tree_setFchSibRoot(etree->tree) ;

{
int      bndJ, count, ierr, I, J, sizeI, sizeJ ;
int      *tmp ;
double   facops, solops, updops ;

facops = solops = updops = 0.0 ;
tmp = IVinit((nfront*(nfront+1))/2, -1) ;
count = 0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   sizeJ = nodwghts[J] ;
   bndJ  = bndwghts[J] ;
   facops += 2*(sizeJ * sizeJ * sizeJ)/3 ;
   solops += 2 * sizeJ * sizeJ * bndJ ;
   tmp[count++] = facops + solops ;
   for ( I = 0 ; I < J ; I++ ) {
      sizeI = nodwghts[I] ;
      updops += 2 * sizeI * sizeJ * (sizeJ + 2*bndJ) ;
      tmp[count++] = updops ;
   }
}
IVqsortUp(count, tmp) ;
/*
IVfp80(stdout, count, tmp, 80, &ierr) ;
*/
IVfree(tmp) ;
fprintf(stdout, 
        "\n factor ops = %.0f, %5.1f per cent of total"
        "\n solve ops  = %.0f, %5.1f per cent of total"
        "\n update ops = %.0f, %5.1f per cent of total",
        facops, 100.*facops/(facops + solops + updops),
        solops, 100.*solops/(facops + solops + updops),
        updops, 100.*updops/(facops + solops + updops)) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   initialize the ETree object
   (1) read ETree object from file
   (2) get the old-to-new permutation
   (3) permute the ETree
   (4) return the old-to-new permutation
 
   created -- 97jul13, cca
   -------------------------------------
*/
IV *
ETree_initFromFile (
   ETree   *frontETree,
   char    *inETreeFileName,
   int     msglvl,
   FILE    *msgFile
) {
double   t1, t2 ;
int      neqns, rc ;
IV       *oldToNewIV ;
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
if ( strcmp(inETreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
ETree_setDefaultFields(frontETree) ;
rc = ETree_readFromFile(frontETree, inETreeFileName) ;
MARKTIME(t2) ;
neqns = frontETree->nvtx   ;
fprintf(msgFile, "\n CPU %8.3f : read in frontETree from file %s",
        t2 - t1, inETreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, frontETree, inETreeFileName) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after reading ETree object from file %s",
           inETreeFileName) ;
   if ( msglvl == 2 ) {
      ETree_writeStats(frontETree, msgFile) ;
   } else {
      ETree_writeForHumanEye(frontETree, msgFile) ;
   }
}
fflush(msgFile) ;
/*
   -----------------------------
   get the permutation IV object
   -----------------------------
*/
MARKTIME(t1) ;
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : get permutation", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n vertex old-to-new IV object") ;
   if ( msglvl == 2 ) {
      IV_writeStats(oldToNewIV, msgFile) ;
   } else {
      IV_writeForHumanEye(oldToNewIV, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   permute the vertices in the ETree object
   ----------------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n before permuting the vertex map") ;
   if ( msglvl == 2 ) {
      ETree_writeStats(frontETree, msgFile) ;
   } else {
      ETree_writeForHumanEye(frontETree, msgFile) ;
   }
   fflush(msgFile) ;
}
MARKTIME(t1) ;
ETree_permuteVertices(frontETree, oldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute ETree", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after permuting the vertex map") ;
   if ( msglvl == 2 ) {
      ETree_writeStats(frontETree, msgFile) ;
   } else {
      ETree_writeForHumanEye(frontETree, msgFile) ;
   }
   fflush(msgFile) ;
}
return(oldToNewIV) ; }

/*--------------------------------------------------------------------*/
