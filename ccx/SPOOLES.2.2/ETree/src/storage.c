/*  storage.c  */

#include "../ETree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the active storage to eliminate J
               using the multifrontal method

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC,
              SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 97may21, cca
   ---------------------------------------------------------------
*/
void
ETree_MFstackProfile (
   ETree    *etree,
   int      symflag,
   double   dvec[]
) {
int    I, J, nDJ, nUI, nUJ, stack ;
int    *bndwghts, *fch, *nodwghts, *sib ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_MFstackProfile(%p,%p)"
           "\n bad input\n", etree, dvec) ;
   exit(-1) ;
}
tree     = ETree_tree(etree) ;
nodwghts = ETree_nodwghts(etree) ;
bndwghts = ETree_bndwghts(etree) ;
fch      = ETree_fch(etree) ;
sib      = ETree_sib(etree) ;
/*
   ---------------------------------------------
   loop over the nodes in a post-order traversal
   ---------------------------------------------
*/
stack = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   nDJ = nodwghts[J] ;
   nUJ = bndwghts[J] ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n working on front %d, nD = %d, nU = %d, stack = %d",
           J, nDJ, nUJ, stack) ;
#endif
   if ( (I = fch[J]) != -1 ) {
/*
      --------------------------------
      remove last child from the stack
      --------------------------------
*/
      while ( sib[I] != -1 ) {
         I = sib[I] ;
      }
      nUI = bndwghts[I] ;
      if ( symflag == SPOOLES_SYMMETRIC 
        || symflag == SPOOLES_HERMITIAN ) {
         stack -= (nUI*(nUI+1))/2 ;
      } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
         stack -= nUI*nUI ;
      }
#if MYDEBUG > 0
      fprintf(stdout, 
           "\n    removing last child %d, stack = %d", I, stack) ;
#endif
   }
/*
   ---------------------------------
   add frontal matrix for J to stack
   and store as the storage for J
   ---------------------------------
*/
   dvec[J] = stack + (nDJ + nUJ)*(nDJ + nUJ) ;
   if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
      dvec[J] = stack + ((nDJ + nUJ)*(nDJ + nUJ+1))/2 ;
   } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
      dvec[J] = stack + (nDJ + nUJ)*(nDJ + nUJ) ;
   }
#if MYDEBUG > 0
      fprintf(stdout, 
           "\n    working storage = %.0f", dvec[J]) ;
#endif
   if ( (I = fch[J]) != -1 ) {
/*
      ------------------------------------
      remove other children from the stack
      ------------------------------------
*/
      while ( sib[I] != -1 ) {
         nUI = bndwghts[I] ;
         if ( symflag == SPOOLES_SYMMETRIC 
           || symflag == SPOOLES_HERMITIAN ) {
            stack -= (nUI*(nUI+1))/2 ;
         } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
            stack -= nUI*nUI ;
         }
#if MYDEBUG > 0
      fprintf(stdout, 
           "\n    removing child %d, stack = %d", I, stack) ;
#endif
         I = sib[I] ;
      }
   }
/*
   --------------------------------
   add update matrix for J to stack
   --------------------------------
*/
   if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
      stack += (nUJ*(nUJ+1))/2 ;
   } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
      stack += nUJ*nUJ ;
   }
#if MYDEBUG > 0
      fprintf(stdout, 
           "\n    adding update matrix, stack = %d", stack) ;
#endif
}
#if MYDEBUG >= 0
fprintf(stdout, "\n    MF: final stack = %d", stack) ;
#endif
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the active storage to eliminate J
               using the left-looking general sparse method

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC,
              SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 97may21, cca
   ---------------------------------------------------------------
*/
void
ETree_GSstorageProfile (
   ETree    *etree,
   int      symflag,
   IVL      *symbfacIVL,
   int      *vwghts,
   double   dvec[]
) {
int    count, ii, I, J, K, nDJ, nfront, nUJ, sizeI, sizeJ, storage, v ;
int    *bndwghts, *head, *indI, *indJ, 
       *link, *nodwghts, *offsets, *vtxToFront ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || symbfacIVL == NULL || dvec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ETree_GSstorageProfile(%p,%p,%p,%p)"
           "\n bad input\n", etree, symbfacIVL, vwghts, dvec) ;
   exit(-1) ;
}
tree       = ETree_tree(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
bndwghts   = ETree_bndwghts(etree) ;
vtxToFront = ETree_vtxToFront(etree) ;
nfront     = ETree_nfront(etree) ;
head       = IVinit(nfront, -1) ;
link       = IVinit(nfront, -1) ;
offsets    = IVinit(nfront,  0) ;
/*
   ---------------------------------------------
   loop over the nodes in a post-order traversal
   ---------------------------------------------
*/
storage = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   nDJ = nodwghts[J] ;
   nUJ = bndwghts[J] ;
   if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
      storage += (nDJ*(nDJ + 1))/2 + nDJ*nUJ ;
   } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
      storage += nDJ*nDJ + 2*nDJ*nUJ ;
   }
   dvec[J] = storage ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n working on front %d, nD = %d, nU = %d, storage = %d",
           J, nDJ, nUJ, storage) ;
#endif
/*
   -----------------------------
   loop over the updating fronts
   -----------------------------
*/
   while ( (I = head[J]) != -1 ) {
      head[J] = link[I] ;
      IVL_listAndSize(symbfacIVL, I, &sizeI, &indI) ;
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n    updating front %d, offset = %d, sizeI = %d", 
              I, offsets[I], sizeI) ;
      IVfprintf(stdout, sizeI, indI) ;
#endif
      for ( ii = offsets[I], count = 0, K = -1 ; ii < sizeI ; ii++ ) {
         v = indI[ii] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n       ii = %d, v = %d, K = %d",
                 ii, v, vtxToFront[v]) ;
         fflush(stdout) ;
#endif
         K = vtxToFront[v] ;
         if ( K < 0 || K >= nfront ) {
            fprintf(stderr, "\n\n fatal error"
                    "\n ii = %d, v = %d, K = %d", ii, v, K) ;
            exit(-1) ;
         }
         if ( (K = vtxToFront[v]) != J ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n       linking to next ancestor %d", K) ;
#endif
            link[I] = head[K] ;
            head[K] = I ;
            offsets[I] = ii ;
            break ;
         }
         count += (vwghts == NULL) ? 1 : vwghts[v] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n       count = %d", count) ;
         fflush(stdout) ;
#endif
      }
      if ( symflag == SPOOLES_SYMMETRIC 
        || symflag == SPOOLES_HERMITIAN ) {
         storage -= count*nodwghts[I] ;
      } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
         storage -= 2*count*nodwghts[I] ;
      }
   }
   if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
      storage -= (nDJ*(nDJ+1))/2 ;
   } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
      storage -= nDJ*nDJ ;
   }
   if ( nUJ > 0 ) {
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indJ) ;
      for ( ii = 0 ; ii < sizeJ ; ii++ ) {
         v = indJ[ii] ;
         if ( (K = vtxToFront[v]) != J ) {
            break ;
         }
      }
      offsets[J] = ii ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    linking to next ancestor %d", K) ;
#endif
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indJ) ;
      link[J] = head[K] ;
      head[K] = J ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n    at end of step %d, storage = %d", 
           J, storage) ;
#endif
}
#if MYDEBUG >= 0
fprintf(stdout, "\n    GS: final storage = %d", storage) ;
#endif
IVfree(head) ;
IVfree(link) ;
IVfree(offsets) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the active storage to eliminate J
               using the right-looking general sparse method

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC,
              SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 98dec19, cca
   ---------------------------------------------------------------
*/
void
ETree_FSstorageProfile (
   ETree    *etree,
   int      symflag,
   IVL      *symbfacIVL,
   double   dvec[]
) {
char   *incore ;
int    ii, J, K, nDJ, nfront, nUJ, sizeJ, storage ;
int    *bndwghts, *indJ, *mark, *nodwghts, *stor, *vtxToFront ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || symbfacIVL == NULL || dvec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ETree_FSstorageProfile(%p,%p,%p)"
           "\n bad input\n", etree, symbfacIVL, dvec) ;
   exit(-1) ;
}
tree       = ETree_tree(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
bndwghts   = ETree_bndwghts(etree) ;
vtxToFront = ETree_vtxToFront(etree) ;
nfront     = ETree_nfront(etree) ;
incore     = CVinit(nfront, 'F') ;
stor       = IVinit(nfront, 0) ;
mark       = IVinit(nfront, -1) ;
/*
   --------------------------------------------
   compute the storage for each front's chevron
   --------------------------------------------
*/
if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
   for ( J = 0 ; J < nfront ; J++ ) {
      nDJ = nodwghts[J] ;
      nUJ = bndwghts[J] ;
      stor[J] = (nDJ*(nDJ+1))/2 + nDJ*nUJ ;
   }
} else {
   for ( J = 0 ; J < nfront ; J++ ) {
      nDJ = nodwghts[J] ;
      nUJ = bndwghts[J] ;
      stor[J] = nDJ*nDJ + 2*nDJ*nUJ ;
   }
}
/*
   ---------------------------------------------
   loop over the nodes in a post-order traversal
   ---------------------------------------------
*/
storage = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( incore[J] == 'F' ) {
      storage += stor[J] ;
      incore[J] = 'T' ;
   }
   IVL_listAndSize(symbfacIVL, J, &sizeJ, &indJ) ;
   mark[J] = J ;
   for ( ii = 0 ; ii < sizeJ ; ii++ ) {
      K = vtxToFront[indJ[ii]] ;
      if ( mark[K] != J ) {
         mark[K] = J ;
         if ( incore[K] == 'F' ) {
            storage += stor[K] ;
            incore[K] = 'T' ;
         }
      }
   }
   dvec[J] = storage ;
   storage -= stor[J] ;
}
IVfree(mark) ;
IVfree(stor) ;
CVfree(incore) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the stack storage to solve for J
               in a forward solve

   created -- 97nov30, cca
   ---------------------------------------------------------------
*/
void
ETree_forwSolveProfile (
   ETree    *etree,
   double   dvec[]
) {
int    I, J, maxstack, nDJ, nfront, nUJ, stack ;
int    *bndwghts, *fch, *nodwghts, *sib ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || dvec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ETree_forwSolveProfile(%p,%p)"
           "\n bad input\n", etree, dvec) ;
   exit(-1) ;
}
tree       = ETree_tree(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
bndwghts   = ETree_bndwghts(etree) ;
nfront     = ETree_nfront(etree) ;
fch        = ETree_fch(etree) ;
sib        = ETree_sib(etree) ;
/*
   ---------------------------------------------
   loop over the nodes in a post-order traversal
   ---------------------------------------------
*/
maxstack = stack = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   nDJ = nodwghts[J] ;
   nUJ = bndwghts[J] ;
   stack += nDJ + nUJ ;
   dvec[J] = stack ;
   if ( maxstack < stack ) {
      maxstack = stack ;
   }
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n working on front %d, nD = %d, nU = %d, stack = %d",
           J, nDJ, nUJ, stack) ;
#endif
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      stack -= bndwghts[I] ;
   }
   stack -= nDJ ;
}
#if MYDEBUG >= 0
fprintf(stdout, 
        "\n    forward solve : final stack = %d, max stack = %d", 
        stack, maxstack) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the stack storage to solve for J
               in a backward solve

   created -- 97nov30, cca
   ---------------------------------------------------------------
*/
void
ETree_backSolveProfile (
   ETree    *etree,
   double   dvec[]
) {
int    J, K, maxstack, nDJ, nfront, nUJ, stack ;
int    *bndwghts, *fch, *nodwghts, *par, *sib ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || dvec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ETree_backSolveProfile(%p,%p)"
           "\n bad input\n", etree, dvec) ;
   exit(-1) ;
}
tree       = ETree_tree(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
bndwghts   = ETree_bndwghts(etree) ;
nfront     = ETree_nfront(etree) ;
par        = ETree_par(etree) ;
fch        = ETree_fch(etree) ;
sib        = ETree_sib(etree) ;
/*
   ---------------------------------------------
   loop over the nodes in a post-order traversal
   ---------------------------------------------
*/
maxstack = stack = 0 ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   nDJ = nodwghts[J] ;
   nUJ = bndwghts[J] ;
   stack += nDJ + nUJ ;
   dvec[J] = stack ;
   if ( maxstack < stack ) {
      maxstack = stack ;
   }
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n working on front %d, nD = %d, nU = %d, stack = %d",
           J, nDJ, nUJ, stack) ;
#endif
   if ( (K = par[J]) != -1 && sib[J] == -1 ) {
      stack -= nodwghts[K] + bndwghts[K] ;
   }
   if ( fch[J] == -1 ) {
      stack -= nodwghts[J] + bndwghts[J] ;
   }
}
#if MYDEBUG >= 0
fprintf(stdout, 
        "\n    forward solve : final stack = %d, max stack = %d", 
        stack, maxstack) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
