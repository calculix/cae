/*  mkNDlinsysQR.c  */

#include "../../FrontMtx.h"
#include "../../EGraph.h"
#include "../../SymbFac.h"
#include "../../Drand.h"
#include "../../misc.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
  ---------------------------------------------------------------------
   purpose -- to create an overdetermined linear system A*X = B
      for a nested dissection ordering of a 2-d or 3-d regular grid

   input --
 
      n1 -- # of nodes in first direction
      n2 -- # of nodes in second direction
      n3 -- # of nodes in third direction
      type -- type of entries
         SPOOLES_REAL or SPOOLES_COMPLEX
      nrhs    -- number of right hand sides
      seed    -- seed for random number generator
      msglvl  -- message level
      msgFile -- message file
 
   output --
 
      pfrontETree -- to be filled with address of front tree 
      psymbfacIVL -- to be filled with address of symbolic factorization
      pmtxA       -- to be filled with address of matrix object A
      pmtxX       -- to be filled with address of matrix object X
      pmtxB       -- to be filled with address of matrix object B
 
   created -- 98may29, cca
   ---------------------------------------------------------------------
*/
void
mkNDlinsysQR (
   int        n1,
   int        n2,
   int        n3,
   int        type,
   int        nrhs,
   int        seed,
   int        msglvl,
   FILE       *msgFile,
   ETree      **pfrontETree,
   IVL        **psymbfacIVL,
   InpMtx     **pmtxA,
   DenseMtx   **pmtxX,
   DenseMtx   **pmtxB
) {
DenseMtx        *mtxB, *mtxX ;
double          one[2] = {1.0, 0.0} ;
Drand           drand ;
double          t1, t2 ;
double          *entries ;
EGraph          *egraph ;
ETree           *etree, *frontETree ;
Graph           *graph ;
InpMtx          *mtxA ;
int             ielem, ii, irow, jrow, mrow, ndim, nelem, nentA, 
                neqns, nfront, nrowA, size ;
int             *indices, *newToOld, *oldToNew ;
IV              *nzerosIV, *oldToNewIV ;
IVL             *adjIVL, *symbfacIVL ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setNormal(&drand, 0.0, 1.0) ;
/*
   -----------------------------------------------------
   get the grid adjacency structure and set up the graph
   -----------------------------------------------------
*/
neqns = n1*n2*n3 ;
MARKTIME(t1) ;
if (  (n1 == 1 && n2 == 1)
   || (n1 == 1 && n3 == 1)
   || (n1 == 2 && n3 == 1) ) {
   fprintf(stderr, "\n fatal error in mkNDlinsysQR"
           "\n n1 = %d, n2 = %d, n3 = %d\n", n1, n2, n3) ;
   exit(-1) ;
}
if ( n1 == 1 ) {
   adjIVL = IVL_make9P(n2, n3, 1) ;
   ndim = 2 ;
} else if ( n2 == 1 ) {
   adjIVL = IVL_make9P(n1, n3, 1) ;
   ndim = 2 ;
} else if ( n3 == 1 ) {
   adjIVL = IVL_make9P(n1, n2, 1) ;
   ndim = 2 ;
} else {
   adjIVL = IVL_make27P(n1, n2, n3, 1) ;
   ndim = 3 ;
}
graph = Graph_new() ;
Graph_init2(graph, 0, neqns, 0, adjIVL->tsize, neqns,
            adjIVL->tsize, adjIVL, NULL, NULL) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : create the grid graph",
        t2 - t1) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n grid graph") ;
   Graph_writeForHumanEye(graph, msgFile) ;
}
/*
   ---------------------------------------------
   make the nested dissection permutation vector
   ---------------------------------------------
*/
MARKTIME(t1) ;
newToOld = IVinit(neqns, -1) ;
oldToNew = IVinit(neqns, -1) ;
mkNDperm(n1, n2, n3, newToOld, 0, n1-1, 0, n2-1, 0, n3-1) ;
for ( ii = 0 ; ii < neqns ; ii++ ) {
   oldToNew[newToOld[ii]] = ii ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : make the nested dissection ordering",
        t2 - t1) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n oldToNew") ;
   IVfprintf(msgFile, neqns, oldToNew) ;
}
/*
   ----------------------------------
   create the elimination tree object
   ----------------------------------
*/
MARKTIME(t1) ;
etree = ETree_new() ;
ETree_initFromGraphWithPerms(etree, graph, newToOld, oldToNew) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n elimination tree") ;
   ETree_writeForHumanEye(etree, msgFile) ;
}
nzerosIV = IV_new() ;
IV_init(nzerosIV, neqns, NULL) ;
IV_fill(nzerosIV, 0) ;
frontETree = ETree_mergeFrontsOne(etree, 0, nzerosIV) ;
IV_free(nzerosIV) ;
IVfree(newToOld) ;
IVfree(oldToNew) ;
ETree_free(etree) ;
nfront = frontETree->nfront ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n front tree") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : create the front tree",
        t2 - t1) ;
/*
   --------------------------------------------
   create the symbolic factorization IVL object
   --------------------------------------------
*/
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromGraph(frontETree, graph) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : compute the symbolic factorization",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
        "\n\n symbolic factorization IVL object in original ordering") ;
   if ( msglvl == 2 ) {
      IVL_writeStats(symbfacIVL, msgFile) ;
   } else {
      IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   --------------------------------------
   permute the vertices in the front tree
   --------------------------------------
*/
MARKTIME(t1) ;
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
oldToNew   = IV_entries(oldToNewIV) ;
ETree_permuteVertices(frontETree, oldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute the front tree",
        t2 - t1) ;
/*
   ------------------------------------------------------
   convert the symbolic factorization to the new ordering
   ------------------------------------------------------
*/
IVL_overwrite(symbfacIVL, oldToNewIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
        "\n\n overwritten symbolic factorization IVL object ") ;
   if ( msglvl == 2 ) {
      IVL_writeStats(symbfacIVL, msgFile) ;
   } else {
      IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   }
   fflush(msgFile) ;
}
IVL_sortUp(symbfacIVL) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
        "\n\n sorted symbolic factorization IVL object ") ;
   if ( msglvl == 2 ) {
      IVL_writeStats(symbfacIVL, msgFile) ;
   } else {
      IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   --------------------------------
   create the natural factor matrix
   --------------------------------
*/
MARKTIME(t1) ;
if ( ndim == 2 ) {
   if ( n1 == 1 ) {
      egraph  = EGraph_make9P(n2, n3, 1) ;
      nelem   = (n2-1)*(n3-1) ;
   } else if ( n2 == 1 ) {
      egraph  = EGraph_make9P(n1, n3, 1) ;
      nelem   = (n1-1)*(n3-1) ;
   } else if ( n3 == 1 ) {
      egraph  = EGraph_make9P(n1, n2, 1) ;
      nelem   = (n1-1)*(n2-1) ;
   }
   mrow = 4 ;
} else {
   egraph = EGraph_make27P(n1, n2, n3, 1) ;
   mrow   = 8 ;
   nelem  = (n1-1)*(n2-1)*(n3-1) ;
}
nrowA = mrow*nelem ;
if ( type == SPOOLES_REAL ) {
   entries = DVinit(mrow, 0.0) ;
} else if ( type == SPOOLES_COMPLEX ) {
   entries = DVinit(2*mrow, 0.0) ;
}
nentA = mrow*nrowA ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : create egraph ", t2 - t1) ;
if ( msglvl > 2 ) {
   EGraph_writeForHumanEye(egraph, msgFile) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
mtxA = InpMtx_new() ;
InpMtx_init(mtxA, INPMTX_BY_ROWS, type, nentA, nrowA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the InpMtx object",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n InpMtx after initialization") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
}
for ( ielem = 0, jrow = 0 ; ielem < nelem ; ielem++ ) {
   IVL_listAndSize(egraph->adjIVL, ielem, &size, &indices) ;
   for ( irow = 0 ; irow < mrow ; irow++, jrow++ ) {
      if ( type == SPOOLES_REAL ) {
         Drand_fillDvector(&drand, size, entries) ;
         InpMtx_inputRealRow(mtxA, jrow, size, indices, entries) ;
      } else if ( type == SPOOLES_COMPLEX ) {
         Drand_fillDvector(&drand, 2*size, entries) ;
         InpMtx_inputComplexRow(mtxA, jrow, size, indices, entries) ;
      }
   }
}
DVfree(entries) ;
EGraph_free(egraph) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n natural factor matrix") ;
   if ( msglvl == 2 ) {
      InpMtx_writeStats(mtxA, msgFile) ;
   } else {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   -------------------------------------------------------------
   permute the InpMtx object into the nested dissection ordering
   -------------------------------------------------------------
*/
MARKTIME(t1) ;
InpMtx_permute(mtxA,  NULL, IV_entries(oldToNewIV)) ;
InpMtx_changeCoordType(mtxA, INPMTX_BY_ROWS) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute inpmtx ", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after permute InpMtx object ") ;
   if ( msglvl == 2 ) {
      InpMtx_writeStats(mtxA, msgFile) ;
   } else {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   }
}
if ( msglvl > 5 ) {
   fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowA, neqns) ;
   InpMtx_writeForMatlab(mtxA, "A", msgFile) ;
}
/*
   --------------------------------------------------------
   generate the linear system
   1. generate solution matrix and fill with random numbers
   2. generate rhs matrix and fill with zeros
   3. compute matrix-matrix multiply
   --------------------------------------------------------
*/
mtxX = DenseMtx_new() ;
DenseMtx_init(mtxX, type, 0, -1, neqns, nrhs, 1, neqns) ;
DenseMtx_fillRandomEntries(mtxX, &drand) ;
mtxB = DenseMtx_new() ;
DenseMtx_init(mtxB, type, 1, -1, nrowA, nrhs, 1, nrowA) ;
DenseMtx_zero(mtxB) ;
InpMtx_nonsym_mmm(mtxA, mtxB, one, mtxX) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n\n mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 5 ) {
   fprintf(msgFile, "\n X = zeros(%d,%d) ;", mtxX->nrow, mtxX->ncol) ;
   DenseMtx_writeForMatlab(mtxX, "X", msgFile) ;
   fprintf(msgFile, "\n B = zeros(%d,%d) ;", mtxB->nrow, mtxB->ncol) ;
   DenseMtx_writeForMatlab(mtxB, "B", msgFile) ;
}
/*
   -----------------------
   set the output pointers
   -----------------------
*/
*pfrontETree = frontETree ;
*psymbfacIVL = symbfacIVL ;
*pmtxA       = mtxA       ;
*pmtxX       = mtxX       ;
*pmtxB       = mtxB       ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
IV_free(oldToNewIV) ;

return ; }

/*--------------------------------------------------------------------*/
