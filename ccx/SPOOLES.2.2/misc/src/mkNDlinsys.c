/*  mkNDlinsys.c  */

#include "../../FrontMtx.h"
#include "../../Drand.h"
#include "../../SymbFac.h"
#include "../../timings.h"
#include "../../misc.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to create a linear system A*X = B
      for a nested dissection ordering of a 2-d or 3-d regular grid
   input --

      n1 -- # of nodes in first direction
      n2 -- # of nodes in second direction
      n3 -- # of nodes in third direction
      maxzeros -- relaxation factor for fronts,
         maximum number of zero entries in a front
      maxsize  -- split parameter for large fronts,
         maximum number of internal vertices in a front
      type -- type of entries
         SPOOLES_REAL or SPOOLES_COMPLEX
      symmetryflag -- symmetry of the matrix
         SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
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

   created -- 98may16, cca
   ---------------------------------------------------------------------
*/
void
mkNDlinsys (
   int        n1,
   int        n2,
   int        n3,
   int        maxzeros,
   int        maxsize,
   int        type,
   int        symmetryflag,
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
DenseMtx   *mtxB, *mtxX ;
InpMtx     *mtxA ;
double     one[2] = {1.0, 0.0} ;
double     *dvec ;
Drand      drand ;
double     t1, t2 ;
ETree      *etree, *etree2, *frontETree ;
Graph      *graph ;
int        ient, ii, nent, neqns, nrow, v, vsize, w ;
int        *ivec1, *ivec2, *newToOld, *oldToNew, *rowind, *vadj ;
IV         *nzerosIV, *oldToNewIV ;
IVL        *adjIVL, *symbfacIVL ;

/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setUniform(&drand, -1.0, 1.0) ;
/*
   --------------------------------
   get the grid adjacency structure
   --------------------------------
*/
neqns = n1 * n2 * n3 ;
MARKTIME(t1) ;
if ( n1 == 1 ) {
   adjIVL = IVL_make9P(n2, n3, 1) ;
} else if ( n2 == 1 ) {
   adjIVL = IVL_make9P(n1, n3, 1) ;
} else if ( n3 == 1 ) {
   adjIVL = IVL_make9P(n1, n2, 1) ;
} else {
   adjIVL = IVL_make27P(n1, n2, n3, 1) ;
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
IVfree(newToOld) ;
IVfree(oldToNew) ;
nzerosIV = IV_new() ;
IV_init(nzerosIV, neqns, NULL) ;
IV_fill(nzerosIV, 0) ;
etree2 = ETree_mergeFrontsOne(etree, 0, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n elimination tree") ;
   ETree_writeForHumanEye(etree, msgFile) ;
}
etree2 = ETree_mergeFrontsOne(etree, maxzeros, nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
etree2 = ETree_mergeFrontsAll(etree, maxzeros, nzerosIV) ;
IV_free(nzerosIV) ;
ETree_free(etree) ;
etree = etree2 ;
frontETree = ETree_splitFronts(etree, NULL, maxsize, 0) ;
ETree_free(etree) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n front tree") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : create the front tree",
        t2 - t1) ;
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
   ------------------------
   set up the InpMtx object
   ------------------------
*/
MARKTIME(t1) ;
mtxA = InpMtx_new() ;
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC    :
case SPOOLES_HERMITIAN    :
   nent = (adjIVL->tsize - neqns)/2 + neqns ;
   break ;
case SPOOLES_NONSYMMETRIC :
   nent = adjIVL->tsize ;
   break ;
default :
   fprintf(stderr, "\n fatal error in mkNDlinsys()"
           "\n invalid symmetryflag %d\n", symmetryflag) ;
   exit(-1) ;
   break ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n neqns = %d, nent = %d", neqns, nent) ;
}
InpMtx_init(mtxA, INPMTX_BY_ROWS, type, nent, 0) ;
ivec1 = InpMtx_ivec1(mtxA) ;
ivec2 = InpMtx_ivec2(mtxA) ;
dvec  = InpMtx_dvec(mtxA) ;
if ( type == SPOOLES_REAL ) {
   Drand_fillDvector(&drand, nent, dvec) ;
} else if ( type == SPOOLES_COMPLEX ) {
   Drand_fillDvector(&drand, 2*nent, dvec) ;
}
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC :
   for ( v = 0, ient = 0 ; v < neqns ; v++ ) {
      IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         if ( vadj[ii] >= v ) {
            ivec1[ient] = v ;
            ivec2[ient] = vadj[ii] ;
            ient++ ;
         }
      }
   }
/*
   -----------------------------
   code for a laplacian operator
   -----------------------------
*/
/*
   if ( n1 == 1 || n2 == 1 || n3 == 1 ) {
      for ( ii = 0 ; ii < ient ; ii++ ) {
         if ( ivec1[ii] == ivec2[ii] ) {
            dvec[ii] = 8.0 ;
         } else {
            dvec[ii] = -1.0 ;
         }
      }
   } else {
      for ( ii = 0 ; ii < ient ; ii++ ) {
         if ( ivec1[ii] == ivec2[ii] ) {
            dvec[ii] = 27.0 ;
         } else {
            dvec[ii] = -1.0 ;
         }
      }
   }
*/
   break ;
case SPOOLES_HERMITIAN :
   for ( v = 0, ient = 0 ; v < neqns ; v++ ) {
      IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) == v ) {
            ivec1[ient] = v ;
            ivec2[ient] = w ;
            dvec[2*ient+1] = 0.0 ;
            ient++ ;
         } else if ( w > v ) {
            ivec1[ient] = v ;
            ivec2[ient] = w ;
            ient++ ;
         }
      }
   }
   break ;
case SPOOLES_NONSYMMETRIC :
   for ( v = 0, ient = 0 ; v < neqns ; v++ ) {
      IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++, ient++ ) {
         ivec1[ient] = v ;
         ivec2[ient] = vadj[ii] ;
      }
   }
   break ;
}
InpMtx_setNent(mtxA, nent) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n raw matrix object") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
}
InpMtx_sortAndCompress(mtxA) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n original mtxA") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fflush(msgFile) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set up the InpMtxA object",
        t2 - t1) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% start MATLAB FILE") ;
   InpMtx_writeForMatlab(mtxA, "A", msgFile) ;
   if (  symmetryflag == SPOOLES_SYMMETRIC ) {
      fprintf(msgFile,
              "\n neqns = %d ; "
              "\n for ii = 1:neqns "
              "\n    for j = ii+1:neqns "
              "\n       A(j,ii) = A(ii,j) ;"
              "\n    end"
              "\n end", neqns) ;
   } else if ( symmetryflag == SPOOLES_HERMITIAN ) {
      fprintf(msgFile,
              "\n neqns = %d ; "
              "\n for ii = 1:neqns "
              "\n    for j = i+1:neqns "
              "\n       A(j,ii) = ctranspose(A(ii,j)) ;"
              "\n    end"
              "\n end", neqns) ;
   }
   fflush(msgFile) ;
   fprintf(msgFile, "\n %% end MATLAB FILE") ;
}
/*
   --------------------------------------------------------
   generate the linear system
   1. generate solution matrix and fill with random numbers
   2. generate rhs matrix and fill with zeros
   3. compute matrix-matrix multiply
   --------------------------------------------------------
*/
MARKTIME(t1) ;
mtxX = DenseMtx_new() ;
DenseMtx_init(mtxX, type, 0, -1, neqns, nrhs, 1, neqns) ;
DenseMtx_fillRandomEntries(mtxX, &drand) ;
mtxB = DenseMtx_new() ;
DenseMtx_init(mtxB, type, 1, -1, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(mtxB) ;
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC : 
   InpMtx_sym_mmm(mtxA, mtxB, one, mtxX) ;
   break ;
case SPOOLES_HERMITIAN :
   InpMtx_herm_mmm(mtxA, mtxB, one, mtxX) ;
   break ;
case SPOOLES_NONSYMMETRIC :
   InpMtx_nonsym_mmm(mtxA, mtxB, one, mtxX) ;
   break ;
default :
   break ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set up the solution and rhs ",
        t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n original mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n\n original mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% start MATLAB FILE") ;
   DenseMtx_writeForMatlab(mtxX, "X", msgFile) ;
   DenseMtx_writeForMatlab(mtxB, "B", msgFile) ;
   fprintf(msgFile, "\n %% end MATLAB FILE") ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------
   permute the matrix into the nested dissection ordering
   ------------------------------------------------------
*/
MARKTIME(t1) ;
InpMtx_permute(mtxA, oldToNew, oldToNew) ;
/*
   ------------------------------------------------
   map entries into the upper triangle if necessary
   ------------------------------------------------
*/
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC : 
   InpMtx_mapToUpperTriangle(mtxA) ;
   break ;
case SPOOLES_HERMITIAN :
   InpMtx_mapToUpperTriangleH(mtxA) ;
   break ;
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   break ;
}
InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute the matrix", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted mtxA") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% start MATLAB FILE") ;
   InpMtx_writeForMatlab(mtxA, "Anew", msgFile) ;
   if (  symmetryflag == SPOOLES_SYMMETRIC ) {
      fprintf(msgFile,
              "\n neqns = %d ; "
              "\n for ii = 1:neqns "
              "\n    for j = ii+1:neqns "
              "\n       Anew(j,ii) = Anew(ii,j) ;"
              "\n    end"
              "\n end", neqns) ;
   } else if ( symmetryflag == SPOOLES_HERMITIAN ) {
      fprintf(msgFile,
              "\n neqns = %d ; "
              "\n for ii = 1:neqns "
              "\n    for j = ii+1:neqns "
              "\n       Anew(j,ii) = ctranspose(Anew(ii,j)) ;"
              "\n    end"
              "\n end", neqns) ;
   }
   fprintf(msgFile, "\n %% end MATLAB FILE") ;
}
/*
   ----------------------------------------
   permute the solution and right hand side
   ----------------------------------------
*/
MARKTIME(t1) ;
DenseMtx_rowIndices(mtxX, &nrow, &rowind) ;
IVcopy(nrow, rowind, oldToNew) ;
DenseMtx_sort(mtxX) ;
DenseMtx_rowIndices(mtxB, &nrow, &rowind) ;
IVcopy(nrow, rowind, oldToNew) ;
DenseMtx_sort(mtxB) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute the solution and rhs",
        t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n\n permuted mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% start MATLAB FILE") ;
   DenseMtx_writeForMatlab(mtxX, "Xnew", msgFile) ;
   DenseMtx_writeForMatlab(mtxB, "Bnew", msgFile) ;
   fprintf(msgFile, "\n %% end MATLAB FILE") ;
}
/*
   --------------------------------------------
   create the symbolic factorization IVL object
   --------------------------------------------
*/
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : compute the symbolic factorization",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n symbolic factorization IVL object") ;
   if ( msglvl == 2 ) {
      IVL_writeStats(symbfacIVL, msgFile) ;
   } else {
      IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   --------------------------------------
   convert the matrix storage to chevrons
   --------------------------------------
*/
MARKTIME(t1) ;
InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : convert to chevron vectors ",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n InpMtx object ") ;
   if ( msglvl == 2 ) {
      InpMtx_writeStats(mtxA, msgFile) ;
   } else if ( msglvl > 3 ) {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   }
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
