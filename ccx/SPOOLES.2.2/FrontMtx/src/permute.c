/*  permute.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
static void FrontMtx_reorderRowIndices ( FrontMtx *frontmtx, int J, 
   int K, int map[], int msglvl, FILE *msgFile ) ;
static void FrontMtx_reorderColumnIndices ( FrontMtx *frontmtx, int J, 
   int K, int map[], int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- permute the upper adjacency structure so that the
      indices in bnd{J} are in ascending order w.r.t. their ancestors

   created -- 98mar05, cca
   ------------------------------------------------------------------
*/
void
FrontMtx_permuteUpperAdj (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) {
int    J, K, neqns ;
int    *map, *par ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_permuteUpperAdj(%p,%d,%p)"
           "\n badn input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
neqns = FrontMtx_neqns(frontmtx) ;
map   = IVinit(neqns, -1) ;
tree  = FrontMtx_frontTree(frontmtx) ;
par   = tree->par ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   if ( (K = par[J]) != -1 ) {
      FrontMtx_reorderColumnIndices(frontmtx, J, K, map,
                                     msglvl, msgFile) ;
   }
}
IVfree(map) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- permute the lower adjacency structure so that the
      indices in bnd{J} are in ascending order w.r.t. their ancestors

   created -- 98mar05, cca
   ------------------------------------------------------------------
*/
void
FrontMtx_permuteLowerAdj (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) {
int    J, K, neqns ;
int    *map, *par ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_permuteLowerAdj(%p,%d,%p)"
           "\n badn input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
neqns = FrontMtx_neqns(frontmtx) ;
map   = IVinit(neqns, -1) ;
tree  = FrontMtx_frontTree(frontmtx) ;
par   = tree->par ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   if ( (K = par[J]) != -1 ) {
      FrontMtx_reorderRowIndices(frontmtx, J, K, map, 
                                  msglvl, msgFile) ;
   }
}
IVfree(map) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- if the columns indices of the front matrix are in
      different order than the column indices of U_{J,bnd{J}}, 
      sort the columns of U_{J,bnd{J}} into ascending order 
      w.r.t the column indices of the front matrix.

   created -- 98mar05, cca
   ------------------------------------------------------------
*/
void
FrontMtx_permuteUpperMatrices (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) {
SubMtx   *mtxUJ ;
int      ii, jj, J, mustdo, neqns, nfront, ncolJ, ncolUJ, nJ ;
int      *map, *colindJ, *colindUJ ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_permuteUpperMatrices(%p,%d,%p)"
           "\n badn input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
neqns  = FrontMtx_neqns(frontmtx) ;
map    = IVinit(neqns, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   mtxUJ = FrontMtx_upperMtx(frontmtx, J, nfront) ;
   if ( mtxUJ != NULL ) {
      nJ = FrontMtx_frontSize(frontmtx, J) ;
      FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
      SubMtx_columnIndices(mtxUJ, &ncolUJ, &colindUJ) ;
      for ( ii = nJ, jj = mustdo = 0 ; ii < ncolJ ; ii++, jj++ ) {
         if ( colindJ[ii] != colindUJ[jj] ) {
            mustdo = 1 ; break ;
         }
      }
      if ( mustdo == 1 ) {
         for ( ii = 0 ; ii < ncolJ ; ii++ ) {
            map[colindJ[ii]] = ii ;
         }
         for ( ii = 0 ; ii < ncolUJ ; ii++ ) {
            colindUJ[ii] = map[colindUJ[ii]] ;
         }
         SubMtx_sortColumnsUp(mtxUJ) ;
         for ( ii = 0 ; ii < ncolUJ ; ii++ ) {
            colindUJ[ii] = colindJ[colindUJ[ii]] ;
         }
      }
   }
}
IVfree(map) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- if the row indices of the front matrix are in
      different order than the row indices of L_{bnd{J},J}, 
      sort the rows of L_{bnd{J},J} into ascending order 
      w.r.t the row indices of the front matrix.

   created -- 98mar05, cca
   --------------------------------------------------------
*/
void
FrontMtx_permuteLowerMatrices (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) {
SubMtx   *mtxLJ ;
int      ii, jj, J, mustdo, neqns, nfront, nJ, nrowJ, nrowUJ ;
int      *map, *rowindJ, *rowindUJ ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_permuteLowerMatrices(%p,%d,%p)"
           "\n badn input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
neqns  = FrontMtx_neqns(frontmtx) ;
map    = IVinit(neqns, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   mtxLJ = FrontMtx_lowerMtx(frontmtx, nfront, J) ;
   if ( mtxLJ != NULL ) {
      nJ = FrontMtx_frontSize(frontmtx, J) ;
      FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
      SubMtx_rowIndices(mtxLJ, &nrowUJ, &rowindUJ) ;
      for ( ii = nJ, jj = mustdo = 0 ; ii < nrowJ ; ii++, jj++ ) {
         if ( rowindJ[ii] != rowindUJ[jj] ) {
            mustdo = 1 ; break ;
         }
      }
      if ( mustdo == 1 ) {
         for ( ii = 0 ; ii < nrowJ ; ii++ ) {
            map[rowindJ[ii]] = ii ;
         }
         for ( ii = 0 ; ii < nrowUJ ; ii++ ) {
            rowindUJ[ii] = map[rowindUJ[ii]] ;
         }
         SubMtx_sortRowsUp(mtxLJ) ;
         for ( ii = 0 ; ii < nrowUJ ; ii++ ) {
            rowindUJ[ii] = rowindJ[rowindUJ[ii]] ;
         }
      }
   }
}
IVfree(map) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to reorder the column indices in bnd{J}
      to be ascending order w.r.t. K cup bnd{K}

   created -- 98mar07, cca
   --------------------------------------------------
*/
static void
FrontMtx_reorderColumnIndices (
   FrontMtx   *frontmtx,
   int         J, 
   int         K,
   int         map[],
   int         msglvl,
   FILE        *msgFile
) {
int   ii, ncolJ, ncolK, nJ ;
int   *colindJ, *colindK ;

if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside reorderColumnIndices(%d,%d)", J, K) ;
   fflush(msgFile) ;
}
nJ = FrontMtx_frontSize(frontmtx, J) ;
FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nJ = %d, ncolJ = %d", nJ, ncolJ) ;
   fflush(msgFile) ;
}
if ( ncolJ == 0 ) {
   return ;
}
FrontMtx_columnIndices(frontmtx, K, &ncolK, &colindK) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, ", ncolK = %d", ncolK) ;
   fflush(msgFile) ;
}
if ( ncolK == 0 ) {
   fprintf(stderr, "\n fatal error FrontMtx_reorderColumnIndices()"
           "\n J = %d, K = %d, nJ = %d, ncolJ = %d, ncolK = %d\n",
           J, K, nJ, ncolJ, ncolK) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n colindJ") ;
   IVfprintf(msgFile, ncolJ, colindJ) ;
   fprintf(msgFile, "\n colindK") ;
   IVfprintf(msgFile, ncolK, colindK) ;
   fflush(msgFile) ;
}
for ( ii = 0 ; ii < ncolK ; ii++ ) {
   map[colindK[ii]] = ii ;
}
for ( ii = nJ ; ii < ncolJ ; ii++ ) {
   colindJ[ii] = map[colindJ[ii]] ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n local colindJ") ;
   IVfprintf(msgFile, ncolJ, colindJ) ;
   fflush(msgFile) ;
}
IVqsortUp(ncolJ - nJ, colindJ + nJ) ;
for ( ii = nJ ; ii < ncolJ ; ii++ ) {
   colindJ[ii] = colindK[colindJ[ii]] ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n global colindJ") ;
   IVfprintf(msgFile, ncolJ, colindJ) ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to reorder the row indices in bnd{J}
      to be ascending order w.r.t. K cup bnd{K}

   created -- 98mar07, cca
   -----------------------------------------------
*/
static void
FrontMtx_reorderRowIndices (
   FrontMtx   *frontmtx,
   int         J, 
   int         K,
   int         map[],
   int         msglvl,
   FILE        *msgFile
) {
int   ii, nrowJ, nrowK, nJ ;
int   *rowindJ, *rowindK ;

if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside reorderRowIndices(%d,%d)", J, K) ;
   fflush(msgFile) ;
}
nJ = FrontMtx_frontSize(frontmtx, J) ;
FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nJ = %d, nrowJ = %d", nJ, nrowJ) ;
   fflush(msgFile) ;
}
if ( nrowJ == 0 ) {
   return ;
}
FrontMtx_rowIndices(frontmtx, K, &nrowK, &rowindK) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, ", nrowK = %d", nrowK) ;
   fflush(msgFile) ;
}
if ( nrowK == 0 ) {
   fprintf(stderr, "\n fatal error FrontMtx_reorderRowIndices()"
           "\n J = %d, K = %d, nJ = %d, nrowJ = %d, nrowK = %d\n",
           J, K, nJ, nrowJ, nrowK) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n rowindJ") ;
   IVfprintf(msgFile, nrowJ, rowindJ) ;
   fprintf(msgFile, "\n rowindK") ;
   IVfprintf(msgFile, nrowK, rowindK) ;
   fflush(msgFile) ;
}
for ( ii = 0 ; ii < nrowK ; ii++ ) {
   map[rowindK[ii]] = ii ;
}
for ( ii = nJ ; ii < nrowJ ; ii++ ) {
   rowindJ[ii] = map[rowindJ[ii]] ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n local rowindJ") ;
   IVfprintf(msgFile, nrowJ, rowindJ) ;
   fflush(msgFile) ;
}
IVqsortUp(nrowJ - nJ, rowindJ + nJ) ;
for ( ii = nJ ; ii < nrowJ ; ii++ ) {
   rowindJ[ii] = rowindK[rowindJ[ii]] ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n global rowindJ") ;
   IVfprintf(msgFile, nrowJ, rowindJ) ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
