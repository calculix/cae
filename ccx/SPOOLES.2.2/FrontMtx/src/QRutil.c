/*  QRutil.c  */

#include "../FrontMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- to setup two data structures for a QR serial
              or multithreaded factorization

   rowsIVL[J]    -- list of rows of A to be assembled into front J
   firstnz[irow] -- column with location of leading nonzero of row in A

   created -- 98may29, cca
   --------------------------------------------------------------------
*/
void
FrontMtx_QR_setup (
   FrontMtx   *frontmtx,
   InpMtx     *mtxA,
   IVL        **prowsIVL,
   int        **pfirstnz,
   int        msglvl,
   FILE       *msgFile
) {
int   count, irow, jcol, J, loc, neqns, nfront, nrowA, rowsize ;
int   *firstnz, *head, *link, *list, *rowind, *vtxToFront ;
IVL   *rowsIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || mtxA == NULL || prowsIVL == NULL
   || pfirstnz == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_QR_setup()"
           "\n bad input\n") ;
   exit(-1) ; 
}
neqns      = FrontMtx_neqns(frontmtx) ;
nfront     = FrontMtx_nfront(frontmtx) ;
vtxToFront = ETree_vtxToFront(frontmtx->frontETree) ;
/*
   ----------------------------------------------------------------
   create the rowsIVL object,
   list(J) = list of rows that are assembled in front J
   firstnz[irowA] = first column with nonzero element in A(irowA,*)
   ----------------------------------------------------------------
*/
InpMtx_changeCoordType(mtxA, INPMTX_BY_ROWS) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
nrowA = 1 + IVmax(InpMtx_nent(mtxA), InpMtx_ivec1(mtxA), &loc) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n nrowA = %d ", nrowA) ;
   fflush(msgFile) ;
}
firstnz = IVinit(nrowA, -1) ;
head    = IVinit(nfront, -1) ;
link    = IVinit(nrowA, -1) ;
for ( irow = nrowA - 1 ; irow >= 0 ; irow-- ) {
   InpMtx_vector(mtxA, irow, &rowsize, &rowind) ;
   if ( rowsize > 0 ) {
      firstnz[irow] = jcol = rowind[0] ;
      J             = vtxToFront[jcol] ;
      link[irow]    = head[J]          ;
      head[J]       = irow             ;
   }
}
rowsIVL = IVL_new() ;
IVL_init2(rowsIVL, IVL_CHUNKED, nfront, nrowA) ;
list = IVinit(neqns, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   count = 0 ;
   for ( irow = head[J] ; irow != -1 ; irow = link[irow] ) {
      list[count++] = irow ;
   }
   if ( count > 0 ) {
      IVL_setList(rowsIVL, J, count, list) ;
   }
}
IVfree(head) ;
IVfree(link) ;
IVfree(list) ;
/*
   ---------------------------
   set the pointers for return
   ---------------------------
*/
*prowsIVL = rowsIVL ;
*pfirstnz = firstnz ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose --  visit front J during a serial 
               or multithreaded QR factorization

   cpus[1] -- initialize and load staircase matrix
   cpus[2] -- factor the matrix
   cpus[3] -- scale and store factor entries
   cpus[4] -- store update entries

   created -- 98may28, cca
   -----------------------------------------------
*/
void
FrontMtx_QR_factorVisit (
   FrontMtx     *frontmtx,
   int          J,
   InpMtx       *mtxA,
   IVL          *rowsIVL,
   int          firstnz[],
   ChvList      *updlist,
   ChvManager   *chvmanager,
   char         status[],
   int          colmap[],
   DV           *workDV,
   double       cpus[],
   double       *pfacops,
   int          msglvl,
   FILE         *msgFile
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || mtxA == NULL || rowsIVL == NULL
     || firstnz == NULL || updlist == NULL || chvmanager == NULL 
     || status == NULL || colmap == NULL || workDV == NULL 
     || cpus == NULL || pfacops == NULL 
     || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(msgFile, "\n fatal error in FrontMtx_QR_factorVisit(%d)"
           "\n bad input\n", J) ;
   exit(-1) ;
}
/*
   ------------------------------------------------------------
   check to see if all incoming updates are present in the list
   ------------------------------------------------------------
*/
if ( ChvList_isCountZero(updlist, J) == 1 ) {
   A2       *frontJ ;
   Chv      *firstchild, *updchv ;
   double   ops, t1, t2 ;
   int      K ;
/*
   ----------------------------------------
   everything is ready to factor this front
   ----------------------------------------
*/
   firstchild = ChvList_getList(updlist, J) ;
/*
   ----------------------------------------
   initialize and load the staircase matrix
   ----------------------------------------
*/
   MARKTIME(t1) ;
   frontJ = FrontMtx_QR_assembleFront(frontmtx, J, mtxA, rowsIVL,
                                      firstnz, colmap, firstchild, 
                                      workDV, msglvl, msgFile) ;
   if ( firstchild != NULL ) {
      ChvManager_releaseListOfObjects(chvmanager, firstchild) ;
   }
   MARKTIME(t2) ;
   cpus[1] += t2 - t1 ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n after assembling front") ;
      A2_writeForHumanEye(frontJ, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ---------------------------
   factor the staircase matrix
   ---------------------------
*/
   MARKTIME(t1) ;
   ops = A2_QRreduce(frontJ, workDV, msglvl, msgFile) ;
   *pfacops += ops ;
   MARKTIME(t2) ;
   cpus[2] += t2 - t1 ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n after factoring") ;
      A2_writeForHumanEye(frontJ, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ----------------------------------
   scale and store the factor entries
   ----------------------------------
*/
   MARKTIME(t1) ;
   FrontMtx_QR_storeFront(frontmtx, J, frontJ, msglvl, msgFile) ;
   MARKTIME(t2) ;
   cpus[3] += t2 - t1 ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n after storing factor entries") ;
      A2_writeForHumanEye(frontJ, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------------
   store the update matrix
   -----------------------
*/
   if ( (K = frontmtx->tree->par[J]) != -1 ) {
      MARKTIME(t1) ;
      updchv = FrontMtx_QR_storeUpdate(frontmtx, J, frontJ, chvmanager,
                                       msglvl, msgFile) ;
      ChvList_addObjectToList(updlist, updchv, K) ;
      MARKTIME(t2) ;
      cpus[4] += t2 - t1 ;
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n after storing update entries") ;
         A2_writeForHumanEye(frontJ, msgFile) ;
         fflush(msgFile) ;
      }
   }
/*
   -------------------------
   free the staircase matrix
   -------------------------
*/
   A2_free(frontJ) ;
/*
   -------------------------------------
   set the status as finished and return
   -------------------------------------
*/
   status[J] = 'F' ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- create and return an A2 object that contains rows
              of A and rows from update matrices of the children.
              the matrix may not be in staircase form

   created -- 98may25, cca
   --------------------------------------------------------------
*/
A2 *
FrontMtx_QR_assembleFront (
   FrontMtx   *frontmtx,
   int        J,
   InpMtx     *mtxA,
   IVL        *rowsIVL,
   int        firstnz[],
   int        colmap[],
   Chv        *firstchild,
   DV         *workDV,
   int        msglvl,
   FILE       *msgFile
) {
A2       *frontJ ;
Chv      *chvI ;
double   *rowI, *rowJ, *rowentA ;
int      ii, irow, irowA, irowI, jcol, jj, jrow, ncolI, ncolJ, 
         nentA, nrowI, nrowJ, nrowFromA ;
int      *colindA, *colindI, *colindJ, *map, *rowids, *rowsFromA ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || mtxA == NULL || rowsIVL == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_QR_assembleFront()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n inside FrontMtx_QR_assembleFront(%d)", J) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------
   set up the map from global to local column indices
   --------------------------------------------------
*/
FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
for ( jcol = 0 ; jcol < ncolJ ; jcol++ ) {
   colmap[colindJ[jcol]] = jcol ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n front %d's column indices", J) ;
   IVfprintf(msgFile, ncolJ, colindJ) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   compute the size of the front and map the global 
   indices of the update matrices into local indices
   -------------------------------------------------
*/
IVL_listAndSize(rowsIVL, J, &nrowFromA, &rowsFromA) ;
nrowJ = nrowFromA ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %d rows from A", nrowFromA) ;
   fflush(msgFile) ;
}
for ( chvI = firstchild ; chvI != NULL ; chvI = chvI->next ) {
   nrowJ += chvI->nD ;
   Chv_columnIndices(chvI, &ncolI, &colindI) ;
   for ( jcol = 0 ; jcol < ncolI ; jcol++ ) {
      colindI[jcol] = colmap[colindI[jcol]] ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n %d rows from child %d", chvI->nD, chvI->id) ;
      fflush(msgFile) ;
   }
}
/*
   ----------------------------------------------------------
   get workspace for the rowids[nrowJ] and map[nrowJ] vectors
   ----------------------------------------------------------
*/
if ( sizeof(int) == sizeof(double) ) {
   DV_setSize(workDV, 2*nrowJ) ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   DV_setSize(workDV, nrowJ) ;
}
rowids = (int *) DV_entries(workDV) ;
map    = rowids + nrowJ ;
IVramp(nrowJ, rowids, 0, 1) ;
IVfill(nrowJ, map, -1) ;
/*
   -----------------------------------------------------------------
   get the map from incoming rows to their place in the front matrix
   -----------------------------------------------------------------
*/
for ( irow = jrow = 0 ; irow < nrowFromA ; irow++, jrow++ ) {
   irowA = rowsFromA[irow] ;
   map[jrow] = colmap[firstnz[irowA]] ;
}
for ( chvI = firstchild ; chvI != NULL ; chvI = chvI->next ) {
   nrowI = chvI->nD ;
   Chv_columnIndices(chvI, &ncolI, &colindI) ;
   for ( irow = 0 ; irow < nrowI ; irow++, jrow++ ) {
      map[jrow] = colindI[irow] ;
   }
}
IV2qsortUp(nrowJ, map, rowids) ;
for ( irow = 0 ; irow < nrowJ ; irow++ ) {
   map[rowids[irow]] = irow ;
}
/*
   ----------------------------
   allocate the A2 front object
   ----------------------------
*/
frontJ = A2_new() ;
A2_init(frontJ, frontmtx->type, nrowJ, ncolJ, ncolJ, 1, NULL) ;
A2_zero(frontJ) ;
/*
   ------------------------------------
   load the original rows of the matrix
   ------------------------------------
*/
for ( jrow = 0 ; jrow < nrowFromA ; jrow++ ) {
   irowA = rowsFromA[jrow] ;
   rowJ  = A2_row(frontJ, map[jrow]) ;
   if ( A2_IS_REAL(frontJ) ) {
      InpMtx_realVector(mtxA, irowA, &nentA, &colindA, &rowentA) ;
   } else if ( A2_IS_COMPLEX(frontJ) ) {
      InpMtx_complexVector(mtxA, irowA, &nentA, &colindA, &rowentA) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n loading row %d", irowA) ;
      fprintf(msgFile, "\n global indices") ;
      IVfprintf(msgFile, nentA, colindA) ;
      fflush(msgFile) ;
   }
   if ( A2_IS_REAL(frontJ) ) {
      for ( ii = 0 ; ii < nentA ; ii++ ) {
         jj = colmap[colindA[ii]] ;
         rowJ[jj] = rowentA[ii] ;
      }
   } else if ( A2_IS_COMPLEX(frontJ) ) {
      for ( ii = 0 ; ii < nentA ; ii++ ) {
         jj = colmap[colindA[ii]] ;
         rowJ[2*jj]   = rowentA[2*ii]   ;
         rowJ[2*jj+1] = rowentA[2*ii+1] ;
      }
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n after assembling rows of A") ;
   A2_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   load the updates from the children 
   ----------------------------------
*/
for ( chvI = firstchild ; chvI != NULL ; chvI = chvI->next ) {
   nrowI = chvI->nD ;
   Chv_columnIndices(chvI, &ncolI, &colindI) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n loading child %d", chvI->id) ;
      fprintf(msgFile, "\n child's column indices") ;
      IVfprintf(msgFile, ncolI, colindI) ;
      Chv_writeForHumanEye(chvI, msgFile) ;
      fflush(msgFile) ;
   }
   rowI = Chv_entries(chvI) ;
   for ( irowI = 0 ; irowI < nrowI ; irowI++, jrow++ ) {
      rowJ = A2_row(frontJ, map[jrow]) ;
      if ( A2_IS_REAL(frontJ) ) {
         for ( ii = irowI ; ii < ncolI ; ii++ ) {
            jj = colindI[ii] ;
            rowJ[jj] = rowI[ii] ;
         }
         rowI += ncolI - irowI - 1 ;
      } else if ( A2_IS_COMPLEX(frontJ) ) {
         for ( ii = irowI ; ii < ncolI ; ii++ ) {
            jj = colindI[ii] ;
            rowJ[2*jj]   = rowI[2*ii]   ;
            rowJ[2*jj+1] = rowI[2*ii+1] ;
         }
         rowI += 2*(ncolI - irowI - 1) ;
      }
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n after assembling child %d", chvI->id) ;
      A2_writeForHumanEye(frontJ, msgFile) ;
      fflush(msgFile) ;
   }
}
return(frontJ) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   store the factor entries of the reduced front matrix

   created -- 98may25, cca
   ----------------------------------------------------
*/
void
FrontMtx_QR_storeFront (
   FrontMtx   *frontmtx,
   int        J,
   A2         *frontJ,
   int        msglvl,
   FILE       *msgFile
) {
A2       tempA2 ;
double   fac, ifac, imag, real, rfac ;
double   *entDJJ, *entUJJ, *entUJN, *row ;
int      inc1, inc2, irow, jcol, ncol, ncolJ, nD, nentD, nentUJJ, 
         nfront, nrow, nU ;
int      *colind, *colindJ, *firstlocs, *sizes ;
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || frontJ == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_QR_storeFront()"
           "\n bad input\n") ;
   exit(-1) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
nrow   = A2_nrow(frontJ) ;
ncol   = A2_ncol(frontJ) ;
A2_setDefaultFields(&tempA2) ;
nD = FrontMtx_frontSize(frontmtx, J) ;
nU = ncol - nD ;
/*
   --------------------------------------
   scale the rows and square the diagonal
   --------------------------------------
*/
row = A2_entries(frontJ) ;
if ( A2_IS_REAL(frontJ) ) {
   for ( irow = 0 ; irow < nD ; irow++ ) {
      if ( row[irow] != 0.0 ) {
         fac = 1./row[irow] ;
         for ( jcol = irow + 1 ; jcol < ncol ; jcol++ ) {
            row[jcol] *= fac ;
         }
         row[irow] = row[irow] * row[irow] ;
      }
      row += ncol ;
   }
} else if ( A2_IS_COMPLEX(frontJ) ) {
   for ( irow = 0 ; irow < nD ; irow++ ) {
      real = row[2*irow] ; imag = row[2*irow+1] ;
      if (  real != 0.0 || imag != 0.0 ) {
         Zrecip(real, imag, &rfac, &ifac) ;
         ZVscale(ncol - irow - 1, & row[2*irow+2], rfac, ifac) ;
         row[2*irow]   = real*real + imag*imag ;
         row[2*irow+1] = 0.0 ;
      }
      row += 2*ncol ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n after scaling rows of A") ;
   A2_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   copy the diagonal entries
   -------------------------
*/
mtx = FrontMtx_diagMtx(frontmtx, J) ;
SubMtx_diagonalInfo(mtx, &nentD, &entDJJ) ;
A2_subA2(&tempA2, frontJ, 0, nD-1, 0, nD-1) ;
A2_copyEntriesToVector(&tempA2, nentD, entDJJ, 
                       A2_DIAGONAL, A2_BY_ROWS) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
IVcopy(nD, colind, colindJ) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n diagonal factor matrix") ;
   SubMtx_writeForHumanEye(mtx, msgFile) ;
   fflush(msgFile) ;
}
if ( (mtx = FrontMtx_upperMtx(frontmtx, J, J)) != NULL ) {
/*
   ------------------------
   copy the U_{J,J} entries
   ------------------------
*/
   SubMtx_denseSubcolumnsInfo(mtx, &nD, &nentUJJ, 
                           &firstlocs, &sizes, &entUJJ) ;
   A2_copyEntriesToVector(&tempA2, nentUJJ, entUJJ, 
                          A2_STRICT_UPPER, A2_BY_COLUMNS) ;
   SubMtx_columnIndices(mtx, &ncol, &colind) ;
   IVcopy(nD, colind, colindJ) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n UJJ factor matrix") ;
      SubMtx_writeForHumanEye(mtx, msgFile) ;
      fflush(msgFile) ;
   }
}
if ( ncolJ > nD ) {
/*
   -----------------------------
   copy the U_{J,bnd{J}} entries
   -----------------------------
*/
   mtx = FrontMtx_upperMtx(frontmtx, J, nfront) ;
   SubMtx_denseInfo(mtx, &nD, &nU, &inc1, &inc2, &entUJN) ;
   A2_subA2(&tempA2, frontJ, 0, nD-1, nD, ncolJ-1) ;
   A2_copyEntriesToVector(&tempA2, nD*nU, entUJN, 
                          A2_ALL_ENTRIES, A2_BY_COLUMNS) ;
   SubMtx_columnIndices(mtx, &ncol, &colind) ;
   IVcopy(nU, colind, colindJ + nD) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n UJN factor matrix") ;
      SubMtx_writeForHumanEye(mtx, msgFile) ;
      fflush(msgFile) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to create and return a Chv object that
              holds the update matrix for front J

   created -- 98may25, cca
   -------------------------------------------------
*/
Chv *
FrontMtx_QR_storeUpdate (
   FrontMtx     *frontmtx,
   int          J,
   A2           *frontJ,
   ChvManager   *chvmanager,
   int          msglvl,
   FILE         *msgFile
) {
A2       tempJ ;
Chv      *chvJ ;
double   *updent ;
int      nbytes, ncolJ, ncolupd, nD, nent, nrowJ, nrowupd ;
int      *colindJ, *updind ;
/*
   -----------------------------------------------
   compute the number of rows in the update matrix
   -----------------------------------------------
*/
nD = FrontMtx_frontSize(frontmtx, J) ;
FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
nrowJ = A2_nrow(frontJ) ;
nrowupd = ((nrowJ >= ncolJ) ? ncolJ : nrowJ) - nD ;
ncolupd = ncolJ - nD ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n inside FrontMtx_QR_storeUpdate(%d)", J) ;
   fprintf(msgFile, "\n nD %d, nrowJ %d, nrowupd %d, ncolupd %d",
           nD, nrowJ, nrowupd, ncolupd) ;
   fflush(msgFile) ;
}
if ( nrowupd > 0 && ncolupd > 0 ) {
   if ( FRONTMTX_IS_REAL(frontmtx) ) {
      nbytes = Chv_nbytesNeeded(nrowupd, 0, ncolupd - nrowupd, 
                                SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
   } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
      nbytes = Chv_nbytesNeeded(nrowupd, 0, ncolupd - nrowupd, 
                                SPOOLES_COMPLEX, SPOOLES_HERMITIAN) ;
   }
   chvJ = ChvManager_newObjectOfSizeNbytes(chvmanager, nbytes) ;
   if ( FRONTMTX_IS_REAL(frontmtx) ) {
       Chv_init(chvJ, J, nrowupd, 0, ncolupd - nrowupd, 
                SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
   } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
       Chv_init(chvJ, J, nrowupd, 0, ncolupd - nrowupd, 
                SPOOLES_COMPLEX, SPOOLES_HERMITIAN) ;
   }
   Chv_columnIndices(chvJ, &ncolupd, &updind) ;
   IVcopy(ncolupd, updind, colindJ + nD) ;
   nent   = Chv_nent(chvJ) ;
   updent = Chv_entries(chvJ) ;
   A2_setDefaultFields(&tempJ) ;
   A2_subA2(&tempJ, frontJ, nD, nrowJ - 1, nD, ncolJ - 1) ;
   A2_copyEntriesToVector(&tempJ, nent, updent, A2_UPPER, A2_BY_ROWS) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n update matrix %d", J) ;
      Chv_writeForHumanEye(chvJ, msgFile) ;
      fflush(msgFile) ;
   }
} else {
   chvJ = NULL ;
}
return(chvJ) ; }

/*--------------------------------------------------------------------*/
