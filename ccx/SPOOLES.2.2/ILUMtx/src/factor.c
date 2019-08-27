/*  factor.c  */

#include "../ILUMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
static void loadDiagonal ( int type, int neqns, InpMtx *mtxA,
   double entD[], int msglvl, FILE *msgFile ) ;
static int getRowStructure ( int jeqn, int indicesU[], InpMtx *mtxA,
  int sizesU[], int *p_indU[], int headL[], int linkL[], int markU[],
  int msglvl, FILE *msgFile ) ;
static void loadOrigRow ( int jeqn, InpMtx *mtxA, int map[], 
  double temp[], int msglvl, FILE *msgFile ) ;
static double updateRow ( int type, int symmetryflag, int jeqn,
   double LJI[], double DII[], int sizeU, int indU[], double entU[],
   int map[], double temp[], int msglvl, FILE *msgFile ) ;
static int storeRow ( int jeqn, int type, double sigma, int countU, 
   int indicesU[], double temp[], double entD[], int sizesU[], 
   int *p_indU[], double *p_entU[], int msglvl, FILE *msgFile ) ;
static int getColumnStructure ( int jeqn, int indicesL[], InpMtx *mtxA,
   int sizesL[], int *p_indL[], int headU[], int linkU[], int markL[],
   int msglvl, FILE *msgFile ) ;
static void loadOrigColumn ( int jeqn, InpMtx *mtxA, int map[],
   double temp[], int msglvl, FILE *msgFile ) ;
static double updateColumn ( int type, int symmetryflag, int jeqn,
   double DII[], double UIJ[], int sizeL, int indL[], double entL[],
   int map[], double temp[], int msglvl, FILE *msgFile ) ;
static int storeColumn ( int jeqn, int type, double sigma, int countL,
   int indicesL[], double temp[], double entD[], int sizesL[], 
   int *p_indL[], double *p_entL[], int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- to factor the linear system
      A = (L + I)D(I + U), A = (U^T + I)D(I + U) or
      A = (U^H + I)D(I + U). 

   if pops is not NULL, then on return *pops has been incremented
   by the number of operations performed in the factorization.

   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtxLDU
     -4  -- bad symmetryflag for mtxLDU
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
     -14 -- mtxA is NULL
     -15 -- types of mtxLDU and mtxA are not the same
     -16 -- mtxA is not in chevron mode
     -17 -- sigma < 0.0
     -18 -- msglvl > 0 and msgFile is NULL
     -19 -- singular pivot found

   created -- 98oct03, cca
   -------------------------------------------------------------------
*/
int
ILUMtx_factor (
   ILUMtx   *mtxLDU,
   InpMtx   *mtxA,
   double   sigma,
   double   *pops,
   int      msglvl,
   FILE     *msgFile
) {
double   ops ;
double   *entD, *entL, *entU, *pDii, *pLji, *pUij, *temp ;
double   **p_entL, **p_entU ;
int      countL, countU, ieqn, ii, ioff, jeqn, keqn, neqns, nexteqn, 
         nkeep, rc, sizeL, sizeU, symmetryflag, type ;
int      *headL, *headU, *indicesL, *indicesU, *indL, *indU, *linkL, 
         *linkU, *map, *markL, *markU, *offsetsL, *offsetsU, *sizesL, 
         *sizesU ;
int      **p_indL, **p_indU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxLDU == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), mtxLDU = NULL\n") ;
   return(-1) ;
}
if ( (neqns = mtxLDU->neqns) <= 0 ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n neqns = %d\n", neqns) ;
   return(-2) ;
}
if ( !(ILUMTX_IS_REAL(mtxLDU) || ILUMTX_IS_COMPLEX(mtxLDU)) ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n type = %d\n", mtxLDU->type) ;
   return(-3) ;
}
if ( !(ILUMTX_IS_SYMMETRIC(mtxLDU) || ILUMTX_IS_HERMITIAN(mtxLDU)
       || ILUMTX_IS_NONSYMMETRIC(mtxLDU)) ) {
   fprintf(stderr, "\n error in ILUMfactor()"
           "\n mtxLDU symmetry = %d\n", mtxLDU->symmetryflag) ;
   return(-4) ;
}
if ( !(ILUMTX_IS_L_BY_ROWS(mtxLDU) 
       || ILUMTX_IS_L_BY_COLUMNS(mtxLDU)) ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n LstorageMode = %d\n", mtxLDU->LstorageMode) ;
   return(-5) ;
}
if ( !(ILUMTX_IS_U_BY_ROWS(mtxLDU) 
       || ILUMTX_IS_U_BY_COLUMNS(mtxLDU)) ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n UstorageMode = %d\n", mtxLDU->UstorageMode) ;
   return(-6) ;
}
type = mtxLDU->type ;
symmetryflag = mtxLDU->symmetryflag ;
sizesU = mtxLDU->sizesU ;
entD   = mtxLDU->entD   ;
p_indU = mtxLDU->p_indU ;
p_entU = mtxLDU->p_entU ;
if ( ILUMTX_IS_SYMMETRIC(mtxLDU) || ILUMTX_IS_HERMITIAN(mtxLDU) ) {
   sizesL = mtxLDU->sizesU ;
   p_indL = mtxLDU->p_indU ;
   p_entL = mtxLDU->p_entU ;
} else {
   sizesL = mtxLDU->sizesL ;
   p_indL = mtxLDU->p_indL ;
   p_entL = mtxLDU->p_entL ;
}
if ( sizesL == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), sizesL = NULL\n") ;
   return(-7) ;
}
if ( sizesU == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), sizesU = NULL\n") ;
   return(-8) ;
}
if ( p_indL == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), p_indL = NULL\n") ;
   return(-9) ;
}
if ( p_indU == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), p_indU = NULL\n") ;
   return(-10) ;
}
if ( entD == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), entD = NULL\n") ;
   return(-11) ;
}
if ( p_entL == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), p_entL = NULL\n") ;
   return(-12) ;
}
if ( p_entU == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), p_entU = NULL\n") ;
   return(-13) ;
}
if ( mtxA == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor(), mtxA = NULL\n") ;
   return(-14) ;
}
if ( mtxA->inputMode != (type = mtxLDU->type) ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n mtxA type = %d, mtxLDU type = %d\n",
           mtxA->inputMode, mtxLDU->type) ;
   return(-15) ;
}
if ( ! INPMTX_IS_BY_CHEVRONS(mtxA) ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n mtxA must be in chevron mode\n") ;
   return(-16) ;
}
if ( sigma < 0.0 ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n sigma = %f\n", sigma) ;
   return(-17) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n error in ILUM_factor()"
           "\n msglvl = %d, msgFile is NULL\n", msglvl) ;
   return(-18) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------
   allocate temporary storage
   --------------------------
*/
indicesU = IVinit(neqns, -1) ;
indicesL = indicesU ;
markU    = IVinit(neqns, -1) ;
headU    = IVinit(neqns, -1) ;
linkU    = IVinit(neqns, -1) ;
offsetsU = IVinit(neqns, -1) ;
map      = IVinit(neqns, -1) ;
if ( ILUMTX_IS_SYMMETRIC(mtxLDU) || ILUMTX_IS_HERMITIAN(mtxLDU) ) {
   headL    = headU ;
   linkL    = linkU ;
   offsetsL = offsetsU ;
   markL    = NULL ;
} else {
   headL    = IVinit(neqns, -1) ;
   linkL    = IVinit(neqns, -1) ;
   offsetsL = IVinit(neqns, -1) ;
   markL    = IVinit(neqns, -1) ;
}
if ( type == SPOOLES_REAL ) {
   temp = DVinit(neqns, 0.0) ;
} else {
   temp = DVinit(2*neqns, 0.0) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------
   load diagonal entries of A
   --------------------------
*/
loadDiagonal(type, neqns, mtxA, entD, msglvl, msgFile) ;
#if MYDEBUG > 0
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n entD after loading diag(A)") ;
   DVfprintf(msgFile, neqns, entD) ;
   fflush(msgFile) ;
}
#endif
/*--------------------------------------------------------------------*/
/*
   -----------------------
   loop over the equations
   -----------------------
*/
rc  = 1 ;
ops = 0.0 ;
for ( jeqn = 0 ; jeqn < neqns ; jeqn++ ) {
#if  MYDEBUG > 0
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n ### working on equation %d", jeqn) ;
      fflush(msgFile) ;
   }
#endif
/*
   --------------------------------------------------
   determine the structure of entries in the row of U
   --------------------------------------------------
*/
   countU = getRowStructure(jeqn, indicesU, mtxA, sizesU, p_indU, 
                            headL, linkL, markU, msglvl, msgFile) ;
#if  MYDEBUG > 0
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n symbolic row structure") ;
      IVfprintf(msgFile, countU, indicesU) ;
   }
#endif
   for ( ii = 0 ; ii < countU ; ii++ ) {
      map[indicesU[ii]] = ii ;
   }
/*
   -------------------------------------------------
   load original entries of row into the temp vector
   -------------------------------------------------
*/
   loadOrigRow(jeqn, mtxA, map, temp, msglvl, msgFile) ;
#if  MYDEBUG > 0
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n original row entries loaded") ;
      IVfprintf(msgFile, countU, indicesU) ;
      if ( type == SPOOLES_REAL ) {
         DVfprintf(msgFile, countU, temp) ;
      } else {
         DVfprintf(msgFile, 2*countU, temp) ;
      }
   }
#endif
/*
   ------------------
   get updates into U
   ------------------
*/
   for ( ieqn = headL[jeqn] ; ieqn != -1 ; ieqn = nexteqn ) {
      sizeL   = sizesL[ieqn] ;
      indL    = p_indL[ieqn] ;
      entL    = p_entL[ieqn] ;
      ioff    = offsetsL[ieqn] ;
      sizeU   = sizesU[ieqn] ;
      indU    = p_indU[ieqn] ;
      entU    = p_entU[ieqn] ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
                 "\n\n  ## update from %d, sizeL %d, ioff %d, sizeU %d",
                 ieqn, sizeL, ioff, sizeU) ;
         fflush(msgFile) ;
      }
#endif
/*
      ------------------------------
      update row jeqn using row ieqn
      ------------------------------
*/
      if ( type == SPOOLES_REAL ) {
         pLji = entL + ioff ;
         pDii = entD + ieqn ;
      } else {
         pLji = entL + 2*ioff ;
         pDii = entD + 2*ieqn ;
      }
      ops += updateRow(type, symmetryflag, jeqn, pLji, pDii, 
                       sizeU, indU, entU, map, temp, msglvl, msgFile) ;
#if  MYDEBUG > 0
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n after update from %d", ieqn) ;
         IVfprintf(msgFile, countU, indicesU) ;
         if ( type == SPOOLES_REAL ) {
            DVfprintf(msgFile, countU, temp) ;
         } else {
            DVfprintf(msgFile, 2*countU, temp) ;
         }
         fflush(msgFile) ;
      }
#endif
/*
      ------------------------------------
      link ieqn to the next row it updates
      ------------------------------------
*/
      nexteqn = linkL[ieqn] ; linkL[ieqn] = -1 ;
      if ( ++ioff < sizeL ) {
         offsetsL[ieqn] = ioff ;
         keqn = indL[ioff] ;
         linkL[ieqn] = headL[keqn] ;
         headL[keqn] = ieqn ;
#if  MYDEBUG > 0
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n linking ieqn %d to keqn %d", 
                    ieqn, keqn) ;
            fflush(msgFile) ;
         }
#endif
      }
   }
/*
   ----------------------
   check for a zero pivot
   ----------------------
*/
   if ( type == SPOOLES_REAL ) {
      if ( temp[0] == 0.0 ) {
         rc = -19 ;
         break ;
      }
   } else {
      if ( temp[0] == 0.0 && temp[1] == 0.0 ) {
         rc = -19 ;
         break ;
      }
   }
/*
   --------------------------
   extract row jeqn and store
   --------------------------
*/
   nkeep = storeRow(jeqn, type, sigma, countU, indicesU, temp, 
                    entD, sizesU, p_indU, p_entU, msglvl, msgFile) ;
#if  MYDEBUG > 0
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n keep %d entries from row %d", nkeep, jeqn) ;
   }
#endif
   if ( nkeep > 0 ) {
/*
      -------------------------------------------------
      link jeqn to the first column of L it must update
      -------------------------------------------------
*/
      keqn = indicesU[0] ;
#if  MYDEBUG > 0
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n U linking jeqn %d to keqn %d", 
                 jeqn, keqn) ;
         fflush(msgFile) ;
      }
#endif
      linkU[jeqn] = headU[keqn] ;
      headU[keqn] = jeqn ;
      offsetsU[jeqn] = 0 ;
   }
   if ( symmetryflag == SPOOLES_NONSYMMETRIC ) {
/*
      -----------------------------------------------------
      determine the structure of entries in the column of U
      -----------------------------------------------------
*/
      countL = getColumnStructure(jeqn, indicesL, mtxA, sizesL, p_indL, 
                                  headU, linkU, markL, msglvl, msgFile);
      for ( ii = 0 ; ii < countL ; ii++ ) {
         map[indicesU[ii]] = ii ;
      }
/*
      ----------------------------------------------------
      load original entries of column into the temp vector
      ----------------------------------------------------
*/
      loadOrigColumn(jeqn, mtxA, map, temp, msglvl, msgFile) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n after loading original entries") ;
         if ( type == SPOOLES_REAL ) {
            DVfprintf(msgFile, countL, temp) ;
         } else {
            DVfprintf(msgFile, 2*countL, temp) ;
         }
         fflush(msgFile) ;
      }
#endif
/*
      ------------------
      get updates into L
      ------------------
*/
      for ( ieqn = headU[jeqn] ; ieqn != -1 ; ieqn = nexteqn ) {
         sizeU   = sizesU[ieqn] ;
         indU    = p_indU[ieqn] ;
         entU    = p_entU[ieqn] ;
         ioff    = offsetsU[ieqn] ;
         sizeL   = sizesL[ieqn] ;
         indL    = p_indL[ieqn] ;
         entL    = p_entL[ieqn] ;
/*
         ------------------------------------
         update column jeqn using column ieqn
         ------------------------------------
*/
         if ( type == SPOOLES_REAL ) {
            pUij = entU + ioff ;
            pDii = entD + ieqn ;
         } else {
            pUij = entU + 2*ioff ;
            pDii = entD + 2*ieqn ;
         }
         ops += updateColumn(type, symmetryflag, jeqn, pUij, pDii, 
                        sizeL, indL, entL, map, temp, msglvl, msgFile) ;
#if  MYDEBUG > 0
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n after update from %d", ieqn) ;
            if ( type == SPOOLES_REAL ) {
               DVfprintf(msgFile, countL, temp) ;
            } else {
               DVfprintf(msgFile, 2*countL, temp) ;
            }
            fflush(msgFile) ;
         }
#endif
/*
         ---------------------------------------
         link ieqn to the next column it updates
         ---------------------------------------
*/
         nexteqn = linkU[ieqn] ; linkU[ieqn] = -1 ;
         if ( ++ioff < sizeU ) {
            offsetsU[ieqn] = ioff ;
            keqn = indU[ioff] ;
            linkU[ieqn] = headU[keqn] ;
            headU[keqn] = ieqn ;
         }
      }
/*
      -----------------------------
      extract column jeqn and store
      -----------------------------
*/
      nkeep = storeColumn(jeqn, type, sigma, countL, indicesL, temp, 
                        entD, sizesL, p_indL, p_entL, msglvl, msgFile) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
                 "\n keep %d entries from column %d", nkeep, jeqn) ;
      }
#endif
      if ( nkeep > 0 ) {
/*
         ----------------------------------------------
         link jeqn to the first row of U it must update
         ----------------------------------------------
*/
         keqn = indicesL[0] ;
#if  MYDEBUG > 0
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n L linking jeqn %d to keqn %d", 
                    jeqn, keqn) ;
            fflush(msgFile) ;
         }
#endif
         linkL[jeqn] = headL[keqn] ;
         headL[keqn] = jeqn ;
         offsetsL[jeqn] = 0 ;
      }
   }
}
/*
   -----------------------------------------
   increment the operation count if not NULL
   -----------------------------------------
*/
if ( pops != NULL ) {
   *pops += ops ;
}
/*
   ----------------------
   free temporary storage
   ----------------------
*/
IVfree(indicesU) ;
IVfree(markU) ;
IVfree(headU) ;
IVfree(linkU) ;
IVfree(offsetsU) ;
if ( ILUMTX_IS_NONSYMMETRIC(mtxLDU) ) {
   IVfree(headL) ;
   IVfree(linkL) ;
   IVfree(offsetsL) ;
   IVfree(markL) ;
}
IVfree(map) ;
DVfree(temp) ;

return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- copy the diagonal entries of A into entD[]

   created -- 98oct08, cca
   -----------------------------------------------------
*/
static void
loadDiagonal (
   int      type,
   int      neqns,
   InpMtx   *mtxA,
   double   entD[],
   int      msglvl,
   FILE     *msgFile
) {
double   *entA ;
int      ii, jeqn, size ;
int      *indA ;

if ( type == SPOOLES_REAL ) {
   for ( jeqn = 0 ; jeqn < neqns ; jeqn++ ) {
      InpMtx_realVector(mtxA, jeqn, &size, &indA, &entA) ;
      for ( ii = 0 ; ii < size ; ii++ ) {
         if ( indA[ii] == 0 ) {
            entD[jeqn] = entA[ii] ;
            break ;
         }
      }
   }
} else {
   for ( jeqn = 0 ; jeqn < neqns ; jeqn++ ) {
      InpMtx_complexVector(mtxA, jeqn, &size, &indA, &entA) ;
      for ( ii = 0 ; ii < size ; ii++ ) {
         if ( indA[ii] == 0 ) {
            entD[2*jeqn]   = entA[2*ii]   ;
            entD[2*jeqn+1] = entA[2*ii+1] ;
            break ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   fill the structure of row jeqn of U into indicesU[].
   note: includes the diagonal entry.
   return value is the number of indices in indicesU.

   created -- 98oct04, cca
   ----------------------------------------------------
*/
static int
getRowStructure (
   int      jeqn,
   int      indicesU[],
   InpMtx   *mtxA,
   int      sizesU[],
   int      *p_indU[],
   int      headL[],
   int      linkL[],
   int      markU[],
   int      msglvl,
   FILE     *msgFile
) {
int   countU, ieqn, ii, keqn, sizeA, sizeU ;
int   *indA, *indU ;
/*
   --------------------------------------------------
   construct the structure of entries in the row of U
   --------------------------------------------------
*/
InpMtx_vector(mtxA, jeqn, &sizeA, &indA) ;
markU[jeqn] = jeqn ;
indicesU[0] = jeqn ;
countU      =   1  ;
#if  MYDEBUG > 0
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n indicesU[%d] = %d", 
           countU-1, indicesU[countU-1]) ;
   fflush(msgFile) ;
}
#endif
for ( ii = sizeA - 1 ; ii >= 0 ; ii-- ) {
   keqn = jeqn + indA[ii] ;
   if ( keqn > jeqn ) {
      markU[keqn] = jeqn ;
      indicesU[countU++] = keqn ;
#if  MYDEBUG > 0
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n 1. indicesU[%d] = %d", 
                 countU-1, indicesU[countU-1]) ;
         fflush(msgFile) ;
      }
#endif
   } else {
      break ;
   }
}
for ( ieqn = headL[jeqn] ; ieqn != -1 ; ieqn = linkL[ieqn] ) {
   sizeU = sizesU[ieqn] ;
   indU  = p_indU[ieqn] ;
#if  MYDEBUG > 0
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n checking out ieqn %d", ieqn) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = sizeU - 1 ; ii >= 0 ; ii-- ) {
#if  MYDEBUG > 0
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n    indU[%d] = %d", ii, indU[ii]) ;
         fflush(msgFile) ;
      }
#endif
      if ( (keqn = indU[ii]) > jeqn ) {
         if ( markU[keqn] != jeqn ) {
            markU[keqn] = jeqn ;
            indicesU[countU++] = keqn ;
#if  MYDEBUG > 0
            if ( msglvl > 3 ) {
               fprintf(msgFile, "\n 2. indicesU[%d] = %d", 
                       countU-1, indicesU[countU-1]) ;
               fflush(msgFile) ;
            }
#endif
         }
      } else {
         break ;
      }
   }
}
IVqsortUp(countU, indicesU) ;

return(countU) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   load entries in upper row jeqn of A into temp vector

   created -- 98oct04, cca
   ----------------------------------------------------
*/
static void
loadOrigRow (
   int      jeqn,
   InpMtx   *mtxA,
   int      map[],
   double   temp[],
   int      msglvl,
   FILE     *msgFile
) {
double   *entA ;
int      ii, keqn, kloc, sizeA ;
int      *indA ;
/*
   -------------------------------------------------
   load original entries of row into the temp vector
   -------------------------------------------------
*/
if ( INPMTX_IS_REAL_ENTRIES(mtxA) ) {
   InpMtx_realVector(mtxA, jeqn, &sizeA, &indA, &entA) ;
   for ( ii = sizeA - 1 ; ii >= 0 ; ii-- ) {
      keqn = jeqn + indA[ii] ;
      if ( keqn >= jeqn ) {
         kloc = map[keqn] ;
         temp[kloc] = entA[ii] ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e", kloc, entA[ii]) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
} else {
   InpMtx_complexVector(mtxA, jeqn, &sizeA, &indA, &entA) ;
   for ( ii = sizeA - 1 ; ii >= 0 ; ii-- ) {
      keqn = jeqn + indA[ii] ;
      if ( keqn >= jeqn ) {
         kloc = map[keqn] ;
         temp[2*kloc]   = entA[2*ii]   ;
         temp[2*kloc+1] = entA[2*ii+1] ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e + i*%12.4e", 
                    kloc, entA[2*ii], entA[2*ii+1]) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- perform an update to row jeqn

   return value ---
     # of ops performed

   created -- 98oct04, cca
   ----------------------------------------
*/
static double
updateRow (
   int      type,
   int      symmetryflag,
   int      jeqn,
   double   LJI[],
   double   DII[],
   int      sizeU,
   int      indU[],
   double   entU[],
   int      map[],
   double   temp[],
   int      msglvl,
   FILE     *msgFile
) {
double   ops = 0.0 ;
int      ii, keqn, kloc ;
/*
   ------------------------------------------
   compute the multiplier and update row jeqn
   ------------------------------------------
*/
if ( type == SPOOLES_REAL ) {
   double   Lji ;

   Lji = LJI[0] * DII[0] ;
#if  MYDEBUG > 0
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n LJI = %12.4e, DII = %12.4e, Lji = %12.4e", 
              LJI[0], DII[0], Lji) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = sizeU - 1 ; ii >= 0 ; ii-- ) {
      if ( (keqn = indU[ii]) >= jeqn ) {
         kloc = map[keqn] ;
         temp[kloc] -= Lji * entU[ii] ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] -= Lji * %12.4e",
                    kloc, entU[ii]) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
   ops += 2*(sizeU - ii) ;
} else {
   double   LjiI, LjiR, UikI, UikR ;

   if ( symmetryflag == SPOOLES_HERMITIAN ) {
      LjiR = LJI[0] * DII[0] + LJI[1] * DII[1] ;
      LjiI = LJI[0] * DII[1] - LJI[1] * DII[0] ;
   } else {
      LjiR = LJI[0] * DII[0] - LJI[1] * DII[1] ;
      LjiI = LJI[0] * DII[1] + LJI[1] * DII[0] ;
   }
#if  MYDEBUG > 0
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n Lji = %12.4e + i*%12.4e", LjiR, LjiI) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = sizeU - 1 ; ii >= 0 ; ii-- ) {
      if ( (keqn = indU[ii]) >= jeqn ) {
         kloc = map[keqn] ;
         UikR = entU[2*ii] ; UikI = entU[2*ii+1] ;
         temp[2*kloc]   -= LjiR*UikR - LjiI*UikI ;
         temp[2*kloc+1] -= LjiR*UikI + LjiI*UikR ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] -= Lji * (%12.4e + i*%12.4e)",
                    kloc, UikR, UikI) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
   ops += 8*(sizeU - ii) ;
}
return(ops) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to extract and store the entries in row jeqn
      note: indicesU[0] = jeqn.

   return value --
      number of entries kept

   created -- 98oct04, cca
   -------------------------------------------------------
*/
static int
storeRow (
   int      jeqn,
   int      type,
   double   sigma,
   int      countU,
   int      indicesU[],
   double   temp[],
   double   entD[],
   int      sizesU[],
   int      *p_indU[],
   double   *p_entU[],
   int      msglvl,
   FILE     *msgFile
) {
double   magDjj, magDkk, magUjk ;
int      ii, keqn, nkeep ;

nkeep = 0 ;
if ( type == SPOOLES_REAL ) {
/*
   ------------------------
   store the diagonal entry
   ------------------------
*/
   entD[jeqn] = temp[0] ;
#if  MYDEBUG > 0
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n entD[%d] = %12.4e", jeqn, entD[jeqn]) ;
      fflush(msgFile) ;
   }
#endif
/*
   --------------------------------------------------
   check out the entries, slide down those to be kept
   --------------------------------------------------
*/
   magDjj = fabs(entD[jeqn]) ;
   for ( ii = 1 ; ii < countU ; ii++ ) {
      keqn   = indicesU[ii] ;
      magUjk = fabs(temp[ii]) ;
      magDkk = fabs(entD[keqn]) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n keqn %d, magUjk %12.4e, magDkk %12.4e",
                 keqn, magUjk, magDkk) ;
         fflush(msgFile) ;
      }
#endif
      if ( magUjk*magUjk > sigma*magDjj*magDkk ) {
         indicesU[nkeep] = keqn ;
         temp[nkeep]     = temp[ii] / entD[jeqn] ;
#if  MYDEBUG > 0
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e",
                    nkeep, temp[nkeep]) ;
            fflush(msgFile) ;
         }
#endif
         nkeep++ ;
      }
   }
   if ( nkeep > 0 ) {
/*
      ------------------------------------
      store indices and entries to be kept
      ------------------------------------
*/
      sizesU[jeqn] = nkeep ;
      p_indU[jeqn] = IVinit(nkeep, -1) ;
      IVcopy(nkeep, p_indU[jeqn], indicesU) ;
      p_entU[jeqn] = DVinit(nkeep, -1) ;
      DVcopy(nkeep, p_entU[jeqn], temp) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n row %d factor indices", jeqn) ;
         IVfprintf(msgFile, nkeep, p_indU[jeqn]) ;
         fprintf(msgFile, "\n row %d factor entries", jeqn) ;
         DVfprintf(msgFile, nkeep, p_entU[jeqn]) ;
         fflush(msgFile) ;
      }
#endif
   }
/*
   --------------------
   zero the temp vector
   --------------------
*/
   DVzero(countU, temp) ;
} else {
   double   rI, rR, tI, tR ;
/*
   ------------------------
   store the diagonal entry
   ------------------------
*/
   entD[2*jeqn]   = temp[0] ;
   entD[2*jeqn+1] = temp[1] ;
#if  MYDEBUG > 0
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n entD[%d] = %12.4e + i*%12.4e", 
              jeqn, entD[2*jeqn], entD[2*jeqn+1]) ;
      fflush(msgFile) ;
   }
#endif
/*
   --------------------------------------------------
   check out the entries, slide down those to be kept
   --------------------------------------------------
*/
   magDjj = Zabs(entD[2*jeqn], entD[2*jeqn+1]) ;
   Zrecip(entD[2*jeqn], entD[2*jeqn+1], &rR, &rI) ;
#if  MYDEBUG > 0
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n rR = %12.4e, rI = %12.4e", rR, rI) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = 1 ; ii < countU ; ii++ ) {
      keqn = indicesU[ii] ;
      magUjk = Zabs(temp[2*ii],temp[2*ii+1]) ;
      magDkk = Zabs(entD[2*keqn], entD[2*keqn+1]) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n keqn %d, magUjk %12.4e, magDkk %12.4e",
                 keqn, magUjk, magDkk) ;
         fflush(msgFile) ;
      }
#endif
      if ( magUjk*magUjk > sigma*magDjj*magDkk ) {
         indicesU[nkeep] = keqn ;
         tR = temp[2*ii] ; tI = temp[2*ii+1] ;
#if  MYDEBUG > 0
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n tR = %12.4e, tI = %12.4e", tR, tI) ;
            fflush(msgFile) ;
         }
#endif
         temp[2*nkeep]   = tR*rR - tI*rI ;
         temp[2*nkeep+1] = tR*rI + tI*rR ;
#if  MYDEBUG > 0
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e + i*%12.4e",
                    nkeep, temp[2*nkeep], temp[2*nkeep+1]) ;
            fflush(msgFile) ;
         }
#endif
         nkeep++ ;
      }
   }
   if ( nkeep > 0 ) {
/*
      ------------------------------------
      store indices and entries to be kept
      ------------------------------------
*/
      sizesU[jeqn] = nkeep ;
      p_indU[jeqn] = IVinit(nkeep, -1) ;
      IVcopy(nkeep, p_indU[jeqn], indicesU) ;
      p_entU[jeqn] = DVinit(2*nkeep, -1) ;
      DVcopy(2*nkeep, p_entU[jeqn], temp) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n row %d factor indices", jeqn) ;
         IVfprintf(msgFile, nkeep, p_indU[jeqn]) ;
         fprintf(msgFile, "\n row %d factor entries", jeqn) ;
         DVfprintf(msgFile, 2*nkeep, p_entU[jeqn]) ;
         fflush(msgFile) ;
      }
#endif
   }
/*
   --------------------
   zero the temp vector
   --------------------
*/
   DVzero(2*countU, temp) ;
}
return(nkeep) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   fill the structure of column jeqn of L into indicesL[].
   note: does not include the diagonal entry.
   return value is the number of indices in indicesL.

   created -- 98oct04, cca
   -------------------------------------------------------
*/
static int
getColumnStructure (
   int      jeqn,
   int      indicesL[],
   InpMtx   *mtxA,
   int      sizesL[],
   int      *p_indL[],
   int      headU[],
   int      linkU[],
   int      markL[],
   int      msglvl,
   FILE     *msgFile
) {
int   countL, ieqn, ii, keqn, sizeA, sizeL ;
int   *indA, *indL ;
/*
   -----------------------------------------------------
   construct the structure of entries in the column of L
   -----------------------------------------------------
*/
InpMtx_vector(mtxA, jeqn, &sizeA, &indA) ;
countL =   0  ;
for ( ii = 0 ; ii < sizeA ; ii++ ) {
   keqn = jeqn - indA[ii] ;
   if ( keqn > jeqn ) {
      markL[keqn] = jeqn ;
      indicesL[countL++] = keqn ;
#if  MYDEBUG > 0
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n 1. indicesL[%d] = %d", 
                 countL-1, indicesL[countL-1]) ;
         fflush(msgFile) ;
      }
#endif
   } else {
      break ;
   }
}
for ( ieqn = headU[jeqn] ; ieqn != -1 ; ieqn = linkU[ieqn] ) {
   sizeL = sizesL[ieqn] ;
   indL  = p_indL[ieqn] ;
   for ( ii = sizeL - 1 ; ii >= 0 ; ii-- ) {
      if ( (keqn = indL[ii]) > jeqn ) {
         if ( markL[keqn] != jeqn ) {
            markL[keqn] = jeqn ;
            indicesL[countL++] = keqn ;
#if  MYDEBUG > 0
            if ( msglvl > 3 ) {
               fprintf(msgFile, "\n 2. indicesL[%d] = %d", 
                       countL-1, indicesL[countL-1]) ;
               fflush(msgFile) ;
            }
#endif
         }
      } else {
         break ;
      }
   }
}
IVqsortUp(countL, indicesL) ;

return(countL) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   load entries in upper column jeqn of A into temp vector

   created -- 98oct04, cca
   -------------------------------------------------------
*/
static void
loadOrigColumn (
   int      jeqn,
   InpMtx   *mtxA,
   int      map[],
   double   temp[],
   int      msglvl,
   FILE     *msgFile
) {
double   *entA ;
int      ii, keqn, kloc, sizeA ;
int      *indA ;
/*
   ----------------------------------------------------
   load original entries of column into the temp vector
   ----------------------------------------------------
*/
if ( INPMTX_IS_REAL_ENTRIES(mtxA) ) {
   InpMtx_realVector(mtxA, jeqn, &sizeA, &indA, &entA) ;
   for ( ii = 0 ; ii < sizeA ; ii++ ) {
      keqn = jeqn - indA[ii] ;
      if ( keqn > jeqn ) {
         kloc = map[keqn] ;
         temp[kloc] = entA[ii] ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e", kloc, entA[ii]) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
} else {
   InpMtx_complexVector(mtxA, jeqn, &sizeA, &indA, &entA) ;
   for ( ii = 0 ; ii < sizeA ; ii++ ) {
      keqn = jeqn - indA[ii] ;
      if ( keqn > jeqn ) {
         kloc = map[keqn] ;
         temp[2*kloc]   = entA[2*ii]   ;
         temp[2*kloc+1] = entA[2*ii+1] ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e + i*%12.4e", 
                    kloc, entA[2*ii], entA[2*ii+1]) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- perform an update to column jeqn

   return value ---
     # of ops performed

   created -- 98oct04, cca
   -------------------------------------------
*/
static double
updateColumn (
   int      type,
   int      symmetryflag,
   int      jeqn,
   double   DII[],
   double   UIJ[],
   int      sizeL,
   int      indL[],
   double   entL[],
   int      map[],
   double   temp[],
   int      msglvl,
   FILE     *msgFile
) {
double   ops = 0.0 ;
int      ii, keqn, kloc ;
/*
   ------------------------------------------
   compute the multiplier and update row jeqn
   ------------------------------------------
*/
if ( type == SPOOLES_REAL ) {
   double   Uij ;

   Uij = DII[0] * UIJ[0] ;
#if  MYDEBUG > 0
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n Uij = %12.4e", Uij) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = sizeL - 1 ; ii >= 0 ; ii-- ) {
      if ( (keqn = indL[ii]) > jeqn ) {
         kloc = map[keqn] ;
         temp[kloc] -= Uij * entL[ii] ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] -= Uij * %12.4e",
                    kloc, entL[ii]) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
   ops += 2*(sizeL - ii) ;
} else {
   double   LkiI, LkiR, UijI, UijR ;

   UijR = UIJ[0] * DII[0] - UIJ[1] * DII[1] ;
   UijI = UIJ[0] * DII[1] + UIJ[1] * DII[0] ;
#if  MYDEBUG > 0
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n Uij = %12.4e + i*%12.4e", UijR, UijI) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = sizeL - 1 ; ii >= 0 ; ii-- ) {
      if ( (keqn = indL[ii]) > jeqn ) {
         kloc = map[keqn] ;
         LkiR = entL[2*ii] ; LkiI = entL[2*ii+1] ;
         temp[2*kloc]   -= LkiR*UijR - LkiI*UijI ;
         temp[2*kloc+1] -= LkiR*UijI + LkiI*UijR ;
#if  MYDEBUG > 0
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n temp[%d] -= Uij * (%12.4e + i*%12.4e)",
                    kloc, LkiR, LkiI) ;
            fflush(msgFile) ;
         }
#endif
      } else {
         break ;
      }
   }
   ops += 8*(sizeL - ii) ;
}
return(ops) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- to extract and store the entries in column jeqn

   return value --
      number of entries kept

   created -- 98oct04, cca
   ----------------------------------------------------------
*/
static int
storeColumn (
   int      jeqn,
   int      type,
   double   sigma,
   int      countL,
   int      indicesL[],
   double   temp[],
   double   entD[],
   int      sizesL[],
   int      *p_indL[],
   double   *p_entL[],
   int      msglvl,
   FILE     *msgFile
) {
double   magDjj, magDkk, magLkj ;
int      ii, keqn, nkeep ;

nkeep = 0 ;
if ( type == SPOOLES_REAL ) {
/*
   --------------------------------------------------
   check out the entries, slide down those to be kept
   --------------------------------------------------
*/
   magDjj = fabs(entD[jeqn]) ;
   for ( ii = 0 ; ii < countL ; ii++ ) {
      keqn   = indicesL[ii] ;
      magLkj = fabs(temp[ii]) ;
      magDkk = fabs(entD[keqn]) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n keqn %d, magLkj %12.4e, magDkk %12.4e",
                 keqn, magLkj, magDkk) ;
         fflush(msgFile) ;
      }
#endif
      if ( magLkj*magLkj > sigma*magDjj*magDkk ) {
         indicesL[nkeep] = keqn ;
         temp[nkeep]     = temp[ii] / entD[jeqn] ;
#if  MYDEBUG > 0
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e",
                    nkeep, temp[nkeep]) ;
            fflush(msgFile) ;
         }
#endif
         nkeep++ ;
      }
   }
   if ( nkeep > 0 ) {
/*
      ------------------------------------
      store indices and entries to be kept
      ------------------------------------
*/
      sizesL[jeqn] = nkeep ;
      p_indL[jeqn] = IVinit(nkeep, -1) ;
      IVcopy(nkeep, p_indL[jeqn], indicesL) ;
      p_entL[jeqn] = DVinit(nkeep, -1) ;
      DVcopy(nkeep, p_entL[jeqn], temp) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n column %d factor indices", jeqn) ;
         IVfprintf(msgFile, nkeep, p_indL[jeqn]) ;
         fprintf(msgFile, "\n column %d factor entries", jeqn) ;
         DVfprintf(msgFile, nkeep, p_entL[jeqn]) ;
         fflush(msgFile) ;
      }
#endif
   }
/*
   --------------------
   zero the temp vector
   --------------------
*/
   DVzero(countL, temp) ;
} else {
   double   rI, rR, tI, tR ;
/*
   --------------------------------------------------
   check out the entries, slide down those to be kept
   --------------------------------------------------
*/
   magDjj = Zabs(entD[2*jeqn], entD[2*jeqn+1]) ;
   Zrecip(entD[2*jeqn], entD[2*jeqn+1], &rR, &rI) ;
#if  MYDEBUG > 0
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n D(%d,%d) = %12.4e + i*%12.4e"
              "\n 1/(D(%d,%d) = %12.4e + i*%12.4e",
              jeqn, jeqn, entD[2*jeqn], entD[2*jeqn+1], 
              jeqn, jeqn, rR, rI) ;
      fflush(msgFile) ;
   }
#endif
   for ( ii = 0 ; ii < countL ; ii++ ) {
      keqn   = indicesL[ii] ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n A(%d,%d) = %12.4e + i*%12.4e",
                 keqn, jeqn, temp[2*ii], temp[2*ii+1]) ;
         fflush(msgFile) ;
      }
#endif
      magLkj = Zabs(temp[2*ii],temp[2*ii+1]) ;
      magDkk = Zabs(entD[2*keqn], entD[2*keqn+1]) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n keqn %d, magLkj %12.4e, magDkk %12.4e",
                 keqn, magLkj, magDkk) ;
         fflush(msgFile) ;
      }
#endif
      if ( magLkj*magLkj > sigma*magDjj*magDkk ) {
         indicesL[nkeep] = keqn ;
         tR = temp[2*ii] ; tI = temp[2*ii+1] ;
         temp[2*nkeep]   = tR*rR - tI*rI ;
         temp[2*nkeep+1] = tR*rI + tI*rR ;
#if  MYDEBUG > 0
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n temp[%d] = %12.4e + i*%12.4e",
                    nkeep, temp[2*nkeep], temp[2*nkeep+1]) ;
            fflush(msgFile) ;
         }
#endif
         nkeep++ ;
      }
   }
   if ( nkeep > 0 ) {
/*
      ------------------------------------
      store indices and entries to be kept
      ------------------------------------
*/
      sizesL[jeqn] = nkeep ;
      p_indL[jeqn] = IVinit(nkeep, -1) ;
      IVcopy(nkeep, p_indL[jeqn], indicesL) ;
      p_entL[jeqn] = DVinit(2*nkeep, -1) ;
      DVcopy(2*nkeep, p_entL[jeqn], temp) ;
#if  MYDEBUG > 0
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n column %d factor indices", jeqn) ;
         IVfprintf(msgFile, nkeep, p_indL[jeqn]) ;
         fprintf(msgFile, "\n column %d factor entries", jeqn) ;
         DVfprintf(msgFile, 2*nkeep, p_entL[jeqn]) ;
         fflush(msgFile) ;
      }
#endif
   }
/*
   --------------------
   zero the temp vector
   --------------------
*/
   DVzero(2*countL, temp) ;
}
return(nkeep) ; }

/*--------------------------------------------------------------------*/
