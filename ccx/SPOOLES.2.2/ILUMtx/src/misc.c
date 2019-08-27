/*  misc.c  */

#include "../ILUMtx.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*
   -------------------------------------------------------------
   purpose -- to fill the indices and entries with random values

   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtx
     -4  -- bad symmetryflag for mtx
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL

   created -- 98oct03, cca
   -------------------------------------------------------------
*/
int
ILUMtx_fillRandom (
   ILUMtx   *mtx,
   int      seed
) {
double   *entD ;
double   **p_entL, **p_entU ;
Drand    *drand ;
int      ieqn, neqns, size ;
int      *list, *sizesL, *sizesU ;
int      **p_indL, **p_indU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), mtx = NULL\n") ;
   return(-1) ;
}
if ( (neqns = mtx->neqns) <= 0 ) {
   fprintf(stderr, "\n error in ILUM_fillRandom()"
           "\n neqns = %d\n", neqns) ;
   return(-2) ;
}
if ( !(ILUMTX_IS_REAL(mtx) || ILUMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n error in ILUM_fillRandom()"
           "\n type = %d\n", mtx->type) ;
   return(-3) ;
}
if ( !(ILUMTX_IS_SYMMETRIC(mtx) || ILUMTX_IS_HERMITIAN(mtx)
       || ILUMTX_IS_NONSYMMETRIC(mtx)) ) {
   fprintf(stderr, "\n error in ILUMfillRandom()"
           "\n mtx symmetry = %d\n", mtx->symmetryflag) ;
   return(-4) ;
}
if ( !(ILUMTX_IS_L_BY_ROWS(mtx) || ILUMTX_IS_L_BY_COLUMNS(mtx)) ) {
  fprintf(stderr, "\n error in ILUM_fillRandom()"
           "\n LstorageMode = %d\n", mtx->LstorageMode) ;
   return(-5) ;
}
if ( !(ILUMTX_IS_U_BY_ROWS(mtx) || ILUMTX_IS_U_BY_COLUMNS(mtx)) ) {
  fprintf(stderr, "\n error in ILUM_fillRandom()"
           "\n UstorageMode = %d\n", mtx->UstorageMode) ;
   return(-6) ;
}
sizesU = mtx->sizesU ;
entD   = mtx->entD   ;
p_indU = mtx->p_indU ;
p_entU = mtx->p_entU ;
if ( ILUMTX_IS_SYMMETRIC(mtx) || ILUMTX_IS_HERMITIAN(mtx) ) {
   sizesL = mtx->sizesU ;
   p_indL = mtx->p_indU ;
   p_entL = mtx->p_entU ;
} else {
   sizesL = mtx->sizesL ;
   p_indL = mtx->p_indL ;
   p_entL = mtx->p_entL ;
}
if ( sizesL == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), sizesL = NULL\n") ;
   return(-7) ;
}
if ( sizesU == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), sizesU = NULL\n") ;
   return(-8) ;
}
if ( p_indL == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), p_indL = NULL\n") ;
   return(-9) ;
}
if ( p_indU == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), p_indU = NULL\n") ;
   return(-10) ;
}
if ( entD == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), entD = NULL\n") ;
   return(-11) ;
}
if ( p_entL == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), p_entL = NULL\n") ;
   return(-12) ;
}
if ( p_entU == NULL ) {
   fprintf(stderr, "\n error in ILUM_fillRandom(), p_entU = NULL\n") ;
   return(-13) ;
}
/*--------------------------------------------------------------------*/
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
/*
   -----------------
   fill entries in D
   -----------------
*/
Drand_setUniform(drand, 0.0, 1.0) ;
if ( ILUMTX_IS_REAL(mtx) ) {
   Drand_fillDvector(drand, neqns, entD) ;
} else {
   Drand_fillDvector(drand, 2*neqns, entD) ;
}
/*
   -----------------------------
   fill indices and entries in U
   -----------------------------
*/
list = IVinit(neqns, -1) ;
for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
   if ( ILUMTX_IS_U_BY_ROWS(mtx) ) {
      size = neqns - ieqn - 1 ;
   } else {
      size = ieqn ;
   }
   Drand_setUniform(drand, 0, size) ;
   sizesU[ieqn] = (int) Drand_value(drand) ;
   if ( sizesU[ieqn] > 0 ) {
      if ( ILUMTX_IS_U_BY_ROWS(mtx) ) {
         IVramp(size, list, ieqn + 1, 1) ;
      } else {
         IVramp(size, list, 0, 1) ;
      }
      IVshuffle(size, list, ++seed) ;
      IVqsortUp(sizesU[ieqn], list) ;
      p_indU[ieqn] = IVinit(sizesU[ieqn], -1) ;
      IVcopy(sizesU[ieqn], p_indU[ieqn], list) ;
      Drand_setUniform(drand, -1.0, 1.0) ;
      if ( ILUMTX_IS_REAL(mtx) ) {
         size = sizesU[ieqn] ;
      } else {
         size = 2*sizesU[ieqn] ;
      }
      p_entU[ieqn] = DVinit(size, 0.0) ;
      Drand_fillDvector(drand, size, p_entU[ieqn]) ;
   }
}
if ( ILUMTX_IS_NONSYMMETRIC(mtx) ) {
/*
   -----------------------------
   fill indices and entries in L
   -----------------------------
*/
   for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
      if ( ILUMTX_IS_L_BY_ROWS(mtx) ) {
         size = ieqn ;
      } else {
         size = neqns - ieqn - 1 ;
      }
      Drand_setUniform(drand, 0, size) ;
      sizesL[ieqn] = (int) Drand_value(drand) ;
      if ( sizesL[ieqn] > 0 ) {
         if ( ILUMTX_IS_L_BY_ROWS(mtx) ) {
            IVramp(size, list, 0, 1) ;
         } else {
            IVramp(size, list, ieqn + 1, 1) ;
         }
         IVshuffle(size, list, ++seed) ;
         IVqsortUp(sizesL[ieqn], list) ;
         p_indL[ieqn] = IVinit(sizesL[ieqn], -1) ;
         IVcopy(sizesL[ieqn], p_indL[ieqn], list) ;
         Drand_setUniform(drand, -1.0, 1.0) ;
         if ( ILUMTX_IS_REAL(mtx) ) {
            size = sizesL[ieqn] ;
         } else {
            size = 2*sizesL[ieqn] ;
         }
         p_entL[ieqn] = DVinit(size, 0.0) ;
         Drand_fillDvector(drand, size, p_entL[ieqn]) ;
      }
   }
}
/*--------------------------------------------------------------------*/
Drand_free(drand) ;
IVfree(list) ;

return(1) ; }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
