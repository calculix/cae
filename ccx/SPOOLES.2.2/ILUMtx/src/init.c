/*  init.c  */

#include "../ILUMtx.h"

#define MYDEBUG 0

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*
   ---------------------------------------------
   purpose -- initialize the ILUMtx object

   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtx
     -4  -- bad symmetryflag for mtx
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- matrix is symmetric or hermitian
            and storage modes are not compatible

   created -- 98oct03, cca
   ---------------------------------------------
*/
int
ILUMtx_init (
   ILUMtx   *mtx,
   int      neqns,
   int      type,
   int      symmetryflag,
   int      LstorageMode,
   int      UstorageMode
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n error in ILUM_init(), mtx = NULL\n") ;
   return(-1) ;
}
if ( neqns <= 0 ) {
   fprintf(stderr, "\n error in ILUM_init()"
           "\n neqns = %d\n", neqns) ;
   return(-2) ;
}
if ( type != SPOOLES_REAL && type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n error in ILUM_init()"
           "\n type = %d\n", type) ;
   return(-3) ;
}
if (   symmetryflag != SPOOLES_SYMMETRIC 
    && symmetryflag != SPOOLES_HERMITIAN 
    && symmetryflag != SPOOLES_NONSYMMETRIC ) {
   fprintf(stderr, "\n error in ILUMinit()"
           "\n symmetry = %d\n", symmetryflag) ;
   return(-4) ;
}
if (   LstorageMode != SPOOLES_BY_ROWS
    && LstorageMode != SPOOLES_BY_COLUMNS ) {
   fprintf(stderr, "\n error in ILUM_init()"
           "\n LstorageMode = %d\n", LstorageMode) ;
   return(-5) ;
}
if (   UstorageMode != SPOOLES_BY_ROWS
    && UstorageMode != SPOOLES_BY_COLUMNS ) {
   fprintf(stderr, "\n error in ILUM_init()"
           "\n UstorageMode = %d\n", UstorageMode) ;
   return(-6) ;
}
if (   (   symmetryflag == SPOOLES_SYMMETRIC 
        || symmetryflag == SPOOLES_HERMITIAN) 
    && (LstorageMode == UstorageMode) ) {
   fprintf(stderr, "\n error in ILUM_init()"
           "\n symmetryflag %d, LstorageMode %d, UstorageMode %d",
           symmetryflag, LstorageMode, UstorageMode) ;
   return(-7) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------
   clear the data
   --------------
*/
ILUMtx_clearData(mtx) ;
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
mtx->neqns        = neqns        ;
mtx->type         = type         ;
mtx->symmetryflag = symmetryflag ;
mtx->LstorageMode = LstorageMode ;
mtx->UstorageMode = UstorageMode ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n mtx->neqns = %d"
        "\n mtx->type = %d"
        "\n mtx->symmetryflag = %d"
        "\n mtx->LstorageMode = %d"
        "\n mtx->UstorageMode = %d",
        mtx->neqns, mtx->type, mtx->symmetryflag,
        mtx->LstorageMode, mtx->UstorageMode) ;
fflush(stdout) ;
#endif
/*
   --------------------
   allocate the vectors
   --------------------
*/
mtx->sizesU = IVinit(neqns, 0) ;
mtx->p_indU = PIVinit(neqns) ;
mtx->p_entU = PDVinit(neqns) ;
if ( type == SPOOLES_REAL ) {
   mtx->entD = DVinit(neqns, 0.0) ;
} else {
   mtx->entD = DVinit(2*neqns, 0.0) ;
}
if ( symmetryflag == SPOOLES_NONSYMMETRIC ) {
   mtx->sizesL = IVinit(neqns, 0) ;
   mtx->p_indL = PIVinit(neqns) ;
   mtx->p_entL = PDVinit(neqns) ;
} else {
   mtx->sizesL = NULL ;
   mtx->p_indL = NULL ;
   mtx->p_entL = NULL ;
}
/*--------------------------------------------------------------------*/
return(1) ; }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
