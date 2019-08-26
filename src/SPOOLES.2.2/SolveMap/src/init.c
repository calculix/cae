/*  init.c  */

#include "../SolveMap.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- set the scalars and allocate the vectors for the object
 
   created -- 98mar19, cca
   ------------------------------------------------------------------
*/
void
SolveMap_init (
   SolveMap   *solvemap,
   int        symmetryflag,
   int        nfront,
   int        nproc,
   int        nblockUpper,
   int        nblockLower
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || symmetryflag < 0 || nfront <= 0 
   || nproc < 0 || nblockUpper < 0 || nblockLower < 0 ) {
   fprintf(stderr, "\n fatal error in SolveMap_init(%p,%d,%d,%d,%d,%d)"
           "\n bad input\n", solvemap, symmetryflag, nfront, 
           nproc, nblockUpper, nblockLower) ;
   exit(-1) ;
}
/*
   ----------------
   clear the object
   ----------------
*/
SolveMap_clearData(solvemap) ;
/*
   ---------------
   set the scalars
   ---------------
*/
solvemap->symmetryflag = symmetryflag ;
solvemap->nfront       = nfront       ;
solvemap->nproc        = nproc        ;
solvemap->nblockUpper  = nblockUpper  ;
solvemap->nblockLower  = nblockLower  ;
/*
   ------------------------
   allocate the data arrays
   ------------------------
*/
solvemap->owners      = IVinit(nfront, -1) ;
solvemap->rowidsUpper = IVinit(nblockUpper, -1) ;
solvemap->colidsUpper = IVinit(nblockUpper, -1) ;
solvemap->mapUpper    = IVinit(nblockUpper, -1) ;
if (  symmetryflag == SPOOLES_NONSYMMETRIC && nblockLower > 0 ) {
   solvemap->rowidsLower = IVinit(nblockLower, -1) ;
   solvemap->colidsLower = IVinit(nblockLower, -1) ;
   solvemap->mapLower    = IVinit(nblockLower, -1) ;
}
return ; }

/*--------------------------------------------------------------------*/
