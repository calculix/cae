/*  info.c  */

#include "../BridgeMT.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose --  generate and return some statistics 
               about the factor and solve

   type -- type of entries
     SPOOLES_REAL or SPOOLES_COMPLEX
   symmetryflag -- symmetry type
     SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   on return ---
      *pnfront     -- # of fronts
      *pnfactorind -- # of factor indices
      *pnfactorent -- # of factor entries
      *pnsolveops  -- # of solve operations 
      *pnfactorops -- # of factor operations 

   return values --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- type is bad, must be SPOOLES_REAL or SPOOLES_COMPLEX
     -3 -- symmetryflag is bad, must be SPOOLES_SYMMETRIC,
           SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
     -4 -- type and symmetryflag mismatch
     -5 -- front tree is not present
     -6 -- pnfront is NULL
     -7 -- pnfactorind is NULL
     -8 -- pnfactorent is NULL
     -9 -- pnsolveops is NULL
    -10 -- pnfactorops is NULL
 
   created -- 98oct01, cca
   --------------------------------------------------------------
*/
int
BridgeMT_factorStats (
   BridgeMT   *bridge,
   int        type,
   int        symmetryflag,
   int        *pnfront,
   int        *pnfactorind,
   int        *pnfactorent,
   int        *pnsolveops,
   double     *pnfactorops
) {
ETree   *etree ;
int     nentD, nentU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n bad type %d\n", type) ;
   return(-3) ;
   break ;
}
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n bad symmetryflag %d\n", symmetryflag) ;
   return(-3) ;
   break ;
}
if ( type == SPOOLES_REAL && symmetryflag == SPOOLES_HERMITIAN ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n type %d, symmetryflag %d, mismatch\n", 
           type, symmetryflag) ;
   return(-4) ;
}
if ( (etree = bridge->frontETree) == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n front tree is not present\n") ;
   return(-5) ;
}
if ( pnfront == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n pnfront is NULL\n") ;
   return(-6) ;
}
if ( pnfactorind == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n pnfactorind is NULL\n") ;
   return(-7) ;
}
if ( pnfactorent == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n pnfactorent is NULL\n") ;
   return(-8) ;
}
if ( pnsolveops == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n pnsolveops is NULL\n") ;
   return(-9) ;
}
if ( pnfactorops == NULL ) {
   fprintf(stderr, "\n error in Bridge_factorStats()"
           "\n pnfactorops is NULL\n") ;
   return(-10) ;
}
*pnfront     = ETree_nfront(etree) ;
*pnfactorind = ETree_nFactorIndices(etree) ;
*pnfactorent = ETree_nFactorEntries(etree, symmetryflag) ;
*pnfactorops = ETree_nFactorOps(etree, type, symmetryflag) ;
nentD = etree->nvtx ;
nentU = *pnfactorent - nentD ;
switch ( type ) {
case SPOOLES_REAL :
   *pnsolveops = 4*nentU + nentD ;
   break ;
case SPOOLES_COMPLEX :
   *pnsolveops = 16*nentU + 8*nentD ;
   break ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
