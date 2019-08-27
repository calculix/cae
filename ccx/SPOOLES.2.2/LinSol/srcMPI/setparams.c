/*  setparams.c  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to set the matrix parameters

   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- neqns <= 0 
    -3 -- type is invalid
    -4 -- symmetryflag is invalid
    -5 -- matrix is hermitian but type is real

   created -- 98sep25, cca
   -------------------------------------------
*/
int
BridgeMPI_setMatrixParams (
   BridgeMPI   *bridge,
   int         neqns, 
   int         type, 
   int         symmetryflag 
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMatrixParams()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( neqns <= 0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMatrixParams()"
           "\n neqns = %d\n", neqns) ;
   return(-2) ;
}
if ( type != SPOOLES_REAL && type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMatrixParams()"
           "\n type = %d\n", type) ;
   return(-3) ;
}
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC :
   break ;
case SPOOLES_HERMITIAN :
   if ( type != SPOOLES_COMPLEX ) {
      fprintf(stderr, "\n\n error in BridgeMPI_setMatrixParams()"
              "\n type = %d\n", type) ;
      return(-5) ;
   }
   break ;
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   if ( type != SPOOLES_COMPLEX ) {
      fprintf(stderr, "\n\n error in BridgeMPI_setMatrixParams()"
              "\n symmetryflag = %d\n", symmetryflag) ;
      return(-4) ;
   }
   break ;
}
bridge->neqns        = neqns        ;
bridge->type         = type         ;
bridge->symmetryflag = symmetryflag ;

return(1) ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to set the MPI parameters

   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- nproc <= 0 
    -3 -- myid < 0 or myid >= nproc

   created -- 98sep25, cca
   -------------------------------------------
*/
int
BridgeMPI_setMPIparams (
   BridgeMPI   *bridge,
   int         nproc, 
   int         myid, 
   MPI_Comm    comm 
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMPIparams()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( nproc <= 0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMPIparams()"
           "\n nproc = %d\n", nproc) ;
   return(-2) ;
}
if ( myid < 0 || myid >= nproc ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMPIparams()"
           "\n myid = %d, nproc = %d\n", myid, nproc) ;
   return(-3) ;
}
bridge->nproc = nproc ;
bridge->myid  = myid  ;
bridge->comm  = comm  ;

return(1) ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to set the ordering parameters

   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- maxdomainsize <= 0 
    -3 -- maxsize <= 0 
    -4 -- compressCutoff > 1.0

   created -- 98sep25, cca
   -------------------------------------------
*/
int
BridgeMPI_setOrderingParams (
   BridgeMPI   *bridge,
   int         maxdomainsize, 
   int         maxnzeros, 
   int         maxsize,
   int         seed,
   int         compressCutoff
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setOrderingParams()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( maxdomainsize <= 0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setOrderingParams()"
           "\n maxdomainsize = %d\n", maxdomainsize) ;
   return(-2) ;
}
if ( maxsize <= 0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setOrderingParams()"
           "\n maxsize = %d\n", maxsize) ;
   return(-2) ;
}
if ( compressCutoff > 1.0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setOrderingParams()"
           "\n compressCutoff = %d\n", compressCutoff) ;
   return(-2) ;
}
bridge->maxdomainsize  = maxdomainsize  ;
bridge->maxnzeros      = maxnzeros      ;
bridge->maxsize        = maxsize        ;
bridge->seed           = seed           ;
bridge->compressCutoff = compressCutoff ;

return(1) ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to set the message info

   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- msglvl > 0 and msgFile is NULL

   created -- 98sep25, cca
   -------------------------------------
*/
int
BridgeMPI_setMessageInfo (
   BridgeMPI   *bridge,
   int         msglvl,
   FILE        *msgFile 
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMessageInfo()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setMessageInfo()"
           "\n msglvl is > 0 and msgFile is NULL\n") ;
   return(-2) ;
}
bridge->msglvl  = msglvl  ;
bridge->msgFile = msgFile ;

return(1) ; }
   
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to set the factorization parameters

   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- sparsityflag is invalid
    -3 -- pivotingflag is invalid
    -4 -- tau < 2.0
    -5 -- droptol < 0.0
    -6 -- lookahead < 0

   created -- 98sep25, cca
   ----------------------------------------------
*/
int
BridgeMPI_setFactorParams (
   BridgeMPI        *bridge,
   int              sparsityflag, 
   int              pivotingflag, 
   double           tau, 
   double           droptol,
   int              lookahead, 
   PatchAndGoInfo   *patchinfo
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setFactorParams()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if (  sparsityflag != FRONTMTX_DENSE_FRONTS 
   && sparsityflag != FRONTMTX_SPARSE_FRONTS ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setFactorParams()"
           "\n sparsityflag = %d\n", sparsityflag) ;
   return(-2) ;
}
if (  pivotingflag != SPOOLES_PIVOTING 
   && pivotingflag != SPOOLES_NO_PIVOTING ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setFactorParams()"
           "\n pivotingflag = %d\n", pivotingflag) ;
   return(-3) ;
}
if (  tau < 2.0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setFactorParams()"
           "\n invalid value %f for tau", tau) ;
   return(-4) ;
}
if (  droptol < 0.0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setFactorParams()"
           "\n invalid value %f for droptol", droptol) ;
   return(-5) ;
}
if (  lookahead < 0 ) {
   fprintf(stderr, "\n\n error in BridgeMPI_setFactorParams()"
           "\n invalid value %d for lookahead", lookahead) ;
   return(-6) ;
}
bridge->sparsityflag = sparsityflag ;
bridge->pivotingflag = pivotingflag ;
bridge->tau          = tau          ;
bridge->droptol      = droptol      ;
bridge->lookahead    = lookahead    ;
bridge->patchinfo    = patchinfo    ;

return(1) ; }
   
/*--------------------------------------------------------------------*/
