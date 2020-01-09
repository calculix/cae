/*  instance.c  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- load *pobj with the address of the
              old-to-new permutation IV object

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL

   created -- 98sep18, cca
   ---------------------------------------------
*/
int
BridgeMPI_oldToNewIV (
   BridgeMPI   *bridge,
   IV          **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_oldToNewIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_oldToNewIV"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->oldToNewIV ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- load *pobj with the address of the
              new-to-old permutation IV object

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL

   created -- 98sep18, cca
   ---------------------------------------------
*/
int
BridgeMPI_newToOldIV (
   BridgeMPI   *bridge,
   IV          **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_newToOldIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_newToOldIV"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->newToOldIV ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the front ETree object

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL

   created -- 98sep18, cca
   --------------------------------------
*/
int
BridgeMPI_frontETree (
   BridgeMPI   *bridge,
   ETree       **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_frontETree"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_frontETree"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->frontETree ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- load *pobj with the address of the
              symbolic factorization IVL object

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL

   created -- 98sep18, cca
   ---------------------------------------------
*/
int
BridgeMPI_symbfacIVL (
   BridgeMPI   *bridge,
   IVL         **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_symbfacIVL"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_symbfacIVL"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->symbfacIVL ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- load *pobj with the address of 
              the submatrix manager object

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL

   created -- 98sep18, cca
   -----------------------------------------
*/
int
BridgeMPI_mtxmanager (
   BridgeMPI       *bridge,
   SubMtxManager   **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_mtxmanager"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_mtxmanager"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->mtxmanager ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the front matrix object

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL

   created -- 98sep18, cca
   --------------------------------------
*/
int
BridgeMPI_frontmtx (
   BridgeMPI   *bridge,
   FrontMtx    **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_frontmtx"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_frontmtx"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->frontmtx ;

return(1) ; }
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the owners IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep25, cca
   --------------------------------------
*/
int
BridgeMPI_ownersIV (
   BridgeMPI   *bridge,
   IV          **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_ownersIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_ownersIV"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->ownersIV ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- load *pobj with the address of
              the solve map SolveMap object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep25, cca
   -----------------------------------------
*/
int
BridgeMPI_solvemap (
   BridgeMPI   *bridge,
   SolveMap    **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_solvemap"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_solvemap"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->solvemap ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the vtxmap IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   --------------------------------------
*/
int
BridgeMPI_vtxmapIV (
   BridgeMPI   *bridge,
   IV          **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_vtxmapIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_vtxmapIV"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->vtxmapIV ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the rowmap IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   --------------------------------------
*/
int
BridgeMPI_rowmapIV (
   BridgeMPI   *bridge,
   IV          **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_rowmapIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_rowmapIV"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->rowmapIV ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- load *pobj with the address
              of the ownedColumns IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   ----------------------------------------
*/
int
BridgeMPI_ownedColumns (
   BridgeMPI   *bridge,
   IV          **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_ownedColumns"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_ownedColumns"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->ownedColumnsIV ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- load *pobj with the address
              of the Xloc DenseMtx object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   ----------------------------------------
*/
int
BridgeMPI_Xloc (
   BridgeMPI   *bridge,
   DenseMtx    **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_Xloc"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_Xloc"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->Xloc ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- load *pobj with the address
              of the Yloc DenseMtx object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   ----------------------------------------
*/
int
BridgeMPI_Yloc (
   BridgeMPI   *bridge,
   DenseMtx    **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_Yloc"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_Yloc"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->Yloc ;
 
return(1) ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- load *pnproc with the number of processors

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pnproc is NULL

   created -- 98sep18, cca
   -----------------------------------------------------
*/
int
BridgeMPI_nproc (
   BridgeMPI   *bridge,
   int         *pnproc
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_nproc()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pnproc == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_nproc()"
           "\n pnproc is NULL\n") ;
   return(-2) ;
}
*pnproc = bridge->nproc ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- load *pmyid with the id of this processor

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pmyid is NULL

   created -- 98sep18, cca
   ----------------------------------------------------
*/
int
BridgeMPI_myid (
   BridgeMPI   *bridge,
   int         *pmyid
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_myid()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pmyid == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_myid()"
           "\n pmyid is NULL\n") ;
   return(-2) ;
}
*pmyid = bridge->myid ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- load *plookahead with the lookahead value 
              for the factorization

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- plookahead is NULL

   created -- 98sep18, cca
   ----------------------------------------------------
*/
int
BridgeMPI_lookahead (
   BridgeMPI   *bridge,
   int         *plookahead
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_lookahead()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( plookahead == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_lookahead()"
           "\n plookahead is NULL\n") ;
   return(-2) ;
}
*plookahead = bridge->lookahead ;

return(1) ; }

/*--------------------------------------------------------------------*/
