/*  instance.c  */

#include "../BridgeMT.h"

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
BridgeMT_oldToNewIV (
   BridgeMT   *bridge,
   IV         **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_oldToNewIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_oldToNewIV"
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
BridgeMT_newToOldIV (
   BridgeMT   *bridge,
   IV         **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_newToOldIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_newToOldIV"
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
BridgeMT_frontETree (
   BridgeMT   *bridge,
   ETree      **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_frontETree"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_frontETree"
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
BridgeMT_symbfacIVL (
   BridgeMT   *bridge,
   IVL        **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_symbfacIVL"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_symbfacIVL"
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
BridgeMT_mtxmanager (
   BridgeMT        *bridge,
   SubMtxManager   **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_mtxmanager"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_mtxmanager"
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
BridgeMT_frontmtx (
   BridgeMT   *bridge,
   FrontMtx   **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_frontmtx"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_frontmtx"
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

   created -- 98sep24, cca
   --------------------------------------
*/
int
BridgeMT_ownersIV (
   BridgeMT   *bridge,
   IV         **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_ownersIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_ownersIV"
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

   created -- 98sep24, cca
   -----------------------------------------
*/
int
BridgeMT_solvemap (
   BridgeMT   *bridge,
   SolveMap   **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solvemap"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solvemap"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->solvemap ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- load *pnthread with the number of threads

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pnthread is NULL

   created -- 98sep24, cca
   ----------------------------------------------------
*/
int
BridgeMT_nthread (
   BridgeMT   *bridge,
   int        *pnthread
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_nthread"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pnthread == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_nthread"
           "\n pnthread is NULL\n") ;
   return(-2) ;
}
*pnthread = bridge->nthread ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- load *plookahead with the lookahead parameter

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- plookahead is NULL

   created -- 98sep24, cca
   --------------------------------------------------------
*/
int
BridgeMT_lookahead (
   BridgeMT   *bridge,
   int        *plookahead
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_lookahead"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( plookahead == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_lookahead"
           "\n plookahead is NULL\n") ;
   return(-2) ;
}
*plookahead = bridge->lookahead ;

return(1) ; }

/*--------------------------------------------------------------------*/
