/*  instance.c  */

#include "../Bridge.h"

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
Bridge_oldToNewIV (
   Bridge   *bridge,
   IV       **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_oldToNewIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in Bridge_oldToNewIV"
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
Bridge_newToOldIV (
   Bridge   *bridge,
   IV       **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_newToOldIV"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in Bridge_newToOldIV"
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
Bridge_frontETree (
   Bridge   *bridge,
   ETree    **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_frontETree"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in Bridge_frontETree"
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
Bridge_symbfacIVL (
   Bridge   *bridge,
   IVL      **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_symbfacIVL"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in Bridge_symbfacIVL"
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
Bridge_mtxmanager (
   Bridge          *bridge,
   SubMtxManager   **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_mtxmanager"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in Bridge_mtxmanager"
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
Bridge_frontmtx (
   Bridge     *bridge,
   FrontMtx   **pobj
) {
/*
   ----------------
   check the output
   ----------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in Bridge_frontmtx"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( pobj == NULL ) {
   fprintf(stderr, "\n error in Bridge_frontmtx"
           "\n pobj is NULL\n") ;
   return(-2) ;
}
*pobj = bridge->frontmtx ;

return(1) ; }

/*--------------------------------------------------------------------*/
