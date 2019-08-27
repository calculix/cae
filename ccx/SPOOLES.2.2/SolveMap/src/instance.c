/*  instance.c  */

#include "../SolveMap.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to return the symmetry flag of the object
  
   created -- 98mar19, cca
   ----------------------------------------------------
*/
int 
SolveMap_symmetryflag (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_symmetryflag(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->symmetryflag) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to return the number of fronts
  
   created -- 98mar19, cca
   -----------------------------------------
*/
int 
SolveMap_nfront (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_nfront(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->nfront) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to return the number of processors
  
   created -- 98mar19, cca
   ---------------------------------------------
*/
int 
SolveMap_nproc (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_nproc(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->nproc) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to return the number of blocks in the upper adjacency
  
   created -- 98mar19, cca
   ----------------------------------------------------------------
*/
int 
SolveMap_nblockUpper (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_nblockUpper(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->nblockUpper) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to return the number of blocks in the lower adjacency
  
   created -- 98mar19, cca
   ----------------------------------------------------------------
*/
int 
SolveMap_nblockLower (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_nblockLower(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->nblockLower) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return a pointer to the owners vector
  
   created -- 98mar19, cca
   ----------------------------------------------------
*/
int *
SolveMap_owners (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_owners(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->owners) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to return a pointer to the row ids vector
              for the upper adjacency structure
  
   created -- 98mar19, cca
   -----------------------------------------------------
*/
int *
SolveMap_rowidsUpper (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_rowidsUpper(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->rowidsUpper) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to return a pointer to the column ids vector
              for the upper adjacency structure
  
   created -- 98mar19, cca
   --------------------------------------------------------
*/
int *
SolveMap_colidsUpper (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_colidsUpper(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->colidsUpper) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to return a pointer to the map vector
              for the upper adjacency structure
  
   created -- 98mar19, cca
   -------------------------------------------------
*/
int *
SolveMap_mapUpper (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_mapUpper(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->mapUpper) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to return a pointer to the row ids vector
              for the upper adjacency structure
  
   created -- 98mar19, cca
   -----------------------------------------------------
*/
int *
SolveMap_rowidsLower (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_rowidsLower(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->rowidsLower) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to return a pointer to the column ids vector
              for the upper adjacency structure
  
   created -- 98mar19, cca
   --------------------------------------------------------
*/
int *
SolveMap_colidsLower (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_colidsLower(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->colidsLower) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to return a pointer to the map vector
              for the upper adjacency structure
  
   created -- 98mar19, cca
   -------------------------------------------------
*/
int *
SolveMap_mapLower (
   SolveMap   *solvemap
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_mapLower(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
return(solvemap->mapLower) ; }

/*--------------------------------------------------------------------*/
