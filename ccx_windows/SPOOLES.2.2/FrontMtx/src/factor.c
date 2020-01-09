/*  factor.c  */

#include "../FrontMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   compute an (U^T + I)D(I + U) or (L + I)D(I + L) factorization of A.
   this is a wrapper method around FrontMtx_factorPencil().

   input --

      frontmtx -- pointer to the FrontMtx object that will hold
                  the factorization
      pencil   -- pointer to the Pencil object that holds A + sigma*B
      tau      -- upper bound on entries in L and U,
                  used only when pivoting enabled
      droptol  -- lower bound on entries in L and U,
                  used only when sparsity enabled
      perror   -- error flag, on return
         *perror >= 0 --> factorization failed at front *perror
         *perror <  0 --> factorization succeeded
      cpus[]   -- timing array
         cpus[0] -- initialize fronts
         cpus[1] -- load original entries
         cpus[2] -- get updates from descendents
         cpus[3] -- assembled postponed data
         cpus[4] -- factor the front
         cpus[5] -- extract postponed data
         cpus[6] -- store factor entries
         cpus[7] -- miscellaneous time
         cpus[8] -- total time
      stats[] -- statistics array
         stats[0] -- # of pivots
         stats[1] -- # of pivot tests
         stats[2] -- # of delayed rows and columns
         stats[3] -- # of entries in D
         stats[4] -- # of entries in L
         stats[5] -- # of entries in U
      msglvl   -- message level
      msgFile  -- message file

   created  -- 98mar27, cca
   modified -- 98mar27, cca
      perror added to argument list
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_factorInpMtx (
   FrontMtx     *frontmtx,
   InpMtx       *inpmtx,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) {
double   zero[2] = {0.0, 0.0} ;
Chv      *rootchv ;
Pencil   pencil ;

Pencil_setDefaultFields(&pencil) ;
Pencil_init(&pencil, frontmtx->type, frontmtx->symmetryflag,
            inpmtx, zero, NULL) ;
rootchv = FrontMtx_factorPencil(frontmtx, &pencil, tau, droptol, 
                     chvmanager, perror, cpus, stats, msglvl, msgFile) ;

return(rootchv) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   compute an (U^T + I)D(I + U) or (L + I)D(I + L) 
   factorization of A + sigma*B. 

   input --

      frontmtx -- pointer to the FrontMtx object that will hold
                  the factorization
      pencil   -- pointer to the Pencil object that holds A + sigma*B
      tau      -- upper bound on entries in L and U,
                  used only when pivoting enabled
      droptol  -- lower bound on entries in L and U,
                  used only when sparsity enabled
      perror   -- error flag, on return
         *perror >= 0 --> factorization failed at front *perror
         *perror <  0 --> factorization succeeded
      cpus[]   -- timing array
         cpus[0] -- initialize fronts
         cpus[1] -- load original entries
         cpus[2] -- get updates from descendents
         cpus[3] -- assembled postponed data
         cpus[4] -- factor the front
         cpus[5] -- extract postponed data
         cpus[6] -- store factor entries
         cpus[7] -- miscellaneous time
         cpus[8] -- total time
      stats[] -- statistics array
         stats[0] -- # of pivots
         stats[1] -- # of pivot tests
         stats[2] -- # of delayed rows and columns
         stats[3] -- # of entries in D
         stats[4] -- # of entries in L
         stats[5] -- # of entries in U
      msglvl   -- message level
      msgFile  -- message file

   created  -- 98mar27, cca
   modified -- 98mar27, cca
      perror added to argument list
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_factorPencil (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) {
char         *status ;
ChvList      *postList ;
Chv          *rootchv ;
Chv          **fronts ;
double       t0, t3 ;
DV           tmpDV ;
ETree        *frontETree ;
int          J, K, ndelayed, nfront, npivots, ntests ;
int          *par ;
IP           **heads ;
IV           pivotsizesIV ;
Tree         *tree ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( frontmtx == NULL || pencil == NULL || cpus == NULL || stats == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_factorPencil()"
           "\n frontmtx = %p, pencil = %p"
           "\n tau = %e, droptol = %e, cpus = %p"
           "\n msglvl = %d, msgFile = %p"
           "\n bad input\n",
           frontmtx, pencil, tau, droptol, cpus, msglvl, msgFile) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n INSIDE FrontMtx_factorPencil()") ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   extract pointers and dimensions
   -------------------------------
*/
frontETree = frontmtx->frontETree ;
nfront     = ETree_nfront(frontETree) ;
tree       = ETree_tree(frontETree) ;
par        = ETree_par(frontETree) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n got pointers and dimensions") ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   allocate and initialize working storage
   ---------------------------------------
*/
heads  = FrontMtx_factorSetup(frontmtx, NULL, 0, msglvl, msgFile) ;
status = CVinit(nfront, 'W') ;
ALLOCATE(fronts, Chv *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   fronts[J] = NULL ;
}
DV_setDefaultFields(&tmpDV) ;
IV_setDefaultFields(&pivotsizesIV) ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   postList = ChvList_new() ;
   ChvList_init(postList, nfront+1, NULL, 0, NULL) ;
} else {
   postList = NULL ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n allocated working storage") ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------
   loop over the tree in a post-order traversal
   --------------------------------------------
*/
*perror = -1 ;
npivots = ndelayed = ntests = 0 ;
rootchv = NULL ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   K = par[J] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n ##### working on front %d, parent %d", 
              J, K) ;
      fflush(msgFile) ;
   }
   FrontMtx_factorVisit(frontmtx, pencil, J, 0, NULL, fronts, 0,
                         tau, droptol, status, heads, &pivotsizesIV, 
                         &tmpDV, par, NULL, postList, 
                         chvmanager, stats, cpus, msglvl, msgFile) ;
   if ( status[J] == 'E' ) {
/*
      -------------------------------
      error found, unable to continue
      -------------------------------
*/
      *perror = J ;
      break ;
   } else if ( status[J] != 'F' ) {
      fprintf(stderr, "\n fatal error, return %c from front %d",
              status[J], J) ;
      exit(-1) ;
   }
}
/*
   ---------------------------------
   get a pointer to the root chevron
   ---------------------------------
*/
if ( postList == NULL ) {
   rootchv = NULL ;
} else {
   rootchv = ChvList_getList(postList, nfront) ;
}
/*
   ------------------
   set the statistics
   ------------------
*/
stats[3] = frontmtx->nentD ;
stats[4] = frontmtx->nentL ;
stats[5] = frontmtx->nentU ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IP_free(heads[nfront+1]) ;
FREE(heads) ;
DV_clearData(&tmpDV) ;
IV_clearData(&pivotsizesIV) ;
CVfree(status) ;
FREE(fronts) ;
if ( postList != NULL ) {
   ChvList_free(postList) ;
}
/*
   --------------------------------
   set final and miscellaneous time
   --------------------------------
*/
MARKTIME(t3) ;
cpus[8] = t3 - t0 ;
cpus[7] = cpus[8] - cpus[0] - cpus[1] - cpus[2]
        - cpus[3] - cpus[4] - cpus[5] - cpus[6] ;

return(rootchv) ; }

/*--------------------------------------------------------------------*/
