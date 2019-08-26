/*  solve.c  */

#include "../FrontMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- to solve a linear system 

   frontmtx -- FrontMtx object that holds the factor matrices
   solmtx   -- DenseMtx that holds the solution
   rhsmtx   -- DenseMtx that holds the right hand side matrix
      note: right hand side entries are copied from rhsmtx,
            and solution entries are copied into solmtx.
            when the row indices of rhsmtx and solmtx are
            identical, rhsmtx and solmtx can be the same object.
   cpus -- vector to hold cpu breakdown time
      cpus[0] -- set up solves
      cpus[1] -- fetch rhs and store solution
      cpus[2] -- forward solve
      cpus[3] -- diagonal solve
      cpus[4] -- backward solve
      cpus[5] -- total time
   mtxmanager -- object that manages working storage
   msglvl  -- message level
   msgFile -- message file

   created -- 98may04, cca
   ------------------------------------------------------------
*/
void
FrontMtx_solve (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   DenseMtx        *rhsmtx,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) {
char     *frontIsDone, *status ;
SubMtx   **p_mtx ;
double   t0, t1, t2 ;
int      J, nfront, nrhs ;
IP       **heads ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( frontmtx == NULL || solmtx == NULL || rhsmtx == NULL 
      || mtxmanager == NULL 
      || cpus == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_solve()"
           "\n bad input :"
           "\n    frontmtx = %p, solmtx = %p, rhsmtx = %p"
           "\n    mtxmanager = %p, cpus = %p"
           "\n    msglvl = %d, msgFile = %p\n",
           frontmtx, solmtx, rhsmtx, mtxmanager, 
           cpus, msglvl, msgFile) ;
   exit(-1) ;
}
nfront  = FrontMtx_nfront(frontmtx) ;
tree    = FrontMtx_frontTree(frontmtx) ;
nrhs    = rhsmtx->ncol ;
/*
   --------------------------------------------------
   set up the update head/links for the forward solve
   --------------------------------------------------
*/
MARKTIME(t1) ;
heads = FrontMtx_forwardSetup(frontmtx, msglvl, msgFile ) ;
frontIsDone = CVinit(nfront, 'N') ;
status      = CVinit(nfront, 'W') ;
MARKTIME(t2) ;
cpus[0] = t2 - t1 ;
/*
   ----------------------------------------------------
   load the right hand side into temporary SubMtx objects
   ----------------------------------------------------
*/
MARKTIME(t1) ;
p_mtx = FrontMtx_loadRightHandSide(frontmtx, rhsmtx, NULL, 0,
                                   mtxmanager, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[1] = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU : load rhs = %8.3f", t2 - t1) ;
}
/*
   ----------------------------
   forward solve: (L + I) y = b
   ----------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n ####### starting forward solve") ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n forward visiting front %d", J) ;
      fflush(msgFile) ;
   }
   FrontMtx_forwardVisit(frontmtx, J, nrhs, NULL, 0, mtxmanager, NULL, 
                          p_mtx, frontIsDone, heads, p_mtx, status, 
                          msglvl, msgFile) ;
}
IP_free(heads[nfront+1]) ;
FREE(heads) ;
MARKTIME(t2) ;
cpus[2] = t2 - t1 ;
/*
   -----------------
   solve D z_J = y_J
   -----------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n ####### starting diagonal solve") ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
CVfill(nfront, frontIsDone, 'N') ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n diagonal visiting front %d", J) ;
      fflush(msgFile) ;
   }
   FrontMtx_diagonalVisit(frontmtx, J, NULL, 0, 
                           p_mtx, frontIsDone, p_mtx, msglvl, msgFile) ;
   frontIsDone[J] = 'D' ;
}
MARKTIME(t2) ;
cpus[3] = t2 - t1 ;
/*
   ---------------------------------------------------
   set up the update head/links for the backward solve
   ---------------------------------------------------
*/
MARKTIME(t1) ;
heads = FrontMtx_backwardSetup(frontmtx, msglvl, msgFile) ;
CVfill(nfront, status, 'W') ;
CVfill(nfront, frontIsDone, 'N') ;
MARKTIME(t2) ;
cpus[0] += t2 - t1 ;
/*
   -----------------------------
   backward solve: (I + U) x = z
   -----------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n ####### starting backward solve") ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n backward visiting front %d", J) ;
      fflush(msgFile) ;
   }
   FrontMtx_backwardVisit(frontmtx, J, nrhs, NULL, 0, mtxmanager, NULL,
                           p_mtx, frontIsDone, heads, p_mtx, status, 
                           msglvl, msgFile) ;
}
MARKTIME(t2) ;
cpus[4] = t2 - t1 ;
/*
   ----------------------------
   store the solution in rhsmtx
   ----------------------------
*/
MARKTIME(t1) ;
FrontMtx_storeSolution(frontmtx, NULL, 0, mtxmanager,
                        p_mtx, solmtx, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[1] += t2 - t1 ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n CPU : store solution = %8.3f", t2 - t1) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IP_free(heads[nfront+1]) ;
FREE(heads) ;
FREE(p_mtx) ;
CVfree(frontIsDone) ;
CVfree(status) ;

MARKTIME(t2) ;
cpus[5] = t2 - t0 ;

return ; }

/*--------------------------------------------------------------------*/
