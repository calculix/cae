/*  QRfactor.c  */

#include "../FrontMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to compute the (U+I)D(I+U) factorization of A^TA,
              where A = QR, R = (D^{1/2} + D^{1/2}U)

   cpus[0] -- setup time
   cpus[1] -- initialize and load staircase matrix
   cpus[2] -- factor the matrix
   cpus[3] -- scale and store factor entries
   cpus[4] -- store update entries
   cpus[5] -- miscellaneous time
   cpus[6] -- total time

   created -- 98may28, cca
   ------------------------------------------------------------
*/
void
FrontMtx_QR_factor (
   FrontMtx     *frontmtx,
   InpMtx       *mtxA,
   ChvManager   *chvmanager,
   double       cpus[],
   double       *pfacops,
   int          msglvl,
   FILE         *msgFile
) {
char      *status ;
ChvList   *updlist ;
double    t0, t1, t2, t3 ;
DV        workDV ;
int       J, neqns, nfront ;
int       *colmap, *firstnz ;
IVL       *rowsIVL ;
Tree      *tree    ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if (  frontmtx == NULL || mtxA == NULL || chvmanager == NULL 
   || cpus == NULL || pfacops == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_QR_factor()"
           "\n bad input\n") ;
   exit(-1) ;
}
neqns  = frontmtx->neqns  ;
nfront = frontmtx->nfront ;
tree   = frontmtx->tree   ;
/*
   ----------------------------------------------------------------
   create the update list object
   create the rowsIVL object, where
      list(J) = list of rows that are assembled in front J
   firstnz[irowA] = first column with nonzero element in A(irowA,*)
   colmap[neqns]  = work vector for mapping
   status[neqns]  = status vector for fronts
   ----------------------------------------------------------------
*/
MARKTIME(t1) ;
updlist = ChvList_new() ;
ChvList_init(updlist, nfront+1, NULL, NO_LOCK, NULL) ;
FrontMtx_QR_setup(frontmtx, mtxA, &rowsIVL, &firstnz, msglvl, msgFile) ;
colmap  = IVinit(neqns, -1) ;
status = CVinit(nfront, 'W') ;
DV_setDefaultFields(&workDV) ;
MARKTIME(t2) ;
cpus[0] += t2 - t1 ;
/*
   --------------------------------------------
   loop over the tree in a post-order traversal
   --------------------------------------------
*/
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   FrontMtx_QR_factorVisit(frontmtx, J, mtxA, rowsIVL, firstnz, updlist,
                           chvmanager, status, colmap, &workDV, cpus, 
                           pfacops, msglvl, msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
CVfree(status) ;
DV_clearData(&workDV) ;
ChvList_free(updlist) ;
IVL_free(rowsIVL) ;
IVfree(colmap) ;
IVfree(firstnz) ;

MARKTIME(t3) ;
cpus[6] = t3 - t0 ;
cpus[5] = t3 - t0 - cpus[0] - cpus[1] 
          - cpus[2] - cpus[3] - cpus[4] - cpus[5] ;

return ; }
   
/*--------------------------------------------------------------------*/
