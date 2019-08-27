/*  solve.c  */

#include "../BridgeMT.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to solve the linear system

   return value ---
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- X is NULL
     -3 -- Y is NULL
     -4 -- frontmtx is NULL
     -5 -- mtxmanager is NULL
     -6 -- oldToNewIV not available
     -7 -- newToOldIV not available

   created -- 98sep18, cca
   -------------------------------------
*/
int
BridgeMT_solve (
   BridgeMT   *bridge,
   int        permuteflag,
   DenseMtx   *X,
   DenseMtx   *Y
) {
double          cputotal, nops, t0, t1, t2 ;
double          cpus[6] ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
int             msglvl ;
SubMtxManager   *mtxmanager ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solve"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( X == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solve"
           "\n X is NULL\n") ;
   return(-2) ;
}
if ( Y == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solve"
           "\n Y is NULL\n") ;
   return(-3) ;
}
if ( (frontmtx = bridge->frontmtx) == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solve"
           "\n frontmtx is NULL\n") ;
   return(-4) ;
}
if ( (mtxmanager = bridge->mtxmanager) == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_solve"
           "\n mtxmanager is NULL\n") ;
   return(-5) ;
}
msglvl  = bridge->msglvl  ;
msgFile = bridge->msgFile ;
/*
   --------------------------
   optionally permute the rhs
   --------------------------
*/
if ( permuteflag == 1 ) {
   int   rc ;
   IV    *oldToNewIV ;

   MARKTIME(t1) ;
   rc = BridgeMT_oldToNewIV(bridge, &oldToNewIV) ;
   if (rc != 1) {
     fprintf(stderr, "\n error in BridgeMT_solve()"
             "\n rc = %d from BridgeMT_oldToNewIV()\n", rc) ;
     return(-6) ;
   }
   DenseMtx_permuteRows(Y, oldToNewIV) ;
   MARKTIME(t2) ;
   bridge->cpus[12] += t2 - t1 ;
}
/*
   ----------------
   solve the system
   ----------------
*/
nops = ETree_nFactorEntries(bridge->frontETree, bridge->symmetryflag) ;
nops *= 2 * X->ncol ;
if ( bridge->type == SPOOLES_COMPLEX ) {
   nops *= 4 ;
}
MARKTIME(t1) ;
DVzero(6, cpus) ;
FrontMtx_MT_solve(frontmtx, X, Y, mtxmanager, 
                  bridge->solvemap, cpus, msglvl, msgFile) ;
MARKTIME(t2) ;
bridge->cpus[13] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n CPU %8.3f : solve the system, %.3f mflops",
           t2 - t1, 1.e-6*nops/(t2 - t1)) ;
}
cputotal = t2 - t1 ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    set up solves               %8.3f %6.2f"
   "\n    load rhs and store solution %8.3f %6.2f"
   "\n    forward solve               %8.3f %6.2f"
   "\n    diagonal solve              %8.3f %6.2f"
   "\n    backward solve              %8.3f %6.2f"
   "\n    total time                  %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal, cputotal) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(X, msgFile) ;
   fflush(stdout) ;
}
/*
   -------------------------------
   optionally permute the solution
   -------------------------------
*/
if ( permuteflag == 1 ) {
   int   rc ;
   IV    *newToOldIV ;

   MARKTIME(t1) ;
   rc = BridgeMT_newToOldIV(bridge, &newToOldIV) ;
   if (rc != 1) {
     fprintf(stderr, "\n error in BridgeMT_solve()"
             "\n rc = %d from BridgeMT_newToOldIV()\n", rc) ;
     return(-7) ;
   }
   DenseMtx_permuteRows(X, newToOldIV) ;
   MARKTIME(t2) ;
   bridge->cpus[14] += t2 - t1 ;
}
MARKTIME(t2) ;
bridge->cpus[15] += t2 - t0 ;

return(1) ; }

/*--------------------------------------------------------------------*/
