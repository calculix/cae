/*  solve.c  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to solve the linear system
      MPI version
      if permuteflag is 1 then
         rhs is permuted into new ordering
         solution is permuted into old ordering

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
   --------------------------------------------
*/
int
BridgeMPI_solve (
   BridgeMPI  *bridge,
   int        permuteflag,
   DenseMtx   *X,
   DenseMtx   *Y
) {
DenseMtx        *Xloc, *Yloc ;
double          cputotal, t0, t1, t2 ;
double          cpus[6] ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
int             firsttag, msglvl, myid, nmycol, nrhs, nrow ;
int             *mycolind, *rowind ;
int             stats[4] ;
IV              *mapIV, *ownersIV ;
MPI_Comm        comm ;
SubMtxManager   *mtxmanager ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_solve"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( (frontmtx = bridge->frontmtx) == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_solve"
           "\n frontmtx is NULL\n") ;
   return(-4) ;
}
if ( (mtxmanager = bridge->mtxmanager) == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_solve"
           "\n mtxmanager is NULL\n") ;
   return(-5) ;
}
myid     = bridge->myid     ;
comm     = bridge->comm     ;
msglvl   = bridge->msglvl   ;
msgFile  = bridge->msgFile  ;
frontmtx = bridge->frontmtx ;
ownersIV = bridge->ownersIV ;
Xloc     = bridge->Xloc     ;
Yloc     = bridge->Yloc     ;
if ( myid != 0 ) {
   X = Y = NULL ;
} else {
   if ( X == NULL ) {
      fprintf(stderr, "\n error in BridgeMPI_solve"
              "\n myid 0, X is NULL\n") ;
      return(-2) ;
   }
   if ( Y == NULL ) {
      fprintf(stderr, "\n error in BridgeMPI_solve"
              "\n myid 0, Y is NULL\n") ;
      return(-3) ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside BridgeMPI_solve()") ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile , "\n\n Xloc") ;
   DenseMtx_writeForHumanEye(Xloc, msgFile) ;
   fprintf(msgFile , "\n\n Yloc") ;
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
if ( myid == 0 ) {
/*
   --------------------------
   optionally permute the rhs
   --------------------------
*/
   if ( permuteflag == 1 ) {
      int   rc ;
      IV    *oldToNewIV ;
   
      MARKTIME(t1) ;
      rc = BridgeMPI_oldToNewIV(bridge, &oldToNewIV) ;
      if (rc != 1) {
        fprintf(stderr, "\n error in BridgeMPI_solve()"
                "\n rc = %d from BridgeMPI_oldToNewIV()\n", rc) ;
        return(-6) ;
      }
      DenseMtx_permuteRows(Y, oldToNewIV) ;
      MARKTIME(t2) ;
      bridge->cpus[15] += t2 - t1 ;
      if ( msglvl > 2 ) {
         fprintf(msgFile , "\n\n permuted Y") ;
         DenseMtx_writeForHumanEye(Y, msgFile) ;
         fflush(msgFile) ;
      }
   }
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   distribute the right hand side matrix
   -------------------------------------
*/
MARKTIME(t1) ;
mapIV = bridge->rowmapIV ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n row map IV object") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
if ( myid == 0 ) {
   nrhs = Y->ncol ;
} else {
   nrhs = 0 ;
}
MPI_Bcast((void *) &nrhs, 1, MPI_INT, 0, comm) ;
firsttag = 0 ;
IVfill(4, stats, 0) ;
DenseMtx_MPI_splitFromGlobalByRows(Y, Yloc, mapIV, 0, stats, 
                                   msglvl, msgFile, firsttag, comm) ;
MARKTIME(t2) ;
bridge->cpus[16] += t2 - t1 ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n local matrix Y after the split") ;
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   initialize the local solution X object
   --------------------------------------
*/
MARKTIME(t1) ;
IV_sizeAndEntries(bridge->ownedColumnsIV, &nmycol, &mycolind) ;
DenseMtx_init(Xloc, bridge->type, -1, -1, nmycol, nrhs, 1, nmycol) ;
if ( nmycol > 0 ) {
   DenseMtx_rowIndices(Xloc, &nrow, &rowind) ;
   IVcopy(nmycol, rowind, mycolind) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n local matrix X") ;
      DenseMtx_writeForHumanEye(Xloc, msgFile) ;
      fflush(msgFile) ;
   }
}
MARKTIME(t2) ;
bridge->cpus[17] += t2 - t1 ;
/*--------------------------------------------------------------------*/
/*
   ----------------
   solve the system
   ----------------
*/
MARKTIME(t1) ;
DVzero(6, cpus) ;
FrontMtx_MPI_solve(frontmtx, Xloc, Yloc, mtxmanager, bridge->solvemap,
                   cpus, stats, msglvl, msgFile, firsttag, comm) ;
MARKTIME(t2) ;
bridge->cpus[18] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n CPU %8.3f : solve the system", t2 - t1) ;
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
   DenseMtx_writeForHumanEye(Xloc, msgFile) ;
   fflush(stdout) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   gather the solution on processor zero
   -------------------------------------
*/
MARKTIME(t1) ;
DenseMtx_MPI_mergeToGlobalByRows(X, Xloc, 0, stats, msglvl, msgFile,
                                 firsttag, comm) ;
MARKTIME(t2) ;
bridge->cpus[19] += t2 - t1 ;
if ( myid == 0 ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n global matrix X in new ordering") ;
      DenseMtx_writeForHumanEye(X, msgFile) ;
      fflush(msgFile) ;
   }
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------
   optionally permute the solution
   -------------------------------
*/
if ( myid == 0 ) {
   if ( permuteflag == 1 ) {
      int   rc ;
      IV    *newToOldIV ;

      rc = BridgeMPI_newToOldIV(bridge, &newToOldIV) ;
      if (rc != 1) {
        fprintf(stderr, "\n error in BridgeMPI_solve()"
                "\n rc = %d from BridgeMPI_newToOldIV()\n", rc) ;
        return(-7) ;
      }
      DenseMtx_permuteRows(X, newToOldIV) ;
   }
   MARKTIME(t2) ;
   bridge->cpus[20] += t2 - t1 ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n global matrix X in old ordering") ;
      DenseMtx_writeForHumanEye(X, msgFile) ;
      fflush(msgFile) ;
   }
}
MARKTIME(t2) ;
bridge->cpus[21] += t2 - t0 ;

return(1) ; }

/*--------------------------------------------------------------------*/
