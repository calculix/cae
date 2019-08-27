/*  JimSolveMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_JimSolve = 0 ;
static double time_JimSolve = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to solve a linear system
     (A - sigma*B) sol[] = rhs[]

   data    -- pointer to bridge data object
   *pnrows -- # of rows in x[] and y[]
   *pncols -- # of columns in x[] and y[]
   rhs[]   -- vector that holds right hand sides
      NOTE: the rhs[] vector is global, not a portion
   sol[]   -- vector to hold solutions
      NOTE: the sol[] vector is global, not a portion

   note: rhs[] and sol[] can be the same array.
  
   on return, *perror holds an error code.

   created -- 98aug28, cca & jcp
   --------------------------------------------------
*/
void 
JimSolveMPI ( 
   int       *pnrows, 
   int       *pncols, 
   double    rhs[], 
   double    sol[],
   void      *data, 
   int       *perror 
) {
BridgeMPI   *bridge = (BridgeMPI *) data ;
DenseMtx    *mtx, *newmtx ;
int         irow, jj, jcol, kk, myid, ncols = *pncols, 
            neqns, nowned, tag = 0 ;
int         *vtxmap ;
int         stats[4] ;
IV          *mapIV ;
#if MYDEBUG > 0
double   t1, t2 ;
count_JimSolve++ ;
MARKTIME(t1) ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) JimSolve() start", count_JimSolve) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) JimSolve() start", count_JimSolve) ;
fflush(bridge->msgFile) ;
#endif
MPI_Barrier(bridge->comm) ;
/*
   ---------------------------------------------
   slide the owned rows of rhs down in the array
   ---------------------------------------------
*/
vtxmap  = IV_entries(bridge->vtxmapIV) ;
neqns   = bridge->neqns ;
myid    = bridge->myid  ;
nowned  = IV_size(bridge->myownedIV) ;
for ( jcol = jj = kk = 0 ; jcol < ncols ; jcol++ ) {
   for ( irow = 0 ; irow < neqns ; irow++, jj++ ) {
      if ( vtxmap[irow] == myid ) {
         sol[kk++] = rhs[jj] ;
      }
   }
}
if ( kk != nowned * ncols ) {
   fprintf(stderr, "\n proc %d : kk %d, nowned %d, ncols %d",
           myid, kk, nowned, ncols) ;
   exit(-1) ;
}
/*
   ----------------------------------------
   call the method that assumes local input
   ----------------------------------------
*/
if ( bridge->msglvl > 1 ) {
   fprintf(bridge->msgFile, "\n calling SolveMPI()") ;
   fflush(bridge->msgFile) ;
}
SolveMPI(&nowned, pncols, sol, sol, data, perror) ;
if ( bridge->msglvl > 1 ) {
   fprintf(bridge->msgFile, "\n return from SolveMPI()") ;
   fflush(bridge->msgFile) ;
}
/*
   ------------------------------------------
   gather all the entries onto processor zero
   ------------------------------------------
*/
mtx = DenseMtx_new() ;
DenseMtx_init(mtx, SPOOLES_REAL, 0, 0, nowned, ncols, 1, nowned) ;
DVcopy (nowned*ncols, DenseMtx_entries(mtx), sol) ;
IVcopy(nowned, mtx->rowind, IV_entries(bridge->myownedIV)) ;
mapIV = IV_new() ;
IV_init(mapIV, neqns, NULL) ;
IV_fill(mapIV, 0) ;
IVfill(4, stats, 0) ;
if ( bridge->msglvl > 1 ) {
   fprintf(bridge->msgFile, "\n calling DenseMtx_split()()") ;
   fflush(bridge->msgFile) ;
}
newmtx = DenseMtx_MPI_splitByRows(mtx, mapIV, stats, bridge->msglvl, 
                                  bridge->msgFile, tag, bridge->comm) ;
if ( bridge->msglvl > 1 ) {
   fprintf(bridge->msgFile, "\n return from DenseMtx_split()()") ;
   fflush(bridge->msgFile) ;
}
DenseMtx_free(mtx) ;
mtx = newmtx ;
IV_free(mapIV) ;
if ( myid == 0 ) {
   DVcopy(neqns*ncols, sol, DenseMtx_entries(mtx)) ;
}
DenseMtx_free(mtx) ;
/*
   ---------------------------------------------
   broadcast the entries to the other processors
   ---------------------------------------------
*/
if ( bridge->msglvl > 1 ) {
   fprintf(bridge->msgFile, "\n calling MPI_Bcast()()") ;
   fflush(bridge->msgFile) ;
}
MPI_Bcast((void *) sol, neqns*ncols, MPI_DOUBLE, 0, bridge->comm) ;
if ( bridge->msglvl > 1 ) {
   fprintf(bridge->msgFile, "\n return from MPI_Bcast()()") ;
   fflush(bridge->msgFile) ;
}
MPI_Barrier(bridge->comm) ;
/*
   ------------------------------------------------------------------
   set the error. (this is simple since when the spooles codes detect
   a fatal error, they print out a message to stderr and exit.)
   ------------------------------------------------------------------
*/
*perror = 0 ;
#if MYDEBUG > 0
MARKTIME(t2) ;
time_JimSolve += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) JimSolve() end", count_JimSolve) ;
   fprintf(stdout, ", %8.3f seconds, %8.3f total time",
           t2 - t1, time_JimSolve) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) JimSolve() end", count_JimSolve) ;
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_JimSolve) ;
fflush(bridge->msgFile) ;
#endif
 
return ; }

/*--------------------------------------------------------------------*/
