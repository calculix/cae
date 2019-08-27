/*  SolveMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_Solve = 0 ;
static double time_Solve = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to solve a linear system
     (A - sigma*B) sol[] = rhs[]

   data    -- pointer to bridge data object
   *pnrows -- # of rows in x[] and y[]
   *pncols -- # of columns in x[] and y[]
   rhs[]   -- vector that holds right hand sides
   sol[]   -- vector to hold solutions

   note: rhs[] and sol[] can be the same array.
  
   on return, *perror holds an error code.
      1 -- normal return
     -1 -- pnrows is NULL
     -2 -- pncols is NULL
     -3 -- rhs is NULL
     -4 -- sol is NULL
     -5 -- data is NULL

   created -- 98aug10, cca & jcp
   ----------------------------------------------
*/
void 
SolveMPI ( 
   int       *pnrows, 
   int       *pncols, 
   double    rhs[], 
   double    sol[],
   void      *data, 
   int       *perror 
) {
BridgeMPI   *bridge = (BridgeMPI *) data ;
DenseMtx    *rhsmtx ;
double      cpus[20] ;
int         nent, ncols, nrows, tag = 0 ;
int         stats[8] ;
#if MYDEBUG > 0
double   t1, t2 ;
count_Solve++ ;
MARKTIME(t1) ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) SolveMPI()", count_Solve) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) SolveMPI()", count_Solve) ;
fflush(bridge->msgFile) ;
#endif
/*
   ---------------
   check the input
   ---------------
*/
if ( perror == NULL ) {
   fprintf(stderr, "\n error in Solve()"
           "\n perror == NULL\n") ;
   return ;
}
if ( pnrows == NULL ) {
   fprintf(stderr, "\n error in Solve()"
           "\n pnrows == NULL\n") ;
   *perror = -1 ; return ;
}
if ( pncols == NULL ) {
   fprintf(stderr, "\n error in Solve()"
           "\n pncols == NULL\n") ;
   *perror = -2 ; return ;
}
if ( rhs == NULL ) {
   fprintf(stderr, "\n error in Solve()"
           "\n rhs == NULL\n") ;
   *perror = -3 ; return ;
}
if ( sol == NULL ) {
   fprintf(stderr, "\n error in Solve()"
           "\n sol == NULL\n") ;
   *perror = -4 ; return ;
}
if ( data == NULL ) {
   fprintf(stderr, "\n error in Solve()"
           "\n data == NULL\n") ;
   *perror = -5 ; return ;
}
/*
   ----------------------------------
   set the number of rows and columns
   ----------------------------------
*/
nrows = *pnrows ;
ncols = *pncols ;
nent  = nrows*ncols ;
/*
   ---------------------------------
   setup rhsmtx as a DenseMtx object
   ---------------------------------
*/
rhsmtx = DenseMtx_new() ;
DenseMtx_init(rhsmtx, SPOOLES_REAL, 0, 0, nrows, ncols, 1, nrows) ;
if ( rhs == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMPI, rhsmtx <-- rhs, rhs is NULL") ;
   exit(-1) ;
}
if ( DenseMtx_entries(rhsmtx) == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SolveMPI, rhsmtx <-- rhs, rhsmtx is NULL") ;
   exit(-1) ;
}
DVcopy (nent, DenseMtx_entries(rhsmtx), rhs) ;
IVcopy(nrows, rhsmtx->rowind, IV_entries(bridge->myownedIV)) ;
if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
   fprintf(bridge->msgFile, "\n\n rhsmtx initialized") ;
   DenseMtx_writeForHumanEye(rhsmtx, bridge->msgFile) ;
   fflush(bridge->msgFile) ;
}
if ( bridge->rowmapIV != NULL ) {
   DenseMtx   *newmtx ;
/*
   --------------------------------
   pivoting may have taken place, 
   redistribute the rows of the rhs
   --------------------------------
*/
   IVzero(4, stats) ;
   newmtx = DenseMtx_MPI_splitByRows(rhsmtx, bridge->rowmapIV, stats,
                   bridge->msglvl, bridge->msgFile, tag, bridge->comm) ;
   DenseMtx_free(rhsmtx) ;
   rhsmtx = newmtx ;
   if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
      fprintf(bridge->msgFile, "\n\n rhsmtx after redistribution") ;
      DenseMtx_writeForHumanEye(rhsmtx, bridge->msgFile) ;
      fflush(bridge->msgFile) ;
   }
}
/*
   -----------------------
   solve the linear system
   -----------------------
*/
IVzero(8, stats) ;
DVzero(20, cpus) ;
FrontMtx_MPI_solve(bridge->frontmtx, rhsmtx, rhsmtx, bridge->mtxmanager,
                   bridge->solvemap, cpus, stats, bridge->msglvl, 
                   bridge->msgFile, tag, bridge->comm) ;
if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
   fprintf(bridge->msgFile, 
           "\n\n solution matrix after redistribution") ;
   DenseMtx_writeForHumanEye(rhsmtx, bridge->msgFile) ;
   fflush(bridge->msgFile) ;
}
if ( bridge->rowmapIV != NULL ) {
   DenseMtx   *newmtx ;
/*
   -------------------------------------
   pivoting may have taken place, 
   redistribute the rows of the solution
   -------------------------------------
*/
   newmtx = DenseMtx_MPI_splitByRows(rhsmtx, bridge->vtxmapIV, stats,
                   bridge->msglvl, bridge->msgFile, tag, bridge->comm) ;
   DenseMtx_free(rhsmtx) ;
   rhsmtx = newmtx ;
   if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
      fprintf(bridge->msgFile, 
              "\n\n solution matrix after redistribution") ;
      DenseMtx_writeForHumanEye(rhsmtx, bridge->msgFile) ;
      fflush(bridge->msgFile) ;
   }
}
/*
   -----------------------------------
   copy solution into output parameter
   -----------------------------------
*/
if ( DenseMtx_entries(rhsmtx) == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SolveMPI, sol <-- rhsmtx, rhsmtx is NULL") ;
   exit(-1) ;
}
if ( sol == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMPI, sol <-- rhsmtx, sol is NULL") ;
   exit(-1) ;
}
DVcopy(nent, sol, DenseMtx_entries(rhsmtx)) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
DenseMtx_free(rhsmtx) ;
/*
   ------------------------------------------------------------------
   set the error. (this is simple since when the spooles codes detect
   a fatal error, they print out a message to stderr and exit.)
   ------------------------------------------------------------------
*/
*perror = 0 ;

#if MYDEBUG > 0
MARKTIME(t2) ;
time_Solve += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, ", %8.3f seconds, %8.3f total time",
           t2 - t1, time_Solve) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_Solve) ;
fflush(bridge->msgFile) ;
#endif
 
return ; }

/*--------------------------------------------------------------------*/
