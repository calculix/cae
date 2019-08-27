/*  Solve.c  */

#include "../Bridge.h"

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
Solve ( 
   int       *pnrows, 
   int       *pncols, 
   double    rhs[], 
   double    sol[],
   void      *data, 
   int       *perror 
) {
Bridge     *bridge = (Bridge *) data ;
DenseMtx   *rhsmtx, *solmtx ;
double     cpus[10] ;
int        nent, ncols, nrows ;
#if MYDEBUG > 0
double   t1, t2 ;
MARKTIME(t1) ;
count_Solve++ ;
fprintf(stdout, "\n (%d) Solve()", count_Solve) ;
fflush(stdout) ;
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
   -------------------------------------------
   setup rhsmtx and solmtx as DenseMtx objects
   ------------------------------------------
*/
rhsmtx = bridge->Y ;
DenseMtx_init(rhsmtx, SPOOLES_REAL, 0, 0, nrows, ncols, 1, nrows) ;
DVcopy (nent, DenseMtx_entries(rhsmtx), rhs) ;
solmtx = bridge->X ;
DenseMtx_init(solmtx, SPOOLES_REAL, 0, 0, nrows, ncols, 1, nrows) ;
DenseMtx_zero(solmtx) ;
/*
   -----------------------
   solve the linear system
   -----------------------
*/
DVzero(10, cpus) ;
FrontMtx_solve(bridge->frontmtx, solmtx, rhsmtx, bridge->mtxmanager,
               cpus, bridge->msglvl, bridge->msgFile) ;
/*
   -----------------------------------
   copy solution into output parameter
   -----------------------------------
*/
DVcopy(nent, sol, DenseMtx_entries(solmtx)) ;
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
fprintf(stdout, ", %8.3f seconds, %8.3f total time", 
        t2 - t1, time_Solve) ;
fflush(stdout) ;
#endif
 
return ; }

/*--------------------------------------------------------------------*/
