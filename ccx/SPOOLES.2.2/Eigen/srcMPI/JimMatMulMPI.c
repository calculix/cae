/*  JimMatMulMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_JimMatMul = 0 ;
static double time_JimMatMul = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose --- to compute a matrix-vector multiply y[] = C * x[]
     where C is the identity, A or B (depending on *pprbtype).

   *pnrows -- # of rows in x[]
   *pncols -- # of columns in x[]
   *pprbtype -- problem type
      *pprbtype = 1 --> vibration problem, matrix is A
      *pprbtype = 2 --> buckling problem, matrix is B
      *pprbtype = 3 --> matrix is identity, y[] = x[]
   x[] -- vector to be multiplied
      NOTE: the x[] vector is global, not a portion
   y[] -- product vector
      NOTE: the y[] vector is global, not a portion

   created -- 98aug28, cca & jcp
   -------------------------------------------------------------
*/
void 
JimMatMulMPI ( 
   int      *pnrows, 
   int      *pncols, 
   double   x[], 
   double   y[],
   int      *pprbtype,
   void     *data
) {
BridgeMPI   *bridge = (BridgeMPI *) data ;
int   ncols, nent, nrows ;
#if MYDEBUG > 0
double   t1, t2 ;
count_JimMatMul++ ;
MARKTIME(t1) ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) JimMatMulMPI() start", count_JimMatMul) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, 
        "\n (%d) JimMatMulMPI() start", count_JimMatMul) ;
fflush(bridge->msgFile) ;
#endif

nrows = *pnrows ;
ncols = *pncols ;
nent  = nrows*ncols ;
if ( *pprbtype == 3 ) {
/*
    --------------------------
    ... matrix is the identity
    --------------------------
*/
   DVcopy(nent, y, x) ;
} else {
   BridgeMPI   *bridge = (BridgeMPI *) data ; 
   DenseMtx    *mtx, *newmtx ;
   int         irow, jcol, jj, kk, myid, neqns, nowned, tag = 0 ;
   int         *vtxmap ;
   int         stats[4] ;
   IV          *mapIV ;
/*
   ---------------------------------------------
   slide the owned rows of x[] down in the array
   ---------------------------------------------
*/
   vtxmap  = IV_entries(bridge->vtxmapIV) ;
   neqns   = bridge->neqns ;
   myid    = bridge->myid  ;
   nowned  = IV_size(bridge->myownedIV) ;
   for ( jcol = jj = kk = 0 ; jcol < ncols ; jcol++ ) {
      for ( irow = 0 ; irow < neqns ; irow++, jj++ ) {
         if ( vtxmap[irow] == myid ) {
            y[kk++] = x[jj] ;
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
   if ( bridge->msglvl > 2 ) {
      fprintf(bridge->msgFile, 
              "\n inside JimMatMulMPI, calling MatMulMpi"
              "\n prbtype %d, nrows %d, ncols %d, nowned %d",
              *pprbtype, *pnrows, *pncols, nowned) ;
      fflush(bridge->msgFile) ;
   }
   MatMulMPI(&nowned, pncols, y, y, pprbtype, data) ;
/*
   -------------------------------------------------
   gather all the entries of y[] onto processor zero
   -------------------------------------------------
*/
   mtx = DenseMtx_new() ;
   DenseMtx_init(mtx, SPOOLES_REAL, 0, 0, nowned, ncols, 1, nowned) ;
   DVcopy (nowned*ncols, DenseMtx_entries(mtx), y) ;
   IVcopy(nowned, mtx->rowind, IV_entries(bridge->myownedIV)) ;
   mapIV = IV_new() ;
   IV_init(mapIV, neqns, NULL) ;
   IV_fill(mapIV, 0) ;
   IVfill(4, stats, 0) ;
   if ( bridge->msglvl > 2 ) {
      fprintf(bridge->msgFile, "\n mtx: %d rows x %d columns",
              mtx->nrow, mtx->ncol) ;
      fflush(bridge->msgFile) ;
   }
   newmtx = DenseMtx_MPI_splitByRows(mtx, mapIV, stats, bridge->msglvl,
                                   bridge->msgFile, tag, bridge->comm) ;
   if ( bridge->msglvl > 2 ) {
      fprintf(bridge->msgFile, "\n newmtx: %d rows x %d columns",
              newmtx->nrow, newmtx->ncol) ;
      fflush(bridge->msgFile) ;
   }
   DenseMtx_free(mtx) ;
   mtx = newmtx ;
   IV_free(mapIV) ;
   if ( myid == 0 ) {
      if ( mtx->nrow != neqns || mtx->ncol != ncols ) {
         fprintf(bridge->msgFile, 
                 "\n\n WHOA: mtx->nrows %d, mtx->ncols %d"
                 ", neqns %d, ncols %d", mtx->nrow, mtx->ncol,
                 neqns, ncols) ;
         exit(-1) ;
      }
      DVcopy(neqns*ncols, y, DenseMtx_entries(mtx)) ;
   }
   DenseMtx_free(mtx) ;
/*
   ---------------------------------------------
   broadcast the entries to the other processors
   ---------------------------------------------
*/
   MPI_Bcast((void *) y, neqns*ncols, MPI_DOUBLE, 0, bridge->comm) ;
   if ( bridge->msglvl > 2 ) {
      fprintf(bridge->msgFile, "\n after the broadcast") ;
      fflush(bridge->msgFile) ;
   }
}
MPI_Barrier(bridge->comm) ;
#if MYDEBUG > 0
MARKTIME(t2) ;
time_JimMatMul += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) JimMatMulMPI() end", count_JimMatMul) ;
   fprintf(stdout, ", %8.3f seconds, %8.3f total time",
           t2 - t1, time_JimMatMul) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, 
        "\n (%d) JimMatMulMPI() end", count_JimMatMul) ;
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_JimMatMul) ;
fflush(bridge->msgFile) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
