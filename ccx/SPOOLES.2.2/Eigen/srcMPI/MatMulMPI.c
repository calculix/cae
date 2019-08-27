/*  MatMulMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_MatMul = 0 ;
static double time_MatMul = 0.0 ;
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

   created -- 98aug11, cca & jcp
   -------------------------------------------------------------
*/
void 
MatMulMPI ( 
   int      *pnrows, 
   int      *pncols, 
   double   x[], 
   double   y[],
   int      *pprbtype,
   void     *data
) {
BridgeMPI   *bridge = (BridgeMPI *) data ; 
int         ncols, nent, nrows ;
#if MYDEBUG > 0
double   t1, t2 ;
count_MatMul++ ;
MARKTIME(t1) ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) MatMulMPI()", count_MatMul) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) MatMulMPI()", count_MatMul) ;
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
   if ( x == NULL ) {
      fprintf(stderr, "\n\n fatal error in MatMulMPI, y <-- x, x is NULL") ;
      exit(-1) ;
   }
   if ( y == NULL ) {
      fprintf(stderr, "\n\n fatal error in MatMulMPI, y <-- x, y is NULL") ;
      exit(-1) ;
   }
   DVcopy(nent, y, x) ;
   return;
} else {
   DenseMtx    *Xloc = bridge->Xloc, *Yloc = bridge->Yloc ;
   double      alpha[2] = {1.0, 0.0} ;
   FILE        *msgFile = bridge->msgFile ;
   int         msglvl = bridge->msglvl, n, nmyowned, tag = 0 ;
   int         *list, *owned ;
   int         stats[4] ;
   MPI_Comm    comm = bridge->comm ;
/*
   -----------------------------------------------
   if the matrices are in global coordinates
   (i.e., this is the first matrix-vector multiply
    following a factorization) then
   map the matrix into local coordinates
   -----------------------------------------------
*/
   if ( bridge->msglvl > 1 ) {
      fprintf(bridge->msgFile, 
              "\n\n inside MatMulMPI, nrow = %d, ncol = %d",
              nrows, ncols) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
              "\n\n bridge->coordFlag = %d", bridge->coordFlag) ;
      fflush(msgFile) ;
   }
   if ( bridge->coordFlag == GLOBAL ) {
      if ( bridge->prbtype == 1 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n ready to permute B") ;
            fflush(msgFile) ;
         }
         MatMul_setLocalIndices(bridge->info, bridge->B) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n matrix B in local coordinates") ;
            InpMtx_writeForHumanEye(bridge->B, msgFile) ;
            fflush(msgFile) ;
         }
      } else if ( bridge->prbtype == 2 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n ready to permute A") ;
            fflush(msgFile) ;
         }
         MatMul_setLocalIndices(bridge->info, bridge->A) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n matrix A in local coordinates") ;
            InpMtx_writeForHumanEye(bridge->A, msgFile) ;
            fflush(msgFile) ;
         }
      }
      bridge->coordFlag = LOCAL ;
   }
/*
   --------------------------------------------------
   check to see that Xloc and Yloc are the right size
   --------------------------------------------------
*/
   if ( Xloc->nrow != nrows ) {
      fprintf(stderr, 
              "\n\n fatal error in MatMulMPI, nrows %d, Xloc->nrow %d",
              nrows, Xloc->nrow) ;
      exit(-1) ;
   }
   if ( Xloc->ncol != ncols ) {
      IV_sizeAndEntries(bridge->myownedIV, &nmyowned, &owned) ;
      DenseMtx_clearData(Xloc) ;
      DenseMtx_init(Xloc, SPOOLES_REAL, 0, 0, 
                    nmyowned, ncols, 1, nmyowned) ;
      DenseMtx_rowIndices(Xloc, &n, &list) ;
      IVcopy(n, list, owned) ;
      DenseMtx_clearData(Yloc) ;
      DenseMtx_init(Yloc, SPOOLES_REAL, 0, 0, 
                    nmyowned, ncols, 1, nmyowned) ;
      DenseMtx_rowIndices(Yloc, &n, &list) ;
      IVcopy(n, list, owned) ;
   }
/*
   ------------------
   copy x[] into Xloc
   ------------------
*/
   DVcopy(nent, DenseMtx_entries(Xloc), x) ;
/*
   ---------
   zero Yloc
   ---------
*/
   DenseMtx_zero(Yloc) ;
/*
   ----------------------------------------
   compute the local matrix-vector multiply
   ----------------------------------------
*/
   if ( *pprbtype == 1 ) {
      IVzero(4, stats) ;
      MatMul_MPI_mmm(bridge->info, Yloc, alpha, bridge->B,
                     Xloc, stats, msglvl, msgFile, tag, comm) ;
   } else if ( *pprbtype == 2 ) {
      IVzero(4, stats) ;
      MatMul_MPI_mmm(bridge->info, Yloc, alpha, bridge->A,
                     Xloc, stats, msglvl, msgFile, tag, comm) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n after mvm, Yloc") ;
      DenseMtx_writeForHumanEye(Yloc, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------------------
   copy entries of Yloc into y[]
   -----------------------------
*/
   if ( DenseMtx_entries(Yloc) == NULL ) {
      fprintf(stderr, 
              "\n\n fatal error in MatMulMPI, y <-- Yloc, Yloc is NULL") ;
      exit(-1) ;
   }
   if ( y == NULL ) {
      fprintf(stderr, "\n\n fatal error in MatMulMPI, y <-- Yloc, y is NULL") ;
      exit(-1) ;
   }
   DVcopy(nent, y, DenseMtx_entries(Yloc)) ;
}
#if MYDEBUG > 0
MARKTIME(t2) ;
time_MatMul += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, ", %8.3f seconds, %8.3f total time",
           t2 - t1, time_MatMul) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_MatMul) ;
fflush(bridge->msgFile) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
