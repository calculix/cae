/*  spoolesMT.h  */

#include "../FrontMtx.h"
#include "../Lock.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   parallel factorization method for A.
   all but two input parameters are the same as the serial method.
   this is a wrapper method around FrontMtx_MT_factorInpMtx().
 
   ownersIV  -- pointer to IV object that holds map from fronts
                to threads
   lookahead -- lookahead parameter that allows computation at
                higher levels of the front tree to proceed when
                lower fronts are not yet finish. use lookahead = 0
                to turn off this option. otherwise lookahead ancestors
                of an active unfinished front can be active.
 
   return value -- pointer to the first Chv object in a list
                   that contains postponed data
 
   created -- 98may16, cca
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_MT_factorInpMtx (
   FrontMtx     *frontmtx,
   InpMtx       *inpmtx,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   IV           *ownersIV,
   int          lookahead,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) ;
/*
   -------------------------------------------------------------------
   parallel factorization method for A + sigma*B.
   all but two input parameters are the same
   as the FrontMtx_factorPencil method.
 
   ownersIV  -- pointer to IV object that holds map from fronts
                to threads
   lookahead -- lookahead parameter that allows computation at
                higher levels of the front tree to proceed when
                lower fronts are not yet finish. use lookahead = 0
                to turn off this option. otherwise lookahead ancestors
                of an active unfinished front can be active.
 
   return value -- pointer to the first Chv object in a list
                   that contains postponed data
 
   created -- 98may16, cca
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_MT_factorPencil (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   IV           *ownersIV,
   int          lookahead,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) ;
/*
   ----------------------------------------------------
   multithreaded solve method for (L + I)D(I + U) X = B
   or (U^T + I)D(I + U) X = B.
 
   created -- 98mar19, cca
   ----------------------------------------------------
*/
void
FrontMtx_MT_solve (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   DenseMtx        *rhsmtx,
   SubMtxManager   *mtxmanager,
   SolveMap        *solvemap,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) ;
/*
   ------------------------------------------------------------
   purpose -- compute a QR factorization using multiple threads
 
   created -- 98may29, cca
   ------------------------------------------------------------
*/
void
FrontMtx_MT_QR_factor (
   FrontMtx     *frontmtx,
   InpMtx       *mtxA,
   ChvManager   *chvmanager,
   IV           *ownersIV,
   double       cpus[],
   double       *pfacops,
   int          msglvl,
   FILE         *msgFile
) ;
/*
   --------------------------------------------------------
   multithreaded version:
 
   minimize ||b - Ax||_2 by
   solving (U^T + I) * D * (I + U) X = A^T * B
      where A = QR = QD(I + U)
   by calling FrontMtx_solve()
 
   note: if A is square, mtxX and mtxB must not be the same
 
   cpus     -- vector of cpu time breakdowns
      cpus[0] -- set up solves
      cpus[1] -- fetch rhs and store solution
      cpus[2] -- forward solve
      cpus[3] -- diagonal solve
      cpus[4] -- backward solve
      cpus[5] -- total solve time
      cpus[6] -- time to compute A^T * B
      cpus[7] -- total time
 
   created -- 98may30, cca
   --------------------------------------------------------
*/
void
FrontMtx_MT_QR_solve (
   FrontMtx        *frontmtx,
   InpMtx          *mtxA,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   SubMtxManager   *mtxmanager,
   SolveMap        *solvemap,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_MT_nonsym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_sym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_herm_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_nonsym_mmm_T (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_nonsym_mmm_H (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
