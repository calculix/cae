#include "../A2.h"
#include "../FrontMtx.h"
#include "../Drand.h"
#include "../SymbFac.h"
#include "../timings.h"
#include "../misc.h"
#include "../DenseMtx.h"
#include "../InpMtx.h"
#include "../Utilities.h"

#define  CONVER_TOL       1.0e-6

 
#define BiCGStabR    0
#define BiCGStabL    1
#define MLBiCGStabR  2
#define MLBiCGStabL  3
#define TFQMRR       4
#define TFQMRL       5
#define PCGR         6
#define PCGL         7
#define BGMRESR      8
#define BGMRESL      9

int
bicgstabr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
bicgstabl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
tfqmrr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
tfqmrl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
pcgr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
pcgl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
mlbicgstabr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxQ,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
mlbicgstabl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxQ,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );

int
bgmresr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             maxnouter,
   int             maxninner,
   int             *pnouter,
   int             *pninner,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );

int
bgmresl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             maxnouter,
   int             maxninner,
   int             *pnouter,
   int             *pninner,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
zbicgstabr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
zbicgstabl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
ztfqmrr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
ztfqmrl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
zpcgr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
zpcgl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
zmlbicgstabr (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxQ,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


int
zmlbicgstabl (
   int             n_matrixSize,
   int             type,
   int             symmetryflag,
   InpMtx          *mtxA,
   FrontMtx        *Precond,
   DenseMtx        *mtxX,
   DenseMtx        *mtxQ,
   DenseMtx        *mtxB,
   int             itermax,
   double          convergetol,
   int             msglvl,
   FILE            *msgFile
 );


 /*************************** NEW UTILITY ROUTINES ********************/
  
double
DenseMtx_frobNorm (
    DenseMtx *mtx
 );
 
double
DenseMtx_twoNormOfColumn (
    DenseMtx *mtx,
    int       jcol
 );
 
/*
   ------------------------------------------------------
   copy column icol of a Dense matrix A to column jcol of a 
   Dense matrix B.  Ai->Bj.
   ------------------------------------------------------
*/
void
DenseMtx_colCopy (
    DenseMtx *mtxB,   int jcol,
    DenseMtx *mtxA,   int icol 
 );
 
/*
   ------------------------------------------------------
   compute dot product of column icol of a Dense matrix A 
   and column jcol of a Dense matrix B
   prod=Ai^H * Bj.
   ------------------------------------------------------
*/
void
DenseMtx_colDotProduct (
    DenseMtx *mtxA,   int icol, 
    DenseMtx *mtxB,   int jcol,
    double *prod 
 );

/*
   ------------------------------------------------------
   compute a general axpy with column icol of a Dense matrix A 
   and column jcol of a Dense matrix B.  Ai=alpha*Ai+beta*Bj.
   ------------------------------------------------------
*/
void
DenseMtx_colGenAxpy (
    double    *alpha,  DenseMtx *mtxA,   int icol,   
    double    *beta,   DenseMtx *mtxB,   int jcol
 );
 
/*
   -----------------------------------------------
   copy col icolA from mtxA into col icolB in mtxB
   -----------------------------------------------
*/
void
DenseMtx_copyCoulmn (
    DenseMtx   *mtxB,
    int        icolB,
    DenseMtx   *mtxA,
    int        icolA
 );
 
/*
   -----------------------------------------------
   FrontMtx_solve with column icol of rhsmtx as right-hand-side.
   -----------------------------------------------
*/ 
void
FrontMtx_solveOneColumn (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   int             jcol,
   DenseMtx        *rhsmtx,
   int             icol,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
);

/*
   -------------------------------------------
   performs the matrix-matrix operations
   C =  beta*C + alpha*(A)*(B)
   A, B and C must be column major.
*/
int 
DenseMtx_mmm(
  char     *A_opt,
  char     *B_opt,
  double   *beta,
  DenseMtx *mtxC, 
  double   *alpha,
  DenseMtx *mtxA, 
  DenseMtx *mtxB
  ); 


/*
   ----------------------------------------------
   return the absolute value of a complex number
   ----------------------------------------------
*/
double 
zabs(double *x);

/*
   -------------------------------------------
   compute the sum of two complex numbers
   -------------------------------------------
*/
void 
zadd
(double *x, double *y, double *u
);

/*
   -------------------------------------------
   divide two complex numbers                  
   -------------------------------------------
*/
void
zdiv
(double *x, double *y, double *u
);

/*
   -------------------------------------------
   multiply two complex numbers
   -------------------------------------------
*/
void
zmul(double *x, double *y,  double *u
);

/*
   ---------------------------------------------
   compute the difference of two complex numbers 
   ---------------------------------------------
*/
void 
zsub
(double *x, double *y, double *u
);
