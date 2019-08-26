/*  pcgr.c  */

#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to solve a symmetric nonnegative  matrix equation

               Ax=b
   using right preconditioned conjugate gradient method 

      x       -- Initial guess
      A       -- Input matrix
      M       -- Front Matrix as the preconditioner
      b       -- Right-hand side
      tol     -- Convergence tolerance
      type -- type of entries
         SPOOLES_REAL or SPOOLES_COMPLEX
      symmetryflag -- symmetry of the matrix
         SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
      nrhs    -- number of right hand sides
      msglvl  -- message level
      msgFile -- message file

      return value  -- final residual norm

   created -- Oct. 27, 1998    Wei-Pai Tang
   ---------------------------------------------------------------------
*/

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
 )
{
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxZ ;
DenseMtx        *vecP, *vecR, *vecQ ;
DenseMtx        *vecX,  *vecZ  ;
double          Alpha, Beta, Rho, Rho0, Init_norm, ratio, Res_norm, Rtmp ;
double          t1, t2,  cpus[9] ;
double          one[2] = {1.0, 0.0} ;
double          Tiny = 0.1e-28;
int             Iter, neqns;
int             stats[6] ;

if (symmetryflag != SPOOLES_SYMMETRIC){
      fprintf(msgFile, "\n\n Fatal Error, \n"
                    " Matrix is not symmetric in PCGR !!") ;
       return(-1);
    };

neqns = n_matrixSize;

/*
   --------------------
   init the vectors in bicgstab
   --------------------
*/
vecP = DenseMtx_new() ;
DenseMtx_init(vecP, type, 0, 0, neqns, 1, 1, neqns) ;

vecR = DenseMtx_new() ;
DenseMtx_init(vecR, type, 0, 0, neqns, 1, 1, neqns) ;

vecX = DenseMtx_new() ;
DenseMtx_init(vecX, type, 0, 0, neqns, 1, 1, neqns) ;

vecQ = DenseMtx_new() ;
DenseMtx_init(vecQ, type, 0, 0, neqns, 1, 1, neqns) ;

vecZ = DenseMtx_new() ;
DenseMtx_init(vecZ, type, 0, 0, neqns, 1, 1, neqns) ;


/*
   --------------------------
   Initialize the iterations
   --------------------------
*/
Init_norm = DenseMtx_twoNormOfColumn(mtxB,0);
if ( Init_norm == 0.0 ){
  Init_norm = 1.0; };
ratio = 1.0;
DenseMtx_zero(vecX) ;

DenseMtx_colCopy(vecR, 0, mtxB, 0);


MARKTIME(t1) ;
Iter = 0;

/*
   ------------------------------
    Main Loop of the iterations
   ------------------------------
*/

while ( ratio > convergetol && Iter <= itermax )
  {
    Iter++;
/*                                                         */
    FrontMtx_solve(Precond, vecZ, vecR, Precond->manager,
               cpus, msglvl, msgFile) ;
    DenseMtx_colDotProduct(vecR, 0, vecZ, 0, &Rho);
/*                                                         */
    if ( Iter == 1 ) {
      DenseMtx_colCopy(vecP, 0, vecZ, 0);
    } else {
      Beta = Rho /(Rho0 + Tiny);
      DenseMtx_colGenAxpy(&Beta, vecP, 0, one, vecZ, 0);
    };

    DenseMtx_zero(vecQ);
    InpMtx_sym_mmm(mtxA, vecQ, one, vecP) ;
    DenseMtx_colDotProduct(vecP, 0, vecQ, 0, &Rtmp);
    Alpha = Rho/(Rtmp+Tiny);
    DenseMtx_colGenAxpy(one, vecX, 0, &Alpha, vecP, 0);
    Rtmp=-Alpha;
    DenseMtx_colGenAxpy(one, vecR, 0, &Rtmp, vecQ, 0);
    Rho0  = Rho;

    /*                                                */
    Res_norm =  DenseMtx_twoNormOfColumn(vecR,0);
    ratio = Res_norm/Init_norm;
    fprintf(msgFile, "\n\n At iteration %d"
	    "  the convergence ratio is  %12.4e", 
	    Iter, ratio) ;
  }
/*            End of while loop              */
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU  : Converges in time: %8.3f ", t2 - t1) ;
fprintf(msgFile, "\n # iterations = %d", Iter) ;

fprintf(msgFile, "\n\n after PCGR") ;
DenseMtx_colCopy(mtxX, 0, vecX, 0);

/*
 
   ------------------------
   free the working storage
   ------------------------
*/

DenseMtx_free(vecP) ;
DenseMtx_free(vecR) ;
DenseMtx_free(vecX) ;
DenseMtx_free(vecQ) ;
DenseMtx_free(vecZ) ;


fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
