/*  TFQMRL.c  */
#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to solve a unsymmetric real matrix equation

               Ax=b

   using Left preconditioned TFQMR method without lookahead 

      x       -- Initial guess as zeros
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

      return value  -- error flag

   created -- Oct. 28, 1998              Wei-Pai Tang
   ---------------------------------------------------------------------
*/

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
 )
{
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *vecD, *vecR, *vecT, *vecU1, *vecU2,  *vecV, *vecW;
DenseMtx        *vecX, *vecY1, *vecY2 ;
double          Alpha, Beta, Cee, Eta, Rho, Rho_new ;
double          Sigma, Tau, Theta;
double          Init_norm,  ratio,  Res_norm;
double          error_trol, m, Rtmp;
double          t1, t2,  cpus[9] ;
double          one[2] = {1.0, 0.0}, zero[2] ={0.0, 0.0} ;
double          Tiny = 0.1e-28;
int             Iter, Imv, neqns;
int             stats[6] ;

neqns = n_matrixSize;


/*
   --------------------
   init the vectors in TFQMRL
   --------------------
*/
vecD = DenseMtx_new() ;
DenseMtx_init(vecD, type, 0, 0, neqns, 1, 1, neqns) ;

vecR = DenseMtx_new() ;
DenseMtx_init(vecR, type, 0, 0, neqns, 1, 1, neqns) ;


vecT = DenseMtx_new() ;
DenseMtx_init(vecT, type, 0, 0, neqns, 1, 1, neqns) ;

vecU1 = DenseMtx_new() ;
DenseMtx_init(vecU1, type, 0, 0, neqns, 1, 1, neqns) ;

vecU2 = DenseMtx_new() ;
DenseMtx_init(vecU2, type, 0, 0, neqns, 1, 1, neqns) ;

vecV = DenseMtx_new() ;
DenseMtx_init(vecV, type, 0, 0, neqns, 1, 1, neqns) ;

vecW = DenseMtx_new() ;
DenseMtx_init(vecW, type, 0, 0, neqns, 1, 1, neqns) ;

vecX = DenseMtx_new() ;
DenseMtx_init(vecX, type, 0, 0, neqns, 1, 1, neqns) ;

vecY1 = DenseMtx_new() ;
DenseMtx_init(vecY1, type, 0, 0, neqns, 1, 1, neqns) ;

vecY2 = DenseMtx_new() ;
DenseMtx_init(vecY2, type, 0, 0, neqns, 1, 1, neqns) ;


/*
   --------------------------
   Initialize the iterations
   --------------------------
*/
/*          ----     Set initial guess as zero  ----               */
DenseMtx_zero(vecX) ;

DenseMtx_colCopy(vecT, 0, mtxB, 0);
/*                                                         */
    FrontMtx_solve(Precond, vecR, vecT, Precond->manager,
               cpus, msglvl, msgFile) ;
/*                                                      */

  
Init_norm = DenseMtx_twoNormOfColumn(vecR,0);
if ( Init_norm == 0.0 ){
  Init_norm = 1.0; 
};
error_trol = Init_norm * convergetol ;

  fprintf(msgFile, "\n TFQMRL Initial norml: %6.2e ", Init_norm ) ;
  fprintf(msgFile, "\n TFQMRL Conveg. Control: %7.3e ", convergetol ) ;
  fprintf(msgFile, "\n TFQMRL Convergen Control: %7.3e ",error_trol ) ;

DenseMtx_zero(vecD) ;
DenseMtx_zero(vecU1) ;
DenseMtx_zero(vecU2) ;
DenseMtx_zero(vecY2) ;

/*       DenseMtx_copy(vecR, mtxB);              */
DenseMtx_colCopy(vecW, 0, vecR, 0);
DenseMtx_colCopy(vecY1, 0, vecR, 0);

Iter = 0;
Imv  = 0;


  switch ( symmetryflag ) {
  case SPOOLES_SYMMETRIC : 
    InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecY1) ;
    break ;
  case SPOOLES_HERMITIAN :
    fprintf(msgFile, "\n TFQMRL Matrix type wrong");
    fprintf(msgFile, "\n Fatal error");
    goto end;
  case SPOOLES_NONSYMMETRIC :
      InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecY1) ;
    break ;
  default :
    fprintf(msgFile, "\n TFQMRL Matrix type wrong");
    fprintf(msgFile, "\n Fatal error");
    goto end;
  }
/*                                                         */
    FrontMtx_solve(Precond, vecV, vecT, Precond->manager,
               cpus, msglvl, msgFile) ;
/*                                                      */
    Imv++;
    DenseMtx_colCopy(vecU1, 0, vecV, 0);
/*

*/
Theta   = 0.0;
Eta     = 0.0;
Tau     = Init_norm ;
Rho     = Tau * Tau ;

/*
   ------------------------------
   TFQMRL   Iteration start
   ------------------------------
*/

MARKTIME(t1) ;


while (  Iter <= itermax )
  {
    Iter++;
    DenseMtx_colDotProduct(vecV, 0, vecR, 0, &Sigma);

    if (Sigma == 0){
      fprintf(msgFile, "\n\n Fatal Error, \n"
	      "  TFQMRL Breakdown, Sigma = 0 !!") ;
      Imv = -1;
      goto end;
    };

    Alpha   = Rho/Sigma;
/*
    ----------------
    Odd step
    ---------------
*/
	
    m      = 2 * Iter - 1;
    Rtmp=-Alpha;
    DenseMtx_colGenAxpy(one, vecW, 0, &Rtmp, vecU1, 0);
    Rtmp   = Theta * Theta * Eta / Alpha ;
    DenseMtx_colGenAxpy(&Rtmp, vecD, 0, one, vecY1, 0);
    Theta  = DenseMtx_twoNormOfColumn(vecW,0)/Tau;
    Cee    = 1.0/sqrt(1.0 + Theta*Theta);
    Tau    = Tau * Theta * Cee ;
    Eta    = Cee * Cee * Alpha ;
    DenseMtx_colGenAxpy(one, vecX, 0, &Eta, vecD, 0);
      fprintf(msgFile, "\n\n Odd step at %d", Imv);
      fprintf(msgFile, " \n Tau is   : %7.3e", Tau) ; 
/*                   
        Debug purpose:  Check the convergence history
	for the true residual norm
*/
/*
      DenseMtx_zero(vecT) ;
      InpMtx_nonsym_mmm(mtxA, vecT, one, vecX) ;
      DenseMtx_sub(vecT, mtxB) ;
      Rtmp = DenseMtx_twoNormOfColumn(vecT,0);
      fprintf(msgFile, "\n TFQMRL Residual norm: %6.2e ", Rtmp) ;
*/
 
/*
    ----------------
    Convergence Test
    ---------------
*/
    if (Tau * sqrt(m + 1)  <= error_trol ) {
/*                                                             */
      DenseMtx_colCopy(mtxX, 0, vecX, 0);
/*
      DenseMtx_zero(vecT) ;
      InpMtx_nonsym_mmm(mtxA, vecT, one, mtxX) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      default :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
      DenseMtx_sub(vecT, mtxB) ;
      Rtmp = DenseMtx_twoNormOfColumn(vecT,0);
      fprintf(msgFile, "\n TFQMRL Residual norm: %6.2e ", Rtmp) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU  : Converges in time: %8.3f ", t2 - t1) ;
      fprintf(msgFile, "\n # iterations = %d", Imv) ;
      fprintf(msgFile, "\n\n after TFQMRL") ;  
      goto end;
    };

/*
    ----------------
    Even step
    ---------------
*/
    DenseMtx_colCopy(vecY2, 0, vecY1, 0);
    Rtmp=-Alpha;
    DenseMtx_colGenAxpy(one, vecY2, 0, &Rtmp, vecV, 0);
/*
    DenseMtx_zero(vecT) ;
    InpMtx_nonsym_mmm(mtxA, vecT, one, vecY2) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecY2) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecY2) ;
	break ;
      default :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
    
    FrontMtx_solve(Precond, vecU2, vecT, Precond->manager,
		   cpus, msglvl, msgFile) ;
    Imv++;
  
    m      = 2 * Iter ;
    Rtmp = -Alpha;
    DenseMtx_colGenAxpy(one, vecW, 0, &Rtmp, vecU2, 0);
    Rtmp   = Theta * Theta * Eta / Alpha ;
    DenseMtx_colGenAxpy(&Rtmp, vecD, 0, one, vecY2, 0);
    Theta  = DenseMtx_twoNormOfColumn(vecW,0)/Tau;
    Cee    = 1.0/sqrt(1.0 + Theta*Theta);
    Tau    = Tau * Theta * Cee ;
    Eta    = Cee * Cee * Alpha ;
    DenseMtx_colGenAxpy(one, vecX, 0, &Eta, vecD, 0);
      fprintf(msgFile, "\n\n Even step at %d", Imv) ;  
    
/*
    ----------------
    Convergence Test for even step
    ---------------
*/
    if (Tau * sqrt(m + 1)  <= error_trol ) {
      DenseMtx_colCopy(mtxX, 0, vecX, 0);
/*
      DenseMtx_zero(vecT) ;
      InpMtx_nonsym_mmm(mtxA, vecT, one, mtxX) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      default :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }

      DenseMtx_sub(vecT, mtxB) ;
      Rtmp = DenseMtx_twoNormOfColumn(vecT,0);
      fprintf(msgFile, "\n TFQMRL Residual norm: %6.2e ", Rtmp) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU  : Converges in time: %8.3f ", t2 - t1) ;
      fprintf(msgFile, "\n # iterations = %d", Imv) ;

      fprintf(msgFile, "\n\n after TFQMRL") ;  
      goto end;
    };



    if (Rho == 0){
      fprintf(msgFile, "\n\n Fatal Error, \n"
	      "  TFQMRL Breakdown, Rho = 0 !!") ;
      Imv = -1;
      goto end;
    };

    DenseMtx_colDotProduct(vecW, 0, vecR, 0, &Rho_new);
    Beta    = Rho_new / Rho;
    Rho     = Rho_new ;

    DenseMtx_colCopy(vecY1, 0, vecY2, 0);
    DenseMtx_colGenAxpy(&Beta, vecY1, 0, one, vecW, 0);
/*
    DenseMtx_zero(vecT) ;
    InpMtx_nonsym_mmm(mtxA, vecT, one, vecY1) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecY1) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecY1) ;
	break ;
      default :
	fprintf(msgFile, "\n TFQMRL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }

    FrontMtx_solve(Precond, vecU1, vecT, Precond->manager,
		   cpus, msglvl, msgFile) ;
    Imv++;

/*                                                         */
    DenseMtx_colCopy(vecT, 0, vecU2, 0);
    DenseMtx_colGenAxpy(one, vecT, 0, &Beta, vecV, 0);
    DenseMtx_colCopy(vecV, 0, vecT, 0);
    DenseMtx_colGenAxpy(&Beta, vecV, 0, one, vecU1, 0);

    Rtmp = Tau*sqrt(m + 1)/Init_norm ;
    fprintf(msgFile, "\n\n At iteration %d"
	    "  the convergence ratio is  %12.4e", 
	    Imv, Rtmp) ;

  }
/*            End of while loop              */
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU  : Total iteration time is : %8.3f ", t2 - t1) ;
fprintf(msgFile, "\n # iterations = %d", Imv) ;
fprintf(msgFile, "\n\n  TFQMRL did not Converge !") ;

fprintf(msgFile, "\n\n after TFQMRL") ;
DenseMtx_colCopy(mtxX, 0, vecX, 0);

/*
 
   ------------------------
   free the working storage
   ------------------------
*/
 end:
DenseMtx_free(vecD) ;
DenseMtx_free(vecR) ;
DenseMtx_free(vecT) ;
DenseMtx_free(vecU1) ;
DenseMtx_free(vecU2) ;
DenseMtx_free(vecV) ;
DenseMtx_free(vecW) ;
DenseMtx_free(vecX) ;
DenseMtx_free(vecY1) ;
DenseMtx_free(vecY2) ;

fprintf(msgFile, "\n") ;

return(Imv) ; }

