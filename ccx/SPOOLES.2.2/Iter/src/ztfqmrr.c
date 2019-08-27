/*  ZTFQMRR.c  */

#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to solve a  complex matrix equation

               Ax=b

   using right preconditioned TFQMR method without lookahead 

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

   created -- Dec. 10, 1998              Wei-Pai Tang
   ---------------------------------------------------------------------
*/

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
 )
{
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *vecD, *vecR, *vecT, *vecU1, *vecU2,  *vecV, *vecW;
DenseMtx        *vecX, *vecY1, *vecY2 ;
double          Alpha[2], Beta[2], Cee, Eta[2], Rho[2], Rho_new[2] ;
double          Sigma[2], Tau, Theta, Rtmp[2], Ttmp[2];
double          Init_norm,  ratio,  Res_norm;
double          error_trol, m ;
double          t1, t2,  cpus[9] ;
double          one[2] = {1.0, 0.0}, zero[2] = {0.0, 0.0} ;
double          Tiny = 0.1e-28;
int             Iter, Imv, neqns;
int             stats[6] ;



neqns = n_matrixSize;


/*
   --------------------
   init the vectors in ZTFQMRR
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

Init_norm = DenseMtx_twoNormOfColumn(mtxB, 0);
if ( Init_norm == 0.0 ){
  Init_norm = 1.0; 
};
error_trol = Init_norm * convergetol ;

  fprintf(msgFile, "\n ZTFQMRR Initial norml: %6.2e ", Init_norm ) ;
  fprintf(msgFile, "\n ZTFQMRR Conveg. Control: %6.2e ", convergetol ) ;
  fprintf(msgFile, "\n ZTFQMRR Convergen Control: %6.2e ",error_trol ) ;

DenseMtx_zero(vecD) ;
DenseMtx_zero(vecU1) ;
DenseMtx_zero(vecU2) ;
DenseMtx_zero(vecY2) ;

Iter = 0;
Imv  = 0;

DenseMtx_colCopy (vecR, 0, mtxB, 0);
DenseMtx_colCopy (vecW, 0, vecR, 0);
DenseMtx_colCopy (vecY1, 0, vecR, 0);


/*                                                         */
    FrontMtx_solve(Precond, vecT, vecY1, Precond->manager,
               cpus, msglvl, msgFile) ;

      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecV, one, vecT) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecV, one, vecT) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecV, one, vecT) ;
	break ;
      default :
	fprintf(msgFile, "\n BiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }


    Imv++;

    DenseMtx_colCopy (vecU1, 0, vecV, 0);
/*

*/
Eta[0]     = 0.0;
Eta[1]     = 0.0;
Tau     = Init_norm ;
Theta   = 0.0;
/*    Rho     = Tau * Tau ;    */
Rho[0]  = Tau * Tau;
Rho[1]  = 0.0;



/*
   ------------------------------
   ZTFQMRR   Iteration start
   ------------------------------
*/

MARKTIME(t1) ;


while (  Iter <= itermax )
  {
    Iter++;
/*       Sigma   = DenseMtx_dot(vecV, vecR);    */
    DenseMtx_colDotProduct (vecV, 0, vecR,0, Sigma);

    if (zabs(Sigma) == 0){
      fprintf(msgFile, "\n\n Fatal Error, \n"
	      "  ZTFQMRR Breakdown, Sigma = 0 !!") ;
      Imv = -1;
      goto end;
    };

/*          Alpha   = Rho/Sigma;    */
    zdiv(Rho, Sigma, Alpha);

/*
    ----------------
    Odd step
    ---------------
*/
	
    m      = 2 * Iter - 1;

    zsub(zero, Alpha, Rtmp);
    DenseMtx_colGenAxpy (one, vecW, 0, Rtmp,  vecU1, 0 );

/*       Rtmp   = Theta * Theta * Eta / Alpha ;  */
    Rtmp[0] = Theta * Theta;
    Rtmp[1] = 0.0;
    zmul(Rtmp, Eta, Ttmp);
    zdiv(Ttmp, Alpha, Rtmp);

    DenseMtx_colGenAxpy (Rtmp, vecD, 0, one,  vecY1, 0 );

/*       Theta  = DenseMtx_fnorm(vecW)/Tau;     */
    Theta =  DenseMtx_twoNormOfColumn(vecW, 0)/Tau;


    Cee    = 1.0/sqrt(1.0 + Theta*Theta);
    Tau    = Tau * Theta * Cee ;
/*       Eta    = Cee * Cee * Alpha ;    */
    Rtmp[0] = Cee * Cee;
    Rtmp[1] = 0.0;
    zmul(Rtmp, Alpha, Eta);

    DenseMtx_colGenAxpy (one, vecX, 0, Eta,  vecD, 0 );

      fprintf(msgFile, "\n\n Odd step at %d", Imv);
      fprintf(msgFile, " \n Tau is   : %6.2e", Tau) ; 
/*   
   ------------------                
     For debug purpose: Check the true residual
   ------------------  
 */
/*
   ---------------
      DenseMtx_colCopy (vecT, 0, vecX, 0);
      FrontMtx_solve(Precond, mtxX, vecT, Precond->manager,
		     cpus, msglvl, msgFile) ;
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      default :
	fprintf(msgFile, "\n ZTFQMRR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
      DenseMtx_sub(vecT, mtxB) ;
      Rtmp[0] = DenseMtx_fnorm(vecT);
      fprintf(msgFile, "\n ZTFQMRR Residual norm: %6.2e ", Rtmp[0]) ;
   ---------------
*/
 
/*
    ----------------
    Convergence Test
    ---------------
*/
    if (Tau * sqrt(m + 1)  <= error_trol ) {
/*                                                             */
      FrontMtx_solve(Precond, mtxX, vecX, Precond->manager,
		     cpus, msglvl, msgFile) ;
/*                                                             */
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      default :
	fprintf(msgFile, "\n ZTFQMRR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
    DenseMtx_sub(vecT, mtxB) ;
    Rtmp[0]  = DenseMtx_twoNormOfColumn(vecT, 0);
    fprintf(msgFile, "\n ZTFQMRR Residual norm: %6.2e ", Rtmp[0]) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU  : Converges in time: %8.3f ", t2 - t1) ;
      fprintf(msgFile, "\n # iterations = %d", Imv) ;
      fprintf(msgFile, "\n\n after ZTFQMRR") ;  
      goto end;
    };

/*
    ----------------
    Even step
    ---------------
*/

    DenseMtx_colCopy (vecY2, 0, vecY1, 0);
    zsub(zero, Alpha, Rtmp);
    DenseMtx_colGenAxpy (one, vecY2, 0, Rtmp,  vecV, 0 );

    FrontMtx_solve(Precond, vecT, vecY2, Precond->manager,
		   cpus, msglvl, msgFile) ;

      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecU2, one, vecT) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecU2, one, vecT) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecU2, one, vecT) ;
	break ;
      default :
	fprintf(msgFile, "\n ZTFQMRR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }

    m      = 2 * Iter ;

    zsub(zero, Alpha, Rtmp);
    DenseMtx_colGenAxpy (one, vecW, 0, Rtmp,  vecU2, 0 );

/*     
    Rtmp   = Theta * Theta * Eta / Alpha ; 
*/
    Rtmp[0] = Theta * Theta;
    Rtmp[1] = 0.0;
    zmul(Rtmp, Eta, Ttmp);
    zdiv(Ttmp, Alpha, Rtmp);
    DenseMtx_colGenAxpy (Rtmp, vecD, 0, one,  vecY2, 0 );

/*      Theta  = DenseMtx_fnorm(vecW)/Tau;    */
    Theta =  DenseMtx_twoNormOfColumn(vecW, 0)/Tau;
   
    Cee    = 1.0/sqrt(1.0 + Theta*Theta);
    Tau    = Tau * Theta * Cee ;
/*       Eta    = Cee * Cee * Alpha ;    */
    Rtmp[0] = Cee * Cee;
    Rtmp[1] = 0.0;
    zmul(Rtmp, Alpha, Eta);

    DenseMtx_colGenAxpy (one, vecX, 0, Eta,  vecD, 0 );
    
      fprintf(msgFile, "\n\n Even step at %d", Imv) ;  
    
/*
    ----------------
    Convergence Test for even step
    ---------------
*/
    if (Tau * sqrt(m + 1)  <= error_trol ) {
/*                                                    */    
      FrontMtx_solve(Precond, mtxX, vecX, Precond->manager,
		     cpus, msglvl, msgFile) ;

      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      default :
	fprintf(msgFile, "\n ZTFQMRR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }


    DenseMtx_sub(vecT, mtxB) ;
    Rtmp[0] = DenseMtx_twoNormOfColumn(vecT, 0);
    fprintf(msgFile, "\n ZTFQMRR Residual norm: %6.2e ", Rtmp[0]) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU  : Converges in time: %8.3f ", t2 - t1) ;
      fprintf(msgFile, "\n # iterations = %d", Imv) ;

      fprintf(msgFile, "\n\n after ZTFQMRR") ; 
      goto end;
    };



    if (zabs(Rho) == 0){
      fprintf(msgFile, "\n\n Fatal Error, \n"
	      "  ZTFQMRR Breakdown, Rho = 0 !!") ;
      Imv = -1;
      goto end;
    };

/*
    Rho_new = DenseMtx_dot(vecW, vecR);
    Beta    = Rho_new / Rho;
    Rho     = Rho_new ;
*/
    DenseMtx_colDotProduct (vecW, 0, vecR,0, Rho_new);
    zdiv(Rho_new, Rho, Beta);
    Rho[0]= Rho_new[0];
    Rho[1]= Rho_new[1];

    DenseMtx_colCopy (vecY1, 0, vecY2, 0);
    DenseMtx_colGenAxpy (Beta, vecY1, 0, one,  vecW, 0 );

/*                                                         */
    FrontMtx_solve(Precond, vecT, vecY1, Precond->manager,
		   cpus, msglvl, msgFile) ;

      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecU1, one, vecT) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecU1, one, vecT) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecU1, one, vecT) ;
	break ;
      default :
	fprintf(msgFile, "\n ZTFQMRR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }


    Imv++;

    DenseMtx_colCopy (vecT, 0, vecU2, 0);
    DenseMtx_colGenAxpy (one, vecT, 0, Beta,  vecV, 0 );
    DenseMtx_colCopy (vecV, 0, vecT, 0);
    DenseMtx_colGenAxpy (Beta, vecV, 0, one,  vecU1, 0 );



    Rtmp[0] = Tau*sqrt(m + 1)/Init_norm ;
    fprintf(msgFile, "\n\n At iteration %d"
	    "  the convergence ratio is  %12.4e", 
	    Imv, Rtmp[0]) ;

  }
/*            End of while loop              */
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU  : Total iteration time is : %8.3f ", t2 - t1) ;
fprintf(msgFile, "\n # iterations = %d", Imv) ;
fprintf(msgFile, "\n\n  ZTFQMRR did not Converge !") ;

      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecX) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecT, one, vecX) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecX) ;
	break ;
      default :
	fprintf(msgFile, "\n ZTFQMRR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
      Imv++;
    DenseMtx_sub(vecT, mtxB) ;

    Rtmp[0]  = DenseMtx_twoNormOfColumn(vecT, 0);
    fprintf(msgFile, "\n ZTFQMRR Residual norm: %6.2e ", Rtmp[0] ) ;

fprintf(msgFile, "\n\n after ZTFQMRR") ;

DenseMtx_colCopy (mtxX, 0, vecX, 0);

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

return(1) ; }

/*--------------------------------------------------------------------*/
