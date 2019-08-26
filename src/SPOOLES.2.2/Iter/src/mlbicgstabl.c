/*  MLBiCGSTABL.c  */


#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to solve a unsymmetric real matrix equation

               Ax=b

   using left preconditioned ML(k)BiCGSTAB method 
   by M-C Yeung and T. Chan
 

      x       -- Initial guess as zeros
      A       -- Input matrix
      Precond -- Front Matrix as the preconditioner
      Q       -- Starting vectors
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

   created -- Nov. 28, 1998              Wei-Pai Tang
   ---------------------------------------------------------------------
*/

int
mlbicgstabl (
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
 )
{
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxD, *mtxG, *mtxU, *mtxW;
DenseMtx        *vecGt, *vecR, *vecT, *vecT1;
DenseMtx        *vecX, *vecZD, *vecZG, *vecZW ;
double          Alpha, Beta, Rho ;
double          c[100] ;
double          Init_norm,  ratio,  Res_norm;
double          error_trol, m, Rtmp;
double          t1, t2,  cpus[9] ;
double          one[2] = {1.0, 0.0} ;
double          zero[2] = {0.0, 0.0} ;
double          minusone[2] = {-1.0, 0.0};
double          Tiny = 0.1e-28;
int             Iter, Imv, neqns, Ik, ii, is;
int             stats[6] ;
int             return_flag;




neqns        = n_matrixSize;
Ik           = mtxQ->ncol;

if (Ik > 100){
      fprintf(msgFile, "\n\n Fatal Error, \n"
                    " Too many starting vectors in Q !!") ;
      return_flag = -1;
      return (return_flag);
    };

return_flag  = 1;

/*
   --------------------
   init the vectors in MLBiCGSTABL
   --------------------
*/
mtxD = DenseMtx_new() ;
DenseMtx_init(mtxD, type, 0, 0, neqns, Ik, 1, neqns) ;

mtxG = DenseMtx_new() ;
DenseMtx_init(mtxG, type, 0, 0, neqns, Ik, 1, neqns) ;

mtxU = DenseMtx_new() ;
DenseMtx_init(mtxU, type, 0, 0, neqns, 2, 1, neqns) ;

mtxW = DenseMtx_new() ;
DenseMtx_init(mtxW, type, 0, 0, neqns, Ik, 1, neqns) ;

vecGt = DenseMtx_new() ;
DenseMtx_init(vecGt, type, 0, 0, neqns, 1, 1, neqns) ;

vecR = DenseMtx_new() ;
DenseMtx_init(vecR, type, 0, 0, neqns, 1, 1, neqns) ;

vecT = DenseMtx_new() ;
DenseMtx_init(vecT, type, 0, 0, neqns, 1, 1, neqns) ;

vecT1 = DenseMtx_new() ;
DenseMtx_init(vecT1, type, 0, 0, neqns, 1, 1, neqns) ;

vecX = DenseMtx_new() ;
DenseMtx_init(vecX, type, 0, 0, neqns, 1, 1, neqns) ;


vecZD = DenseMtx_new() ;
DenseMtx_init(vecZD, type, 0, 0, neqns, 1, 1, neqns) ;

vecZG = DenseMtx_new() ;
DenseMtx_init(vecZG, type, 0, 0, neqns, 1, 1, neqns) ;

vecZW = DenseMtx_new() ;
DenseMtx_init(vecZW, type, 0, 0, neqns, 1, 1, neqns) ;

/*
c     = DV_new();
DV_init(c, Ik, NULL);
*/
for ( ii = 0; ii <100; ii++){
  c[ii] = 0;
}
/*
   --------------------------
   Initialize the iterations
   --------------------------
*/
          /*       ----     Set initial guess as zero  ----       */
DenseMtx_zero(vecX) ;

          /*       ----          If x_0 is not  zero     ----     */
/*
DenseMtx_colCopy(vecX, 0, mtxX, 0);           
*/

/*
InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecX) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecX) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecX) ;
	break ;
      default :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      }

DenseMtx_colCopy(vecR, 0, mtxB, 0);
DenseMtx_sub(vecR, vecT) ;

FrontMtx_solve(Precond, vecT, vecR, Precond->manager,
	       cpus, msglvl, msgFile) ;
DenseMtx_colCopy(vecR, 0, vecT, 0);

Init_norm = DenseMtx_twoNormOfColumn(vecR, 0);
if ( Init_norm == 0.0 ){
  Init_norm = 1.0; 
};
error_trol = Init_norm * convergetol ;

  fprintf(msgFile, "\n MLBiCGSTABL Initial norml: %6.2e ", Init_norm ) ;
  fprintf(msgFile, "\n MLBiCGSTABL Conveg. Control: %7.3e ", convergetol ) ;
  fprintf(msgFile, "\n MLBiCGSTABL Convergen Control: %7.3e ",error_trol ) ;

DenseMtx_zero(mtxG) ;
DenseMtx_zero(mtxD) ;
DenseMtx_zero(mtxW) ;


Iter = 0;
Imv  = 0;

    
DenseMtx_colCopy (mtxG, Ik-1, vecR, 0); 

/*
   ------------------------------
   MLBiCGSTABL   Iteration start
   ------------------------------
*/

MARKTIME(t1) ;


while (  Iter <= itermax ){

    Iter++;
/*
    g_tld = G(:,k);
    W(:,k) = U\(L\( A*g_tld));
*/ 
    DenseMtx_colCopy (vecGt, 0, mtxG, Ik-1);
/*
    InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecGt) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecGt) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecGt) ;
	break ;
      default :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      }

    FrontMtx_solveOneColumn(Precond, mtxW, Ik-1, vecT, 0,
			    Precond->manager, cpus, msglvl, msgFile) ;
    Imv++;

/*    c[Ik] =  *DenseMtx_colDotProduct (mtxQ, 0, mtxG, Ik-1); */
    DenseMtx_colDotProduct (mtxQ, 0, mtxG, Ik-1, c+Ik);

    if (c[Ik] == 0){
      fprintf(msgFile, "\n\n Fatal Error, \n"
	      "  MLBiCGSTABL Breakdown, c[k] = 0 !!") ;
      return_flag   = -1;
      goto end;
    };

/*    Alpha = Q(:,1)'*r/c(k); ; */
    
    DenseMtx_colDotProduct (mtxQ,0, vecR, 0, &Rtmp);
    Alpha = Rtmp/c[Ik];
    DenseMtx_colCopy (mtxU, 0, vecR, 0);
    Rtmp = -Alpha;
    DenseMtx_colGenAxpy (one, mtxU, 0, &Rtmp, mtxW, Ik-1);
/*
    u = r - alpha*W(:,k);
    u_tld = u;
    temp = U\(L\(A*u_tld));
 
    rho =temp'*temp;
*/
    DenseMtx_colCopy (mtxU, 1, mtxU, 0);
    DenseMtx_colCopy (vecT, 0, mtxU, 0);
/*
    InpMtx_nonsym_gmmm(mtxA, zero, vecT1, one, vecT) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT1, one, vecT) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT1, one, vecT) ;
	break ;
      default :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      }

    FrontMtx_solveOneColumn(Precond, vecT, 0, vecT1, 0,
			    Precond->manager, cpus, msglvl, msgFile) ;
    Imv++;

/*    Rho   = temp'*temp;  */
    DenseMtx_colDotProduct (vecT, 0, vecT, 0, &Rho);
    if (Rho == 0){
      fprintf(msgFile, "\n\n Fatal Error, \n"
	      "  MLBiCGSTABL Breakdown, Rho = 0 !!") ;
      return_flag   = -1;
      goto end;
    };

/*    Rho   = -(*DenseMtx_colDotProduct (mtxU, 0, vecT, 0))/Rho; */
    DenseMtx_colDotProduct (mtxU, 0, vecT, 0, &Rtmp);
    Rho = -Rtmp/Rho;

    DenseMtx_colGenAxpy(one, vecX, 0, &Alpha, vecGt, 0);
    Rtmp = -Rho;
    DenseMtx_colGenAxpy (one, vecX, 0, &Rtmp, mtxU, 1) ;
    DenseMtx_colCopy (vecR, 0, mtxU, 0);
    DenseMtx_colGenAxpy (one, vecR, 0, &Rho, vecT, 0) ;
/*     Iter++;    */

/*
    ----------------
    Convergence Test
    ---------------
*/
    Rtmp  =  DenseMtx_twoNormOfColumn ( vecR, 0);
    ratio =  Rtmp/Init_norm;
    fprintf(msgFile, "\n\n At iteration %d"
	    "  the convergence ratio is  %12.4e"
	    "\n Residual norm is %6.2ee", 
	    Imv, ratio, Rtmp) ;
    fflush(msgFile) ;
    if ( Rtmp   <= error_trol ) {
      fprintf(msgFile, "\n # iterations = %d", Imv) ;
      return_flag = Imv;
      goto end;
    };

    for (ii = 1; ii < Ik+1; ii++ ){
      if (Iter > itermax ){
	fprintf(msgFile, "\n # iterations = %d", Imv) ;
        fprintf(msgFile, "\n\n  MLBiCGSTABL did not Converge !") ;
	return_flag = Imv;
	goto end;
      };

      
      DenseMtx_colCopy (vecZD,0, mtxU, 0);
      DenseMtx_colCopy (vecZG,0, vecR, 0);
      DenseMtx_zero(vecZW) ;
     
      if ( Iter > 1 ){
	for ( is = ii ; is < Ik ; is++ ){
/*	  Beta =-( *(DenseMtx_colDotProduct(mtxQ, is, vecZD, 0)))/c[is]; */
	  DenseMtx_colDotProduct(mtxQ, is, vecZD, 0, &Rtmp);
	  Beta = -Rtmp/c[is];
	  DenseMtx_colGenAxpy (one, vecZD, 0, &Beta, mtxD, is-1);
	  DenseMtx_colGenAxpy (one, vecZG, 0, &Beta, mtxG, is-1);
	  DenseMtx_colGenAxpy (one, vecZW, 0, &Beta, mtxW, is-1);
	};
      };

      Beta = Rho * c[Ik];

      if (Beta == 0){
	fprintf(msgFile, "\n\n Fatal Error, \n"
		"  MLBiCGSTABL Breakdown, Beta = 0 !!") ;
	return_flag   = -1;
	goto end;
      };
/*
         Beta = - Q(:,1)'* (r + Rho* zw )/ Beta;
*/
      DenseMtx_colCopy (vecT, 0, vecR, 0);
      DenseMtx_colGenAxpy (one, vecT, 0, &Rho, vecZW, 0);
/*      Beta = - (*DenseMtx_colDotProduct(mtxQ, 0, vecT, 0))/Beta; */
      DenseMtx_colDotProduct(mtxQ, 0, vecT, 0, &Rtmp);
      Beta = - Rtmp/Beta;
/*
      zg = zg + beta*G(:,k);
      zw = rho*(zw + beta*W(:,k));
      zd = r + zw;
*/
      DenseMtx_colGenAxpy (one, vecZG, 0, &Beta, mtxG, Ik-1);
      Rtmp = Rho*Beta;
      DenseMtx_colGenAxpy (&Rho, vecZW, 0, &Rtmp, mtxW, Ik-1);
      DenseMtx_colCopy (vecZD, 0, vecR, 0);
      DenseMtx_colGenAxpy (one, vecZD, 0, one, vecZW, 0);
/*
      
       for s = 1:i-1
          beta = -Q(:,s+1)'*zd/c(s);
          zd = zd + beta*D(:,s);
          zg = zg + beta*G(:,s);
       end
*/
     for ( is = 1; is < ii - 1; is ++){
/*	Beta =  -(*DenseMtx_colDotProduct(mtxQ, is, vecZD, 0))/c[is] ; */
        DenseMtx_colDotProduct(mtxQ, is, vecZD, 0, &Rtmp);
	Beta = - Rtmp/c[is];
	DenseMtx_colGenAxpy (one, vecZD, 0, &Beta, mtxD, is-1);
	DenseMtx_colGenAxpy (one, vecZG, 0, &Beta, mtxG, is-1);
      };

/*
      D(:,i) = zd - u;
      G(:,i) = zg + zw;
*/
      DenseMtx_colCopy (mtxD, ii-1, vecZD, 0);
      DenseMtx_colGenAxpy (one, mtxD, ii-1, minusone, mtxU, 0);
      DenseMtx_colCopy (mtxG, ii-1, vecZG, 0);
      DenseMtx_colGenAxpy (one, mtxG, ii-1, one, vecZW, 0);

/*
      if i < k
        c(i) = Q(:,i+1)'*D(:,i);
*/
      if ( ii < Ik ){
/*	c[ii] = *DenseMtx_colDotProduct(mtxQ, ii, mtxD, ii-1); */
	DenseMtx_colDotProduct(mtxQ, ii, mtxD, ii-1, c+ii);
/*
                            If breakdown ?
*/
	if (c[ii] == 0){
	  fprintf(msgFile, "\n\n Fatal Error, \n"
		  "  MLBiCGSTABL Breakdown, c[ii] = 0 !!") ;
	  return_flag   = -1;
	  goto end;
	};
/*
        alpha = Q(:,i+1)'*u/c(i);
        u = u - alpha*D(:,i);
        g_tld = G(:,i);
*/
	DenseMtx_colDotProduct(mtxQ, ii, mtxU, 0, &Rtmp);
	Alpha = Rtmp/c[ii];
	Rtmp = -Alpha;
	DenseMtx_colGenAxpy (one, mtxU, 0, &Rtmp, mtxD, ii-1);
	DenseMtx_colCopy (vecGt, 0, mtxG, ii-1);
/*
        x = x + rho*alpha*g_tld;
        W(:,i) = U\(L\(A*g_tld));
        r = r - rho*alpha*W(:,i);
*/
	Rtmp = Rho * Alpha;
	DenseMtx_colGenAxpy (one, vecX, 0, &Rtmp, vecGt, 0);
/*
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecGt) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecGt) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecGt) ;
	break ;
      default :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      }

	FrontMtx_solveOneColumn(Precond, mtxW, ii-1, vecT, 0,
				Precond->manager, cpus, msglvl, msgFile) ;
	Imv++;
	Rtmp = -Rtmp;
 	DenseMtx_colGenAxpy (one, vecR, 0, &Rtmp, mtxW, ii-1);

/*
    ----------------
    Convergence Test
    ---------------
*/
	Rtmp =  DenseMtx_twoNormOfColumn ( vecR, 0);
	if ( Rtmp   <= error_trol ) {
	  fprintf(msgFile, "\n # iterations = %d", Imv) ;
	  return_flag = Imv;
	  goto end;
	};
      };

    };

}
/*            End of while loop              */
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU  : Total iteration time is : %8.3f ", t2 - t1) ;
fprintf(msgFile, "\n # iterations = %d", Imv) ;
fprintf(msgFile, "\n\n  MLBiCGSTABL did not Converge !") ;
DenseMtx_colCopy(mtxX, 0, vecX, 0);

/*
 
   ------------------------
   free the working storage
   ------------------------
*/
 end:
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU  : Total iteration time is : %8.3f ", t2 - t1) ;
DenseMtx_colCopy(mtxX, 0, vecX, 0);
/*
InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
*/
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      case SPOOLES_HERMITIAN :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, mtxX) ;
	break ;
      default :
	fprintf(msgFile, "\n MLBiCGSTABL Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	return_flag = -1;
	goto end;
      }

DenseMtx_sub(vecT, mtxB) ;
Rtmp = DenseMtx_twoNormOfColumn(vecT, 0);
fprintf(msgFile, "\n MLBiCGSTABL True Residual norm: %6.2e ", 
	Rtmp) ;
fprintf(msgFile, "\n\n after MLBiCGSTABL") ;
DenseMtx_free(mtxD) ;
DenseMtx_free(mtxG) ;
DenseMtx_free(mtxU) ;
DenseMtx_free(mtxW) ;
DenseMtx_free(vecGt) ;
DenseMtx_free(vecR) ;
DenseMtx_free(vecT) ;
DenseMtx_free(vecT1) ;
DenseMtx_free(vecX) ;
DenseMtx_free(vecZD) ;
DenseMtx_free(vecZG) ;
DenseMtx_free(vecZW) ;


fprintf(msgFile, "\n") ;

return(return_flag) ; }

/*--------------------------------------------------------------------*/
