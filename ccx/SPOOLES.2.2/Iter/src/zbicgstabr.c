/*  zbicgstabr.c  */

#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to solve a complex matrix equation

               Ax=b
   using right preconditioned BiCGSTABR 

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

      return value  -- error flag

   created -- Dec. 10, 1998         Wei-Pai Tang
   ---------------------------------------------------------------------
*/

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
 )
{
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxZ ;
DenseMtx        *vecP, *vecR, *vecR0, *vecT,  *vecV, *vecW;
DenseMtx        *vecX, *vecY, *vecZ ;
double          Alpha[2], Beta[2], Init_norm, Omega[2], Rtmp[2];
double          ratio, Res_norm, Rho[2], Rho_new[2], Ttmp[2];
double          t1, t2,  cpus[9] ;
double          one[2] = {1.0, 0.0}, zero[2] = {0.0, 0.0} ;
double          Tiny = 0.1e-28;
int             Iter, Imv, neqns;
int             stats[6] ;



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

vecR0 = DenseMtx_new() ;
DenseMtx_init(vecR0, type, 0, 0, neqns, 1, 1, neqns) ;

vecT = DenseMtx_new() ;
DenseMtx_init(vecT, type, 0, 0, neqns, 1, 1, neqns) ;

vecV = DenseMtx_new() ;
DenseMtx_init(vecV, type, 0, 0, neqns, 1, 1, neqns) ;

vecW = DenseMtx_new() ;
DenseMtx_init(vecW, type, 0, 0, neqns, 1, 1, neqns) ;

vecX = DenseMtx_new() ;
DenseMtx_init(vecX, type, 0, 0, neqns, 1, 1, neqns) ;

vecY = DenseMtx_new() ;
DenseMtx_init(vecY, type, 0, 0, neqns, 1, 1, neqns) ;

vecZ = DenseMtx_new() ;
DenseMtx_init(vecZ, type, 0, 0, neqns, 1, 1, neqns) ;



/*
   --------------------------
   Initialize the iterations
   --------------------------
*/
Init_norm = DenseMtx_twoNormOfColumn(mtxB, 0);

ratio = 1.0;
DenseMtx_zero(vecX) ;


DenseMtx_colCopy (vecR0, 0, mtxB, 0);
DenseMtx_colCopy (vecR, 0, vecR0, 0);
Iter = 0;
Imv  = 0;

  fprintf(msgFile, "\n ZBiCGSTABR Initial norml: %6.2ef ", Init_norm ) ;
  fprintf(msgFile, "\n ZBiCGSTABR Conveg. Control: %12.8f ", convergetol ) ;

/*

*/
Rho[0]      = 1.0;
Alpha[0]    = 1.0;
Omega[0]    = 1.0;
Rho[1]      = 0.0;
Alpha[1]    = 0.0;
Omega[1]    = 0.0;

DenseMtx_zero(vecV) ;
DenseMtx_zero(vecP) ;
/*
   ------------------------------
   
   ------------------------------
*/

MARKTIME(t1) ;

/*
   -----------------
   factor the matrix
   -----------------
*/
while ( ratio > convergetol && Iter <= itermax )
  {
    Iter++;
/*      Rho_new = DenseMtx_dot(vecR0, vecR);   */
    DenseMtx_colDotProduct (vecR0, 0, vecR,0, Rho_new);
    if ( zabs(Rho_new) == 0.0 || zabs(Omega)== 0.0 ){
      fprintf(stderr, "\n   breakdown in ZBiCGSTABR !! "
              "\n Fatal error   \n");
      exit(-1) ; 
    }

/*      Beta    = (Rho_new / (Rho+Tiny)) * (Alpha / (Omega+Tiny));  */
    zdiv(Rho_new, Rho, Beta);
    zmul(Beta, Alpha, Rtmp);
    zdiv(Rtmp, Omega, Beta);

    Rho[0]     = Rho_new[0];
    Rho[1]     = Rho_new[1];
/*      DenseMtx_axpy(vecP, vecV, -Omega);   */
/*      DenseMtx_aypx(vecP, vecR, Beta);     */
    zsub(zero, Omega, Rtmp);
    DenseMtx_colGenAxpy (one,  vecP, 0, Rtmp, vecV, 0 );
    DenseMtx_colGenAxpy (Beta, vecP, 0, one,  vecR, 0 );
/*                                                         */
    FrontMtx_solve(Precond, vecY, vecP, Precond->manager,
               cpus, msglvl, msgFile) ;
/*                                                         */
/*     InpMtx_nonsym_gmmm(mtxA, zero, vecV, one, vecY) ;   */
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecV, one, vecY) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecV, one, vecY) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecV, one, vecY) ;
	break ;
      default :
	fprintf(msgFile, "\n BiCGSTABR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
    Imv++;

/*   Alpha = Rho / (DenseMtx_dot(vecR0, vecV)+Tiny);           */
    
    DenseMtx_colDotProduct (vecR0, 0, vecV,0, Rtmp);
    if ( zabs(Rtmp) == 0.0 ){
      fprintf(stderr, "\n   breakdown in ZBiCGSTABR !! "
              "\n Fatal error   \n");
      exit(-1) ; 
    }
    zdiv(Rho, Rtmp, Alpha);
    
/*                                                         */
/*    DenseMtx_axpy(vecR, vecV, -Alpha);                   */
    zsub(zero, Alpha, Rtmp);
    DenseMtx_colGenAxpy (one, vecR, 0, Rtmp,  vecV, 0 );

/*                                                         */
    FrontMtx_solve(Precond, vecZ, vecR, Precond->manager,
               cpus, msglvl, msgFile) ;

/*     InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecZ) ;  */
      switch ( symmetryflag ) {
      case SPOOLES_SYMMETRIC : 
	InpMtx_sym_gmmm(mtxA, zero, vecT, one, vecZ) ;
	break ;
      case SPOOLES_HERMITIAN :
	InpMtx_herm_gmmm(mtxA, zero, vecT, one, vecZ) ;
	break ;
      case SPOOLES_NONSYMMETRIC :
	InpMtx_nonsym_gmmm(mtxA, zero, vecT, one, vecZ) ;
	break ;
      default :
	fprintf(msgFile, "\n BiCGSTABR Matrix type wrong");
	fprintf(msgFile, "\n Fatal error");
	goto end;
      }
    Imv++;

/*    Omega = DenseMtx_dot(vecT, vecR)/(DenseMtx_dot(vecT, vecT)+Tiny);  */

    DenseMtx_colDotProduct (vecT, 0, vecT,0, Ttmp);
    if ( zabs(Ttmp) == 0.0 ){
      fprintf(stderr, "\n   breakdown in ZBiCGSTABR !! "
              "\n Fatal error   \n");
      exit(-1) ; 
    };
    DenseMtx_colDotProduct (vecT, 0, vecR, 0, Rtmp);
    zdiv(Rtmp, Ttmp, Omega);


    DenseMtx_colGenAxpy (one, vecX, 0, Alpha, vecY, 0);
    DenseMtx_colGenAxpy (one, vecX, 0, Omega, vecZ, 0);

    /*       DenseMtx_axpy(vecR, vecT, -Omega);                    */
    zsub(zero, Omega, Rtmp);
    DenseMtx_colGenAxpy (one, vecR, 0, Rtmp, vecT, 0);

    Res_norm =  DenseMtx_twoNormOfColumn (vecR, 0);
    ratio = Res_norm/Init_norm;
    fprintf(msgFile, "\n\n At iteration %d"
	    "  the convergence ratio is  %12.4e"
	    "\n Residual norm is %6.2e", 
	    Imv, ratio, Res_norm) ;
  }
/*            End of while loop              */
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU  : Converges in time: %8.3f ", t2 - t1) ;
fprintf(msgFile, "\n # iterations = %d", Imv) ;

fprintf(msgFile, "\n\n after ZBICGSTABR") ;
/*  DenseMtx_copy(mtxX, vecX);  */
DenseMtx_colCopy (mtxX, 0, vecX, 0);
/*
 
   ------------------------
   free the working storage
   ------------------------
*/
 end:
DenseMtx_free(vecP) ;
DenseMtx_free(vecR) ;
DenseMtx_free(vecR0) ;
DenseMtx_free(vecT) ;
DenseMtx_free(vecV) ;
DenseMtx_free(vecW) ;
DenseMtx_free(vecX) ;
DenseMtx_free(vecY) ;
DenseMtx_free(vecZ) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
