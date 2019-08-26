/*  bgmresr.c  */

#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to solve the system A X = B using the block GMRES
   Arnoldi iteration method with right preconditioning

   neqns        -- # of equations
   type         -- must be SPOOLES_REAL or SPOOLES_COMPLEX
   symmetryflag -- must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC
   mtxA         -- pointer to object for A
   mtxM         -- pointer to object for M, may be NULL
   mtxX         -- pointer to object for X, zero on input
   mtxB         -- pointer to object for B
   maxnouter    -- maximum number of outer iterations
   maxninner    -- maximum number of inner iterations
   pninner      -- the last # of inner iterations executed
   pnouter      -- the last # of outer iterations executed
   convergetol  -- tolerance for convergence
      convergence achieved when ||r||_f <= convergetol * ||r_0||_f
   msglvl       -- message level
      1 -- no output
      2 -- iteration #, residual, ratio
      3 -- scalar information
      4 -- scalar and vector information
   msgFile      -- message file

   return value ---
      1 -- normal return
     -1 -- neqns <= 0
     -2 -- type invalid
     -3 -- symmetryflag invalid
     -4 -- mtxA NULL
     -5 -- mtxX NULL
     -6 -- mtxB NULL
     -7 -- convergetol invalid
     -8 -- maxnouter and/or maxninner invalid
     -9 -- pnouter NULL
    -10 -- pninner NULL
    -11 -- msglvl > 0 and msgFile = NULL

   created -- 98dec03, dkw
   ---------------------------------------------------------------------
*/
int
bgmresr (
   int        neqns,
   int        type,
   int        symmetryflag,
   InpMtx     *mtxA,
   FrontMtx   *mtxM,
   DenseMtx   *mtxX,
   DenseMtx   *mtxB,
   int        maxnouter,
   int        maxninner,
   int        *pnouter,
   int        *pninner,
   double     convergetol,
   int        msglvl,
   FILE       *msgFile
) {
A2         *A2G,  *A2Gi, *A2Gj, *A2H, *A2Hij, *A2Hkj, *A2H11, *A2H12,
           *A2Vj, *A2Vk, *A2Z ;
DenseMtx   *G, *Gj, *H, *Hij, *Hkj, *V, *Vi, *Vj, *Vk, *W, *Z;
double     Hkk, Hik, initResNorm, ops, resnorm, t0, t1 ;
double     cpus[9], minusone[2] = {-1.0, 0.0}, one[2] = {1.0, 0.0},
           zero[2] = {0.0, 0.0} ;
double     *col, *val, *Y ;
DV         workDV ;
int        converged, i, ii, iinner, jinner, jouter, kinner,
           nrhs, nstep,
           hijfrow, hijlrow, hijfcol, hijlcol, hkjfcol, hkjlcol,
           iend, inc1, inc2, irow, istart, jcol, ncol, ndiag, nrow,
           rc, vifcol, vilcol, vjfcol, vjlcol, vkfcol, vklcol ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( neqns <= 0 ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n neqns = %d\n", neqns) ;
   return(-1) ;
}
if ( type != SPOOLES_REAL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n type = %d, must be SPOOLES_REAL\n", type) ;
   return(-2) ;
}
if (  symmetryflag != SPOOLES_SYMMETRIC 
   && symmetryflag != SPOOLES_NONSYMMETRIC ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n symmetryflag = %d"
           "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n", 
           symmetryflag) ;
   return(-3) ;
}
if ( mtxA == NULL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n mtxA is NULL\n") ;
   return(-4) ;
}
if ( mtxX == NULL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n mtxX is NULL\n") ;
   return(-5) ;
}
if ( mtxB == NULL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n mtxB is NULL\n") ;
   return(-6) ;
}
if ( convergetol < 0.0 ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n convergetol is %e\n", convergetol) ;
   return(-7) ;
}
if ( maxnouter <= 0 || maxninner <= 0 ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n maxnouter %d, maxninner %d\n", maxnouter, maxninner) ;
   return(-8) ;
}
if ( pnouter == NULL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n pnouter is NULL\n") ;
   return(-9) ;
}
if ( pninner == NULL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n pninner is NULL\n") ;
   return(-10) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n error in bgmresr()"
           "\n msglvl = %d and msgFile is NULL\n", msglvl) ;
   return(-11) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n\n maximum # of outer iterations = %d"
           "\n maximum # of inner iterations = %d",
           maxnouter, maxninner) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   check for zero rhs (if B is zero, then X is zero)
   -------------------------------------------------
*/
initResNorm = DenseMtx_frobNorm(mtxB) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n initial residual = %12.4e", initResNorm) ;
   fflush(msgFile) ;
}
if ( initResNorm == 0.0 ) {
   *pnouter = *pninner = 0 ;
   return(1) ;
}
/*
   ---------------------------
   initialize the working data
   ---------------------------
*/
nrhs = mtxB->ncol ;
Y = DVinit(maxninner+1, 0.0) ;
DV_setDefaultFields(&workDV) ;

H = DenseMtx_new() ;
DenseMtx_init(H, SPOOLES_REAL, 0, 0,
   (maxninner+1)*nrhs, maxninner*nrhs, 1, (maxninner+1)*nrhs) ;
DenseMtx_zero(H) ;

V = DenseMtx_new() ;
DenseMtx_init(V, SPOOLES_REAL, 0, 0, neqns, nrhs*(maxninner+1), 1, neqns) ;
DenseMtx_zero(V) ;

G = DenseMtx_new() ;
DenseMtx_init(G, SPOOLES_REAL, 0, 0,
   (maxninner+1)*nrhs, nrhs, 1, (maxninner+1)*nrhs) ;
DenseMtx_zero(G) ;

Z = DenseMtx_new() ;
DenseMtx_init(Z, SPOOLES_REAL, 0, 0, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(Z) ;

W = DenseMtx_new() ;
DenseMtx_init(W, SPOOLES_REAL, 0, 0, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(W) ;

Vi  = DenseMtx_new() ;
Vj  = DenseMtx_new() ;
Vk  = DenseMtx_new() ;
Gj  = DenseMtx_new() ;
Hij = DenseMtx_new() ;
Hkj = DenseMtx_new() ;
A2G   = A2_new() ; 
A2H   = A2_new() ;
A2Z   = A2_new() ;
DenseMtx_setA2(Z, A2Z) ;
/*
   ---------------------------
   main loop of the iterations
   ---------------------------
*/
converged = 0 ;
for ( jouter = 0 ; jouter < maxnouter ; jouter++ ) {
   A2Gi  = A2_new() ;
   A2Vj  = A2_new() ;
   DenseMtx_zero(H) ;
   DenseMtx_zero(V) ;
   DenseMtx_zero(G) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n outer iteration %d", jouter) ;
      fflush(msgFile) ;
   }
/*
   -----------------------------------------------------
   compute the preconditioned residual and load into V_0
   1. solve M Z = X 
   2. compute Z = A*Z
   -----------------------------------------------------
*/
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n B") ;
      DenseMtx_writeForHumanEye(mtxB, msgFile) ;
      fflush(msgFile) ;
   }
   if ( jouter == 0 ) {
      for ( i = 0 ; i < nrhs ; i++ )
          DenseMtx_colCopy(Z, i, mtxB, i) ;
   }else{
      if ( mtxM != NULL ) {
         FrontMtx_solve(mtxM, W, mtxX, mtxM->manager, 
                        cpus, 0, msgFile) ;
      }
      for ( i = 0 ; i < nrhs ; i++ )
          DenseMtx_colCopy(Z, i, mtxB, i) ;
      if ( symmetryflag == SPOOLES_SYMMETRIC ) {
         InpMtx_sym_mmm(mtxA, Z, minusone, W) ;
      } else {
         InpMtx_nonsym_mmm(mtxA, Z, minusone, W) ;
      }
   }
/*
   ---------------------------------------------------
   Compute the QR factorization of the initial
   residual to get our starting orthogonal vectors V_0
   and our top block in G
   ---------------------------------------------------
*/
   ops = A2_QRreduce(A2Z, &workDV, msglvl, msgFile) ;

   rc = DenseMtx_initAsSubmatrix(Vj, V, 0, neqns-1, 0, nrhs-1) ;
   DenseMtx_setA2(Vj, A2Vj) ;
   A2_computeQ(A2Vj, A2Z, &workDV, msglvl, msgFile) ;

   DenseMtx_initAsSubmatrix(Gj, G, 0, 2*nrhs-1, 0, nrhs-1) ;
   DenseMtx_setA2(Gj, A2Gi) ;

   nrow = A2Z->n1 ;
   ncol = A2Z->n2 ;
   inc1 = A2Z->inc1 ;
   inc2 = A2Z->inc2 ;
   if ( nrow >= ncol ) {
      ndiag = ncol ;
   } else {
      ndiag = nrow ;
   }
   for ( jcol = 0, col = A2Z->entries, val = A2Gi->entries ;
         jcol < ncol ;
         jcol++, col += inc2, val += inc2 ) {
      istart = 0 ;
      iend   = (jcol < ndiag) ? jcol : ndiag - 1 ;
      for ( irow = istart, ii = irow*inc1 ;
            irow <= iend ;
            irow++, ii += inc1 ) {
            val[ii] = col[ii] ;
      }
   }
   resnorm = DenseMtx_frobNorm(G) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n after QRreduce A2Vj =") ;
      fprintf(msgFile, "\n\n V_%d = ", jouter) ;
      A2_writeForHumanEye(A2Vj, msgFile) ;
      fprintf(msgFile, "\n\n G_%d%d = ", jouter, jouter) ;
      A2_writeForHumanEye(A2Gi, msgFile) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n ||G||_f = %12.4e", resnorm) ;
      fflush(msgFile) ;
   }
   if ( resnorm == 0.0 ) {
      break ;
   }
   if ( jouter == 0 ) {
      initResNorm = resnorm ;
   }
/*
   ------------------------------------------------------------
   loop over the inner gmres iterations using Arnoldi iteration
   ------------------------------------------------------------
*/
   nstep = maxninner - 1 ;
   for ( jinner = 0 ; jinner < maxninner ; jinner++ ) {
      A2Gj  = A2_new() ;
      A2Hij = A2_new() ;
      A2Hkj = A2_new() ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n inner iteration %d", jinner) ;
         fflush(msgFile) ;
      }
/*
      ----------------------------
      compute W = A * M^{-1} * V_j
      ----------------------------
*/
      DenseMtx_zero(Z) ;
      vjfcol = jinner*nrhs ;
      vjlcol = vjfcol+nrhs-1 ;
      DenseMtx_initAsSubmatrix(Vj, V, 0, neqns-1, vjfcol, vjlcol) ;
      if ( mtxM != NULL ) {
         FrontMtx_solve(mtxM, W, Vj, mtxM->manager, 
                        cpus, 0, msgFile) ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n ||M^{-1}*V_%d||_f = %12.4e",
                 jinner, DenseMtx_frobNorm(W)) ;
         fflush(msgFile) ;
      }
      if ( symmetryflag == SPOOLES_SYMMETRIC ) {
         InpMtx_sym_mmm(mtxA, Z, one, W) ;
      } else {
         InpMtx_nonsym_mmm(mtxA, Z, one, W) ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n ||A*M^{-1}V_%d||_f = %12.4e",
                 jinner, DenseMtx_frobNorm(Z)) ;
         fflush(msgFile) ;
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n A*M^{-1}*V_%d", jinner) ;
         DenseMtx_writeForHumanEye(Z, msgFile) ;
         fflush(msgFile) ;
      }

      for ( iinner = 0 ; iinner <= jinner ; iinner++ ) {
         vifcol = iinner*nrhs ;
         vilcol = vifcol+nrhs-1 ;
         DenseMtx_initAsSubmatrix(Vi, V, 0, neqns-1, vifcol, vilcol) ;
         hijfrow = vifcol ;
         hijlrow = vilcol ;
         hijfcol = jinner*nrhs ;
         hijlcol = hijfcol+nrhs-1 ;
         DenseMtx_initAsSubmatrix(Hij, H, hijfrow, hijlrow, hijfcol, hijlcol) ;
         DenseMtx_mmm("t", "n", zero, Hij, one, Vi, Z) ;
         DenseMtx_mmm("n", "n", one, Z, minusone, Vi, Hij) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n H_%d%d = ", iinner, jinner) ;
            DenseMtx_writeForHumanEye(Hij, msgFile) ;
            fflush(msgFile) ;
         }
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n Wj after Arnoldi iteration") ;
         DenseMtx_writeForHumanEye(Z, msgFile) ;
         fflush(msgFile) ;
      }
/*
      -----------------------------------------------
      compute V_jinner+1 and H_jinner+1,jinner blocks
      using QR decomposition of W_j then
      copy R from A2Z(1:ncolA,1:ncolA) into H_jinner+1,jinner
      -----------------------------------------------
*/
      ops = A2_QRreduce(A2Z, &workDV, msglvl, msgFile) ;
      vkfcol = vjlcol+1 ;
      vklcol = vkfcol+nrhs-1 ;
      DenseMtx_initAsSubmatrix(Vk, V, 0, neqns-1, vkfcol, vklcol) ;
      A2Vk  = A2_new() ;
      DenseMtx_setA2(Vk, A2Vk) ;
      A2_computeQ(A2Vk, A2Z, &workDV, msglvl, msgFile) ;
      A2_free(A2Vk) ;

      hijfrow = hijlrow+1 ;
      hijlrow = hijfrow+nrhs-1 ;
      DenseMtx_initAsSubmatrix(Hkj, H, hijfrow, hijlrow, hijfcol, hijlcol) ;
      DenseMtx_setA2(Hkj, A2Hkj) ;
      nrow = A2Z->n1 ;
      ncol = A2Z->n2 ;
      inc1 = A2Z->inc1 ;
      inc2 = A2Z->inc2 ;
      if ( nrow >= ncol ) {
         ndiag = ncol ;
      } else {
         ndiag = nrow ;
      }
      for ( jcol = 0, col = A2Z->entries, val = A2Hkj->entries ;
            jcol < ncol ;
            jcol++, col += inc2, val += inc2 ) {
         istart = 0 ;
         iend   = (jcol < ndiag) ? jcol : ndiag - 1 ;
         for ( irow = istart, ii = irow*inc1 ;
               irow <= iend ;
               irow++, ii += inc1 ) {
               val[ii] = col[ii] ;
         }
      }
      if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n V_%d = ", jinner+1) ;
            DenseMtx_writeForHumanEye(Vk, msgFile) ;
            fprintf(msgFile, "\n\n H_%d%d = ", jinner+1, jinner) ;
            DenseMtx_writeForHumanEye(Hij, msgFile) ;
            fflush(msgFile) ;
      }
/*
      -------------------------------------------------
      Before starting the next iteration, triangulate
      the Hessenberg matrix.  To do this we must first
      apply the previous triangulating transformations
      to the newly constructed columns and then
      triangulate the new portion of the Hessenberg
      matrix.  When we compute these new transformations
      we need to also apply them to the corresponding
      portion of G to keep track of the norm of the
      current residual.
 
      First apply the previous transformations
      -------------------------------------------------
*/
      hijfrow = 0 ;
      hijlrow = 2*nrhs-1 ;
      hkjfcol = 0 ;
      hkjlcol = nrhs-1 ;
      for ( iinner = 0 ; iinner <= jinner-1 ; iinner++ ) {
         A2H11 = A2_new() ;
         A2H12 = A2_new() ;
         DenseMtx_initAsSubmatrix(Hij, H, hijfrow, hijlrow, hkjfcol, hkjlcol) ;
         DenseMtx_setA2(Hij, A2H11) ;
         DenseMtx_initAsSubmatrix(Hkj, H, hijfrow, hijlrow, hijfcol, hijlcol) ;
         DenseMtx_setA2(Hkj, A2H12) ;
         A2_applyQT(A2H12, A2H11, A2H12, &workDV, msglvl, msgFile) ;
         hijfrow += nrhs ;
         hijlrow += nrhs ;
         hkjfcol += nrhs ;
         hkjlcol += nrhs ;
         A2_free(A2H11) ;
         A2_free(A2H12) ;
      }
/*
      -------------------------------------------------
      Compute the QR factorization of the diagonal and
      subdiagonal blocks of the block Hessenberg matrix
      -------------------------------------------------
*/
      DenseMtx_initAsSubmatrix(Hij, H, hijfrow, hijlrow, hijfcol, hijlcol) ;
      DenseMtx_setA2(Hij, A2Hij) ;
      ops = A2_QRreduce(A2Hij, &workDV, msglvl, msgFile) ;
/*
      ----------------------------------------------------------
      apply the transformation to the corresponding section of G
      ----------------------------------------------------------
*/
      DenseMtx_initAsSubmatrix(Gj, G, hijfrow, hijlrow, 0, nrhs-1) ;
      DenseMtx_setA2(Gj, A2Gj) ;
      A2_applyQT(A2Gj, A2Hij, A2Gj, &workDV, msglvl, msgFile) ;
      DenseMtx_initAsSubmatrix(Gj, G, hijfrow+nrhs, hijlrow, 0, nrhs-1) ;
      resnorm = DenseMtx_frobNorm(Gj) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n outer %d, inner %d, res = %24.16e",
                 jouter, jinner, resnorm) ;
         fflush(msgFile) ;
      }
/*
      -----------------
      check convergence
      -----------------
*/
      if ( resnorm <= convergetol * initResNorm ) {
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n BGMRESR convergence has been achieved") ;
            fprintf(msgFile, "\n ***convergetol = %24.16e", convergetol) ;
            fprintf(msgFile, "\n ***initResNorm = %24.16e", initResNorm) ;
            fflush(msgFile) ;
         }
/*
         ------------------------------------------------
         convergence has been achieved, break out of loop
         ------------------------------------------------
*/
         nstep = jinner ;
         converged = 1 ;
         A2_free(A2Gj) ;
         A2_free(A2Hij) ;
         A2_free(A2Hkj) ;
         break ;
      }
      A2_free(A2Gj) ;
      A2_free(A2Hij) ;
      A2_free(A2Hkj) ;
   }
/*
   ----------------------------------
   for each rhs copy Gi into vector Y
   solve Y := H^{-1} Y
   and store Y back into Gi
   ----------------------------------
*/
   DenseMtx_setA2(H, A2H) ;
   DenseMtx_setA2(G, A2G) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n H") ;
      A2H->n1 = nstep+1 ;
      A2H->n2 = nstep+1 ;
      A2_writeForHumanEye(A2H, msgFile) ;
      A2G->n1 = nstep+1 ;
      fprintf(msgFile, "\n\n G") ;
      A2_writeForHumanEye(A2G, msgFile) ;
      A2H->n1 = maxninner+1 ;
      A2H->n2 = maxninner+1 ;
      fprintf(msgFile, "\n\n before solve, Y") ;
         fprintf(msgFile, "\n\n before solve ykinner = %24.16e ",Y[7]) ;
      DVfprintf(msgFile, nstep+1, Y) ;
      fflush(msgFile) ;
   }
   for ( i = 0 ; i < nrhs ; i++ ) {
      A2_extractColumn(A2G, Y, i) ;
      for ( kinner = nstep ; kinner >= 0 ; kinner-- ) {
         A2_realEntry(A2H, kinner, kinner, &Hkk) ;
         Y[kinner] = Y[kinner] / Hkk ;
         for ( iinner = 0 ; iinner < kinner ; iinner++ ) {
            A2_realEntry(A2H, iinner, kinner, &Hik) ;
            Y[iinner] -= Y[kinner] * Hik ;
         }
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n after solve, Y") ;
         DVfprintf(msgFile, nstep+1, Y) ;
         fflush(msgFile) ;
      }
      A2_setColumn(A2G, Y, i) ;
   }
/*
   -------------------------------------
   compute the new solution X := X + V*Y
   -------------------------------------
*/
   DenseMtx_mmm("n", "n", one, mtxX, one, V, G) ;
   if ( msglvl > 2 ) {
      resnorm = DenseMtx_frobNorm(mtxX) ;
      fprintf(msgFile, "\n ||X||_f = %24.16e", resnorm) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n X") ;
      DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   }
   if ( converged == 1 ) {
      break ;
   }
   A2_free(A2Gi) ;
   A2_free(A2Vj) ;
}
/*
   ------------------
   final solve MX = Y
   ------------------
*/
if ( mtxM != NULL ) {
   FrontMtx_solve(mtxM, mtxX, mtxX, mtxM->manager, 
                  cpus, 0, msgFile) ;
}
if ( msglvl > 2 ) {
   resnorm = DenseMtx_frobNorm(mtxX) ;
   fprintf(msgFile, "\n ***after solve ||X||_f = %24.16e", resnorm) ;
}
/*
   -----------
   return info
   -----------
*/
*pnouter = jouter + 1 ;
*pninner = jinner + 1 ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n before leaving, *pnouter = %d, *pinner = %d",
           *pnouter, *pninner) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
   A2_free(A2G) ;
   A2_free(A2H) ;
   A2_free(A2Z) ;
   A2_free(A2Gi) ;
   A2_free(A2Vj) ;
   DenseMtx_free(Gj) ;
   DenseMtx_free(G) ;
   DenseMtx_free(Hij) ;
   DenseMtx_free(Hkj) ;
   DenseMtx_free(H) ;
   DenseMtx_free(Vi) ;
   DenseMtx_free(Vj) ;
   DenseMtx_free(Vk) ;
   DenseMtx_free(V) ;
   DenseMtx_free(Z) ;
   DenseMtx_free(W) ;
   DV_clearData(&workDV) ;
   DVfree(Y) ;

return(1) ; }

/*--------------------------------------------------------------------*/
