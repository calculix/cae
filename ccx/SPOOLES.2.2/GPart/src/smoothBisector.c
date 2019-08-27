/*  smoothBisector.c  */

#include "../GPart.h"
#include "../../Ideq.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   static definition of evaluation function
   ----------------------------------------
*/
static float eval ( float alpha, float wS, float wB, float wW) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   smooth a wide separator
 
   nlevel -- number of levels one each side of the separator
             to include into the separator
   alpha  -- partition cost function parameter
   
   created -- 96jun02, cca
   ---------------------------------------------------------
*/
float
GPart_smoothBisector (
   GPart   *gpart,
   int     nlevel,
   float   alpha
) {
FILE      *msgFile ;
float     balance, bestcost, cost, newcost ;
Graph     *g ;
int       ierr, ipass, msglvl ;
int       *compids, *cweights, *vwghts ;
IV        *YCmapIV, *YVmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || nlevel < 0 || alpha < 0.0 ) {
   fprintf(stderr, "\n fatal error in GPart_smoothBisector(%p,%d,%f)"
           "\n bad input\n", gpart, nlevel, alpha) ;
   exit(-1) ;
}
g        = gpart->g        ;
compids  = IV_entries(&gpart->compidsIV)  ;
cweights = IV_entries(&gpart->cweightsIV) ;
vwghts   = g->vwghts ;
msglvl   = gpart->msglvl  ;
msgFile  = gpart->msgFile ;
bestcost = eval(alpha, cweights[0], cweights[1],  cweights[2]) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n smoothBisector : nlevel = %d, alpha = %f", 
           nlevel, alpha) ;
   fprintf(msgFile, "\n old partition cost %.3f, weights = %5d %5d %5d",
           bestcost, cweights[0], cweights[1], cweights[2]) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n compids") ;
   IVfp80(msgFile, g->nvtx, compids, 80, &ierr) ;
}
/*
   ------------------------------------
   loop while improvement is being made
   ------------------------------------
*/
ipass = 0 ;
while ( 1 ) {
   if ( msglvl > 1 ) {
      if ( cweights[1] >= cweights[2] ) {
         balance = ((double) cweights[1]) / cweights[2] ;
      } else {
         balance = ((double) cweights[2]) / cweights[1] ;
      }
      cost = cweights[0] * (1 + alpha * balance) ;
      fprintf(msgFile, 
"\n\n ### pass %d, cweights : %d %d %d, balance = %5.3f, cost = %.1f",
              ipass, cweights[0], cweights[1], cweights[2],
              balance, cost) ;
      fflush(msgFile) ;
   }
/*
   ---------------------------
   identify the wide separator
   ---------------------------
*/
   YVmapIV = GPart_identifyWideSep(gpart, nlevel, nlevel) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n nlevel = %d, |widesep| = %d", 
              nlevel, IV_size(YVmapIV)) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n YVmapIV") ;
      IV_writeForHumanEye(YVmapIV, msgFile) ;
   }
/*
   -------------------------------------------------
   get the Y = Y_0 cup Y_1 cup Y_2 cup Y_3 partition
   -------------------------------------------------
*/
   YCmapIV = GPart_makeYCmap(gpart, YVmapIV) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n YCmapIV found") ;
      fflush(msgFile) ;
   }
/*
   --------------------
   smooth the separator
   --------------------
*/
   newcost = GPart_smoothYSep(gpart, YVmapIV, YCmapIV, alpha) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n newcost = %.3f", newcost) ;
      fflush(msgFile) ;
   }
/*
   ------------------------
   free the two map vectors
   ------------------------
*/
   IV_free(YVmapIV) ;
   IV_free(YCmapIV) ;
/*
   -------------------------------
   check for an improved partition
   -------------------------------
*/
   if ( newcost == bestcost ) {
      break ;
   } else {
      bestcost = newcost ;
   }
   ipass++ ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n\n final partition weights [%d %d %d], cost = %.3f",
           cweights[0], cweights[1], cweights[2], bestcost) ;
      fflush(msgFile) ;
}

return(bestcost) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   partition evaluation function
   -----------------------------
*/
static float
eval (
   float   alpha,
   float   wS,
   float   wB,
   float   wW
) {
float   cost ;
#if MYDEBUG > 2
fprintf(msgFile, "\n alpha = %f, wS = %f, wB = %f, wW = %f",
        alpha, wS, wB, wW) ; 
#endif
if ( wB == 0 || wW == 0 ) {
   cost = (wS + wB + wW) * (wS + wB + wW) ;
} else if ( wB >= wW ) {
   cost = wS*(1. + (alpha*wB)/wW) ;
} else {
   cost = wS*(1. + (alpha*wW)/wB) ;
}
#if MYDEBUG > 2
fprintf(msgFile, ", cost = %f", cost) ;
#endif
return(cost) ; }

/*--------------------------------------------------------------------*/
