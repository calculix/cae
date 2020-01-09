/*  smoothBy2layers.c  */

#include "../GPart.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   static declaration of evaluation function
   -----------------------------------------
*/
static float eval ( float alpha, float wS, float wB, float wW ) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   smooth a bisector using the alternating two-layer algorithm

   option -- network flag
      1 --> make the network bipartite as for the
            Dulmage-Mendelsohn decomposition
      otherwise -- use network induced by the wide separator
   alpha -- cost function parameter

   created -- 96jun08, cca
   -----------------------------------------------------------
*/
void
GPart_smoothBy2layers (
   GPart   *gpart,
   int     option,
   float   alpha
) {
FILE    *msgFile ;
float   bestcost, newcost ;
int     ierr, large, msglvl, nY, pass, small, y ;
int     *cweights, *YCmap ;
IV      *YCmapIV, *YVmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || alpha < 0.0 ) {
   fprintf(stderr, "\n fatal error in GPart_smoothBy2layers(%p,%f)"
           "\n bad input\n", gpart, alpha) ;
   exit(-1) ;
}
pass     = 1 ;
cweights = IV_entries(&gpart->cweightsIV) ;
bestcost = eval(alpha, cweights[0], cweights[1], cweights[2]) ;
msgFile  = gpart->msgFile ;
msglvl   = gpart->msglvl  ;
while ( 1 ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
        "\n\n PASS %d : attempting a move towards the larger component",
        pass) ;
      fflush(msgFile) ;
   }
   if ( cweights[1] >= cweights[2] ) {
      large = 1 ; small = 2 ;
      YVmapIV = GPart_identifyWideSep(gpart, 1, 0) ;
   } else {
      large = 2 ; small = 1 ;
      YVmapIV = GPart_identifyWideSep(gpart, 0, 1) ;
   }
   YCmapIV = GPart_makeYCmap(gpart, YVmapIV) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n YCmapIV") ;
      IV_writeForHumanEye(YCmapIV, msgFile) ;
      fflush(msgFile) ;
   }
   IV_sizeAndEntries(YCmapIV, &nY, &YCmap) ;
   if ( option == 1 ) {
/*
      ----------------------------
      generate a bipartite network
      ----------------------------
*/
      for ( y = 0 ; y < nY ; y++ ) {
         if ( YCmap[y] != small ) {
            YCmap[y] = large ;
         }
      }
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n calling GPartSmoothYSep") ;
      fflush(msgFile) ;
   }
   newcost = GPart_smoothYSep(gpart, YVmapIV, YCmapIV, alpha) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
           "\n\n bestcost = %.2f, newcost = %.2f", bestcost, newcost) ;
      fflush(msgFile) ;
   }
   IV_free(YVmapIV) ;
   IV_free(YCmapIV) ;
   if ( newcost == bestcost ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
       "\n\n PASS %d : attempting a move towards the smaller component",
           pass) ;
         fflush(msgFile) ;
      }
      if ( cweights[1] >= cweights[2] ) {
         large = 1 ; small = 2 ;
         YVmapIV = GPart_identifyWideSep(gpart, 0, 1) ;
      } else {
         large = 2 ; small = 1 ;
         YVmapIV = GPart_identifyWideSep(gpart, 1, 0) ;
      }
      YCmapIV = GPart_makeYCmap(gpart, YVmapIV) ;
      IV_sizeAndEntries(YCmapIV, &nY, &YCmap) ;
      if ( option == 1 ) {
   /*
         ----------------------------
         generate a bipartite network
         ----------------------------
   */
         for ( y = 0 ; y < nY ; y++ ) {
            if ( YCmap[y] != large ) {
               YCmap[y] = small ;
            }
         }
      }
      newcost = GPart_smoothYSep(gpart, YVmapIV, YCmapIV, alpha) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
           "\n\n bestcost = %.2f, newcost = %.2f", bestcost, newcost) ;
         fflush(msgFile) ;
      }
      IV_free(YVmapIV) ;
      IV_free(YCmapIV) ;
   }
   if ( newcost == bestcost ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n no improvement made") ;
         fflush(msgFile) ;
      }
      break ;
   } else {
      bestcost = newcost ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n improvement made") ;
         fflush(msgFile) ;
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n compids") ;
         IVfp80(msgFile, gpart->nvtx, IV_entries(&gpart->compidsIV),
                80, &ierr) ;
         fflush(msgFile) ;
      }
   }
   pass++ ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n leaving smoothBy2layers") ;
   fflush(msgFile) ;
}
return ; }

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
float   bestcost ;
 
if ( wB == 0 || wW == 0 ) {
   bestcost = (wS + wB + wW) * (wS + wB + wW) ;
} else if ( wB >= wW ) {
   bestcost = wS*(1. + (alpha*wB)/wW) ;
} else {
   bestcost = wS*(1. + (alpha*wW)/wB) ;
}
return(bestcost) ; }
 
/*--------------------------------------------------------------------*/
