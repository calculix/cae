/*  init.c  */

#include "../PatchAndGoInfo.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to initialize the patch-and-go information object

   strategy -- type of patch-and-go strategy
      1 -- used with optimization matrices
         if ( |a_{i,i}| <= toosmall ) {
            set a_{i,i} = 1.0
            set offdiagonals = 0.0 
         }
      2 -- used with structural analysis matrices
         if ( |a_{i,i}| <= fudge ) {
            set a_{i,i} = fudge * max (1, |a_{i,*}|, |a_{*,i}|)
         }
   toosmall    -- tolerance for first strategy
   fudge       -- fudge factor for second strategy
   storeids    -- if nonzero, the row and column numbers where
      patches have been applied are stored in an IV object
   storevalues -- if nonzero and strategy 2, the differences
      between the old and new diagonal magnitudes will be 
      stored in a DV object

   created -- 98aug27, cca
   ------------------------------------------------------------
*/
void
PatchAndGoInfo_init (
   PatchAndGoInfo   *info,
   int              strategy,
   double           toosmall,
   double           fudge,
   int              storeids,
   int              storevalues
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( info == NULL || strategy < 1 || strategy > 2
     || toosmall < 0.0 || fudge < 0.0 ) {
   fprintf(stderr, "\n fatal error in PatchAndGoInfo_init()"
           "\n bad input\n") ;
   exit(-1) ;
}
PatchAndGoInfo_clearData(info) ;

info->strategy = strategy ;
info->toosmall = toosmall ;
info->fudge    = fudge    ;
if ( storeids != 0 ) {
   info->fudgeIV = IV_new() ;
}
if ( strategy == 2 && storevalues != 0 ) {
   info->fudgeDV = DV_new() ;
}
return ; }

/*--------------------------------------------------------------------*/
