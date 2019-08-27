/*  PatchAndGoInfo.h  */

#include "../IV.h"
#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   this object is used by the Chv object to implement the 
   "patch-and-go" strategies, used in specialized optimization
   and structural analysis applications.

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
   toosmall -- measure of smallness for diagonal entry
   fudge    -- change factor for diagonal entry
   fudgeIV  -- (optional) stores locations where modifications made
   fudgeDV  -- (optional) stores modifications 
   
   created -- 98aug26, cca
   -----------------------------------------------------------------
*/
typedef struct _PatchAndGoInfo PatchAndGoInfo ;
struct _PatchAndGoInfo {
   int      strategy ;
   double   toosmall ;
   double   fudge    ;
   IV       *fudgeIV ;
   DV       *fudgeDV ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor method
 
   created -- 98aug26, cca
   -----------------------
*/
PatchAndGoInfo *
PatchAndGoInfo_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98aug26, cca
   -----------------------
*/
void
PatchAndGoInfo_setDefaultFields ( 
   PatchAndGoInfo   *info
) ;
/*
   -----------------------
   clear the data fields
 
   created -- 98aug26, cca
   -----------------------
*/
void
PatchAndGoInfo_clearData (
   PatchAndGoInfo   *info
) ;
/*
   -----------------------
   destructor
 
   created -- 98aug26, cca
   -----------------------
*/
void
PatchAndGoInfo_free (
   PatchAndGoInfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
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
) ;
/*--------------------------------------------------------------------*/
