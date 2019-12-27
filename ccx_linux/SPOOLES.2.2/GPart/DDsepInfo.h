/*  DDsepInfo.h  */

#ifndef _DDsepInfo_
#define _DDsepInfo_

#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   this structure is used in the GPart_RBviaDDSEP method

   seed          -- random number seed
   minweight     -- minimum weight for any fishnet domain,
                    default value is 50
   maxweight     -- maximum weight for any fishnet domain,
                    default value is 100
   freeze        -- multiplier used to freeze nodes in the multisector,
                    default value is 4.0
   alpha         -- cost function parameter, default value is 1.0
   maxcompweight -- any leaf component must have 
                    weight less than this value
   ntreeobj      -- accumulator for the number of objects
                    in the domain/separator tree
   DDoption      -- option for finding domain decomposition
      1 --> use Fishnet on each subgraph
      2 --> use Fishnet on graph, then projection on each subgraph
   nlayer        -- number of layers to be used in the smoothing process
      2    --> use two-sided dulmage-mendelsohn smoothing
      2k+1 --> use k layers of wide separator smoothing
   cpuDD     -- cpu time for domain decomposition
   cpuMap    -- cpu time to compute the domain/segment maps
   cpuBPG    -- cpu time to create the domain/segment bipartite graphs
   cpuBKL    -- cpu time to find the initial separators using the
                block kernihan-lin algorithm
   cpuSmooth -- cpu time to smooth the separators
   cpuSplit  -- cpu time to split the subgraphs
   cpuTotal  -- total cpu time 
   msglvl    -- message level, default value is 0
   msgFile   -- message file, default is stdout
   ---------------------------------------------------------------------
*/
typedef 
struct _DDsepInfo {
   int      seed          ;
   int      minweight     ;
   int      maxweight     ;
   double   freeze        ;
   double   alpha         ;
   int      maxcompweight ;
   int      ntreeobj      ;
   int      DDoption      ;
   int      nlayer        ;
   double   cpuDD         ;
   double   cpuMap        ;
   double   cpuBPG        ;
   double   cpuBKL        ;
   double   cpuSmooth     ;
   double   cpuSplit      ;
   double   cpuTotal      ;
   int      msglvl        ;
   FILE     *msgFile      ;
} DDsepInfo ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in DDsepInfo.c   -----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   construct a new instance of the DDsepInfo object

   created -- 96feb24, cca
   ------------------------------------------------
*/
DDsepInfo *
DDsepInfo_new (
   void
) ;
/*
   ---------------------------------------------
   set the default fields of the DDsepInfo object

   created  -- 96feb24, cca
   ---------------------------------------------
*/
void
DDsepInfo_setDefaultFields (
   DDsepInfo   *info
) ;
/*
   ---------------------------------------------
   clear the data fields for a DDsepInfo object

   created  -- 96feb24, cca
   ---------------------------------------------
*/
void
DDsepInfo_clearData (
   DDsepInfo   *info
) ;
/*
   ------------------------
   free the DDsepInfo object

   created  -- 96feb24, cca
   ------------------------
*/
void
DDsepInfo_free (
   DDsepInfo   *info
) ;
/*
   -----------------------
   write the CPU times
 
   created -- 97nov06, cca
   -----------------------
*/
void
DDsepInfo_writeCpuTimes (
   DDsepInfo   *info,
   FILE        *msgFile
) ;
/*--------------------------------------------------------------------*/

#endif
