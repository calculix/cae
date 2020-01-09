/*  BKL.h  */

#include "../BPG.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   the BKL object handles the block kernihan-lin family 
   of algorithms defined on bipartite graphs.

   bpg       -- bipartite graph to work with, 
                not free'd when BKL object is free'd
   ndom      -- number of domains
   nseg      -- number of segments
   nreg      -- number of regions, nreg = ndom + nseg
   totweight -- total weight of the bipartite graph
   npass     -- # of passes made
   npatch    -- # of patches evaluated
   nflips    -- # of domain flips that were executed
   nimprove  -- # of improvement steps
   ngaineval -- # of gain evaluations
   colors    -- map from regions to colors, size nreg
   cweights  -- array to store color weights,
      cweights[0] -- separator weight
      cweights[1] -- black component weight
      cweights[2] -- white component weight
   regwghts  -- vector of region weights. 
      if the Graph g has vertex weights then
         regwght points to its vertex weights
         and is not free'd when the BKL object is free'd
      else
         regwghts is allocated and set to unit weights
         and free'd when the BKL object is free'd
      endif
   alpha -- cost function parameter
      cost = |S|(1 + alpha*max(|B|,|W|)/min(|B|,|W|))

   created  -- 95oct07, cca
   modified -- 95dec07, cca
      directory cleaned up, added some efficiency,
      inserted exhaustive search into fidmat procedure.
   -------------------------------------------------------------------
*/
typedef struct _BKL   BKL ;
struct _BKL {
   BPG     *bpg        ;
   int     ndom        ;
   int     nseg        ;
   int     nreg        ;
   int     totweight   ;
   int     npass       ;
   int     npatch      ;
   int     nflips      ;
   int     nimprove    ;
   int     ngaineval   ;
   int     *colors     ;
   int     cweights[3] ;
   int     *regwghts   ;
   float   alpha       ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor

   created -- 95oct07, cca
   -----------------------
*/
BKL *
BKL_new (
   void
) ;
/*
   -----------------------
   set the default fields

   created -- 95oct07, cca
   -----------------------
*/
void
BKL_setDefaultFields (
   BKL   *bkl
) ;
/*
   -----------------------
   clear the data fields

   created -- 95oct07, cca
   -----------------------
*/
void
BKL_clearData (
   BKL   *bkl
) ;
/*
   -----------------------
   destructor

   created -- 95oct07, cca
   -----------------------
*/
void
BKL_free (
   BKL   *bkl
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   initialize the object
 
   created -- 95oct07, cca
   -----------------------
*/
void
BKL_init (
   BKL     *bkl,
   BPG     *bpg,
   float   alpha
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   set the colors of the domains and segments 
   by coloring the domains black or white randomly.
 
   created -- 95oct07, cca
   ------------------------------------------------
*/
void
BKL_setRandomColors (
   BKL   *bkl,
   int   seed
) ;
/*
   -----------------------------------------
   set the component weights.
   note, we assume the domain colors are set
 
   created -- 95oct07, cca
   -----------------------------------------
*/
void
BKL_setColorWeights (
   BKL   *bkl
) ;
/*
   ---------------------------------------
   return the segment color, a function of
   the colors of its neighboring domains.
  
   created -- 95oct07, cca
   ---------------------------------------
*/
int
BKL_segColor (
   BKL   *bkl,
   int   iseg
) ;
/*
   -----------------------
   flip the domain
 
   created -- 95oct07, cca
   -----------------------
*/
void
BKL_flipDomain (
   BKL   *bkl,
   int   idom
) ;
/*
   ------------------------------
   return the next domain to flip
   in a grey code sequence
 
   created -- 95oct07, cca
   ------------------------------
*/
int
BKL_greyCodeDomain (
   BKL   *bkl,
   int   count
) ;
/*
   ----------------------------------------------------------------
  set the initial partition.
 
   flag -- specifies initial partition type
      flag == 1 --> random coloring of domains
      flag == 2 --> one black domain, (seed % ndom), rest are white
      flag == 3 --> one black pseudoperipheral domain, found using
                    domain (seed % ndom) as root, rest are white
      flag == 4 --> roughly half-half split, breadth first search
                    of domains, (seed % ndom) as root
      flag == 5 --> roughly half-half split, breadth first search
                    of domains, (seed % ndom) as root to find
                    a pseudoperipheral domain as root
      flag == 6 --> use domcolors[] to seed the colors[] array
   seed      -- random number seed, for flag == 1, if seed > 0 then
               we call srand48(seed) to set the random number
                generator.
   domcolors -- vector of domain colors, used when flag == 6
 
   created -- 95oct11, cca
   ----------------------------------------------------------------
*/
float
BKL_setInitPart (
   BKL   *bkl,
   int   flag,
   int   seed,
   int   domcolors[]
) ;
/*
   ---------------------------------------------------
   return 1 if the domain is adjacent to the separator
 
   created -- 95oct11, cca
   ---------------------------------------------------
*/
int
BKL_domAdjToSep ( 
   BKL   *bkl, 
   int   dom
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in exhSearch.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   perform an exhaustive search over a subspace of domains
   
   mdom    -- number of domains in the subspace
   domids  -- vector of domain ids, size mdom
   tcolors -- temporary vector to hold active domain colors, size mdom
 
   note : region colors and component weights of the best
          partition are left in bkl->colors[] and bkl->cweights[].
 
   return value -- cost of best partition
 
   created -- 95oct07, cca
   --------------------------------------------------------------------
*/
float
BKL_exhSearch (
   BKL   *bkl,
   int   mdom,
   int   domids[],
   int   tcolors[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in evalfcn.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   evaluate the partition
 
   created -- 95oct07, cca
   -----------------------
*/
float
BKL_evalfcn (
   BKL   *bkl
) ;
/*
   -----------------------
   evaluate the partition
 
   created -- 95oct07, cca
   -----------------------
*/
float
BKL_eval (
   BKL   *bkl,
   int   Sweight,
   int   Bweight,
   int   Wweight
) ;
/*
   ---------------------------------------------------------
   evaluate the (deltaS, deltaB and deltaW) of a domain flip
 
   created -- 950ct11, cca
   ---------------------------------------------------------
*/
void
BKL_evalgain (
   BKL   *bkl,
   int   dom,
   int   *pdeltaS,
   int   *pdeltaB,
   int   *pdeltaW
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in fidmat.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   improve the partition using the FidMat algorithm
 
   created -- 95oct11, cca
   ------------------------------------------------
*/
float
BKL_fidmat ( 
   BKL   *bkl
) ;
/*--------------------------------------------------------------------*/
