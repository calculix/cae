/*  fidmat.c  */

#include "../BKL.h"

#define MYDEBUG 0
#define MAXNDOM_FOR_EXHAUSTIVE_SEARCH 8

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   structure used in this file
   ---------------------------
*/
typedef struct _cell   Cell ;
struct _cell {
   int    domid  ;
   int    deltaS ;
   int    deltaB ;
   int    deltaW ;
   Cell   *prev  ;
   Cell   *next  ;
} ;

static Cell   Head, *head = &Head ;
static Cell   Undo, *undo = &Undo ;

static float
BKL_fidmatPass ( 
   BKL     *bkl,
   Cell    cells[],
   int     tags[],
   Graph   *DomByDom,
   int     npass
) ;

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   improve the partition using the FidMat algorithm

   created  -- 95oct11, cca
   modified -- 95dec07, cca
      memory leak fixed, comments added, DomByDom inserted
   -------------------------------------------------------
*/
float
BKL_fidmat ( 
   BKL   *bkl
) {
float   cost ;
int     ndom ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_fidmat(%p)"
           "\n bad input\n", bkl) ;
   exit(-1) ;
}
ndom = bkl->ndom ;
/*
   ---------------------------------------------
   if ndom <= MAXNDOM_FOR_EXHAUSTIVE_SEARCH then
      do exhaustive search
   else
      do fidmat sweeps
   endif
   ---------------------------------------------
*/
if ( ndom <= MAXNDOM_FOR_EXHAUSTIVE_SEARCH ) {
   int   mdom, *domids, *tcolors ;
/*
   --------------------
   do exhaustive search
   --------------------
*/
   mdom    = ndom - 1 ;
   domids  = IVinit(mdom, -1) ;
   tcolors = IVinit(mdom, -1) ;
   IVramp(mdom, domids, 1, 1) ;
   BKL_exhSearch(bkl, mdom, domids, tcolors) ;
   IVfree(domids)  ;
   IVfree(tcolors) ;
   cost = BKL_evalfcn(bkl) ;
} else {
   Cell    *cell, *cells ;
   float   bestcost ;
   Graph   *DomByDom ;
   int     idom ;
   int     *tags ;
/*
   ---------------------------------------------------
   initialize the cell objects and the working vectors
   ---------------------------------------------------
*/
   ALLOCATE(cells, struct _cell, ndom) ;
   tags = IVinit(ndom, -1) ;
   for ( idom = 0, cell = cells ; idom < ndom ; idom++, cell++ ) {
      cell->domid = idom ;
      cell->deltaS = cell->deltaB = cell->deltaW = 0 ;
      cell->prev   = cell->next = cell ;
   }
/*
   -------------------------------------
   create the domain-domain Graph object
   -------------------------------------
*/
   DomByDom = BPG_makeGraphXbyX(bkl->bpg) ;
#if MYDEBUG > 1 
   fprintf(stdout, "\n\n domain-domain Graph object") ;
   Graph_writeForHumanEye(DomByDom, stdout) ;
   fflush(stdout) ;
#endif
/*
   --------------------------
   make the first fidmat pass
   --------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n\n pass %d, cost %.2f, < %d %d %d >", 
           bkl->npass, BKL_evalfcn(bkl),
           bkl->cweights[0], bkl->cweights[1], bkl->cweights[2]) ;
   fflush(stdout) ;
#endif
   bkl->npass = 1 ;
   bestcost = BKL_fidmatPass(bkl, cells, tags, DomByDom, bkl->npass) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n pass %d, cost %.2f, < %d %d %d >", 
           bkl->npass, bestcost,
           bkl->cweights[0], bkl->cweights[1], bkl->cweights[2]) ;
   fflush(stdout) ;
#endif
/*
   ---------------------------------------------------
   make additional passes while the partition improves
   ---------------------------------------------------
*/
   while ( 1 ) {
      bkl->npass++ ;
      cost = BKL_fidmatPass(bkl, cells, tags, DomByDom, bkl->npass) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n\n pass %d, cost %.2f, < %d %d %d >", 
              bkl->npass, cost,
              bkl->cweights[0], bkl->cweights[1], bkl->cweights[2]) ;
      fflush(stdout) ;
#endif
      if ( cost < bestcost ) {
         bestcost = cost ;
      } else {
         break ;
      }
   }
/*
   ------------------------
   free the working storage
   ------------------------
*/
   FREE(cells) ;
   IVfree(tags) ;
   Graph_free(DomByDom) ;
}

return(cost) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   make one pass of the FidMat algorithm

   created -- 95oct11, cca
   -------------------------------------
*/
static float
BKL_fidmatPass ( 
   BKL     *bkl,
   Cell    cells[],
   int     tags[],
   Graph   *DomByDom,
   int     npass
) {
Cell    *cell ;
float   bestcost, bettercost, cost ;
int     dom, dom2, ii, ndom, size ;
int     *cweights, *doms ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || cells == NULL || tags == NULL || DomByDom == NULL ){
   fprintf(stderr, "\n fatal error in BKL_fidmatPass(%p,%p,%p,%p,%d)"
           "\n bad input\n", bkl, cells, tags, DomByDom, npass) ;
   exit(-1) ;
}
ndom     = bkl->ndom     ;
cweights = bkl->cweights ;
/*
   ------------------------------
   evaluate the current partition
   ------------------------------
*/
bestcost = BKL_evalfcn(bkl) ;
/*
   -----------------------------------------------------
   fill the cells with domains adjacent to the separator
   -----------------------------------------------------
*/
head->next = head->prev = head ;
undo->next = undo->prev = undo ;
for ( dom = 0 ; dom < ndom ; dom++ ) {
   cell = &cells[dom] ;
   cell->domid = dom ;
   cell->prev  = cell->next = cell ;
   if ( BKL_domAdjToSep(bkl, dom) == 1 ) {
/*
      ----------------------------------------
      domain dom is adjacent to the separator
      evaluate its change and insert into list
      ----------------------------------------
*/
      BKL_evalgain(bkl, dom, &cell->deltaS,
                   &cell->deltaB, &cell->deltaW) ;
#if MYDEBUG > 1
      fprintf(stdout, "\n loading domain %d, <%d %d %d>",
              dom, cell->deltaS, cell->deltaB, cell->deltaW) ;
      fflush(stdout) ;
#endif
      DLIST_TAIL_INSERT(head, cell) ;
   }
}
/*
   ---------------
   loop over moves
   ---------------
*/
while ( head->next != head ) {
/*
   -----------------------------------------------
   find best move to make w.r.t. the cost function
   -----------------------------------------------
*/
   cell = head->next ;
   dom  = cell->domid ;
   bettercost = BKL_eval(bkl, cweights[0] + cell->deltaS,
                              cweights[1] + cell->deltaB,
                              cweights[2] + cell->deltaW) ;
#if MYDEBUG > 1
   fprintf(stdout, "\n domain %d, move cost = %.1f",
           dom, bettercost) ;
   fflush(stdout) ;
#endif
   for ( cell = cell->next ; cell != head ; cell = cell->next ) {
      cost = BKL_eval(bkl, cweights[0] + cell->deltaS,
                           cweights[1] + cell->deltaB,
                           cweights[2] + cell->deltaW) ;
#if MYDEBUG > 1
      fprintf(stdout, "\n domain %d, move cost = %.1f",
              cell->domid, cost) ;
      fflush(stdout) ;
#endif
      if ( cost < bettercost ) {
#if MYDEBUG > 1
         fprintf(stdout, ", better") ;
         fflush(stdout) ;
#endif
         dom = cell->domid ;
         bettercost = cost ;
      }
   }
/*
   -----------------------------
   remove the node from the list
   -----------------------------
*/
   cell = &cells[dom] ;
   DLIST_DELETE(cell) ;
/*
   ---------------
   flip the domain
   ---------------
*/
#if MYDEBUG > 1
   fprintf(stdout, "\n flipping domain %d, cweights <%d %d %d>", 
           dom, cweights[0] + cell->deltaS,
           cweights[1] + cell->deltaB,
           cweights[2] + cell->deltaW) ;
   fflush(stdout) ;
#endif
   BKL_flipDomain(bkl, dom) ;
   cost = BKL_eval(bkl, cweights[0], cweights[1], cweights[2]) ;
#if MYDEBUG > 1
   fprintf(stdout, ", cost = %.1f", cost) ;
   fflush(stdout) ;
#endif
   if ( bestcost > cost ) {
/*
      ---------------------------------------------
      better partition found, set undo list to NULL
      ---------------------------------------------
*/
      bestcost = cost ;
      DLIST_DELETE(undo) ;
      bkl->nimprove++ ;
   } else {
/*
      ----------------------------------------------
      partition is not better, add move to undo list
      ----------------------------------------------
*/
      DLIST_HEAD_INSERT(undo, cell) ;
   }
/*
   --------------------------
   loop over adjacent domains
   --------------------------
*/
   tags[dom] = npass ;
   Graph_adjAndSize(DomByDom, dom, &size, &doms) ;
   for ( ii = 0 ; ii < size ; ii++ ) {
      dom2 = doms[ii] ;
      if ( tags[dom2] < npass && BKL_domAdjToSep(bkl, dom2) == 1 ) {
/*
         ------------------------------------
         domain dom2 has not yet been flipped
         and is adjacent to the separator
         ------------------------------------
*/
#if MYDEBUG > 1
        fprintf(stdout, 
                "\n domain %d, not yet flipped, adj to separator",
                dom2) ;
        fflush(stdout) ;
#endif
         cell = &cells[dom2] ;
         BKL_evalgain(bkl, dom2, &cell->deltaS,
                      &cell->deltaB, &cell->deltaW) ;
#if MYDEBUG > 1
        fprintf(stdout, ", gain < %d %d %d >",
                cell->deltaS, cell->deltaB, cell->deltaW) ;
        fflush(stdout) ;
#endif
         if ( cell->prev == cell ) {
/*
            ---------------------------------------------
            domain is not on the list of domains eligible
            to move, insert on the tail of the list
            ---------------------------------------------
*/
            DLIST_TAIL_INSERT(head, cell) ;
         }
      }
   }
}
/*
   -------------------------------------------------------
   undo the flips of domains since the last best partition
   -------------------------------------------------------
*/
while ( (cell = undo->next) != undo ) {
#if MYDEBUG > 1
   fprintf(stdout, "\n un-flipping domain %d", cell->domid) ;
   fflush(stdout) ;
#endif
   DLIST_DELETE(cell) ;
   BKL_flipDomain(bkl, cell->domid) ;
}
return(bestcost) ; }

/*--------------------------------------------------------------------*/
