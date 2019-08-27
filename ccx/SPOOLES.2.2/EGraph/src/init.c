/*  init.c  */

#include "../EGraph.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- initialize the EGraph object

   created -- 96oct24, cca
   ---------------------------------------
*/
void
EGraph_init (
   EGraph   *egraph,
   int      type,
   int      nelem,
   int      nvtx,
   int      IVL_type
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  egraph == NULL || type < 0 || type > 1 
   || nelem <= 0 || nvtx <= 0 ) {
   fprintf(stderr, "\n fatal error in EGraph_init(%p,%d,%d,%d,%d)"
           "\n bad input\n", egraph, type, nelem, nvtx, IVL_type) ;
   exit(-1) ;
}
/*
   ----------------------------
   clear the data in the object
   ----------------------------
*/
EGraph_clearData(egraph) ;
/*
   ---------------------
   initialize the object
   ---------------------
*/
egraph->type  = type ;
egraph->nelem = nelem ;
egraph->nvtx  = nvtx  ;
egraph->adjIVL = IVL_new() ;
IVL_init1(egraph->adjIVL, IVL_type, nelem) ;
if ( type == 1 ) {
   egraph->vwghts = IVinit(nvtx, 0) ;
}

return ; }

/*--------------------------------------------------------------------*/
