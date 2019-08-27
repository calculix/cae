/*  DDviaProjection.c  */

#include"../GPart.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   set the compids[] vector using a global map from vertices
   to domains and interface nodes.

   DDmapIV -- IV object that contains the map from vertices
              to domains and interface nodes

   created -- 96mar17, cca
   ---------------------------------------------------------
*/
void
GPart_DDviaProjection (
   GPart   *gpart,
   IV      *DDmapIV
) {
int   *compids, *domainMap, *map, *vtxMap ;
int   dom, domloc, ndom, ndomloc, nvtx, vglob, vloc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || DDmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_DDviaProjection(%p,%p)"
           "\n bad input\n", gpart, DDmapIV) ;
   exit(-1) ;
}
nvtx    = gpart->nvtx ;
compids = IV_entries(&gpart->compidsIV) ;
/*
   --------------------------
   find the number of domains
   --------------------------
*/
vtxMap = IV_entries(&gpart->vtxMapIV) ;
map    = IV_entries(DDmapIV) ;
ndom   = IV_max(DDmapIV) ;
/*
   ------------------------
   check for a quick return
   ------------------------
*/
if ( gpart->par == NULL ) {
   IVcopy(nvtx, compids, map) ;
   gpart->ncomp = ndom ;
   return ;
}
/*
   ----------------------------------------
   fill compids[] with the local domain ids
   ----------------------------------------
*/
domainMap = IVinit(ndom+1, -1) ;
ndomloc = 0 ;
for ( vloc = 0 ; vloc < nvtx ; vloc++ ) {
   vglob = vtxMap[vloc] ;
   if ( (dom = map[vglob]) > 0 ) {
      if ( (domloc = domainMap[dom]) == -1 ) {
         domloc = domainMap[dom] = ++ndomloc ;
      }
      compids[vloc] = domloc ;
   } else {
      compids[vloc] = 0 ;
   }
}
gpart->ncomp = ndomloc ;
IVfree(domainMap) ; 

return ; }

/*--------------------------------------------------------------------*/
