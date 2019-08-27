/*  makeYCmap.c  */

#include "../GPart.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   make the map from wide separator vertices Y 
   to components {0, 1, 2, 3}.

   YCmap[y] == 0 --> y is not adjacent to either component
   YCmap[y] == 1 --> y is adjacent to only component 1
   YCmap[y] == 2 --> y is adjacent to only component 2
   YCmap[y] == 3 --> y is adjacent to components 1 and 2

   created -- 96jun09, cca
   -------------------------------------------------------
*/
IV *
GPart_makeYCmap (
   GPart   *gpart,
   IV      *YVmapIV
) {
Graph   *g ;
int     ii, nvtx, nY, v, vsize, w, y ;
int     *compids, *vadj, *VYmap, *YCmap, *YVmap ;
IV      *YCmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (g = gpart->g) == NULL 
   || (nvtx = gpart->nvtx) <= 0
   || YVmapIV == NULL || (nY = IV_size(YVmapIV)) <= 0 
   || (YVmap = IV_entries(YVmapIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_makeYCmap(%p,%p)"
           "\n bad input\n", gpart, YVmapIV) ;
   if ( YVmapIV != NULL ) {
      fprintf(stderr, "\n YVmapIV") ;
      IV_writeForHumanEye(YVmapIV, stderr) ;
   }
   exit(-1) ;
}
compids = IV_entries(&gpart->compidsIV) ;
/*
   --------------------------------
   generate the inverse V --> Y map 
   --------------------------------
*/
VYmap = IVinit(nvtx, -1) ;
for ( y = 0 ; y < nY ; y++ ) {
   v = YVmap[y] ;
   VYmap[v] = y ;
}
/*
   ------------------------------------
   initialize the Y --> C map IV object
   ------------------------------------
*/
YCmapIV = IV_new();
IV_init(YCmapIV, nY, NULL) ;
YCmap = IV_entries(YCmapIV) ;
/*
   ---------------
   fill the fields
   ---------------
*/
for ( y = 0 ; y < nY ; y++ ) {
   YCmap[y] = 0 ;
   v = YVmap[y] ;
   Graph_adjAndSize(g, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      w = vadj[ii] ;
      if ( w < nvtx && VYmap[w] == -1 ) {
/*
         --------------------------------
         w is not in the wide separator Y
         --------------------------------
*/
         if ( compids[w] == 1 ) {
/*
            ---------------------------------------
            v is adjacent to component 1 setminus Y
            ---------------------------------------
*/
            if ( YCmap[y] == 2 ) {
/*
               ------------------------------------
               v is already adjacent to component 2
               so it is adjacent to both components
               ------------------------------------
*/
               YCmap[y] = 3 ;
               break ;
            } else {
/*
               ----------------------------------
               set map value but keep on checking
               ----------------------------------
*/
               YCmap[y] = 1 ;
            }
         } else if ( compids[w] == 2 ) {
/*
            ---------------------------------------
            v is adjacent to component 2 setminus Y
            ---------------------------------------
*/
            if ( YCmap[y] == 1 ) {
/*
               ------------------------------------
               v is already adjacent to component 1
               so it is adjacent to both components
               ------------------------------------
*/
               YCmap[y] = 3 ;
               break ;
            } else {
/*
               ----------------------------------
               set map value but keep on checking
               ----------------------------------
*/
               YCmap[y] = 2 ;
            }
         }
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(VYmap) ;

return(YCmapIV) ; }

/*--------------------------------------------------------------------*/
