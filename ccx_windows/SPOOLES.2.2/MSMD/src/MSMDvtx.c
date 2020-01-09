/*  MSMDvtx.c  */

#include "../MSMD.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   print a vertex  object
 
   created -- 96feb25, cca
   -----------------------
*/
void
MSMDvtx_print (
   MSMDvtx   *v,
   FILE      *fp
) {
int   ierr ;
/*
   ---------------
   check the input
   ---------------
*/
if ( v == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in MSMDvtx_print(%p,%p)"
           "\n bad input\n", v, fp) ;
   exit(-1) ;
}
fprintf(fp, 
        "\n vertex %d, weight %d, mark %c, status %c, stage = %d",
        v->id, v->wght, v->mark, v->status, v->stage) ;
switch ( v->status ) {
case 'O' :
case 'D' :
case 'R' :
case 'B' :
   fprintf(fp, "\n    edges(%d) :", v->nadj) ;
   IVfp80(fp, v->nadj, v->adj, 13, &ierr) ;
   fprintf(fp, "\n    subtrees :") ;
   IP_fp80(fp, v->subtrees, 13) ;
   break ;
case 'L' :
case 'E' :
   fprintf(fp, 
           "\n    parent = %d", (v->par == NULL) ? -1 : v->par->id) ;
   fprintf(fp, "\n    bnd(%d), weight = %d :", v->nadj, v->bndwght) ;
   IVfp80(fp, v->nadj, v->adj, 10, &ierr) ;
   break ;
case 'I' :
   fprintf(fp, 
           "\n    parent = %d", (v->par == NULL) ? -1 : v->par->id) ;
   break ;
default :
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
