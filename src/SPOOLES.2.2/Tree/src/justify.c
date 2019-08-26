/*  justify.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   left-justify a tree by subtree size
   children are linked in ascending order of their subtree size

   created -- 96jun23, cca
   ------------------------------------------------------------
*/
void
Tree_leftJustify (
   Tree   *tree
) {
IV   *tmetricIV, *vmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_leftJustify(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
/*
   ------------------------------------------------------------------
   set the subtree size metric, left justify and free working storage
   ------------------------------------------------------------------
*/
vmetricIV = IV_new() ;
IV_init(vmetricIV, tree->n, NULL) ;
IV_fill(vmetricIV, 1) ;
tmetricIV = Tree_setSubtreeImetric(tree, vmetricIV) ;
Tree_leftJustifyI(tree, tmetricIV) ;
IV_free(vmetricIV) ;
IV_free(tmetricIV) ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   left-justify a tree by a metric
   children are linked in ascending order of their metric

   created -- 96jun23, cca
   ------------------------------------------------------
*/
void
Tree_leftJustifyI (
   Tree   *tree,
   IV     *metricIV
) {
int   i, j, k, n, nexti, prev ;
int   *fch, *metric, *par, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL 
   || (n = tree->n) <= 0 
   || metricIV == NULL 
   || n != IV_size(metricIV)
   || (metric = IV_entries(metricIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_leftJustifyI(%p,%p)"
           "\n bad input\n", tree, metricIV) ;
   exit(-1) ;
}
par = tree->par ;
fch = tree->fch ;
sib = tree->sib ;
/*
   ----------------------------------------------------
   sort all children in decreasing order of metric size
   ----------------------------------------------------
*/
for ( k = 0 ; k < n ; k++ ) {
   for ( i = fch[k], fch[k] = -1 ; i != -1 ; i = nexti ) {
      nexti = sib[i] ;
      for ( j = fch[k], prev = -1 ; j != -1 ; j = sib[j] ) {
         if ( metric[j] < metric[i] ) {
            break ;
         }
         prev = j ;
      }
      if ( prev == -1 ) {
         fch[k] = i ;
      } else {
         sib[prev] = i ;
      }
      sib[i] = j ;
   }
}
/*
   ---------------------------------------------
   sort roots in decreasing order of metric size
   ---------------------------------------------
*/
for ( i = tree->root, tree->root = -1, prev = -1 ; 
      i != -1 ; 
      i = nexti ) {
   nexti = sib[i] ;
   for ( j = tree->root, prev = -1 ; j != -1 ; j = sib[j] ) {
       if ( metric[j] < metric[i] ) {
          break ;
       }
      prev = j ;
   }
   if ( prev == -1 ) {
      tree->root = i ;
   } else {
      sib[prev] = i ;
   }
   sib[i] = j ;
}
 
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   left-justify a tree by a metric
   children are linked in ascending order of their metric

   created -- 96jun23, cca
   ------------------------------------------------------
*/
void
Tree_leftJustifyD (
   Tree   *tree,
   DV     *metricDV
) {
int      i, j, k, n, nexti, prev ;
int      *fch, *par, *sib ;
double   *metric ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL 
   || (n = tree->n) <= 0 
   || metricDV == NULL 
   || n != DV_size(metricDV)
   || (metric = DV_entries(metricDV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_leftJustifyD(%p,%p)"
           "\n bad input\n", tree, metricDV) ;
   exit(-1) ;
}
par = tree->par ;
fch = tree->fch ;
sib = tree->sib ;
/*
   ----------------------------------------------------
   sort all children in decreasing order of metric size
   ----------------------------------------------------
*/
for ( k = 0 ; k < n ; k++ ) {
   for ( i = fch[k], fch[k] = -1 ; i != -1 ; i = nexti ) {
      nexti = sib[i] ;
      for ( j = fch[k], prev = -1 ; j != -1 ; j = sib[j] ) {
         if ( metric[j] < metric[i] ) {
            break ;
         }
         prev = j ;
      }
      if ( prev == -1 ) {
         fch[k] = i ;
      } else {
         sib[prev] = i ;
      }
      sib[i] = j ;
   }
}
/*
   ---------------------------------------------
   sort roots in decreasing order of metric size
   ---------------------------------------------
*/
for ( i = tree->root, tree->root = -1, prev = -1 ; 
      i != -1 ; 
      i = nexti ) {
   nexti = sib[i] ;
   for ( j = tree->root, prev = -1 ; j != -1 ; j = sib[j] ) {
       if ( metric[j] < metric[i] ) {
          break ;
       }
      prev = j ;
   }
   if ( prev == -1 ) {
      tree->root = i ;
   } else {
      sib[prev] = i ;
   }
   sib[i] = j ;
}
 
return ; }

/*--------------------------------------------------------------------*/
