/*  maximizeGain.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- 

   given a gain value assigned to each node,
   find a set of nodes, no two in a child-ancestor
   relationship, that maximizes the total gain.

   this problem arises in finding the optimal domain/schur 
   complement partition for a semi-implicit factorization.

   created -- 98jun20, cca
   -------------------------------------------------------
*/
IV *
Tree_maximizeGainIV (
   Tree   *tree,
   IV     *gainIV,
   int    *ptotalgain,
   int    msglvl,
   FILE   *msgFile
) {
char   *mark ;
int    I, J, K, n, ndom, sum, totalgain ;
int    *compids, *fch, *gain, *par, *sib, *subtreeGain ;
IV     *compidsIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || gainIV == NULL || ptotalgain == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in Tree_maximizeGainIV()"
           "\n bad input\n") ;
   exit(-1) ;
}
n   = tree->n   ;
par = tree->par ;
fch = tree->fch ;
sib = tree->sib ;
if ( n != IV_size(gainIV) ) {
   fprintf(stderr, "\n fatal error in Tree_maximizeGainIV()"
           "\n tree size = %d, gain size = %d", 
           tree->n, IV_size(gainIV)) ;
   exit(-1) ;
}
gain = IV_entries(gainIV) ;
/*
   --------------------------------------------------
   compute the subtree gains and mark candidate roots
   --------------------------------------------------
*/
mark = CVinit(n, 'N') ;
subtreeGain = IVinit(n, 0) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   if ( fch[J] == -1 ) {
      mark[J] = 'R' ;
      subtreeGain[J] = gain[J] ;
   } else {
      for ( I = fch[J], sum = 0 ; I != -1 ; I = sib[I] ) {
         sum += subtreeGain[I] ;
      }
      if ( gain[J] >= sum ) {
         subtreeGain[J] = gain[J] ;
         mark[J] = 'R' ;
      } else {
         subtreeGain[J] = sum ;
      }
   }
}
/*
   ----------------------
   compute the total gain
   ----------------------
*/
for ( J = tree->root, totalgain = 0 ; J != -1 ; J = sib[J] ) {
   totalgain += subtreeGain[J] ;
}
*ptotalgain = totalgain ;
/*
   ----------------------------------------------
   create and initialize the component ids vector
   ----------------------------------------------
*/
compidsIV = IV_new() ;
IV_init(compidsIV, n, NULL) ;
IV_fill(compidsIV, 0) ;
compids = IV_entries(compidsIV) ;
/*
   ----------------------------------------------
   fix the component ids of the nodes in the tree
   ----------------------------------------------
*/
for ( J = Tree_preOTfirst(tree), ndom = 0 ; 
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
   if ( mark[J] == 'R' ) {
      if ( (K = par[J]) != -1 && compids[K] != 0 ) {
         compids[J] = compids[K] ;
      } else {
         compids[J] = ++ndom ;
      }
   } else if ( (K = par[J]) != -1 ) {
      compids[J] = compids[K] ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(subtreeGain) ;
CVfree(mark) ;

return(compidsIV) ; }

/*--------------------------------------------------------------------*/
