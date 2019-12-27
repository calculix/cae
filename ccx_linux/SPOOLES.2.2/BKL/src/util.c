/*  util.c  */

#include "../BKL.h"
#include "../../Drand.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
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
) {
int     ireg, ndom, nreg ;
int     *colors ;
Drand   drand ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || bkl->bpg == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_setRandomColors(%p,%d)"
           "\n bad input\n", bkl, seed) ;
   exit(-1) ;
}
ndom   = bkl->ndom   ;
nreg   = bkl->nreg   ;
colors = bkl->colors ;
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setUniform(&drand, 0.0, 1.0) ;
if ( seed > 0 ) {
/*
   --------------------------
   set the random number seed
   --------------------------
*/
   Drand_setSeed(&drand, seed) ;
}
/*
   -----------------
   color the domains
   -----------------
*/
for ( ireg = 0 ; ireg < ndom ; ireg++ ) {
   colors[ireg] = (Drand_value(&drand) < 0.5) ? 1 : 2 ;
}
/*
   --------------------------------------------
   color the segments and set the color weights
   --------------------------------------------
*/
BKL_setColorWeights(bkl) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   set the component weights. 
   note, we assume the domain colors are set

   created  -- 95oct07, cca
   modified -- 95dec07, cca
      error checking inserted for colors
   -----------------------------------------
*/
void
BKL_setColorWeights (
   BKL   *bkl
) {
int   c, ireg ;
int   *colors ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_setColorsWeights(%p)"
           "\n bad input\n", bkl) ;
   exit(-1) ;
}
colors = bkl->colors ;
bkl->cweights[0] = bkl->cweights[1] = bkl->cweights[2] = 0 ;
/*
   ---------------------
   check out the domains
   ---------------------
*/
for ( ireg = 0 ; ireg < bkl->ndom ; ireg++ ) {
   if ( (c = colors[ireg]) < 1 || 2 < c ) {
      fprintf(stderr, "\n fatal error in BKL_setColorWeights(%p)"
              "\n region %d has color %d", bkl, ireg, c) ;
      exit(-1) ;
   }
   bkl->cweights[c] += bkl->regwghts[ireg] ;
}
/*
   ------------------
   color the segments
   ------------------
*/
for ( ireg = bkl->ndom ; ireg < bkl->nreg ; ireg++ ) {
   if ( (c = BKL_segColor(bkl, ireg)) < 0 || 2 < c ) {
      fprintf(stderr, "\n fatal error in BKL_setColorWeights(%p)"
              "\n region %d has color %d", bkl, ireg, c) ;
      exit(-1) ;
   }
   colors[ireg] = c ;
   bkl->cweights[c] += bkl->regwghts[ireg] ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
int   color, ii, size ;
int   *adj, *colors ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || iseg < bkl->ndom || iseg >= bkl->nreg ) {
   fprintf(stderr, "\n fatal error in BKL_segColor(%p,%d)"
           "\n bad input\n", bkl, iseg) ;
   exit(-1) ;
}
colors = bkl->colors ;
/*
   ----------------------------------------------------
   loop over adjacent domans, 
   break if adjacent to two differently colored domains
   ----------------------------------------------------
*/
Graph_adjAndSize(bkl->bpg->graph, iseg, &size, &adj) ;
color = 0 ;
if ( size > 0 ) {
   color = colors[adj[0]] ;
   for ( ii = 1 ; ii < size ; ii++ ) {
      if ( color != colors[adj[ii]] ) {
         color = 0 ;
         break ;
      }
   }
}
return(color) ; }
   
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   flip the domain

   created  -- 95oct07, cca
   modified -- 95dec07, cca
      simple mod to segment loop made
   ----------------------------------
*/
void
BKL_flipDomain (
   BKL   *bkl,
   int   idom
) {
int   ii, iseg, newcolor, oldcolor, size, wdom, wseg ;
int   *adj, *colors, *regwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || idom < 0 || idom >= bkl->ndom ) {
   fprintf(stderr, "\n fatal error in BKL_flipDomain(%p,%d)"
           "\n bad input\n", bkl, idom) ;
   exit(-1) ;
}
colors   = bkl->colors   ;
regwghts = bkl->regwghts ;
switch ( (oldcolor = colors[idom]) ) {
case 1 :
   newcolor = 2 ;
   break ;
case 2 :
   newcolor = 1 ;
   break ;
default :
   fprintf(stderr, "\n fatal error in BKL_flipDomain(%p,%d)"
           "\n colors[%d] = %d\n", bkl, idom, idom, colors[idom]) ;
   exit(-1) ;
}
colors[idom] = newcolor ;
/*
   --------------------------------------
   adjust color weights for moving domain
   --------------------------------------
*/
wdom = regwghts[idom] ;
bkl->cweights[oldcolor] -= wdom ;
bkl->cweights[newcolor] += wdom ;
/*
   -------------------------------
   loop over the adjacent segments
   -------------------------------
*/
Graph_adjAndSize(bkl->bpg->graph, idom, &size, &adj) ;
for ( ii = 0 ; ii < size ; ii++ ) {
   iseg = adj[ii] ;
   wseg = regwghts[iseg] ;
#if MYDEBUG > 0
   fprintf(stdout, "\n checking out segment %d, weight %d",
           iseg, wseg) ;
#endif
   if ( (oldcolor = colors[iseg]) 
          != (newcolor = BKL_segColor(bkl, iseg)) ) {
      bkl->cweights[oldcolor] -= wseg ;
      bkl->cweights[newcolor] += wseg ;
      colors[iseg] = newcolor ;
   }
}
/*
   -----------------------------
   increment the number of flips
   -----------------------------
*/
bkl->nflips++ ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   return the next domain to flip
   in a grey code sequence

   created  -- 95oct07, cca
   modified -- 95dec07, cca
      error message inserted
   ------------------------------
*/
int
BKL_greyCodeDomain (
   BKL   *bkl,
   int   count
) {
int   chk, idom, res ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_greyCodeDomain(%p)"
           "\n bad input\n", bkl) ;
   exit(-1) ;
}
 
for ( idom = 0, res = 1, chk = 2 ;
      /* no test */ ;
      idom++, res = chk, chk *= 2 ) {
   if ( count % chk == res ) {
      return(idom) ;
   }
}
fprintf(stderr, "\n fatal error in BKL_greyCodeDomain(%p,%d)"
        "\n should never have reached this point\n", bkl, count) ;
exit(-1) ;
 
return(-2) ; }

/*--------------------------------------------------------------------*/
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
                we call set the random number generator.
   domcolors -- vector of domain colors, used when flag == 6

   created  -- 95oct11, cca
   modified -- 95nov30, cca
      switch for one and two domains inserted.
   ----------------------------------------------------------------
*/
float
BKL_setInitPart (
   BKL   *bkl,
   int   flag,
   int   seed,
   int   domcolors[]
) {
BPG     *bpg ;
float   cost ;
int     dom, dom2, dsize, idom, ii, jj, last, ndom, now, root, 
        seg, ssize ;
int     *colors, *cweights, *dadj, *list, *mark, *sadj ;
/*
   ---------------
   check the input
   ---------------
*/
if (  bkl == NULL 
   || flag < 1 || 6 < flag || (flag == 6 && domcolors == NULL) ) { 
   fprintf(stderr, "\n fatal error in BKL_setInitPart(%p,%d,%d,%p)"
           "\n bad input\n", bkl, flag, seed, domcolors) ;
   exit(-1) ;
}
bpg      = bkl->bpg      ;
ndom     = bkl->ndom     ;
colors   = bkl->colors   ;
cweights = bkl->cweights ;
if ( ndom == 1 ) {
   colors[0] = 1 ;
   BKL_setColorWeights(bkl) ;
} else if ( ndom == 2 ) {
   colors[0] = 1 ;
   colors[1] = 2 ;
   BKL_setColorWeights(bkl) ;
} else {
/*
   ---------------------
   switch over the cases
   ---------------------
*/
   switch ( flag ) {
   case 1 : {
      Drand   drand ;
/*
      -------------
      random colors
      -------------
*/
      Drand_setDefaultFields(&drand) ;
      Drand_init(&drand) ;
      Drand_setUniform(&drand, 0.0, 1.0) ;
      if ( seed > 0 ) {
/*
         --------------------------
         set the random number seed
         --------------------------
*/
         Drand_setSeed(&drand, seed) ;
      }
      for ( idom = 0 ; idom < ndom ; idom++ ) {
         colors[idom] = (Drand_value(&drand) < 0.5) ? 1 : 2 ;
      }
      BKL_setColorWeights(bkl) ;
      break ; }
   case 2 :
   case 3 :
/*
      --------------------------------------------------------
      one domain colored black = 1, the rest colored white = 2
      domain is specified (flag = 2) or pseudoperipheral
      (flag = 3)
      --------------------------------------------------------
*/
      IVfill(ndom, colors, 2) ;
      if ( flag == 2 ) {
         colors[seed % ndom] = 1 ;
      } else {
         root = BPG_pseudoperipheralnode(bkl->bpg, seed % ndom) ;
         colors[root] = 1 ;
      }
      BKL_setColorWeights(bkl) ;
      break ;
   case 4 :
   case 5 :
/*
      ----------------------------------------------
      split using a random or pseudoperipheral node 
      as a root of a breadth first traversal
      1. color all domains white
      2. color segments and get color weights
      3. fill the list with the seed domain
      4. do a breadth first traversal of the domains
         5. flip the domain
         6. if |B| >= |W| then break
         7. add unmarked neighboring domains to list
      ----------------------------------------------
*/
      IVfill(ndom, colors, 2) ;
      BKL_setColorWeights(bkl) ;
      list = IVinit(ndom, -1) ;
      mark = IVinit(ndom, -1) ;
      if ( flag == 4 ) {
         list[0] = seed % ndom ;
      } else {
         list[0] = BPG_pseudoperipheralnode(bkl->bpg, seed % ndom) ;
      }
      now = last = 0 ;
      mark[list[0]] = 1 ;
      while ( now <= last ) {
         dom = list[now++] ;
         BKL_flipDomain(bkl, dom) ;
         if ( cweights[1] >= cweights[2] ) {
            break ;
         }
         Graph_adjAndSize(bpg->graph, dom, &dsize, &dadj) ;
         for ( ii = 0 ; ii < dsize ; ii++ ) {
            seg = dadj[ii] ;
            Graph_adjAndSize(bpg->graph, seg, &ssize, &sadj) ;
            for ( jj = 0 ; jj < ssize ; jj++ ) {
               dom2 = sadj[jj] ;
               if ( mark[dom2] == -1 ) {
                  if ( last == ndom - 1 ) {
                     fprintf(stderr, 
                        "\n fatal error in BKL_setInitPart(%p,%d,%d,%p)"
                           "\n list[] size exceeded\n", 
                             bkl, flag, seed, domcolors) ;
                     exit(-1) ;
                  }
                  mark[dom2] = 1 ;
                  list[++last] = dom2 ;
               }
            }
         }
      }
      IVfree(list) ;
      IVfree(mark) ;
      BKL_setColorWeights(bkl) ;
      break ;
   case 6 :
/*
      ------------------
      copy domain colors
      ------------------
*/
      IVcopy(ndom, colors, domcolors) ;
      BKL_setColorWeights(bkl) ;
      break ;
   }
}
cost = BKL_evalfcn(bkl) ;    

return(cost) ; }

/*--------------------------------------------------------------------*/
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
) {
int   ii, size ;
int   *adj, *colors ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || dom < 0 || dom >= bkl->ndom ) {
   fprintf(stderr, "\n fatal error in BKL_domAdjToSep(%p,%d)"
           "\n bad input\n", bkl, dom) ;
   exit(-1) ;
}
colors = bkl->colors ;

Graph_adjAndSize(bkl->bpg->graph, dom, &size, &adj) ;
for ( ii = 0 ; ii < size ; ii++ ) {
   if ( colors[adj[ii]] == 0 ) {
/*
      -------------------------------------
      segment is on the separator, return 1
      -------------------------------------
*/
      return(1) ;
   }
}
return(0) ; }

/*--------------------------------------------------------------------*/
