/*  maps.c  */

#include "../SolveMap.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- map the off diagonal blocks 
      to processes in a random fashion

   created -- 98mar19, cca
   --------------------------------------
*/
void
SolveMap_randomMap (
   SolveMap   *solvemap,
   int        symmetryflag,
   IVL        *upperBlockIVL,
   IVL        *lowerBlockIVL,
   int        nproc,
   IV         *ownersIV,
   int        seed,
   int        msglvl,
   FILE       *msgFile
) {
Drand   drand ;
int     ii, J, K, loc, nadj, nblockLower, nblockUpper, 
        nfront, proc ;
int     *adj, *colids, *map, *owners, *rowids ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || symmetryflag < 0 
   || upperBlockIVL == NULL || ownersIV == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SolveMap_randomMap(%p,%d,%p,%p,%p,%d)"
           "\n bad input\n",
           solvemap, symmetryflag, upperBlockIVL, 
           lowerBlockIVL, ownersIV, seed) ;
   exit(-1) ;
}
nfront = IV_size(ownersIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n\n SolveMap_randomMap(): nfront = %d, nproc = %d",
           nfront, nproc) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   count the number of upper blocks that do not include U(J,J)
   -----------------------------------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n upperBlockIVL = %p", upperBlockIVL) ;
   fflush(msgFile) ;
}
nblockUpper = 0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   IVL_listAndSize(upperBlockIVL, J, &nadj, &adj) ;
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      if ( adj[ii] > J ) {
         nblockUpper++ ;
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nblockUpper = %d", nblockUpper) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   count the number of lower blocks that do not include L(J,J)
   -----------------------------------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n lowerBlockIVL = %p", lowerBlockIVL) ;
   fflush(msgFile) ;
}
nblockLower = 0 ;
if ( lowerBlockIVL != NULL ) {
   for ( J = 0 ; J < nfront ; J++ ) {
      IVL_listAndSize(lowerBlockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         if ( adj[ii] > J ) {
            nblockLower++ ;
         }
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nblockLower = %d", nblockLower) ;
   fflush(msgFile) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
SolveMap_init(solvemap, symmetryflag, nfront, 
              nproc, nblockUpper, nblockLower) ;
owners = SolveMap_owners(solvemap) ;
/*
   ----------------------
   fill the owners vector
   ----------------------
*/
IVcopy(nfront, owners, IV_entries(ownersIV)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n owners") ;
   IVfprintf(msgFile, nfront, owners) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_setUniform(&drand, 0, nproc) ;
/*
   ----------------------------------------
   map the upper blocks in a random fashion
   ----------------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n mapping upper blocks") ;
   fflush(msgFile) ;
}
rowids = SolveMap_rowidsUpper(solvemap) ;
colids = SolveMap_colidsUpper(solvemap) ;
map    = SolveMap_mapUpper(solvemap) ;
for ( J = loc = 0 ; J < nfront ; J++ ) {
   IVL_listAndSize(upperBlockIVL, J, &nadj, &adj) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n J = %d", J) ;
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    K = %d", adj[ii]) ;
         fflush(msgFile) ;
      }
      if ( (K = adj[ii]) > J ) {
         proc = (int) Drand_value(&drand) ;
         rowids[loc] =   J  ;
         colids[loc] =   K  ;
         map[loc]    = proc ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", map[%d] = %d", loc, map[loc]) ;
            fflush(msgFile) ;
         }
         loc++ ;
      }
   }
}
if ( lowerBlockIVL != NULL ) {
/*
   ----------------------------------------
   map the lower blocks in a random fashion
   ----------------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n mapping lower blocks") ;
      fflush(msgFile) ;
   }
   rowids = SolveMap_rowidsLower(solvemap) ;
   colids = SolveMap_colidsLower(solvemap) ;
   map    = SolveMap_mapLower(solvemap) ;
   for ( J = loc = 0 ; J < nfront ; J++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n J = %d", J) ;
         fflush(msgFile) ;
      }
      IVL_listAndSize(lowerBlockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    K = %d", adj[ii]) ;
            fflush(msgFile) ;
         }
         if ( (K = adj[ii]) > J ) {
            proc = (int) Drand_value(&drand) ;
            rowids[loc] =   K  ;
            colids[loc] =   J  ;
            map[loc]    = proc ;
            if ( msglvl > 2 ) {
               fprintf(msgFile, ", map[%d] = %d", loc, map[loc]) ;
               fflush(msgFile) ;
            }
            loc++ ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- map the off diagonal blocks to
      processes in a domain decomposition fashion

   created -- 98mar28, cca
   ----------------------------------------------
*/
void
SolveMap_ddMap (
   SolveMap   *solvemap,
   int        symmetryflag,
   IVL        *upperBlockIVL,
   IVL        *lowerBlockIVL,
   int        nproc,
   IV         *ownersIV,
   Tree       *tree,
   int        seed,
   int        msglvl,
   FILE       *msgFile
) {
char    *mark ;
Drand   drand ;
int     ii, I, J, K, loc, nadj, nblockLower, nblockUpper, 
        nfront, proc ;
int     *adj, *colids, *fch, *map, *owners, *rowids, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || symmetryflag < 0 
   || upperBlockIVL == NULL || ownersIV == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SolveMap_ddMap(%p,%d,%p,%p,%p,%d)"
           "\n bad input\n",
           solvemap, symmetryflag, upperBlockIVL, 
           lowerBlockIVL, ownersIV, seed) ;
   exit(-1) ;
}
nfront = IV_size(ownersIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n\n SolveMap_ddMap(): nfront = %d, nproc = %d",
           nfront, nproc) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   count the number of upper blocks that do not include U(J,J)
   -----------------------------------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n upperBlockIVL = %p", upperBlockIVL) ;
   fflush(msgFile) ;
}
nblockUpper = 0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   IVL_listAndSize(upperBlockIVL, J, &nadj, &adj) ;
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      if ( adj[ii] > J ) {
         nblockUpper++ ;
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nblockUpper = %d", nblockUpper) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   count the number of lower blocks that do not include L(J,J)
   -----------------------------------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n lowerBlockIVL = %p", lowerBlockIVL) ;
   fflush(msgFile) ;
}
nblockLower = 0 ;
if ( lowerBlockIVL != NULL ) {
   for ( J = 0 ; J < nfront ; J++ ) {
      IVL_listAndSize(lowerBlockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         if ( adj[ii] > J ) {
            nblockLower++ ;
         }
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nblockLower = %d", nblockLower) ;
   fflush(msgFile) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
SolveMap_init(solvemap, symmetryflag, nfront, 
              nproc, nblockUpper, nblockLower) ;
owners = SolveMap_owners(solvemap) ;
/*
   ----------------------
   fill the owners vector
   ----------------------
*/
IVcopy(nfront, owners, IV_entries(ownersIV)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n owners") ;
   IVfprintf(msgFile, nfront, owners) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------
   mark a node J in the tree as 'D' if it is in a domain
   (owners[J] = owners[I] for all I a descendent of J)
   and 'S' (for the schur complement) otherwise
   -----------------------------------------------------
*/
mark = CVinit(nfront, 'D') ;
fch  = tree->fch ;
sib  = tree->sib ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      if ( mark[I] != 'D' || owners[I] != owners[J] ) {
         mark[J] = 'S' ; break ;
      }
   }
}
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_setUniform(&drand, 0, nproc) ;
/*
   -------------------------------
   if J is in a domain
      map(J,K) to owners[J]
   else
      map(J,K) to a random process
   -------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n mapping upper blocks") ;
   fflush(msgFile) ;
}
rowids = SolveMap_rowidsUpper(solvemap) ;
colids = SolveMap_colidsUpper(solvemap) ;
map    = SolveMap_mapUpper(solvemap) ;
for ( J = loc = 0 ; J < nfront ; J++ ) {
   IVL_listAndSize(upperBlockIVL, J, &nadj, &adj) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n J = %d", J) ;
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    K = %d", adj[ii]) ;
         fflush(msgFile) ;
      }
      if ( (K = adj[ii]) > J ) {
         if ( mark[J] == 'D' ) {
            proc = owners[J] ;
         } else {
            proc = (int) Drand_value(&drand) ;
         }
         rowids[loc] =   J  ;
         colids[loc] =   K  ;
         map[loc]    = proc ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", map[%d] = %d", loc, map[loc]) ;
            fflush(msgFile) ;
         }
         loc++ ;
      }
   }
}
if ( symmetryflag == SPOOLES_NONSYMMETRIC ) {
/*
   -------------------------------
   if J is in a domain
      map(K,J) to owners[J]
   else
      map(K,J) to a random process
   -------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n mapping lower blocks") ;
      fflush(msgFile) ;
   }
   rowids = SolveMap_rowidsLower(solvemap) ;
   colids = SolveMap_colidsLower(solvemap) ;
   map    = SolveMap_mapLower(solvemap) ;
   for ( J = loc = 0 ; J < nfront ; J++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n J = %d", J) ;
         fflush(msgFile) ;
      }
      IVL_listAndSize(lowerBlockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    K = %d", adj[ii]) ;
            fflush(msgFile) ;
         }
         if ( (K = adj[ii]) > J ) {
            if ( mark[J] == 'D' ) {
               proc = owners[J] ;
            } else {
               proc = (int) Drand_value(&drand) ;
            }
            rowids[loc] =   K  ;
            colids[loc] =   J  ;
            map[loc]    = proc ;
            if ( msglvl > 2 ) {
               fprintf(msgFile, ", map[%d] = %d", loc, map[loc]) ;
               fflush(msgFile) ;
            }
            loc++ ;
         }
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
CVfree(mark) ;

return ; }

/*--------------------------------------------------------------------*/
