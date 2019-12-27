/*  util.c  */

#include "../SolveMap.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- return the owner of block (rowid, colid).

   created -- 98mar19, cca
   ----------------------------------------------------
*/
int
SolveMap_owner (
   SolveMap   *solvemap,
   int        rowid,
   int        colid
) {
int   ii, loc, J, K, nblock ;
int   *colids, *rowids, *map ;
/*
   ---------------
   check the input
   ---------------
*/
if (  solvemap == NULL 
   || rowid < 0 || rowid >= solvemap->nfront
   || colid < 0 || colid >= solvemap->nfront ) {
   fprintf(stderr, "\n fatal error in SolveMap_owner(%p,%d,%d)"
           "\n bad input\n", solvemap, rowid, colid) ;
   exit(-1) ;
}
if ( rowid == colid ) {
/*
   --------------
   diagonal block
   --------------
*/
   return(solvemap->owners[rowid]) ;
} else if ( rowid > colid && solvemap->symmetryflag > 0 ) {
/*
   ---------------------------------
   lower block (K,J) = (rowid,colid)
   ---------------------------------
*/
   nblock = solvemap->nblockLower ;
   rowids = solvemap->rowidsLower ;
   colids = solvemap->colidsLower ;
   map    = solvemap->mapLower    ;
   K      = rowid ;
   J      = colid ;
   loc = IVlocateViaBinarySearch(nblock, colids, J) ;
   if ( loc == -1 ) {
      return(-1) ;
   }
   for ( ii = loc ; ii >= 0 ; ii-- ) {
      if ( colids[ii] == J && rowids[ii] == K ) {
         return(map[ii]) ;
      }
   }
   for ( ii = loc + 1 ; ii < nblock ; ii++ ) {
      if ( colids[ii] == J && rowids[ii] == K ) {
         return(map[ii]) ;
      }
   }
} else {
/*
   -----------------
   upper block (J,K)
   -----------------
*/
   nblock = solvemap->nblockUpper ;
   rowids = solvemap->rowidsUpper ;
   colids = solvemap->colidsUpper ;
   map    = solvemap->mapUpper    ;
   if ( rowid > colid ) {
      J = colid ; K = rowid ;
   } else {
      J = rowid ; K = colid ;
   }
   loc = IVlocateViaBinarySearch(nblock, rowids, J) ;
   if ( loc == -1 ) {
      return(-1) ;
   }
   for ( ii = loc ; ii >= 0 ; ii-- ) {
      if ( rowids[ii] == J && colids[ii] == K ) {
         return(map[ii]) ;
      }
   }
   for ( ii = loc +1 ; ii < nblock ; ii++ ) {
      if ( rowids[ii] == J && colids[ii] == K ) {
         return(map[ii]) ;
      }
   }
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- return an IVL object whose list K contains all
      processes who do not own U(K,K) but own a U(J,K) for some J.

   if myid == -1 then
      the entire IVL object is created and returned
   else 
      only the portion of the IVL object pertinent 
      to myid is created and returned
   endif

   created -- 98may24, cca
   ---------------------------------------------------------------
*/
IVL *
SolveMap_upperSolveIVL (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int   count, K, loc, nblock, nfront, nproc, q ;
int   *colids, *heads, *link, *list, *map, *mark, *owners, *rowids ;
IVL   *solveIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_upperSolveIVL(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
nfront = solvemap->nfront      ;
nproc  = solvemap->nproc       ;
nblock = solvemap->nblockUpper ;
colids = solvemap->colidsUpper ;
rowids = solvemap->rowidsUpper ;
map    = solvemap->mapUpper    ;
owners = solvemap->owners      ;
/*
   ------------------------------------
   link the (J,K,map(J,K)) triples by K
   ------------------------------------
*/
heads = IVinit(nfront, -1) ;
link  = IVinit(nblock, -1) ;
for ( loc = 0 ; loc < nblock ; loc++ ) {
   K         = colids[loc] ;
   link[loc] = heads[K]    ;
   heads[K]  = loc         ;
}
list = IVinit(nproc, -1) ;
mark = IVinit(nproc, -1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n linked triples by columns of U") ;
   for ( K = 0 ; K < nfront ; K++ ) {
      if ( heads[K] != -1 ) {
         fprintf(msgFile, "\n %d :", K) ;
         for ( loc = heads[K] ; loc != -1 ; loc = link[loc] ) {
            fprintf(msgFile, " <%d,%d>", rowids[loc], map[loc]) ;
         }
      }
   }
}
/*
   -------------------------------
   initialize the solve IVL object
   -------------------------------
*/
solveIVL = IVL_new() ;
IVL_init1(solveIVL, IVL_CHUNKED, nfront) ;
/*
   -------------------------------------------------------------
   fill the solve IVL object
   list K contains the process ids that do not own K but use X_K
   -------------------------------------------------------------
*/
for ( K = 0 ; K < nfront ; K++ ) {
   if ( myid == -1 || owners[K] == myid ) {
      mark[owners[K]] = K ;
      count = 0 ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n list for %d :", K) ;
      }
      for ( loc = heads[K] ; loc != -1 ; loc = link[loc] ) {
         q = map[loc] ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, " <%d,%d>", rowids[loc], q) ;
         }
         if ( mark[q] != K ) {
            mark[q]       = K ;
            list[count++] = q ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "*") ;
            }
         }
      } 
      if ( count > 0 ) {
         IVqsortUp(count, list) ;
         IVL_setList(solveIVL, K, count, list) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(heads) ;
IVfree(link) ;
IVfree(list) ;
IVfree(mark) ;

return(solveIVL) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- return an IV object whose entry J contains 
     the number of all processes who do not own U(J,J) 
     but own a U(J,K) for some J.

   if myid == -1 then
      all entries in the vector are filled
   else 
      all those entries pertinent to myid are filled
   endif

   created -- 98mar19, cca
   -----------------------------------------------------
*/
IV *
SolveMap_upperAggregateIV (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int   count, J, loc, nblock, nfront, nproc, q ;
int   *aggcounts, *colids, *heads, *link, *map, 
      *mark, *owners, *rowids ;
IV    *aggIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_upperAggregateIVL(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
nfront = solvemap->nfront      ;
nproc  = solvemap->nproc       ;
nblock = solvemap->nblockUpper ;
colids = solvemap->rowidsUpper ;
rowids = solvemap->rowidsUpper ;
map    = solvemap->mapUpper    ;
owners = solvemap->owners      ;
/*
   ------------------------
   initialize the IV object
   ------------------------
*/
aggIV = IV_new() ;
IV_init(aggIV, nfront, NULL) ;
aggcounts = IV_entries(aggIV) ;
IVzero(nfront, aggcounts) ;
/*
   ------------------------------------
   link the (J,K,map(J,K)) triples by J
   ------------------------------------
*/
heads = IVinit(nfront, -1) ;
link  = IVinit(nblock, -1) ;
for ( loc = 0 ; loc < nblock ; loc++ ) {
   J         = rowids[loc] ;
   link[loc] = heads[J]    ;
   heads[J]  = loc         ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n linked triples by rows of U") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( heads[J] != -1 ) {
         fprintf(msgFile, "\n %d :", J) ;
         for ( loc = heads[J] ; loc != -1 ; loc = link[loc] ) {
            fprintf(msgFile, " <%d,%d>", colids[loc], map[loc]) ;
         }
      }
   }
}
/*
   -------------------------
   fill the aggcounts vector
   -------------------------
*/
mark = IVinit(nproc, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( myid == -1 || owners[J] == myid ) {
      mark[owners[J]] = J ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n list for %d :", J) ;
      }
      for ( loc = heads[J], count = 0 ; loc != -1 ; loc = link[loc] ) {
         q = map[loc] ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, " <%d,%d>", colids[loc], q) ;
         }
         if ( mark[q] != J ) {
            mark[q] = J ;
            count++ ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "*") ;
            }
         }
      }
      aggcounts[J] = count ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(heads) ;
IVfree(link)  ;
IVfree(mark)  ;

return(aggIV) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- return an IV object whose J'th entry contains 
      the number of processes who do not own L(J,J) but own 
      a L(J,I) for some I < J.

   if myid == -1 then
      all entries in the vector are filled
   else 
      all those entries pertinent to myid are filled
   endif

   created -- 98mar20, cca
   --------------------------------------------------------
*/
IV *
SolveMap_lowerAggregateIV (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int   count, J, loc, nblock, nfront, nproc, q ;
int   *aggcounts, *colids, *heads, *link, *map, 
      *mark, *owners, *rowids ;
IV    *aggIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_lowerAggregateIV(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
nfront = solvemap->nfront ;
nproc  = solvemap->nproc  ;
if ( solvemap->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   nblock = solvemap->nblockLower ;
   colids = solvemap->colidsLower ;
   rowids = solvemap->rowidsLower ;
   map    = solvemap->mapLower    ;
} else {
   nblock = solvemap->nblockUpper ;
   colids = solvemap->rowidsUpper ;
   rowids = solvemap->colidsUpper ;
   map    = solvemap->mapUpper    ;
}
owners = solvemap->owners ;
/*
   ------------------------------------
   link the (J,I,map(J,I)) triples by J
   ------------------------------------
*/
heads = IVinit(nfront, -1) ;
link  = IVinit(nblock, -1) ;
for ( loc = 0 ; loc < nblock ; loc++ ) {
   J         = rowids[loc] ;
   link[loc] = heads[J]    ;
   heads[J]  = loc         ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n linked triples by rows of L or U^T") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( heads[J] != -1 ) {
         fprintf(msgFile, "\n %d :", J) ;
         for ( loc = heads[J] ; loc != -1 ; loc = link[loc] ) {
            fprintf(msgFile, " <%d,%d>", colids[loc], map[loc]) ;
         }
      }
   }
}
mark = IVinit(nproc, -1) ;
/*
   ----------------------------------
   initialize the aggregate IV object
   ----------------------------------
*/
aggIV = IV_new() ;
IV_init(aggIV, nfront, NULL) ;
aggcounts = IV_entries(aggIV) ;
IVzero(nfront, aggcounts) ;
/*
   -------------------------------------------------------------
   fill the aggregate IV object
   aggcounts[J] = # of processors besides owner of L_{J,J}
   that own some L_{J,I}
   -------------------------------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( myid == -1 || owners[J] == myid ) {
      mark[owners[J]] = J ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n list for %d :", J) ;
      }
      for ( loc = heads[J], count = 0 ; loc != -1 ; loc = link[loc] ) {
         q = map[loc] ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, " <%d,%d>", colids[loc], q) ;
         }
         if ( mark[q] != J ) {
            mark[q] = J ;
            count++ ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "*") ;
            }
         }
      } 
      aggcounts[J] = count ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(heads) ;
IVfree(link) ;
IVfree(mark) ;

return(aggIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- return an IVL object whose list J contains 
      the processes who do not own L(J,J) but own a L(K,J) 
      for some K > J.

   if myid == -1 then
      the entire IVL object is created and returned
   else 
      only the portion of the IVL object pertinent 
      to myid is created and returned
   endif

   created -- 98may24, cca
   -------------------------------------------------------
*/
IVL *
SolveMap_lowerSolveIVL (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int    count, I, J, loc, nblock, nfront, nproc, q ;
int    *colids, *heads, *link, *list, *map, *mark, *owners, *rowids ;
IVL    *solveIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_lowerSolveIVL(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
nfront = solvemap->nfront ;
nproc  = solvemap->nproc  ;
if ( solvemap->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   nblock = solvemap->nblockLower ;
   rowids = solvemap->rowidsLower ;
   colids = solvemap->colidsLower ;
   map    = solvemap->mapLower    ;
} else {
   nblock = solvemap->nblockUpper ;
   rowids = solvemap->colidsUpper ;
   colids = solvemap->rowidsUpper ;
   map    = solvemap->mapUpper    ;
}
owners = solvemap->owners ;
/*
   ------------------------
   initialize the IV object
   ------------------------
*/
solveIVL = IVL_new() ;
IVL_init1(solveIVL, IVL_CHUNKED, nfront) ;
/*
   ------------------------------------
   link the (J,I,map(J,I)) triples by I
   ------------------------------------
*/
heads = IVinit(nfront, -1) ;
link  = IVinit(nblock, -1) ;
for ( loc = 0 ; loc < nblock ; loc++ ) {
   I         = colids[loc] ;
   link[loc] = heads[I]    ;
   heads[I]  = loc         ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n linked triples by columns of L or U^T") ;
   for ( I = 0 ; I < nfront ; I++ ) {
      if ( heads[I] != -1 ) {
         fprintf(msgFile, "\n %d :", I) ;
         for ( loc = heads[I] ; loc != -1 ; loc = link[loc] ) {
            fprintf(msgFile, " <%d,%d>", rowids[loc], map[loc]) ;
         }
      }
   }
}
/*
   -------------------------
   fill the solve IVL object
   -------------------------
*/
list = IVinit(nproc, -1) ;
mark = IVinit(nproc, -1) ;
for ( I = 0 ; I < nfront ; I++ ) {
   if ( myid == -1 || owners[I] == myid ) {
      mark[owners[I]] = I ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n list for %d :", I) ;
      }
      for ( loc = heads[I], count = 0 ; loc != -1 ; loc = link[loc] ) {
         q = map[loc] ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, " <%d,%d>", rowids[loc], q) ;
         }
         if ( mark[q] != I ) {
            mark[q] = I ;
            list[count++] = q ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "*") ;
            }
         }
      }
      if ( count > 0 ) {
         IVqsortUp(count, list) ;
         IVL_setList(solveIVL, I, count, list) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(heads) ;
IVfree(link)  ;
IVfree(mark)  ;
IVfree(list)  ;

return(solveIVL) ; }

/*--------------------------------------------------------------------*/
