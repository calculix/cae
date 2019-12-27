/*  maxFlow.c  */

#include "../BPG.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   a Cell structure represent an edge (x,y)

   rxy   -- residual for edge
   fxy   -- flow for edge
   xnext -- next Cell in list for x
   ynext -- next Cell in list for y
   ----------------------------------------
*/
typedef struct _Cell   Cell ;
struct _Cell {
   int    x ;
   int    y ;
   int    rxy ;
   int    fxy ;
   Cell   *xnext ;
   Cell   *ynext ;
} ;
/*--------------------------------------------------------------------*/
/*
   -----------------
   static prototypes
   -----------------
*/
static void
setupCells(BPG *bpg, Cell ***pheads, Cell **pcells, 
           int msglvl, FILE *msgFile) ;
static int
findFlowAugmentingPath(BPG *bpg, Cell *heads[], int  xexp, int exp[],
                       int par[], int labels[], int list[], int mark[],
                       int tag, int msglvl, FILE *msgFile) ;
static void
augmentPath(BPG *bpg, int xexp, int yexp, int nvexp[], Cell *heads[], 
            int labels[], int par[], int msglvl, FILE *msgFile) ;
static void
setDMflags(BPG *bpg, Cell *heads[], int nvexp[], int list[], int mark[],
           int dmflags[], int stats[], int msglvl, FILE *msgFile) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   find the DM decomposition of a weighted bipartite graph
   using a simple max flow algorithm

   bpg     -- BPG bipartite graph object
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R

   created -- 96mar08, cca
   -------------------------------------------------------
*/
void
BPG_DMviaMaxFlow (
   BPG    *bpg,
   int    dmflags[],
   int    stats[],
   int    msglvl,
   FILE   *msgFile 
) {
int    ierr, nX, nY, tag, x, xexp, y, yexp ;
int    *labels, *list, *mark, *nvexp, *par, *vwghts ;
Cell   *cell, *cells ;
Cell   **heads ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || dmflags == NULL || stats == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in BPG_DMviaMaxFlow(%p,%p,%p,%d,%p)"
           "\n bad input\n", bpg, dmflags, stats, msglvl, msgFile) ;
   exit(-1) ;
}
nX = bpg->nX ;
nY = bpg->nY ;
/*
   ---------------------------------------------
   set up the Cell structures and working arrays
   ---------------------------------------------
*/
setupCells(bpg, &heads, &cells, msglvl, msgFile) ;
nvexp = IVinit(nX + nY, 1) ;
if ( (vwghts = bpg->graph->vwghts) != NULL ) {
   IVcopy(nX + nY, nvexp, vwghts) ;
}
list   = IVinit(nX + nY, -1) ;
mark   = IVinit(nX + nY, -1) ;
par    = IVinit(nX + nY, -1) ;
labels = IVinit(nX + nY, -1) ;
/*
   ------------------------
   loop over the x vertices
   ------------------------
*/
tag = 0 ;
for ( xexp = 0 ; xexp < nX ; xexp++ ) {
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n checking out x node %d, exposure %d", 
              xexp, nvexp[xexp]) ;
   }
   while ( nvexp[xexp] > 0 ) {
      yexp = findFlowAugmentingPath(bpg, heads, xexp, nvexp,
                                    par, labels, list, mark, tag, 
                                    msglvl, msgFile) ;
      tag++ ;
      if ( yexp == -1 ) {
/*
         -------------------------------------------------------
         no flow augmenting path has been found, node is exposed
         -------------------------------------------------------
*/
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, "\n node %d has exposure %d", 
                    xexp, nvexp[xexp]) ;
         }
         break ;
      } else {
/*
         ----------------------------------------
         flow augmenting path found, augment path
         ----------------------------------------
*/
         augmentPath(bpg, xexp, yexp, nvexp, heads, labels, par, 
                     msglvl, msgFile) ;
      }
   }
}
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n after maximum flow found") ;
   for ( x = 0 ; x < nX ; x++ ) {
      fprintf(msgFile, "\n X node %d : ", x) ;
      for ( cell = heads[x] ; cell != NULL ; cell = cell->xnext ) {
         fprintf(msgFile, "\n (%3d,%3d,%3d)",
                 cell->y, cell->fxy, cell->rxy) ;
      }
   }
   for ( y = nX ; y < nX + nY ; y++ ) {
      fprintf(msgFile, "\n Y node %d : ", y) ;
      for ( cell = heads[y] ; cell != NULL ; cell = cell->ynext ) {
         fprintf(msgFile, "\n (%3d,%3d,%3d)",
                 cell->x, cell->fxy, cell->rxy) ;
      }
   }
}
/*
   ----------------
   set the DM flags
   ----------------
*/
setDMflags(bpg, heads, nvexp, list, mark, 
           dmflags, stats, msglvl, msgFile) ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n X dmflags: ") ;
   IVfp80(msgFile, nX, dmflags, 12, &ierr) ;
   fprintf(msgFile, "\n Y dmflags: ") ;
   IVfp80(msgFile, nY, dmflags + nX, 12, &ierr) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(nvexp)  ;
IVfree(list)   ;
IVfree(mark)   ;
IVfree(par)    ;
IVfree(labels) ;
FREE(heads)    ;
FREE(cells)    ;

return ; }
    
/*--------------------------------------------------------------------*/
/*
   ------------
   set up cells
   ------------
*/
static void
setupCells (
   BPG    *bpg,
   Cell   ***pheads,
   Cell   **pcells,
   int    msglvl,
   FILE   *msgFile
) {
Cell    *cell, *cells ;
Cell    **heads ;
Graph   *graph ;
int     ii, ncell, nX, nY, x, xsize, y ;
int     *vwghts, *xadj ;
/*
   -------------------------
   count the number of cells
   -------------------------
*/
nX    = bpg->nX    ;
nY    = bpg->nY    ;
graph = bpg->graph ;
ncell = 0 ;
for ( x = 0 ; x < nX ; x++ ) {
   Graph_adjAndSize(graph, x, &xsize, &xadj) ;
   ncell += xsize ;
}
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n %d cells required", ncell) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------
   allocate the storage and set pointers for return
   ------------------------------------------------
*/
ALLOCATE(cells, struct _Cell,   ncell) ;
ALLOCATE(heads, struct _Cell *, nX + nY) ;
*pcells = cells ;
*pheads = heads ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n cells = %p", cells) ;
   fflush(msgFile) ;
}
for ( ii = 0 ; ii < nX + nY ; ii++ ) {
   heads[ii] = NULL ;
}
/*
   ----------------------------------
   set the vertex fields of the cells
   ----------------------------------
*/
for ( x = 0, cell = cells ; x < nX ; x++ ) {
   Graph_adjAndSize(graph, x, &xsize, &xadj) ;
   for ( ii = 0 ; ii < xsize ; ii++, cell++ ) {
      y = xadj[ii] ;
      cell->x = x ;
      cell->y = y ;
   }
}
/*
   --------------------------------
   set the link fields of the cells
   --------------------------------
*/
for ( ii = 0, cell = cells ; ii < ncell ; ii++, cell++ ) {
   x = cell->x ;
   cell->xnext = heads[x] ;
   heads[x] = cell ;
   y = cell->y ;
   cell->ynext = heads[y] ;
   heads[y] = cell ;
}
/*
   ---------------------------------------------
   set the residual and flow fields of the cells
   ---------------------------------------------
*/
if ( (vwghts = graph->vwghts) == NULL ) {
   for ( ii = 0, cell = cells ; ii < ncell ; ii++, cell++ ) {
      cell->rxy = 1 ;
      cell->fxy = 0 ;
   }
} else {
   for ( ii = 0, cell = cells ; ii < ncell ; ii++, cell++ ) {
      x = cell->x ;
      y = cell->y ;
      cell->rxy = (vwghts[x] < vwghts[y]) ? vwghts[x] : vwghts[y] ;
      cell->fxy = 0 ;
   }
}
if ( msglvl > 2 && msgFile != NULL ) {
   for ( ii = 0, cell = cells ; ii < ncell ; ii++, cell++ ) {
      fprintf(msgFile, 
           "\n cell %d, address %p : (x,y) = (%3d,%3d), r = %d, f = %d",
           ii, cell, cell->x, cell->y, cell->rxy, cell->fxy) ;
   }
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   search for a flow augmenting path starting at xexp

   created -- 96mar08, cca
   --------------------------------------------------
*/
static int
findFlowAugmentingPath ( 
   BPG    *bpg,
   Cell   *heads[],
   int    xexp,
   int    exp[],
   int    par[],
   int    labels[],
   int    list[],
   int    mark[],
   int    tag,
   int    msglvl,
   FILE   *msgFile
) {
Cell    *cell ;
Graph   *graph ;
int     last, now, nX, x, y ;
/*
   ---------------------------
   get dimensions and pointers
   ---------------------------
*/
nX    = bpg->nX    ;
graph = bpg->graph ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n FIND AUGMENTING PATH FOR NODE %d", xexp) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------
   mark the exposed node with the tag, 
   label with its exposure,
   and put into the list for traversal
   -----------------------------------
*/
mark[xexp]   = tag ;
list[0]      = xexp ;
labels[xexp] = exp[xexp] ;
par[xexp]    = -1 ;
now = last   = 0 ;
while ( now <= last ) {
/*
   -------------------------------
   check out next node on the list
   -------------------------------
*/
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n list[%d] = %d", now, list[now]) ;
      fflush(msgFile) ;
   }
   if ( list[now] < nX ) {
      x = list[now++] ;
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, ", X node %d", x) ;
         fflush(msgFile) ;
      }
      for ( cell = heads[x] ; cell != NULL ; cell = cell->xnext ) {
         y = cell->y ;
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, "\n   adjacent Y node %d, mark %d, r %d",
                    y, mark[y], cell->rxy) ;
            fflush(msgFile) ;
         }
         if ( mark[y] != tag && cell->rxy > 0 ) {
/*
            ----------------------------------------------------------
            y has not yet been visited and edge (x,y) is not saturated
            ----------------------------------------------------------
*/
            mark[y] = tag ;
            par[y]  =  x  ;
            if ( labels[x] < cell->rxy ) {
               labels[y] = labels[x] ;
            } else {
               labels[y] = cell->rxy ;
            }
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, ", label = %d", labels[y]) ;
               fflush(msgFile) ;
            }
            if ( exp[y] > 0 ) {
/*
               ----------------------------------
               y is exposed, set label and return
               ----------------------------------
*/
               if ( labels[y] > exp[y] ) {
                  labels[y] = exp[y] ;
               }
               if ( msglvl > 2 && msgFile != NULL ) {
                  fprintf(msgFile, 
                          "\n      Y node %d is exposed, label = %d",
                          y,labels[y]) ;
                  fflush(msgFile) ;
               }
               return(y) ;
            }
            list[++last] = y ;
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, ", adding to list") ;
               fflush(msgFile) ;
            }
         }
      }
   } else {
      y = list[now++] ;
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, ", Y node %d", y) ;
         fflush(msgFile) ;
      }
      for ( cell = heads[y] ; cell != NULL ; cell = cell->ynext ) {
         x = cell->x ;
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, "\n   adjacent X node %d, mark %d, f %d",
                    x, mark[x], cell->fxy) ;
            fflush(msgFile) ;
         }
         if ( mark[x] != tag && cell->fxy > 0 ) {
/*
            --------------------------------------------------------
            x has not yet been visited and there is flow from x to y
            --------------------------------------------------------
*/
            mark[x] = tag ;
            par[x]  =  y  ;
            if ( labels[y] < cell->fxy ) {
               labels[x] = labels[y] ;
            } else {
               labels[x] = cell->fxy ;
            }
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, ", labels = %d", labels[x]) ;
               fflush(msgFile) ;
            }
            list[++last] = x ;
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, ", adding to list") ;
               fflush(msgFile) ;
            }
         }
      }
   }
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   given an augmenting path starting with xexp and ending with yexp,
   augment the path by the amount in labels[yexp]

   96mar08, cca
   -----------------------------------------------------------------
*/
static void
augmentPath (
   BPG    *bpg,
   int    xexp,
   int    yexp,
   int    nvexp[],
   Cell   *heads[],
   int    labels[],
   int    par[],
   int    msglvl,
   FILE   *msgFile
) {
Cell   *cell ;
int    delta, x, y ;

if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n AUGMENT PATH FOR NODE %d", xexp) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   start at the exposed y node
   ---------------------------
*/
delta = labels[yexp] ;
y     = yexp ;
while ( 1 ) {
   x = par[y] ;
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n (x,y) = (%3d, %3d)", x, y) ;
      fflush(msgFile) ;
   }
   for ( cell = heads[y] ; cell != NULL ; cell = cell->ynext ) {
      if ( cell->x == x ) {
         cell->rxy -= delta ;
         cell->fxy += delta ;
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, ", r = %d, f = %d", cell->rxy, cell->fxy) ;
            fflush(msgFile) ;
         }
         break ;
      }
   }
   if ( cell == NULL ) {
      fprintf(stderr, "\n 1. error, x = %d, y = %d", x, y) ;
      exit(-1) ;
   }
   if ( x == xexp ) {
      break ;
   }
   y = par[x] ;
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n (x,y) = (%3d, %3d)", x, y) ;
      fflush(msgFile) ;
   }
   for ( cell = heads[x] ; cell != NULL ; cell = cell->xnext ) {
      if ( cell->y == y ) {
         cell->rxy += delta ;
         cell->fxy -= delta ;
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, ", r = %d, f = %d", cell->rxy, cell->fxy) ;
            fflush(msgFile) ;
         }
         break ;
      }
   }
   if ( cell == NULL ) {
      fprintf(stderr, "\n 2. error, x = %d, y = %d", x, y) ;
      exit(-1) ;
   }
}
nvexp[xexp] -= delta ;
nvexp[yexp] -= delta ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   set the flags for the DM decomposition

   created -- 96mar08, cca
   --------------------------------------
*/
static void
setDMflags (
   BPG    *bpg,
   Cell   *heads[],
   int    nvexp[],
   int    list[],
   int    mark[],
   int    dmflags[],
   int    stats[],
   int    msglvl,
   FILE   *msgFile
) {
Cell   *cell ;
int    ierr, last, now, nX, nY, x, xwght, y, ywght ;
int    *vwghts ;

if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n SET DM FLAGS") ;
   fflush(msgFile) ;
}

nX = bpg->nX ;
nY = bpg->nY ;
IVfill(nX + nY, dmflags, 0) ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n nvexp : ") ;
   IVfp80(msgFile, nX + nY, nvexp, 20, &ierr) ;
   fflush(msgFile) ;
}
/*
   -----------------------
   load exposed nodes in X
   -----------------------
*/
IVzero(nX + nY, mark) ;
last = -1 ;
for ( x = 0 ; x < nX ; x++ ) {
   if ( nvexp[x] > 0 ) {
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, "\n X node %d is exposed", x) ;
         fflush(msgFile) ;
      }
      mark[x]      = 1 ;
      list[++last] = x ;
      dmflags[x]   = 1 ;
   }
}
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n %d exposed X nodes :", last + 1) ;
   IVfp80(msgFile, 1 + last, list, 20, &ierr) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------------------------
   drop an alternating level structure from the exposed nodes in X
   ---------------------------------------------------------------
*/
now = 0 ;
while ( now <= last ) {
   if ( list[now] < nX ) {
      x = list[now++] ;
      dmflags[x] = 1 ;
      for ( cell = heads[x] ; cell != NULL ; cell = cell->xnext ) {
         y = cell->y ;
/*
         if ( mark[y] != 1 && cell->rxy > 0 ) {
*/
         if ( mark[y] != 1 ) {
            mark[y] = 1 ;
            list[++last] = y ;
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, "\n adding Y node %d", y) ;
               fflush(msgFile) ;
            }
         }
      }
   } else {
      y = list[now++] ;
      dmflags[y] = 2 ;
      for ( cell = heads[y] ; cell != NULL ; cell = cell->ynext ) {
         x = cell->x ;
         if ( mark[x] != 1 && cell->fxy > 0 ) {
            mark[x] = 1 ;
            list[++last] = x ;
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, "\n adding X node %d", x) ;
               fflush(msgFile) ;
            }
         }
      }
   }
}
/*
   -----------------------
   load exposed nodes in Y
   -----------------------
*/
IVzero(nX + nY, mark) ;
last = -1 ;
for ( y = nX ; y < nX + nY ; y++ ) {
   if ( nvexp[y] > 0 ) {
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, "\n Y node %d is exposed", y) ;
         fflush(msgFile) ;
      }
      mark[y]      = 1 ;
      list[++last] = y ;
      dmflags[y]   = 1 ;
   }
}
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n %d exposed Y nodes :", last + 1) ;
   IVfp80(msgFile, 1 + last, list, 20, &ierr) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------------------------
   drop an alternating level structure from the exposed nodes in Y
   ---------------------------------------------------------------
*/
now = 0 ;
while ( now <= last ) {
   if ( list[now] < nX ) {
      x = list[now++] ;
      dmflags[x] = 2 ;
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, "\n checking out X node %d", x) ;
         fflush(msgFile) ;
      }
      for ( cell = heads[x] ; cell != NULL ; cell = cell->xnext ) {
         y = cell->y ;
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, "\n    adjacent Y node %d, mark %d, f %d",
                    y, mark[y], cell->fxy) ;
            fflush(msgFile) ;
         }
         if ( mark[y] != 1 && cell->fxy > 0 ) {
            mark[y] = 1 ;
            list[++last] = y ;
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, ", adding ") ;
               fflush(msgFile) ;
            }
         }
      }
   } else {
      y = list[now++] ;
      dmflags[y] = 1 ;
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, "\n checking out Y node %d", y) ;
         fflush(msgFile) ;
      }
      for ( cell = heads[y] ; cell != NULL ; cell = cell->ynext ) {
         x = cell->x ;
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(msgFile, "\n    adjacent X node %d, mark %d, r %d",
                    x, mark[x], cell->rxy) ;
            fflush(msgFile) ;
         }
/*
         if ( mark[x] != 1 && cell->rxy > 0 ) {
*/
         if ( mark[x] != 1 ) {
            mark[x] = 1 ;
            list[++last] = x ;
            if ( msglvl > 2 && msgFile != NULL ) {
               fprintf(msgFile, ", adding") ;
               fflush(msgFile) ;
            }
         }
      }
   }
}
/*
   --------------------------
   fill the statistics vector
   --------------------------
*/
IVzero(6, stats) ;
vwghts = bpg->graph->vwghts ;
for ( x = 0 ; x < nX ; x++ ) {
   xwght = (vwghts != NULL) ? vwghts[x] : 1 ;
   switch ( dmflags[x] ) {
   case 0 : stats[2] += xwght ; break ;
   case 1 : stats[0] += xwght ; break ;
   case 2 : stats[1] += xwght ; break ;
   default :
      fprintf(stderr, "\n fatal error in BPG_DMviaMaxFlow"
              "\n dmflags[%d] = %d\n", x, dmflags[x]) ;
      exit(-1) ;
   }
}
for ( y = nX ; y < nX + nY ; y++ ) {
   ywght = (vwghts != NULL) ? vwghts[y] : 1 ;
   switch ( dmflags[y] ) {
   case 0 : stats[5] += ywght ; break ;
   case 1 : stats[3] += ywght ; break ;
   case 2 : stats[4] += ywght ; break ;
   default :
      fprintf(stderr, "\n fatal error in BPG_DMviaMaxFlow"
              "\n dmflags[%d] = %d\n", y, dmflags[y]) ;
      exit(-1) ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
