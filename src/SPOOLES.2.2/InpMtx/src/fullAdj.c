/*  fullAdj.c  */

#include "../InpMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to return the full, symmetric adjacency IVL object
              for the graph of A + A^T

   created -- 98jan28, cca
   -------------------------------------------------------------
*/
IVL *
InpMtx_fullAdjacency (
   InpMtx   *inpmtx
) {
int       count, ii, jvtx, jsize, kvtx, nvtx ;
int       *jind, *list, *mark ;
IP        *baseIP, *freeIP, *ip ;
IP        **head ;
IVL       *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_fullAdjacency(%p)"
           "\n NULL input\n", inpmtx) ;
   exit(-1) ;
}
/*
    ----------------------
    check for empty matrix
    ----------------------
*/
if ( inpmtx->nent == 0 ) {
   return(NULL) ;
}
/*
   -------------------------------------------------
   check for invalid coordinate type or storage mode
   -------------------------------------------------
*/
if ( ! (INPMTX_IS_BY_ROWS(inpmtx) || INPMTX_IS_BY_COLUMNS(inpmtx)) ) {
   InpMtx_changeCoordType(inpmtx, INPMTX_BY_ROWS) ;
}
if ( ! INPMTX_IS_BY_VECTORS(inpmtx) ) {
   InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
}
nvtx = 1 + IV_max(&inpmtx->ivec1IV) ;
if ( nvtx < (count = 1 + IV_max(&inpmtx->ivec2IV)) ) {
   nvtx = count ;
}
/*
   -------------------------
   initialize the IVL object
   -------------------------
*/
adjIVL = IVL_new() ;
IVL_init1(adjIVL, IVL_CHUNKED, nvtx) ;
list = IVinit(nvtx, -1) ;
mark = IVinit(nvtx, -1) ;
ALLOCATE(head, struct _IP *, nvtx) ;
baseIP = freeIP = NULL ;
/*
   ---------------------------
   initially link the vertices
   ---------------------------
*/
for ( jvtx = 0 ; jvtx < nvtx ; jvtx++ ) {
   head[jvtx] = NULL ;
}
for ( jvtx = 0 ; jvtx < nvtx ; jvtx++ ) {
   InpMtx_vector(inpmtx, jvtx, &jsize, &jind) ;
   if ( jsize > 0 ) {
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         if ( (kvtx = jind[ii]) < jvtx ) {
/*
            -----------------
            link jvtx to kvtx
            -----------------
*/
            if ( (ip = freeIP) == NULL ) {
               ip = IP_init(nvtx+1, IP_FORWARD) ;
               ip->next = baseIP ; baseIP = ip ;
               freeIP = baseIP + 1 ;
               ip = freeIP ;
            }
            freeIP = ip->next ;
            ip->val    = jvtx ;
            ip->next   = head[kvtx] ;
            head[kvtx] = ip ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    linking %d to %d, %p to %p", 
                    jvtx, kvtx, ip, ip->next) ;
            fflush(stdout) ;
#endif
         }
      }
   }
}   
for ( jvtx = 0 ; jvtx < nvtx ; jvtx++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n working on vertex %d :", jvtx) ;
   fflush(stdout) ;
#endif
   list[0] = jvtx ;
#if MYDEBUG > 0
   fprintf(stdout, "\n    0. adding %d", jvtx) ;
   fflush(stdout) ;
#endif
   mark[jvtx] = jvtx ;
   count = 1 ;
/*
   -------------------------------------
   check previous rows that contain jvtx
   -------------------------------------
*/
   while ( (ip = head[jvtx]) != NULL ) {
      kvtx = ip->val ;
#if MYDEBUG > 0
         fprintf(stdout, "\n    1. kvtx %d", kvtx) ;
         fflush(stdout) ;
#endif
      if ( mark[kvtx] != jvtx ) {
         mark[kvtx] = jvtx ;
         list[count++] = kvtx ;
#if MYDEBUG > 0
         fprintf(stdout, ", adding") ;
         fflush(stdout) ;
#endif
      }
      head[jvtx] = ip->next ;
      ip->next = freeIP ;
      freeIP = ip ;
   }
/*
   ----------------------------
   get the indices for row jvtx
   ----------------------------
*/
   InpMtx_vector(inpmtx, jvtx, &jsize, &jind) ;
   if ( jsize > 0 ) {
/*
      ----------------------------
      add row indices for row jvtx
      ----------------------------
*/
#if MYDEBUG > 0
      fprintf(stdout, "\n    InpMtx row %d :", jvtx) ;
      IVfprintf(stdout, jsize, jind) ;
      fflush(stdout) ;
#endif
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         kvtx = jind[ii] ;
         if ( mark[kvtx] != jvtx ) {
            mark[kvtx] = jvtx ;
            list[count++] = kvtx ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    2. adding %d", kvtx) ;
            fflush(stdout) ;
#endif
         }
         if ( kvtx > jvtx ) {
/*
            -----------------
            link jvtx to kvtx
            -----------------
*/
            if ( (ip = freeIP) == NULL ) {
               ip = IP_init(nvtx+1, IP_FORWARD) ;
               ip->next = baseIP ; baseIP = ip ;
               freeIP = baseIP + 1 ;
               ip = freeIP ;
            }
            freeIP     = ip->next ;
            ip->val    = jvtx ;
            ip->next   = head[kvtx] ;
            head[kvtx] = ip ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    linking %d to %d, %p to %p", 
                    jvtx, kvtx, ip, ip->next) ;
            fflush(stdout) ;
#endif
         }
      }
   }
/*
   ------------------------------------------------
   list is complete, sort and insert into adjacency
   ------------------------------------------------
*/
   IVqsortUp(count, list) ;
   IVL_setList(adjIVL, jvtx, count, list) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list) ;
IVfree(mark) ;
FREE(head) ;
while ( (ip = baseIP) != NULL ) {
   baseIP = baseIP->next ;
   IP_free(ip) ;
}
return(adjIVL) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to return the full, symmetric adjacency IVL object
              for the graph of (A + B) + (A + B)^T

   created -- 97nov05, cca
   -------------------------------------------------------------
*/
IVL *
InpMtx_fullAdjacency2 (
   InpMtx   *inpmtxA,
   InpMtx   *inpmtxB
) {
int       count, ierr, ii, jvtx, jsize, kvtx, nvtx ;
int       *jind, *list, *mark ;
IP        *baseIP, *freeIP, *ip ;
IP        **head ;
IVL       *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtxA == NULL && inpmtxB == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_fullAdjacency2(%p,%p)"
           "\n both input matrices are NULL\n", inpmtxA, inpmtxB) ;
   exit(-1) ;
}
/*
   ------------------------------
   cases where one matrix is NULL
   ------------------------------
*/
if ( inpmtxA == NULL ) {
   adjIVL = InpMtx_fullAdjacency(inpmtxB) ;
   return(adjIVL) ;
} else if ( inpmtxB == NULL ) {
   adjIVL = InpMtx_fullAdjacency(inpmtxA) ;
   return(adjIVL) ;
}
/*
   -------------------------------------------------
   check for invalid coordinate type or storage mode
   -------------------------------------------------
*/
if ( ! (INPMTX_IS_BY_ROWS(inpmtxA) || INPMTX_IS_BY_COLUMNS(inpmtxA)) ) {
   InpMtx_changeCoordType(inpmtxA, INPMTX_BY_ROWS) ;
}
if ( ! INPMTX_IS_BY_VECTORS(inpmtxA) ) {
   InpMtx_changeStorageMode(inpmtxA, INPMTX_BY_VECTORS) ;
}
if ( ! (INPMTX_IS_BY_ROWS(inpmtxB) || INPMTX_IS_BY_COLUMNS(inpmtxB)) ) {
   InpMtx_changeCoordType(inpmtxB, INPMTX_BY_ROWS) ;
}
if ( ! INPMTX_IS_BY_VECTORS(inpmtxB) ) {
   InpMtx_changeStorageMode(inpmtxB, INPMTX_BY_VECTORS) ;
}
nvtx = 1 + IV_max(&inpmtxA->ivec1IV) ;
if ( nvtx < (count = 1 + IV_max(&inpmtxA->ivec2IV)) ) {
   nvtx = count ;
}
if ( nvtx < (count = 1 + IV_max(&inpmtxB->ivec1IV)) ) {
   nvtx = count ;
}
if ( nvtx < (count = 1 + IV_max(&inpmtxB->ivec2IV)) ) {
   nvtx = count ;
}
/*
   -------------------------
   initialize the IVL object
   -------------------------
*/
adjIVL = IVL_new() ;
IVL_init1(adjIVL, IVL_CHUNKED, nvtx) ;
/*
   ------------------------------
   initialize the working storage
   ------------------------------
*/
list = IVinit(nvtx, -1) ;
mark = IVinit(nvtx, -1) ;
ALLOCATE(head, struct _IP *, nvtx) ;
baseIP = freeIP = NULL ;
/*
   ---------------------------
   initially link the vertices
   ---------------------------
*/
for ( jvtx = 0 ; jvtx < nvtx ; jvtx++ ) {
   head[jvtx] = NULL ;
}
for ( jvtx = 0 ; jvtx < nvtx ; jvtx++ ) {
   InpMtx_vector(inpmtxA, jvtx, &jsize, &jind) ;
   if ( jsize > 0 ) {
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         if ( (kvtx = jind[ii]) < jvtx ) {
/*
            -----------------
            link jvtx to kvtx
            -----------------
*/
            if ( (ip = freeIP) == NULL ) {
               ip = IP_init(nvtx+1, IP_FORWARD) ;
               ip->next = baseIP ; baseIP = ip ;
               freeIP = baseIP + 1 ;
               ip = freeIP ;
            }
            freeIP     = ip->next ;
            ip->val    = jvtx ;
            ip->next   = head[kvtx] ;
            head[kvtx] = ip ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    linking %d to %d, %p to %p", 
                    jvtx, kvtx, ip, ip->next) ;
            fflush(stdout) ;
#endif
         }
      }
   }
   InpMtx_vector(inpmtxB, jvtx, &jsize, &jind) ;
   if ( jsize > 0 ) {
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         if ( (kvtx = jind[ii]) < jvtx ) {
/*
            -----------------
            link jvtx to kvtx
            -----------------
*/
            if ( (ip = freeIP) == NULL ) {
               ip = IP_init(nvtx+1, IP_FORWARD) ;
               ip->next = baseIP ; baseIP = ip ;
               freeIP = baseIP + 1 ;
               ip = freeIP ;
            }
            freeIP     = ip->next ;
            ip->val    = jvtx ;
            ip->next   = head[kvtx] ;
            head[kvtx] = ip ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    linking %d to %d, %p to %p", 
                    jvtx, kvtx, ip, ip->next) ;
            fflush(stdout) ;
#endif
         }
      }
   }
}   
for ( jvtx = 0 ; jvtx < nvtx ; jvtx++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n vertex %d :", jvtx) ;
   fflush(stdout) ;
#endif
   list[0] = jvtx ;
#if MYDEBUG > 0
   fprintf(stdout, "\n    0. adding %d", jvtx) ;
   fflush(stdout) ;
#endif
   mark[jvtx] = jvtx ;
   count = 1 ;
/*
   -------------------------------------
   check previous rows that contain jvtx
   -------------------------------------
*/
   while ( (ip = head[jvtx]) != NULL ) {
      kvtx = ip->val ;
#if MYDEBUG > 0
         fprintf(stdout, "\n    1. kvtx %d", kvtx) ;
         fflush(stdout) ;
#endif
      if ( mark[kvtx] != jvtx ) {
         mark[kvtx] = jvtx ;
         list[count++] = kvtx ;
#if MYDEBUG > 0
         fprintf(stdout, ", adding") ;
         fflush(stdout) ;
#endif
      }
      head[jvtx] = ip->next ;
      ip->next = freeIP ;
      freeIP = ip ;
   }
/*
   ----------------------------
   get the indices for row jvtx
   ----------------------------
*/
   InpMtx_vector(inpmtxA, jvtx, &jsize, &jind) ;
   if ( jsize > 0 ) {
/*
      ----------------------------
      add row indices for row jvtx
      ----------------------------
*/
#if MYDEBUG > 0
      fprintf(stdout, "\n    InpMtx row %d :", jvtx) ;
      IVfprintf(stdout, jsize, jind) ;
      fflush(stdout) ;
#endif
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         kvtx = jind[ii] ;
         if ( mark[kvtx] != jvtx ) {
            mark[kvtx] = jvtx ;
            list[count++] = kvtx ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    2. adding %d", kvtx) ;
            fflush(stdout) ;
#endif
         }
         if ( kvtx > jvtx ) {
/*
            -----------------
            link jvtx to kvtx
            -----------------
*/
            if ( (ip = freeIP) == NULL ) {
               ip = IP_init(nvtx+1, IP_FORWARD) ;
               ip->next = baseIP ; baseIP = ip ;
               freeIP = baseIP + 1 ;
               ip = freeIP ;
            }
            freeIP     = ip->next ;
            ip->val    = jvtx ;
            ip->next   = head[kvtx] ;
            head[kvtx] = ip ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    linking %d to %d, %p to %p", 
                    jvtx, kvtx, ip, ip->next) ;
            fflush(stdout) ;
#endif
         }
      }
   }
   InpMtx_vector(inpmtxB, jvtx, &jsize, &jind) ;
   if ( jsize > 0 ) {
/*
      ----------------------------
      add row indices for row jvtx
      ----------------------------
*/
#if MYDEBUG > 0
      fprintf(stdout, "\n    InpMtx row %d :", jvtx) ;
      IVfprintf(stdout, jsize, jind) ;
      fflush(stdout) ;
#endif
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         kvtx = jind[ii] ;
         if ( mark[kvtx] != jvtx ) {
            mark[kvtx] = jvtx ;
            list[count++] = kvtx ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    2. adding %d", kvtx) ;
            fflush(stdout) ;
#endif
         }
         if ( kvtx > jvtx ) {
/*
            -----------------
            link jvtx to kvtx
            -----------------
*/
            if ( (ip = freeIP) == NULL ) {
               ip = IP_init(nvtx+1, IP_FORWARD) ;
               ip->next = baseIP ; baseIP = ip ;
               freeIP = baseIP + 1 ;
               ip = freeIP ;
            }
            freeIP     = ip->next ;
            ip->val    = jvtx ;
            ip->next   = head[kvtx] ;
            head[kvtx] = ip ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    linking %d to %d, %p to %p", 
                    jvtx, kvtx, ip, ip->next) ;
            fflush(stdout) ;
#endif
         }
      }
   }
/*
   ------------------------------------------------
   list is complete, sort and insert into adjacency
   ------------------------------------------------
*/
   IVqsortUp(count, list) ;
   IVL_setList(adjIVL, jvtx, count, list) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list) ;
IVfree(mark) ;
FREE(head) ;
while ( (ip = baseIP) != NULL ) {
   baseIP = baseIP->next ;
   IP_free(ip) ;
}

return(adjIVL) ; }

/*--------------------------------------------------------------------*/
