/*  init.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   return the number of bytes needed to store the chevron

   created -- 98apr30, cca
   ------------------------------------------------------
*/
int
Chv_nbytesNeeded (
   int   nD,
   int   nL,
   int   nU,
   int   type,
   int   symflag
) {
int   nbytes, nent, nint ;
/*
   --------------
   check the data
   --------------
*/
if ( nD < 0 || nL < 0 || nU < 0 ) {
   fprintf(stderr, "\n fatal error in Chv_nbytesNeeded()"
           "\n bad input, nD = %d, nL = %d, nU = %d\n", nD, nL, nU) ;
   exit(-1) ;
}
nbytes = 0 ;
switch ( type ) {
case SPOOLES_REAL :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
      nint = 6 + nD + nU ;
      nent = (nD*(nD+1))/2 + nD*nU ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      nint = 6 + 2*nD + nL + nU ;
      nent = nD*(nD + nL + nU) ;
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in Chv_nbytesNeeded()"
              "\n type = SPOOLES_REAL, invalid symflag = %d"
              "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n",
              symflag) ;
      exit(-1) ;
   }
   if ( 2*sizeof(int) == sizeof(double) ) {
      nbytes = ((nint + 1)/2 + nent)*sizeof(double) ;
   } else if ( sizeof(int) == sizeof(double) ) {
      nbytes = (nint + nent)*sizeof(double) ;
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_nbytesNeeded()"
              "\n sizeof(int) = %d, sizeof(double) = %d",
              sizeof(int), sizeof(double)) ;
      exit(-1) ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN:
      nint = 6 + nD + nU ;
      nent = (nD*(nD+1))/2 + nD*nU ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      nint = 6 + 2*nD + nL + nU ;
      nent = nD*(nD + nL + nU) ;
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in Chv_nbytesNeeded()"
              "\n type = SPOOLES_COMPLEX, invalid symflag = %d"
              "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
              "\n or SPOOLES_NONSYMMETRIC\n",
              symflag) ;
      exit(-1) ;
   }
   if ( 2*sizeof(int) == sizeof(double) ) {
      nbytes = ((nint + 1)/2 + 2*nent)*sizeof(double) ;
   } else if ( sizeof(int) == sizeof(double) ) {
      nbytes = (nint + 2*nent)*sizeof(double) ;
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_nbytesNeeded()"
              "\n sizeof(int) = %d, sizeof(double) = %d",
              sizeof(int), sizeof(double)) ;
      exit(-1) ;
   }
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in Chv_nbytesNeeded()"
           "\n invalid type = %d"
           "\n must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           type) ;
   break ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   return the number of bytes in the workspace owned by this object

   created -- 98apr30, cca
   ----------------------------------------------------------------
*/
int
Chv_nbytesInWorkspace (
   Chv   *chv
) {
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_nbytesInWorkspace(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
return(sizeof(double)*DV_maxsize(&chv->wrkDV)) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   set the number of bytes in the workspace owned by this object

   created -- 98apr30, cca
   ----------------------------------------------------------------
*/
void
Chv_setNbytesInWorkspace (
   Chv   *chv,
   int    nbytes
) {
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_setNbytesInWorkspace(%p,%d)"
           "\n bad input\n", chv, nbytes) ;
   exit(-1) ;
}
DV_setSize(&chv->wrkDV, nbytes/sizeof(double)) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   purpose -- set the fields

   created -- 98apr30, cca
   ----------------------------
*/
void
Chv_setFields (
   Chv     *chv,
   int      id,
   int      nD,
   int      nL,
   int      nU,
   int      type,
   int      symflag
) {
double   *dbuffer ;
int      nint     ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL || nD <= 0 || nL < 0 || nU < 0 ) {
   fprintf(stderr, "\n fatal error in Chv_setFields()"
           "\n bad input, chv %p, nD %d, nL %d, nU %d\n", 
           chv, nD, nL, nU) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
           "\n fatal error in Chv_setFields()"
           "\n type = SPOOLES_REAL, symflag = %d"
           "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n", 
           symflag) ;
      exit(-1) ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in Chv_setFields()"
              "\n type = SPOOLES_COMPLEX, symflag = %d"
              "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
              "\n or SPOOLES_NONSYMMETRIC\n",
              symflag) ;
      exit(-1) ;
   }
   break ;
default :
   fprintf(stderr, 
     "\n fatal error in Chv_setFields()"
     "\n type = %d"
     "\n must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
     type) ;
   exit(-1) ;
} 
dbuffer = DV_entries(&chv->wrkDV) ;
ibuffer = (int *) dbuffer ;
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
chv->id      = ibuffer[0] = id      ;
chv->nD      = ibuffer[1] = nD      ;
chv->nL      = ibuffer[2] = nL      ;
chv->nU      = ibuffer[3] = nU      ;
chv->type    = ibuffer[4] = type    ;
chv->symflag = ibuffer[5] = symflag ;
/*
   -------------------------------------------
   set the colind, rowind and entries pointers
   -------------------------------------------
*/
chv->colind = ibuffer + 6 ;
nint = 6 + nD + nU ;
if ( symflag == SPOOLES_NONSYMMETRIC ) {
   chv->rowind = chv->colind + nD + nU ;
   nint += nD + nL ;
} else {
   chv->rowind = NULL ;
}
if ( sizeof(int) == sizeof(double) ) {
   chv->entries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   chv->entries = dbuffer + (nint + 1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   purpose -- basic initializer

   created -- 98apr30, cca
   ----------------------------
*/
void
Chv_init (
   Chv     *chv,
   int      id,
   int      nD,
   int      nL,
   int      nU,
   int      type,
   int      symflag
) {
int      nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL || nD <= 0 || nL < 0 || nU < 0 ) {
   fprintf(stderr, 
           "\n fatal error in Chv_init()"
           "\n bad input, chv %p, nD %d, nL %d, nU %d\n", 
           chv, nD, nL, nU) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
           "\n fatal error in Chv_init()"
           "\n type = SPOOLES_REAL, symflag = %d"
           "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n", 
           symflag) ;
      exit(-1) ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in Chv_init()"
              "\n type = SPOOLES_COMPLEX, symflag = %d"
              "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
              "\n or SPOOLES_NONSYMMETRIC\n",
              symflag) ;
      exit(-1) ;
   }
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in Chv_init()"
           "\n type = %d"
           "\n must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           type) ;
   exit(-1) ;
} 
/*
   -------------------------------------------------------
   get and set the number of bytes needed in the workspace
   -------------------------------------------------------
*/
nbytes = Chv_nbytesNeeded(nD, nL, nU, type, symflag) ;
Chv_setNbytesInWorkspace(chv, nbytes) ;
/*
   --------------
   set the fields
   --------------
*/
Chv_setFields(chv, id, nD, nL, nU, type, symflag) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- initializer with pointers

   created -- 98apr30, cca
   ------------------------------------
*/
void
Chv_initWithPointers (
   Chv      *chv,
   int      id,
   int      nD,
   int      nL,
   int      nU,
   int      type,
   int      symflag,
   int      *rowind,
   int      *colind,
   double   *entries
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL || nD <= 0 || nL < 0 || nU < 0 ) {
   fprintf(stderr, 
           "\n fatal error in Chv_initWithPointers() "
           "\n chv = %p, nD = %d, nL = %d, nU = %d\n",
           chv, nD, nL, nU) ;
   exit(-1) ;
}
if (  entries == NULL || colind == NULL 
   || (symflag == SPOOLES_NONSYMMETRIC && rowind == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_init()"
           "\n entries = %p, colind = %p, rowind = %p, symflag = %d\n",
           entries, colind, rowind, symflag) ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in Chv_initFromPointers()"
              "\n type = SPOOLES_REAL, symflag = %d"
              "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n", 
              symflag) ;
      exit(-1) ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in Chv_initFromPointers()"
              "\n type = SPOOLES_COMPLEX, symflag = %d"
              "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
              "\n or SPOOLES_NONSYMMETRIC\n",
              symflag) ;
      exit(-1) ;
   }
   break ;
default :
   fprintf(stderr, 
            "\n fatal error in Chv_initFromPointers()"
            "\n type = %d"
            "\n must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
     type) ;
   exit(-1) ;
}
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
chv->id      = id      ;
chv->nD      = nD      ;
chv->nL      = nL      ;
chv->nU      = nU      ;
chv->type    = type    ;
chv->symflag = symflag ;
/*
   --------------------------
   set up the working storage
   --------------------------
*/
chv->entries = entries ;
chv->colind  = colind  ;
if ( symflag == SPOOLES_NONSYMMETRIC ) {
   chv->rowind = rowind ;
} else {
   chv->rowind = NULL ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to initialize the object from its working storage,
              used when the object is an MPI message

   created -- 98apr30
   -------------------------------------------------------------
*/
void
Chv_initFromBuffer (
   Chv   *chv
) {
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_initFromBuffer(%p) "
           "\n bad input\n", chv) ;
   exit(-1) ;
}
ibuffer = (int *) DV_entries(&chv->wrkDV) ;
Chv_setFields(chv, ibuffer[0], ibuffer[1], ibuffer[2], 
               ibuffer[3], ibuffer[4], ibuffer[5]) ;

return ; }

/*--------------------------------------------------------------------*/
