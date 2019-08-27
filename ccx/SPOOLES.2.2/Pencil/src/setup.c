/*  setup.c  */

#include "../Pencil.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   initialize the matrix pencil A + sigma*B

   myid -- id of process, used in MPI implementation
      if myid = 0 then 
         the pencil is loaded with the matrices read from the files
      else
         the pencil is loaded with the empty matrices 
      endif
   symflag -- symmetry flag, 
      SPOOLES_SYMMETRIC    -- symmetric
      SPOOLES_HERMITIAN    -- hermitian
      SPOOLES_NONSYMMETRIC -- nonsymmetric
      if symmetric or hermitian, drop entries in lower triangle 
   inpmtxAfile  -- filename for A
   sigma        -- scaling factor
   inpmtxBfile  -- filename for B
   randomflag   -- random flag, 
      if 1 then fill with random numbers
   drand        -- random number generator
   msglvl       -- message level
   msgFile      -- message file

   return value -- pointer to a Pencil object

   created -- 98may02, cca
   ----------------------------------------------------------------
*/
Pencil *
Pencil_setup (
   int        myid,
   int        symflag,
   char       *inpmtxAfile,
   double     sigma[],
   char       *inpmtxBfile,
   int        randomflag,
   Drand      *drand,
   int        msglvl,
   FILE       *msgFile
) {
InpMtx   *inpmtxA, *inpmtxB ;
double    t1, t2 ;
Pencil   *pencil ;
int       rc ;

switch ( symflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n fatal error in Pencil_setup()"
           "\n bad symflag %d\n", symflag) ;
   exit(-1) ;
   break ;
}

if ( strcmp(inpmtxAfile, "none") != 0 ) {
/*
   --------------------------
   read in the InpMtx object
   --------------------------
*/
   inpmtxA = InpMtx_new() ;
   if ( myid == 0 ) {
      MARKTIME(t1) ;
      rc = InpMtx_readFromFile(inpmtxA, inpmtxAfile) ;
      MARKTIME(t2) ;
      fprintf(msgFile,"\n CPU %8.3f : read in inpmtxA", t2 - t1) ;
      if ( rc != 1 ) {
         fprintf(msgFile, 
                 "\n return value %d from InpMtx_readFromFile(%p,%s)",
                 rc, inpmtxA, inpmtxAfile) ;
         exit(-1) ;
      }
      if ( symflag == SPOOLES_SYMMETRIC 
        || symflag == SPOOLES_HERMITIAN ){
/*
         ----------------------------------------------------
         symmetric matrix, drop entries in the lower triangle
         ----------------------------------------------------
*/
         MARKTIME(t1) ;
         InpMtx_dropLowerTriangle(inpmtxA) ;
         InpMtx_changeStorageMode(inpmtxA, 1) ;
/*
         ------------------------------------
         sort and compress the matrix entries
         ------------------------------------
*/
         InpMtx_sortAndCompress(inpmtxA) ;
         MARKTIME(t2) ;
         fprintf(msgFile,"\n CPU %8.3f : initialize inpmtxA", t2 - t1) ;
      }
      if ( randomflag == 1 ) {
/*
         -----------------------------------
         fill the matrix with random numbers
         -----------------------------------
*/
         MARKTIME(t1) ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) {
            Drand_fillDvector(drand, inpmtxA->nent, 
                              DV_entries(&inpmtxA->dvecDV)) ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) {
            Drand_fillDvector(drand, 2*inpmtxA->nent, 
                              DV_entries(&inpmtxA->dvecDV)) ;
         }
         MARKTIME(t2) ;
         fprintf(msgFile, 
          "\n CPU %8.3f : fill inpmtxA with random numbers", t2 - t1) ;
      } else if ( randomflag == -1 ) {
/*
         -----------------------------
         double one entry in magnitude
         -----------------------------
*/
         double   *dvec = DV_entries(&inpmtxA->dvecDV) ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) {
            dvec[0] *= 2 ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) {
            dvec[0] *= 2 ; dvec[1] *= 2 ;
         }
      }
   }
/*
   --------------------------------------
   change the coordinate type to chevrons
   and the storage mode to sorted triples
   --------------------------------------
*/
   MARKTIME(t1) ;
   InpMtx_changeCoordType(inpmtxA, 3) ;
   InpMtx_changeStorageMode(inpmtxA, 2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %8.3f : change inpmtxA to chevrons", t2 - t1) ;
} else {
   inpmtxA = NULL ;
}
if ( strcmp(inpmtxBfile, "none") != 0 ) {
/*
   --------------------------
   read in the InpMtx object
   --------------------------
*/
   inpmtxB = InpMtx_new() ;
   if ( myid == 0 ) {
      MARKTIME(t1) ;
      rc = InpMtx_readFromFile(inpmtxB, inpmtxBfile) ;
      MARKTIME(t2) ;
      fprintf(msgFile,"\n CPU %8.3f : read in inpmtxB", t2 - t1) ;
      if ( rc != 1 ) {
         fprintf(msgFile, 
                 "\n return value %d from InpMtx_readFromFile(%p,%s)",
                 rc, inpmtxB, inpmtxBfile) ;
         exit(-1) ;
      }
      if ( symflag == SPOOLES_SYMMETRIC 
        || symflag == SPOOLES_HERMITIAN ){
/*
         ----------------------------------------------------
         symmetric matrix, drop entries in the lower triangle
         ----------------------------------------------------
*/
         MARKTIME(t1) ;
         InpMtx_dropLowerTriangle(inpmtxB) ;
         InpMtx_changeStorageMode(inpmtxB, 1) ;
/*
         ------------------------------------
         sort and compress the matrix entries
         ------------------------------------
*/
         InpMtx_sortAndCompress(inpmtxB) ;
         MARKTIME(t2) ;
         fprintf(msgFile,"\n CPU %8.3f : initialize inpmtxB", t2 - t1) ;
      }
      if ( randomflag == 1 ) {
/*
         -----------------------------------
         fill the matrix with random numbers
         -----------------------------------
*/
         MARKTIME(t1) ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) {
            Drand_fillDvector(drand, inpmtxB->nent, 
                              DV_entries(&inpmtxB->dvecDV)) ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) {
            Drand_fillDvector(drand, 2*inpmtxB->nent, 
                              DV_entries(&inpmtxB->dvecDV)) ;
         }
         MARKTIME(t2) ;
         fprintf(msgFile, 
          "\n CPU %8.3f : fill inpmtxB with random numbers", t2 - t1) ;
      } else if ( randomflag == -1 ) {
/*
         -----------------------------
         double one entry in magnitude
         -----------------------------
*/
         double   *dvec = DV_entries(&inpmtxB->dvecDV) ;
         if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) {
            dvec[0] *= 2 ;
         } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) {
            dvec[0] *= 2 ; dvec[1] *= 2 ;
         }
      }
   }
/*
   --------------------------------------
   change the coordinate type to chevrons
   and the storage mode to sorted triples
   --------------------------------------
*/
   MARKTIME(t1) ;
   InpMtx_changeCoordType(inpmtxB, 3) ;
   InpMtx_changeStorageMode(inpmtxB, 2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %8.3f : change inpmtxB to chevrons", t2 - t1) ;
} else {
   inpmtxB = NULL ;
}
/*
   -----------------------------
   initialize the Pencil object
   -----------------------------
*/
if ( inpmtxA != NULL && inpmtxB != NULL ) {
   if ( inpmtxA->inputMode != inpmtxB->inputMode ) {
      fprintf(stderr, "\n fatal error in Pencil_setup()"
              "\n inpmtxA->inputMode = %d, inpmtxB->inputMode = %d\n",
              inpmtxA->inputMode, inpmtxB->inputMode) ;
      exit(-1) ;
   }
   if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) {
      pencil = Pencil_new() ;
      Pencil_init(pencil, SPOOLES_REAL, symflag, 
                  inpmtxA, sigma, inpmtxB) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) {
      pencil = Pencil_new() ;
      Pencil_init(pencil, SPOOLES_COMPLEX, symflag, 
                  inpmtxA, sigma, inpmtxB) ;
   } else {
      fprintf(stderr, "\n fatal error in Pencil_setup()"
              "\n inpmtxA->inputMode = %d\n", inpmtxA->inputMode) ;
      exit(-1) ;
   }
} else if ( inpmtxA != NULL ) {
   if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) {
      pencil = Pencil_new() ;
      Pencil_init(pencil, SPOOLES_REAL, symflag, 
                  inpmtxA, sigma, inpmtxB) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) {
      pencil = Pencil_new() ;
      Pencil_init(pencil, SPOOLES_COMPLEX, symflag, 
                  inpmtxA, sigma, inpmtxB) ;
   } else {
      fprintf(stderr, "\n fatal error in Pencil_setup()"
              "\n inpmtxA->inputMode = %d\n", 
                  inpmtxA->inputMode) ;
      exit(-1) ;
   }
} else if ( inpmtxB != NULL ) {
   if ( INPMTX_IS_REAL_ENTRIES(inpmtxB) ) {
      pencil = Pencil_new() ;
      Pencil_init(pencil, SPOOLES_REAL, symflag, 
                  inpmtxA, sigma, inpmtxB) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxB) ) {
      pencil = Pencil_new() ;
      Pencil_init(pencil, SPOOLES_COMPLEX, symflag, 
                  inpmtxA, sigma, inpmtxB) ;
   } else {
      fprintf(stderr, "\n fatal error in Pencil_setup()"
              "\n inpmtxB->inputMode = %d\n", inpmtxB->inputMode) ;
      exit(-1) ;
   }
}
return(pencil) ; }

/*--------------------------------------------------------------------*/
