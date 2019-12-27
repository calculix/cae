/*  IO.c  */

#include "../ILUMtx.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*
   ---------------------------------------------------------------
   purpose -- to write an ILUMtx object to a file in matlab format

   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtx
     -4  -- bad symmetryflag for mtx
     -5  -- bad LstorageMode for mtx
     -6  -- bad UstorageMode for mtx
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
     -14 -- Lname is NULL
     -15 -- Dname is NULL
     -16 -- Uname is NULL
     -17 -- fp is NULL

   created -- 98oct03, cca
   ---------------------------------------------------------------
*/
int
ILUMtx_writeForMatlab ( 
   ILUMtx   *mtx,
   char     *Lname,
   char     *Dname,
   char     *Uname,
   FILE     *fp
) {
int      ieqn, ii, neqns, size ;
int      *ind, *sizesU, *sizesL ;
int      **p_indL, **p_indU ;
double   *ent, *entD ;
double   **p_entL, **p_entU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, 
           "\n error in ILUMtx_writeForMatlab(), mtx is NULL\n") ;
   return(-1) ;
}
if ( (neqns = mtx->neqns) <= 0 ) {
   fprintf(stderr, "\n error in ILUM_writeForMatlab()"
           "\n neqns = %d\n", neqns) ;
   return(-2) ;
}
if ( !(ILUMTX_IS_REAL(mtx) || ILUMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n error in ILUM_writeForMatlab()"
           "\n mtx type = %d\n", mtx->type) ;
   return(-3) ;
}
if ( !(ILUMTX_IS_SYMMETRIC(mtx) || ILUMTX_IS_HERMITIAN(mtx)
       || ILUMTX_IS_NONSYMMETRIC(mtx)) ) {
   fprintf(stderr, "\n error in ILUM_writeForMatlab()"
           "\n mtx symmetry = %d\n", mtx->symmetryflag) ;
   return(-4) ;
}
if ( !(ILUMTX_IS_L_BY_ROWS(mtx) || ILUMTX_IS_L_BY_COLUMNS(mtx)) ) {
  fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n LstorageMode = %d\n", mtx->LstorageMode) ;
   return(-5) ;
}
if ( !(ILUMTX_IS_U_BY_ROWS(mtx) || ILUMTX_IS_U_BY_COLUMNS(mtx)) ) {
  fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n UstorageMode = %d\n", mtx->UstorageMode) ;
   return(-6) ;
}
sizesU = mtx->sizesU ;
entD   = mtx->entD   ;
p_indU = mtx->p_indU ;
p_entU = mtx->p_entU ;
if ( ILUMTX_IS_SYMMETRIC(mtx) || ILUMTX_IS_HERMITIAN(mtx) ) {
   sizesL = mtx->sizesU ;
   p_indL = mtx->p_indU ;
   p_entL = mtx->p_entU ;
} else {
   sizesL = mtx->sizesL ;
   p_indL = mtx->p_indL ;
   p_entL = mtx->p_entL ;
}
if ( sizesL == NULL ) {
   fprintf(stderr, 
           "\n error in ILUM_writeForMatlab(), sizesL = NULL\n") ;
   return(-7) ;
}
if ( sizesU == NULL ) {
   fprintf(stderr, 
           "\n error in ILUM_writeForMatlab(), sizesU = NULL\n") ;
   return(-8) ;
}
if ( p_indL == NULL ) {
   fprintf(stderr, 
           "\n error in ILUM_writeForMatlab(), p_indL = NULL\n") ;
   return(-9) ;
}
if ( p_indU == NULL ) {
   fprintf(stderr, 
           "\n error in ILUM_writeForMatlab(), p_indU = NULL\n") ;
   return(-10) ;
}
if ( entD == NULL ) {
   fprintf(stderr, "\n error in ILUM_writeForMatlab(), entD = NULL\n") ;
   return(-11) ;
}
if ( p_entL == NULL ) {
   fprintf(stderr, 
           "\n error in ILUM_writeForMatlab(), p_entL = NULL\n") ;
   return(-12) ;
}
if ( p_entU == NULL ) {
   fprintf(stderr, 
           "\n error in ILUM_writeForMatlab(), p_entU = NULL\n") ;
   return(-13) ;
}
if ( Lname == NULL ) {
   fprintf(stderr, 
           "\n error in ILUMtx_writeForMatlab(), Lname is NULL\n") ;
   return(-14) ;
}
if ( Dname == NULL ) {
   fprintf(stderr, 
           "\n error in ILUMtx_writeForMatlab(), Dname is NULL\n") ;
   return(-15) ;
}
if ( Uname == NULL ) {
   fprintf(stderr, 
           "\n error in ILUMtx_writeForMatlab(), Uname is NULL\n") ;
   return(-16) ;
}
if ( fp == NULL ) {
   fprintf(stderr, 
           "\n error in ILUMtx_writeForMatlab(), fp is NULL\n") ;
   return(-17) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------
   write out the entries in D
   --------------------------
*/
fprintf(fp, "\n %s = zeros(%d,%d) ;", Dname, neqns, neqns) ;
if ( ILUMTX_IS_REAL(mtx) ) {
   for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
      fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
              Dname, ieqn + 1, ieqn + 1, entD[ieqn]) ;
   }
} else if ( ILUMTX_IS_COMPLEX(mtx) ) {
   for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
      fprintf(fp, "\n %s(%d,%d) = %24.16e + i*%24.16e ;",
              Dname, ieqn + 1, ieqn + 1, entD[2*ieqn], entD[2*ieqn+1]) ;
   }
}
/*
   --------------------------
   write out the entries in U
   --------------------------
*/
fprintf(fp, "\n %s = eye(%d,%d) ;", Uname, neqns, neqns) ;
if ( ILUMTX_IS_REAL(mtx) ) {
   if ( ILUMTX_IS_U_BY_COLUMNS(mtx) ) {
      for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
         if ( (size = sizesU[ieqn]) > 0 ) {
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = 0 ; ii < size ; ii++ ) {
               fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
                       Uname, ind[ii] + 1, ieqn + 1, ent[ii]) ;
            }
         }
      }
   } else if ( ILUMTX_IS_U_BY_ROWS(mtx) ) {
      for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
         if ( (size = sizesU[ieqn]) > 0 ) {
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = 0 ; ii < size ; ii++ ) {
               fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
                       Uname, ieqn + 1, ind[ii] + 1, ent[ii]) ;
            }
         }
      }
   }
} else if ( ILUMTX_IS_COMPLEX(mtx) ) {
   if ( ILUMTX_IS_U_BY_COLUMNS(mtx) ) {
      for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
         if ( (size = sizesU[ieqn]) > 0 ) {
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = 0 ; ii < size ; ii++ ) {
               fprintf(fp, "\n %s(%d,%d) = %24.16e + i*%24.16e ;",
                       Uname, ind[ii] + 1, ieqn + 1, 
                       ent[2*ii], ent[2*ii+1]) ;
            }
         }
      }
   } else if ( ILUMTX_IS_U_BY_ROWS(mtx) ) {
      for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
         if ( (size = sizesU[ieqn]) > 0 ) {
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = 0 ; ii < size ; ii++ ) {
               fprintf(fp, "\n %s(%d,%d) = %24.16e + i*%24.16e ;",
                       Uname, ieqn + 1, ind[ii] + 1, 
                       ent[2*ii], ent[2*ii+1]) ;
            }
         }
      }
   }
}
if ( ILUMTX_IS_SYMMETRIC(mtx) ) {
   fprintf(fp, "\n %s = transpose(%s) ;", Lname, Uname) ;
} else if ( ILUMTX_IS_HERMITIAN(mtx) ) {
   fprintf(fp, "\n %s = ctranspose(%s) ;", Lname, Uname) ;
} else if ( ILUMTX_IS_NONSYMMETRIC(mtx) ) {
/*
   --------------------------
   write out the entries in L
   --------------------------
*/
   fprintf(fp, "\n %s = eye(%d,%d) ;", Lname, neqns, neqns) ;
   if ( ILUMTX_IS_REAL(mtx) ) {
      if ( ILUMTX_IS_L_BY_COLUMNS(mtx) ) {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (size = sizesL[ieqn]) > 0 ) {
               ind = p_indL[ieqn] ;
               ent = p_entL[ieqn] ;
               for ( ii = 0 ; ii < size ; ii++ ) {
                  fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
                          Lname, ind[ii] + 1, ieqn + 1, ent[ii]) ;
               }
            }
         }
      } else if ( ILUMTX_IS_L_BY_ROWS(mtx) ) {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (size = sizesL[ieqn]) > 0 ) {
               ind = p_indL[ieqn] ;
               ent = p_entL[ieqn] ;
               for ( ii = 0 ; ii < size ; ii++ ) {
                  fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
                          Lname, ieqn + 1, ind[ii] + 1, ent[ii]) ;
               }
            }
         }
      }
   } else if ( ILUMTX_IS_COMPLEX(mtx) ) {
      if ( ILUMTX_IS_L_BY_COLUMNS(mtx) ) {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (size = sizesL[ieqn]) > 0 ) {
               ind = p_indL[ieqn] ;
               ent = p_entL[ieqn] ;
               for ( ii = 0 ; ii < size ; ii++ ) {
                  fprintf(fp, "\n %s(%d,%d) = %24.16e + i*%24.16e ;",
                          Lname, ind[ii] + 1, ieqn + 1, 
                          ent[2*ii], ent[2*ii+1]) ;
               }
            }
         }
      } else if ( ILUMTX_IS_L_BY_ROWS(mtx) ) {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (size = sizesL[ieqn]) > 0 ) {
               ind = p_indL[ieqn] ;
               ent = p_entL[ieqn] ;
               for ( ii = 0 ; ii < size ; ii++ ) {
                  fprintf(fp, "\n %s(%d,%d) = %24.16e + i*%24.16e ;",
                          Lname, ieqn + 1, ind[ii] + 1, 
                          ent[2*ii], ent[2*ii+1]) ;
               }
            }
         }
      }
   }
}
return(1) ; }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
