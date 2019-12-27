/*  copyEntriesToVector.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- copy entries to a vector. 
 
   length     -- length of dvec[]
   npivot     -- number of pivots, may be 0
   pivotsizes -- vector of pivot sizes, may be NULL
   dvec[]     -- vector to receive matrix entries
   copyflag   -- flag to denote what part of the entries to copy
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_DIAGONAL        --> copy diagonal entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
   storeflag  -- flag to denote how to store entries in dvec[]
      CHV_BY_ROWS    --> store by rows
      CHV_BY_COLUMNS --> store by columns
 
   return value -- number of entries copied
 
   created  -- 97jun05, cca
   modified -- 98feb27, cca
     cases 4-7 inserted
   -------------------------------------------------------------------
*/
int
Chv_copyEntriesToVector (
   Chv      *chv,
   int      npivot,
   int      pivotsizes[],
   int      length,
   double   *dvec,
   int      copyflag, 
   int      storeflag
) {
double   *entries ;
int      mm, ncol, nD, nent, nL, nrow, nU, symflag ;
/*
   --------------------------------------------
   check the input, get dimensions and pointers
   and check that length is large enough
   --------------------------------------------
*/
if (  chv == NULL || length < 0 || dvec == NULL ) {
   fprintf(stderr,
           "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,,%d,%d)"
           "\n bad input\n", chv, length, dvec, copyflag, storeflag) ;
   exit(-1) ;
}
switch ( copyflag ) {
case CHV_STRICT_LOWER    :
case CHV_DIAGONAL        :
case CHV_STRICT_UPPER    :
case CHV_STRICT_LOWER_11 :
case CHV_LOWER_21        :
case CHV_STRICT_UPPER_11 :
case CHV_UPPER_12        :
   break ;
default :
   fprintf(stderr,
   "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,%d,%d)"
   "\n bad input\n"
   "\n copyflag = %d, must be\n"
   "\n CHV_STRICT_LOWER    --> copy strict lower entries"
   "\n CHV_DIAGONAL        --> copy diagonal entries"
   "\n CHV_STRICT_UPPER    --> copy strict upper entries"
   "\n CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block"
   "\n CHV_LOWER_21        --> copy lower entries in (2,1) block"
   "\n CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block"
   "\n CHV_UPPER_12        --> copy upper entries in (1,2) block"
   "\n",
   chv, length, dvec, copyflag, storeflag, copyflag) ;
   exit(-1) ;
   break ;
}
switch ( storeflag ) {
case CHV_BY_ROWS    :
case CHV_BY_COLUMNS :
   break ;
default :
   fprintf(stderr,
           "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n bad input\n"
           "\n storeflag = %d, must be\n"
           "\n CHV_BY_ROWS    --> store by rows"
           "\n CHV_BY_COLUMNS --> store by columns"
           "\n",
           chv, length, dvec, copyflag, storeflag, storeflag) ;
   exit(-1) ;
   break ;
}
nD      = chv->nD      ;
nL      = chv->nL      ;
nU      = chv->nU      ;
symflag = chv->symflag ;
nrow    = nD + nL      ;
ncol    = nD + nU      ;
/*
   ------------------------------------------
   compute the number of entries to be copied
   ------------------------------------------
*/
switch ( copyflag ) {
case CHV_STRICT_LOWER : /* strictly lower entries  */
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      fprintf(stderr,
           "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n symflag = %d, copyflag = %d",
           chv, length, dvec, copyflag, storeflag, symflag, copyflag) ;
   exit(-1) ;
   }
   nent = (nD*(nD-1))/2 + nD*nL ;
   break ;
case CHV_DIAGONAL : /* diagonal entries  */
   if ( pivotsizes == NULL ) {
      nent = nD ;
   } else {
      int   ipivot ;
      for ( ipivot = nent = 0 ; ipivot < npivot ; ipivot++ ) {
         nent += (pivotsizes[ipivot]*(pivotsizes[ipivot] + 1))/2 ;
      }
   }
   break ;
case CHV_STRICT_UPPER : /* strictly upper entries  */
   if ( pivotsizes == NULL ) {
      nent = (nD*(nD-1))/2 + nD*nU ;
   } else {
      int   first, ipivot ;
      for ( ipivot = first = nent = 0 ; ipivot < npivot ; ipivot++ ) {
         nent += first * pivotsizes[ipivot] ;
         first += pivotsizes[ipivot] ;
      }
   }
   break ;
case CHV_STRICT_LOWER_11 : /* strictly lower entries in (1,1) block */
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      fprintf(stderr,
           "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n symflag = %d, copyflag = %d",
           chv, length, dvec, copyflag, storeflag, symflag, copyflag) ;
   exit(-1) ;
   }
   nent = (nD*(nD-1))/2 ;
   break ;
case CHV_LOWER_21        : /* lower entries in (2,1) block */
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      fprintf(stderr,
           "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n symflag = %d, copyflag = %d",
           chv, length, dvec, copyflag, storeflag, symflag, copyflag) ;
   exit(-1) ;
   }
   nent = nD*nL ;
   break ;
case CHV_STRICT_UPPER_11 : /* strictly upper entries in (1,1) block */
   if ( pivotsizes == NULL ) {
      nent = (nD*(nD-1))/2 ;
   } else {
      int   first, ipivot ;
      for ( ipivot = first = nent = 0 ; ipivot < npivot ; ipivot++ ) {
         nent += first * pivotsizes[ipivot] ;
         first += pivotsizes[ipivot] ;
      }
   }
   break ;
case CHV_UPPER_12        : /* upper entries in (1,2) block */
   nent = nD*nU ;
   break ;
default :
   break ;
}
if ( nent > length ) {
   fprintf(stderr,
           "\n fatal error in Chv_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n nent = %d, buffer length = %d",
           chv, length, dvec, copyflag, storeflag, nent, length) ;
   exit(-1) ;
}
/*
   --------------------------------------------
   make life simple, unit stride through dvec[]
   --------------------------------------------
*/
entries = Chv_entries(chv) ;
mm = 0 ;
switch ( copyflag ) {
case CHV_STRICT_LOWER :
/*
   ----------------------------------------
   copy lower entries, pivotsizes[] is NULL
   ----------------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS : {
/*
      -------------
      store by rows
      -------------
*/
      int   ii, jj, jjlast, kinc, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = 0 ; ii < nrow ; ii++ ) {
            kk   = kstart ;
            kinc = stride ;
            jjlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
            for ( jj = 0 ; jj <= jjlast ; jj++, mm++ ) {
               dvec[mm] = entries[kk] ;
               kk   += kinc ;
               kinc -=   2  ;
            }
            kstart-- ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = 0 ; ii < nrow ; ii++ ) {
            kk   = kstart ;
            kinc = stride ;
            jjlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
            for ( jj = 0 ; jj <= jjlast ; jj++, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
               kk   += kinc ;
               kinc -=   2  ;
            }
            kstart-- ;
         }
      }
      } break ;
   case CHV_BY_COLUMNS : {
/*
      ----------------
      store by columns
      ----------------
*/
      int   ii, jj, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ; 
            for ( ii = jj + 1 ; ii < nrow ; ii++, kk--, mm++ ) {
               dvec[mm] = entries[kk] ;
            }
            kstart += stride ;
            stride -=    2   ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ; 
            for ( ii = jj + 1 ; ii < nrow ; ii++, kk--, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
            }
            kstart += stride ;
            stride -=    2   ;
         }
      }
      } break ;
   } break ;
case CHV_DIAGONAL :
/*
   ---------------------
   copy diagonal entries
   ---------------------
*/
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      int   first, ii, jj, ipivot, kk, kstart, last, stride ;

      stride = nD + nU ;
      kstart = 0 ;
      if ( pivotsizes == NULL ) {
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0, kk = kstart ; ii < nD ; ii++, mm++ ) {
               dvec[mm] = entries[kk] ;
               kk += stride ;
               stride-- ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0, kk = kstart ; ii < nD ; ii++, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
               kk += stride ;
               stride-- ;
            }
         }
      } else {
         if ( CHV_IS_REAL(chv) ) {
            for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
               last = first + pivotsizes[ipivot] - 1 ;
               for ( ii = first ; ii <= last ; ii++ ) {
                  for ( jj = ii, kk = kstart ; 
                        jj <= last ; 
                        jj++, mm++, kk++ ) {
                     dvec[mm] = entries[kk] ;
                  }
                  kstart += stride ;
                  stride-- ;
               }
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
               last = first + pivotsizes[ipivot] - 1 ;
               for ( ii = first ; ii <= last ; ii++ ) {
                  for ( jj = ii, kk = kstart ; 
                        jj <= last ; 
                        jj++, mm++, kk++ ) {
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                  }
                  kstart += stride ;
                  stride-- ;
               }
            }
         }
      }
   } else {
      int   ii, kk, stride ;

      kk = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
            dvec[mm] = entries[kk] ;
            kk += stride ;
            stride -= 2 ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
            dvec[2*mm]   = entries[2*kk] ;
            dvec[2*mm+1] = entries[2*kk+1] ;
            kk += stride ;
            stride -= 2 ;
         }
      }
   }
   break ;
case CHV_STRICT_UPPER :
/*
   -------------------------
   copy strict upper entries
   -------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS :
/*
      -------------
      store by rows
      -------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   first, ii, ipivot, jj, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1 ; jj < ncol ; jj++, kk++, mm++ ) {
                     dvec[mm] = entries[kk] ;
                  }
                  kstart += stride ;
                  stride-- ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1 ; jj < ncol ; jj++, kk++, mm++ ) {
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                  }
                  kstart += stride ;
                  stride-- ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = first = 0 ;
                     ipivot < npivot ;
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ;
                     for ( jj = last + 1 ; 
                           jj < ncol ; 
                           jj++, kk++, mm++ ) {
                        dvec[mm] = entries[kk] ;
                     }
                     kstart += stride ;
                     stride-- ;
                  }
                  first = last + 1 ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = first = 0 ;
                     ipivot < npivot ;
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ;
                     for ( jj = last + 1 ; 
                           jj < ncol ; 
                           jj++, kk++, mm++ ) {
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                     }
                     kstart += stride ;
                     stride-- ;
                  }
                  first = last + 1 ;
               }
            }
         }
      } else {
         int   ii, jj, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 2 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ; 
               for ( jj = ii + 1 ; jj < ncol ; jj++, kk++, mm++ ) {
                  dvec[mm] = entries[kk] ;
               }
               kstart += stride ;
               stride -= 2 ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ; 
               for ( jj = ii + 1 ; jj < ncol ; jj++, kk++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
               }
               kstart += stride ;
               stride -= 2 ;
            }
         }
      }
      break ;
   case CHV_BY_COLUMNS :
/*
      ----------------
      store by columns
      ----------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   first, ii, iilast, ipivot, jj, 
               kinc, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU - 1 ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( jj = 0 ; jj < ncol ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
                     dvec[mm] = entries[kk] ;
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( jj = 0 ; jj < ncol ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = mm = first = 0 ; 
                     ipivot < npivot ; 
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = 0 ; ii < first ; ii++, mm++ ) {
                        dvec[mm] = entries[kk] ;
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                  }
                  first = last + 1 ;
               }
               for ( jj = nD ; jj < ncol ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
                     dvec[mm] = entries[kk] ;
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = mm = first = 0 ; 
                     ipivot < npivot ; 
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = 0 ; ii < first ; ii++, mm++ ) {
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                  }
                  first = last + 1 ;
               }
               for ( jj = nD ; jj < ncol ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
            }
         }
      } else {
         int   ii, iilast, jj, kinc, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               kk     = kstart ; 
               kinc   = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
                  dvec[mm] = entries[kk] ;
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               kk     = kstart ; 
               kinc   = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
            }
         }
      }
      break ;
   default :
      break ;
   } break ;
case CHV_STRICT_LOWER_11 :
/*
   --------------------------------------------------------------
   copy strict lower entries in (1,1) block, pivotsizes[] is NULL
   --------------------------------------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS : {
/*
      -------------
      store by rows
      -------------
*/
      int   ii, jj, jjlast, kinc, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = 0 ; ii < nD ; ii++ ) {
            kk     = kstart ;
            kinc   = stride ;
            jjlast = ii - 1 ;
            for ( jj = 0 ; jj <= jjlast ; jj++, mm++ ) {
               dvec[mm] = entries[kk] ;
               kk   += kinc ;
               kinc -=   2  ;
            }
            kstart-- ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = 0 ; ii < nD ; ii++ ) {
            kk     = kstart ;
            kinc   = stride ;
            jjlast = ii - 1 ;
            for ( jj = 0 ; jj <= jjlast ; jj++, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
               kk   += kinc ;
               kinc -=   2  ;
            }
            kstart-- ;
         }
      }
      } break ;
   case CHV_BY_COLUMNS : {
/*
      ----------------
      store by columns
      ----------------
*/
      int   ii, jj, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ; 
            for ( ii = jj + 1 ; ii < nD ; ii++, kk--, mm++ ) {
               dvec[mm] = entries[kk] ;
            }
            kstart += stride ;
            stride -=    2   ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ; 
            for ( ii = jj + 1 ; ii < nD ; ii++, kk--, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
            }
            kstart += stride ;
            stride -=    2   ;
         }
      }
      } break ;
   } break ;
case CHV_LOWER_21        :
/*
   ---------------------------------
   copy lower entries in (2,1) block
   ---------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS : {
/*
      -------------
      store by rows
      -------------
*/
      int   ii, jj, kinc, kk, kstart, stride ;

      kstart = nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = nD, mm = 0 ; ii < nrow ; ii++ ) {
            kk   = kstart ;
            kinc = stride ;
            for ( jj = 0 ; jj < nD ; jj++, mm++ ) {
               dvec[mm] = entries[kk] ;
               kk   += kinc ;
               kinc -=   2  ;
            }
            kstart-- ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = nD, mm = 0 ; ii < nrow ; ii++ ) {
            kk   = kstart ;
            kinc = stride ;
            for ( jj = 0 ; jj < nD ; jj++, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
               kk   += kinc ;
               kinc -=   2  ;
            }
            kstart-- ;
         }
      }
      } break ;
   case CHV_BY_COLUMNS : {
/*
      ----------------
      store by columns
      ----------------
*/
      int   ii, jj, kk, kstart, stride ;

      kstart = nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart ; 
            for ( ii = nD ; ii < nrow ; ii++, kk--, mm++ ) {
               dvec[mm] = entries[kk] ;
            }
            kstart += stride ;
            stride -=    2   ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart ; 
            for ( ii = nD ; ii < nrow ; ii++, kk--, mm++ ) {
               dvec[2*mm]   = entries[2*kk] ;
               dvec[2*mm+1] = entries[2*kk+1] ;
            }
            kstart += stride ;
            stride -=    2   ;
         }
      }
      } break ;
   } break ;
case CHV_STRICT_UPPER_11 :
/*
   ----------------------------------------
   copy strict upper entries in (1,1) block
   ----------------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS :
/*
      -------------
      store by rows
      -------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   first, ii, ipivot, jj, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1 ; jj < nD ; jj++, kk++, mm++ ) {
                     dvec[mm] = entries[kk] ;
                  }
                  kstart += stride ;
                  stride-- ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1 ; jj < nD ; jj++, kk++, mm++ ) {
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                  }
                  kstart += stride ;
                  stride-- ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = first = 0 ;
                     ipivot < npivot ;
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ;
                     for ( jj = last + 1 ; jj < nD ; jj++, kk++, mm++ ){
                        dvec[mm] = entries[kk] ;
                     }
                     kstart += stride ;
                     stride-- ;
                  }
                  first = last + 1 ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = first = 0 ;
                     ipivot < npivot ;
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ;
                     for ( jj = last + 1 ; jj < nD ; jj++, kk++, mm++ ){
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                     }
                     kstart += stride ;
                     stride-- ;
                  }
                  first = last + 1 ;
               }
            }
         }
      } else {
         int   ii, jj, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 2 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ; 
               for ( jj = ii + 1 ; jj < nD ; jj++, kk++, mm++ ) {
                  dvec[mm] = entries[kk] ;
               }
               kstart += stride ;
               stride -= 2 ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ; 
               for ( jj = ii + 1 ; jj < nD ; jj++, kk++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
               }
               kstart += stride ;
               stride -= 2 ;
            }
         }
      }
      break ;
   case CHV_BY_COLUMNS :
/*
      ----------------
      store by columns
      ----------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   first, ii, iilast, ipivot, jj, 
               kinc, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU - 1 ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( jj = 0 ; jj < nD ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
                     dvec[mm] = entries[kk] ;
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( jj = 0 ; jj < nD ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = mm = first = 0 ; 
                     ipivot < npivot ; 
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = 0 ; ii < first ; ii++, mm++ ) {
                        dvec[mm] = entries[kk] ;
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                  }
                  first = last + 1 ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = mm = first = 0 ; 
                     ipivot < npivot ; 
                     ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = 0 ; ii < first ; ii++, mm++ ) {
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                  }
                  first = last + 1 ;
               }
            }
         }
      } else {
         int   ii, jj, kinc, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = 0 ; jj < nD ; jj++ ) {
               kk     = kstart ; 
               kinc   = stride ;
               for ( ii = 0 ; ii < jj ; ii++, mm++ ) {
                  dvec[mm] = entries[kk] ;
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = 0 ; jj < nD ; jj++ ) {
               kk     = kstart ; 
               kinc   = stride ;
               for ( ii = 0 ; ii < jj ; ii++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
            }
         }
      }
      break ;
   default :
      break ;
   } break ;
case CHV_UPPER_12        :
/*
   ---------------------------------
   copy upper entries in (1,2) block
   ---------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS :
/*
      -------------
      store by rows
      -------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   ii, jj, kk, kstart, stride ;

         kstart = nD ;
         stride = nD + nU - 1 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart ;
               for ( jj = 0 ; jj < nU ; jj++, kk++, mm++ ) {
                  dvec[mm] = entries[kk] ;
               }
               kstart += stride ;
               stride-- ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart ;
               for ( jj = 0 ; jj < nU ; jj++, kk++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
               }
               kstart += stride ;
               stride-- ;
            }
         }
      } else {
         int   ii, jj, kk, kstart, stride ;

         kstart = 2*nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart ; 
               for ( jj = nD ; jj < ncol ; jj++, kk++, mm++ ) {
                  dvec[mm] = entries[kk] ;
               }
               kstart += stride ;
               stride -= 2 ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart ; 
               for ( jj = nD ; jj < ncol ; jj++, kk++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
               }
               kstart += stride ;
               stride -= 2 ;
            }
         }
      }
      break ;
   case CHV_BY_COLUMNS :
/*
      ----------------
      store by columns
      ----------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   ii, jj, kinc, kk, kstart, stride ;

         kstart = nD ;
         stride = nD + nU - 1 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
                  dvec[mm] = entries[kk] ;
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         }
      } else {
         int   ii, jj, kinc, kk, kstart, stride ;

         kstart = 2*nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk     = kstart ; 
               kinc   = stride ;
               for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
                  dvec[mm] = entries[kk] ;
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk     = kstart ; 
               kinc   = stride ;
               for ( ii = 0 ; ii < nD ; ii++, mm++ ) {
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
            }
         }
      }
      break ;
   default :
      break ;
   } break ;
default :
   break ;
}
return(mm) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- copy large entries to a vector. the portion copied
              can be a union of the strict lower portion,
              the diagonal portion, and the strict upper
              portion. there is one restriction, if the strict
              lower and strict upper are to be copied, the
              diagonal will also be copied.
 
   npivot     -- number of pivots, may be 0
   pivotsizes -- vector of pivot sizes, may be NULL
   sizes[]    -- vector to receive row/column sizes
   ivec[]     -- vector to receive row/column indices
   dvec[]     -- vector to receive matrix entries
   copyflag   -- flag to denote what part of the entries to copy
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
   storeflag  -- flag to denote how to store entries in dvec[]
      CHV_BY_ROWS    --> store by rows
      CHV_BY_COLUMNS --> store by columns
   droptol    -- entry to be copied must be larger than this magnitude
 
   return value -- number of entries copied
 
   created  -- 97jun05, cca
   modified -- 97feb27, cca
      cases 4-7 inserted
   -------------------------------------------------------------------
*/
int
Chv_copyBigEntriesToVector (
   Chv     *chv,
   int      npivot,
   int      pivotsizes[],
   int      sizes[],
   int      ivec[],
   double   dvec[],
   int      copyflag, 
   int      storeflag,
   double   droptol
) {
double   absval ;
double   *entries ;
int      mm, ncol, nD, nL, nrow, nU, symflag ;
/*
   --------------------------------------------
   check the input, get dimensions and pointers
   --------------------------------------------
*/
if (  chv == NULL || sizes == NULL || ivec == NULL || dvec == NULL ) {
   fprintf(stderr,
     "\n fatal error in Chv_copyBigEntriesToVector()"
     "\n chv %p, sizes %p, ivec %p, dvec %p"
     "\n bad input\n", chv, sizes, ivec, dvec) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, 
        "\n\n %% inside Chv_copyBigEntries, copyflag %d, storeflag %d",
        copyflag, storeflag) ;
fflush(stdout) ;
#endif
switch ( copyflag ) {
case CHV_STRICT_LOWER :
case CHV_STRICT_UPPER :
case CHV_STRICT_LOWER_11 :
case CHV_LOWER_21        :
case CHV_STRICT_UPPER_11 :
case CHV_UPPER_12        :
   break ;
default :
   fprintf(stderr,
      "\n fatal error in Chv_copyBigEntriesToVector(%p,%p,%p,%p,%d,%d)"
        "\n bad input\n"
        "\n copyflag = %d, must be\n"
        "\n    1 --> strictly lower entries"
        "\n    3 --> strictly upper entries"
        "\n    4 --> copy strict lower entries in (1,1) block"
        "\n    5 --> copy lower entries in (2,1) block"
        "\n    6 --> copy strict upper entries in (1,1) block"
        "\n    7 --> copy upper entries in (1,2) block",
        chv, sizes, ivec, dvec, copyflag, storeflag, copyflag) ;
   exit(-1) ;
}
if ( storeflag  < 0 || storeflag > 1 ) {
   fprintf(stderr,
      "\n fatal error in Chv_copyBigEntriesToVector(%p,%p,%p,%p,%d,%d)"
      "\n bad input\n"
      "\n storeflag = %d, must be\n"
      "\n    0 --> store by rows"
      "\n    1 --> store by columns",
      chv, sizes, ivec, dvec, copyflag, storeflag, storeflag) ;
   exit(-1) ;
}
nD      = chv->nD      ;
nL      = chv->nL      ;
nU      = chv->nU      ;
symflag = chv->symflag ;
nrow    = nD + nL      ;
ncol    = nD + nU      ;
/*
   --------------------------------------------
   make life simple, unit stride through dvec[]
   --------------------------------------------
*/
entries = Chv_entries(chv) ;
mm = 0 ;
switch ( copyflag ) {
case CHV_STRICT_LOWER :
/*
   ----------------------------------------
   copy lower entries, pivotsizes[] is NULL
   ----------------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS : {
/*
      -------------
      store by rows
      -------------
*/
      int   count, ii, jj, jjlast, kinc, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = 0 ; ii < nrow ; ii++ ) {
            kk     = kstart ;
            kinc   = stride ;
            jjlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
            for ( jj = count = 0 ; jj <= jjlast ; jj++ ) {
               absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  ivec[mm] = jj ;
                  dvec[mm] = entries[kk] ;
                  mm++, count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            sizes[ii] = count ;
            kstart-- ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = 0 ; ii < nrow ; ii++ ) {
            kk     = kstart ;
            kinc   = stride ;
            jjlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
            for ( jj = count = 0 ; jj <= jjlast ; jj++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  ivec[mm] = jj ;
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  mm++, count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            sizes[ii] = count ;
            kstart-- ;
         }
      }
      } break ;
   case CHV_BY_COLUMNS : {
/*
      ----------------
      store by columns
      ----------------
*/
      int   count, ii, jj, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ;
            for ( ii = jj + 1, count = 0 ; ii < nrow ; ii++, kk-- ) {
               absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  ivec[mm] = ii ;
                  dvec[mm] = entries[kk] ;
                  mm++, count++ ;
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[jj] = count ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ;
            for ( ii = jj + 1, count = 0 ; ii < nrow ; ii++, kk-- ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  ivec[mm] = ii ;
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  mm++, count++ ;
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[jj] = count ;
         }
      }
      } break ;
   }
   break ;
case CHV_STRICT_UPPER :
/*
   -------------------------
   copy strict upper entries
   -------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS :
/*
      -------------
      store by rows
      -------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   count, first, ii, ipivot, jj, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1, count = 0 ; 
                        jj < ncol ; 
                        jj++, kk++ ) {
                     absval = fabs(entries[kk]) ;
                     if ( absval >= droptol ) {
                        ivec[mm] = jj ;
                        dvec[mm] = entries[kk] ;
                        mm++, count++ ;
                     }
                  }
                  kstart += stride ;
                  stride-- ;
                  sizes[ii] = count ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1, count = 0 ; 
                        jj < ncol ; 
                        jj++, kk++ ) {
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                     if ( absval >= droptol ) {
                        ivec[mm] = jj ;
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        mm++, count++ ;
                     }
                  }
                  kstart += stride ;
                  stride-- ;
                  sizes[ii] = count ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ; 
                     for ( jj = last + 1, count = 0 ; 
                           jj < ncol ; 
                           jj++, kk++ ) {
                        absval = fabs(entries[kk]) ;
                        if ( absval >= droptol ) {
                           ivec[mm] = jj ;
                           dvec[mm] = entries[kk] ;
                           mm++, count++ ;
                        }
                     }
                     kstart += stride ;
                     stride-- ;
                     sizes[ii] = count ;
                  }
                  first = last + 1 ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ; 
                     for ( jj = last + 1, count = 0 ; 
                           jj < ncol ; 
                           jj++, kk++ ) {
                        absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                        if ( absval >= droptol ) {
                           ivec[mm] = jj ;
                           dvec[2*mm]   = entries[2*kk] ;
                           dvec[2*mm+1] = entries[2*kk+1] ;
                           mm++, count++ ;
                        }
                     }
                     kstart += stride ;
                     stride-- ;
                     sizes[ii] = count ;
                  }
                  first = last + 1 ;
               }
            }
         }
      } else {
         int   count, ii, jj, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 2 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ;
               for ( jj = ii + 1, count = 0 ; jj < ncol ; jj++, kk++ ) {
                  absval = fabs(entries[kk]) ;
                  if ( absval >= droptol ) {
                     ivec[mm] = jj ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
               }
               kstart += stride ;
               stride -= 2 ;
               sizes[ii] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ;
               for ( jj = ii + 1, count = 0 ; jj < ncol ; jj++, kk++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                  if ( absval >= droptol ) {
                     ivec[mm] = jj ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
               }
               kstart += stride ;
               stride -= 2 ;
               sizes[ii] = count ;
            }
         }
      } break ;
   case CHV_BY_COLUMNS :
/*
      ----------------
      store by columns
      ----------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   count, first, ii, iilast, ipivot, jj, 
               kinc, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU - 1 ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( jj = 0 ; jj < ncol ; jj++ ) {
                  kk     = kstart ;
                  kinc   = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                     absval = fabs(entries[kk]) ;
                     if ( absval >= droptol ) {
                        ivec[mm] = ii ;
                        dvec[mm] = entries[kk] ;
                        mm++, count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
                  sizes[jj] = count ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( jj = 0 ; jj < ncol ; jj++ ) {
                  kk     = kstart ;
                  kinc   = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                     if ( absval >= droptol ) {
                        ivec[mm] = ii ;
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        mm++, count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
                  sizes[jj] = count ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = count = 0 ; ii < first ; ii++ ) {
                        absval = fabs(entries[kk]) ;
                        if ( absval >= droptol ) {
                           ivec[mm] = ii ;
                           dvec[mm] = entries[kk] ;
                           mm++, count++ ;
                        }
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                     sizes[jj] = count ;
                  }
                  first = last + 1 ;
               }
               for ( jj = nD ; jj < ncol ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = count = 0 ; ii < nD ; ii++ ) {
                     absval = fabs(entries[kk]) ;
                     if ( absval >= droptol ) {
                        ivec[mm] = ii ;
                        dvec[mm] = entries[kk] ;
                        mm++, count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
                  sizes[jj] = count ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = count = 0 ; ii < first ; ii++ ) {
                        absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                        if ( absval >= droptol ) {
                           ivec[mm] = ii ;
                           dvec[2*mm]   = entries[2*kk] ;
                           dvec[2*mm+1] = entries[2*kk+1] ;
                           mm++, count++ ;
                        }
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                     sizes[jj] = count ;
                  }
                  first = last + 1 ;
               }
               for ( jj = nD ; jj < ncol ; jj++ ) {
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = count = 0 ; ii < nD ; ii++ ) {
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                     if ( absval >= droptol ) {
                        ivec[mm] = ii ;
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        mm++, count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
                  sizes[jj] = count ;
               }
            }
         }
      } else {
         int   count, ii, iilast, jj, kinc, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               kk = kstart ; 
               kinc = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                  absval = fabs(entries[kk]) ;
                  if ( absval >= droptol ) {
                     ivec[mm] = ii ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
               sizes[jj] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               kk = kstart ; 
               kinc = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                  if ( absval >= droptol ) {
                     ivec[mm] = ii ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
               sizes[jj] = count ;
            }
         }
      }
      break ;
   default :
      break ;
   }
   break ;
case CHV_STRICT_LOWER_11 :
/*
   ------------------------------------------------
   copy entries in strict lower part of (1,1) block
   ------------------------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS : {
/*
      -------------
      store by rows
      -------------
*/
      int   count, ii, jj, jjlast, kinc, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = 0 ; ii < nD ; ii++ ) {
            kk     = kstart ;
            kinc   = stride ;
            jjlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
            for ( jj = count = 0 ; jj <= jjlast ; jj++ ) {
               absval = fabs(entries[kk]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = jj ;
                  dvec[mm] = entries[kk] ;
                  mm++, count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            sizes[ii] = count ;
            kstart-- ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = 0 ; ii < nD ; ii++ ) {
            kk     = kstart ;
            kinc   = stride ;
            jjlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
            for ( jj = count = 0 ; jj <= jjlast ; jj++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = jj ;
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  mm++, count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            sizes[ii] = count ;
            kstart-- ;
         }
      }
      } break ;
   case CHV_BY_COLUMNS : {
/*
      ----------------
      store by columns
      ----------------
*/
      int   count, ii, jj, kk, kstart, stride ;

      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ;
            for ( ii = jj + 1, count = 0 ; ii < nD ; ii++, kk-- ) {
               absval = fabs(entries[kk]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = ii ;
                  dvec[mm] = entries[kk] ;
                  mm++, count++ ;
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[jj] = count ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart - 1 ;
            for ( ii = jj + 1, count = 0 ; ii < nD ; ii++, kk-- ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = ii ;
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  mm++, count++ ;
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[jj] = count ;
         }
      }
      } break ;
   }
   break ;
case CHV_LOWER_21        :
/*
   ---------------------------
   copy entries in (2,1) block
   ---------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS : {
/*
      -------------
      store by rows
      -------------
*/
      int   count, ii, jj, kinc, kk, kstart, stride ;

      kstart = nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( ii = nD ; ii < nrow ; ii++ ) {
            kk   = kstart ;
            kinc = stride ;
            for ( jj = count = 0 ; jj < nD ; jj++ ) {
               absval = fabs(entries[kk]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = jj ;
                  dvec[mm] = entries[kk] ;
                  mm++, count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            sizes[ii-nD] = count ;
            kstart-- ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( ii = nD ; ii < nrow ; ii++ ) {
            kk   = kstart ;
            kinc = stride ;
            for ( jj = count = 0 ; jj < nD ; jj++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = jj ;
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  mm++, count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            sizes[ii-nD] = count ;
            kstart-- ;
         }
      }
      } break ;
   case CHV_BY_COLUMNS : {
/*
      ----------------
      store by columns
      ----------------
*/
      int   count, ii, jj, kk, kstart, stride ;

      kstart = nL - 1 ;
      stride = 2*nD + nL + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart ;
            for ( ii = nD, count = 0 ; ii < nrow ; ii++, kk-- ) {
               absval = fabs(entries[kk]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = ii ;
                  dvec[mm] = entries[kk] ;
                  mm++, count++ ;
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[jj] = count ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk = kstart ;
            for ( ii = nD, count = 0 ; ii < nrow ; ii++, kk-- ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
               if ( absval >= droptol ) {
                  ivec[mm] = ii ;
                  dvec[2*mm]   = entries[2*kk] ;
                  dvec[2*mm+1] = entries[2*kk+1] ;
                  mm++, count++ ;
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[jj] = count ;
         }
      }
      } break ;
   }
   break ;
case CHV_STRICT_UPPER_11 :
/*
   ----------------------------------------
   copy strict upper entries in (1,1) block
   ----------------------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS :
/*
      -------------
      store by rows
      -------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   count, first, ii, ipivot, jj, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1, count = 0 ; jj < nD ; jj++, kk++ ){
                     absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                     fprintf(stdout, 
                             "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                             ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                        fprintf(stdout, ", keep") ;
#endif
                        ivec[mm] = jj ;
                        dvec[mm] = entries[kk] ;
                        mm++, count++ ;
                     }
                  }
                  kstart += stride ;
                  stride-- ;
                  sizes[ii] = count ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  kk = kstart + 1 ;
                  for ( jj = ii + 1, count = 0 ; jj < nD ; jj++, kk++ ){
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                     fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                             ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                        fprintf(stdout, ", keep") ;
#endif
                        ivec[mm] = jj ;
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        mm++, count++ ;
                     }
                  }
                  kstart += stride ;
                  stride-- ;
                  sizes[ii] = count ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ; 
                     for ( jj = last + 1, count = 0 ; 
                           jj < nD ; 
                           jj++, kk++ ) {
                        absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                        if ( absval >= droptol ) {
#if MYDEBUG > 0
                           fprintf(stdout, ", keep") ;
#endif
                           ivec[mm] = jj ;
                           dvec[mm] = entries[kk] ;
                           mm++, count++ ;
                        }
                     }
                     kstart += stride ;
                     stride-- ;
                     sizes[ii] = count ;
                  }
                  first = last + 1 ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( ii = first ; ii <= last ; ii++ ) {
                     kk = kstart + last - ii + 1 ; 
                     for ( jj = last + 1, count = 0 ; 
                           jj < nD ; 
                           jj++, kk++ ) {
                        absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                        if ( absval >= droptol ) {
#if MYDEBUG > 0
                           fprintf(stdout, ", keep") ;
#endif
                           ivec[mm] = jj ;
                           dvec[2*mm]   = entries[2*kk] ;
                           dvec[2*mm+1] = entries[2*kk+1] ;
                           mm++, count++ ;
                        }
                     }
                     kstart += stride ;
                     stride-- ;
                     sizes[ii] = count ;
                  }
                  first = last + 1 ;
               }
            }
         }
      } else {
         int   count, ii, jj, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 2 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ;
               for ( jj = ii + 1, count = 0 ; jj < nD ; jj++, kk++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% A(%d,%d) = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = jj ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
               }
               kstart += stride ;
               stride -= 2 ;
               sizes[ii] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart + 1 ;
               for ( jj = ii + 1, count = 0 ; jj < nD ; jj++, kk++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% A(%d,%d) = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = jj ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
               }
               kstart += stride ;
               stride -= 2 ;
               sizes[ii] = count ;
            }
         }
      }
      break ;
   case CHV_BY_COLUMNS :
/*
      ----------------
      store by columns
      ----------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   count, first, ii, iilast, ipivot, jj, 
               kinc, kk, kstart, last, stride ;

         kstart = 0 ;
         stride = nD + nU - 1 ;
         if ( pivotsizes == NULL ) {
            if ( CHV_IS_REAL(chv) ) {
               for ( jj = 0 ; jj < nD ; jj++ ) {
                  kk     = kstart ;
                  kinc   = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                     absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                     fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                             ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                        ivec[mm] = ii ;
                        dvec[mm] = entries[kk] ;
                        mm++, count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
                  sizes[jj] = count ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( jj = 0 ; jj < nD ; jj++ ) {
                  kk     = kstart ;
                  kinc   = stride ;
                  iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
                  for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                     fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                             ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                        ivec[mm] = ii ;
                        dvec[2*mm]   = entries[2*kk] ;
                        dvec[2*mm+1] = entries[2*kk+1] ;
                        mm++, count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
                  sizes[jj] = count ;
               }
            }
         } else {
            if ( CHV_IS_REAL(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = count = 0 ; ii < first ; ii++ ) {
                        absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                        if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                           ivec[mm] = ii ;
                           dvec[mm] = entries[kk] ;
                           mm++, count++ ;
                        }
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                     sizes[jj] = count ;
                  }
                  first = last + 1 ;
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
                  last = first + pivotsizes[ipivot] - 1 ;
                  for ( jj = first ; jj <= last ; jj++ ) {
                     kk   = kstart ;
                     kinc = stride ;
                     for ( ii = count = 0 ; ii < first ; ii++ ) {
                        absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                        if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                           ivec[mm] = ii ;
                           dvec[2*mm]   = entries[2*kk] ;
                           dvec[2*mm+1] = entries[2*kk+1] ;
                           mm++, count++ ;
                        }
                        kk += kinc ;
                        kinc-- ;
                     }
                     kstart++ ;
                     sizes[jj] = count ;
                  }
                  first = last + 1 ;
               }
            }
         }
      } else {
         int   count, ii, iilast, jj, kinc, kk, kstart, stride ;

         kstart = nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = 0 ; jj < nD ; jj++ ) {
               kk = kstart ; 
               kinc = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = ii ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
               sizes[jj] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = 0 ; jj < nD ; jj++ ) {
               kk = kstart ; 
               kinc = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = count = 0 ; ii <= iilast ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = ii ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
               sizes[jj] = count ;
            }
         }
      } break ;
   default :
      break ;
   }
   #if MYDEBUG > 0
   fprintf(stdout, 
           "\n %% break from case CHV_STRICT_UPPER_11, mm = %d", mm) ;
   #endif
   break ;
case CHV_UPPER_12        :
/*
   ---------------------------
   copy entries in (1,2) block
   ---------------------------
*/
   switch ( storeflag ) {
   case CHV_BY_ROWS :
/*
      -------------
      store by rows
      -------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   count, ii, jj, kk, kstart, stride ;

         kstart = nD ;
         stride = nD + nU - 1 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart ;
               for ( jj = nD, count = 0 ; jj < ncol ; jj++, kk++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = jj ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
               }
               kstart += stride ;
               stride-- ;
               sizes[ii] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ii = 0 ; ii < nD ; ii++ ) {
               kk = kstart ;
               for ( jj = nD, count = 0 ; jj < ncol ; jj++, kk++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = jj ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
               }
               kstart += stride ;
               stride-- ;
               sizes[ii] = count ;
            }
         }
      } else {
         int   count, ii, jj, kk, kstart, stride ;

         kstart = 2*nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            kk = kstart ;
            if ( CHV_IS_REAL(chv) ) {
               for ( jj = nD, count = 0 ; jj < ncol ; jj++, kk++ ) {
                  absval = fabs(entries[kk]) ;
                  if ( absval >= droptol ) {
                     ivec[mm] = jj ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
               }
            } else if ( CHV_IS_COMPLEX(chv) ) {
               for ( jj = nD, count = 0 ; jj < ncol ; jj++, kk++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
                  if ( absval >= droptol ) {
                     ivec[mm] = jj ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
               }
            }
            kstart += stride ;
            stride -= 2 ;
            sizes[ii] = count ;
         }
      }
      break ;
   case CHV_BY_COLUMNS :
/*
      ----------------
      store by columns
      ----------------
*/
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         int   count, ii, jj, kinc, kk, kstart, stride ;

         kstart = nD ;
         stride = nD + nU - 1 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk     = kstart ;
               kinc   = stride ;
               for ( ii = count = 0 ; ii < nD ; ii++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = ii ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
               sizes[jj-nD] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk     = kstart ;
               kinc   = stride ;
               for ( ii = count = 0 ; ii < nD ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = ii ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
               sizes[jj-nD] = count ;
            }
         }
      } else {
         int   count, ii, jj, kinc, kk, kstart, stride ;

         kstart = 2*nD + nL - 1 ;
         stride = 2*nD + nL + nU - 3 ;
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk = kstart ; 
               kinc = stride ;
               for ( ii = count = 0 ; ii < nD ; ii++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = ii ;
                     dvec[mm] = entries[kk] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
               sizes[jj-nD] = count ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk = kstart ; 
               kinc = stride ;
               for ( ii = count = 0 ; ii < nD ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     ivec[mm] = ii ;
                     dvec[2*mm]   = entries[2*kk] ;
                     dvec[2*mm+1] = entries[2*kk+1] ;
                     mm++, count++ ;
                  }
                  kk += kinc ;
                  kinc -= 2 ;
               }
               kstart++ ;
               sizes[jj-nD] = count ;
            }
         }
      }
      break ;
   default :
      break ;
   }
   break ;
default :
   break ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n %% return, mm = %d", mm) ;
#endif

return(mm) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- return the number of entries 
              in a portion of the object

   countflag -- which entries to count
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
 
   created -- 98feb27, cca
   -------------------------------------------------------------------
*/
int
Chv_countEntries (
   Chv     *chv,
   int      npivot,
   int      pivotsizes[],
   int      countflag
) {
int      count, nD, nL, nU ;
/*
   --------------------------------------------
   check the input, get dimensions and pointers
   --------------------------------------------
*/
if (  chv == NULL ) {
   fprintf(stderr,
     "\n fatal error in Chv_countEntries(%p,%d,%p,%d)"
     "\n bad input\n", chv, npivot, pivotsizes, countflag) ;
   exit(-1) ;
}
if ( countflag < 1 || countflag > 7 ) {
   fprintf(stderr,
        "\n fatal error in Chv_countEntries(%p,%d,%p,%d)"
        "\n bad input\n"
        "\n countflag = %d, must be\n"
        "\n    1 --> strictly lower entries"
        "\n    2 --> diagonal entries"
        "\n    3 --> strictly upper entries"
        "\n    4 --> strictly lower entries in (1,1) block"
        "\n    5 --> lower entries in (2,1) block"
        "\n    6 --> strictly upper entries in (1,1) block"
        "\n    7 --> upper entries in (1,2) block",
        chv, npivot, pivotsizes, countflag, countflag) ;
   exit(-1) ;
}
if ( (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv))
   && (countflag == 1 || countflag == 4 || countflag == 5 ) ) {
   fprintf(stderr,
        "\n fatal error in Chv_countEntries(%p,%d,%p,%d)"
        "\n countflag = %d --> lower entries but chevron is symmetric",
        chv, npivot, pivotsizes, countflag, countflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
switch ( countflag ) {
case CHV_STRICT_LOWER : /* strictly lower entries */
   count = (nD*(nD-1))/2 + nD*nL ;
   break ;
case CHV_DIAGONAL : /* diagonal entries */
   if ( (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv))
      && pivotsizes != NULL ) {
      count = 2*nD - npivot ;
   } else {
      count = nD ;
   }
   break ;
case CHV_STRICT_UPPER : /* strictly upper entries */
   if ( (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv))
      && pivotsizes != NULL ) {
      count = (nD*(nD-1))/2 - nD + npivot + nD*nU ;
   } else {
      count = (nD*(nD-1))/2 + nD*nU ;
   }
   break ;
case CHV_STRICT_LOWER_11 : /* strictly lower entries in (1,1) block */
   count = (nD*(nD-1))/2 ;
   break ;
case CHV_LOWER_21        : /* lower entries in (2,1) block */
   count = nD*nL ;
   break ;
case CHV_STRICT_UPPER_11 : /* strictly upper entries in (1,1) block */
   if ( (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv))
      && pivotsizes != NULL ) {
      count = (nD*(nD-1))/2 - nD + npivot ;
   } else {
      count = (nD*(nD-1))/2 ;
   }
   break ;
case CHV_UPPER_12        : /* upper entries in (1,2) block */
   count = nD*nU ;
   break ;
}

return(count) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- return the number of entries 
   whose magnitude is larger than droptol.

   countflag -- which entries to count
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
 
   created  -- 97jun07, cca
   modified -- 98feb27, cca
     cases 4-7 inserted
   -------------------------------------------------------------------
*/
int
Chv_countBigEntries (
   Chv     *chv,
   int      npivot,
   int      pivotsizes[],
   int      countflag,
   double   droptol
) {
double   absval ;
double   *entries ;
int      count, nD, nL, nU ;
/*
   --------------------------------------------
   check the input, get dimensions and pointers
   --------------------------------------------
*/
if (  chv == NULL ) {
   fprintf(stderr,
     "\n fatal error in Chv_countBigEntries(%p,%d,%p,%d,%f)"
     "\n bad input\n", chv, npivot, pivotsizes, countflag, droptol) ;
   exit(-1) ;
}
switch ( countflag ) {
case CHV_STRICT_LOWER :
case CHV_STRICT_UPPER :
case CHV_STRICT_LOWER_11 :
case CHV_LOWER_21        :
case CHV_STRICT_UPPER_11 :
case CHV_UPPER_12        :
   break ;
default :
   fprintf(stderr,
        "\n fatal error in Chv_countBigEntries(%p,%d,%p,%d,%f)"
        "\n bad input\n"
        "\n countflag = %d, must be\n"
        "\n    1 --> strictly lower entries"
        "\n    3 --> strictly upper entries"
        "\n    4 --> count strict lower entries in (1,1) block"
        "\n    5 --> count lower entries in (2,1) block"
        "\n    6 --> count strict upper entries in (1,1) block"
        "\n    7 --> count upper entries in (1,2) block",
        chv, npivot, pivotsizes, countflag, droptol, countflag) ;
   exit(-1) ;
}
if ( (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv))
   && (countflag == 1 || countflag == 4 || countflag == 5) ) {
   fprintf(stderr,
        "\n fatal error in Chv_countBigEntries(%p,%d,%p,%d,%f)"
        "\n countflag = %d --> lower entries but chevron is symmetric",
        chv, npivot, pivotsizes, countflag, droptol, countflag) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n %% inside Chv_countBigEntries, countflag %d",
        countflag) ;
fflush(stdout) ;
#endif
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
switch ( countflag ) {
case CHV_STRICT_LOWER : {
   int   ii, jj, jlast, kinc, kk, kstart, nrow, stride ;
/*
      -----------------------------------------------------
      count the number of big entries in the lower triangle
      -----------------------------------------------------
*/
   nrow   = nD + nL ;
   count  = 0 ;
   kstart = nD + nL - 1 ;
   stride = 2*nD + nL + nU - 1 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( ii = 0 ; ii < nrow ; ii++ ) {
         kk   = kstart ;
         kinc = stride ;
         jlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
         for ( jj = 0 ; jj <= jlast ; jj++ ) {
            absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
            if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
               count++ ;
            }
            kk   += kinc ;
            kinc -=   2  ;
         }
         kstart-- ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( ii = 0 ; ii < nrow ; ii++ ) {
         kk   = kstart ;
         kinc = stride ;
         jlast = (ii < nD) ? (ii - 1) : (nD - 1) ;
         for ( jj = 0 ; jj <= jlast ; jj++ ) {
            absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
            if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
               count++ ;
            }
            kk   += kinc ;
            kinc -=   2  ;
         }
         kstart-- ;
      }
   }
   } break ;
case CHV_STRICT_UPPER : {
   int   first, ii, iilast, ipivot, jj, kinc, kk, kstart, 
         last, ncol, stride ;

   ncol = nD + nU ;
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n symmetric chevron, npivot = %d, pivotsizes = %p", 
              npivot, pivotsizes) ;
#endif
/*
      --------------------
      chevron is symmetric
      --------------------
*/
      count  = 0 ;
      kstart = 0 ;
      stride = nD + nU - 1 ;
      if ( pivotsizes == NULL ) {
/*
         -------------
         D is diagonal
         -------------
*/
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = 0 ; ii <= iilast ; ii++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
               for ( ii = 0 ; ii <= iilast ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         }
      } else {
/*
         ---------------------------------
         D can have 1 x 1 and 2 x 2 pivots
         ---------------------------------
*/
         if ( CHV_IS_REAL(chv) ) {
            for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
               last = first + pivotsizes[ipivot] - 1 ;
               for ( jj = first ; jj <= last ; jj++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, 
          "\n ipivot = %d, jj = %d, first = %d, last = %d, kstart = %d",
                       ipivot, jj, first, last, kstart) ;
#endif
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = 0 ; ii < first ; ii++ ) {
                     absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                        fprintf(stdout, ", keep") ;
#endif
                        count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
               first = last + 1 ;
            }
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
               last = first + pivotsizes[ipivot] - 1 ;
               for ( jj = first ; jj <= last ; jj++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, 
          "\n ipivot = %d, jj = %d, first = %d, last = %d, kstart = %d",
                       ipivot, jj, first, last, kstart) ;
#endif
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = 0 ; ii < first ; ii++ ) {
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                        fprintf(stdout, ", keep") ;
#endif
                        count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
               first = last + 1 ;
            }
            for ( jj = nD ; jj < ncol ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               for ( ii = 0 ; ii < nD ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         }
      }
   } else {
/*
      --------------------------------------
      chevron is nonsymmetric, D is diagonal
      --------------------------------------
*/
      count  = 0 ;
      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < ncol ; jj++ ) {
            kk     = kstart ;
            kinc   = stride ;
            iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
            for ( ii = 0 ; ii <= iilast ; ii++ ) {
               absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            kstart++ ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) {
         for ( jj = 0 ; jj < ncol ; jj++ ) {
            kk     = kstart ;
            kinc   = stride ;
            iilast = (jj < nD) ? (jj - 1) : (nD - 1) ;
            for ( ii = 0 ; ii <= iilast ; ii++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            kstart++ ;
         }
      }
   }
   } break ;
case CHV_STRICT_LOWER_11 : {
   int   ii, jj, kinc, kk, kstart, nrow, stride ;
/*
      ----------------------------------------
      count the number of big entries in the 
      strict lower triangle of the (1,1) block
      ----------------------------------------
*/
   nrow   = nD + nL ;
   count  = 0 ;
   kstart = nD + nL - 1 ;
   stride = 2*nD + nL + nU - 1 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( ii = 0 ; ii < nD ; ii++ ) {
         kk   = kstart ;
         kinc = stride ;
         for ( jj = 0 ; jj < ii ; jj++ ) {
            absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
            if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
               count++ ;
            }
            kk   += kinc ;
            kinc -=   2  ;
         }
         kstart-- ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( ii = 0 ; ii < nD ; ii++ ) {
         kk   = kstart ;
         kinc = stride ;
         for ( jj = 0 ; jj < ii ; jj++ ) {
            absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
            if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
               count++ ;
            }
            kk   += kinc ;
            kinc -=   2  ;
         }
         kstart-- ;
      }
   }
   } break ;
case CHV_LOWER_21        : {
   int   ii, jj, kinc, kk, kstart, nrow, stride ;
/*
      --------------------------------------------------
      count the number of big entries in the (2,1) block
      --------------------------------------------------
*/
   nrow   = nD + nL ;
   count  = 0 ;
   kstart = nL - 1 ;
   stride = 2*nD + nL + nU - 1 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( ii = nD ; ii < nrow ; ii++ ) {
         kk   = kstart ;
         kinc = stride ;
         for ( jj = 0 ; jj < nD ; jj++ ) {
            absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
            if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
               count++ ;
            }
            kk   += kinc ;
            kinc -=   2  ;
         }
         kstart-- ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( ii = nD ; ii < nrow ; ii++ ) {
         kk   = kstart ;
         kinc = stride ;
         for ( jj = 0 ; jj < nD ; jj++ ) {
            absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
            if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
               count++ ;
            }
            kk   += kinc ;
            kinc -=   2  ;
         }
         kstart-- ;
      }
   }
   } break ;
case CHV_STRICT_UPPER_11 : {
/*
   ---------------------------------------------------------------
   count the number of big entries in the strict upper (1,1) block
   ---------------------------------------------------------------
*/
   int   first, ii, ipivot, jj, kinc, kk, kstart, 
         last, ncol, stride ;

   ncol = nD + nU ;
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n symmetric chevron, npivot = %d, pivotsizes = %p", 
              npivot, pivotsizes) ;
#endif
/*
      --------------------
      chevron is symmetric
      --------------------
*/
      count  = 0 ;
      kstart = 0 ;
      stride = nD + nU - 1 ;
      if ( pivotsizes == NULL ) {
/*
         -------------
         D is diagonal
         -------------
*/
         if ( CHV_IS_REAL(chv) ) {
            for ( jj = 0 ; jj < nD ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               for ( ii = 0 ; ii < jj ; ii++ ) {
                  absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( jj = 0 ; jj < nD ; jj++ ) {
               kk   = kstart ;
               kinc = stride ;
               for ( ii = 0 ; ii < jj ; ii++ ) {
                  absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
                  if ( absval >= droptol ) {
#if MYDEBUG > 0
                     fprintf(stdout, ", keep") ;
#endif
                     count++ ;
                  }
                  kk += kinc ;
                  kinc-- ;
               }
               kstart++ ;
            }
         }
      } else {
/*
         ---------------------------------
         D can have 1 x 1 and 2 x 2 pivots
         ---------------------------------
*/
         if ( CHV_IS_REAL(chv) ) {
            for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
               last = first + pivotsizes[ipivot] - 1 ;
               for ( jj = first ; jj <= last ; jj++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, 
       "\n %% ipivot = %d, jj = %d, first = %d, last = %d, kstart = %d",
                       ipivot, jj, first, last, kstart) ;
#endif
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = 0 ; ii < first ; ii++ ) {
                     absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                        fprintf(stdout, ", keep") ;
#endif
                        count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
               first = last + 1 ;
            }
         } else if ( CHV_IS_COMPLEX(chv) ) {
            for ( ipivot = first = 0 ; ipivot < npivot ; ipivot++ ) {
               last = first + pivotsizes[ipivot] - 1 ;
               for ( jj = first ; jj <= last ; jj++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, 
       "\n %% ipivot = %d, jj = %d, first = %d, last = %d, kstart = %d",
                       ipivot, jj, first, last, kstart) ;
#endif
                  kk   = kstart ;
                  kinc = stride ;
                  for ( ii = 0 ; ii < first ; ii++ ) {
                     absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                          ii, jj, absval, kk) ;
#endif
                     if ( absval >= droptol ) {
#if MYDEBUG > 0
                        fprintf(stdout, ", keep") ;
#endif
                        count++ ;
                     }
                     kk += kinc ;
                     kinc-- ;
                  }
                  kstart++ ;
               }
               first = last + 1 ;
            }
         }
      }
   } else {
/*
      --------------------------------------
      chevron is nonsymmetric, D is diagonal
      --------------------------------------
*/
      count  = 0 ;
      kstart = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk     = kstart ;
            kinc   = stride ;
            for ( ii = 0 ; ii < jj ; ii++ ) {
               absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            kstart++ ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) { 
         for ( jj = 0 ; jj < nD ; jj++ ) {
            kk     = kstart ;
            kinc   = stride ;
            for ( ii = 0 ; ii < jj ; ii++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            kstart++ ;
         }
      }
   }
   } break ;
case CHV_UPPER_12        : {
   int   ii, jj, kinc, kk, kstart, ncol, stride ;
/*
   ----------------------------------------
   count the big entries in the (1,2) block
   ----------------------------------------
*/
   ncol = nD + nU ;
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n symmetric chevron, npivot = %d, pivotsizes = %p", 
              npivot, pivotsizes) ;
#endif
/*
      --------------------
      chevron is symmetric
      --------------------
*/
      count  = 0 ;
      kstart = nD ;
      stride = nD + nU - 1 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = nD ; jj < ncol ; jj++ ) {
            kk   = kstart ;
            kinc = stride ;
            for ( ii = 0 ; ii < nD ; ii++ ) {
               absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc-- ;
            }
            kstart++ ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) { 
         for ( jj = nD ; jj < ncol ; jj++ ) {
            kk   = kstart ;
            kinc = stride ;
            for ( ii = 0 ; ii < nD ; ii++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                    ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
               fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc-- ;
            }
            kstart++ ;
         }
      }
   } else {
/*
      -----------------------
      chevron is nonsymmetric
      -----------------------
*/
      count  = 0 ;
      kstart = 2*nD + nL - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      if ( CHV_IS_REAL(chv) ) {
         for ( jj = nD ; jj < ncol ; jj++ ) {
            kk     = kstart ;
            kinc   = stride ;
            for ( ii = 0 ; ii < nD ; ii++ ) {
               absval = fabs(entries[kk]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            kstart++ ;
         }
      } else if ( CHV_IS_COMPLEX(chv) ) { 
         for ( jj = nD ; jj < ncol ; jj++ ) {
            kk     = kstart ;
            kinc   = stride ;
            for ( ii = 0 ; ii < nD ; ii++ ) {
               absval = Zabs(entries[2*kk], entries[2*kk+1]) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n %% |A(%d,%d)| = %20.12e, kk = %d",
                       ii, jj, absval, kk) ;
#endif
               if ( absval >= droptol ) {
#if MYDEBUG > 0
                  fprintf(stdout, ", keep") ;
#endif
                  count++ ;
               }
               kk += kinc ;
               kinc -= 2 ;
            }
            kstart++ ;
         }
      }
   }
   } break ;
default :
   break ;
}

return(count) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- copy the trailing chevron that starts at offset
 
   created -- 97may16, cca
   ----------------------------------------------------------
*/
void
Chv_copyTrailingPortion (
   Chv   *chvI,
   Chv   *chvJ,
   int    offset
) {
int   first, ncolI, ncolJ, nDJ, nentToCopy, nLJ, nrowI, nrowJ, nUJ ;
int   *colindI, *colindJ, *rowindI, *rowindJ ;
/*
   --------------
   check the data
   --------------
*/
if ( chvI == NULL || chvJ == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_copyTrailingPortion(%p,%p,%d)"
           "\n bad input\n", chvI, chvJ, offset) ;
   exit(-1) ;
}
Chv_dimensions(chvJ, &nDJ, &nLJ, &nUJ) ;
if ( offset < 0 || offset >= nDJ ) {
   fprintf(stderr, 
           "\n fatal error in Chv_copyTrailingPortion(%p,%p,%d)"
           "\n nDJ = %d, offset = %d\n", 
           chvI, chvJ, offset, nDJ, offset) ;
   exit(-1) ;
}
/*
   ----------------------------------------------
   initialize the chvI object and find the number
   of and first location of the entries to copy
   ----------------------------------------------
*/
Chv_columnIndices(chvJ, &ncolJ, &colindJ) ;
if ( CHV_IS_SYMMETRIC(chvJ) || CHV_IS_HERMITIAN(chvJ) ) {
   Chv_init(chvI, chvJ->id, nDJ - offset, 0, nUJ, 
            chvJ->type, chvJ->symflag) ;
   Chv_columnIndices(chvI, &ncolI, &colindI) ;
   IVcopy(nDJ + nUJ - offset, colindI, colindJ + offset) ;
   first = offset*(nDJ + nUJ) - (offset*(offset-1))/2 ;
   nentToCopy = (nDJ*(nDJ+1))/2 + nDJ*nUJ - first ;
} else {
   Chv_rowIndices(chvJ, &nrowJ, &rowindJ) ;
   Chv_init(chvI, chvJ->id, nDJ - offset, nLJ, nUJ, 
            chvJ->type, chvJ->symflag) ;
   Chv_columnIndices(chvI, &ncolI, &colindI) ;
   IVcopy(nDJ + nUJ - offset, colindI, colindJ + offset) ;
   Chv_rowIndices(chvI, &nrowI, &rowindI) ;
   IVcopy(nDJ + nLJ - offset, rowindI, rowindJ + offset) ;
   first = offset*(2*nDJ + nLJ + nUJ - offset) ;
   nentToCopy = nDJ*(nDJ + nLJ + nUJ) - first ;
}
if ( CHV_IS_REAL(chvJ) ) {
   DVcopy(nentToCopy, Chv_entries(chvI), Chv_entries(chvJ) + first) ;
} else if ( CHV_IS_COMPLEX(chvJ) ) { 
   DVcopy(2*nentToCopy, Chv_entries(chvI), Chv_entries(chvJ) + 2*first);
}
return ; }

/*--------------------------------------------------------------------*/
