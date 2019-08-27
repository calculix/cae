/*  MMM.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- setup the distributed matrix-matrix multiply
              information object

   created -- 98aug21, cca
   -------------------------------------------------------
*/
MatMulInfo *
MatMul_MPI_setup (
   InpMtx     *A,
   int        symflag,
   int        opflag,
   IV         *XownersIV,
   IV         *YownersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        tag,
   MPI_Comm   comm
) {
MatMulInfo   *info ;
int          iproc, myid, nproc, size ;
int          *list, *map, *Xmap, *Ymap ;
IV           *mapIV, *XmapIV, *XownedIV, *XsupIV,
             *YmapIV, *YownedIV, *YsupIV ;
IVL          *XrecvIVL, *XsendIVL, *YrecvIVL, *YsendIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || XownersIV == NULL || YownersIV == NULL
   || stats == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in MatMul_MPI_setup()"
           "\n bad input\n") ;
   exit(-1) ;
}
switch ( symflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n fatal error in MatMul_MPI_setup()"
           "\n bad symflag = %d\n", symflag) ;
   exit(-1) ;
   break ;
}
switch ( opflag ) {
case MMM_WITH_A :
case MMM_WITH_AT :
case MMM_WITH_AH :
   break ;
default :
   fprintf(stderr, "\n fatal error in MatMul_MPI_setup()"
           "\n bad opflag = %d", opflag) ;
   exit(-1) ;
   break ;
}
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
/*
   ------------------------------
   allocate the MatMulInfo object
   ------------------------------
*/
ALLOCATE(info, MatMulInfo, 1) ;
info->symflag  = symflag ;
info->opflag   = opflag  ;
info->YsupIV   = info->XsupIV   = NULL ;
info->YmapIV   = info->XmapIV   = NULL ;
info->XsendIVL = info->XrecvIVL = NULL ;
info->YsendIVL = info->YrecvIVL = NULL ;
info->Xsupp    = info->Ysupp    = NULL ;
/*
   ----------------------------------------------------------------
   get the list of rows of X and Y that are owned by this processor
   ----------------------------------------------------------------
*/
XownedIV = info->XownedIV = IV_targetEntries(XownersIV, myid) ;
YownedIV = info->YownedIV = IV_targetEntries(YownersIV, myid) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n owned rows of X") ;
   IV_writeForHumanEye(XownedIV, msgFile) ;
   fprintf(msgFile, "\n\n owned rows of Y") ;
   IV_writeForHumanEye(YownedIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   find the support of the owned matrix
   ------------------------------------
*/
switch ( symflag ) {
case SPOOLES_SYMMETRIC :
   XsupIV = info->XsupIV = IV_new() ;
   YsupIV = info->YsupIV = XsupIV ;
   InpMtx_supportSym(A, YsupIV) ;
   break ;
case SPOOLES_HERMITIAN :
   XsupIV = info->XsupIV = IV_new() ;
   YsupIV = info->YsupIV = XsupIV ;
   InpMtx_supportHerm(A, YsupIV) ;
   break ;
case SPOOLES_NONSYMMETRIC :
   XsupIV = info->XsupIV = IV_new() ;
   YsupIV = info->YsupIV = IV_new() ;
   switch ( opflag ) {
   case MMM_WITH_A :
      InpMtx_supportNonsym(A, YsupIV, XsupIV) ;
      break ;
   case MMM_WITH_AT :
      InpMtx_supportNonsymT(A, YsupIV, XsupIV) ;
      break ;
   case MMM_WITH_AH :
      InpMtx_supportNonsymH(A, YsupIV, XsupIV) ;
      break ;
   }
   break ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n row support") ;
   IV_writeForHumanEye(YsupIV, msgFile) ;
   fprintf(msgFile, "\n\n column support") ;
   IV_writeForHumanEye(XsupIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   get the maps from global to local coordinates
   ---------------------------------------------
*/
if ( IV_size(XsupIV) > 0 ) {
   XmapIV = info->XmapIV = IV_inverseMap(XsupIV) ;
} else {
   XmapIV = info->XmapIV = IV_new() ;
}
if ( symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN ) {
   YmapIV = info->YmapIV = XmapIV ;
} else {
   if ( IV_size(YsupIV) > 0 ) {
      YmapIV = info->YmapIV = IV_inverseMap(YsupIV) ;
   } else {
      YmapIV = info->YmapIV = IV_new() ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n YsupIV") ;
   IV_writeForHumanEye(YsupIV, msgFile) ;
   fprintf(msgFile, "\n\n YmapIV") ;
   IV_writeForHumanEye(YmapIV, msgFile) ;
   fprintf(msgFile, "\n\n XsupIV") ;
   IV_writeForHumanEye(XsupIV, msgFile) ;
   fprintf(msgFile, "\n\n XmapIV") ;
   IV_writeForHumanEye(XmapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------------------
   create the IVL objects that specify the communication patterns
   --------------------------------------------------------------
*/
XsendIVL = info->XsendIVL = IVL_new() ;
XrecvIVL = info->XrecvIVL = IVL_new() ;
makeSendRecvIVLs(XsupIV, XownersIV, XsendIVL, XrecvIVL,
                 stats, msglvl, msgFile, tag, comm) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n X row send IVL, global") ;
   IVL_writeForHumanEye(XsendIVL, msgFile) ;
   fprintf(msgFile, "\n\n X row receive IVL, global") ;
   IVL_writeForHumanEye(XrecvIVL, msgFile) ;
   fflush(msgFile) ;
}
YsendIVL = info->YsendIVL = IVL_new() ;
YrecvIVL = info->YrecvIVL = IVL_new() ;
makeSendRecvIVLs(YsupIV, YownersIV, YrecvIVL, YsendIVL,
                 stats, msglvl, msgFile, tag, comm) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Y row send IVL") ;
   IVL_writeForHumanEye(YsendIVL, msgFile) ;
   fprintf(msgFile, "\n\n Y row receive IVL") ;
   IVL_writeForHumanEye(YrecvIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------
   make XsendIVL local w.r.t. local X
   make XrecvIVL local w.r.t. supported X
   make YsendIVL local w.r.t. supported Y
   make YrecvIVL local w.r.t. local Y
   --------------------------------------
*/
if ( IV_size(XownedIV) > 0 ) {
   Xmap = IV_entries(XmapIV) ;
   mapIV = IV_inverseMap(XownedIV) ;
   map = IV_entries(mapIV) ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      IVL_listAndSize(XsendIVL, iproc, &size, &list) ;
      IVgather(size, list, map, list) ;
      IVL_listAndSize(XrecvIVL, iproc, &size, &list) ;
      IVgather(size, list, Xmap, list) ;
   }
   IV_free(mapIV) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n X row send IVL, local") ;
   IVL_writeForHumanEye(XsendIVL, msgFile) ;
   fprintf(msgFile, "\n\n X row receive IVL, local") ;
   IVL_writeForHumanEye(XrecvIVL, msgFile) ;
   fflush(msgFile) ;
}
if ( IV_size(YownedIV) > 0 ) {
   Ymap  = IV_entries(YmapIV) ;
   mapIV = IV_inverseMap(YownedIV) ;
   map   = IV_entries(mapIV) ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      IVL_listAndSize(YsendIVL, iproc, &size, &list) ;
      IVgather(size, list, Ymap, list) ;
      IVL_listAndSize(YrecvIVL, iproc, &size, &list) ;
      IVgather(size, list, map, list) ;
   }
   IV_free(mapIV) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Y row send IVL, local") ;
   IVL_writeForHumanEye(YsendIVL, msgFile) ;
   fprintf(msgFile, "\n\n Y row receive IVL, local") ;
   IVL_writeForHumanEye(YrecvIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   set up the Xsupp and Ysupp objects
   ----------------------------------
*/
info->Xsupp = DenseMtx_new() ;
info->Ysupp = DenseMtx_new() ;

return(info) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   set the indices of A to be local with respect to its support

   created -- 98aug21, cca
   ------------------------------------------------------------
*/
void
MatMul_setLocalIndices (
   MatMulInfo   *info,
   InpMtx       *A
) {
if ( info == NULL || A == NULL ) {
   fprintf(stderr, "\n fatal error in MatMul_setLocalIndices()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( A->nent > 0 ) {
/*
   -------------------------------------------
   map the input matrix into local coordinates
   -------------------------------------------
*/
   switch ( info->symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
      InpMtx_mapEntries(A, info->YmapIV, info->XmapIV) ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      switch ( info->opflag ) {
      case MMM_WITH_A :
         InpMtx_mapEntries(A, info->YmapIV, info->XmapIV) ;
         break ;
      case MMM_WITH_AT :
      case MMM_WITH_AH :
         InpMtx_mapEntries(A, info->XmapIV, info->YmapIV) ;
         break ;
      }
      break ;
   default :
      fprintf(stderr, "\n fatal error in MatMul_setLocalIndices()"
              "\n bad symflag = %d\n", info->symflag) ;
      exit(-1) ;
      break ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   set the indices of A to be global with respect to its support

   created -- 98aug21, cca
   -------------------------------------------------------------
*/
void
MatMul_setGlobalIndices (
   MatMulInfo   *info,
   InpMtx       *A
) {
if ( info == NULL || A == NULL ) {
   fprintf(stderr, "\n fatal error in MatMul_setGlobalIndices()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( A->nent > 0 ) {
/*
   --------------------------------------------
   map the input matrix into global coordinates
   --------------------------------------------
*/
   switch ( info->symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
      InpMtx_mapEntries(A, info->YsupIV, info->XsupIV) ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      switch ( info->opflag ) {
      case MMM_WITH_A :
         InpMtx_mapEntries(A, info->YsupIV, info->XsupIV) ;
         break ;
      case MMM_WITH_AT :
      case MMM_WITH_AH :
         InpMtx_mapEntries(A, info->XsupIV, info->YsupIV) ;
         break ;
      }
      break ;
   default :
      fprintf(stderr, "\n fatal error in MatMul_setGlobalIndices()"
              "\n bad symflag = %d\n", info->symflag) ;
      exit(-1) ;
      break ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- compute the distributed matrix-matrix multiply
              Y := Y - alpha * A * X
      where A, Xloc and Yloc are distributed

   created -- 98aug21, cca
   ---------------------------------------------------------
*/
void
MatMul_MPI_mmm (
   MatMulInfo   *info,
   DenseMtx     *Yloc,
   double       alpha[],
   InpMtx       *A,
   DenseMtx     *Xloc,
   int          stats[],
   int          msglvl,
   FILE         *msgFile,
   int          tag,
   MPI_Comm     comm
) {
int   ncolXloc, ncolYloc, nrowXloc, nrowYloc, nrowXsupp, nrowYsupp ;
/*
   ---------------
   check the input
   ---------------
*/
if ( info == NULL || Yloc == NULL || alpha == NULL 
    || A == NULL || Xloc == NULL || stats == NULL 
    || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in MatMul_MPI_mmm()"
           "\n bad input\n") ;
   exit(-1) ;
}
/*
   -------------------------------------------------
   setup Xsupp and gather rows from distributed Xloc
   -------------------------------------------------
*/
DenseMtx_dimensions(Xloc, &nrowXloc, &ncolXloc) ;
nrowXsupp = IV_size(info->XsupIV) ;
DenseMtx_init(info->Xsupp, Xloc->type, 0, 0, 
              nrowXsupp, ncolXloc, 1, nrowXsupp) ;
DenseMtx_zero(info->Xsupp) ;
DenseMtx_MPI_gatherRows(info->Xsupp, Xloc, info->XsendIVL, 
                    info->XrecvIVL, stats, msglvl, msgFile, tag, comm) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n supported matrix Xsupp") ;
   DenseMtx_writeForHumanEye(info->Xsupp, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------
   setup Ysupp
   -----------
*/
DenseMtx_dimensions(Yloc, &nrowYloc, &ncolYloc) ;
nrowYsupp = IV_size(info->YsupIV) ;
DenseMtx_init(info->Ysupp, Yloc->type, 0, 0, 
              nrowYsupp, ncolYloc, 1, nrowYsupp) ;
DenseMtx_zero(info->Ysupp) ;
/*
   ----------------------------------
   compute the matrix-matrix multiply
   ----------------------------------
*/
if ( A->nent > 0 ) {
   switch ( info->symflag ) {
   case SPOOLES_SYMMETRIC :
      InpMtx_sym_mmm(A, info->Ysupp, alpha, info->Xsupp) ;
      break ;
   case SPOOLES_HERMITIAN :
      InpMtx_herm_mmm(A, info->Ysupp, alpha, info->Xsupp) ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      switch ( info->opflag ) {
      case MMM_WITH_A :
         InpMtx_nonsym_mmm(A, info->Ysupp, alpha, info->Xsupp) ;
         break ;
      case MMM_WITH_AT :
         InpMtx_nonsym_mmm_T(A, info->Ysupp, alpha, info->Xsupp) ;
         break ;
      case MMM_WITH_AH :
         InpMtx_nonsym_mmm_H(A, info->Ysupp, alpha, info->Xsupp) ;
         break ;
      }
      break ;
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after mmm, local supported matrix Ysupp") ;
   DenseMtx_writeForHumanEye(info->Ysupp, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   assemble the owned rows of Yloc
   -------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n before scatter/add, local matrix Y") ;
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
   fflush(msgFile) ;
}
DenseMtx_MPI_scatterAddRows(Yloc, info->Ysupp, info->YsendIVL, 
                   info->YrecvIVL, stats, msglvl, msgFile, tag, comm) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after mmm, local matrix Y") ;
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   free the MatMulInfo object and its data structures 

   created -- 98aug21, cca
   --------------------------------------------------
*/
void
MatMul_cleanup (
   MatMulInfo   *info
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in MatMul_cleanup()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( info->XownedIV != NULL ) {
   IV_free(info->XownedIV) ; info->XownedIV = NULL ;
}
if ( info->XsupIV != NULL ) {
   IV_free(info->XsupIV) ; info->XsupIV = NULL ;
}
if ( info->XmapIV != NULL ) {
   IV_free(info->XmapIV) ; info->XmapIV = NULL ;
}
if ( info->XrecvIVL != NULL ) {
   IVL_free(info->XrecvIVL) ; info->XrecvIVL = NULL ;
}
if ( info->XsendIVL != NULL ) {
   IVL_free(info->XsendIVL) ; info->XsendIVL = NULL ;
}
if ( info->Xsupp != NULL ) {
   DenseMtx_free(info->Xsupp) ; info->Xsupp = NULL ;
}
if ( info->YownedIV != NULL ) {
   IV_free(info->YownedIV) ; info->YownedIV = NULL ;
}
if ( info->symflag == SPOOLES_NONSYMMETRIC ) {
   if ( info->YsupIV != NULL ) {
      IV_free(info->YsupIV) ; info->YsupIV = NULL ;
   }
   if ( info->YmapIV != NULL ) {
      IV_free(info->YmapIV) ; info->YmapIV = NULL ;
   }
}
if ( info->YrecvIVL != NULL ) {
   IVL_free(info->YrecvIVL) ; info->YrecvIVL = NULL ;
}
if ( info->YsendIVL != NULL ) {
   IVL_free(info->YsendIVL) ; info->YsendIVL = NULL ;
}
if ( info->Ysupp != NULL ) {
   DenseMtx_free(info->Ysupp) ; info->Ysupp = NULL ;
}
info->symflag = -1 ;
info->opflag  = -1 ;
FREE(info) ;

return ; }

/*--------------------------------------------------------------------*/
