/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef SPOOLES
#include <misc.h>
#include <FrontMtx.h>
#include <SymbFac.h>
#endif

#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void cascade(ITG *ipompc, double **coefmpcp, ITG **nodempcp, ITG *nmpc,
   ITG *mpcfree, ITG *nodeboun, ITG *ndirboun, ITG*nboun, ITG*ikmpc,
   ITG *ilmpc, ITG *ikboun, ITG *ilboun, ITG *mpcend,
   char *labmpc, ITG *nk, ITG *memmpc_, ITG *icascade, ITG *maxlenmpc,
   ITG *callfrommain, ITG *iperturb, ITG *ithermal){

 /*   detects cascaded mpc's and decascades them; checks multiple
     occurrence of the same dependent DOF's in different mpc/spc's

     data structure of ipompc,coefmpc,nodempc:
       for each mpc, e.g. i, 
         -the nodes are stored in nodempc(1,ipompc(i)),
          nodempc(1,nodempc(3,ipompc(i))),
          nodempc(1,nodempc(3,nodempc(3,ipompc(i))))... till
          nodempc(3,nodempc(3,nodempc(3,.......))))))=0;
         -the corresponding directions in nodempc(2,ipompc(i)),
          nodempc(2,nodempc(3,ipompc(i))),.....
         -the corresponding coefficient in coefmpc(ipompc(i)),
          coefmpc(nodempc(3,ipompc(i))),.....
       the mpc is written as a(1)u(i1,j1)+a(2)u(i2,j2)+...
       +....a(k)u(ik,jk)=0, the first term is the dependent term,
       the others are independent, at least after execution of the
       present routine. The mpc's must be homogeneous, otherwise a
       error message is generated and the program stops. */

    ITG i,j,index,id,idof,nterm,idepend,*nodempc=NULL,
	ispooles,iexpand,ichange,indexold,ifluidmpc,
        mpc,indexnew,index1,index2,index1old,index2old,*jmpc=NULL,nl;

    double coef,*coefmpc=NULL;

#ifdef SPOOLES

    ITG irow,icolumn,node,idir,irownl,icolnl,*ipointer=NULL,*icoef=NULL,
	ifree,*indepdof=NULL,nindep;

    double *xcoef=NULL,b;

    DenseMtx        *mtxB, *mtxX ;
    Chv             *rootchv ;
    ChvManager      *chvmanager  ;
    SubMtxManager   *mtxmanager  ;
    FrontMtx        *frontmtx ;
    InpMtx          *mtxA ;
    double          tau = 100.;
    double          cpus[10] ;
    ETree           *frontETree ;
    FILE            *msgFile ;
    Graph           *graph ;
    ITG             jrhs, msglvl=0, nedges,error,
                    nent, neqns, nrhs, pivotingflag=1, seed=389,
                    symmetryflag=2, type=1,maxdomainsize,maxzeros,maxsize;
    ITG             *oldToNew ;
    ITG             stats[20] ;
    IV              *newToOldIV, *oldToNewIV ;
    IVL             *adjIVL, *symbfacIVL ;
#endif

    nodempc=*nodempcp;
    coefmpc=*coefmpcp;
    
    /*      for(i=0;i<*nmpc;i++){
	j=i+1;
	FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
	}*/

    NNEW(jmpc,ITG,*nmpc);
    idepend=0;

/*        check whether a node is used as a dependent node in a MPC
	  and in a SPC */

    for(i=0;i<*nmpc;i++){
	if(*nboun>0){
	    FORTRAN(nident,(ikboun,&ikmpc[i],nboun,&id));}
	else{id=0;}
	if(id>0){
	    if(ikboun[id-1]==ikmpc[i]){
		if(strcmp1(&labmpc[20*i],"FLUID")!=0){
		    printf("*ERROR in cascade: the DOF corresponding to \n node %" ITGFORMAT " in direction %" ITGFORMAT " is detected on the \n dependent side of a MPC and a SPC\n\n",
			   (ikmpc[i])/8+1,ikmpc[i]-8*((ikmpc[i])/8));
		}else{
		    printf("*ERROR in cascade: the DOF corresponding to \n face %" ITGFORMAT " in direction %" ITGFORMAT " is detected on the \n dependent side of a MPC and a SPC\n\n",
			   (-ikmpc[i])/8+1,-ikmpc[i]-8*((-ikmpc[i])/8));
		}
		FORTRAN(stop,());
	    }
	}
    }

/*     check whether there are user mpc's: in user MPC's the
       dependent DOF can change, however, the number of terms
       cannot change   */

    for(i=0;i<*nmpc;i++){

        /* linear mpc */

        /* because of the next line the size of field labmpc
           has to be defined as 20*nmpc+1: without "+1" an
           undefined field is accessed */

	if((strcmp1(&labmpc[20*i],"                    ")==0) ||
	   (strcmp1(&labmpc[20*i],"CYCLIC")==0) ||
	   (strcmp1(&labmpc[20*i],"SUBCYCLIC")==0)||
	   (strcmp1(&labmpc[20*i],"PRETENSION")==0)||
	   (strcmp1(&labmpc[20*i],"THERMALPRET")==0)||
//	   (strcmp1(&labmpc[20*i],"CONTACT")==0)||
	   (strcmp1(&labmpc[20*i],"FLUID")==0)||
           (*iperturb<2)) jmpc[i]=0;

        /* nonlinear mpc */

	else if((strcmp1(&labmpc[20*i],"RIGID")==0) ||
	   (strcmp1(&labmpc[20*i],"KNOT")==0) ||
	   (strcmp1(&labmpc[20*i],"PLANE")==0) ||
	   (strcmp1(&labmpc[20*i],"BEAM")==0) ||
	   (strcmp1(&labmpc[20*i],"STRAIGHT")==0)) jmpc[i]=1;

        /* user mpc */

	else{
	    jmpc[i]=1;
	    if(*icascade==0) *icascade=1;
	}
    }
    
/*     decascading */

    ispooles=0;

    /* decascading using simple substitution */

    do{
        ichange=0;
        for(i=0;i<*nmpc;i++){

	    if(strcmp1(&labmpc[20*i],"FLUID")!=0){
		ifluidmpc=0;
	    }else{
		ifluidmpc=1;
	    }

	    if(jmpc[i]==1) nl=1;
	    else nl=0;
	    iexpand=0;
	    index=nodempc[3*ipompc[i]-1];
	    if(index==0) continue;
	    do{
		if(ifluidmpc==0){
		    /* MPC on node */
		    idof=(nodempc[3*index-3]-1)*8+nodempc[3*index-2];
		}else{
		    if(nodempc[3*index-3]>0){
			/* MPC on face */
			idof=-((nodempc[3*index-3]-1)*8+nodempc[3*index-2]);
		    }else{
			/* MPC on node 
                           SPC number: -nodempc[3*index-3]
                           node: nodeboun[-nodempc[3*index-3]-1] */

			idof=(nodeboun[-nodempc[3*index-3]-1]-1)*8+nodempc[3*index-2];
		    }
		}
		
		FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		if((id>0)&&(ikmpc[id-1]==idof)){
		    
		    /* a term on the independent side of the MPC is
		       detected as dependent node in another MPC */

		    indexold=nodempc[3*index-1];
		    coef=coefmpc[index-1];
		    mpc=ilmpc[id-1];

                    /* no expansion if there is a dependence of a
                       nonlinear MPC on another linear or nonlinear MPC
                       and the call is from main */ 

		    if((jmpc[mpc-1]==1)||(nl==1)){
			*icascade=2;
			if(idepend==0){
			    printf("*INFO in cascade: linear MPCs and\n");
			    printf("       nonlinear MPCs depend on each other\n");
			    printf("       common node: %" ITGFORMAT " in direction %" ITGFORMAT "\n\n",nodempc[3*index-3],nodempc[3*index-2]);
			    idepend=1;}
			if(*callfrommain==1){
			    index=nodempc[3*index-1];
			    if(index!=0) continue;
			    else break;}
		    }

/*		    printf("*INFO in cascade: DOF %" ITGFORMAT " of node %" ITGFORMAT " is expanded\n",
		    nodempc[3*index-2],nodempc[3*index-3]);*/

		    /* collecting terms corresponding to the same DOF */
		    
		    index1=ipompc[i];
		    do{
			index2old=index1;
			index2=nodempc[3*index1-1];
			if(index2==0) break;
			do{
			    if((nodempc[3*index1-3]==nodempc[3*index2-3])&&
			       (nodempc[3*index1-2]==nodempc[3*index2-2])){
				coefmpc[index1-1]+=coefmpc[index2-1];
				nodempc[3*index2old-1]=nodempc[3*index2-1];
				nodempc[3*index2-1]=*mpcfree;
				*mpcfree=index2;
				index2=nodempc[3*index2old-1];
				if(index2==0) break;
			    }
			    else{
				index2old=index2;
				index2=nodempc[3*index2-1];
				if(index2==0) break;
			    }
			}while(1);
			index1=nodempc[3*index1-1];
			if(index1==0) break;
		    }while(1);
		    
		    /* check for zero coefficients on the dependent side */
		    
		    	    index1=ipompc[i];
			    if(fabs(coefmpc[index1-1])<1.e-10){
				printf("*ERROR in cascade: zero coefficient on the\n");
				printf("       dependent side of an equation\n");
				printf("       dependent node: %" ITGFORMAT "",nodempc[3*index1-3]);
				printf("       direction: %" ITGFORMAT "\n\n",nodempc[3*index1-2]);
				FORTRAN(stop,());
			    }
		    
		    ichange=1;iexpand=1;
		    if((strcmp1(&labmpc[20*i],"                    ")==0)&&
		       (strcmp1(&labmpc[20*(mpc-1)],"CYCLIC")==0))
			strcpy1(&labmpc[20*i],"SUBCYCLIC",9);
		    indexnew=ipompc[mpc-1];
		    coef=-coef/coefmpc[indexnew-1];
		    indexnew=nodempc[3*indexnew-1];
		    if(indexnew!=0){
			do{
			    coefmpc[index-1]=coef*coefmpc[indexnew-1];
			    nodempc[3*index-3]=nodempc[3*indexnew-3];
			    nodempc[3*index-2]=nodempc[3*indexnew-2];
			    indexnew=nodempc[3*indexnew-1];
			    if(indexnew!=0){
				nodempc[3*index-1]=*mpcfree;
				index=*mpcfree;
				*mpcfree=nodempc[3**mpcfree-1];
				if(*mpcfree==0){
				    *mpcfree=*memmpc_+1;
				    nodempc[3*index-1]=*mpcfree;
				    *memmpc_=(ITG)(1.1**memmpc_);
				    printf("*INFO in cascade: reallocating nodempc; new size = %" ITGFORMAT "\n\n",*memmpc_);
				    RENEW(nodempc,ITG,3**memmpc_);
				    RENEW(coefmpc,double,*memmpc_);
				    for(j=*mpcfree;j<*memmpc_;j++){
					nodempc[3*j-1]=j+1;
				    }
				    nodempc[3**memmpc_-1]=0;
				}
				continue;
			    }
			    else{
				nodempc[3*index-1]=indexold;
				break;
			    }
			}while(1);
		    }else{
			coefmpc[index-1]=0.;
		    }
		    break;
		}
		else{
		    index=nodempc[3*index-1];
		    if(index!=0) continue;
		    else break;
		}
	    }while(1);
	    if(iexpand==0) continue;
	    
	    /* one term of the mpc was expanded 
	       collecting terms corresponding to the same DOF */
	    
	    index1=ipompc[i];
	    do{
		index2old=index1;
		index2=nodempc[3*index1-1];
		if(index2==0) break;
		do{
		    if((nodempc[3*index1-3]==nodempc[3*index2-3])&&
		       (nodempc[3*index1-2]==nodempc[3*index2-2])){
			coefmpc[index1-1]+=coefmpc[index2-1];
			nodempc[3*index2old-1]=nodempc[3*index2-1];
			nodempc[3*index2-1]=*mpcfree;
			*mpcfree=index2;
			index2=nodempc[3*index2old-1];
			if(index2==0) break;
		    }
		    else{
			index2old=index2;
			index2=nodempc[3*index2-1];
			if(index2==0) break;
		    }
		}while(1);
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	    }while(1);
	    
	    /* check for zero coefficients on the dependent and
	       independent side */
	    
	    index1=ipompc[i];
	    index1old=0;
	    do {
//		printf("index1=%d\n",index1);
//		printf("coefmpc=%e\n",coefmpc[index1-1]);

		if(fabs(coefmpc[index1-1])<1.e-10){
		    if(index1old==0){
			printf("*ERROR in cascade: zero coefficient on the\n");
			printf("       dependent side of an equation\n");
			printf("       dependent node: %" ITGFORMAT "",nodempc[3*index1-3]);
			printf("       direction: %" ITGFORMAT "\n\n",nodempc[3*index1-2]);
			FORTRAN(stop,());
		    }
		    else{
			nodempc[3*index1old-1]=nodempc[3*index1-1];
			nodempc[3*index1-1]=*mpcfree;
			*mpcfree=index1;
			index1=nodempc[3*index1old-1];
		    }
		}
		else{
		    index1old=index1;
		    index1=nodempc[3*index1-1];
		}
		if(index1==0) break;
	    }while(1);
        }
        if(ichange==0) break;
    }while(1);
    
    /* decascading using spooles */
    
#ifdef SPOOLES
    if((*icascade==1)&&(ispooles==1)){
	if ( (msgFile = fopen("spooles.out", "a")) == NULL ) {
	    fprintf(stderr, "\n fatal error in spooles.c"
		    "\n unable to open file spooles.out\n") ;
	}
	NNEW(ipointer,ITG,7**nk);
	NNEW(indepdof,ITG,7**nk);
	NNEW(icoef,ITG,2**memmpc_);
	NNEW(xcoef,double,*memmpc_);
	ifree=0;
	nindep=0;

	for(i=*nmpc-1;i>-1;i--){
	    index=ipompc[i];
	    while(1){
		idof=8*(nodempc[3*index-3]-1)+nodempc[3*index-2]-1;

/* check whether idof is a independent dof which has not yet been
   stored in indepdof */

		FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		if((id==0)||(ikmpc[id-1]!=idof)){
		    FORTRAN(nident,(indepdof,&idof,&nindep,&id));
		    if((id==0)||(indepdof[id-1]!=idof)){
			for(j=nindep;j>id;j--){
			    indepdof[j]=indepdof[j-1];
			}
			indepdof[id]=idof;
			nindep++;
		    }
		}

		icoef[2*ifree]=i+1;
		icoef[2*ifree+1]=ipointer[idof];
		xcoef[ifree]=coefmpc[index-1];
		ipointer[idof]=++ifree;
		index=nodempc[3*index-1];
		if(index==0) break;
	    }
	}

/*     filling the left hand side */

	nent=*memmpc_;
	neqns=*nmpc;
	mtxA = InpMtx_new() ;
	InpMtx_init(mtxA, INPMTX_BY_ROWS, type, nent, neqns) ;
	
	for(i=0;i<*nmpc;i++){
	    idof=ikmpc[i];
	    icolumn=ilmpc[i]-1;
	    if(strcmp1(&labmpc[20*icolumn],"RIGID")==0) icolnl=1;
	    else icolnl=0;
	    index=ipointer[idof-1];
	    while(1){
		irow=icoef[2*index-2]-1;
		if(irow!=icolumn){
		    if(strcmp1(&labmpc[20*irow],"RIGID")==0)irownl=1;
		    else irownl=0;
		    if((irownl==1)||(icolnl==1)){
			*icascade=2;
			InpMtx_free(mtxA);
			printf("*ERROR in cascade: linear and nonlinear MPCs depend on each other\n\n");
			FORTRAN(stop,());
		    }
		}
		if((strcmp1(&labmpc[20*irow],"                    ")==0)&&
		   (strcmp1(&labmpc[20*icolumn],"CYCLIC")==0)){
		    strcpy1(&labmpc[20*irow],"SUBCYCLIC",9);}
		coef=xcoef[index-1];
		InpMtx_inputRealEntry(mtxA,irow,icolumn,coef);
		index=icoef[2*index-1];
		if(index==0) break;
	    }
	    ipointer[idof-1]=0;
	}

	InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
	if ( msglvl > 1 ) {
	    fprintf(msgFile, "\n\n input matrix") ;
	    InpMtx_writeForHumanEye(mtxA, msgFile) ;
	    fflush(msgFile) ;
	}
/*--------------------------------------------------------------------*/
/*
  -------------------------------------------------
  STEP 2 : find a low-fill ordering
  (1) create the Graph object
  (2) order the graph using multiple minimum degree
  -------------------------------------------------
*/
	graph = Graph_new() ;
	adjIVL = InpMtx_fullAdjacency(mtxA) ;
	nedges = IVL_tsize(adjIVL) ;
	Graph_init2(graph, 0, neqns, 0, nedges, neqns, nedges, adjIVL,
		    NULL, NULL) ;
	if ( msglvl > 1 ) {
	    fprintf(msgFile, "\n\n graph of the input matrix") ;
	    Graph_writeForHumanEye(graph, msgFile) ;
	    fflush(msgFile) ;
	}
	maxdomainsize=800;maxzeros=1000;maxsize=64;
	/*maxdomainsize=neqns/100;*/
	/*frontETree = orderViaMMD(graph, seed, msglvl, msgFile) ;*/
	/*frontETree = orderViaND(graph,maxdomainsize,seed,msglvl,msgFile); */
	/*frontETree = orderViaMS(graph,maxdomainsize,seed,msglvl,msgFile);*/
	frontETree=orderViaBestOfNDandMS(graph,maxdomainsize,maxzeros,
	  maxsize,seed,msglvl,msgFile);
	if ( msglvl > 1 ) {
	    fprintf(msgFile, "\n\n front tree from ordering") ;
	    ETree_writeForHumanEye(frontETree, msgFile) ;
	    fflush(msgFile) ;
	}
/*--------------------------------------------------------------------*/
/*
  -----------------------------------------------------
  STEP 3: get the permutation, permute the matrix and 
  front tree and get the symbolic factorization
  -----------------------------------------------------
*/
	oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
	oldToNew = IV_entries(oldToNewIV) ;
	newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
	ETree_permuteVertices(frontETree, oldToNewIV) ;
	InpMtx_permute(mtxA, oldToNew, oldToNew) ;
/*	InpMtx_mapToUpperTriangle(mtxA) ;*/
	InpMtx_changeCoordType(mtxA,INPMTX_BY_CHEVRONS);
	InpMtx_changeStorageMode(mtxA,INPMTX_BY_VECTORS);
	symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
	if ( msglvl > 1 ) {
	    fprintf(msgFile, "\n\n old-to-new permutation vector") ;
	    IV_writeForHumanEye(oldToNewIV, msgFile) ;
	    fprintf(msgFile, "\n\n new-to-old permutation vector") ;
	    IV_writeForHumanEye(newToOldIV, msgFile) ;
	    fprintf(msgFile, "\n\n front tree after permutation") ;
	    ETree_writeForHumanEye(frontETree, msgFile) ;
	    fprintf(msgFile, "\n\n input matrix after permutation") ;
	    InpMtx_writeForHumanEye(mtxA, msgFile) ;
	    fprintf(msgFile, "\n\n symbolic factorization") ;
	    IVL_writeForHumanEye(symbfacIVL, msgFile) ;
	    fflush(msgFile) ;
	}
/*--------------------------------------------------------------------*/
/*
  ------------------------------------------
  STEP 4: initialize the front matrix object
  ------------------------------------------
*/
	frontmtx = FrontMtx_new() ;
	mtxmanager = SubMtxManager_new() ;
	SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
	FrontMtx_init(frontmtx, frontETree, symbfacIVL, type, symmetryflag, 
		      FRONTMTX_DENSE_FRONTS, pivotingflag, NO_LOCK, 0, NULL, 
		      mtxmanager, msglvl, msgFile) ;
/*--------------------------------------------------------------------*/
/*
  -----------------------------------------
  STEP 5: compute the numeric factorization
  -----------------------------------------
*/
	chvmanager = ChvManager_new() ;
	ChvManager_init(chvmanager, NO_LOCK, 1) ;
	DVfill(10, cpus, 0.0) ;
	IVfill(20, stats, 0) ;
	rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, 0.0, chvmanager,
					&error,cpus, stats, msglvl, msgFile) ;
	ChvManager_free(chvmanager) ;
	if ( msglvl > 1 ) {
	    fprintf(msgFile, "\n\n factor matrix") ;
	    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
	    fflush(msgFile) ;
	}
	if ( rootchv != NULL ) {
	    fprintf(msgFile, "\n\n matrix found to be singular\n") ;
	    exit(-1) ;
	}
	if(error>=0){
	    fprintf(msgFile,"\n\nerror encountered at front %" ITGFORMAT "",error);
	    exit(-1);
	}
/*--------------------------------------------------------------------*/
/*
  --------------------------------------
  STEP 6: post-process the factorization
  --------------------------------------
*/
	FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
	if ( msglvl > 1 ) {
	    fprintf(msgFile, "\n\n factor matrix after post-processing") ;
	    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
	    fflush(msgFile) ;
	}

/* reinitialize nodempc */

	*mpcfree=1;
	for(j=0;j<*nmpc;j++){
	    ipompc[j]=0;}

/* filling the RHS */

	jrhs=0;
	nrhs=1;
	mtxB=DenseMtx_new();
	mtxX=DenseMtx_new();

	for(i=nindep;i>0;i--){
	    idof=indepdof[i-1];
	    if(ipointer[idof]>0){
		
/*          new RHS column */

		DenseMtx_init(mtxB, type, 0, 0, neqns, nrhs, 1, neqns) ;
		DenseMtx_zero(mtxB) ;

		index=ipointer[idof];
		while(1){
		    irow=icoef[2*index-2]-1;
		    coef=xcoef[index-1];
		    DenseMtx_setRealEntry(mtxB,irow,jrhs,coef);
		    index=icoef[2*index-1];
		    if(index==0) break;
		}

		if ( msglvl > 1 ) {
		    fprintf(msgFile, "\n\n rhs matrix in original ordering") ;
		    DenseMtx_writeForHumanEye(mtxB, msgFile) ;
		    fflush(msgFile) ;
		}

/*--------------------------------------------------------------------*/
/*
  ---------------------------------------------------------
  STEP 8: permute the right hand side into the new ordering
  ---------------------------------------------------------
*/
		DenseMtx_permuteRows(mtxB, oldToNewIV) ;
		if ( msglvl > 1 ) {
		    fprintf(msgFile, "\n\n right hand side matrix in new ordering") ;
		    DenseMtx_writeForHumanEye(mtxB, msgFile) ;
		    fflush(msgFile) ;
		}
/*--------------------------------------------------------------------*/
/*
  -------------------------------
  STEP 9: solve the linear system
  -------------------------------
*/
		DenseMtx_init(mtxX, type, 0, 0, neqns, nrhs, 1, neqns) ;
		DenseMtx_zero(mtxX) ;
		FrontMtx_solve(frontmtx, mtxX, mtxB, mtxmanager,cpus, msglvl, msgFile) ;
		if ( msglvl > 1 ) {
		    fprintf(msgFile, "\n\n solution matrix in new ordering") ;
		    DenseMtx_writeForHumanEye(mtxX, msgFile) ;
		    fflush(msgFile) ;
		}
/*--------------------------------------------------------------------*/
/*
  --------------------------------------------------------
  STEP 10: permute the solution into the original ordering
  --------------------------------------------------------
*/
		DenseMtx_permuteRows(mtxX, newToOldIV) ;
		if ( msglvl > 1 ) {
		    fprintf(msgFile, "\n\n solution matrix in original ordering") ;
		    DenseMtx_writeForHumanEye(mtxX, msgFile) ;
		    fflush(msgFile) ;
		}
	   

		for(j=0;j<*nmpc;j++){
		    b=DenseMtx_entries(mtxX)[j];
		    if(fabs(b)>1.e-10){
			nodempc[3**mpcfree-1]=ipompc[j];
			node=(ITG)((idof+8)/8);
			idir=idof+1-8*(node-1);
			nodempc[3**mpcfree-3]=node;
			nodempc[3**mpcfree-2]=idir;
			coefmpc[*mpcfree-1]=b;
			ipompc[j]=(*mpcfree)++;
			if(*mpcfree>*memmpc_){
			    *memmpc_=(ITG)(1.1**memmpc_);
			    RENEW(nodempc,ITG,3**memmpc_);
			    RENEW(coefmpc,double,*memmpc_);
			}
		    }
		}
	    }
	}
/*--------------------------------------------------------------------*/
/*
   -----------
   free memory
   -----------
*/
	FrontMtx_free(frontmtx) ;
	DenseMtx_free(mtxB) ;
	DenseMtx_free(mtxX) ;
	IV_free(newToOldIV) ;
	IV_free(oldToNewIV) ;
	InpMtx_free(mtxA) ;
	ETree_free(frontETree) ;
	IVL_free(symbfacIVL) ;
	SubMtxManager_free(mtxmanager) ;
	Graph_free(graph) ;

/* diagonal terms */

	for(i=0;i<*nmpc;i++){
	    j=ilmpc[i]-1;
	    idof=ikmpc[i];
	    node=(ITG)((idof+7)/8);
	    idir=idof-8*(node-1);
	    nodempc[3**mpcfree-1]=ipompc[j];
	    nodempc[3**mpcfree-3]=node;
	    nodempc[3**mpcfree-2]=idir;
	    coefmpc[*mpcfree-1]=1.;
	    ipompc[j]=(*mpcfree)++;
	    if(*mpcfree>*memmpc_){
		*memmpc_=(ITG)(1.1**memmpc_);
		RENEW(nodempc,ITG,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
	    }
	}
	
	SFREE(ipointer);SFREE(indepdof);SFREE(icoef);SFREE(xcoef);

	fclose(msgFile);

    }
#endif

/*     determining the effective size of nodempc and coefmpc for
       the reallocation*/

    *mpcend=0;
    *maxlenmpc=0;
    for(i=0;i<*nmpc;i++){
	index=ipompc[i];
	*mpcend=max(*mpcend,index);
	nterm=1;
	while(1){
	    index=nodempc[3*index-1];
	    if(index==0){
		*maxlenmpc=max(*maxlenmpc,nterm);
		break;
	    }
	    *mpcend=max(*mpcend,index);
	    nterm++;
	}
    }

    SFREE(jmpc);

    *nodempcp=nodempc;
    *coefmpcp=coefmpc;
    
    /*          for(i=0;i<*nmpc;i++){
	j=i+1;
	FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
	}*/
    
    return;
}
