/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "CalculiX.h"

void interpolcycsymcfd(int *nkold, double *cotet, int *neold, int *ipkon,
     int *kon, int **nodempcp, int *ipompc, int *nmpc,
     int *ikmpc, int *ilmpc, double **coefmpcp, char *labmpc,
     int *mpcfree, int *memmpc_, char *lakon,int *nmast,
     int *nslav,int *ithermal,
     double *cs,int *inoslav,int *inomast,int *imast, int *islav){

    /* generate MPC's between new nodes and original nodes:
        (in rectangular coordinates) */

    int j,*nodempc=NULL,idof,id,k,index,number,
        idir,idirmin,idirmax,nodes,nodem;
    
    double *coefmpc=NULL,trabl[7],a[9];

    nodempc=*nodempcp;coefmpc=*coefmpcp;
  
    /* copying the transformation data */

    trabl[0]=cs[5];
    trabl[1]=cs[6];
    trabl[2]=cs[7];
    trabl[3]=cs[8];
    trabl[4]=cs[9];
    trabl[5]=cs[10];
    trabl[6]=-1;

    /* generating MPC's for the master side */

    for(nodes=0;nodes<*nkold;nodes++){
	if(inomast[nodes]==0) continue;
	nodem=inomast[nodes];
	FORTRAN(nident,(imast,&nodem,nmast,&id));
	if(id>0){
	    if(imast[id-1]==nodem) continue;
	}
 		
        /* temperature and pressure degrees of freedom */

	if(*ithermal>1){
	    idirmin=0;idirmax=4;
	}else{
	    idirmin=1;idirmax=4;
	}
	for(idir=idirmin;idir<=idirmax;idir++){
	    if((idir>=1)&&(idir<=3)) continue;
	    idof=8*(nodem-1)+idir;
	    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
	    if(id>0){
		if(ikmpc[id-1]==idof)continue;
	    }
	    (*nmpc)++;
	    ipompc[*nmpc-1]=*mpcfree;
	    strcpy1(&labmpc[20*(*nmpc-1)],"CFDCYCL             ",20);
	    for(k=*nmpc-1;k>id;k--){
		ikmpc[k]=ikmpc[k-1];
		ilmpc[k]=ilmpc[k-1];
	    }
	    ikmpc[id]=idof;
	    ilmpc[id]=*nmpc;
		    
	    /* new node on master side */
	    
	    nodempc[3**mpcfree-3]=nodem;
	    nodempc[3**mpcfree-2]=idir;
	    coefmpc[*mpcfree-1]=1.;
	    index=*mpcfree;
	    *mpcfree=nodempc[3**mpcfree-1];
	    if(*mpcfree==0){
		*mpcfree=*memmpc_+1;
		nodempc[3*index-1]=*mpcfree;
		if(*memmpc_<11)*memmpc_=11;
		*memmpc_=(int)(1.1**memmpc_);
		printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		RENEW(nodempc,int,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
		for(k=*mpcfree;k<*memmpc_;k++){
		    nodempc[3*k-1]=k+1;
		}
		nodempc[3**memmpc_-1]=0;
	    }
		    
	    /* corresponding node in segment */
	    
		nodempc[3**mpcfree-3]=nodes+1;
		nodempc[3**mpcfree-2]=idir;
		coefmpc[*mpcfree-1]=-1.;
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    nodempc[3*index-1]=0;
	}

	/* velocity degrees of freedom: the MPC's are formulated
           in the global rectangular system */

	for(idir=1;idir<4;idir++){

	    FORTRAN(transformatrix,(trabl,&cotet[3*(nodem-1)],a));

	    for(number=1;number<4;number++){
		idof=8*(nodem-1)+number;
		FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		if(id>0){
		    if(ikmpc[id-1]==idof)continue;
		}
		if(fabs(a[3*idir+number-4])<1.e-5) continue;
		(*nmpc)++;
		ipompc[*nmpc-1]=*mpcfree;
		strcpy1(&labmpc[20*(*nmpc-1)],"CFDCYCL             ",20);
		for(k=*nmpc-1;k>id;k--){
		    ikmpc[k]=ikmpc[k-1];
		    ilmpc[k]=ilmpc[k-1];
		}
		ikmpc[id]=idof;
		ilmpc[id]=*nmpc;
		break;
	    }
		    
		/* new master node term in cylindrical coordinates ->
                   first three terms in rectangular coordinates */
	    
	    number--;
	    for(j=1;j<4;j++){
		number++;
		if(number>3) number=1;
		if(fabs(a[3*idir+number-4])<1.e-30) continue;
		nodempc[3**mpcfree-3]=nodem;
		nodempc[3**mpcfree-2]=number;
		coefmpc[*mpcfree-1]=a[3*idir+number-4];
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    }
		    
	    /* corresponding node in segment
               (one term in cylindrical coordinates corresponds
               to three in rectangular coordinates */
	    
	    FORTRAN(transformatrix,(trabl,&cotet[3*nodes],a));
	    
	    for(number=1;number<4;number++){
		if(fabs(a[3*idir+number-4])<1.e-30) continue;
		nodempc[3**mpcfree-3]=nodes+1;
		nodempc[3**mpcfree-2]=number;
		coefmpc[*mpcfree-1]=-a[3*idir+number-4];
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    }
	    nodempc[3*index-1]=0;
	}
    }

    /* generating MPC's for the slave side */

    for(nodem=0;nodem<*nkold;nodem++){
	if(inoslav[nodem]==0) continue;
	nodes=inoslav[nodem];
	FORTRAN(nident,(islav,&nodes,nslav,&id));
	if(id>0){
	    if(islav[id-1]==nodes) continue;
	}
 		
        /* temperature and pressure degrees of freedom */

	if(*ithermal>1){
	    idirmin=0;idirmax=4;
	}else{
	    idirmin=1;idirmax=4;
	}
	for(idir=idirmin;idir<=idirmax;idir++){
	    if((idir>=1)&&(idir<=3)) continue;
	    idof=8*(nodes-1)+idir;
	    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
	    if(id>0){
		if(ikmpc[id-1]==idof)continue;
	    }
	    (*nmpc)++;
	    ipompc[*nmpc-1]=*mpcfree;
	    strcpy1(&labmpc[20*(*nmpc-1)],"CFDCYCL             ",20);
	    for(k=*nmpc-1;k>id;k--){
		ikmpc[k]=ikmpc[k-1];
		ilmpc[k]=ilmpc[k-1];
	    }
	    ikmpc[id]=idof;
	    ilmpc[id]=*nmpc;
		    
	    /* new node on slave side */
	    
	    nodempc[3**mpcfree-3]=nodes;
	    nodempc[3**mpcfree-2]=idir;
	    coefmpc[*mpcfree-1]=1.;
	    index=*mpcfree;
	    *mpcfree=nodempc[3**mpcfree-1];
	    if(*mpcfree==0){
		*mpcfree=*memmpc_+1;
		nodempc[3*index-1]=*mpcfree;
		if(*memmpc_<11)*memmpc_=11;
		*memmpc_=(int)(1.1**memmpc_);
		printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		RENEW(nodempc,int,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
		for(k=*mpcfree;k<*memmpc_;k++){
		    nodempc[3*k-1]=k+1;
		}
		nodempc[3**memmpc_-1]=0;
	    }
		    
	    /* corresponding node in segment */
	    
	    nodempc[3**mpcfree-3]=nodem+1;
	    nodempc[3**mpcfree-2]=idir;
	    coefmpc[*mpcfree-1]=-1.;
	    index=*mpcfree;
	    *mpcfree=nodempc[3**mpcfree-1];
	    if(*mpcfree==0){
		*mpcfree=*memmpc_+1;
		nodempc[3*index-1]=*mpcfree;
		if(*memmpc_<11)*memmpc_=11;
		*memmpc_=(int)(1.1**memmpc_);
		printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		RENEW(nodempc,int,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
		for(k=*mpcfree;k<*memmpc_;k++){
		    nodempc[3*k-1]=k+1;
		}
		nodempc[3**memmpc_-1]=0;
	    }
	    nodempc[3*index-1]=0;
	}

	/* velocity degrees of freedom: the MPC's are formulated
           in the global rectangular system */

	for(idir=1;idir<4;idir++){

	    FORTRAN(transformatrix,(trabl,&cotet[3*(nodes-1)],a));

	    for(number=1;number<4;number++){
		idof=8*(nodes-1)+number;
		FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		if(id>0){
		    if(ikmpc[id-1]==idof)continue;
		}
		if(fabs(a[3*idir+number-4])<1.e-5) continue;
		(*nmpc)++;
		ipompc[*nmpc-1]=*mpcfree;
		strcpy1(&labmpc[20*(*nmpc-1)],"CFDCYCL             ",20);
		for(k=*nmpc-1;k>id;k--){
		    ikmpc[k]=ikmpc[k-1];
		    ilmpc[k]=ilmpc[k-1];
		}
		ikmpc[id]=idof;
		ilmpc[id]=*nmpc;
		break;
	    }
		    
		/* new slave node term in cylindrical coordinates ->
                   first three terms in rectangular coordinates */
	    
	    number--;
	    for(j=1;j<4;j++){
		number++;
		if(number>3) number=1;
		if(fabs(a[3*idir+number-4])<1.e-30) continue;
		nodempc[3**mpcfree-3]=nodes;
		nodempc[3**mpcfree-2]=number;
		coefmpc[*mpcfree-1]=a[3*idir+number-4];
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    }
		    
	    /* subsequent terms (one term in cylindrical coordinates corresponds
               to three in rectangular coordinates */

	    FORTRAN(transformatrix,(trabl,&cotet[3*nodem],a));
	    
	    for(number=1;number<4;number++){
		if(fabs(a[3*idir+number-4])<1.e-30) continue;
		
		/* node is no dependent node of another MPC */ 
		
		nodempc[3**mpcfree-3]=nodem+1;
		nodempc[3**mpcfree-2]=number;
		coefmpc[*mpcfree-1]=-a[3*idir+number-4];
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    }
	    nodempc[3*index-1]=0;
	}
    }

    *nodempcp=nodempc;*coefmpcp=coefmpc;
    
    return;

}
