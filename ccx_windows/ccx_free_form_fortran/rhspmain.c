/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static char *lakonf1;

static ITG num_cpus,*nef1,*ipnei1,*neifa1,*neiel1,*jq1,*irow1,*ielfa1,
    *ifabou1,*neq1,*nzs1,*neij1;

static double *vfa1,*area1,*advfa1,*xlet1,*cosa1,*volume1,*au1,*ad1,
    *ap1,*xle1,*b1,*xxn1,*hfa1,*gradpel1,*bp1,*xxi1,*xlen1,*cosb1,
    *xxicn1,*flux1,*xxnj1,*gradpcfa1;

void rhspmain(ITG *nef,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *advfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,ITG *neq,
	     ITG *nzs,double *hfa,double *gradpel,
	     double *bp,double *xxi,ITG *neij,double *xlen,
	     ITG *iatleastonepressurebc,double *xxicn,
	      double *flux,double *xxnj,double *gradpcfa,double *cosb){

    ITG i,j;
      
    /* variables for multithreading procedure */
    
    ITG sys_cpus,*ithread=NULL;
    char *env,*envloc,*envsys;

    num_cpus = 0;
    sys_cpus=0;

    /* explicit user declaration prevails */

    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
	sys_cpus=atoi(envsys);
	if(sys_cpus<0) sys_cpus=0;
    }

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_CFD");
    if(envloc){
	num_cpus=atoi(envloc);
	if(num_cpus<0){
	    num_cpus=0;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
	
    }

    /* else global declaration, if any, applies */

    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
	if (env)
	    num_cpus = atoi(env);
	if (num_cpus < 1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }

// next line is to be inserted in a similar way for all other paralell parts

    if(*nef<num_cpus) num_cpus=*nef;
    
    pthread_t tid[num_cpus];

    /* calculating the stiffness and/or mass matrix 
       (symmetric part) */

    nef1=nef;lakonf1=lakonf;ipnei1=ipnei;neifa1=neifa;neiel1=neiel;
    vfa1=vfa;area1=area;advfa1=advfa;xlet1=xlet,cosa1=cosa;volume1=volume;
    jq1=jq;irow1=irow;ap1=ap;ielfa1=ielfa;ifabou1=ifabou;xle1=xle;
    xxn1=xxn;neq1=neq;nzs1=nzs;hfa1=hfa;gradpel1=gradpel;bp1=bp;xxi1=xxi;
    neij1=neij;xlen1=xlen;ad1=ad;au1=au;b1=b;xxicn1=xxicn;flux1=flux;
    xxnj1=xxnj;gradpcfa1=gradpcfa;cosb1=cosb;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)rhspmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

/*  at least one pressure bc is needed. If none is applied,
    the last dof is set to 0 

     a pressure bc is only recognized if not all velocity degrees of
     freedom are prescribed on the same face */

    if(*iatleastonepressurebc==0) b[*nef-1]=0.;
  
  return;

}

/* subroutine for multithreading of rhsp */

void *rhspmt(ITG *i){

    ITG indexad,indexb,nefa,nefb,nefdelta;
    long long indexau;

    indexad=*i**neq1;
    indexau=(long long)*i**nzs1;
    indexb=*i**neq1;
    
// ceil -> floor

    nefdelta=(ITG)floor(*nef1/(double)num_cpus);
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(rhsp,(nef1,lakonf1,ipnei1,neifa1,neiel1,vfa1,area1,
		  advfa1,xlet1,cosa1,volume1,au1,ad1,
		  jq1,irow1,ap1,ielfa1,ifabou1,xle1,b1,xxn1,
		  neq1,nzs1,
		  hfa1,gradpel1,bp1,xxi1,neij1,xlen1,&nefa,&nefb,
		  xxicn1,flux1,xxnj1,gradpcfa1,cosb1));

    return NULL;
}
