/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

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

static ITG num_cpus,*nef1,*ipnei1,*neifa1,*neiel1,*jq1,*irow1,*nzs1,
    *ielfa1,*ifabou1,*nbody1,*neq1,*nactdohinv1,*iau61,*iturbulent1;

static double *au1,*ad1,*b1,*vfa1,*xxn1,*area1,*vel1,
    *umfa1,*xlet1,*xle1,*gradtfa1,*xxi1,*body1,*volume1,*dtimef1,*velo1,
    *veloo1,*cvfa1,*hcfa1,*cvel1,*gradvel1,*xload1,*xrlfa1,
    *xxj1,*a11,*a21,*a31,*flux1,*xxni1,*xxnj1,*f11,*of21,*yy1,*umel1,
    *gradkel1,*gradoel1;

void mafillkcompmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
	       double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *xlet,
               double *xle,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cvfa,double *hcfa,double *cvel,
	       double *gradvel,double *xload,double *xrlfa,
	       double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
	       double *flux,ITG *iau6,double *xxni,double *xxnj,
	       ITG *iturbulent,double *f1,double *of2,double *yy,
	       double *umel,double *gradkel,double *gradoel){

    ITG i;
      
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

    nef1=nef;ipnei1=ipnei;neifa1=neifa;neiel1=neiel;vfa1=vfa;xxn1=xxn;
    area1=area;
    jq1=jq;irow1=irow;nzs1=nzs;vel1=vel;umfa1=umfa;xlet1=xlet;xle1=xle;
    gradtfa1=gradtfa;xxi1=xxi;body1=body;volume1=volume;
    ielfa1=ielfa;lakonf1=lakonf;ifabou1=ifabou;
    nbody1=nbody;neq1=neq;dtimef1=dtimef;velo1=velo;veloo1=veloo;
    cvfa1=cvfa;hcfa1=hcfa;cvel1=cvel;gradvel1=gradvel;xload1=xload;
    xrlfa1=xrlfa;xxj1=xxj;nactdohinv1=nactdohinv;a11=a1;
    a21=a2;a31=a3;flux1=flux;iau61=iau6;ad1=ad;au1=au;b1=b;xxni1=xxni;
    xxnj1=xxnj,iturbulent1=iturbulent,f11=f1,of21=of2;yy1=yy;umel1=umel;
    gradkel1=gradkel;gradoel1=gradoel;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillkcompmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of mafillkcomp */

void *mafillkcompmt(ITG *i){

    ITG nefa,nefb,nefdelta;
    
// ceil -> floor

    nefdelta=(ITG)floor(*nef1/(double)num_cpus);
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(nefb<*nef1)) nefb=*nef1;
	      
    FORTRAN(mafillkcomp,(nef1,ipnei1,neifa1,neiel1,vfa1,xxn1,area1,
			 au1,ad1,jq1,irow1,nzs1,
			 b1,vel1,umfa1,xlet1,xle1,gradtfa1,xxi1,
			 body1,volume1,ielfa1,lakonf1,ifabou1,
			 nbody1,neq1,dtimef1,velo1,veloo1,cvfa1,hcfa1,cvel1,
			 gradvel1,xload1,xrlfa1,xxj1,nactdohinv1,
		         a11,a21,a31,flux1,&nefa,&nefb,iau61,xxni1,xxnj1,
		         iturbulent1,f11,of21,yy1,umel1,gradkel1,gradoel1));

    return NULL;
}
