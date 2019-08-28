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

static ITG *ielfa1,*icyclic1,*ifatie1,*ifabou1,*ipnei1,*nef1,*num_cpus1,
    *neifa1,*nface1;

static double *xrlfa1,*vfap1,*vel1,*c1,*xbounact1,*xxn1,*area1,*volume1,
    *gradtel1,*gradtfa1,*xxi1,*vfa1,*rf1,*xle1,*xlet1,*xxj1,*xload1;

void extrapol_telmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradtel,double *gradtfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *xload,double *xlet,double *xxj,double *coefmpc,
             ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
	     ITG *nfaext,ITG *nactdoh,ITG *iitg,double *c,ITG *num_cpus){

    ITG i,is,ie,m;

    double *vfap=NULL;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;

    num_cpus1=num_cpus;
    
    pthread_t tid[*num_cpus];

    NNEW(vfap,double,8**nface);

    /* inter/extrapolation of T at the center of the elements
       to the center of the faces */

    ielfa1=ielfa;xrlfa1=xrlfa;vfap1=vfap;vel1=vel;ifabou1=ifabou;
    xbounact1=xbounact;nef1=nef;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_tel1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* multiple point constraints */

    if(*nmpc>0){
	is=0;
	ie=0;
	FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfap,coefmpc,
			  nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,
			  ifaext,nfaext));
    }

    /* calculate the gradient of the temperatures at the center of
       the elements */

    ipnei1=ipnei;neifa1=neifa;vfap1=vfap;area1=area;xxn1=xxn;volume1=volume;
    gradtel1=gradtel;nef1=nef;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_tel2mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* interpolate/extrapolate the temperature gradient from the
       center of the elements to the center of the faces */

    ielfa1=ielfa;xrlfa1=xrlfa;icyclic1=icyclic;ifatie1=ifatie;gradtfa1=gradtfa;
    gradtel1=gradtel;c1=c;ipnei1=ipnei;xxi1=xxi;ifabou1=ifabou;xxn1=xxn;
    xload1=xload;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_tel3mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    
    /* inter/extrapolation of T at the center of the elements
       to the center of the faces: taking the skewness of the
       elements into account (using the gradient at the center
       of the faces)
       Moukalled et al. p 279 */

    ielfa1=ielfa;vfa1=vfa;vfap1=vfap;gradtfa1=gradtfa;rf1=rf;ifabou1=ifabou;
    ipnei1=ipnei;vel1=vel;xxi1=xxi;xle1=xle;nef1=nef;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_tel4mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* multiple point constraints */

    if(*nmpc>0){
	is=0;
	ie=0;
	FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfap,coefmpc,
			  nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,
			  ifaext,nfaext));
    }

    /* correction loops */

    for(m=0;m<*iitg;m++){

	/* calculate the gradient of the temperatures at the center of
	   the elements */
	
/*	ipnei1=ipnei;neifa1=neifa;vfap1=vfap;area1=area;xxn1=xxn;volume1=volume;
	gradtel1=gradtel;nef1=nef;*/
	ipnei1=ipnei;neifa1=neifa;vfap1=vfa;area1=area;xxn1=xxn;volume1=volume;
	gradtel1=gradtel;nef1=nef;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_tel2mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);
	
	/* interpolate/extrapolate the temperature gradient from the
	   center of the elements to the center of the faces */
	
	ielfa1=ielfa;xrlfa1=xrlfa;icyclic1=icyclic;ifatie1=ifatie;gradtfa1=gradtfa;
	gradtel1=gradtel;c1=c;ipnei1=ipnei;xxi1=xxi;ifabou1=ifabou;xxn1=xxn;
	xload1=xload;nface1=nface;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_tel3mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);
	
	
	/* inter/extrapolation of T at the center of the elements
	   to the center of the faces: taking the skewness of the
	   elements into account (using the gradient at the center
	   of the faces)
	   Moukalled et al. p 279 */

	ielfa1=ielfa;vfa1=vfa;vfap1=vfap;gradtfa1=gradtfa;rf1=rf;ifabou1=ifabou;
	ipnei1=ipnei;vel1=vel;xxi1=xxi;xle1=xle;nef1=nef;nface1=nface;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_tel4mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);
	
	/* multiple point constraints */
	
	if(*nmpc>0){
	    is=0;
	    ie=0;
	    FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfap,coefmpc,
			      nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,
			      ifaext,nfaext));
	}
    }

    SFREE(vfap);
	
    /* correct the facial temperature gradients:
       Moukalled et al. p 289 */

    ielfa1=ielfa;ipnei1=ipnei;vel1=vel;xlet1=xlet;gradtfa1=gradtfa;xxj1=xxj;
    nef1=nef;nface1=nface;
	
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_tel5mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of extrapol_tel1 */

void *extrapol_tel1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_tel1,(ielfa1,xrlfa1,vfap1,vel1,
       ifabou1,xbounact1,nef1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_tel2 */

void *extrapol_tel2mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(extrapol_tel2,(ipnei1,neifa1,vfap1,area1,xxn1,volume1,gradtel1,
			   &nefa,&nefb));

    return NULL;
}

/* subroutine for multithreading of extrapol_tel3 */

void *extrapol_tel3mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_tel3,(ielfa1,xrlfa1,icyclic1,ifatie1,gradtfa1,
			   gradtel1,c1,ipnei1,xxi1,ifabou1,xxn1,xload1,
			   &nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_tel4 */

void *extrapol_tel4mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_tel4,(ielfa1,vfa1,vfap1,gradtfa1,rf1,ifabou1,
       ipnei1,vel1,xxi1,xle1,nef1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_tel5 */

void *extrapol_tel5mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_tel5,(ielfa1,ipnei1,vel1,xlet1,gradtfa1,xxj1,
			   nef1,&nfacea,&nfaceb));

    return NULL;
}
