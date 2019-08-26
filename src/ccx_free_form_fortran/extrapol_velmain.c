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
    *gradvel1,*gradvfa1,*xxi1,*vfa1,*rf1,*xle1,*xlet1,*xxj1;

void extrapol_velmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
		      double *vfa, ITG *ifabou,double *xbounact,ITG *ipnei,
		      ITG *nef,ITG *icyclic,double *c,ITG *ifatie,double *xxn,
		      double *gradvel,double *gradvfa,ITG *neifa,double *rf,
		      double *area,double *volume,double *xle,double *xxi,
		      double *xxj,double *xlet,double *coefmpc,ITG *nmpc,
		      char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
		      ITG *nfaext,ITG *nactdoh,ITG *iitg,ITG *num_cpus){

    ITG i,is,ie,m;

    double *vfap=NULL;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    NNEW(vfap,double,8**nface);

    /* inter/extrapolation of v at the center of the elements
       to the center of the faces */

    ielfa1=ielfa;xrlfa1=xrlfa;icyclic1=icyclic;ifatie1=ifatie;vfap1=vfap;
    vel1=vel;c1=c;ifabou1=ifabou;xbounact1=xbounact;ipnei1=ipnei;
    xxn1=xxn;nef1=nef;num_cpus1=num_cpus;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_vel1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* multiple point constraints */

    if(*nmpc>0){
	is=1;
	ie=3;
	FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfap,coefmpc,
			  nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,
			  ifaext,nfaext));
    }

    /* calculate the gradient of the velocities at the center of
       the elements */

    ipnei1=ipnei;neifa1=neifa;vfap1=vfap;area1=area;xxn1=xxn;volume1=volume;
    gradvel1=gradvel;nef1=nef;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_vel2mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* interpolate/extrapolate the velocity gradient from the
       center of the elements to the center of the faces */

    ielfa1=ielfa;xrlfa1=xrlfa;icyclic1=icyclic;ifatie1=ifatie;gradvfa1=gradvfa;
    gradvel1=gradvel;c1=c;ipnei1=ipnei;xxi1=xxi;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_vel3mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    
    /* inter/extrapolation of v at the center of the elements
       to the center of the faces: taking the skewness of the
       elements into account (using the gradient at the center
       of the faces)
       Moukalled et al. p 279 */

    ielfa1=ielfa;vfa1=vfa;vfap1=vfap;gradvfa1=gradvfa;rf1=rf;ifabou1=ifabou;
    ipnei1=ipnei;xxn1=xxn;vel1=vel;xxi1=xxi;xle1=xle;nef1=nef;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_vel4mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* multiple point constraints */

    if(*nmpc>0){
	is=1;
	ie=3;
	FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfap,coefmpc,
			  nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,
			  ifaext,nfaext));
    }

    /* correction loops */

    for(m=0;m<*iitg;m++){

	/* calculate the gradient of the velocities at the center of
	   the elements */
	
/*	ipnei1=ipnei;neifa1=neifa;vfap1=vfap;area1=area;xxn1=xxn;volume1=volume;
	gradvel1=gradvel;nef1=nef;*/
	ipnei1=ipnei;neifa1=neifa;vfap1=vfa;area1=area;xxn1=xxn;volume1=volume;
	gradvel1=gradvel;nef1=nef;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_vel2mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);
	
	/* interpolate/extrapolate the velocity gradient from the
	   center of the elements to the center of the faces */
	
	ielfa1=ielfa;xrlfa1=xrlfa;icyclic1=icyclic;ifatie1=ifatie;gradvfa1=gradvfa;
	gradvel1=gradvel;c1=c;ipnei1=ipnei;xxi1=xxi;nface1=nface;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_vel3mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);
	
	
	/* inter/extrapolation of v at the center of the elements
	   to the center of the faces: taking the skewness of the
	   elements into account (using the gradient at the center
	   of the faces)
	   Moukalled et al. p 279 */
	
	ielfa1=ielfa;vfa1=vfa;vfap1=vfap;gradvfa1=gradvfa;rf1=rf;ifabou1=ifabou;
	ipnei1=ipnei;xxn1=xxn;vel1=vel;xxi1=xxi;xle1=xle;nef1=nef;nface1=nface;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_vel4mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);
	
	/* multiple point constraints */
	
	if(*nmpc>0){
	    is=1;
	    ie=3;
	    FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfap,coefmpc,
			      nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,
			      ifaext,nfaext));
	}
    }

    SFREE(vfap);
	
    /* correct the facial velocity gradients:
       Moukalled et al. p 289 */

    ielfa1=ielfa;ipnei1=ipnei;vel1=vel;xlet1=xlet;gradvfa1=gradvfa;xxj1=xxj;
    nef1=nef;nface1=nface;
	
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_vel5mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of extrapol_vel1 */

void *extrapol_vel1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;
//    if(*i==1) printf("num_cpus %d %d %d %d",*i,nfacea,nfaceb,*nface1);

    FORTRAN(extrapol_vel1,(ielfa1,xrlfa1,icyclic1,ifatie1,vfap1,vel1,
       c1,ifabou1,xbounact1,ipnei1,xxn1,nef1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_vel2 */

void *extrapol_vel2mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(extrapol_vel2,(ipnei1,neifa1,vfap1,area1,xxn1,volume1,gradvel1,
			   &nefa,&nefb));

    return NULL;
}

/* subroutine for multithreading of extrapol_vel3 */

void *extrapol_vel3mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_vel3,(ielfa1,xrlfa1,icyclic1,ifatie1,gradvfa1,
			   gradvel1,c1,ipnei1,xxi1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_vel4 */

void *extrapol_vel4mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_vel4,(ielfa1,vfa1,vfap1,gradvfa1,rf1,ifabou1,
       ipnei1,xxn1,vel1,xxi1,xle1,nef1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_vel5 */

void *extrapol_vel5mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_vel5,(ielfa1,ipnei1,vel1,xlet1,gradvfa1,xxj1,
			   nef1,&nfacea,&nfaceb));

    return NULL;
}
