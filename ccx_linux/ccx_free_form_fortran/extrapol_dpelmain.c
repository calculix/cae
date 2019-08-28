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

static double *xrlfa1,*vfap1,*vel1,*c1,*xxn1,*area1,*volume1,
    *gradpcel1,*gradpcfa1,*xxi1,*xlet1,*xxj1;

void extrapol_dpelmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
		      double *vfa,
		      ITG *ifabou,double *xbounact,ITG *nef,double *gradpcel,
		      double *gradpcfa,ITG *neifa,double *rf,double *area,
		      double *volume,double *xle,double *xxi,ITG *icyclic,
		      double *xxn,ITG *ipnei,ITG *ifatie,double *coefmpc,
		      ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,
		      ITG *ifaext,ITG *nfaext,ITG *nactdoh,ITG *iflag,
		      double *xxj,double *xlet,double *c,ITG *num_cpus){

    ITG i;

    double *vfap=NULL;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;

    num_cpus1=num_cpus;
    
    pthread_t tid[*num_cpus];

    NNEW(vfap,double,8**nface);

    /* inter/extrapolation of p' at the center of the elements
       to the center of the faces */

    ielfa1=ielfa;xrlfa1=xrlfa;vfap1=vfap;vel1=vel;ifabou1=ifabou;
    nef1=nef;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_dpel1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* multiple point constraints */

    if(*nmpc>0){
	FORTRAN(applympc_dpel,(nface,ielfa,xrlfa,vel,vfa,
		ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
		xle,xxi,icyclic,xxn,ipnei,ifatie,
                coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
		nactdoh,iflag,xxj,xlet));
    }

    /* calculate the gradient of p' at the center of
       the elements */

    ipnei1=ipnei;neifa1=neifa;vfap1=vfap;area1=area;xxn1=xxn;volume1=volume;
    gradpcel1=gradpcel;nef1=nef;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_dpel2mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);


    if(*iflag==1){

	/* interpolate/extrapolate the p' gradient from the
	   center of the elements to the center of the faces */
	
	ielfa1=ielfa;xrlfa1=xrlfa;icyclic1=icyclic;ifatie1=ifatie;gradpcfa1=gradpcfa;
	gradpcel1=gradpcel;c1=c;ipnei1=ipnei;xxi1=xxi;nface1=nface;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)extrapol_dpel3mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
    SFREE(ithread);

    }

    SFREE(vfap);
	
    /* correct the facial p' gradients:
       Moukalled et al. p 289 */

    ielfa1=ielfa;ipnei1=ipnei;vel1=vel;xlet1=xlet;gradpcfa1=gradpcfa;xxj1=xxj;
    nef1=nef;nface1=nface;
	
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapol_dpel5mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
    
    return;
    
}

/* subroutine for multithreading of extrapol_dpel1 */

void *extrapol_dpel1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_dpel1,(ielfa1,xrlfa1,vfap1,vel1,
       ifabou1,nef1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_dpel2 */

void *extrapol_dpel2mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(extrapol_dpel2,(ipnei1,neifa1,vfap1,area1,xxn1,volume1,gradpcel1,
			   &nefa,&nefb));

    return NULL;
}

/* subroutine for multithreading of extrapol_dpel3 */

void *extrapol_dpel3mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_dpel3,(ielfa1,xrlfa1,icyclic1,ifatie1,gradpcfa1,
			   gradpcel1,c1,ipnei1,xxi1,&nfacea,&nfaceb));

    return NULL;
}

/* subroutine for multithreading of extrapol_dpel5 */

void *extrapol_dpel5mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapol_dpel5,(ielfa1,ipnei1,vel1,xlet1,gradpcfa1,xxj1,
			   nef1,&nfacea,&nfaceb));

    return NULL;
}
