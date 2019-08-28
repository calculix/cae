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

static ITG *nef1,*nshcon1,*ielmatf1,*ntmat1_,*mi1,*ithermal1,
    *ncocon1,*ielfa1,*num_cpus1,*nface1;

static double *vel1,*shcon1,*cvel1,*physcon1,*umel1,*vfa1,*cocon1,
    *cvfa1,*umfa1,*hcfa1;

void materialdata_cfd_compmain(ITG *nef,double *vel,double *shcon,ITG *nshcon,
			  ITG *ielmatf,ITG *ntmat_,ITG *mi,double *cvel,
			  double *vfa,double *cocon,ITG *ncocon,
			  double *physcon,double *cvfa,ITG *ithermal,
			  ITG *nface,double *umel,double *umfa,ITG *ielfa,
			  double *hcfa,ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */

    nef1=nef;vel1=vel;shcon1=shcon;nshcon1=nshcon;ielmatf1=ielmatf;
    ntmat1_=ntmat_;mi1=mi;cvel1=cvel;physcon1=physcon;ithermal1=ithermal;
    umel1=umel;num_cpus1=num_cpus;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)materialdata_cfd_comp1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* calculate the density at the face centers */

    shcon1=shcon;nshcon1=nshcon;ielmatf1=ielmatf;ntmat1_=ntmat_;mi1=mi;vfa1=vfa;
    cocon1=cocon;ncocon1=ncocon;physcon1=physcon;cvfa1=cvfa;ithermal1=ithermal;
    umfa1=umfa;ielfa1=ielfa;hcfa1=hcfa;
    num_cpus1=num_cpus;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)materialdata_cfd_comp2mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of materialdatacfd1 */

void *materialdata_cfd_comp1mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(materialdata_cfd_comp1,(nef1,vel1,shcon1,nshcon1,ielmatf1,
			       ntmat1_,mi1,cvel1,physcon1,ithermal1,
			       umel1,&nefa,&nefb));

    return NULL;
}

/* subroutine for multithreading of materialdata_cfd_comp2 */

void *materialdata_cfd_comp2mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(materialdata_cfd_comp2,(shcon1,nshcon1,ielmatf1,ntmat1_,
			       mi1,vfa1,cocon1,ncocon1,physcon1,
			       cvfa1,ithermal1,umfa1,ielfa1,hcfa1,
			       &nfacea,&nfaceb));

    return NULL;
}
