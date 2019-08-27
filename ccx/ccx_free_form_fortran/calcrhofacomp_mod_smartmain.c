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

static ITG *nface1,*ielmatf1,*ntmat1_,*mi1,*ielfa1,*ipnei1,*nef1,*num_cpus1,
    *nface1;

static double *vfa1,*shcon1,*vel1,*flux1,*gradpel1,*gradtel1,*xxj1,*xlet1;

void calcrhofacomp_mod_smartmain(ITG *nface,double *vfa,double *shcon,
				 ITG *ielmatf,ITG *ntmat_,ITG *mi,ITG *ielfa,
				 ITG *ipnei,double *vel,ITG *nef,
				 double *flux,double *gradpel,double *gradtel,
				 double *xxj,double *xlet,ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */
    
    vfa1=vfa;shcon1=shcon;ielmatf1=ielmatf;ntmat1_=ntmat_;mi1=mi;ielfa1=ielfa;
    ipnei1=ipnei;vel1=vel;nef1=nef;flux1=flux;gradpel1=gradpel;gradtel1=gradtel;
    xxj1=xxj;xlet1=xlet;num_cpus1=num_cpus;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcrhofacomp_mod_smart1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of calcgammav1 */

void *calcrhofacomp_mod_smart1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(calcrhofacomp_mod_smart,(vfa1,shcon1,ielmatf1,ntmat1_,mi1,ielfa1,
				      ipnei1,vel1,nef1,flux1,gradpel1,gradtel1,
				      xxj1,xlet1,&nfacea,&nfaceb));

    return NULL;
}
