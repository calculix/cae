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

void resultsini(ITG *nk,double *v,ITG *ithermal,char *filab,ITG *iperturb,
		 double *f,double *fn,ITG *nactdof,ITG *iout,double *qa,
		 double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
		 double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
		 double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
		 double *cam,ITG *neq,double *veold,double *accold,
		 double *bet,double *gam,double *dtime,ITG *mi,double *vini,
		 ITG *nprint,char *prlab,ITG *intpointvarm,ITG *calcul_fn,
		 ITG *calcul_f,ITG *calcul_qa,ITG *calcul_cauchy,ITG *ikin,
		 ITG *intpointvart,char *typeboun,ITG *num_cpus){

    ITG mt,i,j,node,ndir,ist,index,incrementalmpc;

    double bnac,scal1,scal2,fixed_disp;
    
    mt=mi[1]+1;

    if((*iout!=2)&&(*iout>-1)){
	if((*nmethod!=4)||(iperturb[0]<=1)){

            /* steady state: mechanical (static) */

	    if(ithermal[0]!=2){
		for(i=0;i<*nk;i++){
		    for(j=1;j<mt;j++){
			if(nactdof[mt*i+j]>0){
			    bnac=b[nactdof[mt*i+j]-1];
			}else{
			    continue;
			}
			v[mt*i+j]+=bnac;
			if((iperturb[0]!=0)&&(abs(*nmethod)==1)){
			    if(fabs(bnac)>cam[0]){
				cam[0]=fabs(bnac);
				cam[3]=nactdof[mt*i+j]-0.5;
			    }
			}
		    }
		}
	    }

	    /* steady state: thermal */

	    if(ithermal[0]>1){
		for(i=0;i<*nk;i++){
		    if(nactdof[mt*i]>0){
			bnac=b[nactdof[mt*i]-1];
		    }else{
			continue;
		    }
		    v[mt*i]+=bnac;
		    if((iperturb[0]!=0)&&(abs(*nmethod)==1)){
			if(fabs(bnac)>cam[1]){
			    cam[1]=fabs(bnac);
			    cam[4]=nactdof[mt*i]-0.5;
			}
		    }
		}
	    }
	}else{

	    /* direct integration dynamic step
	       b contains the acceleration increment */

	    if(ithermal[0]!=2){
		iniparll(&mt,nactdof,b,v,veold,accold,bet,gam,
			 dtime,cam,nk,num_cpus);
/*		scal1=*bet**dtime**dtime;
		scal2=*gam**dtime;

		for(i=0;i<*nk;i++){
		    for(j=1;j<mt;j++){
			if(nactdof[mt*i+j]>0){
			    bnac=b[nactdof[mt*i+j]-1];
			}else{
			    continue;
			}
			v[mt*i+j]+=scal1*bnac;
			if(fabs(scal1*bnac)>cam[0]){
			    cam[0]=fabs(scal1*bnac);
			    cam[3]=nactdof[mt*i+j]-0.5;
			}
			veold[mt*i+j]+=scal2*bnac;
			accold[mt*i+j]+=bnac;
		    }
		    }*/
	    }
	    
	    /* transient thermal step */

	    if(ithermal[0]>1){
		for(i=0;i<*nk;i++){
		    veold[mt*i]=0.;
		    if(nactdof[mt*i]>0){
			bnac=b[nactdof[mt*i]-1];
		    }else{
			continue;
		    }
		    v[mt*i]+=bnac;
		    if(fabs(bnac)>cam[1]){
			cam[1]=fabs(bnac);
			cam[4]=nactdof[mt*i]-0.5;
		    }
		    if(cam[2]<fabs(v[mt*i]-vini[mt*i])){
			cam[2]=fabs(v[mt*i]-vini[mt*i]);
		    }
		}
	    }
	}
    }

    /* initialization */

    *calcul_fn=0;
    *calcul_f=0;
    *calcul_qa=0;
    *calcul_cauchy=0;

    /* procedure requirements */
    
    if((iperturb[0]>=2)||((iperturb[0]<=0)&&(*iout<0))){
	if((*iout<1)&&(*iout>-2)){
	    *calcul_fn=1;
	    *calcul_f=1;
	    *calcul_qa=1;
	}else if((*iout!=-2)&&(iperturb[1]==1)){
	    *calcul_cauchy=1;
	}
    }

    /* output requests */
    
    if(*iout>0){
	if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[783],"RFL ")==0)){
	    *calcul_fn=1;
	}
	if(*calcul_fn!=1){
	    for(i=0;i<*nprint;i++){
		if((strcmp1(&prlab[6*i],"RF  ")==0)||
		   (strcmp1(&prlab[6*i],"RFL ")==0)){
		    *calcul_fn=1;
		    break;
		}
	    }
	}
    }

    /* initializing fn */

    if(*calcul_fn==1){
	setpardou(fn,0.,mt**nk,*num_cpus);
    }
    
    /* initializing f */

    if(*calcul_f==1){
	setpardou(f,0.,*neq,*num_cpus);
    }

    /* SPC's and MPC's have to be taken into account for 
       iout=0,1 and -1 */

    if(abs(*iout)<2){

        /* single point constraints */

	for(i=0;i<*nboun;i++){
	    if((ndirboun[i]>mi[1])||(typeboun[i]=='F')) continue;
	    fixed_disp=xboun[i];
	    
         /* a discontinuity in the displacements in an initial
            acceleration step (recognized by the "special" time
            increment) should not lead to a change in 
            the acceleration or velocity; actually, such a 
            discontinuity is not allowed since it leads to
            infinite accelerations */
					      
	    if((*nmethod==4)&&(iperturb[0]>1)){
		ndir=ndirboun[i];
		node=nodeboun[i]-1;
		if(ndir>0){

		    /* bnac is the change in acceleration */
		    
		    bnac=(xboun[i]-v[mt*node+ndir])/(*bet**dtime**dtime);
		    if(floor(*dtime*1.e28+0.5)!=123571113){
			veold[mt*node+ndir]+=*gam**dtime*bnac;
			accold[mt*node+ndir]+=bnac;
		    }
		}
	    }
	    v[mt*(nodeboun[i]-1)+ndirboun[i]]=fixed_disp;
	}

	/* multiple point constraints */

	/* the parameter incrementalmpc indicates whether the
	   incremental displacements enter the mpc or the total 
	   displacements (incrementalmpc=0) */
	
	for(i=0;i<*nmpc;i++){
	    if(strcmp1(&labmpc[20*i],"FLUID")==0) continue;

	    if((strcmp1(&labmpc[20*i],"                    ")==0)||
	       (strcmp1(&labmpc[20*i],"CYCLIC")==0)||
	       (strcmp1(&labmpc[20*i],"SUBCYCLIC")==0)){
		incrementalmpc=0;
	    }else{
		if((*nmethod==2)||(*nmethod==3)||
		   ((iperturb[0]==0)&&(abs(*nmethod)==1))){
		    incrementalmpc=0;
		}else{
		    incrementalmpc=1;
		}
	    }

	    ist=ipompc[i]-1;
	    node=nodempc[3*ist]-1;
	    ndir=nodempc[3*ist+1];
	    
	    if(ndir==0){
		if(ithermal[0]<2) continue;
	    }else if(ndir>mi[1]){
		continue;
	    }else{
		if(ithermal[0]==2) continue;
	    }

	    index=nodempc[3*ist+2]-1;
	    fixed_disp=0.;
	    if(index!=-1){
		do{
		    if(incrementalmpc==0){
			fixed_disp-=coefmpc[index]*
			    v[mt*(nodempc[3*index]-1)+nodempc[3*index+1]];
		    }else{
			fixed_disp-=coefmpc[index]*
			    (v[mt*(nodempc[3*index]-1)+nodempc[3*index+1]]-
			     vold[mt*(nodempc[3*index]-1)+nodempc[3*index+1]]);
		    }
		    index=nodempc[3*index+2]-1;
		    if(index==-1) break;
		}while(1);
	    }
	    fixed_disp/=coefmpc[ist];
	    if(incrementalmpc==1){
		fixed_disp+=vold[mt*node+ndir];
	    }

	    /* a discontinuity in the displacements in an initial
               acceleration step (recognized by the "special" time
               increment) should not lead to a change in 
               the acceleration or velocity; actually, such a 
               discontinuity is not allowed since it leads to
               infinite accelerations */

	    if((*nmethod==4)&&(iperturb[0]>1)){
		if(ndir>0){

		    /* bnac is the change in acceleration */

		    bnac=(fixed_disp-v[mt*node+ndir])/(*bet**dtime**dtime);
		    if(floor(*dtime*1.e28+0.5)!=123571113){
			veold[mt*node+ndir]+=*gam**dtime*bnac;
		    }

		    /* different treatment than for xboun: accold is
		       not included in the above if statement */
		    
		    accold[mt*node+ndir]+=bnac;
		}
	    }
	    v[mt*node+ndir]=fixed_disp;
	}
    }

    /* check whether the kinetic energy has to be calculated  */

    *ikin=0;

    /* dynamic calculations */
    
    if((*nmethod==4)&&(iperturb[0]>1)&&(ithermal[0]<=1)){
	*ikin=1;
    }
    
    /* output requests */

    if(*ikin!=1){
	for(i=0;i<*nprint;i++){
	    if(strcmp1(&prlab[6*i],"ELKE")==0){
		*ikin=1;
	    }
	}
    }

    qa[0]=0.;
    qa[1]=0.;

    *intpointvarm=1;
    *intpointvart=1;

    if((*nmethod>=4)&&(*nmethod<=5)&&(iperturb[0]<2)){
	*intpointvarm=0;
	if((strcmp1(&filab[174],"S")==0)||
	   (strcmp1(&filab[261],"E")==0)||
	   (strcmp1(&filab[348],"RF")==0)||
	   (strcmp1(&filab[435],"PEEQ")==0)||
	   (strcmp1(&filab[522],"ENER")==0)||
	   (strcmp1(&filab[609],"SDV")==0)||
	   (strcmp1(&filab[1044],"ZZS")==0)||
	   (strcmp1(&filab[1044],"ERR")==0)||
	   (strcmp1(&filab[1479],"PHS")==0)||
	   (strcmp1(&filab[1653],"MAXS")==0)||
	   (strcmp1(&filab[2175],"CONT")==0)||
	   (strcmp1(&filab[2262],"CELS")==0)){
	    *intpointvarm=1;
	}
	if(*intpointvarm!=1){
	    for(i=0;i<*nprint;i++){
		if((strcmp1(&prlab[6*i],"S")==0)||
		   (strcmp1(&prlab[6*i],"E")==0)||
		   (strcmp1(&prlab[6*i],"PEEQ")==0)||
		   (strcmp1(&prlab[6*i],"ENER")==0)||
		   (strcmp1(&prlab[6*i],"ELKE")==0)||
		   (strcmp1(&prlab[6*i],"CDIS")==0)||
		   (strcmp1(&prlab[6*i],"CSTR")==0)||
		   (strcmp1(&prlab[6*i],"CELS")==0)||
		   (strcmp1(&prlab[6*i],"SDV")==0)||
		   (strcmp1(&prlab[6*i],"RF")==0)){
		    *intpointvarm=1;
		    break;
		}
	    }
	}

	*intpointvart=0;
	if((strcmp1(&filab[696],"HFL")==0)||
	   (strcmp1(&filab[783],"RFL")==0)){
	    *intpointvart=1;
	}
	if(*intpointvart!=1){
	    for(i=0;i<*nprint;i++){
		if((strcmp1(&prlab[6*i],"HFL")==0)||
		   (strcmp1(&prlab[6*i],"RFL")==0)){
		    *intpointvart=1;
		    break;
		}
	    }
	}

	/* if internal forces are requested integration point
	   values have to be calculated */

	if(*calcul_fn==1){
	    *intpointvarm=1;
	    *intpointvart=1;
	}
    }
    
    return;

}

/*!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine resultsini(nk,v,ithermal,filab,iperturb,f,fn,
     &  nactdof,iout,qa,vold,b,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
     &  veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
     &  intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,
     &  ikin,intpointvart,typeboun)
!
!     initialization 
!
!     1. storing the calculated primary variables nodewise
!     2. inserting the boundary conditions nodewise (SPC's and MPC's)
!     3. determining which derived variables (strains, stresses,
!        internal forces...) have to be calculated
!
      implicit none
!
      character*1 typeboun(*)
      character*6 prlab(*)
      character*20 labmpc(*)
      character*87 filab(*)
!
      integer mi(*),nactdof(0:mi(2),*),nodeboun(*),ndirboun(*),
     &  ipompc(*),nodempc(3,*),mt,nk,ithermal(2),i,j,
     &  iperturb(*),iout,nboun,nmpc,nmethod,ist,ndir,node,index,
     &  neq,incrementalmpc,nprint,ikin,calcul_fn,
     &  calcul_f,calcul_cauchy,calcul_qa,intpointvarm,intpointvart,
     &  irefnode,irotnode,iexpnode,irefnodeprev
!
      real*8 v(0:mi(2),*),vini(0:mi(2),*),f(*),fn(0:mi(2),*),
     &  cam(5),vold(0:mi(2),*),b(*),xboun(*),coefmpc(*),
     &  veold(0:mi(2),*),accold(0:mi(2),*),
     &  qa(*),bet,gam,dtime,scal1,scal2,bnac,
     &  fixed_disp
!
      intent(in) nk,ithermal,filab,iperturb,nactdof,iout,vold,b,
     &  nodeboun,ndirboun,xboun,nboun,ipompc,nodempc,coefmpc,labmpc,
     &  nmpc,nmethod,neq,bet,gam,dtime,mi,vini,nprint,prlab
!
      intent(inout) calcul_fn,calcul_cauchy,calcul_qa,calcul_f,
     &  ikin,v,accold,veold,intpointvarm,intpointvart,qa,f,fn,
     &  cam
!
      mt=mi(2)+1
!
      if((iout.ne.2).and.(iout.gt.-1)) then
!     
         if((nmethod.ne.4).or.(iperturb(1).le.1)) then
            if(ithermal(1).ne.2) then
               do i=1,nk
                  do j=1,mi(2)
                     if(nactdof(j,i).gt.0) then
                        bnac=b(nactdof(j,i))
                     else
                        cycle
                     endif
                     v(j,i)=v(j,i)+bnac
                     if((iperturb(1).ne.0).and.(abs(nmethod).eq.1)) then
                        if(dabs(bnac).gt.cam(1)) then
                           cam(1)=dabs(bnac)
                           cam(4)=nactdof(j,i)-0.5d0
                        endif
                     endif
                  enddo
               enddo
            endif
            if(ithermal(1).gt.1) then
               do i=1,nk
                  if(nactdof(0,i).gt.0) then
                     bnac=b(nactdof(0,i))
                  else
                     cycle
                  endif
                  v(0,i)=v(0,i)+bnac
                  if((iperturb(1).ne.0).and.(abs(nmethod).eq.1)) then
                     if(dabs(bnac).gt.cam(2)) then
                        cam(2)=dabs(bnac)
                        cam(5)=nactdof(0,i)-0.5d0
                     endif
                  endif
               enddo
            endif
!     
         else
!     
!     direct integration dynamic step
!     b contains the acceleration increment
!     
            if(ithermal(1).ne.2) then
               scal1=bet*dtime*dtime
               scal2=gam*dtime
               do i=1,nk
                  do j=1,mi(2)
                     if(nactdof(j,i).gt.0) then
                        bnac=b(nactdof(j,i))
                     else
                        cycle
                     endif
                     v(j,i)=v(j,i)+scal1*bnac
                     if(dabs(scal1*bnac).gt.cam(1)) then
                        cam(1)=dabs(scal1*bnac)
                        cam(4)=nactdof(j,i)-0.5d0
                     endif
                     veold(j,i)=veold(j,i)+scal2*bnac
                     accold(j,i)=accold(j,i)+bnac
                  enddo
               enddo
            endif
            if(ithermal(1).gt.1) then
               do i=1,nk
!
!                 no extrapolation is performed for the solution
!                 of parabolic differential equations (e.g. heat
!                 equation)
!
                  veold(0,i)=0.d0
                  if(nactdof(0,i).gt.0) then
                     bnac=b(nactdof(0,i))
                  else
                     cycle
                  endif
                  v(0,i)=v(0,i)+bnac
                  if(dabs(bnac).gt.cam(2)) then
                     cam(2)=dabs(bnac)
                     cam(5)=nactdof(0,i)-0.5d0
                  endif
                  cam(3)=max(cam(3),dabs(v(0,i)-vini(0,i)))
               enddo
            endif
         endif
!     
      endif
!     
!     initialization
!     
      calcul_fn=0
      calcul_f=0
      calcul_qa=0
      calcul_cauchy=0
!     
!     determining which quantities have to be calculated
!     
      if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.0))) 
c      if((iperturb(1).ge.2).or.((iperturb(1).le.1).and.(iout.lt.0))) 
     &     then
         if((iout.lt.1).and.(iout.gt.-2)) then
            calcul_fn=1
            calcul_f=1
            calcul_qa=1
         elseif((iout.ne.-2).and.(iperturb(2).eq.1)) then    
            calcul_cauchy=1
         endif
      endif
!     
      if(iout.gt.0) then
         if((filab(5)(1:4).eq.'RF  ').or.
     &        (filab(10)(1:4).eq.'RFL ')) then
            calcul_fn=1
         else
            do i=1,nprint
               if((prlab(i)(1:4).eq.'RF  ').or.
     &              (prlab(i)(1:4).eq.'RFL ')) then
                  calcul_fn=1
                  exit
               endif
            enddo
         endif
      endif
!     
!     initializing fn
!     
      if(calcul_fn.eq.1) then
         do i=1,nk
            do j=0,mi(2)
               fn(j,i)=0.d0
            enddo
         enddo
      endif
!     
!     initializing f
!     
      if(calcul_f.eq.1) then
         do i=1,neq
            f(i)=0.d0
         enddo
      endif
!     
!     SPC's and MPC's have to be taken into account for 
!     iout=0,1 and -1
!     
      if(abs(iout).lt.2) then
!     
!     inserting the boundary conditions
!     
         do i=1,nboun
            if((ndirboun(i).gt.mi(2)).or.(typeboun(i).eq.'F')) cycle
            fixed_disp=xboun(i)
!
!           a discontinuity in the displacements in an initial
!           acceleration step (recognized by the "special" time
!           increment) should not lead to a change in 
!           the acceleration or velocity; actually, such a 
!           discontinuity is not allowed since it leads to
!           infinite accelerations
!
            if((nmethod.eq.4).and.(iperturb(1).gt.1)) then
               ndir=ndirboun(i)
               node=nodeboun(i)
               if(ndir.gt.0) then
!
!                 bnac is the change of acceleration
!
                  bnac=(xboun(i)-v(ndir,node))/
     &                 (bet*dtime*dtime)
                  if(nint(dtime*1.d28).ne.123571113) then
                     veold(ndir,node)=veold(ndir,node)+
     &                    gam*dtime*bnac                     
                     accold(ndir,node)=accold(ndir,node)+bnac
                  endif
c                  accold(ndir,node)=accold(ndir,node)+bnac
               endif
            endif
            v(ndirboun(i),nodeboun(i))=fixed_disp
         enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
c     
c     to be checked: should replace the lines underneath do i=1,nmpc
c     
c     incrementalmpc=iperturb(2)
!
         do i=1,nmpc
            if(labmpc(i)(1:5).eq.'FLUID') cycle
            if((labmpc(i)(1:20).eq.'                    ').or.
     &           (labmpc(i)(1:6).eq.'CYCLIC').or.
     &           (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
               incrementalmpc=0
            else
               if((nmethod.eq.2).or.(nmethod.eq.3).or.
     &              ((iperturb(1).eq.0).and.(abs(nmethod).eq.1)))
     &              then
                  incrementalmpc=0
               else
                  incrementalmpc=1
               endif
            endif
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.eq.0) then
               if(ithermal(1).lt.2) cycle
            elseif(ndir.gt.mi(2)) then
               cycle
            else
               if(ithermal(1).eq.2) cycle
            endif
            index=nodempc(3,ist)
            fixed_disp=0.d0
            if(index.ne.0) then
               do
                  if(incrementalmpc.eq.0) then
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    v(nodempc(2,index),nodempc(1,index))
                  else
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    (v(nodempc(2,index),nodempc(1,index))-
     &                    vold(nodempc(2,index),nodempc(1,index)))
                  endif
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
            fixed_disp=fixed_disp/coefmpc(ist)
            if(incrementalmpc.eq.1) then
               fixed_disp=fixed_disp+vold(ndir,node)
            endif
!
!           a discontinuity in the displacements in an initial
!           acceleration step (recognized by the "special" time
!           increment) should not lead to a change in 
!           the acceleration or velocity; actually, such a 
!           discontinuity is not allowed since it leads to
!           infinite accelerations
!
            if((nmethod.eq.4).and.(iperturb(1).gt.1)) then
               if(ndir.gt.0) then
!
!                 bnac is the change of acceleration
!
                  bnac=(fixed_disp-v(ndir,node))/
     &                 (bet*dtime*dtime)
                  if(nint(dtime*1.d28).ne.123571113) then
                     veold(ndir,node)=veold(ndir,node)+
     &                    gam*dtime*bnac
                  endif
!
!     different treatment than for xboun: accold is not included in the
!     above if statement
!                  
                  accold(ndir,node)=accold(ndir,node)+bnac
               endif
            endif
            v(ndir,node)=fixed_disp
         enddo
      endif
!
!     storing the knot information in the .dat-file
!
      if(ithermal(1).ne.2) then
         irefnodeprev=0
         do i=1,nmpc
            if(iout.gt.0) then
               if(labmpc(i)(1:4).eq.'KNOT') then
                  irefnode=nodempc(1,nodempc(3,ipompc(i)))
                  if(irefnode.ne.irefnodeprev) then
                     irefnodeprev=irefnode
                     iexpnode=nodempc(1,nodempc(3,nodempc(3,ipompc(i))))
                     if(labmpc(i)(5:5).ne.'2') then
                        irotnode=nodempc(1,nodempc(3,nodempc(3,
     &                       nodempc(3,ipompc(i)))))
                     else
                        irotnode=nodempc(1,nodempc(3,nodempc(3,
     &                     nodempc(3,nodempc(3,nodempc(3,ipompc(i)))))))
                     endif
c                     write(5,*)
c                     write(5,'(a5)') labmpc(i)(1:5)
c                     write(5,'("tra",i10,3(1x,e11.4))')
c     &                    irefnode,(v(j,irefnode),j=1,3)
c                     write(5,'("rot",i10,3(1x,e11.4))')
c     &                    irotnode,(v(j,irotnode),j=1,3)
c                     if(labmpc(i)(5:5).eq.'2') then
c                        write(5,'("exp",i10,3(1x,e11.4))')
c     &                       iexpnode,(v(j,iexpnode),j=1,3)
c                     else
c                        write(5,'("exp",i10,3(1x,e11.4))')
c     &                       iexpnode,v(1,iexpnode)
c                     endif
                  endif
               endif
            endif
         enddo
      endif
!     
!     check whether there are any strain output requests
!     
      ikin=0
!
!     *DYNAMIC
!
      if((nmethod.eq.4).and.(iperturb(1).gt.1).and.
     &   (ithermal(1).le.1)) then
         ikin=1
      endif
!   
!     dat-output
!   
      do i=1,nprint
         if(prlab(i)(1:4).eq.'ELKE') then
            ikin=1
         endif
      enddo
!     
      qa(1)=0.d0
      qa(2)=0.d0
!
!     check whether integration point variables are needed in
!     modal dynamics and steady state dynamics calculations
!
      intpointvarm=1
      intpointvart=1
!
      if((nmethod.ge.4).and.(nmethod.le.5).and.(iperturb(1).lt.2)) then
         intpointvarm=0
         if((filab(3)(1:4).eq.'S   ').or.
     &      (filab(4)(1:4).eq.'E   ').or.
     &      (filab(5)(1:4).eq.'RF  ').or.
     &      (filab(6)(1:4).eq.'PEEQ').or.
     &      (filab(7)(1:4).eq.'ENER').or.
     &      (filab(8)(1:4).eq.'SDV ').or.
     &      (filab(13)(1:4).eq.'ZZS ').or.
     &      (filab(13)(1:4).eq.'ERR ').or.
     &      (filab(18)(1:4).eq.'PHS ').or.
     &      (filab(20)(1:4).eq.'MAXS').or.
     &      (filab(26)(1:4).eq.'CONT').or.
     &      (filab(27)(1:4).eq.'CELS')) intpointvarm=1
         do i=1,nprint
            if((prlab(i)(1:4).eq.'S   ').or.
     &           (prlab(i)(1:4).eq.'E   ').or.
     &           (prlab(i)(1:4).eq.'PEEQ').or.
     &           (prlab(i)(1:4).eq.'ENER').or.
     &           (prlab(i)(1:4).eq.'ELKE').or.
     &           (prlab(i)(1:4).eq.'CDIS').or.
     &           (prlab(i)(1:4).eq.'CSTR').or.
     &           (prlab(i)(1:4).eq.'CELS').or.
     &           (prlab(i)(1:4).eq.'SDV ').or.
     &           (prlab(i)(1:4).eq.'RF  ')) then
               intpointvarm=1
               exit
            endif
         enddo
!
         intpointvart=0
         if((filab(9)(1:4).eq.'HFL ').or.
     &      (filab(10)(1:4).eq.'RFL ')) intpointvart=1
         do i=1,nprint
            if((prlab(i)(1:4).eq.'HFL ').or.
     &         (prlab(i)(1:4).eq.'RFL ')) intpointvart=1
         enddo
!
!        if internal forces are requested integration point
!        values have to be calculated
!
         if(calcul_fn.eq.1) then
            intpointvarm=1
            intpointvart=1
         endif
      endif
!     
      return
      end*/
