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
#include <string.h>
#include <pthread.h>
#include "CalculiX.h"

void stress_sen_2ndorder(double *co,ITG *nk,ITG *kon,ITG *ipkon,
       char *lakon,ITG *ne,
       double *stn,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *emn,
       double *een,ITG *iperturb,double *f,ITG *nactdof,
       double *vold,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,
       double *enern,double *emeini,double *xstaten,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime, ITG *ne0,
       double *xforc, ITG *nforc, double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,ITG *kscale,char *orname,ITG *network,
       ITG *nestart,ITG *neend,ITG *jqs,ITG *irows,ITG *nodedesi,
       double *xdesi,ITG *ndesi,ITG *iobject,ITG *nobject,char *objectset,
       double *g0,double *dgdx,ITG *idesvar1,ITG *idesvar2,ITG *nasym,
       ITG *isolver,double *distmin,ITG *nodeset,double *b1,double *b2,
       double *dgdxdy){

  ITG symmetryflag=0,mt=mi[1]+1,i,j,iactpos,calcul_fn,list,*ilist=NULL,
       calcul_qa,calcul_cauchy,ikin=0,nal,iout=2,icmd=3,nener=0,
      *inum=NULL,nprintl=0;

  double *vnew=NULL,*conew=NULL,*dstn=NULL,*v=NULL,*fn=NULL,
    *stx=NULL,*eei=NULL,qa[4]={0.,0.,-1.,0.},*xstiff=NULL,*ener=NULL,
    *eme=NULL,*dv1=NULL,*dv2=NULL,*dstn1=NULL,*dstn2=NULL,*dstn12=NULL,
    coef1[6]={1,1,0,-1,-1,1},coef2[6]={0,1,1,1,-1,-1};
    
  if(*nasym!=0){symmetryflag=2;}
      
    NNEW(dv1,double,mt**nk);
    NNEW(dv2,double,mt**nk);
    NNEW(vnew,double,mt**nk);
    NNEW(conew,double,3**nk);

    NNEW(eme,double,6*mi[0]**ne);
    NNEW(inum,ITG,*nk);
    
    FORTRAN(resultsnoddir,(nk,dv1,nactdof,b1,ipompc,nodempc,
			   coefmpc,nmpc,mi));
    FORTRAN(resultsnoddir,(nk,dv2,nactdof,b2,ipompc,nodempc,
			   coefmpc,nmpc,mi));
    
    for(j=0;j<6;j++){
      
      /* calculating the perturbed displacements */
      
      for(i=0;i<mt**nk;i++){
	  vnew[i]=vold[i]+coef1[j]*(*distmin)*dv1[i]
	                 +coef2[j]*(*distmin)*dv2[i];
      }
      
      /* copying the unperturbed coordinates */
      
      memcpy(&conew[0],&co[0],sizeof(double)*3**nk);
      
      /* calculating the perturbed coordinates */
      
      iactpos=nodedesi[*idesvar1]-1;
      for(i=0;i<3;i++){
       conew[iactpos*3+i]=co[iactpos*3+i]+xdesi[*idesvar1*3+i];
      }
      iactpos=nodedesi[*idesvar2]-1;
      for(i=0;i<3;i++){
       conew[iactpos*3+i]=co[iactpos*3+i]+xdesi[*idesvar2*3+i];
      }
      
      /* calculating the stress in the perturbed state */
      
      NNEW(v,double,mt**nk);
      NNEW(fn,double,mt**nk);
      NNEW(stx,double,6*mi[0]**ne);
      NNEW(eei,double,6*mi[0]**ne);
      NNEW(dstn,double,(6**nk)*6);
      
      memcpy(&v[0],&vnew[0],sizeof(double)*mt**nk);

      /* setting the output variables */
      
      calcul_fn=0;
      calcul_qa=0;
      calcul_cauchy=1;
      
      list=0;
      FORTRAN(resultsmech,(conew,kon,ipkon,lakon,ne,v,
          stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
          ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
          iprestr,eme,iperturb,fn,&iout,qa,vold,
          nmethod,
          veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
          xstateini,xstiff,xstate,npmat_,matname,mi,ielas,&icmd,
          ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
          springarea,reltime,&calcul_fn,&calcul_qa,&calcul_cauchy,&nener,
          &ikin,&nal,ne0,thicke,emeini,
          pslavsurf,pmastsurf,mortar,clearini,nestart,neend,ielprop,
	  prop,kscale,&list,ilist));

      /* storing results in the .dat file
        extrapolation of integration point values to the nodes
        interpolation of 3d results for 1d/2d elements */
      
      FORTRAN(resultsprint,(conew,nk,kon,ipkon,lakon,ne,v,&dstn[j*6**nk],inum,
              stx,ielorien,norien,orab,t1,ithermal,filab,een,iperturb,fn,
              nactdof,&iout,vold,nodeboun,ndirboun,nboun,nmethod,ttime,xstate,
              epn,mi,
              nstate_,ener,enern,xstaten,eei,set,nset,istartset,iendset,
              ialset,&nprintl,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              nelemload,nload,&ikin,ielmat,thicke,eme,emn,rhcon,nrhcon,shcon,
              nshcon,cocon,ncocon,ntmat_,sideload,icfd,inomat,pslavsurf,islavact,
              cdn,mortar,islavnode,nslavnode,ntie,islavsurf,time,ielprop,prop,
	      veold,ne0,nmpc,ipompc,nodempc,labmpc,energyini,energy,orname,
              xload));
      
      SFREE(v);SFREE(fn);SFREE(stx);SFREE(eei);

    }
      
    /* calculate the first and second order stress derivatives w.r.t.
       idesvar1 and idesvar2 */

    NNEW(dstn1,double,6**nk);
    NNEW(dstn2,double,6**nk);
    NNEW(dstn12,double,6**nk);
    
    for(i=0;i<6**nk;i++){
	dstn1[i]=(dstn[i]-stn[i])/(*distmin);
	dstn2[i]=(dstn[12**nk+i]-stn[i])/(*distmin);
	dstn12[i]=(dstn[6**nk+i]-dstn[30**nk+i]
		   -dstn[18**nk+i]+dstn[24**nk+i])/
	    (4.**distmin**distmin);
    }

    FORTRAN(objective_stress_dx_dy,(nodeset,istartset,iendset,ialset,
       nk,idesvar1,idesvar2,iobject,dgdx,ndesi,nobject,stn,
       dstn1,dstn2,dstn12,objectset,g0,dgdxdy));
    
    SFREE(dstn);SFREE(dstn1);SFREE(dstn2);SFREE(dstn12);

    SFREE(dv1);SFREE(dv2);SFREE(vnew);SFREE(conew);SFREE(eme);SFREE(inum);

}
