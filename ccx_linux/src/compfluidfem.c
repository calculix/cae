/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2020 Guido Dhondt                     */

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
#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif

static char *lakon1,*sideload1, *matname1, *sideface1;

static ITG *nk1,*kon1,*ipkon1,*ne1,*nodeboun1,*ndirboun1,*nboun1,*ipompc1,
  *nodempc1,*nmpc1,*nodeforc1,*ndirforc1,*nforc1,*nelemload1,*nload1,
  *ipobody1,*nbody1,*nactdoh1,*icolv1,*jqv1,*irowv1,neqv1,nzlv1,*nmethod1,
  *ikmpc1,*ilmpc1,*ikboun1,*ilboun1,*nrhcon1,*ielmat1,*ntmat_1,*ithermal1,
  nzsv1,*mi1,*ncmat_1,*nshcon1,*istep1,*iinc1,*ibody1,*iturbulent1,
  *nelemface1,*nface1,compressible1,num_cpus,*icolp1,*jqp1,*irowp1,
  neqp1,nzlp1,nzsp1,iexplicit1,*ncocon1,neqt1,nzst1,*ipvar1,*ipvarf1,
  *nactdok1,neqk1,nzsk1,*isolidsurf1,*nsolidsurf1,*ifreestream1,
  *nfreestream1;

static double *co1,*xboun1,*coefmpc1,*xforc1,*xload1,*xbody1,*rhcon1,*t01,
  *vold1,*vcon1,dtimef1,*physcon1,*shcon1,*ttime1,timef1,*xloadold1,
  *vcontu1,*yy1,*b=NULL,*xbounact1,theta11,*v1,theta21,*cocon1,
  reltimef1,*dt1,*var1,*varf1,*sti1,*bk=NULL,*bt=NULL,*xsolidsurf1,
  *ck=NULL,*ct=NULL;

void compfluidfem(double **cop,ITG *nk,ITG **ipkonp,ITG **konp,char **lakonp,
		  ITG *ne,char **sidefacep,ITG *ifreestream,
		  ITG *nfreestream,ITG *isolidsurf,ITG *neighsolidsurf,
		  ITG *nsolidsurf,ITG **iponoelp,ITG **inoelp,ITG *nshcon,
		  double *shcon,
		  ITG *nrhcon,double *rhcon,double **voldp,ITG *ntmat_,
		  ITG *nodeboun,
		  ITG *ndirboun,ITG *nboun,ITG **ipompcp,ITG **nodempcp,
		  ITG *nmpc,
		  ITG **ikmpcp,ITG **ilmpcp,ITG *ithermal,ITG *ikboun,
		  ITG *ilboun,
		  ITG *iturbulent,ITG *isolver,ITG *iexpl,double *vcontu,
		  double *ttime,
		  double *time,double *dtime,ITG *nodeforc,ITG *ndirforc,
		  double *xforc,
		  ITG *nforc,ITG *nelemload,char *sideload,double *xload,
		  ITG *nload,
		  double *xbody,ITG *ipobody,ITG *nbody,ITG **ielmatp,
		  char *matname,
		  ITG *mi,ITG *ncmat_,double *physcon,ITG *istep,ITG *iinc,
		  ITG *ibody,double *xloadold,double *xboun,
		  double **coefmpcp,ITG *nmethod,double *xforcold,
		  double *xforcact,
		  ITG *iamforc,ITG *iamload,double *xbodyold,double *xbodyact,
		  double *t1old,double *t1,double *t1act,ITG *iamt1,
		  double *amta,
		  ITG *namta,ITG *nam,double *ampli,double *xbounold,
		  double *xbounact,
		  ITG *iamboun,ITG *itg,ITG *ntg,char *amname,double *t0,
		  ITG **nelemfacep,
		  ITG *nface,double *cocon,ITG *ncocon,double *xloadact,
		  double *tper,
		  ITG *jmax,ITG *jout,char *set,ITG *nset,ITG *istartset,
		  ITG *iendset,ITG *ialset,char *prset,char *prlab,ITG *nprint,
		  double *trab,ITG *inotr,ITG *ntrans,char *filab,
		  char **labmpcp,
		  double *sti,ITG *norien,double *orab,char *jobnamef,
		  char *tieset,
		  ITG *ntie,ITG *mcs,ITG *ics,double *cs,ITG *nkon,ITG *mpcfree,
		  ITG *memmpc_,double **fmpcp,ITG *nef,ITG **inomatp,
		  double *qfx){

  /* main computational fluid dynamics routine */

  /* References:

     Zienkiewicz, O.C., Taylor, R.L. and Nithiarasu, P., "The Finite
     Element Method for Fluid Dynamics", 6th Edition, Elsevier (2006)

     Menter, F.R., "Two-Equation Eddy-Viscosity Turbulence Models
     for Engineering Applications", AIAA Journal(1994), 32(8), 
     1598-1605                                                       */
  
  char cflag[1],*labmpc=NULL,*lakon=NULL,*sideface=NULL;

  ITG *ipointer=NULL,*mast1=NULL,*irowt=NULL,*irowv=NULL,*irowp=NULL,
    *irowk=NULL,*icolt=NULL,*icolv=NULL,*icolp=NULL,*icolk=NULL,
    *jqt=NULL,*jqv=NULL,*jqp=NULL,*jqk=NULL,*nactdoh=NULL,i,j,k,
    *nactdok=NULL,*nx=NULL,*ny=NULL,*nz=NULL,nzs,neqt,neqv,neqp,
    neqk,nzst,nzsv,nzsp,nzsk,iexplicit,nzlt,nzlv,nzlp,nzlk,kode,nnstep,
    iconvergence,iout,iit,symmetryflag=0,inputformat=0,compressible,
    nmethodd,nstate_=0,*ielorien=NULL,*inum=NULL,ismooth=0,iqfx=0,isti=0,
    ikin=0,mt=mi[1]+1,*ipvar=NULL,*ipvarf=NULL,nvar_,nvarf_,
    nfield,ndim,iorienglob,icfdout=1,force=0,euler=1,*ithread=NULL,
    *integerglob=NULL,nslav,*islav=NULL,ncs,nmast,*imast=NULL,
    nkref,neref,*nodempc=NULL,*ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,
    *ipkon=NULL,*kon=NULL,*ielmat=NULL,*nelemface=NULL,*inoel=NULL,
    *iponoel=NULL,*ipoface=NULL,*nodface=NULL,inoelfree,*inoslav=NULL,
    *inomast=NULL,*ielslav=NULL,*ielmast=NULL,*nr=NULL,*inomat=NULL,
    nstart=501,memmpcref_,mpcfreeref,nmpcref,nefref,*nodempcref=NULL,
    ithermalref,nrhs=1,ipower=1;

  double *yy=NULL,*xsolidsurf=NULL,*dt=NULL,*vcon=NULL,*x=NULL,
    *y=NULL,*z=NULL,*xo=NULL,*yo=NULL,*zo=NULL,*adbt=NULL,
    *aubt=NULL,*adbv=NULL,*aubv=NULL,*adbp=NULL,*aubp=NULL,
    *adbk=NULL,*aubk=NULL,*v=NULL,*vtu=NULL,timef,ttimef,
    dtimef,*addiv=NULL,*sol=NULL,*aux=NULL,shockscale,*stn=NULL,
    *solk=NULL,*solt=NULL,theta1,theta2,*adb=NULL,*qfn=NULL,
    *aub=NULL,sigma=0.,*dh=NULL,reltimef,*fn=NULL,*thicke=NULL,
    *eme=NULL,*xstate=NULL,*ener=NULL,*adlkk=NULL,*adlkt=NULL,
    csmooth=0.,shockcoef,*sa=NULL,*sav=NULL,*varf=NULL,
    *adlt=NULL,*adlv=NULL,*adlp=NULL,*adlk=NULL,factor=1.,*var=NULL,
    *vconini=NULL,*doubleglob=NULL,*coefmpc=NULL,*fmpc=NULL,
    *co=NULL,*vold=NULL,*rcs=NULL,*zcs=NULL,*rcs0=NULL,*zcs0=NULL,
    sum=0.,sumx=0.,sumxx=0.,sumy[7]={0.,0.,0.,0.,0.,0.,0.},
    sumxy[7]={0.,0.,0.,0.,0.,0.,0.},*del=NULL,reltime,*coefmpcref=NULL,
    *dhel=NULL,*voldo=NULL,dtimefo,temp,ratio,ratio2;

  nodempc=*nodempcp;ipompc=*ipompcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;fmpc=*fmpcp;co=*cop;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielmat=*ielmatp;
  nelemface=*nelemfacep;sideface=*sidefacep;inoel=*inoelp;
  iponoel=*iponoelp;vold=*voldp;inomat=*inomatp;

  /* standard: shockcoef=0 */

#ifdef SGI
  ITG token;
#endif

  /* open frd-file for fluids */

  FORTRAN(openfilefluid,(jobnamef));

  /* variables for multithreading procedure */

  ITG sys_cpus;
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
  
  if(*ne<num_cpus) num_cpus=*ne;
  
  printf(" Using up to %d cpu(s) for CFD.\n", num_cpus);
  
  pthread_t tid[num_cpus];
  
  kode=0;
  
  /*  *iexpl==0:  structure:implicit, fluid:semi-implicit
   *iexpl==1:  structure:implicit, fluid:explicit
   *iexpl==2:  structure:explicit, fluid:semi-implicit
   *iexpl==3:  structure:explicit, fluid:explicit */

  if((*iexpl==1)||(*iexpl==3)){
    iexplicit=1;theta1=0.5;theta2=0.;compressible=1;
  }else{
    iexplicit=0;
    theta1=1.0;theta2=1.0;compressible=0;
  }

  /* if initial conditions are specified for the temperature, 
     it is assumed that the temperature is an unknown */

  ithermalref=*ithermal;
  if(*ithermal==1){
    *ithermal=2;
  }

  /* generating additional elements and MPC's for cyclic symmetric
     structures */

  if(*mcs==1){

    /* find the cyclic symmetry slave nodes */

    /* ncs is the number of fluid and solid master nodes,
       nslav is the number of fluid master nodes only */

    ncs=cs[3];

    NNEW(islav,ITG,ncs);
    NNEW(imast,ITG,ncs);

    NNEW(nr,ITG,ncs);
    NNEW(nz,ITG,ncs);
    NNEW(rcs,double,ncs);
    NNEW(zcs,double,ncs);
    NNEW(rcs0,double,ncs);
    NNEW(zcs0,double,ncs);

    NNEW(inoslav,ITG,*nk);
    NNEW(inomast,ITG,*nk);

    FORTRAN(findslavcfd,(nmpc,labmpc,ipompc,nodempc,islav,&nslav,
			 inoslav,inomast,ics,cs,imast,&nmast,co,inomat,
			 nr,nz,rcs,zcs,rcs0,zcs0,&ncs));
    RENEW(islav,ITG,nslav);
    RENEW(imast,ITG,nmast);

    SFREE(nr);SFREE(nz);SFREE(rcs);SFREE(zcs);SFREE(rcs0);SFREE(zcs0);

    /* generate new nodes and elements in a layer on the slave
       and master side */

    nkref=*nk;neref=*ne;nefref=*nef;

    NNEW(ielslav,ITG,*ne);
    NNEW(ielmast,ITG,*ne);

    RENEW(co,double,3*(3**nk));
    RENEW(vold,double,mt*(3**nk));
    RENEW(ipkon,ITG,3**ne);
    RENEW(lakon,char,8*(3**ne));
    RENEW(kon,ITG,*nkon+8*2**ne);
    RENEW(ielmat,ITG,mi[2]*3**ne);
    RENEW(inomat,ITG,3**nk);

    FORTRAN(gencycsymelemcfd,(cs,islav,
			      &nslav,imast,&nmast,inomat,
			      nk,co,ne,ipkon,lakon,kon,nkon,ielmat,mi,vold,
			      ielslav,ielmast,inoslav,inomast,iponoel,inoel));

    SFREE(ielslav);SFREE(ielmast);

    RENEW(co,double,3**nk);
    RENEW(vold,double,mt**nk);
    RENEW(ipkon,ITG,*ne);
    RENEW(lakon,char,8**ne);
    RENEW(kon,ITG,*nkon);
    RENEW(ielmat,ITG,mi[2]**ne);
    RENEW(inomat,ITG,*nk);

    /* generate new MPC's */

    if(*ithermal>1){
      RENEW(ipompc,ITG,*nmpc+5*(2**nk));
      RENEW(ikmpc,ITG,*nmpc+5*(2**nk));
      RENEW(ilmpc,ITG,*nmpc+5*(2**nk));
      RENEW(labmpc,char,20*(*nmpc+5*(2**nk))+1);
    }else{
      RENEW(ipompc,ITG,*nmpc+4*(2**nk));
      RENEW(ikmpc,ITG,*nmpc+4*(2**nk));
      RENEW(ilmpc,ITG,*nmpc+4*(2**nk));
      RENEW(labmpc,char,20*(*nmpc+4*(2**nk))+1);
    }

    memmpcref_=*memmpc_;mpcfreeref=*mpcfree;nmpcref=*nmpc;
    NNEW(nodempcref,ITG,3**memmpc_);
    for(k=0;k<3**memmpc_;k++){nodempcref[k]=nodempc[k];}
    NNEW(coefmpcref,double,*memmpc_);
    for(k=0;k<*memmpc_;k++){coefmpcref[k]=coefmpc[k];}

    interpolcycsymcfd(&nkref,co,&neref,ipkon,kon,&nodempc,ipompc,nmpc,
		      ikmpc,ilmpc,&coefmpc,labmpc,mpcfree,memmpc_,lakon,
		      &nmast,&nslav,ithermal,cs,inoslav,inomast,imast,islav);

    RENEW(ipompc,ITG,*nmpc);
    RENEW(ikmpc,ITG,*nmpc);
    RENEW(ilmpc,ITG,*nmpc);
    RENEW(labmpc,char,20**nmpc);

    SFREE(inoslav);SFREE(inomast);

    /* decascading the new MPC's */

    ITG mpcend,mpcmult,maxlenmpc;
    ITG icascade=0;
    ITG callfrommain=0;
    ITG iperturb=0;
    cascadefem(ipompc,&coefmpc,&nodempc,nmpc,mpcfree,nodeboun,ndirboun,
	       nboun,ikmpc,ilmpc,ikboun,ilboun,&mpcend,&mpcmult,labmpc,
	       nk,memmpc_,&icascade,&maxlenmpc,&callfrommain,&iperturb,
	       ithermal);

    /* update the field with external faces
       and node-to-element dependence */

    SFREE(iponoel);SFREE(inoel);SFREE(sideface);SFREE(nelemface);
    *nef+=(*ne-neref);
    NNEW(sideface,char,6**nef);
    NNEW(nelemface,ITG,6**nef);
    NNEW(ipoface,ITG,*nk);
    NNEW(nodface,ITG,5*6**nef);
    NNEW(iponoel,ITG,*nk);
    NNEW(inoel,ITG,3*20**nef);

    FORTRAN(precfdcyc,(nelemface,sideface,nface,ipoface,nodface,
		       ne,ipkon,kon,lakon,ikboun,ilboun,xboun,nboun,nk,isolidsurf,
		       nsolidsurf,ifreestream,nfreestream,neighsolidsurf,iponoel,inoel,
		       &inoelfree,nef,co,ipompc,nodempc,ikmpc,ilmpc,nmpc,set,istartset,
		       iendset,ialset,nset,iturbulent));

    RENEW(sideface,char,*nface);
    RENEW(nelemface,ITG,*nface);
    SFREE(ipoface);SFREE(nodface);
    RENEW(inoel,ITG,3*inoelfree);

  }else if(*mcs>1){
    printf(" *ERROR in compfluid: for CFD only one cyclic symmetry\n");
    printf("        conditions is allowed\n");
    FORTRAN(stop,());
  }

  /* determining the matrix structure */

  nzs=1000000;
  
  NNEW(ipointer,ITG,3**nk);
  NNEW(mast1,ITG,nzs);
  NNEW(irowv,ITG,nzs);
  NNEW(irowp,ITG,nzs);
  NNEW(icolv,ITG,3**nk);
  NNEW(icolp,ITG,*nk);
  NNEW(jqv,ITG,3**nk+1);
  NNEW(jqp,ITG,*nk+1);
  NNEW(nactdoh,ITG,mt**nk);

  if(*ithermal>1){
    NNEW(irowt,ITG,nzs);
    NNEW(icolt,ITG,*nk);
    NNEW(jqt,ITG,*nk+1);
  }

  if(*iturbulent!=0){
    NNEW(irowk,ITG,nzs);
    NNEW(icolk,ITG,*nk);
    NNEW(jqk,ITG,*nk+1);
    NNEW(nactdok,ITG,*nk);
  }

  mastructffem(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	       nodempc,nmpc,nactdoh,icolt,icolv,icolp,icolk,jqt,jqv,jqp,
	       jqk,&mast1,&irowt,&irowv,&irowp,&irowk,isolver,&neqt,&neqv,
	       &neqp,&neqk,ikmpc,ilmpc,ipointer,&nzst,&nzsv,&nzsp,&nzsk,
	       ithermal,ikboun,ilboun,iturbulent,nactdok,ifreestream,nfreestream,
	       isolidsurf,nsolidsurf,&nzs,&iexplicit,ielmat,inomat,labmpc);

  SFREE(ipointer);SFREE(mast1);

  /* initialization */

  NNEW(yy,double,*nk);
  NNEW(xsolidsurf,double,*nsolidsurf);
  NNEW(dh,double,*nk);
  NNEW(vcon,double,mt**nk);
  NNEW(vconini,double,mt**nk);
  NNEW(x,double,*nsolidsurf);
  NNEW(y,double,*nsolidsurf);
  NNEW(z,double,*nsolidsurf);
  NNEW(xo,double,*nsolidsurf);
  NNEW(yo,double,*nsolidsurf);
  NNEW(zo,double,*nsolidsurf);
  NNEW(nx,ITG,*nsolidsurf);
  NNEW(ny,ITG,*nsolidsurf);
  NNEW(nz,ITG,*nsolidsurf);
  NNEW(dhel,double,*ne);
  
  FORTRAN(initialcfdfem,(yy,nk,co,ne,ipkon,kon,lakon,x,y,z,xo,yo,zo,nx,ny,nz,
			 isolidsurf,neighsolidsurf,xsolidsurf,dh,nshcon,shcon,
			 nrhcon,rhcon,vold,vcon,ntmat_,iponoel,inoel,
			 &iexplicit,ielmat,nsolidsurf,iturbulent,physcon,
			 &compressible,matname,inomat,vcontu,mi,&euler,
			 ithermal,dhel));
  
  SFREE(x);SFREE(y);SFREE(z);SFREE(xo);SFREE(yo);SFREE(zo);SFREE(nx);SFREE(ny);
  SFREE(nz);SFREE(dhel);

  FORTRAN(normmpc,(nmpc,ipompc,nodempc,coefmpc,inomat));
  
  /* calculating the shape functions, their derivatives and the
     Jacobian determinant in the integration points of the elements */

  nvar_=35**ne;
  NNEW(ipvar,ITG,*ne);
  NNEW(var,double,nvar_);

  nvarf_=8**ne;
  NNEW(ipvarf,ITG,*ne);
  NNEW(varf,double,nvarf_);

  calcshapef(&nvar_,ipvar,&var,ne,lakon,co,ipkon,kon,
             nelemface,sideface,nface,&nvarf_,ipvarf,
             &varf);

  /* composing those left hand sides which do not depend on the increment */

  /* lhs for the energy */

  if(*ithermal>1){
    NNEW(adbt,double,neqt);
    NNEW(aubt,double,nzst);

    FORTRAN(mafilltlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,
			nboun,ipompc,nodempc,coefmpc,nmpc,nactdoh,icolt,jqt,
			irowt,&neqt,&nzlt,ikmpc,ilmpc,ikboun,ilboun,&nzst,
			adbt,aubt,ipvar,var,dhel));

    NNEW(adlt,double,neqt);
    FORTRAN(lump,(adbt,aubt,adlt,irowt,jqt,&neqt));
  
    //spooles_factor_t(adbt,aubt,icolt,irowt,&neqt,&nzst);
  }

  /* lhs for the velocity */

  NNEW(adbv,double,neqv);
  NNEW(aubv,double,nzsv);

  FORTRAN(mafillvlhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
		      xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
		      nactdoh,icolv,jqv,irowv,&neqv,&nzlv,
		      ikmpc,ilmpc,ikboun,ilboun,&nzsv,adbv,aubv,ipvar,var,
		      dhel));

  NNEW(adlv,double,neqv);
  FORTRAN(lump,(adbv,aubv,adlv,irowv,jqv,&neqv));
  
  //spooles_factor_v(adbv,aubv,icolv,irowv,&neqv,&nzsv);

  /* lhs for the pressure  */

  NNEW(adbp,double,neqp);
  NNEW(aubp,double,nzsp);
      
  FORTRAN(mafillplhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
		      xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nactdoh,icolp,jqp,
		      irowp,&neqp,&nzlp,ikmpc,ilmpc,ikboun,ilboun,&nzsp,adbp,
		      aubp,nmethod,&iexplicit,ipvar,var,dhel));

  /* lumping is only applied to compressible fluids */

  if(iexplicit==1){
    NNEW(adlp,double,neqp);
    FORTRAN(lump,(adbp,aubp,adlp,irowp,jqp,&neqp));
    /*spooles_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
      &symmetryflag,&inputformat,&nzsp);*/
  }

  if((iexplicit!=1)&&(neqp>0)){

    /* LU decomposition of the left hand matrix */

    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
		     &symmetryflag,&inputformat,&nzsp);
#else
      printf("*ERROR in compfluid: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,token);
#else
      printf("*ERROR in compfluid: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(adbp,&aubp,adb,aub,&sigma,icolp,&irowp,&neqp,&nzsp);
#else
      printf("*ERROR in compfluid: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
		     &symmetryflag,&inputformat,jqp,&nzsp);
#else
      printf("*ERROR in compfluid: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
      
  }

  /* lhs for the turbulent parameters */

  if(*iturbulent!=0){
    NNEW(adbk,double,neqk);
    NNEW(aubk,double,nzsk);
    FORTRAN(mafillklhs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
			xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
			nactdok,icolk,jqk,irowk,&neqk,&nzlk,
			ikmpc,ilmpc,ikboun,ilboun,&nzsk,adbk,aubk,ipvar,var));

    NNEW(adlk,double,neqk);
    FORTRAN(lump,(adbk,aubk,adlk,irowk,jqk,&neqk));
  }

  /* starting the main loop */

  NNEW(v,double,mt**nk);
  NNEW(vtu,double,2**nk);
      
  /* inserting the velocity and temperature conditions 
     for incompressible materials*/
    
  if(compressible==0){
    FORTRAN(applybounfem,(nodeboun,ndirboun,nboun,xbounact,ithermal,nk,iponoel,
			  inoel,vold,vcontu,t1act,isolidsurf,nsolidsurf,
			  xsolidsurf,nfreestream,ifreestream,iturbulent,vcon,
			  shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
			  &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,
			  inomat,mi,ikboun,ilboun,ilmpc,labmpc));
  }

  /* initialization of voldo: solution from last fluid increment
     and dtimefo: last time step */
  
  NNEW(voldo,double,mt**nk);
  memcpy(voldo,vold,sizeof(double)*mt**nk);
  dtimefo=1.;
  
  /* ttime is the total time up to the start of the present increment
     time is the step time up to the end of the present increment 
     dtime is the present increment size */

  reltime=*time/(*tper);
  //  printf("reltime= %e %e %e\n",*time,*tper,reltime);

  ttimef=*ttime;
  timef=*time-*dtime;
  NNEW(dt,double,*nk);

  if(compressible){
    NNEW(sa,double,neqt);
    NNEW(sav,double,neqv);
    NNEW(del,double,7*jmax[1]);
    shockcoef=physcon[13];
  }

  iit=0;

  do{

    iit++;

    /* determining a new time increment */

    if((*nmethod==4)||((*nmethod==1)&&((iit/ipower)*ipower==iit))){
      FORTRAN(compdt,(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,iponoel,
		      inoel,&dtimef,&iexplicit,ielmat,physcon,dh,cocon,ncocon,
		      ithermal,mi,ipkon,kon,lakon,ne,v,co,iturbulent,vcontu,
		      vcon));
      if(*nmethod==1) ipower*=2;
      printf("iteration %d\n",iit);
    }

    timef+=dtimef;
    if((*time<timef)&&(*nmethod==4)){
      dtimef-=timef-*time;
      timef=*time;
      iconvergence=1;
    }
    reltimef=timef/(*tper);
    if(reltimef>1.) reltimef=1.;

    /* predicting the new velocities by linear extrapolation */

    if(compressible==1){
      ratio=dtimef/(2.*dtimefo);
      for(i=0;i<*nk;i++){
	for(j=1;j<4;j++){
	  temp=vold[i*mt+j];
	  vold[i*mt+j]+=(vold[i*mt+j]-voldo[i*mt+j])*ratio;
	  voldo[i*mt+j]=temp;
	}
      }
      dtimefo=dtimef;
    }

    /* determining the instantaneous load */

    if(*nmethod==1){
      
      /* boundary conditions at end of mechanical increment */

      FORTRAN(temploadfem,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
			   xload,xloadact,iamload,nload,ibody,xbody,nbody,
			   xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			   namta,nam,ampli,time,&reltime,ttime,dtime,ithermal,
			   nmethod,xbounold,xboun,xbounact,iamboun,nboun,
			   nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,co,
			   vold,itg,ntg,amname,ikboun,ilboun,nelemload,
			   sideload,mi,ntrans,trab,inotr,vold,integerglob,
			   doubleglob,tieset,istartset,iendset,ialset,ntie,
			   nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
	  
    }else if(*nmethod==4){

      /* boundary conditions at end of fluid increment */

      FORTRAN(temploadfem,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
			   xload,xloadact,iamload,nload,ibody,xbody,nbody,
			   xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			   namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,
			   ithermal,nmethod,xbounold,xboun,xbounact,iamboun,
			   nboun,nodeboun,ndirboun,nodeforc,ndirforc,istep,
			   iinc,co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,
			   sideload,mi,ntrans,trab,inotr,vold,integerglob,
			   doubleglob,tieset,istartset,iendset,ialset,ntie,
			   nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));
    }

    /*    if((iit/jout[1])*jout[1]==iit){
	  nnstep=6;
	  FORTRAN(frddummy,(co,nk,kon,ipkon,lakon,ne,v,vold,
	  &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon));
	  }*/

    /* STEP 1: velocity correction */

    NNEW(b,double,num_cpus*neqv);

    co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
    nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
    ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
    nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;nforc1=nforc;
    nelemload1=nelemload;sideload1=sideload;xload1=xload;nload1=nload;
    xbody1=xbody;ipobody1=ipobody;nbody1=nbody;nactdoh1=nactdoh;
    icolv1=icolv;jqv1=jqv;irowv1=irowv;neqv1=neqv;nzlv1=nzlv;
    nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
    ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
    ntmat_1=ntmat_;t01=t0;ithermal1=ithermal;vold1=vold;vcon1=vcon;
    nzsv1=nzsv;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
    physcon1=physcon;shcon1=shcon;nshcon1=nshcon;ttime1=ttime;
    timef1=timef;istep1=istep;iinc1=iinc;ibody1=ibody;xloadold1=xloadold;
    iturbulent1=iturbulent;vcontu1=vcontu;yy1=yy;nelemface1=nelemface;
    sideface1=sideface;nface1=nface;compressible1=compressible;
    dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;sti1=sti;
  
    /* create threads and wait */
  
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)mafillv1rhsmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    for(i=0;i<neqv;i++){
      for(j=1;j<num_cpus;j++){
	b[i]+=b[i+j*neqv];
      }
    }
    RENEW(b,double,neqv);SFREE(ithread);

    NNEW(sol,double,neqv);
    NNEW(aux,double,neqv);
    FORTRAN(solveeq,(adbv,aubv,adlv,addiv,b,sol,aux,icolv,irowv,jqv,
		     &neqv,&nzsv,&nzlv));
      
    SFREE(b);SFREE(aux);
      
    //spooles_solve_v(b,&neqv);
      
    /* storing the velocity correction in v */
      
    /*FORTRAN(resultsv1,(nk,nactdoh,v,b,ipompc,nodempc,coefmpc,nmpc,mi));
      SFREE(b);*/
      
    FORTRAN(resultsv1,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
    SFREE(sol);

    /* apply the SPC's and MPC's to Delta V* */

    /*  FORTRAN(applybounv,(nodeboun,ndirboun,nboun,v,nmpc,nodempc,ipompc,
	coefmpc,inomat,mi));*/

    /*  if((iit/jout[1])*jout[1]==iit){
	nnstep=1;
	FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
	&kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
	physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep));
	}*/

    /* STEP 2: pressure correction */

    NNEW(b,double,num_cpus*neqp);

    co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
    nodeboun1=nodeboun;ndirboun1=ndirboun;xbounact1=xbounact;nboun1=nboun;
    ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
    nelemface1=nelemface;sideface1=sideface;nface1=nface;
    nactdoh1=nactdoh;icolp1=icolp;jqp1=jqp;irowp1=irowp;neqp1=neqp;
    nzlp1=nzlp;nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
    ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
    ntmat_1=ntmat_;vold1=vold;vcon1=vcon;nzsp1=nzsp;dtimef1=dtimef;
    matname1=matname;mi1=mi;ncmat_1=ncmat_;shcon1=shcon;nshcon1=nshcon;
    v1=v;theta11=theta1;iexplicit1=iexplicit;physcon1=physcon;
    dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;
  
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)mafillprhsmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    for(i=0;i<neqp;i++){
      for(j=1;j<num_cpus;j++){
	b[i]+=b[i+j*neqp];
      }
    }
    RENEW(b,double,neqp);SFREE(ithread);

    NNEW(sol,double,neqp);
    if((iexplicit==1)&&(neqp>0)){
      NNEW(aux,double,neqp);
      FORTRAN(solveeq,(adbp,aubp,adlp,addiv,b,sol,aux,icolp,irowp,jqp,
		       &neqp,&nzsp,&nzlp));
      SFREE(b);SFREE(aux);
      /*spooles_solve(b,&neqp);

	memcpy(&sol[0],&b[0],sizeof(double)*neqp);
	SFREE(b);*/
    }else if(neqp>0){

      /* solving the system of equations (only for liquids) */

      if(*isolver==0){
#ifdef SPOOLES
	spooles_solve(b,&neqp);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	sgi_solve(b,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_solve(b,&neqp);
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_solve(b,&neqp,&symmetryflag,&inputformat,&nrhs);
#endif
      }

      /* copying the solution into field sol */

      for(i=0;i<neqp;i++){
	sol[i]=b[i]/(theta1*theta2*dtimef*dtimef);
      }
      SFREE(b);

    }

    /* storing the pressure (incompressible) or density (compressible)
       correction in v */

    FORTRAN(resultsp,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,
		      mi));
    SFREE(sol);

    if(iexplicit==0){
      
      /* inserting the pressure boundary conditions for liquids */

      FORTRAN(applybounp,(nodeboun,ndirboun,nboun,xbounact,
			  ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
			  nsolidsurf,xsolidsurf,nfreestream,ifreestream,iturbulent,
			  vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
			  ipompc,nodempc,coefmpc,nmpc,inomat,mi));
    }

    /*      if((iit/jout[1])*jout[1]==iit){
	    nnstep=2;
	    FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
	    &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
	    physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep));
	    }*/
      
    /* STEP 3: velocity correction */

    /*      printf("STEP3: velocity correction\n\n");*/

    NNEW(b,double,num_cpus*neqv);

    co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
    nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
    ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
    nactdoh1=nactdoh;icolv1=icolv;jqv1=jqv;irowv1=irowv;neqv1=neqv;
    nzlv1=nzlv;nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
    ilboun1=ilboun;vold1=vold;nzsv1=nzsv;dtimef1=dtimef;v1=v;
    theta21=theta2;iexplicit1=iexplicit;mi1=mi;
    dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;
  
    /* create threads and wait */
  
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)mafillv2rhsmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    for(i=0;i<neqv;i++){
      for(j=1;j<num_cpus;j++){
	b[i]+=b[i+j*neqv];
      }
    }
    RENEW(b,double,neqv);SFREE(ithread);

    NNEW(sol,double,neqv);
    NNEW(aux,double,neqv);
    FORTRAN(solveeq,(adbv,aubv,adlv,addiv,b,sol,aux,icolv,irowv,jqv,
		     &neqv,&nzsv,&nzlv));
    SFREE(b);SFREE(aux);
      
    //spooles_solve_v(b,&neqv);
      
    /* storing the velocity correction in v */
      
    /*FORTRAN(resultsv2,(nk,nactdoh,v,b,ipompc,nodempc,coefmpc,nmpc,mi));
      SFREE(b);*/

    FORTRAN(resultsv2,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
    SFREE(sol);

    /*  if((iit/jout[1])*jout[1]==iit){
	nnstep=3;
	FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
	&kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
	physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep));
	}*/

    /* STEP 4: energy correction */

    if(*ithermal>1){


      NNEW(b,double,num_cpus*neqt);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
      nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;nforc1=nforc;
      nelemload1=nelemload;sideload1=sideload;xload1=xload;nload1=nload;
      xbody1=xbody;ipobody1=ipobody;nbody1=nbody;nactdoh1=nactdoh;
      neqt1=neqt;
      nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
      ntmat_1=ntmat_;t01=t0;ithermal1=ithermal;vold1=vold;vcon1=vcon;
      nzst1=nzst;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
      physcon1=physcon;shcon1=shcon;nshcon1=nshcon;ttime1=ttime;
      timef1=timef;istep1=istep;iinc1=iinc;ibody1=ibody;xloadold1=xloadold;
      reltimef1=reltimef;cocon1=cocon;ncocon1=ncocon;nelemface1=nelemface;
      sideface1=sideface;nface1=nface;compressible1=compressible;
      vcontu1=vcontu;yy1=yy;iturbulent1=iturbulent;
      dt1=dt;ipvar1=ipvar;var1=var;ipvarf1=ipvarf;varf1=varf;
	  
      NNEW(ithread,ITG,num_cpus);
      for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafilltrhsmt, (void *)&ithread[i]);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	  
      for(i=0;i<neqt;i++){
	for(j=1;j<num_cpus;j++){
	  b[i]+=b[i+j*neqt];
	}
      }
      RENEW(b,double,neqt);SFREE(ithread);
	  
      NNEW(sol,double,neqt);
      NNEW(aux,double,neqt);
      FORTRAN(solveeq,(adbt,aubt,adlt,addiv,b,sol,aux,icolt,irowt,jqt,
		       &neqt,&nzst,&nzlt));
      SFREE(b);SFREE(aux);
	  
      //spooles_solve_t(b,&neqt);
	  
      /* storing the temperature correction in v */
	  
      /*FORTRAN(resultst,(nk,nactdoh,v,b,ipompc,nodempc,coefmpc,nmpc,mi));
	SFREE(b);*/
	  
      FORTRAN(resultst,(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,mi));
      SFREE(sol);
    }

    /*   if((iit/jout[1])*jout[1]==iit){
	 nnstep=4;
	 FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
	 &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
	 physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep));
	 }*/

    /* STEP 5: turbulent correction */

    if(*iturbulent!=0){
      /*	        if((iit/jout[1])*jout[1]==iit){
			nnstep=6;
			FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
			&kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
			physcon,filab,inomat,ntrans,inotr,trab,mi,stn));
			}*/

      NNEW(bk,double,num_cpus*neqk);
      NNEW(bt,double,num_cpus*neqk);
      NNEW(ck,double,num_cpus*neqk);
      NNEW(ct,double,num_cpus*neqk);

      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
      ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
      nelemface1=nelemface;sideface1=sideface;nface1=nface;
      nactdok1=nactdok;neqk1=neqk;
      nmethod1=nmethod;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
      ilboun1=ilboun;rhcon1=rhcon;nrhcon1=nrhcon;ielmat1=ielmat;
      ntmat_1=ntmat_;vold1=vold;vcon1=vcon;
      nzsk1=nzsk;dtimef1=dtimef;matname1=matname;mi1=mi;ncmat_1=ncmat_;
      shcon1=shcon;nshcon1=nshcon;theta11=theta1;
      vcontu1=vcontu;isolidsurf1=isolidsurf;nsolidsurf1=nsolidsurf;
      ifreestream1=ifreestream;nfreestream1=nfreestream;
      xsolidsurf1=xsolidsurf;yy1=yy;compressible1=compressible;
      iturbulent1=iturbulent;ithermal1=ithermal;ipvar1=ipvar;var1=var;
      ipvarf1=ipvarf;varf1=varf;dt1=dt;
	  
      NNEW(ithread,ITG,num_cpus);
      for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillkrhsmt, (void *)&ithread[i]);
      }
      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	  
      for(i=0;i<neqk;i++){
	for(j=1;j<num_cpus;j++){
	  bk[i]+=bk[i+j*neqk];
	  bt[i]+=bt[i+j*neqk];
	  ck[i]+=ck[i+j*neqk];
	  ct[i]+=ct[i+j*neqk];
	}
      }
      RENEW(bk,double,neqk);
      RENEW(bt,double,neqk);
      RENEW(ck,double,neqk);
      RENEW(ct,double,neqk);
      SFREE(ithread);
	  
      NNEW(adlkk,double,neqk);
      //	  for(i=0;i<neqk;i++)adlkk[i]=1./(1./adlk[i]+ck[i]);
      SFREE(ck);

      NNEW(solk,double,neqk);
      NNEW(aux,double,neqk);
      //	  FORTRAN(solveeq,(adbv,aubv,adlkk,addiv,bk,solk,aux,icolk,irowk,jqk,
      //			   &neqk,&nzsk,&nzlk));
      FORTRAN(solveeq,(adbv,aubv,adlk,addiv,bk,solk,aux,icolk,irowk,jqk,
		       &neqk,&nzsk,&nzlk));
      SFREE(bk);SFREE(aux);SFREE(adlkk);
	  
      NNEW(adlkt,double,neqk);
      //	  for(i=0;i<neqk;i++)adlkt[i]=1./(1./adlk[i]+ct[i]);
      SFREE(ct);
	  
      NNEW(solt,double,neqk);
      NNEW(aux,double,neqk);
      //	  FORTRAN(solveeq,(adbv,aubv,adlkt,addiv,bt,solt,aux,icolk,irowk,jqk,
      //			   &neqk,&nzsk,&nzlk));
      FORTRAN(solveeq,(adbv,aubv,adlk,addiv,bt,solt,aux,icolk,irowk,jqk,
		       &neqk,&nzsk,&nzlk));
      SFREE(bt);SFREE(aux);SFREE(adlkt);
	  
      /* storing the turbulence correction in vtu */
	  
      FORTRAN(resultsk,(nk,nactdok,vtu,solk,solt,ipompc,nodempc,
			coefmpc,nmpc));
      SFREE(solk);SFREE(solt);

      /*	         if((iit/jout[1])*jout[1]==iit){
			 nnstep=5;
			 FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
			 &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
			 physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep));
			 }*/
    }
      
    /* update the conservative variables
       (for incompressible fluids: pressure instead of density  */

    FORTRAN(updatecon,(vold,vcon,v,nk,
		       ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
		       nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
		       nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,
		       iturbulent,
		       inomat,nodeboun,ndirboun,nboun,mi,co,&factor));

    /* inserting the boundary conditions for the turbulence
       parameters */

    /*   if(*iturbulent!=0){
	 FORTRAN(applybounk,(nodeboun,ndirboun,nboun,xbounact,
	 iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
	 ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
	 inoel,physcon,&compressible,ielmat,nshcon,shcon,nrhcon,
	 rhcon,vcontu,ntmat_,labmpc,inomat));
	 }*/
 
    /* extrapolating the stresses (for debugging purposes) */

    /*	      nfield=6;
	      ndim=6;
	      cfd=1;
	      if((*norien>0)&&(strcmp1(&filab[179],"L")==0)){
	      iorienglob=1;
	      }else{
	      iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[178],1);
	      NNEW(stn,double,6**nk);
	      NNEW(inum,ITG,*nk);
	      FORTRAN(extrapolatefem,(sti,stn,ipkon,inum,kon,lakon,
	      &nfield,nk,ne,mi,&ndim,orab,ielorien,co,&iorienglob,
	      cflag,nelemload,nload,nodeboun,nboun,ndirboun,
	      vold,ithermal,&force,&cfd));*/
     
    /* smoothing the solution (only for compressible fluids) */

    if((compressible)&&(shockcoef>0.)){

      /*    ismooth=1;*/

      /* shocksmoothing rho * total energy density */

      NNEW(sol,double,neqt);
      NNEW(aux,double,neqt);
      for(i=0;i<neqt;i++){sol[i]=vcon[mt*i];}
      FORTRAN(smoothshock,(adbt,aubt,adlt,addiv,sol,aux,icolt,irowt,jqt,
			   &neqt,&nzlt,sa));
      for(i=0;i<neqt;i++){vcon[mt*i]=sol[i];}
      SFREE(sol);SFREE(aux);

      /* shocksmoothing rho * velocity */

      NNEW(sol,double,neqv);
      NNEW(aux,double,neqv);
      for(i=0;i<neqv/3;i++){
	for(j=0;j<3;j++){
	  sol[3*i+j]=vcon[mt*i+j+1];
	}
      }
      FORTRAN(smoothshock,(adbv,aubv,adlv,addiv,sol,aux,icolv,irowv,jqv,
			   &neqv,&nzlv,sav));
      for(i=0;i<neqv/3;i++){
	for(j=0;j<3;j++){
	  vcon[mt*i+j+1]=sol[3*i+j];
	}
      }
      SFREE(sol);SFREE(aux);

      /* shocksmoothing rho */

      NNEW(sol,double,neqp);
      NNEW(aux,double,neqp);
      for(i=0;i<neqp;i++){sol[i]=vcon[mt*i+4];}
      FORTRAN(smoothshock,(adbp,aubp,adlp,addiv,sol,aux,icolp,irowp,jqp,
			   &neqp,&nzlp,sa));
      for(i=0;i<neqp;i++){vcon[mt*i+4]=sol[i];}
      SFREE(sol);SFREE(aux);
	  
    }

    /* applying the MPC's to the conservative variables */
    
    /*FORTRAN(applympcfem,(nodeboun,ndirboun,nboun,xbounact,ithermal,nk,iponoel,
			 inoel,vold,vcontu,t1act,isolidsurf,nsolidsurf,
			 xsolidsurf,nfreestream,ifreestream,iturbulent,vcon,
			 shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
			 &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,
			 inomat,mi,ikboun,ilboun,ilmpc,labmpc));*/

    /* deriving the physical variables from the conservative
       variables */

    FORTRAN(con2phys,(vold,vcon,v,nk,
		      ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
		      nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
		      nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,iturbulent,
		      inomat,nodeboun,ndirboun,nboun,mi,co,&factor));
      
    /* inserting SPC's to the physical variables */
      
    FORTRAN(applybounfem,(nodeboun,ndirboun,nboun,xbounact,ithermal,nk,iponoel,
			  inoel,vold,vcontu,t1act,isolidsurf,nsolidsurf,
			  xsolidsurf,nfreestream,ifreestream,iturbulent,vcon,
			  shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
			  &compressible,&ismooth,nmpc,nodempc,ipompc,coefmpc,
			  inomat,mi,ikboun,ilboun,ilmpc,labmpc));

    /* calculating the pressure gradient for the shock smoothing in
       the next iteration */

    if((compressible)&&(shockcoef>0.)){
      FORTRAN(presgradient,(iponoel,inoel,sa,sav,&neqt,dt,&shockcoef,
			    &dtimef,ipkon,kon,lakon,vold,mi,
			    &compressible,nmethod,dt,isolidsurf,
			    nsolidsurf,co,&euler));
    }
      
    /* inserting the boundary conditions for the turbulence
       parameters */
      
    /*    if(*iturbulent!=0){
	  FORTRAN(applybounk,(nodeboun,ndirboun,nboun,xbounact,
	  iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
	  ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
	  inoel,physcon,&compressible,ielmat,nshcon,shcon,nrhcon,
	  rhcon,vcontu,ntmat_,labmpc,inomat,mi,ithermal));
	  }*/

    /* check iconvergence */

    FORTRAN(cfdconv,(vold,vcon,v,nk,
		     ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
		     nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
		     nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,iturbulent,
		     inomat,nodeboun,ndirboun,nboun,mi,co,&factor,
		     vconini,&dtimef,del,&sum,&sumx,&sumxx,sumy,sumxy,&nstart,
		     &shockcoef));

    if(compressible==1) iconvergence=0;
      
    if(((iit/jout[1])*jout[1]==iit)||(iconvergence==1)||
       (iit==jmax[1])){

      /* calculating the stress and the heat flow at the
	 integration points, if requested */

      if(strcmp1(&filab[3306],"SF  ")==0)isti=1;
      if(strcmp1(&filab[3393],"HFLF")==0)iqfx=1;
      for(i=0;i<*nprint;i++){
	if(strcmp1(&prlab[6*i],"SF")==0) isti=1;
	if(strcmp1(&prlab[6*i],"HFLF")==0)iqfx=1;
      }
      if((isti==1)||(iqfx==1)){
	FORTRAN(calcstressheatfluxfem,(kon,lakon,ipkon,ielmat,ntmat_,
				       vold,matname,mi,shcon,nshcon,iturbulent,&compressible,
				       ipvar,var,sti,qfx,cocon,ncocon,ne,&isti,&iqfx));
      }
 
      /* extrapolating the stresses */

      if(strcmp1(&filab[3306],"SF  ")==0){
	nfield=6;
	ndim=6;
	if((*norien>0)&&(strcmp1(&filab[2962],"L")==0)){
	  iorienglob=1;
	}else{
	  iorienglob=0;
	}
	strcpy1(&cflag[0],&filab[2962],1);
	NNEW(stn,double,6**nk);
	NNEW(inum,ITG,*nk);
	FORTRAN(extrapolatefem,(sti,stn,ipkon,inum,kon,lakon,
				&nfield,nk,ne,mi,&ndim,orab,ielorien,co,&iorienglob,
				cflag,nelemload,nload,nodeboun,nboun,ndirboun,
				vold,ithermal,&force,&icfdout,ielmat,thicke,filab));
	SFREE(inum);
      }

      /* extrapolating the heat flow */

	  
      if(strcmp1(&filab[3393],"HFLF")==0){
	nfield=3;
	ndim=3;
	if((*norien>0)&&(strcmp1(&filab[3049],"L")==0)){
	  iorienglob=1;
	}else{
	  iorienglob=0;
	}
	strcpy1(&cflag[0],&filab[3049],1);
	NNEW(qfn,double,3**nk);
	NNEW(inum,ITG,*nk);
	FORTRAN(extrapolatefem,(qfx,qfn,ipkon,inum,kon,lakon,
				&nfield,nk,ne,mi,&ndim,orab,ielorien,co,
				&iorienglob,cflag,nelemload,nload,nodeboun,
				nboun,ndirboun,vold,ithermal,&force,&icfdout,
				ielmat,thicke,filab));
	SFREE(inum);
      }
	 
      /* check whether the Mach number is requested */

      if((strcmp1(&filab[1914],"MACH")==0)|| 
	 (strcmp1(&filab[3219],"TTF")==0)||
	 (strcmp1(&filab[3132],"PTF")==0)){
	FORTRAN(calcmach,(vold,vcon,v,nk,
			  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,&iout,
			  nmethod,&iconvergence,physcon,iponoel,inoel,ithermal,
			  nactdoh,&iit,&compressible,&ismooth,vcontu,vtu,
			  iturbulent,inomat,nodeboun,ndirboun,nboun,mi,co,
			  &factor));
      }

      /* print output */

      if(iconvergence==1) timef=*time;
      FORTRAN(printoutfluidfem,(set,nset,istartset,iendset,ialset,nprint,
				prlab,prset,vold,t1,fn,ipkon,lakon,sti,eme,
				xstate,ener,mi,&nstate_,ithermal,co,kon,qfx,
				&timef,trab,inotr,ntrans,orab,ielorien,norien,
				nk,ne,inum,filab,vold,&ikin,ielmat,thicke,eme,
				vcontu,physcon));

      /* lift and drag force */

      FORTRAN(printoutfacefem,(co,rhcon,nrhcon,ntmat_,vold,shcon,nshcon,
	  cocon,ncocon,&compressible,istartset,iendset,ipkon,lakon,kon,
	  ialset,prset,&timef,nset,set,nprint,prlab,ielmat,mi));

      /* frd output */

      nnstep=6;
      FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
			   &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
			   physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,
			   istep));

      if(strcmp1(&filab[3306],"SF  ")==0){SFREE(stn);}
      if(strcmp1(&filab[3393],"HFLF")==0){SFREE(qfn);}

    }
      
    /*  if(((iit/jout[1])*jout[1]==iit)||(iconvergence==1)){
	FORTRAN(printoutfluidfem,(set,nset,istartset,iendset,ialset,nprint,
	prlab,prset,vold,t1,fn,ipkon,lakon,sti,eme,xstate,ener,
	mi,&nstate_,ithermal,co,kon,qfx,&timef,trab,inotr,ntrans,
	orab,ielorien,norien,nk,ne,inum,filab,vold,&ikin,ielmat,thicke,
	eme,vcontu,physcon));*/

    /* lift and drag force */
	  
    /*	  FORTRAN(printoutface,(co,rhcon,nrhcon,ntmat_,vold,shcon,nshcon,
	  cocon,ncocon,&compressible,istartset,iendset,ipkon,lakon,kon,
	  ialset,prset,&timef,nset,set,nprint,prlab,ielmat,mi));
	  }*/
      
    /*   if(iconvergence==1){
	 nnstep=6;
	 FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
	 &kode,&timef,ielmat,matname,&nnstep,vtu,vcontu,vcon,
	 physcon,filab,inomat,ntrans,inotr,trab,mi,stn));
	 break;
	 }*/

    /*      if(iconvergence==1){
	    if(compressible==1){
	    if(shockcoef>0.){
	    shockcoef=shockcoef/2.;
	    if(shockcoef>1.e-2){
	    printf("new chockcoefficient = %e\n",shockcoef);
	    iconvergence=0;
	    sum=0.;sumx=0.;sumxx=0.;iit=500;
	    for(i=0;i<7;i++){sumy[i]=0.;sumxy[i]=0.;}
	    }
	    }
	    }
	    }*/
    if((iit==jmax[1])||(iconvergence==1)) break;
      
    ttimef+=dtimef;
  }while(1);
  
  if((iexplicit!=1)&&(neqp>0)){
    if(*isolver==0){
#ifdef SPOOLES
      spooles_cleanup();
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      sgi_cleanup(token);
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_cleanup();
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_cleanup(&neqp,&symmetryflag,&inputformat);
#endif
    }
  }

  if(compressible){SFREE(sa);SFREE(sav);SFREE(del);}

  SFREE(yy);SFREE(xsolidsurf);SFREE(dt);SFREE(dh);SFREE(vcon);SFREE(vconini);

  SFREE(irowv);SFREE(irowp);
  SFREE(icolv);SFREE(icolp);
  SFREE(jqv);SFREE(jqp);
  SFREE(nactdoh);

  SFREE(adbv);SFREE(adbp);
  SFREE(aubv);SFREE(aubp);
  SFREE(adlv);if(iexplicit==1) SFREE(adlp);

  if(*ithermal>1){
    SFREE(irowt);SFREE(icolt);SFREE(jqt);SFREE(adbt);SFREE(aubt);SFREE(adlt);
  }

  if(*iturbulent!=0){
    SFREE(irowk);SFREE(icolk);SFREE(jqk);SFREE(nactdok);
    SFREE(adbk);SFREE(aubk);SFREE(adlk);
  }

  SFREE(v);SFREE(vtu);SFREE(var);SFREE(ipvar);SFREE(varf);SFREE(ipvarf);
  SFREE(voldo);

  if(*mcs==1){
    *memmpc_=memmpcref_;*mpcfree=mpcfreeref;*nmpc=nmpcref;
    *nk=nkref;*ne=neref;*nef=nefref;

    RENEW(nodempc,ITG,3*memmpcref_);
    for(k=0;k<3*memmpcref_;k++){nodempc[k]=nodempcref[k];}
    RENEW(coefmpc,double,memmpcref_);
    for(k=0;k<memmpcref_;k++){coefmpc[k]=coefmpcref[k];}
    
    RENEW(ipompc,ITG,*nmpc);
    RENEW(ikmpc,ITG,*nmpc);
    RENEW(ilmpc,ITG,*nmpc);
    RENEW(labmpc,char,20**nmpc);
    
    RENEW(co,double,3**nk);
    RENEW(vold,double,mt**nk);
    RENEW(ipkon,ITG,*ne);
    RENEW(lakon,char,8**ne);
    RENEW(kon,ITG,*nkon);
    RENEW(ielmat,ITG,mi[2]**ne);
    RENEW(inomat,ITG,*nk);

    NNEW(sideface,char,6**nef);
    NNEW(nelemface,ITG,6**nef);
    NNEW(iponoel,ITG,*nk);
    NNEW(inoel,ITG,3*20**nef);
  }

  *ithermal=ithermalref;

  *nodempcp=nodempc;*ipompcp=ipompc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*fmpcp=fmpc;*cop=co;
  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielmatp=ielmat;
  *nelemfacep=nelemface;*sidefacep=sideface;*inoelp=inoel;
  *iponoelp=iponoel;*voldp=vold;*inomatp=inomat;
  
  FORTRAN(stop,());
  FORTRAN(closefilefluid,());

  return;
  
} 

/* subroutine for multithreading of mafillv1rhs */

void *mafillv1rhsmt(ITG *i){

  ITG index,nea,neb,nedelta;

  index=*i*neqv1;
    
  nedelta=(ITG)ceil(*ne1/(double)num_cpus);
  nea=*i*nedelta+1;
  neb=(*i+1)*nedelta;
  if(neb>*ne1) neb=*ne1;

  FORTRAN(mafillv1rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
		       xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nodeforc1,
		       ndirforc1,xforc1,
		       nforc1,nelemload1,sideload1,xload1,nload1,xbody1,
		       ipobody1,nbody1,
		       &b[index],nactdoh1,icolv1,jqv1,irowv1,&neqv1,&nzlv1,
		       nmethod1,ikmpc1,ilmpc1,ikboun1,
		       ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,t01,ithermal1,
		       vold1,vcon1,
		       dt1,matname1,mi1,ncmat_1,physcon1,shcon1,nshcon1,ttime1,
		       &timef1,
		       istep1,iinc1,ibody1,xloadold1,iturbulent1,vcontu1,yy1,
		       nelemface1,sideface1,nface1,&compressible1,&nea,&neb,
		       &dtimef1,
		       ipvar1,var1,ipvarf1,varf1,sti1));

  return NULL;
}

/* subroutine for multithreading of mafillprhs */

void *mafillprhsmt(ITG *i){

  ITG index,nea,neb,nedelta;

  index=*i*neqp1;
    
  nedelta=(ITG)ceil(*ne1/(double)num_cpus);
  nea=*i*nedelta+1;
  neb=(*i+1)*nedelta;
  if(neb>*ne1) neb=*ne1;

  FORTRAN(mafillprhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
		      xbounact1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,
		      nelemface1,sideface1,
		      nface1,&b[index],nactdoh1,icolp1,jqp1,irowp1,&neqp1,
		      &nzlp1,nmethod1,ikmpc1,ilmpc1,
		      ikboun1,ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,vold1,
		      vcon1,&nzsp1,
		      dt1,matname1,mi1,ncmat_1,shcon1,nshcon1,v1,&theta11,
		      &iexplicit1,physcon1,&nea,&neb,&dtimef1,ipvar1,var1,
		      ipvarf1,varf1));

  return NULL;
}

/* subroutine for multithreading of mafillv2rhs */

void *mafillv2rhsmt(ITG *i){

  ITG index,nea,neb,nedelta;

  index=*i*neqv1;
    
  nedelta=(ITG)ceil(*ne1/(double)num_cpus);
  nea=*i*nedelta+1;
  neb=(*i+1)*nedelta;
  if(neb>*ne1) neb=*ne1;

  FORTRAN(mafillv2rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
		       xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,
		       &b[index],nactdoh1,icolv1,jqv1,irowv1,&neqv1,&nzlv1,nmethod1,ikmpc1,ilmpc1,ikboun1,
		       ilboun1,vold1,&nzsv1,dt1,v1,&theta21,&iexplicit1,&nea,&neb,mi1,&dtimef1,
		       ipvar1,var1,ipvarf1,varf1));

  return NULL;
}

/* subroutine for multithreading of mafilltrhs */

void *mafilltrhsmt(ITG *i){

  ITG index,nea,neb,nedelta;

  index=*i*neqt1;
    
  nedelta=(ITG)ceil(*ne1/(double)num_cpus);
  nea=*i*nedelta+1;
  neb=(*i+1)*nedelta;
  if(neb>*ne1) neb=*ne1;
	  
  FORTRAN(mafilltrhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
		      xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nodeforc1,
		      ndirforc1,
		      xforc1,
		      nforc1,nelemload1,sideload1,xload1,nload1,xbody1,
		      ipobody1,nbody1,
		      &b[index],nactdoh1,&neqt1,nmethod1,ikmpc1,ilmpc1,ikboun1,
		      ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,t01,ithermal1,
		      vold1,
		      vcon1,&nzst1,
		      dt1,matname1,mi1,ncmat_1,physcon1,shcon1,nshcon1,ttime1,
		      &timef1,
		      istep1,iinc1,ibody1,xloadold1,&reltimef1,cocon1,ncocon1,
		      nelemface1,
		      sideface1,nface1,&compressible1,vcontu1,yy1,iturbulent1,
		      &nea,
		      &neb,&dtimef1,ipvar1,var1,ipvarf1,varf1));

  return NULL;
}

/* subroutine for multithreading of mafillkrhs */

void *mafillkrhsmt(ITG *i){

  ITG index,nea,neb,nedelta;

  index=*i*neqk1;
    
  nedelta=(ITG)ceil(*ne1/(double)num_cpus);
  nea=*i*nedelta+1;
  neb=(*i+1)*nedelta;
  if(neb>*ne1) neb=*ne1;
	  
  FORTRAN(mafillkrhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
		      xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,nelemface1,
		      sideface1,
		      nface1,nactdok1,&neqk1,nmethod1,ikmpc1,ilmpc1,
		      ikboun1,ilboun1,rhcon1,nrhcon1,ielmat1,ntmat_1,vold1,
		      vcon1,
		      &nzsk1,
		      &dtimef1,matname1,mi1,ncmat_1,shcon1,nshcon1,&theta11,
		      &bk[index],&bt[index],vcontu1,isolidsurf1,nsolidsurf1,
		      ifreestream1,nfreestream1,
		      xsolidsurf1,yy1,&compressible1,iturbulent1,ithermal1,
		      ipvar1,var1,
		      ipvarf1,varf1,&nea,&neb,dt1,&ck[index],&ct[index]));

  return NULL;
}
